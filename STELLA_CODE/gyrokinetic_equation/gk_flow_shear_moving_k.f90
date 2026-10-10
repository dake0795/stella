!###############################################################################
!############## PERPENDICULAR FLOW SHEAR WITH A MOVING WAVEVECTOR ##############
!###############################################################################
!
! In the presence of a perpendicular flow shear, gamma_E = g_exb*g_exbfac, the
! radial wavenumber of a mode with fixed ky changes continuously in time,
!     kx(t) = kx(0) - gamma_E ky t.
! The default (Hammett) treatment keeps the kx grid fixed and remaps every mode
! by one radial grid spacing, dkx = akx(2), once the accumulated shift exceeds
! half a grid spacing, i.e. after a time T_ky = |dkx / (gamma_E ky)|.
!
! With <flow_shear_moving_k> = .true., akx(ikx) is instead treated as a fixed
! label and every term of the gyrokinetic equation and of the field equations
! is evaluated at the physical radial wavenumber
!     kx_phys(iky, ikx) = akx(ikx) - gamma_E aky(iky) tau(iky),
! where tau is the residual shear time. The committed value of tau, at the end
! of each accepted time step, is stored in <shift_state>, and |tau| <= T_ky/2.
! At the Runge-Kutta stage with abscissa c, tau = shift_state + c*code_dt is
! stored in <shift_stage>. After each accepted step, rows for which tau has
! moved by a full grid spacing are relabelled: the distribution function is
! moved from label kx to label kx -/+ dkx, which is the same physical mode,
! so the relabelling does not change the solution. Modes that leave the kx grid
! are dropped. kperp2, the Bessel functions and the field denominators are
! recalculated whenever the physical wavenumbers change.
!
! In the ExB nonlinearity the fields and the distribution function are
! transformed to real space in the sheared frame, y' = y - gamma_E tau x, of
! the Runge-Kutta stage, with the derivatives evaluated at kx_phys.
! 
! The algorithm is restricted to local flux tubes with kinetic electrons, phi
! and A_parallel (no B_parallel), the explicit rk4 scheme and no collisions.
! The mirror term and the parallel streaming can be treated implicitly, as a
! split after the explicit step at the committed metric. The response matrix
! of the implicit streaming is calculated at the label wavenumbers, and the
! implicit step is iterated to the solution at the physical wavenumbers (defect
! correction, see gk_implicit_terms). Optionally, additional response matrices
! at sub-cell phases of the low-ky rows are used as better preconditioners.
! All restrictions are checked at initialisation.
!
!###############################################################################
module gk_flow_shear_moving_k

   use response_matrix, only: response_set_type

   implicit none

   ! Make routines available to other modules
   public :: init_flow_shear_moving_k
   public :: prepare_moving_k_stage
   public :: relabel_moving_k
   public :: prepare_moving_k_stream_solve
   public :: rebuild_response_on_label_metric

   private

   ! Only initialise once
   logical :: initialised_moving_k = .false.

   ! Field denominators at the label wavenumbers (phi: denominator_fields, apar: kperp2),
   ! which the response matrix of the implicit streaming is calculated with
   real, dimension(:, :, :), allocatable :: label_denominator, label_kperp2

   ! Additional streaming response matrices at the sub-cell phases tau = m T/M,
   ! m = -M/2, ..., M/2, for the rows iky = 2, ..., set_iky_max (m = 0 is the
   ! label set), with their field denominators; set_installed(iky) is the set in use
   logical :: sets_built = .false.
   integer :: nsets = 0, set_iky_max = 1
   integer, dimension(:), allocatable :: set_installed
   real, dimension(:, :, :, :), allocatable :: set_denominator, set_kperp2
   type(response_set_type) :: label_set
   type(response_set_type), dimension(:), allocatable :: offset_set

contains

   !****************************************************************************
   !                              Initialisation
   !****************************************************************************
   ! Check the input and evaluate the metric at the committed phase, which is
   ! nonzero after a restart, before the fields are solved for the first time.
   !****************************************************************************
   subroutine init_flow_shear_moving_k

      use mp, only: mp_abort
      use arrays, only: shift_state, denominator_fields, apar_denom, kperp2
      use gk_flow_shear, only: moving_k_stream_defect_tol, moving_k_stream_defect_maxit
      use gk_flow_shear, only: moving_k_stream_defect_method
      use gk_implicit_terms, only: stream_defect_tol, stream_defect_maxit, stream_defect_method

      implicit none

      real, dimension(:, :, :), allocatable :: initial_denominator, initial_apar_denom
      real :: error_denominator, error_apar_denom

      !-------------------------------------------------------------------------

      if (initialised_moving_k) return

      ! Field denominators of the original initialisation (label metric)
      allocate (initial_denominator, source=denominator_fields)
      allocate (initial_apar_denom, source=apar_denom)
      allocate (label_denominator, source=denominator_fields)
      allocate (label_kperp2, source=kperp2(:, :, 1, :))

      call check_moving_k_configuration

      ! Defect correction of the implicit streaming
      select case (trim(moving_k_stream_defect_method))
      case ('simple', 'gmres', 'auto')
         stream_defect_method = moving_k_stream_defect_method
      case default
         call mp_abort('moving_k_stream_defect_method should be simple, gmres or auto')
      end select
      stream_defect_tol = moving_k_stream_defect_tol
      stream_defect_maxit = moving_k_stream_defect_maxit
      initialised_moving_k = .true.
      call prepare_moving_k_stage(0.0)

      ! At zero phase the recalculated denominators must equal the original ones
      if (all(abs(shift_state) <= epsilon(0.))) then
         error_denominator = maxval(abs(denominator_fields - initial_denominator)) &
            / max(maxval(abs(initial_denominator)), tiny(0.0))
         error_apar_denom = maxval(abs(apar_denom - initial_apar_denom)) / max(maxval(abs(initial_apar_denom)), tiny(0.0))
         if (max(error_denominator, error_apar_denom) > 1.e-12) &
            call mp_abort('moving-k flow shear: recalculated field denominators differ from the initial ones at zero phase')
      end if

      deallocate (initial_denominator, initial_apar_denom)

   end subroutine init_flow_shear_moving_k

   !****************************************************************************
   !                        Check the input parameters
   !****************************************************************************
   subroutine check_moving_k_configuration

      use mp, only: proc0, mp_abort
      use arrays, only: shift_state
      use geometry, only: q_as_x
      use grids_kxky, only: box, nakx, naky, akx, aky, nalpha, nx
      use grids_z, only: ntubes
      use grids_time, only: code_time, code_dt_max
      use grids_species, only: spec, has_electron_species
      use parallelisation_layouts, only: fields_kxkyz
      use parameters_physics, only: full_flux_surface, radial_variation
      use parameters_physics, only: include_apar, include_bpar, fphi
      use parameters_physics, only: include_nonlinear, include_parallel_nonlinearity
      use parameters_physics, only: include_mirror, include_parallel_streaming
      use parameters_numerical, only: fully_explicit, fully_implicit, flip_flop
      use parameters_numerical, only: stream_implicit, mirror_implicit, drifts_implicit, driftkinetic_implicit
      use parameters_numerical, only: explicit_algorithm_switch, explicit_algorithm_rk4, nitt
      use dissipation_and_collisions, only: include_collisions, hyper_dissipation
      use dissipation_hyper, only: D_hyper
      use gk_sources, only: source_option_switch, source_option_none, include_qn_source
      use parameters_multibox, only: include_multibox_krook
      use neoclassical_terms, only: include_neoclassical_terms
      use file_utils, only: runtype_option_switch, runtype_multibox
      use initialise_distribution_function, only: scale_to_phiinit, init_distribution_switch
      use namelist_initialise_distribution_function, only: init_distribution_option_restart
      use gk_flow_shear, only: hammett_flow_shear, g_exb, g_exbfac, omprimfac
      use gk_flow_shear, only: shift_times, flow_shear_explicit_remap, explicit_junction_sync

      implicit none

      integer :: iky

      !-------------------------------------------------------------------------

      ! Grid and geometry
      if (.not. box .or. nakx < 3 .or. nalpha /= 1 .or. ntubes /= 1 .or. q_as_x) &
         call mp_abort('moving-k flow shear requires a box grid with nakx >= 3, nalpha = ntubes = 1 and q_as_x = F')
      if (full_flux_surface .or. radial_variation .or. runtype_option_switch == runtype_multibox) &
         call mp_abort('moving-k flow shear does not support full flux surface, radial variation or multibox')
      if (akx(2) <= 0.0 .or. any(aky < 0.0)) call mp_abort('moving-k flow shear requires akx(2) > 0 and aky >= 0')

      ! Fields and species: the field denominators are only recalculated for phi and A_parallel
      if (.not. include_apar .or. include_bpar .or. abs(fphi - 1.0) > epsilon(0.)) &
         call mp_abort('moving-k flow shear requires phi and apar (fphi = 1) without bpar')
      if (.not. has_electron_species(spec)) call mp_abort('moving-k flow shear requires kinetic electrons')
      if (fields_kxkyz) call mp_abort('moving-k flow shear requires fields_kxkyz = F')

      ! Time advance: explicit rk4, optionally with the implicit mirror and streaming split
      if (.not. fully_explicit .or. fully_implicit .or. flip_flop &
          .or. (stream_implicit .and. .not. include_parallel_streaming) &
          .or. (mirror_implicit .and. .not. include_mirror) &
          .or. drifts_implicit .or. driftkinetic_implicit .or. explicit_algorithm_switch /= explicit_algorithm_rk4) &
         call mp_abort('moving-k flow shear requires the explicit rk4 scheme; ' // &
            'only the mirror and parallel streaming can be implicit')
      if (trim(explicit_junction_sync) == 'none') &
         call mp_abort('moving-k flow shear requires explicit_junction_sync = step or stage')
      if (flow_shear_explicit_remap) call mp_abort('moving-k flow shear can not be combined with flow_shear_explicit_remap')
      if (stream_implicit .and. nitt /= 1) call mp_abort('moving-k flow shear with implicit streaming requires nitt = 1')

      ! Physics
      if (include_parallel_nonlinearity) call mp_abort('moving-k flow shear does not support the parallel nonlinearity')
      ! With mod(nx, 3) = 1 products of the extreme sheared labels alias into the retained modes
      if (include_nonlinear .and. mod(nx, 3) == 1) &
         call mp_abort('moving-k flow shear with the nonlinearity requires mod(nx, 3) /= 1')
      if (include_collisions) call mp_abort('moving-k flow shear does not support collisions')
      if (hyper_dissipation .and. abs(D_hyper) > 0.0) &
         call mp_abort('moving-k flow shear only supports parallel hyper-dissipation (D_hyper = 0)')
      if (source_option_switch /= source_option_none .or. include_qn_source .or. include_multibox_krook &
         .or. include_neoclassical_terms) call mp_abort('moving-k flow shear does not support sources or neoclassical terms')
      if (omprimfac /= 0.0) call mp_abort('moving-k flow shear requires omprimfac = 0')
      if (.not. hammett_flow_shear) call mp_abort('moving-k flow shear requires hammett_flow_shear = T')
      if (scale_to_phiinit) call mp_abort('moving-k flow shear requires scale_to_phiinit = F')

      ! The committed phase: zero for a new simulation, |tau| <= T_ky/2 after a restart
      if (any(shift_state /= shift_state) .or. any(abs(shift_state) > huge(1.0) / 4)) &
         call mp_abort('moving-k flow shear: the committed phase (shift_state) is not finite')
      if (init_distribution_switch /= init_distribution_option_restart .and. &
          (abs(code_time) > epsilon(0.) .or. any(abs(shift_state) > epsilon(0.)))) &
         call mp_abort('moving-k flow shear requires zero time and phase at the start of a new simulation')
      if (abs(g_exb * g_exbfac) > epsilon(0.)) then
         do iky = 1, naky
            if (aky(iky) == 0.0) then
               if (shift_state(iky) /= 0.0) call mp_abort('moving-k flow shear: nonzero phase of the ky = 0 modes')
            else if (abs(shift_state(iky)) > 0.5 * shift_times(iky) * (1.0 + 1.e-9)) then
               call mp_abort('moving-k flow shear: the committed phase is outside |tau| <= T/2')
            end if
         end do

         ! A time step may never exceed the time needed to shift a mode by one grid spacing
         if (code_dt_max > minval(shift_times)) code_dt_max = minval(shift_times)
      end if

      if (proc0) write (*, *) 'Flow shear: moving-wavevector (continuous) treatment.'

   end subroutine check_moving_k_configuration

   !****************************************************************************
   !                   Physical wavenumbers of a Runge-Kutta stage
   !****************************************************************************
   ! Set kx_phys at the residual shear time tau = shift_state + c*code_dt of the
   ! Runge-Kutta stage with abscissa c, and recalculate kperp2, the Bessel
   ! functions and the field denominators if the wavenumbers changed. This does
   ! not advance the committed phase <shift_state>; c = 0 restores the metric of
   ! the committed state.
   !****************************************************************************
   subroutine prepare_moving_k_stage(c)

      use mp, only: mp_abort
      use arrays, only: flow_shear_moving_k, shift_state, shift_stage, kx_phys
      use grids_kxky, only: naky, nakx, akx, aky
      use grids_time, only: code_dt
      use gk_flow_shear, only: g_exb, g_exbfac, shift_times
      use initialise_arrays, only: refresh_moving_kperp2
      use arrays_gyro_averages, only: refresh_moving_bessel_functions
      use field_equations, only: refresh_moving_two_field_response

      implicit none

      real, intent(in) :: c

      integer :: iky
      real :: gamma
      real, dimension(:), allocatable :: new_stage
      real, dimension(:, :), allocatable :: new_kx

      !-------------------------------------------------------------------------

      if (.not. flow_shear_moving_k) return

      gamma = g_exb * g_exbfac
      if (abs(gamma) > epsilon(0.) .and. (code_dt <= 0.0 .or. code_dt > minval(shift_times))) &
         call mp_abort('moving-k flow shear requires 0 < code_dt <= the minimum shift time')

      allocate (new_stage(naky), new_kx(naky, nakx))
      do iky = 1, naky
         if (aky(iky) == 0.0) then
            new_stage(iky) = 0.0
         else
            new_stage(iky) = shift_state(iky) + c * code_dt
         end if
         new_kx(iky, :) = akx - gamma * aky(iky) * new_stage(iky)
      end do

      ! Identical wavenumbers give identical kernels, so there is nothing to do.
      ! This also leaves <fields_updated> untouched.
      if (allocated(kx_phys)) then
         if (all(new_kx == kx_phys)) then
            shift_stage = new_stage
            deallocate (new_stage, new_kx)
            return
         end if
      end if

      if (.not. allocated(kx_phys)) allocate (kx_phys(naky, nakx))
      if (.not. allocated(shift_stage)) allocate (shift_stage(naky))
      kx_phys = new_kx
      shift_stage = new_stage
      deallocate (new_stage, new_kx)

      call refresh_moving_kperp2
      call refresh_moving_bessel_functions
      call refresh_moving_two_field_response

   end subroutine prepare_moving_k_stage

   !****************************************************************************
   !                    Relabel after an accepted time step
   !****************************************************************************
   ! Advance the committed phase by code_dt and, for every ky row whose phase
   ! has passed half a grid spacing, move the (canonical) distribution function
   ! by an integer number of labels so that |tau| <= T_ky/2 again. Labels that
   ! leave the kx grid are dropped. The physical state is unchanged apart from
   ! the dropped modes. Finally the metric of the new committed state is set.
   !****************************************************************************
   subroutine relabel_moving_k(g)

      use mp, only: mp_abort
      use parallelisation_layouts, only: vmu_lo
      use arrays, only: shift_state
      use grids_kxky, only: naky, nakx, akx, aky
      use grids_z, only: nzgrid, ntubes
      use grids_time, only: code_dt
      use gk_flow_shear, only: g_exb, g_exbfac, shift_times

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      complex, dimension(:, :, :, :), allocatable :: row
      integer :: iky, isrc, idst, nshift
      real :: target, dk

      !-------------------------------------------------------------------------

      ! Without shear nothing moves
      if (abs(g_exb * g_exbfac) <= epsilon(0.)) return

      dk = akx(2)
      allocate (row(nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      do iky = 1, naky
         if (aky(iky) == 0.0) then
            shift_state(iky) = 0.0
            cycle
         end if
         shift_state(iky) = shift_state(iky) + code_dt
         nshift = floor(shift_state(iky) / shift_times(iky) + 0.5)
         if (nshift == 0) cycle
         if (nshift < 0) call mp_abort('moving-k flow shear: negative relabel shift')
         row = g(iky, :, :, :, :)
         g(iky, :, :, :, :) = 0.0
         do isrc = 1, nakx
            target = akx(isrc) - sign(real(nshift) * dk, g_exb * g_exbfac * aky(iky))
            idst = minloc(abs(akx - target), dim=1)
            if (abs(akx(idst) - target) <= 64 * epsilon(0.0) * max(1.0, abs(target))) then
               g(iky, idst, :, :, :) = row(isrc, :, :, :)
            end if
         end do
         shift_state(iky) = shift_state(iky) - real(nshift) * shift_times(iky)
      end do
      deallocate (row)

      call prepare_moving_k_stage(0.0)

   end subroutine relabel_moving_k

   !****************************************************************************
   !             Preconditioner of the implicit streaming solve
   !****************************************************************************
   ! Install, for each low-ky row, the response set nearest to its committed
   ! phase (if additional sets are used), and set the preconditioner scaling.
   !****************************************************************************
   subroutine prepare_moving_k_stream_solve

      use gk_flow_shear, only: g_exb, g_exbfac, moving_k_stream_offset_sets

      implicit none

      !-------------------------------------------------------------------------

      if (moving_k_stream_offset_sets > 0 .and. abs(g_exb * g_exbfac) > epsilon(0.)) call select_response_sets
      call set_stream_defect_scaling

   end subroutine prepare_moving_k_stream_solve

   !****************************************************************************
   !                Scaling of the defect-correction preconditioner
   !****************************************************************************
   ! Ratio of the physical to the label field denominators at the committed
   ! metric, for phi (denominator_fields) and apar (kperp2).
   !****************************************************************************
   subroutine set_stream_defect_scaling

      use arrays, only: kperp2, denominator_fields
      use gk_implicit_terms, only: stream_defect_scale_phi, stream_defect_scale_apar

      implicit none

      real, dimension(:, :, :), allocatable :: ref_denominator, ref_kperp2
      integer :: iky

      !-------------------------------------------------------------------------

      if (.not. allocated(stream_defect_scale_phi)) then
         allocate (stream_defect_scale_phi, mold=label_denominator)
         allocate (stream_defect_scale_apar, mold=label_denominator)
      end if
      allocate (ref_denominator, source=label_denominator)
      allocate (ref_kperp2, source=label_kperp2)

      ! Rows solved with an additional response set use the denominators of that set
      if (sets_built) then
         do iky = 1, size(ref_denominator, 1)
            if (set_installed(iky) /= 0) then
               ref_denominator(iky, :, :) = set_denominator(iky, :, :, set_installed(iky))
               ref_kperp2(iky, :, :) = set_kperp2(iky, :, :, set_installed(iky))
            end if
         end do
      end if

      where (abs(ref_denominator) > epsilon(0.0))
         stream_defect_scale_phi = denominator_fields / ref_denominator
      elsewhere
         stream_defect_scale_phi = 1.0
      end where
      where (ref_kperp2 > epsilon(0.0))
         stream_defect_scale_apar = kperp2(:, :, 1, :) / ref_kperp2
      elsewhere
         stream_defect_scale_apar = 1.0
      end where

      deallocate (ref_denominator, ref_kperp2)

   end subroutine set_stream_defect_scaling

   !****************************************************************************
   !            Response matrix of the implicit streaming at the labels
   !****************************************************************************
   ! A change of the time step recalculates the response matrix with the current
   ! kernels, which may belong to a Runge-Kutta stage. Recalculate it at the
   ! label wavenumbers instead, as at initialisation. The caller restores the
   ! committed metric afterwards.
   !****************************************************************************
   subroutine rebuild_response_on_label_metric

      use arrays, only: kx_phys, shift_stage
      use grids_kxky, only: naky, akx
      use initialise_arrays, only: refresh_moving_kperp2
      use arrays_gyro_averages, only: refresh_moving_bessel_functions
      use field_equations, only: refresh_moving_two_field_response
      use response_matrix, only: init_response_matrix, initialised_response_matrix
      use response_matrix, only: install_response_rows
      use gk_flow_shear, only: moving_k_stream_offset_keep_on_dt_change

      implicit none

      !-------------------------------------------------------------------------

      ! The active low-ky rows may belong to an additional set
      if (sets_built) then
         call install_response_rows(label_set, 1, naky, .true.)
         set_installed = 0
      end if

      kx_phys = spread(akx, 1, naky)
      shift_stage = 0.0
      call refresh_moving_kperp2
      call refresh_moving_bessel_functions
      call refresh_moving_two_field_response
      initialised_response_matrix = .false.
      call init_response_matrix

      ! The additional sets depend on the time step as well. If they are kept, they
      ! are only used as preconditioners, so the converged solution is unchanged.
      if (sets_built .and. .not. moving_k_stream_offset_keep_on_dt_change) call build_response_sets

   end subroutine rebuild_response_on_label_metric

   !****************************************************************************
   !        Streaming response matrices at sub-cell phases of the low-ky rows
   !****************************************************************************
   ! Each set is calculated with the standard construction, at the metric of the
   ! phase tau = m T_ky/M of every row, and stored in its own shared-memory
   ! window. The sets are built at the first use and recalculated in place.
   !****************************************************************************
   subroutine build_response_sets

      use mp, only: proc0, mp_abort
      use mpi, only: mpi_wtime
      use arrays, only: kx_phys, shift_stage, kperp2, denominator_fields
      use grids_kxky, only: naky, akx, aky
      use gk_flow_shear, only: g_exb, g_exbfac, shift_times
      use gk_flow_shear, only: moving_k_stream_offset_sets, moving_k_stream_offset_rows
      use initialise_arrays, only: refresh_moving_kperp2
      use arrays_gyro_averages, only: refresh_moving_bessel_functions
      use field_equations, only: refresh_moving_two_field_response
      use response_matrix, only: save_response_set, build_new_response_set, rebuild_response_set

      implicit none

      integer :: m, iky
      double precision :: start_time
      logical :: first

      !-------------------------------------------------------------------------

      first = .not. sets_built
      nsets = moving_k_stream_offset_sets
      if (mod(nsets, 2) /= 0 .or. nsets < 0) call mp_abort('moving_k_stream_offset_sets should be even and >= 0')
      set_iky_max = min(naky, moving_k_stream_offset_rows + 1)
      start_time = mpi_wtime()

      if (first) then
         call save_response_set(label_set)
         allocate (offset_set(-nsets / 2:nsets / 2))
         allocate (set_denominator(size(label_denominator, 1), size(label_denominator, 2), &
            size(label_denominator, 3), -nsets / 2:nsets / 2))
         allocate (set_kperp2, mold=set_denominator)
         allocate (set_installed(naky))
         set_installed = 0
      end if

      do m = -nsets / 2, nsets / 2
         if (m == 0) cycle
         do iky = 1, naky
            if (aky(iky) == 0.0) then
               kx_phys(iky, :) = akx
            else
               kx_phys(iky, :) = akx - g_exb * g_exbfac * aky(iky) * (real(m) / real(nsets)) * shift_times(iky)
            end if
         end do
         shift_stage = 0.0
         call refresh_moving_kperp2
         call refresh_moving_bessel_functions
         call refresh_moving_two_field_response
         set_denominator(:, :, :, m) = denominator_fields
         set_kperp2(:, :, :, m) = kperp2(:, :, 1, :)
         if (first) then
            call build_new_response_set(offset_set(m), set_iky_max)
         else
            call rebuild_response_set(offset_set(m), set_iky_max)
         end if
      end do
      set_denominator(:, :, :, 0) = label_denominator
      set_kperp2(:, :, :, 0) = label_kperp2
      sets_built = .true.

      ! Restore the metric of the committed state
      call prepare_moving_k_stage(0.0)

      if (proc0) write (*, '(A,I0,A,I0,A,F10.2,A)') ' Moving-k flow shear: built ', nsets, &
         ' additional streaming response sets for ky rows 2 to ', set_iky_max, ' in ', mpi_wtime() - start_time, ' s'

   end subroutine build_response_sets

   !****************************************************************************
   !        Install, per low-ky row, the set nearest to its committed phase
   !****************************************************************************
   subroutine select_response_sets

      use arrays, only: shift_state
      use gk_flow_shear, only: shift_times
      use response_matrix, only: install_response_rows

      implicit none

      integer :: iky, m

      !-------------------------------------------------------------------------

      if (.not. sets_built) call build_response_sets
      do iky = 2, set_iky_max
         m = nint(shift_state(iky) / shift_times(iky) * real(nsets))
         m = max(-nsets / 2, min(nsets / 2, m))
         if (m == set_installed(iky)) cycle
         if (m == 0) then
            call install_response_rows(label_set, iky, iky, .false.)
         else
            call install_response_rows(offset_set(m), iky, iky, .false.)
         end if
         set_installed(iky) = m
      end do

   end subroutine select_response_sets

end module gk_flow_shear_moving_k
