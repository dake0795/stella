!###############################################################################
!################## EXPLICIT TERMS OF THE GYROKINETIC EQUATION #################
!###############################################################################
! 
! This module adds the explicit gyrokinetic terms to the right-hand-side of the 
! gyrokinetic equation, in order to advance the distribution function in time.
! 
! First advance the distribution function <g> in time using the gyrokinetic equation.
! Next, advance the fields (electrostatic potential <phi>, as well as the electromagnetic
! fields <apar> and <bpar>) in time using the quasi-neutrality condition.
! 
!###############################################################################
module gyrokinetic_equation_explicit

   ! Load debug flags
   use debug_flags, only: debug => time_advance_debug

   implicit none
   
   ! Make routines available to stella.f90
   public :: init_gyrokinetic_equation_explicit
   public :: advance_distribution_function_using_explicit_gyrokinetic_terms

   private

   ! Re-synchronise the duplicated twist-and-shift junction points after each
   ! accepted step (junction_sync_step) and before each Runge-Kutta stage
   ! (junction_sync_stage), see <explicit_junction_sync> in the flow_shear namelist
   logical :: junction_sync_step = .false., junction_sync_stage = .false.

   ! Implicit mirror and parallel streaming applied as a split after the accepted
   ! explicit step, instead of in the implicit time advance (see derive_implicit_split)
   logical :: mirror_split = .false., stream_split = .false.

contains

   !****************************************************************************
   !                    INITIALISE THE EXPLICIT TIME ADVANCE
   !****************************************************************************
   subroutine init_gyrokinetic_equation_explicit

      use mp, only: mp_abort
      use gk_flow_shear, only: explicit_junction_sync
      use gk_flow_shear, only: flow_shear_explicit_remap
      use parameters_numerical, only: fully_explicit
      use arrays, only: flow_shear_moving_k
      use gk_flow_shear_moving_k, only: init_flow_shear_moving_k
      use diagnostics_energy_budget, only: init_energy_budget

      implicit none

      !-------------------------------------------------------------------------

      ! Synchronisation of the twist-and-shift junction points
      select case (trim(explicit_junction_sync))
      case ('none')
         junction_sync_step = .false.
         junction_sync_stage = .false.
      case ('step')
         junction_sync_step = .true.
         junction_sync_stage = .false.
      case ('stage')
         junction_sync_step = .true.
         junction_sync_stage = .true.
      case default
         call mp_abort('explicit_junction_sync should be none, step or stage')
      end select

      ! The flow-shear remap of the explicit time advance, and the moving-wavevector
      ! flow shear, are only used when the implicit time advance is never used
      if (flow_shear_explicit_remap .or. flow_shear_moving_k) then
         call derive_implicit_split
         call keep_fully_explicit_with_hyper_dissipation
      end if
      if (flow_shear_explicit_remap .and. .not. fully_explicit) &
         call mp_abort('flow_shear_explicit_remap requires a fully explicit simulation')
      if (flow_shear_moving_k) call init_flow_shear_moving_k

      ! Optional free-energy budget diagnostic
      call init_energy_budget

   end subroutine init_gyrokinetic_equation_explicit

   !****************************************************************************
   !            Fully explicit simulations with parallel hyper-dissipation
   !****************************************************************************
   ! Setting <hyper_dissipation> = .true. switches off <fully_explicit>, because
   ! the perpendicular hyper-dissipation (D_hyper) is applied in the implicit
   ! solve. If D_hyper = 0 and no other term is treated implicitly, the only
   ! active hyper-dissipation is the parallel one (hyp_zed, hyp_vpa), which is
   ! part of the explicit right-hand side, so the simulation is fully explicit.
   !****************************************************************************
   subroutine keep_fully_explicit_with_hyper_dissipation

      use parameters_numerical, only: fully_explicit, fully_implicit, flip_flop
      use parameters_numerical, only: stream_implicit, mirror_implicit
      use parameters_numerical, only: drifts_implicit, driftkinetic_implicit
      use dissipation_and_collisions, only: hyper_dissipation
      use dissipation_hyper, only: D_hyper

      implicit none

      !-------------------------------------------------------------------------

      if (fully_explicit .or. .not. hyper_dissipation .or. abs(D_hyper) > 0.0) return
      if (stream_implicit .or. mirror_implicit .or. drifts_implicit .or. driftkinetic_implicit .or. fully_implicit) return
      fully_explicit = .true.
      flip_flop = .false.

   end subroutine keep_fully_explicit_with_hyper_dissipation

   !****************************************************************************
   !               Implicit mirror and streaming as a split step
   !****************************************************************************
   ! The implicit time advance also holds the flow-shear remap. When the remap is
   ! done by the explicit time advance, the implicit mirror and parallel streaming
   ! are instead applied as a first-order (Lie) split after each accepted explicit
   ! step, in the order of the implicit time advance without flip-flopping:
   ! explicit step -> mirror -> parallel streaming. The implicit time advance is
   ! then not used, so <fully_explicit> is set to true. Only the mirror and the
   ! parallel streaming are supported as implicit terms.
   !****************************************************************************
   subroutine derive_implicit_split

      use mp, only: proc0, mp_abort
      use parameters_numerical, only: fully_explicit, fully_implicit, flip_flop
      use parameters_numerical, only: stream_implicit, mirror_implicit, mirror_semi_lagrange
      use parameters_numerical, only: drifts_implicit, driftkinetic_implicit, nitt
      use parameters_numerical, only: explicit_algorithm_switch, explicit_algorithm_rk4
      use parameters_physics, only: include_mirror, include_parallel_streaming
      use parameters_physics, only: full_flux_surface, radial_variation
      use dissipation_and_collisions, only: include_collisions, hyper_dissipation
      use dissipation_hyper, only: D_hyper
      use neoclassical_terms, only: include_neoclassical_terms

      implicit none

      !-------------------------------------------------------------------------

      mirror_split = mirror_implicit .and. include_mirror
      stream_split = stream_implicit .and. include_parallel_streaming
      if (.not. (mirror_split .or. stream_split)) return

      if (mirror_split .and. mirror_semi_lagrange) &
         call mp_abort('the implicit mirror split requires mirror_semi_lagrange = .false.')
      if (drifts_implicit .or. driftkinetic_implicit .or. fully_implicit) &
         call mp_abort('only the mirror and the parallel streaming can be split from the explicit step')
      if (stream_split .and. nitt /= 1) call mp_abort('the implicit streaming split requires nitt = 1')
      if (flip_flop) call mp_abort('the implicit split requires flip_flop = .false.')
      if (include_collisions) call mp_abort('the implicit split does not support collisions')
      if (hyper_dissipation .and. abs(D_hyper) > 0.0) &
         call mp_abort('the implicit split only supports parallel hyper-dissipation (D_hyper = 0)')
      if (full_flux_surface .or. radial_variation .or. include_neoclassical_terms) &
         call mp_abort('the implicit split requires a local flux tube without neoclassical terms')
      if (explicit_algorithm_switch /= explicit_algorithm_rk4) &
         call mp_abort('the implicit split requires explicit_algorithm = rk4')

      fully_explicit = .true.
      if (proc0) write (*, *) 'Implicit mirror and/or parallel streaming applied as a split after the explicit step.'

   end subroutine derive_implicit_split

   !****************************************************************************
   !                       EXPLICIT TIME ADVANCE SUBROUTINES
   !****************************************************************************
   ! Advance the guiding center distribution equation <g> in k-space in time using 
   ! the gyrokinetic equation. Specifically, the distribution function <g> is updated 
   ! based on all of the terms in the GKE that are advanced explicitly in time.
   !****************************************************************************
   subroutine advance_distribution_function_using_explicit_gyrokinetic_terms(g, restart_time_step, istep)

      ! Parallelisation
      use mp, only: proc0
      use job_manage, only: time_message
      use parallelisation_layouts, only: vmu_lo, iv_idx
      use timers, only: time_gke
      
      ! Fields
      use parameters_physics, only: include_apar
      use arrays_fields, only: phi, apar, bpar
      use field_equations, only: advance_fields, fields_updated

      ! Flow shear
      use parameters_numerical, only: fully_explicit
      use gk_flow_shear, only: prp_shear_enabled, flow_shear_explicit_remap
      use gk_flow_shear, only: advance_perp_flow_shear
      use arrays, only: flow_shear_moving_k
      use gk_flow_shear_moving_k, only: prepare_moving_k_stage, relabel_moving_k
      use gk_flow_shear_moving_k, only: rebuild_response_on_label_metric
      
      ! Free-energy budget
      use diagnostics_energy_budget, only: energy_budget_on, budget_start_step, budget_incoming_sync
      use diagnostics_energy_budget, only: budget_commit_step, budget_sync_start, budget_sync_end
      use diagnostics_energy_budget, only: budget_remap_start, budget_remap_end, budget_end_step
      
      ! Grids
      use grids_z, only: nzgrid
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: periodic
      use grids_extended_zgrid, only: phase_shift
      use gk_parallel_streaming, only: stream_sign
      
      ! Calculations
      use calculations_tofrom_ghf, only: gbar_to_g
      
      ! Numerical time advance schemes
      use parameters_numerical, only: explicit_algorithm_switch
      use parameters_numerical, only: explicit_algorithm_rk3
      use parameters_numerical, only: explicit_algorithm_rk2
      use parameters_numerical, only: explicit_algorithm_rk4
      use parameters_numerical, only: explicit_algorithm_euler

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: ivmu, iv, sgn, iky
      complex, dimension(:, :, :, :, :), allocatable :: gpre

      !-------------------------------------------------------------------------

      ! Start the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

      ! Free-energy budget: start from the free energy of the last accepted state
      call budget_start_step(g)
      
      ! If the fields are not already updated, then update them
      if (include_apar) then
         call advance_fields(g, phi, apar, bpar, dist='g')
      end if

      ! Incoming distribution function is g = h - (Z F0/T) (J0 phi + 4 mu (T/Z) (J1/gamma) bpar)
      ! If <include_apar> = T, convert from g to gbar = g + Z F0/T (2J0 vpa vth apar), as gbar appears in time derivatives
      if (include_apar) then
         call gbar_to_g(g, apar, -1.0)
      end if

      ! Make the junction points of the incoming state single-valued
      if (junction_sync_stage) then
         if (energy_budget_on) allocate (gpre, source=g)
         call synchronise_zgrid_junctions(g)
         if (energy_budget_on) then
            call budget_incoming_sync(gpre, g)
            deallocate (gpre)
         end if
      end if

      ! Use a numerical time advance scheme to advance the distribution function
      ! in time, based on the explicit terms in the gyrokinetic equation. 
      ! Choose between Forward Euler, or 2nd, 3rd or 4rd order Runge-Kutta schemes
      ! The 3rd order Runge-Kutta scheme is the default option.
      select case (explicit_algorithm_switch)
      case (explicit_algorithm_euler)
         call advance_explicit_euler(g, restart_time_step, istep)
      case (explicit_algorithm_rk2)
         call advance_explicit_rk2(g, restart_time_step, istep)
      case (explicit_algorithm_rk3)
         call advance_explicit_rk3(g, restart_time_step, istep)
      case (explicit_algorithm_rk4)
         call advance_explicit_rk4(g, restart_time_step, istep)
      end select

      ! Free-energy budget: commit the stage powers of an accepted step
      if (.not. restart_time_step) call budget_commit_step

      ! Moving-wavevector flow shear: relabel after an accepted step. If the time
      ! step is restarted, restore the metric of the committed state instead.
      if (flow_shear_moving_k) then
         if (restart_time_step) then
            if (stream_split) call rebuild_response_on_label_metric
            call prepare_moving_k_stage(0.0)
            if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')
            return
         end if
         call relabel_moving_k(g)
      end if

      ! Free-energy budget: free energy before the end-of-step projections
      if (.not. restart_time_step) call budget_sync_start(g)

      ! Re-synchronise the junction points of the accepted state. The parallel
      ! vector potential is single-valued at the junctions, so this applies to gbar too.
      if (junction_sync_step .and. .not. restart_time_step) call synchronise_zgrid_junctions(g)

      ! If the fields are not already updated, then update them
      if (include_apar) then
         if (junction_sync_step) fields_updated = .false.
         call advance_fields(g, phi, apar, bpar, dist='gbar')
      end if
      
      ! Later, the implicit solve will use <g> rather than <gbar> to advance the 
      ! distribution function in time. Therefore, convert <gbar> to <g> again
      if (include_apar) then
         call gbar_to_g(g, apar, 1.0)
      end if

      ! Enforce periodicity for periodic (including zonal) modes
      ! <stream_sign> > 0 corresponds to dz/dt < 0
      do iky = 1, naky
         if (periodic(iky)) then
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               iv = iv_idx(vmu_lo, ivmu)
               sgn = stream_sign(iv)
               g(iky, :, sgn * nzgrid, :, ivmu) = g(iky, :, -sgn * nzgrid, :, ivmu) * phase_shift(iky)**(-sgn)
            end do
         end if
      end do

      ! The fields above were solved before the periodic end points were copied;
      ! solve them again so that they are consistent with the final <g>
      if (junction_sync_step) then
         fields_updated = .false.
         call advance_fields(g, phi, apar, bpar, dist='g')
      end if

      ! Free-energy budget: change by the junction synchronisation and the end-point copy
      if (.not. restart_time_step) call budget_sync_end(g)

      ! Remap the perpendicular flow shear in fully explicit simulations. The remap
      ! vacates the edge label of the shifted modes, which changes the pairing at
      ! the ends of the linked chains, so the junction points are synchronised again.
      if (flow_shear_explicit_remap .and. fully_explicit .and. prp_shear_enabled .and. .not. restart_time_step) then
         call budget_remap_start(g)
         call advance_perp_flow_shear(g)
         if (junction_sync_step) call synchronise_zgrid_junctions(g)
         fields_updated = .false.
         call advance_fields(g, phi, apar, bpar, dist='g')
         call budget_remap_end(g)
      end if

      ! Implicit mirror and parallel streaming, split from the explicit step
      if ((mirror_split .or. stream_split) .and. .not. restart_time_step) call advance_implicit_split(g)

      ! Free-energy budget: free energy of the accepted state
      if (.not. restart_time_step) call budget_end_step(g)

      ! Stop the timer for the explicit part of the solve
      if (proc0) call time_message(.false., time_gke(:, 8), ' explicit')

   end subroutine advance_distribution_function_using_explicit_gyrokinetic_terms

   !****************************************************************************
   !                 IMPLICIT MIRROR AND STREAMING SPLIT STEP
   !****************************************************************************
   ! On entry <g> is the accepted state and the fields have been solved from it,
   ! so the implicit mirror uses the stored <apar> to convert between g and gbar,
   ! as in the implicit time advance.
   !****************************************************************************
   subroutine advance_implicit_split(g)

      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo, iv_idx

      ! Grids
      use grids_z, only: nzgrid
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: periodic, phase_shift
      use gk_parallel_streaming, only: stream_sign

      ! Fields
      use arrays_fields, only: phi, apar, bpar
      use field_equations, only: advance_fields, fields_updated

      ! Implicit terms
      use gk_mirror, only: advance_mirror_implicit, refresh_mirror_response_for_metric
      use gk_implicit_terms, only: advance_implicit_terms, stream_defect_active

      ! Moving-wavevector flow shear
      use arrays, only: flow_shear_moving_k
      use gk_flow_shear_moving_k, only: prepare_moving_k_stream_solve

      ! Free-energy budget
      use diagnostics_energy_budget, only: budget_split_start, budget_split_end, budget_split_sync_end
      use diagnostics_energy_budget, only: ib_mi_split, ib_st_split

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      ! Local variables
      integer :: iky, ivmu, iv, sgn

      !-------------------------------------------------------------------------

      if (mirror_split) then
         call budget_split_start(g)
         ! Moving-wavevector flow shear: electromagnetic response at the committed metric
         call refresh_mirror_response_for_metric
         call advance_mirror_implicit(.false., g, apar)
         call budget_split_end(g, ib_mi_split)

         ! The mirror acts point-wise in (kx, ky, z), so restore the single-valued
         ! junction points and the periodic end points
         if (junction_sync_step) then
            call synchronise_zgrid_junctions(g)
         else
            do iky = 1, naky
               if (periodic(iky)) then
                  do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
                     iv = iv_idx(vmu_lo, ivmu)
                     sgn = stream_sign(iv)
                     g(iky, :, sgn * nzgrid, :, ivmu) = g(iky, :, -sgn * nzgrid, :, ivmu) * phase_shift(iky)**(-sgn)
                  end do
               end if
            end do
         end if
         fields_updated = .false.
         call advance_fields(g, phi, apar, bpar, dist='g')
         call budget_split_sync_end(g)
      end if

      if (stream_split) then
         call budget_split_start(g)
         call advance_fields(g, phi, apar, bpar, dist='g')
         ! Moving-wavevector flow shear: iterate to the solution at the physical wavenumbers
         if (flow_shear_moving_k) call prepare_moving_k_stream_solve
         stream_defect_active = flow_shear_moving_k
         call advance_implicit_terms(g, phi, apar, bpar)
         stream_defect_active = .false.
         fields_updated = .false.
         call advance_fields(g, phi, apar, bpar, dist='g')
         if (flow_shear_moving_k) call write_stream_defect_statistics
         call budget_split_end(g, ib_st_split)
      end if

   end subroutine advance_implicit_split

   !****************************************************************************
   !     Iterations of the defect correction of the implicit streaming
   !****************************************************************************
   subroutine write_stream_defect_statistics

      use mp, only: proc0
      use gk_implicit_terms, only: stream_defect_iters

      implicit none

      integer, save :: ncalls = 0, total_iters = 0, max_iters = 0

      !-------------------------------------------------------------------------

      ncalls = ncalls + 1
      total_iters = total_iters + stream_defect_iters
      max_iters = max(max_iters, stream_defect_iters)
      if (proc0 .and. mod(ncalls, 500) == 0) write (*, '(A,I0,A,F8.3,A,I0)') &
         ' Implicit streaming defect correction: steps = ', ncalls, ', mean iterations = ', &
         real(total_iters) / real(ncalls), ', maximum = ', max_iters

   end subroutine write_stream_defect_statistics

   !****************************************************************************
   !                 SYNCHRONISE TWIST-AND-SHIFT JUNCTION POINTS
   !****************************************************************************
   ! In a linked (twist-and-shift) chain the junction point between two
   ! segments is stored twice: at z = +pi of segment i and at z = -pi of
   ! segment i+1. The implicit solve maps <g> onto the extended z-grid and back
   ! every step, which writes a single value to both copies. The explicit
   ! schemes never do this, so the two copies evolve independently and an
   ! antisymmetric junction mode can grow. Here the same projection is applied
   ! (map_to_extended_zgrid keeps the value of the preceding segment and
   ! map_from_extended_zgrid writes it to both copies). For periodic modes the
   ! upwind end point is copied from the other end, as at the end of a step.
   !****************************************************************************
   subroutine synchronise_zgrid_junctions(g)

      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo, iv_idx

      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky
      use grids_extended_zgrid, only: neigen, nsegments, nzed_segment
      use grids_extended_zgrid, only: periodic, phase_shift
      use grids_extended_zgrid, only: map_to_extended_zgrid, map_from_extended_zgrid
      use gk_parallel_streaming, only: stream_sign

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g

      ! Local variables
      complex, dimension(:), allocatable :: gext
      integer :: ivmu, iky, ie, it, ulim, sgn

      !-------------------------------------------------------------------------

      allocate (gext((maxval(nsegments) + 1) * (nzed_segment + 1)))

      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do it = 1, ntubes
            do iky = 1, naky
               if (periodic(iky)) then
                  sgn = stream_sign(iv_idx(vmu_lo, ivmu))
                  g(iky, :, sgn * nzgrid, it, ivmu) = g(iky, :, -sgn * nzgrid, it, ivmu) * phase_shift(iky)**(-sgn)
                  cycle
               end if
               do ie = 1, neigen(iky)
                  if (nsegments(ie, iky) < 2) cycle
                  call map_to_extended_zgrid(it, ie, iky, g(iky, :, :, :, ivmu), gext, ulim)
                  call map_from_extended_zgrid(it, ie, iky, gext(:ulim), g(iky, :, :, :, ivmu))
               end do
            end do
         end do
      end do

      deallocate (gext)

   end subroutine synchronise_zgrid_junctions
   
   !****************************************************************************
   !                NEEDED FOR ALL EXPLICIT TIME ADVANCE SUBROUTINES
   !****************************************************************************
   ! solve_gyrokinetic_equation_explicit accepts as argument pdf, the guiding center distribution function in k-space,
   ! and returns rhs_ky, the right-hand side of the gyrokinetic equation in k-space;
   ! i.e., if dg/dt = r, then rhs_ky = r*dt;
   ! note that if include_apar = T, then the input pdf is actually gbar = g + (Ze/T)*(vpa/c)*<Apar>*F0
   !****************************************************************************
   subroutine add_explicit_gyrokinetic_terms(pdf, rhs_ky, restart_time_step, istep)

      ! Parallelisation
      use job_manage, only: time_message
      use multibox, only: add_multibox_krook
      use parallelisation_layouts, only: vmu_lo
      use calculations_transforms, only: transform_y2ky

      ! Fields
      use arrays_fields, only: phi, apar, bpar
      
      ! Distribution function
      use arrays_distribution_function, only: phi_gyro
      
      ! Calculations
      use arrays_gyro_averages, only: j0_ffs
      use calculations_gyro_averages, only: gyro_average
      use calculations_kxky, only: swap_kxky_back
      use calculations_tofrom_ghf, only: gbar_to_g 

      ! Physics flags
      use parameters_physics, only: include_parallel_nonlinearity
      use parameters_physics, only: include_parallel_streaming
      use parameters_physics, only: include_mirror
      use parameters_physics, only: include_apar
      use parameters_physics, only: include_nonlinear
      use parameters_physics, only: full_flux_surface
      use parameters_physics, only: radial_variation
      use parameters_physics, only: xdriftknob, ydriftknob
      use dissipation_and_collisions, only: include_collisions
      use dissipation_and_collisions, only: hyper_dissipation
      use gk_flow_shear, only: g_exb
      
      ! Krook source
      use parameters_multibox, only: include_multibox_krook
      use gk_sources, only: source_option_switch
      use gk_sources, only: source_option_krook
      use gk_sources, only: add_krook_operator
      
      ! Numerical schemes
      use parameters_numerical, only: stream_implicit
      use parameters_numerical, only: mirror_implicit
      use parameters_numerical, only: drifts_implicit
      use dissipation_and_collisions, only: advance_collisions_explicit
      use dissipation_and_collisions, only: collisions_implicit

      ! Grids
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: zonal_mode, akx
      use grids_kxky, only: ikx_max, ny, naky_all
      
      ! Routines to add terms to the gyrokinetic equation
      use gk_parallel_streaming, only: advance_parallel_streaming_explicit
      use gk_mirror, only: advance_mirror_explicit
      use gk_flow_shear, only: advance_parallel_flow_shear
      use gk_drive, only: advance_wstar_explicit
      use gk_magnetic_drift, only: advance_wdriftx_explicit, advance_wdrifty_explicit
      use gk_parallel_nonlinearity, only: advance_parallel_nonlinearity
      use gk_radial_variation, only: advance_radial_variation
      use gk_ExB_nonlinearity, only: advance_ExB_nonlinearity
      use field_equations_radialvariation, only: get_radial_correction
      use field_equations, only: advance_fields
      use field_equations, only: fields_updated

      ! Free-energy budget
      use diagnostics_energy_budget, only: budget_shear_stage, budget_stage_start, budget_add_term, budget_stage_end
      use diagnostics_energy_budget, only: ib_exb, ib_prl, ib_mirror, ib_wdy, ib_wdx, ib_wstar, ib_stream, ib_hyper, ib_other

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: pdf
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(out), target :: rhs_ky
      logical, intent(out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      complex, dimension(:, :, :, :, :), allocatable, target :: rhs_y
      complex, dimension(:, :, :, :, :), pointer :: rhs
      complex, dimension(:, :), allocatable :: rhs_ky_swap
      integer :: iz, it, ivmu

      !-------------------------------------------------------------------------

      ! Initialise the right-hand-side of the gyrokinetic equation to zero
      rhs_ky = 0.

      ! If full_flux_surface = .true., then initially obtain the RHS of the GKE in alpha-space;
      ! We will later inverse Fourier transform to get RHS in k_alpha-space
      if (full_flux_surface) then
         ! rhs_ky will always be needed as the array returned by the subroutine,
         ! but intermediate array rhs_y (RHS of gke in alpha-space) only needed for full_flux_surface = .true.
         allocate (rhs_y(ny, ikx_max, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
         rhs_y = 0.
         ! rhs is array referred to for both flux tube and full-flux-surface simulations;
         ! for full-flux-surface it should point to rhs_y
         rhs => rhs_y
      else
         ! rhs is array referred to for both flux tube and full-flux-surface simulations;
         ! for flux tube it should point to rhs_ky
         rhs => rhs_ky
      end if

      ! Start with g in k-space and (ky,kx,z) local
      ! obtain fields corresponding to g
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_fields'

      ! If advancing apar, then gbar is evolved in time rather than g
      if (include_apar) then
         call advance_fields(pdf, phi, apar, bpar, dist='gbar')

         ! Free-energy budget: work of the moving metric in this stage
         call budget_shear_stage(pdf)

         ! Convert from gbar to g = h - (Z F0/T)( J0 phi + 4 mu (T/Z) (J1/gamma) bpar),
         ! as all terms on RHS of GKE use g rather than gbar
         call gbar_to_g(pdf, apar, 1.0)
      else
         call advance_fields(pdf, phi, apar, bpar, dist='g')
      end if

      if (radial_variation) call get_radial_correction(pdf, phi, dist='gbar')

      ! Free-energy budget: h of this stage, for the power of each term
      call budget_stage_start(pdf)

      ! Obtain the gyro-average of the electrostatic potential phi and store in phi_gyro;
      ! this can be a particularly costly operation when simulating a full flux surface
      ! due to the coupling of different k-alphas inherent in the gyro-average;
      ! calculate once here to avoid repeated calculation later
      ! TODO-GA : can this be spec up??
      if (full_flux_surface) call gyro_average(phi, phi_gyro, j0_ffs)

      !! INSERT TEST HERE TO SEE IF dg/dy, dg/dx, d<phi>/dy, d<phi>/dx WILL BE NEEDED
      !! IF SO, PRE-COMPUTE ONCE HERE

      ! Default is to continue with same time step size.
      ! if estimated CFL condition for nonlinear terms is violated
      ! then restart_time_step will be set to .true.
      restart_time_step = .false.
      
      ! Calculate and add ExB nonlinearity to RHS of GK eqn
      ! do this first, as the CFL condition may require a change in time step
      ! and thus recomputation of mirror, wdrift, wstar, and parstream
      if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_ExB_nonlinearity'
      if (include_nonlinear) call advance_ExB_nonlinearity(pdf, rhs, restart_time_step, istep)
      if (.not. restart_time_step) call budget_add_term(ib_exb, rhs_ky)

      ! Include contribution from the parallel nonlinearity (aka turbulent acceleration)
      if (include_parallel_nonlinearity .and. .not. restart_time_step) &
         call advance_parallel_nonlinearity(pdf, rhs, restart_time_step)

      if (.not. restart_time_step) then

         ! Include contribution from perp flow shear in the parallel component of the toroidal flow
         if ((g_exb**2) > epsilon(0.0)) call advance_parallel_flow_shear(rhs)
         call budget_add_term(ib_prl, rhs_ky)

         ! Calculate and add mirror term to RHS of GK eqn
         if (include_mirror .and. .not. mirror_implicit) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_mirror_explicit'
            call advance_mirror_explicit(pdf, rhs)
            call budget_add_term(ib_mirror, rhs_ky)
         end if

         if (.not. drifts_implicit) then
            ! Calculate and add alpha-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_wdrifty_explicit'
            if (abs(ydriftknob) > epsilon(0.0)) then
               call advance_wdrifty_explicit(pdf, phi, bpar, rhs)
               call budget_add_term(ib_wdy, rhs_ky)
            end if

            ! Calculate and add psi-component of magnetic drift term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_wdriftx_explicit'
            if (abs(xdriftknob) > epsilon(0.0)) then
               call advance_wdriftx_explicit(pdf, phi, bpar, rhs)
               call budget_add_term(ib_wdx, rhs_ky)
            end if
            
            ! Calculate and add omega_* term to RHS of GK eqn
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_wstar_explicit'
            call advance_wstar_explicit(phi, rhs)
            call budget_add_term(ib_wstar, rhs_ky)
         end if

         ! Calculate and add contribution from collisions to RHS of GK eqn
         if (include_collisions .and. .not. collisions_implicit) call advance_collisions_explicit(pdf, phi, bpar, rhs)

         ! Calculate and add parallel streaming term to RHS of GK eqn
         if (include_parallel_streaming .and. (.not. stream_implicit)) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::advance_parallel_streaming_explicit'
            call advance_parallel_streaming_explicit(pdf, phi, bpar, rhs)
            call budget_add_term(ib_stream, rhs_ky)
         end if
         
         if (hyper_dissipation) then
            call advance_hyper_explicit(pdf, rhs)
            call budget_add_term(ib_hyper, rhs_ky)
         end if

         ! If simulating a full flux surface (flux annulus), all terms to this point have been calculated
         ! in real-space in alpha (y); transform to kalpha (ky) space before adding to RHS of GKE.
         ! NB: it may be that for fully explicit calculation, this transform can be eliminated with additional code changes
         if (full_flux_surface) then
            if (debug) write (*, *) 'time_advance::advance_stella::advance_explicit::solve_gyrokinetic_equation_explicit::transform_y2ky'
            allocate (rhs_ky_swap(naky_all, ikx_max))
            it = 1
            do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
               do iz = -nzgrid, nzgrid
                  call transform_y2ky(rhs_y(:, :, iz, it, ivmu), rhs_ky_swap)
                  call swap_kxky_back(rhs_ky_swap, rhs_ky(:, :, iz, it, ivmu))
               end do
               ! ensure that the kx=ky=0 mode is zeroed out
               if (zonal_mode(1) .and. akx(1) < epsilon(0.)) then
                  rhs_ky(1, 1, :, it, ivmu) = 0.0
               end if
            end do
            deallocate (rhs_ky_swap)
         end if

         if (radial_variation) call advance_radial_variation(pdf, rhs)

         if (source_option_switch == source_option_krook) call add_krook_operator(pdf, rhs)

         if (include_multibox_krook) call add_multibox_krook(pdf, rhs)

         ! Free-energy budget: any remaining terms
         call budget_add_term(ib_other, rhs_ky)

      end if
      call budget_stage_end

      ! if advancing apar, need to convert input pdf back from g to gbar
      if (include_apar) call gbar_to_g(pdf, apar, -1.0)

      fields_updated = .false.

      if (allocated(rhs_y)) deallocate (rhs_y)
      nullify (rhs)

   end subroutine add_explicit_gyrokinetic_terms

   !****************************************************************************
   !                                      Title
   !****************************************************************************
   subroutine advance_hyper_explicit(gin, gout)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      use dissipation_hyper, only: advance_hyper_vpa, advance_hyper_zed
      use dissipation_hyper, only: hyp_zed, hyp_vpa

      implicit none
      
      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gin
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: gout

      ! Local variables
      complex, dimension(:, :, :, :, :), allocatable :: dg

      !-------------------------------------------------------------------------

      allocate (dg(naky, nakx, -nzgrid:nzgrid, ntubes, vmu_lo%llim_proc:vmu_lo%ulim_alloc)); dg = 0.0

      if (hyp_zed) then
         call advance_hyper_zed(gin, dg)
         gout = gout + dg
      end if
      if (hyp_vpa) then
         call advance_hyper_vpa(gin, dg)
      end if
      deallocate (dg)
      
   end subroutine advance_hyper_explicit
   
   
!###############################################################################
!######################## EXPLICIT TIME ADVANCE SCHEMES ########################
!###############################################################################

   !****************************************************************************
   !                        EXPLICIT EULER TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses forward Euler to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_euler(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use gk_radial_variation, only: mb_communicate
      use parameters_multibox, only: rk_step

      ! Distribution function
      use arrays_distribution_function, only: g0
      
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      !-------------------------------------------------------------------------

      ! <rk_step> is only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g

      ! Calculate the right-hand-side of the gyrokinetic equation and store it in <g>
      call add_explicit_gyrokinetic_terms(g0, g, restart_time_step, istep)

      ! Add the right-hand-side of the gyrokinetic equation to <g0>
      g = g0 + g

   end subroutine advance_explicit_euler

   !****************************************************************************
   !                         EXPLICIT RK2 TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses strong stability-preserving rk2 to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_rk2(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate
      
      ! Distribution function
      use arrays_distribution_function, only: g0, g1

      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: counter

      !-------------------------------------------------------------------------

      ! <rk_step> is only true if running in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g
      
      ! Count the iterations
      counter = 1

      ! SSP rk2 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! add_explicit_gyrokinetic_terms returns rhs*dt
      do while (counter <= 2)
      
         ! Second order Runge-Kutta Scheme
         select case (counter)
         case (1)
            call add_explicit_gyrokinetic_terms(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            if (junction_sync_stage) call synchronise_zgrid_junctions(g1)
            call add_explicit_gyrokinetic_terms(g1, g, restart_time_step, istep)
         end select
         
         ! If the code_dt is reset, we need to quit this loop and restart the timestep again
         ! Otherwise increase the <counter> that keeps track of the iterations
         if (restart_time_step) then 
            counter = 10
         else
            counter = counter + 1
         end if
         
      end do

      ! This is g at intermediate time level
      g = 0.5 * g0 + 0.5 * (g1 + g)

   end subroutine advance_explicit_rk2

   !****************************************************************************
   !                         EXPLICIT RK3 TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses strong stability-preserving rk3 to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_rk3(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate

      ! Distribution function
      use arrays_distribution_function, only: g0, g1, g2
      
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: counter

      !-------------------------------------------------------------------------

      ! rk_STEP = false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g
      
      ! Count the iterations
      counter = 1

      ! SSP rk3 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! add_explicit_gyrokinetic_terms returns rhs*dt
      do while (counter <= 3)
      
         ! Third order Runge-Kutta Scheme
         select case (counter)
         case (1)
            call add_explicit_gyrokinetic_terms(g0, g1, restart_time_step, istep)
         case (2)
            g1 = g0 + g1
            if (rk_step) call mb_communicate(g1)
            if (junction_sync_stage) call synchronise_zgrid_junctions(g1)
            call add_explicit_gyrokinetic_terms(g1, g2, restart_time_step, istep)
         case (3)
            g2 = g1 + g2
            if (rk_step) call mb_communicate(g2)
            if (junction_sync_stage) call synchronise_zgrid_junctions(g2)
            call add_explicit_gyrokinetic_terms(g2, g, restart_time_step, istep)
         end select
         
         ! If the code_dt is reset, we need to quit this loop and restart the timestep again
         ! Otherwise increase the <counter> that keeps track of the iterations
         if (restart_time_step) then
            counter = 10
         else
            counter = counter + 1
         end if
      end do

      ! This is g at intermediate time level
      g = g0 / 3.+0.5 * g1 + (g2 + g) / 6.

   end subroutine advance_explicit_rk3

   !****************************************************************************
   !                         EXPLICIT RK4 TIME ADVANCE SUBROUTINE
   !****************************************************************************
   ! Uses rk4 to advance one time step.
   !****************************************************************************
   subroutine advance_explicit_rk4(g, restart_time_step, istep)
   
      ! Parallelisation
      use parallelisation_layouts, only: vmu_lo
      use parameters_multibox, only: rk_step
      use gk_radial_variation, only: mb_communicate

      ! Distribution function
      use arrays_distribution_function, only: g0, g1, g2, g3

      ! Moving-wavevector flow shear: metric of each Runge-Kutta stage
      use gk_flow_shear_moving_k, only: prepare_moving_k_stage
      use diagnostics_energy_budget, only: budget_set_stage_weight
      
      ! Grids
      use grids_z, only: nzgrid

      implicit none

      ! Arguments
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in out) :: g
      logical, intent(in out) :: restart_time_step
      integer, intent(in) :: istep

      ! Local variables
      integer :: counter

      !-------------------------------------------------------------------------

      ! rk_step is false unless in multibox mode
      if (rk_step) call mb_communicate(g)

      ! Store the distribution function
      g0 = g
      
      ! Count the iterations
      counter = 1

      ! RK4 algorithm to advance explicit part of code
      ! if GK equation written as dg/dt = rhs - vpar . grad h,
      ! add_explicit_gyrokinetic_terms returns rhs*dt
      do while (counter <= 4)
      
         ! Fourth order Runge-Kutta Scheme
         select case (counter)
         case (1)
            call prepare_moving_k_stage(0.0)
            call budget_set_stage_weight(1.0 / 6.0)
            call add_explicit_gyrokinetic_terms(g0, g1, restart_time_step, istep)
         case (2)
            g3 = g0 + 0.5 * g1 ! g1 is h*k1
            if (rk_step) call mb_communicate(g3)
            if (junction_sync_stage) call synchronise_zgrid_junctions(g3)
            call prepare_moving_k_stage(0.5)
            call budget_set_stage_weight(1.0 / 3.0)
            call add_explicit_gyrokinetic_terms(g3, g2, restart_time_step, istep)
            g1 = g1 + 2.*g2
         case (3)
            g2 = g0 + 0.5 * g2 ! g2 is h*k2
            if (rk_step) call mb_communicate(g2)
            if (junction_sync_stage) call synchronise_zgrid_junctions(g2)
            call prepare_moving_k_stage(0.5)
            call budget_set_stage_weight(1.0 / 3.0)
            call add_explicit_gyrokinetic_terms(g2, g3, restart_time_step, istep)
            g1 = g1 + 2.*g3
         case (4)
            g3 = g0 + g3 ! g3 is h*k3
            if (rk_step) call mb_communicate(g3)
            if (junction_sync_stage) call synchronise_zgrid_junctions(g3)
            call prepare_moving_k_stage(1.0)
            call budget_set_stage_weight(1.0 / 6.0)
            call add_explicit_gyrokinetic_terms(g3, g, restart_time_step, istep)
            g1 = g1 + g
         end select
         
         ! If the code_dt is reset, we need to quit this loop and restart the timestep again
         ! Otherwise increase the <counter> that keeps track of the iterations
         if (restart_time_step) then
            counter = 10
         else
            counter = counter + 1
         end if
         
      end do

      ! This is g at intermediate time level
      g = g0 + g1 / 6.

   end subroutine advance_explicit_rk4

end module gyrokinetic_equation_explicit
