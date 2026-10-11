!###############################################################################
!######################## FREE-ENERGY BUDGET DIAGNOSTIC #########################
!###############################################################################
!
! Optional diagnostic (<energy_budget_diagnostics> in the flow_shear namelist)
! for the fully explicit rk4 time advance, used with the moving-wavevector flow
! shear and with the explicit flow-shear remap. It measures the change of the
! free energy
!     W = 1/2 sum m zw H |g|^2 + 1/2 sum m zw P |phi|^2 + sum m zw K |apar|^2 / beta,
! with H = n_s T_s w_vpa w_mu / F0_s, zw the normalised field-line weights
! dl/B (upper z end point excluded), m = 1 for ky = 0 and 2 for ky > 0,
! P = denominator_fields and K = kperp2 at the current metric, and attributes it
!   - to every term N of the explicit right-hand side, through
!     dW/dt|_N = sum m zw H Re[conj(h) N], h = g + (Z/T) F0 J0 phi, accumulated
!     with the rk4 weights of the stages;
!   - to the change of the metric at fixed canonical distribution (moving
!     wavevector), dW/dt = Re[conj(phi) rdot] - 1/2 Pdot |phi|^2
!     + [(Kdot + Ddot) |apar|^2 - 2 Re[conj(apar) cdot]] / beta;
!   - to every projection (relabelling, remap, junction synchronisation and the
!     implicit mirror and streaming split steps), as W(after) - W(before).
! Energies are cumulative since the start of the simulation (or the restart),
! and W - W0 = sum(E_*) + R_RK, with R_RK the remainder (time discretisation).
! The diagnostic only reads g and the fields: the fields used for W are solved
! into local arrays and <fields_updated> is restored.
!
! Every accumulator holds the local kinetic part on every processor and the
! field part on proc0 only; they are summed over processors when written.
!
!###############################################################################
module diagnostics_energy_budget

   implicit none

   ! Make routines available to other modules
   public :: init_energy_budget, write_energy_budget
   public :: budget_start_step, budget_incoming_sync, budget_commit_step, budget_end_step
   public :: budget_set_stage_weight, budget_shear_stage
   public :: budget_stage_start, budget_add_term, budget_stage_end
   public :: budget_relabel_start, budget_relabel_end
   public :: budget_sync_start, budget_sync_end
   public :: budget_remap_start, budget_remap_end
   public :: budget_split_start, budget_split_end, budget_split_sync_end
   public :: energy_budget_on

   ! Contributions
   public :: ib_exb, ib_prl, ib_mirror, ib_wdy, ib_wdx, ib_wstar, ib_stream, ib_hyper, ib_other
   public :: ib_mi_split, ib_st_split

   private

   integer, parameter :: nbud = 18
   integer, parameter :: ib_exb = 1, ib_prl = 2, ib_mirror = 3, ib_wdy = 4, ib_wdx = 5, ib_wstar = 6, &
      ib_stream = 7, ib_hyper = 8, ib_other = 9, ib_shear = 10, ib_rl_out = 11, ib_rl_err = 12, &
      ib_rm_jump = 13, ib_sync = 14, ib_mi_split = 15, ib_sync_in = 16, ib_rm_out_info = 17, ib_st_split = 18
   character(len=16), parameter :: bud_name(nbud) = [character(len=16) :: 'E_exb', 'E_prl_shear', 'E_mirror_rhs', &
      'E_wdrifty', 'E_wdriftx', 'E_wstar', 'E_stream', 'E_hyper', 'E_other_rhs', 'E_shear_metric', 'E_relabel_out', &
      'E_relabel_reterr', 'E_remap_jump', 'E_sync_step', 'E_mirror_split', 'E_sync_incoming', 'I_remap_edge', &
      'E_stream_split']

   logical :: energy_budget_on = .false.
   logical :: bud_started = .false., bud_wlast_ok = .false., bud_relabel_any = .false., bud_remap = .false.
   real :: bud_step(nbud) = 0.0, bud_cum(nbud) = 0.0, bud_b = 0.0, bud_sync_in_maxchange = 0.0, bud_s = 0.0
   real, dimension(6) :: bud_W0 = 0.0, bud_wlast = 0.0, bud_wstep_start = 0.0
   real, dimension(6) :: bud_wsync0 = 0.0, bud_wpre = 0.0, bud_wout = 0.0, bud_wsplit0 = 0.0
   integer :: bud_nacc = 0, bud_nrelabel = 0, bud_nremap = 0, bud_unit = -1
   real, dimension(:, :), allocatable :: bud_H
   real, dimension(:), allocatable :: bud_zw, bud_m
   complex, dimension(:, :, :, :, :), allocatable :: bud_hstage

contains

   !****************************************************************************
   !                               Initialisation
   !****************************************************************************
   subroutine init_energy_budget

      use mp, only: proc0, mp_abort
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, aky
      use grids_species, only: spec, has_electron_species
      use grids_velocity, only: wgts_vpa, wgts_mu, maxwell_vpa, maxwell_mu, maxwell_fac
      use geometry, only: dl_over_b
      use parameters_physics, only: include_bpar, full_flux_surface, radial_variation, fphi
      use parameters_numerical, only: explicit_algorithm_switch, explicit_algorithm_rk4, fully_explicit
      use gk_flow_shear, only: energy_budget_diagnostics
      use file_utils, only: run_name

      implicit none

      integer :: ivmu, iv, imu, is, iz, iky

      !-------------------------------------------------------------------------

      energy_budget_on = energy_budget_diagnostics
      if (.not. energy_budget_on) return

      if (.not. fully_explicit) call mp_abort('energy_budget_diagnostics requires a fully explicit time advance')
      if (explicit_algorithm_switch /= explicit_algorithm_rk4) call mp_abort('energy_budget_diagnostics requires rk4')
      if (include_bpar .or. full_flux_surface .or. radial_variation .or. ntubes /= 1 .or. abs(fphi - 1.0) > epsilon(0.)) &
         call mp_abort('energy_budget_diagnostics requires a local flux tube with ntubes = 1, fphi = 1 and no bpar')
      if (.not. has_electron_species(spec)) call mp_abort('energy_budget_diagnostics requires kinetic electrons')

      allocate (bud_zw(-nzgrid:nzgrid), bud_m(naky))
      bud_zw = dl_over_b(1, :) / sum(dl_over_b(1, :))
      do iky = 1, naky
         bud_m(iky) = merge(1.0, 2.0, aky(iky) == 0.0)
      end do
      allocate (bud_H(-nzgrid:nzgrid, vmu_lo%llim_proc:vmu_lo%ulim_alloc))
      bud_H = 0.0
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do iz = -nzgrid, nzgrid
            bud_H(iz, ivmu) = bud_zw(iz) * spec(is)%dens * spec(is)%temp * wgts_vpa(iv) * wgts_mu(1, iz, imu) &
               / (maxwell_vpa(iv, is) * maxwell_mu(1, iz, imu, is) * maxwell_fac(is))
         end do
      end do
      bud_cum = 0.0
      bud_step = 0.0

      if (proc0) then
         open (newunit=bud_unit, file=trim(run_name)//'.budget', status='unknown', position='append', action='write')
         write (bud_unit, '(A)') '# Free-energy budget. Energies are cumulative since the start of this run (W0 = W of the'
         write (bud_unit, '(A)') '# initial state). W - W0 = sum(E_*) + R_RK; I_remap_edge is for information only (already'
         write (bud_unit, '(A)') '# included in E_remap_jump). E_* of explicit terms are rk4-weighted stage powers; projections'
         write (bud_unit, '(A)') '# are W(after) - W(before). Junction synchronisation of the stages is included in R_RK.'
         write (bud_unit, '(A)', advance='no') '# 1:t 2:istep 3:code_dt 4:W 5:W_kin 6:W_phi 7:W_apar 8:W_nonzonal 9:W0'
         do iv = 1, nbud
            write (bud_unit, '(1X,I0,A,A)', advance='no') 9 + iv, ':', trim(bud_name(iv))
         end do
         write (bud_unit, '(1X,I0,A,1X,I0,A,1X,I0,A,1X,I0,A,1X,I0,A,1X,I0,A)') 10 + nbud, ':R_RK', 11 + nbud, ':n_accepted', &
            12 + nbud, ':n_relabel_events', 13 + nbud, ':n_remap_events', 14 + nbud, ':max_incoming_sync_change', &
            15 + nbud, ':max_abs_shift_state'
         flush (bud_unit)
      end if

   end subroutine init_energy_budget

   !****************************************************************************
   !                       Free energy of a distribution
   !****************************************************************************
   ! W of g (dist = 'g') or gbar (dist = 'gbar') at the current metric, split as
   ! w = [kinetic zonal, kinetic non-zonal, phi zonal, phi non-zonal, apar zonal,
   ! apar non-zonal]; the kinetic parts are local, the field parts on proc0 only.
   ! With mask(naky, nakx), wmask is the same quantity restricted to the masked modes.
   !****************************************************************************
   subroutine budget_W(g, dist, w, mask, wmask)

      use mp, only: proc0
      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid, ntubes
      use grids_kxky, only: naky, nakx
      use arrays, only: kperp2, denominator_fields
      use parameters_physics, only: beta, include_apar
      use field_equations, only: advance_fields, fields_updated
      use calculations_tofrom_ghf, only: gbar_to_g

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      character(*), intent(in) :: dist
      real, dimension(6), intent(out) :: w
      logical, dimension(:, :), intent(in), optional :: mask
      real, dimension(6), intent(out), optional :: wmask

      complex, dimension(:, :, :, :), allocatable :: phl, apl, bpl
      complex, dimension(:, :, :, :, :), allocatable :: gg
      logical :: saved

      !-------------------------------------------------------------------------

      allocate (phl(naky, nakx, -nzgrid:nzgrid, ntubes), apl(naky, nakx, -nzgrid:nzgrid, ntubes), &
         bpl(naky, nakx, -nzgrid:nzgrid, ntubes))
      saved = fields_updated
      fields_updated = .false.
      if (include_apar) then
         call advance_fields(g, phl, apl, bpl, dist=dist)
      else
         call advance_fields(g, phl, apl, bpl, dist='g')
      end if
      fields_updated = saved
      if (dist == 'gbar' .and. include_apar) then
         allocate (gg, source=g)
         call gbar_to_g(gg, apl, 1.0)
         call add_parts(gg)
         deallocate (gg)
      else
         call add_parts(g)
      end if
      deallocate (phl, apl, bpl)

   contains

      subroutine add_parts(gp)

         complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gp
         integer :: ivmu, iz, ikx, iky, iw
         real :: e, wt

         w = 0.0
         if (present(wmask)) wmask = 0.0
         do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
            do iz = -nzgrid, nzgrid
               if (bud_H(iz, ivmu) == 0.0) cycle
               do ikx = 1, nakx
                  do iky = 1, naky
                     e = 0.5 * bud_m(iky) * bud_H(iz, ivmu) * abs(gp(iky, ikx, iz, 1, ivmu))**2
                     iw = merge(1, 2, bud_m(iky) == 1.0)
                     w(iw) = w(iw) + e
                     if (present(mask)) then
                        if (mask(iky, ikx)) wmask(iw) = wmask(iw) + e
                     end if
                  end do
               end do
            end do
         end do
         if (proc0) then
            do iz = -nzgrid, nzgrid
               do ikx = 1, nakx
                  do iky = 1, naky
                     wt = bud_m(iky) * bud_zw(iz)
                     iw = merge(0, 1, bud_m(iky) == 1.0)
                     e = 0.5 * wt * denominator_fields(iky, ikx, iz) * abs(phl(iky, ikx, iz, 1))**2
                     w(3 + iw) = w(3 + iw) + e
                     if (present(mask)) then
                        if (mask(iky, ikx)) wmask(3 + iw) = wmask(3 + iw) + e
                     end if
                     if (include_apar .and. beta > 0.0) then
                        e = wt * kperp2(iky, ikx, 1, iz) * abs(apl(iky, ikx, iz, 1))**2 / beta
                        w(5 + iw) = w(5 + iw) + e
                        if (present(mask)) then
                           if (mask(iky, ikx)) wmask(5 + iw) = wmask(5 + iw) + e
                        end if
                     end if
                  end do
               end do
            end do
         end if

      end subroutine add_parts

   end subroutine budget_W

   !****************************************************************************
   !                       Bookkeeping of a time step
   !****************************************************************************
   ! Every attempt (including a restarted time step) starts from the free energy
   ! of the last accepted state; W0 is the free energy of the very first state.
   subroutine budget_start_step(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      if (.not. energy_budget_on) return
      bud_step = 0.0
      if (.not. bud_started) then
         call budget_W(g, 'g', bud_W0)
         bud_wstep_start = bud_W0
         bud_started = .true.
      end if
      bud_wlast = bud_wstep_start
      bud_wlast_ok = .true.

   end subroutine budget_start_step

   ! Junction synchronisation of the incoming state (gpre -> g, canonical gbar);
   ! measured only if any value changed
   subroutine budget_incoming_sync(gpre, g)

      use mp, only: max_allreduce
      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: gpre, g
      real :: chg
      real, dimension(6) :: wb, wa

      chg = 0.0
      if (vmu_lo%ulim_proc >= vmu_lo%llim_proc) chg = maxval(abs(g(:, :, :, :, vmu_lo%llim_proc:vmu_lo%ulim_proc) &
         - gpre(:, :, :, :, vmu_lo%llim_proc:vmu_lo%ulim_proc)))
      call max_allreduce(chg)
      bud_sync_in_maxchange = max(bud_sync_in_maxchange, chg)
      if (chg > 0.0) then
         call budget_W(gpre, 'gbar', wb)
         call budget_W(g, 'gbar', wa)
         bud_step(ib_sync_in) = bud_step(ib_sync_in) + sum(wa) - sum(wb)
         bud_wlast = bud_wlast + (wa - wb)
      end if

   end subroutine budget_incoming_sync

   ! Commit the stage powers of an accepted step
   subroutine budget_commit_step

      if (.not. energy_budget_on) return
      bud_cum = bud_cum + bud_step
      bud_nacc = bud_nacc + 1
      bud_wlast_ok = .false.

   end subroutine budget_commit_step

   ! Free energy of the accepted state at the end of the step
   subroutine budget_end_step(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      if (.not. energy_budget_on) return
      if (.not. bud_wlast_ok) call budget_W(g, 'g', bud_wlast)
      bud_wstep_start = bud_wlast

   end subroutine budget_end_step

   !****************************************************************************
   !                     Terms of the explicit right-hand side
   !****************************************************************************
   subroutine budget_set_stage_weight(b)

      real, intent(in) :: b

      bud_b = b

   end subroutine budget_set_stage_weight

   ! h = g + (Z/T) F0 J0 phi of the stage, for the term powers
   subroutine budget_stage_start(pdf)

      use arrays_fields, only: phi, bpar
      use calculations_tofrom_ghf, only: g_to_h
      use parameters_physics, only: fphi
      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: pdf

      if (.not. energy_budget_on) return
      allocate (bud_hstage, source=pdf)
      call g_to_h(bud_hstage, phi, bpar, fphi)
      bud_s = 0.0

   end subroutine budget_stage_start

   ! Add the rk4-weighted power of the term just added to rhs (s = <h, rhs>)
   subroutine budget_add_term(ib, rhs)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx

      implicit none

      integer, intent(in) :: ib
      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: rhs
      integer :: ivmu, iz, ikx, iky
      real :: s

      if (.not. energy_budget_on) return
      s = 0.0
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         do iz = -nzgrid, nzgrid
            if (bud_H(iz, ivmu) == 0.0) cycle
            do ikx = 1, nakx
               do iky = 1, naky
                  s = s + bud_m(iky) * bud_H(iz, ivmu) * real(conjg(bud_hstage(iky, ikx, iz, 1, ivmu)) * rhs(iky, ikx, iz, 1, ivmu))
               end do
            end do
         end do
      end do
      bud_step(ib) = bud_step(ib) + bud_b * (s - bud_s)
      bud_s = s

   end subroutine budget_add_term

   subroutine budget_stage_end

      if (allocated(bud_hstage)) deallocate (bud_hstage)

   end subroutine budget_stage_end

   !****************************************************************************
   !        Work done by the moving metric at fixed canonical distribution
   !****************************************************************************
   ! For the current Runge-Kutta stage; phi and apar are the fields solved from
   ! x (dist = 'gbar') at the stage metric. Adds b * code_dt * power.
   !****************************************************************************
   subroutine budget_shear_stage(x)

      use mp, only: proc0, sum_allreduce
      use parallelisation_layouts, only: vmu_lo, iv_idx, imu_idx, is_idx
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx, aky
      use grids_time, only: code_dt
      use grids_species, only: spec
      use grids_velocity, only: wgts_vpa, wgts_mu, maxwell_vpa, maxwell_mu, maxwell_fac, vpa, vperp2
      use grids_extended_zgrid, only: neigen, nsegments, ikxmod
      use geometry, only: bmag, gradx_dot_gradx, gradx_dot_grady
      use arrays, only: flow_shear_moving_k, kx_phys
      use arrays_fields, only: phi, apar
      use arrays_gyro_averages, only: aj0x, aj1x
      use parameters_physics, only: beta, include_apar
      use gk_flow_shear, only: g_exb, g_exbfac

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: x

      real, dimension(:, :, :), allocatable :: kdot, pdot, ddot
      complex, dimension(:, :, :), allocatable :: rdx, cdx
      integer :: ivmu, iv, imu, is, iz, ikx, iky, ie, iseg
      real :: gam, fac, wv, j0d, f0, p, tmp

      !-------------------------------------------------------------------------

      if (.not. (energy_budget_on .and. flow_shear_moving_k)) return
      gam = g_exb * g_exbfac
      if (abs(gam) <= epsilon(0.)) return

      allocate (kdot(naky, nakx, -nzgrid:nzgrid), pdot(naky, nakx, -nzgrid:nzgrid), ddot(naky, nakx, -nzgrid:nzgrid))
      allocate (rdx(naky, nakx, -nzgrid:nzgrid), cdx(naky, nakx, -nzgrid:nzgrid))

      ! d(kperp2)/dt, averaged over connected points as for kperp2
      do iz = -nzgrid, nzgrid
         do ikx = 1, nakx
            do iky = 1, naky
               kdot(iky, ikx, iz) = -gam * aky(iky) * (2.0 * kx_phys(iky, ikx) * gradx_dot_gradx(1, iz) &
                  + 2.0 * aky(iky) * gradx_dot_grady(1, iz))
            end do
         end do
      end do
      do iky = 1, naky
         do ie = 1, neigen(iky)
            if (nsegments(ie, iky) > 1) then
               do iseg = 2, nsegments(ie, iky)
                  tmp = 0.5 * (kdot(iky, ikxmod(iseg - 1, ie, iky), nzgrid) + kdot(iky, ikxmod(iseg, ie, iky), -nzgrid))
                  kdot(iky, ikxmod(iseg, ie, iky), -nzgrid) = tmp
                  kdot(iky, ikxmod(iseg - 1, ie, iky), nzgrid) = tmp
               end do
            end if
         end do
      end do

      ! dJ0/dt = -(J1(a)/a) a adot, with a^2 = (bess_fac smz / B)^2 vperp2 kperp2
      rdx = 0.0
      cdx = 0.0
      pdot = 0.0
      ddot = 0.0
      do ivmu = vmu_lo%llim_proc, vmu_lo%ulim_proc
         iv = iv_idx(vmu_lo, ivmu)
         imu = imu_idx(vmu_lo, ivmu)
         is = is_idx(vmu_lo, ivmu)
         do iz = -nzgrid, nzgrid
            fac = 0.5 * (spec(is)%bess_fac * spec(is)%smz_psi0 / bmag(1, iz))**2 * vperp2(1, iz, imu)
            wv = wgts_mu(1, iz, imu) * wgts_vpa(iv)
            f0 = maxwell_vpa(iv, is) * maxwell_mu(1, iz, imu, is) * maxwell_fac(is)
            do ikx = 1, nakx
               do iky = 1, naky
                  if (aky(iky) == 0.0) cycle
                  j0d = -aj1x(iky, ikx, iz, ivmu) * fac * kdot(iky, ikx, iz)
                  rdx(iky, ikx, iz) = rdx(iky, ikx, iz) + wv * spec(is)%z * spec(is)%dens_psi0 * j0d * x(iky, ikx, iz, 1, ivmu)
                  cdx(iky, ikx, iz) = cdx(iky, ikx, iz) + wv * beta * spec(is)%z * spec(is)%dens_psi0 * spec(is)%stm_psi0 &
                     * vpa(iv) * j0d * x(iky, ikx, iz, 1, ivmu)
                  pdot(iky, ikx, iz) = pdot(iky, ikx, iz) - wv * spec(is)%z**2 * spec(is)%dens_psi0 / spec(is)%temp &
                     * 2.0 * aj0x(iky, ikx, iz, ivmu) * j0d * f0
                  ddot(iky, ikx, iz) = ddot(iky, ikx, iz) + wv * 2.0 * beta * spec(is)%z**2 * spec(is)%dens / spec(is)%mass &
                     * vpa(iv)**2 * f0 * 2.0 * aj0x(iky, ikx, iz, ivmu) * j0d
               end do
            end do
         end do
      end do
      call sum_allreduce(rdx)
      call sum_allreduce(cdx)
      call sum_allreduce(pdot)
      call sum_allreduce(ddot)

      p = 0.0
      if (proc0) then
         do iz = -nzgrid, nzgrid
            do ikx = 1, nakx
               do iky = 1, naky
                  if (aky(iky) == 0.0) cycle
                  tmp = real(conjg(phi(iky, ikx, iz, 1)) * rdx(iky, ikx, iz)) &
                     - 0.5 * pdot(iky, ikx, iz) * abs(phi(iky, ikx, iz, 1))**2
                  if (include_apar .and. beta > 0.0) tmp = tmp + ((kdot(iky, ikx, iz) + ddot(iky, ikx, iz)) &
                     * abs(apar(iky, ikx, iz, 1))**2 - 2.0 * real(conjg(apar(iky, ikx, iz, 1)) * cdx(iky, ikx, iz))) / beta
                  p = p + bud_m(iky) * bud_zw(iz) * tmp
               end do
            end do
         end do
         bud_step(ib_shear) = bud_step(ib_shear) + bud_b * code_dt * p
      end if

      deallocate (kdot, pdot, ddot, rdx, cdx)

   end subroutine budget_shear_stage

   !****************************************************************************
   !                                 Projections
   !****************************************************************************
   ! Moving-wavevector relabelling of the canonical gbar: before the relabelling,
   ! W at the stage c = 1 metric and the modal W of the labels that leave the grid
   ! (mask); after it, W at the new committed metric
   subroutine budget_relabel_start(g, mask)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      logical, dimension(:, :), intent(in) :: mask

      bud_relabel_any = energy_budget_on .and. any(mask)
      if (bud_relabel_any) call budget_W(g, 'gbar', bud_wpre, mask, bud_wout)

   end subroutine budget_relabel_start

   subroutine budget_relabel_end(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      real, dimension(6) :: wpost

      if (.not. bud_relabel_any) return
      call budget_W(g, 'gbar', wpost)
      bud_cum(ib_rl_out) = bud_cum(ib_rl_out) - sum(bud_wout)
      bud_cum(ib_rl_err) = bud_cum(ib_rl_err) + sum(wpost) - sum(bud_wpre) + sum(bud_wout)
      bud_wlast = wpost
      bud_wlast_ok = .true.
      bud_nrelabel = bud_nrelabel + 1
      bud_relabel_any = .false.

   end subroutine budget_relabel_end

   ! Junction synchronisation and periodic end-point copy of the accepted state
   subroutine budget_sync_start(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      if (.not. energy_budget_on) return
      if (.not. bud_wlast_ok) call budget_W(g, 'gbar', bud_wlast)
      bud_wsync0 = bud_wlast

   end subroutine budget_sync_start

   subroutine budget_sync_end(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      if (.not. energy_budget_on) return
      call budget_W(g, 'g', bud_wlast)
      bud_wlast_ok = .true.
      bud_cum(ib_sync) = bud_cum(ib_sync) + sum(bud_wlast) - sum(bud_wsync0)

   end subroutine budget_sync_end

   ! Flow-shear remap of the explicit time advance; the modal W of the vacated
   ! edge label of each shifted row is recorded for information
   subroutine budget_remap_start(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid
      use grids_kxky, only: naky, nakx, ikx_max, zonal_mode
      use arrays, only: shift_state
      use gk_flow_shear, only: perp_flow_shear_pending, hammett_flow_shear, shift_times, shift_sign

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      logical, dimension(:, :), allocatable :: mask
      integer :: iky

      bud_remap = .false.
      if (.not. energy_budget_on) return
      bud_remap = perp_flow_shear_pending()
      if (.not. bud_remap) return
      allocate (mask(naky, nakx))
      mask = .false.
      if (hammett_flow_shear) then
         do iky = 1, naky
            if (zonal_mode(iky)) cycle
            if (shift_state(iky) > 0.5 * shift_times(iky)) mask(iky, merge(ikx_max + 1, ikx_max, shift_sign < 0)) = .true.
         end do
      end if
      call budget_W(g, 'g', bud_wpre, mask, bud_wout)
      deallocate (mask)

   end subroutine budget_remap_start

   subroutine budget_remap_end(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      if (.not. bud_remap) return
      call budget_W(g, 'g', bud_wlast)
      bud_wlast_ok = .true.
      bud_cum(ib_rm_jump) = bud_cum(ib_rm_jump) + sum(bud_wlast) - sum(bud_wpre)
      bud_cum(ib_rm_out_info) = bud_cum(ib_rm_out_info) + sum(bud_wout)
      bud_nremap = bud_nremap + 1

   end subroutine budget_remap_end

   ! Implicit mirror or streaming split step (ib = ib_mi_split or ib_st_split).
   ! For the mirror, the subsequent end-point synchronisation is measured separately.
   subroutine budget_split_start(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      if (.not. energy_budget_on) return
      if (.not. bud_wlast_ok) call budget_W(g, 'g', bud_wlast)
      bud_wsplit0 = bud_wlast

   end subroutine budget_split_start

   subroutine budget_split_end(g, ib)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g
      integer, intent(in) :: ib

      if (.not. energy_budget_on) return
      call budget_W(g, 'g', bud_wlast)
      bud_wlast_ok = .true.
      bud_cum(ib) = bud_cum(ib) + sum(bud_wlast) - sum(bud_wsplit0)
      bud_wsplit0 = bud_wlast

   end subroutine budget_split_end

   ! End-point synchronisation after the mirror split step
   subroutine budget_split_sync_end(g)

      use parallelisation_layouts, only: vmu_lo
      use grids_z, only: nzgrid

      implicit none

      complex, dimension(:, :, -nzgrid:, :, vmu_lo%llim_proc:), intent(in) :: g

      if (.not. energy_budget_on) return
      call budget_W(g, 'g', bud_wlast)
      bud_wlast_ok = .true.
      bud_cum(ib_sync) = bud_cum(ib_sync) + sum(bud_wlast) - sum(bud_wsplit0)

   end subroutine budget_split_sync_end

   !****************************************************************************
   !                    Write one line to <run_name>.budget
   !****************************************************************************
   subroutine write_energy_budget(istep)

      use mp, only: proc0, sum_allreduce, max_allreduce
      use grids_time, only: code_time, code_dt
      use arrays, only: shift_state

      implicit none

      integer, intent(in) :: istep
      real, dimension(12 + nbud) :: v
      real :: wn, w0, closure, mx
      integer :: i

      !-------------------------------------------------------------------------

      if (.not. energy_budget_on) return
      v(1:6) = bud_wstep_start
      v(7:12) = bud_W0
      v(13:12 + nbud) = bud_cum
      call sum_allreduce(v)
      mx = bud_sync_in_maxchange
      call max_allreduce(mx)
      if (.not. proc0) return
      wn = sum(v(1:6))
      w0 = sum(v(7:12))
      closure = (wn - w0)
      do i = 1, nbud
         if (i /= ib_rm_out_info) closure = closure - v(12 + i)
      end do
      write (bud_unit, '(ES24.16,1X,I0,7(1X,ES24.16))', advance='no') code_time, istep, code_dt, wn, v(1) + v(2), &
         v(3) + v(4), v(5) + v(6), v(2) + v(4) + v(6), w0
      do i = 1, nbud
         write (bud_unit, '(1X,ES24.16)', advance='no') v(12 + i)
      end do
      write (bud_unit, '(1X,ES24.16,3(1X,I0),2(1X,ES24.16))') closure, bud_nacc, bud_nrelabel, bud_nremap, mx, &
         maxval(abs(shift_state))
      flush (bud_unit)

   end subroutine write_energy_budget

end module diagnostics_energy_budget
