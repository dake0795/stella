!###############################################################################
!############ READ STELLA NAMELISTS FOR RADIAL VARIATION (MULTIBOX) ############
!###############################################################################
! 
! This module will read the namelists associated with shear flow:
!   
!   flow_shear
!     prp_shear_enabled = .false.
!     hammett_flow_shear = .true.
!     explicit_junction_sync = 'none'
!     flow_shear_explicit_remap = .false.
!     flow_shear_moving_k = .false.
!     g_exb = 0.0
!     g_exbfac = 1.0
!     omprimfac = 1.0
! 
! For each namelists two (or three) routines exist:
!    - set_default_parameters_<namelist>
!    - read_namelist_<namelist>
!    - check_inputs_<namelist>
! 
! First the default input parameters are set, then the default options are
! overwritten with those specified in the input file. Optionally, it is
! checked whether any input variables are clashing with each other.
! 
!###############################################################################
module namelist_flow_shear

   implicit none

   ! Make reading routines accesible to other modules
   public :: read_namelist_flow_shear

   private

   ! These variables are used in every single subroutine, so make them global
   integer :: in_file
   logical :: dexist

contains

   
   !****************************************************************************
   !                               FLOW SHEAR TERMS                            !
   !****************************************************************************
   subroutine read_namelist_flow_shear(prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, &
      omprimfac, explicit_junction_sync, flow_shear_explicit_remap, flow_shear_moving_k)

      use mp, only: proc0

      implicit none

      ! Variables that are read from the input file
      logical, intent(out) :: prp_shear_enabled, hammett_flow_shear 
      character(len=8), intent(out) :: explicit_junction_sync
      logical, intent(out) :: flow_shear_explicit_remap
      logical, intent(out) :: flow_shear_moving_k
      real :: g_exb, g_exbfac, omprimfac
      
      !-------------------------------------------------------------------------

      if (.not. proc0) return
      call set_default_parameters_flow_shear
      call read_input_file_flow_shear
      call write_parameters_to_input_file

   contains
      
      !------------------------ Default input parameters -----------------------
      subroutine set_default_parameters_flow_shear

         implicit none

         prp_shear_enabled = .false.
         hammett_flow_shear = .true.
         ! Explicit schemes: re-synchronise the duplicated twist-and-shift junction points
         ! of linked modes ('none' = no synchronisation, 'step' = after each accepted step,
         ! 'stage' = also before each Runge-Kutta stage). The implicit solve always does this.
         explicit_junction_sync = 'none'
         ! Fully explicit simulations: apply the perpendicular flow-shear remap after each
         ! accepted explicit step. Otherwise the remap is only applied in the implicit solve.
         flow_shear_explicit_remap = .false.
         ! Treat the perpendicular flow shear continuously in time: each mode keeps its
         ! label akx and its physical kx moves as kx - g_exb*ky*t (moving wavevector)
         ! instead of being remapped by one radial grid spacing at a time.
         flow_shear_moving_k = .false.
         g_exb = 0.0
         g_exbfac = 1.0
         omprimfac = 1.0

      end subroutine set_default_parameters_flow_shear

      !---------------------------- Read input file ----------------------------
      subroutine read_input_file_flow_shear

         use file_utils, only: input_unit_exist
         implicit none

         namelist /flow_shear/ prp_shear_enabled, hammett_flow_shear, g_exb, g_exbfac, omprimfac, &
            explicit_junction_sync, &
            flow_shear_explicit_remap, &
            flow_shear_moving_k
         in_file = input_unit_exist('flow_shear', dexist)
         if (dexist) read (unit=in_file, nml=flow_shear)

      end subroutine read_input_file_flow_shear
      
      !------------------------- Write input parameters ------------------------
      subroutine write_parameters_to_input_file

         use file_units, only: unit => unit_input_file_with_defaults

         implicit none

         !-------------------------------------------------------------------------

         write (unit, '(A)') '&flow_shear'
         write (unit, '(A, L1)') '  prp_shear_enabled = ', prp_shear_enabled
         write (unit, '(A, L1)') '  hammett_flow_shear = ', hammett_flow_shear
         write (unit, '(A, A, A)') '  explicit_junction_sync = "', trim(explicit_junction_sync), '"'
         write (unit, '(A, L1)') '  flow_shear_explicit_remap = ', flow_shear_explicit_remap
         write (unit, '(A, L1)') '  flow_shear_moving_k = ', flow_shear_moving_k
         write (unit, '(A, ES0.4)') '  g_exb = ', g_exb
         write (unit, '(A, ES0.4)') '  g_exbfac = ', g_exbfac
         write (unit, '(A, ES0.4)') '  omprimfac = ', omprimfac
         write (unit, '(A)') '/'
         write (unit, '(A)') ''

      end subroutine write_parameters_to_input_file

   end subroutine read_namelist_flow_shear
   
end module namelist_flow_shear
