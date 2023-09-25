!===================================================================================================
!  Created by LiXiaohan on 20/08/10
!
!  interface of the prognostic cloud macrophysics (Park scheme , Lin scheme)
!
!===================================================================================================

 module grist_macrop

     use grist_handle_error,                 only: endrun
     use grist_constants,                    only: r8, i4

     implicit none
     private
     save

     public          :: macrop_driver_init,        &
                        read_nml_macrop,           &
                        macrop_driver_tend,        &
                        end_of_macrop
 ! Private:
     character(len=16)  :: macrop_scheme   ! Macrophysics scheme

 contains

     subroutine read_nml_macrop(nlfile)

     use phys_control,                       only: phys_getopts
     use grist_macrop_park,                  only: read_nml_macrop_park
     use grist_macrop_lin,                   only: read_nml_macrop_lin

     character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

     call phys_getopts(macrop_scheme_out=macrop_scheme)

     select case (macrop_scheme)
     case ('park')
         call read_nml_macrop_park(nlfile)
     case ('lin')
         call read_nml_macrop_lin(nlfile)
     case default
        call endrun('read_nml_macrop:: unrecognized macrop_scheme')
     end select

     end subroutine read_nml_macrop


     subroutine macrop_driver_init(ncol)
     use grist_macrop_park,                  only: macrop_driver_init_park
     use grist_macrop_lin,                   only: macrop_driver_init_lin
! io
    integer, intent(in)  :: ncol

     select case (macrop_scheme)
     case ('park')
         call macrop_driver_init_park(ncol)
     case ('lin')
         call macrop_driver_init_lin(ncol)
     case default
        call endrun('read_nml_macrop:: unrecognized macrop_scheme')
     end select

     end subroutine macrop_driver_init


     subroutine end_of_macrop
     use grist_macrop_park,                  only: end_of_macrop_park
     use grist_macrop_lin,                   only: end_of_macrop_lin

     select case (macrop_scheme)
     case ('park')
         call end_of_macrop_park
     case ('lin')
         call end_of_macrop_lin
     case default
        call endrun('read_nml_macrop:: unrecognized macrop_scheme')
     end select

     end subroutine end_of_macrop


     subroutine macrop_driver_tend(ncol  , dtime , istep   ,     &
                                   dlf   , dlf2  ,  cmfmc2 ,     &
                                   det_s , det_ice)

     use grist_nml_module,                   only: nlev, nlevp
     use grist_macrop_park,                  only: macrop_driver_tend_park
     use grist_macrop_lin,                   only: macrop_driver_tend_lin

! io
     integer,  intent(in)  :: ncol
     integer,  intent(in)  :: istep
     real(r8), intent(in)  :: dtime                   ! Timestep
     real(r8), intent(in)  :: dlf(nlev, ncol)         ! Detrained water from convection schemes
     real(r8), intent(in)  :: dlf2(nlev, ncol)        ! Detrained water from shallow convection scheme
     real(r8), intent(in)  :: cmfmc2(nlevp, ncol)     ! Shallow convective mass flux [ kg/s/m^2 ]
 
     ! These two variables are needed for energy check    
     real(r8), intent(out) :: det_s(ncol)             ! Integral of detrained static energy from ice
     real(r8), intent(out) :: det_ice(ncol)           ! Integral of detrained ice for energy check

     select case (macrop_scheme)
     case ('park')
         call macrop_driver_tend_park(ncol  , dtime , istep   ,   &
                                    dlf   , dlf2  ,  cmfmc2 ,     &
                                    det_s , det_ice)
     case ('lin')
         call macrop_driver_tend_lin(ncol  , dtime , istep   ,    &
                                    dlf   , dlf2  ,  cmfmc2 ,     &
                                    det_s , det_ice)
     case default
        call endrun('read_nml_macrop:: unrecognized macrop_scheme')
     end select

     end subroutine macrop_driver_tend

 end module grist_macrop
