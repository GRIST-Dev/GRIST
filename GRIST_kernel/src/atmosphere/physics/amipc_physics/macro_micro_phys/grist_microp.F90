!===================================================================================================
!
!  Created by LiXiaohan on 19/07/10, adopted from CAM5
!
!     interface of the prognostic cloud macrophysics
!
!===================================================================================================

 module grist_microp

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8, i4
    use phys_control,                       only: phys_getopts

    implicit none
    private
    save

    public          :: microp_driver_init,        &
                       microp_driver_tend,        &
                       end_of_microp,             &
                       read_nml_microp

! Private:
    character(len=16)  :: microp_scheme   ! Microphysics scheme

 contains

    subroutine read_nml_microp(nlfile)

     use micro_mg,                           only: read_nml_micro_mg
     use micro_mg_lin,                       only: read_nml_micro_mg_lin

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input


    call phys_getopts(microp_scheme_out=microp_scheme)

    select case (microp_scheme)
    case ('MG')
       call read_nml_micro_mg(nlfile)
    case ('MG_LIN')
       call read_nml_micro_mg_lin(nlfile)
    case default
       call endrun('read_nml_microp:: unrecognized microp_scheme')
    end select

    end subroutine read_nml_microp


    subroutine microp_driver_init(ncol)
    use grist_nml_module,               only: nlev, nlevp, ntracer
    use cldwat2m_macro,                 only: ini_macro
    use cldwat2m_macro_lin,             only: ini_macro_lin
    use micro_mg,                       only: micro_mg_init 
    use micro_mg_lin,                   only: micro_mg_init_lin 

! io
    integer, intent(in)  :: ncol
! local
    character(len=16)  :: macrop_scheme   ! Macrophysics scheme

    call phys_getopts(macrop_scheme_out=macrop_scheme)

    select case (macrop_scheme)
    case ('park')
       call ini_macro
    case ('lin')
       call ini_macro_lin 
    case default
       call endrun('read_nml_macrop:: unrecognized macrop_scheme')
    end select

    select case (microp_scheme)
    case ('MG')
       call micro_mg_init(ncol)
    case ('MG_LIN')
       call micro_mg_init_lin(ncol)
    case default
       call endrun('microp_driver_init:: unrecognized microp_scheme')
    end select

    end subroutine microp_driver_init


    subroutine end_of_microp

    use micro_mg,                       only: micro_mg_end 
    use micro_mg_lin,                   only: micro_mg_end_lin 

    select case (microp_scheme)
    case ('MG')
       call micro_mg_end
    case ('MG_LIN')
       call micro_mg_end_lin
    case default
        call endrun("end_of_microp:: unrecognized microp_scheme")
    end select

    end subroutine end_of_microp


    subroutine microp_driver_tend(ncol, dtime)

    use micro_mg,       only: micro_mg_tend
    use micro_mg_lin,       only: micro_mg_tend_lin
! io
    integer , intent(in)  :: ncol
    real(r8), intent(in)  :: dtime                    ! Timestep

    select case (microp_scheme)
    case ('MG')
       call micro_mg_tend(ncol, dtime)
    case ('MG_LIN')
       call micro_mg_tend_lin(ncol, dtime)
    case default
       call endrun('microp_driver_tend: unrecognized microp_scheme')
    end select

    end subroutine microp_driver_tend


 end module grist_microp
