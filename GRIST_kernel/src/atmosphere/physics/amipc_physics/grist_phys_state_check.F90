!======================================================
!
!  Created by LiXiaohan on 20/7/25.
!  
!  check physics variables
!
!======================================================

module grist_phys_state_check

    use grist_constants,                only: i4, r8, cp, rad2deg
    use grist_nml_module,               only: nlev, nlevp
    use grist_PhysC_data_structure,     only: pstate_cam
    use grist_physics_data_structure,   only: pstate
    use phys_control,                   only: phys_getopts
    use grist_handle_error,             only: endrun
    use grist_mpi

    implicit none
    private
    public :: phys_state_check
          
    logical, parameter :: print_check = .true.

    interface check_non
        module procedure check_nan_1d
        module procedure check_nan_2d
        module procedure check_nan_3d
    end interface

contains
    
    subroutine phys_state_check(name, ncol, lat, lon)
! io
    integer,             intent(in)           :: ncol  
    real(r8),            intent(in)           :: lat(ncol)             ! longitude in radian
    real(r8),            intent(in)           :: lon(ncol)             ! latitude in radian  
    character(len=*),    intent(in), optional :: name
! local
    logical :: state_debug_checks
    character(len=1080) :: msg
    integer :: m

    call phys_getopts(state_debug_checks_out=state_debug_checks)
    if(.not.state_debug_checks)return
 
    ! Check for NaN first to avoid any IEEE exceptions.
    if(present(name))then
        msg = 'NaN produced in physics_state by package '//trim(name)
    else
        msg = 'NaN produced in physics_state'
    end if

    call check_non(pstate%u_wind_at_pc_full_level%f, ncol, 'pstate%u_wind_at_pc_full_level', msg,lat,lon)
    call check_non(pstate%v_wind_at_pc_full_level%f, ncol, 'pstate%v_wind_at_pc_full_level', msg,lat,lon)
    call check_non(pstate%omega_at_pc_full_level%f, ncol, 'pstate%omega_at_pc_full_level', msg,lat,lon)
    call check_non(pstate%temp_at_pc_full_level%f, ncol, 'pstate%temp_at_pc_full_level', msg,lat,lon, lt=100._r8)
    call check_non(pstate%tracer_mxrt_at_pc_full_level%f, ncol, 'pstate%tracer_mxrt_at_pc_full_level', msg,lat,lon)
    call check_non(pstate%z_at_pc_full_level%f, ncol, 'pstate%z_at_pc_full_level', msg,lat,lon, lt=0._r8)
    call check_non(pstate%pressure_at_pc_surface%f, ncol, 'pstate%pressure_at_pc_surface', msg,lat,lon, lt=0._r8)
    call check_non(pstate%pressure_at_pc_full_level%f, ncol, 'pstate%pressure_at_pc_full_level', msg,lat,lon, lt=0._r8)
    call check_non(pstate%ts_at_pc_surface%f, ncol, 'pstate%ts_at_pc_surface', msg,lat,lon, lt=100._r8)
    call check_non(pstate%ustar_at_pc_surface%f, ncol, 'pstate%ustar_at_pc_surface', msg,lat,lon)
    call check_non(pstate%atm_in_shflx_at_pc_surface%f, ncol, 'pstate%atm_in_shflx_at_pc_surface', msg,lat,lon)
    call check_non(pstate%atm_in_qflx_at_pc_surface%f, ncol, 'pstate%atm_in_qflx_at_pc_surface', msg,lat,lon)
    call check_non(pstate%atm_in_lwup_at_pc_surface%f, ncol, 'pstate%atm_in_lwup_at_pc_surface', msg,lat,lon)
    call check_non(pstate%scalar_prect_surface%f, ncol, 'pstate%scalar_prect_surface', msg,lat,lon)


    if(trim(name) .eq. 'vertical diffusion')then
    if(allocated(pstate_cam%pbl_tke_at_pc_face_level%f))then
        call check_non(pstate_cam%pbl_tke_at_pc_face_level%f, ncol, 'pstate_cam%pbl_tke_at_pc_face_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%pbl_kvm_at_pc_face_level%f))then
        call check_non(pstate_cam%pbl_kvm_at_pc_face_level%f, ncol, 'pstate_cam%pbl_kvm_at_pc_face_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%pbl_kvh_at_pc_face_level%f))then
        call check_non(pstate_cam%pbl_kvh_at_pc_face_level%f, ncol, 'pstate_cam%pbl_kvh_at_pc_face_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%pbl_kvt_at_pc_face_level%f))then
        call check_non(pstate_cam%pbl_kvt_at_pc_face_level%f, ncol, 'pstate_cam%pbl_kvt_at_pc_face_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%pbl_pblh_at_pc_surface%f))then
        call check_non(pstate_cam%pbl_pblh_at_pc_surface%f, ncol, 'pstate_cam%pbl_pblh_at_pc_surface', msg,lat,lon)
    end if
    end if

    if(trim(name) .eq. 'shallow convection' .or. trim(msg) .eq. 'deep convection')then
    if(allocated(pstate_cam%sh_shfrc_at_pc_full_level%f))then
        call check_non(pstate_cam%sh_shfrc_at_pc_full_level%f, ncol, 'pstate_cam%sh_shfrc_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%sh_icwmr_at_pc_full_level%f))then
        call check_non(pstate_cam%sh_icwmr_at_pc_full_level%f, ncol, 'pstate_cam%sh_icwmr_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%dp_icwmr_at_pc_full_level%f))then
        call check_non(pstate_cam%dp_icwmr_at_pc_full_level%f, ncol, 'pstate_cam%dp_icwmr_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%updraft_mass_flux%f))then
        call check_non(pstate_cam%updraft_mass_flux%f, ncol, 'pstate_cam%updraft_mass_flux', msg,lat,lon)
    end if
    if(allocated(pstate_cam%cld_sh_frac_at_pc_full_level%f))then
        call check_non(pstate_cam%cld_sh_frac_at_pc_full_level%f, ncol, 'pstate_cam%cld_sh_frac_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%cld_dp_frac_at_pc_full_level%f))then
        call check_non(pstate_cam%cld_dp_frac_at_pc_full_level%f, ncol, 'pstate_cam%cld_dp_frac_at_pc_full_level', msg,lat,lon)
    end if
    end if

    if(trim(name) .eq. 'macrophysics')then
    if(allocated(pstate_cam%macrop_cmeliq_at_pc_full_level%f))then
        call check_non(pstate_cam%macrop_cmeliq_at_pc_full_level%f, ncol, 'pstate_cam%macrop_cmeliq_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%macrop_cld_at_pc_full_level%f))then
        call check_non(pstate_cam%macrop_cld_at_pc_full_level%f, ncol, 'pstate_cam%macrop_cld_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%macrop_concld_at_pc_full_level%f))then
        call check_non(pstate_cam%macrop_concld_at_pc_full_level%f, ncol, 'pstate_cam%macrop_concld_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%macrop_ast_at_pc_full_level%f))then
        call check_non(pstate_cam%macrop_ast_at_pc_full_level%f, ncol, 'pstate_cam%macrop_ast_at_pc_full_level', msg,lat,lon)
    end if
    end if

    if(trim(name) .eq. 'microphysics')then
    if(allocated(pstate_cam%microp_wsedl_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_wsedl_at_pc_full_level%f, ncol, 'pstate_cam%microp_wsedl_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_dei_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_dei_at_pc_full_level%f, ncol, 'pstate_cam%microp_dei_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_des_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_des_at_pc_full_level%f, ncol, 'pstate_cam%microp_des_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_res_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_res_at_pc_full_level%f, ncol, 'pstate_cam%microp_res_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_mu_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_mu_at_pc_full_level%f, ncol, 'pstate_cam%microp_mu_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_lambdac_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_lambdac_at_pc_full_level%f, ncol, 'pstate_cam%microp_lambdac_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_rei_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_rei_at_pc_full_level%f, ncol, 'pstate_cam%microp_rei_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_rel_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_rel_at_pc_full_level%f, ncol, 'pstate_cam%microp_rel_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_iciwp_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_iciwp_at_pc_full_level%f, ncol, 'pstate_cam%microp_iciwp_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_iclwp_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_iclwp_at_pc_full_level%f, ncol, 'pstate_cam%microp_iclwp_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_icswp_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_icswp_at_pc_full_level%f, ncol, 'pstate_cam%microp_icswp_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_cc_t_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_cc_t_at_pc_full_level%f, ncol, 'pstate_cam%microp_cc_t_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_cc_qv_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_cc_qv_at_pc_full_level%f, ncol, 'pstate_cam%microp_cc_qv_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_cc_ql_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_cc_ql_at_pc_full_level%f, ncol, 'pstate_cam%microp_cc_ql_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_cc_qi_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_cc_qi_at_pc_full_level%f, ncol, 'pstate_cam%microp_cc_qi_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_cc_nl_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_cc_nl_at_pc_full_level%f, ncol, 'pstate_cam%microp_cc_nl_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_cc_ni_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_cc_ni_at_pc_full_level%f, ncol, 'pstate_cam%microp_cc_ni_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%microp_cc_qlst_at_pc_full_level%f))then
        call check_non(pstate_cam%microp_cc_qlst_at_pc_full_level%f, ncol, 'pstate_cam%microp_cc_qlst_at_pc_full_level', msg,lat,lon)
    end if
    end if

    if(trim(name) .eq. 'radiation')then
    if(allocated(pstate_cam%lw_qrl_at_pc_full_level%f))then
        call check_non(pstate_cam%lw_qrl_at_pc_full_level%f, ncol, 'pstate_cam%lw_qrl_at_pc_full_level', msg,lat,lon)
    end if
    if(allocated(pstate_cam%sw_qrs_at_pc_full_level%f))then
        call check_non(pstate_cam%sw_qrs_at_pc_full_level%f, ncol, 'pstate_cam%sw_qrs_at_pc_full_level', msg,lat,lon)
    end if
    end if

    end subroutine phys_state_check


    subroutine check_nan_1d(var, ncol, varname, msg, lat, lon, gt, lt)
! io
    real(r8),         intent(in)    :: var(:)
    integer,          intent(in)    :: ncol      
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: msg
    real(r8),         intent(in)    :: lat(ncol)             ! longitude in radian
    real(r8),         intent(in)    :: lon(ncol)             ! latitude in radian  
 
    real(r8),         intent(in), optional  :: gt
    real(r8),         intent(in), optional  :: lt

    if(any(isnan(var(1:ncol))))then
        print*,trim(msg),':',trim(varname)
        print*,'rank=',mpi_rank()
        call endrun('check_non')
    endif

    if(present(gt))then
        if(any(var(1:ncol) .gt. gt))then
            print*,trim(msg),':'
            print*,trim(varname),' is greater than',gt
            print*,maxval(var(1:ncol))
            print*,'rank=',mpi_rank()
       !     call endrun('check_non')
        endif
    endif

    if(present(lt))then
        if(any(var(1:ncol) .lt. lt))then
            print*,trim(msg),':'
            print*,trim(varname),' is less than',lt
            print*,minval(var(1:ncol))
            print*,'rank=',mpi_rank()
       !     call endrun('check_non')
        endif
    endif

    end subroutine check_nan_1d

    subroutine check_nan_2d(var, ncol, varname, msg, lat, lon, gt, lt)
! io
    real(r8),         intent(in)    :: var(:,:)
    integer,          intent(in)    :: ncol      
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: msg
    real(r8),         intent(in)    :: lat(ncol)             ! longitude in radian
    real(r8),         intent(in)    :: lon(ncol)             ! latitude in radian  
 
 
    real(r8),         intent(in), optional  :: gt
    real(r8),         intent(in), optional  :: lt
! local
    integer                         :: i, k
    integer                         :: var_size1, var_size2

    if(any(isnan(var(:,1:ncol))))then
        print*,trim(msg),':',trim(varname)
        print*,'rank=',mpi_rank()

        if(print_check)then
            var_size1 = size(var,1)    
            !var_size2 = size(var,2)    
            var_size2 = ncol
  
            do k = 1, var_size1
            do i = 1, var_size2
                if(isnan(var(k,i)))then
                    print*,'k=',k,'i=',i,'lon=',lon(i)*rad2deg,'lat=',lat(i)*rad2deg
                end if
            end do
            end do
        end if

        call endrun('check_non')

    endif

    if(present(gt))then
        if(any(var(:,1:ncol) .gt. gt))then
            print*,trim(msg),':'
            print*,trim(varname),' is greater than',gt
            print*,maxval(var(:,1:ncol))
            print*,'rank=',mpi_rank()

            if(print_check)then
                var_size1 = size(var,1)    
                !var_size2 = size(var,2)    
                var_size2 = ncol
                do k = 1, var_size1
                do i = 1, var_size2
                    if(var(k,i) .gt. gt)then
                        print*,'k=',k,'i=',i,'var=',var(k,i),'lon=',lon(i)*rad2deg,'lat=',lat(i)*rad2deg
                    end if
                end do
                end do
            end if

       !     call endrun('check_non')
        endif
    endif

    if(present(lt))then
        if(any(var(:,1:ncol) .lt. lt))then
            print*,trim(msg),':'
            print*,trim(varname),' is less than',lt
            print*,minval(var(:,1:ncol))
            print*,'rank=',mpi_rank()

            if(print_check)then
                var_size1 = size(var,1)    
                !var_size2 = size(var,2)    
                var_size2 = ncol
                do k = 1, var_size1
                do i = 1, var_size2
                    if(var(k,i) .lt. lt)then
                        print*,'k=',k,'i=',i,'var=',var(k,i),'lon=',lon(i)*rad2deg,'lat=',lat(i)*rad2deg
                    end if
                end do
                end do
            end if

       !     call endrun('check_non')
        endif
    endif

    end subroutine check_nan_2d

    subroutine check_nan_3d(var, ncol, varname, msg, lat, lon, gt, lt)
! io
    real(r8),         intent(in)    :: var(:,:,:)
    integer,          intent(in)    :: ncol      
    character(len=*), intent(in)    :: varname
    character(len=*), intent(in)    :: msg
    real(r8),         intent(in)    :: lat(ncol)             ! longitude in radian
    real(r8),         intent(in)    :: lon(ncol)             ! latitude in radian  
    real(r8),         intent(in), optional  :: gt
    real(r8),         intent(in), optional  :: lt
! local
    integer                         :: i, k, m
    integer                         :: var_size1, var_size2, var_size3



    if(any(isnan(var(:,:,1:ncol))))then
        print*,trim(msg),':',trim(varname)
        print*,'rank=',mpi_rank()

        if(print_check)then
            var_size1 = size(var,1)    
            var_size2 = size(var,2)    
            !var_size3 = size(var,3)    
            var_size3 = ncol
            do m = 1, var_size1
            do k = 1, var_size2
            do i = 1, var_size3
                if(isnan(var(m,k,i)))then
                    print*,'m=',m,'k=',k,'i=',i,'lon=',lon(i)*rad2deg,'lat=',lat(i)*rad2deg
                end if
            end do
            end do
            end do
        end if

        call endrun('check_non')
    endif

    if(present(gt))then
        if(any(var(:,:,1:ncol) .gt. gt))then
            print*,trim(msg),':'
            print*,trim(varname),' is greater than',gt
            print*,maxval(var(:,:,1:ncol))
            print*,'rank=',mpi_rank()

            if(print_check)then
                var_size1 = size(var,1)    
                var_size2 = size(var,2)    
                !var_size3 = size(var,3)    
                var_size3 = ncol
                do m = 1, var_size1
                do k = 1, var_size2
                do i = 1, var_size3
                    if(var(m,k,i) .gt. gt)then
                        print*,'m=',m,'k=',k,'i=',i,'var=',var(m,k,i),'lon=',lon(i)*rad2deg,'lat=',lat(i)*rad2deg
                    end if
                end do
                end do
                end do
            end if

      !      call endrun('check_non')
        endif
    endif

    if(present(lt))then
        if(any(var(:,:,1:ncol) .lt. lt))then
            print*,trim(msg),':'
            print*,trim(varname),' is less than',lt
            print*,minval(var(:,:,1:ncol))
            print*,'rank=',mpi_rank()

            if(print_check)then
                var_size1 = size(var,1)    
                var_size2 = size(var,2)    
                !var_size3 = size(var,3)    
                var_size3 = ncol
                do m = 1, var_size1
                do k = 1, var_size2
                do i = 1, var_size3
                    if(var(m,k,i) .lt. lt)then
                        print*,'m=',m,'k=',k,'i=',i,'var=',var(m,k,i),'lon=',lon(i)*rad2deg,'lat=',lat(i)*rad2deg
                    end if
                end do
                end do
                end do
            end if

      !      call endrun('check_non')
        endif
    endif

    end subroutine check_nan_3d


end module grist_phys_state_check

