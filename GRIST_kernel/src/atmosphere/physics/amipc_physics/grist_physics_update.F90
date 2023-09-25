!======================================================
!
!  Created by LiXiaohan on 19/6/30.
!  contains:
!  coupling of each physical process in the physpkg
!  time level update in physics variables
!======================================================

 module grist_physics_update


    use grist_constants,                    only: i4, r8, gravity, rdry, cp, zvir
    use grist_nml_module,                   only: nlev, nlevp, ntracer
 
    implicit none
    private
    save

    public  :: ptimelevels,                 &
               old_time_level,              &
               phys_time_level_init,        &
               phys_time_level_update,      &
               phys_update,                 &
               geopotential_hydrostatic,    &
               geopotential_dse

    integer :: ptimelevels
    integer :: old_time_level

contains

    subroutine phys_time_level_init
    use grist_nml_module,                   only: working_mode

    if(trim(working_mode) .eq. 'scm')then
#ifdef SCAM
        ! for three time level integration
        ptimelevels = 2
#else
        ptimelevels = 1
#endif
    else
        ! for two time level integration
        ptimelevels = 1
    end if

    old_time_level = 1

    end subroutine phys_time_level_init


    subroutine phys_time_level_update

    old_time_level = mod(old_time_level, ptimelevels) + 1

    end subroutine phys_time_level_update


    subroutine phys_update(subname, dt, ncol, ptend1, ptend2)
    use grist_physics_data_structure,       only: pstate, ptend_f3,                 &
                                                  phy_tracer_info
    use grist_PhysC_data_structure,         only: tend_in_cam_physics
    use grist_mpi

! io
    character(*),          intent(in)              :: subname
    real(r8),              intent(in)              :: dt
    integer(i4),           intent(in)              :: ncol
    type(tend_in_cam_physics), intent(inout)           :: ptend1
    type(tend_in_cam_physics), intent(inout), optional :: ptend2
! local
    integer(i4) :: n, k, i
    integer(i4) :: ixcldice, ixcldliq, ixnumice, ixnumliq
    real(r8)    :: q_old(ntracer,nlev,ncol)

    ! inquire cloud index
    ixcldice = -1; ixcldliq = -1; ixnumice = -1; ixnumliq = -1
    do n = 1, ntracer
        if(phy_tracer_info(n)%longname .eq. 'cloud_liquid')        ixcldliq = n
        if(phy_tracer_info(n)%longname .eq. 'cloud_ice')           ixcldice = n
        if(phy_tracer_info(n)%longname .eq. 'cloud_liquid_number') ixnumliq = n
        if(phy_tracer_info(n)%longname .eq. 'cloud_ice_number')    ixnumice = n
    end do

     q_old(:,:,1:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)

    if(allocated(ptend1%tend_u%f))then
         pstate%u_wind_at_pc_full_level%f(:,1:ncol)        = pstate%u_wind_at_pc_full_level%f(:,1:ncol)         &
                                                           + ptend1%tend_u%f(:,1:ncol)*dt
         ptend_f3%tend_u_wind_at_pc_full_level%f(:,1:ncol) = ptend_f3%tend_u_wind_at_pc_full_level%f(:,1:ncol)  &
                                                           + ptend1%tend_u%f(:,1:ncol)
    !     ptend1%tend_u%f = 0._r8
    end if

    if(allocated(ptend1%tend_v%f))then
         pstate%v_wind_at_pc_full_level%f(:,1:ncol)        = pstate%v_wind_at_pc_full_level%f(:,1:ncol)         &
                                                           + ptend1%tend_v%f(:,1:ncol)*dt
         ptend_f3%tend_v_wind_at_pc_full_level%f(:,1:ncol) = ptend_f3%tend_v_wind_at_pc_full_level%f(:,1:ncol)  &
                                                           + ptend1%tend_v%f(:,1:ncol)
    !     ptend1%tend_v%f = 0._r8
    end if

    if(allocated(ptend1%tend_s%f))then
         pstate%static_energy_at_pc_full_level%f(:,1:ncol) = pstate%static_energy_at_pc_full_level%f(:,1:ncol)       &
                                                           + ptend1%tend_s%f(:,1:ncol)*dt
         ptend_f3%tend_temp_at_pc_full_level%f(:,1:ncol)   = ptend_f3%tend_temp_at_pc_full_level%f(:,1:ncol)  &
                                                           + ptend1%tend_s%f(:,1:ncol)/cp
    !     ptend1%tend_s%f = 0._r8
    end if

    if(allocated(ptend1%tend_q%f))then
         pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)        = pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)       &
                                                                  + ptend1%tend_q%f(:,:,1:ncol)*dt
    !     ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(:,:,1:ncol) = ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)  &
    !                                                              + ptend1%tend_q%f(:,:,1:ncol)
    !     ptend1%tend_q%f = 0._r8
    end if

    if(present(ptend2))then
    if(allocated(ptend2%tend_u%f))then
        pstate%u_wind_at_pc_full_level%f(:,1:ncol)        = pstate%u_wind_at_pc_full_level%f(:,1:ncol) &
                                                          + ptend2%tend_u%f(:,1:ncol)*dt
        ptend_f3%tend_u_wind_at_pc_full_level%f(:,1:ncol) = ptend_f3%tend_u_wind_at_pc_full_level%f(:,1:ncol)  &
                                                          + ptend2%tend_u%f(:,1:ncol)
    !    ptend2%tend_u%f = 0._r8
    end if

    if(allocated(ptend2%tend_v%f))then
        pstate%v_wind_at_pc_full_level%f(:,1:ncol)        = pstate%v_wind_at_pc_full_level%f(:,1:ncol) &
                                                          + ptend2%tend_v%f(:,1:ncol)*dt
        ptend_f3%tend_v_wind_at_pc_full_level%f(:,1:ncol) = ptend_f3%tend_v_wind_at_pc_full_level%f(:,1:ncol)  &
                                                          + ptend2%tend_v%f(:,1:ncol)
    !    ptend2%tend_v%f = 0._r8
    end if

    if(allocated(ptend2%tend_s%f))then
        pstate%static_energy_at_pc_full_level%f(:,1:ncol) = pstate%static_energy_at_pc_full_level%f(:,1:ncol)   &
                                                          + ptend2%tend_s%f(:,1:ncol)*dt
        ptend_f3%tend_temp_at_pc_full_level%f(:,1:ncol)   = ptend_f3%tend_temp_at_pc_full_level%f(:,1:ncol)  &
                                                          + ptend2%tend_s%f(:,1:ncol)/cp
    !    ptend2%tend_s%f = 0._r8
    end if

    if(allocated(ptend2%tend_q%f))then
        pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)        = pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)   &
                                                                 + ptend2%tend_q%f(:,:,1:ncol)*dt
    !    ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(:,:,1:ncol) = ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)  &
    !                                                             + ptend2%tend_q%f(:,:,1:ncol)
    !    ptend2%tend_q%f = 0._r8
    end if
    end if

    do n = 1, ntracer
        if(n .eq. ixnumliq .or. n .eq. ixnumice)then
             do k = 1, nlev
                do i = 1, ncol
                   pstate%tracer_mxrt_at_pc_full_level%f(n,k,i) = max(1.e-12_r8,pstate%tracer_mxrt_at_pc_full_level%f(n,k,i))
                   pstate%tracer_mxrt_at_pc_full_level%f(n,k,i) = min(1.e10_r8,pstate%tracer_mxrt_at_pc_full_level%f(n,k,i))
                 end do
             end do
        else
             call qneg3(subname, ncol, nlev, n, n, phy_tracer_info(n)%qmin, &
                        pstate%tracer_mxrt_at_pc_full_level%f(n,:,1:ncol))
        end if
    end do

    ! Special tests for cloud liquid and ice: Enforce a minimum non-zero value.
    if(subname .eq. 'deep convection')then
        call state_cnst_min_nz(1.e-36_r8, ixcldliq, ixnumliq)
        call state_cnst_min_nz(1.e-36_r8, ixcldice, ixnumice)
    end if

    if(allocated(ptend1%tend_s%f) .or. allocated(ptend1%tend_q%f) .or.                     &
      (present(ptend2) .and. (allocated(ptend2%tend_s%f) .or. allocated(ptend2%tend_q%f))) )then

    call geopotential_dse(ncol, pstate%pressure_at_pc_face_level%f(:,1:ncol),              &
                          pstate%pressure_at_pc_full_level%f(:,1:ncol),                    &
                          pstate%delp_at_pc_full_level%f(:,1:ncol),                        &
                          pstate%geop_at_pc_surface%f(1:ncol),                             &
                          pstate%static_energy_at_pc_full_level%f(:,1:ncol),               &
                          pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol),               &
                          pstate%temp_at_pc_full_level%f(:,1:ncol),                        &
                          pstate%z_at_pc_full_level%f(:,1:ncol),                           &
                          pstate%z_at_pc_face_level%f(:,1:ncol))
    end if

    if(allocated(ptend1%tend_q%f))then
         ptend1%tend_q%f(:,:,1:ncol)                              = (pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol) - q_old(:,:,1:ncol))/dt
         ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(:,:,1:ncol) = ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)  &
                                                                  + ptend1%tend_q%f(:,:,1:ncol)
 
         if(present(ptend2) .and. allocated(ptend2%tend_q%f))then
             ptend2%tend_q%f = 0._r8
         end if
    end if

    contains

    subroutine state_cnst_min_nz(lim, qix, numix)
        ! Small utility function for setting minimum nonzero
        ! constituent concentrations.

        ! Lower limit and constituent index
        real(r8), intent(in) :: lim
        integer,  intent(in) :: qix
        integer,  intent(in) :: numix

        if (numix > 0) then
            ! Where q is too small, zero mass and number concentration.
            where (pstate%tracer_mxrt_at_pc_full_level%f(qix,:,:ncol) < lim)
                pstate%tracer_mxrt_at_pc_full_level%f(qix,:,:ncol) = 0._r8
                pstate%tracer_mxrt_at_pc_full_level%f(numix,:,:ncol) = 0._r8
            end where
        else
            ! If no number index, just do mass.
            where (pstate%tracer_mxrt_at_pc_full_level%f(qix,:,:ncol) < lim)
                pstate%tracer_mxrt_at_pc_full_level%f(qix,:,:ncol) = 0._r8
            end where
        end if
    end subroutine state_cnst_min_nz


    end subroutine phys_update

    subroutine geopotential_hydrostatic(ncol, temp, pface, pfull, dpfull, qv, z_full, z_face)
! io
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: temp(nlev, ncol)
    real(r8), intent(in)  :: pface(nlevp, ncol)
    real(r8), intent(in)  :: pfull(nlev, ncol)
    real(r8), intent(in)  :: dpfull(nlev, ncol)
    real(r8), intent(in)  :: qv(nlev, ncol)

    real(r8), intent(out) :: z_full(nlev, ncol)
    real(r8), intent(out) :: z_face(nlevp, ncol)

! local
    integer  :: k, i
    real(r8) :: hkl(ncol), hkk(ncol)
    real(r8) :: geo_full(nlev, ncol), geo_face(nlevp, ncol)
    real(r8) :: tvfac, tv

    ! The surface height is zero by definition ?
    geo_face(nlevp,:) = 0._r8

    do k = nlev, 1, -1
        do i = 1, ncol
            hkl(i) = log(pface(k+1,i))-log(pface(k,i))
            hkk(i) = 1._r8-pface(k,i)*hkl(i)/dpfull(k,i)
        end do
        
        do i = 1, ncol
            tvfac = 1._r8+zvir*qv(k,i)
            tv    = temp(k,i)*tvfac
            geo_face(k,i) = geo_face(k+1,i) + rdry*tv*hkl(i)
            geo_full(k,i) = geo_face(k+1,i) + rdry*tv*hkk(i)

            z_full(k,i)   = geo_full(k,i)/gravity
            z_face(k,i)   = geo_face(k,i)/gravity
        end do
    end do

   end subroutine geopotential_hydrostatic


   subroutine geopotential_dse(ncol, pface, pfull, dpfull, phis, dse, qv,     &
                               temp, z_full, z_face)
! io
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: pface(nlevp, ncol)
    real(r8), intent(in)  :: pfull(nlev, ncol)
    real(r8), intent(in)  :: dpfull(nlev, ncol)
    real(r8), intent(in)  :: qv(nlev, ncol)
    real(r8), intent(in)  :: phis(ncol)
    real(r8), intent(in)  :: dse(nlev, ncol)

    real(r8), intent(out) :: temp(nlev, ncol)
    real(r8), intent(out) :: z_full(nlev, ncol)
    real(r8), intent(out) :: z_face(nlevp, ncol)
! local
    integer  :: k, i
    real(r8) :: hkl(ncol), hkk(ncol)
    real(r8) :: tvfac, tv, rog
    
    rog = rdry/gravity

    do i = 1, ncol
        z_face(nlevp, i) = 0._r8
    end do

    do k = nlev, 1, -1
        do i = 1, ncol
            hkl(i) = log(pface(k+1,i))-log(pface(k,i))
            hkk(i) = 1._r8-pface(k,i)*hkl(i)/dpfull(k,i)
        end do
        
        do i = 1, ncol
            tvfac       = 1._r8+zvir*qv(k,i)
            tv          = (dse(k,i)-phis(i)-gravity*z_face(k+1,i))/(cp/tvfac+rdry*hkk(i))
            temp(k,i)   = tv/tvfac
            
            z_full(k,i) = z_face(k+1,i)+rog*tv*hkk(i)
            z_face(k,i) = z_face(k+1,i)+rog*tv*hkl(i)
        end do
    end do

   end subroutine geopotential_dse

 end module grist_physics_update
