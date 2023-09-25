!======================================================
!
!  Created by LiXiaohan on 19/5/20.
!  Including:
!  dyn-core Euralian update using the var_obs;
!  dyn-core and physics coupling in SCM;
!  update var_obs if needed.
!======================================================

 module grist_scm_coupling_module
    use grist_constants,                    only: i4, r8, latvap, gravity, rdry, cp, stebol,    &
                                                  zvir, p00, fillvalue, zero
    use grist_nml_module,                   only: nlev, nlevp, ntracer, model_timestep, physpkg
    use grist_physics_data_structure,       only: pstate, ptend_f3, phy_tracer_info
#ifdef AMIPC_PHYSICS 
    use grist_physics_update,               only: geopotential_hydrostatic
#endif

#ifdef AMIPW_PHYSICS
    use grist_PhysW_data_structure,         only: pstate_wrf, psurf_wrf, ptend_wrf
#endif

    use grist_scm_dyn_vars
    use grist_scm_comm_module
    use grist_handle_error,                 only: endrun
    use grist_mpi

    implicit none

    private
    public  :: grist_scm_dtp_coupling,      &
               grist_scm_ptd_coupling,      &
               grist_scm_coupling_init 

contains

    subroutine grist_scm_coupling_init
    use grist_scm_pressure_update,          only: time_integration_renew_mass_state
    use grist_hpe_constants,                only: deta_face, deta_full

! local
    integer(i4)                         :: k, m 

! geopotential at surface and albedo, same as SCAM test cases
    if(trim(scm_test_name) .eq. 'arm95' .or. trim(scm_test_name) .eq. 'ARM95' .or.        &
       trim(scm_test_name) .eq. 'arm97' .or. trim(scm_test_name) .eq. 'ARM97')then
        pstate%geop_at_pc_surface%f(1)   = 4930.40595298484
        pstate%landfrac_at_pc_surface%f  = 1._r8
        pstate%icefrac_at_pc_surface%f   = 0._r8
        pstate%ocnfrac_at_pc_surface%f   = 0._r8

    elseif(trim(scm_test_name) .eq. 'bomex' .or. trim(scm_test_name) .eq. 'BOMEX' .or.      &
           trim(scm_test_name) .eq. 'gateIII' .or. trim(scm_test_name) .eq. 'GATEIII' .or.  &
           trim(scm_test_name) .eq. 'dycoms' .or. trim(scm_test_name) .eq. 'DYCOMS' .or.  &
           trim(scm_test_name) .eq. 'cgilsS11' .or. trim(scm_test_name) .eq. 'CGILSS11'.or.  &
           trim(scm_test_name) .eq. 'cgilsS6' .or. trim(scm_test_name) .eq. 'CGILSS6' .or.  &
           trim(scm_test_name) .eq. 'cgilsS12' .or. trim(scm_test_name) .eq. 'CGILSS12')then 
       pstate%geop_at_pc_surface%f(1)   = 0._r8
       pstate%landfrac_at_pc_surface%f  = 0._r8
       pstate%icefrac_at_pc_surface%f   = 0._r8
       pstate%ocnfrac_at_pc_surface%f   = 1._r8

    elseif(trim(scm_test_name) .eq. 'twp06' .or. trim(scm_test_name) .eq. 'TWP06')then
       pstate%geop_at_pc_surface%f(1)   = 0._r8
       pstate%landfrac_at_pc_surface%f  = 0._r8
       pstate%icefrac_at_pc_surface%f   = 0._r8
       pstate%ocnfrac_at_pc_surface%f   = 1._r8

       !pstate%geop_at_pc_surface%f(1)   = 1213.06045693356
       !pstate%landfrac_at_pc_surface%f  = 1._r8
       !pstate%icefrac_at_pc_surface%f   = 0._r8
       !pstate%ocnfrac_at_pc_surface%f   = 0._r8

    else
        if(mpi_rank()==0)then
            print*,'The SCM test has no geopotential height, please set it in grist_scm_coupling_module.F90'
            call endrun('stop at subroutine grist_scm_coupling_init')
        end if
    end if
    end subroutine grist_scm_coupling_init


    subroutine grist_scm_dtp_coupling
    use grist_scm_pressure_update,          only: time_integration_renew_mass_state
    use grist_hpe_constants,                only: deta_face, deta_full

! local
    integer(i4)                         :: k, m 

    pstate%u_wind_at_pc_full_level%f(:,1)         = scm_u(:)
    pstate%v_wind_at_pc_full_level%f(:,1)         = scm_v(:)
    pstate%omega_at_pc_full_level%f(:,1)          = scm_omega(:)
    pstate%temp_at_pc_full_level%f(:,1)           = scm_t3(:,n3m2)
    pstate%pressure_at_pc_surface%f(1)            = scm_ps3(n3m2)
    pstate%tracer_mxrt_at_pc_full_level%f(1:ntracer,:,1)  = scm_q3(1:ntracer,:,n3m2)

#ifdef AMIPC_PHYSICS 
    do m = 1, ntracer
    where(pstate%tracer_mxrt_at_pc_full_level%f(m,:,1) .lt. phy_tracer_info(m)%qmin)   &
          pstate%tracer_mxrt_at_pc_full_level%f(m,:,1) = phy_tracer_info(m)%qmin
    end do
#endif
 
!    if(trim(scm_test_name) .eq. 'bomex' .or. trim(scm_test_name) .eq.'BOMEX')then
!        !use the LHF,SHF,SST,ps in the UNICON paper (Shin and Park, 2020), LiXH
!        pstate%sst_at_pc_surface%f(1) = 300.4_r8
!        pstate%ts_at_pc_surface%f(1)  = 300.4_r8
!        pstate%pressure_at_pc_surface%f(1) = 101500._r8
!    end if

! update model level pressure
    call time_integration_renew_mass_state(1, pstate%pressure_at_pc_surface%f(1)     , &
                                              pstate%pressure_at_pc_face_level%f(:,1), &
                                              pstate%delp_at_pc_full_level%f(:,1)    , &
                                              pstate%pressure_at_pc_full_level%f(:,1), &
                                              pstate%delp_at_pc_face_level%f(:,1))

    do k = 1, nlev
#ifdef EXNER_P00 
        pstate%exner_at_pc_full_level%f(k,1)      = (p00/pstate%pressure_at_pc_full_level%f(k,1))**(rdry/cp)
#else
        pstate%exner_at_pc_full_level%f(k,1)      = (pstate%pressure_at_pc_surface%f(1)/pstate%pressure_at_pc_full_level%f(k,1))**(rdry/cp)
#endif
    end do

! update geopotential and model level height
! in physpkg, z_full and z_face do not contain geopotential surface height.   LiXH
    call geopotential_hydrostatic(1, pstate%temp_at_pc_full_level%f(:,1),             &
                                  pstate%pressure_at_pc_face_level%f(:,1),            &
                                  pstate%pressure_at_pc_full_level%f(:,1),            &
                                  pstate%delp_at_pc_full_level%f(:,1),                &
                                  pstate%tracer_mxrt_at_pc_full_level%f(1,:,1),       &
                                  pstate%z_at_pc_full_level%f(:,1),                   &
                                  pstate%z_at_pc_face_level%f(:,1) )

    do k = 1, nlev
        pstate%static_energy_at_pc_full_level%f(k,1) = cp*pstate%temp_at_pc_full_level%f(k,1)       &
                                                      +gravity*pstate%z_at_pc_full_level%f(k,1)     & 
                                                      +pstate%geop_at_pc_surface%f(1)
    end do

! geopotential at surface and albedo, same as SCAM test cases
    if(use_phys_vars)then
        if(have_tg)then
            pstate%ts_at_pc_surface%f(1)                  = tground
        else
            pstate%ts_at_pc_surface%f(1)                  = tobs(nlev)
        end if

        if(pstate%ocnfrac_at_pc_surface%f(1) .eq. 1._r8)then
            pstate%sst_at_pc_surface%f(1) = pstate%ts_at_pc_surface%f(1)
        end if

#if (defined AMIPC_PHYSICS || defined WRFCAM_RRTMG)
        pstate%atm_in_lwup_at_pc_surface%f(1)             = stebol * pstate%ts_at_pc_surface%f(1)**4
#endif
    end if

#ifdef AMIPW_PHYSICS
    do k = 1, nlev
        pstate_wrf%u_phy (1,nlev+1-k,1)  = pstate%u_wind_at_pc_full_level%f(k,1)
        pstate_wrf%v_phy (1,nlev+1-k,1)  = pstate%v_wind_at_pc_full_level%f(k,1)
        pstate_wrf%t_phy (1,nlev+1-k,1)  = pstate%temp_at_pc_full_level%f(k,1)
        pstate_wrf%p_phy (1,nlev+1-k,1)  = pstate%pressure_at_pc_full_level%f(k,1)
        pstate_wrf%th_phy(1,nlev+1-k,1)  = pstate_wrf%t_phy(1,nlev+1-k,1)*((p00/pstate_wrf%p_phy(1,nlev+1-k,1))**(rdry/cp))
        pstate_wrf%zzz   (1,nlev+1-k,1)  = pstate%z_at_pc_full_level%f(k,1)+pstate%geop_at_pc_surface%f(1)/gravity
        !Note: exner function in CAM = (p00/p)**(Rdry/cp)
        !      exner function in WRF = (p/p00)**(Rdry/cp)
        pstate_wrf%pi_phy(1,nlev+1-k,1)  = 1._r8/pstate%exner_at_pc_full_level%f(k,1) 
        pstate_wrf%rhom  (1,nlev+1-k,1)  = (pstate%pressure_at_pc_face_level%f(k+1,1)-pstate%pressure_at_pc_face_level%f(k,1))/&
                                           ((pstate%z_at_pc_face_level%f(k,1)-pstate%z_at_pc_face_level%f(k+1,1))*gravity) 
        pstate_wrf%moist (1,nlev+1-k,1,1:ntracer)  = pstate%tracer_mxrt_at_pc_full_level%f(1:ntracer,k,1)
        pstate_wrf%dz8w(  1,nlev+1-k,1)  = pstate%z_at_pc_face_level%f(k,1)-pstate%z_at_pc_face_level%f(k+1,1)
        pstate_wrf%omega( 1,nlev+1-k,1)  = pstate%omega_at_pc_full_level%f(k,1)
        ptend_wrf%rthdyten(1,nlev+1-k,1) = t_tend_dycore(k)*((p00/pstate_wrf%p_phy(1,nlev+1-k,1))**(rdry/cp))
        ptend_wrf%rqvdyten(1,nlev+1-k,1) = q_tend_dycore(k) 
    end do

    pstate_wrf%u_phy (1,nlevp,1)  = fillvalue
    pstate_wrf%v_phy (1,nlevp,1)  = fillvalue
    pstate_wrf%t_phy (1,nlevp,1)  = fillvalue
    pstate_wrf%th_phy(1,nlevp,1)  = fillvalue
    pstate_wrf%p_phy (1,nlevp,1)  = fillvalue
    pstate_wrf%zzz   (1,nlevp,1)  = fillvalue
    pstate_wrf%pi_phy(1,nlevp,1)  = fillvalue
    pstate_wrf%rhom  (1,nlevp,1)  = fillvalue
    pstate_wrf%moist (1,nlevp,1,1:ntracer) = fillvalue
    pstate_wrf%dz8w  (1,nlevp,1)  = fillvalue
    pstate_wrf%omega (1,nlevp,1)  = fillvalue
    pstate_wrf%w0avg (1,nlevp,1)  = fillvalue

    pstate_wrf%rhod=pstate_wrf%rhom         !LiXH, temperorily set, not used in WSM6
! w-level var
    call www_at_face(nlevp, pstate_wrf%p_phy(1,:,1), pstate_wrf%zzz(1,:,1), pstate_wrf%omega(1,:,1), &
                     pstate_wrf%www(1,:,1))

    pstate_wrf%w0avg(1,1:nlev,1)  = -pstate_wrf%omega(1,1:nlev,1)/gravity/pstate_wrf%rhom(1,1:nlev,1)       !for hydro-dynamics
    pstate_wrf%p8w(1,1:nlevp,1)   = pstate%pressure_at_pc_face_level%f(nlevp:1:-1,1)
    pstate_wrf%z8w(1,1:nlevp,1)   = pstate%z_at_pc_face_level%f(nlevp:1:-1,1) 
! interpolate w-level t based on m-level t
    do k = 2, nlev ! ilev level of WRF > nlev+2-ilev face of GRIST, the surrounding full levels are nlev+2-ilev, nlev+1-ilev
        pstate_wrf%t8w(1,k,1) = 0.5*(deta_full(nlev+1-k)/deta_face(nlev+2-k)*pstate%temp_at_pc_full_level%f(nlev+2-k,1)+&
                                     deta_full(nlev+2-k)/deta_face(nlev+2-k)*pstate%temp_at_pc_full_level%f(nlev+1-k,1))
    end do
! we can use some extrapolate here but simply use full-level value now
    pstate_wrf%t8w(1,nlevp,1)   = pstate%temp_at_pc_full_level%f(1,1)

! surface
    psurf_wrf%ht  (1,1)  = pstate%geop_at_pc_surface%f(1)/gravity
    psurf_wrf%psfc(1,1)  = pstate%pressure_at_pc_surface%f(1)
    psurf_wrf%tsk (1,1)  = pstate%ts_at_pc_surface%f(1)
    pstate_wrf%t8w(1,1,1)= psurf_wrf%tsk (1,1) 


    contains

    subroutine www_at_face(nlevp, wrf_p_phys, wrf_z_phys, wrf_omega, wrf_w_face)
!io 
    integer,  intent(in)    :: nlevp
    real(r8), intent(in)    :: wrf_p_phys(:)
    real(r8), intent(in)    :: wrf_z_phys(:)
    real(r8), intent(in)    :: wrf_omega(:)
    real(r8), intent(out)   :: wrf_w_face(:)
!local
    real(r8)                :: rho_face(nlevp)
    real(r8)                :: omega_face(nlevp)

    omega_face(1:nlevp) = wfldint(nlevp:1:-1)
    do k = 2, nlevp-1
        rho_face(k)   = (wrf_p_phys(k-1)-wrf_p_phys(k))/((wrf_z_phys(k)-wrf_z_phys(k-1))*gravity)
        wrf_w_face(k) = -omega_face(k)/(rho_face(k)*gravity)
    end do

    wrf_w_face(1)     = 0._r8
    wrf_w_face(nlevp) = 0._r8

    end subroutine www_at_face

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

#endif

    end subroutine grist_scm_dtp_coupling


    subroutine grist_scm_ptd_coupling(dtime)
! io
    real(r8), intent(in)    :: dtime
! local
    integer :: m, k

! physics tendency
#ifdef AMIPC_PHYSICS 
    scm_u_tend(:) = ptend_f3%tend_u_wind_at_pc_full_level%f(:,1)
    scm_v_tend(:) = ptend_f3%tend_v_wind_at_pc_full_level%f(:,1)
    scm_t_tend(:) = ptend_f3%tend_temp_at_pc_full_level%f(:,1)

#ifdef SCAM
    scm_qminus(1:ntracer,:) = pstate%tracer_mxrt_at_pc_full_level%f(1:ntracer,:,1)
    do m = 1, ntracer
        where(scm_qminus(m,:) .lt. phy_tracer_info(m)%qmin)   &
              scm_qminus(m,:) = phy_tracer_info(m)%qmin
    end do

#else
    scm_q_tend(1:ntracer,:) = ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1:ntracer,:,1)
#endif

    ptend_f3%tend_u_wind_at_pc_full_level%f      = 0._r8
    ptend_f3%tend_v_wind_at_pc_full_level%f      = 0._r8
    ptend_f3%tend_temp_at_pc_full_level%f        = 0._r8
    ptend_f3%tend_tracer_mxrt_at_pc_full_level%f = 0._r8
#endif

#ifdef AMIPW_PHYSICS
    do k = 1, nlev
    scm_u_tend(k) = ptend_f3%tend_u_wind_at_pc_full_level%f(nlev+1-k,1)
    scm_v_tend(k) = ptend_f3%tend_v_wind_at_pc_full_level%f(nlev+1-k,1)
    scm_t_tend(k) = ptend_f3%tend_potential_temp_at_pc_full_level%f(nlev+1-k,1)*((pstate%pressure_at_pc_full_level%f(k,1)/p00)**(rdry/cp))

#ifdef SCAM

    scm_qminus(1:ntracer,k) = pstate_wrf%moist(1,nlev+1-k,1,1:ntracer)
    do m = 1, ntracer
        if(scm_qminus(m,k) .lt. phy_tracer_info(m)%qmin)   &
              scm_qminus(m,k) = phy_tracer_info(m)%qmin
    end do

#else
    scm_q_tend(1:ntracer,k) = ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1:ntracer,nlev+1-k,1)
#endif

    end do

    ptend_f3%tend_u_wind_at_pc_full_level%f         = 0._r8
    ptend_f3%tend_v_wind_at_pc_full_level%f         = 0._r8
    ptend_f3%tend_potential_temp_at_pc_full_level%f = 0._r8
    ptend_f3%tend_tracer_mxrt_at_pc_full_level%f    = 0._r8

#endif

    end subroutine grist_scm_ptd_coupling

 end module grist_scm_coupling_module
