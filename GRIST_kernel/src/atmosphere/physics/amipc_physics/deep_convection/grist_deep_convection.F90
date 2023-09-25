!===================================================================================================
!
!  Created by LiXiaohan on 19/07/10, adopted from CAM5
!
!     interface of the deep convection scheme 
!     only the Zhang&Mcfarlane scheme is avaiable now
!===================================================================================================

 module grist_deep_convection

    use grist_handle_error,                 only: endrun
    use grist_nml_module,                   only: ntracer, nlev, nlevp
    use grist_constants,                    only: r8, i4, p00,            &
                                                  cp, cpliq, cpwv, tmelt, &
                                                  epsilo, latvap, latice, &
                                                  gravity, rdry, rvap
    use grist_mpi
    implicit none
    private
    public          :: convect_deep_init,        &
                       end_of_convect_deep,      &
                       convect_deep_tend,        &
                       convect_deep_tend_2

! Private:
    character(len=16) :: deep_scheme
! if the variables below are used by other modules (pbuf in CAM), change to pstate%~, Lixh
    real(r8),    allocatable     :: evapcdp(:,:)    ! 'NEVAPR_DPCU' Evaporation of deep convection precipitation [ kg/kg/s ]
    real(r8),    allocatable     :: flxprec(:,:)    ! 'DP_FLXPRC'   Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8),    allocatable     :: flxsnow(:,:)    ! 'DP_FLXSNW'   Convective-scale flux of snow   at interfaces (kg/m2/s)
! other pbuf includes:
!   'DP_CLDLIQ' , 'DP_CLDICE'

 contains

    subroutine convect_deep_init(ncol)
    use grist_hpe_constants,                only: eta_face
    use grist_zm_conv,                      only: zm_conv_init
    use phys_control,                       only: phys_deepconv_pbl, phys_getopts
    use grist_PhysC_data_structure,         only: pstate_cam
! io
    integer(i4), intent(in)     :: ncol
! local
    real(r8)                    :: ref_pressure_face(nlevp)
    logical                     :: no_deep_pbl

    ! Initialize pstate:
    allocate(pstate_cam%scalar_precc_dp_surface%f(ncol))
    allocate(pstate_cam%dp_icwmr_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%rprddp_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%updraft_mass_flux%f(nlevp,ncol))

    pstate_cam%scalar_precc_dp_surface%pos       = 0
    pstate_cam%dp_icwmr_at_pc_full_level%pos     = 0
    pstate_cam%rprddp_at_pc_full_level%pos       = 0
    pstate_cam%updraft_mass_flux%pos             = 0
 
    pstate_cam%dp_icwmr_at_pc_full_level%f       = 0._r8
    pstate_cam%scalar_precc_dp_surface%f         = 0._r8
    pstate_cam%rprddp_at_pc_full_level%f         = 0._r8
    pstate_cam%updraft_mass_flux%f               = 0._r8
 
    call phys_getopts(deep_scheme_out = deep_scheme)

    select case ( deep_scheme)
    case('off')
        if(mpi_rank()==0)print*,'convect_deep: no deep convection selected'
    case('ZM')
        allocate(evapcdp(nlev,  ncol))
        allocate(flxprec(nlevp, ncol))
        allocate(flxsnow(nlevp, ncol))

        ref_pressure_face(:) = eta_face(:)*p00
        no_deep_pbl          = phys_deepconv_pbl()

        call zm_conv_init(ncol   , nlev    , nlevp   , ntracer   ,  &
                          ref_pressure_face, tmelt   , epsilo    ,  &
                          latvap , latice  , cp      , cpliq     ,  &
                          cpwv   , gravity , rdry    , rvap      ,  &
                          no_deep_pbl)
    end select

    end subroutine convect_deep_init


    subroutine end_of_convect_deep
    use grist_PhysC_data_structure,      only: pstate_cam
    use grist_zm_conv,                   only: zm_conv_end

    select case ( deep_scheme)
    case('off')
    case('ZM')
        if(allocated(evapcdp)) deallocate(evapcdp)
        if(allocated(flxprec)) deallocate(flxprec)
        if(allocated(flxsnow)) deallocate(flxsnow)
#ifndef CMAPI
        call zm_conv_end
#endif
    end select
#ifndef CMAPI
    if(allocated(pstate_cam%scalar_precc_dp_surface%f))          & 
      deallocate(pstate_cam%scalar_precc_dp_surface%f)

    if(allocated(pstate_cam%dp_icwmr_at_pc_full_level%f))        &
      deallocate(pstate_cam%dp_icwmr_at_pc_full_level%f)

    if(allocated(pstate_cam%rprddp_at_pc_full_level%f))          &
      deallocate(pstate_cam%rprddp_at_pc_full_level%f)

    if(allocated(pstate_cam%updraft_mass_flux%f))                &
      deallocate(pstate_cam%updraft_mass_flux%f)
#endif

    end subroutine end_of_convect_deep


    subroutine convect_deep_tend(ncol, cme, dlf, pflx, zdu, rliq, dt)
    use grist_zm_conv,                      only: zm_conv_tend
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam, ptend_deep_convection
    use grist_physics_update,               only: old_time_level
 
! io
    integer(i4), intent(in) :: ncol
    real(r8),    intent(in) :: dt                  ! delta t (model time increment)
      

    real(r8),    intent(out) :: dlf(nlev, ncol)    ! scattrd version of the detraining cld h2o tend
    real(r8),    intent(out) :: pflx(nlevp, ncol)  ! scattered precip flux at each level
    real(r8),    intent(out) :: cme(nlev, ncol)    ! cmf condensation - evaporation
    real(r8),    intent(out) :: zdu(nlev, ncol)    ! detraining mass flux
    real(r8),    intent(out) :: rliq(ncol)         ! reserved liquid (not yet in cldliq) for energy integrals
! local
    real(r8)                 :: landfrac(ncol)      ! Land fraction
    real(r8)                 :: pblh(ncol)         ! Planetary boundary layer height 
    real(r8)                 :: tpert(ncol)        ! Thermal temperature excess
    real(r8)                 :: snow(ncol)
    real(r8)                 :: prec(ncol)
    real(r8)                 :: t(nlev, ncol)
    real(r8)                 :: s(nlev, ncol)
    real(r8)                 :: u(nlev, ncol)
    real(r8)                 :: v(nlev, ncol)
    real(r8)                 :: q(ntracer, nlev, ncol)
    real(r8)                 :: ps(ncol)
    real(r8)                 :: phis(ncol)
    real(r8)                 :: pdel(nlev, ncol)
    real(r8)                 :: pmid(nlev, ncol)
    real(r8)                 :: pint(nlevp, ncol)
    real(r8)                 :: zm(nlev, ncol)
    real(r8)                 :: zi(nlevp, ncol)
    real(r8)                 :: ptend_u(nlev, ncol)
    real(r8)                 :: ptend_v(nlev, ncol)
    real(r8)                 :: ptend_s(nlev, ncol)
    real(r8)                 :: ptend_q(ntracer, nlev, ncol)
    real(r8)                 :: mcon(nlevp, ncol)

    select case ( deep_scheme )
    case('off')

    ptend_deep_convection%tend_u%f = 0._r8
    ptend_deep_convection%tend_v%f = 0._r8
    ptend_deep_convection%tend_s%f = 0._r8
    ptend_deep_convection%tend_q%f = 0._r8

    dlf  = 0._r8
    pflx = 0._r8
    cme  = 0._r8
    zdu  = 0._r8
    rliq = 0._r8
    pstate_cam%updraft_mass_flux%f(:,1:ncol)     = 0._r8


    case('ZM')
 
    landfrac(1:ncol)   = pstate%landfrac_at_pc_surface%f(1:ncol)
    pblh(1:ncol)       = pstate_cam%pbl_pblh_at_pc_surface%f(1:ncol)
    tpert(1:ncol)      = pstate_cam%pbl_tpert_at_pc_surface%f(1:ncol) 
    u(:,1:ncol)        = pstate%u_wind_at_pc_full_level%f(:,1:ncol)
    v(:,1:ncol)        = pstate%v_wind_at_pc_full_level%f(:,1:ncol) 
    t(:,1:ncol)        = pstate%temp_at_pc_full_level%f(:,1:ncol)
    s(:,1:ncol)        = pstate%static_energy_at_pc_full_level%f(:,1:ncol)
    q(:,:,1:ncol)      = pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)
    pmid(:,1:ncol)     = pstate%pressure_at_pc_full_level%f(:,1:ncol)
    pint(:,1:ncol)     = pstate%pressure_at_pc_face_level%f(:,1:ncol)
    pdel(:,1:ncol)     = pstate%delp_at_pc_full_level%f(:,1:ncol)
    ps(1:ncol)         = pstate%pressure_at_pc_surface%f(1:ncol)
    zm(:,1:ncol)       = pstate%z_at_pc_full_level%f(:,1:ncol)
    zi(:,1:ncol)       = pstate%z_at_pc_face_level%f(:,1:ncol)
    phis(1:ncol)       = pstate%geop_at_pc_surface%f(1:ncol)

    call zm_conv_tend( pblh    , mcon    , cme     , tpert   , dlf    , &
                       pflx    , zdu     , landfrac, rliq    , dt     , &
                       t       , q       , u       , v       , ps     , &
                       s       , pmid    , pint    , pdel    ,          &
                       zm      , zi      , phis    ,                    &
                       pstate_cam%cumulus_cldtop%f(1:ncol)   ,                    &
                       pstate_cam%cumulus_cldbot%f(1:ncol)   ,                    &
                       pstate_cam%dp_icwmr_at_pc_full_level%f(:,1:ncol)   , &
                       pstate_cam%aerosol_fracis%f(:,:,1:ncol)   ,          &
                       pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol) , &
                       prec    , snow    , pstate_cam%rprddp_at_pc_full_level%f(:,1:ncol), &
                       evapcdp , ptend_u , ptend_v , ptend_q , ptend_s, &
                       flxprec , flxsnow)


    pstate_cam%scalar_snowc_surface%f(1:ncol)    = pstate_cam%scalar_snowc_surface%f(1:ncol) + snow(1:ncol)
    pstate_cam%scalar_precc_dp_surface%f(1:ncol) = pstate_cam%scalar_precc_dp_surface%f(1:ncol) + prec(1:ncol)
    pstate_cam%updraft_mass_flux%f(:,1:ncol)     = mcon(:,1:ncol)

    ! re-initialize ptend
    ptend_deep_convection%tend_u%f = 0._r8
    ptend_deep_convection%tend_v%f = 0._r8
    ptend_deep_convection%tend_s%f = 0._r8
    ptend_deep_convection%tend_q%f = 0._r8

    ptend_deep_convection%tend_u%f(:,1:ncol)            = ptend_u(:,1:ncol)
    ptend_deep_convection%tend_v%f(:,1:ncol)            = ptend_v(:,1:ncol)
    ptend_deep_convection%tend_s%f(:,1:ncol)            = ptend_s(:,1:ncol)
    ptend_deep_convection%tend_q%f(:,:,1:ncol)          = ptend_q(:,:,1:ncol)

    end select
    
    end subroutine convect_deep_tend


    subroutine convect_deep_tend_2(ncol)
    use grist_zm_conv,                      only: zm_conv_tend_2
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam, ptend_deep_convection_2
 
! io
    integer(i4), intent(in) :: ncol
! local
    real(r8)                 :: pdeldry(nlev, ncol)
    real(r8)                 :: ptend_q(ntracer, nlev, ncol)
 
    if ( deep_scheme .eq. 'ZM' ) then
    !---------------LiXH Test----------------
    ! This part may be error ! Should use dry pressure ??
    pdeldry(:,1:ncol) = pstate%delp_at_pc_full_level%f(:,1:ncol)
    !---------------LiXH Test----------------

    call zm_conv_tend_2(pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol),  &
                        pstate_cam%aerosol_fracis%f(:,:,1:ncol),            &
                        pdeldry, ptend_q)

    ptend_deep_convection_2%tend_q%f(:,:,1:ncol)          = ptend_q(:,:,1:ncol)

    end if

    end subroutine convect_deep_tend_2

 end module grist_deep_convection
