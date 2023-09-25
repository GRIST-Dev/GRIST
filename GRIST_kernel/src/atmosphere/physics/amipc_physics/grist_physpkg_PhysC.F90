!======================================================
!
!  Created by LiXiaohan on 19/5/13.
!  Provides the interface to GRIST physics package
!======================================================

 module grist_physpkg_PhysC

    use grist_domain_types,                 only: global_domain
    use grist_constants,                    only: i4, r8, deg2rad, one
    use grist_nml_module,                   only: physpkg, working_mode, sub_physpkg, &
                                                  nlev, nlevp, ntracer, test_real_case
    use grist_hpe_constants,                only: eta_full
    use grist_physics_data_structure,       only: pstate, ptend_f3
    use grist_PhysC_data_structure,          only: pstate_cam
    use grist_physics_idealized_package,    only: grist_idealized_physics_dcmip2016

    implicit none
    private
    public  :: grist_physpkg_PhysC_init,    &
               grist_physpkgbc_PhysC_run,   &
               grist_physpkgac_PhysC_run,   &
               grist_physpkg_PhysC_final

contains
    subroutine grist_physpkg_PhysC_init(ncell, dtime, istep_begin, dxmean, lat, lon)
    use grist_wv_saturation,                only: wv_sat_init
    use grist_vertical_diffusion,           only: vertical_diffusion_init
    use grist_shallow_convection,           only: convect_shallow_init
    use grist_deep_convection,              only: convect_deep_init
    use cloud_fraction,                     only: cldfrc_init
    use grist_macrop,                       only: macrop_driver_init
    use grist_microp,                       only: microp_driver_init
    use grist_rad_constituents,             only: rad_cnst_init
    use grist_physics_update,               only: phys_time_level_init
    use grist_microp_aero,                  only: microp_aero_init
    use grist_radiation,                    only: radiation_defaultopts,        &
                                                  radiation_setopts,            &
                                                  radiation_init
    use cloud_rad_props,                    only: cloud_rad_props_init
    use grist_cloud_diagnostics,            only: cloud_diagnostics_init
    use rad_solar_var,                      only: rad_solar_var_init
    use grist_nml_module,                   only: start_ymd
    use grist_zenith,                       only: orb_params
    use solar_data,                         only: solar_data_init
    use ghg_data,                           only: ghg_data_init
    use chem_surfvals,                      only: chem_surfvals_init
    use aer_rad_props,                      only: aer_rad_props_init
    use grist_aerosol_intr,                 only: aerosol_init
    use grist_prescribed_aero,              only: prescribed_aero_init
    use grist_prescribed_ozone,             only: prescribed_ozone_init
    use grist_nml_module,                   only: test_real_case,               &
                                                  levsoil, nlev            !--cheyz
#ifdef USE_NOAHMP
    use grist_lsm_noahmp_init,              only: lsm_noahmp_init
    use grist_lsm_noahmp_resVars,           only: grist_lsm_resVars_construct
#endif

! io
    integer(i4), intent(in)            :: ncell
    real(r8)   , intent(in)            :: dtime
    integer    , intent(in)            :: istep_begin 
    real(r8)   , intent(in)            :: dxmean(ncell)          ! mesh scale
    real(r8)   , intent(in)            :: lat(ncell)             ! longitude in radian
    real(r8)   , intent(in)            :: lon(ncell)             ! latitude in radian  
 
! local
    integer(i4)     :: ncol 
! Radiative heating rate calculation options
    integer :: iradsw        ! freq. of shortwave radiation calc in time steps (positive)
                             ! or hours (negative).  Default: -1
    integer :: iradlw        ! frequency of longwave rad. calc. in time steps (positive)
                             ! or hours (negative).  Default: -1
    integer :: iradae        ! frequency of absorp/emis calc in time steps (positive)
                             ! or hours (negative).  Default: -12
    integer :: irad_always   ! Specifies length of time in timesteps (positive)
                             ! or hours (negative) SW/LW radiation will be run continuously
                             ! from the start of an initial run.  Default: 0
    logical :: spectralflux  ! calculate fluxes (up and down) per band. Default: FALSE
    integer :: start_step
    integer :: year


    ncol = ncell

    !--------initialize orbital parameters----------->
    ! move outside the atm model ? LiXH
    year  = start_ymd/10000
    call orb_params( year, .true.)
    !<-------initialize orbital parameters------------


    if(trim(physpkg) .eq. 'AMIPC_PHYSICS')then
    call radiation_defaultopts(iradsw_out      = iradsw,      &
                               iradlw_out      = iradlw,      &
                               iradae_out      = iradae,      &
                               irad_always_out = irad_always, &
                               spectralflux_out = spectralflux )

    call radiation_setopts( dtime,                            &
                            start_in       = istep_begin,     &
                            iradsw_in      = iradsw,          &
                            iradlw_in      = iradlw,          &
                            iradae_in      = iradae,          &
                            irad_always_in = irad_always,     &
                            spectralflux_in = spectralflux )

    ! read namelist
    call read_phys_nml
    call phys_time_level_init
    call chem_surfvals_init
    call phys_inidat(ncol)
    ! wv_saturation is relatively independent of everything else and
    ! low level, so init it early. Must at least do this before radiation.
    call wv_sat_init

    ! Prescribed tracers
    call ghg_data_init(ncol)
    call prescribed_ozone_init(ncol, dtime, lat, lon)
    call prescribed_aero_init(ncol, dtime, lat, lon)

    ! Initialize rad constituents and their properties
    call rad_cnst_init(ncol)
    call aer_rad_props_init(ncol)
    call cloud_rad_props_init

    ! Initialize some aerosol code
    call aerosol_init(ncol)

    ! solar irradiance data modules
    call solar_data_init(0, dtime)

    call vertical_diffusion_init(ncol)
    call radiation_init(ncol,dtime)          
    call rad_solar_var_init
    call cloud_diagnostics_init(ncol)
    call convect_shallow_init(ncol, dxmean)
    call cldfrc_init(nlev, nlevp, ncol)
    call convect_deep_init(ncol)
    call macrop_driver_init(ncol)
    call microp_aero_init(ncol)
    call microp_driver_init(ncol)
#ifdef USE_NOAHMP
    if(test_real_case)then
    call lsm_noahmp_init (ncol, nlev, levsoil,lat, lon)  !--cheyz
    call grist_lsm_resVars_construct(ncol)
    end if
#endif

    ! CAM: call phys_inidat         *    
    ! CAM: call wv_sat_init         *
    ! CAM: if (cam3_aero_data_on) call cam3_aero_data_init
    ! CAM: call rad_cnst_init       *
    ! CAM: call aer_rad_props_init  *
    ! CAM: call cloud_rad_props_init   * 
    ! CAM: call aerosol_init
    ! CAM: call carma_init
    ! CAM: call solar_data_init     *
    ! CAM: call chem_init
    ! CAM: call prescribed tracers:prescribed_ozone_init, prescribed_ghg_init, &
    !                              prescribed_aero_init, aerodep_flx_init,     &
    !                              aircraft_emit_init, prescribed_volcaero_init
    ! CAM: if (co2_transport()) call co2_init
    ! CAM: if (cam3_ozone_data_on) call cam3_ozone_data_init
    ! CAM: call gw_inti
    ! CAM: call rayleigh_friction_init
    ! CAM: call pbl_utils_init
    ! CAM: call vertical_diffusion_init *
    ! CAM: call tsinti
    ! CAM: call radiation_init          *
    ! CAM: call rad_solar_var_init      *
    ! CAM: call cloud_diagnostics_init  *
    ! CAM: call radheat_init            *
    ! CAM: call convect_shallow_init    *
    ! CAM: call cldfrc_init             *
    ! CAM: call convect_deep_init       *
    ! CAM: call macrop_driver_init      *
    ! CAM: call microp_aero_init        *
    ! CAM: microp_driver_init,          *
    ! CAM: conv_water_init              *
    ! CAM: sslt_rebin_init
    ! CAM: tropopause_init

    else

    ! Other physical parameterization packages will be added here, Lixh
    end if

    end subroutine grist_physpkg_PhysC_init


    subroutine grist_physpkg_PhysC_final
    use grist_vertical_diffusion,           only: end_of_vertical_diffusion
    use grist_wv_saturation,                only: end_of_wv_sat
    use grist_shallow_convection,           only: end_of_convect_shallow
    use grist_deep_convection,              only: end_of_convect_deep
    use grist_macrop,                       only: end_of_macrop
    use grist_microp,                       only: end_of_microp
    use cloud_fraction,                     only: end_of_cldfrc
    use grist_microp_aero,                  only: end_of_microp_aero
    use grist_radiation,                    only: end_of_radiation
    use cloud_rad_props,                    only: end_of_cloud_rad_props
    use solar_data,                         only: end_of_solar_data
    use rad_solar_var,                      only: end_of_rad_solar_var
    use grist_cloud_diagnostics,            only: end_of_cloud_diagnostics
    use ghg_data,                           only: end_of_ghg_data
    use grist_aerosol_intr,                 only: end_of_grist_aerosol

    ! destruct variables in each physical parameterization
    call end_of_vertical_diffusion
    call end_of_wv_sat
    call end_of_convect_shallow
    call end_of_convect_deep
    call end_of_macrop
    call end_of_microp
    call end_of_microp_aero
    call end_of_cldfrc
    call end_of_solar_data
    call end_of_rad_solar_var
    call end_of_radiation
    call end_of_cloud_rad_props
    call end_of_cloud_diagnostics
    call end_of_ghg_data
    call end_of_grist_aerosol
    
    end subroutine grist_physpkg_PhysC_final


    subroutine read_phys_nml
    use grist_vertical_diffusion,           only: read_nml_vertical_diffusion
    use phys_control,                       only: read_nml_phys_ctl
    use grist_wv_saturation,                only: read_nml_wv_sat
    use cloud_fraction,                     only: read_nml_cldfrc
    use grist_zm_conv,                      only: read_nml_zmconv
    use grist_uwshcu,                       only: read_nml_uwshcu
    use grist_physics_data_structure,       only: get_tracer_information
    use grist_macrop,                       only: read_nml_macrop
    use grist_microp,                       only: read_nml_microp
    use grist_rad_constituents,             only: read_nml_rad_cnst
    use modal_aer_opt,                      only: read_nml_modal_aer_opt
    use grist_microp_aero,                  only: read_nml_microp_aero
    use solar_data,                         only: read_nml_solar_data
    use chem_surfvals,                      only: read_nml_chem_surfvals
    use grist_radiation,                    only: read_nml_radiation
    use grist_aerosol_intr,                 only: read_nml_aerosol
    use grist_prescribed_ozone,             only: read_nml_prescribed_ozone
    use grist_prescribed_aero,              only: read_nml_prescribed_aero

! local
    character(len=300) :: filename

    filename = 'grist_amipc_phys.nml'
    call get_tracer_information(physpkg,sub_physpkg)
    call read_nml_chem_surfvals(filename)
    call read_nml_phys_ctl(filename)
    call read_nml_wv_sat(filename)
    call read_nml_radiation(filename)
    call read_nml_macrop(filename)
    call read_nml_microp(filename)
    call read_nml_microp_aero(filename)
    call read_nml_cldfrc(filename)
    call read_nml_zmconv(filename)
    call read_nml_uwshcu(filename)
    call read_nml_rad_cnst(filename)
    call read_nml_modal_aer_opt(filename)
    call read_nml_solar_data(filename)
    call read_nml_vertical_diffusion(filename)
    call read_nml_aerosol(filename)
    call read_nml_prescribed_ozone(filename)
    call read_nml_prescribed_aero(filename)


!   call spmd_utils_readnl(nlfilename)
!   call physconst_readnl(nlfilename)
!   call chem_surfvals_readnl(nlfilename)       *
!   call phys_ctl_readnl(nlfilename)            *
!   call wv_sat_readnl(nlfilename)              *
!   call ref_pres_readnl(nlfilename)
!   call cam3_aero_data_readnl(nlfilename)
!   call cam3_ozone_data_readnl(nlfilename)
!   call macrop_driver_readnl(nlfilename)       *
!   call microp_driver_readnl(nlfilename)       *
!   call microp_aero_readnl(nlfilename)         *
!   call cldfrc_readnl(nlfilename)              *
!   call zmconv_readnl(nlfilename)              *
!   call cldwat_readnl(nlfilename)
!   call hkconv_readnl(nlfilename)
!   call uwshcu_readnl(nlfilename)              *
!   call cld_sediment_readnl(nlfilename)
!   call gw_drag_readnl(nlfilename)
!   call phys_debug_readnl(nlfilename)
!   call rad_cnst_readnl(nlfilename)            *
!   call rad_data_readnl(nlfilename)
!   call modal_aer_opt_readnl(nlfilename)       *
!   call chem_readnl(nlfilename)
!   call prescribed_volcaero_readnl(nlfilename)
!   call solar_data_readnl(nlfilename)          *
!   call carma_readnl(nlfilename)
!   call tropopause_readnl(nlfilename)
!   call aoa_tracers_readnl(nlfilename)
!   call aerodep_flx_readnl(nlfilename)
!   call prescribed_ozone_readnl(nlfilename)    *
!   call prescribed_aero_readnl(nlfilename)     *
!   call prescribed_ghg_readnl(nlfilename)
!   call co2_cycle_readnl(nlfilename)
!   call aircraft_emit_readnl(nlfilename)
!   call cospsimulator_intr_readnl(nlfilename)
!   call sat_hist_readnl(nlfilename, hfilename_spec, mfilt, fincl, nhtfrq, avgflag_pertape)
!   call diag_readnl(nlfilename)
!#if (defined WACCM_PHYS)
!   call waccm_forcing_readnl(nlfilename)
!#endif
!   call vd_readnl(nlfilename)                  *
!#if ( defined OFFLINE_DYN )
!   call metdata_readnl(nlfilename)
!#endif

    end subroutine read_phys_nml


    subroutine grist_physpkgbc_PhysC_run(ncell, lat, lon, istep, dtime, model_timestep)
    !use grist_lnd_static_vars_module,       only: scalar_static_landmask_at_pc
    use lixh_test_simple_physics,           only: large_scale_prcp,             &
                                                  surface_flux,                 &
                                                  grist_idealized_physics_mitc
! test
    use grist_amp_surface_flux_module , only: grist_amp_surface_flux_run
    use grist_PhysC_data_structure,          only: ptend_vertical_diffusion,     &
                                                  ptend_deep_convection,        &
                                                  ptend_deep_convection_2,      &
                                                  ptend_shallow_convection,     &
                                                  ptend_macrophysics,           &
                                                  ptend_microphysics,           &
                                                  ptend_radiation,              &
                                                  ptend_dry_adjustment
    use grist_physics_update,               only: phys_update,                  &
                                                  phys_time_level_update
    use grist_vertical_diffusion,           only: vertical_diffusion_tend
    use grist_shallow_convection,           only: convect_shallow_tend
    use grist_deep_convection,              only: convect_deep_tend,            &
                                                  convect_deep_tend_2
    use grist_dry_adjustment,               only: dadadj
    use grist_macrop,                       only: macrop_driver_tend
    use grist_microp,                       only: microp_driver_tend
    use grist_microp_aero,                  only: microp_aero_driver
    use grist_radiation,                    only: radiation_tend
    use grist_aerosol_intr,                 only: aerosol_wet_intr
    use grist_mpi
    use grist_PhysC_diag,                   only: grist_PhysC_diag_run
    use grist_cloud_diagnostics,            only: cloud_diagnostics_calc,       &
                                                  phys_diagnostics_calc
    use grist_phys_state_check,             only: phys_state_check

! io
    integer(i4), intent(in)            :: ncell 
    real(r8)   , intent(in)            :: lat(:)             ! longitude in radian
    real(r8)   , intent(in)            :: lon(:)             ! latitude in radian  
    integer(i4), intent(in)            :: istep
    real(r8)   , intent(in)            :: dtime
    real(r8)   , intent(in)            :: model_timestep
! local
    integer(i4)                        :: i
    integer(i4)                        :: ncol
    real(r8), allocatable              :: net_flx(:)  
    real(r8), allocatable              :: zdu(:,:)            ! detraining mass flux from deep convection
    real(r8), allocatable              :: cmfcme(:,:)         ! cmf condensation - evaporation
    real(r8), allocatable              :: cmfmc2(:,:)         ! Updraft mass flux by shallow convection [ kg/s/m2 ]
    real(r8), allocatable              :: dlf(:,:)            ! Detraining cld H20 from shallow + deep convections
    real(r8), allocatable              :: dlf2(:,:)           ! Detraining cld H20 from shallow convections
    real(r8), allocatable              :: pflx(:,:)           ! Conv rain flux thru out btm of lev
    real(r8), allocatable              :: rliq(:)             ! vertical integral of liquid not yet in q(ixcldliq)
    real(r8), allocatable              :: rliq2(:)            ! vertical integral of liquid from shallow scheme
    real(r8), allocatable              :: det_s(:)            ! vertical integral of detrained static energy from ice
    real(r8), allocatable              :: det_ice(:)          ! vertical integral of detrained ice
!
! for testing simple physics   here
!
    real(r8)               :: local_uuu(ncell,nlev)
    real(r8)               :: local_vvv(ncell,nlev)
    real(r8)               :: local_ttt(ncell,nlev)
    real(r8)               :: local_qqq(ncell,nlev)
    real(r8)               :: local_pmid(ncell,nlev)
    real(r8)               :: local_pint(ncell,nlev+1)
    real(r8)               :: local_pdel(ncell,nlev)
    real(r8)               :: local_rpdel(ncell,nlev)
    real(r8)               :: local_dudt(ncell,nlev)
    real(r8)               :: local_dvdt(ncell,nlev)
    real(r8)               :: local_dtdt(ncell,nlev)
    real(r8)               :: local_dqdt(ncell,nlev)

    ncol = ncell

    if(physpkg .eq. 'simple_physics' )then
        call grist_idealized_physics_mitc(ncol, nlev, dtime, lat,                 &
                                pstate%temp_at_pc_full_level%f(:,1:ncol),                  &
                                pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol),         &
                                pstate%u_wind_at_pc_full_level%f(:,1:ncol),                &
                                pstate%v_wind_at_pc_full_level%f(:,1:ncol),                &
                                pstate%pressure_at_pc_full_level%f(:,1:ncol),              &
                                pstate%pressure_at_pc_face_level%f(:,1:ncol),              &
                                eta_full,                                                  &
                                pstate%pressure_at_pc_surface%f(1:ncol),                   &
                                pstate%scalar_precl_surface%f(1:ncol),                     &
                                working_mode)

    else if(trim(physpkg) .eq. 'AMIPC_PHYSICS'.and. trim(sub_physpkg).eq.'DCMIP2016-TC' )then
        ! This produces exactly the same results as tested  outside given the same config
        if(mpi_rank().eq.0) print*,"test  dcmip2016-tc inside CAM5PHYSICS"
        local_uuu(1:ncol,1:nlev)   = transpose(pstate%u_wind_at_pc_full_level%f(1:nlev,1:ncol))
        local_vvv(1:ncol,1:nlev)   = transpose(pstate%v_wind_at_pc_full_level%f(1:nlev,1:ncol))
        local_ttt(1:ncol,1:nlev)   = transpose(pstate%temp_at_pc_full_level%f(1:nlev,1:ncol))
        local_qqq(1:ncol,1:nlev)   = transpose(pstate%tracer_mxrt_at_pc_full_level%f(1,1:nlev,1:ncol))
        local_pmid(1:ncol,1:nlev)  = transpose(pstate%pressure_at_pc_full_level%f(1:nlev,1:ncol))
        local_pint(1:ncol,1:nlevp) = transpose(pstate%pressure_at_pc_face_level%f(1:nlevp,1:ncol))
        local_pdel(1:ncol,1:nlev)  = transpose(pstate%delp_at_pc_full_level%f(1:nlev,1:ncol))
        local_rpdel(1:ncol,1:nlev) = one/local_pdel(1:ncol,1:nlev)
   
        call grist_idealized_physics_dcmip2016(ncol, nlev, dtime, lat , &
                                               local_ttt, local_qqq, local_uuu, local_vvv, &
                                               local_pmid, local_pint, local_pdel, local_rpdel, pstate%pressure_at_pc_surface%f(1:ncol), &
                                               pstate%scalar_precl_surface%f(1:ncol), 0, .true., .false.              , &
                                               local_dudt, local_dvdt, local_dtdt, local_dqdt)

        ptend_f3%tend_u_wind_at_pc_full_level%f(1:nlev,1:ncol)        = transpose(local_dudt(1:ncol,1:nlev))
        ptend_f3%tend_v_wind_at_pc_full_level%f(1:nlev,1:ncol)        = transpose(local_dvdt(1:ncol,1:nlev))
        ptend_f3%tend_temp_at_pc_full_level%f(1:nlev,1:ncol)          = transpose(local_dtdt(1:ncol,1:nlev))
        ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1,1:nlev,1:ncol) = transpose(local_dqdt(1:ncol,1:nlev))
                           
    else if(trim(physpkg) .eq. 'AMIPC_PHYSICS'.and.trim(sub_physpkg) .eq. 'none')then
        allocate(zdu     (nlev,  ncol));zdu    = 0._r8
        allocate(cmfcme  (nlev,  ncol));cmfcme = 0._r8
        allocate(cmfmc2  (nlevp, ncol));cmfmc2 = 0._r8
        allocate(pflx    (nlevp, ncol));pflx   = 0._r8
        allocate(dlf     (nlev,  ncol));dlf    = 0._r8
        allocate(dlf2    (nlev,  ncol));dlf2   = 0._r8
        allocate(rliq    (       ncol));rliq   = 0._r8
        allocate(rliq2   (       ncol));rliq2  = 0._r8
        allocate(det_s   (       ncol));det_s  = 0._r8
        allocate(det_ice (       ncol));det_ice= 0._r8
        allocate(net_flx (       ncol));net_flx= 0._r8

        pstate_cam%aerosol_fracis%f            = 1._r8
        ! Advance time information
        call phys_timestep_init(ncol, istep, model_timestep, lat, lon)
        call phys_state_check('phys_timestep_init', ncol, lat, lon)

        ! Dry adiabatic adjustment
        call dadadj(ncol, dtime)
        call phys_update('dadadj',dtime, ncol, ptend_dry_adjustment)
        call phys_state_check('dadadj', ncol, lat, lon)

        ! Moist convection
        call convect_deep_tend(ncol, cmfcme, dlf, pflx, zdu, rliq, dtime)
        call phys_update('deep convection',dtime, ncol, ptend_deep_convection)
        call phys_state_check('deep convection', ncol, lat, lon)

        call convect_shallow_tend(ncol, istep, dtime, cmfmc2,      &
                                  dlf, dlf2, rliq, rliq2 )
        call phys_update('shallow convection',dtime, ncol, ptend_shallow_convection)
        call phys_state_check('shallow convection', ncol, lat, lon)

        ! Micro- and Macrophysics
        call macrop_driver_tend(ncol, dtime, istep, dlf, dlf2,     &
                                cmfmc2, det_s, det_ice)
        call phys_update('macrophysics',dtime, ncol, ptend_macrophysics)
        call phys_state_check('macrophysics', ncol, lat, lon)

        call microp_aero_driver(ncol, dtime)

        call microp_driver_tend(ncol, dtime)
        call phys_update('microphysics',dtime, ncol, ptend_microphysics)
        call phys_state_check('microphysics', ncol, lat, lon)

        !  Aerosol wet chemistry determines scavenging fractions, and transformations
        !  Then do convective transport of all trace species except water vapor and
        !     cloud liquid and ice (we needed to do the scavenging first
        !     to determine the interstitial fraction) 
        call aerosol_wet_intr(ncol)

        call convect_deep_tend_2(ncol)
        call phys_update('deep convection 2',dtime, ncol, ptend_deep_convection_2)
        call phys_state_check('deep convection 2', ncol, lat, lon)
        
        ! Moist physical parameteriztions complete
        ! cloud and other physics diagnostics
        call cloud_diagnostics_calc(ncol)

        call phys_diagnostics_calc(ncol)

        call grist_PhysC_diag_run(ncol)
        
        ! Raidation computations
        call radiation_tend(ncol, istep, model_timestep, lon, lat, net_flx)
        call phys_update('radiation',dtime, ncol, ptend_radiation)
        call phys_state_check('radiation', ncol, lat, lon)

        deallocate(zdu    )
        deallocate(cmfcme )
        deallocate(cmfmc2 )
        deallocate(pflx   )
        deallocate(dlf    )
        deallocate(dlf2   )
        deallocate(rliq   )
        deallocate(rliq2  )
        deallocate(det_s  )
        deallocate(det_ice)
        deallocate(net_flx)

    end if

    return
    end subroutine grist_physpkgbc_PhysC_run

   subroutine grist_physpkgac_PhysC_run(ncell, lat, lon, istep, dtime, model_timestep)
    use grist_PhysC_data_structure,         only: ptend_vertical_diffusion
    use grist_physics_update,               only: phys_update,                  &
                                                  phys_time_level_update
    use grist_vertical_diffusion,           only: vertical_diffusion_tend
    use grist_aerosol_intr,                 only: aerosol_drydep_intr
    use grist_shallow_convection,           only: trigger_function_xie
    use grist_phys_state_check,             only: phys_state_check
    use grist_mpi

! io
    integer(i4), intent(in)            :: ncell 
    real(r8)   , intent(in)            :: lat(:)             ! longitude in radian
    real(r8)   , intent(in)            :: lon(:)             ! latitude in radian  
    integer(i4), intent(in)            :: istep
    real(r8)   , intent(in)            :: dtime
    real(r8)   , intent(in)            :: model_timestep
! local
    integer(i4)                        :: i
    integer(i4)                        :: ncol
 
    ncol = ncell
   
    if(trim(physpkg) .eq. 'AMIPC_PHYSICS'.and.trim(sub_physpkg) .eq. 'none')then
    ! Vertical diffusion/pbl calculation
    call vertical_diffusion_tend(ncol, dtime)
    call phys_update('vertical diffusion',dtime, ncol, ptend_vertical_diffusion)
    call phys_state_check('vertical diffusion', ncol, lat, lon)
    !  aerosol dry deposition processes
    call aerosol_drydep_intr(ncol)

    call phys_time_level_update

    !----------------->
    ! trigger function for deep convection, LiXH, 2021/10/14.
    call trigger_function_xie(ncol)
    !<-----------------
    end if

    return
    end subroutine grist_physpkgac_PhysC_run


    subroutine phys_timestep_init(ncol, istep, dtime, lat, lon)
    use ghg_data,                           only: ghg_data_timestep_init
    use grist_prescribed_ozone,             only: prescribed_ozone_adv
    use solar_data,                         only: solar_data_advance
    use grist_prescribed_aero,              only: prescribed_aero_adv

! io
    integer(i4), intent(in) :: ncol
    integer(i4), intent(in) :: istep
    real(r8)   , intent(in) :: dtime
    real(r8)   , intent(in) :: lat(ncol)             ! longitude in radian
    real(r8)   , intent(in) :: lon(ncol)             ! latitude in radian  
 
    ! Chemistry surface values
    !call chem_surfvals_set()

    ! Solar irradiance
    call solar_data_advance(istep, dtime)

    ! Time interpolate for chemistry
    !call chem_timestep_init()

    ! Prescribed tracers
    call prescribed_ozone_adv(ncol, istep, dtime, lon, lat)
    !call prescribed_ghg_adv()
    call prescribed_aero_adv(ncol, istep, dtime, lon, lat)
    !call aircraft_emit_adv()
    !call prescribed_volcaero_adv()

    ! prescribed aerosol deposition fluxes
    !call aerodep_flx_adv()                          !need?

    ! Time interpolate data models of gasses
    call ghg_data_timestep_init(ncol, lat)

    ! Time interpolate for vertical diffusion upper boundary condition
    !call vertical_diffusion_ts_init()

    ! age of air tracers
    !call aoa_tracers_timestep_init()

    end subroutine phys_timestep_init


    subroutine phys_inidat(ncol)
    use grist_physics_update,               only: ptimelevels
! io
    integer(i4), intent(in) :: ncol

! surface and land model
    allocate(pstate%snowhice_at_pc_surface%f(ncol))
    if(.not.allocated(pstate%icefrac_at_pc_surface%f))  allocate(pstate%icefrac_at_pc_surface%f(ncol))


! surface flux of precipitation and snow
    allocate(pstate_cam%str_prec_surface%f(ncol))
    allocate(pstate_cam%str_snow_surface%f(ncol))
    allocate(pstate_cam%sed_prec_surface%f(ncol))
    allocate(pstate_cam%sed_snow_surface%f(ncol))
    allocate(pstate_cam%pcw_prec_surface%f(ncol))
    allocate(pstate_cam%pcw_snow_surface%f(ncol))

! atm out
    allocate(pstate%atm_out_sols_at_pc_surface%f(ncol))
    allocate(pstate%atm_out_soll_at_pc_surface%f(ncol))
    allocate(pstate%atm_out_solsd_at_pc_surface%f(ncol))
    allocate(pstate%atm_out_solld_at_pc_surface%f(ncol))

! gwd
    allocate(pstate_cam%sgh30_at_pc_surface%f(ncol))


! initialize position for output
    pstate%snowhice_at_pc_surface%pos        = 0
    pstate%icefrac_at_pc_surface%pos         = 0

    pstate_cam%str_prec_surface%pos              = 0
    pstate_cam%str_snow_surface%pos              = 0
    pstate_cam%sed_prec_surface%pos              = 0
    pstate_cam%sed_snow_surface%pos              = 0
    pstate_cam%pcw_prec_surface%pos              = 0
    pstate_cam%pcw_snow_surface%pos              = 0

    pstate%atm_out_sols_at_pc_surface%pos    = 0
    pstate%atm_out_soll_at_pc_surface%pos    = 0
    pstate%atm_out_solsd_at_pc_surface%pos   = 0
    pstate%atm_out_solld_at_pc_surface%pos   = 0

    pstate_cam%sgh30_at_pc_surface%pos   = 0

! initialize to 0, if we have initial data, set to initial value.
    pstate%snowhice_at_pc_surface%f          = 0._r8
    pstate%icefrac_at_pc_surface%f           = 0._r8

    pstate_cam%str_prec_surface%f                = 0._r8
    pstate_cam%str_snow_surface%f                = 0._r8
    pstate_cam%sed_prec_surface%f                = 0._r8
    pstate_cam%sed_snow_surface%f                = 0._r8
    pstate_cam%pcw_prec_surface%f                = 0._r8
    pstate_cam%pcw_snow_surface%f                = 0._r8

    pstate%atm_out_sols_at_pc_surface%f      = 0._r8
    pstate%atm_out_soll_at_pc_surface%f      = 0._r8
    pstate%atm_out_solsd_at_pc_surface%f     = 0._r8
    pstate%atm_out_solld_at_pc_surface%f     = 0._r8

    pstate_cam%sgh30_at_pc_surface%f     = 0._r8

    end subroutine phys_inidat

 end module grist_physpkg_PhysC
