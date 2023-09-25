!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: This is the calling interface to CRM-style full physics, as
!              a template for any of such type (e.g., WRF). The
!              full physpkg depends on the pstate/ptend_f3 structure while
!              simple physics does not. But when we call simple-physics inside 
!              a full physics package, it must (tested) produce the same 
!              results as we directly call the simple-physics interface given 
!              the same configuration. This ensures the MPI-behavior of full
!              physpkg is as same as simple-physics, so seq-grist does not need
!              full physpkg currently, and we have no need to refactor the simple 
!              physics interface.
!              IMP to know: WRF's nlev is our nlevp, we should use a different
!              name to indicate the vertical level inside wrf physics: nLevel
!
! Revision history:
!            1. AMIPW physics (PhysW) is entering from here
!----------------------------------------------------------------------------

 module grist_dtp_interface_full_physpkg2

! Infrastructure
   use grist_mpi
   use grist_constants,                  only: r8, i4, pi, rdry, rvap, cp, p00, omega, gravity, rearth, zero, one, half, t00, fillvalue, ptfactor,  rad2deg
   use grist_nml_module,                 only: physpkg, sub_physpkg, working_mode, nh_dynamics, doAquaPlanet, levsoil, test_real_case, start_ymd, start_tod, &
                                               ptendSubDiabPhys, physics_use_som, nspecies
   use grist_domain_types,               only: global_domain
   use grist_math_module,                only: convert_vector_sph2cart
   use grist_time_manager,               only: get_curr_calday
   use grist_zenith,                     only: zenith, orb_params
! Dynamics
   use grist_hpe_constants,              only: deta_face, deta_full
   use grist_dycore_gcd_recon_module_2d, only: vector_recon_perot_edge2cell_uv_2d
! Dtp, just use some data
   use grist_dtp_interface_module,       only: scalar_U_wind_at_pc_full_level, &
                                               scalar_V_wind_at_pc_full_level, &
                                               scalar_qqq_at_pc_full_level
! Physics
   use grist_physics_data_structure,     only: pstate
   use grist_physics_idealized_package,  only: grist_idealized_physics_kessler_klemp15
#ifdef AMIPW_PHYSICS
   use grist_PhysW_pkg,                only: grist_PhysW_pkg_run_bc, &
                                               grist_PhysW_pkg_run_ac
   use grist_PhysW_data_structure,       only: pstate_wrf, psurf_wrf, ptend_wrf
   use grist_PhysW_nml_module,           only: lsm_scheme=>PhysW_lm_scheme
#endif
! Datam
   use grist_datam_static_data_module,   only: staticData_sst_at_pc_surface  , &
                                               staticData_sic_at_pc_surface  , &
                                               staticData_var2d_at_pc_surface, &
                                               staticData_oc12d_at_pc_surface, &
                                               staticData_oa2d1_at_pc_surface, &
                                               staticData_oa2d2_at_pc_surface, &
                                               staticData_oa2d3_at_pc_surface, &
                                               staticData_oa2d4_at_pc_surface, &
                                               staticData_ol2d1_at_pc_surface, &
                                               staticData_ol2d2_at_pc_surface, &
                                               staticData_ol2d3_at_pc_surface, &
                                               staticData_ol2d4_at_pc_surface
   use grist_datam_initial_data_module,  only: initialData_skintemp_at_pc_surface
! Amp
   use grist_prescribe_seaice_module,    only: sicetemp_at_pc_surface, grist_seaice_init, grist_seaice_run
   use grist_amp_shr_flux_mod,           only: shr_flux_atmOcn

   implicit none
   private
   
   public :: grist_full_physpkg2_init, &
             grist_full_physpkg2_final
#ifdef AMIPW_PHYSICS
   public :: grist_full_physpkg2_PhysW_run
#endif

  contains
!---------------------------------------------------------
! initialization/finalization associated with physpkg
!---------------------------------------------------------

  subroutine  grist_full_physpkg2_init(local_block, dtime, ncell, nlevel, nspecies)
   use grist_domain_types,       only: block_structure
#ifdef AMIPW_PHYSICS
   use grist_PhysW_nml_module,   only: set_PhysW_nml
   use grist_PhysW_pkg,        only: grist_PhysW_pkg_init
   use grist_lnd_static_vars_module,     only: grist_lnd_static_vars_construct, &
                                               scalar_static_landmask_at_pc
#endif
   use grist_slab_ocean_module,  only: grist_prt_oml_init

   type(block_structure) ,target, intent(inout) :: local_block
   real(r8), intent(in)                         :: dtime
   integer(i4),  intent(in)                     :: ncell, nlevel, nspecies
   integer(i4)                                  :: year, icell

!-----------------------------------------------------------------------
! normally, a column model only needs halo(1)
!-----------------------------------------------------------------------

#ifdef AMIPW_PHYSICS

    call set_PhysW_nml

!-----------------------------------------------------------------------
! land data
!-----------------------------------------------------------------------

    local_block%full_domain%nv = local_block%full_domain%nv_halo(1)
    call grist_lnd_static_vars_construct(local_block%full_domain)
    local_block%full_domain%nv = local_block%full_domain%nv_compute

    pstate%landfrac_at_pc_surface%f(1:ncell) = scalar_static_landmask_at_pc%f(1:ncell)

!-----------------------------------------------------------------------
! physpkg init
!-----------------------------------------------------------------------

!    if(trim(physpkg).eq.'AMIPW_PHYSICS')then
     call grist_PhysW_pkg_init(ncell, local_block%full_domain%nv_full, nlevel, nspecies, dtime,&
                                                        local_block%full_domain%vtx_lat(1:ncell),&
                                                        local_block%full_domain%vtx_lon(1:ncell))
     if(mpi_rank()==0) print*,"Sucessfully initiate AMIPW model physics"
!    end  if

    pstate_wrf%xlat(1:ncell,1)  =  local_block%full_domain%vtx_lat(1:ncell)*rad2deg
    pstate_wrf%xlong(1:ncell,1) =  local_block%full_domain%vtx_lon(1:ncell)*rad2deg

    do icell = 1, ncell
       pstate_wrf%dxmean(icell)   = local_block%full_domain%vtxCellLeng(icell)*rearth
       pstate_wrf%coriolis(icell) = 2._r8*omega*sin(local_block%full_domain%vtx_lat(icell))
    end do

!-----------------------------------------------------------------------
! added for noahmp
!-----------------------------------------------------------------------

    where(scalar_static_landmask_at_pc%f(1:ncell).le.0.5_r8) pstate%ocnfrac_at_pc_surface%f(1:ncell)  = 1._r8
    where(scalar_static_landmask_at_pc%f(1:ncell).gt.0.5_r8) pstate%ocnfrac_at_pc_surface%f(1:ncell)  = 0._r8

!-----------------------------------------------------------------------
! init all ts and sst as sst data (real sst data + analytic sst function)
!-----------------------------------------------------------------------

     pstate%ts_at_pc_surface%f(1:ncell)       = staticData_sst_at_pc_surface%f(1:ncell) ! default as SST
#ifdef USE_NOAHMP
    if(test_real_case)then ! default as that in lsm init data
       pstate%ts_at_pc_surface%f(1:ncell)     = initialData_skintemp_at_pc_surface%f(1:ncell)
!-----------------------------------------------------------------------
! static gwdo
!-----------------------------------------------------------------------
       psurf_wrf%var2d(1:ncell,1) = staticData_var2d_at_pc_surface%f(1:ncell)
       psurf_wrf%oc12d(1:ncell,1) = staticData_oc12d_at_pc_surface%f(1:ncell)
       psurf_wrf%oa2d1(1:ncell,1) = staticData_oa2d1_at_pc_surface%f(1:ncell)
       psurf_wrf%oa2d2(1:ncell,1) = staticData_oa2d2_at_pc_surface%f(1:ncell)
       psurf_wrf%oa2d3(1:ncell,1) = staticData_oa2d3_at_pc_surface%f(1:ncell)
       psurf_wrf%oa2d4(1:ncell,1) = staticData_oa2d4_at_pc_surface%f(1:ncell)
       psurf_wrf%ol2d1(1:ncell,1) = staticData_ol2d1_at_pc_surface%f(1:ncell)
       psurf_wrf%ol2d2(1:ncell,1) = staticData_ol2d2_at_pc_surface%f(1:ncell)
       psurf_wrf%ol2d3(1:ncell,1) = staticData_ol2d3_at_pc_surface%f(1:ncell)
       psurf_wrf%ol2d4(1:ncell,1) = staticData_ol2d4_at_pc_surface%f(1:ncell)
    end if
#endif

    where(pstate%landfrac_at_pc_surface%f(1:ncell).le.0.5_r8) pstate%ts_at_pc_surface%f(1:ncell) = staticData_sst_at_pc_surface%f(1:ncell)
    pstate%sst_at_pc_surface%f(1:ncell)  = staticData_sst_at_pc_surface%f(1:ncell)

    if(trim(working_mode).eq.'amipw') psurf_wrf%tsk (1:ncell,1)  = pstate%ts_at_pc_surface%f(1:ncell)

    year  = start_ymd/10000
    call orb_params( year, .true.)

!-----------------------------------------------------------------------
! slab ocean model, init tsk using SST data
!-----------------------------------------------------------------------
    if(physics_use_som) call grist_prt_oml_init(50._r8, pstate%ts_at_pc_surface%f(1:ncell), local_block%full_domain%nv_halo(1))

    if(mpi_rank()==0) print*, "Sucessfully initiate AMIPW Physics (PhysW) in GRIST"

#endif

   return
  end subroutine  grist_full_physpkg2_init

  subroutine  grist_full_physpkg2_final
! Add physics package compiling option here

#ifdef AMIPW_PHYSICS
   use grist_PhysW_pkg,      only: grist_PhysW_pkg_final
   use grist_slab_ocean_module,only: grist_prt_oml_final
#endif

#ifdef AMIPW_PHYSICS
    if(trim(physpkg).eq.'AMIPW_PHYSICS')then
      call grist_PhysW_pkg_final
      if(physics_use_som) call grist_prt_oml_final
      if(mpi_rank()==0) print*,"Sucessfully delete AMIPW model physics"
    end  if
#endif

    return
  end subroutine  grist_full_physpkg2_final

!-----------------------------------------------------------------------
!  AMIPW PHYSICS (PhysW)
!-----------------------------------------------------------------------

#ifdef AMIPW_PHYSICS

  subroutine grist_full_physpkg2_PhysW_run(mesh,ntracer,nmif,nlev,ncell,istep,dtime      , &
                                               scalar_mif_at_pc_full_level              , &
                                               tend_ptdyn_at_pc_full_level              , &
                                               tend_qvdyn_at_pc_full_level              , &
                                               scalar_geopotential_at_pc_full_level     , &
                                               scalar_geopotential_at_pc_face_level     , &
                                               scalar_mpressure_at_pc_full_level        , &
                                               scalar_mpressure_at_pc_face_level        , &
                                               scalar_normal_velocity_at_edge_full_level, &
                                               scalar_potential_temp_at_pc_full_level   , & !thetam
                                               scalar_temp_at_pc_full_level             , &
                                               scalar_tracer_mxrt_at_pc_full_level      , &
                                               scalar_www_at_pc_face_level              , & ! time-averaged
                                               scalar_omega_at_pc_full_level            , & ! time-averaged
                                               scalar_delhp_at_pc_full_level       ) !    , &
                                               !tend_normal_velocity_at_edge_full_level  , &
                                               !tend_potential_temp_at_pc_full_level     , &
                                               !tend_tracer_mxrt_at_pc_full_level        )

    use grist_physics_data_structure,   only: pstate, ptend_f1, ptend_f2, ptend_f3, ptend_rk
    use grist_datam_static_data_module, only: staticData_sst_at_pc_surface
#ifdef USE_NOAHMP
    use grist_lsm_noahmp_driver,        only: grist_noahmp_run
#endif
!    use grist_amp_surface_flux_module , only: grist_amp_surface_flux_run
! io
    use omp_lib
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: ntracer
    integer(i4)          , intent(in)   :: nmif
    integer(i4)          , intent(in)   :: nlev
    integer(i4)          , intent(in)   :: ncell
    integer(i4)          , intent(in)   :: istep
    real(r8)             , intent(in)   :: dtime
    real(r8),allocatable , intent(in)   :: scalar_mif_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: tend_ptdyn_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: tend_qvdyn_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(in)   :: scalar_www_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_omega_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_delhp_at_pc_full_level(:,:)
! local
    real(r8)                            :: lats(ncell),lons(ncell)
    real(r8)                            :: scalar_pt_at_pc_full_level(nlev,mesh%nv_full) ! PT
    integer(i4)                         :: itracer, ie, iv, ilev, icell, icell1, icell2, nlevp
    real(r8)                            :: hfx(ncell,1), qfx(ncell,1), tsk(ncell,1), &
                                           qsfc_lsm(ncell), emiss_lsm(ncell), snowh_lsm(ncell), xice_lsm(ncell)
    real(r8)                            :: lwup(ncell), taux(ncell), tauy(ncell)
    real(r8)                            :: asdir(ncell,1), asdif(ncell,1), aldir(ncell,1), aldif(ncell,1)
    real(r8)                            :: julday, coszrs(ncell)
    real(r8)                            :: scalar_zbot_at_pc_surface(ncell)
    real(r8)                            :: scalar_us_at_pc_surface  (ncell)
    real(r8)                            :: scalar_vs_at_pc_surface  (ncell)

        nlevp = nlev+1
!
! prepare data
!
        if(ncell.ne.mesh%nv_halo(1)) then
           print*, "dim stuff wrong in grist_dtp_interface_physics, model aborts"
           call mpi_abort()
        end if

        do iv = 1, ncell
           lats(iv) = mesh%vtx_lat(iv) ! radian
           lons(iv) = mesh%vtx_lon(iv)
        end do

        call get_curr_calday(start_ymd, start_tod, istep, dtime, julday)
        call zenith(julday, coszrs, ncell, lons, lats,doAquaPlanet)

!----------------------------------------------------------------------
! obtain U&V
!----------------------------------------------------------------------

        mesh%nv = mesh%nv_halo(1)
        call vector_recon_perot_edge2cell_uv_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                                      scalar_U_wind_at_pc_full_level, &
                                                      scalar_V_wind_at_pc_full_level, &
                                                      nlev)

!----------------------------------------------------------------------
! for mxrt-like quantity, from dry to moist
! I am not sure, but I feel that the qv in WRF-physics is dry rather
! than moist mixing ratio; so no conversion here
! If number concentration follows mixing ratio, this line can still
! be generalized
! water vapor is always No.1
!----------------------------------------------------------------------

        scalar_qqq_at_pc_full_level(1:ntracer,1:nlev,1:ncell) = scalar_tracer_mxrt_at_pc_full_level(1:ntracer,1:nlev,1:ncell)
!
! obtain pt from ptm
!
        do iv = 1, ncell
           do ilev = 1, nlev
              scalar_pt_at_pc_full_level(ilev,iv) = scalar_potential_temp_at_pc_full_level(ilev,iv)/(one+ptfactor*scalar_qqq_at_pc_full_level(1,ilev,iv))
           end do
        end do

!----------------------------------------------------------------------
! fill pstate with dynamics state
! note, ptate_wrf's vertical ranges from 1->nlevp,
! while other from  1->nlev
!----------------------------------------------------------------------

        call dtp_interface_fill_pstate(mesh,ntracer,nlev,nlevp,ncell,dtime, lats, &
                                       scalar_U_wind_at_pc_full_level, &
                                       scalar_V_wind_at_pc_full_level, &
                                       scalar_pt_at_pc_full_level    , &
                                       scalar_temp_at_pc_full_level  , &
                                       scalar_qqq_at_pc_full_level   , & ! dry
                                       scalar_mpressure_at_pc_face_level, &
                                       scalar_mpressure_at_pc_full_level, &
                                       scalar_geopotential_at_pc_face_level, &
                                       scalar_geopotential_at_pc_full_level, &
                                       scalar_www_at_pc_face_level,&
                                       scalar_omega_at_pc_full_level,&
                                       scalar_delhp_at_pc_full_level,&
                                       tend_ptdyn_at_pc_full_level,&
                                       tend_qvdyn_at_pc_full_level)

!----------------------------------------------------------------------
! call ocean surface scheme, which relies on pstate of model physics
!----------------------------------------------------------------------

!        call grist_amp_surface_flux_run(ncell,nlev)

!----------------------------------------------------------------------
! call interface for WRF-physics tendency, relies on pstate_wrf and produce 
! ptend_f3; this will can sfclay for surface
!----------------------------------------------------------------------

       call grist_PhysW_pkg_run_bc(ncell,nlevp,nspecies,istep,dtime,rearth*mesh%mean_edt_dist)

!-----------------------------------------------------------------------------------------
! sfclay can do ocean point, so only modify albedo over ocean
!
!IF(doAquaPlanet.or.trim(lsm_scheme).eq.'noahmp')then

!$omp parallel  private(icell) 
!$omp do schedule(dynamic,10) 
       do icell=1, ncell
         if (psurf_wrf%xland(icell,1).ge.1.5_r8)then ! ocean, include ice

!-----------------------------------------------------------------------
! this reset is a must for each step, because seaice will merge sea albedo with ice albedo, 
! we need a pseudo-"fixed" sea albedo here; otherwise, asdir computed at night or not will 
! affect bit solution; because a zero value at this time will be used for asdir input value at next
! time the current code does not have this issue for sea and seaice!
!-----------------------------------------------------------------------

              if(coszrs(icell) .gt. zero)then ! daytime
              ! as in ampflux
              !  psurf_wrf%asdir(icell,1)  = 0.037_r8/(1.1_r8*coszrs(icell)**1.4_r8 + 0.15_r8)
              !  psurf_wrf%asdif(icell,1)  = psurf_wrf%asdir(icell,1)
              !  psurf_wrf%aldir(icell,1)  = 0.06_r8
              !  psurf_wrf%aldif(icell,1)  = 0.06_r8
! CAM3 code, default
                psurf_wrf%aldir(icell,1)  = (.026_r8/(coszrs(icell)**1.7_r8 + .065_r8)) + &
                                          (.15_r8*(coszrs(icell)- 0.10_r8)*(coszrs(icell) - 0.50_r8)*(coszrs(icell) - 1._r8))
                psurf_wrf%asdir(icell,1)  = psurf_wrf%aldir(icell,1)
                psurf_wrf%asdif(icell,1)  = 0.06_r8
                psurf_wrf%aldif(icell,1)  = 0.06_r8
!               psurf_wrf%asdir(icell,1)  = 0.08_r8 
!               psurf_wrf%asdif(icell,1)  = 0.08_r8
!               psurf_wrf%aldir(icell,1)  = 0.08_r8
!               psurf_wrf%aldif(icell,1)  = 0.08_r8
              else
                psurf_wrf%asdir(icell,1)  = zero
                psurf_wrf%asdif(icell,1)  = zero
                psurf_wrf%aldir(icell,1)  = zero
                psurf_wrf%aldif(icell,1)  = zero
             end if
          end if
       end do
!$omp end do nowait
!$omp end parallel 
!ENDIF

!---------------------------------------------------------------------------------
! Use shr_flux_atmOcn OFFLINE over OCEAN for diagnostic purpose for taux, tauy.
! For other vars (lat,lwup,evap,tref,qref,etc), it has been checked
! that they are very similar to the values produced by WRF-physics over ocean.
! still exactly the same (checked); if called after seaice, exactly the same
! results as previously called inside gcm_diagnose_h1 (verified)
!---------------------------------------------------------------------------------

        scalar_zbot_at_pc_surface(1:ncell) = pstate_wrf%zzz(1:ncell,1,1)-scalar_geopotential_at_pc_face_level(nlevp,1:ncell)/gravity
        scalar_us_at_pc_surface(1:ncell)   = zero
        scalar_vs_at_pc_surface(1:ncell)   = zero

        call shr_flux_atmOcn(nmax    = ncell,  & ! in 
                    zbot    = scalar_zbot_at_pc_surface        , & ! in
                    ubot    = pstate_wrf%u_phy (1:ncell,1,1)   , & ! in
                    vbot    = pstate_wrf%v_phy (1:ncell,1,1)   , & ! in
                    thbot   = pstate_wrf%th_phy(1:ncell,1,1)   , & ! in
                    qbot    = pstate_wrf%moist (1:ncell,1,1,1) , & ! in
                    rbot    = pstate_wrf%rhom  (1:ncell,1,1)   , & ! in
                    tbot    = pstate_wrf%t_phy (1:ncell,1,1)   , & ! in
                    us      = scalar_us_at_pc_surface          , & ! in
                    vs      = scalar_vs_at_pc_surface          , & ! in
                    ts      = psurf_wrf%tsk (1:ncell,1)        , & ! in
                    mask    = int(pstate%ocnfrac_at_pc_surface%f(1:ncell),i4)  , & ! in
                    sen     = psurf_wrf%hfx_atmOcn%f(1:ncell)  , & ! out
                    lat     = psurf_wrf%lh_atmOcn %f(1:ncell)  , & ! out
                    lwup    = psurf_wrf%lwup_atmOcn%f(1:ncell) , & ! out
                    evap    = psurf_wrf%qfx_atmOcn %f(1:ncell) , & ! out
                    taux    = psurf_wrf%taux_atmOcn%f(1:ncell) , & ! out
                    tauy    = psurf_wrf%tauy_atmOcn%f(1:ncell) , & ! out
                    tref    = psurf_wrf%t2m_atmOcn %f(1:ncell) , & ! out
                    qref    = psurf_wrf%q2m_atmOcn %f(1:ncell) , & ! out
                    duu10n  = psurf_wrf%uu10m_atmOcn%f(1:ncell), & ! out
                    ustar_sv= psurf_wrf%ustar_atmOcn%f(1:ncell))

        psurf_wrf%hfx_atmOcn%f(1:ncell)  = -psurf_wrf%hfx_atmOcn%f(1:ncell)
        psurf_wrf%lh_atmOcn%f(1:ncell)   = -psurf_wrf%lh_atmOcn%f(1:ncell)
        psurf_wrf%qfx_atmOcn%f(1:ncell)  = -psurf_wrf%qfx_atmOcn%f(1:ncell)
        psurf_wrf%lwup_atmOcn%f(1:ncell) = -psurf_wrf%lwup_atmOcn%f(1:ncell)
        psurf_wrf%taux_atmOcn%f(1:ncell) = -psurf_wrf%taux_atmOcn%f(1:ncell)
        psurf_wrf%tauy_atmOcn%f(1:ncell) = -psurf_wrf%tauy_atmOcn%f(1:ncell)

        where(int(pstate%ocnfrac_at_pc_surface%f(1:ncell),i4)/=0) psurf_wrf%taux%f(1:ncell) = psurf_wrf%taux_atmOcn%f(1:ncell)
        where(int(pstate%ocnfrac_at_pc_surface%f(1:ncell),i4)/=0) psurf_wrf%tauy%f(1:ncell) = psurf_wrf%tauy_atmOcn%f(1:ncell)

!-----------------------------------------------------------------------------------------
! Update seaice flux, albedo and ts here, combine these ocean values with seaice values
! following rongxy's implementation for an old GRIST version at 2020-07/08
! this will overwrite taux over ocean diagnosed above, for pure-water points,
! still exactly the same (checked )
!
       if(.not.doAquaPlanet.and.test_real_case)then
          call grist_seaice_run(ncell, nlev, lats, coszrs, pstate%ocnfrac_at_pc_surface%f(1:ncell) , & ! in
                                pseaice  = staticData_sic_at_pc_surface%f(1:ncell), & ! in
                                snowh    = psurf_wrf%snow  (1:ncell,1), & ! in
                                plwds    = psurf_wrf%glw   (1:ncell,1), & ! in
                                pswds    = psurf_wrf%swdown(1:ncell,1), & ! in
                                pahfs    = psurf_wrf%hfx   (1:ncell,1), & ! inout
                                pahfl    = psurf_wrf%lh    (1:ncell,1), & ! inout
                                pevap    = psurf_wrf%qfx   (1:ncell,1), & ! inout
                                plwup    = lwup(1:ncell), & ! not used WRF_phys just dum
                                ptaux    = psurf_wrf%taux%f(1:ncell), & ! not used WRF_phys just diag 
                                ptauy    = psurf_wrf%tauy%f(1:ncell), & ! not used WRF_phys just diag
                                ptsi     = sicetemp_at_pc_surface%f(1:ncell),& ! out
                                pts      = psurf_wrf%tsk  (1:ncell,1), &       ! inout
                                palbdvis = psurf_wrf%asdir(1:ncell,1), &       ! inout
                                palbivis = psurf_wrf%asdif(1:ncell,1), &       ! inout
                                palbdnir = psurf_wrf%aldir(1:ncell,1), &       ! inout
                                palbinir = psurf_wrf%aldif(1:ncell,1), &       ! inout
                                zbot     = pstate_wrf%zzz   (1:ncell,1,1)  , & ! in
                                ubot     = pstate_wrf%u_phy (1:ncell,1,1)  , & ! in
                                vbot     = pstate_wrf%v_phy (1:ncell,1,1)  , & ! in
                                thbot    = pstate_wrf%th_phy(1:ncell,1,1)  , & ! in
                                tbot     = pstate_wrf%t_phy (1:ncell,1,1)  , & ! in
                                qbot     = pstate_wrf%moist (1:ncell,1,1,1), & ! in
                                rbot     = pstate_wrf%rhom  (1:ncell,1,1) )
! overwrite xice with sic
           psurf_wrf%xice(1:ncell,1) = staticData_sic_at_pc_surface%f(1:ncell)
       endif

#ifdef USE_NOAHMP
if(.not.doAquaPlanet)then

#ifdef LAM_DOMAIN

! Acid test (QJ)
#ifdef REG_QJLAM
        if (trim(lsm_scheme)=='noahmp') then
#else
! LAM in reanalysis-data may cause LSM problem, so only call after one model step
        if (istep.gt.1.and.trim(lsm_scheme)=='noahmp') then
#endif

#else
        if (trim(lsm_scheme)=='noahmp') then !---cheyz
#endif

            snowh_lsm(1:ncell) = psurf_wrf%snow(1:ncell,1)   ! mm
            xice_lsm(1:ncell)  = psurf_wrf%xice(1:ncell,1)
            emiss_lsm(1:ncell) = psurf_wrf%emiss(1:ncell,1)
            qsfc_lsm(1:ncell)  = psurf_wrf%qsfc(1:ncell,1)
            tsk(1:ncell,1)     = psurf_wrf%tsk (1:ncell,1)
            hfx(1:ncell,1)     = psurf_wrf%hfx (1:ncell,1)
! assign or not does not affect solutions
            !qfx(1:ncell,1)     = psurf_wrf%qfx (1:ncell,1)
!
! yiz: we've reset this driver with updated iostream
!
                call grist_noahmp_run ( ncell,istep,nlevp,levsoil,              &  !in
                                        pstate_wrf%u_phy (1:ncell,1:nlevp,1)  , &  !in
                                        pstate_wrf%v_phy (1:ncell,1:nlevp,1)  , &  !in
                                        pstate_wrf%t_phy (1:ncell,1:nlevp,1)  , &  !in
                                        pstate_wrf%moist (1:ncell,1:nlevp,1,1), &  !in
                                        pstate_wrf%p8w   (1:ncell,1:nlevp,1  ), &  !in
                                        pstate_wrf%dz8w  (1:ncell,1:nlevp,1)  , &  !in
                                        qsfc_lsm(1:ncell)                     , &  !inout
                                        snowh_lsm(1:ncell)                    , &  !inout
                                        xice_lsm(1:ncell)                     , &  !inout
                                        emiss_lsm(1:ncell)                    , &  !inout
                                        coszrs(1:ncell)                       , &  !in
                                        julday, & !in
                                        psurf_wrf%swdown (1:ncell,1)          , &  !in , SWDOWN  [W m-2]
                                        psurf_wrf%glw    (1:ncell,1)          , &  !in , longwave down at surface [W m-2]
                                        psurf_wrf%raint  (1:ncell,1)          , &  !in , total rain within one step [ assuem m]
                                        !--out later
                                        lwup(1:ncell)  , & ! out, not used by WRF_phys
                                        hfx (1:ncell,1), & ! out, hfx
                                        qfx (1:ncell,1), & ! out, qfx
                                        taux(1:ncell)  , & ! out, not used by WRF_phys
                                        tauy(1:ncell)  , & ! out, not used by WRF_phys
                                        asdir  (1:ncell, 1), & ! 0.2-0.7 micro-meter srfc alb: direct rad
                                        asdif  (1:ncell, 1), & ! 0.2-0.7 micro-meter srfc alb: diffuse rad
                                        aldir  (1:ncell, 1), & ! 0.7-5.0 micro-meter srfc alb: direct rad
                                        aldif  (1:ncell, 1), & ! 0.7-5.0 micro-meter srfc alb: diffuse rad   ! out
                                        tsk (1:ncell,1)  )     ! inout

! only modify land points
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%tsk (1:ncell,1)   = tsk(1:ncell,1)    ! No-land value changed by lsm
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%hfx (1:ncell,1)   = hfx(1:ncell,1)    ! ...
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%qfx (1:ncell,1)   = qfx(1:ncell,1)    ! ...
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%taux%f(1:ncell)   = taux(1:ncell) ! for diagnostics
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%tauy%f(1:ncell)   = tauy(1:ncell) ! for diagnostics

!-----------------------------------------------------------------------
! V381 sfclay needs qsfc from last-step lsm as input
! as surface is called before lsm,
! the lsm restart-var qsfc should overwrite psurf_wrf%qsfc after reading;
! This land-qsfc is bad for some vr-grid applications (g8x16), check why (lt. gt??)! 
! we still use those in sfclay for all ocean/land cells (further
! tests show similar solutions)
!-----------------------------------------------------------------------

#ifndef BEF_PATCH
           where(psurf_wrf%xland(1:ncell,1).le.1.5_r8) psurf_wrf%qsfc(1:ncell,1)  = qsfc_lsm(1:ncell)
#endif
! the folloing vars will produce itendical solution, even if no "where statement" was used: identical solutions
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%emiss(1:ncell,1) = emiss_lsm(1:ncell)
! noahmp will modify seaice part albedo, do not use that for now
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%asdir(1:ncell,1) = asdir(1:ncell,1)
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%asdif(1:ncell,1) = asdif(1:ncell,1)
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%aldir(1:ncell,1) = aldir(1:ncell,1)
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%aldif(1:ncell,1) = aldif(1:ncell,1)
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%snow(1:ncell,1)  = snowh_lsm(1:ncell) ! no-land value not changed by lsm
           where(psurf_wrf%xland(1:ncell,1).lt.1.5_r8) psurf_wrf%xice(1:ncell,1)  = xice_lsm(1:ncell)

        end if
end if
#endif

!-----------------------------------------------------------------------
! set albedo, not used in fact
!-----------------------------------------------------------------------

       psurf_wrf%albedo = psurf_wrf%asdir

       call grist_PhysW_pkg_run_ac(ncell,nlevp,nspecies,istep,dtime,rearth*mesh%mean_edt_dist,coszrs)

       mesh%nv = mesh%nv_compute

!----------------------------------------------------------------------
!  Transform to model-dynamics compatible tendencies
!  mind  the vertical sequence
!----------------------------------------------------------------------

!$omp parallel  private(iv,ilev,itracer)  
!$omp do schedule(dynamic,10) 
       do iv  =  1, ncell
          do ilev = 1, nlev
             do itracer  = 1, ntracer
                ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,ilev,iv) = ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(itracer,nlev+1-ilev,iv)
             end do
          end  do
       end  do
!$omp end do nowait
!$omp end parallel 

!-----------------------------------------------------------------------
! use f3 for f2 tend
! ptm has moisture and pt contribution, but when one tend is zero, donot have
!-----------------------------------------------------------------------
      call dtp_interface_convert_tend_phy2dyn(mesh, ntracer, nlev, ncell,.true. , &
                                  scalar_mpressure_at_pc_full_level      , &
                                  scalar_tracer_mxrt_at_pc_full_level    , &
                                  scalar_pt_at_pc_full_level             , &
                                  scalar_mif_at_pc_full_level            , &
                                  ptend_f3%tend_u_wind_at_pc_full_level%f, &
                                  ptend_f3%tend_v_wind_at_pc_full_level%f, &
                                  ptend_f3%tend_potential_temp_at_pc_full_level%f   , &
                                  ptend_f3%tend_tracer_mxrt_at_pc_full_level%f(1:ntracer,:,:)    , &
                                  ptend_f2%tend_normal_velocity_at_edge_full_level%f, & ! send back to model dynamics
                                  ptend_f2%tend_potential_temp_at_pc_full_level%f  )    ! send back to model dynamics

!-----------------------------------------------------------------------
! use f1 for rk tend of U and PT, F1 also contain moisture but should not be
! added as dycore_f1, only for tracer
!-----------------------------------------------------------------------
      call dtp_interface_convert_tend_phy2dyn(mesh, ntracer, nlev, ncell,.true., &
                                  scalar_mpressure_at_pc_full_level      , &
                                  scalar_tracer_mxrt_at_pc_full_level    , &
                                  scalar_pt_at_pc_full_level             , &
                                  scalar_mif_at_pc_full_level            , &
                                  ptend_f1%tend_u_wind_at_pc_full_level%f, &
                                  ptend_f1%tend_v_wind_at_pc_full_level%f, &
                                  ptend_f1%tend_potential_temp_at_pc_full_level%f   , &
                                  ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(1:ntracer,:,:)      , &
                                  ptend_rk%tend_normal_velocity_at_edge_full_level%f, & ! send back to model dynamics
                                  ptend_rk%tend_potential_temp_at_pc_full_level%f )     ! send back to model dynamics
      !if(any(ptend_rk%tend_normal_velocity_at_edge_full_level%f.ne.0)) print*,"plz double check interface2, rk wind"

! t1: let pt contribution all in rk, let qv contribution all in f2 (with (no) sub for
! ndc?)

      if(nh_dynamics.or.ptendSubDiabPhys)then ! f2's ptm tend also has moisture contribution
! add pt's f2 tendency to rk, and f2 will be substracted
         ptend_rk%tend_potential_temp_at_pc_full_level%f = ptend_rk%tend_potential_temp_at_pc_full_level%f+ptend_f2%tend_potential_temp_at_pc_full_level%f
      end if

    return
  end subroutine grist_full_physpkg2_PhysW_run

#endif

!************************************************************
!   PRIVATE
!************************************************************

!------------------------------------------------------------
!           Fill pstate
!------------------------------------------------------------
#if (defined AMIPW_PHYSICS )

  subroutine dtp_interface_fill_pstate(mesh,ntracer,nlev,nlevp,ncell, dtime, lats, &
                                       scalar_Uwind_at_pc_full_level, &
                                       scalar_Vwind_at_pc_full_level, &
                                       scalar_potential_temp_at_pc_full_level, &
                                       scalar_temp_at_pc_full_level , &
                                       scalar_qqq_at_pc_full_level  , & ! moist
                                       scalar_mpressure_at_pc_face_level, &
                                       scalar_mpressure_at_pc_full_level, &
                                       scalar_geopotential_at_pc_face_level, &
                                       scalar_geopotential_at_pc_full_level, &
                                       scalar_www_at_pc_face_level  ,&
                                       scalar_omega_at_pc_full_level,&
                                       scalar_delhp_at_pc_full_level,&
                                       tend_ptdyn_at_pc_full_level  ,&
                                       tend_qvdyn_at_pc_full_level)

  use omp_lib
  use grist_physics_data_structure,   only: pstate
  use grist_datam_static_data_module, only: staticData_sst_at_pc_surface
  use grist_slab_ocean_module,        only: grist_prt_oml_run
! io
    type(global_domain),  intent(in)  :: mesh
    integer(i4)        ,  intent(in)  :: ntracer
    integer(i4)        ,  intent(in)  :: nlev
    integer(i4)        ,  intent(in)  :: nlevp
    integer(i4)        ,  intent(in)  :: ncell
    real(r8)           ,  intent(in)  :: dtime, lats(ncell)
    real(r8)           ,  intent(in)  :: scalar_Uwind_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_Vwind_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_potential_temp_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_temp_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_qqq_at_pc_full_level(ntracer,nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_mpressure_at_pc_face_level(nlevp,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_mpressure_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_geopotential_at_pc_face_level(nlevp,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_geopotential_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_www_at_pc_face_level(nlevp,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_omega_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_delhp_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: tend_ptdyn_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: tend_qvdyn_at_pc_full_level(nlev,mesh%nv_full)
! local
    integer(i4)       :: ilev, icell

!-----------------------------------------------------------------------
! Note that all wrf state from down to top, while GRIST model from top to down
! dz8w is dz between our face-levels (WRF's full level), t8w and p8w
! are t and p at our face levels (WRF's full), and the top-level value is not given
! note that in WRF, the "top-face" (GRIST's imagniary full) has nothing
! also note that WRF's nlev is our nlevp actually, so I should use nLevel
! inside all WRF physics to distinguish this
!-----------------------------------------------------------------------

! 3D 
!$omp parallel  private(icell,ilev) 
!$omp do schedule(dynamic,10) 
      DO icell = 1, ncell
! m-level var
         do ilev = 1, nlev ! for full-level of grist
            pstate_wrf%u_phy (icell,nlev+1-ilev,1)  = scalar_Uwind_at_pc_full_level(ilev,icell)
            pstate_wrf%v_phy (icell,nlev+1-ilev,1)  = scalar_Vwind_at_pc_full_level(ilev,icell)
            pstate_wrf%t_phy (icell,nlev+1-ilev,1)  = scalar_temp_at_pc_full_level (ilev,icell)
            pstate_wrf%th_phy(icell,nlev+1-ilev,1)  = scalar_potential_temp_at_pc_full_level(ilev,icell)
            pstate_wrf%p_phy (icell,nlev+1-ilev,1)  = scalar_mpressure_at_pc_full_level (ilev,icell)
            pstate_wrf%zzz   (icell,nlev+1-ilev,1)  = scalar_geopotential_at_pc_full_level(ilev,icell)/gravity
            pstate_wrf%pi_phy(icell,nlev+1-ilev,1)  = (scalar_mpressure_at_pc_full_level(ilev,icell)/p00)**(rdry/cp)
            pstate_wrf%rhom  (icell,nlev+1-ilev,1)  = (scalar_mpressure_at_pc_face_level(ilev+1,icell)-scalar_mpressure_at_pc_face_level(ilev,icell))/&
                                                      (scalar_geopotential_at_pc_face_level(ilev,icell)-scalar_geopotential_at_pc_face_level(ilev+1,icell))
            pstate_wrf%rhod  (icell,nlev+1-ilev,1)  = scalar_delhp_at_pc_full_level(ilev,icell)/&
                                                     (scalar_geopotential_at_pc_face_level(ilev,icell)-scalar_geopotential_at_pc_face_level(ilev+1,icell))
            pstate_wrf%moist (icell,nlev+1-ilev,1,1:ntracer)  = scalar_qqq_at_pc_full_level(1:ntracer,ilev,icell)

            pstate_wrf%dz8w(  icell,nlev+1-ilev,1)  = (scalar_geopotential_at_pc_face_level(ilev,icell)-scalar_geopotential_at_pc_face_level(ilev+1,icell))/gravity
            pstate_wrf%omega( icell,nlev+1-ilev,1)  = scalar_omega_at_pc_full_level(ilev,icell) ! only used by CU for hdc
            ptend_wrf%rthdyten(icell,nlev+1-ilev,1) = tend_ptdyn_at_pc_full_level(ilev,icell)   ! only used by CU
            ptend_wrf%rqvdyten(icell,nlev+1-ilev,1) = tend_qvdyn_at_pc_full_level(ilev,icell)   ! only used by CU
         end do
         pstate_wrf%u_phy (icell,nlevp,1)  = fillvalue
         pstate_wrf%v_phy (icell,nlevp,1)  = fillvalue
         pstate_wrf%t_phy (icell,nlevp,1)  = fillvalue
         pstate_wrf%th_phy(icell,nlevp,1)  = fillvalue
         pstate_wrf%p_phy (icell,nlevp,1)  = fillvalue
         pstate_wrf%zzz   (icell,nlevp,1)  = fillvalue
         pstate_wrf%pi_phy(icell,nlevp,1)  = fillvalue
         pstate_wrf%rhom  (icell,nlevp,1)  = fillvalue
         pstate_wrf%rhod  (icell,nlevp,1)  = fillvalue
         pstate_wrf%moist (icell,nlevp,1,:)= fillvalue
         pstate_wrf%dz8w  (icell,nlevp,1)  = fillvalue
         pstate_wrf%omega (icell,nlevp,1)  = fillvalue
         pstate_wrf%w0avg (icell,nlevp,1)  = fillvalue
! face-level www and already time-averaged, only nh_dynamics has meaning for this; do a reverse
         pstate_wrf%www(icell,1:nlevp,1)   = scalar_www_at_pc_face_level(nlevp:1:-1,icell)
! Used for CAMUW PBL scheme, LiXH add        
         pstate_wrf%z8w(icell,1:nlevp,1)   = scalar_geopotential_at_pc_face_level(nlevp:1:-1,icell)/gravity
         if(nh_dynamics)then
! set w0avg as a full level var, as input in wrf-physics
           pstate_wrf%w0avg(icell,1:nlev,1)= 0.5_r8*(pstate_wrf%www(icell,1:nlev,1)+pstate_wrf%www(icell,2:nlevp,1))
         else
! for hydrostatic, www convered from omega based on their definition
           pstate_wrf%w0avg(icell,1:nlev,1)= -pstate_wrf%omega(icell,1:nlev,1)/gravity/pstate_wrf%rhom(icell,1:nlev,1)
         end if

         pstate_wrf%p8w(icell,1:nlevp,1)   = scalar_mpressure_at_pc_face_level(nlevp:1:-1,icell)
! interpolate w-level t based on m-level t
         do ilev = 2, nlev ! ilev level of WRF > nlev+2-ilev face of GRIST, the surrounding full levels are nlev+2-ilev, nlev+1-ilev
             pstate_wrf%t8w(icell,ilev,1)= 0.5*(deta_full(nlev+1-ilev)/deta_face(nlev+2-ilev)*scalar_temp_at_pc_full_level(nlev+2-ilev,icell)+&
                                                deta_full(nlev+2-ilev)/deta_face(nlev+2-ilev)*scalar_temp_at_pc_full_level(nlev+1-ilev,icell))
         end do
! we can use some extrapolate here but simply use full-level value now
         pstate_wrf%t8w(icell,nlevp,1)   = scalar_temp_at_pc_full_level(1,icell)
      END DO
!$omp end do nowait
!$omp end parallel 

!-----------------------------------------------------------------------
! once called, slab ocean will modify staticData_sst_at_pc_surface 
! directly at each fast-physics step
!-----------------------------------------------------------------------

      pstate%sst_at_pc_surface%f   = staticData_sst_at_pc_surface%f ! sst will change each month

      if(physics_use_som)then
! every input must be initialized, or will be bugged
         call grist_prt_oml_run(ncell,ncell,nlev,dtime,lats(1:ncell)   , &
                             pstate%sst_at_pc_surface%f(1:ncell),& ! inout
                             psurf_wrf%ust(1:ncell,1)        , &
                             pstate_wrf%u_phy (1:ncell,1,1)  , &
                             pstate_wrf%v_phy (1:ncell,1,1)  , &
                             psurf_wrf%xland  (1:ncell,1)    , &
                             psurf_wrf%hfx    (1:ncell,1)    , &
                             psurf_wrf%lh     (1:ncell,1)    , &
                             psurf_wrf%gsw    (1:ncell,1)    , &
                             psurf_wrf%glw    (1:ncell,1)    , &
                             psurf_wrf%emiss  (1:ncell,1))
      end if
! this line should in front of t8w, check after regression
! where(psurf_wrf%xland(1:ncell,1).ge.1.5_r8) psurf_wrf%tsk (1:ncell,1) = staticData_sst_at_pc_surface%f(1:ncell)

!-----------------------------------------------------------------------
! Modify ts over sea
! Do not consider to include seaice for ts now
!-----------------------------------------------------------------------

      do icell = 1, ncell
         if(psurf_wrf%xland(icell,1).ge.1.5_r8)then
            psurf_wrf%tsk(icell,1) = pstate%sst_at_pc_surface%f(icell)
#ifdef DTP_TSICE
! For seaice conc >0.01 (threshold), combine ts with sst and sicetemp following Rongxy's suggestion
            if(staticData_sic_at_pc_surface%f(icell) .ge. 0.01_r8) then
                psurf_wrf%tsk(icell,1) = staticData_sst_at_pc_surface%f(icell)*(one-staticData_sic_at_pc_surface%f(icell))+&
                                               sicetemp_at_pc_surface%f(icell)*     staticData_sic_at_pc_surface%f(icell)
            endif
#endif
         end if
      end do

!-----------------------------------------------------------------------
! Surface
!-----------------------------------------------------------------------

      psurf_wrf%ht  (1:ncell,1)   = scalar_geopotential_at_pc_face_level(nlevp,1:ncell)/gravity ! used by lin microphysics
      psurf_wrf%psfc(1:ncell,1)   = scalar_mpressure_at_pc_face_level(nlevp,1:ncell)
      pstate_wrf%t8w(1:ncell,1,1) = psurf_wrf%tsk(1:ncell,1)

      return
  end subroutine dtp_interface_fill_pstate

#endif

!------------------------------------------------------------
! Transform physics tendencies to model-dynamics compatible
! style, this is for general purposese, should be regressed
! to the before state given the same config (2019-Nov)
!------------------------------------------------------------

  subroutine dtp_interface_convert_tend_phy2dyn(mesh, ntracer, nlev, ncell , flag,  &
                                  scalar_mpressure_at_pc_full_level,   &
                                  scalar_tracer_mxrt_at_pc_full_level, & ! dry mxrt
                                  scalar_pt_at_pc_full_level  ,&
                                  scalar_mif_at_pc_full_level, &
                                  tend_duDt_at_pc_full_level,  &
                                  tend_dvDt_at_pc_full_level,  &
                                  tend_dptDt_at_pc_full_level, &
                                  tend_dqDt_at_pc_full_level,  &
                                  tend_dunDt_at_edge_full_level, &
                                  tend_dptmDt_at_pc_full_level )
! io
    use omp_lib
    type(global_domain)  , intent(in)     :: mesh
    integer(i4),           intent(in)     :: ntracer
    integer(i4),           intent(in)     :: nlev
    integer(i4),           intent(in)     :: ncell
    logical    ,           intent(in)     :: flag
    real(r8)   ,           intent(in)     :: scalar_mpressure_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: scalar_tracer_mxrt_at_pc_full_level(ntracer,nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: scalar_pt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: scalar_mif_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_duDt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_dvDt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_dptDt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_dqDt_at_pc_full_level(ntracer,nlev,mesh%nv_full)
    real(r8)   ,           intent(inout)  :: tend_dunDt_at_edge_full_level(nlev,mesh%ne_full)
    real(r8)   ,           intent(inout)  :: tend_dptmDt_at_pc_full_level(nlev,mesh%nv_full)
! local
    integer(i4)   :: icell1, icell2, ilev, iv, ie
    real(r8)      :: point_cart_coort(3), vector_velocity(3)
    real(r8)      :: tend_U_wind_at_edge_full_level 
    real(r8)      :: tend_V_wind_at_edge_full_level

!-----------------------------------------------------------------------
! Transform u,v tendency from cell to normal 
! wind tendency at edge
!-----------------------------------------------------------------------

!$omp parallel private(ie,icell1,icell2,point_cart_coort,ilev,tend_U_wind_at_edge_full_level,tend_V_wind_at_edge_full_level,vector_velocity)  
!$omp do schedule(dynamic,10) 
        do ie = 1, mesh%ne
           icell1 = mesh%edt_v(1,ie)
           icell2 = mesh%edt_v(2,ie)
           point_cart_coort(:) = mesh%edt_c_p(1:3,ie)
           do ilev = 1, nlev
              tend_U_wind_at_edge_full_level = half*(tend_duDt_at_pc_full_level(ilev,icell1)+tend_duDt_at_pc_full_level(ilev,icell2))
              tend_V_wind_at_edge_full_level = half*(tend_dvDt_at_pc_full_level(ilev,icell1)+tend_dvDt_at_pc_full_level(ilev,icell2))
              call convert_vector_sph2cart(tend_U_wind_at_edge_full_level,&
                                           tend_V_wind_at_edge_full_level,&
                                           point_cart_coort, vector_velocity)
              tend_dunDt_at_edge_full_level(nlev+1-ilev,ie) = dot_product(vector_velocity, mesh%edp_nr(1:3,ie))
           end do
        end do
!$omp end do nowait
!$omp end parallel 

!-----------------------------------------------------------------------
! find mxrt and thetam tendency
! partial, reverse vertical sequence
!-----------------------------------------------------------------------

      if(flag)then ! default
        do iv = 1, ncell
           do ilev = 1, nlev
                 tend_dptmDt_at_pc_full_level(ilev,iv) = (one+rvap/rdry*scalar_tracer_mxrt_at_pc_full_level(1,ilev,iv))*&
                                                          tend_dptDt_at_pc_full_level(nlev+1-ilev,iv)+&
                                                         (rvap/rdry*scalar_pt_at_pc_full_level(ilev,iv)*&
                                                          tend_dqDt_at_pc_full_level(1,nlev+1-ilev,iv))
           end do
        end do
      else ! ignore dqDt in dptm/dt
        do iv = 1, ncell
           do ilev = 1, nlev
                 tend_dptmDt_at_pc_full_level(ilev,iv) = (one+rvap/rdry*scalar_tracer_mxrt_at_pc_full_level(1,ilev,iv))*&
                                                          tend_dptDt_at_pc_full_level(nlev+1-ilev,iv)
           end do
        end do
      end if

     return
  end subroutine dtp_interface_convert_tend_phy2dyn

 end module grist_dtp_interface_full_physpkg2
