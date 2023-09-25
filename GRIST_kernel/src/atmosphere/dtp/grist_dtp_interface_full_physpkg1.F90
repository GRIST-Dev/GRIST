
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: This is the calling interface to UVTQ-style full physics, as 
!              a template for any of such type (e.g., CAM, ECHAM). The
!              full physpkg depends on the pstate/ptend_f3 structure while 
!              simple physics does not. But when we call simple-physics inside 
!              a full physics package, it must (tested) produce the same 
!              results as we directly call the simple-physics interface given 
!              the same configuration. This ensures the MPI-behavior of full 
!              physpkg is as same as simple-physics, so seq-grist does not need
!              full physpkg currently, and we have no need to refactor the simple 
!              physics interface.
!
! Revision history:
!            1. CAM5 package is added
!            1. Now is PhysC for AMIPC
!----------------------------------------------------------------------------

 module grist_dtp_interface_full_physpkg1

   use grist_mpi
   use grist_constants,                  only: r8, i4, pi, rdry, rvap, cp, p00, omega, gravity, rearth, zero, one, half, ptfactor
   use grist_nml_module,                 only: physpkg, doAquaPlanet, test_real_case, start_ymd, start_tod
   use grist_domain_types,               only: global_domain
   use grist_dycore_gcd_recon_module_2d, only: vector_recon_perot_edge2cell_uv_2d
   use grist_math_module,                only: convert_vector_sph2cart
   use grist_dtp_interface_module,       only: scalar_U_wind_at_pc_full_level, &
                                               scalar_V_wind_at_pc_full_level, &
                                               scalar_qqq_at_pc_full_level
   use grist_time_manager,               only: get_curr_calday
   use grist_zenith,                     only: zenith
#ifdef AMIPC_PHYSICS
   use grist_physpkg_PhysC,          only: grist_physpkgbc_PhysC_run,&
                                               grist_physpkgac_PhysC_run
   use phys_control,                     only: lsm_scheme  !--cheyz  20200425
#endif
    use grist_physics_data_structure,    only: pstate, ptend_f3
    use grist_datam_static_data_module,  only: staticData_sst_at_pc_surface, &
                                               staticData_sic_at_pc_surface, &
                                               staticData_var2d_at_pc_surface
    use grist_datam_initial_data_module, only: initialData_skintemp_at_pc_surface
    use grist_prescribe_seaice_module,   only: sicetemp_at_pc_surface, grist_seaice_init, grist_seaice_run


   implicit none
   private
   
   public :: grist_full_physpkg1_init, &
             grist_full_physpkg1_final

#ifdef AMIPC_PHYSICS
   public :: grist_full_physpkg1_PhysC_run
#endif

  contains
!---------------------------------------------------------
! initialization/finalization associated with physpkg
!---------------------------------------------------------

  subroutine  grist_full_physpkg1_init(local_block, model_timestep)
   use grist_domain_types,               only: block_structure

#ifdef AMIPC_PHYSICS
   use grist_PhysC_data_structure,       only: grist_PhysC_data_structure_construct
   use grist_physpkg_PhysC,          only: grist_physpkg_PhysC_init
   use grist_lnd_static_vars_module,     only: grist_lnd_static_vars_construct, &
                                               scalar_static_landmask_at_pc
   use grist_PhysC_data_structure,       only: pstate_cam
#endif
    type(block_structure) ,target, intent(inout)  :: local_block
    real(r8),                      intent(in)     :: model_timestep

    integer(i4) :: ncell
    real(r8), allocatable    :: dxmean(:)
#ifndef NO_PHYSHALO
    ncell = local_block%full_domain%nv_halo(1)
#else
    ncell = local_block%full_domain%nv_compute
#endif
!
! model physics package init
! normally, a column model only needs the halo(1) circle
!
#ifdef AMIPC_PHYSICS
     if(trim(physpkg).eq.'AMIPC_PHYSICS')then

#ifndef NO_PHYSHALO
        local_block%full_domain%nv = local_block%full_domain%nv_halo(1)
#else
        local_block%full_domain%nv = local_block%full_domain%nv_compute
#endif
        call grist_PhysC_data_structure_construct(local_block%full_domain)

! LiXH add mesh dx for scale_factor of convection.
        allocate(dxmean(ncell));dxmean=0._r8
        dxmean(1:ncell)  = local_block%full_domain%vtxCellLeng(1:ncell)*rearth

#ifndef NO_PHYSHALO
        call grist_physpkg_PhysC_init(local_block%full_domain%nv_halo(1), model_timestep, 0, dxmean,    &
                                          local_block%full_domain%vtx_lat(1:ncell), local_block%full_domain%vtx_lon(1:ncell))
#else
        call grist_physpkg_PhysC_init(local_block%full_domain%nv_compute, model_timestep, 0, dxmean,    &
                                          local_block%full_domain%vtx_lat(1:ncell), local_block%full_domain%vtx_lon(1:ncell))
#endif

!        local_block%full_domain%nv = local_block%full_domain%nv_halo(1)
        call grist_lnd_static_vars_construct(local_block%full_domain)
        local_block%full_domain%nv = local_block%full_domain%nv_compute
!
! init land surface
!

        pstate%landfrac_at_pc_surface%f(1:ncell) = scalar_static_landmask_at_pc%f(1:ncell)

        where(scalar_static_landmask_at_pc%f(1:ncell).le.0.5_r8) pstate%ocnfrac_at_pc_surface%f(1:ncell)  = 1._r8
        where(scalar_static_landmask_at_pc%f(1:ncell).gt.0.5_r8) pstate%ocnfrac_at_pc_surface%f(1:ncell)  = 0._r8

!
! LiXH add sgh30 (var2d) for PBL scheme.
! Note that CAM uses sgh30 (variation from 30s to 10min), var2d maybe not fitted !! 
!
        pstate_cam%sgh30_at_pc_surface%f(1:ncell) = staticData_var2d_at_pc_surface%f(1:ncell)


!
! We have no seaice model, set pstate%icefrac_at_pc_surface%f(i)  = 0._r8
! It should be modified after the seaice model being coupled in.
!        pstate%icefrac_at_pc_surface%f  = 0._r8

!
! init all ts and sst as sst data(real sst data+ analytic sst function)
!
#ifdef USE_NOAHMP
        if(test_real_case) pstate%ts_at_pc_surface%f(1:ncell)       = initialData_skintemp_at_pc_surface%f(1:ncell)
#endif
!
! init all ts and sst as sst data(real sst data+ analytic sst function)
!
        where(pstate%landfrac_at_pc_surface%f(1:ncell).le.0.5_r8) pstate%ts_at_pc_surface%f(1:ncell) = staticData_sst_at_pc_surface%f(1:ncell)

        pstate%sst_at_pc_surface%f(1:ncell)  = staticData_sst_at_pc_surface%f(1:ncell)

        if(mpi_rank()==0) print*, "Sucessfully initiate AMIPC physics (PhysC) in GRIST"

        deallocate(dxmean)
     end if
#endif

   return
  end subroutine  grist_full_physpkg1_init

  subroutine  grist_full_physpkg1_final
! Add physics package compiling option here

#ifdef AMIPC_PHYSICS
   use grist_PhysC_data_structure,     only: grist_PhysC_data_structure_destruct
   use grist_physpkg_PhysC,        only: grist_physpkg_PhysC_final
   use grist_lnd_static_vars_module,   only: grist_lnd_static_vars_destruct
#endif

#ifdef AMIPC_PHYSICS
      call grist_PhysC_data_structure_destruct
      call grist_physpkg_PhysC_final
      call grist_lnd_static_vars_destruct
#endif
    return
  end subroutine  grist_full_physpkg1_final

#ifdef AMIPC_PHYSICS

  subroutine grist_full_physpkg1_PhysC_run(mesh,ntracer,nmif,nlev,ncell,istep,dtime , &
                                                   scalar_mif_at_pc_full_level              , &
                                                   scalar_geopotential_at_pc_full_level     , &
                                                   scalar_geopotential_at_pc_face_level     , &
                                                   scalar_mpressure_at_pc_full_level        , &
                                                   scalar_mpressure_at_pc_face_level        , &
                                                   scalar_normal_velocity_at_edge_full_level, &
                                                   scalar_potential_temp_at_pc_full_level   , & !thetam
                                                   scalar_temp_at_pc_full_level             , &
                                                   scalar_tracer_mxrt_at_pc_full_level      , &
                                                   scalar_omega_at_pc_full_level            , &
                                                   tend_normal_velocity_at_edge_full_level  , &
                                                   tend_potential_temp_at_pc_full_level     , &
                                                   tend_tracer_mxrt_at_pc_full_level        )

    use grist_amp_surface_flux_module , only: grist_amp_surface_flux_run
    use grist_radiation,                only: call_radiation
    use grist_PhysC_data_structure,     only: pstate_cam
    use grist_data_types,               only: exchange_field_list_2d
    use grist_config_partition,         only: exchange_data_2d_add, &
                                              exchange_data_2d
    use grist_lsm_noahmp_driver,        only: grist_noahmp_run   !--cheyz
    use grist_nml_module,               only: levsoil            !--cheyz
! io
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: ntracer
    integer(i4)          , intent(in)   :: nmif
    integer(i4)          , intent(in)   :: nlev
    integer(i4)          , intent(in)   :: ncell
    integer(i4)          , intent(in)   :: istep
    real(r8)             , intent(in)   :: dtime
    real(r8),allocatable , intent(in)   :: scalar_mif_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(in)   :: scalar_omega_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_tracer_mxrt_at_pc_full_level(:,:,:)
! local
    real(r8)                            :: lats(ncell),lons(ncell)
    integer(i4)                         :: itracer, ie, iv, ilev, icell1, icell2, nlevp
    real(r8)                            :: dzbot_1(ncell),dzbot_2(ncell), rhobot(ncell)
    type(exchange_field_list_2d),pointer:: field_head_2d
    real(r8)                            :: julday, coszrs(ncell)
    real(r8)                            :: lhbot(ncell), thbot(ncell), zbot(ncell)

    real(r8)                            :: snowh_lsm(ncell), xice_lsm(ncell),                   &
                                           lwup_lsm(ncell),  shflx_lsm(ncell), qflx_lsm(ncell), &
                                           taux_lsm(ncell),  tauy_lsm(ncell),                   &
                                           asdir_lsm(ncell), asdif_lsm(ncell),                  &
                                           aldir_lsm(ncell), aldif_lsm(ncell), ts_lsm(ncell) 
 

        field_head_2d =>null()
        nlevp = nlev+1
!
! prepare data
!
#ifndef NO_PHYSHALO
        if(ncell.ne.mesh%nv_halo(1)) then
#else
        if(ncell.ne.mesh%nv_compute) then
#endif
           print*, "dim stuff wrong in grist_full_physpkg1_PhysCPhysics_run, model aborts"
           call mpi_abort()
        end if

        do iv = 1, ncell
           lats(iv) = mesh%vtx_lat(iv) ! radian
           lons(iv) = mesh%vtx_lon(iv)
        end do

        call get_curr_calday(start_ymd, start_tod, istep, dtime, julday)
        call zenith(julday, coszrs, ncell, lons, lats)

!----------------------------------------------------------------------
! obtain U&V
!----------------------------------------------------------------------

#ifndef NO_PHYSHALO
        mesh%nv = mesh%nv_halo(1)
#else
        mesh%nv = mesh%nv_compute
#endif
        call vector_recon_perot_edge2cell_uv_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                                      scalar_U_wind_at_pc_full_level, &
                                                      scalar_V_wind_at_pc_full_level, &
                                                      nlev)
        scalar_qqq_at_pc_full_level = scalar_tracer_mxrt_at_pc_full_level

!----------------------------------------------------------------------
! for mxrt-like quantity, from dry to moist
!----------------------------------------------------------------------
! for cam physics, it is just coincidence that nmif also represent number of mxrt-like var from 1st
! so use it here
        do itracer  = 1, nmif
           scalar_qqq_at_pc_full_level(itracer,1:nlev,1:ncell) = scalar_tracer_mxrt_at_pc_full_level(itracer,1:nlev,1:ncell)*&
                                                                 scalar_mif_at_pc_full_level(1:nlev,1:ncell)
        end do
! from dry to moist; for regressing tc_a (an special case)
        !scalar_qqq_at_pc_full_level = scalar_tracer_mxrt_at_pc_full_level(:,:,:)/(one+scalar_tracer_mxrt_at_pc_full_level(:,:,:))

!----------------------------------------------------------------------
! fill pstate with dynamics state
!----------------------------------------------------------------------

        call dtp_interface_fill_pstate(mesh,ntracer,nlev,nlevp,ncell,  &
                                       scalar_U_wind_at_pc_full_level, &
                                       scalar_V_wind_at_pc_full_level, &
                                       scalar_temp_at_pc_full_level  , &
                                       scalar_qqq_at_pc_full_level   , & ! moist
                                       scalar_omega_at_pc_full_level , &
                                       scalar_mpressure_at_pc_face_level, &
                                       scalar_mpressure_at_pc_full_level, &
                                       scalar_geopotential_at_pc_face_level, &
                                       scalar_geopotential_at_pc_full_level)

!----------------------------------------------------------------------
! call ocean surface scheme, which relies on pstate of model physics
!----------------------------------------------------------------------

!!        call grist_amp_surface_flux_run(ncell,nlev)

!----------------------------------------------------------------------
! call interface for CAM5 tendency, relies on pstate and produce ptend_f3
! same as call bc, surface, ac sequentially
!----------------------------------------------------------------------

!        call grist_physpkg_PhysC_run(ncell, lats, lons, istep, dtime, dtime)

!----------------------------------------------------------------------
!  seperate physac, physbc, surface here, not in physpkg
!  for testing purpose
!----------------------------------------------------------------------

        call grist_physpkgbc_PhysC_run(ncell, lats, lons, istep, dtime, dtime)


!*******************THIS part is for testing Real-World Climate Modeling****************
!
! ocean surface and ac physics, ocean flux will overwrite if wrfldn module also
! modify ocean points; tested ok
!
        if(doAquaPlanet.or.trim(lsm_scheme).eq.'noahmp')then
                call grist_amp_surface_flux_run(ncell,nlev,istep,dtime,lats,lons,call_radiation)
        endif
 
! update seaice flux, albedo and ts here, combine these ocean values with seaice values
! following rongxy's implementation for an old GRIST version at 2020-07/08
!
       if(.not.doAquaPlanet.and.test_real_case)then
           rhobot(1:ncell)   = pstate%delp_at_pc_full_level%f(nlev,1:ncell)/(pstate%geop_at_pc_face_level%f(nlev,1:ncell)-pstate%geop_at_pc_face_level%f(nlev+1,1:ncell))

           do iv = 1, ncell
               thbot(iv) = scalar_potential_temp_at_pc_full_level(nlev,iv)/(one+ptfactor*scalar_qqq_at_pc_full_level(1,nlev,iv))

               !Note: zbot includes surface height, different with that used in CAMphys, LiXH
               zbot(iv)  = pstate%z_at_pc_full_level%f(nlev,iv)+pstate%geop_at_pc_surface%f(iv)/gravity
           end do

           call grist_seaice_run(ncell, nlev, lats, coszrs, pstate%ocnfrac_at_pc_surface%f(1:ncell) , & ! in
                                pseaice  = staticData_sic_at_pc_surface%f(1:ncell),                   & ! in
                                snowh    = pstate%snowhland_at_pc_surface%f(1:ncell),                 & ! in 
                                plwds    = pstate_cam%flwds_at_pc_surface%f(1:ncell),                 & ! in
                                pswds    = pstate_cam%fswds_at_pc_surface%f(1:ncell),                 & ! in
                                pahfs    = pstate%atm_in_shflx_at_pc_surface%f(1:ncell),              & ! inout
                                pahfl    = lhbot(1:ncell),                                            & ! inout, not used CAM_phys just dum
                                pevap    = pstate%atm_in_qflx_at_pc_surface%f(1,1:ncell),             & ! inout
                                plwup    = pstate%atm_in_lwup_at_pc_surface%f(  1:ncell),             & ! inout
                                ptaux    = pstate%atm_in_taux_at_pc_surface%f(  1:ncell),             & ! inout
                                ptauy    = pstate%atm_in_tauy_at_pc_surface%f(  1:ncell),             & ! inout
                                ptsi     = sicetemp_at_pc_surface%f(1:ncell),                         & ! out
                                pts      = pstate%ts_at_pc_surface%f(1:ncell),                        & ! inout
                                palbdvis = pstate%atm_in_asdir_at_pc_surface%f( 1:ncell),             & ! inout
                                palbivis = pstate%atm_in_asdif_at_pc_surface%f( 1:ncell),             & ! inout
                                palbdnir = pstate%atm_in_aldir_at_pc_surface%f( 1:ncell),             & ! inout
                                palbinir = pstate%atm_in_aldif_at_pc_surface%f( 1:ncell),             & ! inout
                                zbot     = zbot(1:ncell) ,                                            & ! in
                                ubot     = pstate%u_wind_at_pc_full_level%f(nlev,1:ncell),            & ! in
                                vbot     = pstate%v_wind_at_pc_full_level%f(nlev,1:ncell),            & ! in
                                thbot    = thbot(1:ncell),                                            & ! in
                                tbot     = pstate%temp_at_pc_full_level%f(nlev,1:ncell),              & ! in
                                qbot     = pstate%tracer_mxrt_at_pc_full_level%f(1,nlev,1:ncell),     & ! in
                                rbot     = rhobot  (1:ncell)   )                                        ! in
! overwrite icefrac with sic
! Note: icefrac_at_pc_surface only contains sea-ice, LiXH.
           pstate%icefrac_at_pc_surface%f(1:ncell) = staticData_sic_at_pc_surface%f(1:ncell)
       endif

!
! land surface
!

if(.not.doAquaPlanet)then
        if(lsm_scheme=='wrflnd')then   !---cheyz

        else if (trim(lsm_scheme)=='noahmp') then !---cheyz
                snowh_lsm = pstate%snowhland_at_pc_surface%f(1:ncell)
                xice_lsm  = pstate%icefrac_at_pc_surface%f(1:ncell)
                ts_lsm    = pstate%ts_at_pc_surface%f(1:ncell) 
                
                call grist_noahmp_run ( ncell,istep,nlev,levsoil,                          &  !in  
                                        pstate%u_wind_at_pc_full_level%f(1:nlev,1:ncell) , &  !in  
                                        pstate%v_wind_at_pc_full_level%f(1:nlev,1:ncell) , &  !in  
                                        pstate%temp_at_pc_full_level%f(  1:nlev,1:ncell) , &  !in  
                                        pstate%tracer_mxrt_at_pc_full_level%f(1,1:nlev,1:ncell), &   !in
                                        pstate%pressure_at_pc_face_level%f (    1:nlevp,1:ncell),&   !in 
                                        pstate_cam%fswds_at_pc_surface%f(1:ncell)    , &  !in ! SWDOWN  [W m-2]
                                        pstate_cam%flwds_at_pc_surface%f(1:ncell)    , &  !in ! longwave down at surface [W m-2]
                                        pstate%scalar_prect_surface%f       (1:ncell), &  !in
                                        !--out later
                                        snowh_lsm, xice_lsm,                                     &
                                        lwup_lsm,  shflx_lsm, qflx_lsm,  taux_lsm,  tauy_lsm,    &
                                        asdir_lsm, asdif_lsm, aldir_lsm, aldif_lsm, ts_lsm   )
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%snowhland_at_pc_surface%f(1:ncell)    = snowh_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%icefrac_at_pc_surface%f(1:ncell)      = xice_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_lwup_at_pc_surface%f(1:ncell)  = lwup_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_shflx_at_pc_surface%f(1:ncell) = shflx_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_qflx_at_pc_surface%f(1,1:ncell)= qflx_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_taux_at_pc_surface%f(1:ncell)  = taux_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_tauy_at_pc_surface%f(1:ncell)  = tauy_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_asdir_at_pc_surface%f(1:ncell) = asdir_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_asdif_at_pc_surface%f(1:ncell) = asdif_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_aldir_at_pc_surface%f(1:ncell) = aldir_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%atm_in_aldif_at_pc_surface%f(1:ncell) = aldif_lsm(1:ncell)
                where(pstate%landfrac_at_pc_surface%f(1:ncell) .gt. 0.5_r8) pstate%ts_at_pc_surface%f(1:ncell)           = ts_lsm(1:ncell)
        end if
end if
        

        call grist_physpkgac_PhysC_run(ncell, lats, lons, istep, dtime, dtime)

#ifdef NO_PHYSHALO
        call exchange_data_2d_add(mesh,field_head_2d,ptend_f3%tend_u_wind_at_pc_full_level)
        call exchange_data_2d_add(mesh,field_head_2d,ptend_f3%tend_v_wind_at_pc_full_level)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

!*******************THIS part is for testing Real-World Climate Modeling****************

        mesh%nv = mesh%nv_compute

!----------------------------------------------------------------------
!  Transform to model-dynamics compatible tendencies
!----------------------------------------------------------------------

        call dtp_interface_convert_tend_phy2dyn(mesh, ntracer, nlev, ncell ,         &
                                  scalar_mpressure_at_pc_full_level,   &
                                  scalar_tracer_mxrt_at_pc_full_level, &
                                  scalar_qqq_at_pc_full_level        , &       
                                  scalar_temp_at_pc_full_level       , &
                                  scalar_mif_at_pc_full_level        , &
                                  ptend_f3%tend_u_wind_at_pc_full_level%f, &
                                  ptend_f3%tend_v_wind_at_pc_full_level%f, &
                                  ptend_f3%tend_temp_at_pc_full_level%f,   &
                                  ptend_f3%tend_tracer_mxrt_at_pc_full_level%f, &
                                  tend_normal_velocity_at_edge_full_level, &
                                  tend_potential_temp_at_pc_full_level   , &
                                  tend_tracer_mxrt_at_pc_full_level      )
!
! reset tend to zero for use of Next step
!
       ptend_f3%tend_u_wind_at_pc_full_level%f       = zero
       ptend_f3%tend_v_wind_at_pc_full_level%f       = zero
       ptend_f3%tend_temp_at_pc_full_level%f         = zero
       ptend_f3%tend_tracer_mxrt_at_pc_full_level%f  = zero


      return
  end subroutine grist_full_physpkg1_PhysC_run

#endif

!************************************************************
!   PRIVATE
!************************************************************

!------------------------------------------------------------
!           Fill pstate
!------------------------------------------------------------
#if (defined  AMIPC_PHYSICS)

  subroutine dtp_interface_fill_pstate(mesh,ntracer,nlev,nlevp,ncell,  &
                                       scalar_Uwind_at_pc_full_level, &
                                       scalar_Vwind_at_pc_full_level, &
                                       scalar_temp_at_pc_full_level  , &
                                       scalar_qqq_at_pc_full_level   , & ! moist
                                       scalar_omega_at_pc_full_level , &
                                       scalar_mpressure_at_pc_face_level, &
                                       scalar_mpressure_at_pc_full_level, &
                                       scalar_geopotential_at_pc_face_level, &
                                       scalar_geopotential_at_pc_full_level)

  use grist_physics_data_structure,   only: pstate
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,     only: pstate_cam
#endif
! io
    type(global_domain),  intent(in)  :: mesh
    integer(i4)        ,  intent(in)  :: ntracer
    integer(i4)        ,  intent(in)  :: nlev
    integer(i4)        ,  intent(in)  :: nlevp
    integer(i4)        ,  intent(in)  :: ncell
    real(r8)           ,  intent(in)  :: scalar_Uwind_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_Vwind_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_temp_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_omega_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_qqq_at_pc_full_level(ntracer,nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_mpressure_at_pc_face_level(nlevp,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_mpressure_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_geopotential_at_pc_face_level(nlevp,mesh%nv_full)
    real(r8)           ,  intent(in)  :: scalar_geopotential_at_pc_full_level(nlev,mesh%nv_full)
! local
    integer(i4)       :: ilev, icell

!
! u,v,t,q
!
      pstate%u_wind_at_pc_full_level%f(1:nlev,1:ncell)                =        scalar_Uwind_at_pc_full_level(1:nlev,1:ncell)
      pstate%v_wind_at_pc_full_level%f(1:nlev,1:ncell)                =        scalar_Vwind_at_pc_full_level(1:nlev,1:ncell)
      pstate%omega_at_pc_full_level%f(1:nlev,1:ncell)                 =        scalar_omega_at_pc_full_level(1:nlev,1:ncell)
      pstate%temp_at_pc_full_level%f(1:nlev,1:ncell)                  =          scalar_temp_at_pc_full_level(1:nlev,1:ncell)
      pstate%tracer_mxrt_at_pc_full_level%f(1:ntracer,1:nlev,1:ncell) = scalar_qqq_at_pc_full_level(1:ntracer,1:nlev,1:ncell)
!
! evaluate pressure related quantities
!
      pstate%pressure_at_pc_surface%f(1:ncell)            = scalar_mpressure_at_pc_face_level(nlevp,1:ncell)
      pstate%pressure_at_pc_face_level%f(1:nlevp,1:ncell) = scalar_mpressure_at_pc_face_level(1:nlevp,1:ncell)
      pstate%pressure_at_pc_full_level%f(1:nlev ,1:ncell) = scalar_mpressure_at_pc_full_level(1:nlev ,1:ncell)
      pstate%delp_at_pc_full_level%f(1:nlev,1:ncell)      = scalar_mpressure_at_pc_face_level(2:nlevp,1:ncell)-scalar_mpressure_at_pc_face_level(1:nlev,1:ncell)
      do ilev = 1, nlev
         pstate%exner_at_pc_full_level%f(ilev,1:ncell)    = (p00/pstate%pressure_at_pc_full_level%f(ilev,1:ncell))**(rdry/cp)
      end do
!
! surface
!
!
      pstate%geop_at_pc_surface%f(1:ncell)          = scalar_geopotential_at_pc_face_level(nlevp,1:ncell)
!
! z and geop is zero in CAM5 by definition
!
      pstate%z_at_pc_face_level%f(nlevp,1:ncell)    = zero
      pstate%geop_at_pc_face_level%f(nlevp,1:ncell) = zero
!
! height with surface substracted
!
      do ilev = 1, nlev
         pstate%z_at_pc_full_level%f(ilev,1:ncell)    =(scalar_geopotential_at_pc_full_level(ilev,1:ncell)-pstate%geop_at_pc_surface%f(1:ncell))/gravity
         pstate%z_at_pc_face_level%f(ilev,1:ncell)    =(scalar_geopotential_at_pc_face_level(ilev,1:ncell)-pstate%geop_at_pc_surface%f(1:ncell))/gravity
         pstate%geop_at_pc_full_level%f(ilev,1:ncell) = scalar_geopotential_at_pc_full_level(ilev,1:ncell)-pstate%geop_at_pc_surface%f(1:ncell)
         pstate%geop_at_pc_face_level%f(ilev,1:ncell) = scalar_geopotential_at_pc_face_level(ilev,1:ncell)-pstate%geop_at_pc_surface%f(1:ncell)
      end do
!
! static energy (for cam5)
!
      do ilev  = 1, nlev
         pstate%static_energy_at_pc_full_level%f(ilev,1:ncell) = cp*pstate%temp_at_pc_full_level%f(ilev,1:ncell) &
                                                                   +pstate%geop_at_pc_full_level%f(ilev,1:ncell) &
                                                                   !+scalar_geopotential_at_pc_full_level(ilev,1:ncell) &
                                                                   +pstate%geop_at_pc_surface%f(1:ncell)
      end do
! uncomment for regressing AP
      !pstate_cam%pbl_tpert_at_pc_surface%f       = zero
      !pstate_cam%pbl_qpert_at_pc_surface%f       = zero
! only zero is ocean
      pstate%sst_at_pc_surface%f(1:ncell)   = staticData_sst_at_pc_surface%f(1:ncell)

      if(test_real_case)then 
      do icell = 1, ncell
        if(pstate%landfrac_at_pc_surface%f(icell).le.0.5_r8)then
! for seaice conc >0.01 (threshold), combine ts with sst and sicetemp following Rongxy's implementation
            if(staticData_sic_at_pc_surface%f(icell) .ge. 0.01_r8) then
                pstate%ts_at_pc_surface%f(icell) = staticData_sst_at_pc_surface%f(icell)*(one-staticData_sic_at_pc_surface%f(icell))+&
                                                   sicetemp_at_pc_surface%f(icell)           *staticData_sic_at_pc_surface%f(icell)
            else
                pstate%ts_at_pc_surface%f(icell) = staticData_sst_at_pc_surface%f(icell)
            end if
         end if
      end do

      else  
      where(pstate%landfrac_at_pc_surface%f .le. 0.5_r8)pstate%ts_at_pc_surface%f = staticData_sst_at_pc_surface%f
      end if  

      return
  end subroutine dtp_interface_fill_pstate

#endif

!------------------------------------------------------------
! Transform physics tendencies to model-dynamics compatible
! style, this is for general purposese, should be regressed
! to the before state given the same config (2019-Nov)
!------------------------------------------------------------

  subroutine dtp_interface_convert_tend_phy2dyn(mesh, ntracer, nlev, ncell ,         &
                                  scalar_mpressure_at_pc_full_level,   &
                                  scalar_tracer_mxrt_at_pc_full_level, & ! dry mxrt
                                  scalar_tracer_sphm_at_pc_full_level, & ! moist mxrt, for regression
                                  scalar_temp_at_pc_full_level,&
                                  scalar_mif_at_pc_full_level, &
                                  tend_duDt_at_pc_full_level,  &
                                  tend_dvDt_at_pc_full_level,  &
                                  tend_dtDt_at_pc_full_level,  &
                                  tend_dsDt_at_pc_full_level,  &
                                  tend_dunDt_at_edge_full_level, &
                                  tend_dptmDt_at_pc_full_level,&
                                  tend_dqDt_at_pc_full_level )
! io
    type(global_domain)  , intent(inout)  :: mesh
    integer(i4),           intent(in)     :: ntracer
    integer(i4),           intent(in)     :: nlev
    integer(i4),           intent(in)     :: ncell
    real(r8)   ,           intent(in)     :: scalar_mpressure_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: scalar_tracer_mxrt_at_pc_full_level(ntracer,nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: scalar_tracer_sphm_at_pc_full_level(ntracer,nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: scalar_temp_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: scalar_mif_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_duDt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_dvDt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_dtDt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(in)     :: tend_dsDt_at_pc_full_level(ntracer,nlev,mesh%nv_full)
    real(r8)   ,           intent(inout)  :: tend_dunDt_at_edge_full_level(nlev,mesh%ne_full)
    real(r8)   ,           intent(inout)  :: tend_dptmDt_at_pc_full_level(nlev,mesh%nv_full)
    real(r8)   ,           intent(inout)  :: tend_dqDt_at_pc_full_level(ntracer,nlev,mesh%nv_full)
! local
    integer(i4)   :: icell1, icell2, ilev, iv, ie
    real(r8)      :: point_cart_coort(3), vector_velocity(3)
    real(r8)      :: tend_U_wind_at_edge_full_level 
    real(r8)      :: tend_V_wind_at_edge_full_level

!
! transform u,v tendency from cell to normal wind tendency at edge
!
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
              tend_dunDt_at_edge_full_level(ilev,ie) = dot_product(vector_velocity, mesh%edp_nr(1:3,ie))
           end do
        end do
!
! find mxrt and thetam tendency
!
        do iv = 1, ncell
           do ilev = 1, nlev
!
! do not change for the time being
!             !------------ for regress tc_a, which is a special case)
              !tend_dqDt_at_pc_full_level(1:ntracer,ilev,iv) = tend_dsDt_at_pc_full_level(1:ntracer,ilev,iv)/&
              !                                               ((one-scalar_tracer_sphm_at_pc_full_level(1,ilev,iv))**2)
              !--------------------------------------------------------------

              tend_dqDt_at_pc_full_level(1:ntracer,ilev,iv) = tend_dsDt_at_pc_full_level(1:ntracer,ilev,iv)/scalar_mif_at_pc_full_level(ilev,iv) ! for general purpose
!
! partial
!
              tend_dptmDt_at_pc_full_level(ilev,iv) = ((p00/scalar_mpressure_at_pc_full_level(ilev,iv))**(rdry/cp))*&
                                                     (tend_dtDt_at_pc_full_level(ilev,iv)*(one+(rvap/rdry)*scalar_tracer_mxrt_at_pc_full_level(1,ilev,iv))+&
                                                      scalar_temp_at_pc_full_level(ilev,iv)*(rvap/rdry)*tend_dqDt_at_pc_full_level(1,ilev,iv))
           end do
        end do

     return
  end subroutine dtp_interface_convert_tend_phy2dyn

 end module grist_dtp_interface_full_physpkg1
