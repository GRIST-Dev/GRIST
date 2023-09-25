
!----------------------------------------------------------------------------
! Created on 2019
! Version 1.0
! Description: This module contains major data structure for the physics driver,
!              which helps to ensure that different physics (dry or moist) impacts
!              the model in a similar manner.
! Revision history:
!              Only data have an I/O requirement (needs type definition) are 
!              contained here; physpkg related data are mostly in their own 
!              data module. YiZ: try to minimize data array in this module.
!   202108: ntracer is the number for dynamics and PDC, nspecies are the
!   number for all species in model physics; ntracer<=nspecies, not all need
!   to be advected and PDC. Array allocated for nspecies
!----------------------------------------------------------------------------

 module grist_physics_data_structure

   use grist_constants,     only: r8, i4, rdry, cp, p00, zero, mwh2o, mwdry
   use grist_domain_types,  only: global_domain
   use grist_data_types,    only: scalar_1d_field, scalar_2d_field, scalar_3d_field, &
                                  wrap_allocate_data1d,   wrap_allocate_data2d,   wrap_allocate_data3d, &
                                  wrap_deallocate_data1d, wrap_deallocate_data2d, wrap_deallocate_data3d
   use grist_nml_module,    only: nlev, nlevp, ntracer, nspecies
   use grist_handle_error,  only: endrun
 
   implicit none
   private

   public :: ptend_f1,         & ! used by dtp
             ptend_f2,         & ! used by dtp
             ptend_rk,         & ! used by dtp
             ptend_f3,         & ! used inside physics
             pstate  ,         &
             phy_tracer_info,                        &
             get_tracer_information,                 &
             grist_physics_data_structure_construct, &
             grist_physics_data_structure_destruct
!
! these common data were in scm_comm_, moved here for a clean add
!

    real(r8),    public        :: tground
    logical,     public        :: have_shflx
    logical,     public        :: have_lhflx
    logical,     public        :: have_tg
!
! all physics data are defined at primal cell (voronoi)
!

    type physics_tend
! from model physics
         type(scalar_2d_field) :: tend_u_wind_at_pc_full_level
         type(scalar_2d_field) :: tend_v_wind_at_pc_full_level
         type(scalar_2d_field) :: tend_potential_temp_at_pc_full_level ! pt or ptm
         type(scalar_2d_field) :: tend_temp_at_pc_full_level
         type(scalar_3d_field) :: tend_tracer_mxrt_at_pc_full_level    ! q or s
! to model dynamics
         type(scalar_2d_field) :: tend_normal_velocity_at_edge_full_level
         type(scalar_2d_field) :: tend_mass_pt_at_pc_full_level
         type(scalar_3d_field) :: tend_tracer_mass_at_pc_full_level
         type(scalar_2d_field) :: tend_www_at_pc_face_level
    end type physics_tend

!=========================================================
! Variables only for physics
! Note: all physics variables are located at primal cell
! LiXiaohan, 2019
! This has been seperated to physics_state and physics_buff
! physics state contains the type that can be shared by
! different physics modules (yizhang 2019)
!=========================================================

    type physics_state
! topo
         type(scalar_1d_field) :: geop_at_pc_surface 
! 3d dynamics to physics
#if (defined AMIPC_PHYSICS || defined SCM_PHYSICS)
         type(scalar_2d_field) :: u_wind_at_pc_full_level
         type(scalar_2d_field) :: v_wind_at_pc_full_level
         type(scalar_2d_field) :: omega_at_pc_full_level
         type(scalar_2d_field) :: temp_at_pc_full_level
         type(scalar_2d_field) :: static_energy_at_pc_full_level
         type(scalar_3d_field) :: tracer_mxrt_at_pc_full_level
         type(scalar_2d_field) :: z_at_pc_full_level
         type(scalar_2d_field) :: z_at_pc_face_level
         type(scalar_2d_field) :: geop_at_pc_full_level
         type(scalar_2d_field) :: geop_at_pc_face_level
         type(scalar_2d_field) :: exner_at_pc_full_level
         type(scalar_1d_field) :: pressure_at_pc_surface
         type(scalar_2d_field) :: pressure_at_pc_face_level
         type(scalar_2d_field) :: pressure_at_pc_full_level
         type(scalar_2d_field) :: delp_at_pc_full_level
         type(scalar_2d_field) :: delp_at_pc_face_level
#endif

! land model
         type(scalar_1d_field) :: snowhland_at_pc_surface   ! snow depth (liquid water equivalent) over land 
         type(scalar_1d_field) :: snowhice_at_pc_surface    ! snow depth (liquid water equivalent) over ice
         type(scalar_1d_field) :: landfrac_at_pc_surface    ! land area fraction
         type(scalar_1d_field) :: icefrac_at_pc_surface     ! sea-ice areal fraction
         type(scalar_1d_field) :: ocnfrac_at_pc_surface     ! ocean areal fraction
         type(scalar_1d_field) :: ts_at_pc_surface          ! merged surface temperature [K]
         type(scalar_1d_field) :: sst_at_pc_surface         ! sea surface temperature    [K]
         type(scalar_1d_field) :: ustar_at_pc_surface       ! ustar

#ifdef AMIPW_PHYSICS
! zhangyi: below used for storing psurf_wrf output and restart
         type(scalar_1d_field) :: qfx_at_pc_surface         ! znt
         type(scalar_1d_field) :: znt_at_pc_surface         ! znt
         type(scalar_1d_field) :: mol_at_pc_surface         ! mol
         type(scalar_1d_field) :: pblh_at_pc_surface        ! pblh
         type(scalar_1d_field) :: mavail_at_pc_surface      ! pblh
         type(scalar_1d_field) :: scalar_rainc_surface      ! accumulated convective rain [mm]
         type(scalar_1d_field) :: scalar_rainnc_surface     ! accumulated non-convective rain[mm]
         type(scalar_1d_field) :: scalar_snownc_surface     ! accumulated non-conv snow [mm]
         type(scalar_1d_field) :: scalar_grapnc_surface     ! accumulated non-conv graupel [mm]
         type(scalar_1d_field) :: scalar_raincv_surface     ! convective rain [mm] per model_step (fast-physics step)
#endif

! ---------------------------------------------------------------------------

! atm in
         type(scalar_1d_field) :: atm_in_taux_at_pc_surface ! surface momentum flux in the longitudinal direction [N/m2]
         type(scalar_1d_field) :: atm_in_tauy_at_pc_surface ! surface momentum flux in the latitudinal direction [N/m2]
         type(scalar_1d_field) :: atm_in_shflx_at_pc_surface! sensible heat flux at the surface [w/m2]
         type(scalar_2d_field) :: atm_in_qflx_at_pc_surface ! surface constituent flux [kg/m2/s]
         type(scalar_1d_field) :: atm_in_lwup_at_pc_surface ! upward longwave heat flux 
         type(scalar_1d_field) :: atm_in_asdir_at_pc_surface! 0.2-0.7 micro-meter srfc alb: direct rad 
         type(scalar_1d_field) :: atm_in_asdif_at_pc_surface! 0.2-0.7 micro-meter srfc alb: diffuse rad 
         type(scalar_1d_field) :: atm_in_aldir_at_pc_surface! 0.7-5.0 micro-meter srfc alb: direct rad 
         type(scalar_1d_field) :: atm_in_aldif_at_pc_surface! 0.7-5.0 micro-meter srfc alb: diffuse rad 

! atm out
         type(scalar_1d_field) :: atm_out_sols_at_pc_surface  ! Direct solar rad on surface (< 0.7)
         type(scalar_1d_field) :: atm_out_soll_at_pc_surface  ! Direct solar rad on surface (>= 0.7)
         type(scalar_1d_field) :: atm_out_solsd_at_pc_surface ! Diffuse solar rad on surface (< 0.7)
         type(scalar_1d_field) :: atm_out_solld_at_pc_surface ! Diffuse solar rad on surface (>= 0.7)
         type(scalar_1d_field) :: atm_out_flwds_at_pc_surface ! Down longwave flux at surface
         type(scalar_1d_field) :: atm_out_fswds_at_pc_surface ! Down shortwave flux at surface
         type(scalar_1d_field) :: atm_out_netsw_at_pc_surface ! Surface solar absorbed flux

! precipitation & snow
         type(scalar_1d_field) :: scalar_precl_surface      ! large scale precipitation flux at surface [m/s]
         type(scalar_1d_field) :: scalar_precc_surface      ! convective (shallow+deep) precipitation flux at surface [m/s]
         type(scalar_1d_field) :: scalar_prect_surface      ! total precipitation flux at surface [m/s]
         type(scalar_1d_field) :: scalar_snowl_surface      ! large scale snow flux at surface [m/s]
         type(scalar_1d_field) :: scalar_snowc_surface      ! (wrf-phys has zero)'SNOW_SH' + 'SNOW_DP' convective (shallow+deep) snow flux at surface [m/s]
         type(scalar_1d_field) :: scalar_grapl_surface      ! grid-scale grap 

    end type physics_state

    type(physics_tend)         :: ptend_rk
    type(physics_tend)         :: ptend_f1
    type(physics_tend)         :: ptend_f2
    type(physics_tend)         :: ptend_f3
    type(physics_state)        :: pstate

!================================================
! Information of tracers
! LiXiaohan, 2019
!================================================

    type tracer_info
         character(256)        :: longname
         real(r8)              :: molec_weight
         real(r8)              :: qmin
         logical               :: cnst_fixed_ubc
         logical               :: cnst_fixed_ubflx
    end type tracer_info

    type(tracer_info), allocatable  :: phy_tracer_info(:)

  contains

  subroutine grist_physics_data_structure_construct(mesh)

    type(global_domain), intent(in) :: mesh
!
! all physics tendency is collocated at cell
!
      call wrap_allocate_data2d(mesh%nv, nlev, ptend_f1%tend_u_wind_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, ptend_f1%tend_v_wind_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, ptend_f1%tend_potential_temp_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, ptend_f1%tend_temp_at_pc_full_level)
      call wrap_allocate_data3d(mesh%nv, nlev, nspecies, ptend_f1%tend_tracer_mxrt_at_pc_full_level)
!
! for dycore
!
      ! only this is edge dim
      if(.not.allocated(ptend_f1%tend_normal_velocity_at_edge_full_level%f))then
         allocate(ptend_f1%tend_normal_velocity_at_edge_full_level%f(nlev,mesh%ne))
         ptend_f1%tend_normal_velocity_at_edge_full_level%f   = zero
         ptend_f1%tend_normal_velocity_at_edge_full_level%pos = 6
      end if

      call wrap_allocate_data2d(mesh%nv, nlev,          ptend_f1%tend_mass_pt_at_pc_full_level)
      call wrap_allocate_data3d(mesh%nv, nlev, nspecies,ptend_f1%tend_tracer_mass_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlevp,         ptend_f1%tend_www_at_pc_face_level)
!
! statś
!
      call wrap_allocate_data1d(mesh%nv,       pstate%geop_at_pc_surface)
#if (defined AMIPC_PHYSICS || defined SCM_PHYSICS)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%u_wind_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%v_wind_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%omega_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%temp_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%static_energy_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%z_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlevp,pstate%z_at_pc_face_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%geop_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlevp,pstate%geop_at_pc_face_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%exner_at_pc_full_level)
      call wrap_allocate_data1d(mesh%nv,       pstate%pressure_at_pc_surface)
      call wrap_allocate_data2d(mesh%nv, nlevp,pstate%pressure_at_pc_face_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%pressure_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlev, pstate%delp_at_pc_full_level)
      call wrap_allocate_data2d(mesh%nv, nlevp,pstate%delp_at_pc_face_level)
      call wrap_allocate_data3d(mesh%nv, nlev, ntracer, pstate%tracer_mxrt_at_pc_full_level)
#endif
!
! surface
!
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_precl_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_precc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_prect_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_snowl_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_snowc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_grapl_surface)

      call wrap_allocate_data1d(mesh%nv, pstate%snowhland_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%ts_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%sst_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%ustar_at_pc_surface)

#ifdef AMIPW_PHYSICS
      call wrap_allocate_data1d(mesh%nv, pstate%qfx_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%znt_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%mol_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%pblh_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%mavail_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_rainc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_rainnc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_snownc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_grapnc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%scalar_raincv_surface)
#endif
      call wrap_allocate_data1d(mesh%nv, pstate%landfrac_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%ocnfrac_at_pc_surface)
#if (defined AMIPC_PHYSICS || defined AMIPW_PHYSICS)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_taux_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_tauy_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_lwup_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_shflx_at_pc_surface)
      call wrap_allocate_data2d(mesh%nv, ntracer, pstate%atm_in_qflx_at_pc_surface)

      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_asdir_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_asdif_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_aldir_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_in_aldif_at_pc_surface)

      call wrap_allocate_data1d(mesh%nv, pstate%atm_out_flwds_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_out_fswds_at_pc_surface)
      call wrap_allocate_data1d(mesh%nv, pstate%atm_out_netsw_at_pc_surface)
#endif

       ptend_f2 = ptend_f1
       ptend_rk = ptend_f1
       ptend_f3 = ptend_f1

!==================================================
! Information of travers (used by SCM/CAM_physics)
! LiXH
!==================================================

      allocate(phy_tracer_info(ntracer))

      return
    end subroutine grist_physics_data_structure_construct

    subroutine grist_physics_data_structure_destruct
!
! all physics tendency is collocated at cell
!
      call wrap_deallocate_data2d(ptend_f1%tend_u_wind_at_pc_full_level)
      call wrap_deallocate_data2d(ptend_f1%tend_v_wind_at_pc_full_level)
      call wrap_deallocate_data2d(ptend_f1%tend_potential_temp_at_pc_full_level)
      call wrap_deallocate_data2d(ptend_f1%tend_temp_at_pc_full_level)
      call wrap_deallocate_data3d(ptend_f1%tend_tracer_mxrt_at_pc_full_level)
!
! for dycore usage
!
      call wrap_deallocate_data2d(ptend_f1%tend_normal_velocity_at_edge_full_level)
      call wrap_deallocate_data2d(ptend_f1%tend_mass_pt_at_pc_full_level)
      call wrap_deallocate_data3d(ptend_f1%tend_tracer_mass_at_pc_full_level)
      call wrap_deallocate_data2d(ptend_f1%tend_www_at_pc_face_level)

! state
      call wrap_deallocate_data1d(pstate%geop_at_pc_surface)
#if (defined AMIPC_PHYSICS || defined SCM_PHYSICS)
      call wrap_deallocate_data2d(pstate%u_wind_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%v_wind_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%omega_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%temp_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%static_energy_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%z_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%z_at_pc_face_level)
      call wrap_deallocate_data2d(pstate%geop_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%geop_at_pc_face_level)
      call wrap_deallocate_data2d(pstate%exner_at_pc_full_level)
      call wrap_deallocate_data1d(pstate%pressure_at_pc_surface)
      call wrap_deallocate_data2d(pstate%pressure_at_pc_face_level)
      call wrap_deallocate_data2d(pstate%pressure_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%delp_at_pc_full_level)
      call wrap_deallocate_data2d(pstate%delp_at_pc_face_level)
      call wrap_deallocate_data3d(pstate%tracer_mxrt_at_pc_full_level)
#endif

      call wrap_deallocate_data1d(pstate%scalar_precl_surface)
      call wrap_deallocate_data1d(pstate%scalar_precc_surface)
      call wrap_deallocate_data1d(pstate%scalar_prect_surface)
      call wrap_deallocate_data1d(pstate%scalar_snowc_surface)
      call wrap_deallocate_data1d(pstate%scalar_snowl_surface)
      call wrap_deallocate_data1d(pstate%scalar_grapl_surface)
!
! for tracer information, LiXH
!
      deallocate(phy_tracer_info)
!
! for physical parameterizations, LiXH, these are init inside physpkg
!
      call wrap_deallocate_data1d(pstate%atm_in_taux_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_in_tauy_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_in_shflx_at_pc_surface)
      call wrap_deallocate_data2d(pstate%atm_in_qflx_at_pc_surface)

      call wrap_deallocate_data1d(pstate%snowhland_at_pc_surface)
      call wrap_deallocate_data1d(pstate%snowhice_at_pc_surface)
      call wrap_deallocate_data1d(pstate%landfrac_at_pc_surface)
      call wrap_deallocate_data1d(pstate%icefrac_at_pc_surface)
      call wrap_deallocate_data1d(pstate%ocnfrac_at_pc_surface)

      call wrap_deallocate_data1d(pstate%ts_at_pc_surface)
      call wrap_deallocate_data1d(pstate%sst_at_pc_surface)
      call wrap_deallocate_data1d(pstate%ustar_at_pc_surface)
#ifdef AMIPW_PHYSICS
      call wrap_deallocate_data1d(pstate%qfx_at_pc_surface)
      call wrap_deallocate_data1d(pstate%znt_at_pc_surface)
      call wrap_deallocate_data1d(pstate%mol_at_pc_surface)
      call wrap_deallocate_data1d(pstate%pblh_at_pc_surface)
      call wrap_deallocate_data1d(pstate%mavail_at_pc_surface)
      call wrap_deallocate_data1d(pstate%scalar_rainc_surface)
      call wrap_deallocate_data1d(pstate%scalar_rainnc_surface)
      call wrap_deallocate_data1d(pstate%scalar_snownc_surface)
      call wrap_deallocate_data1d(pstate%scalar_grapnc_surface)
      call wrap_deallocate_data1d(pstate%scalar_raincv_surface)
#endif

      call wrap_deallocate_data1d(pstate%atm_in_lwup_at_pc_surface)

      call wrap_deallocate_data1d(pstate%atm_in_asdir_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_in_asdif_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_in_aldir_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_in_aldif_at_pc_surface)

      call wrap_deallocate_data1d(pstate%atm_out_sols_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_out_soll_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_out_solsd_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_out_solld_at_pc_surface)

      call wrap_deallocate_data1d(pstate%atm_out_flwds_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_out_fswds_at_pc_surface)
      call wrap_deallocate_data1d(pstate%atm_out_netsw_at_pc_surface)

      return
    end subroutine grist_physics_data_structure_destruct

    subroutine get_tracer_information(physpkg,sub_physpkg)
#ifdef AMIPW_PHYSICS
    use grist_PhysW_data_structure, only: p_qv, p_qc, p_qr, p_qi, p_qs, p_qg
#endif
    ! io
    character(len=*),intent(in)              :: physpkg
    character(len=*),intent(in), optional    :: sub_physpkg
    
    if(trim(physpkg) .eq. 'AMIPC_PHYSICS' .and. trim(sub_physpkg) .eq. 'none') then
        ! CAM physical_pkg contains 5 tracers:
        ! Specific humidity, cloud liquid, cloud ice, cloud liquid number, cloud ice number
        if(ntracer .ne. 5) call endrun("ntracer should be equal to 5")
        call add_tracer(1, 'specific_humidity',   mwh2o, 1.E-12_r8)
        call add_tracer(2, 'cloud_liquid',        mwdry, 0._r8    )
        call add_tracer(3, 'cloud_ice',           mwdry, 0._r8    )
        call add_tracer(4, 'cloud_liquid_number', mwdry, 0._r8    )
        call add_tracer(5, 'cloud_ice_number',    mwdry, 0._r8    )

    elseif(trim(physpkg) .eq. 'AMIPC_PHYSICS' .and. trim(sub_physpkg) .eq. 'DCMIP2016-TC') then
        if(ntracer .ne. 1) call endrun("ntracer should be equal to 1")
        call add_tracer(1, 'specific_humidity',   mwh2o, 1.E-12_r8)

#ifdef AMIPW_PHYSICS
    elseif(trim(physpkg) .eq. 'AMIPW_PHYSICS')then
        if(p_qv.le.ntracer)call add_tracer(p_qv, 'water_vapor', mwh2o, 1.E-20_r8     )
        if(p_qc.le.ntracer)call add_tracer(p_qc, 'cloud_water', mwdry, 0._r8         )
        if(p_qr.le.ntracer)call add_tracer(p_qr, 'rain_water',  mwdry, 0._r8         )
        if(p_qi.le.ntracer)call add_tracer(p_qi, 'cloud_ice',   mwdry, 0._r8         )
        if(p_qs.le.ntracer)call add_tracer(p_qs, 'snow',        mwdry, 0._r8         )
        if(p_qg.le.ntracer)call add_tracer(p_qg, 'graupel',     mwdry, 0._r8         )
#endif

    else
        ! other phy_pkg tracers will be added here.
        ! not completed, Lixh
    end if

    end subroutine get_tracer_information

    subroutine add_tracer(ind, longname, molec_weight, qmin,        &
                          fixed_ubc, fixed_ubflx)
      integer(i4),intent(in)           :: ind
      character(len=*),intent(in)      :: longname
      real(r8),intent(in)              :: molec_weight
      real(r8),intent(in)              :: qmin
      logical, optional                :: fixed_ubc
      logical, optional                :: fixed_ubflx
 
      phy_tracer_info(ind)%longname     = longname
      phy_tracer_info(ind)%molec_weight = molec_weight
      phy_tracer_info(ind)%qmin         = qmin

      if ( present(fixed_ubc) ) then
          phy_tracer_info(ind)%cnst_fixed_ubc   = fixed_ubc
      else
          phy_tracer_info(ind)%cnst_fixed_ubc   = .false.
      end if

      if ( present(fixed_ubflx) ) then
          phy_tracer_info(ind)%cnst_fixed_ubflx = fixed_ubflx
      else
          phy_tracer_info(ind)%cnst_fixed_ubflx = .false.
      end if
 
    end subroutine add_tracer

 end module grist_physics_data_structure
