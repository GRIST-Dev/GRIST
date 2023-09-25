
!==========================================================================
! This module is separated from grist_physics_data_structure, and only
! contains physics info specific to PhysC model physics
! The original grist_physics_data_structure can be shared by different
! model physics packages
! Yi Zhang 2019
!==========================================================================

 module grist_PhysC_data_structure

   use grist_constants,     only: r8, i4, rdry, cp, p00, zero, mwh2o, mwdry
   use grist_domain_types,  only: global_domain
   use grist_data_types,    only: scalar_1d_field, scalar_2d_field, scalar_3d_field
   use grist_nml_module,    only: nlev, nlevp, ntracer
   use grist_handle_error,  only: endrun
 
   implicit none
   private

   public :: pstate_cam,               &
             tend_in_cam_physics,      &
             ptend_vertical_diffusion, &
             ptend_deep_convection,    &
             ptend_deep_convection_2,  &
             ptend_shallow_convection, &
             ptend_radiation,          &
             ptend_microphysics,       &
             ptend_microp_aero,        &
             ptend_macrophysics,       &
             ptend_dry_adjustment,     &
             grist_PhysC_data_structure_construct, &
             grist_PhysC_data_structure_destruct             

!
! This varaibles are inited inside specific physics module
! currently it is related to CAM-style Physics, similar to
! the original buffer style vars; other physpkg can also
! define such type for their own usage
!
    type tracer_data
        character(16)          :: name
        integer                :: idx 
        real(r8),dimension(:,:), allocatable :: f
    end type tracer_data
 
 
    type physics_state_cam
        
! vertical diffusion
         type(scalar_2d_field) :: pbl_tke_at_pc_face_level  ! 'tke'   turbulent kinetic energy [ m2/s2 ]
         type(scalar_2d_field) :: pbl_kvm_at_pc_face_level  ! 'kvm'   eddy diffusivity for momentum [ m2/s ]
         type(scalar_2d_field) :: pbl_kvh_at_pc_face_level  ! 'kvh'   eddy diffusivity for heat [m2/s]
         type(scalar_2d_field) :: pbl_kvt_at_pc_face_level  ! 'kvt'   Kinematic molecular conductivity
         type(scalar_1d_field) :: pbl_pblh_at_pc_surface    ! 'pblh'  planetary boundary layer height [m]
         type(scalar_1d_field) :: pbl_tpert_at_pc_surface   ! 'tpert' perturbation temperature [K]
         type(scalar_1d_field) :: pbl_qpert_at_pc_surface   ! 'qpert' perturbation specific humidity [kg/kg]
         type(scalar_1d_field) :: pbl_tauresx_at_pc_surface ! 'tauresx'
         type(scalar_1d_field) :: pbl_tauresy_at_pc_surface ! 'tauresy'
         ! for Lin Macro ------------------------>
         type(scalar_2d_field) :: pbl_lengi_at_pc_face_level
         type(scalar_2d_field) :: pbl_shi_at_pc_face_level
         type(scalar_2d_field) :: pbl_smi_at_pc_face_level
         type(scalar_1d_field) :: pbl_wstarPBL
         ! for Lin Macro <------------------------

! gwd
         type(scalar_1d_field) :: sgh30_at_pc_surface       ! 'sgh30'    

! shallow and deep convection
         type(scalar_2d_field) :: sh_shfrc_at_pc_full_level ! 'shfrc'   UW Shallow cumulus cloud fraction
         type(scalar_2d_field) :: sh_icwmr_at_pc_full_level ! 'ICWMRSH' In cloud water + ice mixing ratio
         type(scalar_2d_field) :: dp_icwmr_at_pc_full_level ! 'ICWMRDP' In cloud water + ice mixing ratio
         type(scalar_1d_field) :: cumulus_cldtop            ! 'CLDTOP'  Cumulus cloud top interface index
         type(scalar_1d_field) :: cumulus_cldbot            ! 'CLDBOT'  Cumulus cloud botom interface index
         type(scalar_1d_field) :: cumulus_cldtop2           ! 'CLDTOP2' Shallow Cumulus cloud top interface index
         type(scalar_1d_field) :: cumulus_cldbot2           ! 'CLDBOT2' Shallow Cumulus cloud botom interface index
 
         type(scalar_1d_field) :: pcape_before_dycore       ! LiXH add, CAPE before dynamic core, following Xie et al. (2019) 
         type(scalar_1d_field) :: pcape_previous            ! LiXH add, CAPE at the last time step
         type(scalar_1d_field) :: cumulus_cldtop_pressure
         type(scalar_1d_field) :: cumulus_cldbot_pressure
         type(scalar_2d_field) :: rprddp_at_pc_full_level   ! 'RPRDDP' rain production rate in deep convection
         type(scalar_2d_field) :: rprdsh_at_pc_full_level   ! 'RPRDSH' rain production rate in shallow convection
         type(scalar_2d_field) :: rprdtot_at_pc_full_level  ! 'RPRDTOT' rain production rate
         type(scalar_2d_field) :: updraft_mass_flux         ! 'cmfmc' Moist deep + shallow convection cloud mass flux [ kg/s/m2 ]
         type(scalar_2d_field) :: cumulus_cush              ! 'cush'  Convective scale height [ m ]

         ! for Lin Macro ------------------------>
         type(scalar_2d_field) :: sh_qtu_shal_at_pc_face_level
         type(scalar_2d_field) :: sh_thlu_shal_at_pc_face_level
         type(scalar_2d_field) :: sh_umf_shal_at_pc_face_level
         ! for Lin Macro <------------------------

! macrophysics & stratiform
         type(scalar_2d_field) :: macrop_cmeliq_at_pc_full_level ! 'CMELIQ' Net condensation rate of stratiform liquid [kg/kg/s] 
         type(scalar_2d_field) :: macrop_sgm_at_pc_full_level    ! sub-grid scale standard deviation for GaussPDF macrophysics 

         !multi time levels 
         type(scalar_3d_field) :: macrop_cld_at_pc_full_level    ! 'CLD' Total cloud fraction
         type(scalar_3d_field) :: macrop_concld_at_pc_full_level ! 'CONCLD' Convective (deep+shallow) cloud fraction
         type(scalar_3d_field) :: macrop_ast_at_pc_full_level    ! 'AST' Stratus cloud fraction

         type(scalar_3d_field) :: macrop_alst_at_pc_full_level   ! 'ALST' Physical liquid stratus fraction 
         type(scalar_3d_field) :: macrop_aist_at_pc_full_level   ! 'AIST' Physical ice stratus fraction
         type(scalar_3d_field) :: macrop_qlst_at_pc_full_level   ! 'QLST' Physical in-cloud LWC
         type(scalar_3d_field) :: macrop_qist_at_pc_full_level   ! 'QIST' Physical in-cloud IWC

         type(scalar_3d_field) :: macrop_qcwat_at_pc_full_level  ! 'QCWAT' Cloud water old q
         type(scalar_3d_field) :: macrop_lcwat_at_pc_full_level  ! 'LCWAT' Cloud liquid water old q
         type(scalar_3d_field) :: macrop_iccwat_at_pc_full_level ! 'ICCWAT' Cloud ice water old q
         type(scalar_3d_field) :: macrop_nlwat_at_pc_full_level  ! 'NLWAT' Cloud liquid droplet number condentration. old.
         type(scalar_3d_field) :: macrop_niwat_at_pc_full_level  ! 'NIWAT' Cloud ice droplet number condentration. old.
         type(scalar_3d_field) :: macrop_tcwat_at_pc_full_level  ! 'TCWAT' Cloud water old temperature

! microphysics
         type(scalar_2d_field) :: microp_wsedl_at_pc_full_level  ! 'WSEDL' Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
         type(scalar_2d_field) :: microp_dei_at_pc_full_level    ! 'DEI' Ice effective diameter (meters) (AG: microns?)
         type(scalar_2d_field) :: microp_des_at_pc_full_level    ! 'DES' Snow effective diameter
         type(scalar_2d_field) :: microp_res_at_pc_full_level    !       Snow effective radius    Linhan: for RRTMG 4DDA
         type(scalar_2d_field) :: microp_mu_at_pc_full_level     ! 'MU' Size distribution shape parameter for radiation
         type(scalar_2d_field) :: microp_lambdac_at_pc_full_level! 'LAMBDAC' Size distribution slope parameter for radiation
         type(scalar_2d_field) :: microp_rei_at_pc_full_level    ! 'REI' Ice effective drop size (microns)
         type(scalar_2d_field) :: microp_rel_at_pc_full_level    ! 'REL' Liquid effective drop size (microns)
         type(scalar_2d_field) :: microp_rel_fn_at_pc_full_level ! 'REL_FN' Liquid effective drop size (microns)
         type(scalar_2d_field) :: microp_iciwp_at_pc_full_level  ! 'ICIWP' In-cloud ice water path for radiation
         type(scalar_2d_field) :: microp_iclwp_at_pc_full_level  ! 'ICLWP' In-cloud liquid water path for radiation
         type(scalar_2d_field) :: microp_icswp_at_pc_full_level  ! 'ICSWP' In-cloud snow water path

         type(scalar_2d_field) :: microp_qrout_at_pc_full_level
         type(scalar_2d_field) :: microp_nrout_at_pc_full_level

         type(scalar_2d_field) :: microp_lwc_at_pc_full_level    ! 'LWC' Grid box average liquid water content 
         type(scalar_2d_field) :: microp_iwc_at_pc_full_level    ! 'IWC+SWC' Grid box average ice water and snow content, only for io

         !multi time levels 
         type(scalar_3d_field) :: microp_cc_t_at_pc_full_level     ! 'CC_T' Grid-mean microphysical tendency
         type(scalar_3d_field) :: microp_cc_qv_at_pc_full_level    ! 'CC_qv' Grid-mean microphysical tendency
         type(scalar_3d_field) :: microp_cc_ql_at_pc_full_level    ! 'CC_ql' Grid-mean microphysical tendency
         type(scalar_3d_field) :: microp_cc_qi_at_pc_full_level    ! 'CC_qi' Grid-mean microphysical tendency
         type(scalar_3d_field) :: microp_cc_nl_at_pc_full_level    ! 'CC_nl' Grid-mean microphysical tendency
         type(scalar_3d_field) :: microp_cc_ni_at_pc_full_level    ! 'CC_ni' Grid-mean microphysical tendency
         type(scalar_3d_field) :: microp_cc_qlst_at_pc_full_level  ! 'CC_qlst' In-liquid stratus microphysical tendency
         type(scalar_3d_field) :: microp_cldo_at_pc_full_level     ! 'CLDO' Old cloud fraction
         type(scalar_3d_field) :: microp_cldfsnow_at_pc_full_level ! 'CLDFSNOW' Cloud fraction for liquid+snow

! microphysics aero
         type(scalar_2d_field) :: microp_areo_naai                 ! 'NAAI' number of activated aerosol for ice nucleation
         type(scalar_2d_field) :: microp_areo_naai_hom             ! 'NAAI_HOM' number of activated aerosol for ice nucleation 
                                                                   !            (homogeneous freezing only)
         type(scalar_2d_field) :: microp_areo_npccn                ! 'NPCCN' number of CCN (liquid activated)
         type(scalar_3d_field) :: microp_areo_rndst                ! 'RNDST' radius of 4 dust bins for contact freezing
         type(scalar_3d_field) :: microp_areo_nacon                ! 'NACON' number in 4 dust bins for contact freezing 

! cloud
         type(scalar_2d_field) :: cld_sh_frac_at_pc_full_level   ! 'SH_FRAC' Shallow convective cloud fraction 
         type(scalar_2d_field) :: cld_dp_frac_at_pc_full_level   ! 'DP_FRAC' Deep convective cloud fraction

! radiation
         type(scalar_2d_field) :: lw_qrl_at_pc_full_level   ! 'QRL' lw radiative heating rate
         type(scalar_2d_field) :: lw_qrlc_at_pc_full_level  ! 'QRLC' clearsky lw radiative heating rate
         type(scalar_1d_field) :: flwdsc_at_pc_surface      ! 'FSDSC' Clear sky surface downward lw flux
         type(scalar_1d_field) :: flwusc_at_pc_surface      ! 'FSUSC' Clear sky surface upward lw flux
         type(scalar_1d_field) :: flwdtc_at_pc_top          ! 'FSDTC' Clear sky downward lw flux at TOA
         type(scalar_1d_field) :: flwutc_at_pc_top          ! 'FSUTC' Clear sky upweward lw flux at TOA
         type(scalar_1d_field) :: flwds_at_pc_surface       ! 'FSDS'  Flux lw downward surface  
         type(scalar_1d_field) :: flwus_at_pc_surface       ! 'FSUS'  Flux lw upward surface
         type(scalar_1d_field) :: flwdt_at_pc_top           ! 'FSDT'  Flux lw downward at TOA
         type(scalar_1d_field) :: flwut_at_pc_top           ! 'FSUT'  Flux lw upward at TOA

         type(scalar_2d_field) :: sw_qrs_at_pc_full_level   ! 'QRS' sw radiative heating rate
         type(scalar_2d_field) :: sw_qrsc_at_pc_full_level  ! 'QRSC' clearsky sw radiative heating rate
         type(scalar_1d_field) :: fswdsc_at_pc_surface      ! 'FSDSC' Clear sky surface downwelling solar flux
         type(scalar_1d_field) :: fswusc_at_pc_surface      ! 'FSUSC' Clear sky surface upwelling solar flux
         type(scalar_1d_field) :: fswdtc_at_pc_top          ! 'FSDTC' Clear sky downwelling solar flux at TOA
         type(scalar_1d_field) :: fswutc_at_pc_top          ! 'FSUTC' Clear sky upwelling solar flux at TOA
         type(scalar_1d_field) :: fswds_at_pc_surface       ! 'FSDS'  Flux shortwave downwelling surface  
         type(scalar_1d_field) :: fswus_at_pc_surface       ! 'FSUS'  Flux shortwave upwelling surface
         type(scalar_1d_field) :: fswdt_at_pc_top           ! 'FSDT'  Flux shortwave downwelling at TOA
         type(scalar_1d_field) :: fswut_at_pc_top           ! 'FSUT'  Flux shortwave upwelling at TOA

         type(scalar_1d_field) :: fsntoa_at_pc_top          ! 'FSNTOA' Net solar flux at TOA
         type(scalar_1d_field) :: fsntoac_at_pc_top         ! 'FSNTOAC'Clear sky Net solar flux at TOA

         type(scalar_1d_field) :: lwcf_at_pc_top            ! 'LWCF' longwave cloud forcing  
         type(scalar_1d_field) :: swcf_at_pc_top            ! 'SWCF' shortwave cloud forcing

! precipitation
         type(scalar_1d_field) :: str_prec_surface          ! 'PREC_STR' Surface flux of precipitation from stratiform [m/s]
         type(scalar_1d_field) :: str_snow_surface          ! 'SNOW_STR' Surface flux of snow from stratiform [m/s]
         type(scalar_1d_field) :: sed_prec_surface          ! 'PREC_SED' Surface flux of total cloud water from sedimentation [m/s]
         type(scalar_1d_field) :: sed_snow_surface          ! 'SNOW_SED' Surface flux of total cloud ice from sedimentation [m/s]
         type(scalar_1d_field) :: pcw_prec_surface          ! 'PREC_PCW' Surface flux of precipitation from microphysics [m/s]
         type(scalar_1d_field) :: pcw_snow_surface          ! 'SNOW_PCW' Surface flux of snow from microphysics [m/s]

! prec
         type(scalar_1d_field) :: scalar_precc_sh_surface   ! 'PREC_SH' shallow convective precipitation flux at surface [m/s]
         type(scalar_1d_field) :: scalar_precc_dp_surface   ! 'PREC_DP' deep convective precipitation flux at surface [m/s]
         type(scalar_1d_field) :: scalar_snowc_surface      ! 'SNOW_SH+DP' deep convective precipitation flux at surface [m/s]

! aerosol
         type(scalar_3d_field) :: aerosol_fracis            ! 'FRACIS'  fraction of transported species that are insoluble
         type(scalar_3d_field) :: aerosol_dgnum             ! 'DGNUM'
         type(scalar_3d_field) :: aerosol_dgnumwet          ! 'DGNUMWET'
         type(scalar_3d_field) :: aerosol_wetdens_ap        ! 'WETDENS_AP'
         type(scalar_3d_field) :: aerosol_qaerwat           ! 'QAERWAT'


! ghg_data: N2O, CH4, CHC11, CFC12, CO2, O2, O3
         type(tracer_data),dimension(:),allocatable :: ghg_at_pc_full_level
         integer  :: total_ghg_num

! aerosol
         type(tracer_data),dimension(:),allocatable :: aerosol_at_pc_full_level
         integer  :: total_aerosol_num


! diagnostics
         type(scalar_1d_field) :: diag_cloud_tot        
         type(scalar_1d_field) :: diag_cloud_low
         type(scalar_1d_field) :: diag_cloud_med
         type(scalar_1d_field) :: diag_cloud_hgh
         type(scalar_1d_field) :: diag_z_at_500hpa          ! height at 500 hPa (m)
         type(scalar_1d_field) :: diag_u_at_850hpa          ! zonal wind at 850 hPa (m/s)
         type(scalar_1d_field) :: diag_u_at_200hpa          ! zonal wind at 200 hPa (m/s)
         type(scalar_1d_field) :: diag_v_at_850hpa          ! Meridional wind at 850 hPa (m/s)
         type(scalar_2d_field) :: diag_relhum               ! 'RELHUM' relativite humidity
         type(scalar_1d_field) :: diag_tmq                  ! 'TMQ' Total (vertically integrated) precipitable water
         type(scalar_1d_field) :: diag_tgliqwp              ! 'TGCLDLWP' Vertically integrated liquid water path
         type(scalar_1d_field) :: diag_tgicewp              ! 'TGCLDIWP' Vertically integrated ice water path

    end type physics_state_cam
 
    type tend_in_cam_physics
        type(scalar_2d_field)  :: tend_u
        type(scalar_2d_field)  :: tend_v
        type(scalar_2d_field)  :: tend_s
        type(scalar_3d_field)  :: tend_q
    end type tend_in_cam_physics

    type(physics_state_cam)    :: pstate_cam

    type(tend_in_cam_physics)      :: ptend_vertical_diffusion
    type(tend_in_cam_physics)      :: ptend_deep_convection
    type(tend_in_cam_physics)      :: ptend_deep_convection_2
    type(tend_in_cam_physics)      :: ptend_shallow_convection
    type(tend_in_cam_physics)      :: ptend_radiation
    type(tend_in_cam_physics)      :: ptend_macrophysics
    type(tend_in_cam_physics)      :: ptend_microphysics
    type(tend_in_cam_physics)      :: ptend_microp_aero
    type(tend_in_cam_physics)      :: ptend_dry_adjustment


  contains

  subroutine grist_PhysC_data_structure_construct(mesh)

    type(global_domain), intent(in) :: mesh

!
! physics tendency variables 
!
      allocate(ptend_vertical_diffusion%tend_u%f(nlev,mesh%nv))
      allocate(ptend_vertical_diffusion%tend_v%f(nlev,mesh%nv))
      allocate(ptend_vertical_diffusion%tend_s%f(nlev,mesh%nv))
      allocate(ptend_vertical_diffusion%tend_q%f(ntracer,nlev,mesh%nv))

      allocate(ptend_deep_convection%tend_u%f(nlev,mesh%nv))
      allocate(ptend_deep_convection%tend_v%f(nlev,mesh%nv))
      allocate(ptend_deep_convection%tend_s%f(nlev,mesh%nv))
      allocate(ptend_deep_convection%tend_q%f(ntracer,nlev,mesh%nv))

      allocate(ptend_deep_convection_2%tend_q%f(ntracer,nlev,mesh%nv))

      allocate(ptend_shallow_convection%tend_u%f(nlev,mesh%nv))
      allocate(ptend_shallow_convection%tend_v%f(nlev,mesh%nv))
      allocate(ptend_shallow_convection%tend_s%f(nlev,mesh%nv))
      allocate(ptend_shallow_convection%tend_q%f(ntracer,nlev,mesh%nv))

      allocate(ptend_radiation%tend_s%f(nlev,mesh%nv))

      allocate(ptend_microphysics%tend_s%f(nlev,mesh%nv))
      allocate(ptend_microphysics%tend_q%f(ntracer,nlev,mesh%nv))

      allocate(ptend_microp_aero%tend_q%f(ntracer,nlev,mesh%nv))

      allocate(ptend_macrophysics%tend_s%f(nlev,mesh%nv))
      allocate(ptend_macrophysics%tend_q%f(ntracer,nlev,mesh%nv))

      allocate(ptend_dry_adjustment%tend_s%f(nlev,mesh%nv))
      allocate(ptend_dry_adjustment%tend_q%f(ntracer,nlev,mesh%nv))

! init
      ptend_vertical_diffusion%tend_u%f                    = zero
      ptend_vertical_diffusion%tend_v%f                    = zero
      ptend_vertical_diffusion%tend_s%f                    = zero
      ptend_vertical_diffusion%tend_q%f                    = zero
      ptend_deep_convection%tend_u%f                       = zero
      ptend_deep_convection%tend_v%f                       = zero
      ptend_deep_convection%tend_s%f                       = zero
      ptend_deep_convection%tend_q%f                       = zero
      ptend_deep_convection_2%tend_q%f                     = zero
      ptend_shallow_convection%tend_u%f                    = zero
      ptend_shallow_convection%tend_v%f                    = zero
      ptend_shallow_convection%tend_s%f                    = zero
      ptend_shallow_convection%tend_q%f                    = zero
      ptend_radiation%tend_s%f                             = zero
      ptend_microphysics%tend_s%f                          = zero
      ptend_microphysics%tend_q%f                          = zero
      ptend_microp_aero%tend_q%f                           = zero
      ptend_macrophysics%tend_s%f                          = zero
      ptend_macrophysics%tend_q%f                          = zero
      ptend_dry_adjustment%tend_s%f                        = zero
      ptend_dry_adjustment%tend_q%f                        = zero

! location
      ptend_vertical_diffusion%tend_u%pos                  = 0
      ptend_vertical_diffusion%tend_v%pos                  = 0
      ptend_vertical_diffusion%tend_s%pos                  = 0
      ptend_vertical_diffusion%tend_q%pos                  = 0
      ptend_deep_convection%tend_u%pos                     = 0
      ptend_deep_convection%tend_v%pos                     = 0
      ptend_deep_convection%tend_s%pos                     = 0
      ptend_deep_convection%tend_q%pos                     = 0
      ptend_deep_convection_2%tend_q%pos                   = 0
      ptend_shallow_convection%tend_u%pos                  = 0
      ptend_shallow_convection%tend_v%pos                  = 0
      ptend_shallow_convection%tend_s%pos                  = 0
      ptend_shallow_convection%tend_q%pos                  = 0
      ptend_radiation%tend_s%pos                           = 0
      ptend_microphysics%tend_s%pos                        = 0
      ptend_microphysics%tend_q%pos                        = 0
      ptend_microp_aero%tend_q%pos                         = 0
      ptend_macrophysics%tend_s%pos                        = 0
      ptend_macrophysics%tend_q%pos                        = 0
      ptend_dry_adjustment%tend_s%pos                      = 0
      ptend_dry_adjustment%tend_q%pos                      = 0

      return
    end subroutine grist_PhysC_data_structure_construct

    subroutine grist_PhysC_data_structure_destruct()

!
! for tendency inside physics package, LiXH
!
      deallocate(ptend_vertical_diffusion%tend_u%f)
      deallocate(ptend_vertical_diffusion%tend_v%f)
      deallocate(ptend_vertical_diffusion%tend_s%f)
      deallocate(ptend_vertical_diffusion%tend_q%f)
      deallocate(ptend_deep_convection%tend_u%f)
      deallocate(ptend_deep_convection%tend_v%f)
      deallocate(ptend_deep_convection%tend_s%f)
      deallocate(ptend_deep_convection%tend_q%f)
      deallocate(ptend_deep_convection_2%tend_q%f)
      deallocate(ptend_shallow_convection%tend_u%f)
      deallocate(ptend_shallow_convection%tend_v%f)
      deallocate(ptend_shallow_convection%tend_s%f)
      deallocate(ptend_shallow_convection%tend_q%f)
      deallocate(ptend_radiation%tend_s%f)
      deallocate(ptend_microphysics%tend_s%f)
      deallocate(ptend_microphysics%tend_q%f)
      deallocate(ptend_microp_aero%tend_q%f)
      deallocate(ptend_macrophysics%tend_s%f)
      deallocate(ptend_macrophysics%tend_q%f)
      deallocate(ptend_dry_adjustment%tend_s%f)
      deallocate(ptend_dry_adjustment%tend_q%f)

      if(allocated(pstate_cam%str_prec_surface%f))                 &
        deallocate(pstate_cam%str_prec_surface%f)

      if(allocated(pstate_cam%str_snow_surface%f))                 &
        deallocate(pstate_cam%str_snow_surface%f)

      if(allocated(pstate_cam%sed_prec_surface%f))                 &
        deallocate(pstate_cam%sed_prec_surface%f)

      if(allocated(pstate_cam%sed_snow_surface%f))                 &
        deallocate(pstate_cam%sed_snow_surface%f)

      if(allocated(pstate_cam%pcw_prec_surface%f))                 &
        deallocate(pstate_cam%pcw_prec_surface%f)

      if(allocated(pstate_cam%pcw_snow_surface%f))                 &
        deallocate(pstate_cam%pcw_snow_surface%f)

      if(allocated(pstate_cam%sgh30_at_pc_surface%f))              &
        deallocate(pstate_cam%sgh30_at_pc_surface%f)

      return
    end subroutine grist_PhysC_data_structure_destruct

 end module grist_PhysC_data_structure
