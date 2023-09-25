
module grist_lsm_noahmp_vars
    use grist_constants,                    only: i4, r8, deg2rad

    implicit none
    public
    save

    integer,public:: ids,ide,jds,jde,kds,kde
    integer,public:: ims,ime,jms,jme,kms,kme
    integer,public:: its,ite,jts,jte,kts,kte
    !=================================================================================================================
    !.. variables and arrays related to land-surface parameterization:
    !=================================================================================================================

    logical,parameter:: &
    ua_phys=.false.    !option to activate UA Noah changes: a different snow-cover physics in the land-surface
                        !scheme. That option is not currently implemented in MPAS.

    integer,parameter:: &
    opt_thcnd = 1      !option to treat thermal conductivity in NoahLSM (new option implemented in WRF 3.8.0).
                        != 1, original (default).
                        != 2, McCumber and Pielke for silt loam and sandy loam.

    integer,parameter:: &
    fasdas = 0         !for WRF surface data assimilation system (not used in MPAS).

    integer,public:: &
    sf_surface_physics !used to define the land surface scheme by a number instead of name. It
                        !is only needed in module_ra_rrtmg_sw.F to define the spectral surface
                        !albedos as functions of the land surface scheme.

    integer,public:: &
    num_soils          !number of soil layers                                                                  [-]
    
    integer,dimension(:,:),allocatable:: &
    isltyp_2d,         &!dominant soil type category                                                            [-]
    ivgtyp_2d           !dominant vegetation category                                                           [-]

    real(r8),dimension(:),allocatable:: &
    dzs_2d              !thickness of soil layers                                                               [m]
    real(r8),dimension(:,:,:),allocatable:: &
    smcrel_2d,         &!soil moisture threshold below which transpiration starts to stress                     [-]
    sh2o_2d,           &!unfrozen soil moisture content                                       [volumetric fraction]
    smois_2d,          &!soil moisture                                                        [volumetric fraction]
    tslb_2d             !soil temperature                                                                       [K]

    real(r8),dimension(:,:),allocatable:: &
    acsnom_2d,         &!accumulated melted snow                                                           [kg m-2]
    acsnow_2d,         &!accumulated snow                                                                  [kg m-2]
    canwat_2d,         &!canopy water                                                                      [kg m-2]
    chklowq_2d,        &!surface saturation flag                                                                [-]
    grdflx_2d,         &!ground heat flux                                                                   [W m-2]
    lai_2d,            &!leaf area index                                                                        [-]
    noahres_2d,        &!residual of the noah land-surface scheme energy budget                             [W m-2]
    potevp_2d,         &!potential evaporation                                                              [W m-2]
    qz0_2d,            &!specific humidity at znt                                                         [kg kg-1]
    rainbl_2d,         &!
    sfcrunoff_2d,      &!surface runoff                                                                     [m s-1]
    shdmin_2d,         &!minimum areal fractional coverage of annual green vegetation                           [-]
    shdmax_2d,         &!maximum areal fractional coverage of annual green vegetation                           [-]
    smstav_2d,         &!moisture availability                                                                  [-]
    smstot_2d,         &!total moisture                                                                    [m3 m-3]
    snopcx_2d,         &!snow phase change heat flux                                                        [W m-2]
    snotime_2d,        &!
    snowc_2d,          &!snow water equivalent                                                             [kg m-2]
    snowh_2d,          &!physical snow depth                                                                    [m]
    swdown_2d,         &!downward shortwave flux at the surface                                             [W m-2]
    udrunoff_2d,       &!sub-surface runoff                                                                 [m s-1]
    tmn_2d,            &!soil temperature at lower boundary                                                     [K]
    vegfra_2d,         &!vegetation fraction                                                                    [-]
    z0_2d               !background roughness length                                                            [m]

    real(r8),dimension(:,:),allocatable:: &
    alswvisdir_2d,     &!direct-beam surface albedo in visible spectrum                                         [-]
    alswvisdif_2d,     &!diffuse-beam surface albedo in visible spectrum                                        [-]
    alswnirdir_2d,     &!direct-beam surface albedo in near-IR spectrum                                         [-]
    alswnirdif_2d       !diffuse-beam surface albedo in near-IR spectrum                                        [-]

    !.. arrays needed to run UA Noah changes (different snow-cover physics):
    real(r8),dimension(:,:),allocatable:: &
    flxsnow_2d,        &!energy added to sensible heat flux when ua_phys=true                               [W m-2]
    fvbsnow_2d,        &!fraction of vegetation with snow beneath when ua_phys=true                             [-]
    fbursnow_2d,       &!fraction of canopy buried when ua_phys=true                                            [-]
    fgsnsnow_2d         !fraction of ground snow cover when ua_phys=true                                        [-]

    !.. arrays needed in the argument list in the call to the Noah LSM urban parameterization: note that these arrays
    !.. are initialized to zero since we do not run an urban model:
    integer,dimension(:,:),allocatable:: &
    utype_urb_2d        !urban type                                                                             [-]

    real(r8),dimension(:,:),allocatable:: &
    frc_urb_2d,        &!urban fraction                                                                         [-]
    ust_urb_2d          !urban u* in similarity theory                                                        [m/s]

    !.. arrays needed in the argument list in the call to the Noah LSM hydrology model: note that these arrays are
    !.. initialized to zero since we do not run a hydrology model:
    real(r8),dimension(:,:),allocatable:: &
    infxsrt_2d,        &!timestep infiltration excess                                                          [mm]
    sfcheadrt_2d,      &!surface water detph                                                                   [mm]
    soldrain_2d         !soil column drainage                                                                  [mm]

    !=================================================================================================================
    !.. variables and arrays related to surface characteristics:
    !=================================================================================================================

    real(r8),dimension(:,:),allocatable:: &
    xlat_2d,           &!longitude, west is negative                                                      [degrees]
    xlon_2d             !latitude, south is negative                                                      [degrees]

    real(r8),dimension(:,:),allocatable:: &
    sfc_albedo_2d,     &!surface albedo                                                                         [-]
    sfc_albbck_2d,     &!surface background albedo                                                              [-]
    sfc_emibck_2d,     &!land surface background emissivity                                                     [-]
    sfc_emiss_2d,      &!land surface emissivity                                                                [-]
    snoalb_2d,         &!annual max snow albedo                                                                 [-]
    snow_2d,           &!snow water equivalent                                                             [kg m-2]
    tsk_2d,            &!surface-skin temperature                                                               [K]
    sst_2d,            &!sea-surface temperature                                                                [K]
    xice_2d,           &!ice mask                                                                               [-]
    xland_2d            !land mask    (1 for land; 2 for water)                                                 [-]

    !=================================================================================================================
    !.. variables needed for the surface layer scheme and land surface scheme when config_frac_seaice
    !   is set to true. the arrays below have the same definition as the corresponding "_2d" arrays:
    !=================================================================================================================

    real(r8),dimension(:,:),allocatable:: br_sea,ch_sea,chs_sea,chs2_sea,cpm_sea,cqs2_sea,      &
                                          flhc_sea,flqc_sea,gz1oz0_sea,hfx_sea,lh_sea,mavail_sea,mol_sea,      &
                                          psih_sea,psim_sea,fh_sea,fm_sea,qfx_sea,qgh_sea,qsfc_sea,regime_sea, &
                                          rmol_sea,ust_sea,wspd_sea,znt_sea,zol_sea,tsk_sea,xland_sea
    real(r8),dimension(:,:),allocatable:: t2m_sea,th2m_sea,q2_sea,u10_sea,v10_sea
    real(r8),dimension(:,:),allocatable:: cd_sea,cda_sea,ck_sea,cka_sea,ustm_sea

    real(r8),dimension(:,:),allocatable:: regime_hold
    real(r8),dimension(:,:),allocatable:: tsk_ice







    































end module grist_lsm_noahmp_vars