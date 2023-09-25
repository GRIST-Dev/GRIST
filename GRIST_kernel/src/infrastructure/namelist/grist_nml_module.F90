 !======================================================
 !
 !  Created by zhangyi on 16/8/5.
 !  Latest Updated at DEC, 2018
 !
 !  NAMELIST module for GRIST, also included are some 
 !  global parameters.
 !======================================================

 module grist_nml_module

   use grist_constants,   only: i4, r8, day2sec
   use grist_time_manager,only: get_curr_calday
   use grist_mpi
   use grist_handle_error,only: endrun
   use grist_handle_error,only: delete_file_sub

   implicit none
   
   private

   public  :: nml_set_global_vars,   &
!
! ctl_para: setup related to time, history, run_mode, etc
!
              day_duration,      & ! how many days to run, overwritten if final_ymd final_tod are set
              model_timestep,    & ! for 3d model, it should be equal to the largest timestep of that mode, e.g., dycore_timestep for rkfb
              run_type,          & ! 'initial', 'restart'
              h1_restart_freq,   & ! frequency of writing h1 restart file
              h0_restart_freq,   & ! frequency of writing h0 restart file
              h1_history_freq,   & ! frequency of writing h1 history
              h1_1d_history_freq,& ! frequency of writing h1's 1d var history
              h1_2d_history_freq,& ! frequency of writing h1's 2d var history
              h1_3d_history_freq,& ! frequency of writing h1's 3d var history
              model_instance,    & ! select of model instance, default 'gcm',
              working_mode,      & ! running mode of 3d model, 'see grist_gcm_control', only follow model_instance 'gcm'
              write_history,     & ! write history or not (an old safeguard)
              write_history_h0,  & ! write h0 history file
              write_history_h1,  & ! write h1 history file
              write_history_h1_separate,  & ! write h1 history file
              write_restart_h0,  & ! write restart (n months) or not
              write_restart_h1,  & ! write restart (n-step) or not
              write_verbose,     & ! write verbose debug info
              write_stepinfo,    & ! write step info
              test_real_case,    & ! whether run with real-world initial and bdy conditions, amp mode
              start_ymd,         & ! e.g., 20120526
              start_tod,         & ! 0
              final_ymd,         & ! e.g., 20120528
              final_tod,         & ! 0
              doAquaPlanet,      & ! whether run AquaPlanet, default: .false.
              grid_info,         & ! specify unique grid name for output data, e.g., G6, G5B3X4, X20EA, ...
              comm_group_size,   & ! group size for I/O, see liuz's GMDD paper for info
!
! data_para: setup related to external datasets
!
              outdir,            & ! data ourput directory
              gridFilePath,      & ! grid file directory
              gridFileNameHead,  & ! head name of grid file
              read_partition,    & ! Whether to read partition or not, set to true for high-resolution simulation
              pardir,            & ! The directory containing the partition files (liuz's offline partition)
              initialAtmFilePath,& ! full path+filename
              initialAtmUFilePath,&! in case of large-scale, read seperately 
              initialAtmVFilePath,&! in case of large-scale, read seperately
              initialAtmTFilePath,&! in case of large-scale, read seperately
              initialAtmQFilePath,&! in case of large-scale, read seperately
              large_atm_file_on  ,&
              initialLndFilePath,& ! full path+filename
              initialDataSorc,   & ! initial data source, ERAIP, ERAIM, GFS..
               staticFilePath,   & ! static file: full path+filename
               sftopoFilePath,   & ! a seperate topo file, can be the same as static File, full path+filename, unit: m*gravity
                  sstFilePath,   & ! filepath 
              aqua_sst_style,    & ! style of aqua planet sst
              numMonSST,         & ! number of sst time in a file
              sstFileNameHead,   & ! for climo-12 sst: it is $head.cycle.$tail; for amip sst, it is $head.19x0/20x0.$tail
              sstFileNameTail,   & ! as above
              sstFile_year_beg,  & ! the begining year of this decade in which start_ymd falls, e,g, 1976 shoud have sstFile_year_beg=1970
              real_sst_style,    & ! NEW to replace cycle_sst: 'CYCLE','AMIP','DAILY'
              firstClipSic,      & ! clip sic after just read sic, default: false
              glamAtmFilePath,   & ! only filepath, filename is determined inside code
              glamDataSorc,      & ! bdy data source: 'ERA_ml', 'GFS_ml', 'ERA_Pres', 'MERRA_Pres'
              glamDataInterval,  &
              glamLbcDiffusion,  & ! add F2 term in LBC forcing 
              glamLbcCoef1,      &
              glamLbcCoef2,      &
              glamForcePhi,      &
              glam_use_relaxation,&
!
! swe_para: set up related to shallow water model, only a few shared with 3D model
!
              swe_testcase,          & ! sw testcase
              swe_initialfield,      & ! sw initial field
              swe_timestep,          & ! to control the swm
              swe_time_scheme,       & ! as is
              swe_advection_scheme,  & ! adv
              swe_isolate_advection_test, & ! in isolation mode
              swe_test_dual_cell,    &
              swe_use_limiter_rk3,   &
              swe_use_tspas,         &
              swe_use_streamf,       &
              swe_scalar_type,       &
              swe_limiter_type,      &
              conserve_scheme,       & ! shared with 3D model
              pv_order,              & ! shared with 3D model, default: 2,2,3
              use_tr_mbs,            & ! shared with 3D model
              ke_method,             & ! shared with 3D model
              ke_alpha,              & ! shared with 3D model
              mass_wtd_vor,          & ! shared with 3D model
              modon_radius,          & ! shared with 3D model
              two_polyfit,           &
              one_polyfit,           &
              half_polyfit,          &
!
! dycore: setup related to (dry) dynamics core (main solver)
!
              gcm_testcase,      & ! initial conditions: ''
              mas_adv_flag,      & ! mass advection, default: 2,2,2 (hard coded)
              pot_adv_flag,      & ! pt advection,   default: 33,33,33
              phi_adv_flag,      & ! phi advection,  default: 33,33,33
              www_adv_flag,      & ! www advecetion, default: 33,33,33
              ver_adv_flag,      & ! vertical flux order, default: 3
              hor_pgf_flag,      & ! hori pgd,       default: 6
              nrk,               & ! rk2 or rk3, default: 3
              nh_dynamics,       & ! nh_dynamics flag
              www_damping_coef,  & ! w-damping coef
              d3d_damping_coef,  & ! 3d div damping
              zd              ,  & ! w-damping depth
              imbeta,            & ! implicit weighting coef, default: 0.55
              imkesi,            & ! implicit weighting coef, default: 0.55
              hpgf_pfull_alpha,  &
              hpgf_pface_alpha,  &
              dycore_pgf_kesi,   &
              smg_coef,          & ! square of Smg coefficient
              ref_leng,          & ! for VR
              ko4_coef,          & ! 4th hyper coef
              ko6_coef,          & ! 6rh hyper coef (not needed)
              topo_ko_coef,      & !
              do_dycore_laplacian_2nd , ad_dycore_laplacian_2nd ,&
              do_dycore_laplacian_4th , ad_dycore_laplacian_4th ,&
              do_dycore_laplacian_6th , ad_dycore_laplacian_6th ,&
              do_div2d_o4_damping, div2d_o4_damping_coef, &
              eqs_vert_diff,     & ! default yes
              nexpdif_level,     & ! how many levels to add explicit diffusion, default, as nlev
              nsponge,           & ! how many top sponge layer
              spcoef,            & ! sponge layer coefficient
              tend_nct_once,     & ! only calc nct tend once per rk step, obsolete
              use_laplacian_vi,  & ! use vector identity version of velocity Laplacian, do not recommend
              dycore_timestep,   & ! dycore time step
              profile_type,      & ! profile used in PGF, 'tp1'
              vr_mode,           & ! use vr_mode or not, for dynamics, only dycore is VR sensitive
              smooth_topo,       & ! smooth topo or not
              smooth_type,       & ! smooth type
              nsmooth_topo,      & ! how many passes of smoothing topo
              topo_type,         & ! data for topo height source: static, init, usgs
              use_expdif_tstep,  & ! use explicit diffusion in tracer step rather than in model_step (more stable when t-p is split)
              www_k2coef_type,   & ! what type 2nd Laplacian for www: SPLIKE (sponge-like), SMG
              use_www_hyperDiffusion,      & ! whether to use Laplacian 4th for www, with same coef as ko4_coef and scaled in vr_mode
              restore_hydro,     &
              restore_hydro_minsteps,&
              restore_hydro_intsteps,&
              adjphi,      &
              doNotDiagnose,     & ! do not do diagnose (remove MPI_REDUCE) for dycore, tracer, etc
!
! tracer transport: setup related to passive transport of moist species (2nd 3D solver)
!
              ntracer,              & ! how many tracer
              nmif,                 & ! how many index in mif
              mif_index,            & ! to store indexes in mif, max=6 now
              nrk_hori,             & ! how many rk steps in horizontal (1,2,3)
              nrk_vert,             & ! how many rk steps in vertical
              tracer_hadv_flag,     & ! hori tracer select
              tracer_vadv_flag,     & ! vert tracer select
              tracer_hori_limiter,  & ! hori limiter use or not
              tracer_vert_limiter,  & ! vert limiter use or not
              tracer_timestep,      & ! it should be equal the highest of tracer hori/vert step
              tracer_hori_timestep, & ! hori tracer time step
              tracer_vert_timestep, & ! vert tracer time step
              tracer_hvsplit_method,& ! h->v or v->h, default h->v: 1
              tracer_vert_aimp,     & ! use adaptively implicit or not (only final step)
              tracer_vert_fimp,     & ! use fully implicit when using adaptively implicit (final step) 
              tracer_transport_hori_scheme, & ! name for hori transport
              tracer_transport_vert_scheme, & ! name for vert transport
              tracer_mxrt_fixer,    &   ! use mass fixer or not, default not use
              do_tracer_laplacian_2nd,& ! only retained for DCMIP-superecell
              isolate_tracer_test,  &   ! for isolated tracer test (DCMIP1)
!
! physics (common physics config besides those specifically related to each physpkg)
!
              use_phys,          &  ! true or false
              physpkg,           &  ! CAM5based, AMIPW_PHYSICS, or any simple suite
              sub_physpkg,       &  ! sub-physpkg inside the first physpkg
              ptend_wind_rk_on,  &  ! add rk wind tend in RK integration
              ptend_heat_rk_on,  &  ! add rk heat tend in RK integration
              ptend_dycore_f1_on,&      ! if used, please add "USE_PTENDF1" compiling options
              ptend_dycore_heat_f1_on,& ! if used, please add "USE_PTENDF1" compiling options
              ptend_dycore_wind_f1_on,& ! if used, please add "USE_PTENDF1" compiling options
              ptend_tracer_f1_on,&      ! if used, please add "USE_PTENDF1" compiling options
              ptend_f2_on,       &      ! ptend_f2 operator
              ptend_f2drbl_on,   &      ! dribling f2
              ptendSubDiabPhys,  &      ! only for ndc to substract additional heating within rk heat but for f2
              physics_coupling,  &      ! only for AMIPW physics now
              physics_use_som ,  &      ! when true, to use slab ocean model
              seaice_albedo_scheme   ,  &
              TC_pbl_flag     ,  &      ! used by Reed-Jablonowski simple physics
              nspecies        ,  &      ! how many species
!
! mesh
!
              mesh_nv          , & ! needed for running 3d model
              mesh_glv         , & ! needed for running 3d model
              nlev             , & ! needed for running 3d model
              nlevp            , & ! needed for running 3d model
              nlev_inidata     , & ! needed for running 3d model in real-case
              nlev_lamdata     , & ! needed for running LAM model in real-case 
              ntracer_lamdata  , & ! needed for running LAM model in real-case 
              levsoil          , & !---cheyz, needed for running 3d model
!
! below shared with seq_grist and gg
!
              mesh_directly    , &
              mesh_kind        , &
              mesh_node        , &
              mesh_opt         , &
              mesh_readfile_name,&
              mesh_file_source,  &
              mesh_rotate      , &
              mesh_ref_lat     , &
              mesh_ref_lon     , &
              mesh_tgt_lat     , &
              mesh_tgt_lon     , &
              mesh_bird_ccw_deg, &
              mesh_write_topo  , &
!
              config_edt       , &
              config_edp       , &
              config_plg       , &
              config_tri       , &
              write_cdo_file   , &
              write_metis_file , &
              write_grid_file  , &
!
! mesh variable resolution: shared with gg, not used by ParGRIST runtime
!
              ccvt_alpha       , &
              ccvt_beta        , &
              ccvt_gama        , &
              center_lon       , &
              center_lat       , &

!
! derived and runtime-dependent, no need to set
!
              sec_duration,      & 
              nsteps,            &
              fname_output,      &
              fname_restart,     &
              stencil_width,         &
              stencil_exchange_flag, &
              lbdry_flag,            &
              index_flag              ! Reordering flag for vtx-index of inner domain, to improve the cache hits
!
! CTL
!
   real(r8)           :: day_duration              ! days of duration
   real(r8)           :: model_timestep            ! timestep for integration
   character(100)     :: run_type                  ! 'init' or ' 'restart'
   integer(i4)        :: h1_restart_freq           ! how many steps for outputing restart files
   integer(i4)        :: h0_restart_freq           ! how many INTEGER months for h0 monthly restart
   integer(i4)        :: h1_history_freq           ! how many steps for outputing output files (h1 is default)
   integer(i4)        :: h1_1d_history_freq        ! how many steps for outputing output files (h1's 1d is default)
   integer(i4)        :: h1_2d_history_freq        ! how many steps for outputing output files (h1's 2d is default)
   integer(i4)        :: h1_3d_history_freq        ! how many steps for outputing output files (h1's 3d is default)
   character(100)     :: model_instance            !
   character(100)     :: working_mode              !
   logical            :: write_history
   logical            :: write_restart_h0
   logical            :: write_restart_h1
   logical            :: write_verbose
   logical            :: write_stepinfo
   logical            :: write_history_h0
   logical            :: write_history_h1
   logical            :: write_history_h1_separate
   integer            :: comm_group_size = 1200
   logical            :: test_real_case 
   integer(i4)        :: start_ymd   = 20000601  ! start date in yyyymmdd format
   integer(i4)        :: start_tod   = 0         ! time of day relative to start date
   integer(i4)        :: final_ymd   = 0
   integer(i4)        :: final_tod   = 0         !
   logical            :: doAquaPlanet 
   character(100)     :: grid_info   = 'null'    ! default as null
!
! data
!
   character(1024)    :: outdir                  ! output dierctory
   character(1024)    :: gridFilePath            ! grid dierctory
   character(1024)    :: gridFileNameHead
   character(1024)    :: initialAtmFilePath
   character(1024)    :: initialAtmUFilePath
   character(1024)    :: initialAtmVFilePath
   character(1024)    :: initialAtmTFilePath
   character(1024)    :: initialAtmQFilePath
   logical            :: large_atm_file_on
   character(1024)    :: initialLndFilePath
   character(1024)    :: initialDataSorc
   character(1024)    ::  staticFilePath 
   character(1024)    ::  sftopoFilePath 
   character(1024)    ::     sstFilePath
   character(1024)    ::     sstFileNameHead
   character(1024)    ::     sstFileNameTail
   character(1024)    :: aqua_sst_style            ! RJ, NK_CNTL, NK_PEAK, NK_FLAT, NK_CT5N, NK_1K3K, NK_3KW1
   character(1024)    :: real_sst_style            !
   integer(i4)        :: numMonSST  = 12           ! default, climatological SST 12 months
   integer(i4)        :: sstFile_year_beg          ! for amip sst, 1st year of SST data
   logical            :: firstClipSic  = .false.   !
   logical            :: read_partition = .false.
   character(1024)    :: pardir               ! partition dierctory; only necessary for partition large mesh files such as G9B3 (QU 5km)
   character(1024)    :: glamAtmFilePath
   character(1024)    :: glamDataSorc
   real(r8)           :: glamDataInterval
   logical            :: glamLbcDiffusion  = .true.
   real(r8)           :: glamLbcCoef1 = 0.005
   real(r8)           :: glamLbcCoef2 = 0.005
   logical            :: glamForcePhi = .true.
   logical            :: glam_use_relaxation = .true.
!
! SWM
!
   integer(i4)        :: swe_testcase              ! velocity case
   integer(i4)        :: swe_initialfield          ! initial scalar
   real(r8)           :: swe_timestep              ! timestep for integration
   character(100)     :: swe_time_scheme           ! 'rk3','rk4','ssprk','ab3','abmpc','asselin','lftz'
   integer(i4)        :: swe_advection_scheme      !
   logical            :: swe_isolate_advection_test! only test advection
   logical            :: swe_test_dual_cell        ! only test adv on a dual cell
   logical            :: swe_use_limiter_rk3       ! only test adv on a dual cell
   logical            :: swe_use_tspas             ! use tspas or not
   logical            :: swe_use_streamf
   character(100)     :: swe_scalar_type           ! ways of initials 
   character(100)     :: swe_limiter_type          ! ways of initials 
   character(100)     :: conserve_scheme           ! total-energy conserving in R10
   integer(i4)        :: pv_order(3)               !
   logical            :: use_tr_mbs                ! tr stencil use MBS
   integer(i4)        :: ke_method                 ! 0: r10, 1: s12/g13
   real(r8)           :: ke_alpha                  ! alpha in the s12/g13 formulation
   logical            :: mass_wtd_vor              ! true: advect pvor, false: advect avor
   real(r8)           :: modon_radius              ! initial radius in modon test, unit: m
   logical            :: two_polyfit
   logical            :: one_polyfit
   logical            :: half_polyfit
!
! dycore
!
   character(100)     :: gcm_testcase
   integer(i4)        :: mas_adv_flag
   integer(i4)        :: pot_adv_flag(3) ! for each pass of RK3/2
   integer(i4)        :: phi_adv_flag(3)
   integer(i4)        :: www_adv_flag(3)
   integer(i4)        :: ver_adv_flag
   integer(i4)        :: hor_pgf_flag
   integer(i4)        :: nrk             ! how many pass of RK
   logical            :: nh_dynamics
   real(r8)           :: www_damping_coef
   real(r8)           :: d3d_damping_coef
   real(r8)           :: zd
   real(r8)           :: imbeta
   real(r8)           :: imkesi
   real(r8)           :: hpgf_pfull_alpha = 1._r8 
   real(r8)           :: hpgf_pface_alpha = 1._r8
   real(r8)           :: dycore_pgf_kesi  = 1._r8
   real(r8)           :: smg_coef
   real(r8)           :: ref_leng
   real(r8)           :: ko4_coef
   real(r8)           :: ko6_coef
   real(r8)           :: topo_ko_coef
   logical            :: do_dycore_laplacian_2nd, ad_dycore_laplacian_2nd
   logical            :: do_dycore_laplacian_4th, ad_dycore_laplacian_4th
   logical            :: do_dycore_laplacian_6th, ad_dycore_laplacian_6th
   logical            :: do_div2d_o4_damping
   real(r8)           :: div2d_o4_damping_coef
   logical            :: eqs_vert_diff
   integer(i4)        :: nexpdif_level
   integer(i4)        :: nsponge = 0
   real(r8)           :: spcoef  = 0._r8
   logical            :: tend_nct_once
   logical            :: use_laplacian_vi 
   real(r8)           :: dycore_timestep
   character(100)     :: profile_type
   logical            :: vr_mode
   logical            :: smooth_topo
   character(100)     :: smooth_type  = ''
   integer(i4)        :: nsmooth_topo = 1 ! default 1 pass
   character(100)     :: topo_type  = 'static' ! default from static data
   logical            :: use_expdif_tstep
   character(100)     :: www_k2coef_type
   logical            :: use_www_hyperDiffusion
   logical            :: restore_hydro  = .false.
   logical            :: adjphi   = .false.
   integer(i4)        :: restore_hydro_minsteps = 999999999
   integer(i4)        :: restore_hydro_intsteps = 999999999
   logical            :: doNotDiagnose 
!
! tracer transport
!
   integer(i4)        :: ntracer
   integer(i4)        :: nmif
   integer(i4)        :: mif_index(6)
   integer(i4)        :: nrk_hori
   integer(i4)        :: nrk_vert
   integer(i4)        :: tracer_hadv_flag
   integer(i4)        :: tracer_vadv_flag
   logical            :: tracer_hori_limiter
   logical            :: tracer_vert_limiter
   real(r8)           :: tracer_timestep
   real(r8)           :: tracer_hori_timestep
   real(r8)           :: tracer_vert_timestep
   integer(i4)        :: tracer_hvsplit_method
   logical            :: tracer_vert_aimp
   logical            :: tracer_vert_fimp
   character(100)     :: tracer_transport_hori_scheme
   character(100)     :: tracer_transport_vert_scheme
   logical            :: tracer_mxrt_fixer
   logical            :: do_tracer_laplacian_2nd 
   logical            :: isolate_tracer_test
!
! physics cntl
!
   logical            :: use_phys
   character(100)     :: physpkg
   character(100)     :: sub_physpkg = 'none' ! default is none
   logical            :: ptend_wind_rk_on
   logical            :: ptend_heat_rk_on
   logical            :: ptend_dycore_f1_on
   logical            :: ptend_dycore_heat_f1_on
   logical            :: ptend_dycore_wind_f1_on
   logical            :: ptend_tracer_f1_on
   logical            :: ptend_f2_on
   logical            :: ptend_f2drbl_on
   logical            :: ptendSubDiabPhys
   character(100)     :: physics_coupling    ! coupling inside physics, for wrf_physics now
   logical            :: physics_use_som     ! slab ocean flag 
   character(100)     :: seaice_albedo_scheme = 'ECHAM5'
   logical            :: TC_pbl_flag
   integer(i4)        :: nspecies = 0
!
! MESH
!
   integer(i4)        :: mesh_nv
   integer(i4)        :: mesh_glv
   integer(i4)        :: nlev
   integer(i4)        :: nlevp
   integer(i4)        :: nlev_inidata
   integer(i4)        :: nlev_lamdata
   integer(i4)        :: ntracer_lamdata
   integer(i4)        :: levsoil  !--cheyz

   integer(i4)        :: mesh_directly
   character(100)     :: mesh_kind
   character(100)     :: mesh_node
   character(100)     :: mesh_opt
   character(100)     :: mesh_readfile_name
   character(100)     :: mesh_file_source   ! 'xyz (ASCII)', 'GRIST', 'MPAS'
   logical            :: mesh_rotate
   logical            :: mesh_write_topo
   real(r8)           :: mesh_ref_lat
   real(r8)           :: mesh_ref_lon
   real(r8)           :: mesh_tgt_lat
   real(r8)           :: mesh_tgt_lon
   real(r8)           :: mesh_bird_ccw_deg

   logical            :: config_edp = .true.
   logical            :: config_edt = .true.
   logical            :: config_plg = .true.
   logical            :: config_tri = .true.
   logical            :: write_cdo_file   = .true.
   logical            :: write_metis_file = .false.
   logical            :: write_grid_file  = .true.
!
! CCVT
!
   real(r8)           :: ccvt_alpha
   real(r8)           :: ccvt_beta
   real(r8)           :: ccvt_gama
   real(r8)           :: center_lon
   real(r8)           :: center_lat
!
! derived and runtime-dependent
!
   real(r8)           :: sec_duration              ! duration of integration
   real(r8)           :: nsteps                    ! number of steps
   character(200)     :: fname_output              ! filename for output 
   character(200)     :: fname_restart             ! filename for restart
! although USE_HALO2 is default, it needs to be explicitly define by Macro
#ifdef USE_HALO2
   integer, parameter :: stencil_width = 2
#else
   integer, parameter :: stencil_width = 3
#endif
   logical            :: stencil_exchange_flag = .true.
   logical            :: lbdry_flag
   character(100)     :: index_flag   ! 'bfs', 'hilbert', 'hilbert2', 'morton'

  contains

  subroutine nml_set_global_vars()

!================================================
! global vars have been defined in the header
!================================================

! local
  character(len=300) :: filename
  character(len=300) :: buffer
  integer (i4)       :: fileunit
  integer (i4)       :: leap_year_count, yr
  real(r8)           :: juld2, juld1
  logical            :: file_exists

  namelist /ctl_para/ day_duration,           &
                      model_timestep,         &
                      run_type,               &
                      h1_restart_freq,        &
                      h0_restart_freq,        &
                      h1_history_freq,        &
                      h1_1d_history_freq,     &
                      h1_2d_history_freq,     &
                      h1_3d_history_freq,     &
                      model_instance,         &
                      working_mode,           &
                      write_history,          &
                      write_history_h0,       &
                      write_history_h1,       &
                      write_history_h1_separate,  &
                      write_restart_h0,       &
                      write_restart_h1,       &
                      write_verbose,          &
                      write_stepinfo,         &
                      comm_group_size,        &
                      test_real_case,         &
                      start_ymd,              &
                      start_tod,              &
                      final_ymd,              &
                      final_tod,              &
                      doAquaPlanet,           &
                      grid_info

  namelist /data_para/ outdir,                &
                       gridFilePath,          &
                       gridFileNameHead,      &
                       read_partition,        &
                       pardir,                &
                       initialAtmFilePath,    &
                       initialAtmUFilePath,   &
                       initialAtmVFilePath,   &
                       initialAtmTFilePath,   &
                       initialAtmQFilePath,   &
                       large_atm_file_on,     &
                       initialLndFilePath,    &
                       initialDataSorc,       &
                        staticFilePath,       &
                        sftopoFilePath,       &
                           sstFilePath,       &
                           sstFileNameHead,   &
                           sstFileNameTail,   &
                       sstFile_year_beg,      &
                       aqua_sst_style,        &
                       real_sst_style,        &
                       numMonSST,             &
                       firstClipSic,          &
                       glamAtmFilePath,       &
                       glamDataSorc,          &
                       glamDataInterval,      &
                       glamLbcDiffusion,      &
                       glamLbcCoef1,          &
                       glamLbcCoef2,          &
                       glamForcePhi,          &
                       glam_use_relaxation

  namelist /swe_para/ swe_testcase,           & ! for swe&adv tests
                      swe_initialfield,       &
                      swe_timestep,           &
                      swe_time_scheme,        &
                      swe_advection_scheme,   &
                      swe_isolate_advection_test, &
                      swe_test_dual_cell,     &
                      swe_use_limiter_rk3,    &
                      swe_use_tspas,          &
                      swe_use_streamf,        &
                      swe_scalar_type,        &
                      swe_limiter_type,       &
                      conserve_scheme,        &
                      pv_order      ,         &
                      use_tr_mbs,             &
                      ke_method,              &
                      ke_alpha,               &
                      mass_wtd_vor,           &
                      modon_radius,           &
                      two_polyfit,            &
                      one_polyfit,            &
                      half_polyfit

   namelist /dycore_para/ gcm_testcase   , &
                       mas_adv_flag   , &
                       pot_adv_flag   , &
                       phi_adv_flag   , &
                       www_adv_flag   , &
                       ver_adv_flag   , &
                       hor_pgf_flag   , &
                       nrk            , &
                       nh_dynamics    , &
                       www_damping_coef,&
                       d3d_damping_coef,&
                       zd             , &
                       imbeta         , &
                       imkesi         , &
                       hpgf_pfull_alpha,&
                       hpgf_pface_alpha,&
                       dycore_pgf_kesi, &
                       smg_coef       , &
                       ref_leng       , &
                       ko4_coef       , &
                       ko6_coef       , &
                       topo_ko_coef   , &
                       do_dycore_laplacian_2nd, ad_dycore_laplacian_2nd, &
                       do_dycore_laplacian_4th, ad_dycore_laplacian_4th, &
                       do_dycore_laplacian_6th, ad_dycore_laplacian_6th, &
                       do_div2d_o4_damping, div2d_o4_damping_coef, &
                       eqs_vert_diff  , &
                       nexpdif_level  , &
                       nsponge        , &
                       spcoef         , &
                       tend_nct_once  , &
                       use_laplacian_vi,&
                       dycore_timestep, &
                       profile_type   , &
                       vr_mode        , &
                       smooth_topo    , &
                       smooth_type    , &
                       nsmooth_topo   , &
                       topo_type      , &
                       use_expdif_tstep,&
                       www_k2coef_type, &
                       use_www_hyperDiffusion, &
                       restore_hydro  , &
                       adjphi         , &
                       restore_hydro_minsteps, &
                       restore_hydro_intsteps, &
                       lbdry_flag            , &
                       index_flag            , &
                       doNotDiagnose

   namelist /tracer_para/ ntracer            , &
                          nmif               , &
                          mif_index          , &
                          nrk_hori           , & ! can be set by tracer_transport_hori/vert_scheme
                          nrk_vert           , & ! can be set by tracer_transport_hori/vert_scheme 
                          tracer_hadv_flag   , & ! can be set by tracer_transport_hori/vert_scheme
                          tracer_vadv_flag   , & ! can be set by tracer_transport_hori/vert_scheme
                          tracer_hori_limiter, & ! can be set by tracer_transport_hori/vert_scheme
                          tracer_vert_limiter, & ! can be set by tracer_transport_hori/vert_scheme
                          tracer_timestep    , &
                          tracer_hori_timestep,&
                          tracer_vert_timestep,&
                          tracer_hvsplit_method,&
                          tracer_vert_aimp,    & !
                          tracer_vert_fimp,    &
                          tracer_transport_hori_scheme, &
                          tracer_transport_vert_scheme, &
                          tracer_mxrt_fixer,   &
                          do_tracer_laplacian_2nd,&
                          isolate_tracer_test

   namelist /physics_para/use_phys,           &
                          physpkg,            &
                          sub_physpkg,        &
                          ptend_wind_rk_on,   &
                          ptend_heat_rk_on,   &
                          ptend_dycore_f1_on, &
                          ptend_dycore_heat_f1_on, &
                          ptend_dycore_wind_f1_on, &
                          ptend_tracer_f1_on, &
                          ptend_f2_on,        &
                          ptend_f2drbl_on,    &
                          ptendSubDiabPhys,   &
                          physics_coupling,   &
                          physics_use_som,    & ! use slab ocean or not
                          seaice_albedo_scheme,&
                          TC_pbl_flag,        &
                          nspecies

   namelist /mesh_para/ mesh_nv,              &
                        mesh_glv,             &
                        nlev,                 &
                        nlevp,                &
                        nlev_inidata,         &
                        nlev_lamdata,         &
                        ntracer_lamdata,      &
                        levsoil,              & !--cheyz

! below is not required by ParGRIST
                        mesh_directly,        &
                        mesh_kind,            &
                        mesh_node,            &
                        mesh_opt,             &
                        mesh_readfile_name,   &
                        mesh_file_source,     &
                        mesh_rotate,          &
                        mesh_write_topo,      &
                        mesh_ref_lat,         &
                        mesh_ref_lon,         &
                        mesh_tgt_lat,         &
                        mesh_tgt_lon,         &
                        mesh_bird_ccw_deg,    &
!
                        config_edt,           &
                        config_edp,           &
                        config_tri,           &
                        config_plg,           &
                        write_cdo_file,       &
                        write_grid_file,      &
                        write_metis_file

   namelist /ccvt_para/ ccvt_alpha,           &
                        ccvt_beta,            &
                        ccvt_gama,            &
                        center_lon,           &
                        center_lat
! above is not required by ParGRIST

    filename = "grist.nml"

    if(mpi_rank() .eq. 0)then

       print*,"----------------------------------------------------------"
       print*,"   GRIST Framework's parameter file is: ", trim(filename)
       print*,"----------------------------------------------------------"

    end if

!================================================
! pre-setup in case of w/o reading
!================================================
! Default should regress old dtp results
! confirmed at 20210509, cma-pi, dcmip-tc, f1_rk and f2_sudden coupled NDC with simple-physics

    nmif         = 1
    mif_index(1) = 1

#ifdef AMIPW_PHYSICS
! set based on G6 model climate (ice-phase mxrt (WSM6) unreal)
! (known behaviors of its impact, small but do exist)
    nmif         = 3
    mif_index(1) = 1
    mif_index(2) = 2
    mif_index(3) = 3
    mif_index(4) = 999 ! if used, blow it up
    mif_index(5) = 999
    mif_index(6) = 999
#endif

#ifdef AMIPC_PHYSICS
! qv,qc,qi,nc,ni (MG-2m), just a typical setup for CAM5 climate runs
    nmif         = 3
    mif_index(1) = 1 
    mif_index(2) = 2
    mif_index(3) = 3
    mif_index(4) = 999
    mif_index(5) = 999 
#endif

!
! default set these to standard O3
!
    pv_order(1)     = 3;  pv_order(2)     = 3;  pv_order(3)     = 3
    pot_adv_flag(1) = 33; pot_adv_flag(2) = 33; pot_adv_flag(3) = 33
    phi_adv_flag(1) = 33; phi_adv_flag(2) = 33; phi_adv_flag(3) = 33
    www_adv_flag(1) = 33; www_adv_flag(2) = 33; www_adv_flag(3) = 33

!================================================
!           A parameters file must exist
!================================================

    fileunit = 1

    open  (fileunit, status='old',file=filename)
    read  (fileunit, nml=ctl_para)
    read  (fileunit, nml=data_para)
    read  (fileunit, nml=swe_para)
    read  (fileunit, nml=dycore_para)
    read  (fileunit, nml=tracer_para)
    read  (fileunit, nml=physics_para)
    read  (fileunit, nml=mesh_para)
    !read  (fileunit, nml=ccvt_para)
    close (fileunit)

!================================================
!                  calculated
!================================================
! if set final_ymd, this will overwrite day_duration set in namelist
! or, just use day_duration from namelist
    if(final_ymd.ne.0)then
       call get_curr_calday(start_ymd,start_tod,0,0._r8,juld1)
       call get_curr_calday(final_ymd,final_tod,0,0._r8,juld2)
#ifdef USE_LEAP_YEAR
       leap_year_count = 0
       do yr = start_ymd/10000,final_ymd/10000-1
         if((mod(yr,4).eq.0.and.mod(yr,100).ne.0).or.(mod(yr,400).eq.0))then
            leap_year_count = leap_year_count+1
         end if   
       end do
#endif
       day_duration = (final_ymd/10000-start_ymd/10000)*365+leap_year_count+juld2-juld1
    end if

    sec_duration = day_duration*day2sec
    nsteps       = sec_duration/model_timestep

! default, nspecies is equal to ntracer; only for AMIPW_Physics now
    if(nspecies.eq.0) nspecies = ntracer
! now force this to be NLEV
    nexpdif_level = nlev

!==============================================================
!                      For NH-specific
!==============================================================

    if(nh_dynamics)then
       if(trim(physpkg).eq.'AMIPW_PHYSICS')then
          physics_coupling  = 'P3'
          ptendSubDiabPhys  = .true.
       end if 
       if(trim(physpkg).eq.'AMIPC_PHYSICS')then
          print*,"The present AMIPC_PHYSICS can not be coupled to NH-Dynamics yet"
          stop
       end if
       if(test_real_case) adjphi  = .true.
#ifdef LAM_DOMAIN
! not a must but for our DCMIPBW regression
       adjphi  = .true.
#endif
    end if

!==============================================================
! For tracer transport, we use string name to define some 
! available options for user friendly
!==============================================================

    select case (tracer_transport_hori_scheme)
    case('rk3o2')
        nrk_hori = 3;   tracer_hadv_flag = 2;  tracer_hori_limiter = .true.
    case('rk3o3')
        nrk_hori = 3;   tracer_hadv_flag = 33;  tracer_hori_limiter= .true.
    case('rk3o4')
        nrk_hori = 3;   tracer_hadv_flag = 4;  tracer_hori_limiter = .true.
    case('rk3o5')
        nrk_hori = 3;   tracer_hadv_flag = 8;  tracer_hori_limiter = .true.
        if(stencil_width.eq.2) call endrun('rk3o5 can not be used with 2 halo layers')
    case('rk3o6')
        nrk_hori = 3;   tracer_hadv_flag = 9;  tracer_hori_limiter = .true.
        if(stencil_width.eq.2) call endrun('rk3o6 can not be used with 2 halo layers')
    case('ffsl')
        nrk_hori = 1;   tracer_hadv_flag = 28;  tracer_hori_limiter = .true.
    case('tspas')
        nrk_hori = 1;   tracer_hadv_flag = 99;  tracer_hori_limiter = .false.
    case('rk2o2')
        nrk_hori = 2;   tracer_hadv_flag = 2;  tracer_hori_limiter = .true.
    case('rk2o3')
        nrk_hori = 2;   tracer_hadv_flag = 33;  tracer_hori_limiter= .true.
    case('rk2o4')
        nrk_hori = 2;   tracer_hadv_flag = 4;  tracer_hori_limiter = .true.
    case('rk2o5')
        nrk_hori = 2;   tracer_hadv_flag = 8;  tracer_hori_limiter = .true.
        if(stencil_width.eq.2) call endrun('rk3o5 can not be used with 2 halo layers')
    case('rk2o6')
        nrk_hori = 2;   tracer_hadv_flag = 9;  tracer_hori_limiter = .true.
        if(stencil_width.eq.2) call endrun('rk3o6 can not be used with 2 halo layers')
    case default
        nrk_hori = 3;   tracer_hadv_flag = 33;  tracer_hori_limiter = .true.
    end select

    select case (tracer_transport_vert_scheme)
    case('rk3o3')
        nrk_vert = 3;   tracer_vadv_flag = 3;  tracer_vert_limiter = .true.; tracer_vert_aimp = .false.
    case('ppm')
        nrk_vert = 1;   tracer_vadv_flag = 9;  tracer_vert_limiter = .true.; tracer_vert_aimp = .false.
    case('rk3o3-aimp')
        nrk_vert = 3;   tracer_vadv_flag = 3;  tracer_vert_limiter = .true.; tracer_vert_aimp = .true.
    case('rk3o3-fimp')
        nrk_vert = 3;   tracer_vadv_flag = 3;  tracer_vert_limiter = .true.; tracer_vert_aimp = .true. ; tracer_vert_fimp = .true.
    case('ppm-aimp')
        nrk_vert = 1;   tracer_vadv_flag = 9;  tracer_vert_limiter = .true.; tracer_vert_aimp = .true.
    case default
! rk3o3
        nrk_vert = 3;   tracer_vadv_flag = 3;  tracer_vert_limiter = .true.; tracer_vert_aimp = .false.
    end select


    if(mpi_rank() .eq. 0)then
       print*,"--------------------------------------------------------"
       print*,"                                                        "
       print*,"              Key Parameters Input                      "
       print*,"                                                        "
       print*,"        Day_duration:"    , day_duration
       print*,"        Timestep:"        , model_timestep
       print*,"        PV:"              , pv_order(1:3)
       print*,"        masAdv:"          , mas_adv_flag
       print*,"        potAdv:"          , pot_adv_flag(1:3)
       print*,"        wwwAdv:"          , www_adv_flag(1:3)
       print*,"        phiAdv:"          , phi_adv_flag(1:3)
       print*,"        swe_use_limiter:" , swe_use_limiter_rk3
       print*,"        Use_tr_mbs:"      , use_tr_mbs
       print*,"        SWE Use_tspas:"   , swe_use_tspas
       print*,"        Use_nh:"          , nh_dynamics
       print*,"        mass_wtd_vor:"    , mass_wtd_vor
       print*,"        Test_dual_cell:"  , swe_test_dual_cell
       print*,"        mif_index is:"    , mif_index(1:nmif)
       print*,"        ntracer is:"      , ntracer
       print*,"        nspecies is:"     , nspecies
       print*,"                                                        "
       print*,"--------------------------------------------------------"
    end if

    fileunit = 1
    if(mpi_rank().eq.0) then

       inquire(file='GRIST_runtime_nml_config.txt', exist=file_exists)
       if(file_exists) then
          print *, 'File exists. Deleting...'
          call delete_file_sub('GRIST_runtime_nml_config.txt')
       end if

       open (fileunit, status='new',file='GRIST_runtime_nml_config.txt')
       write(fileunit, nml= ctl_para)
       write(fileunit, nml= data_para)
       write(fileunit, nml= swe_para)
       write(fileunit, nml= dycore_para)
       write(fileunit, nml= tracer_para)
       write(fileunit, nml= physics_para)
       write(fileunit, nml= mesh_para)
       write(fileunit, nml= ccvt_para)
       close(fileunit)

    end if
   return

  end subroutine nml_set_global_vars

  end module grist_nml_module
