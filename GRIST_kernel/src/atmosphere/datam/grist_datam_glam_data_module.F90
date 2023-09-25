!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Data and method module for GRIST_LAM (GLAM)
!              All added data and method related to LAM should be defined here,
!              insead of modifying data, method and workflow outside this
!              module, to the highest extent
!
!              Name Convention:
!              glam_data_${varname}_at_$(location)
!
!        (1) LAM DATA have the same format as initial data, and this module
!            is also rewritten based on initilization method in datam_initial_data_module
!
!        (2) This module also provides methods for transforming the raw bdy
!            data to GRIST-applicable prognostic variables, following dtp_initial
!
!            currently we support: (i)   General pressure-based data (do not allow missing values)
!                                  (ii)  Hybri-(moist) pressure coordinate like EC-IFS/ERA, NCEP(FV3)'s hyai * p0 + hybi * ps
!                                  (iii) GRIST native
!
!        (3) The logic is: each time when lam bdy is needed, the data are exported
!            from external source (can be ERAIP, GFS, ERAIM, GRIST), 
!            then converted for GRIST importing; then a linear interp in time
!            is done (as that done for SST/SIC) to obtain instantaneous state for GRIST;
!            next, these states are used, together with the current model states, to derive
!            relaxation forcing tendency, and this forcing tendency is finally added to
!            model dynamics (dycore and tracer) as the bdy forcing at the relax region
!            for bdy data, we assume they are in hydrostatic balance, thus, in this module
!            (full) air pressure = mpressure in hydrostatic balance
!----------------------------------------------------------------------------

  module grist_datam_glam_data_module
#ifdef LAM_DOMAIN
! constant
    use grist_constants,                 only: rearth, i4, i8, r8, pi, one, zero, half, ptfactor, p00,rvap,rdry,cp, gravity
    use grist_nml_module,                only: model_timestep,glamAtmFilePath, glamDataSorc, ntracer_lamData, nlev_lamData, start_ymd, start_tod, working_mode,&
                                               nmif, mif_index, glamDataInterval, ntracer, glamLbcDiffusion, nh_dynamics, glamLbcCoef1,glamLbcCoef2,glamForcePhi
    use grist_hpe_constants,             only: eta_face_a, eta_face_b, eta_full_a, eta_full_b, p0
! data
    use grist_domain_types,              only: global_domain
    use grist_data_types,                only: scalar_1d_field, scalar_2d_field, scalar_3d_field, exchange_field_list_2d, exchange_field_list_1d
    use grist_mpi
    use grist_datam_static_data_module,  only: staticData_phis_at_pc_surface
    use grist_dycore_vars_module,        only: dycoreVarCellFull, dycoreVarEdgeFull, dycoreVarCellFace, dycoreVarSurface
    use grist_tracer_transport_vars_module,only: tracerVarCellFull
! method
    use grist_math_module,               only: lininterp
    use grist_math_module,               only: convert_vector_sph2cart
    use grist_time_manager,              only: get_curr_date
    use grist_hpe_hydro_pgf,             only: calc_hpe_hydro, calc_hpe_get_full_mass
    use grist_config_partition,          only: exchange_data_2d_add, exchange_data_2d, exchange_data_1d_add, exchange_data_1d
    use grist_fileio_list_1d_module_par, only: wrap_read_1d_group_rst
    use grist_fileio_list_2d_module_par, only: wrap_read_2d_group_rst
    use grist_fileio_list_3d_module_par, only: wrap_read_3d_group_rst

    implicit none
    private

    public :: grist_glam_data_construct, &
              grist_glam_data_destruct,  &
              grist_glam_data_read,      &
              grist_glam_data_time_interp,&
              grist_glam_data_idealized_dcmip2016_tc,&
              grist_glam_data_idealized_jwss,        &
              update_assignValueArea_1d,             &
              update_assignValueArea_2d,             &
              update_assignValueArea_3d,             &
              grist_glam_add_dynamics_bdy_forcing_term
              !grist_glam_add_dycore_bdy_forcing_term, &
              !grist_glam_add_tracer_bdy_forcing_term

!================================================
! primal cell, full level
!================================================

    type glamDataExternal     ! data exported from external source: ERA_P, GFS_P, GRIST...
       type(scalar_2d_field)     ::  uuu_at_pc_full_level
       type(scalar_2d_field)     ::  vvv_at_pc_full_level
       type(scalar_2d_field)     ::  ttt_at_pc_full_level
       type(scalar_3d_field)     ::  qqq_at_pc_full_level    ! in most reanalysis, this is specific humidity (moist mixing ratio), and we treat it so..
                                                             ! in GRIST, qqq is dry mixing ratio, no need to convert
       type(scalar_1d_field)     ::  ps_at_pc_surface        ! if from GRIST, this var=>hps
    end type glamDataExternal

    type(scalar_1d_field) :: glamData_plev
    type(scalar_1d_field) :: glamData_hyai
    type(scalar_1d_field) :: glamData_hybi
    type(scalar_1d_field) :: glamData_hyam
    type(scalar_1d_field) :: glamData_hybm

! this type contain those compatible with model dynamics
    type glamDataDynamics
! for glamData, we assume data is on hydrostatic balance, so mpressure is also full air pressure
        real(r8), allocatable    :: scalar_hpressure_at_pc_face_level_n(:,:)
        real(r8), allocatable    :: scalar_hpressure_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_delhp_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_mif_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_mpressure_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_mpressure_at_pc_face_level_n(:,:)

        real(r8), allocatable    :: scalar_delp_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_geopotential_at_pc_full_level_n(:,:) ! diagnostic

        real(r8), allocatable    :: scalar_U_wind_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_V_wind_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_temp_at_pc_full_level_n(:,:)
        real(r8), allocatable    :: scalar_tracer_mxrt_at_pc_full_level_n(:,:,:)! only has water vapor dimension for ntracer
    end type glamDataDynamics

    type glamDataPrognose    ! prognostic var used outside for adding bdy forcing
        type(scalar_2d_field)    :: scalar_normal_velocity_at_edge_full_level_n ! prognostic; will be used by exchange function
        real(r8), allocatable    :: scalar_mass_pt_at_pc_full_level_n(:,:)      ! prognostic
        real(r8), allocatable    :: scalar_tracer_mxrt_at_pc_full_level_n(:,:,:)! prognostic
        real(r8), allocatable    :: scalar_hpressure_at_pc_surface_n(:)         ! prognostic
        real(r8), allocatable    :: scalar_geopotential_at_pc_face_level_n(:,:) ! prognostic (ndc only)
    end type glamDataPrognose

    type(glamDataExternal)       :: t1_glamData, t2_glamData
    type(glamDataPrognose), public       :: t1_glamProg, t2_glamProg, t_glamProg
    type(glamDataDynamics)       :: glamData_dynamics

! common time used by interpolation function
    integer(i4)           :: t1_year,t1_mon,t1_day,t1_sec
    integer(i4)           :: t2_year,t2_mon,t2_day,t2_sec
    character(len=4)      :: char_t1_year, char_t2_year
    character(len=2)      :: char_t1_mon,  char_t2_mon
    character(len=2)      :: char_t1_day,  char_t2_day
    character(len=5)      :: char_t1_sec,  char_t2_sec

  CONTAINS

   subroutine grist_glam_data_construct(mesh,nlev,nlevp)
! data will be constructed at halo2, for bdy region, the two specific region layers are represented by halo1&2
! io
    type(global_domain), intent(in), target :: mesh
    integer(i4) , intent(in) :: nlev, nlevp
! local
    integer(i4)  :: it, ie, iv

if(.not.allocated(t1_glamProg%scalar_hpressure_at_pc_surface_n))       allocate(t1_glamProg%scalar_hpressure_at_pc_surface_n(mesh%nv_full))
if(.not.allocated(t1_glamProg%scalar_geopotential_at_pc_face_level_n)) allocate(t1_glamProg%scalar_geopotential_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n))  allocate(t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(ntracer,nlev,mesh%nv_full))
if(.not.allocated(t1_glamProg%scalar_mass_pt_at_pc_full_level_n))      allocate(t1_glamProg%scalar_mass_pt_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(t1_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)) allocate(t1_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(nlev,mesh%ne_full))

if(.not.allocated(t2_glamProg%scalar_hpressure_at_pc_surface_n))       allocate(t2_glamProg%scalar_hpressure_at_pc_surface_n(mesh%nv_full))
if(.not.allocated(t2_glamProg%scalar_geopotential_at_pc_face_level_n)) allocate(t2_glamProg%scalar_geopotential_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(t2_glamProg%scalar_tracer_mxrt_at_pc_full_level_n))  allocate(t2_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(ntracer,nlev,mesh%nv_full))
if(.not.allocated(t2_glamProg%scalar_mass_pt_at_pc_full_level_n))      allocate(t2_glamProg%scalar_mass_pt_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(t2_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)) allocate(t2_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(nlev,mesh%ne_full))

if(.not.allocated(t_glamProg%scalar_hpressure_at_pc_surface_n))       allocate(t_glamProg%scalar_hpressure_at_pc_surface_n(mesh%nv_full))
if(.not.allocated(t_glamProg%scalar_geopotential_at_pc_face_level_n)) allocate(t_glamProg%scalar_geopotential_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n))  allocate(t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(ntracer,nlev,mesh%nv_full))
if(.not.allocated(t_glamProg%scalar_mass_pt_at_pc_full_level_n))      allocate(t_glamProg%scalar_mass_pt_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)) allocate(t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(nlev,mesh%ne_full))

    if(mpi_rank().eq.0) print*,"GRIST_LAM data constructed sucessfully in grist_glam_data_construct"
    return
   end subroutine grist_glam_data_construct

   subroutine grist_glam_data_destruct

if(allocated(t1_glamProg%scalar_hpressure_at_pc_surface_n))       deallocate(t1_glamProg%scalar_hpressure_at_pc_surface_n)
if(allocated(t1_glamProg%scalar_geopotential_at_pc_face_level_n)) deallocate(t1_glamProg%scalar_geopotential_at_pc_face_level_n)
if(allocated(t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n))  deallocate(t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n)
if(allocated(t1_glamProg%scalar_mass_pt_at_pc_full_level_n))      deallocate(t1_glamProg%scalar_mass_pt_at_pc_full_level_n)
if(allocated(t1_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)) deallocate(t1_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)

if(allocated(t2_glamProg%scalar_hpressure_at_pc_surface_n))       deallocate(t2_glamProg%scalar_hpressure_at_pc_surface_n)
if(allocated(t2_glamProg%scalar_geopotential_at_pc_face_level_n)) deallocate(t2_glamProg%scalar_geopotential_at_pc_face_level_n)
if(allocated(t2_glamProg%scalar_tracer_mxrt_at_pc_full_level_n))  deallocate(t2_glamProg%scalar_tracer_mxrt_at_pc_full_level_n)
if(allocated(t2_glamProg%scalar_mass_pt_at_pc_full_level_n))      deallocate(t2_glamProg%scalar_mass_pt_at_pc_full_level_n)
if(allocated(t2_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)) deallocate(t2_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)

if(allocated(t_glamProg%scalar_hpressure_at_pc_surface_n))       deallocate(t_glamProg%scalar_hpressure_at_pc_surface_n)
if(allocated(t_glamProg%scalar_geopotential_at_pc_face_level_n)) deallocate(t_glamProg%scalar_geopotential_at_pc_face_level_n)
if(allocated(t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n))  deallocate(t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n)
if(allocated(t_glamProg%scalar_mass_pt_at_pc_full_level_n))      deallocate(t_glamProg%scalar_mass_pt_at_pc_full_level_n)
if(allocated(t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)) deallocate(t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)

      return
   end subroutine grist_glam_data_destruct

!================================================
! glam data read is done every xx-hour
! this time frequency is controlled here or within
! gcm_control? I think should be within gcm_control
!================================================

   subroutine grist_glam_data_read(mesh, dtime, itimestep,nlev,nlevp)
! io
     type(global_domain), intent(in), target :: mesh
     real(r8)     :: dtime
     integer(i4)  :: itimestep, nlev, nlevp
! local
     integer(i8)  :: dim2_len, dim3_len
     character(len=10)  :: level_type

     call get_curr_date(start_ymd, start_tod, itimestep , dtime, t1_year, t1_mon, t1_day, t1_sec)
     call handle_time
     if(mpi_rank().eq.0) print*,"t1_time:", char_t1_year, char_t1_mon, char_t1_day, char_t1_sec
     if(mpi_rank().eq.0) print*,"t2_time:", char_t2_year, char_t2_mon, char_t2_day, char_t2_sec 

     dim2_len = int(nlev_lamData,i8)
     dim3_len = int(ntracer_lamData,i8)

     if(.not.allocated(t1_glamData%uuu_at_pc_full_level%f))then
        allocate(t1_glamData%uuu_at_pc_full_level%f(nlev_lamData,mesh%nv_full)) ;t1_glamData%uuu_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t1_glamData%vvv_at_pc_full_level%f))then
        allocate(t1_glamData%vvv_at_pc_full_level%f(nlev_lamData,mesh%nv_full)) ;t1_glamData%vvv_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t1_glamData%ttt_at_pc_full_level%f))then
        allocate(t1_glamData%ttt_at_pc_full_level%f(nlev_lamData,mesh%nv_full)) ;t1_glamData%ttt_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t1_glamData%qqq_at_pc_full_level%f))then
        allocate(t1_glamData%qqq_at_pc_full_level%f(ntracer_lamData,nlev_lamData,mesh%nv_full),source=zero) ;t1_glamData%qqq_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t1_glamData%ps_at_pc_surface%f))then
        allocate(t1_glamData%ps_at_pc_surface%f(mesh%nv_full)) ;t1_glamData%ps_at_pc_surface%pos    = 0
     end if

     if(.not.allocated(t2_glamData%uuu_at_pc_full_level%f))then
        allocate(t2_glamData%uuu_at_pc_full_level%f(nlev_lamData,mesh%nv_full)) ;t2_glamData%uuu_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t2_glamData%vvv_at_pc_full_level%f))then
        allocate(t2_glamData%vvv_at_pc_full_level%f(nlev_lamData,mesh%nv_full)) ;t2_glamData%vvv_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t2_glamData%ttt_at_pc_full_level%f))then
        allocate(t2_glamData%ttt_at_pc_full_level%f(nlev_lamData,mesh%nv_full)) ;t2_glamData%ttt_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t2_glamData%qqq_at_pc_full_level%f))then
        allocate(t2_glamData%qqq_at_pc_full_level%f(ntracer_lamData,nlev_lamData,mesh%nv_full),source=zero) ;t2_glamData%qqq_at_pc_full_level%pos = 0
     end if
     if(.not.allocated(t2_glamData%ps_at_pc_surface%f))then
        allocate(t2_glamData%ps_at_pc_surface%f(mesh%nv_full))  ;t2_glamData%ps_at_pc_surface%pos    = 0
     end if

!================================================
! Time 1
!================================================

     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t1_year)//trim(char_t1_mon)//trim(char_t1_day)//trim(char_t1_sec)//'.nc',&
          'U', dim2_len, 0, t1_glamData%uuu_at_pc_full_level%f)
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t1_year)//trim(char_t1_mon)//trim(char_t1_day)//trim(char_t1_sec)//'.nc',&
          'V', dim2_len, 0, t1_glamData%vvv_at_pc_full_level%f)
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),& 
          'GRIST.lamData.'//trim(char_t1_year)//trim(char_t1_mon)//trim(char_t1_day)//trim(char_t1_sec)//'.nc',&
          'T', dim2_len, 0, t1_glamData%ttt_at_pc_full_level%f)


     if(trim(working_mode).ne.'dycore')&
     call wrap_read_3d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t1_year)//trim(char_t1_mon)//trim(char_t1_day)//trim(char_t1_sec)//'.nc',&
          'Q', dim2_len, 0, t1_glamData%qqq_at_pc_full_level%f,dim3_len)

!================================================
! Time 2
!================================================

     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
          'U', dim2_len, 0, t2_glamData%uuu_at_pc_full_level%f)
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
          'V', dim2_len, 0, t2_glamData%vvv_at_pc_full_level%f)
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),& 
          'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
          'T', dim2_len, 0, t2_glamData%ttt_at_pc_full_level%f)
     if(trim(working_mode).ne.'dycore')&
     call wrap_read_3d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
          'Q', dim2_len, 0, t2_glamData%qqq_at_pc_full_level%f,dim3_len)

     select case(trim(glamDataSorc))
!
! hybrid model level based on moisture pressure coordinate, e.g., merra(fv3), era(ifs)
!
     case('ERA_ml','MERRA_ml','ERA_Pres','GFS_Pres')

        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t1_year)//trim(char_t1_mon)//trim(char_t1_day)//trim(char_t1_sec)//'.nc',&
          'PS',          0, t1_glamData%ps_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
          'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
          'PS',          0, t2_glamData%ps_at_pc_surface%f)

        if(trim(glamDataSorc).eq.'ERA_ml'.or.trim(glamDataSorc).eq.'MERRA_ml')then

           level_type = 'mlev'
! nlev_lamData is for full-layer, and so nlev_lamData+1 is for face-level
           if(.not.allocated(glamData_hyai%f))          allocate(glamData_hyai%f(nlev_lamData+1))
           if(.not.allocated(glamData_hybi%f))          allocate(glamData_hybi%f(nlev_lamData+1))
           if(.not.allocated(glamData_hyam%f))          allocate(glamData_hyam%f(nlev_lamData))
           if(.not.allocated(glamData_hybm%f))          allocate(glamData_hybm%f(nlev_lamData))

           call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
                'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
                'hyai', 8, glamData_hyai%f, nlev_lamData+1)

           call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
                'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
                'hybi', 8, glamData_hybi%f, nlev_lamData+1)

           call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
                'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
                'hyam', 8, glamData_hyam%f, nlev_lamData)

           call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
                'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
                'hybm', 8, glamData_hybm%f, nlev_lamData)
        end if

        if(trim(glamDataSorc).eq.'ERA_Pres'.or.trim(glamDataSorc).eq.'GFS_Pres')then
           level_type = 'plev'

           if(.not.allocated(glamData_plev%f))          allocate(glamData_plev%f(nlev_lamData))

           call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
                'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
                'plev', 8, glamData_plev%f, nlev_lamData)
        end if


        call grist_glam_data_transform_pmlev(mesh, nlev, nlevp, &
                                          t1_glamData%uuu_at_pc_full_level%f,&
                                          t1_glamData%vvv_at_pc_full_level%f,&
                                          t1_glamData%ttt_at_pc_full_level%f,&
                                          t1_glamData%qqq_at_pc_full_level%f,&
                                          t1_glamData%ps_at_pc_surface%f    ,&
                                          t1_glamProg%scalar_normal_velocity_at_edge_full_level_n,& ! we need for LAM relaxation
                                          t1_glamProg%scalar_mass_pt_at_pc_full_level_n          ,& ! we need for LAM relaxation
                                          t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n      ,& ! we need for LAM relaxation
                                          t1_glamProg%scalar_hpressure_at_pc_surface_n           ,& ! we need for LAM relaxation
                                          t1_glamProg%scalar_geopotential_at_pc_face_level_n, trim(level_type) )      ! we need for LAM relaxation (optional)

        call grist_glam_data_transform_pmlev(mesh, nlev, nlevp, &
                                          t2_glamData%uuu_at_pc_full_level%f,&
                                          t2_glamData%vvv_at_pc_full_level%f,&
                                          t2_glamData%ttt_at_pc_full_level%f,&
                                          t2_glamData%qqq_at_pc_full_level%f,&
                                          t2_glamData%ps_at_pc_surface%f    ,&
                                          t2_glamProg%scalar_normal_velocity_at_edge_full_level_n,&
                                          t2_glamProg%scalar_mass_pt_at_pc_full_level_n          ,&
                                          t2_glamProg%scalar_tracer_mxrt_at_pc_full_level_n      ,&
                                          t2_glamProg%scalar_hpressure_at_pc_surface_n           ,&
                                          t2_glamProg%scalar_geopotential_at_pc_face_level_n, trim(level_type) )

        if(allocated(glamData_hyai%f))     deallocate(glamData_hyai%f)
        if(allocated(glamData_hybi%f))     deallocate(glamData_hybi%f)
        if(allocated(glamData_hyam%f))     deallocate(glamData_hyam%f)
        if(allocated(glamData_hybm%f))     deallocate(glamData_hybm%f)
        if(allocated(glamData_plev%f))     deallocate(glamData_plev%f)

     case('GRIST')
!
! GRIST's own model level
! hps to replace ps
!
        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
           'GRIST.lamData.'//trim(char_t1_year)//trim(char_t1_mon)//trim(char_t1_day)//trim(char_t1_sec)//'.nc',&
           'HPS',          0, t1_glamData%ps_at_pc_surface%f)

        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(glamAtmFilePath),&
           'GRIST.lamData.'//trim(char_t2_year)//trim(char_t2_mon)//trim(char_t2_day)//trim(char_t2_sec)//'.nc',&
           'HPS',          0, t2_glamData%ps_at_pc_surface%f)

        call grist_glam_data_transform_grist_mlev(mesh, nlev, nlevp, &
                                          t1_glamData%uuu_at_pc_full_level%f,&
                                          t1_glamData%vvv_at_pc_full_level%f,&
                                          t1_glamData%ttt_at_pc_full_level%f,&
                                          t1_glamData%qqq_at_pc_full_level%f,&
                                          t1_glamData%ps_at_pc_surface%f    ,&
                                          t1_glamProg%scalar_normal_velocity_at_edge_full_level_n, &
                                          t1_glamProg%scalar_mass_pt_at_pc_full_level_n          , &
                                          t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n      , &
                                          t1_glamProg%scalar_hpressure_at_pc_surface_n           , &
                                          t1_glamProg%scalar_geopotential_at_pc_face_level_n )

        call grist_glam_data_transform_grist_mlev(mesh, nlev, nlevp, &
                                          t2_glamData%uuu_at_pc_full_level%f,&
                                          t2_glamData%vvv_at_pc_full_level%f,&
                                          t2_glamData%ttt_at_pc_full_level%f,&
                                          t2_glamData%qqq_at_pc_full_level%f,&
                                          t2_glamData%ps_at_pc_surface%f    ,&
                                          t2_glamProg%scalar_normal_velocity_at_edge_full_level_n, &
                                          t2_glamProg%scalar_mass_pt_at_pc_full_level_n          , &
                                          t2_glamProg%scalar_tracer_mxrt_at_pc_full_level_n      , &
                                          t2_glamProg%scalar_hpressure_at_pc_surface_n           , &
                                          t2_glamProg%scalar_geopotential_at_pc_face_level_n)

     case default
         if (mpi_rank().eq.0) print*,"glamDataSorc is undefined or defined not right, current support: ERAIP, GFS, GRIST"
     end select
!
! clean data
!
     if(allocated(t1_glamData%uuu_at_pc_full_level%f)) deallocate(t1_glamData%uuu_at_pc_full_level%f)
     if(allocated(t1_glamData%vvv_at_pc_full_level%f)) deallocate(t1_glamData%vvv_at_pc_full_level%f)
     if(allocated(t1_glamData%ttt_at_pc_full_level%f)) deallocate(t1_glamData%ttt_at_pc_full_level%f)
     if(allocated(t1_glamData%qqq_at_pc_full_level%f)) deallocate(t1_glamData%qqq_at_pc_full_level%f)
     if(allocated(t1_glamData%ps_at_pc_surface%f))     deallocate(t1_glamData%ps_at_pc_surface%f)

     if(allocated(t2_glamData%uuu_at_pc_full_level%f)) deallocate(t2_glamData%uuu_at_pc_full_level%f)
     if(allocated(t2_glamData%vvv_at_pc_full_level%f)) deallocate(t2_glamData%vvv_at_pc_full_level%f)
     if(allocated(t2_glamData%ttt_at_pc_full_level%f)) deallocate(t2_glamData%ttt_at_pc_full_level%f)
     if(allocated(t2_glamData%qqq_at_pc_full_level%f)) deallocate(t2_glamData%qqq_at_pc_full_level%f)
     if(allocated(t2_glamData%ps_at_pc_surface%f))     deallocate(t2_glamData%ps_at_pc_surface%f)

     return
   end subroutine grist_glam_data_read

   subroutine grist_glam_data_time_interp(dtime, itimestep)
!
! the interpolated var is Prog var directly used for updating the relaxation tendency so that
! transform only done once after each data reading
!
!io
     real(r8)     :: dtime
     integer(i4)  :: itimestep
!local
     integer(i4)  :: t_year,t_mon,t_day,t_sec
     integer(i4)  :: tn_sec

! interpolte in time based on t1_data and t2_data
     call get_curr_date(start_ymd, start_tod, itimestep , dtime, t_year, t_mon, t_day, t_sec)
     if(mpi_rank().eq.0) print*,"current time is: ",t_year,t_mon,t_day,t_sec,t1_sec,t2_sec

! t1----t-----tn; we use linear interpolation in time for obtain t_glamProg based on t1_glamProg and t2_glamProg
     tn_sec = t2_sec
     if(t2_sec.eq.0)  tn_sec = 86400

     t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f = t1_glamProg%scalar_normal_velocity_at_edge_full_level_n%f +&
     real(t_sec-t1_sec)/real(tn_sec-t1_sec)*(t2_glamProg%scalar_normal_velocity_at_edge_full_level_n%f-t1_glamProg%scalar_normal_velocity_at_edge_full_level_n%f)

     t_glamProg%scalar_mass_pt_at_pc_full_level_n     = t1_glamProg%scalar_mass_pt_at_pc_full_level_n +&
     real(t_sec-t1_sec)/real(tn_sec-t1_sec)*(t2_glamProg%scalar_mass_pt_at_pc_full_level_n-t1_glamProg%scalar_mass_pt_at_pc_full_level_n)

     t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n = t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n +&
     real(t_sec-t1_sec)/real(tn_sec-t1_sec)*(t2_glamProg%scalar_tracer_mxrt_at_pc_full_level_n-t1_glamProg%scalar_tracer_mxrt_at_pc_full_level_n)

     t_glamProg%scalar_hpressure_at_pc_surface_n      = t1_glamProg%scalar_hpressure_at_pc_surface_n +&
     real(t_sec-t1_sec)/real(tn_sec-t1_sec)*(t2_glamProg%scalar_hpressure_at_pc_surface_n-t1_glamProg%scalar_hpressure_at_pc_surface_n)

     t_glamProg%scalar_geopotential_at_pc_face_level_n= t1_glamProg%scalar_geopotential_at_pc_face_level_n +&
     real(t_sec-t1_sec)/real(tn_sec-t1_sec)*(t2_glamProg%scalar_geopotential_at_pc_face_level_n-t1_glamProg%scalar_geopotential_at_pc_face_level_n)

     return
   end subroutine grist_glam_data_time_interp

   subroutine grist_glam_add_dynamics_bdy_forcing_term(mesh,nlev)
! io
    type(global_domain),  intent(in) :: mesh
    integer(i4) :: nlev
! local
    integer(i4) :: ie, iv, icell1, icell2, cell_mark, nlevp, ilev, itracer
    real(r8)    :: f1_coef, delhp(1:nlev)
    type(exchange_field_list_2d),pointer :: field_head_2d

    field_head_2d=>null()
  
    nlevp = nlev+1
!
! We regard the edges of the outmost relaxation-zone cell form the physical boundary of the LAM;
! thus, for spec region, only cell-point data matter (for obtaining edge-point normal-v at the outmost relaxation-zone edge;
! for relax region, we nudge normal velocity, mass-pt, mass-qv, towards their large-scale state;
!
    !if(mesh%hasBdyDomain)then

       do ie = 1, mesh%ne_full
          if(mesh%edt_mark(ie).eq.188)then
             f1_coef   = one/5_r8/(glamLbcCoef1*mesh%edt_leng(ie)*rearth)
#ifdef LAM_VCOEF
             icell1 = mesh%edt_v(1,ie)
             icell2 = mesh%edt_v(2,ie)
             cell_mark = min(mesh%plg_mark(icell1),mesh%plg_mark(icell2))
             if(cell_mark.eq.199)then
                cell_mark=0
             elseif(cell_mark.eq.177)then
                cell_mark=5
             else
                cell_mark=cell_mark-1880 ! from 1 to 5
             end if
             if(cell_mark.lt.0.or.cell_mark.gt.5)then
                print*,"cell_mark for edge is wrong",cell_mark
             end if
             f1_coef   = cell_mark*one/25_r8/(glamLbcCoef1*mesh%edt_leng(ie)*rearth)
#endif
             dycoreVarEdgeFull%scalar_normal_velocity_n%f(1:nlev,ie) = dycoreVarEdgeFull%scalar_normal_velocity_n%f(1:nlev,ie)+model_timestep*f1_coef*&
                                                                   (t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(1:nlev,ie)-&
                                                             dycoreVarEdgeFull%scalar_normal_velocity_n%f(1:nlev,ie))
          end if
       end do

       if(glamLbcDiffusion)then
          call add_laplacian_forcing_velocity(mesh, nlev, model_timestep, &
                                              t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f,&
                                              dycoreVarEdgeFull%scalar_normal_velocity_n)
          !call add_laplacian_forcing_scalar_1d(mesh, model_timestep, t_glamProg%scalar_hpressure_at_pc_surface_n, dycoreVarSurface%scalar_hpressure_n%f)
          !call add_laplacian_forcing_scalar_2d(mesh, nlev, model_timestep, t_glamProg%scalar_mass_pt_at_pc_full_level_n, dycoreVarCellFull%scalar_mass_pt_n%f)
          call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
          call exchange_data_2d(mesh%local_block,field_head_2d)
       end if

       do iv = 1, mesh%nv_full
          if(mesh%plg_mark(iv).eq.1881.or.mesh%plg_mark(iv).eq.1882.or.mesh%plg_mark(iv).eq.1883.or.mesh%plg_mark(iv).eq.1884.or.mesh%plg_mark(iv).eq.1885)then
             f1_coef   = one/5_r8/(glamLbcCoef1*mesh%vtxCellLeng(iv)*rearth)
#ifdef LAM_VCOEF
             f1_coef   = (mesh%plg_mark(iv)-1880)/25_r8/(glamLbcCoef1*mesh%vtxCellLeng(iv)*rearth)
#endif

             dycoreVarCellFull%scalar_mass_pt_n%f(1:nlev,iv) = dycoreVarCellFull%scalar_mass_pt_n%f(1:nlev,iv)+model_timestep*f1_coef*&
                                                       (t_glamProg%scalar_mass_pt_at_pc_full_level_n(1:nlev,iv)-&
                                                        dycoreVarCellFull%scalar_mass_pt_n%f(1:nlev,iv))

             dycoreVarSurface%scalar_hpressure_n%f(iv) = dycoreVarSurface%scalar_hpressure_n%f(iv)+model_timestep*f1_coef*&
                                                              (t_glamProg%scalar_hpressure_at_pc_surface_n(iv)-dycoreVarSurface%scalar_hpressure_n%f(iv))
#ifdef DYCORE_NHD
             if(nh_dynamics.and.glamForcePhi)then
                dycoreVarCellFace%scalar_phi_n%f(1:nlev+1,iv) = dycoreVarCellFace%scalar_phi_n%f(1:nlev+1,iv)+model_timestep*f1_coef*&
                                                              (t_glamProg%scalar_geopotential_at_pc_face_level_n(1:nlev+1,iv)-&
                                                               dycoreVarCellFace%scalar_phi_n%f(1:nlev+1,iv))
             end if
#endif
          end if
       end do

       if(trim(working_mode).ne.'dycore')then
          do iv = 1, mesh%nv_full
             delhp(1:nlev) = (eta_face_a(2:nlevp)-eta_face_a(1:nlev))*p0+(eta_face_b(2:nlevp)-eta_face_b(1:nlev))*t_glamProg%scalar_hpressure_at_pc_surface_n(iv)
             if(mesh%plg_mark(iv).eq.1881.or.mesh%plg_mark(iv).eq.1882.or.mesh%plg_mark(iv).eq.1883.or.mesh%plg_mark(iv).eq.1884.or.mesh%plg_mark(iv).eq.1885)then
                f1_coef   = one/5_r8/(glamLbcCoef1*mesh%vtxCellLeng(iv)*rearth)
#ifdef LAM_VCOEF
                f1_coef   = (mesh%plg_mark(iv)-1880)/25_r8/(glamLbcCoef1*mesh%vtxCellLeng(iv)*rearth)
#endif
                do itracer = 1, ntracer
                   tracerVarCellFull%scalar_tracer_mass_n%f(itracer,1:nlev,iv) = tracerVarCellFull%scalar_tracer_mass_n%f(itracer,1:nlev,iv)+model_timestep*f1_coef*&
                                                              (delhp(1:nlev)*t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(itracer,1:nlev,iv)-&
                                                        tracerVarCellFull%scalar_tracer_mass_n%f(itracer,1:nlev,iv))

                   tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,1:nlev,iv) = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,1:nlev,iv)+model_timestep*f1_coef*&
                                                (t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(itracer,1:nlev,iv)-&
                                          tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,1:nlev,iv))
                end do
             end if
          end do
       end if

    return
   end subroutine grist_glam_add_dynamics_bdy_forcing_term

   subroutine update_assignValueArea_1d(mesh,varName,varValue)
!io
    type(global_domain),  intent(in)      :: mesh
    character(len=*), intent(in)          :: varName
    type(scalar_1d_field), intent(inout)  :: varValue
!local
    integer(i4)                           :: iv

      if(mesh%hasBdyDomain)then
         select case(trim(varName))
         case('hpressure_surface')
           do iv = mesh%nv_compute+1, mesh%nv_full
              if(mesh%plg_mark(iv).eq.177)then
                 varValue%f(iv) = t_glamProg%scalar_hpressure_at_pc_surface_n(iv)
              end if
           end do
         case default
           if(mpi_rank().eq.0) print*, "update_assignValueArea_1d bad behavior"
         end select
      end if

      return
   end subroutine update_assignValueArea_1d

   subroutine update_assignValueArea_2d(mesh, nlev, varName,varValue)
!io
    type(global_domain),  intent(in)     :: mesh
    integer(i4),          intent(in)     :: nlev
    character(len=*),     intent(in)     :: varName
    type(scalar_2d_field),intent(inout)  :: varValue
!local
    integer(i4)                           :: iv, ilev, ie, icell1, icell2

      if(mesh%hasBdyDomain)then

         select case(trim(varName))
         case('mass_pt')

           do iv = mesh%nv_compute+1, mesh%nv_full
              if(mesh%plg_mark(iv).eq.177)then
                 varValue%f(1:nlev,iv) = t_glamProg%scalar_mass_pt_at_pc_full_level_n(1:nlev,iv)
              end if
           end do

         case('normal_velocity')
           do ie = mesh%ne_compute+1, mesh%ne_full
              icell1 = mesh%edt_v(1,ie)
              icell2 = mesh%edt_v(2,ie)
              if(mesh%edt_mark(ie).eq.1771.or.mesh%edt_mark(ie).eq.1772)then
                 varValue%f(1:nlev,ie) = t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(1:nlev,ie)
              end if
           end do

         case('wwwFace','tend_et_www')
! not actually used for evolving www, only to provide zero dw/dt for a part of hpgf at the domain
! boundary
           do iv = mesh%nv_compute+1, mesh%nv_full
              if(mesh%plg_mark(iv).eq.177)then
                 varValue%f(1:nlev,iv) = zero
              end if
           end do
         case('phiFace')

           do iv = mesh%nv_compute+1, mesh%nv_full
              if(mesh%plg_mark(iv).eq.177)then
                 varValue%f(1:nlev,iv) = t_glamProg%scalar_geopotential_at_pc_face_level_n(1:nlev,iv)
              end if
           end do

         case('phiFull')
           do iv = mesh%nv_compute+1, mesh%nv_full
              if(mesh%plg_mark(iv).eq.177)then
                 varValue%f(1:nlev,iv) = 0.5_r8*(t_glamProg%scalar_geopotential_at_pc_face_level_n(1:nlev,iv)+&
                                                 t_glamProg%scalar_geopotential_at_pc_face_level_n(2:nlev+1,iv))
              end if
           end do
        case default
           if(mpi_rank().eq.0) print*, "update_assignValueArea_2d bad behavior"
        end select
      end if

      return
   end subroutine update_assignValueArea_2d

   subroutine update_assignValueArea_3d(mesh,nlev,ntracer,varName,delhp_avg,tracer_mass,tracer_mxrt)
!io
    type(global_domain),  intent(in)     :: mesh
    integer(i4),          intent(in)     :: nlev, ntracer
    character(len=*),     intent(in)     :: varName
    real(r8),allocatable ,intent(in)     :: delhp_avg(:,:)
    type(scalar_3d_field),intent(inout)  :: tracer_mass,tracer_mxrt
!local
    integer(i4)                          :: iv, ilev, itracer
    real(r8)                             :: delhp(nlev)

      if(mesh%hasBdyDomain)then
         do iv = mesh%nv_compute+1, mesh%nv_full
            if(mesh%plg_mark(iv).eq.177)then
               delhp(1:nlev) = (eta_face_a(2:nlev+1)-eta_face_a(1:nlev))*p0+(eta_face_b(2:nlev+1)-eta_face_b(1:nlev))*t_glamProg%scalar_hpressure_at_pc_surface_n(iv)
               do ilev = 1, nlev
                  do itracer = 1, ntracer
                     tracer_mxrt%f(itracer,ilev,iv) = t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(itracer,ilev,iv)
                     tracer_mass%f(itracer,ilev,iv) = delhp(ilev)*t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(itracer,ilev,iv)
                  end do
               end do
            end if
         end do
      end if

      return
   end subroutine update_assignValueArea_3d

!================================================
! Private
!================================================

   subroutine grist_glam_data_idealized_dcmip2016_tc(mesh,nlev,nlevp)

!================================================
! we use this routine for deriving bdy conditions
! for dcmip2016, used directly to force model
!================================================
! this routine will repalce grist_glam_data_read+grist_glam_data_transform
! to derive constant t_glamProg does not change over time
!
    use grist_dtp_dcmip2016_tc,      only: tropical_cyclone_test, dp, rp, exppr, cen_lat, cen_lon,&
                                           ztrop,Ttrop,ptrop, Rd, T0, q0, exponent, gamma
    use grist_dtp_initial_module,    only: evaluate_dry_mass_height

!io
    type(global_domain),  intent(inout) :: mesh
    integer(i4)        ,  intent(in)    :: nlev, nlevp
!local
    real(r8)                            :: utmp, vtmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
    integer(i4)                         :: iv, ie, ilev, icell1, icell2
    real(r8)                            :: vector_velocity(3)
    type(exchange_field_list_2d),pointer:: field_head_2d

if(.not.allocated(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)) allocate(glamData_dynamics%scalar_hpressure_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)) allocate(glamData_dynamics%scalar_hpressure_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)) allocate(glamData_dynamics%scalar_mpressure_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)) allocate(glamData_dynamics%scalar_mpressure_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_geopotential_at_pc_full_level_n)) allocate(glamData_dynamics%scalar_geopotential_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_U_wind_at_pc_full_level_n))    allocate(glamData_dynamics%scalar_U_wind_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_V_wind_at_pc_full_level_n))    allocate(glamData_dynamics%scalar_V_wind_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_temp_at_pc_full_level_n))      allocate(glamData_dynamics%scalar_temp_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n))allocate(glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(ntracer,nlev,mesh%nv_full))

!
! evaluate dry mass and geometric height at each level
!
      if(.not.allocated(mesh%vtx)) allocate(mesh%vtx(mesh%nv_full))
      do iv = 1,mesh%nv_full
         mesh%vtx(iv)%lon = mesh%vtx_lon(iv)
         mesh%vtx(iv)%lat = mesh%vtx_lat(iv)
      end do
      call evaluate_dry_mass_height(mesh, mesh%nv_full, nlev, 'DCMIP2016-TC' , &
                                           t_glamProg%scalar_hpressure_at_pc_surface_n      , &
                                    glamData_dynamics%scalar_hpressure_at_pc_face_level_n   , &
                                    glamData_dynamics%scalar_hpressure_at_pc_full_level_n   , &
                                           t_glamProg%scalar_geopotential_at_pc_face_level_n, &
                                    glamData_dynamics%scalar_geopotential_at_pc_full_level_n)
      if(allocated(mesh%vtx)) deallocate(mesh%vtx)
!
! final evaluation based on geop full
!
      do iv = 1, mesh%nv_full
        do ilev = 1, nlev
           call tropical_cyclone_test(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8)  ,&
                                      glamData_dynamics%scalar_mpressure_at_pc_full_level_n(ilev,iv)       ,&
                                      glamData_dynamics%scalar_geopotential_at_pc_full_level_n(ilev,iv), 1 ,&
                                      glamData_dynamics%scalar_U_wind_at_pc_full_level_n(ilev,iv)          ,&
                                      glamData_dynamics%scalar_V_wind_at_pc_full_level_n(ilev,iv)          ,&
                                      glamData_dynamics%scalar_temp_at_pc_full_level_n(ilev,iv)            ,&
                                      tmp1                                             ,&  ! not used
                                      t_glamProg%scalar_geopotential_at_pc_face_level_n(nlevp,iv)          ,&
                                      glamData_dynamics%scalar_mpressure_at_pc_face_level_n(nlevp,iv)      ,&  ! not used
                                      tmp2      ,&  ! not used
                                      glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1,ilev,iv))
! moist to dry
        t_glamProg%scalar_tracer_mxrt_at_pc_full_level_n(1,ilev,iv) = glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1,ilev,iv)/&
                                                                 (one-glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1,ilev,iv))
 
        end do
      end do

      t_glamProg%scalar_geopotential_at_pc_face_level_n        = gravity*t_glamProg%scalar_geopotential_at_pc_face_level_n
      glamData_dynamics%scalar_geopotential_at_pc_full_level_n = gravity*glamData_dynamics%scalar_geopotential_at_pc_full_level_n

      do iv = 1, mesh%nv_full
         do ilev = 1, nlev
            t_glamProg%scalar_mass_pt_at_pc_full_level_n(ilev,iv)  = (glamData_dynamics%scalar_hpressure_at_pc_face_level_n(ilev+1,iv)-&
                                                                      glamData_dynamics%scalar_hpressure_at_pc_face_level_n(ilev,iv))*&
                                                          glamData_dynamics%scalar_temp_at_pc_full_level_n(ilev,iv)/&
                                                  ( (glamData_dynamics%scalar_mpressure_at_pc_full_level_n(ilev,iv)/p00)**(rdry/cp))*&
                                                  (one+ptfactor*glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1,ilev,iv)/&
                                                           (one-glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1,ilev,iv)))
         end do
      end do

      field_head_2d => null()
! until halo1
      do ie = 1, mesh%ne_halo(1)
         icell1 = mesh%edt_v(1,ie)
         icell2 = mesh%edt_v(2,ie)
         do ilev = 1, nlev
            utmp = half*(glamData_dynamics%scalar_U_wind_at_pc_full_level_n(ilev,icell1)+glamData_dynamics%scalar_U_wind_at_pc_full_level_n(ilev,icell2))
            vtmp = half*(glamData_dynamics%scalar_V_wind_at_pc_full_level_n(ilev,icell1)+glamData_dynamics%scalar_V_wind_at_pc_full_level_n(ilev,icell2))
            call convert_vector_sph2cart(utmp, vtmp, real(mesh%edt_c_p(1:3,ie),r8), vector_velocity)
            t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(ilev,ie) = dot_product(vector_velocity, real(mesh%edp_nr(1:3,ie),r8))
         end do
      end do
! exchange data
      call exchange_data_2d_add(mesh,field_head_2d,t_glamProg%scalar_normal_velocity_at_edge_full_level_n)
      call exchange_data_2d(mesh%local_block,field_head_2d)

      if(mesh%hasBdyDomain)then
         do ie = mesh%ne_halo(1)+1,mesh%ne_halo(2)
            t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(1:nlev,ie) = zero
         end do
      end if 

if(allocated(glamData_dynamics%scalar_hpressure_at_pc_face_level_n))    deallocate(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)
if(allocated(glamData_dynamics%scalar_hpressure_at_pc_full_level_n))    deallocate(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)
if(allocated(glamData_dynamics%scalar_mpressure_at_pc_face_level_n))    deallocate(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)
if(allocated(glamData_dynamics%scalar_mpressure_at_pc_full_level_n))    deallocate(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)
if(allocated(glamData_dynamics%scalar_geopotential_at_pc_full_level_n)) deallocate(glamData_dynamics%scalar_geopotential_at_pc_full_level_n)
if(allocated(glamData_dynamics%scalar_U_wind_at_pc_full_level_n))       deallocate(glamData_dynamics%scalar_U_wind_at_pc_full_level_n)
if(allocated(glamData_dynamics%scalar_V_wind_at_pc_full_level_n))       deallocate(glamData_dynamics%scalar_V_wind_at_pc_full_level_n)
if(allocated(glamData_dynamics%scalar_temp_at_pc_full_level_n))         deallocate(glamData_dynamics%scalar_temp_at_pc_full_level_n)
if(allocated(glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n))  deallocate(glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n)

     return
   end subroutine grist_glam_data_idealized_dcmip2016_tc

   subroutine grist_glam_data_idealized_jwss(mesh,nlev,nlevp, itimestep)
!================================================
! we use this routine for deriving bdy conditions
! for dcmip2016, used directly to force model
!================================================
! this routine will repalce grist_glam_data_read+grist_glam_data_transform
! to derive constant t_glamProg does not change over time
!
    use grist_dtp_dcmip2016_tc,      only: tropical_cyclone_test, dp, rp, exppr, cen_lat, cen_lon,&
                                           ztrop,Ttrop,ptrop, Rd, T0, q0, exponent, gamma
    use grist_dtp_initial_module,    only: evaluate_dry_mass_height

!io
    type(global_domain),  intent(inout) :: mesh
    integer(i4)        ,  intent(in)    :: nlev, nlevp, itimestep
!local
    real(r8)                            :: utmp, vtmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
    integer(i4)                         :: iv, ie, ilev, icell1, icell2
    real(r8)                            :: vector_velocity(3)
    type(exchange_field_list_2d),pointer:: field_head_2d
     
    if(itimestep.le.1)then
       t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(1:nlev,:) = dycoreVarEdgeFull%scalar_normal_velocity_n%f(1:nlev,:)
       t_glamProg%scalar_mass_pt_at_pc_full_level_n(1:nlev,:)             = dycoreVarCellFull%scalar_mass_pt_n%f(1:nlev,:)
       t_glamProg%scalar_hpressure_at_pc_surface_n(:)                     = dycoreVarSurface%scalar_hpressure_n%f(:)
    end if
#ifdef LAM_DOMAIN
    if(mesh%hasBdyDomain)then
       !do ie = mesh%ne_halo(1)+1,mesh%ne_halo(2)
       !   t_glamProg%scalar_normal_velocity_at_edge_full_level_n%f(1:nlev,ie) = zero
       !end do
    end if
#endif
    return
   end subroutine grist_glam_data_idealized_jwss

!
! when using this routine, we assume to use pressure-level data at a fixed number of vertical levels
! note that qqq is defined as specific humidity as typically provided by most 
! analysis/reanalysis/climate-model datasets
! 2024, Oct. generalized for era/merra model level data
!

   subroutine grist_glam_data_transform_pmlev(mesh, nlev, nlevp, &
                                              glamData_uuu_at_pc_full_level,&
                                              glamData_vvv_at_pc_full_level,&
                                              glamData_ttt_at_pc_full_level,&
                                              glamData_qqq_at_pc_full_level,&
                                              glamData_ps_at_pc_surface    ,&
                                              scalar_normal_velocity_at_edge_full_level_n, &
                                              scalar_mass_pt_at_pc_full_level_n          , &
                                              scalar_tracer_mxrt_at_pc_full_level_n      , &
                                              scalar_hpressure_at_pc_surface_n           , &
                                              scalar_geopotential_at_pc_face_level_n     , &
                                              level_type)
! io
   type(global_domain),   intent(in)    :: mesh
   integer(i4)        ,   intent(in)    :: nlev, nlevp
   real(r8), allocatable, intent(in)    :: glamData_uuu_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: glamData_vvv_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: glamData_ttt_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: glamData_qqq_at_pc_full_level(:,:,:)
   real(r8), allocatable, intent(in)    :: glamData_ps_at_pc_surface(:)
   character(len=*),      intent(in)    :: level_type
   !inout
   type(scalar_2d_field), intent(inout) :: SCALAR_NORMAL_VELOCITY_AT_EDGE_FULL_LEVEL_N
   real(r8), allocatable, intent(inout) :: SCALAR_MASS_PT_AT_PC_FULL_LEVEL_N(:,:)
   real(r8), allocatable, intent(inout) :: SCALAR_TRACER_MXRT_AT_PC_FULL_LEVEL_N(:,:,:)
   real(r8), allocatable, intent(inout) :: SCALAR_HPRESSURE_AT_PC_SURFACE_N(:)
   real(r8), allocatable, intent(inout) :: SCALAR_GEOPOTENTIAL_AT_PC_FACE_LEVEL_N(:,:)
!local
   real(r8)                             :: pface(nlev_lamdata+1), hpface(nlev_lamdata+1)
   real(r8)                             :: pfull(nlev_lamdata),   hpfull(nlev_lamdata)
   real(r8)                             :: dpfull(nlev_lamdata),  dhpfull(nlev_lamdata)
   real(r8)                             :: utmp, vtmp
   real(r8)                             :: vector_velocity(3)
! for real-case
   real(r8)                             :: scalar_template_1d_nlev_a(nlev)
   real(r8)                             :: scalar_template_1d_nlev_b(nlev)
   real(r8)                             :: scalar_template_1d_nlev_c(nlev)
   real(r8)                             :: scalar_template_1d_nlevp_a(nlevp)
   real(r8)                             :: scalar_template_1d_nlevp_b(nlevp)
   real(r8)                             :: scalar_template_a
   integer(i4)                          :: nlev_eff
   integer(i4)                          :: iv, ie, icell1, icell2, ilev
   type(exchange_field_list_2d),pointer :: field_head_2d

!==============================================================
! This entry only takes care of ERA/GFS pressure-based data
!==============================================================

if(.not.allocated(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)) allocate(glamData_dynamics%scalar_hpressure_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)) allocate(glamData_dynamics%scalar_hpressure_at_pc_full_level_n(nlev ,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_delhp_at_pc_full_level_n))     allocate(glamData_dynamics%scalar_delhp_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mif_at_pc_full_level_n))       allocate(glamData_dynamics%scalar_mif_at_pc_full_level_n(nlev  ,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)) allocate(glamData_dynamics%scalar_mpressure_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)) allocate(glamData_dynamics%scalar_mpressure_at_pc_full_level_n(nlev ,mesh%nv_full))

if(.not.allocated(glamData_dynamics%scalar_U_wind_at_pc_full_level_n))    allocate(glamData_dynamics%scalar_U_wind_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_V_wind_at_pc_full_level_n))    allocate(glamData_dynamics%scalar_V_wind_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_temp_at_pc_full_level_n))      allocate(glamData_dynamics%scalar_temp_at_pc_full_level_n(nlev,mesh%nv_full))

     Do iv = 1, mesh%nv_full
        
        if(trim(level_type).eq.'plev')then
! nominal pressure at ERA pressure levels
           pfull(1:nlev_lamData) = glamData_plev%f(1:nlev_lamData)  ! Pa
! check effective full level number
           nlev_eff = 99999
           do ilev = 1, nlev_lamdata
              if(pfull(ilev).ge.glamData_ps_at_pc_surface(iv))then
                 nlev_eff = ilev
                 exit
              end if
           end do
           if(nlev_eff.eq.99999) nlev_eff = nlev_lamdata

           pface(2:nlev_eff)     = (pfull(1:nlev_eff-1)+pfull(2:nlev_eff))*half
           pface(1)              = zero
           pface(nlev_eff+1)     = glamData_ps_at_pc_surface(iv)      ! if ps<pface(nlev_eff), let it be, good assumption! no problem till now...
        end if

        if(level_type.eq.'mlev')then
           nlev_eff              = nlev_lamData
! Reconstruct their moist (hydro) pressure at their coordinates
           pfull(1:nlev_eff)     = glamData_hyam%f(1:nlev_eff)*1e5   + glamData_hybm%f(1:nlev_eff)  *glamData_ps_at_pc_surface(iv)  ! Pa
           pface(1:nlev_eff+1)   = glamData_hyai%f(1:nlev_eff+1)*1e5 + glamData_hybi%f(1:nlev_eff+1)*glamData_ps_at_pc_surface(iv)  ! Pa
        end if

        dpfull(1:nlev_eff)    = pface(2:nlev_eff+1)-pface(1:nlev_eff)
        hpface(1)             = pface(1)
!
! evaluate dry mass at ERA/GFS pressure level and surface, assume only water vapor and dry air exist in the LAM atmosphere
! that is, q(c,i,...) = 0
!
        do ilev = 2, nlev_eff+1
           hpface(ilev) = hpface(ilev-1)+dpfull(ilev-1)*(one-glamData_qqq_at_pc_full_level(1,ilev-1,iv))
        end do
        hpfull(1:nlev_eff) = half*(hpface(1:nlev_eff)+hpface(2:nlev_eff+1))
!
! Milestone 1: surface dry air mass pressure in LamData
!
        SCALAR_HPRESSURE_AT_PC_SURFACE_N(iv)                               = hpface(nlev_eff+1)
        glamData_dynamics%scalar_hpressure_at_pc_face_level_n(1:nlevp,iv)  = eta_face_a(1:nlevp)*p0+eta_face_b(1:nlevp)*hpface(nlev_eff+1) ! Pa
        glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv)   = eta_full_a(1:nlev) *p0+eta_full_b(1:nlev) *hpface(nlev_eff+1)
        glamData_dynamics%scalar_delhp_at_pc_full_level_n(1:nlev,iv)       = glamData_dynamics%scalar_hpressure_at_pc_face_level_n(2:nlevp,iv)-&
                                                                             glamData_dynamics%scalar_hpressure_at_pc_face_level_n(1:nlev,iv)
!
! interpolate from ERA pressure level to GRIST model level based on dry air mass
!

        call lininterp(glamData_uuu_at_pc_full_level(1:nlev_eff,iv),   hpfull(1:nlev_eff), 1, nlev_eff, &
                       glamData_dynamics%scalar_U_wind_at_pc_full_level_n(1:nlev,iv), glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)
        call lininterp(glamData_vvv_at_pc_full_level(1:nlev_eff,iv),   hpfull(1:nlev_eff), 1, nlev_eff, &
                       glamData_dynamics%scalar_V_wind_at_pc_full_level_n(1:nlev,iv), glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)
        call lininterp(glamData_ttt_at_pc_full_level(1:nlev_eff,iv),   hpfull(1:nlev_eff), 1, nlev_eff, &
                       glamData_dynamics%scalar_temp_at_pc_full_level_n(1:nlev,iv)  , glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)
        call lininterp(glamData_qqq_at_pc_full_level(1,1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
                       scalar_tracer_mxrt_at_pc_full_level_n(1,1:nlev,iv), glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)

        ! moist to dry mxrt, assume air only has dry-part and water vapor
        ! i.e., q = s/(1-s)
        SCALAR_TRACER_MXRT_AT_PC_FULL_LEVEL_N(1,1:nlev,iv) = scalar_tracer_mxrt_at_pc_full_level_n(1,1:nlev,iv)/&
                                                        (one-scalar_tracer_mxrt_at_pc_full_level_n(1,1:nlev,iv))
        SCALAR_TRACER_MXRT_AT_PC_FULL_LEVEL_N(2:ntracer,1:nlev,iv) = zero

        do ilev = 1, nlev
           glamData_dynamics%scalar_mif_at_pc_full_level_n(ilev,iv) = one/(one+scalar_tracer_mxrt_at_pc_full_level_n(1,ilev,iv))
        end do

     end do

! find mpressure at model levels

     call calc_hpe_get_full_mass(nlev, mesh%nv_full, 'dtp'                                    , & ! in
                                       glamData_dynamics%scalar_hpressure_at_pc_face_level_n  , & ! in
                                       glamData_dynamics%scalar_hpressure_at_pc_full_level_n  , & ! in
                                       glamData_dynamics%scalar_delhp_at_pc_full_level_n      , & ! in
                                       glamData_dynamics%scalar_mif_at_pc_full_level_n        , & ! in
                                       glamData_dynamics%scalar_mpressure_at_pc_full_level_n  , & ! out
                                       glamData_dynamics%scalar_mpressure_at_pc_face_level_n)     ! out

!
! evaluate geop based on the moist hydrostatic equation
!
     do iv = 1, mesh%nv_full

        scalar_template_1d_nlev_a  = glamData_dynamics%scalar_temp_at_pc_full_level_n(1:nlev,iv)*&
                                     (one+(rvap-rdry)/rdry*(scalar_tracer_mxrt_at_pc_full_level_n(1,1:nlev,iv)/&
                                                (one+scalar_tracer_mxrt_at_pc_full_level_n(1,1:nlev,iv))))      ! Tv
        scalar_template_1d_nlevp_a = glamData_dynamics%scalar_mpressure_at_pc_face_level_n(1:nlevp,iv) ! pressure face
        scalar_template_1d_nlev_b  = glamData_dynamics%scalar_mpressure_at_pc_face_level_n(2:nlevp,iv)-&
                                     glamData_dynamics%scalar_mpressure_at_pc_face_level_n(1:nlev,iv)  ! delp full
        scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)                      ! geopotential at surface

        call calc_hpe_hydro(nlev,nlevp,scalar_template_1d_nlev_a  ,& ! Tv or T
                            scalar_template_1d_nlevp_a ,& ! p or hp
                            scalar_template_1d_nlev_b  ,& ! delp or delhp
                            scalar_template_a          ,& ! phis
                            scalar_template_1d_nlevp_b ,& ! face geop
                            scalar_template_1d_nlev_c )   ! full geop

        SCALAR_GEOPOTENTIAL_AT_PC_FACE_LEVEL_N(1:nlevp,iv)  =  scalar_template_1d_nlevp_b

     end do
!
! common procedures, horizontal momentum on the model mesh
! evaluate normal wind at edge
!
     field_head_2d => null()

     do ie = 1, mesh%ne_halo(1)
        icell1 = mesh%edt_v(1,ie)
        icell2 = mesh%edt_v(2,ie)
        do ilev = 1, nlev
           utmp = half*(glamData_dynamics%scalar_U_wind_at_pc_full_level_n(ilev,icell1)+glamData_dynamics%scalar_U_wind_at_pc_full_level_n(ilev,icell2))
           vtmp = half*(glamData_dynamics%scalar_V_wind_at_pc_full_level_n(ilev,icell1)+glamData_dynamics%scalar_V_wind_at_pc_full_level_n(ilev,icell2))
           call convert_vector_sph2cart(utmp, vtmp, real(mesh%edt_c_p(1:3,ie),r8), vector_velocity)
           SCALAR_NORMAL_VELOCITY_AT_EDGE_FULL_LEVEL_N%f(ilev,ie) = dot_product(vector_velocity, real(mesh%edp_nr(1:3,ie),r8))
        end do
     end do
! exchange data
     call exchange_data_2d_add(mesh,field_head_2d,SCALAR_NORMAL_VELOCITY_AT_EDGE_FULL_LEVEL_N)
     call exchange_data_2d(mesh%local_block,field_head_2d)

! will not be used
     do ie = 1, mesh%ne_full
        if(mesh%edt_mark(ie).eq.1772)then
           SCALAR_NORMAL_VELOCITY_AT_EDGE_FULL_LEVEL_N%f(1:nlev,ie) = zero
        end if
     end do

! obtain mass-pt
     do iv = 1, mesh%nv_full
           SCALAR_MASS_PT_AT_PC_FULL_LEVEL_N(1:nlev,iv)  =  &
                                                  glamData_dynamics%scalar_temp_at_pc_full_level_n(1:nlev,iv)/&
                                                ((glamData_dynamics%scalar_mpressure_at_pc_full_level_n(1:nlev,iv)/p00)**(rdry/cp))*&
                                    (one+ptfactor*SCALAR_TRACER_MXRT_AT_PC_FULL_LEVEL_N(1,1:nlev,iv))*glamData_dynamics%scalar_delhp_at_pc_full_level_n(1:nlev,iv)
     end do

     if(allocated(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)) deallocate(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)
     if(allocated(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)) deallocate(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)
     if(allocated(glamData_dynamics%scalar_delhp_at_pc_full_level_n))     deallocate(glamData_dynamics%scalar_delhp_at_pc_full_level_n)
     if(allocated(glamData_dynamics%scalar_mif_at_pc_full_level_n))       deallocate(glamData_dynamics%scalar_mif_at_pc_full_level_n)
     if(allocated(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)) deallocate(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)
     if(allocated(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)) deallocate(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)

     if(allocated(glamData_dynamics%scalar_U_wind_at_pc_full_level_n))    deallocate(glamData_dynamics%scalar_U_wind_at_pc_full_level_n)
     if(allocated(glamData_dynamics%scalar_V_wind_at_pc_full_level_n))    deallocate(glamData_dynamics%scalar_V_wind_at_pc_full_level_n)
     if(allocated(glamData_dynamics%scalar_temp_at_pc_full_level_n))      deallocate(glamData_dynamics%scalar_temp_at_pc_full_level_n)

    return
   end subroutine grist_glam_data_transform_pmlev

!
! when use this routine, assumption is that data is split from a GRIST global model run that contains the same vertical
! layer setup as this lam run, i.e., global-to-regional downscaling, nlev_lamdata=nlev; so no need for vertical interpolation
! for some model-level vars
!

   subroutine grist_glam_data_transform_grist_mlev(mesh, nlev, nlevp, &
                                              glamData_uuu_at_pc_full_level,&
                                              glamData_vvv_at_pc_full_level,&
                                              glamData_ttt_at_pc_full_level,&
                                              glamData_qqq_at_pc_full_level,&
                                              glamData_hps_at_pc_surface   ,&
                                              scalar_normal_velocity_at_edge_full_level_n, &
                                              scalar_mass_pt_at_pc_full_level_n          , &
                                              scalar_tracer_mxrt_at_pc_full_level_n      , &
                                              scalar_hpressure_at_pc_surface_n           , &
                                              scalar_geopotential_at_pc_face_level_n)
! io
   type(global_domain),   intent(in)    :: mesh
   integer(i4)        ,   intent(in)    :: nlev, nlevp
   real(r8), allocatable, intent(in)    :: glamData_uuu_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: glamData_vvv_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: glamData_ttt_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: glamData_qqq_at_pc_full_level(:,:,:)
   real(r8), allocatable, intent(in)    :: glamData_hps_at_pc_surface(:)
   !inout
   type(scalar_2d_field), intent(inout) :: scalar_normal_velocity_at_edge_full_level_n
   real(r8), allocatable, intent(inout) :: scalar_mass_pt_at_pc_full_level_n(:,:)
   real(r8), allocatable, intent(inout) :: scalar_tracer_mxrt_at_pc_full_level_n(:,:,:)
   real(r8), allocatable, intent(inout) :: scalar_hpressure_at_pc_surface_n(:)
   real(r8), allocatable, intent(inout) :: scalar_geopotential_at_pc_face_level_n(:,:)
!local
   real(r8)                             :: utmp, vtmp
   real(r8)                             :: vector_velocity(3)
! for real-case
   real(r8)                             :: scalar_template_1d_nlev_a(nlev)
   real(r8)                             :: scalar_template_1d_nlev_b(nlev)
   real(r8)                             :: scalar_template_1d_nlev_c(nlev)
   real(r8)                             :: scalar_template_1d_nlevp_a(nlevp)
   real(r8)                             :: scalar_template_1d_nlevp_b(nlevp)
   real(r8)                             :: scalar_template_a
   integer(i4)                          :: iv, ie, icell1, icell2, ilev, itracer
   type(exchange_field_list_2d),pointer :: field_head_2d

!==============================================================
! this entry only takes care of same grist model level
!==============================================================

if(.not.allocated(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)) allocate(glamData_dynamics%scalar_hpressure_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)) allocate(glamData_dynamics%scalar_hpressure_at_pc_full_level_n(nlev ,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_delhp_at_pc_full_level_n))     allocate(glamData_dynamics%scalar_delhp_at_pc_full_level_n(nlev,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mif_at_pc_full_level_n))       allocate(glamData_dynamics%scalar_mif_at_pc_full_level_n(nlev  ,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)) allocate(glamData_dynamics%scalar_mpressure_at_pc_face_level_n(nlevp,mesh%nv_full))
if(.not.allocated(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)) allocate(glamData_dynamics%scalar_mpressure_at_pc_full_level_n(nlev ,mesh%nv_full))

     do iv = 1, mesh%nv_full
!
! evaluate dry mass at GRIST dry-air model level
!
        scalar_hpressure_at_pc_surface_n(iv)                               = glamData_hps_at_pc_surface(iv)
        glamData_dynamics%scalar_hpressure_at_pc_face_level_n(1:nlevp,iv)  = eta_face_a(1:nlevp)*p0+eta_face_b(1:nlevp)*glamData_hps_at_pc_surface(iv) ! Pa
        glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev ,iv)  = eta_full_a(1:nlev) *p0+eta_full_b(1:nlev) *glamData_hps_at_pc_surface(iv)
        glamData_dynamics%scalar_delhp_at_pc_full_level_n(1:nlev,iv)       = glamData_dynamics%scalar_hpressure_at_pc_face_level_n(2:nlevp,iv)-&
                                                                             glamData_dynamics%scalar_hpressure_at_pc_face_level_n(1:nlev ,iv)
!
! interpolate from ERA pressure level to GRIST model level based on dry air mass
!
        !call lininterp (glamData_uuu_at_pc_full_level(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
        !                glamData_dynamics%scalar_U_wind_at_pc_full_level_n(1:nlev,iv), glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)
        !call lininterp (glamData_vvv_at_pc_full_level(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
        !                glamData_dynamics%scalar_V_wind_at_pc_full_level_n(1:nlev,iv), glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)
        !call lininterp (glamData_ttt_at_pc_full_level(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
        !                glamData_dynamics%scalar_temp_at_pc_full_level_n(1:nlev,iv)  , glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)
        !call lininterp (glamData_qqq_at_pc_full_level(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
        !                glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1:nlev,iv), glamData_dynamics%scalar_hpressure_at_pc_full_level_n(1:nlev,iv), nlev)
        ! moist to dry
        !glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1:nlev,iv) = glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1:nlev,iv)/&
        !                                                                (one-glamData_dynamics%scalar_tracer_mxrt_at_pc_full_level_n(1:nlev,iv))

        do itracer = 1, ntracer
           scalar_tracer_mxrt_at_pc_full_level_n(itracer,1:nlev,iv) = glamData_qqq_at_pc_full_level(itracer,1:nlev,iv)
        end do

        do ilev = 1, nlev
           glamData_dynamics%scalar_mif_at_pc_full_level_n(ilev,iv) = one/(one+sum(glamData_qqq_at_pc_full_level(mif_index(1:nmif),ilev,iv)))
        end do

     end do

! find mpressure at model levels
     call calc_hpe_get_full_mass(nlev, mesh%nv_full, 'dtp'                                    , & ! in
                                       glamData_dynamics%scalar_hpressure_at_pc_face_level_n  , & ! in
                                       glamData_dynamics%scalar_hpressure_at_pc_full_level_n  , & ! in
                                       glamData_dynamics%scalar_delhp_at_pc_full_level_n      , & ! in
                                       glamData_dynamics%scalar_mif_at_pc_full_level_n        , & ! in
                                       glamData_dynamics%scalar_mpressure_at_pc_full_level_n  , & ! out
                                       glamData_dynamics%scalar_mpressure_at_pc_face_level_n)     ! out

     !glamData_dynamics%scalar_delp_at_pc_full_level_n(1:nlev,:) = glamData_dynamics%scalar_mpressure_at_pc_face_level_n(2:nlevp,:)-&
     !                                                             glamData_dynamics%scalar_mpressure_at_pc_face_level_n(1:nlev,:)
!
! evaluate geop based on the moist hydrostatic equation, VIP
!
     do iv = 1, mesh%nv_full

        scalar_template_1d_nlev_a  = glamData_ttt_at_pc_full_level(1:nlev,iv)*&
                                     (one+(rvap-rdry)/rdry*(glamData_qqq_at_pc_full_level(1,1:nlev,iv)/&
                                                       (one+glamData_qqq_at_pc_full_level(1,1:nlev,iv))))  ! Tv

        scalar_template_1d_nlevp_a = glamData_dynamics%scalar_mpressure_at_pc_face_level_n(1:nlevp,iv) ! pressure face
        scalar_template_1d_nlev_b  = glamData_dynamics%scalar_mpressure_at_pc_face_level_n(2:nlevp,iv)-&
                                     glamData_dynamics%scalar_mpressure_at_pc_face_level_n(1:nlev ,iv) ! delp full
        scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)                      ! geopotential at surface

        call calc_hpe_hydro(nlev,nlevp,scalar_template_1d_nlev_a  ,& ! Tv or T
                            scalar_template_1d_nlevp_a ,& ! p or hp
                            scalar_template_1d_nlev_b  ,& ! delp or delhp
                            scalar_template_a          ,& ! phis
                            scalar_template_1d_nlevp_b ,& ! face geop
                            scalar_template_1d_nlev_c )   ! full geop

        scalar_geopotential_at_pc_face_level_n(1:nlevp,iv)  =  scalar_template_1d_nlevp_b

     end do
!
! Normal wind
!
     field_head_2d => null()

     do ie = 1, mesh%ne_halo(1)
        icell1 = mesh%edt_v(1,ie)
        icell2 = mesh%edt_v(2,ie)
        do ilev = 1, nlev
           utmp = half*(glamData_uuu_at_pc_full_level(ilev,icell1)+glamData_uuu_at_pc_full_level(ilev,icell2))
           vtmp = half*(glamData_vvv_at_pc_full_level(ilev,icell1)+glamData_vvv_at_pc_full_level(ilev,icell2))
           call convert_vector_sph2cart(utmp, vtmp, real(mesh%edt_c_p(1:3,ie),r8), vector_velocity)
           scalar_normal_velocity_at_edge_full_level_n%f(ilev,ie) = dot_product(vector_velocity, real(mesh%edp_nr(1:3,ie),r8))
        end do
     end do

     call exchange_data_2d_add(mesh,field_head_2d,scalar_normal_velocity_at_edge_full_level_n)
     call exchange_data_2d(mesh%local_block,field_head_2d)

! exchange data
     !if(mesh%hasBdyDomain)then
        !do ie = 1, mesh%ne_full
     do ie = 1, mesh%ne_full
        if(mesh%edt_mark(ie).eq.1772)then
           scalar_normal_velocity_at_edge_full_level_n%f(1:nlev,ie) = zero
        end if
     end do
     !end if

     !call exchange_data_2d_add(mesh,field_head_2d,scalar_normal_velocity_at_edge_full_level_n)
     !call exchange_data_2d(mesh%local_block,field_head_2d)

! MASS-PT
    do iv = 1, mesh%nv_full
          scalar_mass_pt_at_pc_full_level_n(1:nlev,iv)  =  glamData_ttt_at_pc_full_level(1:nlev,iv)/&
                                                ((glamData_dynamics%scalar_mpressure_at_pc_full_level_n(1:nlev,iv)/p00)**(rdry/cp))*&
                                    (one+ptfactor*glamData_qqq_at_pc_full_level(1,1:nlev,iv))*glamData_dynamics%scalar_delhp_at_pc_full_level_n(1:nlev,iv)
    end do

    if(allocated(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)) deallocate(glamData_dynamics%scalar_hpressure_at_pc_face_level_n)
    if(allocated(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)) deallocate(glamData_dynamics%scalar_hpressure_at_pc_full_level_n)
    if(allocated(glamData_dynamics%scalar_delhp_at_pc_full_level_n))     deallocate(glamData_dynamics%scalar_delhp_at_pc_full_level_n)
    if(allocated(glamData_dynamics%scalar_mif_at_pc_full_level_n))       deallocate(glamData_dynamics%scalar_mif_at_pc_full_level_n)
    if(allocated(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)) deallocate(glamData_dynamics%scalar_mpressure_at_pc_face_level_n)
    if(allocated(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)) deallocate(glamData_dynamics%scalar_mpressure_at_pc_full_level_n)

    return
   end subroutine grist_glam_data_transform_grist_mlev

   subroutine handle_time
     integer(i4)  :: endInMon(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/) ! 1-12, 365 days
#ifdef USE_LEAP_YEAR
     if((mod(t1_year,4).eq.0.and.mod(t1_year,100).ne.0).or.(mod(t1_year,400).eq.0))then
        endInMon=(/31,29,31,30,31,30,31,31,30,31,30,31/)
     else
        endInMon=(/31,28,31,30,31,30,31,31,30,31,30,31/)
     end if
#endif

    t2_sec = t1_sec+glamDataInterval

     if(t2_sec.ge.86400)then
        t2_sec = t2_sec-86400
        t2_day = t1_day+1
     else
        t2_day = t1_day
     end if
     if(t2_day.gt.endInMon(t1_mon))then
        t2_day = t2_day-endInMon(t1_mon)
        t2_mon = t1_mon+1
     else
        t2_mon = t1_mon
     end if
     if(t2_mon.gt.12)then
        t2_mon  = t2_mon-12
        t2_year = t1_year+1
     else
        t2_year = t1_year
     end if

!================================================
! obtain char
!================================================
     write(char_t1_year,'(i4)') t1_year
     write(char_t2_year,'(i4)') t2_year

     if(t1_mon.ge.10)then
        write(char_t1_mon,'(i2)') t1_mon
     else
        write(char_t1_mon,'(i1)') t1_mon
        char_t1_mon="0"//char_t1_mon
     end if

     if(t2_mon.ge.10)then
        write(char_t2_mon,'(i2)') t2_mon
     else
        write(char_t2_mon,'(i1)') t2_mon
        char_t2_mon="0"//char_t2_mon
     end if

      if(t1_day.ge.10)then
         write(char_t1_day,'(i2)') t1_day
      else
         write(char_t1_day,'(i1)') t1_day
         char_t1_day="0"//char_t1_day
      end if

      if(t2_day.ge.10)then
         write(char_t2_day,'(i2)') t2_day
      else
         write(char_t2_day,'(i1)') t2_day
         char_t2_day="0"//char_t2_day
      end if
 
      if(t1_sec.eq.0)then ! t1_sec only can be 00000, 10800, 21600, every 3-hour
         char_t1_sec="00000"
      elseif(t1_sec.eq.3600.or.t1_sec.eq.7200)then
         write(char_t1_sec,'(i4)') t1_sec
         char_t1_sec = "0"//char_t1_sec
      else
         write(char_t1_sec,'(i5)') t1_sec
      end if

      if(t2_sec.eq.0)then ! t2_sec only can be 00000, 10800, 21600, every 3-hour
         char_t2_sec="00000"
      elseif(t2_sec.eq.3600.or.t2_sec.eq.7200)then
         write(char_t2_sec,'(i4)') t2_sec
         char_t2_sec = "0"//char_t2_sec
      else
         write(char_t2_sec,'(i5)') t2_sec
      end if

!================================================
! obtain char
!================================================

   end subroutine handle_time

!================================================
! Add secondary laplacian-style forcing
!================================================

   subroutine add_laplacian_forcing_velocity(mesh, nLevel, dtime, &
                                                   scalar_lgs_normal_velocity_at_edge,&
                                                   scalar_mod_normal_velocity_at_edge)

    use grist_dycore_diffusion_module, only: perot_weight_at_pc, perot_weight_at_dc
! io
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: nLevel
     real(r8),              intent(in)    :: dtime
     real(r8), allocatable, intent(in)    :: scalar_lgs_normal_velocity_at_edge(:,:)
     type(scalar_2d_field), intent(inout) :: scalar_mod_normal_velocity_at_edge
! local
     real(r8)                             :: velocity_vector_at_pc_full_level(3,nLevel,mesh%nv_halo(1))
     real(r8)                             :: velocity_vector_at_dc_full_level(3,nLevel,mesh%nt_halo(1))
     integer(i4)                          :: iv, it, ie, ilev, inb, cell_mark
     integer(i4)                          :: pc1, pc2, dc1, dc2              ! cell index
     real(r8)                             :: un_pc1, un_pc2, un_dc1, un_dc2  ! normal wind follow edp's nr
     real(r8)                             :: ut_pc1, ut_pc2, ut_dc1, ut_dc2  ! tangen wind follow edp's tg
     real(r8)                             :: pc1pc2(3), dc1dc2(3), flag1, flag2

!$omp parallel  private(iv,inb,ie,ilev)    
!$omp do schedule(dynamic,50) 
      do iv  = 1, mesh%nv_halo(1)
         velocity_vector_at_pc_full_level(:,1:nLevel,iv) = zero
         do inb  = 1, mesh%vtx_nnb(iv)
            ie = mesh%vtx_ed(inb,iv)
            do ilev = 1, nLevel
                velocity_vector_at_pc_full_level(:,ilev,iv) = velocity_vector_at_pc_full_level(:,ilev,iv)+perot_weight_at_pc(:,inb,iv)*&
                                                   (scalar_lgs_normal_velocity_at_edge(ilev,ie)-scalar_mod_normal_velocity_at_edge%f(ilev,ie))
            end do
         end do
      end do
!$omp end do nowait
!$omp end parallel 

!$omp parallel  private(it,inb,ilev,ie) 
!$omp do schedule(dynamic,50) 
      do it  = 1, mesh%nt_halo(1)
         velocity_vector_at_dc_full_level(:,:,it) = zero
         do inb = 1, mesh%tri_nnb(it)
            ie = mesh%tri_ed(inb,it)
            do ilev = 1, nLevel
               velocity_vector_at_dc_full_level(:,ilev,it) = velocity_vector_at_dc_full_level(:,ilev,it)+perot_weight_at_dc(:,inb,it)*&
                                                   (scalar_lgs_normal_velocity_at_edge(ilev,ie)-scalar_mod_normal_velocity_at_edge%f(ilev,ie))
            end do
         end do
      end do
!$omp end do nowait
!$omp end parallel 

!-----------------------------------------------------------------------
!         pc2
!          |
!    dc1---e---dc2
!          |
!         pc1
!-----------------------------------------------------------------------

!$omp parallel  private(ie,pc1,pc2,dc1,dc2,pc1pc2,dc1dc2,flag1,flag2,ilev,un_pc1,un_pc2,un_dc1,un_dc2,ut_pc1,ut_pc2,ut_dc1,ut_dc2,un_pc1pc2,ut_pc1pc2,un_dc1dc2,ut_dc1dc2) 
!$omp do schedule(dynamic,20) 
      DO ie = 1, mesh%ne_compute

         IF(mesh%edt_mark(ie).eq.188)then

         pc1     = mesh%edt_v(1,ie)
         pc2     = mesh%edt_v(2,ie)
         dc1     = mesh%edp_v(1,ie)
         dc2     = mesh%edp_v(2,ie)
         cell_mark = min(mesh%plg_mark(pc1),mesh%plg_mark(pc2))

         do ilev = 1, nLevel
            un_pc1  = dot_product(velocity_vector_at_pc_full_level(1:3,ilev,pc1),mesh%edp_nr(1:3,ie))
            un_pc2  = dot_product(velocity_vector_at_pc_full_level(1:3,ilev,pc2),mesh%edp_nr(1:3,ie))
            un_dc1  = dot_product(velocity_vector_at_dc_full_level(1:3,ilev,dc1),mesh%edp_nr(1:3,ie))
            un_dc2  = dot_product(velocity_vector_at_dc_full_level(1:3,ilev,dc2),mesh%edp_nr(1:3,ie))
                                                                                                
            ut_pc1  = dot_product(velocity_vector_at_pc_full_level(1:3,ilev,pc1),mesh%edp_tg(1:3,ie))
            ut_pc2  = dot_product(velocity_vector_at_pc_full_level(1:3,ilev,pc2),mesh%edp_tg(1:3,ie))
            ut_dc1  = dot_product(velocity_vector_at_dc_full_level(1:3,ilev,dc1),mesh%edp_tg(1:3,ie))
            ut_dc2  = dot_product(velocity_vector_at_dc_full_level(1:3,ilev,dc2),mesh%edp_tg(1:3,ie))

     ! second-order laplacian of normal wind, because we should use 1/2(edt&edp's length) for double center differences but not,
     ! so this laplacian needs a *4 when using outside, for laplacian vi, we also normalize it like this;
     ! why not *4 here? currently just for bit-bit regression!
     ! math: ((upc1-un)/(0.5*edt_leng)-(un-upc2)/(0.5*edt_leng*rearth))/(0.5*edt_leng*rearth)
     ! math: ((udc1-un)/(0.5*edp_leng)-(un-udc2)/(0.5*edp_leng*rearth))/(0.5*edp_leng*rearth)
     ! rearth^2 is selected outside and canceld by earth in 3DX
     ! 12 = f1*4.

            scalar_mod_normal_velocity_at_edge%f(ilev,ie) = scalar_mod_normal_velocity_at_edge%f(ilev,ie)-dtime*12._r8*mesh%edt_leng(ie)/rearth*&
             ((un_pc1+un_pc2-2._r8*(scalar_lgs_normal_velocity_at_edge(ilev,ie)-scalar_mod_normal_velocity_at_edge%f(ilev,ie)))/((mesh%edt_leng(ie))**2)+&
              (un_dc1+un_dc2-2._r8*(scalar_lgs_normal_velocity_at_edge(ilev,ie)-scalar_mod_normal_velocity_at_edge%f(ilev,ie)))/((mesh%edp_leng(ie))**2))
         end do

         END IF

       END DO
!$omp end do nowait
!$omp end parallel

       return
    end subroutine add_laplacian_forcing_velocity

    subroutine add_laplacian_forcing_scalar_1d(mesh, dtime, &
                                                     scalar_lgs_at_pc , &
                                                     scalar_mod_at_pc )
! io
      type(global_domain),   intent(in)    :: mesh
      real(r8),              intent(in)    :: dtime
      real(r8), allocatable, intent(in)    :: scalar_lgs_at_pc(:)
      real(r8), allocatable, intent(inout) :: scalar_mod_at_pc(:)
! local
      type(scalar_1d_field)  :: gradient_at_prime_edge
      integer(i4)            :: ie, iv, inb,ilev, v1, v2
      real(r8)               :: v1v2(3), flag, div_sum
      type(exchange_field_list_1d),pointer :: field_head_1d

      field_head_1d=>null()

      if(.not.allocated(gradient_at_prime_edge%f)) allocate(gradient_at_prime_edge%f(mesh%ne_full))
!
! gradient at edge, counter edge's normal direction
!
      do ie = 1, mesh%ne_halo(1)     ! global index
         v1      = mesh%edt_v(1,ie)
         v2      = mesh%edt_v(2,ie)
         v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
         flag    = sign(1._r8,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
         gradient_at_prime_edge%f(ie) = flag*(scalar_lgs_at_pc(v2)-scalar_mod_at_pc(v2)-&
                                              scalar_lgs_at_pc(v1)+scalar_mod_at_pc(v1))/(rearth*mesh%edt_leng(ie))
      end do
      call exchange_data_1d_add(mesh,field_head_1d,gradient_at_prime_edge)
      call exchange_data_1d(mesh%local_block,field_head_1d)
!
! divergence at cell
!
        do iv = 1, mesh%nv_full
          IF(mesh%plg_mark(iv).eq.1881.or.mesh%plg_mark(iv).eq.1882.or.mesh%plg_mark(iv).eq.1883.or.mesh%plg_mark(iv).eq.1884.or.mesh%plg_mark(iv).eq.1885)then
           div_sum = zero
           do inb = 1, mesh%vtx_nnb(iv)
              ie       = mesh%vtx_ed(inb,iv)
              div_sum  = div_sum+gradient_at_prime_edge%f(ie)*mesh%plg_nr(inb,iv)*mesh%edp_leng(ie)
           end do
           scalar_mod_at_pc(iv) = scalar_mod_at_pc(iv)-dtime*9._r8*div_sum*mesh%vtxCellLeng(iv)/mesh%plg_areag(iv)
          END IF
        end do

      if(allocated(gradient_at_prime_edge%f)) deallocate(gradient_at_prime_edge%f)

      return
    end subroutine add_laplacian_forcing_scalar_1d

    subroutine add_laplacian_forcing_scalar_2d(mesh, nLevel, dtime , &
                                                  scalar_lgs_at_pc , &
                                                  scalar_mod_at_pc )
! io
      type(global_domain),   intent(in)    :: mesh
      integer(i4)          , intent(in)    :: nLevel
      real(r8),              intent(in)    :: dtime
      real(r8), allocatable, intent(in)    :: scalar_lgs_at_pc(:,:)
      real(r8), allocatable, intent(inout) :: scalar_mod_at_pc(:,:)
! local
      type(scalar_2d_field)  :: gradient_at_prime_edge
      integer(i4)            :: ie, iv, inb,ilev, v1, v2
      real(r8)               :: v1v2(3), flag, div_sum
      type(exchange_field_list_2d),pointer :: field_head_2d

      field_head_2d=>null()

      if(.not.allocated(gradient_at_prime_edge%f)) allocate(gradient_at_prime_edge%f(nLevel,mesh%ne_full))
!
! gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1)     ! global index
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(1._r8,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           do ilev = 1, nLevel
              gradient_at_prime_edge%f(ilev,ie) = flag*(scalar_lgs_at_pc(ilev,v2)-scalar_mod_at_pc(ilev,v2)-&
                                                        scalar_lgs_at_pc(ilev,v1)+scalar_mod_at_pc(ilev,v1))/(rearth*mesh%edt_leng(ie))
           end do
        end do
        call exchange_data_2d_add(mesh,field_head_2d,gradient_at_prime_edge)
        call exchange_data_2d(mesh%local_block,field_head_2d)
!
! divergence at cell
!
        do iv = 1, mesh%nv_full
          IF(mesh%plg_mark(iv).eq.1881.or.mesh%plg_mark(iv).eq.1882.or.mesh%plg_mark(iv).eq.1883.or.mesh%plg_mark(iv).eq.1884.or.mesh%plg_mark(iv).eq.1885)then
           do ilev = 1, nLevel
              div_sum = zero
              do inb = 1, mesh%vtx_nnb(iv)
                 ie  = mesh%vtx_ed(inb,iv)
                 div_sum  = div_sum+gradient_at_prime_edge%f(ilev,ie)*mesh%plg_nr(inb,iv)*mesh%edp_leng(ie)
              end do
              !tend_laplacian_2nd_at_pc(ilev,iv) = div_sum/((rearth**2)*mesh%plg_areag(iv))
              ! this is a reduced version..
              scalar_mod_at_pc(ilev,iv) = scalar_mod_at_pc(ilev,iv)-dtime*9._r8*div_sum*mesh%vtxCellLeng(iv)/mesh%plg_areag(iv)
           end do
         end if
        end do

        if(allocated(gradient_at_prime_edge%f)) deallocate(gradient_at_prime_edge%f)

        return
    end subroutine add_laplacian_forcing_scalar_2d
#endif

  end module grist_datam_glam_data_module
