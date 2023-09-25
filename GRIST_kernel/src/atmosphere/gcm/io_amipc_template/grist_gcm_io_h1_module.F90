
 module grist_gcm_io_h1_module

  use grist_lib
  use grist_constants,        only: i4, i8
  use grist_hpe_constants,    only: eta_full_a, eta_full_b, eta_face_a, eta_face_b
  use grist_domain_types,     only: global_domain, group_comm, block_structure
  use grist_data_types,       only: scalar_1d_field

  use grist_nml_module,       only: fname_output,          &
                                    working_mode,          &
                                    advection_scheme,      &
                                    conserve_scheme,       &
                                    model_timestep,        &
                                    testcase,              &
                                    outdir, nlev, ntracer, &
                                    comm_group_size,       &
                                    start_ymd, start_tod, grid_info

  use grist_time_manager,     only: get_current_date, get_curr_cdate ,&
                                    get_curr_date

  use grist_fileio_list_3d_module_par, only: wrap_output_init_3d     ,&
                                             wrap_add_field_3d       ,&
                                             wrap_output_3d_group    ,&
                                             wrap_output_3d_group_sp ,&
                                             wrap_output_clean_3d

  use grist_fileio_list_2d_module_par, only: wrap_output_init_2d      ,&
                                             wrap_add_field_2d        ,&
                                             wrap_output_2d_group     ,&
                                             wrap_output_2d_group_sp  ,&
                                             wrap_output_clean_2d

  use grist_fileio_list_1d_module_par, only: wrap_output_init_1d    ,&
                                             wrap_add_field_1d      ,&
                                             wrap_output_1d_group   ,&
                                             wrap_output_1d_group_sp,&
                                             wrap_output_clean_1d   ,&
                                             wrap_read_1d_group_rst
! data module
  use grist_dycore_vars_module
  use grist_tracer_transport_vars_module
  use grist_dtp_vars_module
  use grist_physics_data_structure, only: pstate
  use grist_datam_static_data_module,  only: staticData_phis_at_pc_surface, &
                                             staticData_sst_at_pc_surface
  use grist_datam_initial_data_module, only: initialData_uuu_at_pc_full_level, &
                                             initialData_vvv_at_pc_full_level, &
                                             initialData_ttt_at_pc_full_level, &
                                             initialData_ps_at_pc_surface
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,  only: pstate_cam
  use grist_physics_update,        only: old_time_level
#endif
!
! gcm-diagnose
!
  use grist_gcm_diagnose_h1_module,only: gcm_h1_accu_physics_variables, &
                                         gcm_h1_dump_physics_variables, &
                                         gcm_h1_rest_physics_variables, &
                                         diag_physics_vars_h1
  use grist_dycore_diagnose_module_2d,  only: dycore_diagnose_variables

  implicit none

  private

  public  :: gcm_output_history_h1, &
             gcm_output_history_h1_1d, & ! just to plz compilation
             gcm_output_history_h1_2d, &
             gcm_output_history_h1_3d, &
             gcm_output_atm_h1_file

  contains
!================================================
! [1] Construct filename according to itimestep
! [2] Calculate Diag Vars using Prog vars
! [3] Do an Output Flow
!================================================

  subroutine gcm_output_history_h1(mesh,itimestep,nsteps,history_freq)
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps 
   integer(i4)        , intent(in)     :: history_freq 
! local
   integer(i4) :: yr,mn,dy,sc 

    call get_curr_date(start_ymd,start_tod,itimestep,model_timestep,yr,mn,dy,sc)
    !-------we only need JJA data after the spinup period, study of diurnal cycle, LiXH---------
    !if(yr .gt. start_ymd/10000 .and. mn .ge. 6 .and. mn .le. 8)then    
    !if(yr .gt. start_ymd/10000 )then                !Omit the spin-up period, LiXH
    call gcm_h1_accu_physics_variables
    if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,history_freq)==0))then
        ! only for those vars not in Diag-h1, otherwise we should put this line before h1_accu as in h0-diag;
        ! currently we do not need to accumulate variables from dycore_diagnoise to h1's Diag vars
        call dycore_diagnose_variables(mesh)
        call gcm_h1_dump_physics_variables
        call gcm_output_atm_h1_file(mesh, itimestep)
        call gcm_h1_rest_physics_variables
    end if
    !------------------------------------------------------------------
    !end if
    return
  end subroutine gcm_output_history_h1

  subroutine gcm_output_history_h1_1d(mesh,itimestep,nsteps,history_freq)
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps
   integer(i4)        , intent(in)     :: history_freq
  end subroutine gcm_output_history_h1_1d

  subroutine gcm_output_history_h1_2d(mesh,itimestep,nsteps,history_freq)
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps
   integer(i4)        , intent(in)     :: history_freq
  end subroutine gcm_output_history_h1_2d

  subroutine gcm_output_history_h1_3d(mesh,itimestep,nsteps,history_freq)
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps
   integer(i4)        , intent(in)     :: history_freq
  end subroutine gcm_output_history_h1_3d

  subroutine gcm_output_atm_h1_file(mesh, itimestep)
    use grist_util_module, only: write_string
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
! local
   character(128)                      :: c_glevel
   character(128)                      :: c_testcase
   character(len=4)                    :: cyear
   character(len=2)                    :: cmon
   character(len=2)                    :: cday
   character(len=5)                    :: csec
   character(len=5)                    :: day
   character(len=5)                    :: sec

!================================================
! Info for constructing filename
!================================================

   call write_string(mesh%glevel     , c_glevel)
   call write_string(testcase        , c_testcase)
! continous date from 00000-00000
   call get_current_date(itimestep,model_timestep,day,sec)
! calendar date
   call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, &
                       cyear, cmon, cday, csec)
!
! ATM 3d file, not output when using dycore mode
!
    IF(.not.trim(working_mode).eq.'dycore')then
       call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme, &
                                    cyear,cmon,cday,csec,day,sec,"3d.h1",fname_output)
   !    if(.false.)then
       call write_atm_file3d(mesh,fname_output)
   !    end if
    Endif
!
! ATM 2d file
!
    call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"2d.h1",fname_output)
  !  if(.false.)then
    call write_atm_file2d(mesh,fname_output)
  !  end if
!
! ATM 1d file
!
    call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"1d.h1",fname_output)
    call write_atm_file1d(mesh,fname_output)
!
! LND 1d file
!
!    call obtain_lnd_fname_output(c_glevel,working_mode,conserve_scheme,&
!                                 cyear,cmon,cday,csec,day,sec,"1d.h1",fname_output)
!    call write_lnd_file1d(mesh,fname_output)

    return
  end subroutine gcm_output_atm_h1_file

  subroutine write_atm_file3d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output
! local
   type(scalar_3d_field)               :: tracer_mxrt
   integer(i4)                         :: ilev, itracer

      
    call wrap_output_init_3d(mesh)
    call wrap_add_field_3d(tracerVarCellFull%scalar_tracer_mxrt_n,"tracerMxrt")
!    call wrap_add_field_3d(scalar_tracer_mass_at_pc_full_level_n,"tracerMass")

#ifdef SPIO
    call wrap_output_3d_group_sp(mesh, outdir,fname_output)
#else
    call wrap_output_3d_group(   mesh, outdir,fname_output)
#endif
    call wrap_output_clean_3d()

    return
  end subroutine write_atm_file3d

  subroutine write_atm_file2d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output
! local
   type(scalar_2d_field)               :: density
   integer(i4)                         :: ilev

!
! output 2d variables
!
!   density   = scalar_delphi_at_pc_full_level_n
!   density%f = -scalar_delphi_at_pc_full_level_n%f/scalar_delhp_at_pc_full_level_n%f

! start output
    call wrap_output_init_2d(mesh)
! add variables
    IF(.not.trim(working_mode).eq.'tracer')then
    call wrap_add_field_2d(dycoreVarCellFull%scalar_U_wind_n           ,"uPC"      )
    call wrap_add_field_2d(dycoreVarCellFull%scalar_V_wind_n           ,"vPC"      )
!    call wrap_add_field_2d(scalar_www_at_pc_face_level_n              ,"wwwFace"  )
!    call wrap_add_field_2d(scalar_www_at_pc_full_level_n              ,"wwwFull"  )
!    call wrap_add_field_2d(scalar_abs_vor_at_dc_full_level_n          ,"aVor" )
!    call wrap_add_field_2d(scalar_rel_vor_at_dc_full_level_n          ,"rVor" )
!    call wrap_add_field_2d(scalar_divergence_at_pc_full_level_n       ,"div"  )
    call wrap_add_field_2d(dycoreVarCellFull%scalar_temp_n             ,"temperature" )
!    call wrap_add_field_2d(scalar_potential_temp_at_pc_full_level_n   ,"pt")
!    call wrap_add_field_2d(scalar_potential_temp_at_pc_full_level_iniupt,"pt0")
    call wrap_add_field_2d(dycoreVarCellFull%scalar_geopotential_n     ,"geopotentialFull")
!    call wrap_add_field_2d(scalar_geopotential_at_pc_face_level_n     ,"geopotentialFace")
    call wrap_add_field_2d(dycoreVarCellFull%scalar_pressure_n         ,"pressureFull"  )
!    call wrap_add_field_2d(scalar_pressure_at_pc_face_level_n         ,"pressureFace"  )
!    call wrap_add_field_2d(scalar_delp_at_pc_full_level_np1           ,"delp"  )
!    call wrap_add_field_2d(scalar_delhp_at_pc_full_level_n            ,"delhp" )
!    call wrap_add_field_2d(scalar_omega_at_pc_full_level_n            ,"omega" )
    !call wrap_add_field_2d(scalar_U_wind_at_edge_full_level_n         ,"uEdge"      )
    !call wrap_add_field_2d(scalar_V_wind_at_edge_full_level_n         ,"vEdge"      )
    !call wrap_add_field_2d(scalar_hpressure_at_pc_full_level_n        ,"hPressureFull" )
    !call wrap_add_field_2d(scalar_hpressure_at_pc_face_level_n        ,"hPressureFace" )
!    call wrap_add_field_2d(scalar_delphi_at_pc_full_level_n           ,"delphi")
!    call wrap_add_field_2d(density                                    ,"density")


    if(trim(working_mode).ne.'dycore'.and.trim(working_mode).ne.'tracer')then
!    call wrap_add_field_2d(scalar_thetav_at_pc_full_level_n      ,"ptv")
    end if
! static
!    call wrap_add_field_2d(scalar_static_albedo_at_pc, "albedo2d")
    END IF
!
! outputing
!

#ifdef SPIO
    call wrap_output_2d_group_sp(mesh, outdir, fname_output)
#else
    call wrap_output_2d_group(mesh, outdir, fname_output)
#endif
!
! clean
!
    call wrap_output_clean_2d()

    return
  end subroutine write_atm_file2d

  subroutine write_atm_file1d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output
! local
   type(scalar_1d_field)               :: vert_level, vertp_level
   type(scalar_1d_field)               :: hyam, hybm, hyai, hybi
   type(scalar_1d_field)               :: scalar_io_1d_template
   integer(i4)                         :: ilev

    if(.not.allocated(vert_level%f))  allocate(vert_level%f(nlev))
    if(.not.allocated(vertp_level%f)) allocate(vertp_level%f(nlev+1))
    if(.not.allocated(hyam%f))        allocate(hyam%f(nlev))
    if(.not.allocated(hybm%f))        allocate(hybm%f(nlev))
    if(.not.allocated(hyai%f))        allocate(hyai%f(nlev+1))
    if(.not.allocated(hybi%f))        allocate(hybi%f(nlev+1))
    if(.not.allocated(scalar_io_1d_template%f)) allocate(scalar_io_1d_template%f(mesh%nv_halo(1)))

    vert_level%pos  = 8
    vertp_level%pos = 9
    hyam%pos        = 8
    hybm%pos        = 8
    hyai%pos        = 9
    hybi%pos        = 9

    do ilev = 1, nlev
       vert_level%f(ilev)  = ilev
       vertp_level%f(ilev) = ilev
    end do
       vertp_level%f(nlev+1) = nlev+1

    hyam%f(:)  = eta_full_a
    hybm%f(:)  = eta_full_b
    hyai%f(:)  = eta_face_a
    hybi%f(:)  = eta_face_b

! start output
    call wrap_output_init_1d(mesh)
! add variables
! geometric
    call wrap_add_field_1d(dycoreVarGeography%scalar_lon_at_pc                         , "lon_nv")
    call wrap_add_field_1d(dycoreVarGeography%scalar_lat_at_pc                         , "lat_nv")
!    call wrap_add_field_1d(scalar_lon_at_dc                         , "lon_nt")
!    call wrap_add_field_1d(scalar_lat_at_dc                         , "lat_nt")
!    call wrap_add_field_1d(scalar_lon_at_edge                       , "lon_ne")
!    call wrap_add_field_1d(scalar_lat_at_edge                       , "lat_ne")
!    call wrap_add_field_1d(scalar_area_at_pc                        , "areaPC")
!    call wrap_add_field_1d(scalar_area_at_dc                        , "areaDC")
!    call wrap_add_field_1d(scalar_leng_at_edp                       , "lengEdp")
!    call wrap_add_field_1d(scalar_leng_at_edt                       , "lengEdt")
! vertical
!    call wrap_add_field_1d(vert_level                               , "nlev")
!    call wrap_add_field_1d(vertp_level                              , "nlevp")
!    call wrap_add_field_1d(hyam                                     , "hyam")
!    call wrap_add_field_1d(hybm                                     , "hybm")
!    call wrap_add_field_1d(hyai                                     , "hyai")
!    call wrap_add_field_1d(hybi                                     , "hybi")
! vars
!    call wrap_add_field_1d(scalar_geopotential_at_pc_surface_n      , "phis")
!    call wrap_add_field_1d(scalar_hpressure_at_pc_surface_n         , "hps")
!    call wrap_add_field_1d(scalar_pressure_at_pc_surface_n          , "ps")

    if(trim(working_mode).ne.'dycore'.and.trim(working_mode).ne.'tracer')then
!       call wrap_add_field_1d(diag_physics_vars_h1%precc           , "precc")
!       call wrap_add_field_1d(diag_physics_vars_h1%precl           , "precl")
!       call wrap_add_field_1d(diag_physics_vars_h1%prect           , "prect")
!       if(trim(working_mode).eq.'amipc'.or.trim(working_mode).eq.'amipw')then
!         call wrap_add_field_1d(pstate%scalar_precc_surface           , "precc")
         call wrap_add_field_1d(pstate%scalar_prect_surface           , "prect")
!         call wrap_add_field_1d(pstate%ts_at_pc_surface               , "ts")
!       end if
       if(trim(working_mode).eq.'amipc')then
!          call wrap_add_field_1d(diag_physics_vars_h1%flwut           , "flwut")
#ifdef AMIPC_PHYSICS
          call wrap_add_field_1d(pstate_cam%flwut_at_pc_top            , "fluwt")
#endif
!         call wrap_add_field_1d(pstate%sst_at_pc_surface              , "sstPstate")
!         call wrap_add_field_1d(staticData_sst_at_pc_surface          , "sstRaw")
!         call wrap_add_field_1d(pstate%landfrac_at_pc_surface         , "landfrac")
!         call wrap_add_field_1d(pstate%ocnfrac_at_pc_surface          , "ocnfrac")
!         call wrap_add_field_1d(pstate%icefrac_at_pc_surface          , "icefrac")
!         call wrap_add_field_1d(pstate%snowhland_at_pc_surface        , "snowdepth")

#ifdef USE_NOAHMP
!         call wrap_add_field_1d(pstate%atm_in_lwup_at_pc_surface      , "lwup")
!         call wrap_add_field_1d(pstate%atm_in_shflx_at_pc_surface     , "shflx")
!         scalar_io_1d_template%f(:mesh%nv_halo(1)) = pstate%atm_in_qflx_at_pc_surface%f(1,:mesh%nv_halo(1))
!         scalar_io_1d_template%pos              = pstate%atm_in_qflx_at_pc_surface%pos
!         call wrap_add_field_1d(scalar_io_1d_template                 , "qflx")
!         call wrap_add_field_1d(pstate%atm_in_taux_at_pc_surface      , "taux")
!         call wrap_add_field_1d(pstate%atm_in_tauy_at_pc_surface      , "tauy")
!         call wrap_add_field_1d(pstate%atm_in_asdir_at_pc_surface     , "asdir")
!         call wrap_add_field_1d(pstate%atm_in_asdif_at_pc_surface     , "asdif")
!         call wrap_add_field_1d(pstate%atm_in_aldir_at_pc_surface     , "aldir")
!         call wrap_add_field_1d(pstate%atm_in_aldif_at_pc_surface     , "aldif")
!         call wrap_add_field_1d(pstate%atm_out_fswds_at_pc_surface    , "swdown")
!         call wrap_add_field_1d(pstate%atm_out_flwds_at_pc_surface    , "lwdown")
#endif
       end if
    end if
! outputing
#ifdef SPIO
    call wrap_output_1d_group_sp(mesh, outdir, fname_output)
#else
    call wrap_output_1d_group(mesh, outdir, fname_output)
#endif
! clean
    call wrap_output_clean_1d()
    deallocate(vert_level%f,vertp_level%f,hyam%f,hybm%f,hyai%f,hybi%f)
    deallocate(scalar_io_1d_template%f)

   return
  end subroutine write_atm_file1d

  subroutine write_lnd_file1d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output

    call wrap_output_init_1d(mesh)

#ifdef SPIO
    call wrap_output_1d_group_sp(mesh, outdir, fname_output)
#else
    call wrap_output_1d_group(mesh, outdir, fname_output)
#endif
! clean
    call wrap_output_clean_1d()
    return
  end subroutine write_lnd_file1d

  subroutine obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,dim,fname_output)
! io
  character(len=*),  intent(in)    :: c_glevel, working_mode, conserve_scheme
  character(len=*),  intent(in)    :: cyear, cmon, cday, csec, day, sec, dim
  character(len=*),  intent(inout) :: fname_output

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_output="GRIST.ATM."//trim(grid_info)//"."//trim(working_mode)//"."//&
#ifdef CDATE
                               trim(cyear)//"-"//&
                               trim(cmon)//"-"//&
                               trim(cday)//"-"//&
                               trim(csec)//&
#else
                               trim(day)//"-"//&
                               trim(sec)//&
#endif
                               "."//trim(dim)//".nc"
    return
  end subroutine obtain_atm_fname_output

  subroutine obtain_lnd_fname_output(c_glevel,working_mode,conserve_scheme,&
                                     cyear,cmon,cday,csec,day,sec,dim,fname_output)
! io
  character(len=*),  intent(in)    :: c_glevel, working_mode, conserve_scheme
  character(len=*),  intent(in)    :: cyear, cmon, cday, csec, day, sec, dim
  character(len=*),  intent(inout) :: fname_output

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_output="GRIST.LND."//trim(grid_info)//"."// trim(working_mode)//"."//&
#ifdef CDATE
                               trim(cyear)//"-"//&
                               trim(cmon)//"-"//&
                               trim(cday)//"-"//&
                               trim(csec)//&
#else
                               trim(day)//"-"//&
                               trim(sec)//&
#endif
                                "."//trim(dim)//".nc"
    return
  end subroutine obtain_lnd_fname_output

  end module grist_gcm_io_h1_module
