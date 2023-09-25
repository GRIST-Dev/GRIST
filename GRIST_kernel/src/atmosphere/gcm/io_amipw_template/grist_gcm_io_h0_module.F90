
module grist_gcm_io_h0_module

  use grist_constants,                  only: i4
  use grist_hpe_constants,              only: eta_full_a, eta_full_b, eta_face_a, eta_face_b
  use grist_data_types,                 only: scalar_1d_field, scalar_2d_field
  use grist_domain_types,               only: global_domain
  use grist_dycore_vars_module,         only: dycoreVarGeography 

  use grist_nml_module,                 only: outdir, nlev, ntracer,               &
                                              start_ymd, start_tod,                &
                                              working_mode,                        & 
                                              fname_output,                        & 
                                              model_timestep, grid_info
  
  use grist_time_manager,               only: get_curr_cdate

  use grist_gcm_diagnose_h0_module,     only: gcm_h0_accu_physics_variables, &
                                              gcm_h0_dump_physics_variables, &
                                              gcm_h0_rest_physics_variables, &
                                              diag_physics_vars_1d_h0,       &
                                              diag_physics_vars_2d_h0
  use grist_dycore_diagnose_module_2d,  only: dycore_diagnose_variables

  implicit none
  private
  public  :: gcm_output_history_h0

! private
  integer(i4)      :: restart_timestep
  character(len=2) :: old_mon
  character(len=4) :: old_year

  contains

  subroutine gcm_output_history_h0(mesh, itimestep)
! io
  type(global_domain), intent(in)     :: mesh
  integer(i4)        , intent(in)     :: itimestep
! local
  character(len=4)                    :: cyear
  character(len=2)                    :: cmon
  character(len=2)                    :: cday
  character(len=5)                    :: csec
  logical                             :: does_file_exist

  call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, &
                      cyear, cmon, cday, csec) 

    inquire(file=trim(outdir)//'step_restart.txt',exist=does_file_exist)
    if(does_file_exist) then
       open(1,file=trim(outdir)//'step_restart.txt',status='old')
       read(1,*) restart_timestep
       close(1)
    else
       restart_timestep = 1
    endif  !if(does_file_exist)

  if(itimestep .eq. 1.or.itimestep.eq.restart_timestep)then
     old_mon = cmon
     old_year= cyear
  end if  

  if(cmon .eq. old_mon)then
    call dycore_diagnose_variables(mesh)
    call gcm_h0_accu_physics_variables
    return
  end if
  
  if(cmon .ne. old_mon)then
     call gcm_h0_dump_physics_variables
     call gcm_mon_output_file(mesh, old_year, old_mon)
     call gcm_h0_rest_physics_variables

     old_mon = cmon
     if(old_year .ne. cyear) old_year= cyear
     
     call dycore_diagnose_variables(mesh)
     call gcm_h0_accu_physics_variables
  end if

  end subroutine gcm_output_history_h0

  subroutine gcm_mon_output_file(mesh, cyear, cmon)
  use grist_util_module,               only: write_string
  use grist_fileio_list_2d_module_par, only: wrap_output_init_2d      ,&
                                             wrap_add_field_2d        ,&
                                             wrap_output_2d_group     ,&
                                             wrap_output_2d_group_sp  ,&
                                             wrap_output_clean_2d
  use grist_fileio_list_1d_module_par, only: wrap_output_init_1d      ,&
                                             wrap_add_field_1d        ,&
                                             wrap_output_1d_group     ,&
                                             wrap_output_1d_group_sp  ,&
                                             wrap_output_clean_1d 


! io
  type(global_domain), intent(in)     :: mesh
  character(len=4),    intent(in)     :: cyear
  character(len=2),    intent(in)     :: cmon
! local
  character(128)                      :: c_glevel
  integer(i4)                         :: ilev
  type(scalar_1d_field)               :: vert_level, vertp_level
  type(scalar_1d_field)               :: hyam, hybm, hyai, hybi
  type(scalar_1d_field)               :: scalar_io_1d_template
  type(scalar_2d_field)               :: scalar_io_2d_template

  if(.not.allocated(vert_level%f)) allocate(vert_level%f(nlev))
  if(.not.allocated(vertp_level%f))allocate(vertp_level%f(nlev+1))
  if(.not.allocated(hyam%f))       allocate(hyam%f(nlev))
  if(.not.allocated(hybm%f))       allocate(hybm%f(nlev))
  if(.not.allocated(hyai%f))       allocate(hyai%f(nlev+1))
  if(.not.allocated(hybi%f))       allocate(hybi%f(nlev+1))

  vert_level%pos = 8
  vertp_level%pos = 9
  hyam%pos = 8
  hybm%pos = 8
  hyai%pos = 9
  hybi%pos = 9

  do ilev = 1, nlev
    vert_level%f(ilev)  = ilev
    vertp_level%f(ilev) = ilev
  end do
  vertp_level%f(nlev+1) = nlev+1

  hyam%f(:)  = eta_full_a
  hybm%f(:)  = eta_full_b
  hyai%f(:)  = eta_face_a
  hybi%f(:)  = eta_face_b

  call write_string(mesh%glevel, c_glevel)
!
! 2d file
!
  if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

  fname_output="GRIST.ATM."//trim(grid_info)//"."//trim(working_mode)//".MonAvg."//trim(cyear)//"-"//trim(cmon)//".2d.h0.nc"

! output
  call wrap_output_init_2d(mesh)

  call wrap_add_field_2d(diag_physics_vars_2d_h0%uwind           ,"uPC"     )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%vwind           ,"vPC"     )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%omegaFull       ,"omegaFull")
#ifdef DYCORE_NHD
  call wrap_add_field_2d(diag_physics_vars_2d_h0%wwwFace         ,"wwwFace" )
#endif
  call wrap_add_field_2d(diag_physics_vars_2d_h0%phiFace         ,"geopotentialFace")
  call wrap_add_field_2d(diag_physics_vars_2d_h0%mpressureFace   ,"mpressureFace" )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%temp            ,"temp"    )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qv              ,"qv"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%cloud           ,"cloudF"  )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%relhum          ,"relhum"  )
#ifdef AMIPW_PHYSICS
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qc              ,"qc"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qr              ,"qr"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qi              ,"qi"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qs              ,"qs"      )
  call wrap_add_field_2d(diag_physics_vars_2d_h0%qg              ,"qg"      )
#endif

#ifdef SPIO
  call wrap_output_2d_group_sp(mesh, outdir, fname_output)  
#else
  call wrap_output_2d_group(mesh, outdir, fname_output)
#endif

! clean
  call wrap_output_clean_2d()  

!
! 1d file
!
  if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

  fname_output="GRIST.ATM."//trim(grid_info)//"."//trim(working_mode)//".MonAvg."//trim(cyear)//"-"//trim(cmon)//".1d.h0.nc"

! output
  call wrap_output_init_1d(mesh)

  call wrap_add_field_1d(dycoreVarGeography%scalar_lon_at_pc      , "lon")
  call wrap_add_field_1d(dycoreVarGeography%scalar_lat_at_pc      , "lat")
  call wrap_add_field_1d(vert_level                               , "nlev"  )
  call wrap_add_field_1d(vertp_level                              , "nlevp" )
  call wrap_add_field_1d(hyam                                     , "hyam"  )
  call wrap_add_field_1d(hybm                                     , "hybm"  )
  call wrap_add_field_1d(hyai                                     , "hyai"  )
  call wrap_add_field_1d(hybi                                     , "hybi"  )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%ps               , "ps"    )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%ts               , "ts"    )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%shflx            , "shflx" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%qflx             , "qflx"  )

  call wrap_add_field_1d(diag_physics_vars_1d_h0%prect            , "prect" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%precc            , "precc" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%precl            , "precl" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%snowl            , "precSnowl" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%grapl            , "precGrapl" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldtot           , "cldtot" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldlow           , "cldlow" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldmed           , "cldmed" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%cldhgh           , "cldhgh" )

  call wrap_add_field_1d(diag_physics_vars_1d_h0%swcf             , "swcf"  )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%lwcf             , "lwcf"  )

  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwut            , "flwut" ) 
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwdt            , "flwdt" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwus            , "flwus" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwds            , "flwds" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwutc           , "flwutc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwdtc           , "flwdtc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwusc           , "flwusc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%flwdsc           , "flwdsc")

  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswut            , "fswut" ) 
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswdt            , "fswdt" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswus            , "fswus" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswds            , "fswds" )
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswutc           , "fswutc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswdtc           , "fswdtc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswusc           , "fswusc")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%fswdsc           , "fswdsc")
#ifdef AMIPW_PHYSICS
  call wrap_add_field_1d(diag_physics_vars_1d_h0%asdir            , "asdir")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%asdif            , "asdif")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%aldir            , "aldir")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%aldif            , "aldif")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%emiss            , "emiss")
  call wrap_add_field_1d(diag_physics_vars_1d_h0%snowhland        , "snowhland")
#endif
  !call wrap_add_field_1d(diag_physics_vars_h0%z500                   , "z500"  )
  !call wrap_add_field_1d(diag_physics_vars_h0%u200                   , "u200"  )
  !call wrap_add_field_1d(diag_physics_vars_h0%u850                   , "u850"  )
  !call wrap_add_field_1d(diag_physics_vars_h0%v850                   , "v850"  )
#ifdef SPIO
  call wrap_output_1d_group_sp(mesh, outdir, fname_output)  
#else
  call wrap_output_1d_group(mesh, outdir, fname_output)
#endif

! clean
  call wrap_output_clean_1d()  

  if(allocated(vert_level%f))  deallocate(vert_level%f)
  if(allocated(vertp_level%f)) deallocate(vertp_level%f)
  if(allocated(hyam%f))        deallocate(hyam%f)
  if(allocated(hybm%f))        deallocate(hybm%f)
  if(allocated(hyai%f))        deallocate(hyai%f)
  if(allocated(hybi%f))        deallocate(hybi%f)

  end subroutine gcm_mon_output_file

 end module grist_gcm_io_h0_module
