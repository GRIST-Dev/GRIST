
 !================================================
 !  Created by zhangyi on 16/8/5.
 !
 !  Restart module
 !  Rev: 2020Dec, yiz adds restart vars for phy2&lsm
 !  tested ok for current configuration
 !================================================

 module grist_gcm_restart_module

  use grist_lib
  use grist_domain_types, only: global_domain
  use grist_constants
  use grist_nml_module,   only: fname_restart,    &
                                working_mode,     &
                                model_timestep,   &
                                dycore_timestep,  &
                                tracer_timestep,  &
                                outdir,           &
                                conserve_scheme,  &
                                nlev,             &
                                nlevp,            &
                                ntracer,nspecies, &
                                ptend_heat_rk_on, ptend_wind_rk_on, ptendSubDiabPhys, &
                                ptend_dycore_heat_f1_on, ptend_dycore_wind_f1_on, ptend_dycore_f1_on, test_real_case, &
                                start_ymd, start_tod, model_timestep, grid_info, large_atm_file_on

  use grist_fileio_list_1d_module_par, only: wrap_output_init_1d , &
                                             wrap_add_field_1d   , &
                                             wrap_output_1d_group, &
                                             wrap_output_clean_1d, &
                                             wrap_read_1d_group_rst

  use grist_fileio_list_2d_module_par, only: wrap_output_init_2d ,   &
                                             wrap_add_field_2d   ,   &
                                             wrap_output_2d_group,   &
                                             wrap_output_clean_2d,   &
                                             wrap_read_2d_group_rst

  use grist_fileio_list_3d_module_par, only: wrap_output_init_3d ,   &
                                             wrap_add_field_3d   ,   &
                                             wrap_output_3d_group,   &
                                             wrap_output_clean_3d,   &
                                             wrap_read_3d_group_rst
  !use grist_list_array
  use grist_time_manager,  only: get_current_date, get_curr_cdate
  use grist_dycore_vars_module, only: dycoreVarSurface,  &
                                      dycoreVarEdgeFull, &
                                      dycoreVarCellFull, &
                                      dycoreVarCellFace

  use grist_tracer_transport_vars_module, only: tracerVarCellFull, &
                                                tracerVarCellFace
  use grist_physics_data_structure,  only: pstate, ptend_f1, ptend_f2, ptend_rk
#ifdef AMIPW_PHYSICS
  use grist_PhysW_data_structure,    only: psurf_wrf, ptend_wrf, pstate_wrf
  use grist_PhysW_nml_module,        only: step_cu, step_bl
  use grist_PhysW_nml_module,        only: PhysW_mp_scheme
#endif

#ifdef AMIPC_PHYSICS
  use grist_gcm_restart_phys1_module, only: gcm_write_restart_file_phy1,    &
                                            gcm_read_restart_file_phy1
#endif

  use grist_util_module,    only: write_string
  use grist_data_types,     only: scalar_3d_field, scalar_2d_field, scalar_1d_field, data2d_temp
  implicit none

  private

  public        :: gcm_read_restart_file    ,&
                   gcm_write_restart_file_h0,&
                   gcm_write_restart_file_h1,&
                   restart_timestep,         &
                   itimestep_pre

    integer(i4) :: restart_timestep
    integer(i4) :: itimestep_pre = 0
    character(len=2) :: old_mon
    character(len=4) :: old_year

  contains

!================================================
!  If Called, read from restart file. 
!  Overwrite data in swe_vars_module
!  Although we write so many data into restart,
!  not all of them are used.
!================================================

  subroutine gcm_read_restart_file(mesh)
! io
    type(global_domain), intent(inout), target :: mesh
    character(len=4)    :: cyear
    character(len=2)    :: cmon
    character(len=2)    :: cday

    call gcm_read_restart_file_dyn(mesh,cyear,cmon,cday)  ! obtain cyear, cmon, cday

#ifdef AMIPW_PHYSICS
    call gcm_read_restart_file_phy2(mesh,cyear,cmon,cday) ! input cyear, cmon, cday
#endif

#ifdef AMIPC_PHYSICS
    call gcm_read_restart_file_phy1(mesh,cyear,cmon,cday)
#endif

#ifdef USE_NOAHMP
    if(test_real_case) call gcm_read_restart_file_lnd(mesh,cyear,cmon,cday)
#endif

    return
  end subroutine gcm_read_restart_file

  subroutine gcm_read_restart_file_dyn(mesh,cyear,cmon,cday)
    use grist_dycore_gcd_recon_module_2d, only: vector_recon_perot_edge2cell_uv_2d
    ! io
    type(global_domain), intent(inout), target :: mesh
    character(len=4),    intent(out)           :: cyear
    character(len=2),    intent(out)           :: cmon
    character(len=2),    intent(out)           :: cday
    ! local
    character(128)            :: c_glevel
    character(128)            :: fname_restart1,fname_restart2,fname_restart3
    character(len=5)          :: day
    character(len=5)          :: sec
    integer                   :: i,j
    logical                   :: does_file_exist
    integer(i8)               :: lev,levp
    real(r8)                  :: timestep

    character(len=5)          :: csec

!================================================
! Construct Filename
!================================================

! read the timestep at which we generate this restart file

    lev  = nlev
    levp = nlevp
    call write_string(mesh%glevel, c_glevel)

!================================================
!  if file exist, read the restart_timestep
!================================================
    inquire(file=trim(outdir)//'step_restart.txt',exist=does_file_exist)
    if(does_file_exist) then
       open(1,file=trim(outdir)//'step_restart.txt',status='old')
       read(1,*) restart_timestep
       close(1)

       open(2,file=trim(outdir)//'step_stop.txt',status='old')
       read(2,*) itimestep_pre
       close(2)
    else
       restart_timestep = 1
    endif  !if(does_file_exist)
    
    SELECT CASE(trim(working_mode))
    CASE('dycore')
      timestep = dycore_timestep
    CASE('tracer')
      timestep = tracer_timestep
    CASE('dtp','amipc','amipw')
      timestep = model_timestep
    END SELECT

    call get_current_date(restart_timestep,timestep,day,sec)
! restart is written at itimestep+1, so here use  minus 1 to recover cyear, cmon, cday
    call get_curr_cdate(start_ymd, start_tod, restart_timestep-1, model_timestep, cyear, cmon, cday, csec)

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)
     
    fname_restart1 = "GRIST.RST.Dyn."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2 = "GRIST.RST.Dyn."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart3 = "GRIST.RST.Dyn."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"

    IF(does_file_exist)THEN

    call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'hps', 0,              dycoreVarSurface%scalar_hpressure_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'uNormal',     lev, 6, dycoreVarEdgeFull%scalar_normal_velocity_n%f)
#ifdef DYCORE_NHD
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'wwwFace',     levp,0, dycoreVarCellFace%scalar_www_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'phi',         levp,0, dycoreVarCellFace%scalar_phi_n%f)
#endif
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'massPt',      lev, 0, dycoreVarCellFull%scalar_mass_pt_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ptemp',       lev, 0, dycoreVarCellFull%scalar_potential_temp_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'pressureFull',lev, 0, dycoreVarCellFull%scalar_pressure_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'mpressureFull',lev,0, dycoreVarCellFull%scalar_mpressure_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'pressureFace',levp,0, dycoreVarCellFace%scalar_pressure_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'delhpFull',   lev, 0, dycoreVarCellFull%scalar_delhp_n%f)
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'delhpFace',   levp,0, dycoreVarCellFace%scalar_delhp_n%f)
!    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'hpressureFace',levp,0,dycoreVarCellFace%scalar_hpressure_n%f)
!    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'mifEdge',      lev,6, scalar_mif_at_edge_full_level_n%f)
!    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'mifFull',      lev,0, tracerVarCellFull%scalar_mif_n%f)
!    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'mifFace',      levp,0,scalar_mif_at_pc_face_level_n%f)

! as the following vars are accumulated during gcm-diag, reconstruct them here
!
     call vector_recon_perot_edge2cell_uv_2d(mesh, dycoreVarEdgeFull%scalar_normal_velocity_n%f,&
                                                   dycoreVarCellFull%scalar_U_wind_n%f,&
                                                   dycoreVarCellFull%scalar_V_wind_n%f,nlev)
     dycoreVarSurface%scalar_pressure_n%f = dycoreVarCellFace%scalar_pressure_n%f(nlevp,:)

! tendency restarts
    if(ptend_heat_rk_on)then
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ptendRKMassPt', lev, 0, ptend_rk%tend_mass_pt_at_pc_full_level%f)
    end if
    if(ptend_wind_rk_on)then
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ptendRKNormalV',lev, 6, ptend_rk%tend_normal_velocity_at_edge_full_level%f)
    end if
    if(ptendSubDiabPhys)then
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ptendF2MassPt', lev, 0, ptend_f2%tend_mass_pt_at_pc_full_level%f)
    end if
    if(ptend_dycore_heat_f1_on.or.ptend_dycore_wind_f1_on.or.ptend_dycore_f1_on)then
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ptendF1MassPt' ,lev, 0, ptend_f1%tend_mass_pt_at_pc_full_level%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ptendF1NormalV',lev, 6, ptend_f1%tend_normal_velocity_at_edge_full_level%f)
    end if

    if(trim(working_mode) .ne. 'dycore') then
       if(.false.)then ! this is only for ad-hoc reading G11 with ntracer=6
          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMxrt.1.nc", 'tracerMxrt' ,lev, 0, tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:))
          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMass.1.nc", 'tracerMass' ,lev, 0, tracerVarCellFull%scalar_tracer_mass_n%f(1,:,:))

          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMxrt.2.nc", 'tracerMxrt' ,lev, 0, tracerVarCellFull%scalar_tracer_mxrt_n%f(2,:,:))
          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMass.2.nc", 'tracerMass' ,lev, 0, tracerVarCellFull%scalar_tracer_mass_n%f(2,:,:))

          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMxrt.3.nc", 'tracerMxrt' ,lev, 0, tracerVarCellFull%scalar_tracer_mxrt_n%f(3,:,:))
          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMass.3.nc", 'tracerMass' ,lev, 0, tracerVarCellFull%scalar_tracer_mass_n%f(3,:,:))

          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMxrt.4.nc", 'tracerMxrt' ,lev, 0, tracerVarCellFull%scalar_tracer_mxrt_n%f(4,:,:))
          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMass.4.nc", 'tracerMass' ,lev, 0, tracerVarCellFull%scalar_tracer_mass_n%f(4,:,:))

          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMxrt.5.nc", 'tracerMxrt' ,lev, 0, tracerVarCellFull%scalar_tracer_mxrt_n%f(5,:,:))
          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMass.5.nc", 'tracerMass' ,lev, 0, tracerVarCellFull%scalar_tracer_mass_n%f(5,:,:))

          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMxrt.6.nc", 'tracerMxrt' ,lev, 0, tracerVarCellFull%scalar_tracer_mxrt_n%f(6,:,:))
          call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, "tracerMass.6.nc", 'tracerMass' ,lev, 0, tracerVarCellFull%scalar_tracer_mass_n%f(6,:,:))
       else
          call wrap_read_3d_group_rst(mesh%gcomm_read, outdir, fname_restart3, 'tracerMxrt', lev, 0, tracerVarCellFull%scalar_tracer_mxrt_n%f)
          call wrap_read_3d_group_rst(mesh%gcomm_read, outdir, fname_restart3, 'tracerMass', lev, 0, tracerVarCellFull%scalar_tracer_mass_n%f)
       end if
    end if

    END IF

    return
  end subroutine gcm_read_restart_file_dyn

!================================================
! write restart as h1 file (N-step)
!================================================

  subroutine gcm_write_restart_file_h1(mesh,itimestep,nsteps,restart_freq)
! io
    type(global_domain), intent(in)  :: mesh
    integer(i4),         intent(in)  :: itimestep
    real(r8),            intent(in)  :: nsteps 
    integer(i4),         intent(in)  :: restart_freq

    character(len=4)                 :: cyear
    character(len=2)                 :: cmon
    character(len=2)                 :: cday
    character(len=5)                 :: csec

    call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, cyear, cmon, cday, csec)

    if((itimestep==nsteps .or.mod(itimestep,restart_freq)==0))then
       call gcm_write_restart_file(mesh,(itimestep+1),cyear,cmon,cday)
       if(mpi_rank()==0)then
           print*,"=========================================================="
           print*,"    gcm_write_restart_file_h1, write step =", itimestep, restart_freq
           print*,"=========================================================="
       endif
    end if
    return
  end subroutine gcm_write_restart_file_h1

!================================================
! write restart as h0 file (instant state at every
! new month)
!================================================

  subroutine gcm_write_restart_file_h0(mesh,itimestep,nsteps)
! io
    type(global_domain), intent(in)  :: mesh
    integer(i4),         intent(in)  :: itimestep
    real(r8),            intent(in)  :: nsteps
! local
    character(len=4)                 :: cyear
    character(len=2)                 :: cmon
    character(len=2)                 :: cday
    character(len=5)                 :: csec

     call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, cyear, cmon, cday, csec)

     if(itimestep .eq. 1 .or. itimestep.eq.restart_timestep)then
        old_mon = cmon
        old_year= cyear
     end if

     if(cmon .eq. old_mon)return
     if(cmon .ne. old_mon)then
        call gcm_write_restart_file(mesh,(itimestep+1),cyear,cmon,cday)
        if(mpi_rank()==0)then
           print*,"=========================================================="
           print*,"    gcm_write_restart_file_h0, write step =", itimestep
           print*,"=========================================================="
        endif
        old_mon = cmon
        if(old_year .ne. cyear) old_year = cyear
     end if

    return
  end subroutine gcm_write_restart_file_h0

!================================================
! [1] Construct Filename
! [2] Write Restart File
!================================================

  subroutine gcm_write_restart_file(mesh,itimestep,cyear,cmon,cday)
! io
    type(global_domain), intent(in)  :: mesh
    integer(i4),         intent(in)  :: itimestep
    character(len=4)                 :: cyear
    character(len=2)                 :: cmon
    character(len=2)                 :: cday


    call gcm_write_restart_file_dyn(mesh,itimestep,cyear,cmon,cday)

#ifdef AMIPW_PHYSICS
    call gcm_write_restart_file_phy2(mesh,itimestep,cyear,cmon,cday)
#endif

#ifdef AMIPC_PHYSICS
    call gcm_write_restart_file_phy1(mesh,itimestep,cyear,cmon,cday)
#endif

#ifdef USE_NOAHMP
    ! Aqua-Planet also can compile NOAHMP, so avoid here
    if(test_real_case) call gcm_write_restart_file_lnd(mesh,itimestep,cyear,cmon,cday)
#endif

! write the timestep at which we generate last restart file
    if(mpi_rank() == 0) then
      open(1,file=trim(outdir)//'step_restart.txt',status='unknown')
      write(1,*) itimestep
      close(1)
    endif

    return

  end subroutine gcm_write_restart_file

  subroutine gcm_write_restart_file_dyn(mesh,itimestep,cyear,cmon,cday)
! io
    type(global_domain), intent(in)  :: mesh
    integer(i4),         intent(in)  :: itimestep
    character(len=4)                 :: cyear
    character(len=2)                 :: cmon
    character(len=2)                 :: cday
! local
    character(128)     :: c_glevel
    character(128)     :: fname_restart1, fname_restart2, fname_restart3
    !character(len=5)   :: myid_c
    character(len=5)   :: day
    character(len=5)   :: sec
    !integer(i4)        :: myid, i, j
    !logical            :: timecheck
    real(r8)           :: timestep

!================================================
! [1] Construct Filename
!================================================

    SELECT CASE(trim(working_mode))
    CASE('dycore')
      timestep = dycore_timestep
    CASE('tracer')
      timestep = tracer_timestep
    CASE('dtp','amipc','amipw')
      timestep = model_timestep
    END SELECT

    call write_string(mesh%glevel     , c_glevel)

    call get_current_date(itimestep,timestep,day,sec)

! default, use glevel
    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_restart1 = "GRIST.RST.Dyn."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2 = "GRIST.RST.Dyn."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart3 = "GRIST.RST.Dyn."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"

!================================================
! [2] Write Restart File
!================================================

    call wrap_output_init_1d(mesh)
    call wrap_add_field_1d(dycoreVarSurface%scalar_hpressure_n,"hps")
    call wrap_output_1d_group(mesh, outdir,fname_restart1)
    call wrap_output_clean_1d()

    call wrap_output_init_2d(mesh)
    call wrap_add_field_2d(dycoreVarEdgeFull%scalar_normal_velocity_n,"uNormal")
#ifdef DYCORE_NHD
    call wrap_add_field_2d(dycoreVarCellFace%scalar_www_n              ,"wwwFace")
    call wrap_add_field_2d(dycoreVarCellFace%scalar_phi_n              ,"phi")
#endif
    call wrap_add_field_2d(dycoreVarCellFull%scalar_mass_pt_n          ,"massPt")
    call wrap_add_field_2d(dycoreVarCellFull%scalar_potential_temp_n   ,"ptemp")
    call wrap_add_field_2d(dycoreVarCellFull%scalar_delhp_n            ,"delhpFull")
    call wrap_add_field_2d(dycoreVarCellFace%scalar_delhp_n            ,"delhpFace")
    call wrap_add_field_2d(dycoreVarCellFull%scalar_pressure_n         ,"pressureFull")
    call wrap_add_field_2d(dycoreVarCellFace%scalar_pressure_n         ,"pressureFace")
    call wrap_add_field_2d(dycoreVarCellFull%scalar_mpressure_n        ,"mpressureFull") ! for hdc-cutdk-omega
   ! call wrap_add_field_2d(dycoreVarCellFace%scalar_hpressure_n        ,"hpressureFace")
   ! call wrap_add_field_2d(scalar_mif_at_edge_full_level_n            ,"mifEdge")
   ! call wrap_add_field_2d(tracerVarCellFull%scalar_mif_n              ,"mifFull")
   ! call wrap_add_field_2d(scalar_mif_at_pc_face_level_n              ,"mifFace")
    if(ptend_heat_rk_on) call wrap_add_field_2d(ptend_rk%tend_mass_pt_at_pc_full_level            ,"ptendRKMassPt")
    if(ptend_wind_rk_on) call wrap_add_field_2d(ptend_rk%tend_normal_velocity_at_edge_full_level  ,"ptendRKNormalV")
    if(ptendSubDiabPhys) call wrap_add_field_2d(ptend_f2%tend_mass_pt_at_pc_full_level            ,"ptendF2MassPt")
    if(ptend_dycore_heat_f1_on.or.ptend_dycore_wind_f1_on.or.ptend_dycore_f1_on)then
       call wrap_add_field_2d(ptend_f1%tend_mass_pt_at_pc_full_level            ,"ptendF1MassPt")
       call wrap_add_field_2d(ptend_f1%tend_normal_velocity_at_edge_full_level  ,"ptendF1NormalV")
    end if
    call wrap_output_2d_group(mesh, outdir,fname_restart2)
    call wrap_output_clean_2d()

    if(trim(working_mode) .ne. 'dycore') then
      call wrap_output_init_3d(mesh)
      call wrap_add_field_3d(tracerVarCellFull%scalar_tracer_mxrt_n    ,"tracerMxrt")
      call wrap_add_field_3d(tracerVarCellFull%scalar_tracer_mass_n    ,"tracerMass")
      call wrap_output_3d_group(mesh, outdir,fname_restart3)
      call wrap_output_clean_3d()
    end if

    return
  end subroutine gcm_write_restart_file_dyn

#ifdef AMIPW_PHYSICS

!================================================
! Restart info for using WRF physics (Phy2)
!================================================

  subroutine gcm_write_restart_file_phy2(mesh,itimestep,cyear,cmon,cday)
! io
    type(global_domain), intent(in)  :: mesh
    integer(i4),         intent(in)  :: itimestep
    character(len=4)                 :: cyear
    character(len=2)                 :: cmon
    character(len=2)                 :: cday
! local
    character(128)     :: c_glevel
    character(128)     :: fname_restart1, fname_restart2, fname_restart3
    !character(len=5)   :: myid_c
    character(len=5)   :: day
    character(len=5)   :: sec
    integer(i4)        :: ncell, ndim3, itracer, icell, ilev
    !logical            :: timecheck
    real(r8)           :: timestep
    type(scalar_3d_field) :: output3d_temp1
    type(scalar_2d_field) :: output2d_temp1, output2d_temp2,  output2d_temp3,  output2d_temp4,  output2d_temp5, output2d_temp6, output2d_temp7
    type(scalar_2d_field) :: output2d_temp8, output2d_temp9,  output2d_temp10, output2d_temp11, output2d_temp12,output2d_temp13
    type(scalar_1d_field) :: output1d_temp1, output1d_temp2,  output1d_temp3,  output1d_temp4,  output1d_temp5, output1d_temp6, output1d_temp7, output1d_temp8
    type(scalar_1d_field) :: output1d_temp9, output1d_temp10, output1d_temp11, output1d_temp12, output1d_temp13,output1d_temp14,output1d_temp15,output1d_temp16

!================================================
! [1] Construct Filename
!================================================

    timestep = model_timestep
    ncell    = mesh%nv_halo(1)

    call write_string(mesh%glevel     , c_glevel)
    call get_current_date(itimestep,timestep,day,sec)

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_restart1 = "GRIST.RST.Phy."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2 = "GRIST.RST.Phy."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
!
! write flux that affecets diagnose, because rad is not called each step, so
! restart will zero some flux states
!
    call wrap_allocate_data1d(mesh,output1d_temp1)
    call wrap_allocate_data1d(mesh,output1d_temp2)
    call wrap_allocate_data1d(mesh,output1d_temp3)
    call wrap_allocate_data1d(mesh,output1d_temp4)
    call wrap_allocate_data1d(mesh,output1d_temp5)
    call wrap_allocate_data1d(mesh,output1d_temp6)
    call wrap_allocate_data1d(mesh,output1d_temp7)
    call wrap_allocate_data1d(mesh,output1d_temp8)
    call wrap_allocate_data1d(mesh,output1d_temp9)
    call wrap_allocate_data1d(mesh,output1d_temp10)
    call wrap_allocate_data1d(mesh,output1d_temp11)
    call wrap_allocate_data1d(mesh,output1d_temp12)
    call wrap_allocate_data1d(mesh,output1d_temp13)
    call wrap_allocate_data1d(mesh,output1d_temp14)
    call wrap_allocate_data1d(mesh,output1d_temp15)
    call wrap_allocate_data1d(mesh,output1d_temp16)

    output1d_temp1%f(1:ncell) = pstate_wrf%lwupt(1:ncell,1)
    output1d_temp2%f(1:ncell) = pstate_wrf%lwupb(1:ncell,1)
    output1d_temp3%f(1:ncell) = pstate_wrf%lwdnb(1:ncell,1)
    output1d_temp4%f(1:ncell) = pstate_wrf%lwuptc(1:ncell,1)
    output1d_temp5%f(1:ncell) = pstate_wrf%lwupbc(1:ncell,1)
    output1d_temp6%f(1:ncell) = pstate_wrf%lwdnbc(1:ncell,1)
    output1d_temp7%f(1:ncell) = pstate_wrf%lwcf(1:ncell,1)
    output1d_temp8%f(1:ncell) = pstate_wrf%swupt(1:ncell,1)
    output1d_temp9%f(1:ncell) = pstate_wrf%swdnt(1:ncell,1)
    output1d_temp10%f(1:ncell)= pstate_wrf%swupb(1:ncell,1)
    output1d_temp11%f(1:ncell)= pstate_wrf%swdnb(1:ncell,1)
    output1d_temp12%f(1:ncell)= pstate_wrf%swuptc(1:ncell,1)
    output1d_temp13%f(1:ncell)= pstate_wrf%swdntc(1:ncell,1)
    output1d_temp14%f(1:ncell)= pstate_wrf%swupbc(1:ncell,1)
    output1d_temp15%f(1:ncell)= pstate_wrf%swdnbc(1:ncell,1)
    output1d_temp16%f(1:ncell)= pstate_wrf%swcf(1:ncell,1)

!================================================
! [2] Write Restart File
!================================================		

    call wrap_output_init_1d(mesh)

    call wrap_add_field_1d(pstate%ts_at_pc_surface,           "ts")
    call wrap_add_field_1d(pstate%scalar_rainc_surface,       "rainc")
    call wrap_add_field_1d(pstate%scalar_rainnc_surface,      "rainnc")
    call wrap_add_field_1d(pstate%scalar_snownc_surface,      "snownc")
    call wrap_add_field_1d(pstate%scalar_grapnc_surface,      "grapnc")
    call wrap_add_field_1d(pstate%scalar_raincv_surface,      "raincv")
    call wrap_add_field_1d(pstate%atm_in_shflx_at_pc_surface, "hfx")
    call wrap_add_field_1d(pstate%qfx_at_pc_surface,          "qfx")
! rad flux
    call wrap_add_field_1d(pstate%atm_out_flwds_at_pc_surface,"glw")
    call wrap_add_field_1d(pstate%atm_out_netsw_at_pc_surface,"gsw")
    call wrap_add_field_1d(pstate%atm_out_fswds_at_pc_surface,"swdown")
! surface layer+pbl
    call wrap_add_field_1d(pstate%ustar_at_pc_surface,        "ust")
    call wrap_add_field_1d(pstate%znt_at_pc_surface,          "znt")
    call wrap_add_field_1d(pstate%mol_at_pc_surface,          "mol")
    call wrap_add_field_1d(pstate%pblh_at_pc_surface,         "pblh")
    call wrap_add_field_1d(pstate%mavail_at_pc_surface,       "mavail")
! albedo
    call wrap_add_field_1d(pstate%snowhland_at_pc_surface,    "snowh")
    call wrap_add_field_1d(pstate%atm_in_asdir_at_pc_surface, "asdir")
    call wrap_add_field_1d(pstate%atm_in_asdif_at_pc_surface, "asdif")
    call wrap_add_field_1d(pstate%atm_in_aldir_at_pc_surface, "aldir")
    call wrap_add_field_1d(pstate%atm_in_aldif_at_pc_surface, "aldif")

    call wrap_add_field_1d(output1d_temp1, "lwupt")
    call wrap_add_field_1d(output1d_temp2, "lwupb")
    call wrap_add_field_1d(output1d_temp3, "lwdnb")
    call wrap_add_field_1d(output1d_temp4, "lwuptc")
    call wrap_add_field_1d(output1d_temp5, "lwupbc")
    call wrap_add_field_1d(output1d_temp6, "lwdnbc")
    call wrap_add_field_1d(output1d_temp7, "lwcf")
    call wrap_add_field_1d(output1d_temp8, "swupt")
    call wrap_add_field_1d(output1d_temp9, "swdnt")
    call wrap_add_field_1d(output1d_temp10,"swupb")
    call wrap_add_field_1d(output1d_temp11,"swdnb")
    call wrap_add_field_1d(output1d_temp12,"swuptc")
    call wrap_add_field_1d(output1d_temp13,"swdntc")
    call wrap_add_field_1d(output1d_temp14,"swupbc")
    call wrap_add_field_1d(output1d_temp15,"swdnbc")
    call wrap_add_field_1d(output1d_temp16,"swcf")

    call wrap_output_1d_group(mesh, outdir,fname_restart1)
    call wrap_output_clean_1d()
! deallocate 1d
    call wrap_deallocate_data1d(output1d_temp1)
    call wrap_deallocate_data1d(output1d_temp2)
    call wrap_deallocate_data1d(output1d_temp3)
    call wrap_deallocate_data1d(output1d_temp4)
    call wrap_deallocate_data1d(output1d_temp5)
    call wrap_deallocate_data1d(output1d_temp6)
    call wrap_deallocate_data1d(output1d_temp7)
    call wrap_deallocate_data1d(output1d_temp8)
    call wrap_deallocate_data1d(output1d_temp9)
    call wrap_deallocate_data1d(output1d_temp10)
    call wrap_deallocate_data1d(output1d_temp11)
    call wrap_deallocate_data1d(output1d_temp12)
    call wrap_deallocate_data1d(output1d_temp13)
    call wrap_deallocate_data1d(output1d_temp14)
    call wrap_deallocate_data1d(output1d_temp15)
    call wrap_deallocate_data1d(output1d_temp16)
!
! write slow physics tendency, currently CU, RAD, BL
!
    if(step_cu.gt.1)then ! if cu is slow-phys
       call wrap_allocate_data2d(mesh,nlev,output2d_temp1)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp2)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp3)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp4)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp5)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp6)
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       call wrap_allocate_data2d(mesh,nlev,output2d_temp6)
    end if
    
    ! rad always slow
    call wrap_allocate_data2d(mesh,nlev,output2d_temp7)
    if(step_bl.gt.1)then ! if bl is slow-phys
       call wrap_allocate_data2d(mesh,nlev,output2d_temp8)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp9)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp10)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp11)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp12)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp13)
    end if

    if(step_cu.gt.1)then ! if cu is slow-phys
       output2d_temp1%f(1:nlev, 1:ncell) = transpose(ptend_wrf%ruucuten(1:ncell,1:nlev,1))
       output2d_temp2%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rvvcuten(1:ncell,1:nlev,1))
       output2d_temp3%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rthcuten(1:ncell,1:nlev,1))
       output2d_temp4%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rqvcuten(1:ncell,1:nlev,1))
       output2d_temp5%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rqccuten(1:ncell,1:nlev,1))
       output2d_temp6%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       output2d_temp6%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rqicuten(1:ncell,1:nlev,1))
    end if
    ! rad is definitely slow-physics for any scale, so hard-coded here
    output2d_temp7%f(1:nlev, 1:ncell)     = transpose(ptend_wrf%rthraten(1:ncell,1:nlev,1))
    if(step_bl.gt.1)then ! if bl is slow-phys
       output2d_temp8%f (1:nlev, 1:ncell) = transpose(ptend_wrf%ruublten(1:ncell,1:nlev,1))
       output2d_temp9%f (1:nlev, 1:ncell) = transpose(ptend_wrf%rvvblten(1:ncell,1:nlev,1))
       output2d_temp10%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rthblten(1:ncell,1:nlev,1))
       output2d_temp11%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rqvblten(1:ncell,1:nlev,1))
       output2d_temp12%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rqcblten(1:ncell,1:nlev,1))
       output2d_temp13%f(1:nlev, 1:ncell) = transpose(ptend_wrf%rqiblten(1:ncell,1:nlev,1))
    end if

    call wrap_output_init_2d(mesh)
    if(step_cu.gt.1)then ! if cu is slow-phys
       call wrap_add_field_2d(output2d_temp1          ,"ruucuten")
       call wrap_add_field_2d(output2d_temp2          ,"rvvcuten")
       call wrap_add_field_2d(output2d_temp3          ,"rthcuten")
       call wrap_add_field_2d(output2d_temp4          ,"rqvcuten")
       call wrap_add_field_2d(output2d_temp5          ,"rqccuten")
       call wrap_add_field_2d(output2d_temp6          ,"rqicuten")
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       call wrap_add_field_2d(output2d_temp6          ,"rqicuten")
    end if
    ! rad is definitely slow-physics for any scale, so hard-coded here
    call wrap_add_field_2d(output2d_temp7             ,"rthraten")
    if(step_bl.gt.1)then ! if bl is slow-phys
       call wrap_add_field_2d(output2d_temp8          ,"ruublten")
       call wrap_add_field_2d(output2d_temp9          ,"rvvblten")
       call wrap_add_field_2d(output2d_temp10         ,"rthblten")
       call wrap_add_field_2d(output2d_temp11         ,"rqvblten")
       call wrap_add_field_2d(output2d_temp12         ,"rqcblten")
       call wrap_add_field_2d(output2d_temp13         ,"rqiblten")
    end if
    call wrap_output_2d_group(mesh, outdir,fname_restart2)
    call wrap_output_clean_2d()

!
! deallocate 2d
!
    if(step_cu.gt.1)then ! if cu is slow-phys
       call wrap_deallocate_data2d(output2d_temp1)
       call wrap_deallocate_data2d(output2d_temp2)
       call wrap_deallocate_data2d(output2d_temp3)
       call wrap_deallocate_data2d(output2d_temp4)
       call wrap_deallocate_data2d(output2d_temp5)
       call wrap_deallocate_data2d(output2d_temp6)
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       call wrap_deallocate_data2d(output2d_temp6)
    end if
    call wrap_deallocate_data2d(output2d_temp7)
    if(step_bl.gt.1)then ! if bl is slow-phys
       call wrap_deallocate_data2d(output2d_temp8)
       call wrap_deallocate_data2d(output2d_temp9)
       call wrap_deallocate_data2d(output2d_temp10)
       call wrap_deallocate_data2d(output2d_temp11)
       call wrap_deallocate_data2d(output2d_temp12)
       call wrap_deallocate_data2d(output2d_temp13)
    end if

!================================================
! 3D var for restart, only when nspecies>ntracer
!================================================
    if(nspecies.gt.ntracer)then
       fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
       ndim3 = nspecies-ntracer
       call wrap_allocate_data3d(mesh,nlev,ndim3,output3d_temp1)
       do itracer = 1, ndim3
          output3d_temp1%f(itracer,1:nlev,1:ncell) = transpose(pstate_wrf%moist(1:ncell,1:nlev,1,ntracer+itracer))
       end do

       call wrap_output_init_3d(mesh)
       call wrap_add_field_3d(output3d_temp1,"PhysStateMoistSpecies")
       call wrap_output_3d_group(mesh, outdir, fname_restart3, ndim3)
       call wrap_output_clean_3d
       call wrap_deallocate_data3d(output3d_temp1)
    end if

    return
  end subroutine gcm_write_restart_file_phy2

  subroutine gcm_read_restart_file_phy2(mesh,cyear,cmon,cday)
! io
    type(global_domain), intent(inout), target :: mesh
    character(len=4)          :: cyear
    character(len=2)          :: cmon
    character(len=2)          :: cday
! local
    character(128)            :: c_glevel
    character(128)            :: fname_restart1,fname_restart2,fname_restart3
    character(len=5)          :: day
    character(len=5)          :: sec
    !integer                   :: i,j
    logical                   :: does_file_exist
    integer(i8)               :: lev,levp, ncell, ndim3, itracer
    real(r8)                  :: timestep
    type(scalar_3d_field) :: output3d_temp1
    type(scalar_2d_field) :: output2d_temp1, output2d_temp2,  output2d_temp3,  output2d_temp4,  output2d_temp5, output2d_temp6, output2d_temp7
    type(scalar_2d_field) :: output2d_temp8, output2d_temp9,  output2d_temp10, output2d_temp11, output2d_temp12,output2d_temp13
    type(scalar_1d_field) :: output1d_temp1, output1d_temp2,  output1d_temp3,  output1d_temp4,  output1d_temp5, output1d_temp6, output1d_temp7, output1d_temp8
    type(scalar_1d_field) :: output1d_temp9, output1d_temp10, output1d_temp11, output1d_temp12, output1d_temp13,output1d_temp14,output1d_temp15,output1d_temp16

!================================================
! Construct Filename
!================================================

! read the timestep at which we generate this restart file
    lev  = nlev
    levp = nlevp
    ncell= mesh%nv_halo(1)
    call write_string(mesh%glevel     , c_glevel)

!================================================
!  if file exist, read the restart_timestep
!================================================
    inquire(file=trim(outdir)//'step_restart.txt',exist=does_file_exist)
    if(does_file_exist) then
       open(1,file=trim(outdir)//'step_restart.txt',status='old')
       read(1,*) restart_timestep
       close(1)

       open(2,file=trim(outdir)//'step_stop.txt',status='old')
       read(2,*) itimestep_pre
       close(2)
    else
       restart_timestep = 1
    endif  !if(does_file_exist)
    
    SELECT CASE(trim(working_mode))
    CASE('amipw')
      timestep = model_timestep
    END SELECT

    call get_current_date(restart_timestep,timestep,day,sec)

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_restart1="GRIST.RST.Phy."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2="GRIST.RST.Phy."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"

    call wrap_allocate_data1d(mesh,output1d_temp1)
    call wrap_allocate_data1d(mesh,output1d_temp2)
    call wrap_allocate_data1d(mesh,output1d_temp3)
    call wrap_allocate_data1d(mesh,output1d_temp4)
    call wrap_allocate_data1d(mesh,output1d_temp5)
    call wrap_allocate_data1d(mesh,output1d_temp6)
    call wrap_allocate_data1d(mesh,output1d_temp7)
    call wrap_allocate_data1d(mesh,output1d_temp8)
    call wrap_allocate_data1d(mesh,output1d_temp9)
    call wrap_allocate_data1d(mesh,output1d_temp10)
    call wrap_allocate_data1d(mesh,output1d_temp11)
    call wrap_allocate_data1d(mesh,output1d_temp12)
    call wrap_allocate_data1d(mesh,output1d_temp13)
    call wrap_allocate_data1d(mesh,output1d_temp14)
    call wrap_allocate_data1d(mesh,output1d_temp15)
    call wrap_allocate_data1d(mesh,output1d_temp16)

    if(does_file_exist) then
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'ts',    0, pstate%ts_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'rainc', 0, pstate%scalar_rainc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'rainnc',0, pstate%scalar_rainnc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'snownc',0, pstate%scalar_snownc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'grapnc',0, pstate%scalar_grapnc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'raincv',0, pstate%scalar_raincv_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'hfx',   0, pstate%atm_in_shflx_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'qfx',   0, pstate%qfx_at_pc_surface%f)
! rad flux
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'glw',   0, pstate%atm_out_flwds_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'gsw',   0, pstate%atm_out_netsw_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swdown',0, pstate%atm_out_fswds_at_pc_surface%f)
! surface layer+pbl
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'ust',   0, pstate%ustar_at_pc_surface%f  )
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'znt',   0, pstate%znt_at_pc_surface%f    )
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'mol',   0, pstate%mol_at_pc_surface%f    )
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'pblh',  0, pstate%pblh_at_pc_surface%f   )
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'mavail',0, pstate%mavail_at_pc_surface%f )
! albedo
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'snowh', 0, pstate%snowhland_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'asdir', 0, pstate%atm_in_asdir_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'asdif', 0, pstate%atm_in_asdif_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'aldir', 0, pstate%atm_in_aldir_at_pc_surface%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'aldif', 0, pstate%atm_in_aldif_at_pc_surface%f)
! rad flux diag
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwupt' , 0, output1d_temp1%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwupb' , 0, output1d_temp2%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwdnb' , 0, output1d_temp3%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwuptc', 0, output1d_temp4%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwupbc', 0, output1d_temp5%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwdnbc', 0, output1d_temp6%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwcf ' , 0, output1d_temp7%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swupt' , 0, output1d_temp8%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swdnt' , 0, output1d_temp9%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swupb' , 0, output1d_temp10%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swdnb' , 0, output1d_temp11%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swuptc', 0, output1d_temp12%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swdntc', 0, output1d_temp13%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swupbc', 0, output1d_temp14%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swdnbc', 0, output1d_temp15%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'swcf'  , 0, output1d_temp16%f)
    endif
!
! Fill physics internal state
!
    psurf_wrf%tsk (1:ncell,1)   = pstate%ts_at_pc_surface%f(1:ncell)
    psurf_wrf%rainc (1:ncell,1) = pstate%scalar_rainc_surface%f(1:ncell)
    psurf_wrf%rainnc(1:ncell,1) = pstate%scalar_rainnc_surface%f(1:ncell)
    psurf_wrf%snownc(1:ncell,1) = pstate%scalar_snownc_surface%f(1:ncell)
    psurf_wrf%graupelnc(1:ncell,1) = pstate%scalar_grapnc_surface%f(1:ncell)
    psurf_wrf%raincv(1:ncell,1) = pstate%scalar_raincv_surface%f(1:ncell)
    psurf_wrf%hfx (1:ncell,1)   = pstate%atm_in_shflx_at_pc_surface%f(1:ncell)
    psurf_wrf%qfx (1:ncell,1)   = pstate%qfx_at_pc_surface%f (1:ncell)
! rad flux
    psurf_wrf%glw (1:ncell,1)   = pstate%atm_out_flwds_at_pc_surface%f(1:ncell)
    psurf_wrf%gsw (1:ncell,1)   = pstate%atm_out_netsw_at_pc_surface%f(1:ncell)
    psurf_wrf%swdown(1:ncell,1) = pstate%atm_out_fswds_at_pc_surface%f(1:ncell)
! surface layer+pbl
    psurf_wrf%ust(1:ncell,1)    = pstate%ustar_at_pc_surface%f(1:ncell)
    psurf_wrf%znt(1:ncell,1)    = pstate%znt_at_pc_surface%f(1:ncell)
    psurf_wrf%mol(1:ncell,1)    = pstate%mol_at_pc_surface%f(1:ncell)
    psurf_wrf%pblh(1:ncell,1)   = pstate%pblh_at_pc_surface%f(1:ncell)
    psurf_wrf%mavail(1:ncell,1) = pstate%mavail_at_pc_surface%f(1:ncell)
! albedo
    psurf_wrf%snow(1:ncell,1)   = pstate%snowhland_at_pc_surface%f(1:ncell) ! m
    psurf_wrf%asdir(1:ncell,1)  = pstate%atm_in_asdir_at_pc_surface%f(1:ncell)
    psurf_wrf%asdif(1:ncell,1)  = pstate%atm_in_asdif_at_pc_surface%f(1:ncell)
    psurf_wrf%aldir(1:ncell,1)  = pstate%atm_in_aldir_at_pc_surface%f(1:ncell)
    psurf_wrf%aldif(1:ncell,1)  = pstate%atm_in_aldif_at_pc_surface%f(1:ncell)
! rad flux for diag
    pstate_wrf%lwupt(1:ncell,1) = output1d_temp1%f(1:ncell) 
    pstate_wrf%lwupb(1:ncell,1) = output1d_temp2%f(1:ncell) 
    pstate_wrf%lwdnb(1:ncell,1) = output1d_temp3%f(1:ncell) 
    pstate_wrf%lwuptc(1:ncell,1)= output1d_temp4%f(1:ncell) 
    pstate_wrf%lwupbc(1:ncell,1)= output1d_temp5%f(1:ncell) 
    pstate_wrf%lwdnbc(1:ncell,1)= output1d_temp6%f(1:ncell) 
    pstate_wrf%lwcf(1:ncell,1)  = output1d_temp7%f(1:ncell) 
    pstate_wrf%swupt(1:ncell,1) = output1d_temp8%f(1:ncell) 
    pstate_wrf%swdnt(1:ncell,1) = output1d_temp9%f(1:ncell) 
    pstate_wrf%swupb(1:ncell,1) = output1d_temp10%f(1:ncell)
    pstate_wrf%swdnb(1:ncell,1) = output1d_temp11%f(1:ncell)
    pstate_wrf%swuptc(1:ncell,1)= output1d_temp12%f(1:ncell)
    pstate_wrf%swdntc(1:ncell,1)= output1d_temp13%f(1:ncell)
    pstate_wrf%swupbc(1:ncell,1)= output1d_temp14%f(1:ncell)
    pstate_wrf%swdnbc(1:ncell,1)= output1d_temp15%f(1:ncell)
    pstate_wrf%swcf(1:ncell,1)  = output1d_temp16%f(1:ncell)

!================================================
!2D: read slow physics tendency, currently CU,RAD
!================================================

    if(step_cu.gt.1)then
       call wrap_allocate_data2d(mesh,nlev,output2d_temp1)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp2)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp3)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp4)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp5)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp6)
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       call wrap_allocate_data2d(mesh,nlev,output2d_temp6)
    end if
    call wrap_allocate_data2d(mesh,nlev,output2d_temp7)
    if(step_bl.gt.1)then
       call wrap_allocate_data2d(mesh,nlev,output2d_temp8)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp9)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp10)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp11)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp12)
       call wrap_allocate_data2d(mesh,nlev,output2d_temp13)
    end if

    if(step_cu.gt.1)then
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ruucuten' , lev, 0, output2d_temp1%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rvvcuten' , lev, 0, output2d_temp2%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rthcuten' , lev, 0, output2d_temp3%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rqvcuten' , lev, 0, output2d_temp4%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rqccuten' , lev, 0, output2d_temp5%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rqicuten' , lev, 0, output2d_temp6%f)
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rqicuten' , lev, 0, output2d_temp6%f)
    end if
    call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rthraten' , lev, 0, output2d_temp7%f)
    if(step_bl.gt.1)then
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'ruublten' , lev, 0, output2d_temp8%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rvvblten' , lev, 0, output2d_temp9%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rthblten' , lev, 0, output2d_temp10%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rqvblten' , lev, 0, output2d_temp11%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rqcblten' , lev, 0, output2d_temp12%f)
       call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'rqiblten' , lev, 0, output2d_temp13%f)
    end if

    if(step_cu.gt.1)then
       ptend_wrf%ruucuten(1:ncell,1:nlev,1) = transpose(output2d_temp1%f(1:nlev, 1:ncell))
       ptend_wrf%rvvcuten(1:ncell,1:nlev,1) = transpose(output2d_temp2%f(1:nlev, 1:ncell))
       ptend_wrf%rthcuten(1:ncell,1:nlev,1) = transpose(output2d_temp3%f(1:nlev, 1:ncell))
       ptend_wrf%rqvcuten(1:ncell,1:nlev,1) = transpose(output2d_temp4%f(1:nlev, 1:ncell))
       ptend_wrf%rqccuten(1:ncell,1:nlev,1) = transpose(output2d_temp5%f(1:nlev, 1:ncell))
       ptend_wrf%rqicuten(1:ncell,1:nlev,1) = transpose(output2d_temp6%f(1:nlev, 1:ncell))
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       ptend_wrf%rqicuten(1:ncell,1:nlev,1) = transpose(output2d_temp6%f(1:nlev, 1:ncell))
    end if
    ptend_wrf%rthraten(1:ncell,1:nlev,1) = transpose(output2d_temp7%f(1:nlev, 1:ncell))
    if(step_bl.gt.1)then
       ptend_wrf%ruublten(1:ncell,1:nlev,1) = transpose(output2d_temp8%f(1:nlev, 1:ncell))
       ptend_wrf%rvvblten(1:ncell,1:nlev,1) = transpose(output2d_temp9%f(1:nlev, 1:ncell))
       ptend_wrf%rthblten(1:ncell,1:nlev,1) = transpose(output2d_temp10%f(1:nlev, 1:ncell))
       ptend_wrf%rqvblten(1:ncell,1:nlev,1) = transpose(output2d_temp11%f(1:nlev, 1:ncell))
       ptend_wrf%rqcblten(1:ncell,1:nlev,1) = transpose(output2d_temp12%f(1:nlev, 1:ncell))
       ptend_wrf%rqiblten(1:ncell,1:nlev,1) = transpose(output2d_temp13%f(1:nlev, 1:ncell))
    end if

    if(step_cu.gt.1)then   ! if cu is slow-phys
       call wrap_deallocate_data2d(output2d_temp1)
       call wrap_deallocate_data2d(output2d_temp2)
       call wrap_deallocate_data2d(output2d_temp3)
       call wrap_deallocate_data2d(output2d_temp4)
       call wrap_deallocate_data2d(output2d_temp5)
       call wrap_deallocate_data2d(output2d_temp6)
    end if
    if(step_cu.eq.1.and.trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381')then
       call wrap_deallocate_data2d(output2d_temp6)
    end if
    call wrap_deallocate_data2d(output2d_temp7)
    if(step_bl.gt.1)then   ! if bl is slow-phys
       call wrap_deallocate_data2d(output2d_temp8)
       call wrap_deallocate_data2d(output2d_temp9)
       call wrap_deallocate_data2d(output2d_temp10)
       call wrap_deallocate_data2d(output2d_temp11)
       call wrap_deallocate_data2d(output2d_temp12)
       call wrap_deallocate_data2d(output2d_temp13)
    end if

!================================================
! 3D var for read, only when nspecies>ntracer
!================================================
    if(nspecies.gt.ntracer)then
       fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
       ndim3 = nspecies-ntracer
       call wrap_allocate_data3d(mesh,nlev,int(ndim3,i4),output3d_temp1)
       call wrap_read_3d_group_rst(mesh%gcomm_read, outdir, fname_restart3, 'PhysStateMoistSpecies', lev, 0, output3d_temp1%f, ndim3)
       do itracer = 1, int(ndim3,i4)
! index <=ntracer are handled by dynamics
          pstate_wrf%moist(1:ncell,1:nlev,1,ntracer+itracer) = transpose(output3d_temp1%f(itracer,1:nlev,1:ncell))
       end do
       call wrap_deallocate_data3d(output3d_temp1)
    end if
    return
  end subroutine gcm_read_restart_file_phy2
#endif

#ifdef USE_NOAHMP

  subroutine gcm_write_restart_file_lnd(mesh,itimestep,cyear,cmon,cday)
    use grist_lsm_noahmp_resVars, only: pstate_lsm, fill_lsm_resVars_for_write
    use grist_prescribe_seaice_module, only: sicetemp_at_pc_surface
! io
    type(global_domain), intent(in)  :: mesh
    integer(i4),         intent(in)  :: itimestep
    character(len=4)                 :: cyear
    character(len=2)                 :: cmon
    character(len=2)                 :: cday
! local
    character(128)     :: c_glevel
    character(128)     :: fname_restart1, fname_restart2
    character(len=5)   :: day
    character(len=5)   :: sec
    logical            :: timecheck
    real(r8)           :: timestep

!================================================
! [1] Construct Filename
!================================================

    timestep = model_timestep

    call write_string(mesh%glevel     , c_glevel)
    call get_current_date(itimestep,timestep,day,sec)

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_restart1 = "GRIST.RST.Lnd."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2 = "GRIST.RST.Lnd."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"

!================================================
! [2] Write Restart File
!================================================		
    call fill_lsm_resVars_for_write(mesh%nv_halo(1))

    call wrap_output_init_1d(mesh)
    call wrap_add_field_1d(pstate_lsm%tsk,       "tsk")
    call wrap_add_field_1d(pstate_lsm%hfx,       "hfx")
    call wrap_add_field_1d(pstate_lsm%qfx,       "qfx")
    call wrap_add_field_1d(pstate_lsm%lh,        "lh" )
    call wrap_add_field_1d(pstate_lsm%grdflx,    "grdflx")
    call wrap_add_field_1d(pstate_lsm%smstav,    "smstav")
    call wrap_add_field_1d(pstate_lsm%smstot,    "smstot")
    call wrap_add_field_1d(pstate_lsm%sfcrunoff, "sfcrunoff")
    call wrap_add_field_1d(pstate_lsm%udrunoff,  "udrunoff")
    call wrap_add_field_1d(pstate_lsm%albedo,    "albedo")
    call wrap_add_field_1d(pstate_lsm%snowc,     "snowc")
    call wrap_add_field_1d(pstate_lsm%snow,      "snow")
    call wrap_add_field_1d(pstate_lsm%snowh,     "snowh")
    call wrap_add_field_1d(pstate_lsm%canwat,    "canwat")
    call wrap_add_field_1d(pstate_lsm%acsnom,    "acsnom")
    call wrap_add_field_1d(pstate_lsm%acsnow,    "acsnow")
    call wrap_add_field_1d(pstate_lsm%emiss,     "emiss")
    call wrap_add_field_1d(pstate_lsm%qsfc,      "qsfc")
    call wrap_add_field_1d(pstate_lsm%z0,        "z0")
    call wrap_add_field_1d(pstate_lsm%znt,       "znt")
    call wrap_add_field_1d(pstate_lsm%isnowxy,   "isnowxy")
    call wrap_add_field_1d(pstate_lsm%tvxy,      "tvxy")
    call wrap_add_field_1d(pstate_lsm%tgxy,      "tgxy")
    call wrap_add_field_1d(pstate_lsm%canicexy,  "canicexy")
    call wrap_add_field_1d(pstate_lsm%canliqxy,  "canliqxy")
    call wrap_add_field_1d(pstate_lsm%eahxy,     "eahxy")
    call wrap_add_field_1d(pstate_lsm%tahxy,     "tahxy")
    call wrap_add_field_1d(pstate_lsm%cmxy,      "cmxy")
    call wrap_add_field_1d(pstate_lsm%chxy,      "chxy")
    call wrap_add_field_1d(pstate_lsm%fwetxy,    "fwetxy")
    call wrap_add_field_1d(pstate_lsm%sneqvoxy,  "sneqvoxy")
    call wrap_add_field_1d(pstate_lsm%alboldxy,  "alboldxy")
    call wrap_add_field_1d(pstate_lsm%qsnowxy,   "qsnowxy")
    call wrap_add_field_1d(pstate_lsm%wslakexy,  "wslakexy")
    call wrap_add_field_1d(pstate_lsm%zwtxy,     "zwtxy")
    call wrap_add_field_1d(pstate_lsm%waxy,      "waxy")
    call wrap_add_field_1d(pstate_lsm%wtxy,      "wtxy")
    call wrap_add_field_1d(pstate_lsm%lfmassxy,  "lfmassxy")
    call wrap_add_field_1d(pstate_lsm%rtmassxy,  "rtmassxy")
    call wrap_add_field_1d(pstate_lsm%stmassxy,  "stmassxy")
    call wrap_add_field_1d(pstate_lsm%woodxy,    "woodxy")
    call wrap_add_field_1d(pstate_lsm%grainxy,   "grainxy")
    call wrap_add_field_1d(pstate_lsm%gddxy,     "gddxy")
    call wrap_add_field_1d(pstate_lsm%stblcpxy,  "stblcpxy")
    call wrap_add_field_1d(pstate_lsm%fastcpxy,  "fastcpxy")
    call wrap_add_field_1d(pstate_lsm%xlaixy,    "xlaixy")
    call wrap_add_field_1d(pstate_lsm%xsaixy,    "xsaixy")
    call wrap_add_field_1d(pstate_lsm%taussxy,   "taussxy")
    call wrap_add_field_1d(pstate_lsm%smcwtdxy,  "smcwtdxy")
    call wrap_add_field_1d(pstate_lsm%deeprechxy,"deeprechxy")
    call wrap_add_field_1d(pstate_lsm%rechxy,    "rechzy")
! write seaice temp here !!
    call wrap_add_field_1d(sicetemp_at_pc_surface,"sicetemp")
    call wrap_output_1d_group(mesh, outdir,fname_restart1)
    call wrap_output_clean_1d()

! 2d
    call wrap_output_init_2d(mesh)
    call wrap_add_field_2d(pstate_lsm%tslb,       "tslb")
    call wrap_add_field_2d(pstate_lsm%smois,      "smois")
    call wrap_add_field_2d(pstate_lsm%smoiseq,    "smoiseq")
    call wrap_add_field_2d(pstate_lsm%sh2o,       "sh2o")
    call wrap_add_field_2d(pstate_lsm%tsnoxy,     "tsnoxy")
    call wrap_add_field_2d(pstate_lsm%zsnsoxy,    "zsnsoxy")
    call wrap_add_field_2d(pstate_lsm%snicexy,    "snicexy")
    call wrap_add_field_2d(pstate_lsm%snliqxy,    "snliqxy")
    call wrap_output_2d_group(mesh, outdir,fname_restart2)
    call wrap_output_clean_2d()

    return
  end subroutine gcm_write_restart_file_lnd

  subroutine gcm_read_restart_file_lnd(mesh,cyear,cmon,cday)
    use grist_lsm_noahmp_resVars, only: pstate_lsm, fill_lsm_resVars_after_read
    use grist_prescribe_seaice_module, only: sicetemp_at_pc_surface
! io
    type(global_domain), intent(inout), target :: mesh
    character(len=4)          :: cyear
    character(len=2)          :: cmon
    character(len=2)          :: cday
! local
    character(128)            :: c_glevel
    character(128)            :: fname_restart1, fname_restart2
    character(len=5)          :: day
    character(len=5)          :: sec
    logical                   :: does_file_exist
    real(r8)                  :: timestep

!================================================
! Construct Filename
!================================================

! read the timestep at which we generate this restart file
    call write_string(mesh%glevel     , c_glevel)

!================================================
!  if file exist, read the restart_timestep
!================================================
    inquire(file=trim(outdir)//'step_restart.txt',exist=does_file_exist)
    if(does_file_exist) then
       open(1,file=trim(outdir)//'step_restart.txt',status='old')
       read(1,*) restart_timestep
       close(1)

       open(2,file=trim(outdir)//'step_stop.txt',status='old')
       read(2,*) itimestep_pre
       close(2)
    else
       restart_timestep = 1
    endif
    
    timestep = model_timestep

    call get_current_date(restart_timestep,timestep,day,sec)

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_restart1="GRIST.RST.Lnd."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2="GRIST.RST.Lnd."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"

    if(does_file_exist) then
! 1d
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "tsk",      0,pstate_lsm%tsk%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "hfx",      0,pstate_lsm%hfx%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "qfx",      0,pstate_lsm%qfx%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "lh" ,      0,pstate_lsm%lh%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "grdflx",   0,pstate_lsm%grdflx%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "smstav",   0,pstate_lsm%smstav%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "smstot",   0,pstate_lsm%smstot%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "sfcrunoff",0,pstate_lsm%sfcrunoff%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "udrunoff" ,0,pstate_lsm%udrunoff%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "albedo"   ,0,pstate_lsm%albedo%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "snowc"    ,0,pstate_lsm%snowc%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "snow"     ,0,pstate_lsm%snow%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "snowh"    ,0,pstate_lsm%snowh%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "canwat"   ,0,pstate_lsm%canwat%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "acsnom"   ,0,pstate_lsm%acsnom%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "acsnow"   ,0,pstate_lsm%acsnow%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "emiss"    ,0,pstate_lsm%emiss%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "qsfc"     ,0,pstate_lsm%qsfc%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "z0"       ,0,pstate_lsm%z0%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "znt"      ,0,pstate_lsm%znt%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "isnowxy"  ,0,pstate_lsm%isnowxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "tvxy"     ,0,pstate_lsm%tvxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "tgxy"     ,0,pstate_lsm%tgxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "canicexy" ,0,pstate_lsm%canicexy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "canliqxy" ,0,pstate_lsm%canliqxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "eahxy"    ,0,pstate_lsm%eahxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "tahxy"    ,0,pstate_lsm%tahxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "cmxy"     ,0,pstate_lsm%cmxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "chxy"     ,0,pstate_lsm%chxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "fwetxy"   ,0,pstate_lsm%fwetxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "sneqvoxy" ,0,pstate_lsm%sneqvoxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "alboldxy" ,0,pstate_lsm%alboldxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "qsnowxy"  ,0,pstate_lsm%qsnowxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "wslakexy" ,0,pstate_lsm%wslakexy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "zwtxy"    ,0,pstate_lsm%zwtxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "waxy"     ,0,pstate_lsm%waxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "wtxy"     ,0,pstate_lsm%wtxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "lfmassxy" ,0,pstate_lsm%lfmassxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "rtmassxy" ,0,pstate_lsm%rtmassxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "stmassxy" ,0,pstate_lsm%stmassxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "woodxy"   ,0,pstate_lsm%woodxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "grainxy"  ,0,pstate_lsm%grainxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "gddxy"    ,0,pstate_lsm%gddxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "stblcpxy" ,0,pstate_lsm%stblcpxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "fastcpxy" ,0,pstate_lsm%fastcpxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "xlaixy"   ,0,pstate_lsm%xlaixy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "xsaixy"   ,0,pstate_lsm%xsaixy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "taussxy"  ,0,pstate_lsm%taussxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "smcwtdxy" ,0,pstate_lsm%smcwtdxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "deeprechxy",0,pstate_lsm%deeprechxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "rechzy",   0,pstate_lsm%rechxy%f)
      call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, "sicetemp", 0,sicetemp_at_pc_surface%f)
! 2d
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'tslb',   int8(4), 0, pstate_lsm%tslb%f)
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'smois',  int8(4), 0, pstate_lsm%smois%f)
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'smoiseq',int8(4), 0, pstate_lsm%smoiseq%f)
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'sh2o',   int8(4), 0, pstate_lsm%sh2o%f)
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'tsnoxy', int8(3), 0, pstate_lsm%tsnoxy%f)
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'zsnsoxy',int8(7), 0, pstate_lsm%zsnsoxy%f)
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'snicexy',int8(3), 0, pstate_lsm%snicexy%f)
      call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'snliqxy',int8(3), 0, pstate_lsm%snliqxy%f)

      call fill_lsm_resVars_after_read(mesh%nv_halo(1))
#ifdef AMIPW_PHYSICS
! for sfclay-v381
      psurf_wrf%qsfc(1:mesh%nv_halo(1),1) = pstate_lsm%qsfc%f(1:mesh%nv_halo(1)) 
#endif

    endif
    return
  end subroutine gcm_read_restart_file_lnd
#endif

!================================================
! purely private
!================================================

  subroutine wrap_allocate_data1d(mesh,var)
     type(global_domain),   intent(in)    :: mesh
     type(scalar_1d_field), intent(inout) :: var
     if(.not.allocated(var%f)) allocate(var%f(mesh%nv_full))
! field and position
     var%f    = zero
     var%pos  = 0
     return
  end subroutine wrap_allocate_data1d

  subroutine wrap_allocate_data2d(mesh,nLevel,var)
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: nLevel
     type(scalar_2d_field), intent(inout) :: var
     if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%nv_full))
     var%f    = zero
     var%pos  = 0
     return
  end subroutine wrap_allocate_data2d

  subroutine wrap_allocate_data3d(mesh,nLevel,ndim3,var)
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: nLevel, ndim3
     type(scalar_3d_field), intent(inout) :: var
     if(.not.allocated(var%f)) allocate(var%f(ndim3,nLevel,mesh%nv_full))
     var%f    = zero
     var%pos  = 0
     return
  end subroutine wrap_allocate_data3d

  subroutine wrap_deallocate_data3d(var)
     type(scalar_3d_field), intent(inout) :: var
     if(allocated(var%f)) deallocate(var%f)
     return
  end subroutine wrap_deallocate_data3d

  subroutine wrap_deallocate_data2d(var)
     type(scalar_2d_field), intent(inout) :: var
     if(allocated(var%f)) deallocate(var%f)
     return
  end subroutine wrap_deallocate_data2d

  subroutine wrap_deallocate_data1d(var)
     type(scalar_1d_field), intent(inout) :: var
     if(allocated(var%f)) deallocate(var%f)
     return
  end subroutine wrap_deallocate_data1d

end module grist_gcm_restart_module
