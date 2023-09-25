
  module grist_gcm_restart_phys1_module

#ifdef AMIPC_PHYSICS
  use grist_lib
  use grist_domain_types, only: global_domain
  use grist_constants
  use grist_nml_module,   only: working_mode,     &
                                outdir,           &
                                nlev,             &
                                nlevp,            &
                                ntracer,          &
                                start_ymd, start_tod, model_timestep, grid_info
   use grist_time_manager,only: get_current_date
   use grist_physics_data_structure,  only: pstate
   use grist_PhysC_data_structure,    only: pstate_cam 
   use phys_control,                  only: phys_getopts
   use grist_physics_update,          only: old_time_level
   use grist_rad_constituents,        only: rad_cnst_get_info

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
 
   use grist_util_module,    only: write_string
   use grist_data_types,     only: scalar_3d_field, scalar_2d_field, scalar_1d_field

   implicit none

   private

   public        :: gcm_write_restart_file_phy1,    &
                    gcm_read_restart_file_phy1

contains

    subroutine gcm_write_restart_file_phy1(mesh,itimestep,cyear,cmon,cday)
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
    character(len=16)  :: shallow_scheme, macrop_scheme

    integer(i4)        :: ncell, ndim3, nmodes, itracer, icell, ilev
    !logical            :: timecheck
    real(r8)           :: timestep
    type(scalar_2d_field) :: output2d_temp1, output2d_temp2,  output2d_temp3,  output2d_temp4,  output2d_temp5, output2d_temp6, output2d_temp7
    type(scalar_2d_field) :: output2d_temp8, output2d_temp9,  output2d_temp10, output2d_temp11, output2d_temp12,output2d_temp13,output2d_temp14
    type(scalar_2d_field) :: output2d_temp15, output2d_temp16
    type(scalar_1d_field) :: output1d_temp1
    type(scalar_3d_field) :: output3d_temp1

!================================================
! [1] Construct Filename
!================================================

    timestep = model_timestep

#ifndef NO_PHYSHALO
    ncell    = mesh%nv_halo(1)
#else
    ncell    = mesh%nv_compute
#endif

    call write_string(mesh%glevel     , c_glevel)
    call get_current_date(itimestep,timestep,day,sec)

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_restart1 = "GRIST.RST.Phy."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2 = "GRIST.RST.Phy."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
!
! write flux that affecets diagnose, because rad is not called each step, so
! restart will zero some flux states
!
!================================================
! [2] Write Restart File
!================================================		


    call phys_getopts(shallow_scheme_out = shallow_scheme,      &
                      macrop_scheme_out  = macrop_scheme )
! 1D
    call wrap_output_init_1d(mesh)
    call wrap_add_field_1d(pstate_cam%pbl_pblh_at_pc_surface,       "pblh"          )
    call wrap_add_field_1d(pstate_cam%pbl_tauresx_at_pc_surface,    "tauresx"       )
    call wrap_add_field_1d(pstate_cam%pbl_tauresy_at_pc_surface,    "tauresy"       )
    call wrap_add_field_1d(pstate_cam%pbl_tpert_at_pc_surface,      "tpert"         )
    call wrap_add_field_1d(pstate%snowhland_at_pc_surface,          "snowh"         )
    call wrap_add_field_1d(pstate%atm_in_asdir_at_pc_surface,       "asdir"         )
    call wrap_add_field_1d(pstate%atm_in_asdif_at_pc_surface,       "asdif"         )
    call wrap_add_field_1d(pstate%atm_in_aldir_at_pc_surface,       "aldir"         )
    call wrap_add_field_1d(pstate%atm_in_aldif_at_pc_surface,       "aldif"         )
    call wrap_add_field_1d(pstate%atm_in_lwup_at_pc_surface,        "lwup"          )
    if(trim(shallow_scheme) .eq. 'double_plume')then
    call wrap_add_field_1d(pstate_cam%pcape_previous,               "pcape_previous")
    call wrap_add_field_1d(pstate_cam%pcape_before_dycore,          "pcape_before_d")
    end if
    if(trim(macrop_scheme) .eq. 'lin')then
    call wrap_add_field_1d(pstate_cam%pbl_wstarPBL,                 "wstarPBL"      )
    end if
    call wrap_allocate_data1d(mesh,output1d_temp1)
    output1d_temp1%f(1:ncell) = pstate_cam%cumulus_cush%f(old_time_level,1:ncell)
    call wrap_add_field_1d(output1d_temp1, "cush")

    call wrap_output_1d_group(mesh, outdir,fname_restart1)
    call wrap_output_clean_1d()

    call wrap_deallocate_data1d(output1d_temp1)

! 2D
    call wrap_output_init_2d(mesh)
    call wrap_add_field_2d(pstate_cam%pbl_tke_at_pc_face_level,    "tke")
    call wrap_add_field_2d(pstate_cam%pbl_kvm_at_pc_face_level,    "kvm") 
    call wrap_add_field_2d(pstate_cam%pbl_kvh_at_pc_face_level,    "kvh")
    call wrap_add_field_2d(pstate_cam%sw_qrs_at_pc_full_level,     "qrs")
    call wrap_add_field_2d(pstate_cam%lw_qrl_at_pc_full_level,     "qrl")

    call wrap_allocate_data2d(mesh,nlev,output2d_temp1)
    call wrap_allocate_data2d(mesh,nlev,output2d_temp15)
    output2d_temp1%f(1:nlev,  1:ncell) = pstate_cam%microp_cldo_at_pc_full_level%f(old_time_level, 1:nlev, 1:ncell)
    output2d_temp15%f(1:nlev, 1:ncell) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level, 1:nlev, 1:ncell)

    call wrap_add_field_2d(output2d_temp1,      "cldo")
    call wrap_add_field_2d(output2d_temp15,     "cld")

    if(trim(macrop_scheme) .eq. 'park')then
        call wrap_allocate_data2d(mesh,nlev,output2d_temp2)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp3)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp4)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp5)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp6)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp7)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp8)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp9)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp10)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp11)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp12)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp13)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp14)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp16)

        output2d_temp2%f(1:nlev, 1:ncell)  = pstate_cam%macrop_qcwat_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp3%f(1:nlev, 1:ncell)  = pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp4%f(1:nlev, 1:ncell)  = pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp5%f(1:nlev, 1:ncell)  = pstate_cam%macrop_nlwat_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp6%f(1:nlev, 1:ncell)  = pstate_cam%macrop_niwat_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp7%f(1:nlev, 1:ncell)  = pstate_cam%macrop_tcwat_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp8%f(1:nlev, 1:ncell)  = pstate_cam%microp_cc_t_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp9%f(1:nlev, 1:ncell)  = pstate_cam%microp_cc_qv_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp10%f(1:nlev, 1:ncell) = pstate_cam%microp_cc_ql_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp11%f(1:nlev, 1:ncell) = pstate_cam%microp_cc_qi_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp12%f(1:nlev, 1:ncell) = pstate_cam%microp_cc_nl_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp13%f(1:nlev, 1:ncell) = pstate_cam%microp_cc_ni_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp14%f(1:nlev, 1:ncell) = pstate_cam%microp_cc_qlst_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)
        output2d_temp16%f(1:nlev, 1:ncell) = pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,1:nlev,1:ncell)

        call wrap_add_field_2d(output2d_temp2,       "qcwat"  )
        call wrap_add_field_2d(output2d_temp3,       "lcwat"  )
        call wrap_add_field_2d(output2d_temp4,       "iccwat" )
        call wrap_add_field_2d(output2d_temp5,       "nlwat"  )
        call wrap_add_field_2d(output2d_temp6,       "niwat"  )
        call wrap_add_field_2d(output2d_temp7,       "tcwat"  )
        call wrap_add_field_2d(output2d_temp8,       "cc_t"   )
        call wrap_add_field_2d(output2d_temp9,       "cc_qv"  )
        call wrap_add_field_2d(output2d_temp10,      "cc_ql"  )
        call wrap_add_field_2d(output2d_temp11,      "cc_qi"  )
        call wrap_add_field_2d(output2d_temp12,      "cc_nl"  )
        call wrap_add_field_2d(output2d_temp13,      "cc_ni"  )
        call wrap_add_field_2d(output2d_temp14,      "cc_qlst")
        call wrap_add_field_2d(output2d_temp16,      "concld" )
    end if

    if(trim(macrop_scheme) .eq. 'lin')then
        call wrap_add_field_2d(pstate_cam%pbl_lengi_at_pc_face_level,    "lengi")
        call wrap_add_field_2d(pstate_cam%pbl_shi_at_pc_face_level,      "shi"  )
    end if

    call wrap_output_2d_group(mesh, outdir,fname_restart2)
    call wrap_output_clean_2d()

    call wrap_deallocate_data2d(output2d_temp1)
    call wrap_deallocate_data2d(output2d_temp15)
    if(trim(macrop_scheme) .eq. 'park')then
        call wrap_deallocate_data2d(output2d_temp2)
        call wrap_deallocate_data2d(output2d_temp3)
        call wrap_deallocate_data2d(output2d_temp4)
        call wrap_deallocate_data2d(output2d_temp5)
        call wrap_deallocate_data2d(output2d_temp6)
        call wrap_deallocate_data2d(output2d_temp7)
        call wrap_deallocate_data2d(output2d_temp8)
        call wrap_deallocate_data2d(output2d_temp9)
        call wrap_deallocate_data2d(output2d_temp10)
        call wrap_deallocate_data2d(output2d_temp11)
        call wrap_deallocate_data2d(output2d_temp12)
        call wrap_deallocate_data2d(output2d_temp13)
        call wrap_deallocate_data2d(output2d_temp14)
        call wrap_deallocate_data2d(output2d_temp16)
    end if

! 3D
    fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".ghg.nc"
    ndim3 = pstate_cam%total_ghg_num
    call wrap_allocate_data3d(mesh,nlev,ndim3,output3d_temp1)
    do itracer = 1, ndim3
         output3d_temp1%f(itracer,1:nlev,1:ncell) = pstate_cam%ghg_at_pc_full_level(itracer)%f(1:nlev,1:ncell)
    end do
    call wrap_output_init_3d(mesh)
    call wrap_add_field_3d(output3d_temp1, "ghg")
    call wrap_output_3d_group(mesh, outdir, fname_restart3, ndim3)
    call wrap_output_clean_3d
    call wrap_deallocate_data3d(output3d_temp1)

    fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".aero.nc"
    ndim3 = pstate_cam%total_aerosol_num
    call wrap_allocate_data3d(mesh,nlev,ndim3,output3d_temp1)
    do itracer = 1, ndim3
        output3d_temp1%f(itracer,1:nlev,1:ncell) =  pstate_cam%aerosol_at_pc_full_level(itracer)%f(1:nlev,1:ncell)
    end do
    call wrap_output_init_3d(mesh)
    call wrap_add_field_3d(output3d_temp1, "aero")
    call wrap_output_3d_group(mesh, outdir, fname_restart3, ndim3)
    call wrap_output_clean_3d
    call wrap_deallocate_data3d(output3d_temp1)


    fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".dgnum.nc"
    call rad_cnst_get_info(0, nmodes=nmodes)
    ndim3 = nmodes*2
    call wrap_allocate_data3d(mesh,nlev,ndim3,output3d_temp1)
    do itracer = 1,nmodes 
        output3d_temp1%f(itracer,1:nlev,1:ncell) = pstate_cam%aerosol_dgnum%f(itracer,1:nlev,1:ncell)
        output3d_temp1%f(nmodes+itracer,1:nlev,1:ncell) = pstate_cam%aerosol_dgnumwet%f(itracer,1:nlev,1:ncell)
    end do
    call wrap_output_init_3d(mesh)
    call wrap_add_field_3d(output3d_temp1, "dgnum")
    call wrap_output_3d_group(mesh, outdir, fname_restart3, ndim3)
    call wrap_output_clean_3d
    call wrap_deallocate_data3d(output3d_temp1)

    if(mpi_rank() == 0) then
        open(1,file=trim(outdir)//'restart_oldtimelevel_4camphys.txt',status='unknown')
        write(1,*) old_time_level
        close(1)
    end if

    return
    end subroutine gcm_write_restart_file_phy1

    subroutine gcm_read_restart_file_phy1(mesh,cyear,cmon,cday)
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
    character(len=16)         :: shallow_scheme, macrop_scheme
    logical                   :: does_file_exist
    integer(i8)               :: lev,levp, ncell, ndim3, itracer
    real(r8)                  :: timestep
    integer(i4)               :: restart_timestep,nmodes
    type(scalar_3d_field) :: output3d_temp1
    type(scalar_1d_field) :: output1d_temp1
    type(scalar_2d_field) :: output2d_temp1, output2d_temp2,  output2d_temp3,  output2d_temp4,  output2d_temp5, output2d_temp6, output2d_temp7
    type(scalar_2d_field) :: output2d_temp8, output2d_temp9,  output2d_temp10, output2d_temp11, output2d_temp12,output2d_temp13,output2d_temp14
    type(scalar_2d_field) :: output2d_temp15,output2d_temp16
 
!================================================
! Construct Filename
!================================================

! read the timestep at which we generate this restart file
    lev  = nlev
    levp = nlevp
#ifndef NO_PHYSHALO
    ncell    = mesh%nv_halo(1)
#else
    ncell    = mesh%nv_compute
#endif

    call write_string(mesh%glevel     , c_glevel)

!================================================
!  if file exist, read the restart_timestep
!================================================
    inquire(file=trim(outdir)//'step_restart.txt',exist=does_file_exist)
    if(does_file_exist) then
       open(1,file=trim(outdir)//'step_restart.txt',status='old')
       read(1,*) restart_timestep
       close(1)

       open(2,file=trim(outdir)//'restart_oldtimelevel_4camphys.txt',status='old')
       read(2,*) old_time_level
       close(2)
    else
       restart_timestep = 1
    endif  !if(does_file_exist)
 
    timestep = model_timestep

    call get_current_date(restart_timestep,timestep,day,sec)
    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)
    fname_restart1="GRIST.RST.Phy."//trim(grid_info)//".1d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"
    fname_restart2="GRIST.RST.Phy."//trim(grid_info)//".2d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".nc"

    call phys_getopts(shallow_scheme_out = shallow_scheme,      &
                      macrop_scheme_out  = macrop_scheme )


! 1D
    if(does_file_exist) then
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'pblh',     0, pstate_cam%pbl_pblh_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'tauresx',  0, pstate_cam%pbl_tauresx_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'tauresy',  0, pstate_cam%pbl_tauresy_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'tpert',    0, pstate_cam%pbl_tpert_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'snowh',    0, pstate%snowhland_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'asdir',    0, pstate%atm_in_asdir_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'asdif',    0, pstate%atm_in_asdif_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'aldir',    0, pstate%atm_in_aldir_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'aldif',    0, pstate%atm_in_aldif_at_pc_surface%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'lwup',     0, pstate%atm_in_lwup_at_pc_surface%f)
        if(trim(shallow_scheme) .eq. 'double_plume')then
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'pcape_previous',     0, pstate_cam%pcape_previous%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'pcape_before_d',     0, pstate_cam%pcape_before_dycore%f)
        end if
        if(trim(macrop_scheme) .eq. 'lin')then
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'wstarPBL', 0, pstate_cam%pbl_wstarPBL%f)
        end if
        call wrap_allocate_data1d(mesh,output1d_temp1)
        call wrap_read_1d_group_rst(mesh%gcomm_read, outdir, fname_restart1, 'cush' ,    0, output1d_temp1%f)
        pstate_cam%cumulus_cush%f(old_time_level,1:ncell) = output1d_temp1%f(1:ncell)
        call wrap_deallocate_data1d(output1d_temp1)
 
!2D
        call wrap_allocate_data2d(mesh,nlev,output2d_temp1)
        call wrap_allocate_data2d(mesh,nlev,output2d_temp15)
        if(trim(macrop_scheme) .eq. 'park')then
            call wrap_allocate_data2d(mesh,nlev,output2d_temp2)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp3)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp4)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp5)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp6)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp7)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp8)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp9)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp10)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp11)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp12)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp13)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp14)
            call wrap_allocate_data2d(mesh,nlev,output2d_temp16)
        end if
        call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'tke' ,    levp, 0, pstate_cam%pbl_tke_at_pc_face_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'kvm' ,    levp, 0, pstate_cam%pbl_kvm_at_pc_face_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'kvh' ,    levp, 0, pstate_cam%pbl_kvh_at_pc_face_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cldo',    lev,  0, output2d_temp1%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cld' ,    lev,  0, output2d_temp15%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'qrs' ,    lev,  0, pstate_cam%sw_qrs_at_pc_full_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'qrl' ,    lev,  0, pstate_cam%lw_qrl_at_pc_full_level%f)
        if(trim(macrop_scheme) .eq. 'park')then
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'qcwat',    lev,  0, output2d_temp2%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'lcwat',    lev,  0, output2d_temp3%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'iccwat',   lev,  0, output2d_temp4%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'nlwat',    lev,  0, output2d_temp5%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'niwat',    lev,  0, output2d_temp6%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'tcwat',    lev,  0, output2d_temp7%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cc_t',     lev,  0, output2d_temp8%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cc_qv',    lev,  0, output2d_temp9%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cc_ql',    lev,  0, output2d_temp10%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cc_qi',    lev,  0, output2d_temp11%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cc_nl',    lev,  0, output2d_temp12%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cc_ni',    lev,  0, output2d_temp13%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'cc_qlst',  lev,  0, output2d_temp14%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'concld',   lev,  0, output2d_temp16%f)
        end if
        if(trim(macrop_scheme) .eq. 'lin')then
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'lengi',    levp, 0, pstate_cam%pbl_lengi_at_pc_face_level%f)
            call wrap_read_2d_group_rst(mesh%gcomm_read, outdir, fname_restart2, 'shi',      levp, 0, pstate_cam%pbl_shi_at_pc_face_level%f)
        end if

        pstate_cam%microp_cldo_at_pc_full_level%f(old_time_level, 1:lev, 1:ncell) = output2d_temp1%f(1:lev,  1:ncell)
        pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level, 1:lev, 1:ncell)  = output2d_temp15%f(1:lev, 1:ncell)
        if(trim(macrop_scheme) .eq. 'park')then
            pstate_cam%macrop_qcwat_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp2%f(1:lev,  1:ncell)
            pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp3%f(1:lev,  1:ncell)
            pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,1:lev,1:ncell)= output2d_temp4%f(1:lev,  1:ncell)
            pstate_cam%macrop_nlwat_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp5%f(1:lev,  1:ncell)
            pstate_cam%macrop_niwat_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp6%f(1:lev,  1:ncell)
            pstate_cam%macrop_tcwat_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp7%f(1:lev,  1:ncell)
            pstate_cam%microp_cc_t_at_pc_full_level%f(old_time_level,1:lev,1:ncell)  = output2d_temp8%f(1:lev,  1:ncell)
            pstate_cam%microp_cc_qv_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp9%f(1:lev,  1:ncell)
            pstate_cam%microp_cc_ql_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp10%f(1:lev, 1:ncell)
            pstate_cam%microp_cc_qi_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp11%f(1:lev, 1:ncell)
            pstate_cam%microp_cc_nl_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp12%f(1:lev, 1:ncell)
            pstate_cam%microp_cc_ni_at_pc_full_level%f(old_time_level,1:lev,1:ncell) = output2d_temp13%f(1:lev, 1:ncell)
            pstate_cam%microp_cc_qlst_at_pc_full_level%f(old_time_level,1:lev,1:ncell)=output2d_temp14%f(1:lev, 1:ncell)
            pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,1:lev,1:ncell)= output2d_temp16%f(1:lev, 1:ncell)
        end if

        call wrap_deallocate_data2d(output2d_temp1)
        call wrap_deallocate_data2d(output2d_temp15)
        if(trim(macrop_scheme) .eq. 'park')then
            call wrap_deallocate_data2d(output2d_temp2)
            call wrap_deallocate_data2d(output2d_temp3)
            call wrap_deallocate_data2d(output2d_temp4)
            call wrap_deallocate_data2d(output2d_temp5)
            call wrap_deallocate_data2d(output2d_temp6)
            call wrap_deallocate_data2d(output2d_temp7)
            call wrap_deallocate_data2d(output2d_temp8)
            call wrap_deallocate_data2d(output2d_temp9)
            call wrap_deallocate_data2d(output2d_temp10)
            call wrap_deallocate_data2d(output2d_temp11)
            call wrap_deallocate_data2d(output2d_temp12)
            call wrap_deallocate_data2d(output2d_temp13)
            call wrap_deallocate_data2d(output2d_temp14)
            call wrap_deallocate_data2d(output2d_temp16)
        end if

! 3D    
        fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".ghg.nc"
        ndim3 = pstate_cam%total_ghg_num
        call wrap_allocate_data3d(mesh,nlev,int(ndim3,i4),output3d_temp1)
        call wrap_read_3d_group_rst(mesh%gcomm_read, outdir, fname_restart3, 'ghg', lev, 0, output3d_temp1%f, ndim3)
        do itracer = 1, ndim3
            pstate_cam%ghg_at_pc_full_level(itracer)%f(1:nlev,1:ncell) = output3d_temp1%f(itracer,1:nlev,1:ncell)
        end do
        call wrap_deallocate_data3d(output3d_temp1)

        fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".aero.nc"
        ndim3 = pstate_cam%total_aerosol_num
        call wrap_allocate_data3d(mesh,nlev,int(ndim3,i4),output3d_temp1)
        call wrap_read_3d_group_rst(mesh%gcomm_read, outdir, fname_restart3, 'aero', lev, 0, output3d_temp1%f, ndim3)
        do itracer = 1, ndim3
            pstate_cam%aerosol_at_pc_full_level(itracer)%f(1:nlev,1:ncell) = output3d_temp1%f(itracer,1:nlev,1:ncell)
        end do
        call wrap_deallocate_data3d(output3d_temp1)

        fname_restart3 = "GRIST.RST.Phy."//trim(grid_info)//".3d."//trim(working_mode)//"."//trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"."//trim(day)//"-"//trim(sec)//".dgnum.nc"
        call rad_cnst_get_info(0, nmodes=nmodes)
        ndim3 = nmodes*2
        call wrap_allocate_data3d(mesh,nlev,int(ndim3,i4),output3d_temp1)
        call wrap_read_3d_group_rst(mesh%gcomm_read, outdir, fname_restart3, 'dgnum', lev, 0, output3d_temp1%f, ndim3)
        do itracer = 1,nmodes 
            pstate_cam%aerosol_dgnum%f(itracer,1:nlev,1:ncell)      = output3d_temp1%f(itracer,1:nlev,1:ncell)
            pstate_cam%aerosol_dgnumwet%f(itracer,1:nlev,1:ncell)   = output3d_temp1%f(nmodes+itracer,1:nlev,1:ncell)
        end do
        call wrap_deallocate_data3d(output3d_temp1)

    end if
    end subroutine gcm_read_restart_file_phy1



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

#endif
  end module grist_gcm_restart_phys1_module
