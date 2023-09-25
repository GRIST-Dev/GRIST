
!----------------------------------------------------------------------------
! Created on 2018
! Author: Yi Zhang
! Version 1.0
! Description: Update at 2019 with the working mode
!              Main Driver for the GCM model
!          (1) Set global vars
!          (2) Initial conditions
!          (3) Time_integration
!          (4) Diagnose
!          (5) File I/O
!
! Revision history: 1. 201908, seperate working modes with isolated subroutines
!                   2. scm control is isolated from gcm as a clean add
!                      which requires mesh info and some parameters from
!                      the model part
!----------------------------------------------------------------------------

 module grist_scm_control_driver
! setting
   use grist_constants,      only: i4, r8, rearth, deg2rad
   use grist_domain_types,   only: global_domain
   use grist_nml_module,     only: model_timestep, nsteps, h1_history_freq, run_type, &
                                   write_history, outdir, nlevp, ntracer, doAquaPlanet
! mpi 
   use grist_mpi
! scm
   use grist_scm_comm_module,            only: grist_scm_vars_construct,           &
                                               grist_scm_vars_destruct,            &
                                               scm_lat, scm_lon,                   &
                                               start_ymd, start_tod  
   use grist_scm_nml_module,             only: grist_scm_read_nml
   use grist_scm_io_module,              only: grist_scm_io_output
   use grist_scm_coupling_module,        only: grist_scm_dtp_coupling,             &
                                               grist_scm_ptd_coupling,             &
                                               grist_scm_coupling_init 
   use grist_domain_types,               only: block_structure
#ifdef SCAM
   use grist_scm_dyn_core,               only: grist_scm_dyn_init,                 &
                                               grist_scm_dyn_end,                  &
                                               grist_scm_dyn_pkg
#else
   use grist_scm_dyn_core2,              only: grist_scm_dyn_init2,                &
                                               grist_scm_dyn_end2,                 &
                                               grist_scm_dyn_pkg2
#endif
   use grist_scm_surface_module,         only: grist_scm_surface 
#ifdef AMIPC_PHYSICS
! PhysC-physics interface
   use grist_physpkg_PhysC,          only: grist_physpkg_PhysC_init,       &
                                               grist_physpkgbc_PhysC_run,      &
                                               grist_physpkgac_PhysC_run,      &
                                               grist_physpkg_PhysC_final
#endif

#ifdef AMIPW_PHYSICS
   use grist_PhysW_nml_module,         only: set_PhysW_nml
   use grist_PhysW_pkg,                only: grist_PhysW_pkg_init,                &
                                               grist_PhysW_pkg_run_ac,              &
                                               grist_PhysW_pkg_run_bc,              &
                                               grist_PhysW_pkg_final
   use grist_zenith,                     only: zenith,orb_params
   use grist_PhysW_data_structure,       only: pstate_wrf
   use grist_time_manager,               only: get_curr_calday 
#endif

! scm-diagnose
   use grist_scm_diagnose_module,        only: scm_accumulate_physics_variables,  &
                                               scm_dump_physics_variables,        &
                                               scm_reset_physics_variables,       &
                                               scm_diagnose_init

   implicit none

   private
   public   :: grist_scm_init,&
               grist_scm_run, &
               grist_scm_final
 
! module local
   integer(i4)          :: itimestep
   integer(i4)          :: istep_beg
   integer(i4)          :: istep_end

   contains

  subroutine grist_scm_init(local_block)
    use grist_physics_data_structure,     only: grist_physics_data_structure_construct
#ifdef AMIPC_PHYSICS
    use grist_PhysC_data_structure,        only: grist_PhysC_data_structure_construct
#endif
    use grist_hpe_constants,              only: init_hpe_constants

     type(block_structure) ,target, intent(inout)  :: local_block

     call grist_physics_data_structure_construct(local_block%full_domain)
#ifdef AMIPC_PHYSICS
     call grist_PhysC_data_structure_construct(local_block%full_domain)
#endif
! read static file, read checked ok, currently not used for dev
! now has been combined to static data construct
!     call gcm_read_file(local_block%full_domain,comm)
! hpe initial
     call init_hpe_constants

     call scm_diagnose_init(local_block%full_domain)

     local_block%full_domain%nv = local_block%full_domain%nv_compute
     local_block%full_domain%ne = local_block%full_domain%ne_compute
     local_block%full_domain%nt = local_block%full_domain%nt_compute

     return
  end subroutine grist_scm_init

!================================================
!          Major stepping of the model
!================================================

  subroutine grist_scm_run(mesh)
! io
   type(global_domain), intent(inout) :: mesh
   real(r8)  :: dtime
   real(r8)  :: clat(1), clon(1)
   real(r8)  :: time_beg, step_beg
   real(r8)  :: time_end, step_end
   real(r8)  :: time_elapse, step_elapse
   real(r8)  :: sum_time,max_time 
   real(r8)  :: julday
   integer(i4) :: year 
   integer(i4) :: ierr

   real(r8), allocatable    :: dxmean(:), coszrs(:)


      if(mpi_rank() == 0) print*, ""
      if(mpi_rank() == 0) print*, "Start to run GRIST Single Column Model"
#ifdef SCAM
      istep_beg = 0
      istep_end  = nsteps ! model step
#else
      istep_beg = 1
      istep_end  = nsteps ! model step
#endif
      !reading nml and initial data.
      call grist_scm_read_nml
      clat(1) = scm_lat*deg2rad
      clon(1) = scm_lon*deg2rad

      !use 120km resolution in SCM, LiXH
      allocate(dxmean(1));dxmean=120._r8*1000.

      call grist_scm_vars_construct
      call grist_scm_coupling_init
#ifdef AMIPC_PHYSICS
      call grist_physpkg_PhysC_init(1, model_timestep, istep_beg, dxmean, clat, clon)
      if(mpi_rank()==0) print*,"Sucessfully initiate CAM5 model physics"
#endif
#ifdef AMIPW_PHYSICS
      allocate(coszrs(1));coszrs=0.
      call set_PhysW_nml
      call grist_PhysW_pkg_init(1, 1, nlevp, ntracer, model_timestep, clat, clon)
      if(mpi_rank()==0) print*,"Sucessfully initiate WRF2 model physics"
      pstate_wrf%dxmean(1) = dxmean(1)
      
      year  = start_ymd/10000
      call orb_params( year, .true.)  
#endif

#ifdef SCAM
      call grist_scm_dyn_init
#else
      call grist_scm_dyn_init2
#endif
      if(trim(run_type) .eq. 'init') call grist_scm_io_output(mesh,0)
!
!  begin
!
      call cpu_time(time_beg)

      DO itimestep = istep_beg, istep_end
          if(mpi_rank() == 0) print*, "itimestep=", itimestep
#ifdef SCAM
          if(itimestep .ge. 1)then
              dtime = 2.*model_timestep
          else
              dtime = model_timestep
          end if
#else
          dtime = model_timestep
#endif
          call grist_scm_dtp_coupling

#ifdef AMIPC_PHYSICS
          call grist_physpkgbc_PhysC_run(1, clat, clon, itimestep, dtime, model_timestep)

          call scm_accumulate_physics_variables(1)
          if((itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0).and.write_history)then
              call scm_dump_physics_variables
              call grist_scm_io_output(mesh, itimestep)
              call scm_reset_physics_variables
          end if

          call grist_scm_surface(itimestep, dtime, clat, clon)
          call grist_physpkgac_PhysC_run(1, clat, clon, itimestep, dtime, model_timestep)
#endif

#ifdef AMIPW_PHYSICS
          call grist_PhysW_pkg_run_bc(1, nlevp, ntracer, itimestep, dtime,rearth*mesh%mean_edt_dist)
          call get_curr_calday(start_ymd, start_tod, itimestep, dtime, julday)
          call zenith(julday, coszrs, 1, clon, clat,doAquaPlanet)
          call grist_scm_surface(itimestep, dtime, clat, clon, coszrs)
          call grist_PhysW_pkg_run_ac(1, nlevp, ntracer, itimestep, dtime,rearth*mesh%mean_edt_dist, coszrs)

          call scm_accumulate_physics_variables(1)
          if((itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0).and.write_history)then
              call scm_dump_physics_variables
              call grist_scm_io_output(mesh, itimestep)
              call scm_reset_physics_variables
          end if
#endif

          call grist_scm_ptd_coupling(dtime)
#ifdef SCAM
          call grist_scm_dyn_pkg(itimestep, dtime)
#else
          call grist_scm_dyn_pkg2(itimestep, dtime)
#endif

      END DO
#ifdef SCAM
      call grist_scm_dyn_end
#else
      call grist_scm_dyn_end2
#endif

#ifdef AMIPC_PHYSICS
      call grist_physpkg_PhysC_final
#endif

#ifdef AMIPW_PHYSICS
      call grist_PhysW_pkg_final  
#endif

      call grist_scm_vars_destruct

      if(mpi_rank() == 0) print*, ""
      if(mpi_rank() == 0) print*, "End of GRIST Single Column Model"

!
! Final elapsed time statistics
! 
     call cpu_time(time_end)
     time_elapse = time_end-time_beg
     sum_time = 0._r8
     call reduce(time_elapse, sum_time, 'sum')

     sum_time = sum_time/mpi_size()
     max_time = 0._r8
     call reduce(time_elapse, max_time, 'max')

     if(mpi_rank() == 0)then
        open(1,file=trim(outdir)//"/elapsed_time.txt",status='unknown')
        write(1,*) "average elapsed time:", sum_time,"  secs"
        write(1,*) "max elapsed time:"    , max_time,"  secs"
        write(1,*) "number of processors:", mpi_size()
        write(1,*) "max time_elapsed="    , max_time,"  secs"
        write(1,*) "average time_elapsed=", sum_time,"  secs"
        write(1,*) "number of processors=", mpi_size()
        close(1)
     end if

     deallocate(dxmean)

  end subroutine grist_scm_run

  subroutine grist_scm_final(local_block)
    ! physics data structure
  use grist_physics_data_structure,     only: grist_physics_data_structure_destruct
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,        only: grist_PhysC_data_structure_destruct
#endif
  use grist_scm_io_module,              only: scm_output_final

     type(block_structure) ,target, intent(inout)  :: local_block

    call grist_physics_data_structure_destruct
#ifdef AMIPC_PHYSICS
    call grist_PhysC_data_structure_destruct
#endif

    call scm_output_final(local_block%full_domain)
  end subroutine grist_scm_final

  end module grist_scm_control_driver
