
!----------------------------------------------------------------------------------
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
!                   2. Bring grist_gcm_tracer_init into this module
!
!    (This is not a final mark)
!    A working_mode contains its own data, computational logic, interface, and 
! experimental protocols. Note that SCM and SWM contain independent drivers, 
! not treated as a working_mode here, because they are not 3D model;
!
!    SWM is used to test the model infrastructure/framework in a 2D mode;
!    SCM is used to test the model physics   in a 1D mode;
!
!   (1) dycore: dry dynamical core [DCMIP2008-2012]
!   (2) tracer: 3D passive tracer transport [DCMIP1]
!   (3) dtp: dycore-tracer-physics coupling [DCMIP2016]
!   (4) amp: inherit dtp but with extension to water surface
!          [real/simple physics aqua-(mountain) planet]; this
!          mode should regress dtp when using the same config;
!
!      These above are prerequiste before a real-world model; when a real-world
!   model is ready, we may construct various model instances (copy this
!   gcm_control_driver,  name a new file, and call it in a branch in grist_atmos),
!   depending on different target applications, provided that these applications 
!   would require substantially different model configurations, for instance:
!
!   (5) atm: Typical weather-climate scale (10-15 km to 100 km) atmosphere model [AMIP]
!   (6) gcrm: GCRM (<=10-15 km) configuration (AMIP, DyaMond?)
!
!     If these applications only differ from, e.g., model physics, or
!  some small changes, this should be done by a nml config while following
!  an existing model instance and working_mode
!
!----------------------------------------------------------------------------------

 module grist_gcm_tracer_control_driver

   use grist_constants,      only: i4, r8, rearth
   use grist_domain_types,   only: block_structure, global_domain

   use grist_nml_module,     only: model_timestep,dycore_timestep,tracer_timestep,&
                                   nlev, ntracer, nmif, nexpdif_level, do_dycore_laplacian_2nd, do_tracer_laplacian_2nd, &
                                   nsteps, working_mode, outdir,&
                                   run_type,h1_restart_freq,       &
                                   h1_history_freq, fname_output,  &
                                   test_dual_cell, use_tr_mbs,  &
                                   write_history, write_restart_h1,&
                                   gcm_testcase, use_phys, physpkg, ptend_f2_on, &
                                   test_real_case, isolate_tracer_test
! mpi
   use grist_mpi
! io
   use grist_gcm_io_h1_module,              only: gcm_output_atm_h1_file
 
! dycore
!   use grist_dycore_time_integration_2d, only: dycore_time_integration_run

   use grist_dycore_diagnose_module_2d,  only: dycore_diagnose_total_energy,&
                                               dycore_diagnose_mass,        &
                                               dycore_diagnose_mass_pt,     &
                                               dycore_diagnose_variables
! static data
!   use grist_datam_static_data_module,   only: grist_static_data_generate_sst
! diffusion
!   use grist_dycore_diffusion_module,    only: grist_dycore_diffusion_run, grist_tracer_diffusion_run
! tracer transport
   use grist_tracer_transport_time_integration_hvsplit, only: tracer_transport_time_integration_run
   use grist_tracer_transport_diagnose_module,          only: tracer_transport_diagnose_global_tracer_mass, &
                                                              tracer_transport_diagnose_global_crnum
   use grist_tracer_transport_vars_module,              only: tracer_transport_evaluate_mif
! dtp
!   use grist_dtp_coupling_module,                       only: grist_dtp_coupling_driver_tend,&
!                                                              grist_dtp_coupling_driver_coup,&
!                                                              grist_d2t_coupling_driver
! gcm-diagnose
   !use grist_gcm_diagnose_h1_module,      only: gcm_h1_accu_physics_variables, &
   !                                             gcm_h1_dump_physics_variables, &
   !                                             gcm_h1_rest_physics_variables
! restart
   use grist_gcm_restart_module,  only: restart_timestep,&
                                        itimestep_pre,&
                                        gcm_read_restart_file,&
                                        gcm_write_restart_file_h1

   implicit none

   private
   public   :: grist_gcm_tracer_init,&
               grist_gcm_tracer_run, &
               grist_gcm_tracer_final
 
! module local
   real(8)              :: time_beg, step_beg
   real(8)              :: time_end, step_end
   real(8)              :: time_elapse, step_elapse
   real(r8)             :: potential_enstropy
   integer(i4)          :: itimestep
   integer(i4)          :: it
   integer(i4)          :: ie
   integer(i4)          :: iv
   integer(i4)          :: istep_beg
   integer(i4)          :: istep_end
   logical              :: input_tend
   integer(i4)          :: tmplevl,ierr
   integer(i4)          :: idstep, dstep_in_mstep, dstep_in_tstep    ! used for dtp
   integer(i4)          :: itstep, tstep_in_mstep    ! used for dtp
   real(r8)             :: sum_time,max_time

   contains

!=========================================================
! Name convention for 'init'
! construct: init data structure for a mode (e.g., dycore vars)
! initial:   init initial condition
! init:      init data used across the target module
! these above are used in a loose way, if some usages
! deviate from the definition above, we should pull it 
! all init now moved here
!=========================================================

   subroutine grist_gcm_tracer_init(local_block)
! dycore
   use grist_dycore_vars_module,         only: grist_dycore_vars_construct
!   use grist_dycore_initial_module,      only: grist_dycore_initial
!   use grist_dycore_time_integration_2d, only: dycore_time_integration_init
! dycore-diffusion
!   use grist_dycore_diffusion_module,    only: grist_diffusion_init
! dycore-hpe
   use grist_hpe_constants,              only: init_hpe_constants
! Tracer Transport
   use grist_tracer_transport_vars_module,              only: tracer_transport_vars_construct
   use grist_tracer_transport_time_integration_hvsplit, only: tracer_transport_time_integration_init
   use grist_tracer_transport_prescribe_module,         only: tracer_transport_prescribe_ambient, &
                                                              tracer_transport_prescribe_initial
   use grist_tracer_transport_ffsl_module,              only: tracer_transport_ffsl_flux_sm10_init
! dtp
!   use grist_dtp_vars_module,            only: grist_dtp_vars_construct
!   use grist_dtp_initial_module,         only: grist_dtp_initial
!   use grist_dtp_interface_module,       only: grist_dtp_interface_init
! datam
!   use grist_datam_static_data_module,   only: grist_static_data_construct
!   use grist_datam_initial_data_module,  only: grist_initial_data_construct
! diagnose
   use grist_gcm_diagnose_h1_module,        only: gcm_h1_diagnose_init
! Physics
!   use grist_physics_data_structure,     only: grist_physics_data_structure_construct
! full physpkg
!   use grist_dtp_interface_full_physpkg1,only: grist_full_physpkg1_init

    type(block_structure) ,target, intent(inout)  :: local_block

!-------------------------------------------------------------------
!                 Component data construction
!-------------------------------------------------------------------

!
! dycore data
!
     call grist_dycore_vars_construct(local_block%full_domain)
!
! tracer data
!
     call tracer_transport_vars_construct(local_block%full_domain)
!
! physics data, only common data that are not related to specific physpkg done here
!
!     call grist_physics_data_structure_construct(local_block%full_domain)
!
! dtp data
!
!     call grist_dtp_vars_construct(local_block%full_domain)
!
! datam static data
!
!     if(test_real_case) call grist_static_data_construct(local_block%full_domain)
!     if(test_real_case) call grist_initial_data_construct(local_block%full_domain)

!-------------------------------------------------------------------
!                          State initial
!-------------------------------------------------------------------

!
! hpe initial
!
     call init_hpe_constants
!
! dycore initial
!
!     select case(trim(working_mode)) 
!     case('dycore')
!        call grist_dycore_initial(local_block%full_domain, gcm_testcase)
!     case('dtp','amp')
!        call grist_dycore_initial(local_block%full_domain, gcm_testcase) ! for reg
!        ! overwrite above if activated
!        call grist_dtp_initial(local_block%full_domain, gcm_testcase)
!     case('tracer')
!
! tracer initial
!
     if(isolate_tracer_test)then
        local_block%full_domain%ne = local_block%full_domain%ne_halo(2)
#ifdef USE_HALO2
        local_block%full_domain%ne = local_block%full_domain%ne_halo(1)
#endif
        call tracer_transport_prescribe_ambient(local_block%full_domain,0._r8,gcm_testcase)
        call tracer_transport_prescribe_initial(local_block%full_domain,gcm_testcase)
        local_block%full_domain%ne = local_block%full_domain%ne_full
     else
        print*,"isolate_tracer_test must be true if running tracer mode"
        call mpi_abort()
     end if
!     case default
!        print*,"GRIST_atmos: you must select a working mode in namelist"
!        stop
!     end select

!-------------------------------------------------------------------
!                     Specific module init
!-------------------------------------------------------------------

!
! dycore-time-integration
!
!     call dycore_time_integration_init(local_block%full_domain)
!
! diffusion
!
!     call grist_diffusion_init(local_block%full_domain)
!
! tracer transport
!
     call tracer_transport_time_integration_init(local_block%full_domain)
!
! dtp coupling
!
!     call grist_dtp_interface_init(local_block%full_domain,nlev,ntracer)
!
! gcm io and diagnose
!
     call gcm_h1_diagnose_init(local_block%full_domain)

#ifdef USE_HALO2
     local_block%full_domain%ne = local_block%full_domain%ne_halo(1)
#else
     local_block%full_domain%ne = local_block%full_domain%ne_halo(2)
#endif
     call tracer_transport_ffsl_flux_sm10_init(local_block%full_domain)
     local_block%full_domain%ne = local_block%full_domain%ne_full

!
! full physics package
!
!     call grist_full_physpkg1_init(local_block, model_timestep)

     local_block%full_domain%nv = local_block%full_domain%nv_compute
     local_block%full_domain%ne = local_block%full_domain%ne_compute
     local_block%full_domain%nt = local_block%full_domain%nt_compute

     return
  end subroutine grist_gcm_tracer_init

!================================================
!          Major stepping of the model
!================================================

  subroutine grist_gcm_tracer_run(mesh)
! io
   type(global_domain), intent(inout) :: mesh

!================================================
!      initial conditions and topography
!================================================

     istep_beg  = 1
     istep_end  = nsteps ! model step
!
! how many dycore and tracer steps withing one model step
! normally tracer_timestep is 3 to 5 times than dycore_timestep
!
     dstep_in_mstep = int(model_timestep/dycore_timestep) ! typically, this should be 3-5
     dstep_in_tstep = int(tracer_timestep/dycore_timestep)
     tstep_in_mstep = int(model_timestep/tracer_timestep) ! typically, this should be 1

     if(mpi_rank().eq.0)then
       print*,"dstep_in_mstep=",dstep_in_mstep
       print*,"dstep_in_tstep=",dstep_in_tstep
       print*,"tstep_in_mstep=",tstep_in_mstep
     end if

     if(trim(run_type).eq.'restart')then ! overwrite data and istep_beg
       call gcm_read_restart_file(mesh)
       istep_beg = restart_timestep
       if(mpi_rank()==0) then
         print*,"=========================================================="
         print*,"    gcm_read_restart_file, restart step =", istep_beg
         print*,"=========================================================="
       endif
     end if

!================================================
!           WRITE INITIAL CONDITIONS
!================================================

     if(istep_beg == 1) then
       call dycore_diagnose_variables(mesh)
       call dycore_diagnose_total_energy(mesh)
       call dycore_diagnose_mass(mesh)
       call dycore_diagnose_mass_pt(mesh)

       call tracer_transport_diagnose_global_tracer_mass(mesh)
       call tracer_transport_diagnose_global_crnum(mesh,tracer_timestep)
     end if
     IF(trim(run_type).eq.'init'.and.write_history) call gcm_output_atm_h1_file(mesh,0)

!================================================
!      major working mode and time marching
! 1) dycore: dry dycore test (need special config to reproduce JAMES2019 results)
! 2) tracer: passive tracer transport test
! 3) dtp: moist model test with/without phys
!================================================

     call cpu_time(time_beg)

!    SELECT CASE(trim(working_mode))
!    CASE('dycore') ! pure dry core test
!      call grist_gcm_control_run_dycore(mesh)
!    CASE('tracer') ! pure tracer test
      call grist_gcm_control_run_tracer(mesh)
!    CASE('dtp')    ! moist model test, inside each component, we only have RK sub-cycle
!      call grist_gcm_control_run_dtp(mesh)
!    CASE('amp')    ! inherit dtp with extensions to amp 
!      call grist_gcm_control_run_amp(mesh)
!    CASE DEFAULT
!         if(mpi_rank() == 0) print*," you must select a working mode"
!         stop
!    END SELECT
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

     return
  end subroutine grist_gcm_tracer_run

!---------------------------------------------------
!                 Tracer Transport 
!---------------------------------------------------

  subroutine grist_gcm_control_run_tracer(mesh)
   type(global_domain), intent(inout) :: mesh
    DO itimestep = istep_beg, istep_end
       if(mpi_rank() == 0) print*,"CURRENT TIMESTEP=",ITIMESTEP,":",istep_end
       call cpu_time(step_beg)
       !integration
       call tracer_transport_time_integration_run(mesh,tracer_timestep,itimestep,gcm_testcase)
       !diagnose
       if(itimestep .gt. itimestep_pre)then
         call tracer_transport_diagnose_global_tracer_mass(mesh)
         call tracer_transport_diagnose_global_crnum(mesh,tracer_timestep)
         call dycore_diagnose_total_energy(mesh)
         call dycore_diagnose_mass(mesh)
         call dycore_diagnose_mass_pt(mesh)
       end if
       !I/O
       if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0).and.write_history)then
          call dycore_diagnose_variables(mesh)
          call gcm_output_atm_h1_file(mesh, itimestep)
       end if

       if(mpi_rank() == 0 .and. trim(run_type).eq.'restart') then
         open(1,file=trim(outdir)//'step_stop.txt',status='unknown')
         write(1,*) itimestep
         close(1)
       endif
       !if((itimestep==nsteps .or. mod(itimestep,restart_freq)==0).and.write_restart )then
       !   call gcm_write_restart_file(mesh,(itimestep+1))
       !   if(mpi_rank()==0) then
       !       print*,"=========================================================="
       !       print*,"    gcm_write_restart_file, write step =", itimestep,restart_freq
       !       print*,"=========================================================="
       !   endif
       !end if
       call cpu_time(step_end)
       step_elapse = step_end-step_beg
       if(mpi_rank() == 0) print*, "elapsed time of this step is ", step_elapse
    END DO

   return
  end subroutine grist_gcm_control_run_tracer

  subroutine grist_gcm_tracer_final(mesh)

!   use grist_datam_static_data_module,     only: grist_static_data_destruct
!   use grist_datam_initial_data_module,    only: grist_initial_data_destruct
!   use grist_dycore_time_integration_2d,   only: dycore_time_integration_final
   use grist_dycore_vars_module,           only: grist_dycore_vars_destruct
!   use grist_dycore_diffusion_module,      only: grist_diffusion_destruct 
!   use grist_dycore_ref_atmos,             only: grist_ref_atmos_final
!   use grist_physics_data_structure,       only: grist_physics_data_structure_destruct
!   use grist_dtp_vars_module,              only: grist_dtp_vars_destruct
!   use grist_dtp_interface_module,         only: grist_dtp_interface_destruct
   use grist_gcm_diagnose_h1_module,       only: gcm_h1_diagnose_final

   use grist_tracer_transport_vars_module, only: tracer_transport_vars_destruct
   use grist_tracer_transport_time_integration_hvsplit, only: tracer_transport_time_integration_final
   use grist_tracer_transport_ffsl_module, only: tracer_transport_ffsl_flux_sm10_final

!   use grist_dtp_interface_full_physpkg1,  only: grist_full_physpkg1_final

   type(global_domain), intent(inout) :: mesh

! data clean and destruct
!     call grist_static_data_destruct
!     call grist_initial_data_destruct
!     call dycore_time_integration_final
     call grist_dycore_vars_destruct
!     call grist_diffusion_destruct
!     call grist_ref_atmos_final
!     call grist_physics_data_structure_destruct
!     call grist_dtp_vars_destruct
!     call grist_dtp_interface_destruct
     call gcm_h1_diagnose_final
! tracer clean
     call tracer_transport_time_integration_final
     call tracer_transport_ffsl_flux_sm10_final
     call tracer_transport_vars_destruct
! full physpkg
!     call grist_full_physpkg1_final
    return
  end subroutine grist_gcm_tracer_final

  end module grist_gcm_tracer_control_driver
