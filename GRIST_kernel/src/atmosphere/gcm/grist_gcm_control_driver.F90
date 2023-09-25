
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
!                   2. Bring grist_gcm_init into this module
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

 module grist_gcm_control_driver

   use grist_constants,      only: i4, r8, rearth
   use grist_domain_types,   only: block_structure, global_domain

   use grist_nml_module,     only: model_timestep,dycore_timestep,tracer_timestep,&
                                   nlev, ntracer, nspecies, nmif, nexpdif_level, use_expdif_tstep, &
                                   do_dycore_laplacian_2nd, do_tracer_laplacian_2nd, do_dycore_laplacian_4th, &
                                   nsteps, working_mode, outdir,&
                                   run_type, h1_restart_freq, nlevp,&
                                   h1_history_freq, fname_output,  &
                                   write_history, write_restart_h1, write_restart_h0, &
                                   gcm_testcase, use_phys, physpkg, ptend_f2_on, ptend_f2drbl_on, &
                                   test_real_case, isolate_tracer_test, write_history_h1, write_history_h0, &
                                   write_history_h1_separate, h1_1d_history_freq, h1_2d_history_freq, h1_3d_history_freq
 
! dycore
   use grist_dycore_time_integration_2d, only: dycore_time_integration_run

   use grist_dycore_diagnose_module_2d,  only: dycore_diagnose_total_energy,&
                                               dycore_diagnose_total_energy_moist_e2, &
                                               dycore_diagnose_total_energy_moist_e3, &
                                               dycore_diagnose_mass,        &
                                               dycore_diagnose_mass_pt,     &
                                               dycore_diagnose_variables
! static data
   use grist_datam_static_data_module,   only: grist_static_data_generate_analytic_sst, &
                                               grist_static_fill_current_sst
   use grist_prescribe_seaice_module,    only: grist_seaice_init
! io
   use grist_gcm_io_h1_module,           only: gcm_output_atm_h1_file, gcm_output_history_h1, &
                                               gcm_output_history_h1_1d, &
                                               gcm_output_history_h1_2d, &
                                               gcm_output_history_h1_3d
                                               
   use grist_gcm_io_h0_module,           only: gcm_output_history_h0
! diffusion
   use grist_dycore_diffusion_module,    only: grist_dycore_diffusion_run, grist_tracer_diffusion_run
! mpi
   use grist_mpi
! tracer transport
   use grist_tracer_transport_time_integration_hvsplit, only: tracer_transport_time_integration_run
   use grist_tracer_transport_diagnose_module,          only: tracer_transport_diagnose_global_tracer_mass, &
                                                              tracer_transport_diagnose_global_crnum
   use grist_tracer_transport_vars_module,              only: tracer_transport_evaluate_mif
! dtp
   use grist_dtp_coupling_module,                       only: grist_dtp_coupling_driver_tend,&
                                                              grist_dtp_coupling_driver_coup,&
                                                              grist_d2t_coupling_driver     ,&
                                                              grist_dtp_time_average
! restart
   use grist_gcm_restart_module,  only: restart_timestep,&
                                        itimestep_pre,&
                                        gcm_read_restart_file,&
                                        gcm_write_restart_file_h1, &
                                        gcm_write_restart_file_h0
   use grist_clocks,             only: clock_id, clock_begin, clock_end

#ifdef LAM_DOMAIN
   use grist_nml_module,     only: glamDataInterval, glam_use_relaxation
   use grist_constants,      only: zero
   use grist_datam_glam_data_module, only: grist_glam_add_dynamics_bdy_forcing_term, &
                                           grist_glam_data_read, &
                                           grist_glam_data_time_interp, &
                                           grist_glam_data_idealized_jwss, &
                                           grist_glam_data_idealized_dcmip2016_tc
#endif

   implicit none

   private
   public   :: grist_gcm_init,&
               grist_gcm_run, &
               grist_gcm_final
 
! module local
   real(8)              :: time_beg, step_beg
   real(8)              :: time_end, step_end
   real(8)              :: time_elapse, step_elapse
   integer(i4)          :: itimestep
   !integer(i4)          :: it
   !integer(i4)          :: ie
   !integer(i4)          :: iv
   integer(i4)          :: istep_beg
   integer(i4)          :: istep_end
   !logical              :: input_tend
   !integer(i4)          :: tmplevl,ierr
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

   subroutine grist_gcm_init(local_block)

! gcm diagnose
   use grist_gcm_diagnose_h0_module,     only: gcm_h0_diagnose_init
   use grist_gcm_diagnose_h1_module,     only: gcm_h1_diagnose_init
! datam
   use grist_datam_static_data_module,   only: grist_static_data_construct
   use grist_datam_initial_data_module,  only: grist_initial_data_construct
! dycore
   use grist_dycore_vars_module,         only: grist_dycore_vars_construct
   use grist_dycore_initial_module,      only: grist_dycore_initial
   use grist_dycore_time_integration_2d, only: dycore_time_integration_init
! dycore-diffusion
   use grist_dycore_diffusion_module,    only: grist_diffusion_init
! dycore-hpe
   use grist_hpe_constants,              only: init_hpe_constants
! Tracer Transport
   use grist_tracer_transport_vars_module,              only: tracer_transport_vars_construct
   use grist_tracer_transport_time_integration_hvsplit, only: tracer_transport_time_integration_init
   use grist_tracer_transport_prescribe_module,         only: tracer_transport_prescribe_ambient, tracer_transport_prescribe_initial
! dtp
   use grist_dtp_vars_module,            only: grist_dtp_vars_construct
   use grist_dtp_initial_module,         only: grist_dtp_initial
   use grist_dtp_interface_module,       only: grist_dtp_interface_init
   use grist_dtp_interface_full_physpkg1,only: grist_full_physpkg1_init
   use grist_dtp_interface_full_physpkg2,only: grist_full_physpkg2_init
! Physics
   use grist_physics_data_structure,     only: grist_physics_data_structure_construct

#ifdef LAM_DOMAIN
   use grist_datam_glam_data_module,     only: grist_glam_data_construct
#endif

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
     call grist_physics_data_structure_construct(local_block%full_domain)
!
! dtp data
!
     call grist_dtp_vars_construct(local_block%full_domain)
!
! datam static data
!
     if(test_real_case)then
        call grist_initial_data_construct(local_block%full_domain)! initial data
        call grist_static_data_construct(local_block%full_domain) ! sst from here
        call grist_seaice_init(local_block%full_domain%nv_full)   ! seaice temp init
     else
        call grist_static_data_generate_analytic_sst(local_block%full_domain) ! sst from here, for ap-planet
     end if

#ifdef LAM_DOMAIN
     call grist_glam_data_construct(local_block%full_domain,nlev,nlevp)
#endif
!-------------------------------------------------------------------
!                          State initial
!-------------------------------------------------------------------

!
! hpe initial, necessary for all 3d models
!
     call init_hpe_constants(nlev,nlevp)
!
! dycore initial
!
     select case(trim(working_mode)) 
     case('dycore')
        call grist_dycore_initial(local_block%full_domain, gcm_testcase, nlev, nlevp)
     case('tracer')
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

     case('dtp','amipc','amipw')
        call grist_dycore_initial(local_block%full_domain, gcm_testcase, nlev, nlevp) ! for reg
        ! overwrite above if activated
        call grist_dtp_initial(local_block%full_domain, gcm_testcase)
     case default
        print*,"GRIST_atmos: you must select a working mode in namelist"
        stop
     end select

!-------------------------------------------------------------------
!                     Specific module init
!-------------------------------------------------------------------

!
! dycore-time-integration
!
     call dycore_time_integration_init(local_block%full_domain,nlev)
!
! diffusion
!
     call grist_diffusion_init(local_block%full_domain)
!
! tracer transport
!
     call tracer_transport_time_integration_init(local_block%full_domain)
     local_block%full_domain%ne = local_block%full_domain%ne_full
!
! dtp coupling
!
     call grist_dtp_interface_init(local_block%full_domain,nlev,ntracer)
!
! full physics package
!
     call grist_full_physpkg1_init(local_block, model_timestep)
     call grist_full_physpkg2_init(local_block, model_timestep, local_block%full_domain%nv_halo(1), nlevp, nspecies)
     if(mpi_rank().eq.0) print*,"physpkg1&2 init done"

!
! gcm io and diagnose
!
     local_block%full_domain%nv = local_block%full_domain%nv_full   
     if(write_history_h0) call gcm_h0_diagnose_init(local_block%full_domain)
     if(write_history_h1) call gcm_h1_diagnose_init(local_block%full_domain)
     local_block%full_domain%nv = local_block%full_domain%nv_compute
     local_block%full_domain%ne = local_block%full_domain%ne_compute
     local_block%full_domain%nt = local_block%full_domain%nt_compute

     return
  end subroutine grist_gcm_init

!================================================
!          Major stepping of the model
!================================================

  subroutine grist_gcm_run(mesh)
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
       call dycore_diagnose_total_energy(mesh) ! raw
       if(trim(working_mode).ne.'dycore' .and.trim(working_mode).ne.'tracer')then
          call dycore_diagnose_total_energy_moist_e2(mesh)
          call dycore_diagnose_total_energy_moist_e3(mesh)
       end if
       call dycore_diagnose_mass(mesh)
       call dycore_diagnose_mass_pt(mesh)

       call tracer_transport_diagnose_global_tracer_mass(mesh)
       call tracer_transport_diagnose_global_crnum(mesh,tracer_timestep)
     end if

     IF(trim(run_type).eq.'init'.and.write_history_h1) call gcm_output_history_h1(mesh,0,nsteps,h1_history_freq)

!================================================
!      major working mode and time marching
! 1) dycore: dry dycore test (need special config to reproduce JAMES2019 results)
! 2) tracer: passive tracer transport test
! 3) dtp: moist model test with/without phys
!================================================

    call cpu_time(time_beg)

    SELECT CASE(trim(working_mode))
    CASE('dycore') ! pure dry core test
      call grist_gcm_control_run_dycore(mesh)
    CASE('tracer') ! pure tracer test
      call grist_gcm_control_run_tracer(mesh)
    CASE('dtp')    ! moist model test, inside each component, we only have RK sub-cycle
      call grist_gcm_control_run_dtp(mesh)
    CASE('amipc') ! inherit dtp with extensions to amp
      call grist_gcm_control_run_amipc(mesh)
    CASE('amipw') ! inherit dtp with extensions to amp
      call grist_gcm_control_run_amipw(mesh)
    CASE DEFAULT
         if(mpi_rank() == 0) print*," you must select a working mode"
         stop
    END SELECT
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
  end subroutine grist_gcm_run

!---------------------------------------------------
!             Dry Dynamical Core
!---------------------------------------------------

  subroutine grist_gcm_control_run_dycore(mesh)
   type(global_domain), intent(inout) :: mesh

#ifdef LAM_DOMAIN
    call grist_glam_data_read(mesh,  model_timestep, 0, nlev, nlevp)
    !call grist_glam_data_idealized_jwss(mesh,nlev,nlevp,1)
#endif

    DO itimestep = istep_beg, istep_end

#ifdef LAM_DOMAIN
       if(mod(model_timestep*itimestep,glamDataInterval).eq.zero)then
          call grist_glam_data_read(mesh,  model_timestep, itimestep, nlev, nlevp)
       end if
       call grist_glam_data_time_interp(model_timestep, itimestep)
       call grist_glam_add_dynamics_bdy_forcing_term(mesh,nlev)
#endif

       if(mpi_rank() == 0) print*,"CURRENT TIMESTEP=",ITIMESTEP,":",istep_end
       call cpu_time(step_beg)
       ! call physics forcing here
       if(use_phys.and.trim(physpkg).eq.'HSDRY_EXP') call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep, itimestep, physpkg)
       ! call explicit diffusion tendencies
       if(do_dycore_laplacian_2nd.or.do_dycore_laplacian_4th) call grist_dycore_diffusion_run(mesh, nexpdif_level,ntracer)

       call tracer_transport_evaluate_mif(mesh,nlev) ! actually works nothing but make code work
       !integration
       call dycore_time_integration_run(mesh,dycore_timestep,itimestep)
       !diagnose
       if(itimestep .gt. itimestep_pre)then
         call dycore_diagnose_total_energy(mesh)
         call dycore_diagnose_mass(mesh)
         call dycore_diagnose_mass_pt(mesh)
       end if

       if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0).and.write_history_h1)then
          call dycore_diagnose_variables(mesh)
          call gcm_output_atm_h1_file(mesh, itimestep)
       end if
       !if(write_history_h1) call gcm_output_history_h1(mesh, itimestep,nsteps,h1_history_freq) ! most vars are instant state; default

       if(mpi_rank() == 0 .and. trim(run_type).eq.'restart') then
         open(1,file=trim(outdir)//'step_stop.txt',status='unknown')
         write(1,*) itimestep
         close(1)
       endif
       if(write_restart_h1) call gcm_write_restart_file_h1(mesh,itimestep,nsteps,h1_restart_freq)
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
  end subroutine grist_gcm_control_run_dycore

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
       !if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0).and.write_history)then
       !   call dycore_diagnose_variables(mesh)
       !   call gcm_output_atm_h1_file(mesh, itimestep)
       !end if
       if(write_history_h1) call gcm_output_history_h1(mesh, itimestep,nsteps,h1_history_freq) ! most vars are instant state; default

       if(mpi_rank() == 0 .and. trim(run_type).eq.'restart') then
         open(1,file=trim(outdir)//'step_stop.txt',status='unknown')
         write(1,*) itimestep
         close(1)
       endif
       if(write_restart_h1) call gcm_write_restart_file_h1(mesh,itimestep,nsteps,h1_restart_freq)
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

!---------------------------------------------------
!               Dycore-Tracer-Physics
!---------------------------------------------------

  subroutine grist_gcm_control_run_dtp(mesh)
   type(global_domain), intent(inout) :: mesh
! local
   integer(i4)  :: clock_dtploop

    clock_dtploop = clock_id('DtpMainLoop')

     ! this should recover the previous hsdry_exp case
     !if(use_phys.and.ptend_rk_on) call grist_dtp_coupling_driver(mesh, nlev, ntracer, model_timestep,physpkg)
     ! call diffusion operator here
     !if(do_laplacian_2nd) call grist_diffusion_run(mesh, nexpdif_level,ntracer)
#ifdef LAM_DOMAIN
    call grist_glam_data_read(mesh,  model_timestep, 0, nlev, nlevp)
#endif

    if(use_phys)then
       call tracer_transport_evaluate_mif(mesh,nlev) ! newly added for real-physics, does not affect regression old results,
       call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep,max(1,istep_beg),physpkg)
    end if
    call clock_begin(clock_dtploop)
    DO itimestep = istep_beg, istep_end

#ifdef LAM_DOMAIN
       if(mod(model_timestep*itimestep,glamDataInterval).eq.zero)then
          call grist_glam_data_read(mesh,  model_timestep, itimestep, nlev, nlevp)
       end if
       call grist_glam_data_time_interp(model_timestep, itimestep)
       if(glam_use_relaxation) call grist_glam_add_dynamics_bdy_forcing_term(mesh,nlev)
#endif

       if(mpi_rank() == 0) print*,"CURRENT TIMESTEP=",ITIMESTEP,":",istep_end
       call cpu_time(step_beg)
!
! test calling phys before of dyn and do update
!
       if(do_dycore_laplacian_2nd.or.do_dycore_laplacian_4th) call grist_dycore_diffusion_run(mesh, nexpdif_level,ntracer)
       if(do_tracer_laplacian_2nd) call grist_tracer_diffusion_run(mesh, nexpdif_level,ntracer)
!
! dynamics integration
!
       do itstep = 1, tstep_in_mstep
!
! ptend_f2drbl is similar to f2, but when dynamics-physics are split,
! it is like dribbling as in CAM, should not be used with ptend_f2_on
!
          if(use_phys.and.ptend_f2drbl_on) call grist_dtp_coupling_driver_coup(mesh, nlev, ntracer, tracer_timestep, physpkg)

          call tracer_transport_evaluate_mif(mesh,nlev)
          do idstep = 1, dstep_in_tstep
             call dycore_time_integration_run(mesh,dycore_timestep,itimestep)
             call grist_d2t_coupling_driver(mesh, nlev, idstep, dstep_in_tstep)
             call grist_dtp_time_average(mesh, nlev, nlevp, idstep, dstep_in_tstep,itstep, tstep_in_mstep) ! time-average www/omega
          end do
          call tracer_transport_time_integration_run(mesh,tracer_timestep,itimestep,gcm_testcase)
       end do
!
! diagnostics at model step
!
       if(itimestep .gt. itimestep_pre)then
         call dycore_diagnose_total_energy(mesh)
         call dycore_diagnose_total_energy_moist_e2(mesh)
         call dycore_diagnose_total_energy_moist_e3(mesh)
         call dycore_diagnose_mass(mesh)
         call dycore_diagnose_mass_pt(mesh)
         call tracer_transport_diagnose_global_tracer_mass(mesh)
         call tracer_transport_diagnose_global_crnum(mesh,tracer_timestep)
       end if
!
! call physics forcing here
!
       if(use_phys)then
          call tracer_transport_evaluate_mif(mesh,nlev) ! newly added for real-physics, does not affect regression old results,
                                                        ! which do not depend on this
          call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep,itimestep, physpkg)
          ! update f2 physics
          if(ptend_f2_on) call grist_dtp_coupling_driver_coup(mesh, nlev, ntracer, model_timestep,physpkg)
       end if
!       
! diag accumulated vars each model step; call gcm_accu_physics_variables; I/O
!
       !if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0).and.write_history)then
       !   call dycore_diagnose_variables(mesh)
       !  ! call gcm_dump_physics_variables
       !   call gcm_output_atm_h1_file(mesh, itimestep)
       !  ! call gcm_rest_physics_variables
       !end if
       if(write_history_h1) call gcm_output_history_h1(mesh, itimestep,nsteps,h1_history_freq) ! most vars are instant state; default

       if(mpi_rank() == 0 .and. trim(run_type).eq.'restart') then
         open(1,file=trim(outdir)//'step_stop.txt',status='unknown')
         write(1,*) itimestep
         close(1)
       endif

       if(write_restart_h1) call gcm_write_restart_file_h1(mesh,itimestep,nsteps,h1_restart_freq)

       call cpu_time(step_end)
       step_elapse = step_end-step_beg
       if(mpi_rank() == 0) print*, "elapsed time of this step is ", step_elapse
    END DO
    call clock_end(clock_dtploop)
    return
  end subroutine grist_gcm_control_run_dtp

!---------------------------------------------------
!      Aqua-Mountain-Planet (CAM5-based physpkg)
!---------------------------------------------------

  subroutine grist_gcm_control_run_amipc(mesh)
   type(global_domain), intent(inout) :: mesh
   integer(i4)  :: clock_dynamics, clock_camphys, clock_diagio, clock_diagte, clock_tracer, clock_d2t
   logical      :: is_first_restart

    clock_dynamics= clock_id('dynamics')
    clock_camphys = clock_id('camPhys')
    clock_diagio  = clock_id('diagio')
    clock_diagte  = clock_id('diagte')
    clock_tracer  = clock_id('tracer')
    clock_d2t     = clock_id('d2t')
! local

    if(use_phys.and.istep_beg.eq.1)then
       call tracer_transport_evaluate_mif(mesh,nlev)
       call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep,max(1,istep_beg),physpkg)
    end if
    if(test_real_case)then
        if(istep_beg-1.eq.0)then
           is_first_restart = .false. ! whether this first step is a restart initial step
        else
           is_first_restart = .true.
        end if
        call grist_static_fill_current_sst(mesh,istep_beg-1,mesh%nv_full,model_timestep,is_first_restart) ! only update SST when day is updated
    end if

    DO itimestep = istep_beg, istep_end

       if(mpi_rank() == 0) print*,"CURRENT TIMESTEP=",ITIMESTEP,":",istep_end
       call cpu_time(step_beg)
       if(test_real_case.and.mod(model_timestep*itimestep,86400.).eq.0)then
          if(mpi_rank().eq.0) print*, "update SST at this step"
          call grist_static_fill_current_sst(mesh,itimestep,mesh%nv_full,model_timestep,.false.) ! only update SST when day is updated
       end if
!
! test calling phys before of dyn and do update
!
       if(do_dycore_laplacian_2nd.or.do_dycore_laplacian_4th) call grist_dycore_diffusion_run(mesh, nexpdif_level,ntracer)
       if(do_tracer_laplacian_2nd) call grist_tracer_diffusion_run(mesh, nexpdif_level,ntracer)
!
! integration
!
       call clock_begin(clock_dynamics)
       do itstep = 1, tstep_in_mstep
          call tracer_transport_evaluate_mif(mesh,nlev)
          do idstep = 1, dstep_in_tstep
             call dycore_time_integration_run(mesh,dycore_timestep,itimestep)
             call clock_begin(clock_d2t)
             call grist_d2t_coupling_driver(mesh, nlev, idstep, dstep_in_tstep)
             call grist_dtp_time_average(mesh, nlev, nlevp, idstep, dstep_in_tstep,itstep, tstep_in_mstep) ! time-average www/omega
             call clock_end(clock_d2t)
          end do
          call clock_begin(clock_tracer)
          call tracer_transport_time_integration_run(mesh,tracer_timestep,itimestep,gcm_testcase)
          call clock_end(clock_tracer)
       end do
       call clock_end(clock_dynamics)
!
! diagnostics at model step
!
       call clock_begin(clock_diagte)
       if(itimestep .gt. itimestep_pre)then
         call dycore_diagnose_total_energy(mesh)
         call dycore_diagnose_total_energy_moist_e2(mesh)
         call dycore_diagnose_total_energy_moist_e3(mesh)
         call dycore_diagnose_mass(mesh)
         call dycore_diagnose_mass_pt(mesh)
         call tracer_transport_diagnose_global_tracer_mass(mesh)
         call tracer_transport_diagnose_global_crnum(mesh,tracer_timestep)
       end if
       call clock_end(clock_diagte)
!       
! call physics forcing here
!
       call clock_begin(clock_camphys)
       if(use_phys)then
          call tracer_transport_evaluate_mif(mesh,nlev) ! newly added for real-physics, does not affect regression old results,
                                                        ! which do not depend on this
          call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep,itimestep, physpkg)
          ! update f2 physics
          if(ptend_f2_on) call grist_dtp_coupling_driver_coup(mesh, nlev, ntracer, model_timestep,physpkg)
       end if
       call clock_end(clock_camphys)
!
! I/O: h0: monthly mean, h1: export according to h1_history_freq (default)
!
       ! currently h0 and h1 are mutually exclusive
       call clock_begin(clock_diagio)
       if(write_history_h0) call gcm_output_history_h0(mesh, itimestep) ! monthly-averaged
       if(write_history_h1) call gcm_output_history_h1(mesh, itimestep,nsteps,h1_history_freq) ! most vars are instant state; default
       call clock_end(clock_diagio)

       if(mpi_rank() == 0) then
         open(1,file=trim(outdir)//'step_stop.txt',status='unknown')
         write(1,*) itimestep
         close(1)
       endif

       if(write_restart_h1) call gcm_write_restart_file_h1(mesh,itimestep,nsteps,h1_restart_freq) ! write n-step
       if(write_restart_h0) call gcm_write_restart_file_h0(mesh,itimestep,nsteps)                 ! write each month

       call cpu_time(step_end)
       step_elapse = step_end-step_beg
       if(mpi_rank() == 0) print*, "elapsed time of this step is ", step_elapse
    END DO

    return
  end subroutine grist_gcm_control_run_amipc

!---------------------------------------------------
!      Aqua-Mountain-Planet (WRF-based physpkg)
!---------------------------------------------------
  subroutine grist_gcm_control_run_amipw(mesh)
   type(global_domain), intent(inout) :: mesh
! local
   integer(i4)  :: clock_dynamics, clock_physw, clock_diagio, clock_diagte, clock_tracer, clock_d2t
   integer(i4)  :: clock_outloop, clock_fillsst, clock_write_rst
   logical      :: is_first_restart

    clock_outloop   = clock_id('Outloop')
    clock_fillsst   = clock_id('Fillsst_sic')
    clock_dynamics  = clock_id('Dynamics')
    clock_physw     = clock_id('PhysW')
    clock_diagio    = clock_id('write_history')
    clock_diagte    = clock_id('diagte')
    clock_tracer    = clock_id('tracer_transport')
    clock_d2t       = clock_id('d2t')
    clock_write_rst = clock_id('write_rst')
!
! SST
!
!    call grist_static_data_generate_sst(mesh)

#ifdef LAM_DOMAIN
    if(istep_beg-1.ne.0)then ! restart step
       call grist_glam_data_read(mesh,  model_timestep, istep_beg-1, nlev, nlevp)
    else
       call grist_glam_data_read(mesh,  model_timestep, 0, nlev, nlevp)
    end if

#endif

    if(use_phys.and.istep_beg.eq.1)then
       call tracer_transport_evaluate_mif(mesh,nlev)
       ! this affect amipw regression due to WRF radiation
       call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep,max(1,istep_beg),physpkg) ! old
       !call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep,0,physpkg)  ! new
    end if

    if(test_real_case)then
       if(istep_beg-1.eq.0)then
          is_first_restart = .false. ! whether this first step is a restart initial step
       else
          is_first_restart = .true.
       end if
       call grist_static_fill_current_sst(mesh,istep_beg-1,mesh%nv_full,model_timestep,is_first_restart) ! only update SST when day is updated
    end if

    call clock_begin(clock_outloop)

    DO itimestep = istep_beg, istep_end

#ifdef LAM_DOMAIN
       if(mod(model_timestep*itimestep,glamDataInterval).eq.zero)then
          call grist_glam_data_read(mesh,  model_timestep, itimestep, nlev, nlevp)
       end if
       call grist_glam_data_time_interp(model_timestep, itimestep)
       call grist_glam_add_dynamics_bdy_forcing_term(mesh,nlev)
#endif

       if(mpi_rank() == 0) print*,"CURRENT TIMESTEP=",ITIMESTEP,":",istep_end
       call cpu_time(step_beg)

       call clock_begin(clock_fillsst)
       if(test_real_case.and.mod(model_timestep*itimestep,86400.).eq.0)then
          if(mpi_rank().eq.0) print*, "update SST at this step"
          call grist_static_fill_current_sst(mesh,itimestep,mesh%nv_full,model_timestep,.false.) ! only update SST when day is updated
       end if
       call clock_end(clock_fillsst)
!
! test calling phys before of dyn and do update
! this is the default seting for all previous tests and papers till 2021
!

       if(.not.use_expdif_tstep)then
          if(do_dycore_laplacian_2nd.or.do_dycore_laplacian_4th) call grist_dycore_diffusion_run(mesh, nexpdif_level,ntracer)
          if(do_tracer_laplacian_2nd) call grist_tracer_diffusion_run(mesh, nexpdif_level,ntracer)
       end if
!
! integration
!

       call clock_begin(clock_dynamics)
       do itstep = 1, tstep_in_mstep
!
! new option in 2021-feb-22
! call expdif before each dynamics (n dycore + 1 tracer) step
!
          if(use_expdif_tstep)then
             if(do_dycore_laplacian_2nd.or.do_dycore_laplacian_4th) call grist_dycore_diffusion_run(mesh, nexpdif_level,ntracer)
             if(do_tracer_laplacian_2nd) call grist_tracer_diffusion_run(mesh, nexpdif_level,ntracer)
          end if
          call tracer_transport_evaluate_mif(mesh,nlev)
          do idstep = 1, dstep_in_tstep
             call dycore_time_integration_run(mesh,dycore_timestep,itimestep,idstep,dstep_in_tstep,itstep,tstep_in_mstep)
             call clock_begin(clock_d2t)
             call grist_d2t_coupling_driver(mesh, nlev, idstep, dstep_in_tstep)
             call grist_dtp_time_average(mesh, nlev, nlevp, idstep, dstep_in_tstep,itstep,tstep_in_mstep) ! time-average www/omega
             call clock_end(clock_d2t)
          end do
          call clock_begin(clock_tracer)
          call tracer_transport_time_integration_run(mesh,tracer_timestep,itimestep,gcm_testcase,itstep,tstep_in_mstep)
          call clock_end(clock_tracer)
       end do
       call clock_end(clock_dynamics)
!
! diagnostics at model step
!
       call clock_begin(clock_diagte)
       if(itimestep .gt. itimestep_pre)then
         call dycore_diagnose_total_energy(mesh)
         call dycore_diagnose_total_energy_moist_e2(mesh)
         call dycore_diagnose_total_energy_moist_e3(mesh)
         call dycore_diagnose_mass(mesh)
         call dycore_diagnose_mass_pt(mesh)
         call tracer_transport_diagnose_global_tracer_mass(mesh)
         call tracer_transport_diagnose_global_crnum(mesh,tracer_timestep)
       end if
       call clock_end(clock_diagte)
!       
! call physics forcing here
!
       call clock_begin(clock_physw)
       if(use_phys)then
          call tracer_transport_evaluate_mif(mesh,nlev) ! newly added for real-physics, does not affect regression old results,
                                                        ! which do not depend on this
          call grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, model_timestep,itimestep, physpkg)
          ! update f2 physics
          if(ptend_f2_on) call grist_dtp_coupling_driver_coup(mesh, nlev, ntracer, model_timestep,physpkg)
       end if
       call clock_end(clock_physw)
!
! diag accumulated vars each model step call gcm_accu_physics_variables
! I/O
!
       call clock_begin(clock_diagio)
       if(write_history_h0) call gcm_output_history_h0(mesh, itimestep) ! monthly-averaged

       if(write_history_h1.and..not.write_history_h1_separate) call gcm_output_history_h1(mesh, itimestep,nsteps,h1_history_freq)
       if(write_history_h1_separate)then
          call gcm_output_history_h1_1d(mesh, itimestep,nsteps,h1_1d_history_freq)
          call gcm_output_history_h1_2d(mesh, itimestep,nsteps,h1_2d_history_freq)
          call gcm_output_history_h1_3d(mesh, itimestep,nsteps,h1_3d_history_freq)
       end if
       call clock_end(clock_diagio)

       if(mpi_rank() == 0) then
         open(1,file=trim(outdir)//'step_stop.txt',status='unknown')
         write(1,*) itimestep
         close(1)
       endif
!
! how to write restart (but the rst files are the same)
!
       call clock_begin(clock_write_rst)
       if(write_restart_h1) call gcm_write_restart_file_h1(mesh,itimestep,nsteps,h1_restart_freq) ! write n-step
       if(write_restart_h0) call gcm_write_restart_file_h0(mesh,itimestep,nsteps)                 ! write each month
       call clock_end(clock_write_rst)

       call cpu_time(step_end)
       step_elapse = step_end-step_beg
       if(mpi_rank() == 0) print*, "elapsed time of this step is ", step_elapse

    END DO

    call clock_end(clock_outloop)

    return
  end subroutine grist_gcm_control_run_amipw

  subroutine grist_gcm_final(mesh)

!gcm
   use grist_gcm_diagnose_h0_module,       only: gcm_h0_diagnose_final
   use grist_gcm_diagnose_h1_module,       only: gcm_h1_diagnose_final
! datam
   use grist_datam_static_data_module,     only: grist_static_data_destruct
! dycore
   use grist_dycore_time_integration_2d,   only: dycore_time_integration_final
   use grist_dycore_vars_module,           only: grist_dycore_vars_destruct
   use grist_dycore_diffusion_module,      only: grist_diffusion_destruct 
! tracer
   use grist_tracer_transport_vars_module, only: tracer_transport_vars_destruct
   use grist_tracer_transport_time_integration_hvsplit, only: tracer_transport_time_integration_final
! dtp
   use grist_dtp_vars_module,              only: grist_dtp_vars_destruct
   use grist_dtp_interface_module,         only: grist_dtp_interface_destruct
   use grist_dtp_interface_full_physpkg1,  only: grist_full_physpkg1_final
   use grist_dtp_interface_full_physpkg2,  only: grist_full_physpkg2_final
! physics
   use grist_physics_data_structure,       only: grist_physics_data_structure_destruct
! infrastructure
   use grist_domain_types,                 only: grist_domain_clean_runtime_allocated
#ifdef LAM_DOMAIN
   use grist_datam_glam_data_module,       only: grist_glam_data_destruct 
#endif

   type(global_domain), intent(inout) :: mesh

! data clean and destruct
     if(write_history_h0) call gcm_h0_diagnose_final
     if(write_history_h1) call gcm_h1_diagnose_final

     call grist_static_data_destruct
! dycore
     call dycore_time_integration_final
     call grist_dycore_vars_destruct
     call grist_diffusion_destruct
! tracer clean
     call tracer_transport_time_integration_final
     call tracer_transport_vars_destruct
! dtp
     call grist_dtp_vars_destruct
     call grist_dtp_interface_destruct
     call grist_full_physpkg1_final
     call grist_full_physpkg2_final
! physics
     call grist_physics_data_structure_destruct
! infrastructure 
     call grist_domain_clean_runtime_allocated(mesh)
#ifdef LAM_DOMAIN
     call grist_glam_data_destruct
#endif
    return
  end subroutine grist_gcm_final

  end module grist_gcm_control_driver
