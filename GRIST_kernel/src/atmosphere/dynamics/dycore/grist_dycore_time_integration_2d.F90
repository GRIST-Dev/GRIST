
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: This is the major time stepping driver for the 3d model in
!              either, hdc or ndc configuration
! 
! Revision history:
!
!              1. DYCORE_OPT is optimized code will change bit-regression,
!              retained due to we want regresssion some old solutions
!
!              2. Refactor to allow optimal implementation for both G-R
!
! References: Zhang et al. (2019), JAMES, doi:10.1029/2018MS001539 (Dry core)
!             Zhang et al. (2020), MWR,   doi:10.1175/MWR-D-19-0305.1 (Moist core)
!             Zhang et al. (2024), QJRMS, doi:https://doi.org/10.1002/qj.4804 (LAM)
!
!----------------------------------------------------------------------------

  module grist_dycore_time_integration_2d

!-----------------------------------------------------------------------
! infrastructure
!-----------------------------------------------------------------------

! mpi
    use grist_mpi
    use grist_domain_types,    only: global_domain
    use grist_data_types,      only: scalar_1d_field, scalar_2d_field, exchange_field_list_2d, exchange_field_list_1d
    use grist_constants,       only: rvap,rdry,ptfactor,cp,p00,zero, i4, r8, gravity, one, half, rearth
    use grist_nml_module,      only: model_timestep, nlev, nlevp, mas_adv_flag, pot_adv_flag, ver_adv_flag, nrk, hor_pgf_flag, &
                                     nh_dynamics, gcm_testcase,  nsteps, &
                                     ad_dycore_laplacian_2nd, ad_dycore_laplacian_4th, ad_dycore_laplacian_6th, &
                                     tend_nct_once, working_mode, use_phys, physpkg, &
                                     ptendSubDiabPhys, ptend_wind_rk_on, ptend_heat_rk_on, ptend_dycore_f1_on, ptend_dycore_heat_f1_on, ptend_dycore_wind_f1_on, &
                                     write_stepinfo, write_verbose, use_www_hyperDiffusion, &
                                     restore_hydro, restore_hydro_minsteps,restore_hydro_intsteps, adjphi
#ifndef SEQ_GRIST
    use grist_config_partition,       only: debug_data_1d, debug_data_2d,&
                                            exchange_data_2d_add, exchange_data_2d,&
                                            exchange_data_1d_add, exchange_data_1d
    use grist_clocks,                 only: clock_id, clock_begin, clock_end
#endif

!-----------------------------------------------------------------------
! atmosphere
!-----------------------------------------------------------------------

! dycore&tracer data
    use grist_dycore_vars_module,              only: dycoreVarCellFace, dycoreVarCellFull, dycoreVarEdgeFull, dycoreVarEdgeFace, dycoreVarSurface
    use grist_tracer_transport_vars_module,    only: tracerVarCellFace, tracerVarEdgeFull, tracerVarCellFull
! horizontal module
    use grist_dycore_primal_flux_operators_2d, only: calc_primal_normal_flux_at_edge
    use grist_dycore_hori_swe_module_2d,       only: calc_tend_nct_at_edge_full_level, calc_grad_kinetic_energy
    use grist_dycore_ref_atmos,                only: grist_dycore_ref_atmos_create_ptb, &
                                                     scalar_alphad_at_pc_full_level_ptb, &
                                                     scalar_pressure_at_pc_full_level_ptb, &
                                                     scalar_geop_at_pc_full_level_ptb, &
                                                     scalar_delp_at_pc_full_level_ptb, &
                                                     scalar_grad_hpres_at_edge_full_level_bar
! hpe 
    use grist_hpe_constants,          only: deta_full, deta_face,eta_full,eta_full_a,eta_full_b,eta_face_a,eta_face_b
    use grist_hpe_continuity,         only: calc_hpe_tend_continuity_2d
    use grist_hpe_vertical_advection, only: calc_hpe_vert_advection, calc_hpe_tend_vert_mass_flux_2d
    use grist_hpe_hydro_pgf,          only: calc_hpe_hydro, calc_hpe_get_full_mass
#ifdef DYCORE_NHD
! now nhd is add-on
    use grist_nh_driver_module,       only: grist_nh_dynamics_run, ndc_restore_flag
#endif
! physics
    use grist_physics_data_structure, only: ptend_f1, ptend_rk, ptend_f2
#ifdef AMIPW_PHYSICS
    use grist_PhysW_nml_module,     only: step_cu
#endif
#ifdef LAM_DOMAIN
    use grist_datam_glam_data_module, only: update_assignValueArea_1d, update_assignValueArea_2d
#endif

    implicit none

    private
    public :: dycore_time_integration_init,&
              dycore_time_integration_run ,&
              dycore_time_integration_final
!
! local np1 var (prognostic) 
!
     type(scalar_2d_field)               :: scalar_www_at_pc_face_level_np1               ! halo data exchange
     type(scalar_2d_field)               :: scalar_phi_at_pc_face_level_np1               ! halo data exchange
     type(scalar_1d_field)               :: scalar_hpressure_at_pc_surface_np1            ! halo data exchange
     type(scalar_2d_field)               :: scalar_normal_velocity_at_edge_full_level_np1 ! halo data exchange
     type(scalar_2d_field)               :: scalar_mass_pt_at_pc_full_level_np1           ! halo data exchange
!
! local np1 var (diagnostic)
!
     type(scalar_2d_field)               :: scalar_delhp_at_pc_face_level_np1     ! used by nhd
     type(scalar_2d_field)               :: scalar_hpressure_at_pc_full_level_np1
     type(scalar_2d_field)               :: scalar_hpressure_at_pc_face_level_np1
!
! clock
!
     integer(i4)                         :: clock_rkl, clock_mas, clock_nct, clock_vau, clock_ket, clock_diagom, &
                                            clock_nhd, clock_ptm, clock_pgf, clock_mainexch, clock_fnlrk

     real(r8), parameter :: damp_kesi = -1 ! this is a regression

   CONTAINS

!-----------------------------------------------------------------------------
! Name convention:
! for prognostic variables
! q(n+1) = q(n)+F(rk)
! n+1, new state after each sub stage
! n, state at begining of each time step
! rk itermediate state used for explicit tendency computin
!
! for other variable, _n denotes an intermediate state
! np1 denotes next (if needed)
!------------------------------------------------------------------------------

   subroutine dycore_time_integration_run(mesh, dtime, itimestep,idstep,dstep_in_tstep,itstep,tstep_in_mstep)
!
! io
!
     use omp_lib
     type(global_domain),  intent(inout) :: mesh
     real(r8)          ,   intent(in)    :: dtime
     integer(i4)       ,   intent(in)    :: itimestep
     integer(i4), optional,intent(in)    :: idstep, dstep_in_tstep, itstep, tstep_in_mstep
!
! variables used within RK step, local rk
!
     type(scalar_2d_field)               :: scalar_delhp_at_pc_full_level_rk
     type(scalar_2d_field)               :: scalar_normal_velocity_at_edge_full_level_rk
     type(scalar_2d_field)               :: scalar_www_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_phi_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_delhp_at_pc_face_level_rk      ! used by nhd
!
! local
!
     real(r8)                            :: rk_substep, cr_num, cr_act, cr_sum, eta, kv 
     real(r8)                            :: div_sum(nlev), tmp1(nlev), tmp2(nlev), tmp3(nlev), tmp4(nlev)
     real(r8)                            :: scalar_mpressure_at_pc_full_level_old(nlev,mesh%nv_full)
     real(r8)                            :: scalar_pressure_at_pc_face_level_old(nlevp,mesh%nv_full)
     real(r8)                            :: scalar_pressure_at_pc_full_level_old(nlev,mesh%nv_full)
#ifdef AMIPW_PHYSICS
     real(r8)                            :: scalar_pt_at_pc_full_level_old(nlev,mesh%nv_full)
#endif
     integer(i4)                         :: irk_step, rk_number ! number to divide DT in RK step
     integer(i4)                         :: it, ie, iv, ilev, inb, icell1, icell2, index_edge, ncell_do
     integer(i4)                         :: iblock
! template
     real(r8)                            :: scalar_template_1d_nlev_a(nlev)
     real(r8)                            :: scalar_template_1d_nlev_b(nlev)
     real(r8)                            :: scalar_template_1d_nlev_c(nlev)
     real(r8)                            :: scalar_template_1d_nlevp_a(nlev+1)
     real(r8)                            :: scalar_template_1d_nlevp_b(nlev+1)
     real(r8)                            :: scalar_template_a
#ifndef SEQ_GRIST
     type(exchange_field_list_2d),pointer:: field_head_2d
     type(exchange_field_list_1d),pointer:: field_head_1d
#endif
      
      field_head_2d=>null()
      field_head_1d=>null()
      iblock = mpi_rank()

!-----------------------------------------------------------------------
! added for cumulus
!-----------------------------------------------------------------------

#ifdef AMIPW_PHYSICS
! only at the first dycore and tracer step
      IF(present(idstep))then
         if(idstep.eq.1.and.itstep.eq.1.and.mod(itimestep-1,step_cu).eq.0) dycoreVarCellFull%tend_pt_n%f = zero
      END IF
#endif
!
! Prepare for dycore integration
!
#ifdef DYCORE_NHD
      if(nh_dynamics)then
         scalar_www_at_pc_face_level_rk            = dycoreVarCellFace%scalar_www_n
         scalar_phi_at_pc_face_level_rk            = dycoreVarCellFace%scalar_phi_n
         dycoreVarCellFull%scalar_pressure_rk      = dycoreVarCellFull%scalar_pressure_n
         scalar_delhp_at_pc_face_level_rk          = dycoreVarCellFace%scalar_delhp_n
      end if
#endif
      scalar_normal_velocity_at_edge_full_level_rk = dycoreVarEdgeFull%scalar_normal_velocity_n
      scalar_delhp_at_pc_full_level_rk             = dycoreVarCellFull%scalar_delhp_n
      scalar_mpressure_at_pc_full_level_old        = dycoreVarCellFull%scalar_mpressure_n%f

#ifdef AMIPW_PHYSICS
! raw pt
      scalar_pt_at_pc_full_level_old = dycoreVarCellFull%scalar_potential_temp_n%f/(one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:))
#endif

!-----------------------------------------------------------------------
! we have commented if statement assocaited stencil_exchange_flag,
! stencil_exchange_flag=.true. is default now (-DUSE_HALO2)
!-----------------------------------------------------------------------

#ifndef SEQ_GRIST
     call clock_begin(clock_rkl)
#endif

!----------------------------------------------------------------------
! The overall workflow of dycore solver (Z19, JAMES):
! predictor-corrector (rk2/3) plus forward-backward (mesinger-arakawa)
!
! (1) Solve mass continuity equation forward in time and do some necessary diagnostics.
! (2) Solve mass-pt equation forward in time and do some necessary update (e.g., diffusion, physics).
!     In future, the above two procedures can be further combined to be more unified.
!
!     mass_pt needed by ndc, pt needed by hdc, so halo exchange at least for halo1 are required.
!     in the old implementation (OLD_MASS_PT), we first exchange pt flux, then do update until nv_halo1. This is now changed
!     to directly do update until nv_compute and do exchange. This removes later mass_pt exchange for this substep.
!     regresssion same for global modeling. For LAM, doing so affects some halo-177 cell at the data zone, as we have
!     put lam_update_halo before exchange_halo, so not exactly regression, but self-par consistency ok.
!
! (3) Calculate the forward-part tendency of velocity equation (nct, grad-ke, vadv, utill ne_compute).
!     As all advection terms are forward in time, they are consistent.
!
! (3.1) Diagnose OMEGA (mpressure-based vertial speed).
!
! (4) Solve ndc equation (prognose www&phi) or solve hdc equation (diagnose) to obtain renewed geopotential.
!  old (NHD_NOOPT): www&phi are solved until nv_halo1 by exchange edge-based flux data, then evaluting hpgf until ne_compute for Un.
!  new: solve them until nv_compute, then do cell-based exchange for them, then evaluting hpgf until ne_compute for Un.
!  for ndc, full air pressure is from gas law; for hdc, full air pressure is moist mass pressure
!
! (5) Calculate the backward-part tendency of velocity equation (hpgf, ne_compute) and do an update.
!
! Note:
!
! (1) Diffusion, physics tend entry can be slightly flexible, but we fix the current
! update sequence for bitwise regressing old-code solutions (a check for code refactoring)
!
! (2) Some code refactoring (OLD_MASS_PT, NHD_NOOPT), while regression (bitwise) ok for global, is not so for LAM (e.g., 
!   garbage numbers may occur at data zone for some uninit vars). Be carefull and understand whether changes are acceptable.
!   Overall, global modeling is our 1st guiding principle (design, implementation) for LAM, i.e., LAM has to follow, as long as 
!   some desiable consistencies are guarranteed (e.g., G-R, par-self, restart).
!   Bespoke, optimization and added value of LAM, they are another story (see words in Z24, QJRMS).
!----------------------------------------------------------------------

     DO irk_step = 1, nrk

        rk_number  = nrk+1-irk_step
        rk_substep = dtime/rk_number
!(1)
#ifndef SEQ_GRIST
        call clock_begin(clock_mas)
#endif
        call dycore_solve_mass_eq_diag
#ifndef SEQ_GRIST
        call clock_end(clock_mas)
#endif

!(2)
#ifndef SEQ_GRIST
        call clock_begin(clock_ptm)
#endif
        call dycore_solve_mass_pt_eq
#ifndef SEQ_GRIST
        call clock_end(clock_ptm)
#endif

!(3)
        call dycore_velocityeq_forward

#ifndef SEQ_GRIST
        call clock_begin(clock_diagom)
#endif
        call dycore_diagnose_omega
#ifndef SEQ_GRIST
        call clock_end(clock_diagom)
#endif

!----------------------------------------------------------------------
! Until this line, HYDROSTATIC and NONHYDROSTATIC is same
!----------------------------------------------------------------------

!(4)
#ifndef SEQ_GRIST
        call clock_begin(clock_nhd) ! either nhd or hpe
#endif

#ifdef DYCORE_NHD
        IF(nh_dynamics)then

           if(iblock .eq. 0.and.write_verbose) print*,"before nh dynamics"

           dycoreVarCellFace%scalar_delp_n%f(2:nlev,:)  = dycoreVarCellFull%scalar_pressure_rk%f(2:nlev,:)-dycoreVarCellFull%scalar_pressure_rk%f(1:nlev-1,:)
           dycoreVarCellFace%scalar_delp_n%f(1,:)       = dycoreVarCellFull%scalar_pressure_rk%f(1,:)     -dycoreVarCellFace%scalar_pressure_n%f(1,:)
           dycoreVarCellFace%scalar_delp_n%f(nlev+1,:)  = dycoreVarCellFace%scalar_pressure_n%f(nlev+1,:) -dycoreVarCellFull%scalar_pressure_rk%f(nlev,:)
   
#ifdef GRIST_DEBUG
           call debug_data_2d(irk_step,1,mesh%e_index,mesh%ne,"o1o",dycoreVarEdgeFace%scalar_normal_mass_flux_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o2o",dycoreVarCellFace%tend_mass_hori%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o3o",dycoreVarCellFull%scalar_eta_mass_flux_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o4o",scalar_phi_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o5o",scalar_www_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o6o",dycoreVarCellFace%scalar_delp_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o7o",dycoreVarCellFace%scalar_delhp_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o8o",dycoreVarCellFull%scalar_pressure_rk%f)
           call debug_data_1d(irk_step,  mesh%v_index,mesh%nv,"o9o",dycoreVarSurface%scalar_geopotential_n%f)!1d
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o10o",dycoreVarCellFull%scalar_pressure_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o11o",dycoreVarCellFace%scalar_phi_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o12o",dycoreVarCellFull%scalar_delhp_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o13o",dycoreVarCellFace%scalar_www_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o14o",dycoreVarCellFull%scalar_mass_pt_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o15o",dycoreVarCellFull%scalar_delhp_np1%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o16o",scalar_delhp_at_pc_face_level_np1%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o17o",scalar_mass_pt_at_pc_full_level_np1%f)
#endif

           call grist_nh_dynamics_run(mesh,nlev,nlevp,rk_substep, itimestep, irk_step, idstep, itstep, &
                                         dycoreVarSurface%tend_hpressure_cnty,          & ! time tendency of hps
                                         scalar_hpressure_at_pc_face_level_np1,         & ! only need hp top
                                         dycoreVarEdgeFace%scalar_normal_mass_flux_n,   & ! rk, follow comment, donot fooled by vars' name which is for simplicity
                                         dycoreVarCellFace%tend_mass_hori,              & ! rk
                                         dycoreVarCellFull%scalar_eta_mass_flux_n,      & ! rk
                                         dycoreVarCellFace%scalar_eta_mass_flux_n,      & ! rk
                                         scalar_phi_at_pc_face_level_rk,                & ! rk
                                         scalar_www_at_pc_face_level_rk,                & ! rk
                                         dycoreVarCellFace%scalar_delp_n,               & ! rk
                                         scalar_delhp_at_pc_face_level_rk,              & ! rk
                                         dycoreVarCellFace%scalar_delhp_n,              & ! n
                                         dycoreVarCellFull%scalar_pressure_rk,          & ! rk
                                         dycoreVarSurface%scalar_geopotential_n,        & ! n
                                         dycoreVarCellFull%scalar_pressure_n,           & ! n
                                         dycoreVarCellFace%scalar_phi_n,                & ! n
                                         dycoreVarCellFull%scalar_delhp_n,              & ! n
                                         dycoreVarCellFace%scalar_www_n,                & ! n
                                         dycoreVarCellFull%scalar_mass_pt_n,            & ! n
                                         dycoreVarCellFull%scalar_delhp_np1,            & ! np1
                                         scalar_delhp_at_pc_face_level_np1,             & ! np1
                                         scalar_mass_pt_at_pc_full_level_np1,           & ! np1
                                         tracerVarCellFace%scalar_mif_n,                & ! tracer's n
                                         tracerVarCellFull%scalar_mif_n,                & ! tracer's n
                                         tracerVarCellFull%scalar_tracer_mxrt_n,        & ! tracer's n
                                         scalar_phi_at_pc_face_level_np1,               & ! np1, out
                                         scalar_www_at_pc_face_level_np1,               & ! np1, out
                                         dycoreVarCellFull%scalar_temp_n,               & ! np1, out!
                                         dycoreVarCellFull%scalar_geopotential_n,       & ! np1, out!
                                         dycoreVarCellFace%scalar_geopotential_n,       & ! np1, out
                                         dycoreVarCellFull%scalar_pressure_np1,         & ! np1, out
                                         dycoreVarCellFace%scalar_pressure_n,           & ! np1, out!
                                         dycoreVarCellFull%scalar_delp_np1,             & ! np1, out
                                         dycoreVarCellFull%scalar_alpha_np1)

#ifdef GRIST_DEBUG
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x1x",scalar_phi_at_pc_face_level_np1%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x2x",scalar_www_at_pc_face_level_np1%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x3x",dycoreVarCellFull%scalar_temp_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv_halo(1),"x4x",dycoreVarCellFull%scalar_geopotential_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x5x",dycoreVarCellFace%scalar_geopotential_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x6x",dycoreVarCellFull%scalar_pressure_np1%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv_halo(1),"x7x",dycoreVarCellFace%scalar_pressure_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv_halo(1),"x8x",dycoreVarCellFull%scalar_delp_np1%f)
#endif
          if(iblock .eq. 0.and.write_verbose) PRINT*,"finish NDC equation"

!-----------------------------------------------------------------------
! Optionally, to obtain a hydrostatic balance state after minsteps
! for every-intsteps, or based on dynamic control
! note for ndc, this will only modify phi, so do this for nv_full to
! avoid re-exchange phi. This can be sure because mass_pt, delhp have
! been exchanged, and mif remains unchanged within dycore. If only compute to
! halo1, then we need exchange phi after this adjust, which is undesiable
!
! Bitwise regression ok as an old implementation for global
!-----------------------------------------------------------------------

           IF((restore_hydro.and.itimestep.gt.restore_hydro_minsteps.and.mod(itimestep,restore_hydro_intsteps).eq.0).or.& ! mannually control
               ndc_restore_flag.eq.1)then ! dynamic control
              if(iblock.eq.0.and.write_stepinfo) print*,"instant restored at step,", itimestep

!-----------------------------------------------------------------------
! verified: if activate this under ndc every timestep, it produces 
! exactly identical results as use a HDC, except some diagnostic
! variables and www (i.e., code refactoring ok for ndc also works for hdc)
!-----------------------------------------------------------------------

              call hydrostatic_adjust
           ENDIF

           IF(adjphi.and.itimestep.lt.16) call hydrostatic_adjust   ! legacy for regresssion (random set)
#ifdef REAL_NDC
! if strong init imbalance occurs, adjphi for less than 10 steps is typically more than enough
           IF(adjphi.and.itimestep.lt.int(3600/model_timestep)) call hydrostatic_adjust ! may formal set to 1h in real case
#endif
        ELSE
#endif
           call hydrostatic_adjust
           if(iblock .eq. 0.and.write_verbose) PRINT*,"finish HDC equation"
#ifdef DYCORE_NHD
        END IF
#endif

#ifndef SEQ_GRIST
        call clock_end(clock_nhd)
#endif

!(5)
#ifndef SEQ_GRIST
      call clock_begin(clock_pgf)
#endif
      call dycore_velocityeq_hpgf
#ifndef SEQ_GRIST
      call clock_end(clock_pgf)
      call clock_begin(clock_fnlrk)
#endif

        scalar_normal_velocity_at_edge_full_level_np1%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f + dtime/rk_number*&
       ( dycoreVarEdgeFull%tend_normal_velocity_vadv%f + dycoreVarEdgeFull%tend_normal_velocity_pgf%f  +&
         dycoreVarEdgeFull%tend_normal_velocity_ke%f   + dycoreVarEdgeFull%tend_normal_velocity_nct%f )

         if(use_phys.and.ptend_wind_rk_on)then
            scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime/rk_number*ptend_rk%tend_normal_velocity_at_edge_full_level%f
         end if

!-----------------------------------------------------------------------
! diffusion terms added here as phys, must be written in this way for
! bitwise regression (not a correctness requirement)
!-----------------------------------------------------------------------

         if(ad_dycore_laplacian_2nd.and. .not.ad_dycore_laplacian_4th)then
            scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_2nd%f
         else if(ad_dycore_laplacian_4th .and. .not. ad_dycore_laplacian_2nd)then
            scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_4th%f
         else if(ad_dycore_laplacian_2nd.and.ad_dycore_laplacian_4th)then
            scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime/rk_number*(dycoreVarEdgeFull%tend_hwind_laplacian_4th%f+&
                                                                    dycoreVarEdgeFull%tend_hwind_laplacian_2nd%f)
         end if

#ifdef DYCORE_NHD

!-----------------------------------------------------------------------
! add Lap4th for www equation in nh_dynamics
!-----------------------------------------------------------------------

         if(nh_dynamics.and.ad_dycore_laplacian_4th.and.use_www_hyperDiffusion)then
            scalar_www_at_pc_face_level_np1%f = scalar_www_at_pc_face_level_np1%f+dtime/rk_number*dycoreVarCellFace%tend_www_laplacian_4th%f
         end if

#if (!defined DCMIP21)
#ifndef REG_AMP21
         if(nh_dynamics.and.ad_dycore_laplacian_2nd)then
             scalar_www_at_pc_face_level_np1%f = scalar_www_at_pc_face_level_np1%f+dtime/rk_number*dycoreVarCellFace%tend_www_laplacian_2nd%f
         end if
#endif
#endif

#endif

#ifdef GRIST_DYCORE_TEST
! diff 6th is only for testing, not used for all cases
         if(ad_dycore_laplacian_6th)then
             scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_6th%f
         end if
#endif

!-----------------------------------------------------------------------
! RK-sub-step update
!-----------------------------------------------------------------------

#ifdef GRIST_DEBUG
       if(irk_step .eq. 1 .and. itimestep .eq. int(nsteps)) then
       call debug_data_1d(irk_step,      mesh%v_index,mesh%nv,"shaps_n" ,scalar_hpressure_at_pc_surface_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"smpapf_n",scalar_mass_pt_at_pc_full_level_np1%f)
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"swapf_n" ,scalar_www_at_pc_face_level_np1%f)
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"sdpapf_n",scalar_phi_at_pc_face_level_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%e_index,mesh%ne,"snvaef_n",scalar_normal_velocity_at_edge_full_level_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"sdhapf_n",dycoreVarCellFull%scalar_delhp_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"spapf_n" ,dycoreVarCellFull%scalar_pressure_np1%f)      
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"diag_1"  ,scalar_delhp_at_pc_face_level_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"diag_2"  ,dycoreVarCellFull%scalar_potential_temp_n%f) !step 2 error
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"diag_3"  ,scalar_hpressure_at_pc_face_level_np1%f)
       end if
#endif

#ifdef LAM_DOMAIN
       call update_assignValueArea_2d(mesh,nlev,"normal_velocity",scalar_normal_velocity_at_edge_full_level_np1)
#ifdef OLD_MASS_PT
       call update_assignValueArea_2d(mesh,nlev,"mass_pt",scalar_mass_pt_at_pc_full_level_np1)
#endif
#endif

#ifndef SEQ_GRIST
       call clock_begin(clock_mainexch)
#ifdef OLD_MASS_PT
       call exchange_data_2d_add(mesh,field_head_2d,scalar_mass_pt_at_pc_full_level_np1)
#endif
       call exchange_data_2d_add(mesh,field_head_2d,scalar_normal_velocity_at_edge_full_level_np1)
#endif

#ifdef DYCORE_NHD
       if(nh_dynamics)then
#ifdef LAM_DOMAIN
!-----------------------------------------------------------------------
! as long as halo region being modified, we should reset halo bdy, thou 
! no need to exchange, below are due to www diffusion and hydro_adjust phi
!-----------------------------------------------------------------------
          call update_assignValueArea_2d(mesh,nlevp,"wwwFace", scalar_www_at_pc_face_level_np1)
          call update_assignValueArea_2d(mesh,nlevp,"phiFace", scalar_phi_at_pc_face_level_np1)
#endif

#ifndef SEQ_GRIST

!-----------------------------------------------------------------------
! exchange due to www 2nd diffusion added before, if do opt (put diffusion
! tend before the previous www data exchange, we can remove this line 
! but bit-solutions will change (physically still correct), as www from ndc
! will affect wind
!-----------------------------------------------------------------------

          call exchange_data_2d_add(mesh,field_head_2d,scalar_www_at_pc_face_level_np1)

!-----------------------------------------------------------------------
! in this new cell-based NDC data exchange, see(nhd), if hydro_adjust donot
! compute until nv_full, we may need this as it modifies halo, but now to
! nv_full, so no need for exchange note that LAM's bitwise regression is
! more sensitive to exchange pattern because bdy-halo does not exchange data
!-----------------------------------------------------------------------

#ifdef NHD_NOOPT
          call exchange_data_2d_add(mesh,field_head_2d,scalar_phi_at_pc_face_level_np1)
#endif

#endif
       end if
#endif

#ifndef SEQ_GRIST
       call exchange_data_2d(mesh%local_block,field_head_2d)
       call clock_end(clock_mainexch)
#endif

#ifdef OLD_MASS_PT
       dycoreVarCellFull%scalar_potential_temp_n%f = scalar_mass_pt_at_pc_full_level_np1%f/&  ! Overwritten each RK step
                                                     dycoreVarCellFull%scalar_delhp_np1%f
#endif

! prog
#ifdef DYCORE_NHD
      if(nh_dynamics)then
         scalar_www_at_pc_face_level_rk            = scalar_www_at_pc_face_level_np1
         scalar_phi_at_pc_face_level_rk            = scalar_phi_at_pc_face_level_np1
         dycoreVarCellFull%scalar_pressure_rk      = dycoreVarCellFull%scalar_pressure_np1
! diag
         scalar_delhp_at_pc_face_level_rk          = scalar_delhp_at_pc_face_level_np1
      end if
#endif
      scalar_normal_velocity_at_edge_full_level_rk = scalar_normal_velocity_at_edge_full_level_np1
      scalar_delhp_at_pc_full_level_rk             = dycoreVarCellFull%scalar_delhp_np1

      if(iblock .eq. 0 .and. write_stepinfo) print*,"---------- rkfb ",irk_step,"----------"
#ifndef SEQ_GRIST
      call clock_end(clock_fnlrk)
#endif

   END DO
#ifndef SEQ_GRIST
   call clock_end(clock_rkl)
#endif

!----------------------------------------------
! final update for the next dycore step
!----------------------------------------------

#ifdef DYCORE_NHD
      if(nh_dynamics)then
         dycoreVarCellFace%scalar_www_n      = scalar_www_at_pc_face_level_np1
         dycoreVarCellFace%scalar_phi_n      = scalar_phi_at_pc_face_level_np1
         dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_pressure_np1
         dycoreVarCellFace%scalar_delhp_n    = scalar_delhp_at_pc_face_level_np1
      end if
#endif
      dycoreVarEdgeFull%scalar_normal_velocity_n  = scalar_normal_velocity_at_edge_full_level_np1
      dycoreVarCellFull%scalar_mass_pt_n          = scalar_mass_pt_at_pc_full_level_np1
      dycoreVarSurface%scalar_hpressure_n         = scalar_hpressure_at_pc_surface_np1
      dycoreVarCellFull%scalar_delhp_n            = dycoreVarCellFull%scalar_delhp_np1    ! needed by physics and dtp

!-----------------------------------------------------------------------
! this is only for diagnose PS, not put to dycore_diag because of 
! historry-h1 does not dycore_diag at each step
!-----------------------------------------------------------------------

      dycoreVarSurface%scalar_pressure_n%f = dycoreVarCellFace%scalar_pressure_n%f(nlevp,:)

!-----------------------------------------------------------------------
! substract diabatic physical influence, used with F2, so we should substract F2
! here; for idealized test, f2 has the same tendency, while for full-physics
! test, rk has its own total tendency, only the portion that f2 covers should be
! substracted
!-----------------------------------------------------------------------

      if(ptendSubDiabPhys)then
          dycoreVarCellFull%scalar_mass_pt_n%f =  dycoreVarCellFull%scalar_mass_pt_n%f-&
                                                 !dtime*ptend_rk%tend_mass_pt_at_pc_full_level%f
                                                 dtime*ptend_f2%tend_mass_pt_at_pc_full_level%f
#ifndef SEQ_GRIST
          call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
          call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
          dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/& ! overwritten each RK step
                                                       dycoreVarCellFull%scalar_delhp_n%f
      end if

!-----------------------------------------------------------------------
! pt now only affected by adv and diffusion, get this tend and accumulate until this itimestep
! reachs step_cu and all dycore+tracer steps done this tend is only for cumulus tiedtke
!-----------------------------------------------------------------------

#ifdef AMIPW_PHYSICS
      IF(present(idstep))then
      dycoreVarCellFull%tend_pt_n%f = dycoreVarCellFull%tend_pt_n%f+&
                                    (dycoreVarCellFull%scalar_potential_temp_n%f/(one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:))-&
                                     scalar_pt_at_pc_full_level_old)/dtime
      if(idstep.eq.dstep_in_tstep.and.itstep.eq.tstep_in_mstep.and.mod(itimestep,step_cu).eq.0)then
         dycoreVarCellFull%tend_pt_n%f = dycoreVarCellFull%tend_pt_n%f/(dstep_in_tstep*tstep_in_mstep*step_cu)
      end if
      END IF
#endif


#ifdef USE_PTENDF1

!-----------------------------------------------------------------------
!                          Coupling of F1 physics forcing
! after each dycore step, _n is the updated state, on which model physics reside
! Recent tests suggest that F1 physics for dycore is less useful,
! For HDC, just use ptend_f2_sudden;
! For NDC, use ptend_rk for wind, heat (slow-physics), use ptend_f2_sudden
! for tracer update, and some OS heat and wind (fast-physics)
! These codes are still here for reference purpose, but we added a compipling option.
! So as in tracer transport' F1
! 2021-01-05: Based on many full-physics tests, we can basically ignore F1-tend
! can be replaced by ptend_rk, ptend_f2 and f-s combination
! leave here for futher update and regression
!-----------------------------------------------------------------------

      if(use_phys.and.ptend_dycore_heat_f1_on)then
         dycoreVarCellFull%scalar_mass_pt_n%f           = dycoreVarCellFull%scalar_mass_pt_n%f+&
                                                         ptend_f1%tend_mass_pt_at_pc_full_level%f*dtime
#ifndef SEQ_GRIST
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
         call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
         dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
      end if

      if(use_phys.and.ptend_dycore_wind_f1_on)then
         dycoreVarEdgeFull%scalar_normal_velocity_n%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f+&
                                                         ptend_f1%tend_normal_velocity_at_edge_full_level%f*dtime
#ifndef SEQ_GRIST
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
         call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
      end if

      if(use_phys.and.ptend_dycore_f1_on)then
         dycoreVarCellFull%scalar_mass_pt_n%f           = dycoreVarCellFull%scalar_mass_pt_n%f+&
                                                         ptend_f1%tend_mass_pt_at_pc_full_level%f*dtime
         dycoreVarEdgeFull%scalar_normal_velocity_n%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f+&
                                                         ptend_f1%tend_normal_velocity_at_edge_full_level%f*dtime
#ifndef SEQ_GRIST
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
         call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
         dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
      end if
#endif

      RETURN

CONTAINS

!-----------------------------------------------------------------------
!SEGMENT1
!-----------------------------------------------------------------------

    subroutine dycore_solve_mass_eq_diag

!-----------------------------------------------------------------------
! compute normal mass flux {DELP*V} and mass flux div GRAD.(DELP*V)} 
! at each full level
!-----------------------------------------------------------------------

        mesh%ne = mesh%ne_halo(1)
! inline_opt_here
! normal mass flux is new only until ne_halo1
! scalar_normal_mass_flux_n is also latest for a dycore step (i.e., a complte RK cycle)
! it does not use exchange, this why we need an exchange within dtp

!$omp parallel  private(ie,ilev)
!$omp do schedule(static,100)
        do ie = 1, mesh%ne_halo(1)
           do ilev = 1, nlev
              dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev, ie) = 0.5*(scalar_delhp_at_pc_full_level_rk%f(ilev,mesh%edt_v(1,ie))+&
                                                                             scalar_delhp_at_pc_full_level_rk%f(ilev,mesh%edt_v(2,ie)))
           end do
        end do
!$omp end do nowait
!$omp end parallel 
! inline_opt_here
        dycoreVarEdgeFull%scalar_normal_mass_flux_n%f = dycoreVarEdgeFull%scalar_normal_mass_flux_n%f*scalar_normal_velocity_at_edge_full_level_rk%f

!-----------------------------------------------------------------------
! exchange scheme, needed if o3 is used, but o2 is used so hard coded here!
! call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_mass_flux_n)
! call exchange_data_2d(mesh%local_block,field_head_2d)
! mesh%ne = mesh%ne_halo(1)
!-----------------------------------------------------------------------
        mesh%nv = mesh%nv_halo(1)

!-----------------------------------------------------------------------
! GRAD.(delhp.V) new until cell halo1
!-----------------------------------------------------------------------

!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do schedule(dynamic,5)
        do iv = 1, mesh%nv
           div_sum(1:nlev) = zero
           do inb = 1, mesh%vtx_nnb(iv)
              index_edge  = mesh%vtx_ed(inb,iv)
              do ilev = 1, nlev
#ifndef DYCORE_OPT
                 div_sum(ilev)  =  div_sum(ilev)+dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
#else
                 div_sum(ilev)  =  div_sum(ilev)+dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*mesh%edp_leng(index_edge)
#endif
              end do
           end do
#ifndef DYCORE_OPT
           dycoreVarCellFull%tend_mass_hori%f(1:nlev,iv) = div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))
#else
           dycoreVarCellFull%tend_mass_hori%f(1:nlev,iv) = div_sum(1:nlev)/(rearth*mesh%plg_areag(iv))
#endif
        end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute normal mass flux",mas_adv_flag

        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute

!-----------------------------------------------------------------------
! compute PS tendency and face-level vertical mass flux (m*etadot)
! Z19 (JAMES, Eq. 26,27)
! new until hps tend and eta_mass_flux are new until cell halo1
!-----------------------------------------------------------------------

        call calc_hpe_tend_continuity_2d(mesh%nv_halo(1), nlev                      , &
                                         dycoreVarCellFull%tend_mass_hori%f         , & ! in
                                         dycoreVarSurface%tend_hpressure_cnty%f     , & ! out
                                         dycoreVarCellFace%scalar_eta_mass_flux_n%f)    ! out

        dycoreVarCellFull%scalar_eta_mass_flux_n%f(1:nlev,:) = 0.5_r8*dycoreVarCellFace%scalar_eta_mass_flux_n%f(1:nlev,:)+&
                                                               0.5_r8*dycoreVarCellFace%scalar_eta_mass_flux_n%f(2:nlev+1,:)

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute PS and face level V mass flux",nh_dynamics

!-----------------------------------------------------------------------
! diag FACE LEVEL MASS FLUX VARS
! Z19 (JAMES, Eq.28)
! new until cell halo1
!-----------------------------------------------------------------------

#ifdef DYCORE_NHD

        IF(nh_dynamics)then

!$omp parallel  private(ilev)
!$omp do schedule(static,5)
        do ilev = 2, nlev
             dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(ilev,:) = &
                                     0.5*(deta_full(ilev-1)/deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,:)+&
                                          deta_full(ilev)  /deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev-1,:))
        end do
!$omp end do nowait
!$omp end parallel

        dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(1,:)     = 0.5_r8*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(1,:)
        dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(nlevp,:) = 0.5_r8*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(nlev,:)

        mesh%nv = mesh%nv_halo(1)  ! only used by nh_dynamics originally, now also used for diagnosing dpidry/dt
        
        !call divergence_operator_2d(mesh,dycoreVarEdgeFace%scalar_normal_mass_flux_n%f, &
        !                                 dycoreVarCellFace%tend_mass_hori%f, & ! face-level adv
        !                                 nlev+1)
!
! inline_opt of divergence_operator_2d below
!

!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do schedule(dynamic,5)
        do iv = 1, mesh%nv
           div_sum(1:nlevp) = zero
           do inb = 1, mesh%vtx_nnb(iv)
              index_edge  = mesh%vtx_ed(inb, iv)
              do ilev = 1, nlevp
#ifndef DYCORE_OPT
                 div_sum(ilev)  = div_sum(ilev) + dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(index_edge)
#else
                 div_sum(ilev)  = div_sum(ilev) + dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb,iv)*mesh%edp_leng(index_edge)
#endif
              end do
           end do
#ifndef DYCORE_OPT
           dycoreVarCellFace%tend_mass_hori%f(1:nlevp,iv) = div_sum(1:nlevp)/((rearth**2)*mesh%plg_areag(iv))
#else
           dycoreVarCellFace%tend_mass_hori%f(1:nlevp,iv) = div_sum(1:nlevp)/(rearth*mesh%plg_areag(iv))
#endif
        end do
!$omp end do nowait
!$omp end parallel

        END IF
#endif

!-----------------------------------------------------------------------
! end of face-level continuty eq, needed by ndc
!-----------------------------------------------------------------------
        
        mesh%nv = mesh%nv_compute ! reset

!-----------------------------------------------------------------------
! update hpressure_surface here, new-time value till halo1, but still 
! exchange because we want all
!-----------------------------------------------------------------------

        scalar_hpressure_at_pc_surface_np1%f = dycoreVarSurface%scalar_hpressure_n%f + dtime/rk_number*dycoreVarSurface%tend_hpressure_cnty%f

#ifndef SEQ_GRIST
        call exchange_data_1d_add(mesh,field_head_1d,scalar_hpressure_at_pc_surface_np1)
        call exchange_data_1d(mesh%local_block,field_head_1d)
#endif

!-----------------------------------------------------------------------
! sequence exchange_data and update_assignValueArea does not matter for LAM (as it should be)
!-----------------------------------------------------------------------

#ifdef LAM_DOMAIN
        call update_assignValueArea_1d(mesh,"hpressure_surface",scalar_hpressure_at_pc_surface_np1)
#endif

!-----------------------------------------------------------------------
! renew mass/hpressure state with new surface hpressure
!-----------------------------------------------------------------------

        call time_integration_renew_mass_state(mesh, mesh%nv_full, nlev               ,&
                                               scalar_hpressure_at_pc_surface_np1%f   ,&
                                               scalar_hpressure_at_pc_face_level_np1%f,& ! overwritten each RK step
                                               dycoreVarCellFull%scalar_delhp_np1%f   ,& ! overwritten each RK step
                                               scalar_hpressure_at_pc_full_level_np1%f,& ! overwritten each RK step
                                               scalar_delhp_at_pc_face_level_np1%f)      ! overwritten each RK step

!-----------------------------------------------------------------------
! get full (moist) mass, if dycore, just hpres (dry mass)
!-----------------------------------------------------------------------

#ifdef NHD_NOOPT
!-----------------------------------------------------------------------
! for regression old LAM; global is not sensitive to this
! note that halo2 is not computed, so POTENTIALLY create trash number for data
! zone in de LAM mode, but not affect model solutions
!-----------------------------------------------------------------------
        call calc_hpe_get_full_mass(nlev, mesh%nv_halo(1), working_mode            , & ! in
#else
        call calc_hpe_get_full_mass(nlev, mesh%nv_full,    working_mode            , & ! in
#endif
                                          scalar_hpressure_at_pc_face_level_np1%f  , & ! in
                                          scalar_hpressure_at_pc_full_level_np1%f  , & ! in
                                          dycoreVarCellFull%scalar_delhp_np1%f     , & ! in
                                          tracerVarCellFull%scalar_mif_n%f         , & ! in
                                          dycoreVarCellFull%scalar_mpressure_n%f   , & ! out
                                          dycoreVarCellFace%scalar_mpressure_n%f)      ! out

    end subroutine dycore_solve_mass_eq_diag

!-----------------------------------------------------------------------
! SEGMENT2
!-----------------------------------------------------------------------

    subroutine dycore_solve_mass_pt_eq

!-----------------------------------------------------------------------
! compute vertical potential temperature mass flux tendency 
! [partial (m*etadot*theta)]
!-----------------------------------------------------------------------
 
        call calc_hpe_tend_vert_mass_flux_2d(mesh%nv_full,mesh%nv_halo(1),nlev,nlev+1,     & ! in
                                             dycoreVarCellFull%scalar_potential_temp_n%f , & ! in
                                             dycoreVarCellFace%scalar_eta_mass_flux_n%f  , & ! in
                                             ver_adv_flag                                , & ! in
                                             dycoreVarCellFull%tend_mass_pt_vert%f) ! out

!--------------------------------------------------------------------------
! compute horizontal potential temperature (pt) flux and pt mass flux div
! it computes pt-flux until ne_compute, but we need ne_halo(1), so do excg
! for LAM, because at bdy we use 2nd-order flux, so we can compute until
! ne_halo(1), but only inside this routine, and will be replaced
!--------------------------------------------------------------------------

        call calc_primal_normal_flux_at_edge(mesh, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f   , &
                                                   dycoreVarEdgeFull%scalar_normal_mass_flux_n%f   , &
                                                   dycoreVarCellFull%scalar_potential_temp_n%f     , &
                                                   dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f, &
                                                   pot_adv_flag(irk_step), rk_substep, nlev)

        dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f = dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f*&
                                                           dycoreVarEdgeFull%scalar_normal_mass_flux_n%f
!
! exchange scheme
!
#ifdef OLD_MASS_PT
#ifndef SEQ_GRIST
        call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
#endif

        ! data init
        scalar_mass_pt_at_pc_full_level_np1%f  = dycoreVarCellFull%scalar_mass_pt_n%f

#ifdef OLD_MASS_PT
        mesh%ne = mesh%ne_halo(1)
        mesh%nv = mesh%nv_halo(1)
#else
        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute
#endif

!inline_opt_here
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do schedule(dynamic,5)
        do iv = 1, mesh%nv
           div_sum(:) = 0._r8
           do inb = 1, mesh%vtx_nnb(iv)
              index_edge  = mesh%vtx_ed(inb,iv)
              do ilev = 1, nlev
                 div_sum(ilev)  = div_sum(ilev)+&
#ifndef DYCORE_OPT
                 dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
#else
                 dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*mesh%edp_leng(index_edge)
#endif
              end do
           end do

!-----------------------------------------------------------------------
! direct (prog) mass*pt update here
!-----------------------------------------------------------------------

#ifndef DYCORE_OPT
        !   dycoreVarCellFull%tend_mass_pt_hori%f(1:nlev,iv) = -div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))
        scalar_mass_pt_at_pc_full_level_np1%f(1:nlev,iv)  = dycoreVarCellFull%scalar_mass_pt_n%f(1:nlev,iv)+&
                             dtime/rk_number*(-div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))+dycoreVarCellFull%tend_mass_pt_vert%f(1:nlev,iv))
#else
        !   dycoreVarCellFull%tend_mass_pt_hori%f(1:nlev,iv) = -div_sum(1:nlev)/(rearth*mesh%plg_areag(iv))

        scalar_mass_pt_at_pc_full_level_np1%f(1:nlev,iv)  = dycoreVarCellFull%scalar_mass_pt_n%f(1:nlev,iv)+&
                             dtime/rk_number*(-div_sum(1:nlev)/(rearth*mesh%plg_areag(iv))+dycoreVarCellFull%tend_mass_pt_vert%f(1:nlev,iv))
#endif
        end do
!$omp end do nowait
!$omp end parallel

!inline_opt_here
        !dycoreVarCellFull%tend_mass_pt_hori%f = -1._r8*dycoreVarCellFull%tend_mass_pt_hori%f

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute horizontal pt flux"

        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute

!-----------------------------------------------------------------------
!old update code
        ! mass-pt opt
        !scalar_mass_pt_at_pc_full_level_np1%f  = dycoreVarCellFull%scalar_mass_pt_n%f+&
        !                                         dtime/rk_number*&
        !                                        (dycoreVarCellFull%tend_mass_pt_hori%f+&
        !                                         dycoreVarCellFull%tend_mass_pt_vert%f)
        ! mass-pt opt
!-----------------------------------------------------------------------

! phys tend add
        if(use_phys.and.ptend_heat_rk_on)then
            scalar_mass_pt_at_pc_full_level_np1%f = scalar_mass_pt_at_pc_full_level_np1%f+&
                                                  dtime/rk_number*ptend_rk%tend_mass_pt_at_pc_full_level%f
        end if
! 2nd lap tend(smg)

        if(ad_dycore_laplacian_2nd)then
            scalar_mass_pt_at_pc_full_level_np1%f = scalar_mass_pt_at_pc_full_level_np1%f+&
                                                    dtime/rk_number*dycoreVarCellFull%tend_mass_pt_laplacian_2nd%f
        end if

#ifndef OLD_MASS_PT

#ifdef LAM_DOMAIN
       call update_assignValueArea_2d(mesh,nlev,"mass_pt",scalar_mass_pt_at_pc_full_level_np1)
#endif
!Exchange data
#ifndef SEQ_GRIST
       call exchange_data_2d_add(mesh,field_head_2d,scalar_mass_pt_at_pc_full_level_np1)
       call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

#endif

!-----------------------------------------------------------------------
! HDC needs dycoreVarCellFull%scalar_potential_temp_n%f for updating geop,
! while NDC needs scalar_mass_pt_at_pc_full_level_np1, so update for all
! halos here for once. The old_mass_pt procedure requires twice update,
! one here only for halo1, and one later for all halos
!-----------------------------------------------------------------------

        dycoreVarCellFull%scalar_potential_temp_n%f = scalar_mass_pt_at_pc_full_level_np1%f/& ! overwritten each RK step
                                                      dycoreVarCellFull%scalar_delhp_np1%f

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute vertical pt flux"

    end subroutine dycore_solve_mass_pt_eq

    subroutine dycore_velocityeq_forward

!----------------------------------------------------------------------
! compute velocity tendency from Coriolis and Kinetic Energy gradient
!----------------------------------------------------------------------

#ifndef SEQ_GRIST
        call clock_begin(clock_nct)
#endif
        !================================================
        ! no longer use tend_nct_once even for testing
        ! IF((tend_nct_once.and.irk_step.eq.1).or.(.not.tend_nct_once)) &  ! control the frequency of NCT evaluation
        !================================================

! donot inline_opt for this module
        call calc_tend_nct_at_edge_full_level(mesh, nlev, irk_step, dtime                         ,&
                                                    scalar_delhp_at_pc_full_level_rk%f            ,&  ! in
                                                    scalar_normal_velocity_at_edge_full_level_rk%f,&  ! in
                                                    dycoreVarEdgeFull%scalar_normal_mass_flux_n%f ,&  ! in
                                                    dycoreVarEdgeFull%tend_normal_velocity_nct%f )   ! ou
#ifndef SEQ_GRIST
        call clock_end(clock_nct)
#endif

        mesh%nv = mesh%nv_halo(1)
        mesh%nt = mesh%nt_halo(1)

!----------------------------------------------------------------------
! KE definition of its gradient
!----------------------------------------------------------------------

#ifndef SEQ_GRIST
        call clock_begin(clock_ket)
#endif
! donot inline_opt for this module 
        call calc_grad_kinetic_energy(mesh,nlev,&
                                           scalar_normal_velocity_at_edge_full_level_rk%f,&
                                           dycoreVarEdgeFull%tend_normal_velocity_ke%f)
#ifndef SEQ_GRIST
        call clock_end(clock_ket)
#endif
        if(iblock .eq. 0.and.write_verbose) print*,"finish NCT and KE"

#ifndef SEQ_GRIST
        call clock_begin(clock_vau)
#endif
!----------------------------------------------------------
! compute vertical advection of normal velocity at edge
!----------------------------------------------------------
!$omp parallel private(ie,icell1,icell2,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_a,scalar_template_1d_nlev_b,scalar_template_1d_nlev_c)
!$omp do schedule(dynamic,5)
        do ie = 1, mesh%ne
!
! center intepolate delhp and mass eta velocity from pc to edge
!
          icell1 = mesh%edt_v(1,ie)
          icell2 = mesh%edt_v(2,ie)

          scalar_template_1d_nlevp_a = &
          0.5_r8*(dycoreVarCellFace%scalar_eta_mass_flux_n%f(1:nlevp,iCell1)+&
                  dycoreVarCellFace%scalar_eta_mass_flux_n%f(1:nlevp,iCell2))
          
          scalar_template_1d_nlev_a  = &
          0.5_r8*(scalar_delhp_at_pc_full_level_rk%f(1:nlev,iCell1)+&
                  scalar_delhp_at_pc_full_level_rk%f(1:nlev,iCell2))

          scalar_template_1d_nlev_b(1:nlev) = scalar_normal_velocity_at_edge_full_level_rk%f(1:nlev,ie)

          call calc_hpe_vert_advection(nlev,nlevp,&
                                       scalar_template_1d_nlevp_a, &
                                       scalar_template_1d_nlev_a , &
                                       scalar_template_1d_nlev_b , &
                                       scalar_template_1d_nlev_c )

          dycoreVarEdgeFull%tend_normal_velocity_vadv%f(1:nlev,ie) = scalar_template_1d_nlev_c

        end do
!$omp end do nowait
!$omp end parallel
        if(iblock .eq. 0.and.write_verbose) print*,"finish vertical advection of normal wind"
#ifndef SEQ_GRIST
        call clock_end(clock_vau)
#endif

    end subroutine dycore_velocityeq_forward

!-----------------------------------------------------------------------
! SEGMENT3
!-----------------------------------------------------------------------

    subroutine dycore_diagnose_omega

     real(r8)                            :: scalar_mpressure_at_edge_full_level_n(nlev,mesh%ne_full)
     real(r8)                            :: tend_pidpiv_at_pc_full_level_hori(nlev,mesh%nv_full)

!===========================================================================================
! Diagnose dpim/dt (omegaMoist)
! Besides being a diagnostics:
! for HDC this diag is necessary for WRF physics; for NDC this diag produces similar results
! as we use www, g and rhom for diag, but not itentical even in logic;
! 2021-01: to use this for ndc's cu scheme make tropical rainfall looks better real-world 
!===========================================================================================

       if(irk_step.eq.nrk)then ! diagnose d(mpressure)/dt at the final rk stage of one dycore_step
           mesh%ne = mesh%ne_compute
! pi flux
! inline_opt_here
           !call calc_primal_normal_flux_at_edge(mesh, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, & ! dum
           !                                           dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, & ! dum
           !                                           scalar_mpressure_at_pc_full_level_old        , & ! halo(1)->ne
           !                                           scalar_mpressure_at_edge_full_level_n        , &
           !                                           2, rk_substep, nlev)
!$omp parallel  private(ie,ilev,icell1,icell2)
!$omp do schedule(static,100) 
          do ie = 1, mesh%ne
             iCell1      = mesh%edt_v(1, ie)
             iCell2      = mesh%edt_v(2, ie)
             do ilev = 1, nlev
                scalar_mpressure_at_edge_full_level_n(ilev, ie) = 0.5*(scalar_mpressure_at_pc_full_level_old(ilev,iCell1)+&
                                                                       scalar_mpressure_at_pc_full_level_old(ilev,iCell2))
             end do
          end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here

! pi*dpi*V flux (old-time)
           scalar_mpressure_at_edge_full_level_n = scalar_mpressure_at_edge_full_level_n*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f

! inline_opt_here
           !call divergence_operator_2d(mesh,scalar_mpressure_at_edge_full_level_n, &
           !                                 tend_pidpiv_at_pc_full_level_hori    , &
           !                                 nlev)
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!$omp do schedule(dynamic,5)
           do iv = 1, mesh%nv
              div_sum(:) = zero
              do inb = 1, mesh%vtx_nnb(iv)
                 index_edge  = mesh%vtx_ed(inb, iv)
                 do ilev = 1, nlev
#ifndef DYCORE_OPT
                    div_sum(ilev)  = div_sum(ilev)+scalar_mpressure_at_edge_full_level_n(ilev,index_edge)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(index_edge)
#else
                    div_sum(ilev)  = div_sum(ilev)+scalar_mpressure_at_edge_full_level_n(ilev,index_edge)*mesh%plg_nr(inb,iv)*mesh%edp_leng(index_edge)
#endif
                 end do
              end do
#ifndef DYCORE_OPT
              tend_pidpiv_at_pc_full_level_hori(1:nlev,iv) = div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))
#else
              tend_pidpiv_at_pc_full_level_hori(1:nlev,iv) = div_sum(1:nlev)/(rearth*mesh%plg_areag(iv))
#endif
           end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here

           mesh%nv = mesh%nv_compute
           do iv = 1, mesh%nv
              do ilev = 1, nlev
                 dycoreVarCellFull%scalar_omega_n%f(ilev,iv)=(dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)-scalar_mpressure_at_pc_full_level_old(ilev,iv))/dtime+&
                                                              dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)+&
                                                            ((tend_pidpiv_at_pc_full_level_hori(ilev,iv)-dycoreVarCellFull%tend_mass_hori%f(ilev,iv)*&
                                                              scalar_mpressure_at_pc_full_level_old(ilev,iv))/scalar_delhp_at_pc_full_level_rk%f(ilev,iv))
                 if(eta_full_b(ilev).eq.0)then
                    dycoreVarCellFull%scalar_omega_n%f(ilev,iv) = dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)
                 end if
              end do
           end do
        mesh%nv = mesh%nv_halo(1)
       end if

    end subroutine dycore_diagnose_omega

!-----------------------------------------------------------------------
! SEGMENT4
!-----------------------------------------------------------------------

    subroutine hydrostatic_adjust

!-----------------------------------------------------------------------
! The modifiction of this part for mhdc may affect bit reproduce of dry hdc,
! although physically identical
!-----------------------------------------------------------------------
#ifdef DYCORE_NHD
#ifndef REG_OLD_NHDAMP
        if(nh_dynamics)then
! hold old value
           scalar_pressure_at_pc_face_level_old    = dycoreVarCellFace%scalar_pressure_n%f
           scalar_pressure_at_pc_full_level_old    = dycoreVarCellFull%scalar_pressure_np1%f
        end if
#endif
#endif
        dycoreVarCellFull%scalar_pressure_np1%f = dycoreVarCellFull%scalar_mpressure_n%f
        dycoreVarCellFace%scalar_pressure_n%f   = dycoreVarCellFace%scalar_mpressure_n%f

!-----------------------------------------------------------------------
! Diagnose geopotential and evaluate its gradient
!-----------------------------------------------------------------------

!$omp parallel  private(iv,scalar_template_1d_nlev_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_b,scalar_template_a,scalar_template_1d_nlevp_b,scalar_template_1d_nlev_c)
!$omp do schedule(dynamic,5)

!-----------------------------------------------------------------------
! necessary, as halo1 phi is used for updating ne_compute hpgf
!-----------------------------------------------------------------------

        ncell_do = mesh%nv_full
#ifdef NHD_NOOPT
! for regression LAM; global is not sensitive to this
        ncell_do = mesh%nv_halo(1)
#endif

        do iv = 1, ncell_do

!-----------------------------------------------------------------------
! recover temperature, and store
!-----------------------------------------------------------------------

           !scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_np1%f(:,iv)*&
           scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)*&
                                      ((dycoreVarCellFull%scalar_pressure_np1%f(:,iv)/p00)**(rdry/cp))/&
                                       (one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv))

           dycoreVarCellFull%scalar_temp_n%f(:,iv) = scalar_template_1d_nlev_a ! output, overwritten each RK step
           scalar_template_1d_nlev_a  = scalar_template_1d_nlev_a*(one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)))) ! Tv

!-----------------------------------------------------------------------
! new pressure face
!-----------------------------------------------------------------------

           !scalar_template_1d_nlevp_a = scalar_hpressure_at_pc_face_level_n%f(:,iv)
           scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f(:,iv)

!-----------------------------------------------------------------------
! new delp full, not in effect actually
!-----------------------------------------------------------------------
           !scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delhp_np1%f(:,iv)
           scalar_template_1d_nlev_b  = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,iv)-dycoreVarCellFace%scalar_pressure_n%f(1:nlev,iv)
!-----------------------------------------------------------------------
! geopotential at surface
!-----------------------------------------------------------------------
           scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)
!-----------------------------------------------------------------------
! diagnose geop at face based on pressure, which is related to
! the definition of model layer
!-----------------------------------------------------------------------
           call calc_hpe_hydro(nlev,nlevp,&
                               scalar_template_1d_nlev_a  ,& ! T or Tv
                               scalar_template_1d_nlevp_a ,& ! hp or p
                               scalar_template_1d_nlev_b  ,& ! delhp or delp
                               scalar_template_a          ,& ! phis
                               scalar_template_1d_nlevp_b ,& ! phi-face
                               scalar_template_1d_nlev_c )   ! phi-full

           dycoreVarCellFace%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlevp_b
           dycoreVarCellFull%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlev_c  ! output
#ifdef REG_OLD_NHDAMP
           scalar_phi_at_pc_face_level_np1%f(:,iv)         = dycoreVarCellFace%scalar_geopotential_n%f(:,iv)
#endif
        end do

!$omp end do nowait
!$omp end parallel

#ifdef DYCORE_NHD
#ifndef REG_OLD_NHDAMP
!-----------------------------------------------------------------------
! REG_OLD_NHDAMP is directly restore a hydro state as damping, below one
! is more subtle; this routine is called under nh_dynamics only for damping
!-----------------------------------------------------------------------
        if(nh_dynamics)then
           scalar_phi_at_pc_face_level_np1%f = scalar_phi_at_pc_face_level_np1%f+&
                                              damp_kesi*(scalar_phi_at_pc_face_level_np1%f-dycoreVarCellFace%scalar_geopotential_n%f)

           dycoreVarCellFace%scalar_geopotential_n%f = scalar_phi_at_pc_face_level_np1%f ! used for hpgf
           dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,:) = 0.5_r8*(dycoreVarCellFace%scalar_geopotential_n%f(1:nlev ,:)+&
                                                                         dycoreVarCellFace%scalar_geopotential_n%f(2:nlevp,:))

           dycoreVarCellFull%scalar_pressure_np1%f = scalar_pressure_at_pc_full_level_old+&
                                          damp_kesi*(scalar_pressure_at_pc_full_level_old-dycoreVarCellFull%scalar_mpressure_n%f)

           dycoreVarCellFace%scalar_pressure_n%f   = scalar_pressure_at_pc_face_level_old+&
                                          damp_kesi*(scalar_pressure_at_pc_face_level_old-dycoreVarCellFace%scalar_mpressure_n%f)
        end if
#endif
#endif

          ! since in hdc, delp is only used for pgf, and we use (6), so delp should be delhp, to cancel the last term in pgf6, not so !
          ! dycoreVarCellFull%scalar_delp_np1%f   = dycoreVarCellFull%scalar_delhp_np1%f ! only valid for dry hdc
          ! old dry-hydrostatic model only uses this
          ! dycoreVarCellFull%scalar_alpha_np1%f  = rdry*dycoreVarCellFull%scalar_temp_n%f/scalar_hpressure_at_pc_full_level_np1%f

!-----------------------------------------------------------------------
! real layer-averaged
!-----------------------------------------------------------------------
! affect hpgf
          dycoreVarCellFull%scalar_delp_np1%f   = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,:)-dycoreVarCellFace%scalar_pressure_n%f(1:nlev,:) 
          dycoreVarCellFull%scalar_alpha_np1%f  =(dycoreVarCellFace%scalar_geopotential_n%f(1:nlev,:)-dycoreVarCellFace%scalar_geopotential_n%f(2:nlevp,:))/&
                                                  dycoreVarCellFull%scalar_delhp_np1%f

          return
    end subroutine hydrostatic_adjust

!-----------------------------------------------------------------------
! SEGMENT5
!-----------------------------------------------------------------------

    subroutine dycore_velocityeq_hpgf

!----------------------------------------------------------------------
! use below subroutines with smallest modification to compute PGF for
! normal velocity. The below should not change in the case of ndc or hdc
!----------------------------------------------------------------------

      if(iblock .eq. 0.and.write_verbose) print*,"finish geop, grad_delp, grad_pressure gradient"

!-----------------------------------------------------------------------
! compute real PGF tendency
!-----------------------------------------------------------------------

      select case(hor_pgf_flag)

      case(6)  ! ptb, default with tp1 profile

          call grist_dycore_ref_atmos_create_ptb
          !dycoreVarEdgeFull%tend_normal_velocity_gz%f  = zero

!$omp parallel  private(ie,icell1,icell2,cr_num,tmp1,tmp2,tmp3,tmp4)
!$omp do schedule(dynamic,5)
          do ie = 1, mesh%ne ! ne_compute

             icell1 = mesh%edt_v(1,ie)
             icell2 = mesh%edt_v(2,ie)
             cr_num = rearth*mesh%edt_leng(ie)

             tmp1   = scalar_grad_hpres_at_edge_full_level_bar(1:nlev,ie)*&
                      0.5_r8*(scalar_alphad_at_pc_full_level_ptb(1:nlev,icell1)+scalar_alphad_at_pc_full_level_ptb(1:nlev,icell2))

             tmp2   = 0.5_r8*(dycoreVarCellFull%scalar_alpha_np1%f(1:nlev,icell1)+dycoreVarCellFull%scalar_alpha_np1%f(1:nlev,icell2))*&
                             (scalar_pressure_at_pc_full_level_ptb(1:nlev,icell2)-scalar_pressure_at_pc_full_level_ptb(1:nlev,icell1))/cr_num

             tmp3   = (scalar_geop_at_pc_full_level_ptb(1:nlev,icell2)-scalar_geop_at_pc_full_level_ptb(1:nlev,icell1))/cr_num

             tmp4   = 0.5_r8*(scalar_delp_at_pc_full_level_ptb(1:nlev,icell1)/dycoreVarCellFull%scalar_delhp_np1%f(1:nlev,icell1) +&
                              scalar_delp_at_pc_full_level_ptb(1:nlev,icell2)/dycoreVarCellFull%scalar_delhp_np1%f(1:nlev,icell2))*&
                             (dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,icell2)-dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,icell1))/cr_num

             dycoreVarEdgeFull%tend_normal_velocity_pgf%f(1:nlev,ie) = -1._r8*(tmp1+tmp2+tmp3+tmp4)*tracerVarEdgeFull%scalar_mif_n%f(1:nlev,ie)
          end do
!$omp end do nowait
!$omp end parallel 
      case default
          if(iblock .eq. 0) print*, "you must select a PGF term in hor_pgf_flag"
          stop
      end select

      if(iblock .eq. 0.and.write_verbose) print*,"finish HPGF"

    end subroutine dycore_velocityeq_hpgf

   end subroutine dycore_time_integration_run

   subroutine dycore_time_integration_init(mesh,nlev)

#ifdef DYCORE_NHD
! dycore-nh
     use grist_nh_driver_module,           only: grist_nh_dynamics_init
     use grist_nh_explicit_tend_module_2d, only: grist_nh_et_init
#endif
     use grist_dycore_ref_atmos,           only: grist_dycore_ref_atmos_init
! io
     type(global_domain),  intent(inout) :: mesh
     integer(i4)        ,  intent(in)    :: nlev
! local
     integer(i4)      :: iv
     real(r8)         :: scalar_template_1d_nlev_a(nlev)
     real(r8)         :: scalar_template_1d_nlev_b(nlev)
     real(r8)         :: scalar_template_1d_nlev_c(nlev)
     real(r8)         :: scalar_template_1d_nlevp_a(nlev+1)
     real(r8)         :: scalar_template_1d_nlevp_b(nlev+1)
     real(r8)         :: scalar_template_a

#ifdef DYCORE_NHD
     call grist_nh_dynamics_init(mesh)
     call grist_nh_et_init(mesh,nlev)
#endif
!
! initially compute hpressure (full&face), delhp (full&face), based on surface hpressure
!
        call time_integration_renew_mass_state(mesh, mesh%nv_full, nlev, &
                                                     dycoreVarSurface%scalar_hpressure_n%f  ,&
                                                     dycoreVarCellFace%scalar_hpressure_n%f ,&
                                                     dycoreVarCellFull%scalar_delhp_n%f     ,&
                                                     dycoreVarCellFull%scalar_hpressure_n%f ,&
                                                     dycoreVarCellFace%scalar_delhp_n%f    )

        if(trim(working_mode).eq.'dycore')then
           dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_hpressure_n ! init as hpressure
           dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_hpressure_n ! init as hpressure
        end if
! above for dycore mode

        if(trim(working_mode).ne.'dycore')then
           ! already defined in dtp
           !dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_pressure_n
           !dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_pressure_n
        end if

        do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)  =  &
                                                          dycoreVarCellFull%scalar_temp_n%f(:,iv)/&
                                                  ((dycoreVarCellFull%scalar_pressure_n%f(:,iv)/p00)**(rdry/cp))*&
                                                  (one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv))

            dycoreVarCellFull%scalar_mass_pt_n%f(:,iv)  =  dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)*&
                                                          dycoreVarCellFull%scalar_delhp_n%f(:,iv)
        end do

        IF(trim(gcm_testcase).ne.'SCHAR'     .and.&
           trim(gcm_testcase).ne.'DCMIP2-1'  .and.&
           trim(gcm_testcase).ne.'real-ERAIM'.and.&
           trim(gcm_testcase).ne.'real-ERAIP'.and.&
           trim(gcm_testcase).ne.'real-GFS'  .and.&
           trim(gcm_testcase).ne.'real-WRFDA'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-BW'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-TC'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC-A'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC-B'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC1' .and.&
           trim(gcm_testcase).ne.'DCMIP2016-SCXX' .and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC')then
        
        do iv = 1, mesh%nv

           scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_temp_n%f(:,iv)      ! temperature
           scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_hpressure_n%f(:,iv) ! hpressure face
           scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delhp_n%f(:,iv)     ! delhp full
           scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv) ! geopotential at surface

           call calc_hpe_hydro(nlev,nlevp,scalar_template_1d_nlev_a  ,&
                               scalar_template_1d_nlevp_a ,&
                               scalar_template_1d_nlev_b  ,&
                               scalar_template_a          ,&
                               scalar_template_1d_nlevp_b ,&
                               scalar_template_1d_nlev_c )

           dycoreVarCellFace%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlevp_b
           dycoreVarCellFull%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlev_c ! output
#ifdef DYCORE_NHD
           if(nh_dynamics) dycoreVarCellFace%scalar_phi_n%f(:,iv)          = dycoreVarCellFace%scalar_geopotential_n%f(:,iv)
#endif
        end do

#ifdef DYCORE_NHD
        if(nh_dynamics)then
           dycoreVarCellFace%scalar_www_n%f            = zero 
           dycoreVarCellFace%scalar_phi_n%f(nlev+1,:)  = dycoreVarSurface%scalar_geopotential_n%f ! surface is given in ini condition
        end if
#endif
        END IF

        mesh%ne = mesh%ne_compute
        call grist_dycore_ref_atmos_init(mesh)
        mesh%ne = mesh%ne_full

!-----------------------------------------------------------------------
! init mpressure
!-----------------------------------------------------------------------

        call calc_hpe_get_full_mass(nlev, mesh%nv_halo(1), working_mode          , & ! in
                                          dycoreVarCellFace%scalar_hpressure_n%f , & ! in
                                          dycoreVarCellFull%scalar_hpressure_n%f , & ! in
                                          dycoreVarCellFull%scalar_delhp_n%f     , & ! in
                                          tracerVarCellFull%scalar_mif_n%f       , & ! in
                                          dycoreVarCellFull%scalar_mpressure_n%f , & ! out
                                          dycoreVarCellFace%scalar_mpressure_n%f)    ! out
! init type
      scalar_delhp_at_pc_face_level_np1     = dycoreVarCellFace%scalar_geopotential_n
      scalar_hpressure_at_pc_full_level_np1 = dycoreVarCellFull%scalar_hpressure_n
      scalar_hpressure_at_pc_face_level_np1 = dycoreVarCellFace%scalar_geopotential_n
      scalar_phi_at_pc_face_level_np1       = dycoreVarCellFace%scalar_geopotential_n

#ifdef DYCORE_NHD
      if(nh_dynamics) scalar_www_at_pc_face_level_np1  = dycoreVarCellFace%scalar_www_n 
#endif

#ifndef SEQ_GRIST
      clock_mas      = clock_id('mas')
      clock_nct      = clock_id('nct') 
      clock_nhd      = clock_id('nhd') 
      clock_ptm      = clock_id('ptm')
      clock_pgf      = clock_id('pgf')
      clock_vau      = clock_id('vau')
      clock_ket      = clock_id('ket')
      clock_diagom   = clock_id('diagom')
      clock_mainexch = clock_id('mainexch')
      clock_rkl      = clock_id('rkl')
      clock_fnlrk    = clock_id('fnlrk')
#endif


      return
   end subroutine dycore_time_integration_init

   subroutine dycore_time_integration_final

#ifdef DYCORE_NHD
     use grist_nh_driver_module,             only: grist_nh_dynamics_final
     use grist_nh_explicit_tend_module_2d,   only: grist_nh_et_final
#endif
     use grist_dycore_ref_atmos,             only: grist_dycore_ref_atmos_final

#ifdef DYCORE_NHD
      call grist_nh_dynamics_final
      call grist_nh_et_final
#endif
      call grist_dycore_ref_atmos_final

   end subroutine dycore_time_integration_final

!-----------------------------------------------------------------------
! PRIVATE BELOW
!-----------------------------------------------------------------------

   subroutine time_integration_renew_mass_state(mesh,ncell_do, nlev, &
                                                scalar_hpressure_at_pc_surface    ,&
                                                scalar_hpressure_at_pc_face_level ,&
                                                scalar_delhp_at_pc_full_level     ,&
                                                scalar_hpressure_at_pc_full_level ,&
                                                scalar_delhp_at_pc_face_level     )
! io
     use omp_lib
     type(global_domain)  , intent(in)    :: mesh
     integer(i4)          , intent(in)    :: ncell_do, nlev
     real(r8)             , intent(in)    :: scalar_hpressure_at_pc_surface(:)
     real(r8)             , intent(inout) :: scalar_hpressure_at_pc_face_level(:,:)
     real(r8)             , intent(inout) :: scalar_delhp_at_pc_full_level(:,:)
     real(r8)             , intent(inout) :: scalar_hpressure_at_pc_full_level(:,:)
     real(r8)             , intent(inout) :: scalar_delhp_at_pc_face_level(:,:)
! local
     integer(i4)                          :: iv, ilev

!-----------------------------------------------------------------------
! renew mass state, compute hpressure at face level, full level, 
! delhp, based on input surface hpressure
!-----------------------------------------------------------------------

!$omp parallel  private(iv,ilev)
!$omp do schedule(dynamic,5)
        do iv = 1, ncell_do
!-----------------------------------------------------------------------
! inline_opt
!-----------------------------------------------------------------------
           scalar_hpressure_at_pc_face_level(1:nlevp,iv) = eta_face_a(1:nlevp)*p00+eta_face_b(1:nlevp)*scalar_hpressure_at_pc_surface(iv)
           scalar_delhp_at_pc_full_level(1:nlev,iv)    = scalar_hpressure_at_pc_face_level(2:nlevp,iv)-scalar_hpressure_at_pc_face_level(1:nlev,iv)
           scalar_hpressure_at_pc_full_level(1:nlev,iv)  = eta_full_a(1:nlev)*p00 +eta_full_b(1:nlev) *scalar_hpressure_at_pc_surface(iv)

           do ilev = 2, nlev
              scalar_delhp_at_pc_face_level(ilev,iv)= scalar_hpressure_at_pc_full_level(ilev,iv)-scalar_hpressure_at_pc_full_level(ilev-1,iv)
           end do

           scalar_delhp_at_pc_face_level(1,iv)      = 0.5_r8*scalar_delhp_at_pc_full_level(1,iv)
           scalar_delhp_at_pc_face_level(nlev+1,iv) = 0.5_r8*scalar_delhp_at_pc_full_level(nlev,iv)
!-----------------------------------------------------------------------
! inline_opt
!-----------------------------------------------------------------------
        end do
!$omp end do nowait
!$omp end parallel

        return
   end subroutine time_integration_renew_mass_state

  end module grist_dycore_time_integration_2d
