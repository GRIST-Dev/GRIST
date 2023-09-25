
!----------------------------------------------------------------------------
! Created on Aug 18 2023
! Author: Siyuan Chen
! Version 1.0
! Description: Building upon the grist_dycore_time_integration_2d module, 
!              replace selected variables with single precision while maintaining similarity 
!              to the results of double precision variables.
!----------------------------------------------------------------------------
module grist_dycore_time_integration_2d_mixed

    use grist_constants,       only: rvap,rdry,ptfactor,cp,p00,zero, i4, r8, gravity, one, rearth
    use grist_constants,       only: r4 => ns
    use grist_nml_module,      only: nlev, nlevp, mas_adv_flag, pot_adv_flag, ver_adv_flag, nrk, hor_pgf_flag, &
                                     nh_dynamics, gcm_testcase,  nsteps, &
                                     ad_dycore_laplacian_2nd, ad_dycore_laplacian_4th, ad_dycore_laplacian_6th, &
                                     tend_nct_once, working_mode, use_phys, physpkg, &
                                     ptendSubDiabPhys, ptend_wind_rk_on, ptend_heat_rk_on, ptend_dycore_f1_on, ptend_dycore_heat_f1_on, &
                                     ptend_dycore_wind_f1_on, write_stepinfo, write_verbose, adjphi, &
                                     restore_hydro, restore_hydro_minsteps,restore_hydro_intsteps
#ifdef AMIPW_PHYSICS
    use grist_physw_nml_module, only: step_cu
#endif
    use grist_domain_types,    only: global_domain
    use grist_data_types,      only: scalar_1d_field, scalar_2d_field, exchange_field_list_2d, exchange_field_list_1d
!
! horizontal module
!
    use grist_dycore_primal_flux_operators_2d_mixed, only: calc_primal_normal_flux_at_edge
    use grist_dycore_hori_swe_module_2d_mixed,       only: calc_tend_nct_at_edge_full_level, calc_grad_kinetic_energy, grist_dycore_hori_swe_clock_init
    use grist_dycore_gcd_recon_module_2d_mixed,      only: divergence_operator_2d
    use grist_dycore_ref_atmos_mixed
!
! vertical qh
!
    use grist_hpe_hydro_pgf_mixed,          only: calc_hpe_hydro               , &
                                            	  calc_hpe_hpressure_face_level , &
                                            	  calc_hpe_hpressure_full_level , &
                                            	  calc_hpe_delhp, &
                                            	  calc_hpe_get_full_mass, &
                                                  calc_hpe_get_full_mass1

    use grist_hpe_vertical_advection_mixed, only: calc_hpe_vert_advection, &
                                            	  calc_hpe_tend_vert_mass_flux_2d
    use grist_hpe_continuity_mixed        , only: calc_hpe_tend_continuity_2d
! nh_dynamics
    use grist_nh_driver_module_mixed,       only: grist_nh_dynamics_run, ndc_restore_flag
! data
    use grist_dycore_module_vars_control_mixed, only: deta_full, deta_face,eta_full,eta_full_b,edp_leng
! dycore data
    use grist_dycore_vars_module,     only: dycoreVarCellFace, dycoreVarCellFull, dycoreVarEdgeFull, dycoreVarEdgeFace, dycoreVarSurface
! tracer
    use grist_tracer_transport_vars_module, only: tracerVarCellFace, tracerVarEdgeFull, tracerVarCellFull
!
! physics
! 
    use grist_physics_data_structure, only: ptend_f1, ptend_rk, ptend_f2
! mpi
    use grist_mpi
#ifndef SEQ_GRIST
    use grist_config_partition,       only: debug_data_1d,&
                                            debug_data_2d,&
                                            exchange_data_2d_add , &
                                            exchange_data_2d,&
                                            exchange_data_1d_add,&
                                            exchange_data_1d,&
                                            exchange_data_1d_r4,&
                                            exchange_data_2d_r4
    use grist_clocks,                 only: clock_id, clock_begin, clock_end
#endif

    implicit none

    private

    public :: dycore_time_integration_init,&
              dycore_time_integration_run ,&
              dycore_time_integration_final
!
! local varitable for template io of subroutine
!
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_a
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_b
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_c
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_d
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlevp_a
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlevp_b
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlevp_c
     real(r8)                            :: scalar_template_a
     type(scalar_2d_field)               :: scalar_template_2d_ne_b
     type(scalar_2d_field)               :: scalar_pressure_at_pc_face_level_rk  ! nh-diag-var
     type(scalar_2d_field)               :: scalar_www_at_pc_face_level_np1      ! nh-prog-var
     type(scalar_2d_field)               :: scalar_phi_at_pc_face_level_np1      ! nh-prog-var
     type(scalar_2d_field)               :: scalar_hpressure_at_pc_face_level_np1
     type(scalar_2d_field)               :: scalar_delhp_at_pc_face_level_np1
     type(scalar_2d_field)               :: scalar_delhp_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_hpressure_at_pc_full_level_np1
     integer(i4) :: clock_rkl, clock_mas, clock_nct, clock_vau, clock_ket, clock_diagom, clock_nhd, clock_ptm, clock_pgf, clock_mainexch, clock_fnlrk, clock_update, clock_prep
     integer(i4) :: clock_fnlrk1, clock_fnlrk2, clock_fnlrk3

#ifdef NEW_KESI
     real(r8), parameter :: damp_kesi = -0.5 ! this is a regression
#else
     real(r8), parameter :: damp_kesi = -1 ! this is a regression
#endif      

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

!     type(scalar_1d_field)               :: scalar_hpressure_at_pc_surface_rk
!     type(scalar_2d_field)               :: scalar_mass_pt_at_pc_full_level_rk
     type(scalar_2d_field)               :: scalar_www_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_phi_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_normal_velocity_at_edge_full_level_rk
     type(scalar_2d_field)               :: scalar_delhp_at_pc_full_level_rk
! local np1 vars
     type(scalar_1d_field)               :: scalar_hpressure_at_pc_surface_np1
     type(scalar_2d_field)               :: scalar_normal_velocity_at_edge_full_level_np1
     type(scalar_2d_field)               :: scalar_mass_pt_at_pc_full_level_np1
! local
     real(r4)                            :: dtime_r4
     real(r8)                            :: rk_substep, fb_step
     real(r8)                            :: cr_num, cr_act, cr_sum, eta, kv 
     real(r4)                            :: div_sum_r4(nlev)
     real(r8)                            :: div_sum(nlev), tmp1(nlev), tmp2(nlev), tmp3(nlev), tmp4(nlev)

     real(r8)                            :: tend_pidpiv_at_pc_full_level_hori(nlev,mesh%nv_full)
     real(r4)                            :: scalar_mpressure_at_edge_full_level_n(nlev,mesh%ne_full)

     real(r4)                            :: scalar_mpressure_at_pc_full_level_old(nlev,mesh%nv_full)
     real(r4)                            :: scalar_pressure_at_pc_face_level_old(nlevp,mesh%nv_full)
     real(r8)                            :: scalar_pressure_at_pc_full_level_old(nlev,mesh%nv_full)

     real(r8)                            :: scalar_pt_at_pc_full_level_old(nlev,mesh%nv_full)
     integer(i4)                         :: nsplit                 ! how many small steps within DT
     integer(i4)                         :: rk_number              ! number to divide DT in RK step
     integer(i4)                         :: irk_step
     integer(i4)                         :: ifb_step
     integer(i4)                         :: it, ie, iv, ilev, inb
     integer(i4)                         :: icell1, icell2, index_edge
     integer(i4)                         :: iblock
! local
     real(r8) :: scalar_template_1d_nlevp_d(nlevp)
     real(r8) :: scalar_template_1d_nlevp_e(nlevp)
     real(r8) :: scalar_template_1d_nlevp_f(nlevp)

      type(exchange_field_list_2d),pointer :: field_head_2d
      type(exchange_field_list_1d),pointer :: field_head_1d
      
      dtime_r4 = dtime

      field_head_2d=>null()
      field_head_1d=>null()
      iblock = mpi_rank()
      !iblock = -1
!================================================
! added for cumulus wrfphys
!================================================

#ifdef AMIPW_PHYSICS
      IF(present(idstep))then
! only at the first dycore and tracer step
         if(idstep.eq.1.and.itstep.eq.1.and.mod(itimestep-1,step_cu).eq.0) dycoreVarCellFull%tend_pt_n%f = zero
      END IF
#endif
!================================================
! added for cumulus wrfphys
!================================================

!
#ifndef SEQ_GRIST
     call clock_begin(clock_prep)
#endif
!
! prepare for RK dyn integration
!
!      scalar_hpressure_at_pc_surface_rk            = dycoreVarSurface%scalar_hpressure_n
!      scalar_mass_pt_at_pc_full_level_rk           = dycoreVarCellFull%scalar_mass_pt_n
      if(nh_dynamics)then
         scalar_www_at_pc_face_level_rk%f_r4      = dycoreVarCellFace%scalar_www_n%f
         scalar_www_at_pc_face_level_rk%pos       = dycoreVarCellFace%scalar_www_n%pos
         scalar_phi_at_pc_face_level_rk%f_r4      = dycoreVarCellFace%scalar_phi_n%f
         scalar_phi_at_pc_face_level_rk%pos       = dycoreVarCellFace%scalar_phi_n%pos
         dycoreVarCellFull%scalar_pressure_rk= dycoreVarCellFull%scalar_pressure_n
      end if
      scalar_normal_velocity_at_edge_full_level_rk%f_r4 = dycoreVarEdgeFull%scalar_normal_velocity_n%f
      scalar_normal_velocity_at_edge_full_level_rk%pos  = dycoreVarEdgeFull%scalar_normal_velocity_n%pos
      scalar_delhp_at_pc_full_level_rk%f_r4             = dycoreVarCellFull%scalar_delhp_n%f
      scalar_delhp_at_pc_full_level_rk%pos              = dycoreVarCellFull%scalar_delhp_n%pos
      scalar_delhp_at_pc_face_level_rk             = dycoreVarCellFace%scalar_delhp_n
      scalar_mpressure_at_pc_full_level_old        = dycoreVarCellFull%scalar_mpressure_n%f

! raw pt
      tracerVarCellFull%scalar_tracer_mxrt_n%f = tracerVarCellFull%scalar_tracer_mxrt_n%f_r4  ! copy

      scalar_pt_at_pc_full_level_old = dycoreVarCellFull%scalar_potential_temp_n%f/(one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:))
#ifndef SEQ_GRIST
     call clock_end(clock_prep)
#endif
!
! BEGIN THIS RK STEP, use zzz to represent code that should be sleep foreever
!
! we have commented if statement assocaited stencil_exchange_flag,
! stencil_exchange_flag=.true. is default now
!
#ifndef SEQ_GRIST
     call clock_begin(clock_rkl)
#endif

     DO irk_step = 1, nrk

        rk_number  = nrk+1-irk_step
        rk_substep = dtime/rk_number

!-----------------------------------------------------------------------
! compute normal mass flux {DELP*V} and mass flux div {GRAD.(DELP*V)} 
! at each full level
!-----------------------------------------------------------------------
#ifndef SEQ_GRIST
        call clock_begin(clock_mas)
#endif

        mesh%ne = mesh%ne_halo(1)
! inline_opt_here
!        call calc_primal_normal_flux_at_edge(mesh, scalar_normal_velocity_at_edge_full_level_rk%f, &
!                                                   scalar_normal_velocity_at_edge_full_level_rk%f, & ! dum
!                                                   scalar_delhp_at_pc_full_level_rk%f            , &
!                                                   dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, &
!                                                   mas_adv_flag, rk_substep, nlev)

!$omp parallel  private(ie,ilev)
!$omp do schedule(static,100) 
        do ie = 1, mesh%ne_halo(1)
           do ilev = 1, nlev
              dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4(ilev, ie) = 0.5_r4*(scalar_delhp_at_pc_full_level_rk%f_r4(ilev,mesh%edt_v(1,ie))+&
                                                                             scalar_delhp_at_pc_full_level_rk%f_r4(ilev,mesh%edt_v(2,ie)))
           end do
        end do
!$omp end do nowait
!$omp end parallel 
! inline_opt_here
        dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4 = dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4*&
                                                        scalar_normal_velocity_at_edge_full_level_rk%f_r4
!================================================
! exchange scheme, needed if o3 is used, but o2
! is used so hard coded here!
!================================================

!        call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_mass_flux_n)
!        call exchange_data_2d(mesh%local_block,field_head_2d)
        !mesh%ne = mesh%ne_halo(1)

        ! from ne_halo(1) to nv_halo(1)
        mesh%nv = mesh%nv_halo(1)
! inline_opt_here
!        call divergence_operator_2d(mesh, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, &
!                                          dycoreVarCellFull%tend_mass_hori%f, nlev)
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do schedule(dynamic,5)
        do iv = 1, mesh%nv_halo(1)
           div_sum(1:nlev) = zero
           do inb = 1, mesh%vtx_nnb(iv)
              index_edge  = mesh%vtx_ed(inb,iv)
              do ilev = 1, nlev
                 div_sum(ilev)  =  div_sum(ilev)+dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
              end do
           end do
#ifdef DBL_MASS
           dycoreVarCellFull%tend_mass_hori%f(1:nlev,iv) = div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))
#else
           dycoreVarCellFull%tend_mass_hori%f_r4(1:nlev,iv) = div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))
#endif
        end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute normal mass flux",mas_adv_flag

        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute
!
! compute PS tendency and face level vertical mass flux (m*etadot)
!
        call calc_hpe_tend_continuity_2d(mesh%nv_halo(1), nlev                     , &
#ifdef DBL_MASS
                                         dycoreVarCellFull%tend_mass_hori%f         , & ! in
                                         dycoreVarSurface%tend_hpressure_cnty%f       , & ! out
                                         dycoreVarCellFace%scalar_eta_mass_flux_n%f)    ! out
#else
                                         dycoreVarCellFull%tend_mass_hori%f_r4         , & ! in
                                         dycoreVarSurface%tend_hpressure_cnty%f_r4       , & ! out
                                         dycoreVarCellFace%scalar_eta_mass_flux_n%f_r4)    ! out
#endif

#ifdef DBL_MASS
        dycoreVarCellFull%scalar_eta_mass_flux_n%f(1:nlev,:) = 0.5_r8*dycoreVarCellFace%scalar_eta_mass_flux_n%f(1:nlev,:)+&
                                                               0.5_r8*dycoreVarCellFace%scalar_eta_mass_flux_n%f(2:nlev+1,:)
#else
        dycoreVarCellFull%scalar_eta_mass_flux_n%f_r4(1:nlev,:) = 0.5_r4*dycoreVarCellFace%scalar_eta_mass_flux_n%f_r4(1:nlev,:)+&
                                                               0.5_r4*dycoreVarCellFace%scalar_eta_mass_flux_n%f_r4(2:nlev+1,:)
#endif

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute PS and face level V mass flux",nh_dynamics

!================================================
! diag FACE LEVEL VARS
!================================================

!$omp parallel  private(ilev)
!$omp do schedule(static,5)
        do ilev = 2, nlev
             dycoreVarEdgeFace%scalar_normal_mass_flux_n%f_r4(ilev,:) = &
                                     0.5*(deta_full(ilev-1)/deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4(ilev,:)+&
                                          deta_full(ilev)  /deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4(ilev-1,:))
        end do
!$omp end do nowait
!$omp end parallel

        dycoreVarEdgeFace%scalar_normal_mass_flux_n%f_r4(1,:)     = 0.5_r4*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4(1,:)
        dycoreVarEdgeFace%scalar_normal_mass_flux_n%f_r4(nlevp,:) = 0.5_r4*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4(nlev,:)

        mesh%nv = mesh%nv_halo(1)  ! only used by nh_dynamics originally, now also used for diagnosing dpidry/dt
        call divergence_operator_2d(mesh,dycoreVarEdgeFace%scalar_normal_mass_flux_n%f_r4, &
                                         dycoreVarCellFace%tend_mass_hori%f_r4, & ! face-level adv
                                         nlev+1)
        mesh%nv = mesh%nv_compute

        scalar_hpressure_at_pc_surface_np1%f = dycoreVarSurface%scalar_hpressure_n%f+&
#ifdef DBL_MASS
                                               dtime/rk_number*dycoreVarSurface%tend_hpressure_cnty%f
#else
                                               dtime/rk_number*dycoreVarSurface%tend_hpressure_cnty%f_r4
#endif
#ifndef SEQ_GRIST
        call exchange_data_1d_add(mesh,field_head_1d,scalar_hpressure_at_pc_surface_np1)
        call exchange_data_1d(mesh%local_block,field_head_1d)
#endif
!
! renew mass/hpressure state
!
        !mesh%nv = mesh%nv_halo(1)
        mesh%nv = mesh%nv_full
        call time_integration_renew_mass_state_new(mesh, scalar_hpressure_at_pc_surface_np1   ,&
                                                     scalar_hpressure_at_pc_face_level_np1,& ! overwritten each RK step
                                                     dycoreVarCellFull%scalar_delhp_np1    ,& ! overwritten each RK step
                                                     scalar_hpressure_at_pc_full_level_np1,& ! overwritten each RK step
                                                     scalar_delhp_at_pc_face_level_np1  )    ! overwritten each RK step
! get full mass, if dycore, just hpres
        call calc_hpe_get_full_mass(nlev, mesh%nv_halo(1), working_mode            , & ! in
                                          scalar_hpressure_at_pc_face_level_np1%f_r4  , & ! in
                                          scalar_hpressure_at_pc_full_level_np1%f  , & ! in
                                          dycoreVarCellFull%scalar_delhp_np1%f     , & ! in
                                          tracerVarCellFull%scalar_mif_n%f         , & ! in
                                          dycoreVarCellFull%scalar_mpressure_n%f   , & ! out
                                          dycoreVarCellFace%scalar_mpressure_n%f_r4)      ! out
#ifndef SEQ_GRIST
        call clock_end(clock_mas)
#endif

!--------------------------------------------------------------------------
! compute horizontal potential temperature (pt) flux and pt mass flux div
!--------------------------------------------------------------------------

!
! begin pt(m)-flux computing
!
#ifndef SEQ_GRIST
        call clock_begin(clock_ptm)
#endif
        call calc_primal_normal_flux_at_edge(mesh, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4, &
                                                   dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4, &
                                                   dycoreVarCellFull%scalar_potential_temp_n%f    , &
                                                   dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f, &
                                                   pot_adv_flag(irk_step), rk_substep, nlev)

        scalar_template_2d_ne_b%f = dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f*&
                                    dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4
!
! exchange scheme
!
#ifndef SEQ_GRIST
        call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_b)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

        mesh%ne = mesh%ne_halo(1)
        mesh%nv = mesh%nv_halo(1)
!inline_opt_here
      !  call divergence_operator_2d(mesh,scalar_template_2d_ne_b%f, &
      !                              dycoreVarCellFull%tend_mass_pt_hori%f, &
      !                              nlev)
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!$omp do schedule(dynamic,5)
        do iv = 1, mesh%nv
           div_sum(:) = 0._r8
           do inb = 1, mesh%vtx_nnb(iv)
              index_edge  = mesh%vtx_ed(inb,iv)
              do ilev = 1, nlev
                 div_sum(ilev)  = div_sum(ilev)+scalar_template_2d_ne_b%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
              end do
           end do
           dycoreVarCellFull%tend_mass_pt_hori%f(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
        end do
!$omp end do nowait
!$omp end parallel
!inline_opt_here
        dycoreVarCellFull%tend_mass_pt_hori%f = -1._r8*dycoreVarCellFull%tend_mass_pt_hori%f

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute horizontal pt flux"
        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute
!
! compute vertical potential temperature mass flux tendency ( partial (m*etadot*theta) )
! 
        call calc_hpe_tend_vert_mass_flux_2d(mesh%nv_full,mesh%nv_halo(1),nlev, & ! in
                                             dycoreVarCellFull%scalar_potential_temp_n%f , & ! in
#ifdef DBL_MASS
                                             dycoreVarCellFace%scalar_eta_mass_flux_n%f  , & ! in
#else
                                             dycoreVarCellFace%scalar_eta_mass_flux_n%f_r4  , & ! in
#endif
                                             ver_adv_flag                               , & ! in
                                             dycoreVarCellFull%tend_mass_pt_vert%f) ! out
!
! mass*pt update
!
        scalar_mass_pt_at_pc_full_level_np1%f  = dycoreVarCellFull%scalar_mass_pt_n%f+&
                                                 dtime/rk_number*&
                                                (dycoreVarCellFull%tend_mass_pt_hori%f+&
                                                 dycoreVarCellFull%tend_mass_pt_vert%f)
        if(use_phys.and.ptend_heat_rk_on)then
            scalar_mass_pt_at_pc_full_level_np1%f = scalar_mass_pt_at_pc_full_level_np1%f+&
                                                  dtime/rk_number*ptend_rk%tend_mass_pt_at_pc_full_level%f
        end if
        if(ad_dycore_laplacian_2nd)then
            scalar_mass_pt_at_pc_full_level_np1%f = scalar_mass_pt_at_pc_full_level_np1%f+&
                                                    dtime/rk_number*dycoreVarCellFull%tend_mass_pt_laplacian_2nd%f
        end if
        dycoreVarCellFull%scalar_potential_temp_n%f = scalar_mass_pt_at_pc_full_level_np1%f/& ! overwritten each RK step
                                                      dycoreVarCellFull%scalar_delhp_np1%f
#ifndef SEQ_GRIST
        call clock_end(clock_ptm)
#endif

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute vertical pt flux"

!----------------------------------------------------------------------
! compute velocity tendency from Coriolis and Kinetic Energy gradient
!----------------------------------------------------------------------
#ifndef SEQ_GRIST
        call clock_begin(clock_nct)
#endif

        IF((tend_nct_once.and.irk_step.eq.1).or.(.not.tend_nct_once)) &   ! control the frequency of NCT evaluation

        call calc_tend_nct_at_edge_full_level(mesh, scalar_delhp_at_pc_full_level_rk            ,&  ! in
                                                    scalar_normal_velocity_at_edge_full_level_rk,&  ! in
                                                    dycoreVarEdgeFull%scalar_normal_mass_flux_n ,&  ! in
                                                    dycoreVarEdgeFull%tend_normal_velocity_nct  ,&  ! out
                                                    dtime, nlev , irk_step)
#ifndef SEQ_GRIST
        call clock_end(clock_nct)
#endif

        mesh%nv = mesh%nv_halo(1)
        mesh%nt = mesh%nt_halo(1)
!
! KE definition of its gradient
!
#ifndef SEQ_GRIST
        call clock_begin(clock_ket)
#endif
        call calc_grad_kinetic_energy(mesh,scalar_normal_velocity_at_edge_full_level_rk%f_r4,&
                                           dycoreVarEdgeFull%tend_normal_velocity_ke%f_r4  ,&
                                           nlev)
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
#ifdef DBL_MASS
          0.5_r8*(dycoreVarCellFace%scalar_eta_mass_flux_n%f(:,icell1)+&
                  dycoreVarCellFace%scalar_eta_mass_flux_n%f(:,icell2))
#else
          0.5_r4*(dycoreVarCellFace%scalar_eta_mass_flux_n%f_r4(:,icell1)+&
                  dycoreVarCellFace%scalar_eta_mass_flux_n%f_r4(:,icell2))
#endif
          
          scalar_template_1d_nlev_a  = &
          0.5_r4*(scalar_delhp_at_pc_full_level_rk%f_r4(:,icell1)+&
                  scalar_delhp_at_pc_full_level_rk%f_r4(:,icell2))

          scalar_template_1d_nlev_b(:) = scalar_normal_velocity_at_edge_full_level_rk%f_r4(:,ie)

          call calc_hpe_vert_advection(scalar_template_1d_nlevp_a, &
                                       scalar_template_1d_nlev_a , &
                                       scalar_template_1d_nlev_b , &
                                       scalar_template_1d_nlev_c )

          dycoreVarEdgeFull%tend_normal_velocity_vadv%f_r4(:,ie) = scalar_template_1d_nlev_c
      
        end do
!$omp end do nowait
!$omp end parallel 

#ifndef SEQ_GRIST
        call clock_end(clock_vau)
#endif

        if(iblock .eq. 0.and.write_verbose) print*,"finish vertical advection of normal wind"

!===========================================================================================
! diagnose dpim/dt (omegaMoist)
! Besides being a diagnostics:
! for HDC this diag is necessary for WRF physics; for NDC this diag produces similar results
! as we use www, g and rhom for diag, but not itentical even in logic;
! 2021-01: to use this for ndc's cu scheme make tropical rainfall looks better real-world 
!===========================================================================================
#ifndef SEQ_GRIST
        call clock_begin(clock_diagom)
#endif
        if(irk_step.eq.nrk)then ! diagnose d(mpressure)/dt for this dycore_step (final rk step)
           dycoreVarCellFace%scalar_mpressure_n%f = dycoreVarCellFace%scalar_mpressure_n%f_r4

           mesh%ne = mesh%ne_compute
! pi flux
! inline_opt_here
           !call calc_primal_normal_flux_at_edge(mesh, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, & ! dum
           !                                           dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, & ! dum
           !                                           scalar_mpressure_at_pc_full_level_old         , & ! halo(1)->ne
           !                                           scalar_mpressure_at_edge_full_level_n         , &
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
           scalar_mpressure_at_edge_full_level_n = scalar_mpressure_at_edge_full_level_n*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f_r4

! inline_opt_here
           !call divergence_operator_2d(mesh,scalar_mpressure_at_edge_full_level_n, &
           !                                 tend_pidpiv_at_pc_full_level_hori    , &
           !                                 nlev)
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!$omp do schedule(dynamic,5)
           do iv = 1, mesh%nv
              div_sum(:) = 0._r8
              do inb = 1, mesh%vtx_nnb(iv)
                 index_edge  = mesh%vtx_ed(inb, iv)
                 do ilev = 1, nlev
                    div_sum(ilev)  = div_sum(ilev)+scalar_mpressure_at_edge_full_level_n(ilev,index_edge)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(index_edge)
                 end do
              end do
              tend_pidpiv_at_pc_full_level_hori(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
           end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here

           mesh%nv = mesh%nv_compute
            do iv = 1, mesh%nv
               do ilev = 1, nlev
                  dycoreVarCellFull%scalar_omega_n%f(ilev,iv) = (dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)-scalar_mpressure_at_pc_full_level_old(ilev,iv))/dtime+&
#ifdef DBL_MASS
                                                                 dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)+&
                                                                 ((tend_pidpiv_at_pc_full_level_hori(ilev,iv)-dycoreVarCellFull%tend_mass_hori%f(ilev,iv)*&
#else
                                                                 dycoreVarCellFull%scalar_eta_mass_flux_n%f_r4(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)+&
                                                                 ((tend_pidpiv_at_pc_full_level_hori(ilev,iv)-dycoreVarCellFull%tend_mass_hori%f_r4(ilev,iv)*&
#endif
                                                              scalar_mpressure_at_pc_full_level_old(ilev,iv))/scalar_delhp_at_pc_full_level_rk%f_r4(ilev,iv))
                  if(eta_full_b(ilev).eq.0)then
#ifdef DBL_MASS
                     dycoreVarCellFull%scalar_omega_n%f(ilev,iv) = dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)
#else
                     dycoreVarCellFull%scalar_omega_n%f(ilev,iv) = dycoreVarCellFull%scalar_eta_mass_flux_n%f_r4(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)
#endif
                  end if
               end do
            end do
        mesh%nv = mesh%nv_halo(1)
        end if
#ifndef SEQ_GRIST
        call clock_end(clock_diagom)
#endif
!
! Until this line, HYDROSTATIC and NONHYDROSTATIC should be the same
!
#ifndef SEQ_GRIST
        call clock_begin(clock_nhd) ! either nhd or hpe
#endif
        IF(nh_dynamics)then

           if(iblock .eq. 0.and.write_verbose) print*,"before nh dynamics"

           dycoreVarCellFace%scalar_delp_n%f(2:nlev,:)  = dycoreVarCellFull%scalar_pressure_rk%f(2:nlev,:)-dycoreVarCellFull%scalar_pressure_rk%f(1:nlev-1,:)
           dycoreVarCellFace%scalar_delp_n%f(1,:)       = dycoreVarCellFull%scalar_pressure_rk%f(1,:)     -scalar_pressure_at_pc_face_level_rk%f_r4(1,:)
           dycoreVarCellFace%scalar_delp_n%f(nlev+1,:)  = scalar_pressure_at_pc_face_level_rk%f_r4(nlev+1,:) -dycoreVarCellFull%scalar_pressure_rk%f(nlev,:)
   
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
        !    dycoreVarCellFull%scalar_eta_mass_flux_n%f_r4 = dycoreVarCellFull%scalar_eta_mass_flux_n%f
        !    dycoreVarEdgeFace%scalar_normal_mass_flux_n%f_r4 = dycoreVarEdgeFace%scalar_normal_mass_flux_n%f
        !    dycoreVarCellFace%tend_mass_hori%f_r4 = dycoreVarCellFace%tend_mass_hori%f
           dycoreVarCellFace%scalar_www_n%f_r4 = dycoreVarCellFace%scalar_www_n%f
           tracerVarCellFace%scalar_mif_n%f_r4 = tracerVarCellFace%scalar_mif_n%f
           dycoreVarCellFull%scalar_delhp_np1%f_r4 = dycoreVarCellFull%scalar_delhp_np1%f
        !    dycoreVarCellFace%scalar_phi_n%f_r4 = dycoreVarCellFace%scalar_phi_n%f
           call grist_nh_dynamics_run(mesh,rk_substep, itimestep, irk_step, idstep, itstep, &
                                         scalar_hpressure_at_pc_face_level_np1,         & !
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


          if(iblock .eq. 0.and.write_verbose) print*,"end nh dynamics"
!
! optionally to obtain a hydrostatic balance state after minsteps for every-intsteps, or based on dynamic control
!
           IF((restore_hydro.and.itimestep.gt.restore_hydro_minsteps.and.mod(itimestep,restore_hydro_intsteps).eq.0).or.& ! mannually control
               ndc_restore_flag.eq.1)then ! dynamic control
! verified: if activate this under ndc every timestep, it produces exactly identical results as use a HDC,
! except some diagnostic variables and www
              call hydrostatic_adjust
              if(iblock.eq.0.and.write_stepinfo) print*,"instant restored at step,", itimestep
           ENDIF

           IF(adjphi.and.itimestep.lt.16)then
              call hydrostatic_adjust
              !if(iblock.eq.0.and.write_stepinfo) print*,"restore hydro during model spinup", itimestep
           ENDIF

        ELSE
           call hydrostatic_adjust
           if(iblock .eq. 0.and.write_verbose) PRINT*,"finish diagnosing geopotential in HDC"

        END IF
#ifndef SEQ_GRIST
        call clock_end(clock_nhd)
#endif
      dycoreVarCellFace%scalar_pressure_n%f = dycoreVarCellFace%scalar_pressure_n%f_r4

!=================================================================================
! use below subroutines with smallest modification to compute PGF for normal
! velocity; The below part should not change too much in the case of nhdc or hdc
!=================================================================================
#ifndef SEQ_GRIST
      call clock_begin(clock_pgf)
#endif
      if(iblock .eq. 0.and.write_verbose) print*,"finish geop, grad_delp, grad_pressure gradient"
!
! compute real PGF tendency
!
      select case(hor_pgf_flag)

      case(6)  ! ptb, default with tp1 profile

          call grist_dycore_ref_atmos_create_ptb

!$omp parallel  private(ie,icell1,icell2,cr_num,tmp1,tmp2,tmp3,tmp4)    
!$omp do schedule(dynamic,5)

          do ie = 1, mesh%ne
             icell1 = mesh%edt_v(1,ie)
             icell2 = mesh%edt_v(2,ie)
             cr_num = rearth*mesh%edt_leng(ie)

             tmp1   = scalar_grad_hpres_at_edge_full_level_bar(:,ie)*&
                      0.5_r8*(scalar_alphad_at_pc_full_level_ptb(:,icell1)+scalar_alphad_at_pc_full_level_ptb(:,icell2))

             tmp2   = 0.5_r8*(dycoreVarCellFull%scalar_alpha_np1%f(:,icell1)+dycoreVarCellFull%scalar_alpha_np1%f(:,icell2))*&
                             (scalar_pressure_at_pc_full_level_ptb(:,icell2)-scalar_pressure_at_pc_full_level_ptb(:,icell1))/cr_num

             tmp3   = (scalar_geop_at_pc_full_level_ptb(:,icell2)-scalar_geop_at_pc_full_level_ptb(:,icell1))/cr_num

             tmp4   = 0.5_r8*(scalar_delp_at_pc_full_level_ptb(:,icell1)/dycoreVarCellFull%scalar_delhp_np1%f(:,icell1)+&
                              scalar_delp_at_pc_full_level_ptb(:,icell2)/dycoreVarCellFull%scalar_delhp_np1%f(:,icell2))*&
                             (dycoreVarCellFull%scalar_geopotential_n%f(:,icell2)-dycoreVarCellFull%scalar_geopotential_n%f(:,icell1))/cr_num

             dycoreVarEdgeFull%tend_normal_velocity_pgf%f_r4(:,ie) = -1._r8*(tmp1+tmp2+tmp3+tmp4)*tracerVarEdgeFull%scalar_mif_n%f(:,ie)
          end do
!$omp end do nowait
!$omp end parallel 
      case default
          if(iblock .eq. 0) print*, "you must select a PGF term in hor_pgf_flag"
          stop
      end select

      if(iblock .eq. 0.and.write_verbose) print*,"finish PGF"
#ifndef SEQ_GRIST
      call clock_end(clock_pgf)

      call clock_begin(clock_fnlrk)
#endif
!
! normal velocity updatee
!
        call clock_begin(clock_fnlrk1)
         scalar_normal_velocity_at_edge_full_level_np1%f = &
         dycoreVarEdgeFull%scalar_normal_velocity_n%f +dtime_r4/rk_number*&
       ( dycoreVarEdgeFull%tend_normal_velocity_vadv%f_r4 +&
         dycoreVarEdgeFull%tend_normal_velocity_pgf%f_r4  +&
         dycoreVarEdgeFull%tend_normal_velocity_ke%f_r4   +&
         dycoreVarEdgeFull%tend_normal_velocity_nct%f_r4)

         if(use_phys.and.ptend_wind_rk_on)then
             scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime/rk_number*ptend_rk%tend_normal_velocity_at_edge_full_level%f
         end if
!
! diffusion terms added here as phys, must be written in this way for bit regression
!
         if(ad_dycore_laplacian_2nd.and. .not.ad_dycore_laplacian_4th)then
             scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime_r4/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_2nd%f_r4
         else if(ad_dycore_laplacian_4th .and. .not. ad_dycore_laplacian_2nd)then
             scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                   dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_4th%f
         else if(ad_dycore_laplacian_2nd.and.ad_dycore_laplacian_4th)then
             ! True
             scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
                                                dtime_r4/rk_number*(dycoreVarEdgeFull%tend_hwind_laplacian_4th%f_r4+&
                                                                    dycoreVarEdgeFull%tend_hwind_laplacian_2nd%f_r4)
         end if

#if (!defined DCMIP21)
         if(nh_dynamics.and.ad_dycore_laplacian_2nd)then
             scalar_www_at_pc_face_level_np1%f = scalar_www_at_pc_face_level_np1%f+&
                                                 dtime_r4/rk_number*dycoreVarCellFace%tend_www_laplacian_2nd%f_r4
         end if
#endif
        call clock_end(clock_fnlrk1)
        call clock_begin(clock_fnlrk2)
!
! diff 6th is only for testing, not used for all cases
!
        !  if(ad_dycore_laplacian_6th)then
        !      scalar_normal_velocity_at_edge_full_level_np1%f = scalar_normal_velocity_at_edge_full_level_np1%f+&
        !                                            dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_6th%f
        !  end if
!
! RK-sub-step update
!

!exchange data
#ifndef SEQ_GRIST
       call clock_begin(clock_mainexch)
       call exchange_data_2d_add(mesh,field_head_2d,scalar_mass_pt_at_pc_full_level_np1)
       call exchange_data_2d_add(mesh,field_head_2d,scalar_normal_velocity_at_edge_full_level_np1)
#endif

       if(nh_dynamics)then
#ifndef SEQ_GRIST
          call exchange_data_2d_add(mesh,field_head_2d,scalar_www_at_pc_face_level_np1)
          call exchange_data_2d_add(mesh,field_head_2d,scalar_phi_at_pc_face_level_np1)
#endif
       end if
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!       call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_delhp_np1)
!       call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_pressure_np1)
!       call exchange_data_2d_add(mesh,field_head_2d,scalar_delhp_at_pc_face_level_np1)
!       call exchange_data_2d_add(mesh,field_head_2d,scalar_hpressure_at_pc_face_level_np1)
!zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
#ifndef SEQ_GRIST
       call exchange_data_2d(mesh%local_block,field_head_2d)
       call clock_end(clock_mainexch)
#endif
        call clock_end(clock_fnlrk2)
        call clock_begin(clock_fnlrk3)

       dycoreVarCellFull%scalar_potential_temp_n%f = scalar_mass_pt_at_pc_full_level_np1%f/& ! overwritten each RK step
                                                     dycoreVarCellFull%scalar_delhp_np1%f

! prog
      !scalar_hpressure_at_pc_surface_rk            = scalar_hpressure_at_pc_surface_np1
      !scalar_mass_pt_at_pc_full_level_rk           = scalar_mass_pt_at_pc_full_level_np1
      if(nh_dynamics)then
         scalar_www_at_pc_face_level_rk%f_r4            = scalar_www_at_pc_face_level_np1%f
         scalar_phi_at_pc_face_level_rk%f_r4            = scalar_phi_at_pc_face_level_np1%f
         dycoreVarCellFull%scalar_pressure_rk      = dycoreVarCellFull%scalar_pressure_np1
         scalar_pressure_at_pc_face_level_rk%f_r4       = dycoreVarCellFace%scalar_pressure_n%f_r4
      end if
      scalar_normal_velocity_at_edge_full_level_rk%f_r4 = scalar_normal_velocity_at_edge_full_level_np1%f
!       scalar_delhp_at_pc_full_level_rk             = dycoreVarCellFull%scalar_delhp_np1
      scalar_delhp_at_pc_full_level_rk%f_r4             = dycoreVarCellFull%scalar_delhp_np1%f
! diag
      scalar_delhp_at_pc_face_level_rk             = scalar_delhp_at_pc_face_level_np1
      !dycoreVarCellFull%scalar_potential_temp_n     = dycoreVarCellFull%scalar_potential_temp_np1
      !scalar_hpressure_at_pc_face_level_n          = scalar_hpressure_at_pc_face_level_np1

      call clock_end(clock_fnlrk3)
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
#ifndef SEQ_GRIST
      call clock_begin(clock_update)
#endif

      if(nh_dynamics)then
         dycoreVarCellFace%scalar_www_n     = scalar_www_at_pc_face_level_np1
         dycoreVarCellFace%scalar_phi_n      = scalar_phi_at_pc_face_level_np1
         dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_pressure_np1
      end if
      dycoreVarEdgeFull%scalar_normal_velocity_n  = scalar_normal_velocity_at_edge_full_level_np1
      dycoreVarCellFull%scalar_mass_pt_n          = scalar_mass_pt_at_pc_full_level_np1
      dycoreVarSurface%scalar_hpressure_n         = scalar_hpressure_at_pc_surface_np1
      dycoreVarCellFull%scalar_delhp_n            = dycoreVarCellFull%scalar_delhp_np1    ! needed by physics and dtp
      dycoreVarCellFace%scalar_delhp_n            = scalar_delhp_at_pc_face_level_np1
!
! this is only for diagnose PS, not put to dycore_diag because of historry-h1 does not dycore_diag at each step
      dycoreVarSurface%scalar_pressure_n%f = dycoreVarCellFace%scalar_pressure_n%f_r4(nlevp,:)
!
! substract diabatic physical influence, used with F2, so we should substract F2
! here; for idealized test, f2 has the same tendency, while for full-physics
! test, rk has its own total tendency, only the portion that f2 covers should be
! substracted
!
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

!================================================
! pt now only affected by adv and diffusion, get 
! this tend and accumulate until this itimestep
! reachs step_cu and all dycore+tracer steps done
! this tend is only for cumulus tiedtke
!================================================

#ifdef AMIPW_PHYSICS
      IF(present(idstep))then
      dycoreVarCellFull%tend_pt_n%f = dycoreVarCellFull%tend_pt_n%f+&
                                     (dycoreVarCellFull%scalar_potential_temp_n%f/(one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(1,:,:))-&
                                     scalar_pt_at_pc_full_level_old)/dtime
      if(idstep.eq.dstep_in_tstep.and.itstep.eq.tstep_in_mstep.and.mod(itimestep,step_cu).eq.0)then
         dycoreVarCellFull%tend_pt_n%f = dycoreVarCellFull%tend_pt_n%f/(dstep_in_tstep*tstep_in_mstep*step_cu)
      end if
      END IF
#endif

!=====================================================================================
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
!=====================================================================================

! #ifdef USE_PTENDF1
!       if(use_phys.and.ptend_dycore_heat_f1_on)then
!          dycoreVarCellFull%scalar_mass_pt_n%f           = dycoreVarCellFull%scalar_mass_pt_n%f+&
!                                                          ptend_f1%tend_mass_pt_at_pc_full_level%f*dtime
! #ifndef SEQ_GRIST
!          call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
!          call exchange_data_2d(mesh%local_block,field_head_2d)
! #endif
!          dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
!       end if

!       if(use_phys.and.ptend_dycore_wind_f1_on)then
!          dycoreVarEdgeFull%scalar_normal_velocity_n%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f+&
!                                                          ptend_f1%tend_normal_velocity_at_edge_full_level%f*dtime
! #ifndef SEQ_GRIST
!          call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
!          call exchange_data_2d(mesh%local_block,field_head_2d)
! #endif
!       end if

!       if(use_phys.and.ptend_dycore_f1_on)then
!          dycoreVarCellFull%scalar_mass_pt_n%f           = dycoreVarCellFull%scalar_mass_pt_n%f+&
!                                                          ptend_f1%tend_mass_pt_at_pc_full_level%f*dtime
!          dycoreVarEdgeFull%scalar_normal_velocity_n%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f+&
!                                                          ptend_f1%tend_normal_velocity_at_edge_full_level%f*dtime
! #ifndef SEQ_GRIST
!          call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
!          call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
!          call exchange_data_2d(mesh%local_block,field_head_2d)
! #endif
!          dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
!       end if
! #endif
#ifndef SEQ_GRIST
      call clock_end(clock_update)
#endif
      return

contains

    subroutine hydrostatic_adjust
!
! the modifiction of this part for mhdc may affect bit reproduce of dry hdc,
! although physically identical
!

#ifndef REG_OLD_NHDAMP
        if(nh_dynamics)then
! hold old value
           scalar_pressure_at_pc_face_level_old    = dycoreVarCellFace%scalar_pressure_n%f_r4
           scalar_pressure_at_pc_full_level_old    = dycoreVarCellFull%scalar_pressure_np1%f
        end if
#endif

        dycoreVarCellFull%scalar_pressure_np1%f    = dycoreVarCellFull%scalar_mpressure_n%f
        dycoreVarCellFace%scalar_pressure_n%f_r4   = dycoreVarCellFace%scalar_mpressure_n%f_r4
!
! diagnose geopotential and evaluate its gradient
!
!$omp parallel  private(iv,scalar_template_1d_nlev_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_b,scalar_template_a,scalar_template_1d_nlevp_b,scalar_template_1d_nlev_c)
!$omp do schedule(dynamic,5) 
        do iv = 1, mesh%nv_halo(1)
          !
          ! recover temperature, and store
          !
           !scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_np1%f(:,iv)*&
           scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)*&
                                      ((dycoreVarCellFull%scalar_pressure_np1%f(:,iv)/p00)**(rdry/cp))/&
                                       (one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv))

           dycoreVarCellFull%scalar_temp_n%f(:,iv) = scalar_template_1d_nlev_a ! output, overwritten each RK step
           scalar_template_1d_nlev_a  = scalar_template_1d_nlev_a*(one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)))) ! Tv
          !
          ! new pressure face
          !
           !scalar_template_1d_nlevp_a = scalar_hpressure_at_pc_face_level_n%f(:,iv)
           scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f_r4(:,iv)
          !
          ! new delp full, not in effect actually
          !
           !scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delhp_np1%f(:,iv)
           scalar_template_1d_nlev_b  = dycoreVarCellFace%scalar_pressure_n%f_r4(2:nlevp,iv)-dycoreVarCellFace%scalar_pressure_n%f_r4(1:nlev,iv)
          !
          ! geopotential at surface
          !
           scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)
          ! diagnose geop at face based on pressure, which is related to
          ! the definition of model layer
           call calc_hpe_hydro(scalar_template_1d_nlev_a  ,& ! T or Tv
                               scalar_template_1d_nlevp_a ,& ! hp or p
                               scalar_template_1d_nlev_b  ,& ! delhp or delp
                               scalar_template_a          ,& ! phis
                               scalar_template_1d_nlevp_b ,& ! phi-face
                               scalar_template_1d_nlev_c )   ! phi-full

           dycoreVarCellFace%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlevp_b
           dycoreVarCellFull%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlev_c  ! output

#ifdef REG_OLD_NHDAMP
           scalar_phi_at_pc_face_level_np1%f(:,iv)        = dycoreVarCellFace%scalar_geopotential_n%f(:,iv)
#endif
        end do
!$omp end do nowait
!$omp end parallel

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

           dycoreVarCellFace%scalar_pressure_n%f_r4= scalar_pressure_at_pc_face_level_old+&
                                          damp_kesi*(scalar_pressure_at_pc_face_level_old-dycoreVarCellFace%scalar_mpressure_n%f_r4)
        end if
#endif


          ! since in hdc, delp is only used for pgf, and we use (6), so delp should be delhp, to cancel the last term in pgf6, not so !
          ! dycoreVarCellFull%scalar_delp_np1%f   = dycoreVarCellFull%scalar_delhp_np1%f
          ! new
          dycoreVarCellFull%scalar_delp_np1%f_r4   = dycoreVarCellFace%scalar_pressure_n%f_r4(2:nlevp,:)-dycoreVarCellFace%scalar_pressure_n%f_r4(1:nlev,:) 
          !scalar_kesi_at_pc_full_level_np1%f   = 1._r8
          ! old dry-hydrostatic model only uses this
          !dycoreVarCellFull%scalar_alpha_np1%f  = rdry*dycoreVarCellFull%scalar_temp_n%f/scalar_hpressure_at_pc_full_level_np1%f
          ! real layer-averaged
          dycoreVarCellFull%scalar_alpha_np1%f  = (dycoreVarCellFace%scalar_geopotential_n%f(1:nlev,:)-dycoreVarCellFace%scalar_geopotential_n%f(2:nlevp,:))/&
                                                  dycoreVarCellFull%scalar_delhp_np1%f
          return
    end subroutine hydrostatic_adjust

   end subroutine dycore_time_integration_run

   subroutine dycore_time_integration_init(mesh)
! dycore-nh
     use grist_nh_driver_module_mixed,           only: grist_nh_dynamics_init
     use grist_nh_explicit_tend_module_2d_mixed, only: grist_nh_et_init
! io
     type(global_domain),  intent(inout) :: mesh
! local
     integer(i4)      :: iv

     call grist_nh_dynamics_init(mesh)
     call grist_nh_et_init(mesh)
     call grist_dycore_hori_swe_clock_init

        if(.not.allocated(scalar_template_1d_nlev_a))  allocate(scalar_template_1d_nlev_a(nlev))
        if(.not.allocated(scalar_template_1d_nlev_b))  allocate(scalar_template_1d_nlev_b(nlev))
        if(.not.allocated(scalar_template_1d_nlev_c))  allocate(scalar_template_1d_nlev_c(nlev))
        if(.not.allocated(scalar_template_1d_nlev_d))  allocate(scalar_template_1d_nlev_d(nlev))
        if(.not.allocated(scalar_template_1d_nlevp_a)) allocate(scalar_template_1d_nlevp_a(nlevp))
        if(.not.allocated(scalar_template_1d_nlevp_b)) allocate(scalar_template_1d_nlevp_b(nlevp))
        if(.not.allocated(scalar_template_1d_nlevp_c)) allocate(scalar_template_1d_nlevp_c(nlevp))
        if(.not.allocated(scalar_template_2d_ne_b%f))  allocate(scalar_template_2d_ne_b%f(nlev,mesh%ne))
!
! initially compute hpressure (full&face), delhp (full&face), based on surface hpressure
!
        call time_integration_renew_mass_state(mesh, dycoreVarSurface%scalar_hpressure_n    ,&
                                                     dycoreVarCellFace%scalar_hpressure_n ,&
                                                     dycoreVarCellFull%scalar_delhp_n     ,&
                                                     dycoreVarCellFull%scalar_hpressure_n ,&
                                                     dycoreVarCellFace%scalar_delhp_n    )

        if(trim(working_mode).eq.'dycore')then
           dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_hpressure_n ! init as hpressure
           dycoreVarCellFace%scalar_pressure_n%f_r4 = dycoreVarCellFace%scalar_hpressure_n%f ! init as hpressure
           scalar_pressure_at_pc_face_level_rk%f_r4= dycoreVarCellFace%scalar_hpressure_n%f ! init as hpressure
        end if
! above for dycore mode

        if(trim(working_mode).ne.'dycore')then
           ! already defined in dtp
           !dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_pressure_n
           !dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_pressure_n
           dycoreVarCellFace%scalar_pressure_n%f_r4= dycoreVarCellFace%scalar_pressure_n%f
           scalar_pressure_at_pc_face_level_rk%f_r4= dycoreVarCellFace%scalar_pressure_n%f_r4
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
           scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)   ! geopotential at surface

           call calc_hpe_hydro(scalar_template_1d_nlev_a  ,&
                               scalar_template_1d_nlevp_a ,&
                               scalar_template_1d_nlev_b  ,&
                               scalar_template_a          ,&
                               scalar_template_1d_nlevp_b ,&
                               scalar_template_1d_nlev_c )

           dycoreVarCellFace%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlevp_b
           dycoreVarCellFull%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlev_c ! output
           dycoreVarCellFace%scalar_phi_n%f(:,iv)          = dycoreVarCellFace%scalar_geopotential_n%f(:,iv)
        end do
        dycoreVarCellFace%scalar_www_n%f   = 0._r8
        dycoreVarCellFace%scalar_phi_n%f(nlev+1,:)  = dycoreVarSurface%scalar_geopotential_n%f ! surface is given in ini condition
        END IF

        mesh%ne = mesh%ne_compute
        call grist_dycore_ref_atmos_init(mesh)
        mesh%ne = mesh%ne_full

!
! init mpressure
!
        call calc_hpe_get_full_mass1(nlev, mesh%nv_halo(1), working_mode          , & ! in
                                          dycoreVarCellFace%scalar_hpressure_n%f , & ! in
                                          dycoreVarCellFull%scalar_hpressure_n%f , & ! in
                                          dycoreVarCellFull%scalar_delhp_n%f     , & ! in
                                          tracerVarCellFull%scalar_mif_n%f       , & ! in
                                          dycoreVarCellFull%scalar_mpressure_n%f , & ! out
                                          dycoreVarCellFace%scalar_mpressure_n%f_r4)     ! out
! init array
      scalar_www_at_pc_face_level_np1        = dycoreVarCellFace%scalar_www_n 
      scalar_phi_at_pc_face_level_np1        = dycoreVarCellFace%scalar_www_n
      scalar_delhp_at_pc_face_level_np1      = dycoreVarCellFace%scalar_www_n
      scalar_hpressure_at_pc_face_level_np1%f_r4  = dycoreVarCellFace%scalar_www_n%f
      scalar_hpressure_at_pc_face_level_np1%pos   = dycoreVarCellFace%scalar_www_n%pos
      scalar_hpressure_at_pc_full_level_np1  = dycoreVarCellFull%scalar_hpressure_n
#ifndef SEQ_GRIST
      clock_prep     = clock_id('prep')
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
      clock_update   = clock_id('update')

      clock_fnlrk1   = clock_id('fnlrk1')
      clock_fnlrk2   = clock_id('fnlrk2')
      clock_fnlrk3   = clock_id('fnlrk3')
#endif

      return
   end subroutine dycore_time_integration_init

   subroutine dycore_time_integration_final

     use grist_nh_driver_module_mixed,             only: grist_nh_dynamics_final
     use grist_nh_explicit_tend_module_2d_mixed,   only: grist_nh_et_final

      if(allocated(scalar_template_1d_nlev_a))  deallocate(scalar_template_1d_nlev_a)
      if(allocated(scalar_template_1d_nlev_b))  deallocate(scalar_template_1d_nlev_b)
      if(allocated(scalar_template_1d_nlev_c))  deallocate(scalar_template_1d_nlev_c)
      if(allocated(scalar_template_1d_nlev_d))  deallocate(scalar_template_1d_nlev_d)
      if(allocated(scalar_template_1d_nlevp_a)) deallocate(scalar_template_1d_nlevp_a)
      if(allocated(scalar_template_1d_nlevp_b)) deallocate(scalar_template_1d_nlevp_b)
      if(allocated(scalar_template_1d_nlevp_c)) deallocate(scalar_template_1d_nlevp_c)
      if(allocated(scalar_template_2d_ne_b%f))  deallocate(scalar_template_2d_ne_b%f)

      call grist_nh_dynamics_final
      call grist_nh_et_final
      call grist_dycore_ref_atmos_final

   end subroutine dycore_time_integration_final

!--------------------
!
! PRIVATE BELOW
!
!--------------------

   subroutine time_integration_renew_mass_state(mesh, scalar_hpressure_at_pc_surface    ,&
                                                      scalar_hpressure_at_pc_face_level ,&
                                                      scalar_delhp_at_pc_full_level     ,&
                                                      scalar_hpressure_at_pc_full_level ,&
                                                      scalar_delhp_at_pc_face_level     )
! io
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     type(scalar_1d_field), intent(in)    :: scalar_hpressure_at_pc_surface
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_face_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_face_level
! local
     integer(i4)                          :: iv, ilev

!
! renew mass state, compute hpressure at face level, full level, 
! delhp, based on input surface hpressure
!

!$omp parallel  private(iv,scalar_template_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_a,scalar_template_1d_nlev_b,ilev) 
!$omp do schedule(dynamic,5) 
        do iv = 1, mesh%nv

           scalar_template_a   =   scalar_hpressure_at_pc_surface%f(iv)

           call calc_hpe_hpressure_face_level(scalar_template_a, &       ! in
                                             scalar_template_1d_nlevp_a) ! out

           call calc_hpe_delhp(scalar_template_1d_nlevp_a, &  ! in
                               scalar_template_1d_nlev_a)     ! out

           call calc_hpe_hpressure_full_level(scalar_template_a         , & ! in
                                              scalar_template_1d_nlevp_a, & ! in
                                              scalar_template_1d_nlev_a , & ! in
                                              scalar_template_1d_nlev_b)    ! out

           scalar_hpressure_at_pc_face_level%f(:,iv)  = scalar_template_1d_nlevp_a
               scalar_delhp_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_a
           scalar_hpressure_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_b

           do ilev = 2, nlev
              scalar_delhp_at_pc_face_level%f(ilev,iv)= scalar_hpressure_at_pc_full_level%f(ilev,iv)-&
                                                        scalar_hpressure_at_pc_full_level%f(ilev-1,iv)
           end do

           scalar_delhp_at_pc_face_level%f(1,iv)      = 0.5_r8*scalar_delhp_at_pc_full_level%f(1,iv)
           scalar_delhp_at_pc_face_level%f(nlev+1,iv) = 0.5_r8*scalar_delhp_at_pc_full_level%f(nlev,iv)

        end do
!$omp end do nowait
!$omp end parallel

        return
   end subroutine time_integration_renew_mass_state

   subroutine time_integration_renew_mass_state_new(mesh, scalar_hpressure_at_pc_surface    ,&
                                                          scalar_hpressure_at_pc_face_level ,&
                                                          scalar_delhp_at_pc_full_level     ,&
                                                          scalar_hpressure_at_pc_full_level ,&
                                                          scalar_delhp_at_pc_face_level     )
! io
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     type(scalar_1d_field), intent(in)    :: scalar_hpressure_at_pc_surface
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_face_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_face_level
! local
     integer(i4)                          :: iv, ilev

!
! renew mass state, compute hpressure at face level, full level, 
! delhp, based on input surface hpressure
!

!$omp parallel  private(iv,scalar_template_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_a,scalar_template_1d_nlev_b,ilev) 
!$omp do schedule(dynamic,5) 
        do iv = 1, mesh%nv

           scalar_template_a   =   scalar_hpressure_at_pc_surface%f(iv)

           call calc_hpe_hpressure_face_level(scalar_template_a, &       ! in
                                             scalar_template_1d_nlevp_a) ! out

           call calc_hpe_delhp(scalar_template_1d_nlevp_a, &  ! in
                               scalar_template_1d_nlev_a)     ! out

           call calc_hpe_hpressure_full_level(scalar_template_a         , & ! in
                                              scalar_template_1d_nlevp_a, & ! in
                                              scalar_template_1d_nlev_a , & ! in
                                              scalar_template_1d_nlev_b)    ! out

           scalar_hpressure_at_pc_face_level%f_r4(:,iv)  = scalar_template_1d_nlevp_a
               scalar_delhp_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_a
           scalar_hpressure_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_b

           do ilev = 2, nlev
              scalar_delhp_at_pc_face_level%f(ilev,iv)= scalar_hpressure_at_pc_full_level%f(ilev,iv)-&
                                                        scalar_hpressure_at_pc_full_level%f(ilev-1,iv)
           end do

           scalar_delhp_at_pc_face_level%f(1,iv)      = 0.5_r8*scalar_delhp_at_pc_full_level%f(1,iv)
           scalar_delhp_at_pc_face_level%f(nlev+1,iv) = 0.5_r8*scalar_delhp_at_pc_full_level%f(nlev,iv)

        end do
!$omp end do nowait
!$omp end parallel

        return
   end subroutine time_integration_renew_mass_state_new

end module grist_dycore_time_integration_2d_mixed
