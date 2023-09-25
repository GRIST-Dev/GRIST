
!----------------------------------------------------------------------------
! Created on 2018
! Author: Yi Zhang
! Version 1.0
! Description: The main driver module for the nh-dynamics; The explicit and
!              implicit part are seperated, and only the explicit part needs
!              horizontal data structure (i.e., parallel).
! Revision history:
!              1. add exploop and imploop options for regression usage
!              2. default use exploop
!----------------------------------------------------------------------------

   module grist_nh_driver_module_mixed
! constant
    use grist_constants,               only: rearth,i4, r8, fillvalue, gravity, p00, rdry, rvap, ptfactor, prfactor, one, cp, cv, pi, pi_r4 => pi_ns, zero, r4 => ns, one_r4 => one_ns, gravity_r4 => gravity_ns
    use grist_nml_module,              only: nlev, nlevp, phi_adv_flag, www_adv_flag, d3d_damping_coef, www_damping_coef,imbeta,zd, &
                                             write_verbose, write_stepinfo
! data definition
    use grist_data_types,              only: scalar_1d_field, scalar_2d_field, scalar_3d_field
    use grist_domain_types,            only: global_domain
    use grist_math_module,             only: extrapolate_bdy
! hori operator
! nh module
    use grist_nh_explicit_tend_module_2d_mixed, only: grist_nh_et_adv_face
    use grist_nh_implicit_module_2d_mixed,      only: grist_nh_solve_www_equation
! hpe
!     use grist_hpe_constants,              only: deta_full, deta_face,eta_full,eta_face
    use grist_dycore_module_vars_control_mixed, only: deta_full, deta_face,eta_full,eta_face,www_damping_coef_r4 => www_damping_coef, zd_r4 => zd
! mpi
#ifndef SEQ_GRIST
    use grist_mpi
#endif
    use grist_lib
#ifndef SEQ_GRIST
    use grist_config_partition,           only: debug_data_2d,debug_data_1d
#else
    use debug_module,                     only: debug_data_2d,debug_data_1d
#endif

#ifndef SEQ_GRIST
    use grist_clocks,                 only: clock_id, clock_begin, clock_end
#endif
    implicit none

      private
      public :: grist_nh_dynamics_init, &
                grist_nh_dynamics_final,&
                grist_nh_dynamics_run
!
! local 1d, hori
!
        real(r8), allocatable    :: scalar_www_top(:)
        real(r8), allocatable    :: scalar_www_bot(:)
        real(r8), allocatable    :: alpha_n(:)
        real(r8), allocatable    :: alpha_np1(:)         
        real(r8), parameter, public :: ndc_restore_num = 10._r8
        integer(i4), public         :: ndc_restore_flag
        integer(i4)                 :: itimestep_old

        integer(i4)     :: clock_nhd_exp, clock_nhd_imp, clock_nhd_diag

  CONTAINS
!
! main interface to time integration
!
     subroutine grist_nh_dynamics_run(mesh,dtime, itimestep, irk_step, idstep, itstep, &
                                      scalar_hpressure_at_pc_face_level_np1,         &
                                      scalar_normal_mass_flux_at_edge_face_level_rk, &
                                      tend_mass_at_pc_face_level_hori_rk,            &
                                      scalar_eta_mass_flux_at_pc_full_level_rk,      &
                                      scalar_eta_mass_flux_at_pc_face_level_rk,      &
                                      scalar_phi_at_pc_face_level_rk,                &
                                      scalar_www_at_pc_face_level_rk,                &
                                      scalar_delp_at_pc_face_level_rk,               &
                                      scalar_delhp_at_pc_face_level_rk,              &
                                      scalar_delhp_at_pc_face_level_n,               &
                                      scalar_pressure_at_pc_full_level_rk,           &
                                      scalar_geopotential_at_pc_surface_n,           &
                                      scalar_pressure_at_pc_full_level_n,            &
                                      scalar_phi_at_pc_face_level_n,                 &
                                      scalar_delhp_at_pc_full_level_n,               &
                                      scalar_www_at_pc_face_level_n,                 &
                                      scalar_potential_temp_at_pc_full_level_n,      &
                                      scalar_delhp_at_pc_full_level_np1,             &
                                      scalar_delhp_at_pc_face_level_np1,             &
                                      scalar_potential_temp_at_pc_full_level_np1,    &
                                      scalar_mif_at_pc_face_level_n,                 &
                                      scalar_mif_at_pc_full_level_n,                 &
                                      scalar_tracer_mxrt_at_pc_full_level_n,         &
                                      scalar_phi_at_pc_face_level_np1,               &  ! out
                                      scalar_www_at_pc_face_level_np1,               &  ! out
                                      scalar_temp_at_pc_full_level_np1,              &  ! out
                                      scalar_geopotential_at_pc_full_level_np1,      &  ! out
                                      scalar_geopotential_at_pc_face_level_np1,      &  ! out
                                      scalar_pressure_at_pc_full_level_np1,          &  ! out
                                      scalar_pressure_at_pc_face_level_np1,          &  ! out
                                      scalar_delp_at_pc_full_level_np1,              &  ! out
                                      scalar_alpha_at_pc_full_level_np1)                ! out
! io
        type(global_domain),   intent(inout) :: mesh
        real(r8)           ,   intent(in)    :: dtime
        integer(i4)        ,   intent(in)    :: itimestep, irk_step, idstep, itstep
        type(scalar_2d_field), intent(in)    :: scalar_hpressure_at_pc_face_level_np1
        type(scalar_2d_field), intent(in)    :: scalar_normal_mass_flux_at_edge_face_level_rk
        type(scalar_2d_field), intent(in)    :: tend_mass_at_pc_face_level_hori_rk
        type(scalar_2d_field), intent(in)    :: scalar_eta_mass_flux_at_pc_full_level_rk
        type(scalar_2d_field), intent(in)    :: scalar_eta_mass_flux_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)    :: scalar_phi_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)    :: scalar_www_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)    :: scalar_delp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_n
        type(scalar_2d_field), intent(in)    :: scalar_pressure_at_pc_full_level_rk
        type(scalar_1d_field), intent(in)    :: scalar_geopotential_at_pc_surface_n
        type(scalar_2d_field), intent(in)    :: scalar_pressure_at_pc_full_level_n
        type(scalar_2d_field), intent(in)    :: scalar_phi_at_pc_face_level_n
        type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_full_level_n
        type(scalar_2d_field), intent(in)    :: scalar_www_at_pc_face_level_n
        type(scalar_2d_field), intent(in)    :: scalar_potential_temp_at_pc_full_level_n
        type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_full_level_np1
        type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_np1
        type(scalar_2d_field), intent(in)    :: scalar_potential_temp_at_pc_full_level_np1
        type(scalar_2d_field), intent(in)    :: scalar_mif_at_pc_face_level_n
        type(scalar_2d_field), intent(in)    :: scalar_mif_at_pc_full_level_n
        type(scalar_3d_field), intent(in)    :: scalar_tracer_mxrt_at_pc_full_level_n
        type(scalar_2d_field), intent(inout) :: scalar_phi_at_pc_face_level_np1            ! out
        type(scalar_2d_field), intent(inout) :: scalar_www_at_pc_face_level_np1            ! out
        type(scalar_2d_field), intent(inout) :: scalar_temp_at_pc_full_level_np1           ! out
        type(scalar_2d_field), intent(inout) :: scalar_geopotential_at_pc_full_level_np1   ! out
        type(scalar_2d_field), intent(inout) :: scalar_geopotential_at_pc_face_level_np1   ! out
        type(scalar_2d_field), intent(inout) :: scalar_pressure_at_pc_full_level_np1       ! out
        type(scalar_2d_field), intent(inout) :: scalar_pressure_at_pc_face_level_np1       ! out
        type(scalar_2d_field), intent(inout) :: scalar_delp_at_pc_full_level_np1           ! out
        !type(scalar_2d_field), intent(inout) :: scalar_kesi_at_pc_full_level_np1           ! out
        type(scalar_2d_field), intent(inout) :: scalar_alpha_at_pc_full_level_np1          ! out
! local
        type(scalar_2d_field)                :: tend_et_phi_at_pc_face_level_rk
        type(scalar_2d_field)                :: tend_et_www_at_pc_face_level_rk
        !integer, save                        :: irk_step = 0

        !irk_step = irk_step + 1  ! yiz: already as input now

! just init
        tend_et_phi_at_pc_face_level_rk%f_r4   = scalar_pressure_at_pc_face_level_np1%f_r4
        tend_et_www_at_pc_face_level_rk%f_r4   = scalar_pressure_at_pc_face_level_np1%f_r4

#ifndef SEQ_GRIST
        call clock_begin(clock_nhd_exp)
#endif
        if(mpi_rank() == 0.and.write_verbose) print*,"grist_nh_dynamics_run_explicit"
        call grist_nh_dynamics_run_explicit(mesh, dtime, itimestep, irk_step             , &
                                            scalar_normal_mass_flux_at_edge_face_level_rk, &
                                            tend_mass_at_pc_face_level_hori_rk           , &
                                            scalar_eta_mass_flux_at_pc_full_level_rk     , &
                                            scalar_eta_mass_flux_at_pc_face_level_rk     , &
                                            scalar_phi_at_pc_face_level_rk               , &
                                            scalar_www_at_pc_face_level_n                , &
                                            scalar_www_at_pc_face_level_rk               , &
                                            scalar_delhp_at_pc_face_level_rk             , &
                                            scalar_delhp_at_pc_face_level_n              , &
                                            scalar_delhp_at_pc_face_level_np1            , &
                                            scalar_geopotential_at_pc_surface_n          , &
                                            tend_et_phi_at_pc_face_level_rk              , & ! need
                                            tend_et_www_at_pc_face_level_rk)                 ! need
#ifndef SEQ_GRIST
        call clock_end(clock_nhd_exp)
#endif
#ifndef SEQ_GRIST
        if(irk_step .eq. 1. .and. .false.)then
           call debug_data_2d(irk_step,19,mesh%e_index,mesh%ne_halo(1),"x1x",scalar_normal_mass_flux_at_edge_face_level_rk%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"x2x",tend_mass_at_pc_face_level_hori_rk%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"x3x",scalar_eta_mass_flux_at_pc_full_level_rk%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"x4x",scalar_phi_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"x5x",scalar_www_at_pc_face_level_n%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"x6x",scalar_www_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"x7x",scalar_delhp_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"x8x",tend_et_phi_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"x9x",tend_et_www_at_pc_face_level_rk%f)
           call debug_data_1d(irk_step,mesh%v_index,mesh%nv_halo(1),"x10x",scalar_geopotential_at_pc_surface_n%f)
        end if
#endif

#ifndef SEQ_GRIST
        call clock_begin(clock_nhd_imp)
#endif
        if(mpi_rank() == 0.and.write_verbose) print*,"grist_nh_dynamics_run_implicit"
        !need change halo(1)
        call grist_nh_dynamics_run_implicit(mesh, dtime,                               &
                                            scalar_pressure_at_pc_full_level_n,        &
                                            scalar_phi_at_pc_face_level_n,             &
                                            scalar_delhp_at_pc_full_level_n,           &
                                            scalar_www_at_pc_face_level_n,             &
                                            scalar_potential_temp_at_pc_full_level_n,  &
                                            scalar_eta_mass_flux_at_pc_full_level_rk,  &
                                            scalar_www_at_pc_face_level_rk,            &
                                            tend_et_www_at_pc_face_level_rk,           &
                                            tend_et_phi_at_pc_face_level_rk,           &
                                            scalar_potential_temp_at_pc_full_level_np1,&
                                            scalar_delhp_at_pc_full_level_np1,         &
                                            scalar_delhp_at_pc_face_level_np1,         &
                                            scalar_delp_at_pc_face_level_rk,           &
                                            scalar_delhp_at_pc_face_level_rk,          &
                                            scalar_mif_at_pc_face_level_n,             &
                                            scalar_mif_at_pc_full_level_n,             &
                                            scalar_www_at_pc_face_level_np1,           &
                                            scalar_phi_at_pc_face_level_np1)
                                            !scalar_kesi_at_pc_full_level_np1,         &
        scalar_www_at_pc_face_level_np1%f_r4 = scalar_www_at_pc_face_level_np1%f

        scalar_phi_at_pc_face_level_np1%f(nlev+1,:) = scalar_geopotential_at_pc_surface_n%f
#ifndef SEQ_GRIST
        call clock_end(clock_nhd_imp)
#endif
#ifndef SEQ_GRIST
        if(irk_step .eq. 1. .and. .false.)then
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"y1y",scalar_pressure_at_pc_full_level_rk%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"y2y",scalar_pressure_at_pc_full_level_n%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"y3y",scalar_phi_at_pc_face_level_n%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"y4y",scalar_phi_at_pc_face_level_np1%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"y5y",scalar_www_at_pc_face_level_np1%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"y6y",scalar_delhp_at_pc_full_level_n%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"y7y",scalar_delhp_at_pc_full_level_np1%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"y8y",scalar_potential_temp_at_pc_full_level_n%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"y9y",scalar_potential_temp_at_pc_full_level_np1%f)
           call debug_data_1d(irk_step,mesh%v_index,mesh%nv_halo(1),"y10y",scalar_geopotential_at_pc_surface_n%f)
        end if
#endif

#ifndef SEQ_GRIST
        call clock_begin(clock_nhd_diag)
#endif
        if(mpi_rank() == 0.and.write_verbose) print*," grist_nh_dynamics_run_diagvars"
        call grist_nh_dynamics_run_diagvars(mesh, itimestep, irk_step,idstep, itstep, dtime, &
                                                  scalar_delp_at_pc_face_level_rk,           &  ! in
                                                  scalar_delhp_at_pc_face_level_rk,          &  ! in
                                                  scalar_hpressure_at_pc_face_level_np1,     &  ! in
                                                  scalar_delhp_at_pc_face_level_np1,         &  ! in
                                                  scalar_pressure_at_pc_full_level_rk,       &  ! in
                                                  scalar_pressure_at_pc_full_level_n,        &  ! in
                                                  scalar_phi_at_pc_face_level_n,             &  ! in
                                                  scalar_phi_at_pc_face_level_np1,           &  ! in
                                                  scalar_www_at_pc_face_level_n,             &  ! in
                                                  scalar_www_at_pc_face_level_np1,           &  ! in
                                                  tend_et_www_at_pc_face_level_rk,           &  ! in
                                                  scalar_delhp_at_pc_full_level_n,           &  ! in
                                                  scalar_delhp_at_pc_full_level_np1,         &  ! in
                                                  scalar_potential_temp_at_pc_full_level_n  ,&  ! in
                                                  scalar_potential_temp_at_pc_full_level_np1,&  ! in
                                                  scalar_geopotential_at_pc_surface_n,       &  ! in
                                                  scalar_mif_at_pc_face_level_n,             &  ! in
                                                  scalar_mif_at_pc_full_level_n,             &  ! in
                                                  scalar_tracer_mxrt_at_pc_full_level_n,     &  ! in
                                                  scalar_temp_at_pc_full_level_np1,          &  ! out
                                                  scalar_geopotential_at_pc_full_level_np1,  &  ! out
                                                  scalar_geopotential_at_pc_face_level_np1,  &  ! out
                                                  scalar_pressure_at_pc_full_level_np1    ,  &  ! out
                                                  scalar_pressure_at_pc_face_level_np1    ,  &  ! out
                                                  scalar_delp_at_pc_full_level_np1,          &
                                                  !scalar_kesi_at_pc_full_level_np1,          &
                                                  scalar_alpha_at_pc_full_level_np1)
#ifndef SEQ_GRIST
        call clock_end(clock_nhd_diag)
#endif
#ifndef SEQ_GRIST
        if(irk_step .eq. 1 .and. .false.)then
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"k1k",scalar_temp_at_pc_full_level_np1%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"k2k",scalar_geopotential_at_pc_full_level_np1%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"k3k",scalar_geopotential_at_pc_face_level_np1%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"k4k",scalar_pressure_at_pc_full_level_np1%f)
           call debug_data_2d(irk_step,19,mesh%v_index,mesh%nv_halo(1),"k5k",scalar_pressure_at_pc_face_level_np1%f)
           call debug_data_2d(irk_step,18,mesh%v_index,mesh%nv_halo(1),"k6k",scalar_delp_at_pc_full_level_np1%f)
        end if
#endif

       return
     end subroutine grist_nh_dynamics_run

!================================================
!  BELOW IS PRIVATE
!================================================

     subroutine grist_nh_dynamics_run_explicit(mesh, dtime, itimestep, irk                  , &
                                               scalar_normal_mass_flux_at_edge_face_level_rk, &
                                               tend_mass_at_pc_face_level_hori_rk           , &
                                               scalar_eta_mass_flux_at_pc_full_level_rk     , &
                                               scalar_eta_mass_flux_at_pc_face_level_rk     , &
                                               scalar_phi_at_pc_face_level_rk               , &
                                               scalar_www_at_pc_face_level_n                , &
                                               scalar_www_at_pc_face_level_rk               , &
                                               scalar_delhp_at_pc_face_level_rk             , &
                                               scalar_delhp_at_pc_face_level_n              , &
                                               scalar_delhp_at_pc_face_level_np1            , &
                                               scalar_geopotential_at_pc_surface_n          , &
                                               tend_et_phi_at_pc_face_level_rk              , &
                                               tend_et_www_at_pc_face_level_rk)
! io
        use omp_lib
        type(global_domain),   intent(inout)  :: mesh
        real(r8)           ,   intent(in)     :: dtime
        integer(i4)        ,   intent(in)     :: itimestep, irk
        type(scalar_2d_field), intent(in)     :: scalar_normal_mass_flux_at_edge_face_level_rk
        type(scalar_2d_field), intent(in)     :: tend_mass_at_pc_face_level_hori_rk
        type(scalar_2d_field), intent(in)     :: scalar_eta_mass_flux_at_pc_full_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_eta_mass_flux_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_phi_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_np1
        type(scalar_1d_field), intent(in)     :: scalar_geopotential_at_pc_surface_n
        type(scalar_2d_field), intent(inout)  :: tend_et_phi_at_pc_face_level_rk
        type(scalar_2d_field), intent(inout)  :: tend_et_www_at_pc_face_level_rk
! local
        integer(i4)                           :: iv, ilev, ie, v1, v2, index_edge, ierr
        real(r8)                              :: tmp_sum, www_max, www_max_global
        real(r8)                              :: www_edge_bot(mesh%ne_full)

           if(mpi_rank() == 0.and.write_verbose) print*,"grist_nh_et_phi_face"
        !need halo(1)
! old-code
!           call grist_nh_et_www_face(mesh , &
!                                     scalar_normal_mass_flux_at_edge_face_level_rk , & ! in, DELHP*V, OVERWRITE EACH RK STEP
!                                     tend_mass_at_pc_face_level_hori_rk            , & ! in, DIV.(delhp*v) at face
!                                     scalar_eta_mass_flux_at_pc_full_level_rk      , & ! in, M*ETADOT
!                                     scalar_phi_at_pc_face_level_rk                , & ! in
!                                     scalar_delhp_at_pc_face_level_rk              , & ! in
!                                     tend_et_phi_at_pc_face_level_rk, phi_adv_flag)    ! out, no minus sign
   
   
           if(mpi_rank() == 0.and.write_verbose) print*,"grist_nh_et_www_face"
           !need halo(1)
! old-code
!           call grist_nh_et_www_face(mesh , &
!                                     scalar_normal_mass_flux_at_edge_face_level_rk , & ! in, DELHP*V
!                                     tend_mass_at_pc_face_level_hori_rk            , & ! in,
!                                     scalar_eta_mass_flux_at_pc_full_level_rk      , & ! in, m*etadot
!                                     scalar_www_at_pc_face_level_rk                , & ! in
!                                     scalar_delhp_at_pc_face_level_rk              , & ! in
!                                     tend_et_www_at_pc_face_level_rk, www_adv_flag)                  ! out, no minus sign

           call grist_nh_et_adv_face(mesh, dtime, &
                                     scalar_normal_mass_flux_at_edge_face_level_rk ,& ! DELHP*V
                                     tend_mass_at_pc_face_level_hori_rk            ,& ! div.(mass)
                                     scalar_eta_mass_flux_at_pc_full_level_rk      ,& ! m*etadot
                                     scalar_eta_mass_flux_at_pc_face_level_rk      ,& ! m*etadot
                                     scalar_delhp_at_pc_face_level_rk              ,&
                                     scalar_delhp_at_pc_face_level_n               ,&
                                     scalar_delhp_at_pc_face_level_np1             ,&
                                     scalar_www_at_pc_face_level_rk                ,&
                                     scalar_phi_at_pc_face_level_rk                ,&
                                     tend_et_www_at_pc_face_level_rk               ,&
                                     tend_et_phi_at_pc_face_level_rk               ,&
                                     www_adv_flag(irk))
           
           scalar_www_top(:) = 0._r8
! #ifdef NHDWPIER
!            scalar_www_bot(:) = tend_et_phi_at_pc_face_level_rk%f(nlev+1,:)/gravity-scalar_www_at_pc_face_level_rk%f_r4(nlev+1,:) 
! #elif defined (NHDWPIEN)
!            scalar_www_bot(:) = tend_et_phi_at_pc_face_level_rk%f(nlev+1,:)/gravity-scalar_www_at_pc_face_level_n%f(nlev+1,:) 
! #else
           scalar_www_bot(:) = tend_et_phi_at_pc_face_level_rk%f_r4(nlev+1,:)/gravity
! #endif

           www_max = 0._r8

!$omp parallel  private(iv) 
!$omp do schedule(static,5) reduction(max:www_max)
           do iv = 1, mesh%nv
              www_max = max(www_max,abs(scalar_www_bot(iv)))
           end do
!$omp end do nowait
!$omp end parallel 

#ifndef SEQ_GRIST
           call reduce(www_max, www_max_global, 'max')
#else
           www_max_global = www_max
#endif
           if(mpi_rank() == 0.and.write_verbose) print*,"max www_bot at this step is", www_max_global
           if(mpi_rank() == 0.and.write_verbose) print*,"out of grist_nh_dynamics_run_explicit"

       return
    end subroutine grist_nh_dynamics_run_explicit

    subroutine grist_nh_dynamics_run_implicit(mesh, dtime,                               &
                                              scalar_pressure_at_pc_full_level_n,        &
                                              scalar_phi_at_pc_face_level_n,             &
                                              scalar_delhp_at_pc_full_level_n,           &
                                              scalar_www_at_pc_face_level_n,             &
                                              scalar_potential_temp_at_pc_full_level_n,  &
                                              scalar_eta_mass_flux_at_pc_face_level_rk,  &
                                              scalar_www_at_pc_face_level_rk,            &
                                              tend_et_www_at_pc_face_level_rk,           &
                                              tend_et_phi_at_pc_face_level_rk,           &
                                              scalar_potential_temp_at_pc_full_level_np1,&
                                              scalar_delhp_at_pc_full_level_np1,         &
                                              scalar_delhp_at_pc_face_level_np1,         &
                                              scalar_delp_at_pc_face_level_rk,           &
                                              scalar_delhp_at_pc_face_level_rk,          &
                                              scalar_mif_at_pc_face_level_n,             &
                                              scalar_mif_at_pc_full_level_n,             &
                                              scalar_www_at_pc_face_level_np1,           &
                                              scalar_phi_at_pc_face_level_np1)
                                              !scalar_kesi_at_pc_full_level_np1,          &
! io
        use omp_lib
        type(global_domain),   intent(in)     :: mesh
        real(r8)           ,   intent(in)     :: dtime
        type(scalar_2d_field), intent(in)     :: scalar_pressure_at_pc_full_level_n
        type(scalar_2d_field), intent(in)     :: scalar_phi_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_full_level_n
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_potential_temp_at_pc_full_level_n
        type(scalar_2d_field), intent(in)     :: scalar_eta_mass_flux_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: tend_et_www_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: tend_et_phi_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_potential_temp_at_pc_full_level_np1
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_full_level_np1
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_np1
        type(scalar_2d_field), intent(in)     :: scalar_delp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_mif_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_mif_at_pc_full_level_n
        type(scalar_2d_field), intent(inout)  :: scalar_www_at_pc_face_level_np1
        type(scalar_2d_field), intent(inout)  :: scalar_phi_at_pc_face_level_np1
        !type(scalar_2d_field), intent(inout)  :: scalar_kesi_at_pc_full_level_np1
! local
        !type(scalar_2d_field)                 :: scalar_kesi_at_pc_face_level_np1
        integer(i4)                           :: iv, ilev
        real(r4)                              :: dtime_r4
        ! real(r8)                              :: tmp_coef,cr_num,cr_act
        real(r8)                              :: tmp_coef
        real(r8)                              :: cr_num,cr_act
        integer                               :: ii

        dtime_r4 = dtime
          ! scalar_kesi_at_pc_face_level_np1 = scalar_www_at_pc_face_level_np1 ! init

           call grist_nh_solve_www_equation(mesh, nlev, nlevp, &
                                            scalar_pressure_at_pc_full_level_n%f,        &
                                            scalar_delhp_at_pc_full_level_n%f,           &
                                            scalar_potential_temp_at_pc_full_level_n%f,  &
                                            scalar_potential_temp_at_pc_full_level_np1%f,&
                                            scalar_delhp_at_pc_full_level_np1%f,         &
                                            scalar_phi_at_pc_face_level_n%f,             &
                                            scalar_www_at_pc_face_level_n%f,             &
                                            scalar_www_at_pc_face_level_rk%f_r4,            &
                                            tend_et_www_at_pc_face_level_rk%f_r4,           &
                                            tend_et_phi_at_pc_face_level_rk%f_r4,           &
                                            scalar_delhp_at_pc_face_level_np1%f,         &
                                            scalar_delp_at_pc_face_level_rk%f,           &
                                            scalar_delhp_at_pc_face_level_rk%f,          &
                                            scalar_www_at_pc_face_level_np1%f,           &
                                            scalar_mif_at_pc_face_level_n%f,             &
                                            scalar_mif_at_pc_full_level_n%f,             &
                                            dtime,                               &
                                            scalar_www_top, scalar_www_bot)
! #ifdef NHDWPIER
!            scalar_www_at_pc_face_level_np1%f =scalar_www_at_pc_face_level_np1%f+scalar_www_at_pc_face_level_rk%f
! #endif
! #ifdef NHDWPIEN
!            scalar_www_at_pc_face_level_np1%f =scalar_www_at_pc_face_level_np1%f+scalar_www_at_pc_face_level_n%f
! #endif
!
! optional w-damping as in KDH08 for upper levels
!
!$omp parallel  private(ii,iv,ilev,tmp_coef) 
!$omp do schedule(dynamic,100) 
     do ii = 1, mesh%nv_halo(1)*(nlev+1),1
        iv=ceiling(ii/real((nlev+1),r4))
        ilev=ii-(iv-1)*(nlev+1)

!           do iv = 1, mesh%nv_halo(1)
!              do ilev = 1, nlev+1
                 if(zd.gt.0._r8.and.www_damping_coef.gt.0._r8.and.&
                    scalar_phi_at_pc_face_level_n%f(ilev,iv).ge.(scalar_phi_at_pc_face_level_n%f(1,iv)-zd*gravity))then
                    tmp_coef = www_damping_coef*(sin((pi*0.5_r8)*(1._r8-(scalar_phi_at_pc_face_level_n%f(1,iv)-scalar_phi_at_pc_face_level_n%f(ilev,iv))/(zd*gravity)))**2)
                    scalar_www_at_pc_face_level_np1%f(ilev,iv) = scalar_www_at_pc_face_level_np1%f(ilev,iv)/(1._r4+dtime*tmp_coef)
                 end if
!              end do
!           end do
     end do
!$omp end do nowait
!$omp end parallel 
!
!  do a w-damping in case of violation of CFL
!  if this is activated, in most cases, the model will blow
!
!$omp parallel  private(ii,iv,ilev,cr_num,cr_act)  
!$omp do schedule(dynamic,100) 
     do ii = 1,mesh%nv_halo(1)*(nlev-1),1
        iv=ceiling(ii/real((nlev-1),r8))
        ilev=1+ii-(iv-1)*(nlev-1)

!           do iv = 1, mesh%nv_halo(1)
!              do ilev = 2, nlev
#ifdef DBL_MASS
                 cr_num = scalar_eta_mass_flux_at_pc_face_level_rk%f(ilev,iv)*dtime/scalar_delhp_at_pc_face_level_rk%f(ilev,iv)
#else
                 cr_num = scalar_eta_mass_flux_at_pc_face_level_rk%f_r4(ilev,iv)*dtime/scalar_delhp_at_pc_face_level_rk%f(ilev,iv)
#endif
                 cr_act = max(0._r8,abs(cr_num)-1._r8)
                 if(cr_act.gt.0._r8) print*,"violation of vertical cfl at, ilev=",ilev
                 cr_act = cr_act*sign(1._r8,scalar_www_at_pc_face_level_rk%f_r4(ilev,iv))
                 scalar_www_at_pc_face_level_np1%f(ilev,iv) = scalar_www_at_pc_face_level_np1%f(ilev,iv)-0.3_r8*cr_act*dtime
!              end do
!           end do
     end do
!$omp end do nowait
!$omp end parallel 
!
! diagnose kesi
!
        !   scalar_kesi_at_pc_face_level_np1%f     = (scalar_www_at_pc_face_level_np1%f-scalar_www_at_pc_face_level_n%f)/dtime+&
        !                                             tend_et_www_at_pc_face_level_rk%f+gravity-gravity*(1._r8-imbeta)*&
        !                                             scalar_delp_at_pc_face_level_rk%f/scalar_delhp_at_pc_face_level_rk%f
        !   scalar_kesi_at_pc_face_level_np1%f     = scalar_kesi_at_pc_face_level_np1%f/(imbeta*gravity)
        !   scalar_kesi_at_pc_face_level_np1%f(1,:)= 1._r8
        !   do ilev = 1, nlev
        !      scalar_kesi_at_pc_full_level_np1%f(ilev,:) = 0.5_r8*(scalar_kesi_at_pc_face_level_np1%f(ilev,:)+&
        !                                                           scalar_kesi_at_pc_face_level_np1%f(ilev+1,:))
        !   end do
        !   
        !   scalar_kesi_at_pc_full_level_np1%f(nlev,:) = 1._r8
!
! update phi equation, not updated one time due to an intel bug
!
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_n%f(1:nlev,:)
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_np1%f(1:nlev,:)-dtime*tend_et_phi_at_pc_face_level_rk%f_r4(1:nlev,:)
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_np1%f(1:nlev,:)+&
                                                         (1._r8-imbeta)*gravity*dtime*scalar_www_at_pc_face_level_rk%f_r4(1:nlev,:)
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_np1%f(1:nlev,:)+&
                                                                 imbeta*gravity*dtime*scalar_www_at_pc_face_level_np1%f(1:nlev,:)
!
! enforce dphi to be positve
!
!           do iv = 1, mesh%nv_halo(1)
!              do ilev = nlev, 1
!                 scalar_phi_at_pc_face_level_np1%f(ilev,iv) = max(scalar_phi_at_pc_face_level_np1%f(ilev,iv),&
!                                                                  scalar_phi_at_pc_face_level_np1%f(ilev+1,iv)+20._r8*gravity)
!                if((scalar_phi_at_pc_face_level_np1%f(ilev,iv)-scalar_phi_at_pc_face_level_np1%f(ilev+1,iv)).lt.20._r8*gravity)then
!                  print*,"dz<0"
!                  scalar_phi_at_pc_face_level_np1%f(ilev,iv) = scalar_phi_at_pc_face_level_np1%f(ilev+1,iv)+20._r8*gravity
!                end if
           !     if((scalar_phi_at_pc_face_level_np1%f(ilev,iv)-scalar_phi_at_pc_face_level_np1%f(ilev+1,iv)).lt.2._r8*gravity)then
           !       print*,"dz<0, again?"
           !     end if
!              end do
!           end do

       return
     end subroutine grist_nh_dynamics_run_implicit

     subroutine grist_nh_dynamics_run_diagvars(mesh, itimestep, irk_step, idstep, itstep, &
                                                     dtime,    &
                                                     scalar_delp_at_pc_face_level_rk,         &
                                                     scalar_delhp_at_pc_face_level_rk,        &
                                                     scalar_hpressure_at_pc_face_level_np1,   &
                                                     scalar_delhp_at_pc_face_level_np1,       &
                                                     scalar_pressure_at_pc_full_level_rk,     &
                                                     scalar_pressure_at_pc_full_level_n,      &
                                                     scalar_phi_at_pc_face_level_n,           &
                                                     scalar_phi_at_pc_face_level_np1,         &
                                                     scalar_www_at_pc_face_level_n,           &
                                                     scalar_www_at_pc_face_level_np1,         &
                                                     tend_et_www_at_pc_face_level_rk,         &
                                                     scalar_delhp_at_pc_full_level_n,         &
                                                     scalar_delhp_at_pc_full_level_np1,       &
                                                     scalar_potential_temp_at_pc_full_level_n,&
                                                     scalar_potential_temp_at_pc_full_level_np1,&
                                                     scalar_geopotential_at_pc_surface_n,     &
                                                     scalar_mif_at_pc_face_level_n,           &
                                                     scalar_mif_at_pc_full_level_n,           &
                                                     scalar_tracer_mxrt_at_pc_full_level_n,   &
                                                     scalar_temp_at_pc_full_level_np1,        &
                                                     scalar_geopotential_at_pc_full_level_np1,&
                                                     scalar_geopotential_at_pc_face_level_np1,&
                                                     scalar_pressure_at_pc_full_level_np1    ,&
                                                     scalar_pressure_at_pc_face_level_np1    ,&
                                                     scalar_delp_at_pc_full_level_np1,        &
                                                     !scalar_kesi_at_pc_full_level_np1,        &
                                                     scalar_alpha_at_pc_full_level_np1)
! io
        use omp_lib
        type(global_domain),   intent(in)     :: mesh
        integer(i4),           intent(in)     :: itimestep, irk_step, idstep, itstep 
        real(r8),              intent(in)     :: dtime
        type(scalar_2d_field), intent(in)     :: scalar_delp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_hpressure_at_pc_face_level_np1
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_np1 
        type(scalar_2d_field), intent(in)     :: scalar_pressure_at_pc_full_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_pressure_at_pc_full_level_n
        type(scalar_2d_field), intent(in)     :: scalar_phi_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_phi_at_pc_face_level_np1
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_np1
        type(scalar_2d_field), intent(in)     :: tend_et_www_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_full_level_n
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_full_level_np1
        type(scalar_2d_field), intent(in)     :: scalar_potential_temp_at_pc_full_level_n
        type(scalar_2d_field), intent(in)     :: scalar_potential_temp_at_pc_full_level_np1
        type(scalar_1d_field), intent(in)     :: scalar_geopotential_at_pc_surface_n
        type(scalar_2d_field), intent(in)     :: scalar_mif_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_mif_at_pc_full_level_n
        type(scalar_3d_field), intent(in)     :: scalar_tracer_mxrt_at_pc_full_level_n
        type(scalar_2d_field), intent(inout)  :: scalar_temp_at_pc_full_level_np1
        type(scalar_2d_field), intent(inout)  :: scalar_geopotential_at_pc_full_level_np1
        type(scalar_2d_field), intent(inout)  :: scalar_geopotential_at_pc_face_level_np1
        type(scalar_2d_field), intent(inout)  :: scalar_pressure_at_pc_full_level_np1
        type(scalar_2d_field), intent(inout)  :: scalar_pressure_at_pc_face_level_np1
        type(scalar_2d_field), intent(inout)  :: scalar_delp_at_pc_full_level_np1
       ! type(scalar_2d_field), intent(inout)  :: scalar_kesi_at_pc_full_level_np1
        type(scalar_2d_field), intent(inout)  :: scalar_alpha_at_pc_full_level_np1
! local
       ! type(scalar_2d_field)  :: scalar_kesi_at_pc_face_level_np1
        real(r4)       :: dtime_r4
        real(r4)       :: sound_speed,  sound_speed_max,  sound_speed_max_global
        real(r4)       :: sound_speed1, sound_speed1_max, sound_speed1_max_global
        real(r4)       :: max_kesi,min_kesi, tmp1(mesh%nv_full), tmp2(mesh%nv_full)
        integer(i4)    :: ie,ilev, iv, v1, v2, ierr
        integer(i4)    :: ii

        dtime_r4 = dtime
       ! scalar_kesi_at_pc_face_level_np1 = scalar_phi_at_pc_face_level_n
!
! compute density with updated phi and delhp, then compute new pressure and temp
!
!
! par-seq consistenty under intel-O2 sensitive to loop implementation on some machines
! default not use this imploop, leading to nan at cma-pi for par, but not seq
!

! #ifdef NHDRV_IMPLOOP
!            do ilev = 1, nlev

!               alpha_np1(:) = scalar_phi_at_pc_face_level_np1%f(ilev+1,:)-scalar_phi_at_pc_face_level_np1%f(ilev,:)
!               alpha_n(:)   = scalar_phi_at_pc_face_level_n%f(ilev+1,:)  -scalar_phi_at_pc_face_level_n%f(ilev,:)
!               scalar_alpha_at_pc_full_level_np1%f(ilev,:) = -alpha_np1(:)/scalar_delhp_at_pc_full_level_np1%f(ilev,:) ! ensure >0
! #if !defined NEWLGE && !defined LGE1
!               scalar_pressure_at_pc_full_level_np1%f(ilev,:) = scalar_pressure_at_pc_full_level_n%f(ilev,:)+&
!                                                                (cp/cv)*scalar_pressure_at_pc_full_level_n%f(ilev,:)*&
!                                                               ((scalar_potential_temp_at_pc_full_level_np1%f(ilev,:)/&
!                                                                 scalar_potential_temp_at_pc_full_level_n%f(ilev,:))-&
!                                                                 alpha_np1(:)/alpha_n(:))
! #endif
! #ifdef LGE1
!               alpha_np1(:) = (scalar_phi_at_pc_face_level_np1%f(ilev+1,:)-scalar_phi_at_pc_face_level_np1%f(ilev,:))/&
!                               scalar_delhp_at_pc_full_level_np1%f(ilev,:)
!               alpha_n(:)   = (scalar_phi_at_pc_face_level_n%f(ilev+1,:)  -scalar_phi_at_pc_face_level_n%f(ilev,:))/&
!                               scalar_delhp_at_pc_full_level_n%f(ilev,:)
!               scalar_pressure_at_pc_full_level_np1%f(ilev,:) = scalar_pressure_at_pc_full_level_n%f(ilev,:)+&
!                                                                (cp/cv)*scalar_pressure_at_pc_full_level_n%f(ilev,:)*&
!                                                               ((scalar_potential_temp_at_pc_full_level_np1%f(ilev,:)/&
!                                                                 scalar_potential_temp_at_pc_full_level_n%f(ilev,:))-&
!                                                                 alpha_np1(:)/alpha_n(:)-&
!                                                                 scalar_delhp_at_pc_full_level_np1%f(ilev,:)/scalar_delhp_at_pc_full_level_n%f(ilev,:)+one)
! #endif

! ! do a 3d div damping
!               scalar_pressure_at_pc_full_level_np1%f(ilev,:) = scalar_pressure_at_pc_full_level_np1%f(ilev,:)+&
!                                                                d3d_damping_coef*&
!                                                               (scalar_pressure_at_pc_full_level_np1%f(ilev,:)-&
!                                                                scalar_pressure_at_pc_full_level_n%f(ilev,:))
!                                                                !scalar_pressure_at_pc_full_level_rk%f(ilev,:)) ! testing

!               scalar_temp_at_pc_full_level_np1%f(ilev,:) = scalar_potential_temp_at_pc_full_level_np1%f(ilev,:)/scalar_delhp_at_pc_full_level_np1%f(ilev,:)*&
!                                                           ((scalar_pressure_at_pc_full_level_np1%f(ilev,:)/p00)**(rdry/cp))/&
!                                                           (one+ptfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,ilev,:))
            
!           end do
! #else
!$omp parallel  private(iv,ilev)   
!$omp do schedule(dynamic,20) 
          do iv = 1, mesh%nv_halo(1)
             do ilev = 1, nlev

              alpha_np1(iv) = scalar_phi_at_pc_face_level_np1%f(ilev+1,iv)-scalar_phi_at_pc_face_level_np1%f(ilev,iv)
              alpha_n(iv)   = scalar_phi_at_pc_face_level_n%f(ilev+1,iv)  -scalar_phi_at_pc_face_level_n%f(ilev,iv)
              scalar_alpha_at_pc_full_level_np1%f(ilev,iv) = -alpha_np1(iv)/scalar_delhp_at_pc_full_level_np1%f(ilev,iv) ! ensure >0
#if !defined NEWLGE && !defined LGE1
              scalar_pressure_at_pc_full_level_np1%f(ilev,iv) = scalar_pressure_at_pc_full_level_n%f(ilev,iv)+&
                                                               (cp/cv)*scalar_pressure_at_pc_full_level_n%f(ilev,iv)*&
                                                              ((scalar_potential_temp_at_pc_full_level_np1%f(ilev,iv)/&
                                                                scalar_potential_temp_at_pc_full_level_n%f(ilev,iv))-&
                                                                alpha_np1(iv)/alpha_n(iv))
#endif
! #ifdef LGE1
!               alpha_np1(iv) = (scalar_phi_at_pc_face_level_np1%f(ilev+1,iv)-scalar_phi_at_pc_face_level_np1%f(ilev,iv))/&
!                               scalar_delhp_at_pc_full_level_np1%f(ilev,iv)
!               alpha_n(iv)   = (scalar_phi_at_pc_face_level_n%f(ilev+1,iv)  -scalar_phi_at_pc_face_level_n%f(ilev,iv))/&
!                               scalar_delhp_at_pc_full_level_n%f(ilev,iv)
!               scalar_pressure_at_pc_full_level_np1%f(ilev,iv) = scalar_pressure_at_pc_full_level_n%f(ilev,iv)+&
!                                                                (cp/cv)*scalar_pressure_at_pc_full_level_n%f(ilev,iv)*&
!                                                               ((scalar_potential_temp_at_pc_full_level_np1%f(ilev,iv)/&
!                                                                 scalar_potential_temp_at_pc_full_level_n%f(ilev,iv))-&
!                                                                 alpha_np1(iv)/alpha_n(iv)-&
!                                                                 scalar_delhp_at_pc_full_level_np1%f(ilev,iv)/scalar_delhp_at_pc_full_level_n%f(ilev,iv)+one)
! #endif

! do a 3d div damping
              scalar_pressure_at_pc_full_level_np1%f(ilev,iv) = scalar_pressure_at_pc_full_level_np1%f(ilev,iv)+&
                                                               d3d_damping_coef*&
                                                              (scalar_pressure_at_pc_full_level_np1%f(ilev,iv)-&
                                                               scalar_pressure_at_pc_full_level_n%f(ilev,iv))
                                                               !scalar_pressure_at_pc_full_level_rk%f(ilev,iv)) ! testing

              scalar_temp_at_pc_full_level_np1%f(ilev,iv) = scalar_potential_temp_at_pc_full_level_np1%f(ilev,iv)/scalar_delhp_at_pc_full_level_np1%f(ilev,iv)*&
                                                          ((scalar_pressure_at_pc_full_level_np1%f(ilev,iv)/p00)**(rdry/cp))/&
                                                          (one+ptfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,ilev,iv))
             end do
           end do
!$omp end do nowait
!$omp end parallel 

! #endif

           sound_speed_max  = 0._r4
           sound_speed1_max = 0._r4
!$omp parallel  private(iv,ii,ilev,sound_speed,sound_speed1)  
!$omp do schedule(dynamic,100) reduction(max:sound_speed_max,sound_speed1_max)
     do ii = 1, mesh%nv_compute*nlev,1
        iv=ceiling(ii/real(nlev,r8))
        ilev=ii-(iv-1)*nlev
!           do iv = 1, mesh%nv_compute
!              do ilev = 1, nlev
                 sound_speed  = (cp/cv)*scalar_pressure_at_pc_full_level_np1%f(ilev,iv)*&
                                       (scalar_phi_at_pc_face_level_np1%f(ilev,iv)-scalar_phi_at_pc_face_level_np1%f(ilev+1,iv))/scalar_delhp_at_pc_full_level_np1%f(ilev,iv)
                 sound_speed  = sound_speed/(one+ptfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,ilev,iv))
                 sound_speed1 = (cp/cv)*rdry*scalar_temp_at_pc_full_level_np1%f(ilev,iv)
                 sound_speed_max  = max(sound_speed_max ,sqrt(sound_speed))
                 sound_speed1_max = max(sound_speed1_max,sqrt(sound_speed1))
!              end do
!           end do
     end do
!$omp end do nowait
!$omp end parallel

#ifndef SEQ_GRIST
           call reduce(sound_speed_max, sound_speed_max_global, 'max')
           call reduce(sound_speed1_max, sound_speed1_max_global, 'max')
! #else
!            sound_speed_max_global = sound_speed_max
!            sound_speed1_max_global= sound_speed1_max
#endif
           if(mpi_rank() == 0.and.write_stepinfo) print*,"max sound speed at this step is", sound_speed_max_global, sound_speed1_max_global
!
! this is to dynamic control
!
           if(mpi_rank().eq.0)then
              if(abs(sound_speed_max_global-sound_speed1_max_global).gt.ndc_restore_num)then
                 if(irk_step.eq.1.and.idstep.eq.1.and.itstep.eq.1.and.ndc_restore_flag.ne.1)then
                    ndc_restore_flag = 1
                    itimestep_old    = itimestep
                 end if
              end if
              if(itimestep.ne.itimestep_old) ndc_restore_flag = 0
           end if
#ifndef SEQ_GRIST
           CALL mpi_bcast(ndc_restore_flag,1,MPI_INTEGER,0,grist_comm(),ierr)
           CALL mpi_bcast(itimestep_old   ,1,MPI_INTEGER,0,grist_comm(),ierr)
#endif

! #ifdef NHDRV_IMPLOOP2
! !
! ! compute geopotential at face and full levels, just use another name
! !
!            scalar_geopotential_at_pc_face_level_np1%f  = scalar_phi_at_pc_face_level_np1%f
!            do ilev = 1, nlev
!               scalar_geopotential_at_pc_full_level_np1%f(ilev,:) = 0.5_r8*(scalar_geopotential_at_pc_face_level_np1%f(ilev+1,:)+&
!                                                                            scalar_geopotential_at_pc_face_level_np1%f(ilev,:))
!            end do
! ! face level pressure, as in hdc
!            scalar_pressure_at_pc_face_level_np1%f(1,:) = scalar_hpressure_at_pc_face_level_np1%f(1,:)!*&
!                                         !                 (one+prfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,1,:))

!            do ilev = 2, nlev+1
! ! #ifdef FACEP1
! !               scalar_pressure_at_pc_face_level_np1%f(ilev,:) = scalar_pressure_at_pc_face_level_np1%f(ilev-1,:)+scalar_delhp_at_pc_full_level_np1%f(ilev-1,:)*&
! !                                                                (0.5_r8/gravity*&
! !                                                                (((scalar_www_at_pc_face_level_np1%f(ilev-1,:)-scalar_www_at_pc_face_level_n%f(ilev-1,:))/dtime+&
! !                                                                   tend_et_www_at_pc_face_level_rk%f(ilev-1,:))+&
! !                                                                 ((scalar_www_at_pc_face_level_np1%f(ilev,:)-scalar_www_at_pc_face_level_n%f(ilev,:))/dtime+&
! !                                                                   tend_et_www_at_pc_face_level_rk%f(ilev,:)))+1._r8)/scalar_mif_at_pc_full_level_n%f(ilev-1,:)
! ! #else 
!               tmp1(:) =(((scalar_www_at_pc_face_level_np1%f(ilev-1,:)-scalar_www_at_pc_face_level_n%f(ilev-1,:))/dtime+tend_et_www_at_pc_face_level_rk%f(ilev-1,:))/gravity+one)/&
!                        scalar_mif_at_pc_face_level_n%f(ilev-1,:)
!               tmp2(:) =(((scalar_www_at_pc_face_level_np1%f(ilev,:)-scalar_www_at_pc_face_level_n%f(ilev,:))/dtime+tend_et_www_at_pc_face_level_rk%f(ilev,:))/gravity+one)/&
!                        scalar_mif_at_pc_face_level_n%f(ilev,:)
! ! #ifdef FBDP
! !               tmp1(:) = (tmp1(:)-(one-imbeta)*scalar_delp_at_pc_face_level_rk%f(ilev-1,:)/scalar_delhp_at_pc_face_level_rk%f(ilev-1,:))/imbeta
! !               tmp2(:) = (tmp2(:)-(one-imbeta)*scalar_delp_at_pc_face_level_rk%f(ilev,:)  /scalar_delhp_at_pc_face_level_rk%f(ilev  ,:))/imbeta

! ! #endif
!               scalar_pressure_at_pc_face_level_np1%f(ilev,:) = scalar_pressure_at_pc_face_level_np1%f(ilev-1,:)+scalar_delhp_at_pc_full_level_np1%f(ilev-1,:)*&
!                                                                (tmp1(:)+tmp2(:))*0.5_r8
! ! #endif
!            end do

! ! #ifdef REGRESSION
! !            do ilev = 2, nlev
! !               scalar_pressure_at_pc_face_level_np1%f(ilev,:) = 0.5_r8*(deta_full(ilev-1)/deta_face(ilev)*scalar_pressure_at_pc_full_level_np1%f(ilev,:)+&
! !                                                                        deta_full(ilev)  /deta_face(ilev)*scalar_pressure_at_pc_full_level_np1%f(ilev-1,:))
! !            end do
! !            scalar_pressure_at_pc_face_level_np1%f(1,:)      = 2._r8*scalar_pressure_at_pc_full_level_np1%f(1,:)-scalar_pressure_at_pc_face_level_np1%f(2,:)
! !            scalar_pressure_at_pc_face_level_np1%f(nlev+1,:) = 2._r8*scalar_pressure_at_pc_full_level_np1%f(nlev,:)-scalar_pressure_at_pc_face_level_np1%f(nlev,:)
! ! #endif
! !
! ! compute delp at full level
! !
!         do ilev = 1, nlev
!            scalar_delp_at_pc_full_level_np1%f(ilev,:) = scalar_pressure_at_pc_face_level_np1%f(ilev+1,:)-&
!                                                         scalar_pressure_at_pc_face_level_np1%f(ilev,:)
!         end do

! ! NHDRV_IMPLOOP2
! #else
        
!
! compute geopotential at face and full levels, just use another name
!
        scalar_geopotential_at_pc_face_level_np1%f  = scalar_phi_at_pc_face_level_np1%f  ! sensitive

!$omp parallel  private(iv,ilev)  
!$omp do schedule(dynamic,20) 
        do iv = 1, mesh%nv_halo(1)
           do ilev = 1, nlev
              scalar_geopotential_at_pc_full_level_np1%f(ilev,iv) = 0.5_r8*(scalar_geopotential_at_pc_face_level_np1%f(ilev+1,iv)+&
                                                                            scalar_geopotential_at_pc_face_level_np1%f(ilev,iv))
           end do
! face level pressure, as in hdc
           scalar_pressure_at_pc_face_level_np1%f_r4(1,iv) = scalar_hpressure_at_pc_face_level_np1%f_r4(1,iv) !*&
                                      !                 (one+prfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,1,:))
        end do
!$omp end do nowait
!$omp end parallel 

!$omp parallel default(shared) private(iv,ilev)  
!$omp do schedule(dynamic,20) 
        do iv = 1, mesh%nv_halo(1)
           do ilev = 2, nlev+1
! #ifdef FACEP1
!               scalar_pressure_at_pc_face_level_np1%f(ilev,iv) = scalar_pressure_at_pc_face_level_np1%f(ilev-1,iv)+scalar_delhp_at_pc_full_level_np1%f(ilev-1,iv)*&
!                                                                (0.5_r8/gravity*&
!                                                                (((scalar_www_at_pc_face_level_np1%f(ilev-1,iv)-scalar_www_at_pc_face_level_n%f(ilev-1,iv))/dtime+&
!                                                                   tend_et_www_at_pc_face_level_rk%f(ilev-1,iv))+&
!                                                                 ((scalar_www_at_pc_face_level_np1%f(ilev,iv)-scalar_www_at_pc_face_level_n%f(ilev,iv))/dtime+&
!                                                                   tend_et_www_at_pc_face_level_rk%f(ilev,iv)))+1._r8)/scalar_mif_at_pc_full_level_n%f(ilev-1,iv)
! #else 
              tmp1(iv) =(((scalar_www_at_pc_face_level_np1%f_r4(ilev-1,iv)-scalar_www_at_pc_face_level_n%f_r4(ilev-1,iv))/dtime_r4+tend_et_www_at_pc_face_level_rk%f_r4(ilev-1,iv))/gravity_r4+one_r4)/&
                           scalar_mif_at_pc_face_level_n%f_r4(ilev-1,iv)
              tmp2(iv) =(((scalar_www_at_pc_face_level_np1%f_r4(ilev,iv)-scalar_www_at_pc_face_level_n%f_r4(ilev,iv))/dtime_r4+tend_et_www_at_pc_face_level_rk%f_r4(ilev,iv))/gravity_r4+one_r4)/&
                           scalar_mif_at_pc_face_level_n%f_r4(ilev,iv)
! #ifdef FBDP
!               tmp1(iv) = (tmp1(iv)-(one-imbeta)*scalar_delp_at_pc_face_level_rk%f(ilev-1,iv)/scalar_delhp_at_pc_face_level_rk%f(ilev-1,iv))/imbeta
!               tmp2(iv) = (tmp2(iv)-(one-imbeta)*scalar_delp_at_pc_face_level_rk%f(ilev,iv)  /scalar_delhp_at_pc_face_level_rk%f(ilev  ,iv))/imbeta

! #endif
              scalar_pressure_at_pc_face_level_np1%f_r4(ilev,iv) = scalar_pressure_at_pc_face_level_np1%f_r4(ilev-1,iv)+scalar_delhp_at_pc_full_level_np1%f_r4(ilev-1,iv)*&
                                                               (tmp1(iv)+tmp2(iv))*0.5_r4
! #endif
           end do
        end do
!$omp end do nowait
!$omp end parallel 
        ! scalar_pressure_at_pc_face_level_np1%f = scalar_pressure_at_pc_face_level_np1%f_r4
! #ifdef REGRESSION
! !$omp parallel  private(iv,ilev) 
! !$omp do schedule(dynamic,20) 
!         do iv = 1, mesh%nv_halo(1)
!            do ilev = 2, nlev
!               scalar_pressure_at_pc_face_level_np1%f(ilev,iv) = 0.5_r8*(deta_full(ilev-1)/deta_face(ilev)*scalar_pressure_at_pc_full_level_np1%f(ilev  ,iv)+&
!                                                                         deta_full(ilev)  /deta_face(ilev)*scalar_pressure_at_pc_full_level_np1%f(ilev-1,iv))
!            end do
!            scalar_pressure_at_pc_face_level_np1%f(1,iv)      = 2._r8*scalar_pressure_at_pc_full_level_np1%f(1   ,iv)-scalar_pressure_at_pc_face_level_np1%f(2,iv)
!            scalar_pressure_at_pc_face_level_np1%f(nlev+1,iv) = 2._r8*scalar_pressure_at_pc_full_level_np1%f(nlev,iv)-scalar_pressure_at_pc_face_level_np1%f(nlev,iv)
!         end do
! !$omp end do nowait
! !$omp end parallel 
! #endif

!
! compute delp at full level
!
!$omp parallel  private(iv,ii,ilev)  
!$omp do schedule(dynamic,100) 
     do ii = 1, mesh%nv_halo(1)*nlev,1
        iv=ceiling(ii/real(nlev,r8))
        ilev=ii-(iv-1)*nlev
!        do iv = 1, mesh%nv_halo(1)
!           do ilev = 1, nlev
              scalar_delp_at_pc_full_level_np1%f_r4(ilev,iv) = scalar_pressure_at_pc_face_level_np1%f_r4(ilev+1,iv)-&
                                                            scalar_pressure_at_pc_face_level_np1%f_r4(ilev  ,iv)
!           end do
!        end do
     end do
!$omp end do nowait
!$omp end parallel 

! #endif

        return
     end subroutine grist_nh_dynamics_run_diagvars

     subroutine grist_nh_dynamics_init(mesh)
! io
        type(global_domain),   intent(in)     :: mesh
! hori
        if(.not.allocated(scalar_www_top))    allocate(scalar_www_top(mesh%nv))
        if(.not.allocated(scalar_www_bot))    allocate(scalar_www_bot(mesh%nv))
        if(.not.allocated(alpha_n))           allocate(alpha_n(mesh%nv))
        if(.not.allocated(alpha_np1))         allocate(alpha_np1(mesh%nv))

        ndc_restore_flag= 0
        itimestep_old   = 1

#ifndef SEQ_GRIST
        clock_nhd_exp     = clock_id('nhd_exp')
        clock_nhd_imp     = clock_id('nhd_imp')
        clock_nhd_diag    = clock_id('nhd_diag')
#endif

        return
     end subroutine grist_nh_dynamics_init

     subroutine grist_nh_dynamics_final

        deallocate(scalar_www_top)
        deallocate(scalar_www_bot)
        deallocate(alpha_n)
        deallocate(alpha_np1)

        return
     end subroutine grist_nh_dynamics_final

   end module grist_nh_driver_module_mixed
