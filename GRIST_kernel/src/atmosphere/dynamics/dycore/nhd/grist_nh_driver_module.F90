
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: The main driver module for the nh-dynamics; The explicit and
!              implicit part are seperated, and only the explicit part needs
!              horizontal data structure (i.e., parallel).
! Revision history:
!              1. add exploop and imploop options for regression usage
!              2. default use exploop
!              3. Add cell-based data exchange pattern for NHD to skip
!                 edge-flux based data exchange pattern. With this, dycore
!                 has no edge-flux data exchange now. For global models, the code
!                 still EXACTLY regress before. For LAM, some changes for bdy-halo update
!                 sequence, and phi's computing to nv_full to ensure self-par consistency
!                 and correct restart. To EXACTLY regress old LAM solutions, we need
!                 -DOLD_MASS_PT & -DNHD_NOOPT (see dycore_time_integration for our
!                 principle for code refactoring).
!----------------------------------------------------------------------------

   module grist_nh_driver_module
! constant
    use grist_constants,               only: i4, r8, fillvalue, gravity, p00, rdry, rvap, ptfactor, prfactor, one, cp, cv, pi, zero
    use grist_nml_module,              only: phi_adv_flag, www_adv_flag, d3d_damping_coef, www_damping_coef,imbeta,zd, write_verbose, write_stepinfo
! data definition
    use grist_data_types,              only: scalar_1d_field, scalar_2d_field, scalar_3d_field
    use grist_domain_types,            only: global_domain
    use grist_math_module,             only: extrapolate_bdy
! hori operator
! nh module
    use grist_nh_explicit_tend_module_2d, only: grist_nh_et_adv_face
    use grist_nh_implicit_module_2d,      only: grist_nh_solve_www_equation
! hpe
    use grist_hpe_constants,              only: deta_full, deta_face,eta_full,eta_face
! mpi
#ifndef SEQ_GRIST
    use grist_mpi
#endif
    use grist_lib
#ifndef SEQ_GRIST
    use grist_config_partition,           only: debug_data_2d,debug_data_1d
    use grist_config_partition,           only: exchange_data_2d_add , exchange_data_2d
    use grist_data_types,                 only: exchange_field_list_2d
#else
    use debug_module,                     only: debug_data_2d,debug_data_1d
#endif
#ifdef LAM_DOMAIN
    use grist_datam_glam_data_module,     only: update_assignValueArea_2d
#endif

    implicit none

      private
      public :: grist_nh_dynamics_init, &
                grist_nh_dynamics_final,&
                grist_nh_dynamics_run
!
! local 1d, hori
!
        real(r8), allocatable       :: scalar_www_top(:)
        real(r8), allocatable       :: scalar_www_bot(:)
        real(r8), allocatable       :: alpha_n(:)
        real(r8), allocatable       :: alpha_np1(:)         
        real(r8), parameter, public :: ndc_restore_num = 10._r8
        integer(i4), public         :: ndc_restore_flag
        integer(i4)                 :: itimestep_old

  CONTAINS
!
! main interface to time integration
!
     subroutine grist_nh_dynamics_run(mesh,nlev,nlevp,dtime, itimestep, irk_step, idstep, itstep, &
                                      tend_hpressure_at_pc_surface_cnty,             &
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
        integer(i4),           intent(in)    :: nlev, nlevp
        real(r8)           ,   intent(in)    :: dtime
        integer(i4)        ,   intent(in)    :: itimestep, irk_step, idstep, itstep
        type(scalar_1d_field), intent(in)    :: tend_hpressure_at_pc_surface_cnty
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
        type(scalar_2d_field), intent(inout) :: scalar_alpha_at_pc_full_level_np1          ! out
! local
        type(scalar_2d_field)                :: tend_et_phi_at_pc_face_level_rk
        type(scalar_2d_field)                :: tend_et_www_at_pc_face_level_rk
        !integer, save                        :: irk_step = 0
#ifndef NHD_NOOPT
        type(exchange_field_list_2d),pointer :: field_head_2d
        field_head_2d=>null()
#endif

        !irk_step = irk_step + 1  ! yiz: already as input now

! just init
        tend_et_phi_at_pc_face_level_rk   = scalar_pressure_at_pc_face_level_np1
        tend_et_www_at_pc_face_level_rk   = scalar_pressure_at_pc_face_level_np1

        if(mpi_rank() == 0.and.write_verbose) print*,"grist_nh_dynamics_run_explicit"
        call grist_nh_dynamics_run_explicit(mesh, nlevp, dtime, itimestep, irk_step             , &
                                            tend_hpressure_at_pc_surface_cnty            , &
                                            scalar_normal_mass_flux_at_edge_face_level_rk, &
                                            tend_mass_at_pc_face_level_hori_rk           , &
                                            scalar_eta_mass_flux_at_pc_full_level_rk     , &
                                            scalar_phi_at_pc_face_level_rk               , &
                                            scalar_www_at_pc_face_level_n                , &
                                            scalar_www_at_pc_face_level_rk               , &
                                            scalar_delhp_at_pc_face_level_rk             , &
                                            scalar_geopotential_at_pc_surface_n          , &
                                            tend_et_phi_at_pc_face_level_rk              , & ! need
                                            tend_et_www_at_pc_face_level_rk)                 ! need
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

        if(mpi_rank() == 0.and.write_verbose) print*,"grist_nh_dynamics_run_implicit"
        !need change halo(1)
        call grist_nh_dynamics_run_implicit(mesh, nlevp, dtime,                               &
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

        scalar_phi_at_pc_face_level_np1%f(nlev+1,:) = scalar_geopotential_at_pc_surface_n%f

#ifndef NHD_NOOPT

#ifdef LAM_DOMAIN
        call update_assignValueArea_2d(mesh,nlevp,"wwwFace", scalar_www_at_pc_face_level_np1)    ! zero dw/dt at data zone
        call update_assignValueArea_2d(mesh,nlevp,"phiFace", scalar_phi_at_pc_face_level_np1)
        call update_assignValueArea_2d(mesh,nlevp,"tend_et_www",tend_et_www_at_pc_face_level_rk) ! zero dw/dt at data zone
#endif
        call exchange_data_2d_add(mesh,field_head_2d,tend_et_www_at_pc_face_level_rk)
        call exchange_data_2d_add(mesh,field_head_2d,scalar_www_at_pc_face_level_np1)
        call exchange_data_2d_add(mesh,field_head_2d,scalar_phi_at_pc_face_level_np1)
        call exchange_data_2d(mesh%local_block,field_head_2d)
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

        if(mpi_rank() == 0.and.write_verbose) print*," grist_nh_dynamics_run_diagvars"
        call grist_nh_dynamics_run_diagvars(mesh, nlevp, itimestep, irk_step,idstep, itstep, dtime, &
                                                  scalar_delp_at_pc_face_level_rk,           &  ! in
                                                  scalar_hpressure_at_pc_face_level_np1,     &  ! in
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
                                                  scalar_alpha_at_pc_full_level_np1)
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

     subroutine grist_nh_dynamics_run_explicit(mesh, nlevp, dtime, itimestep, irk                  , &
                                               tend_hpressure_at_pc_surface_cnty            , &
                                               scalar_normal_mass_flux_at_edge_face_level_rk, &
                                               tend_mass_at_pc_face_level_hori_rk           , &
                                               scalar_eta_mass_flux_at_pc_full_level_rk     , &
                                               scalar_phi_at_pc_face_level_rk               , &
                                               scalar_www_at_pc_face_level_n                , &
                                               scalar_www_at_pc_face_level_rk               , &
                                               scalar_delhp_at_pc_face_level_rk             , &
                                               scalar_geopotential_at_pc_surface_n          , &
                                               tend_et_phi_at_pc_face_level_rk              , &
                                               tend_et_www_at_pc_face_level_rk)
! io
        use omp_lib
        type(global_domain),   intent(inout)  :: mesh
        integer(i4)        ,   intent(in)     :: nlevp
        real(r8)           ,   intent(in)     :: dtime
        integer(i4)        ,   intent(in)     :: itimestep, irk
        type(scalar_1d_field), intent(in)     :: tend_hpressure_at_pc_surface_cnty
        type(scalar_2d_field), intent(in)     :: scalar_normal_mass_flux_at_edge_face_level_rk
        type(scalar_2d_field), intent(in)     :: tend_mass_at_pc_face_level_hori_rk
        type(scalar_2d_field), intent(in)     :: scalar_eta_mass_flux_at_pc_full_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_phi_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_n
        type(scalar_2d_field), intent(in)     :: scalar_www_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_delhp_at_pc_face_level_rk
        type(scalar_1d_field), intent(in)     :: scalar_geopotential_at_pc_surface_n
        type(scalar_2d_field), intent(inout)  :: tend_et_phi_at_pc_face_level_rk
        type(scalar_2d_field), intent(inout)  :: tend_et_www_at_pc_face_level_rk
! local
        integer(i4)                           :: iv, ilev, ie, v1, v2, index_edge, ierr
        real(r8)                              :: tmp_sum, www_max, www_max_global
        real(r8)                              :: www_edge_bot(mesh%ne_full)

           call grist_nh_et_adv_face(mesh, nlevp-1, nlevp, dtime, &
                                     scalar_normal_mass_flux_at_edge_face_level_rk ,& ! DELHP*V
                                     tend_mass_at_pc_face_level_hori_rk%f          ,& ! div.(mass)
                                     scalar_eta_mass_flux_at_pc_full_level_rk%f    ,& ! m*etadot
                                     scalar_delhp_at_pc_face_level_rk%f            ,&
                                     scalar_www_at_pc_face_level_rk%f              ,&
                                     scalar_phi_at_pc_face_level_rk%f              ,&
                                     tend_et_www_at_pc_face_level_rk%f             ,&
                                     tend_et_phi_at_pc_face_level_rk%f             ,&
                                     www_adv_flag(irk))
           
           scalar_www_top(:) = 0._r8
           scalar_www_bot(:) = tend_et_phi_at_pc_face_level_rk%f(nlevp,:)/gravity

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

    subroutine grist_nh_dynamics_run_implicit(mesh, nlevp, dtime,                               &
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
! io
        use omp_lib
        type(global_domain),   intent(in)     :: mesh
        integer(i4)        ,   intent(in)     :: nlevp
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
! local
        integer(i4)                           :: iv, ilev
        real(r8)                              :: tmp_coef,cr_num,cr_act
        integer                               :: ii, ncell_do, nlev

         nlev=nlevp-1


           call grist_nh_solve_www_equation(mesh, nlev, nlevp, dtime, &
                                            scalar_pressure_at_pc_full_level_n%f,        &
                                            scalar_delhp_at_pc_full_level_n%f,           &
                                            scalar_potential_temp_at_pc_full_level_n%f,  &
                                            scalar_potential_temp_at_pc_full_level_np1%f,&
                                            scalar_delhp_at_pc_full_level_np1%f,         &
                                            scalar_phi_at_pc_face_level_n%f,             &
                                            scalar_www_at_pc_face_level_n%f,             &
                                            scalar_www_at_pc_face_level_rk%f,            &
                                            tend_et_www_at_pc_face_level_rk%f,           &
                                            tend_et_phi_at_pc_face_level_rk%f,           &
                                            scalar_delhp_at_pc_face_level_np1%f,         &
                                            scalar_delp_at_pc_face_level_rk%f,           &
                                            scalar_delhp_at_pc_face_level_rk%f,          &
                                            scalar_www_at_pc_face_level_np1%f,           &
                                            scalar_mif_at_pc_face_level_n%f,             &
                                            scalar_mif_at_pc_full_level_n%f,             &
                                            scalar_www_top, scalar_www_bot)

        ncell_do = mesh%nv_halo(1)
!
! optional w-damping as in KDH08 for upper levels
!
!$omp parallel  private(ii,iv,ilev,tmp_coef) 
!$omp do schedule(dynamic,100) 
     do ii = 1, ncell_do*(nlev+1),1
        iv=ceiling(ii/real((nlev+1),r8))
        ilev=ii-(iv-1)*(nlev+1)

!           do iv = 1, ncell_do
!              do ilev = 1, nlev+1
                 if(zd.gt.0._r8.and.www_damping_coef.gt.0._r8.and.&
                    scalar_phi_at_pc_face_level_n%f(ilev,iv).ge.(scalar_phi_at_pc_face_level_n%f(1,iv)-zd*gravity))then
                    tmp_coef = www_damping_coef*(sin((pi*0.5_r8)*(1._r8-(scalar_phi_at_pc_face_level_n%f(1,iv)-scalar_phi_at_pc_face_level_n%f(ilev,iv))/(zd*gravity)))**2)
                    scalar_www_at_pc_face_level_np1%f(ilev,iv) = scalar_www_at_pc_face_level_np1%f(ilev,iv)/(1._r8+dtime*tmp_coef)
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
     do ii = 1,ncell_do*(nlev-1),1
        iv=ceiling(ii/real((nlev-1),r8))
        ilev=1+ii-(iv-1)*(nlev-1)

!           do iv = 1, ncell_do
!              do ilev = 2, nlev
                 cr_num = scalar_eta_mass_flux_at_pc_face_level_rk%f(ilev,iv)*dtime/scalar_delhp_at_pc_face_level_rk%f(ilev,iv)
                 cr_act = max(0._r8,abs(cr_num)-1._r8)
                 if(cr_act.gt.0._r8) print*,"violation of vertical cfl at, ilev=",ilev
                 cr_act = cr_act*sign(1._r8,scalar_www_at_pc_face_level_rk%f(ilev,iv))
                 scalar_www_at_pc_face_level_np1%f(ilev,iv) = scalar_www_at_pc_face_level_np1%f(ilev,iv)-0.3_r8*cr_act*dtime
!              end do
!           end do
     end do
!$omp end do nowait
!$omp end parallel 

!
! update phi equation, not updated one time due to an intel bug
!
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_n%f(1:nlev,:)
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_np1%f(1:nlev,:)-dtime*tend_et_phi_at_pc_face_level_rk%f(1:nlev,:)
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_np1%f(1:nlev,:)+&
                                                         (1._r8-imbeta)*gravity*dtime*scalar_www_at_pc_face_level_rk%f(1:nlev,:)
           scalar_phi_at_pc_face_level_np1%f(1:nlev,:) = scalar_phi_at_pc_face_level_np1%f(1:nlev,:)+&
                                                                 imbeta*gravity*dtime*scalar_www_at_pc_face_level_np1%f(1:nlev,:)

       return
     end subroutine grist_nh_dynamics_run_implicit

     subroutine grist_nh_dynamics_run_diagvars(mesh, nlevp, itimestep, irk_step, idstep, itstep, &
                                                     dtime,    &
                                                     scalar_delp_at_pc_face_level_rk,         &
                                                     scalar_hpressure_at_pc_face_level_np1,   &
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
                                                     scalar_alpha_at_pc_full_level_np1)
! io
        use omp_lib
        type(global_domain),   intent(in)     :: mesh
        integer(i4),           intent(in)     :: nlevp, itimestep, irk_step, idstep, itstep 
        real(r8),              intent(in)     :: dtime
        type(scalar_2d_field), intent(in)     :: scalar_delp_at_pc_face_level_rk
        type(scalar_2d_field), intent(in)     :: scalar_hpressure_at_pc_face_level_np1
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
        type(scalar_2d_field), intent(inout)  :: scalar_alpha_at_pc_full_level_np1
! local
        real(r8)       :: sound_speed,  sound_speed_max,  sound_speed_max_global
        real(r8)       :: sound_speed1, sound_speed1_max, sound_speed1_max_global
        real(r8)       :: max_kesi,min_kesi, tmp1(mesh%nv_full), tmp2(mesh%nv_full)
        integer(i4)    :: ie,ilev, iv, v1, v2, ierr
        integer(i4)    :: ii, ncell_do, nlev

        ncell_do = mesh%nv_halo(1)

        nlev = nlevp-1

!
! compute density with updated phi and delhp, then compute new pressure and temp
!
!
! par-seq consistenty under intel-O2 sensitive to loop implementation on some machines
! default not use this imploop, leading to nan at cma-pi for par, but not seq
! remember this experience when similar situations come out in other machines
!

#ifdef NHDRV_IMPLOOP
           do ilev = 1, nlev

              alpha_np1(:) = scalar_phi_at_pc_face_level_np1%f(ilev+1,:)-scalar_phi_at_pc_face_level_np1%f(ilev,:)
              alpha_n(:)   = scalar_phi_at_pc_face_level_n%f(ilev+1,:)  -scalar_phi_at_pc_face_level_n%f(ilev,:)
              scalar_alpha_at_pc_full_level_np1%f(ilev,:) = -alpha_np1(:)/scalar_delhp_at_pc_full_level_np1%f(ilev,:) ! ensure >0
              scalar_pressure_at_pc_full_level_np1%f(ilev,:) = scalar_pressure_at_pc_full_level_n%f(ilev,:)+&
                                                               (cp/cv)*scalar_pressure_at_pc_full_level_n%f(ilev,:)*&
                                                              ((scalar_potential_temp_at_pc_full_level_np1%f(ilev,:)/&
                                                                scalar_potential_temp_at_pc_full_level_n%f(ilev,:))-&
                                                                alpha_np1(:)/alpha_n(:))

! do a 3d div damping
              scalar_pressure_at_pc_full_level_np1%f(ilev,:) = scalar_pressure_at_pc_full_level_np1%f(ilev,:)+&
                                                               d3d_damping_coef*&
                                                              (scalar_pressure_at_pc_full_level_np1%f(ilev,:)-&
                                                               scalar_pressure_at_pc_full_level_n%f(ilev,:))

              scalar_temp_at_pc_full_level_np1%f(ilev,:) = scalar_potential_temp_at_pc_full_level_np1%f(ilev,:)/scalar_delhp_at_pc_full_level_np1%f(ilev,:)*&
                                                          ((scalar_pressure_at_pc_full_level_np1%f(ilev,:)/p00)**(rdry/cp))/&
                                                          (one+ptfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,ilev,:))
            
          end do
#else
!$omp parallel  private(iv,ilev)   
!$omp do schedule(dynamic,20) 
          do iv = 1, ncell_do 
             do ilev = 1, nlev

              alpha_np1(iv) = scalar_phi_at_pc_face_level_np1%f(ilev+1,iv)-scalar_phi_at_pc_face_level_np1%f(ilev,iv)
              alpha_n(iv)   = scalar_phi_at_pc_face_level_n%f(ilev+1,iv)  -scalar_phi_at_pc_face_level_n%f(ilev,iv)
              scalar_alpha_at_pc_full_level_np1%f(ilev,iv) = -alpha_np1(iv)/scalar_delhp_at_pc_full_level_np1%f(ilev,iv) ! ensure >0
              scalar_pressure_at_pc_full_level_np1%f(ilev,iv) = scalar_pressure_at_pc_full_level_n%f(ilev,iv)+&
                                                               (cp/cv)*scalar_pressure_at_pc_full_level_n%f(ilev,iv)*&
                                                              ((scalar_potential_temp_at_pc_full_level_np1%f(ilev,iv)/&
                                                                scalar_potential_temp_at_pc_full_level_n%f(ilev,iv))-&
                                                                alpha_np1(iv)/alpha_n(iv))
! do a 3d div damping
              scalar_pressure_at_pc_full_level_np1%f(ilev,iv) = scalar_pressure_at_pc_full_level_np1%f(ilev,iv)+&
                                                               d3d_damping_coef*&
                                                              (scalar_pressure_at_pc_full_level_np1%f(ilev,iv)-&
                                                               scalar_pressure_at_pc_full_level_n%f(ilev,iv))

              scalar_temp_at_pc_full_level_np1%f(ilev,iv) = scalar_potential_temp_at_pc_full_level_np1%f(ilev,iv)/scalar_delhp_at_pc_full_level_np1%f(ilev,iv)*&
                                                          ((scalar_pressure_at_pc_full_level_np1%f(ilev,iv)/p00)**(rdry/cp))/&
                                                          (one+ptfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,ilev,iv))
             end do
           end do
!$omp end do nowait
!$omp end parallel 

#endif

           sound_speed_max  = 0._r8
           sound_speed1_max = 0._r8
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
#else
           sound_speed_max_global = sound_speed_max
           sound_speed1_max_global= sound_speed1_max
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

#ifdef NHDRV_IMPLOOP2
!
! compute geopotential at face and full levels, just use another name
!
           scalar_geopotential_at_pc_face_level_np1%f  = scalar_phi_at_pc_face_level_np1%f
           do ilev = 1, nlev
              scalar_geopotential_at_pc_full_level_np1%f(ilev,:) = 0.5_r8*(scalar_geopotential_at_pc_face_level_np1%f(ilev+1,:)+&
                                                                           scalar_geopotential_at_pc_face_level_np1%f(ilev,:))
           end do
! face level pressure, as in hdc
           scalar_pressure_at_pc_face_level_np1%f(1,:) = scalar_hpressure_at_pc_face_level_np1%f(1,:)

           do ilev = 2, nlev+1
#ifdef FACEP1
              scalar_pressure_at_pc_face_level_np1%f(ilev,:) = scalar_pressure_at_pc_face_level_np1%f(ilev-1,:)+scalar_delhp_at_pc_full_level_np1%f(ilev-1,:)*&
                                                               (0.5_r8/gravity*&
                                                               (((scalar_www_at_pc_face_level_np1%f(ilev-1,:)-scalar_www_at_pc_face_level_n%f(ilev-1,:))/dtime+&
                                                                  tend_et_www_at_pc_face_level_rk%f(ilev-1,:))+&
                                                                ((scalar_www_at_pc_face_level_np1%f(ilev,:)-scalar_www_at_pc_face_level_n%f(ilev,:))/dtime+&
                                                                  tend_et_www_at_pc_face_level_rk%f(ilev,:)))+1._r8)/scalar_mif_at_pc_full_level_n%f(ilev-1,:)
#else 
              tmp1(:) =(((scalar_www_at_pc_face_level_np1%f(ilev-1,:)-scalar_www_at_pc_face_level_n%f(ilev-1,:))/dtime+tend_et_www_at_pc_face_level_rk%f(ilev-1,:))/gravity+one)/&
                       scalar_mif_at_pc_face_level_n%f(ilev-1,:)
              tmp2(:) =(((scalar_www_at_pc_face_level_np1%f(ilev,:)-scalar_www_at_pc_face_level_n%f(ilev,:))/dtime+tend_et_www_at_pc_face_level_rk%f(ilev,:))/gravity+one)/&
                       scalar_mif_at_pc_face_level_n%f(ilev,:)
              scalar_pressure_at_pc_face_level_np1%f(ilev,:) = scalar_pressure_at_pc_face_level_np1%f(ilev-1,:)+scalar_delhp_at_pc_full_level_np1%f(ilev-1,:)*&
                                                               (tmp1(:)+tmp2(:))*0.5_r8
#endif
           end do

!
! compute delp at full level
!
        do ilev = 1, nlev
           scalar_delp_at_pc_full_level_np1%f(ilev,:) = scalar_pressure_at_pc_face_level_np1%f(ilev+1,:)-&
                                                        scalar_pressure_at_pc_face_level_np1%f(ilev,:)
        end do

! NHDRV_IMPLOOP2
#else
        
!
! compute geopotential at face and full levels, just use another name
!
        scalar_geopotential_at_pc_face_level_np1%f  = scalar_phi_at_pc_face_level_np1%f

!$omp parallel  private(iv,ilev)
!$omp do schedule(dynamic,20)
        do iv = 1, ncell_do
           do ilev = 1, nlev
              scalar_geopotential_at_pc_full_level_np1%f(ilev,iv) = 0.5_r8*(scalar_geopotential_at_pc_face_level_np1%f(ilev+1,iv)+&
                                                                            scalar_geopotential_at_pc_face_level_np1%f(ilev,iv))
           end do
! face level pressure, as in hdc
           scalar_pressure_at_pc_face_level_np1%f(1,iv) = scalar_hpressure_at_pc_face_level_np1%f(1,iv) !*&
                                      !                 (one+prfactor*scalar_tracer_mxrt_at_pc_full_level_n%f(1,1,:))
        end do
!$omp end do nowait
!$omp end parallel 

!$omp parallel default(shared) private(iv,ilev)  
!$omp do schedule(dynamic,20) 
        do iv = 1, ncell_do 
           do ilev = 2, nlev+1
#ifdef FACEP1
              scalar_pressure_at_pc_face_level_np1%f(ilev,iv) = scalar_pressure_at_pc_face_level_np1%f(ilev-1,iv)+scalar_delhp_at_pc_full_level_np1%f(ilev-1,iv)*&
                                                               (0.5_r8/gravity*&
                                                               (((scalar_www_at_pc_face_level_np1%f(ilev-1,iv)-scalar_www_at_pc_face_level_n%f(ilev-1,iv))/dtime+&
                                                                  tend_et_www_at_pc_face_level_rk%f(ilev-1,iv))+&
                                                                ((scalar_www_at_pc_face_level_np1%f(ilev,iv)-scalar_www_at_pc_face_level_n%f(ilev,iv))/dtime+&
                                                                  tend_et_www_at_pc_face_level_rk%f(ilev,iv)))+1._r8)/scalar_mif_at_pc_full_level_n%f(ilev-1,iv)
#else 
              tmp1(iv) =(((scalar_www_at_pc_face_level_np1%f(ilev-1,iv)-scalar_www_at_pc_face_level_n%f(ilev-1,iv))/dtime+tend_et_www_at_pc_face_level_rk%f(ilev-1,iv))/gravity+one)/&
                           scalar_mif_at_pc_face_level_n%f(ilev-1,iv)
              tmp2(iv) =(((scalar_www_at_pc_face_level_np1%f(ilev,iv)-scalar_www_at_pc_face_level_n%f(ilev,iv))/dtime+tend_et_www_at_pc_face_level_rk%f(ilev,iv))/gravity+one)/&
                           scalar_mif_at_pc_face_level_n%f(ilev,iv)
              scalar_pressure_at_pc_face_level_np1%f(ilev,iv) = scalar_pressure_at_pc_face_level_np1%f(ilev-1,iv)+scalar_delhp_at_pc_full_level_np1%f(ilev-1,iv)*&
                                                               (tmp1(iv)+tmp2(iv))*0.5_r8
#endif
           end do
        end do
!$omp end do nowait
!$omp end parallel 

!
! compute delp at full level
!
!$omp parallel  private(iv,ii,ilev)  
!$omp do schedule(dynamic,100) 
     do ii = 1, ncell_do*nlev,1
        iv=ceiling(ii/real(nlev,r8))
        ilev=ii-(iv-1)*nlev
!        do iv = 1, ncell_do 
!           do ilev = 1, nlev
              scalar_delp_at_pc_full_level_np1%f(ilev,iv) = scalar_pressure_at_pc_face_level_np1%f(ilev+1,iv)-&
                                                            scalar_pressure_at_pc_face_level_np1%f(ilev  ,iv)
!           end do
!        end do
     end do
!$omp end do nowait
!$omp end parallel 

#endif

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

        return
     end subroutine grist_nh_dynamics_init

     subroutine grist_nh_dynamics_final

        deallocate(scalar_www_top)
        deallocate(scalar_www_bot)
        deallocate(alpha_n)
        deallocate(alpha_np1)

        return
     end subroutine grist_nh_dynamics_final

   end module grist_nh_driver_module
