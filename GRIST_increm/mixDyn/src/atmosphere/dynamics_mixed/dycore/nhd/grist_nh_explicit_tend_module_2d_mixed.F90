module grist_nh_explicit_tend_module_2d_mixed

    use grist_constants,                                only: i4, r8, zero, half, rearth, r4 => ns, zero_r4 => zero_ns, half_r4 => half_ns
    use grist_domain_types,                             only: global_domain
    use grist_data_types,                               only: scalar_1d_field, scalar_2d_field, exchange_field_list_1d, exchange_field_list_2d
    use grist_nml_module,                               only: nlev, nlevp, eqs_vert_diff, write_verbose, ver_adv_flag
    use grist_lib,                                      only: mpi_rank
    use grist_dycore_gcd_recon_module_2d_mixed,         only: divergence_operator_2d_var2
    use grist_dycore_primal_flux_operators_2d_mixed,    only: calc_primal_normal_flux_at_edge_var2
! mpi
#ifndef SEQ_GRIST
    use grist_config_partition,       only: exchange_data_2d_add , exchange_data_2d, exchange_data_1d_add , exchange_data_1d, debug_data_2d, &
                                            exchange_data_2d_r4
#endif

! data
   !  use grist_hpe_constants,            only: deta_face, deta_full
    use grist_dycore_module_vars_control_mixed, only: deta_face, deta_full
    use grist_dycore_ref_atmos_mixed,   only: scalar_geop_at_pc_face_level_bar
    use grist_dycore_diffusion_module_mixed,  only: perot_weight_at_pc
    ! use grist_dycore_module_vars_control_mixed, only: perot_weight_at_pc

    implicit none

    real(r8), allocatable   :: scalar_template_1d_nlevp_c(:)

contains

!===================================================
! same as grist_nh_et_www_face but evaluate 2 vars
!===================================================

   subroutine grist_nh_et_adv_face(mesh, dtime,  &
                                   scalar_normal_mass_flux_at_edge_face_level ,& ! DELHP*V
                                   tend_mass_at_pc_face_level_hori_rk         ,& ! div.(mass)
                                   scalar_eta_mass_flux_at_pc_full_level      ,& ! m*etadot
                                   scalar_eta_mass_flux_at_pc_face_level      ,& ! m*etadot
                                   scalar_delhp_at_pc_face_level_rk           ,&
                                   scalar_delhp_at_pc_face_level_n            ,&
                                   scalar_delhp_at_pc_face_level_np1          ,&
                                   scalar_www_at_pc_face_level                ,&
                                   scalar_phi_at_pc_face_level                ,&
                                   tend_et_www_face_level                     ,& 
                                   tend_et_phi_face_level                     ,&
                                   adv_face_flag)
! io
      use omp_lib
      type(global_domain)  , intent(inout) :: mesh
      real(r8)             , intent(in)    :: dtime
      type(scalar_2d_field), intent(in)    :: scalar_normal_mass_flux_at_edge_face_level
      type(scalar_2d_field), intent(in)    :: tend_mass_at_pc_face_level_hori_rk
      type(scalar_2d_field), intent(in)    :: scalar_eta_mass_flux_at_pc_full_level
      type(scalar_2d_field), intent(in)    :: scalar_eta_mass_flux_at_pc_face_level
      type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_rk
      type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_n
      type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_face_level_np1
      type(scalar_2d_field), intent(in)    :: scalar_www_at_pc_face_level
      type(scalar_2d_field), intent(in)    :: scalar_phi_at_pc_face_level
      type(scalar_2d_field), intent(inout) :: tend_et_www_face_level
      type(scalar_2d_field), intent(inout) :: tend_et_phi_face_level
      integer(i4)          , intent(in)    :: adv_face_flag
! local
      type(scalar_2d_field)                :: scalar_phipie_at_pc_face_level
      type(scalar_2d_field)                :: tend1_hori_final, tend2_hori_final
      type(scalar_2d_field)                :: tend1_vert_tmp,   tend2_vert_tmp
      type(scalar_2d_field)                :: tend1_vert_final, tend2_vert_final
      type(scalar_2d_field)                :: scalar_template_2d_ne_a
      type(scalar_2d_field)                :: scalar_template_2d_ne_b
      type(scalar_1d_field)                :: scalar_template_1d_ne_a
      type(scalar_1d_field)                :: scalar_template_1d_ne_b
#ifndef SEQ_GRIST
      type(exchange_field_list_1d),pointer :: field_head_1d
      type(exchange_field_list_2d),pointer :: field_head_2d
#endif
      integer(i4)                          :: iv, ilev, icell1, icell2, ie, inb
      real(r4)                             :: vector_sum1(3), vector_sum2(3), vector_sum3(3)

#ifndef SEQ_GRIST 
        field_head_1d => null()
        field_head_2d => null()
#endif
! init
        if(.not.allocated(scalar_template_2d_ne_a%f_r4)) allocate(scalar_template_2d_ne_a%f_r4(nlevp,mesh%ne_full))
        if(.not.allocated(scalar_template_2d_ne_b%f_r4)) allocate(scalar_template_2d_ne_b%f_r4(nlevp,mesh%ne_full))
        if(.not.allocated(scalar_template_1d_ne_a%f)) allocate(scalar_template_1d_ne_a%f(mesh%ne_full))
        if(.not.allocated(scalar_template_1d_ne_b%f)) allocate(scalar_template_1d_ne_b%f(mesh%ne_full))

        scalar_template_2d_ne_a%pos  = 6; scalar_template_2d_ne_a%f_r4    = zero_r4
        scalar_template_2d_ne_b%pos  = 6; scalar_template_2d_ne_b%f_r4    = zero_r4
 
        tend1_hori_final%f_r4              = scalar_www_at_pc_face_level%f_r4   ! just init
        tend1_vert_final%f_r4              = scalar_www_at_pc_face_level%f_r4   ! just init
        tend1_vert_tmp%f_r4             = scalar_www_at_pc_face_level%f_r4   ! just init
        tend2_hori_final%f_r4              = scalar_phi_at_pc_face_level%f_r4   ! just init
        tend2_vert_final%f_r4              = scalar_phi_at_pc_face_level%f_r4   ! just init
        tend2_vert_tmp%f_r4           = scalar_phi_at_pc_face_level%f_r4   ! just init
        scalar_phipie_at_pc_face_level%f_r4= scalar_phi_at_pc_face_level%f_r4

! info
        if(mpi_rank() == 0.and.write_verbose) print*,"before compute horizontal face level mass flux diff"

        !if(.not.stencil_exchange_flag) mesh%ne = mesh%ne_halo(1)
        !if(.not.stencil_exchange_flag) mesh%nv = mesh%nv_halo(1)

!        select case(wwwphi_hori_adv_method)
!        case(1)
!
! hori evaluation
!
        call calc_primal_normal_flux_at_edge_var2(mesh, scalar_normal_mass_flux_at_edge_face_level%f_r4, &
                                                        scalar_normal_mass_flux_at_edge_face_level%f_r4, &
                                                        scalar_www_at_pc_face_level%f_r4               , &
                                                        scalar_phipie_at_pc_face_level%f_r4            , &
                                                        scalar_template_2d_ne_a%f_r4                   , &
                                                        scalar_template_2d_ne_b%f_r4                   , &
                                                        adv_face_flag, dtime, nlev+1, .true.)

        scalar_template_2d_ne_a%f_r4 = scalar_template_2d_ne_a%f_r4*scalar_normal_mass_flux_at_edge_face_level%f_r4
        scalar_template_2d_ne_b%f_r4 = scalar_template_2d_ne_b%f_r4*scalar_normal_mass_flux_at_edge_face_level%f_r4

        !if(.not.stencil_exchange_flag)then
        !  call divergence_operator_2d_var2(mesh, scalar_template_2d_ne_a%f, &
        !                                         scalar_template_2d_ne_b%f, &
        !                                         tend1_hori_final%f       , &
        !                                         tend2_hori_final%f       , &
        !                                         nlev+1)
        !  tend1_hori_final%f = tend1_hori_final%f-scalar_www_at_pc_face_level%f*tend_mass_at_pc_face_level_hori_rk%f
        !  tend2_hori_final%f = tend2_hori_final%f-scalar_phi_at_pc_face_level%f*tend_mass_at_pc_face_level_hori_rk%f
        !end if

        !if(stencil_exchange_flag)then
#ifndef SEQ_GRIST
        call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_a,"r4")
        call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_b,"r4")
        call exchange_data_2d_r4(mesh%local_block,field_head_2d)
#endif

        mesh%ne = mesh%ne_halo(1)
        mesh%nv = mesh%nv_halo(1)
        call divergence_operator_2d_var2(mesh, scalar_template_2d_ne_a%f_r4, &
                                               scalar_template_2d_ne_b%f_r4, &
                                               tend1_hori_final%f_r4       , &
                                               tend2_hori_final%f_r4       , &
                                               nlev+1)
        tend1_hori_final%f_r4 = tend1_hori_final%f_r4-scalar_www_at_pc_face_level%f_r4   *tend_mass_at_pc_face_level_hori_rk%f_r4
        tend2_hori_final%f_r4 = tend2_hori_final%f_r4-scalar_phipie_at_pc_face_level%f_r4*tend_mass_at_pc_face_level_hori_rk%f_r4
        !end if

        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute

!         case(2)
! ! yizhang, added newly 2021 Nov, faster direct method; works for idealized modeling, but not very stable for real-world VR

!         do ie = 1, mesh%ne_halo(1)
!            icell1 = mesh%edt_v(1,ie)
!            icell2 = mesh%edt_v(2,ie)
!            do ilev = 1, nlevp
!               scalar_template_2d_ne_a%f(ilev,ie) = (scalar_www_at_pc_face_level%f(ilev   ,icell2)-scalar_www_at_pc_face_level%f(ilev   ,icell1))/(rearth*mesh%edt_leng(ie))
!               scalar_template_2d_ne_b%f(ilev,ie) = (scalar_phipie_at_pc_face_level%f(ilev,icell2)-scalar_phipie_at_pc_face_level%f(ilev,icell1))/(rearth*mesh%edt_leng(ie))
!            end do
!         end do

!         do iv  = 1, mesh%nv_halo(1)
!            do ilev = 1, nlevp
!              vector_sum1 = zero
!              vector_sum2 = zero
!              vector_sum3 = zero
!              do inb  = 1, mesh%vtx_nnb(iv)
!                 ie = mesh%vtx_ed(inb,iv)
!                 vector_sum1(:) = vector_sum1(:) + perot_weight_at_pc(:,inb,iv)*scalar_normal_mass_flux_at_edge_face_level%f(ilev,ie)
!                 vector_sum2(:) = vector_sum2(:) + perot_weight_at_pc(:,inb,iv)*scalar_template_2d_ne_a%f(ilev,ie)
!                 vector_sum3(:) = vector_sum3(:) + perot_weight_at_pc(:,inb,iv)*scalar_template_2d_ne_b%f(ilev,ie)
!            end do
!            tend1_hori_final%f(ilev,iv) = dot_product(vector_sum1,vector_sum2) ! www
!            tend2_hori_final%f(ilev,iv) = dot_product(vector_sum1,vector_sum3) ! phi
!           end do
!         end do

!         mesh%ne = mesh%ne_compute
!         mesh%nv = mesh%nv_compute

!         case(3)
! ! as case 1, but do not use residule directly, use it indirectly
! ! must be used together with vert_adv_method=3
! ! This is for testing

! !
! ! hori evaluation
! !
!         call calc_primal_normal_flux_at_edge_var2(mesh, scalar_normal_mass_flux_at_edge_face_level%f, &
!                                                         scalar_normal_mass_flux_at_edge_face_level%f, &
!                                                         scalar_www_at_pc_face_level%f               , &
!                                                         scalar_phipie_at_pc_face_level%f            , &
!                                                         scalar_template_2d_ne_a%f                   , &
!                                                         scalar_template_2d_ne_b%f                   , &
!                                                         adv_face_flag, dtime, nlev+1, .true.)

!         scalar_template_2d_ne_a%f = scalar_template_2d_ne_a%f*scalar_normal_mass_flux_at_edge_face_level%f
!         scalar_template_2d_ne_b%f = scalar_template_2d_ne_b%f*scalar_normal_mass_flux_at_edge_face_level%f
! #ifndef SEQ_GRIST
!         call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_a)
!         call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_b)
!         call exchange_data_2d(mesh%local_block,field_head_2d)
! #endif

!         mesh%ne = mesh%ne_halo(1)
!         mesh%nv = mesh%nv_halo(1)
!         call divergence_operator_2d_var2(mesh, scalar_template_2d_ne_a%f, &
!                                                scalar_template_2d_ne_b%f, &
!                                                tend1_hori_final%f       , &
!                                                tend2_hori_final%f       , &
!                                                nlev+1)
!         do iv = 1, mesh%nv_halo(1)
!            do ilev = 1, nlev
!               tend1_hori_final%f(ilev,iv) = tend1_hori_final%f(ilev,iv)+scalar_www_at_pc_face_level%f(ilev,iv) *&
!                                            (scalar_delhp_at_pc_face_level_np1%f(ilev,iv)-scalar_delhp_at_pc_face_level_n%f(ilev,iv))/dtime
!               tend2_hori_final%f(ilev,iv) = tend2_hori_final%f(ilev,iv)+scalar_phipie_at_pc_face_level%f(ilev,iv)*&
!                                            (scalar_delhp_at_pc_face_level_np1%f(ilev,iv)-scalar_delhp_at_pc_face_level_n%f(ilev,iv))/dtime
!            end do
!            tend1_hori_final%f(nlevp,iv) = tend1_hori_final%f(nlevp,iv)-scalar_www_at_pc_face_level%f(nlevp,iv)   *tend_mass_at_pc_face_level_hori_rk%f(nlevp,iv)
!            tend2_hori_final%f(nlevp,iv) = tend2_hori_final%f(nlevp,iv)-scalar_phipie_at_pc_face_level%f(nlevp,iv)*tend_mass_at_pc_face_level_hori_rk%f(nlevp,iv)
!         end do
!         !end if

!         mesh%ne = mesh%ne_compute
!         mesh%nv = mesh%nv_compute

!        case default
!            print*, "one must set wwwphi_hori_adv_method in dycore_para"
!            call mpi_abort()
!        end select
! info
        if(mpi_rank() == 0.and.write_verbose) print*,"finish compute horizontal face level mass flux diff"
!
! vert evaluation
!
!      select case(wwwphi_vert_adv_method)
!      case(1)
           call  calc_hpe_tend_vert_mass_flux_face_2d_var2(mesh,&
                                                           scalar_www_at_pc_face_level%f_r4    , & ! face level scalar1
                                                           scalar_phipie_at_pc_face_level%f_r4 , & ! face level scalar2
#ifdef DBL_MASS
                                                           scalar_eta_mass_flux_at_pc_full_level%f , & ! full level etadot
#else
                                                           scalar_eta_mass_flux_at_pc_full_level%f_r4 , & ! full level etadot
#endif
                                                           ver_adv_flag                  , &
                                                           tend1_vert_tmp%f_r4, tend2_vert_tmp%f_r4)    ! face level tend

           scalar_template_1d_nlevp_c(1)      = zero_r4
           scalar_template_1d_nlevp_c(nlev+1) = zero_r4
!$omp parallel  private(iv) firstprivate(scalar_template_1d_nlevp_c)
!$omp do schedule(static,5)
           do iv = 1, mesh%nv_halo(1)
#ifdef DBL_MASS
              scalar_template_1d_nlevp_c(2:nlev) = scalar_eta_mass_flux_at_pc_full_level%f(2:nlev,iv)-scalar_eta_mass_flux_at_pc_full_level%f(1:nlev-1,iv)
#else
              scalar_template_1d_nlevp_c(2:nlev) = scalar_eta_mass_flux_at_pc_full_level%f_r4(2:nlev,iv)-scalar_eta_mass_flux_at_pc_full_level%f_r4(1:nlev-1,iv)
#endif
              tend1_vert_final%f_r4(:,iv)    = tend1_vert_tmp%f_r4(:,iv)-scalar_www_at_pc_face_level%f_r4(:,iv)*scalar_template_1d_nlevp_c(:)
              tend2_vert_final%f_r4(:,iv)    = tend2_vert_tmp%f_r4(:,iv)-scalar_phipie_at_pc_face_level%f_r4(:,iv)*scalar_template_1d_nlevp_c(:)
           end do
!$omp end do nowait
!$omp end parallel
!       case(2)
! ! yizhang added newly, 2021 nov

!           do iv = 1, mesh%nv_halo(1)
!              do ilev = 2, nlev
!                 tend1_vert_final%f(ilev,iv) = 0.5*(scalar_www_at_pc_face_level%f(ilev+1,iv)-scalar_www_at_pc_face_level%f(ilev-1,iv))*&
!                                                    scalar_eta_mass_flux_at_pc_face_level%f(ilev,iv)
!                 tend2_vert_final%f(ilev,iv) = 0.5*(scalar_phipie_at_pc_face_level%f(ilev+1,iv)-scalar_phipie_at_pc_face_level%f(ilev-1,iv))*&
!                                                    scalar_eta_mass_flux_at_pc_face_level%f(ilev,iv)
!              end do
!                 tend1_vert_final%f(1,iv)     = zero
!                 tend2_vert_final%f(1,iv)     = zero
!                 tend1_vert_final%f(nlevp,iv) = zero
!                 tend2_vert_final%f(nlevp,iv) = zero
!           end do

!       case(3)
! ! as case 1, but using indirectly residule, 

!            call  calc_hpe_tend_vert_mass_flux_face_2d_var2(mesh,&
!                                                            scalar_www_at_pc_face_level%f    , & ! face level scalar1
!                                                            scalar_phipie_at_pc_face_level%f , & ! face level scalar2
!                                                            scalar_eta_mass_flux_at_pc_full_level%f , & ! full level etadot
!                                                            ver_adv_flag                  , &
!                                                            tend1_vert_tmp%f, tend2_vert_tmp%f)  ! face level tend

! !$omp parallel  private(iv) firstprivate(scalar_template_1d_nlevp_c)
! !$omp do schedule(static,5)
!            do iv = 1, mesh%nv_halo(1)
!               tend1_vert_final%f(1:nlevp,iv)    = tend1_vert_tmp%f(1:nlevp,iv)
!               tend2_vert_final%f(1:nlevp,iv)    = tend2_vert_tmp%f(1:nlevp,iv)
!            end do
! !$omp end do nowait
! !$omp end parallel
!      case default
!          print*,"one must set wwwphi_vert_adv_method in dycore_para"
!          call mpi_abort()
!    end select

! ASSIGN; still divided by local delhp at face leve
        tend_et_www_face_level%f_r4   = (tend1_hori_final%f_r4+tend1_vert_final%f_r4)/scalar_delhp_at_pc_face_level_rk%f
        tend_et_phi_face_level%f_r4   = (tend2_hori_final%f_r4+tend2_vert_final%f_r4)/scalar_delhp_at_pc_face_level_rk%f
! info
        if(mpi_rank() == 0.and.write_verbose) print*,"finish compute vertical www&phi mass flux diff"
! clean

!================================================
! FOR TESTING PURPOSE
!================================================
! #ifdef USE_PEROTWS
! ! grad component at edge
! !$omp parallel  private(ie,icell1,icell2) 
! !$omp do schedule(static,5)
!         do ie = 1, mesh%ne_halo(1)
!            icell1 = mesh%edt_v(1,ie)
!            icell2 = mesh%edt_v(2,ie)
!            scalar_template_1d_ne_a%f(ie) = (scalar_phipie_at_pc_face_level%f_r4(nlev+1,icell2)-scalar_phipie_at_pc_face_level%f_r4(nlev+1,icell1))/(rearth*mesh%edt_leng(ie))
!            scalar_template_1d_ne_b%f(ie) = 2._r8*scalar_normal_mass_flux_at_edge_face_level%f_r4(nlev+1,ie)/&
!                                                 (scalar_delhp_at_pc_face_level_rk%f(nlev+1,icell1)+scalar_delhp_at_pc_face_level_rk%f(nlev+1,icell2))
!         end do
! !$omp end do nowait
! !$omp end parallel 
!         !call exchange_data_1d_add(mesh,field_head_1d,scalar_template_1d_ne_a)
!         !call exchange_data_1d_add(mesh,field_head_1d,scalar_template_1d_ne_b)
!         !call exchange_data_1d(mesh%local_block,field_head_1d)

! !$omp parallel  private(iv,vector_sum1,vector_sum2,inb,ie)
! !$omp do schedule(dynamic,5)
!         do iv  = 1, mesh%nv_halo(1)
!            vector_sum1 = 0._r4
!            vector_sum2 = 0._r4
!            do inb  = 1, mesh%vtx_nnb(iv)
!               ie = mesh%vtx_ed(inb,iv)
!               vector_sum1(:) = vector_sum1(:) + perot_weight_at_pc(:,inb,iv)*scalar_template_1d_ne_a%f(ie)
!               vector_sum2(:) = vector_sum2(:) + perot_weight_at_pc(:,inb,iv)*scalar_template_1d_ne_b%f(ie)
!            end do
!            tend_et_phi_face_level%f(nlev+1,iv) = dot_product(vector_sum1,vector_sum2)
!         end do
! !$omp end do nowait
! !$omp end parallel 
! #endif
!================================================
! FOR TESTING PURPOSE
!================================================

        deallocate(scalar_template_2d_ne_a%f_r4,scalar_template_2d_ne_b%f_r4,scalar_template_1d_ne_a%f,scalar_template_1d_ne_b%f)
        deallocate(tend1_hori_final%f_r4,tend2_hori_final%f_r4,tend1_vert_final%f_r4,tend2_vert_final%f_r4)

     return
   end subroutine grist_nh_et_adv_face

!
! same as above, but for two vars
!

   subroutine calc_hpe_tend_vert_mass_flux_face_2d_var2(mesh,&
                                                   scalar1_at_pc_face_level            , &
                                                   scalar2_at_pc_face_level            , &
                                                   scalar_eta_mass_velocity_at_pc_full , &
                                                   order                               , &
                                                   tend1_vert_mass_flux_at_pc_face     , &
                                                   tend2_vert_mass_flux_at_pc_face     )
! io
    use omp_lib
    type(global_domain),    intent(in)   :: mesh
    real(r4), allocatable,  intent(in)   :: scalar1_at_pc_face_level(:,:)
    real(r4), allocatable,  intent(in)   :: scalar2_at_pc_face_level(:,:)
#ifdef DBL_MASS
    real(r8), allocatable,  intent(in)   :: scalar_eta_mass_velocity_at_pc_full(:,:)
#else
    real(r4), allocatable,  intent(in)   :: scalar_eta_mass_velocity_at_pc_full(:,:)
#endif
    integer(i4)          ,  intent(in)   :: order
    real(r4), allocatable,  intent(inout):: tend1_vert_mass_flux_at_pc_face(:,:)  ! with minus sign
    real(r4), allocatable,  intent(inout):: tend2_vert_mass_flux_at_pc_face(:,:)  ! with minus sign
! local
    real(r4),  allocatable               :: scalar1_at_pc_full_level(:,:)
    real(r4),  allocatable               :: scalar2_at_pc_full_level(:,:)
    integer(i4)                          :: iv,ilev, ii

    if(.not.allocated(scalar1_at_pc_full_level)) allocate(scalar1_at_pc_full_level(nlev,mesh%nv_full))
    if(.not.allocated(scalar2_at_pc_full_level)) allocate(scalar2_at_pc_full_level(nlev,mesh%nv_full))

    call calc_hpe_vert_flux_operator_face_2d_var2(mesh,&
                                             scalar1_at_pc_face_level   ,&
                                             scalar2_at_pc_face_level   ,&
                                             scalar_eta_mass_velocity_at_pc_full ,&
                                             order                      ,&
                                             scalar1_at_pc_full_level   ,&
                                             scalar2_at_pc_full_level )

!$omp parallel  private(ii,iv,ilev) 
!$omp do schedule(dynamic,100) 
    do ii=1,mesh%nv_halo(1)*(nlev-1),1
        iv=ceiling(ii/real(nlev-1,r8))
        ilev=1+ii-(iv-1)*(nlev-1)
!     do iv = 1, mesh%nv_halo(1)
!     do ilev =  2, nlev
! no minus sign here
        tend1_vert_mass_flux_at_pc_face(ilev,iv) = (scalar_eta_mass_velocity_at_pc_full(ilev  ,iv)*scalar1_at_pc_full_level(ilev  ,iv)-&
                                                    scalar_eta_mass_velocity_at_pc_full(ilev-1,iv)*scalar1_at_pc_full_level(ilev-1,iv))
        tend2_vert_mass_flux_at_pc_face(ilev,iv) = (scalar_eta_mass_velocity_at_pc_full(ilev  ,iv)*scalar2_at_pc_full_level(ilev  ,iv)-&
                                                    scalar_eta_mass_velocity_at_pc_full(ilev-1,iv)*scalar2_at_pc_full_level(ilev-1,iv))

!     end do
!     end do
    end do
!$omp end do nowait
!$omp end parallel

     tend1_vert_mass_flux_at_pc_face(1,:)      = 0._r4; tend2_vert_mass_flux_at_pc_face(1,:)      = 0._r4
     tend1_vert_mass_flux_at_pc_face(nlev+1,:) = 0._r4; tend2_vert_mass_flux_at_pc_face(nlev+1,:) = 0._r4

     deallocate(scalar1_at_pc_full_level)
     deallocate(scalar2_at_pc_full_level)

     return
   end subroutine calc_hpe_tend_vert_mass_flux_face_2d_var2

!
! private, from face level to full level
!

   subroutine calc_hpe_vert_flux_operator_face_2d_var2(mesh,&
                                                  scalar1_at_pc_face_level    ,&
                                                  scalar2_at_pc_face_level    ,&
                                                  scalar_eta_mass_velocity_at_pc_full,&
                                                  order                       ,&
                                                  scalar1_at_pc_full_level    ,&
                                                  scalar2_at_pc_full_level)
!
    use omp_lib
    type(global_domain),     intent(in)     :: mesh
    real(r4), allocatable,   intent(in)     :: scalar1_at_pc_face_level(:,:)
    real(r4), allocatable,   intent(in)     :: scalar2_at_pc_face_level(:,:)
#ifdef DBL_MASS
    real(r8), allocatable,   intent(in)     :: scalar_eta_mass_velocity_at_pc_full(:,:)
#else
    real(r4), allocatable,   intent(in)     :: scalar_eta_mass_velocity_at_pc_full(:,:)
#endif
    integer               ,  intent(in)     :: order
    real(r4), allocatable,   intent(inout)  :: scalar1_at_pc_full_level(:,:)
    real(r4), allocatable,   intent(inout)  :: scalar2_at_pc_full_level(:,:)
! local
    integer(i4)    :: ilev, iv, ii
    real(r4)       :: part1, der_k, der_kp1
    

!$omp parallel  private(iv) 
!$omp do schedule(static,5) 
     do iv = 1, mesh%nv_halo(1)
        scalar1_at_pc_full_level(1,iv)      = half_r4*(scalar1_at_pc_face_level(1,iv)   +scalar1_at_pc_face_level(2,iv))
        scalar1_at_pc_full_level(nlev,iv)   = half_r4*(scalar1_at_pc_face_level(nlev,iv)+scalar1_at_pc_face_level(nlev+1,iv))
        scalar2_at_pc_full_level(1,iv)      = half_r4*(scalar2_at_pc_face_level(1,iv)   +scalar2_at_pc_face_level(2,iv))
        scalar2_at_pc_full_level(nlev,iv)   = half_r4*(scalar2_at_pc_face_level(nlev,iv)+scalar2_at_pc_face_level(nlev+1,iv))
     end do
!$omp end do nowait
!$omp end parallel

    IF(EQS_VERT_DIFF)THEN
    select case(order)
   !  case(2)
   !     do iv = 1, mesh%nv_halo(1)
   !       do ilev = 2, nlev-1
   !          scalar1_at_pc_full_level(ilev,iv) = half*(scalar1_at_pc_face_level(ilev,iv)+scalar1_at_pc_face_level(ilev+1,iv))
   !          scalar2_at_pc_full_level(ilev,iv) = half*(scalar2_at_pc_face_level(ilev,iv)+scalar2_at_pc_face_level(ilev+1,iv))
   !       end do
   !     end do

    case(3) 
!$omp parallel  private(ii,iv,ilev)      
!$omp do schedule(dynamic,100) 
    do ii=1,mesh%nv_halo(1)*(nlev-2),1
        iv=ceiling(ii/real(nlev-2,r8))
        ilev=1+ii-(iv-1)*(nlev-2)
!       do iv = 1, mesh%nv_halo(1)
!         do ilev = 2, nlev-1
            scalar1_at_pc_full_level(ilev,iv) = (scalar1_at_pc_face_level(ilev,iv)  +scalar1_at_pc_face_level(ilev+1,iv))*(7._r4/12._r4)-&
                                                (scalar1_at_pc_face_level(ilev+2,iv)+scalar1_at_pc_face_level(ilev-1,iv))*(1._r4/12._r4)+&
#ifdef DBL_MASS
                                 sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*&
#else
                                 sign(1._r4,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*&
#endif
                                          ((scalar1_at_pc_face_level(ilev+2,iv)-scalar1_at_pc_face_level(ilev-1,iv))-&
                                     3._r4*(scalar1_at_pc_face_level(ilev+1,iv)-scalar1_at_pc_face_level(ilev,iv)))/12._r4

            scalar2_at_pc_full_level(ilev,iv) = (scalar2_at_pc_face_level(ilev,iv)  +scalar2_at_pc_face_level(ilev+1,iv))*(7._r4/12._r4)-&
                                                (scalar2_at_pc_face_level(ilev+2,iv)+scalar2_at_pc_face_level(ilev-1,iv))*(1._r4/12._r4)+&
#ifdef DBL_MASS
                                 sign(1._r8,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*&
#else
                                 sign(1._r4,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*&
#endif
                                          ((scalar2_at_pc_face_level(ilev+2,iv)-scalar2_at_pc_face_level(ilev-1,iv))-&
                                     3._r4*(scalar2_at_pc_face_level(ilev+1,iv)-scalar2_at_pc_face_level(ilev,iv)))/12._r4

!         end do
       !end do
    end do
!$omp end do nowait
!$omp end parallel 
     case default
        print*," you must select a vert order in calc_hpe_vert_flux_operator_face"
        stop
     end select

!      ELSE

!      select case(order)
!    !   case(2) 
!    !     do iv = 1, mesh%nv_halo(1)
!    !       do ilev = 2, nlev-1
!    !          scalar1_at_pc_full_level(ilev,iv) = half*(scalar1_at_pc_face_level(ilev,iv)+scalar1_at_pc_face_level(ilev+1,iv))
!    !          scalar2_at_pc_full_level(ilev,iv) = half*(scalar2_at_pc_face_level(ilev,iv)+scalar2_at_pc_face_level(ilev+1,iv))
!    !       end do
!    !     end do
!      case(3)
!       do iv = 1, mesh%nv_halo(1)
!          do ilev = 2, nlev-1
!             ! var1
!             part1 = half_r4*(scalar1_at_pc_face_level(ilev,iv)+scalar1_at_pc_face_level(ilev+1,iv))

!             der_k = 2._r4*(deta_face(ilev)**2)*((scalar1_at_pc_face_level(ilev+1,iv)-scalar1_at_pc_face_level(ilev,iv))/deta_full(ilev)-&
!                                                 (scalar1_at_pc_face_level(ilev,iv)-scalar1_at_pc_face_level(ilev-1,iv))/deta_full(ilev-1))/&
!                                                 (deta_full(ilev)+deta_full(ilev-1))

!             der_kp1 = 2._r4*(deta_face(ilev+1)**2)*((scalar1_at_pc_face_level(ilev+2,iv)-scalar1_at_pc_face_level(ilev+1,iv))/deta_full(ilev+1)-&
!                                                     (scalar1_at_pc_face_level(ilev+1,iv)-scalar1_at_pc_face_level(ilev,iv))/deta_full(ilev))/&
!                                                     (deta_full(ilev+1)+deta_full(ilev))

!             scalar1_at_pc_full_level(ilev,iv)= part1-(der_k+der_kp1)/12._r4+sign(1._r4,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*(der_kp1-der_k)
         
!             ! var2
!             part1 = half_r4*(scalar2_at_pc_face_level(ilev,iv)+scalar2_at_pc_face_level(ilev+1,iv))

!             der_k = 2._r8*(deta_face(ilev)**2)*((scalar2_at_pc_face_level(ilev+1,iv)-scalar2_at_pc_face_level(ilev,iv))/deta_full(ilev)-&
!                                                 (scalar2_at_pc_face_level(ilev,iv)  -scalar2_at_pc_face_level(ilev-1,iv))/deta_full(ilev-1))/&
!                                                 (deta_full(ilev)+deta_full(ilev-1))

!             der_kp1 = 2._r4*(deta_face(ilev+1)**2)*((scalar2_at_pc_face_level(ilev+2,iv)-scalar2_at_pc_face_level(ilev+1,iv))/deta_full(ilev+1)-&
!                                                     (scalar2_at_pc_face_level(ilev+1,iv)-scalar2_at_pc_face_level(ilev,iv))  /deta_full(ilev))/&
!                                                     (deta_full(ilev+1)+deta_full(ilev))

!             scalar2_at_pc_full_level(ilev,iv)= part1-(der_k+der_kp1)/12._r4+sign(1._r4,scalar_eta_mass_velocity_at_pc_full(ilev,iv))*(der_kp1-der_k)

!          end do
!       end do
!      case default
!         print*," you must select a vert order in calc_hpe_vert_flux_operator_face_2d"
!         stop
!      end select

     END IF
     
     return
   end subroutine calc_hpe_vert_flux_operator_face_2d_var2

   subroutine grist_nh_et_init(mesh)
! io
     type(global_domain),  intent(in) :: mesh
! local
      if(.not.allocated(scalar_template_1d_nlevp_c)) allocate(scalar_template_1d_nlevp_c(nlev+1))

     return
   end subroutine grist_nh_et_init
   
   subroutine grist_nh_et_final

      deallocate(scalar_template_1d_nlevp_c)

   end subroutine grist_nh_et_final

end module grist_nh_explicit_tend_module_2d_mixed
