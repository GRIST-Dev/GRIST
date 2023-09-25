
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: The model's explicit diffusion module
! Revision history:
!           1. For second-order diffusion, SMG and Comp only runs either in 
!           one runtime; In real-run, SMG is enough for now
!           2. No fourth-order or higher-order hyperviscorsity is used for now
!           3. vertical direction also needs no diffusion for now
!           4. DCMIP-SC will choose to use constant coef
!           5. 2020, add velocity Laplacian based on vector intentity (vi)
!           6. 2020, hyperdiffusion is added for hori wind optional, 5&6 
!              par-seq consistency checked ok
!           7. in old terms, nsponge_layer denotes layer for adding diffusion, not
!             really "sponge-layer"; now as we added real sponge layer,
!             this should be called as nexpdif_layer. NSPONGE is only added for TOP-level 
!             normal velocity currently
!           8. use_laplacian_vi means use vector identity for computiny velocity
!              Laplacian, this is less tested and not recommenened, only retained
!              for testing purpose (produce worse supercell shape).
!           9. unify scalar_laplacian_2nd at full and face, by inputing K2 coef
!              and target variable levels (Regression OK); www can choose to use
!              SPLIKE K2 coef
!----------------------------------------------------------------------------

 module grist_dycore_diffusion_module

  use grist_constants,        only: i4, r8, rearth, zero, one
  use grist_domain_types,     only: global_domain
#ifndef SEQ_GRIST
  use grist_data_types,       only: scalar_2d_field, exchange_field_list_2d
#else
  use grist_data_types,       only: scalar_2d_field
#endif
  use grist_math_module,      only: convert_vector_cart2sph, cross_product, norm, &
                                    extrapolate_bdy
  use grist_nml_module,       only: nlev, nlevp, nsponge, spcoef, smg_coef, ref_leng, ko4_coef, ko6_coef, do_div2d_o4_damping, div2d_o4_damping_coef, &
                                    do_dycore_laplacian_2nd, do_dycore_laplacian_4th, do_dycore_laplacian_6th, use_laplacian_vi, &
                                    www_k2coef_type, nh_dynamics
  use grist_nml_module,       only: dtime=>model_timestep, vr_mode,use_www_hyperDiffusion
  use grist_mesh_weight_icosh,only: project_sphere2tplane
  use grist_hpe_constants,    only: deta_full, deta_face, eta_full, eta_face
#ifndef SEQ_GRIST
  use grist_config_partition, only: debug_data_1d, debug_data_2d
  use grist_config_partition, only: exchange_data_2d_add, exchange_data_2d
#endif
! data
  use grist_dycore_vars_module,           only: dycoreVarEdgeFull, dycoreVarCellFull, dycoreVarCellFace
  use grist_tracer_transport_vars_module, only: tracerVarCellFull
! operator
  use grist_dycore_gcd_recon_module_2d,   only: divergence_operator_2d, curl_operator_2d 
  use grist_mpi
#ifndef SEQ_GRIST
  use grist_clocks,                       only: clock_id, clock_begin, clock_end
#endif

  implicit none

  private
  public :: grist_diffusion_init          , &
            grist_diffusion_destruct      , &
            grist_dycore_diffusion_run    , &
            grist_tracer_diffusion_run    , &
            perot_weight_at_pc            , &
            perot_weight_at_dc
!
! internal data
!
     real(r8), allocatable  :: perot_weight_at_pc(:,:,:) ! 3, max nnb, nv
     real(r8), allocatable  :: perot_weight_at_dc(:,:,:) ! 3, max nnb, nt
     real(r8), allocatable  :: scalar_smg_eddy_coef_at_edge_full_level(:,:)  !nlev,  ne
     real(r8), allocatable  :: scalar_smg_eddy_coef_at_edge_face_level(:,:)  !nlevp, ne
     real(r8), allocatable  :: scalar_hwind_laplace_at_edge_full_level(:,:)  !nlev,  ne, only grad2(), no coef and sign
     type(scalar_2d_field)  :: scalar_hwind_laplace4th_at_edge_full_level    !nlev,  ne, only grad4(), no coef and sign, needs exchange
     integer(i4)            :: clock_diffusion

  CONTAINS

!=======================================================================================
! Major diffusion driver, once called, providing ptend_rk style tendency based on the 
! present model state, added like a ptend_rk physics
! Diffusion options: (all on model level)
! 1. NormalWind: 2nd SMG; 4th/6th-order hyperdiffusion; 2nd top-level sponge
! 2. potential temp: 2nd SMG
! 3. www: 2nd SMG or K2 SPLIKE const (controlled by spcoef)
! 4. tracer: 2nd SMG
! NOTE: SMG coef will opt to const coef in DCMIP-SC tests
!====================================================================

   subroutine grist_dycore_diffusion_run(mesh,nexpdif_level,ntracer)
! io
     type(global_domain), intent(inout) :: mesh
     integer(i4),         intent(in)    :: nexpdif_level
     integer(i4),         intent(in)    :: ntracer
! local
     integer(i4)                        :: itracer
#ifndef SEQ_GRIST
     type(exchange_field_list_2d),pointer :: field_head_2d
#endif
     real(r8)        :: www_k2coef(nexpdif_level+1,mesh%ne_halo(1))
     real(r8)        :: div2d_o4_damping_tendency(nlev,mesh%ne_halo(2))
     integer(i4)     :: ie

#ifndef SEQ_GRIST
      field_head_2d=>null()
      clock_diffusion = clock_id('diffusion')
      call clock_begin(clock_diffusion)
#endif
!
! for 2nd laplacian and smg coef
!

      if(do_dycore_laplacian_2nd)then ! smg for velocities and PT
! coef and vector Laplacian
       call calc_smg_eddy_coef(mesh,dycoreVarEdgeFull%scalar_normal_velocity_n%f)
! un
       call calc_hrwind_laplacian_2nd_full(mesh,dycoreVarEdgeFull%scalar_normal_velocity_n%f,&
                                                dycoreVarEdgeFull%tend_hwind_laplacian_2nd%f,&
                                                nexpdif_level)
! pt
       call calc_scalar_laplacian_2nd_FullFace(mesh,scalar_smg_eddy_coef_at_edge_full_level      ,&
                                                    dycoreVarCellFull%scalar_potential_temp_n%f   ,&
                                                    dycoreVarCellFull%tend_mass_pt_laplacian_2nd%f,&
                                                    nexpdif_level)

       dycoreVarCellFull%tend_mass_pt_laplacian_2nd%f = dycoreVarCellFull%tend_mass_pt_laplacian_2nd%f*&
                                                        dycoreVarCellFull%scalar_delhp_n%f
!
! www diffusion, we may control the diffusion coef here
!
       select case(www_k2coef_type)
       case('SPLIKE') ! sponge-like
           do ie = 1, mesh%ne_halo(1)  ! global index
               www_k2coef(:,ie) = 0.5_r8*(mesh%vtxCellLeng(mesh%edt_v(1,ie))+mesh%vtxCellLeng(mesh%edt_v(2,ie)))*rearth*spcoef
           end do
       case('SMG')    ! Smagorinsky
           www_k2coef(1:nexpdif_level+1,1:mesh%ne_halo(1)) = scalar_smg_eddy_coef_at_edge_face_level(1:nexpdif_level+1,1:mesh%ne_halo(1))
       case default   ! SMG
           www_k2coef(1:nexpdif_level+1,1:mesh%ne_halo(1)) = scalar_smg_eddy_coef_at_edge_face_level(1:nexpdif_level+1,1:mesh%ne_halo(1))
       end select

#ifdef DYCORE_NHD
       if(nh_dynamics) call calc_scalar_laplacian_2nd_FullFace(mesh,www_k2coef                  ,&
                                                    dycoreVarCellFace%scalar_www_n%f            ,&
                                                    dycoreVarCellFace%tend_www_laplacian_2nd%f  ,&
                                                    nexpdif_level+1)
#endif
       end if

       if(do_dycore_laplacian_4th)then
! only un is diffused using 4th if needed
          call calc_hrwind_laplacian_4th_full(mesh,2, ko4_coef, real(mesh%edt_scale_4th,r8), .true.,& 
                                                   scalar_hwind_laplace_at_edge_full_level,&
                                                   dycoreVarEdgeFull%tend_hwind_laplacian_4th%f)

! the VI version now support VR at 20210612, needs test
          if(use_laplacian_vi) call calc_hrwind_laplacian_4th_full_vi(mesh,2, ko4_coef, real(mesh%edt_scale_4th,r8), .true.,&
                                                   scalar_hwind_laplace_at_edge_full_level,&
                                                   dycoreVarEdgeFull%tend_hwind_laplacian_4th%f)
#ifdef DYCORE_NHD
!
! 20210827, add Laplacian 4th for phi, with sign "minus"
! 20220923, change Laplacian4th for www instead of phi
!
          if(nh_dynamics.and.use_www_hyperDiffusion)then
             call calc_scalar_laplacian_4th_face(mesh, dycoreVarCellFace%scalar_www_n%f          , &
                                                       dycoreVarCellFace%tend_www_laplacian_4th%f, &
                                                       nexpdif_level+1)
          end if
#endif

       end if

! if use div2d damping, add its tendency into laplacian_4th, so that only needs to compute it until ne_compute,
! and can share MPI pattern of hori_diffusion 4th

       if(do_div2d_o4_damping)then
          call calc_hori_div2d_o4_damping(mesh, nlev, div2d_o4_damping_coef, &
                                          dycoreVarEdgeFull%scalar_normal_velocity_n%f, &
                                          div2d_o4_damping_tendency)

          dycoreVarEdgeFull%tend_hwind_laplacian_4th%f(1:nlev,1:mesh%ne_compute) = dycoreVarEdgeFull%tend_hwind_laplacian_4th%f(1:nlev,1:mesh%ne_compute)+&
                                         div2d_o4_damping_tendency(1:nlev,1:mesh%ne_compute)
       end if

#ifdef GRIST_DYCORE_TEST
       if(do_dycore_laplacian_6th)then
!LAM does not support 6th even for testing purpose!
#ifndef SEQ_GRIST
! exchange data, because the computation of 4th ends at ne_compute/ne
          call exchange_data_2d_add(mesh,field_head_2d,scalar_hwind_laplace4th_at_edge_full_level)
          call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
! do computation now, not VI version
          call calc_hrwind_laplacian_4th_full(mesh,3, ko6_coef, real(mesh%edt_scale_6th,r8), .false., &
                                                   scalar_hwind_laplace4th_at_edge_full_level%f,&
                                                   dycoreVarEdgeFull%tend_hwind_laplacian_6th%f)

       end if
#endif

#ifndef SEQ_GRIST
      call clock_end(clock_diffusion)
#endif

      return
   end subroutine grist_dycore_diffusion_run
 
   subroutine grist_tracer_diffusion_run(mesh,nexpdif_level,ntracer)
! io
     type(global_domain), intent(in) :: mesh
     integer(i4),         intent(in) :: nexpdif_level
     integer(i4),         intent(in) :: ntracer
     integer(i4)        :: itracer

! tracer
       call calc_scalar3d_laplacian_2nd_full(mesh,dycoreVarEdgeFull%scalar_normal_velocity_n%f,&
                                                  tracerVarCellFull%scalar_tracer_mxrt_n%f      ,&
                                                  tracerVarCellFull%tend_tracer_mass_laplacian_2nd%f,&
                                                  nexpdif_level,ntracer)
       do itracer = 1, ntracer
          tracerVarCellFull%tend_tracer_mass_laplacian_2nd%f(itracer,:,:) = tracerVarCellFull%tend_tracer_mass_laplacian_2nd%f(itracer,:,:)*&
                                                                            dycoreVarCellFull%scalar_delhp_n%f
       end do

      return
   end subroutine grist_tracer_diffusion_run

!================================================
! init constant for diffusion calculation
!================================================

    subroutine grist_diffusion_init(mesh)
! io
     type(global_domain), intent(in) :: mesh
! local
     integer(i4)       :: iv, ie, it, inb
     real(r8)          :: local_basis(3),tmp,vector_tmp(3)

        if(.not.allocated(   perot_weight_at_pc)) allocate(   perot_weight_at_pc(3,8,mesh%nv_halo(2)))  ! 8 means max number of nnb 
        if(.not.allocated(   perot_weight_at_dc)) allocate(   perot_weight_at_dc(3,3,mesh%nt_halo(1)))
        if(.not.allocated(scalar_smg_eddy_coef_at_edge_full_level)) allocate(scalar_smg_eddy_coef_at_edge_full_level(nlev,mesh%ne_halo(1)))
        if(.not.allocated(scalar_smg_eddy_coef_at_edge_face_level)) allocate(scalar_smg_eddy_coef_at_edge_face_level(nlevp,mesh%ne_halo(1)))
        if(.not.allocated(scalar_hwind_laplace_at_edge_full_level)) allocate(scalar_hwind_laplace_at_edge_full_level(nlev,mesh%ne_halo(1)))
        if(do_dycore_laplacian_6th.and.(.not.allocated(scalar_hwind_laplace4th_at_edge_full_level%f))) allocate(scalar_hwind_laplace4th_at_edge_full_level%f(nlev,mesh%ne_full))
!
! reconstruction weight used for reconstruct velocity vector at primal&dual cells 
! using perot method, this weight is time-invariant, double check
!
        do iv  = 1, mesh%nv_halo(2)
           do inb  = 1, mesh%vtx_nnb(iv)
              ie = mesh%vtx_ed(inb,iv)
              call project_sphere2tplane(real(mesh%edt_c_p(1:3,ie),r8), real(mesh%vtx_p(1:3,iv),r8),local_basis) ! unit sphere
              perot_weight_at_pc(:,inb,iv) = mesh%edp_leng(ie)*dot_product(local_basis/norm(local_basis),mesh%edp_nr(1:3,ie))*local_basis/mesh%plg_areag(iv)
           end do
        end do

        do it  = 1, mesh%nt_halo(1)
           do inb  = 1, mesh%tri_nnb(it)
              ie = mesh%tri_ed(inb,it)
              call project_sphere2tplane(real(mesh%edt_c_p(1:3,ie),r8),real(mesh%tri_c_p(1:3,it),r8),local_basis) ! unit sphere
              tmp       = dot_product(mesh%edp_nr(1:3,ie),cross_product(real(mesh%edt_c_p(1:3,ie),r8),local_basis)/norm(cross_product(real(mesh%edt_c_p(1:3,ie),r8),local_basis)))
              vector_tmp= cross_product(real(mesh%edt_c_p(1:3,ie),r8),local_basis)
              perot_weight_at_dc(:,inb,it) = vector_tmp*tmp*mesh%edt_leng(ie)/mesh%tri_areag(it)
           end do
        end do

      return
    end subroutine grist_diffusion_init

    subroutine calc_scalar_laplacian_4th_face(mesh,scalar_at_pc_face_level,             &
                                                   tend_laplacian_4th_at_pc_face_level, &
                                                   nLevel)
! io
      type(global_domain),   intent(in)    :: mesh
      real(r8), allocatable, intent(in)    :: scalar_at_pc_face_level(:,:)
      real(r8), allocatable, intent(inout) :: tend_laplacian_4th_at_pc_face_level(:,:)
      integer(i4)          , intent(in)    :: nLevel
! local
      real(r8)              :: gradient_at_prime_edge(nLevel,mesh%ne_halo(1))
      integer(i4)           :: ie, iv, inb, ilev, v1, v2, order
      real(r8)              :: v1v2(3), flag, div_sum, expp

!
! 1st pass: gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1)
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           do ilev = 1, nLevel
              gradient_at_prime_edge(ilev,ie) = flag*(scalar_at_pc_face_level(ilev,v2)-scalar_at_pc_face_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
           end do
        end do
!
! 1st pass: divergence at cell
!
        do iv = 1, mesh%nv_halo(1)
           tend_laplacian_4th_at_pc_face_level(:,iv) = zero
           do ilev = 1, nLevel
              div_sum = zero
              do inb = 1, mesh%vtx_nnb(iv)
                 ie       = mesh%vtx_ed(inb,iv)
                 div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(ie)
              end do
              tend_laplacian_4th_at_pc_face_level(ilev,iv) = div_sum/((rearth**2)*mesh%plg_areag(iv))
           end do
        end do
!
! 2nd pass: gradient at edge, counter edge's normal direction, untill ne_compute
!
        if(mesh%ne.ne.mesh%ne_compute)then
           if(mpi_rank().eq.0) print*,"ne .not equal ne_compute in calc_scalar_laplacian_4th_face"
        end if
        if(mesh%nv.ne.mesh%nv_compute)then
           if(mpi_rank().eq.0) print*,"nv .not equal nv_compute in calc_scalar_laplacian_4th_face"
        end if

        do ie = 1, mesh%ne_compute
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           do ilev = 1, nLevel
              gradient_at_prime_edge(ilev,ie) = flag*(tend_laplacian_4th_at_pc_face_level(ilev,v2)-tend_laplacian_4th_at_pc_face_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
           end do
        end do
!
! 2nd pass: divergence at cell
!
        expp  = 3.3219
        order = 2

        IF(vr_mode)THEN
        do iv = 1, mesh%nv_compute
           tend_laplacian_4th_at_pc_face_level(:,iv) = zero
           do ilev = 1, nLevel
              div_sum = zero
              do inb = 1, mesh%vtx_nnb(iv)
                 ie  = mesh%vtx_ed(inb,iv)
                 div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(ie)
              end do
! multiplied coef here with scaling factor
              tend_laplacian_4th_at_pc_face_level(ilev,iv) = -ko4_coef*div_sum/((rearth**2)*mesh%plg_areag(iv))*(mesh%vtxCellLeng(iv)/(ref_leng/rearth))**expp 
           end do
        end do
   
        ELSE

        do iv = 1, mesh%nv_compute
           tend_laplacian_4th_at_pc_face_level(:,iv) = zero
           do ilev = 1, nLevel
              div_sum = zero
              do inb = 1, mesh%vtx_nnb(iv)
                 ie  = mesh%vtx_ed(inb,iv)
                 div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(ie)
              end do
! multiplied coef here with scaling factor
              tend_laplacian_4th_at_pc_face_level(ilev,iv) = -ko4_coef*div_sum/((rearth**2)*mesh%plg_areag(iv))
           end do
        end do

        END IF

      return
    end subroutine calc_scalar_laplacian_4th_face

    subroutine calc_scalar_laplacian_2nd_face(mesh,scalar_at_pc_face_level, &
                                                   tend_laplacian_2nd_at_pc_face_level, &
                                                   nLevel)
! io
      type(global_domain),   intent(in)    :: mesh
      real(r8), allocatable, intent(in)    :: scalar_at_pc_face_level(:,:)
      real(r8), allocatable, intent(inout) :: tend_laplacian_2nd_at_pc_face_level(:,:)
      integer(i4)          , intent(in)    :: nLevel
! local
      real(r8), allocatable  :: gradient_at_prime_edge(:,:)
      integer(i4)            :: ie, iv, inb,ilev, v1, v2
      real(r8)               :: v1v2(3), flag, div_sum

        if(.not.allocated(gradient_at_prime_edge)) allocate(gradient_at_prime_edge(nLevel,mesh%ne_halo(1)))
!
! gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1) ! global index
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(1._r8,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           do ilev = 1, nLevel
              gradient_at_prime_edge(ilev,ie) = flag*(scalar_at_pc_face_level(ilev,v2)-scalar_at_pc_face_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
              !gradient_at_prime_edge(ilev,ie) = scalar_smg_eddy_coef_at_edge_face_level(ilev,ie)*gradient_at_prime_edge(ilev,ie)
! tests this for www in real-world
              gradient_at_prime_edge(ilev,ie) = 0.5_r8*(mesh%vtxCellLeng(v1)+mesh%vtxCellLeng(v2))*rearth*spcoef*gradient_at_prime_edge(ilev,ie)
           end do
        end do
!
! divergence at cell
!
        do iv = 1, mesh%nv_halo(1)
           tend_laplacian_2nd_at_pc_face_level(:,iv) = 0._r8
           do ilev = 1, nLevel
              div_sum = 0._r8
              do inb = 1, mesh%vtx_nnb(iv)
                 ie  = mesh%vtx_ed(inb,iv)
                 !div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg(iv)%nr(inb)*rearth*mesh%edp(ie)%leng
                 div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(ie)
              end do
              tend_laplacian_2nd_at_pc_face_level(ilev,iv) = div_sum/((rearth**2)*mesh%plg_areag(iv))
           end do
        end do

        if(allocated(gradient_at_prime_edge)) deallocate(gradient_at_prime_edge)

      return
    end subroutine calc_scalar_laplacian_2nd_face

! default 2d data scalar

    subroutine calc_scalar_laplacian_2nd_full(mesh,scalar_normal_velocity_at_edge_full_level, &
                                                   scalar_at_pc_full_level                  , &
                                                   tend_laplacian_2nd_at_pc_full_level,  &
                                                   nLevel)
! io
      type(global_domain),   intent(in)    :: mesh
      real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
      real(r8), allocatable, intent(in)    :: scalar_at_pc_full_level(:,:)
      real(r8), allocatable, intent(inout) :: tend_laplacian_2nd_at_pc_full_level(:,:)
      integer(i4)          , intent(in)    :: nLevel
! local
      real(r8), allocatable  :: gradient_at_prime_edge(:,:)
      integer(i4)            :: ie, iv, inb,ilev, v1, v2
      real(r8)               :: v1v2(3), flag, div_sum

        if(.not.allocated(gradient_at_prime_edge)) allocate(gradient_at_prime_edge(nLevel,mesh%ne_halo(1)))
!
! gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1)     ! global index
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(1._r8,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           do ilev = 1, nLevel
              gradient_at_prime_edge(ilev,ie) = flag*(scalar_at_pc_full_level(ilev,v2)-scalar_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
              gradient_at_prime_edge(ilev,ie) = scalar_smg_eddy_coef_at_edge_full_level(ilev,ie)*gradient_at_prime_edge(ilev,ie)
           end do
        end do
!
! divergence at cell
!
        do iv = 1, mesh%nv_halo(1)
           tend_laplacian_2nd_at_pc_full_level(:,iv) = 0._r8 ! default zero
           do ilev = 1, nLevel
              div_sum = 0._r8
              do inb = 1, mesh%vtx_nnb(iv)
                 ie  = mesh%vtx_ed(inb,iv)
                 !div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg(iv)%nr(inb)*rearth*mesh%edp(ie)%leng
                 div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(ie)
              end do
              tend_laplacian_2nd_at_pc_full_level(ilev,iv) = div_sum/((rearth**2)*mesh%plg_areag(iv))
           end do
        end do

        if(allocated(gradient_at_prime_edge)) deallocate(gradient_at_prime_edge)

        return
    end subroutine calc_scalar_laplacian_2nd_full

!
! This routine generalizes the two above, for both face and full levels,
! but only for SMG-style 2nd Laplacian diffusion because smg-coef is inserted
!

    subroutine calc_scalar_laplacian_2nd_FullFace(mesh,scalar_smg_eddy_coef_at_edge  , &
                                                       scalar_at_pc                  , &
                                                       tend_laplacian_2nd_at_pc      , &
                                                       nLevel)
! io
      type(global_domain),   intent(in)    :: mesh
      real(r8),              intent(in)    :: scalar_smg_eddy_coef_at_edge(1:nLevel,1:mesh%ne_halo(1))
      real(r8), allocatable, intent(in)    :: scalar_at_pc(:,:)
      real(r8), allocatable, intent(inout) :: tend_laplacian_2nd_at_pc(:,:)
      integer(i4)          , intent(in)    :: nLevel
! local
      real(r8)               :: gradient_at_prime_edge(nLevel,mesh%ne_halo(1))
      integer(i4)            :: ie, iv, inb,ilev, v1, v2
      real(r8)               :: v1v2(3), flag, div_sum

!
! gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1)     ! global index
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(1._r8,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           do ilev = 1, nLevel
              gradient_at_prime_edge(ilev,ie) = flag*(scalar_at_pc(ilev,v2)-scalar_at_pc(ilev,v1))/(rearth*mesh%edt_leng(ie))
              gradient_at_prime_edge(ilev,ie) = scalar_smg_eddy_coef_at_edge(ilev,ie)*gradient_at_prime_edge(ilev,ie)
           end do
        end do
!
! divergence at cell
!
        do iv = 1, mesh%nv_halo(1)
           tend_laplacian_2nd_at_pc(:,iv) = zero ! default zero
           do ilev = 1, nLevel
              div_sum = zero
              do inb = 1, mesh%vtx_nnb(iv)
                 ie  = mesh%vtx_ed(inb,iv)
                 !div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg(iv)%nr(inb)*rearth*mesh%edp(ie)%leng
                 div_sum  = div_sum+gradient_at_prime_edge(ilev,ie)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(ie)
              end do
              tend_laplacian_2nd_at_pc(ilev,iv) = div_sum/((rearth**2)*mesh%plg_areag(iv))
           end do
        end do

        return
    end subroutine calc_scalar_laplacian_2nd_FullFace

    subroutine calc_scalar3d_laplacian_2nd_full(mesh,scalar_normal_velocity_at_edge_full_level, &
                                                     scalar_at_pc_full_level                  , &
                                                     tend_laplacian_2nd_at_pc_full_level      , &
                                                     nexpdif, ntracer)
! io
      type(global_domain),   intent(in)    :: mesh
      real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
      real(r8), allocatable, intent(in)    :: scalar_at_pc_full_level(:,:,:)
      real(r8), allocatable, intent(inout) :: tend_laplacian_2nd_at_pc_full_level(:,:,:)
      integer(i4)          , intent(in)    :: nexpdif
      integer(i4)          , intent(in)    :: ntracer
! local
      real(r8), allocatable  :: gradient_at_prime_edge(:,:,:)
      integer(i4)            :: ie, iv, inb,ilev, v1, v2, itracer
      real(r8)               :: v1v2(3), flag, div_sum

        if(.not.allocated(gradient_at_prime_edge)) allocate(gradient_at_prime_edge(ntracer,nexpdif,mesh%ne_halo(1)))
!
! gradient at edge, counter edge's normal direction
!
        do ie = 1, mesh%ne_halo(1)     ! global index
           v1      = mesh%edt_v(1,ie)
           v2      = mesh%edt_v(2,ie)
           v1v2    = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
           flag    = sign(1._r8,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
           do ilev = 1, nexpdif
              do itracer = 1, ntracer
                 gradient_at_prime_edge(itracer,ilev,ie) = flag*(scalar_at_pc_full_level(itracer,ilev,v2)-scalar_at_pc_full_level(itracer,ilev,v1))/(rearth*mesh%edt_leng(ie))
                 gradient_at_prime_edge(itracer,ilev,ie) = scalar_smg_eddy_coef_at_edge_full_level(ilev,ie)*gradient_at_prime_edge(itracer,ilev,ie)
              end do
           end do
        end do
!
! divergence at cell
!
        do iv = 1, mesh%nv_halo(1)
           tend_laplacian_2nd_at_pc_full_level(:,:,iv) = 0._r8 ! default zero
! need optimization of the loop sequence in future
           do ilev = 1, nexpdif
              do itracer = 1, ntracer
                 div_sum = 0._r8
                 do inb = 1, mesh%vtx_nnb(iv)
                    ie  = mesh%vtx_ed(inb,iv)
                    !div_sum  = div_sum+gradient_at_prime_edge(itracer,ilev,ie)*mesh%plg(iv)%nr(inb)*rearth*mesh%edp(ie)%leng
                    div_sum  = div_sum+gradient_at_prime_edge(itracer,ilev,ie)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(ie)
                 end do
                 tend_laplacian_2nd_at_pc_full_level(itracer,ilev,iv) = div_sum/((rearth**2)*mesh%plg_areag(iv))
              end do
           end do
        end do

        if(allocated(gradient_at_prime_edge)) deallocate(gradient_at_prime_edge)

        return
    end subroutine calc_scalar3d_laplacian_2nd_full

    subroutine calc_hrwind_laplacian_2nd_full(mesh,scalar_normal_velocity_at_edge_full_level, &
                                                   tend_laplacian_2nd_at_pc_full_level, &
                                                   nexpdif)
! io
      type(global_domain),   intent(inout) :: mesh
      real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
      real(r8), allocatable, intent(inout) :: tend_laplacian_2nd_at_pc_full_level(:,:)
      integer(i4)          , intent(in)    :: nexpdif
! local
      integer(i4)            :: ie, ilev,cell1, cell2
!
! if use laplacian vi, will overwrite that from smg_coef
!
      if(use_laplacian_vi) call calc_hrwind_laplacian_2nd_full_vi(mesh,scalar_normal_velocity_at_edge_full_level)
  
      do ie = 1, mesh%ne_compute
         tend_laplacian_2nd_at_pc_full_level(:,ie) = 0._r8  ! default zero
         do ilev = 1, nexpdif
             tend_laplacian_2nd_at_pc_full_level(ilev,ie) = 4._r8*scalar_smg_eddy_coef_at_edge_full_level(ilev,ie)*&
                                                                  scalar_hwind_laplace_at_edge_full_level(ilev,ie)
#ifdef DCMIP_SUPERCELL
             tend_laplacian_2nd_at_pc_full_level(ilev,ie) = 4._r8*500._r8*&
                                                                  scalar_hwind_laplace_at_edge_full_level(ilev,ie)
#endif
         end do
      end do
!
! add top-nsponge layer with real 2nd diffusion
!
      if(nsponge.gt.0)then ! default should be 5
        do ie = 1, mesh%ne_compute
           cell1 = mesh%edt_v(1,ie)
           cell2 = mesh%edt_v(2,ie)
! from top to bottom
           do ilev = 1, nsponge
               tend_laplacian_2nd_at_pc_full_level(ilev,ie) = tend_laplacian_2nd_at_pc_full_level(ilev,ie)+ &
                                            4._r8*scalar_hwind_laplace_at_edge_full_level(ilev,ie)*&
                                       0.5_r8*(mesh%vtxCellLeng(cell1)+mesh%vtxCellLeng(cell2))*rearth*(7-ilev)*spcoef/6._r8
           end do
        end do
      end if

      return
    end subroutine calc_hrwind_laplacian_2nd_full

    subroutine grist_diffusion_destruct

        if(allocated(perot_weight_at_pc))                           deallocate(perot_weight_at_pc)
        if(allocated(perot_weight_at_dc))                           deallocate(perot_weight_at_dc)
        if(allocated(scalar_smg_eddy_coef_at_edge_full_level))      deallocate(scalar_smg_eddy_coef_at_edge_full_level)
        if(allocated(scalar_smg_eddy_coef_at_edge_face_level))      deallocate(scalar_smg_eddy_coef_at_edge_face_level)
        if(allocated(scalar_hwind_laplace_at_edge_full_level))      deallocate(scalar_hwind_laplace_at_edge_full_level)
        if(allocated(scalar_hwind_laplace4th_at_edge_full_level%f)) deallocate(scalar_hwind_laplace4th_at_edge_full_level%f)

    end subroutine grist_diffusion_destruct
!
! This routine also compute grad2 of normal velocity because of efficiency to do them all
!
    subroutine calc_smg_eddy_coef(mesh, &
                                  scalar_normal_velocity_at_edge_full_level)
! io
     use omp_lib
     type(global_domain),   intent(in) :: mesh
     real(r8), allocatable, intent(in) :: scalar_normal_velocity_at_edge_full_level(:,:)
! local
     real(r8), allocatable             :: velocity_vector_at_pc_full_level(:,:,:)
     real(r8), allocatable             :: velocity_vector_at_dc_full_level(:,:,:)
     integer(i4)                       :: iv, it, ie, ilev, inb
     integer(i4)                       :: pc1, pc2, dc1, dc2              ! cell index
     real(r8)                          :: un_pc1, un_pc2, un_dc1, un_dc2  ! normal wind follow edp's nr
     real(r8)                          :: ut_pc1, ut_pc2, ut_dc1, ut_dc2  ! tangent wind follow edp's tg
     real(r8)                          :: pc1pc2(3), dc1dc2(3), flag1, flag2
     real(r8)                          :: un_pc1pc2, ut_pc1pc2, un_dc1dc2, ut_dc1dc2

! calculate velocity vector at primal and dual cells

      if(.not.allocated(velocity_vector_at_pc_full_level)) allocate(velocity_vector_at_pc_full_level(3,nlev,mesh%nv_halo(2))) ! 3 means 3 components of a cartesian vector
      if(.not.allocated(velocity_vector_at_dc_full_level)) allocate(velocity_vector_at_dc_full_level(3,nlev,mesh%nt_halo(1)))

!$omp parallel  private(iv,inb,ie,ilev)    
!$omp do schedule(dynamic,50) 
      do iv  = 1, mesh%nv_halo(2)
         velocity_vector_at_pc_full_level(:,:,iv) = 0._r8
         do inb  = 1, mesh%vtx_nnb(iv)
            ie = mesh%vtx_ed(inb,iv)
            do ilev = 1, nlev
                velocity_vector_at_pc_full_level(:,ilev,iv) = velocity_vector_at_pc_full_level(:,ilev,iv)+perot_weight_at_pc(:,inb,iv)*&
                                                   scalar_normal_velocity_at_edge_full_level(ilev,ie)
            end do
         end do
      end do
!$omp end do nowait
!$omp end parallel 

!$omp parallel  private(it,inb,ilev,ie) 
!$omp do schedule(dynamic,50) 
      do it  = 1, mesh%nt_halo(1)
         velocity_vector_at_dc_full_level(:,:,it) = 0._r8
         do inb = 1, mesh%tri_nnb(it)
            ie = mesh%tri_ed(inb,it)
            do ilev = 1, nlev
               velocity_vector_at_dc_full_level(:,ilev,it) = velocity_vector_at_dc_full_level(:,ilev,it)+perot_weight_at_dc(:,inb,it)*&
                                                   scalar_normal_velocity_at_edge_full_level(ilev,ie)
            end do
         end do
      end do
!$omp end do nowait
!$omp end parallel 

      !call debug_data_2d(1,1,mesh%v_index,mesh%nv_halo(1),"x1x", velocity_vector_at_pc_full_level(1,:,:))
      !call debug_data_2d(1,1,mesh%v_index,mesh%nv_halo(1),"x2x", perot_weight_at_pc(1,:,:))
      !call debug_data_2d(1,1,mesh%t_index,mesh%nt_halo(1),"y1y", velocity_vector_at_dc_full_level(1,:,:))
      !call debug_data_2d(1,1,mesh%t_index,mesh%nt_halo(1),"y2y", perot_weight_at_dc(1,:,:))

!-----------------------------------------------------------------------
!         pc2
!          |
!    dc1---e---dc2
!          |
!         pc1
!-----------------------------------------------------------------------

!$omp parallel  private(ie,pc1,pc2,dc1,dc2,pc1pc2,dc1dc2,flag1,flag2,ilev,un_pc1,un_pc2,un_dc1,un_dc2,ut_pc1,ut_pc2,ut_dc1,ut_dc2,un_pc1pc2,ut_pc1pc2,un_dc1dc2,ut_dc1dc2) 
!$omp do schedule(dynamic,20) 
      DO ie = 1, mesh%ne_halo(1)

         pc1     = mesh%edt_v(1,ie)
         pc2     = mesh%edt_v(2,ie)
         dc1     = mesh%edp_v(1,ie)
         dc2     = mesh%edp_v(2,ie)
         pc1pc2  = mesh%vtx_p(1:3,pc2)-mesh%vtx_p(1:3,pc1)
         dc1dc2  = mesh%tri_c_p(1:3,dc2)-mesh%tri_c_p(1:3,dc1)
         flag1   = sign(1._r8,dot_product(pc1pc2,real(mesh%edp_nr(1:3,ie),r8)))
         flag2   = sign(1._r8,dot_product(dc1dc2,real(mesh%edp_tg(1:3,ie),r8)))

         do ilev = 1, nlev
            un_pc1  = dot_product(velocity_vector_at_pc_full_level(:,ilev,pc1),mesh%edp_nr(1:3,ie))
            un_pc2  = dot_product(velocity_vector_at_pc_full_level(:,ilev,pc2),mesh%edp_nr(1:3,ie))
            un_dc1  = dot_product(velocity_vector_at_dc_full_level(:,ilev,dc1),mesh%edp_nr(1:3,ie))
            un_dc2  = dot_product(velocity_vector_at_dc_full_level(:,ilev,dc2),mesh%edp_nr(1:3,ie))
                                                                                              
            ut_pc1  = dot_product(velocity_vector_at_pc_full_level(:,ilev,pc1),mesh%edp_tg(1:3,ie))
            ut_pc2  = dot_product(velocity_vector_at_pc_full_level(:,ilev,pc2),mesh%edp_tg(1:3,ie))
            ut_dc1  = dot_product(velocity_vector_at_dc_full_level(:,ilev,dc1),mesh%edp_tg(1:3,ie))
            ut_dc2  = dot_product(velocity_vector_at_dc_full_level(:,ilev,dc2),mesh%edp_tg(1:3,ie))

            un_pc1pc2 = flag1*(un_pc2-un_pc1)/(rearth*mesh%edt_leng(ie))
            ut_pc1pc2 = flag1*(ut_pc2-ut_pc1)/(rearth*mesh%edt_leng(ie))
            un_dc1dc2 = flag2*(un_dc2-un_dc1)/(rearth*mesh%edp_leng(ie))
            ut_dc1dc2 = flag2*(ut_dc2-ut_dc1)/(rearth*mesh%edp_leng(ie))

            if(vr_mode)then
               scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*(rearth**2)*mesh%edt_leng(ie)*mesh%edp_leng(ie)*sqrt((un_pc1pc2-ut_dc1dc2)**2+(un_dc1dc2+ut_pc1pc2)**2)
            else ! QU
               scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*((rearth*mesh%mean_edt_dist)**2)*sqrt((un_pc1pc2-ut_dc1dc2)**2+(un_dc1dc2+ut_pc1pc2)**2)
            end if

            !----------------------------------------depreciated------------------------------------------------
!#if (defined VRVIS)
            !scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*(ref_leng**2)*mesh%edt_scale_2nd(ie)*&
            !scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*((rearth*mesh%min_edt_dist)**2)*mesh%edt_scale_2nd(ie)*&
            !scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*((rearth*mesh%edt_leng(ie))**2)*&
            !scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*(ref_leng**2)*mesh%edt_scale_2nd(ie)*&
            !scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*(rearth**2)*mesh%edt_leng(ie)*mesh%edp_leng(ie)*&
!#elif (defined VRSMG)
        ! old setup for exact regression
            !scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*((rearth*mesh%edt_leng(ie))**2)*&
!#else
            !scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = smg_coef*((rearth*mesh%mean_edt_dist)**2)*&
!#endif
            !----------------------------------------------------------------------------------------------------------

#ifdef DCMIP_SUPERCELL
            scalar_smg_eddy_coef_at_edge_full_level(ilev,ie) = 1500._r8
#endif
     ! second-order laplacian of normal wind, because we should use 1/2(edt&edp's length) for double center differences but not,
     ! so this laplacian needs a *4 when using outside, for laplacian vi, we also normalize it like this;
     ! why not *4 here? currently just for bit-bit regression!
     ! math: ((upc1-un)/(0.5*edt_leng)-(un-upc2)/(0.5*edt_leng))/(0.5*edt_leng)
            scalar_hwind_laplace_at_edge_full_level(ilev,ie) = (un_pc1+un_pc2-2._r8*scalar_normal_velocity_at_edge_full_level(ilev,ie))/((rearth*mesh%edt_leng(ie))**2)+&
                                                               (un_dc1+un_dc2-2._r8*scalar_normal_velocity_at_edge_full_level(ilev,ie))/((rearth*mesh%edp_leng(ie))**2)
         end do
!
! face level value is a combination of full level
! 
         do ilev = 2, nlev
             scalar_smg_eddy_coef_at_edge_face_level(ilev,ie) = &
                                    0.5*(deta_full(ilev-1)/deta_face(ilev)*scalar_smg_eddy_coef_at_edge_full_level(ilev,ie)+&
                                         deta_full(ilev)  /deta_face(ilev)*scalar_smg_eddy_coef_at_edge_full_level(ilev-1,ie))
         end do
!#ifndef REGRESSION
! bdy extrapolation
!         scalar_smg_eddy_coef_at_edge_face_level(1,ie)     = extrapolate_bdy(scalar_smg_eddy_coef_at_edge_full_level(1,ie),&
!                                                                             scalar_smg_eddy_coef_at_edge_full_level(2,ie),&
!                                                                             scalar_smg_eddy_coef_at_edge_full_level(3,ie),&
!                                                                             eta_full(1), eta_full(2), eta_full(3), eta_face(1))
 
!         scalar_smg_eddy_coef_at_edge_face_level(nlevp,ie)= extrapolate_bdy(scalar_smg_eddy_coef_at_edge_full_level(nlev-2,ie),&
!                                                                             scalar_smg_eddy_coef_at_edge_full_level(nlev-1,ie),&
!                                                                             scalar_smg_eddy_coef_at_edge_full_level(nlev,ie)  ,&
!                                                                             eta_full(nlev-2),eta_full(nlev-1), eta_full(nlev), eta_face(nlevp))
!#else
         scalar_smg_eddy_coef_at_edge_face_level(1,ie)     = 2.*scalar_smg_eddy_coef_at_edge_full_level(1,ie)-scalar_smg_eddy_coef_at_edge_face_level(2,ie)
         scalar_smg_eddy_coef_at_edge_face_level(nlevp,ie) = 2.*scalar_smg_eddy_coef_at_edge_full_level(nlev,ie)-scalar_smg_eddy_coef_at_edge_face_level(nlev,ie)
!#endif
#ifdef DCMIP_SUPERCELL
         scalar_smg_eddy_coef_at_edge_face_level(:,:) = 500._r8
#endif

       END DO
!$omp end do nowait
!$omp end parallel 

       if(allocated(velocity_vector_at_pc_full_level)) deallocate(velocity_vector_at_pc_full_level)
       if(allocated(velocity_vector_at_dc_full_level)) deallocate(velocity_vector_at_dc_full_level)

       return
    end subroutine calc_smg_eddy_coef

!
! This calc_hrwind_laplacian_4th_full subroutine is feasible 
! for 4th, 6th..., note: 2nd is done within SMG
!
    subroutine calc_hrwind_laplacian_4th_full(mesh, order, coef, scaling_factor, flag, &
                                              scalar_hwind_lap2nd_at_edge_full_level, &
                                              scalar_hwind_lap4th_at_edge_full_level)
! io
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: order
     real(r8),              intent(in)    :: coef
     real(r8),              intent(in)    :: scaling_factor(:)
     logical ,              intent(in)    :: flag
     real(r8), allocatable, intent(in)    :: scalar_hwind_lap2nd_at_edge_full_level(:,:)
     real(r8), allocatable, intent(inout) :: scalar_hwind_lap4th_at_edge_full_level(:,:)
! local
     real(r8), allocatable             :: velocity_vector_at_pc_full_level(:,:,:)
     real(r8), allocatable             :: velocity_vector_at_dc_full_level(:,:,:)
     integer(i4)                       :: iv, it, ie, ilev, inb
     integer(i4)                       :: pc1, pc2, dc1, dc2              ! cell index
     real(r8)                          :: un_pc1, un_pc2, un_dc1, un_dc2  ! normal wind follow edp's nr
     real(r8)                          :: ut_pc1, ut_pc2, ut_dc1, ut_dc2  ! tangent wind follow edp's tg
!
! calculate Laplacian vector at primal and dual cells
!        
      if(.not.allocated(velocity_vector_at_pc_full_level)) allocate(velocity_vector_at_pc_full_level(3,nlev,mesh%nv_halo(1))) ! 3 means 3 components of a cartesian vector
      if(.not.allocated(velocity_vector_at_dc_full_level)) allocate(velocity_vector_at_dc_full_level(3,nlev,mesh%nt))

      do iv  = 1, mesh%nv_halo(1)
         velocity_vector_at_pc_full_level(:,:,iv) = 0._r8
         do inb  = 1, mesh%vtx_nnb(iv)
            ie = mesh%vtx_ed(inb,iv)
            do ilev = 1, nlev
                velocity_vector_at_pc_full_level(:,ilev,iv) = velocity_vector_at_pc_full_level(:,ilev,iv)+perot_weight_at_pc(:,inb,iv)*&
                                                   scalar_hwind_lap2nd_at_edge_full_level(ilev,ie)
            end do
         end do
      end do

      do it  = 1, mesh%nt
         velocity_vector_at_dc_full_level(:,:,it) = 0._r8
         do inb = 1, mesh%tri_nnb(it)
            ie = mesh%tri_ed(inb,it)
            do ilev = 1, nlev
               velocity_vector_at_dc_full_level(:,ilev,it) = velocity_vector_at_dc_full_level(:,ilev,it)+perot_weight_at_dc(:,inb,it)*&
                                                   scalar_hwind_lap2nd_at_edge_full_level(ilev,ie)
            end do
         end do
      end do

!-----------------------------------------------------------------------
!         pc2
!          |
!    dc1---e---dc2
!          |
!         pc1
!-----------------------------------------------------------------------

      DO ie = 1, mesh%ne

         pc1     = mesh%edt_v(1,ie)
         pc2     = mesh%edt_v(2,ie)
         dc1     = mesh%edp_v(1,ie)
         dc2     = mesh%edp_v(2,ie)

         do ilev = 1, nlev
            un_pc1  = dot_product(velocity_vector_at_pc_full_level(:,ilev,pc1),mesh%edp_nr(1:3,ie))
            un_pc2  = dot_product(velocity_vector_at_pc_full_level(:,ilev,pc2),mesh%edp_nr(1:3,ie))
            un_dc1  = dot_product(velocity_vector_at_dc_full_level(:,ilev,dc1),mesh%edp_nr(1:3,ie))
            un_dc2  = dot_product(velocity_vector_at_dc_full_level(:,ilev,dc2),mesh%edp_nr(1:3,ie))

!
! high-order laplacian of normal wind for hyperviscosity
!
            scalar_hwind_lap4th_at_edge_full_level(ilev,ie) = (-1)**(order+1)*(4._r8**order)*coef*scaling_factor(ie)*&  ! coef and constant

                                                             (un_pc1+un_pc2-2._r8*scalar_hwind_lap2nd_at_edge_full_level(ilev,ie))/((rearth*mesh%edt_leng(ie))**2)+&
                                                             (un_dc1+un_dc2-2._r8*scalar_hwind_lap2nd_at_edge_full_level(ilev,ie))/((rearth*mesh%edp_leng(ie))**2)
             if(flag.and.do_dycore_laplacian_6th)then ! store 4th for later 6th
                scalar_hwind_laplace4th_at_edge_full_level%f(ilev,ie) = (un_pc1+un_pc2-2._r8*scalar_hwind_lap2nd_at_edge_full_level(ilev,ie))/((rearth*mesh%edt_leng(ie))**2)+&
                                                                        (un_dc1+un_dc2-2._r8*scalar_hwind_lap2nd_at_edge_full_level(ilev,ie))/((rearth*mesh%edp_leng(ie))**2)
             end if
         end do

      End do

      if(allocated(velocity_vector_at_pc_full_level)) deallocate(velocity_vector_at_pc_full_level)
      if(allocated(velocity_vector_at_dc_full_level)) deallocate(velocity_vector_at_dc_full_level)

      return

    end subroutine calc_hrwind_laplacian_4th_full

    subroutine calc_hori_div2d_o4_damping(mesh, nlev, coef, &
                                          scalar_normal_velocity_at_edge_full_level, &
                                          div_damping_4th_tendency)
! io
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: nlev
     real(r8),              intent(in)    :: coef
     real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
     real(r8), intent(inout) :: div_damping_4th_tendency(nlev,mesh%ne_halo(2))
! local     
     real(r8)                             :: div_sum(nlev)
     real(r8)                             :: v1v2(3), flag
     integer(i4)                          :: iv, ilev, inb, index_edge, ie
     real(r8)                             :: div_like_tendency(nlev,mesh%nv_halo(2))
     real(r8)                             :: grad_of_div(nlev,mesh%ne_halo(1))

! pass1  
!!!!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!!!!$omp do schedule(dynamic,5)
    do iv = 1, mesh%nv_halo(2)
       div_sum(1:nlev) = 0._r8
       do inb = 1, mesh%vtx_nnb(iv)
          index_edge  = mesh%vtx_ed(inb, iv)
          do ilev = 1, nlev
             div_sum(ilev)  = div_sum(ilev)+scalar_normal_velocity_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(inb, iv)*&
                                            mesh%edp_leng(index_edge)
          end do
       end do
! divergence
       div_like_tendency(1:nlev,iv) = div_sum(1:nlev)/(rearth*mesh%plg_areag(iv))
    end do
!!!!$omp end do nowait
!!!!$omp end parallel

! pass2
    do ie = 1, mesh%ne_halo(1)
       do ilev = 1, nlev
          grad_of_div(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(div_like_tendency(ilev,mesh%edt_v(2,ie))-&
                                                           div_like_tendency(ilev,mesh%edt_v(1,ie)))/(rearth*mesh%edt_leng(ie))
       end do
    end do

! pass3
    do iv = 1, mesh%nv_halo(1)
       div_sum(1:nlev) = zero
       do inb = 1, mesh%vtx_nnb(iv)
          index_edge  = mesh%vtx_ed(inb, iv)
          do ilev = 1, nlev
             div_sum(ilev)  = div_sum(ilev)+grad_of_div(ilev,index_edge)*mesh%plg_nr(inb, iv)*mesh%edp_leng(index_edge)
          end do
       end do
       div_like_tendency(1:nlev,iv) = div_sum(1:nlev)/(rearth*mesh%plg_areag(iv))
    end do

! pass4
    do ie = 1, mesh%ne_compute
       do ilev = 1, nlev
          div_damping_4th_tendency(ilev,ie) = -coef*mesh%edt_edpNr_edtTg(ie)*(div_like_tendency(ilev,mesh%edt_v(2,ie))-&
                                                                        div_like_tendency(ilev,mesh%edt_v(1,ie)))/(rearth*mesh%edt_leng(ie))
       end do
    end do

    return
    end subroutine calc_hori_div2d_o4_damping

!-------------------------------------------------------------------------
! These routines use the vector identity to calculate
! laplacian operator of hori normal wind, instead of center difference
! it uses the same compute pattern as in calc_smg_eddy_coef
! so no new parallel issue
! NOT DEFAULT
!-------------------------------------------------------------------------

    subroutine calc_hrwind_laplacian_2nd_full_vi(mesh, &
                                  scalar_normal_velocity_at_edge_full_level)
! io
     type(global_domain),   intent(inout) :: mesh
     real(r8), allocatable, intent(in) :: scalar_normal_velocity_at_edge_full_level(:,:)
! local
     real(r8), allocatable             :: div_at_pc_full_level(:,:)
     real(r8), allocatable             :: vor_at_dc_full_level(:,:)
     integer(i4)                       :: iv, it, ie, ilev, inb
     integer(i4)                       :: pc1, pc2, dc1, dc2              ! cell index
     real(r8)                          :: un_pc1, un_pc2, un_dc1, un_dc2  ! normal wind follow edp's nr
!     real(r8)                          :: ut_pc1, ut_pc2, ut_dc1, ut_dc2  ! tangent wind follow edp's tg
     real(r8)                          :: pc1pc2(3), dc1dc2(3), flag1, flag2
!     real(r8)                          :: un_pc1pc2, ut_pc1pc2, un_dc1dc2, ut_dc1dc2

!
! calculate div and relvor at primal and dual cells
!  
      if(.not.allocated(div_at_pc_full_level)) allocate(div_at_pc_full_level(nlev,mesh%nv_halo(2))) ! 3 means 3 components of a cartesian vector
      if(.not.allocated(vor_at_dc_full_level)) allocate(vor_at_dc_full_level(nlev,mesh%nt_halo(1)))

      mesh%nv = mesh%nv_halo(2)
      call divergence_operator_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                        div_at_pc_full_level, &
                                        nlev)
      mesh%nv = mesh%nv_compute

      mesh%nt = mesh%nt_halo(1)
      call curl_operator_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                  vor_at_dc_full_level,&
                                  nlev)
      mesh%nt = mesh%nt_compute
!
! use vector identicy to evaluate laplacian of normal wind, note that
! we should use plus tangent direvative as our tg is in opposite direction
! as in the typical definition of this formula
!

!-----------------------------------------------------------------------
!         pc2
!          |
!    dc1---e---dc2
!          |
!         pc1
!-----------------------------------------------------------------------

      DO ie = 1, mesh%ne_halo(1)

         pc1     = mesh%edt_v(1,ie)
         pc2     = mesh%edt_v(2,ie)
         dc1     = mesh%edp_v(1,ie)
         dc2     = mesh%edp_v(2,ie)
         pc1pc2  = mesh%vtx_p(1:3,pc2)-mesh%vtx_p(1:3,pc1)
         dc1dc2  = mesh%tri_c_p(1:3,dc2)-mesh%tri_c_p(1:3,dc1)
         flag1   = sign(1._r8,dot_product(pc1pc2,mesh%edp_nr(1:3,ie)))
         flag2   = sign(1._r8,dot_product(dc1dc2,mesh%edp_tg(1:3,ie)))

         do ilev = 1, nlev
! second-order laplacian of normal wind
            scalar_hwind_laplace_at_edge_full_level(ilev,ie) = flag1*(div_at_pc_full_level(ilev,pc2)-div_at_pc_full_level(ilev,pc1))/(rearth*mesh%edt_leng(ie))+&
                                                               flag2*(vor_at_dc_full_level(ilev,dc2)-vor_at_dc_full_level(ilev,dc1))/(rearth*mesh%edp_leng(ie))
          ! divided by 4 to give the same magnitude as the center-difference style Laplacian for final SMG 2nd
            scalar_hwind_laplace_at_edge_full_level(ilev,ie) = scalar_hwind_laplace_at_edge_full_level(ilev,ie)/4._r8
         end do

       END DO

       if(allocated(div_at_pc_full_level)) deallocate(div_at_pc_full_level)
       if(allocated(vor_at_dc_full_level)) deallocate(vor_at_dc_full_level)

       return
    end subroutine calc_hrwind_laplacian_2nd_full_vi

!
! alternative version based on the vector identity 
!

    subroutine calc_hrwind_laplacian_4th_full_vi(mesh, order, coef, scaling_factor, flag, &
                                            scalar_hwind_lap2nd_at_edge_full_level, &
                                            scalar_hwind_lap4th_at_edge_full_level)
! io
     type(global_domain),   intent(inout) :: mesh
     integer(i4),           intent(in)    :: order
     real(r8),              intent(in)    :: coef
     real(r8),              intent(in)    :: scaling_factor(:)
     logical ,              intent(in)    :: flag 
     real(r8), allocatable, intent(in)    :: scalar_hwind_lap2nd_at_edge_full_level(:,:)
     real(r8), allocatable, intent(inout) :: scalar_hwind_lap4th_at_edge_full_level(:,:)
! local
     real(r8), allocatable             :: div_at_pc_full_level(:,:)
     real(r8), allocatable             :: vor_at_dc_full_level(:,:)
     integer(i4)                       :: iv, it, ie, ilev, inb
     integer(i4)                       :: pc1, pc2, dc1, dc2              ! cell index
     real(r8)                          :: un_pc1, un_pc2, un_dc1, un_dc2  ! normal wind follow edp's nr
     real(r8)                          :: ut_pc1, ut_pc2, ut_dc1, ut_dc2  ! tangent wind follow edp's tg
     real(r8)                          :: pc1pc2(3), dc1dc2(3), flag1, flag2
     real(r8)                          :: un_pc1pc2, ut_pc1pc2, un_dc1dc2, ut_dc1dc2
!
! calculate Laplacian vector at primal and dual cells
!        
      if(.not.allocated(div_at_pc_full_level)) allocate(div_at_pc_full_level(nlev,mesh%nv_halo(1))) ! 3 means 3 components of a cartesian vector
      if(.not.allocated(vor_at_dc_full_level)) allocate(vor_at_dc_full_level(nlev,mesh%nt))


      mesh%nv = mesh%nv_halo(1)
      call divergence_operator_2d(mesh, scalar_hwind_lap2nd_at_edge_full_level, &
                                        div_at_pc_full_level, &
                                        nlev)
      mesh%nv = mesh%nv_compute

      call curl_operator_2d(mesh, scalar_hwind_lap2nd_at_edge_full_level, &
                                  vor_at_dc_full_level,&
                                  nlev)

!-----------------------------------------------------------------------
!         pc2
!          |
!    dc1---e---dc2
!          |
!         pc1
!-----------------------------------------------------------------------

      DO ie = 1, mesh%ne

         pc1     = mesh%edt_v(1,ie)
         pc2     = mesh%edt_v(2,ie)
         dc1     = mesh%edp_v(1,ie)
         dc2     = mesh%edp_v(2,ie)
         pc1pc2  = mesh%vtx_p(1:3,pc2)-mesh%vtx_p(1:3,pc1)
         dc1dc2  = mesh%tri_c_p(1:3,dc2)-mesh%tri_c_p(1:3,dc1)
         flag1   = sign(1._r8,dot_product(pc1pc2,mesh%edp_nr(1:3,ie)))
         flag2   = sign(1._r8,dot_product(dc1dc2,mesh%edp_tg(1:3,ie)))

         do ilev = 1, nlev
! second-order laplacian of normal wind
            scalar_hwind_lap4th_at_edge_full_level(ilev,ie) = flag1*(div_at_pc_full_level(ilev,pc2)-div_at_pc_full_level(ilev,pc1))/(rearth*mesh%edt_leng(ie))+&
                                                              flag2*(vor_at_dc_full_level(ilev,dc2)-vor_at_dc_full_level(ilev,dc1))/(rearth*mesh%edp_leng(ie))
          ! this is the "real" Laplacian-4th, just * -ko4_coef and do a scaling for coef, not need to divided by 4**power as
          ! center-version
            scalar_hwind_lap4th_at_edge_full_level(ilev,ie) = (-1)**(order+1)*coef*scaling_factor(ie)*scalar_hwind_lap4th_at_edge_full_level(ilev,ie)
         end do

      END DO

      if(allocated(div_at_pc_full_level)) deallocate(div_at_pc_full_level)
      if(allocated(vor_at_dc_full_level)) deallocate(vor_at_dc_full_level)

      return

    end subroutine calc_hrwind_laplacian_4th_full_vi

  end module grist_dycore_diffusion_module
