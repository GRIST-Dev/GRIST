
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: 2D version of some selected operators in gcd and recon module
! Revision history:
!----------------------------------------------------------------------------

 module grist_dycore_gcd_recon_module_2d

  use grist_constants,       only: gravity, i4, r8, rearth, omega, zero
  use grist_domain_types  ,  only: global_domain
  use grist_math_module,     only: convert_vector_cart2sph, norm
  use grist_mesh_weight_icosh,only: project_sphere2tplane
  use grist_dycore_hori_swe_module_2d, only: calc_coriolis_term

  implicit none

  private

  public :: gradient_operator_2d       , &
            gradient_operator_2d_var2  , &
            gradient_operator_2d_var3  , &
            gradient_operator_2d_var4  , &
            divergence_operator_2d     , &
            curl_operator_2d           , &
            divergence_operator_2d_var2, &
            project_uv_2d              , &
            calc_vorticity_at_dual_cell_2d, &
            vector_recon_perot_edge2cell_uv_2d

  contains

  subroutine gradient_operator_2d(mesh, scalar_at_pc_full_level, &
                                        tend_grad_at_edge_full_level, &
                                        nlev)
!io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_at_pc_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend_grad_at_edge_full_level(:,:)
   integer(i4)          , intent(in)    :: nlev
!local
   real(r8)                             :: v1v2(3)
   real(r8)                             :: flag
   integer(i4)                          :: ilev, ie, v1, v2

   do ie = 1, mesh%ne
      v1                         = mesh%edt_v(1, ie)
      v2                         = mesh%edt_v(2, ie)
      do ilev = 1, nlev
         tend_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar_at_pc_full_level(ilev,v2)-scalar_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
      end do
   end do

   return
  end subroutine gradient_operator_2d

!=======================================================================
! same as gradient_operator_2d but simultanesouly return two tendencies
!=======================================================================

  subroutine gradient_operator_2d_var2(mesh, scalar1_at_pc_full_level, &
                                             scalar2_at_pc_full_level, &
                                             tend1_grad_at_edge_full_level, &
                                             tend2_grad_at_edge_full_level, &
                                             nlev)
!io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar1_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar2_at_pc_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend1_grad_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend2_grad_at_edge_full_level(:,:)
   integer(i4)          , intent(in)    :: nlev
!local
   real(r8)                             :: v1v2(3)
   real(r8)                             :: flag
   integer(i4)                          :: ilev, ie, v1, v2

   do ie = 1, mesh%ne
      v1                         = mesh%edt_v(1, ie)
      v2                         = mesh%edt_v(2, ie)
      do ilev = 1, nlev
         tend1_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar1_at_pc_full_level(ilev,v2)-scalar1_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
         tend2_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar2_at_pc_full_level(ilev,v2)-scalar2_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
      end do
   end do

   return
  end subroutine gradient_operator_2d_var2

!=======================================================================
! same as gradient_operator_2d but simultanesouly return 3 tendencies
!=======================================================================

  subroutine gradient_operator_2d_var3(mesh, scalar1_at_pc_full_level     , &
                                             scalar2_at_pc_full_level     , &
                                             scalar3_at_pc_face_level     , &
                                             tend1_grad_at_edge_full_level, &
                                             tend2_grad_at_edge_full_level, &
                                             tend3_grad_at_edge_face_level, &
                                             nlev)
!io
   use omp_lib
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar1_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar2_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar3_at_pc_face_level(:,:)
   real(r8), allocatable, intent(inout) :: tend1_grad_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend2_grad_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend3_grad_at_edge_face_level(:,:)
   integer(i4)          , intent(in)    :: nlev
!local
   real(r8)                             :: v1v2(3)
   real(r8)                             :: flag
   integer(i4)                          :: ilev, ie, v1, v2

!$omp parallel  private(ie,v1,v2,ilev)     
!$omp do schedule(dynamic,50) 
   do ie = 1, mesh%ne
      v1                         = mesh%edt_v(1, ie)
      v2                         = mesh%edt_v(2, ie)
      do ilev = 1, nlev
         tend1_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar1_at_pc_full_level(ilev,v2)-scalar1_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
         tend2_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar2_at_pc_full_level(ilev,v2)-scalar2_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
      end do
      do ilev = 1, nlev+1
         tend3_grad_at_edge_face_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar3_at_pc_face_level(ilev,v2)-scalar3_at_pc_face_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
      end do
   end do
!$omp end do nowait
!$omp end parallel 

   return
  end subroutine gradient_operator_2d_var3

  subroutine gradient_operator_2d_var4(mesh, scalar1_at_pc_full_level     , &
                                             scalar2_at_pc_full_level     , &
                                             scalar3_at_pc_full_level     , &
                                             scalar4_at_pc_face_level     , &
                                             tend1_grad_at_edge_full_level, &
                                             tend2_grad_at_edge_full_level, &
                                             tend3_grad_at_edge_full_level, &
                                             tend4_grad_at_edge_face_level, &
                                             nlev)
!io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar1_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar2_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar3_at_pc_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar4_at_pc_face_level(:,:)
   real(r8), allocatable, intent(inout) :: tend1_grad_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend2_grad_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend3_grad_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: tend4_grad_at_edge_face_level(:,:)
   integer(i4)          , intent(in)    :: nlev
!local
   real(r8)                             :: v1v2(3)
   real(r8)                             :: flag
   integer(i4)                          :: ilev, ie, v1, v2

   do ie = 1, mesh%ne
      v1                         = mesh%edt_v(1, ie)
      v2                         = mesh%edt_v(2, ie)
      do ilev = 1, nlev
         tend1_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar1_at_pc_full_level(ilev,v2)-scalar1_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
         tend2_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar2_at_pc_full_level(ilev,v2)-scalar2_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
         tend3_grad_at_edge_full_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar3_at_pc_full_level(ilev,v2)-scalar3_at_pc_full_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
      end do
      do ilev = 1, nlev+1
         tend4_grad_at_edge_face_level(ilev,ie) = mesh%edt_edpNr_edtTg(ie)*(scalar4_at_pc_face_level(ilev,v2)-scalar4_at_pc_face_level(ilev,v1))/(rearth*mesh%edt_leng(ie))
      end do
   end do

   return
  end subroutine gradient_operator_2d_var4

  subroutine curl_operator_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                    scalar_relative_vorticity_at_dual_cell_full_level,&
                                    nlev)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: scalar_relative_vorticity_at_dual_cell_full_level(:,:)
   integer(i4),           intent(in)    :: nlev
! local
   integer(i4)                          :: it
   integer(i4)                          :: ie, ilev
   integer(i4)                          :: index_of_edge
   integer(i4)                          :: t_edge_v
   real(r8)                             :: length_of_edge
   real(r8)                             :: scalar_normal_velocity(nlev)
   real(r8)                             :: tmp(nlev)

   do it = 1, mesh%nt
      tmp = 0._r8
      do ie = 1, mesh%tri_nnb(it)
         index_of_edge  = mesh%tri_ed(ie, it)
#ifndef DYCORE_OPT
         length_of_edge = rearth*mesh%edt_leng(index_of_edge)
#else
         length_of_edge = mesh%edt_leng(index_of_edge)
#endif
         t_edge_v       = mesh%edt_edpTg_edtNr(index_of_edge)*mesh%tri_nr(ie, it) ! R10: kXnr(-hxtg) relative to dual cell center
         do ilev = 1, nlev
            scalar_normal_velocity(ilev) = scalar_normal_velocity_at_edge_full_level(ilev,index_of_edge)
            tmp(ilev)  = tmp(ilev)+t_edge_v*length_of_edge*scalar_normal_velocity(ilev)
         end do
      end do
#ifndef DYCORE_OPT
      scalar_relative_vorticity_at_dual_cell_full_level(:,it) = tmp(:)/((rearth**2)*mesh%tri_areag(it))
#else
      scalar_relative_vorticity_at_dual_cell_full_level(:,it) = tmp(:)/(rearth*mesh%tri_areag(it))
#endif
   end do

   return
  end subroutine curl_operator_2d

  subroutine divergence_operator_2d(mesh, scalar_normal_flux_at_edge_full_level, &
                                          tend_div_at_pc_full_level, &
                                          nlev)
! io
   use omp_lib
   type(global_domain)  ,  intent(in)   :: mesh
   real(r8),               intent(in)   :: scalar_normal_flux_at_edge_full_level(nlev,mesh%ne_full)
   real(r8),               intent(inout):: tend_div_at_pc_full_level(nlev,mesh%nv_full)
   integer(i4)          ,  intent(in)   :: nlev
! local
   real(r8)                             :: div_sum(nlev)
   integer(i4)                          :: iv, ilev, inb
   integer(i4)                          :: index_edge

!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do schedule(dynamic,5)
    do iv = 1, mesh%nv
       div_sum(:) = 0._r8
       do inb = 1, mesh%vtx_nnb(iv)
          index_edge  = mesh%vtx_ed(inb, iv)
          do ilev = 1, nlev
#ifndef DYCORE_OPT
             div_sum(ilev)  = div_sum(ilev)+scalar_normal_flux_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
#else
             div_sum(ilev)  = div_sum(ilev)+scalar_normal_flux_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(inb, iv)*mesh%edp_leng(index_edge)
#endif
          end do
       end do
#ifndef DYCORE_OPT
       tend_div_at_pc_full_level(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
#else
       tend_div_at_pc_full_level(:,iv) = div_sum(:)/(rearth*mesh%plg_areag(iv))
#endif
    end do
!$omp end do nowait
!$omp end parallel 

   return
  end subroutine divergence_operator_2d

!-----------------------------------------------------------------------
! bespoke for nhd flux divergence
!-----------------------------------------------------------------------

  subroutine divergence_operator_2d_var2(mesh, scalar1_normal_flux_at_edge_full_level,&
                                               scalar2_normal_flux_at_edge_full_level,&
                                               tend1_div_at_pc_full_level            ,&
                                               tend2_div_at_pc_full_level            ,&
                                               nlev)
! io
   use omp_lib
   type(global_domain)  ,  intent(in)   :: mesh
   real(r8), allocatable,  intent(in)   :: scalar1_normal_flux_at_edge_full_level(:,:)
   real(r8), allocatable,  intent(in)   :: scalar2_normal_flux_at_edge_full_level(:,:)
   real(r8),               intent(inout):: tend1_div_at_pc_full_level(nlev,mesh%nv_full)
   real(r8),               intent(inout):: tend2_div_at_pc_full_level(nlev,mesh%nv_full)
   integer(i4)          ,  intent(in)   :: nlev
! local
   real(r8)                             :: div_sum(2,nlev)
   integer(i4)                          :: iv, ilev, inb
   integer(i4)                          :: index_edge

!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!$omp do schedule(dynamic,5)
    do iv = 1, mesh%nv
       div_sum(:,:) = 0._r8
       do inb = 1, mesh%vtx_nnb(iv)
          index_edge  = mesh%vtx_ed(inb, iv)
          do ilev = 1, nlev
#ifndef DYCORE_OPT
             div_sum(1,ilev)  = div_sum(1,ilev)+scalar1_normal_flux_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
             div_sum(2,ilev)  = div_sum(2,ilev)+scalar2_normal_flux_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
#else
             div_sum(1,ilev)  = div_sum(1,ilev)+scalar1_normal_flux_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(inb, iv)*mesh%edp_leng(index_edge)
             div_sum(2,ilev)  = div_sum(2,ilev)+scalar2_normal_flux_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(inb, iv)*mesh%edp_leng(index_edge)
#endif
          end do
       end do
#ifndef DYCORE_OPT
       tend1_div_at_pc_full_level(:,iv) = div_sum(1,:)/((rearth**2)*mesh%plg_areag(iv))
       tend2_div_at_pc_full_level(:,iv) = div_sum(2,:)/((rearth**2)*mesh%plg_areag(iv))
#else
       tend1_div_at_pc_full_level(:,iv) = div_sum(1,:)/(rearth*mesh%plg_areag(iv))
       tend2_div_at_pc_full_level(:,iv) = div_sum(2,:)/(rearth*mesh%plg_areag(iv))
#endif
    end do
!$omp end do nowait
!$omp end parallel 

   return
  end subroutine divergence_operator_2d_var2

  subroutine calc_vorticity_at_dual_cell_2d(mesh, scalar_normal_velocity_at_edge_full_level,&
                                                  scalar_absolute_vorticity_at_dual_cell_full_level,&
                                                  scalar_relative_vorticity_at_dual_cell_full_level,&
                                                  nlev)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: scalar_absolute_vorticity_at_dual_cell_full_level(:,:)
   real(r8), allocatable, intent(inout) :: scalar_relative_vorticity_at_dual_cell_full_level(:,:)
   integer(i4)          , intent(in)    :: nlev
! local
   integer(i4)                          :: ilev, ie
   integer(i4)                          :: it
   integer(i4)                          :: v1
   integer(i4)                          :: v2
   real(r8)                             :: coriolis

   call curl_operator_2d(mesh,scalar_normal_velocity_at_edge_full_level,&
                              scalar_relative_vorticity_at_dual_cell_full_level,&
                              nlev)
   do it = 1, mesh%nt
      do ilev = 1, nlev
         coriolis = 2._r8*omega*sin(mesh%tri_c_lat(it))
         scalar_absolute_vorticity_at_dual_cell_full_level(ilev,it) = scalar_relative_vorticity_at_dual_cell_full_level(ilev,it)+coriolis
      end do
   end do

   return
  end subroutine calc_vorticity_at_dual_cell_2d

!==================================================
!  Given a normal component of vector, reconstruct
! U and V components
!==================================================

  subroutine project_uv_2d(mesh,scalar_normal_velocity_at_edge_full_level,&
                                scalar_u_full_level, scalar_v_full_level, &
                                nlev)

! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: scalar_u_full_level(:,:)
   real(r8), allocatable, intent(inout) :: scalar_v_full_level(:,:)
   integer(i4),           intent(in)    :: nlev
! local
   real(r8), allocatable                :: scalar_tangent_velocity_at_edge_full_level(:,:)
   real(r8), allocatable                :: tmp(:,:)
   real(r8)                             :: cos_theta_east
   real(r8)                             :: cos_theta_north
   real(r8)                             :: vector(1:3)
   real(r8)                             :: utmp1
   real(r8)                             :: vtmp1
   integer(i4)                          :: ie, ilev

   allocate(scalar_tangent_velocity_at_edge_full_level(nlev,mesh%ne_full))
   allocate(tmp(nlev,mesh%ne_full))

   tmp(:,:) = 1._r8
!
!    We let pv_at_edge(tmp)=1, and input Normal 
!  V, routine calc_coriolis_term can be directly
!  used to reconstruct tangent wind.
!
    call calc_coriolis_term(mesh,tmp,&
                                 scalar_normal_velocity_at_edge_full_level,&
                                 scalar_tangent_velocity_at_edge_full_level,&
                                 nlev)

    tmp = scalar_normal_velocity_at_edge_full_level

    do ie  = 1, mesh%ne
       do ilev = 1, nlev
       vector(:) = scalar_normal_velocity_at_edge_full_level(ilev,ie)*mesh%edp_nr(1:3,ie)+&
                  scalar_tangent_velocity_at_edge_full_level(ilev,ie)*mesh%edp_tg(1:3,ie)

       !call convert_vector_cart2sph(real(mesh%edp_c_p(1:3,ie),r8), vector , utmp1, vtmp1)
       ! change here to edt_c_p such that edp_c_p can be deleted
       call convert_vector_cart2sph(real(mesh%edt_c_p(1:3,ie),r8), vector , utmp1, vtmp1)

       scalar_u_full_level(ilev,ie) = utmp1
       scalar_v_full_level(ilev,ie) = vtmp1

       end do
    end do

    deallocate(scalar_tangent_velocity_at_edge_full_level)
    deallocate(tmp)

    return
   end subroutine project_uv_2d

!
! use perot reconstrucion to recover UV at centroid
! based on normal wind at edges.
!
   subroutine vector_recon_perot_edge2cell_uv_2d(mesh,scalar_normal_velocity_at_edge_full_level,&
                                                      scalar_u_full_level, scalar_v_full_level,&
                                                      nlev)

! io
   use omp_lib
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: scalar_u_full_level(:,:)  ! at Voronoi centroid
   real(r8), allocatable, intent(inout) :: scalar_v_full_level(:,:)  ! at Voronoi centroid
   integer(i4)          , intent(in)    :: nlev
! local
   real(r8)                             :: vector_sum(1:3,nlev),vector_tmp(1:3),local_basis(1:3)
   real(r8)                             :: utmp, vtmp
   integer(i4)                          :: iv ,inb, ie, ilev

!$omp parallel private(iv,inb,vector_sum,ie,vector_tmp,local_basis,ilev,utmp,vtmp) 
!$omp do schedule(dynamic,10) 
      do iv  = 1, mesh%nv
         vector_sum = 0._r8
         do inb  = 1, mesh%vtx_nnb(iv)
            ie = mesh%vtx_ed(inb,iv)
            call project_sphere2tplane(real(mesh%edt_c_p(1:3,ie),r8),real(mesh%vtx_p(1:3,iv),r8),local_basis) ! unit sphere
            vector_tmp = dot_product(local_basis/norm(local_basis),mesh%edp_nr(1:3,ie))*local_basis*mesh%edp_leng(ie)/(mesh%plg_areag(iv))
            do ilev = 1, nlev
               vector_sum(:,ilev) = vector_sum(:,ilev) + vector_tmp*scalar_normal_velocity_at_edge_full_level(ilev,ie)
               !vector_sum(:,ilev) = vector_sum(:,ilev) +
               !perot_weight_at_pc(:,inb,iv)*scalar_normal_velocity_at_edge_full_level(ilev,ie) !affect bit repro due to var sequence of
               !calculating perot weight is different
            end do
         end do
         do ilev = 1, nlev
            call convert_vector_cart2sph(real(mesh%vtx_p(1:3,iv),r8),real(vector_sum(:,ilev),r8),utmp,vtmp)
            scalar_u_full_level(ilev,iv)  = utmp
            scalar_v_full_level(ilev,iv)  = vtmp
         end do
      end do
!$omp end do nowait
!$omp end parallel 

      return
  
  end subroutine vector_recon_perot_edge2cell_uv_2d

 end module grist_dycore_gcd_recon_module_2d
