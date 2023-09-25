
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description:  For 3D tracer transport with FFSL flux in the sequantial code, 
!               1d version with ntracer outside  is faster than 2d version 
!               with ntracer inside, this needs to be checked in the dm version
! Revision history:
!----------------------------------------------------------------------------

 module grist_tracer_transport_ffsl_module

  use grist_domain_types    , only: global_domain
  use grist_constants       , only: i4, r8, rearth
  use grist_math_module     , only: norm

  implicit none

  private
  public ::  tracer_transport_ffsl_flux_sm10, &
             tracer_transport_ffsl_flux_sm10_init, &
             tracer_transport_ffsl_flux_sm10_final

    real(r8),    parameter     :: one = 1._r8
    real(r8),    allocatable   :: v1v2_global(:,:)
    real(r8),    allocatable   :: v3v4_global(:,:)
    real(r8),    allocatable   :: t1_point(:,:)
    real(r8),    allocatable   :: t2_point(:,:)
    real(r8),    allocatable   :: ed_point(:,:)
    integer(i4), allocatable   :: v1v2_edp_nr(:)
    integer(i4), allocatable   :: v3v4_edp_tg(:)

  contains

  subroutine tracer_transport_ffsl_flux_sm10_init(mesh)
    type(global_domain),   intent(in)    :: mesh
! local
    integer(i4)  :: ie, v1, v2,v3,v4

     if(.not.allocated(v1v2_edp_nr)) allocate(v1v2_edp_nr(mesh%ne))
     if(.not.allocated(v3v4_edp_tg)) allocate(v3v4_edp_tg(mesh%ne))
     if(.not.allocated(v1v2_global)) allocate(v1v2_global(3,mesh%ne))
     if(.not.allocated(v3v4_global)) allocate(v3v4_global(3,mesh%ne))
     if(.not.allocated(t1_point))    allocate(t1_point(3,mesh%ne))
     if(.not.allocated(t2_point))    allocate(t2_point(3,mesh%ne))
     if(.not.allocated(ed_point))    allocate(ed_point(3,mesh%ne))

     do ie = 1, mesh%ne

        !v1             = mesh%edt(ie)%v(1)
        !v2             = mesh%edt(ie)%v(2)
        v1             = mesh%edt_v(1,ie)
        v2             = mesh%edt_v(2,ie)
        v1v2_global(:,ie) = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)

        !v3             = mesh%edp(ie)%v(1)
        !v4             = mesh%edp(ie)%v(2)
        v3             = mesh%edp_v(1,ie)
        v4             = mesh%edp_v(2,ie)
        v3v4_global(:,ie) = mesh%tri_c_p(1:3,v4)-mesh%tri_c_p(1:3,v3)
        v1v2_edp_nr(ie)= sign(1._r8,dot_product(v1v2_global(:,ie)/norm(v1v2_global(:,ie)),real(mesh%edp_nr(1:3,ie),r8))) ! v1v2 direction relative to edp's nr
        v3v4_edp_tg(ie)= sign(1._r8,dot_product(v3v4_global(:,ie)/norm(v3v4_global(:,ie)),real(mesh%edp_tg(1:3,ie),r8))) ! v3v4 direction relative to edp's tg

        !t1_point(:,ie) = mesh%tri(mesh%edp(ie)%v(1))%c%p
        !t2_point(:,ie) = mesh%tri(mesh%edp(ie)%v(2))%c%p
        t1_point(:,ie) = mesh%tri_c_p(1:3,mesh%edp_v(1,ie))
        t2_point(:,ie) = mesh%tri_c_p(1:3,mesh%edp_v(2,ie))
        ed_point(:,ie) =(t1_point(:,ie)+t2_point(:,ie))*0.5_r8
        ed_point(:,ie) = ed_point(:,ie)/norm(ed_point(:,ie))

     end do
    
    return
  end subroutine tracer_transport_ffsl_flux_sm10_init

  subroutine tracer_transport_ffsl_flux_sm10_final

    deallocate(v1v2_edp_nr)
    deallocate(v3v4_edp_tg)
    deallocate(v1v2_global)
    deallocate(v3v4_global)
    deallocate(t1_point)
    deallocate(t2_point)
    deallocate(ed_point)

  end subroutine tracer_transport_ffsl_flux_sm10_final

!================================================
!       SM10 FLUX AT THE HEXAGONAL EDGE
!================================================

  subroutine tracer_transport_ffsl_flux_sm10(mesh,scalar_normal_velocity_at_edge_full_level , &
                                                  scalar_tangen_velocity_at_edge_full_level , &
                                                  scalar_tracer_mxrt_at_pc_full_level       , &
                                                  scalar_normal_mxrt_at_edge_full_level     , &
                                                  dtime, nlev, ntracer)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar_tangen_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
   real(r8), allocatable, intent(inout) :: scalar_normal_mxrt_at_edge_full_level(:,:,:)
   real(r8)   ,           intent(in)    :: dtime
   integer(i4),           intent(in)    :: nlev
   integer(i4),           intent(in)    :: ntracer
! local
   real(r8)                             :: v1v2(3), v3v4(3) ! 1->2, 3->4
   real(r8)                             :: g1(3),g2(3),g3(3),g4(3),g5(3)
   real(r8)                             :: nr_wind_edge
   real(r8)                             :: tg_wind_edge
   real(r8)                             :: g1_lob_x, g2_lob_x, g3_lob_x, g4_lob_x, g5_lob_x
   real(r8)                             :: g1_lob_y, g2_lob_y, g3_lob_y, g4_lob_y, g5_lob_y
   real(r8)                             :: scalar_upwind(ntracer)
   real(r8)                             :: scalar_tmp(8)  ! assume maxnnb = 8
   real(r8)                             :: alpha
   integer(i4)                          :: v_upwind, mark
   integer(i4)                          :: iv,ie, ilev, itracer
   integer(i4)                          :: flag, flag0, flag1, flag2

!================================================
! normal wind follows edp's nr positive
! tangen wind follows edp's tg positive
!================================================

     do ie = 1, mesh%ne

        v1v2       = v1v2_global(:,ie)
        v3v4       = v3v4_global(:,ie)
        flag       = v1v2_edp_nr(ie)  ! v1v2 direction relative to nr
        flag0      = v3v4_edp_tg(ie)  ! v3v4 direction relative to tg
     !   flag       = sign(1._r8,dot_product(v1v2/norm(v1v2),mesh%edp(ie)%nr)) ! v1v2 direction relative to edp's nr
     !   flag0      = sign(1._r8,dot_product(v3v4/norm(v3v4),mesh%edp(ie)%tg)) ! v3v4 direction relative to edp's tg

        do ilev = 1, nlev

           flag1         = sign(1._r8,scalar_normal_velocity_at_edge_full_level(ilev,ie))  ! sign of normal wind
           flag2         = sign(1._r8,scalar_tangen_velocity_at_edge_full_level(ilev,ie))  ! sign of tangen wind
           nr_wind_edge  = scalar_normal_velocity_at_edge_full_level(ilev,ie)
           tg_wind_edge  = scalar_tangen_velocity_at_edge_full_level(ilev,ie)

           mark      = int(-0.5_r8*sign(1,flag*flag1)+1.5_r8)
           !v_upwind  = mesh%edt(ie)%v(mark)
           v_upwind  = mesh%edt_v(mark,ie)

!================================================
! compute Gauss quadrature points
! alpha=1 recovers de original formulation
! some optimization can be made here
! to further pre-store time-invariant geometic 
! pattern in future
!================================================
           alpha = 1._r8

           g2   = ed_point(:,ie) !+abs(nr_wind_edge*dtime)*(1.-alpha)/rearth*v1v2/norm(v1v2)*flag*flag1
           g2   =         g2 !+abs(tg_wind_edge*dtime)*(1.-alpha)/rearth*v3v4/norm(v3v4)*flag0*flag2

           g3   = ed_point(:,ie)-abs(nr_wind_edge*dtime)*alpha/rearth*v1v2/norm(v1v2)*flag*flag1
           g3   =               g3-abs(tg_wind_edge*dtime)*alpha/rearth*v3v4/norm(v3v4)*flag0*flag2

           g1   = ed_point(:,ie)+0.5_r8*abs(nr_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v1v2/norm(v1v2)*flag*flag1
           g1   =         g1+0.5_r8*abs(tg_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v3v4/norm(v3v4)*flag0*flag2

           g4   = t1_point(:,ie)+0.5_r8*abs(nr_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v1v2/norm(v1v2)*flag*flag1
           g4   =             g4+0.5_r8*abs(tg_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v3v4/norm(v3v4)*flag0*flag2

           g5   = t2_point(:,ie)+0.5_r8*abs(nr_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v1v2/norm(v1v2)*flag*flag1
           g5   =             g5+0.5_r8*abs(tg_wind_edge*dtime)*(1._r8-2._r8*alpha)/rearth*v3v4/norm(v3v4)*flag0*flag2



! only upwind set
           g1   = g1-mesh%vtx_p(1:3,v_upwind)! vector points from upwind cell to centroid
           g2   = g2-mesh%vtx_p(1:3,v_upwind)! vector points from upwind cell to centroid
           g3   = g3-mesh%vtx_p(1:3,v_upwind)! vector points from upwind cell to centroid
           g4   = g4-mesh%vtx_p(1:3,v_upwind)! vector points from upwind cell to centroid
           g5   = g5-mesh%vtx_p(1:3,v_upwind)! vector points from upwind cell to centroid

           g1_lob_x = dot_product(g1,mesh%vtx_lob_nx(1,:,v_upwind))
           g1_lob_y = dot_product(g1,mesh%vtx_lob_ny(1,:,v_upwind))
           g2_lob_x = dot_product(g2,mesh%vtx_lob_nx(1,:,v_upwind))
           g2_lob_y = dot_product(g2,mesh%vtx_lob_ny(1,:,v_upwind))
           g3_lob_x = dot_product(g3,mesh%vtx_lob_nx(1,:,v_upwind))
           g3_lob_y = dot_product(g3,mesh%vtx_lob_ny(1,:,v_upwind))
           g4_lob_x = dot_product(g4,mesh%vtx_lob_nx(1,:,v_upwind))
           g4_lob_y = dot_product(g4,mesh%vtx_lob_ny(1,:,v_upwind))
           g5_lob_x = dot_product(g5,mesh%vtx_lob_nx(1,:,v_upwind))
           g5_lob_y = dot_product(g5,mesh%vtx_lob_ny(1,:,v_upwind))
           call calc_sm10_area_at_hx(mesh,scalar_tracer_mxrt_at_pc_full_level(:,ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind),v_upwind)), &
                                             scalar_tracer_mxrt_at_pc_full_level(:,ilev,v_upwind), v_upwind, &
                                             g1_lob_x, g1_lob_y, &
                                             g2_lob_x, g2_lob_y, &
                                             g3_lob_x, g3_lob_y, &
                                             g4_lob_x, g4_lob_y, &
                                             g5_lob_x, g5_lob_y, &
                                             scalar_upwind     , &
                                             mesh%plg_stencil_number_2nd(v_upwind),ntracer)

          scalar_normal_mxrt_at_edge_full_level(:,ilev,ie) = scalar_upwind

       end do
     end do

   return
  end subroutine tracer_transport_ffsl_flux_sm10

!================================================
! second-order reconstruction here
!================================================

  subroutine calc_sm10_area_at_hx(mesh, s_array, s0, iv   ,& 
                                        g1_lob_x, g1_lob_y,&
                                        g2_lob_x, g2_lob_y,&
                                        g3_lob_x, g3_lob_y,&
                                        g4_lob_x, g4_lob_y,&
                                        g5_lob_x, g5_lob_y,&
                                        scalar_upwind     ,&
                                        nlength, ntracer)
! io
   type(global_domain),   intent(in) :: mesh
   real(r8)            ,  intent(in) :: s_array(:,:) ! ntracer, nlength
   real(r8)            ,  intent(in) :: s0(:)        ! ntracer
   integer(i4)         ,  intent(in) :: iv
   real(r8)            ,  intent(in) :: g1_lob_x, g1_lob_y
   real(r8)            ,  intent(in) :: g2_lob_x, g2_lob_y
   real(r8)            ,  intent(in) :: g3_lob_x, g3_lob_y
   real(r8)            ,  intent(in) :: g4_lob_x, g4_lob_y
   real(r8)            ,  intent(in) :: g5_lob_x, g5_lob_y
   real(r8)            ,  intent(out):: scalar_upwind(:) ! ntracer
   integer(i4)         ,  intent(in) :: nlength
   integer(i4)         ,  intent(in) :: ntracer
! local
   real(r8)                          :: b_array(5,nlength)
   real(r8)                          :: m_array(nlength)
   real(r8)                          :: f_array(ntracer,5)
   real(r8)                          :: sub_gauss_inte(ntracer)
   real(r8)                          :: total_sum(ntracer)
   real(r8)                          :: correction(ntracer)
   integer(i4)                       :: ii,inb,itracer

      b_array      = mesh%plg_blocal(:,:,iv)
      !b_array      = mesh%plg(iv)%blocal(1,:,:)
      do itracer   = 1, ntracer
         m_array   = s_array(itracer,:)-s0(itracer)
         f_array(itracer,:) = matmul(b_array,m_array(1:nlength))
      end do
!
! compute correction
!
      total_sum(:) = 0._r8
      do inb = 1, mesh%vtx_nnb(iv)
         sub_gauss_inte(:) = 0._r8
         do ii = 1, 3
           do itracer = 1, ntracer
              sub_gauss_inte(itracer) = sub_gauss_inte(itracer)+polyfit2d(s0(itracer), &
                                                       f_array(itracer,:),&
                                                       real(mesh%plg_sub_triangle_midp_lob(inb,ii,1,iv),r8),&
                                                       real(mesh%plg_sub_triangle_midp_lob(inb,ii,2,iv),r8))
           end do
         end do
         sub_gauss_inte(:) = sub_gauss_inte(:)*mesh%plg_sub_triangle_area(inb,iv)/3._r8
         total_sum         = total_sum+sub_gauss_inte
      end do
      total_sum     = total_sum/mesh%plg_areag(iv)
      correction    = total_sum-s0

      do itracer = 1, ntracer
         scalar_upwind(itracer) = (2*polyfit2d(s0(itracer), f_array(itracer,:),g1_lob_x,g1_lob_y)+&
                                     polyfit2d(s0(itracer), f_array(itracer,:),g2_lob_x,g2_lob_y)+& 
                                     polyfit2d(s0(itracer), f_array(itracer,:),g3_lob_x,g3_lob_y)+& 
                                     polyfit2d(s0(itracer), f_array(itracer,:),g4_lob_x,g4_lob_y)+& 
                                     polyfit2d(s0(itracer), f_array(itracer,:),g5_lob_x,g5_lob_y))/6._r8-correction(itracer)
      end do

      return
  end subroutine calc_sm10_area_at_hx

  real(r8) function polyfit2d(scalar0,f_array,lob_x,lob_y)
! io
   real(r8),  intent(in) :: scalar0
   real(r8),  intent(in) :: f_array(:)
   real(r8),  intent(in) :: lob_x, lob_y

     polyfit2d = scalar0+f_array(1)*lob_x+&
                         f_array(2)*lob_y+&
                         f_array(3)*lob_x**2+&
                         f_array(4)*lob_x*lob_y+&
                         f_array(5)*lob_y**2
    return
   end function polyfit2d

  end module grist_tracer_transport_ffsl_module
