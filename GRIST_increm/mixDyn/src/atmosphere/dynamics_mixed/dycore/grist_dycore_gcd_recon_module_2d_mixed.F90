module grist_dycore_gcd_recon_module_2d_mixed

    use grist_constants,        only: i4, r8, rearth, r4 => ns, rearth_r4 => rearth_ns
    use grist_domain_types,     only: global_domain
    use grist_dycore_module_vars_control_mixed, only: edp_leng, plg_areag

    implicit none

contains

  subroutine divergence_operator_2d(mesh, scalar_normal_flux_at_edge_full_level, &
                                          tend_div_at_pc_full_level, &
                                          nlev)
! io
   use omp_lib
   type(global_domain)  ,  intent(in)   :: mesh
   real(r4),               intent(in)   :: scalar_normal_flux_at_edge_full_level(nlev,mesh%ne_full)
   real(r4),               intent(inout):: tend_div_at_pc_full_level(nlev,mesh%nv_full)
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
             div_sum(ilev)  = div_sum(ilev)+scalar_normal_flux_at_edge_full_level(ilev,index_edge)*&
                                            mesh%plg_nr(inb, iv)*&
                                            rearth*mesh%edp_leng(index_edge)
          end do
       end do
       tend_div_at_pc_full_level(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
    end do
!$omp end do nowait
!$omp end parallel 

   return
  end subroutine divergence_operator_2d

  subroutine divergence_operator_2d_var2(mesh, scalar1_normal_flux_at_edge_full_level,&
                                               scalar2_normal_flux_at_edge_full_level,&
                                               tend1_div_at_pc_full_level            ,&
                                               tend2_div_at_pc_full_level            ,&
                                               nlev)
! io
   use omp_lib
   type(global_domain)  ,  intent(in)   :: mesh
   real(r4), allocatable,  intent(in)   :: scalar1_normal_flux_at_edge_full_level(:,:)
   real(r4), allocatable,  intent(in)   :: scalar2_normal_flux_at_edge_full_level(:,:)
   real(r4), allocatable,  intent(inout):: tend1_div_at_pc_full_level(:,:)
   real(r4), allocatable,  intent(inout):: tend2_div_at_pc_full_level(:,:)
   integer(i4)          ,  intent(in)   :: nlev
! local
   real(r4)                             :: div_sum(2,nlev)
   integer(i4)                          :: iv, ilev, inb
   integer(i4)                          :: index_edge

!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!$omp do schedule(dynamic,5)
    do iv = 1, mesh%nv
       div_sum(:,:) = 0._r4
       do inb = 1, mesh%vtx_nnb(iv)
          index_edge  = mesh%vtx_ed(inb, iv)
          do ilev = 1, nlev
             div_sum(1,ilev)  = div_sum(1,ilev)+scalar1_normal_flux_at_edge_full_level(ilev,index_edge)*&
                                                mesh%plg_nr(inb, iv)*&
                                                rearth_r4*edp_leng(index_edge)
             div_sum(2,ilev)  = div_sum(2,ilev)+scalar2_normal_flux_at_edge_full_level(ilev,index_edge)*&
                                                mesh%plg_nr(inb, iv)*&
                                                rearth_r4*edp_leng(index_edge)
          end do
       end do
       tend1_div_at_pc_full_level(:,iv) = div_sum(1,:)/((rearth**2)*mesh%plg_areag(iv))
       tend2_div_at_pc_full_level(:,iv) = div_sum(2,:)/((rearth**2)*mesh%plg_areag(iv))
    end do
!$omp end do nowait
!$omp end parallel 

   return
  end subroutine divergence_operator_2d_var2

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

end module grist_dycore_gcd_recon_module_2d_mixed