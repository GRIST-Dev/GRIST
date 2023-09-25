
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: initilize some additional mesh info not in
! mesh_weight and grid generator.
! Revision history:
!       O58 similar to O8, but using double reconstruction
!----------------------------------------------------------------------------

  module grist_geometric_info_icosh
  
  use grist_data_types  ,    only: exchange_field_list_2d, scalar_2d_field
  use grist_domain_types,    only: global_domain
  use grist_constants_dbl,   only: r8, i4, eps
  use grist_mpi
  use grist_math_module,     only: det, tri_kite_areas
  use grist_config_partition,only: exchange_data_2d_add, exchange_data_2d
  use grist_handle_error,    only: endrun
  use grist_element_types_icosh, only: maxnnb, nDual
  use grist_nml_module,      only: stencil_width
  
  implicit none

    private
    public :: init_geometric_info,&
              config_plg_tri_nrtg_corrector, &
              init_geometric_info_global ! old seq code

  type(exchange_field_list_2d), pointer :: field_head_2d
  contains

   subroutine init_geometric_info_global(mesh)
! io
    type(global_domain), intent(inout) :: mesh
! local
    integer(i4)     :: ie, tr1, tr2, kk,i
    integer(i4)     :: iv, v(6)
    integer(i4)     :: icell, cell_index, inb, nnb
    integer(i4)     :: my_edge_index, ur_cell_index

    if(mpi_rank() == 0) print*,"global init_geometric_info"
    do ie = 1, mesh%ne
       tr1   = mesh%edp(ie)%v(1)
       tr2   = mesh%edp(ie)%v(2)
       v(1)  = mesh%tri(tr1)%v(1)
       v(2)  = mesh%tri(tr1)%v(2)
       v(3)  = mesh%tri(tr1)%v(3)
       v(4)  = mesh%tri(tr2)%v(1)
       v(5)  = mesh%tri(tr2)%v(2)
       v(6)  = mesh%tri(tr2)%v(3)
       kk    = 1
       do iv = 1, 6
          if(v(iv).ne.mesh%edt(ie)%v(1).and.v(iv).ne.mesh%edt(ie)%v(2))then
             mesh%edp(ie)%vtx(kk) = v(iv)
             kk    = kk +1
          end if
       end do
       if(kk.ne.3)then
          print*,"error in config geometric info I, stop with kk=", kk,ie,v
          stop
       end if
    end do

    ! call config_index_for_o58(mesh)

    DO ie = 1, mesh%ne
       DO icell = 1, 2
          cell_index = mesh%edt(ie)%v(icell)
          nnb        = mesh%vtx(cell_index)%nnb
! this sequence follows that in computing the B matrix
          mesh%edt(ie)%my_edge_on_edge(icell,1:nnb) = mesh%vtx(cell_index)%ed(1:nnb)
          mesh%edt(ie)%ur_cell_on_edge(icell,1:nnb) = mesh%vtx(cell_index)%nb(1:nnb)
          mesh%edt(ie)%my_edge_on_edge_num(icell)   = nnb
          mesh%edt(ie)%ur_cell_on_edge_num(icell)   = nnb
! for each my_edge-ur_cell pair, change ur_cell_on_edge index from global
! to local in my_edge's two ends
          do inb = 1, nnb
             my_edge_index = mesh%edt(ie)%my_edge_on_edge(icell,inb)
             ur_cell_index = mesh%edt(ie)%ur_cell_on_edge(icell,inb)
             if(mesh%edt(my_edge_index)%v(1).eq.ur_cell_index)then
                mesh%edt(ie)%ur_cell_on_edge(icell,inb) = 1
             else if(mesh%edt(my_edge_index)%v(2).eq.ur_cell_index)then
                mesh%edt(ie)%ur_cell_on_edge(icell,inb) = 2
             else
                print*,"this should not happen, stop"
                stop
             end if
          end do
       END DO
    END DO

    return
   end subroutine init_geometric_info_global
  
  
   subroutine init_geometric_info(mesh)
! io
    type(global_domain), intent(inout) :: mesh
    if(mpi_rank() == 0) print*,"local init_geometric_info"
    call config_edp_vtx(mesh)  !unuse and edp%vtx is global index
#ifndef USE_HALO2
    call config_index_for_o58(mesh)
#endif
    call config_kite_area(mesh)
    return
   end subroutine init_geometric_info

!================================================
!  PRIVATE
!================================================
 
  subroutine config_edp_vtx(mesh)
! io
   type(global_domain), intent(inout) :: mesh
! local
   integer(i4)           :: ie, tr1, tr2, kk, i, g_index
   integer(i4)           :: iv, v(6)
   type(scalar_2d_field) :: tmp_data_exchange

   mesh%ne = mesh%ne_compute
   do ie = 1, mesh%ne
      tr1   = mesh%edp(ie)%v(1)
      tr2   = mesh%edp(ie)%v(2)
      v(1)  = mesh%tri(tr1)%v(1)
      v(2)  = mesh%tri(tr1)%v(2)
      v(3)  = mesh%tri(tr1)%v(3)
      v(4)  = mesh%tri(tr2)%v(1)
      v(5)  = mesh%tri(tr2)%v(2)
      v(6)  = mesh%tri(tr2)%v(3)
      kk    = 1
      do iv = 1, 6
         if(v(iv).ne.mesh%edt(ie)%v(1).and.v(iv).ne.mesh%edt(ie)%v(2))then
            mesh%edp(ie)%vtx(kk) = v(iv)
            kk = kk +1
         end if
      end do
      if(kk.ne.3)then
         print*,"error in config geometric info I, stop with kk=", kk,ie,v
         stop
      end if
   end do
   mesh%ne = mesh%ne_full

   allocate(tmp_data_exchange%f(2,mesh%ne_full))
   do ie=1,mesh%ne_compute
      do i=1,2
        g_index = mesh%v_index(mesh%edp(ie)%vtx(i))
        tmp_data_exchange%f(i,ie)=dble(g_index)
      end do
   end do

   field_head_2d =>null()

   call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
   call exchange_data_2d(mesh%local_block,field_head_2d)

   do ie=mesh%ne_compute+1,mesh%ne_full
     do i=1,2
       g_index = int(tmp_data_exchange%f(i,ie))
       mesh%edp(ie)%vtx(i) = mesh%map_g2l_v%find(g_index)
     end do
   end do
   deallocate(tmp_data_exchange%f)

   return
  end subroutine config_edp_vtx

! o58 needs two special cell-edge relations to speed up runtime effciency

   subroutine config_index_for_o58(mesh)
! io
   type(global_domain),  intent(inout) :: mesh
! local
   integer(i4)                         :: ie, icell, cell_index, inb, nnb,i,j,k
   integer(i4)                         :: my_edge_index, ur_cell_index,g_index
   type(scalar_2d_field)               :: tmp_data_exchange_my,tmp_data_exchange_ur,&
                                           tmp_data_exchange_my_n,tmp_data_exchange_ur_n
! added
   integer(i4)                         :: iedge, edge_gindex, kk
   real(r8)                            :: p1(3), p2(3), p3(3), ccw_vector(3), det_tmp

   mesh%ne = mesh%ne_compute
   DO ie = 1, mesh%ne
      DO icell = 1, 2
         cell_index = mesh%edt(ie)%v(icell)
         nnb        = mesh%vtx(cell_index)%nnb
! this sequence follows that in computing the B matrix
         mesh%edt(ie)%my_edge_on_edge(icell,1:nnb) = mesh%vtx(cell_index)%ed(1:nnb)
         mesh%edt(ie)%my_edge_on_edge(icell,nnb+1:8) = -1
         mesh%edt(ie)%ur_cell_on_edge(icell,1:nnb) = mesh%vtx(cell_index)%nb(1:nnb)
         mesh%edt(ie)%ur_cell_on_edge(icell,nnb+1:8) = -1
         mesh%edt(ie)%my_edge_on_edge_num(icell)   = nnb
         mesh%edt(ie)%ur_cell_on_edge_num(icell)   = nnb
! for each my_edge-ur_cell pair, change ur_cell_on_edge index from global
! to local in my_edge's two ends
         do inb = 1, nnb
            my_edge_index = mesh%edt(ie)%my_edge_on_edge(icell,inb)
            ur_cell_index = mesh%edt(ie)%ur_cell_on_edge(icell,inb)
            if(mesh%edt(my_edge_index)%v(1).eq.ur_cell_index)then
               mesh%edt(ie)%ur_cell_on_edge(icell,inb) = 1
            else if(mesh%edt(my_edge_index)%v(2).eq.ur_cell_index)then
               mesh%edt(ie)%ur_cell_on_edge(icell,inb) = 2
            else
               print*,"this should not happen, stop"
               stop
            end if
         end do
      END DO
   END DO
   mesh%ne = mesh%ne_full
    
    !exchange halo data
   allocate(tmp_data_exchange_my%f(2*8,mesh%ne_full))
   allocate(tmp_data_exchange_ur%f(2*8,mesh%ne_full))
   allocate(tmp_data_exchange_ur_n%f(2,mesh%ne_full))
   allocate(tmp_data_exchange_my_n%f(2,mesh%ne_full))
   tmp_data_exchange_my%f = -1
   tmp_data_exchange_ur%f = -1
   do ie=1,mesh%ne_compute
      do k=1,8
         do j=1,2
            if(mesh%edt(ie)%my_edge_on_edge(j,k) <= mesh%ne_full .and. &
               mesh%edt(ie)%my_edge_on_edge(j,k) > 0) then
               g_index = mesh%e_index(mesh%edt(ie)%my_edge_on_edge(j,k))
               tmp_data_exchange_my%f(j+(k-1)*2,ie)=dble(g_index)
            end if
            tmp_data_exchange_ur%f(j+(k-1)*2,ie)=dble(mesh%edt(ie)%ur_cell_on_edge(j,k))
         end do
      end do
      do j=1,2
         tmp_data_exchange_my_n%f(j,ie)=dble(mesh%edt(ie)%my_edge_on_edge_num(j))
         tmp_data_exchange_ur_n%f(j,ie)=dble(mesh%edt(ie)%ur_cell_on_edge_num(j))
      end do
   end do

   call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_my)
   call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_ur)
   call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_my_n)
   call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_ur_n)
   call exchange_data_2d(mesh%local_block,field_head_2d)
    
   do ie=mesh%ne_compute+1,mesh%ne_full
      do k=1,8
         do j=1,2
            g_index = int(tmp_data_exchange_my%f(j+(k-1)*2,ie))
            mesh%edt(ie)%my_edge_on_edge(j,k) = mesh%map_g2l_e%find(g_index)
            mesh%edt(ie)%ur_cell_on_edge(j,k) = int(tmp_data_exchange_ur%f(j+(k-1)*2,ie))
         end do
      end do
      do j=1,2
         mesh%edt(ie)%my_edge_on_edge_num(j) = int(tmp_data_exchange_my_n%f(j,ie))
         mesh%edt(ie)%ur_cell_on_edge_num(j) = int(tmp_data_exchange_ur_n%f(j,ie))
      end do
   end do
   deallocate(tmp_data_exchange_my%f)
   deallocate(tmp_data_exchange_ur%f)
   deallocate(tmp_data_exchange_my_n%f)
   deallocate(tmp_data_exchange_ur_n%f)

   allocate(mesh%edt_my_edge_on_edge(2,8,mesh%ne_full))
   allocate(mesh%edt_ur_cell_on_edge(2,8,mesh%ne_full))
   do ie = 1, mesh%ne_full
      mesh%edt_my_edge_on_edge(1:2,1:8,ie) = mesh%edt(ie)%my_edge_on_edge(1:2,1:8)
      mesh%edt_ur_cell_on_edge(1:2,1:8,ie) = mesh%edt(ie)%ur_cell_on_edge(1:2,1:8)
   end do

   return
   end subroutine config_index_for_o58

   subroutine config_kite_area(mesh)
! io
    type(global_domain), intent(inout) :: mesh
! local
    type(scalar_2d_field)              :: tmp_tri_kite_area
    integer(i4)  :: it

#ifdef LAM_DOMAIN
    mesh%nt = mesh%nt_halo(stencil_width-1)
#else
    mesh%nt = mesh%nt_compute
#endif
    call tri_kite_areas(mesh)
    mesh%nt = mesh%nt_full

#ifndef LAM_DOMAIN
    allocate(tmp_tri_kite_area%f(3,mesh%nt_full))
    do it=1,mesh%nt_compute
       tmp_tri_kite_area%f(1:3,it) = mesh%tri(it)%kite_area(1:3)
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_tri_kite_area)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do it=mesh%nt_compute+1,mesh%nt_full
       mesh%tri(it)%kite_area(1:3) = tmp_tri_kite_area%f(1:3,it)
    end do
    deallocate(tmp_tri_kite_area%f)
#endif
    end subroutine config_kite_area

    subroutine config_plg_tri_nrtg_corrector(mesh)
!
! Because this part is currently used for directly replace 'read mesh file', it better to call
! it before mesh weight calculation, imediately after reading mesh, to please all codes; in future,
! some further optimization can be done
!

! io
    type(global_domain), intent(inout) :: mesh
! local
    type(scalar_2d_field)              :: tmp1, tmp2
    real(r8)                           :: pp0(3), pp1(3), normal1(3), real_flag
    integer(i4)                        :: iv, it, ie, ie_local, flag, inb

!
! Since tg corrector is not used, we need only to define nr here

!
!                  local_v1
!!                 /\  
!            ed1  /  \
!           p1   /    \  p1
!               /  p0  \
!              /        \
!            v2 ---------v3

       allocate(tmp1%f(nDual, mesh%nt_full))
       tmp1%f  =  -999

#ifdef LAM_DOMAIN
       do it = 1, mesh%nt_full
#else
       do it = 1, mesh%nt_compute
#endif
          pp0     = mesh%tri(it)%c%p(1:3)
          if(mesh%tri(it)%nnb.ne.nDual)then
            call endrun('tri_nnb is not equal to nDual')
          end if
          do ie_local = 1, mesh%tri(it)%nnb       ! loop over nb tr, ith nbtr coonects ith local edge
             ie      = mesh%tri(it)%ed(ie_local)
#ifdef LAM_DOMAIN
             if (ie .gt. 0 .and. mesh%tri(it)%nb(ie_local) .gt. 0) then
#endif
             normal1 = mesh%edt(ie)%nr(1:3)
             pp1     = mesh%tri(mesh%tri(it)%nb(ie_local))%c%p(1:3)
!
! because p0->p1 is not always perpendicular to mesh%tri_ed (e.g., cubed-sphere), so require dot_product.ge.0
!
             if(abs(dot_product((pp1-pp0),(pp0+pp1))).gt.eps*10)then ! must be close to perp
                print*,"vector ill in calc mesh%tri_nr, check",dot_product((pp1-pp0),(pp0+pp1)),eps
                call mpi_abort()
             end if
             flag = 999
             if(dot_product((pp1-pp0),normal1).gt.eps*100)then
                mesh%tri(it)%nr(ie_local) = 1
                flag = 1
             end if
             if(dot_product((pp1-pp0),normal1).lt.-eps)then
                mesh%tri(it)%nr(ie_local) = -1
                flag = 1
             end if
             if(flag.eq.999)then
                call endrun("bad in config_plg_tri_nrtg_corrector")
             end if
             tmp1%f(ie_local,it) = mesh%tri(it)%nr(ie_local)
#ifdef LAM_DOMAIN
             end if
#endif
          end do
       end do

!
!            nb
!
!       tr2--ed2--tr1
!        |       |
!    nb  |  vtx  | ed1  nb
!        |       |
!        tr3------tr4
!
!            nb
       allocate(tmp2%f(maxnnb,mesh%nv_full))
       tmp2%f = -999

#ifdef LAM_DOMAIN
       do iv = 1, mesh%nv_full
#else
       do iv = 1, mesh%nv_compute
#endif
          allocate(mesh%plg(iv)%nr(1:mesh%vtx(iv)%nnb))
          do inb = 1, mesh%vtx(iv)%nnb
! configure mesh%plg_nr
#ifdef LAM_DOMAIN
             if (mesh%vtx(iv)%nb(inb) .gt. 0) then
#endif
             real_flag = dot_product(mesh%edp(mesh%vtx(iv)%ed(inb))%nr(1:3),(mesh%vtx(mesh%vtx(iv)%nb(inb))%p(1:3)-mesh%vtx(iv)%p(1:3)))
             flag = 999
             if(real_flag.gt.eps*100)then
                mesh%plg(iv)%nr(inb) =  1
                flag = 1
             end if
             if(real_flag.lt.-eps*100)then
                mesh%plg(iv)%nr(inb) = -1
                flag = 1
             end if
             if(flag.eq.999)then
                call endrun("mesh%edp_nr is nearly perpendicular to iv->inb")
                stop
             end if
             tmp2%f(inb,iv) = mesh%plg(iv)%nr(inb)
#ifdef LAM_DOMAIN
             end if
#endif
          end do
       end do

    field_head_2d =>null()
    call exchange_data_2d_add(mesh,field_head_2d,tmp1)
    call exchange_data_2d_add(mesh,field_head_2d,tmp2)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do it = mesh%nt_compute+1, mesh%nt_full
      do inb =1, nDual
         mesh%tri(it)%nr(inb) = int(tmp1%f(inb,it))
      end do
    end do

    do iv = mesh%nv_compute+1, mesh%nv_full
      if(.not. allocated(mesh%plg(iv)%nr)) allocate(mesh%plg(iv)%nr(1:mesh%vtx(iv)%nnb))
      do inb =1, mesh%vtx(iv)%nnb
         mesh%plg(iv)%nr(inb) = int(tmp2%f(inb,iv))
      end do
    end do

    deallocate(tmp1%f)
    deallocate(tmp2%f)

    return
   end subroutine config_plg_tri_nrtg_corrector

 end module grist_geometric_info_icosh
