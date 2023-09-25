
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: Calculate the mesh weight for the entire mesh
!              Kinds of weights:
!              1) mapping weights for vector recons
!              2) The B weight for polynomial recons
!              3) scaling factor for VR (hyperdiffusion)
! Revision history: 
!              1) Modified for LAM_Domain (same results for global modeling)
!----------------------------------------------------------------------------

 module grist_mesh_weight_icosh

  !use grist_constants,      only: i4, r8, rearth
  use grist_constants_dbl,   only: i4, r8, rearth, zero, one, two, pi, deg2rad, half, r4
  use grist_nml_module,      only: ref_leng, vr_mode, tracer_hadv_flag, do_dycore_laplacian_4th, do_dycore_laplacian_6th
  use grist_domain_types  ,  only: global_domain
  use grist_math_module,     only: plg_kite_areas, norm, cross_product, sphtriarea, arcdistll
  use grist_svd_module,      only: pseudo_inverse, svd_lapack, migs
  use grist_handle_error,    only: endrun
  use grist_list,            only: index_list, check_and_add
  use grist_mpi
  use grist_data_types,      only: exchange_field_list_2d, scalar_2d_field
  use grist_config_partition,only: exchange_data_2d_add, exchange_data_2d
  use grist_geometric_info_icosh, only: config_plg_tri_nrtg_corrector

  implicit none

  private

  public :: calc_trsk_weights_plg      ,&   ! trsk weight at plg cell, TRSK09
            init_blocal_prime_cell     ,&   ! 2nd order recons, hx
            init_blocal_prime_cell_4th ,&   ! 4th order recons, hx
            init_blocal_dual_cell      ,&   ! 2nd order recons, tr (default now)
            init_blocal_dual_cell_4th  ,&   ! 4th order recons, tr, bad not use
            init_blocal_dual_cell_mbs  ,&   ! 2nd order recons, tr, more balanced stencil
            init_blocal_dual_cell_4th_mbs, &! 4th order recons, tr, more balanced stencil
            init_plg_sub_triangle_info ,&   ! sub-triangle info, for ffsl
            project_sphere2tplane      ,&
            fill_key_compute_pattern2  ,&
! Raw global seq-version
            calc_trsk_weights_plg_global     ,&
            init_blocal_prime_cell_global    ,&
            init_plg_sub_triangle_info_global,&
            init_blocal_prime_cell_4th_global,&
            init_blocal_dual_cell_mbs_global ,&
            init_blocal_dual_cell_global     ,&
            init_blocal_dual_cell_4th_mbs_global,&
            init_blocal_dual_cell_4th_global

  type(exchange_field_list_2d), pointer :: field_head_2d
  interface project_sphere2tplane 
#ifdef SPCODE
      module procedure project_sphere2tplane_r4
#else
      module procedure project_sphere2tplane_r4, project_sphere2tplane_r8
#endif
  end interface project_sphere2tplane

  contains

!================================================================
!   Calculates weights for each edge of prime cell i, weight is 
! stored in a 2D array, in which the 1st D denotes e (pushed), 
! and the 2nd D denotes e'(push)
!================================================================

  subroutine calc_trsk_weights_plg_global(mesh)
! io
    type(global_domain), intent(inout) :: mesh
! local
!Index for edges
    integer(i4)                         :: iv, j, kk
    integer(i4)                         :: ie
    integer(i4)                         :: iee      ! e  in wee'
    integer(i4)                         :: iep      ! e' in wee'
    integer(i4)                         :: tev2     ! see TRSK's EQ(33), left side
    integer(i4)                         :: nei      ! see TRSK's EQ(33), right side
    integer(i4)                         :: signcor  ! signcor=nei*tev2
    integer(i4)                         :: global_index_of_trv2
    integer(i4)                         :: global_index_of_iee
    integer(i4)                         :: local_index_iee_trv2
!Area ratios
    real(r8), allocatable               :: riv(:)
    real(r8)                            :: riv_sum
    real(r8)                            :: tmp_sum
    real(r8), parameter                 :: alpha=0.5
    real(r8)                            :: tmp
    integer(i4)                         :: icell, v_index, e_index
!
! calculate intersecting areas between dual and prime cells
!
    if(mpi_rank() == 0) print *,"global calc_trsk_weights_plg"
    do iv=1, mesh%nv
       if(.not. allocated(mesh%plg(iv)%kite_area)) allocate(mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb))
       mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb) = plg_kite_areas(iv, mesh%vtx(iv)%nnb, mesh)
       mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb) = mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb)/mesh%plg(iv)%areag
    end do

    do ie = 1, mesh%ne

       tmp                      = sign(1._r8,dot_product(mesh%edp(ie)%tg,mesh%edt(ie)%nr))
       mesh%edt(ie)%edpTg_edtNr   = int(tmp)
       mesh%edp(ie)%edpTg_edtNr   = int(tmp)

       if(abs(mesh%edt(ie)%edpTg_edtNr).ne.1) print*,&
         "bad news #1 in calc_trsk_weights, mesh%edt(ie)%edpTg_edtNr=",mesh%edt(ie)%edpTg_edtNr

    end do

    ! TRSK weights
    do iv=1, mesh%nv

       if(.not. allocated(mesh%plg(iv)%trsk_weight)) allocate(mesh%plg(iv)%trsk_weight(1:mesh%vtx(iv)%nnb,1:mesh%vtx(iv)%nnb))

       if(.not. allocated(riv)) allocate(riv(1:mesh%vtx(iv)%nnb))

       riv(1:mesh%vtx(iv)%nnb)= mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb)
       tmp_sum              = sum(riv(1:mesh%vtx(iv)%nnb))

!================================================
!               Calculating TRSKW
!================================================

!
! for each edge e
!
       do iee  = 1, mesh%vtx(iv)%nnb
!
! for each edge e'
!
         do iep  = 1, mesh%vtx(iv)%nnb
            if(iep.eq.iee)then ! zero weight for itself
               mesh%plg(iv)%trsk_weight(iee,iep)=0._r8
               cycle
            end if

            if(iep.gt.iee)then  ! e.g., w(2,5)
               riv_sum = 0._r8
               do j = iep, mesh%vtx(iv)%nnb
                  riv_sum = riv_sum + riv(j)
               end do
               do j = 1, iee-1  ! if iee-1=0, this loop works nothing.
                  riv_sum = riv_sum + riv(j)
               end do
            end if

            if(iep.lt.iee)then  ! e.g., w(5,2)
               riv_sum = 0._r8
               do j=iep, iee-1
                  riv_sum = riv_sum+riv(j)
               end do
            end if
            
            mesh%plg(iv)%trsk_weight(iee,iep) = riv_sum-alpha ! without sign yet
!
! the binary flag of hx(iv)'s e'
!
            nei = mesh%plg(iv)%nr(iep)

!================================================
! tev2 belongs to triangle at dual cell v2
! so we need to know:
!       (1) dcell_structure at v2
!       (2) local index of iee within triangle v2
!           the number of v2=iee-1
!================================================

            if(iee-1.ne.0)then
               global_index_of_trv2 = mesh%vtx(iv)%tr(iee-1)
            else  ! = 0
               global_index_of_trv2 = mesh%vtx(iv)%tr(mesh%vtx(iv)%nnb)
            end if

            !triangle_v2_structure = mesh%tri(global_index_of_trv2) ! (1) gotta 

            global_index_of_iee   = mesh%vtx(iv)%ed(iee)
!
! find local index within the triangle_v2 for iee
!
            local_index_iee_trv2  = -999

            do kk = 1, 3
               if(mesh%tri(global_index_of_trv2)%ed(kk).eq.global_index_of_iee)then
                  local_index_iee_trv2 = kk
               end if
            end do
!
! Thanks, we don't have bad news :)
!
            if(local_index_iee_trv2.eq.-999) call endrun("bad news #2 in calc_trsk")

            tev2    = mesh%edt(global_index_of_iee)%edpTg_edtNr*mesh%tri(global_index_of_trv2)%nr(local_index_iee_trv2)
            signcor = tev2*nei

            if(abs(signcor).ne.1) call endrun("bad news #3 in calc_trsk")

            mesh%plg(iv)%trsk_weight(iee,iep)=signcor*(riv_sum-alpha)

         end do
       end do
       deallocate(riv)
    end do

!
! for runtime efficiency, we restore trsk weight based on edge,
! also stored is the global edge index, their sequence should be same
! this is similar to that we store the blocal weight on an edge basis
! both ie and e_index belongs to this cell
! e_index's local index relative to this cell is inb
! then we need to know ie's local index relative to this cell
! then what we need is a trsk weight that e_index->ie
! i.e., iep=>e_index & iee=>ie
!
    do ie = 1, mesh%ne
       kk = 1
       do icell = 1, 2  ! two cells share edge ie
          v_index = mesh%edt(ie)%v(icell)
          do iee = 1, mesh%vtx(v_index)%nnb
             if(mesh%vtx(v_index)%ed(iee).eq.ie)then  ! only once
                do iep = 1, mesh%vtx(v_index)%nnb     ! loop all nb of this cell
                   e_index = mesh%vtx(v_index)%ed(iep)
                   mesh%edp(ie)%edge_on_edge(kk) = e_index
                   mesh%edp(ie)%trsk_on_edge(kk) = mesh%plg(v_index)%trsk_weight(iee,iep)
                   mesh%edp(ie)%edpl_on_edge(kk) = rearth*mesh%edp(e_index)%leng
                   kk = kk+1
                end do
             end if
          end do
       end do
       mesh%edp(ie)%nedge = kk-1
    end do

   return
  end subroutine calc_trsk_weights_plg_global
  
  subroutine calc_trsk_weights_plg(mesh)
! io
    type(global_domain), intent(inout) :: mesh
! local
!Index for edges
    integer(i4)                         :: iv, j, kk
    integer(i4)                         :: ie, it
    integer(i4)                         :: iee      ! e  in wee'
    integer(i4)                         :: iep      ! e' in wee'
    integer(i4)                         :: tev2     ! see TRSK's EQ(33), left side
    integer(i4)                         :: nei      ! see TRSK's EQ(33), right side
    integer(i4)                         :: signcor  ! signcor=nei*tev2
    integer(i4)                         :: global_index_of_trv2
    integer(i4)                         :: global_index_of_iee
    integer(i4)                         :: local_index_iee_trv2
!Area ratios
    real(r8), allocatable               :: riv(:)
    real(r8)                            :: riv_sum
    real(r8)                            :: tmp_sum
    real(r8), parameter                 :: alpha=0.5
    real(r8)                            :: tmp
    integer(i4)                         :: icell, v_index, e_index

    type(scalar_2d_field)               :: tmp_data_exchange_edge,tmp_data_exchange_edpl,&
                                           tmp_data_exchange_trsk,tmp_data_exchange_nedge
    integer(i4)                         :: g_index,k
    integer(i4)                         :: v1, v2
    real(r8),dimension(1:3)             :: v1v2            ! direction 1->2

    call config_plg_tri_nrtg_corrector(mesh)
!
! calculate intersecting areas between dual and prime cells
!
    if(mpi_rank() == 0) print *,"local calc_trsk_weights_plg"
    do iv=1, mesh%nv
       if(.not. allocated(mesh%plg(iv)%kite_area)) allocate(mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb))
       mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb) = plg_kite_areas(iv, mesh%vtx(iv)%nnb, mesh)
       mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb) = mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb)/mesh%plg(iv)%areag
    end do

    do ie = 1, mesh%ne

       tmp                      = sign(1._r8,dot_product(mesh%edp(ie)%tg,mesh%edt(ie)%nr))
       mesh%edt(ie)%edpTg_edtNr   = int(tmp)
       mesh%edp(ie)%edpTg_edtNr   = int(tmp)

       if(abs(mesh%edt(ie)%edpTg_edtNr).ne.1) print*,&
         "bad news #1 in calc_trsk_weights, mesh%edt(ie)%edpTg_edtNr=",mesh%edt(ie)%edpTg_edtNr

       tmp                      = sign(1._r8,dot_product(mesh%edp(ie)%nr,mesh%edt(ie)%tg))
       mesh%edt(ie)%edpNr_edtTg   = int(tmp)
       mesh%edp(ie)%edpNr_edtTg   = int(tmp)

       if(abs(mesh%edt(ie)%edpNr_edtTg).ne.1) print*,&
         "bad news #1 in calc_trsk_weights, mesh%edt(ie)%edpNr_edtTg=",mesh%edt(ie)%edpNr_edtTg
    end do

    do ie = 1, mesh%ne

       v1        = mesh%edp(ie)%v(1)
       v2        = mesh%edp(ie)%v(2)

       v1v2      = mesh%tri(v2)%c%p-mesh%tri(v1)%c%p

       tmp       = sign(1._r8,dot_product(v1v2,mesh%edp(ie)%tg))
       mesh%edt(ie)%edp12_edp_tg = int(tmp)
       mesh%edp(ie)%edp12_edp_tg = int(tmp)

       if(abs(mesh%edt(ie)%edp12_edp_tg).ne.1) print*,&
         "bad news #1 in calc_trsk_weights, mesh%edt(ie)%edp12_edp_tg=",mesh%edt(ie)%edp12_edp_tg

       tmp       = sign(1._r8,dot_product(v1v2,mesh%edt(ie)%nr))
       mesh%edt(ie)%edp12_edt_nr = int(tmp)
       mesh%edp(ie)%edp12_edt_nr = int(tmp)

       if(abs(mesh%edt(ie)%edp12_edt_nr).ne.1) print*,&
         "bad news #1 in calc_trsk_weights, mesh%edt(ie)%edp12_edt_nr=",mesh%edt(ie)%edp12_edt_nr
    end do

    ! TRSK weights
    do iv=1, mesh%nv
       if(.not. allocated(mesh%plg(iv)%trsk_weight)) allocate(mesh%plg(iv)%trsk_weight(1:mesh%vtx(iv)%nnb,1:mesh%vtx(iv)%nnb))
       if(.not. allocated(riv)) allocate(riv(1:mesh%vtx(iv)%nnb))
#ifdef LAM_DOMAIN
#ifdef REG_QJLAM
       if(mesh%plg_mark(iv).eq.177)then
#else
       if(mesh%plg_mark(iv).eq.177 .and.iv.gt.mesh%nv_halo(1))then
#endif
          mesh%plg(iv)%trsk_weight(:,:) = zero
          cycle
       end if
#endif
       riv(1:mesh%vtx(iv)%nnb)= mesh%plg(iv)%kite_area(1:mesh%vtx(iv)%nnb)
       tmp_sum              = sum(riv(1:mesh%vtx(iv)%nnb))

!================================================
!             Calculating TRSKW
!================================================

!
! for each edge e
!
       do iee  = 1, mesh%vtx(iv)%nnb
!
! for each edge e'
!
         do iep  = 1, mesh%vtx(iv)%nnb
            if(iep.eq.iee)then ! zero weight for itself
               mesh%plg(iv)%trsk_weight(iee,iep)=0._r8
               cycle
            end if
!
! key to remember, who-1, who is e
!
            if(iep.gt.iee)then  ! e.g., w(2,5)
               riv_sum = 0._r8
               do j = iep, mesh%vtx(iv)%nnb
                  riv_sum = riv_sum + riv(j)
               end do
               do j = 1, iee-1  ! if iee-1=0, this loop works nothing.
                  riv_sum = riv_sum + riv(j)
               end do
            end if

            if(iep.lt.iee)then  ! e.g., w(5,2)
               riv_sum = 0._r8
               do j=iep, iee-1
                  riv_sum = riv_sum+riv(j)
               end do
            end if
            
            mesh%plg(iv)%trsk_weight(iee,iep) = riv_sum-alpha ! without sign yet
!
! the binary flag of hx(iv)'s e'
!
            nei     = mesh%plg(iv)%nr(iep)

!================================================
! tev2 belongs to triangle at dual cell v2
! so we need to know:
! (1) dcell_structure at v2
! (2) local index of iee within triangle v2; the number of v2=iee-1
!================================================

            if(iee-1.ne.0)then
               global_index_of_trv2 = mesh%vtx(iv)%tr(iee-1)
            else  ! = 0
               global_index_of_trv2 = mesh%vtx(iv)%tr(mesh%vtx(iv)%nnb)
            end if

            !triangle_v2_structure = mesh%tri(global_index_of_trv2) ! (1) gotta 

            global_index_of_iee   = mesh%vtx(iv)%ed(iee)
!
! find local index within the triangle_v2 for iee
!
            local_index_iee_trv2  = -999

            do kk = 1, 3
               if(mesh%tri(global_index_of_trv2)%ed(kk).eq.global_index_of_iee)then
                  local_index_iee_trv2 = kk
               end if
            end do
!
! Thanks, we don't have bad news :)
!
            if(local_index_iee_trv2.eq.-999) call endrun("bad news #2 in calc_trsk")

            tev2    = mesh%edt(global_index_of_iee)%edpTg_edtNr*mesh%tri(global_index_of_trv2)%nr(local_index_iee_trv2)

            signcor = tev2*nei
            if(abs(signcor).ne.1) call endrun("bad news #3 in calc_trsk")

            mesh%plg(iv)%trsk_weight(iee,iep) = signcor*(riv_sum-alpha)

         end do
       end do
       deallocate(riv)
    end do
!
! for runtime efficiency, we restore trsk weight based on edge,
! also stored is the global edge index, their sequence should be same
! this is similar to that we store the blocal weight on an edge basis
! both ie and e_index belongs to this cell
! e_index's local index relative to this cell is inb
! then we need to know ie's local index relative to this cell
! then what we need is a trsk weight that e_index->ie
! i.e., iep=>e_index & iee=>ie
!
    allocate(tmp_data_exchange_edge%f(16,mesh%ne_full), source=-1._r8)
    allocate(tmp_data_exchange_edpl%f(16,mesh%ne_full))
    allocate(tmp_data_exchange_trsk%f(16,mesh%ne_full))
    allocate(tmp_data_exchange_nedge%f(1,mesh%ne_full))

#ifndef LAM_DOMAIN
    mesh%ne = mesh%ne_compute
#endif
    do ie = 1, mesh%ne
       mesh%edp(ie)%edge_on_edge = -1
       mesh%edp(ie)%trsk_on_edge = 0._r8
       mesh%edp(ie)%edpl_on_edge = 0._r8
       kk = 1
       do icell = 1, 2  ! two cells share edge ie
          v_index = mesh%edt(ie)%v(icell)
#ifdef LAM_DOMAIN
          if (v_index .gt. 0) then
#endif
          do iee = 1, mesh%vtx(v_index)%nnb
             if(mesh%vtx(v_index)%ed(iee).eq.ie)then  ! only once
                do iep = 1, mesh%vtx(v_index)%nnb     ! loop all nb of this cell
                   e_index = mesh%vtx(v_index)%ed(iep)
                   mesh%edp(ie)%edge_on_edge(kk) = e_index
                   mesh%edp(ie)%trsk_on_edge(kk) = mesh%plg(v_index)%trsk_weight(iee,iep)
                   mesh%edp(ie)%edpl_on_edge(kk) = rearth*mesh%edp(e_index)%leng
                   kk = kk+1
                end do
             end if
          end do
#ifdef LAM_DOMAIN
          end if
#endif
       end do
       ! for halo exchange
       mesh%edp(ie)%nedge = kk-1
       do k = 1, 16
         if ((mesh%edp(ie)%edge_on_edge(k) .gt. 0) .and. (mesh%edp(ie)%edge_on_edge(k) .le. mesh%ne_full)) then
           g_index = mesh%e_index(mesh%edp(ie)%edge_on_edge(k))
           tmp_data_exchange_edge%f(k,ie) = g_index
         end if
         tmp_data_exchange_edpl%f(k,ie) = mesh%edp(ie)%edpl_on_edge(k)
         tmp_data_exchange_trsk%f(k,ie) = mesh%edp(ie)%trsk_on_edge(k)
       end do
       tmp_data_exchange_nedge%f(1,ie) = mesh%edp(ie)%nedge
!
! for gassmann nct info; see implemented code (not verified) at some time during 2018-2019
!
    end do

#ifndef LAM_DOMAIN
    mesh%ne = mesh%ne_full
#endif

    !exchange halo data
    field_head_2d=>null()
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_edge)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_edpl)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_trsk)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_nedge)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do ie=mesh%ne_compute+1,mesh%ne_full
      do k=1,16
        g_index = int(tmp_data_exchange_edge%f(k,ie))
        mesh%edp(ie)%edge_on_edge(k) = mesh%map_g2l_e%find(g_index)
        mesh%edp(ie)%edpl_on_edge(k) = tmp_data_exchange_edpl%f(k,ie)
        mesh%edp(ie)%trsk_on_edge(k) = tmp_data_exchange_trsk%f(k,ie)
      end do
      mesh%edp(ie)%nedge = int(tmp_data_exchange_nedge%f(1,ie))
    end do

    deallocate(tmp_data_exchange_edge%f)
    deallocate(tmp_data_exchange_edpl%f)
    deallocate(tmp_data_exchange_trsk%f)
    deallocate(tmp_data_exchange_nedge%f)
!
!exchange gassnct's halo data info
!
   call fill_key_compute_pattern1(mesh)
   
   return
  end subroutine calc_trsk_weights_plg

  subroutine fill_key_compute_pattern1(mesh)
! io
  type(global_domain), intent(inout) :: mesh
! local
  integer(i4)          :: it, ie, iv
  real(r8)             :: dist1, dist2, gama, lamda, rho1, rho2, expp
  real(r8)             :: minDens, minDens_global, tmp_sum
  integer(i4)          :: ierr,inb
!
! Coriolis
!
    !allocate(mesh%edp_trsk_on_edge(16,mesh%ne_full))
    !allocate(mesh%edp_edpl_on_edge(16,mesh%ne_full))
    allocate(mesh%edp_edge_on_edge(16,mesh%ne_full))
    allocate(mesh%edp_nedge(mesh%ne_full))
    allocate(mesh%edp_trsk_edpl_on_edge(16,mesh%ne_full))

    do ie = 1, mesh%ne_full
      !mesh%edp_trsk_on_edge(1:16, ie) = mesh%edp(ie)%trsk_on_edge(1:16)
      !mesh%edp_edpl_on_edge(1:16, ie) = mesh%edp(ie)%edpl_on_edge(1:16)
      mesh%edp_edge_on_edge(1:16, ie) = mesh%edp(ie)%edge_on_edge(1:16)
                   mesh%edp_nedge(ie) = mesh%edp(ie)%nedge
      mesh%edp_trsk_edpl_on_edge(1:16,ie) = mesh%edp(ie)%trsk_on_edge(1:16)*mesh%edp(ie)%edpl_on_edge(1:16)
    end do
!
! divergence
!
    allocate(mesh%vtx_nnb(mesh%nv_full))
    allocate(mesh%vtx_ed(mesh%maxvnb,mesh%nv_full))
    allocate(mesh%plg_nr(mesh%maxvnb,mesh%nv_full))
    allocate(mesh%plg_areag(mesh%nv_full))
    allocate(mesh%edp_leng(mesh%ne_full))
    do iv = 1, mesh%nv_full
      mesh%vtx_nnb(iv) = mesh%vtx(iv)%nnb
      mesh%vtx_ed(1:mesh%vtx_nnb(iv), iv) = mesh%vtx(iv)%ed(1:mesh%vtx_nnb(iv))
      mesh%plg_nr(1:mesh%vtx_nnb(iv), iv) = mesh%plg(iv)%nr(1:mesh%vtx_nnb(iv))
      mesh%plg_areag(iv) = mesh%plg(iv)%areag
    end do
    do ie = 1, mesh%ne_full
      mesh%edp_leng(ie) = mesh%edp(ie)%leng
    end do
!
! plg flux
!
    allocate(mesh%edt_v(2, mesh%ne_full))
    allocate(mesh%edt_edpNr_edtTg(mesh%ne_full))
    allocate(mesh%edt_edpTg_edtNr(mesh%ne_full))
    allocate(mesh%edt_leng(mesh%ne_full))
    do ie = 1, mesh%ne_full
      mesh%edt_v(1:2,ie)= mesh%edt(ie)%v(1:2)
      mesh%edt_edpNr_edtTg(ie) = mesh%edt(ie)%edpNr_edtTg
      mesh%edt_edpTg_edtNr(ie) = mesh%edt(ie)%edpTg_edtNr
      mesh%edt_leng(ie)      = mesh%edt(ie)%leng
    end do
!
! tri flux
!
    allocate(mesh%edp_v(2, mesh%ne_full))
    allocate(mesh%edp12_edp_tg(mesh%ne_full))
    allocate(mesh%edp12_edt_nr(mesh%ne_full))
    do ie = 1, mesh%ne_full
      mesh%edp_v(1:2,ie)    = mesh%edp(ie)%v(1:2)
      mesh%edp12_edp_tg(ie) = mesh%edp(ie)%edp12_edp_tg
      mesh%edp12_edt_nr(ie) = mesh%edp(ie)%edp12_edt_nr
    end do
!
! vorticity(curl)
!
    allocate(mesh%tri_nnb(mesh%nt_full))
    allocate(mesh%tri_ed(3,mesh%nt_full))
    allocate(mesh%tri_nr(3,mesh%nt_full))
    allocate(mesh%tri_areag(mesh%nt_full))
    allocate(mesh%tri_c_lat(mesh%nt_full))
    do it = 1, mesh%nt_full
      mesh%tri_nnb(it)     = mesh%tri(it)%nnb ! for scvt, this is always 3
      mesh%tri_ed(1:3, it) = mesh%tri(it)%ed(1:3)
      mesh%tri_nr(1:3, it) = mesh%tri(it)%nr(1:3)
      mesh%tri_areag(it)   = mesh%tri(it)%areag
      mesh%tri_c_lat(it)   = mesh%tri(it)%c%lat
    end do
!
! kinetic_energy
!
    allocate(mesh%vtx_tr(mesh%maxvnb,mesh%nv_full))
    allocate(mesh%plg_kite_area(mesh%maxvnb,mesh%nv_full))
    do iv = 1, mesh%nv_full
      mesh%vtx_tr(1:mesh%vtx_nnb(iv), iv)        = mesh%vtx(iv)%tr(1:mesh%vtx_nnb(iv))
      mesh%plg_kite_area(1:mesh%vtx_nnb(iv), iv) = mesh%plg(iv)%kite_area(1:mesh%vtx_nnb(iv))
    end do
!
!----------------zhangyi added for VR-----------------------------
!
!    gama  = (one/ccvt_gama)**4
!    lamda = (one/ccvt_lamda)**4

    if(.not.allocated(mesh%vtxCellLeng))allocate(mesh%vtxCellLeng(mesh%nv_full))
!    select case (density_function_type)
!    case (0)  ! QU mesh
!       mesh%vtxCellLeng = one
!       if(mpi_rank().eq.0) print*,"we use the uniform density"
!    case (1)  ! single-refinment mesh, Eq. (3)
!       if(mpi_rank().eq.0) print*,"we use the single refinment style centered at",center1_lon," ",center1_lat
!       do iv = 1, mesh%nv_full
!          dist1 = arcdistll(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,center1_lon*deg2rad,center1_lat*deg2rad)
!          mesh%vtxCellLeng(iv) = half/(one-gama)*(tanh((ccvt_beta*deg2rad-dist1)/(ccvt_alpha*deg2rad))+one)+gama
!       end do
!    case (2)  ! hiarachy-gradual style, Eq. (4)
!       if(mpi_rank().eq.0) print*,"we use the hiarachy-refinment mesh centered at",center1_lon," ",center1_lat
       do iv = 1, mesh%nv_full
          !dist1 = arcdistll(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,center1_lon*deg2rad,center1_lat*deg2rad)
          !mesh%vtxCellLeng(iv) = half/(one-gama)*((one -lamda)/(one-gama)*tanh((ccvt_beta1*deg2rad-dist1)/(ccvt_alpha1*deg2rad))+&
          !                                        (lamda-gama)/(one-gama)*tanh((ccvt_beta2*deg2rad-dist1)/(ccvt_alpha2*deg2rad))+one)+gama
          tmp_sum = zero
          do inb = 1, mesh%vtx(iv)%nnb
             tmp_sum = tmp_sum+mesh%edt(mesh%vtx(iv)%ed(inb))%leng
          end do
          mesh%vtxCellLeng(iv) = tmp_sum/mesh%vtx(iv)%nnb
       end do

!    case(3)   ! multi-center style, Eq. (7)
!       if(mpi_rank().eq.0) print*,"we use multi-center style centered at",center1_lon," ",center1_lat, center2_lon," ",center2_lat
!       do iv = 1, mesh%nv_full
!          dist1 = arcdistll(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,center1_lon*deg2rad,center1_lat*deg2rad)
!          dist2 = arcdistll(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,center2_lon*deg2rad,center2_lat*deg2rad)
!          mesh%vtxCellLeng(iv) = half/(one-gama)*(tanh((ccvt_beta*deg2rad-dist1)/(ccvt_alpha*deg2rad))+&
!                                                  tanh((ccvt_beta*deg2rad-dist2)/(ccvt_alpha*deg2rad))+two)+gama
!       end do

!    case(4)   ! special case for G8X4L2 in the VR paper
!       if(mpi_rank().eq.0) print*,"special case for G8X4L2 in the paper"
!       do iv = 1, mesh%nv_full
!          dist1 = arcdistll(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,center1_lon*deg2rad,center1_lat*deg2rad)
!          mesh%vtxCellLeng(iv) = half/(one-gama)*(0.94_r8*tanh((ccvt_beta1*deg2rad-dist1)/(ccvt_alpha1*deg2rad))+&
!                                                  0.06_r8*tanh((ccvt_beta2*deg2rad-dist2)/(ccvt_alpha2*deg2rad))+one)+gama
!       end do

!    case (5)  ! channel style
!       if(mpi_rank().eq.0) print*,"we use the channel refinment mesh centered at",center1_lon," ",center1_lat
!       do iv = 1, mesh%nv_full
!          dist1 = arcdistll(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,center1_lon*deg2rad,center1_lat*deg2rad)
!          rho1  = half/(one-gama)*(tanh((-ccvt_beta1*deg2rad+dist1)/(ccvt_alpha*deg2rad))+one)+gama
!          rho2  = half/(one-gama)*(tanh(( ccvt_beta2*deg2rad-dist1)/(ccvt_alpha*deg2rad))+one)+gama
!          mesh%vtxCellLeng(iv) = rho1*rho2
!       end do

!    case default
!       mesh%vtxCellLeng = one
!       if(mpi_rank().eq.0) print*,"default we use the uniform-density for scaling"
!    end select
    minDens = minval(mesh%vtxCellLeng(1:mesh%nv_compute))
!
! use mpi_reduce to obtain the global minimum mesh%vtxCellLeng
!
    call reduce(minDens, minDens_global, 'min')
    mesh%minCellLeng = minDens_global
    if(mpi_rank().eq.0) print*,"minCellLeng is ", mesh%minCellLeng*rearth, mesh%minCellLeng
!
! zhangyi added for VR: when we have this, we may deallocate mesh%vtxCellLeng
!
    !if(.not.allocated(mesh%edt_scale_2nd))allocate(mesh%edt_scale_2nd(mesh%ne_full))
    if(.not.allocated(mesh%edt_scale_4th).and.do_dycore_laplacian_4th) allocate(mesh%edt_scale_4th(mesh%ne_full))
    if(.not.allocated(mesh%edt_scale_6th).and.do_dycore_laplacian_6th) allocate(mesh%edt_scale_6th(mesh%ne_full))

    expp = 3.3219 ! empirically following Zarzychi et al. (2014)
!    select case (density_function_type)
!    case (0)  ! QU mesh
!       mesh%edt_scale_2nd = one
!       mesh%edt_scale_4th = one
!    case (1,2,3) ! VR mesh
     if(do_dycore_laplacian_4th)then
        do ie = 1, mesh%ne_halo(1)  ! halo1 is enough
          !mesh%edt_scale_2nd(ie) = one/(((mesh%vtxCellLeng(mesh%edt_v(1,ie))+mesh%vtxCellLeng(mesh%edt_v(2,ie)))*half )**0.25_r8)
          !mesh%edt_scale_4th(ie) = one/(((mesh%vtxCellLeng(mesh%edt_v(1,ie))+mesh%vtxCellLeng(mesh%edt_v(2,ie)))*half )**(0.25_r8*expp))
          !if(mesh%minCellLeng.ne.0)then
             !mesh%edt_scale_2nd(ie) =  ((mesh%vtxCellLeng(mesh%edt_v(1,ie))+mesh%vtxCellLeng(mesh%edt_v(2,ie)))*half )/(ref_leng/rearth)  ! not used
             mesh%edt_scale_4th(ie) = (((mesh%vtxCellLeng(mesh%edt_v(1,ie))+mesh%vtxCellLeng(mesh%edt_v(2,ie)))*half )/(ref_leng/rearth))**expp
          !else
          !  if(mpi_rank().eq.0) print*,"impossible case in minCellLeng"
          !end if
        end do
     end if
     if(do_dycore_laplacian_6th)then
        do ie = 1, mesh%ne_halo(1)  ! halo1 is enough
           mesh%edt_scale_6th(ie) = (((mesh%vtxCellLeng(mesh%edt_v(1,ie))+mesh%vtxCellLeng(mesh%edt_v(2,ie)))*half )/(ref_leng/rearth))**expp  ! same as for 4th
        end do
     end if
!    case default
!       mesh%edt_scale_2nd = one
!       mesh%edt_scale_4th = one
!    end select

    if(.not.vr_mode)then ! reset to one for QU
       !mesh%edt_scale_2nd = one
       if(do_dycore_laplacian_4th) mesh%edt_scale_4th = one
       if(do_dycore_laplacian_6th) mesh%edt_scale_6th = one
    end if

    return
  end subroutine fill_key_compute_pattern1
!
! added to replace all five elements, called after all mesh set
! once these are added, all five mesh elements can be deallocated;
!

  subroutine fill_key_compute_pattern2(mesh)
! io
  type(global_domain), intent(inout) :: mesh
! local
  integer(i4)          :: it, ie, iv

!#if (defined CP2_VTX) || (defined CP2_PLG) || (defined CP2_TRI) || (defined CP2_EDT) || (defined CP2_EDP)
!
! vtx
!
    if(.not.allocated(mesh%vtx_lon)) allocate(mesh%vtx_lon(mesh%nv_full))
    if(.not.allocated(mesh%vtx_lat)) allocate(mesh%vtx_lat(mesh%nv_full))
    if(.not.allocated(mesh%vtx_p))   allocate(mesh%vtx_p(3,mesh%nv_full))
    if(.not.allocated(mesh%vtx_nb))  allocate(mesh%vtx_nb(mesh%maxvnb,mesh%nv_full))
    do iv = 1, mesh%nv_full
       mesh%vtx_lon(iv)   = mesh%vtx(iv)%lon
       mesh%vtx_lat(iv)   = mesh%vtx(iv)%lat
       mesh%vtx_p(1:3,iv) = mesh%vtx(iv)%p(1:3)
       mesh%vtx_nb(1:mesh%vtx_nnb(iv),iv) = mesh%vtx(iv)%nb(1:mesh%vtx_nnb(iv))
    end do
!
! vtx, plg for ffsl
!
    IF(tracer_hadv_flag.eq.28)THEN ! only allocate for FFSL
    if(.not.allocated(mesh%vtx_lob_nx)) allocate(mesh%vtx_lob_nx(mesh%maxvnb,3,mesh%nv_full))
    if(.not.allocated(mesh%vtx_lob_ny)) allocate(mesh%vtx_lob_ny(mesh%maxvnb,3,mesh%nv_full))
    if(.not.allocated(mesh%plg_blocal))                allocate(mesh%plg_blocal(5,mesh%maxvnb,mesh%nv_full))                  ! (1)*n*stencil_number X nv
    if(.not.allocated(mesh%plg_sub_triangle_midp_lob)) allocate(mesh%plg_sub_triangle_midp_lob(mesh%maxvnb,3,2,mesh%nv_full)) ! nnb*3edges*2(lob_x,lob_y) X nplg
    if(.not.allocated(mesh%plg_sub_triangle_area))     allocate(mesh%plg_sub_triangle_area(    mesh%maxvnb,    mesh%nv_full)) ! nnb X nv
    do iv = 1, mesh%nv_full
       mesh%vtx_lob_nx(    1:mesh%vtx_nnb(iv),1:3,iv)                = mesh%vtx(iv)%lob_nx(1:mesh%vtx_nnb(iv),1:3)
       mesh%vtx_lob_ny(    1:mesh%vtx_nnb(iv),1:3,iv)                = mesh%vtx(iv)%lob_ny(1:mesh%vtx_nnb(iv),1:3)
       mesh%plg_blocal(1:5,1:mesh%vtx_nnb(iv),iv)                    = mesh%plg(iv)%blocal(1,1:5,1:mesh%vtx_nnb(iv))
       mesh%plg_sub_triangle_midp_lob(1:mesh%vtx_nnb(iv),1:3,1:2,iv) = mesh%plg(iv)%sub_triangle_midp_lob(1:mesh%vtx_nnb(iv),1:3,1:2)
       mesh%plg_sub_triangle_area(1:mesh%vtx_nnb(iv),iv)             = mesh%plg(iv)%sub_triangle_area    (1:mesh%vtx_nnb(iv))
    end do
    END IF
!
! tri
!
    if(.not.allocated(mesh%tri_c_p))      allocate(mesh%tri_c_p(3,mesh%nt_full))       ! 3, nt
    if(.not.allocated(mesh%tri_c_lon))    allocate(mesh%tri_c_lon(mesh%nt_full))       ! nt
    if(.not.allocated(mesh%tri_v))        allocate(mesh%tri_v(3,mesh%nt_full))         ! nb=3, nt
    if(.not.allocated(mesh%tri_kite_area))allocate(mesh%tri_kite_area(3,mesh%nt_full)) ! nb=3, nv
    do it = 1, mesh%nt_full
       mesh%tri_c_p(1:3,it)         = mesh%tri(it)%c%p(1:3)
       mesh%tri_c_lon(it)           = mesh%tri(it)%c%lon
       mesh%tri_v(1:3,it)           = mesh%tri(it)%v(1:3)
       mesh%tri_kite_area(1:3,it)   = mesh%tri(it)%kite_area(1:3)
    end do
!
! edt
!
    if(.not.allocated(mesh%edt_c_p))   allocate(mesh%edt_c_p(3,mesh%ne_full))
    if(.not.allocated(mesh%edt_c_lon)) allocate(mesh%edt_c_lon(mesh%ne_full))
    if(.not.allocated(mesh%edt_c_lat)) allocate(mesh%edt_c_lat(mesh%ne_full))
    if(.not.allocated(mesh%edt_nr))    allocate(mesh%edt_nr(3, mesh%ne_full))
    do ie = 1, mesh%ne_full
       mesh%edt_c_p(1:3,ie)  = mesh%edt(ie)%c%p(1:3)
       mesh%edt_c_lon(ie)    = mesh%edt(ie)%c%lon
       mesh%edt_c_lat(ie)    = mesh%edt(ie)%c%lat
       mesh%edt_nr( 1:3,ie)  = mesh%edt(ie)%nr
    end do
!
! edp
!
    if(.not.allocated(mesh%edp_c_p))   allocate(mesh%edp_c_p(3,mesh%ne_full))
    if(.not.allocated(mesh%edp_nr))    allocate(mesh%edp_nr(3, mesh%ne_full))
    if(.not.allocated(mesh%edp_tg))    allocate(mesh%edp_tg(3, mesh%ne_full))
    do ie = 1, mesh%ne_full
       mesh%edp_c_p(1:3,ie)  = mesh%edp(ie)%c%p(1:3)
       mesh%edp_nr( 1:3,ie)  = mesh%edp(ie)%nr
       mesh%edp_tg( 1:3,ie)  = mesh%edp(ie)%tg
    end do
!#endif

    return
  end subroutine fill_key_compute_pattern2
  
!================================================
!       B array for 2D polynomial on HX
!================================================

  subroutine init_blocal_prime_cell_global(mesh)

! io
    type(global_domain), intent(inout)  :: mesh

!================================================
!            vars for local basis
!================================================

    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp
    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny
    real(r8)                            :: xy_list(5)  ! x, y, x2, xy, y2

    integer                             :: p1_index
    integer                             :: iv
    integer                             :: ie
    integer                             :: icell
    integer                             :: iv_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:)
    real(r8), allocatable               :: b_local(:,:)
    integer(i4)                         :: m   ! number of LOB direction
    integer(i4)                         :: n   ! 5 for 2d poly
    integer(i4)                         :: tmp,max_num_2nd,i,j,k,g_index

    if(mpi_rank() == 0) print *,"global init_blocal_prime_cell"
    n   = 5                   ! Xi, Yi, Xi^2, XiYi, Yi^2
    max_num_2nd = mesh%maxvnb !mesh%plg(iv)%stencil_number_2nd = mesh%vtx(iv)%nnb

    DO iv = 1, mesh%nv ! traverse each cell

!================================================
!             common for this cell
!================================================
        
       op0                             = mesh%vtx(iv)%p     ! OP0 vector
       m                               = mesh%vtx(iv)%nnb
       mesh%plg(iv)%stencil_number_2nd = mesh%vtx(iv)%nnb

       if(.not.allocated(p_local))                       allocate(p_local(mesh%plg(iv)%stencil_number_2nd,n))
       if(.not.allocated(b_local))                       allocate(b_local(n,mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%plg(iv)%blocal))           allocate(mesh%plg(iv)%blocal(m,n,mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%plg(iv)%mid_point_lob))    allocate(mesh%plg(iv)%mid_point_lob(m,m,2))
       if(.not.allocated(mesh%plg(iv)%stencil_index_2nd))allocate(mesh%plg(iv)%stencil_index_2nd(1:mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%vtx(iv)%lob_nx))           allocate(mesh%vtx(iv)%lob_nx(m,3))
       if(.not.allocated(mesh%vtx(iv)%lob_ny))           allocate(mesh%vtx(iv)%lob_ny(m,3))

       mesh%plg(iv)%stencil_index_2nd(:) = mesh%vtx(iv)%nb(:)

!================================================
!  Local orthogonal basis using each nbrs icell
!================================================

       do icell = 1, m   ! loop each nb of home cell

          p1_index   = mesh%vtx(iv)%nb(icell)    ! using P1 to form a local basis
          op1(1:3)   = mesh%vtx(p1_index)%p(1:3)   ! OP1 vector
          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)  ! we want is nx, ny

          mesh%vtx(iv)%lob_nx(icell,:) = local_nx
          mesh%vtx(iv)%lob_ny(icell,:) = local_ny

          do iv_local = 1, mesh%plg(iv)%stencil_number_2nd ! loop each cell in the stencil 

             p1_index = mesh%vtx(iv)%nb(iv_local)
             op1(1:3) = mesh%vtx(p1_index)%p(1:3)     ! OP1 vector
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p2nd(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(iv_local,:) = xy_list

!================================================
!  Store each edge's midpoint XY in LOB
!================================================

             p1_index       = mesh%vtx(iv)%ed(iv_local)
             op1(1:3)       = mesh%edp(p1_index)%c%p(1:3)    ! ed1's midpoint vector
             call project_sphere2tplane(op1,op0,p0p1_tp)

             ! the 1st dim stores the direction of EACH LOB
             ! the 2nd dim stores each neighbouring edge of hx(iv)
             ! the 3rd dim stores x and y locations in LOB
             ! e.g., (2,3,1) means the X location of 3rd neighbouring edge of
             ! plg(iv) under the LOB formed by its 2th neighouring point

             mesh%plg(iv)%mid_point_lob(icell,iv_local,1) = dot_product(p0p1_tp,local_nx) ! x
             mesh%plg(iv)%mid_point_lob(icell,iv_local,2) = dot_product(p0p1_tp,local_ny) ! y

          end do

!================================================
!      calculate b_local for current cell iv
!================================================

          call calc_array_b(m,n,p_local,b_local)
          mesh%plg(iv)%blocal(icell,:,:) = b_local

       end do

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing vertice

    DO ie = 1, mesh%ne  ! store for each edge

       v1_index = mesh%edt(ie)%v(1)
       v2_index = mesh%edt(ie)%v(2)

!       allocate(mesh%edt(ie)%v0_weight_prime_cell(mesh%vtx(v1_index)%nnb))
!       allocate(mesh%edt(ie)%v1_weight_prime_cell(mesh%vtx(v2_index)%nnb))
!       allocate(mesh%edt(ie)%v0_weight_1st_prime_cell(mesh%vtx(v1_index)%nnb))
!       allocate(mesh%edt(ie)%v1_weight_1st_prime_cell(mesh%vtx(v2_index)%nnb))

!
! Transver V1
!
       tmp = -999
       do iv_local = 1, mesh%vtx(v1_index)%nnb
          if(mesh%vtx(v1_index)%ed(iv_local) .eq. ie)then
                 mesh%edt(ie)%v_weight_prime_cell(1,1:mesh%vtx(v1_index)%nnb) = mesh%plg(v1_index)%blocal(iv_local,3,:)
             mesh%edt(ie)%v_weight_1st_prime_cell(1,1:mesh%vtx(v1_index)%nnb) = mesh%plg(v1_index)%blocal(iv_local,1,:)
             mesh%edp(ie)%local_index_prime_cell(1)                           = iv_local
             mesh%edp(ie)%v_mid_point_lob_prime_cell(1,:)                     = mesh%plg(v1_index)%mid_point_lob(1,iv_local,:)
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_prime_cell, v1, stop")
!
! Transver V2
!
       tmp = -999
       do iv_local = 1, mesh%vtx(v2_index)%nnb
          if(mesh%vtx(v2_index)%ed(iv_local) .eq. ie)then
                 mesh%edt(ie)%v_weight_prime_cell(2,1:mesh%vtx(v2_index)%nnb) = mesh%plg(v2_index)%blocal(iv_local,3,:)
             mesh%edt(ie)%v_weight_1st_prime_cell(2,1:mesh%vtx(v2_index)%nnb) = mesh%plg(v2_index)%blocal(iv_local,1,:)
             mesh%edp(ie)%local_index_prime_cell(2)                           = iv_local
             mesh%edp(ie)%v_mid_point_lob_prime_cell(2,:)                     = mesh%plg(v2_index)%mid_point_lob(1,iv_local,:)
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_prime_cell, v2, stop")

    END DO

    DO iv = 1, mesh%nv ! Only blocal(1,:,:) needs to be stored
       allocate(b_local(n,mesh%plg(iv)%stencil_number_2nd))
       b_local(:,:) = mesh%plg(iv)%blocal(1,:,:)
       deallocate(mesh%plg(iv)%blocal)
       allocate(mesh%plg(iv)%blocal(1,n,mesh%plg(iv)%stencil_number_2nd))
       mesh%plg(iv)%blocal(1,:,:) = b_local(:,:)
       deallocate(b_local)

       deallocate(mesh%plg(iv)%mid_point_lob)
    END DO

  end subroutine init_blocal_prime_cell_global
  
  subroutine init_blocal_prime_cell(mesh)

! io
    type(global_domain), intent(inout)  :: mesh

!================================================
!            vars for local basis
!================================================

    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp
    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny
    real(r8)                            :: xy_list(5)  ! x, y, x2, xy, y2

    integer                             :: p1_index
    integer                             :: iv
    integer                             :: ie
    integer                             :: icell
    integer                             :: iv_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:)
    real(r8), allocatable               :: b_local(:,:)
    integer(i4)                         :: m   ! number of LOB direction
    integer(i4)                         :: n   ! 5 for 2d poly
    integer(i4)                         :: tmp,max_num_2nd,i,j,k,g_index

    type(scalar_2d_field)               :: tmp_data_exchange_bl,tmp_data_exchange_lob, tmp_data_exchange_2nd, &
                                           tmp_data_exchange_x,tmp_data_exchange_y,&
                                           tmp_data_exchange_1,tmp_data_exchange_2,&
                                           tmp_data_exchange_3,tmp_data_exchange_4

    if(mpi_rank() == 0) print *,"local init_blocal_prime_cell"
    n   = 5                   ! Xi, Yi, Xi^2, XiYi, Yi^2
    max_num_2nd = mesh%maxvnb !mesh%plg(iv)%stencil_number_2nd = mesh%vtx(iv)%nnb

    allocate(tmp_data_exchange_bl%f(mesh%maxvnb*n*max_num_2nd,mesh%nv_full), source=0._r8)
    allocate(tmp_data_exchange_lob%f(mesh%maxvnb*mesh%maxvnb*2,mesh%nv_full),source=0._r8)
    allocate(tmp_data_exchange_2nd%f(max_num_2nd,mesh%nv_full), source=-1._r8)
    allocate(tmp_data_exchange_x%f(mesh%maxvnb*3,mesh%nv_full), source=0._r8)
    allocate(tmp_data_exchange_y%f(mesh%maxvnb*3,mesh%nv_full), source=0._r8)
 
#ifdef LAM_DOMAIN
    mesh%nv = mesh%nv_halo(1)
#else
    mesh%nv = mesh%nv_compute
#endif
    DO iv = 1, mesh%nv ! traverse each cell

!================================================
!             common for this cell
!================================================

       op0                             = mesh%vtx(iv)%p     ! OP0 vector
       m                               = mesh%vtx(iv)%nnb
       mesh%plg(iv)%stencil_number_2nd = mesh%vtx(iv)%nnb

       if(.not.allocated(p_local))                       allocate(p_local(mesh%plg(iv)%stencil_number_2nd,n))
       if(.not.allocated(b_local))                       allocate(b_local(n,mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%plg(iv)%blocal))           allocate(mesh%plg(iv)%blocal(m,n,mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%plg(iv)%mid_point_lob))    allocate(mesh%plg(iv)%mid_point_lob(m,m,2))
       if(.not.allocated(mesh%plg(iv)%stencil_index_2nd))allocate(mesh%plg(iv)%stencil_index_2nd(1:mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%vtx(iv)%lob_nx))           allocate(mesh%vtx(iv)%lob_nx(m,3))
       if(.not.allocated(mesh%vtx(iv)%lob_ny))           allocate(mesh%vtx(iv)%lob_ny(m,3))

#ifdef LAM_DOMAIN
#ifdef REG_QJLAM
       if(mesh%plg_mark(iv).eq.177)then
          mesh%plg(iv)%blocal(:,:,:) = zero
          mesh%plg(iv)%stencil_index_2nd(:) = mesh%vtx(iv)%nb(:)
          deallocate(p_local)
          deallocate(b_local)
          cycle
       end if
#endif
#endif

       mesh%plg(iv)%stencil_index_2nd(:) = mesh%vtx(iv)%nb(:)

!================================================
!  Local orthogonal basis using each nbrs icell
!================================================

       do icell = 1, m   ! loop each nb of home cell

          p1_index   = mesh%vtx(iv)%nb(icell)    ! using P1 to form a local basis
          op1(1:3)   = mesh%vtx(p1_index)%p(1:3)   ! OP1 vector
          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)  ! we want is nx, ny

          mesh%vtx(iv)%lob_nx(icell,:) = local_nx
          mesh%vtx(iv)%lob_ny(icell,:) = local_ny

          do iv_local = 1, mesh%plg(iv)%stencil_number_2nd ! loop each cell in the stencil 

             p1_index = mesh%vtx(iv)%nb(iv_local)
             op1(1:3) = mesh%vtx(p1_index)%p(1:3)     ! OP1 vector
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p2nd(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(iv_local,:) = xy_list

!================================================
!  Store each edge's midpoint XY in LOB
!================================================

             p1_index       = mesh%vtx(iv)%ed(iv_local)
             op1(1:3)       = mesh%edp(p1_index)%c%p(1:3)    ! ed1's midpoint vector
             call project_sphere2tplane(op1,op0,p0p1_tp)

             ! the 1st dim stores the direction of EACH LOB
             ! the 2nd dim stores each neighbouring edge of hx(iv)
             ! the 3rd dim stores x and y locations in LOB
             ! e.g., (2,3,1) means the X location of 3rd neighbouring edge of
             ! plg(iv) under the LOB formed by its 2th neighouring point

             mesh%plg(iv)%mid_point_lob(icell,iv_local,1) = dot_product(p0p1_tp,local_nx) ! x
             mesh%plg(iv)%mid_point_lob(icell,iv_local,2) = dot_product(p0p1_tp,local_ny) ! y
          end do

!================================================
!      calculate b_local for current cell iv
!================================================

          call calc_array_b(m,n,p_local,b_local)
          mesh%plg(iv)%blocal(icell,:,:) = b_local
       end do

       deallocate(p_local)
       deallocate(b_local)

       ! for halo exchange
       do k=1,mesh%plg(iv)%stencil_number_2nd
         do j=1,n
            do i=1,mesh%vtx(iv)%nnb
              tmp_data_exchange_bl%f(i+(j-1)*mesh%maxvnb+(k-1)*mesh%maxvnb*n,iv) = mesh%plg(iv)%blocal(i,j,k)
            end do
         end do
         g_index = mesh%v_index(mesh%plg(iv)%stencil_index_2nd(k))
         tmp_data_exchange_2nd%f(k,iv) = dble(g_index)
       end do
       do k=1,3
         do j=1,mesh%vtx(iv)%nnb
            do i=1,mesh%vtx(iv)%nnb
              if(k <= 2)tmp_data_exchange_lob%f(i+(j-1)*mesh%maxvnb+(k-1)*mesh%maxvnb*mesh%maxvnb,iv) = &
                mesh%plg(iv)%mid_point_lob(i,j,k)
            end do
            tmp_data_exchange_x%f(j+(k-1)*mesh%maxvnb,iv) = mesh%vtx(iv)%lob_nx(j,k)
            tmp_data_exchange_y%f(j+(k-1)*mesh%maxvnb,iv) = mesh%vtx(iv)%lob_ny(j,k)
         end do
       end do
    END DO ! end of traversing vertice
    mesh%nv = mesh%nv_full

    !exchange halo data
    do iv=mesh%nv_compute+1,mesh%nv_full
       m                               = mesh%vtx(iv)%nnb
       mesh%plg(iv)%stencil_number_2nd = mesh%vtx(iv)%nnb

       if(.not.allocated(mesh%plg(iv)%blocal))           allocate(mesh%plg(iv)%blocal(m,n,mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%plg(iv)%mid_point_lob))    allocate(mesh%plg(iv)%mid_point_lob(m,m,2))
       if(.not.allocated(mesh%plg(iv)%stencil_index_2nd))allocate(mesh%plg(iv)%stencil_index_2nd(1:mesh%plg(iv)%stencil_number_2nd))
       if(.not.allocated(mesh%vtx(iv)%lob_nx))           allocate(mesh%vtx(iv)%lob_nx(m,3))
       if(.not.allocated(mesh%vtx(iv)%lob_ny))           allocate(mesh%vtx(iv)%lob_ny(m,3))
    end do
 
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_bl)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_2nd)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_lob)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_x)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_y)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do iv=mesh%nv_compute+1,mesh%nv_full
#ifdef LAM_DOMAIN
      if (mesh%plg_mark(iv) .eq. 177) cycle
#endif
      do k=1,mesh%plg(iv)%stencil_number_2nd
        do j=1,n
          do i=1,mesh%vtx(iv)%nnb
            mesh%plg(iv)%blocal(i,j,k) = tmp_data_exchange_bl%f(i+(j-1)*mesh%maxvnb+(k-1)*mesh%maxvnb*n,iv)
          end do
        end do
        g_index = int(tmp_data_exchange_2nd%f(k,iv))
        mesh%plg(iv)%stencil_index_2nd(k) = mesh%map_g2l_v%find(g_index)
      end do
      do k=1,3
        do j=1,mesh%vtx(iv)%nnb
          do i=1,mesh%vtx(iv)%nnb
            if(k <= 2)mesh%plg(iv)%mid_point_lob(i,j,k) = &
              tmp_data_exchange_lob%f(i+(j-1)*mesh%maxvnb+(k-1)*mesh%maxvnb*mesh%maxvnb,iv)
          end do
          mesh%vtx(iv)%lob_nx(j,k) = tmp_data_exchange_x%f(j+(k-1)*mesh%maxvnb,iv)
          mesh%vtx(iv)%lob_ny(j,k) = tmp_data_exchange_y%f(j+(k-1)*mesh%maxvnb,iv)
        end do
      end do
    end do

    deallocate(tmp_data_exchange_bl%f)
    deallocate(tmp_data_exchange_lob%f)
    deallocate(tmp_data_exchange_2nd%f)
    deallocate(tmp_data_exchange_x%f)
    deallocate(tmp_data_exchange_y%f)
 
    allocate(tmp_data_exchange_1%f(mesh%maxvnb*2,mesh%ne_full))!v_weight_prime_cell
    allocate(tmp_data_exchange_2%f(mesh%maxvnb*2,mesh%ne_full))!v_weight_1st_prime_cell
    allocate(tmp_data_exchange_3%f(2,mesh%ne_full))!local_index_prime_cell
    allocate(tmp_data_exchange_4%f(2*2,mesh%ne_full))!v_mid_point_lob_prime_cell

#ifndef LAM_DOMAIN
    mesh%ne = mesh%ne_compute
#endif
    if(mpi_rank().eq.0) print*,"mesh%ne:",mesh%ne,"mesh%ne_compute:",mesh%ne_compute,"mesh%ne_halo1:",mesh%ne_halo(1),"mesh%ne_full:",mesh%ne_full
    
    DO ie = 1, mesh%ne  ! store for each edge

       v1_index = mesh%edt(ie)%v(1)
       v2_index = mesh%edt(ie)%v(2)

!       allocate(mesh%edt(ie)%v0_weight_prime_cell(mesh%vtx(v1_index)%nnb))
!       allocate(mesh%edt(ie)%v1_weight_prime_cell(mesh%vtx(v2_index)%nnb))
!       allocate(mesh%edt(ie)%v0_weight_1st_prime_cell(mesh%vtx(v1_index)%nnb))
!       allocate(mesh%edt(ie)%v1_weight_1st_prime_cell(mesh%vtx(v2_index)%nnb))

!
! Transver V1
!
#ifdef LAM_DOMAIN
       if (v1_index .le. 0) then
          mesh%edt(ie)%v_weight_prime_cell(1,:)        = zero
          mesh%edt(ie)%v_weight_1st_prime_cell(1,:)    = zero
          mesh%edp(ie)%local_index_prime_cell(1)       = 1
          mesh%edp(ie)%v_mid_point_lob_prime_cell(1,:) = zero
       else
#endif
       tmp = -999
       do iv_local = 1, mesh%vtx(v1_index)%nnb
          if(mesh%vtx(v1_index)%ed(iv_local) .eq. ie)then
                 mesh%edt(ie)%v_weight_prime_cell(1,1:mesh%vtx(v1_index)%nnb) = mesh%plg(v1_index)%blocal(iv_local,3,:)
                 mesh%edt(ie)%v_weight_prime_cell(1,mesh%vtx(v1_index)%nnb+1:mesh%maxvnb) = 0._r8
             mesh%edt(ie)%v_weight_1st_prime_cell(1,1:mesh%vtx(v1_index)%nnb) = mesh%plg(v1_index)%blocal(iv_local,1,:)
             mesh%edt(ie)%v_weight_1st_prime_cell(1,mesh%vtx(v1_index)%nnb+1:mesh%maxvnb) = 0._r8
             mesh%edp(ie)%local_index_prime_cell(1)                           = iv_local
             mesh%edp(ie)%v_mid_point_lob_prime_cell(1,:)                     = mesh%plg(v1_index)%mid_point_lob(1,iv_local,:)
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_prime_cell, v1, stop")
#ifdef LAM_DOMAIN
       end if
#endif
!
! Transver V2
!
#ifdef LAM_DOMAIN
       if(v2_index .le. 0) then
          mesh%edt(ie)%v_weight_prime_cell(2,:)        = zero
          mesh%edt(ie)%v_weight_1st_prime_cell(2,:)    = zero
          mesh%edp(ie)%local_index_prime_cell(2)       = 1
          mesh%edp(ie)%v_mid_point_lob_prime_cell(2,:) = zero
       else
#endif
       tmp = -999
       do iv_local = 1, mesh%vtx(v2_index)%nnb
          if(mesh%vtx(v2_index)%ed(iv_local) .eq. ie)then
                 mesh%edt(ie)%v_weight_prime_cell(2,1:mesh%vtx(v2_index)%nnb) = mesh%plg(v2_index)%blocal(iv_local,3,:)
                 mesh%edt(ie)%v_weight_prime_cell(2,mesh%vtx(v2_index)%nnb+1:mesh%maxvnb) = 0._r8
             mesh%edt(ie)%v_weight_1st_prime_cell(2,1:mesh%vtx(v2_index)%nnb) = mesh%plg(v2_index)%blocal(iv_local,1,:)
             mesh%edt(ie)%v_weight_1st_prime_cell(2,mesh%vtx(v2_index)%nnb+1:mesh%maxvnb) = 0._r8
             mesh%edp(ie)%local_index_prime_cell(2)                           = iv_local
             mesh%edp(ie)%v_mid_point_lob_prime_cell(2,:)                     = mesh%plg(v2_index)%mid_point_lob(1,iv_local,:)
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_prime_cell, v2, stop")
#ifdef LAM_DOMAIN
       end if
#endif

       do k=1,mesh%maxvnb
         do j=1,2
           tmp_data_exchange_1%f(j+(k-1)*2,ie) = mesh%edt(ie)%v_weight_prime_cell(j,k)
           tmp_data_exchange_2%f(j+(k-1)*2,ie) = mesh%edt(ie)%v_weight_1st_prime_cell(j,k)
         end do
       end do
       do k=1,2
         do j=1,2
           tmp_data_exchange_4%f(j+(k-1)*2,ie) = mesh%edp(ie)%v_mid_point_lob_prime_cell(j,k)
         end do
         tmp_data_exchange_3%f(k,ie) = dble(mesh%edp(ie)%local_index_prime_cell(k))
       end do
    END DO
#ifndef LAM_DOMAIN
    mesh%ne = mesh%ne_full
#endif

    !exchange halo data
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_1)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_2)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_3)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_4)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do ie=mesh%ne_compute+1,mesh%ne_full
      do k=1,mesh%maxvnb
        do j=1,2
          mesh%edt(ie)%v_weight_prime_cell(j,k) = tmp_data_exchange_1%f(j+(k-1)*2,ie)
          mesh%edt(ie)%v_weight_1st_prime_cell(j,k) = tmp_data_exchange_2%f(j+(k-1)*2,ie)
        end do
      end do
      do k=1,2
        do j=1,2
          mesh%edp(ie)%v_mid_point_lob_prime_cell(j,k) = tmp_data_exchange_4%f(j+(k-1)*2,ie)
        end do
        mesh%edp(ie)%local_index_prime_cell(k) = int(tmp_data_exchange_3%f(k,ie))
      end do
    end do
 
    deallocate(tmp_data_exchange_1%f)
    deallocate(tmp_data_exchange_2%f)
    deallocate(tmp_data_exchange_3%f)
    deallocate(tmp_data_exchange_4%f)

    allocate(b_local(n,max_num_2nd))
    DO iv = 1, mesh%nv ! Only blocal(1,:,:) needs to be stored
       b_local(:,1:mesh%plg(iv)%stencil_number_2nd) = mesh%plg(iv)%blocal(1,:,:)
       deallocate(mesh%plg(iv)%blocal)
       allocate(mesh%plg(iv)%blocal(1,n,mesh%plg(iv)%stencil_number_2nd))
       mesh%plg(iv)%blocal(1,:,:) = b_local(:,1:mesh%plg(iv)%stencil_number_2nd)

       deallocate(mesh%plg(iv)%mid_point_lob)
    END DO
    deallocate(b_local)
 
    allocate(mesh%edt_v_weight_prime_cell(mesh%maxvnb,2,mesh%ne_full))!v_weight_prime_cell
    allocate(mesh%plg_stencil_number_2nd(mesh%nv_full))
    allocate(mesh%plg_stencil_index_2nd(mesh%maxvnb,mesh%nv_full))
    do iv = 1, mesh%nv_full
       mesh%plg_stencil_number_2nd(iv) = mesh%plg(iv)%stencil_number_2nd
       mesh%plg_stencil_index_2nd(1:mesh%plg(iv)%stencil_number_2nd, iv) = &
           mesh%plg(iv)%stencil_index_2nd(1:mesh%plg(iv)%stencil_number_2nd)
    end do
    do ie = 1, mesh%ne_full
      do k = 1, 2
        mesh%edt_v_weight_prime_cell(1:mesh%maxvnb, k, ie) = &
          mesh%edt(ie)%v_weight_prime_cell(k,1:mesh%maxvnb)
      end do
    end do
  end subroutine init_blocal_prime_cell

!================================================
! Create B array for use in 4th polynomial of HX
!================================================

  subroutine init_blocal_prime_cell_4th_global(mesh)

! io
    type(global_domain), intent(inout)  :: mesh

!================================================
!            Vars for local basis
!================================================

    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp
    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny
    real(r8), dimension(1:14)           :: xy_list

    integer(i4)                         :: p1_index
    integer(i4)                         :: iv        ! global index
    integer(i4)                         :: ie        ! global index
    integer(i4)                         :: icell
    integer(i4)                         :: iv_local
    integer(i4)                         :: v1_index
    integer(i4)                         :: v2_index

!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions of LOB
    integer(i4)                         :: n            ! 14 for 4d
    integer(i4)                         :: tmp
  
    n   =  14

    if(mpi_rank() == 0) print*,"global init_blocal_prime_cell_4th"
    DO iv = 1, mesh%nv

! create local stencil surrounding this cell
       call create_index_list_hx4th(mesh,iv)

!================================================
!            Common for this cell iv
!================================================

       op0     = mesh%vtx(iv)%p        ! OP0 vector
       m       = mesh%vtx(iv)%nnb      ! how many direction of LOB axis

       if(.not.allocated(p_local))                 allocate(p_local(mesh%plg(iv)%stencil_number_4th,n))
       if(.not.allocated(b_local))                 allocate(b_local(n,mesh%plg(iv)%stencil_number_4th))
       if(.not.allocated(mesh%plg(iv)%blocal_4th)) allocate(mesh%plg(iv)%blocal_4th(m,n,mesh%plg(iv)%stencil_number_4th))

!================================================
!      Local Orthogonal Basis using each nbr
!================================================

       do icell = 1, m

          p1_index   = mesh%vtx(iv)%nb(icell)    ! using P1 to form a local basis
          op1(1:3)   = mesh%vtx(p1_index)%p(1:3)   ! OP1 vector
          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          do iv_local = 1, mesh%plg(iv)%stencil_number_4th    ! for cell in the stencil

             p1_index   = mesh%plg(iv)%stencil_index_4th(iv_local)
             op1(1:3)   = mesh%vtx(p1_index)%p(1:3)             ! OP1 vector
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p4th(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(iv_local,1:14) = xy_list

          end do

!================================================
!      calculate b_local for current cell iv
!================================================

          call calc_array_b(mesh%plg(iv)%stencil_number_4th,n,p_local,b_local)
          mesh%plg(iv)%blocal_4th(icell,1:n,1:mesh%plg(iv)%stencil_number_4th) = b_local(1:n,1:mesh%plg(iv)%stencil_number_4th)
       end do

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing vertice
 
! store for each edge
    DO ie = 1, mesh%ne

       v1_index = mesh%edt(ie)%v(1)
       v2_index = mesh%edt(ie)%v(2)

!       allocate(mesh%edt(ie)%v0_weight_4th_prime_cell(mesh%plg(v1_index)%stencil_number_4th))
!       allocate(mesh%edt(ie)%v1_weight_4th_prime_cell(mesh%plg(v2_index)%stencil_number_4th))
!       allocate(mesh%edt(ie)%v0_weight_4th_2nd_prime_cell(mesh%plg(v1_index)%stencil_number_4th))
!       allocate(mesh%edt(ie)%v1_weight_4th_2nd_prime_cell(mesh%plg(v2_index)%stencil_number_4th))

       ! transver v1
       tmp = -999
       do iv_local = 1, mesh%vtx(v1_index)%nnb
          if(mesh%vtx(v1_index)%ed(iv_local) .eq. ie)then
                 mesh%edt(ie)%v_weight_4th_prime_cell(1,1:mesh%plg(v1_index)%stencil_number_4th) &
                         = mesh%plg(v1_index)%blocal_4th(iv_local,10,:)
             mesh%edt(ie)%v_weight_4th_2nd_prime_cell(1,1:mesh%plg(v1_index)%stencil_number_4th) &
                     = mesh%plg(v1_index)%blocal_4th(iv_local,3,:)
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_4th_prime_cell, v1, stop")

       ! Transver v2
       tmp = -999
       do iv_local = 1, mesh%vtx(v2_index)%nnb
          if(mesh%vtx(v2_index)%ed(iv_local) .eq. ie)then
                 mesh%edt(ie)%v_weight_4th_prime_cell(2,1:mesh%plg(v2_index)%stencil_number_4th) &
                         = mesh%plg(v2_index)%blocal_4th(iv_local,10,:)
             mesh%edt(ie)%v_weight_4th_2nd_prime_cell(2,1:mesh%plg(v2_index)%stencil_number_4th) &
                     = mesh%plg(v2_index)%blocal_4th(iv_local,3,:)
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_4th_prime_cell, v2, stop")

    END DO

    DO iv = 1, mesh%nv
       deallocate(mesh%plg(iv)%blocal_4th)
    END DO

  end subroutine init_blocal_prime_cell_4th_global
  
  subroutine init_blocal_prime_cell_4th(mesh)

! io
    type(global_domain), intent(inout)  :: mesh

!================================================
!            Vars for local basis
!================================================

    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp
    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny
    real(r8), dimension(1:14)           :: xy_list

    integer(i4)                         :: p1_index
    integer(i4)                         :: iv        ! global index
    integer(i4)                         :: ie        ! global index
    integer(i4)                         :: icell
    integer(i4)                         :: iv_local
    integer(i4)                         :: v1_index
    integer(i4)                         :: v2_index

!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions of LOB
    integer(i4)                         :: n            ! 14 for 4d
    integer(i4)                         :: tmp
    
    type(scalar_2d_field)               :: tmp_data_exchange,tmp_data_exchange_bl
    integer(i4)                         :: max_num_4d,i,j,k,g_index, max_num_4d_global, err

    n   =  14
    !max_num_4d = 18

    if(mpi_rank() == 0) print*,"local init_blocal_prime_cell_4th"
    mesh%nv = mesh%nv_compute
    DO iv = 1, mesh%nv

! create local stencil surrounding this cell
       call create_index_list_hx4th(mesh,iv)

!================================================
!            Common for this cell iv
!================================================

       op0     = mesh%vtx(iv)%p        ! OP0 vector
       m       = mesh%vtx(iv)%nnb      ! how many direction of LOB axis

       if(.not.allocated(p_local))                 allocate(p_local(mesh%plg(iv)%stencil_number_4th,n))
       if(.not.allocated(b_local))                 allocate(b_local(n,mesh%plg(iv)%stencil_number_4th))
       if(.not.allocated(mesh%plg(iv)%blocal_4th)) allocate(mesh%plg(iv)%blocal_4th(m,n,mesh%plg(iv)%stencil_number_4th))

!================================================
!      Local Orthogonal Basis using each nbr
!================================================

       do icell = 1, m

          p1_index   = mesh%vtx(iv)%nb(icell)    ! using P1 to form a local basis
          op1(1:3)   = mesh%vtx(p1_index)%p(1:3)   ! OP1 vector
          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          do iv_local = 1, mesh%plg(iv)%stencil_number_4th    ! for cell in the stencil

             p1_index   = mesh%plg(iv)%stencil_index_4th(iv_local)
             op1(1:3)   = mesh%vtx(p1_index)%p(1:3)             ! OP1 vector
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p4th(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(iv_local,1:14) = xy_list

          end do

!================================================
!      calculate b_local for current cell iv
!================================================

          call calc_array_b(mesh%plg(iv)%stencil_number_4th,n,p_local,b_local)
          mesh%plg(iv)%blocal_4th(icell,1:n,1:mesh%plg(iv)%stencil_number_4th) = b_local(1:n,1:mesh%plg(iv)%stencil_number_4th)
       end do

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing vertice
    mesh%nv = mesh%nv_full

    !exchange halo data
    allocate(tmp_data_exchange%f(1,mesh%nv_full))
    do iv=1,mesh%nv_compute
       tmp_data_exchange%f(1,iv)=dble(mesh%plg(iv)%stencil_number_4th)
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    max_num_4d = 0
    do iv=1,mesh%nv_full
       if(iv > mesh%nv_compute) mesh%plg(iv)%stencil_number_4th = int(tmp_data_exchange%f(1,iv))
       if(mesh%plg(iv)%stencil_number_4th > max_num_4d) max_num_4d = mesh%plg(iv)%stencil_number_4th
    end do
    deallocate(tmp_data_exchange%f)
    call MPI_ALLREDUCE(max_num_4d, max_num_4d_global, 1, MPI_INT, MPI_MAX, grist_comm(), err)

    do iv=mesh%nv_compute+1,mesh%nv_full
       m = mesh%vtx(iv)%nnb
       if(.not.allocated(mesh%plg(iv)%stencil_index_4th)) &
          allocate(mesh%plg(iv)%stencil_index_4th(mesh%plg(iv)%stencil_number_4th))
       if(.not.allocated(mesh%plg(iv)%blocal_4th)) allocate(mesh%plg(iv)%blocal_4th(m,n,mesh%plg(iv)%stencil_number_4th))
    end do
 
    allocate(tmp_data_exchange%f(max_num_4d_global,mesh%nv_full), source=0._r8)!stencil_index_4th
    allocate(tmp_data_exchange_bl%f(mesh%maxvnb*n*max_num_4d_global,mesh%nv_full), source=0._r8)!blocal_4th
    do iv=1,mesh%nv_compute
      m = mesh%vtx(iv)%nnb
      do k=1,mesh%plg(iv)%stencil_number_4th
        do j=1,n
          do i=1,m
            tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,iv) = mesh%plg(iv)%blocal_4th(i,j,k)
          end do
        end do
        g_index = mesh%v_index(mesh%plg(iv)%stencil_index_4th(k))
        tmp_data_exchange%f(k,iv)=dble(g_index)
      end do
    end do
 
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_bl)
    call exchange_data_2d(mesh%local_block,field_head_2d)
 
    do iv=mesh%nv_compute+1,mesh%nv_full!iv=nv_compute+1 better
      m = mesh%vtx(iv)%nnb
      do k=1,mesh%plg(iv)%stencil_number_4th
        do j=1,n
          do i=1,m
            mesh%plg(iv)%blocal_4th(i,j,k) = tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,iv)
          end do
        end do
        g_index = int(tmp_data_exchange%f(k,iv))
        mesh%plg(iv)%stencil_index_4th(k) = mesh%map_g2l_v%find(g_index)
      end do
    end do

    deallocate(tmp_data_exchange%f)!stencil_index_4th
    deallocate(tmp_data_exchange_bl%f)!blocal_4th
! store for each edge

    mesh%ne = mesh%ne_compute
    DO ie = 1, mesh%ne

       v1_index = mesh%edt(ie)%v(1)
       v2_index = mesh%edt(ie)%v(2)

!       allocate(mesh%edt(ie)%v0_weight_4th_prime_cell(mesh%plg(v1_index)%stencil_number_4th))
!       allocate(mesh%edt(ie)%v1_weight_4th_prime_cell(mesh%plg(v2_index)%stencil_number_4th))
!       allocate(mesh%edt(ie)%v0_weight_4th_2nd_prime_cell(mesh%plg(v1_index)%stencil_number_4th))
!       allocate(mesh%edt(ie)%v1_weight_4th_2nd_prime_cell(mesh%plg(v2_index)%stencil_number_4th))

       ! transver v1
       tmp = -999
       do iv_local = 1, mesh%vtx(v1_index)%nnb
          if(mesh%vtx(v1_index)%ed(iv_local) .eq. ie)then
             mesh%edt(ie)%v_weight_4th_prime_cell(1,1:mesh%plg(v1_index)%stencil_number_4th) &
                     = mesh%plg(v1_index)%blocal_4th(iv_local,10,:)
             mesh%edt(ie)%v_weight_4th_prime_cell(1,mesh%plg(v1_index)%stencil_number_4th+1:23) = 0._r8
             mesh%edt(ie)%v_weight_4th_2nd_prime_cell(1,1:mesh%plg(v1_index)%stencil_number_4th) &
                     = mesh%plg(v1_index)%blocal_4th(iv_local,3,:)
             mesh%edt(ie)%v_weight_4th_2nd_prime_cell(1,mesh%plg(v1_index)%stencil_number_4th+1:23) = 0._r8
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_4th_prime_cell, v1, stop")

       ! Transver v2
       tmp = -999
       do iv_local = 1, mesh%vtx(v2_index)%nnb
          if(mesh%vtx(v2_index)%ed(iv_local) .eq. ie)then
             mesh%edt(ie)%v_weight_4th_prime_cell(2,1:mesh%plg(v2_index)%stencil_number_4th) &
                     = mesh%plg(v2_index)%blocal_4th(iv_local,10,:)
             mesh%edt(ie)%v_weight_4th_prime_cell(2,mesh%plg(v2_index)%stencil_number_4th+1:23) = 0._r8
             mesh%edt(ie)%v_weight_4th_2nd_prime_cell(2,1:mesh%plg(v2_index)%stencil_number_4th) &
                     = mesh%plg(v2_index)%blocal_4th(iv_local,3,:)
             mesh%edt(ie)%v_weight_4th_2nd_prime_cell(2,mesh%plg(v2_index)%stencil_number_4th+1:23) = 0._r8
             tmp = 0
             exit
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_4th_prime_cell, v2, stop")

    END DO
    mesh%ne = mesh%ne_full

    !exchange halo data
    allocate(tmp_data_exchange%f(2*23,mesh%ne_full))!v_weight_4th_prime_cell
    allocate(tmp_data_exchange_bl%f(2*23,mesh%ne_full))!v_weight_4th_2nd_prime_cell
    do ie=1,mesh%ne_compute
      do j=1,23
        do i=1,2
          tmp_data_exchange%f(i+(j-1)*2,ie) = mesh%edt(ie)%v_weight_4th_prime_cell(i,j)
          tmp_data_exchange_bl%f(i+(j-1)*2,ie) = mesh%edt(ie)%v_weight_4th_2nd_prime_cell(i,j)
        end do
      end do
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_bl)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do ie=mesh%ne_compute+1,mesh%ne_full
      do j=1,23
        do i=1,2
          mesh%edt(ie)%v_weight_4th_prime_cell(i,j) = tmp_data_exchange%f(i+(j-1)*2,ie)
          mesh%edt(ie)%v_weight_4th_2nd_prime_cell(i,j) = tmp_data_exchange_bl%f(i+(j-1)*2,ie)
        end do
      end do
    end do
    deallocate(tmp_data_exchange%f)
    deallocate(tmp_data_exchange_bl%f)
 
    DO iv = 1, mesh%nv
      deallocate(mesh%plg(iv)%blocal_4th)
    END DO

    allocate(mesh%edt_v_weight_4th_prime_cell(23,2,mesh%ne_full))!v_weight_prime_cell
    allocate(mesh%edt_v_weight_4th_2nd_prime_cell(23,2,mesh%ne_full))!v_weight_prime_cell
    allocate(mesh%plg_stencil_number_4th(mesh%nv_full))
    allocate(mesh%plg_stencil_index_4th(23,mesh%nv_full))
    do iv = 1, mesh%nv_full
      mesh%plg_stencil_number_4th(iv) = mesh%plg(iv)%stencil_number_4th
      mesh%plg_stencil_index_4th(1:mesh%plg(iv)%stencil_number_4th, iv) = &
         mesh%plg(iv)%stencil_index_4th(1:mesh%plg(iv)%stencil_number_4th)
    end do
    do ie = 1, mesh%ne_full
      do k = 1, 2
        mesh%edt_v_weight_4th_prime_cell(1:23, k, ie) = &
          mesh%edt(ie)%v_weight_4th_prime_cell(k,1:23)
        mesh%edt_v_weight_4th_2nd_prime_cell(1:23, k, ie) = &
          mesh%edt(ie)%v_weight_4th_2nd_prime_cell(k,1:23)
      end do
    end do
  end subroutine init_blocal_prime_cell_4th

!================================================
!     B array for 2D polynomial on TR
!================================================

  subroutine init_blocal_dual_cell_global(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp
    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny
    real(r8), dimension(1:5)            :: xy_list   ! x, y, x2, xy, y2

    integer                             :: p1_index
    integer                             :: it        ! global index 
    integer                             :: ie        ! global index 
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin
 
    n = 5   ! Xi, Yi, Xi^2, XiYi, Yi^2

    if(mpi_rank() == 0) print *,"global init_blocal_dual_cell"
    DO it = 1, mesh%nt   ! traverse each dual cell

!================================================
!            common for this tr cell
!================================================
        call create_index_list_tr(mesh,it)

        num_point_of_stencil = mesh%tri(it)%stencil_number_2nd

        op0     = mesh%tri(it)%c%p      ! OP0 vector (home cell)
        m       = 3                    ! 3 directions for LOB

        if(.not.allocated(mesh%tri(it)%blocal)) allocate(mesh%tri(it)%blocal(m,n,num_point_of_stencil))
        if(.not.allocated(b_local))             allocate(b_local(n,num_point_of_stencil))
        if(.not.allocated(p_local))             allocate(p_local(num_point_of_stencil,n))
        if(.not.allocated(mesh%tri(it)%lob_nx)) allocate(mesh%tri(it)%lob_nx(m,3))
        if(.not.allocated(mesh%tri(it)%lob_ny)) allocate(mesh%tri(it)%lob_ny(m,3))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)    ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          mesh%tri(it)%lob_nx(icell,:) = local_nx
          mesh%tri(it)%lob_ny(icell,:) = local_ny

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_2nd(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p2nd(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

    !================================================
    !      calculate B using P for current TR
    !================================================
          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal(icell,:,:) = b_local

          ! calculate each neighbouring edge's midpoint location in LOB
          do it_local = 1,3 

              p1_index       = mesh%tri(it)%ed(it_local)
              op1(1:3)       = mesh%edt(p1_index)%c%p(1:3)      ! ed1's midpoint vector
              call project_sphere2tplane(op1,op0,p0p1_tp)

            ! the 1st dim stores the direction of EACH LOB
            ! the 2nd dim stores each neighbouring edge of tr(it)
            ! the 3rd dim stores x and y locations in LOB
            ! e.g., (2,3,1) means the x location of 3rd neighbouring edge of
            ! tr(it) under the LOB formed by its 2th neighouring point
              mesh%tri(it)%mid_point_lob(icell,it_local,1) = dot_product(p0p1_tp,local_nx) !x
              mesh%tri(it)%mid_point_lob(icell,it_local,2) = dot_product(p0p1_tp,local_ny) !y
          end do

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing triangle 
 
    DO ie = 1, mesh%ne  ! store for each hx'S edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_dual_cell(mesh%tri(v1_index)%stencil_number_2nd))
!       allocate(mesh%edp(ie)%v1_weight_dual_cell(mesh%tri(v2_index)%stencil_number_2nd))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_2nd)  &
                  = mesh%tri(v1_index)%blocal(it_local,3,:)
             mesh%edt(ie)%local_index_dual_cell(1)       = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(1,:) = mesh%tri(v1_index)%mid_point_lob(1,it_local,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v1, stop")

       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_2nd)  &
                   = mesh%tri(v2_index)%blocal(it_local,3,:)
             mesh%edt(ie)%local_index_dual_cell(2)        = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(2,:) = mesh%tri(v2_index)%mid_point_lob(1,it_local,:) 
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v2, stop")

    END DO

    DO it = 1, mesh%nt ! Only blocal(1,:,:) needs to be stored
       allocate(b_local(n,mesh%tri(it)%stencil_number_2nd))
       b_local(:,:) = mesh%tri(it)%blocal(1,:,:)
       deallocate(mesh%tri(it)%blocal)
       allocate(mesh%tri(it)%blocal(1,n,mesh%tri(it)%stencil_number_2nd))
       mesh%tri(it)%blocal(1,:,:) = b_local(:,:)
       deallocate(b_local)
    END DO
  end subroutine init_blocal_dual_cell_global
  
  subroutine init_blocal_dual_cell(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp
    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny
    real(r8), dimension(1:5)            :: xy_list   ! x, y, x2, xy, y2

    integer                             :: p1_index
    integer                             :: it        ! global index 
    integer                             :: ie        ! global index 
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin
    
    type(scalar_2d_field)               :: tmp_data_exchange,tmp_data_exchange_bl,&
                                           tmp_data_exchange_x,tmp_data_exchange_y,&
                                           tmp_data_exchange1
    integer(i4)                         :: max_num_2nd,i,j,k,g_index  !max_num_2nd = 16?
 
    n = 5   ! Xi, Yi, Xi^2, XiYi, Yi^2
    m = 3   ! 3 directions for LOB
    max_num_2nd = 16
    num_point_of_stencil = 9

    if(mpi_rank() == 0) print *,"local init_blocal_dual_cell"

    allocate(tmp_data_exchange%f(1,mesh%nt_full), source=9._r8)
    allocate(tmp_data_exchange1%f(max_num_2nd,mesh%nt_full), source=0._r8)!stencil_index_2nd
    allocate(tmp_data_exchange_bl%f(max_num_2nd*m*n,mesh%nt_full), source=0._r8)!mesh%tri(it)%blocal
    allocate(tmp_data_exchange_x%f(3*m,mesh%nt_full), source=0._r8)!mesh%tri(it)%lob_nx
    allocate(tmp_data_exchange_y%f(3*m,mesh%nt_full), source=0._r8)!mesh%tri(it)%lob_ny
#ifdef LAM_DOMAIN
    mesh%nt = mesh%nt_halo(1)
#else
    mesh%nt = mesh%nt_compute
#endif
    DO it = 1, mesh%nt   ! traverse each dual cell

#ifdef LAM_DOMAIN
        if(mesh%tri_mark(it).eq.177)then
           if(.not.allocated(mesh%tri(it)%blocal))           allocate(mesh%tri(it)%blocal(m,n,num_point_of_stencil))
           mesh%tri(it)%blocal          = zero
           mesh%tri(it)%mid_point_lob   = zero
           mesh%tri(it)%stencil_number_2nd = num_point_of_stencil
           if(.not.allocated(mesh%tri(it)%stencil_index_2nd))allocate(mesh%tri(it)%stencil_index_2nd(num_point_of_stencil))
           cycle
        end if
#endif
!================================================
!            common for this tr cell
!================================================
       call create_index_list_tr(mesh,it)

       num_point_of_stencil = mesh%tri(it)%stencil_number_2nd

       op0     = mesh%tri(it)%c%p     ! OP0 vector (home cell)

       if(.not.allocated(mesh%tri(it)%blocal)) allocate(mesh%tri(it)%blocal(m,n,num_point_of_stencil))
       if(.not.allocated(b_local))             allocate(b_local(n,num_point_of_stencil))
       if(.not.allocated(p_local))             allocate(p_local(num_point_of_stencil,n))
       if(.not.allocated(mesh%tri(it)%lob_nx)) allocate(mesh%tri(it)%lob_nx(m,3))
       if(.not.allocated(mesh%tri(it)%lob_ny)) allocate(mesh%tri(it)%lob_ny(m,3))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)    ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          mesh%tri(it)%lob_nx(icell,:) = local_nx
          mesh%tri(it)%lob_ny(icell,:) = local_ny

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_2nd(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p2nd(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

    !================================================
    !      calculate B using P for current TR
    !================================================
          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal(icell,:,:) = b_local

          ! calculate each neighbouring edge's midpoint location in LOB, for ffsl
          do it_local = 1,3

              p1_index       = mesh%tri(it)%ed(it_local)
              op1(1:3)       = mesh%edt(p1_index)%c%p(1:3)      ! ed1's midpoint vector
              call project_sphere2tplane(op1,op0,p0p1_tp)

            ! the 1st dim stores the direction of EACH LOB
            ! the 2nd dim stores each neighbouring edge of tr(it)
            ! the 3rd dim stores x and y locations in LOB
            ! e.g., (2,3,1) means the x location of 3rd neighbouring edge of
            ! tr(it) under the LOB formed by its 2th neighouring point
              mesh%tri(it)%mid_point_lob(icell,it_local,1) = dot_product(p0p1_tp,local_nx) !x
              mesh%tri(it)%mid_point_lob(icell,it_local,2) = dot_product(p0p1_tp,local_ny) !y
          end do

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)
       ! for halo exchange
       tmp_data_exchange%f(1,it)=mesh%tri(it)%stencil_number_2nd
       do k=1,mesh%tri(it)%stencil_number_2nd
         do j=1,n
           do i=1,m
             tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it) = mesh%tri(it)%blocal(i,j,k)
           end do
         end do
         g_index = mesh%t_index(mesh%tri(it)%stencil_index_2nd(k))
         tmp_data_exchange1%f(k,it)=g_index
       end do
       do j=1,3
         do i=1,m
           tmp_data_exchange_x%f(i+(j-1)*m,it) = mesh%tri(it)%lob_nx(i,j)
           tmp_data_exchange_y%f(i+(j-1)*m,it) = mesh%tri(it)%lob_ny(i,j)
         end do
       end do 
    END DO ! end of traversing triangle 
    mesh%nt = mesh%nt_full

    !exchange halo data

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d(mesh%local_block,field_head_2d)
    do it=mesh%nt_compute+1,mesh%nt_full
       mesh%tri(it)%stencil_number_2nd = int(tmp_data_exchange%f(1,it))
    end do
    deallocate(tmp_data_exchange%f)

    do it = mesh%nt_compute+1,mesh%nt_full
       if(.not.allocated(mesh%tri(it)%blocal)) allocate(mesh%tri(it)%blocal(m,n,mesh%tri(it)%stencil_number_2nd))
       if(.not.allocated(mesh%tri(it)%lob_nx)) allocate(mesh%tri(it)%lob_nx(m,3))
       if(.not.allocated(mesh%tri(it)%lob_ny)) allocate(mesh%tri(it)%lob_ny(m,3))
       if(.not.allocated(mesh%tri(it)%stencil_index_2nd)) &
          allocate(mesh%tri(it)%stencil_index_2nd(mesh%tri(it)%stencil_number_2nd))
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange1)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_bl)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_x)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_y)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do it=mesh%nt_compute+1,mesh%nt_full
#ifdef LAM_DOMAIN
      if (mesh%tri_mark(it) .eq. 177) cycle
#endif
      do k=1,mesh%tri(it)%stencil_number_2nd
        do j=1,n
          do i=1,m
            mesh%tri(it)%blocal(i,j,k) = tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it)
          end do
        end do 
        g_index = int(tmp_data_exchange1%f(k,it))
        mesh%tri(it)%stencil_index_2nd(k) = mesh%map_g2l_t%find(g_index)
      end do
      do j=1,3
        do i=1,m
          mesh%tri(it)%lob_nx(i,j) = tmp_data_exchange_x%f(i+(j-1)*m,it)
          mesh%tri(it)%lob_ny(i,j) = tmp_data_exchange_y%f(i+(j-1)*m,it)
        end do
      end do 
    end do
    deallocate(tmp_data_exchange1%f)
    deallocate(tmp_data_exchange_bl%f)
    deallocate(tmp_data_exchange_x%f)
    deallocate(tmp_data_exchange_y%f)

    allocate(tmp_data_exchange%f(2,mesh%ne_full))!local_index_dual_cell
    allocate(tmp_data_exchange_x%f(max_num_2nd*2,mesh%ne_full))!v_weight_dual_cell
    allocate(tmp_data_exchange_y%f(2*2,mesh%ne_full))!v_mid_point_lob_dual_cell
#ifdef LAM_DOMAIN
    mesh%ne = mesh%ne_halo(1)
#else
    mesh%ne = mesh%ne_compute
#endif
    DO ie = 1, mesh%ne  ! store for each hx'S edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_dual_cell(mesh%tri(v1_index)%stencil_number_2nd))
!       allocate(mesh%edp(ie)%v1_weight_dual_cell(mesh%tri(v2_index)%stencil_number_2nd))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_2nd)  &
                  = mesh%tri(v1_index)%blocal(it_local,3,:)
             mesh%edp(ie)%v_weight_dual_cell(1,mesh%tri(v1_index)%stencil_number_2nd+1:16) = 0._r8
             mesh%edt(ie)%local_index_dual_cell(1)       = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(1,:) = mesh%tri(v1_index)%mid_point_lob(1,it_local,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v1, stop")

       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_2nd)  &
                   = mesh%tri(v2_index)%blocal(it_local,3,:)
             mesh%edp(ie)%v_weight_dual_cell(2,mesh%tri(v2_index)%stencil_number_2nd+1:16) = 0._r8
             mesh%edt(ie)%local_index_dual_cell(2)        = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(2,:) = mesh%tri(v2_index)%mid_point_lob(1,it_local,:) 
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v2, stop")

       do j=1,2
         do i=1,2
           tmp_data_exchange_y%f(i+(j-1)*2,ie) = mesh%edt(ie)%v_mid_point_lob_dual_cell(i,j)
         end do
         tmp_data_exchange%f(j,ie) = dble(mesh%edt(ie)%local_index_dual_cell(j))
       end do
       do j=1, max_num_2nd
         do i=1,2
           tmp_data_exchange_x%f(i+(j-1)*2,ie) = mesh%edp(ie)%v_weight_dual_cell(i,j)
         end do
       end do
    END DO
    mesh%ne = mesh%ne_full

    !exchange halo data
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_x)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_y)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do ie=mesh%ne_compute+1,mesh%ne_full
      do j=1,2
        do i=1,2
          mesh%edt(ie)%v_mid_point_lob_dual_cell(i,j) = tmp_data_exchange_y%f(i+(j-1)*2,ie)
        end do
        mesh%edt(ie)%local_index_dual_cell(j) = int(tmp_data_exchange%f(j,ie))
      end do
      do j=1, max_num_2nd
        do i=1,2
          mesh%edp(ie)%v_weight_dual_cell(i,j) = tmp_data_exchange_x%f(i+(j-1)*2,ie)
        end do
      end do
    end do

    deallocate(tmp_data_exchange%f)
    deallocate(tmp_data_exchange_x%f)
    deallocate(tmp_data_exchange_y%f)

    allocate(b_local(n,max_num_2nd))
    DO it = 1, mesh%nt ! Only blocal(1,:,:) needs to be stored
       b_local(:,1:mesh%tri(it)%stencil_number_2nd) = mesh%tri(it)%blocal(1,:,:)
       deallocate(mesh%tri(it)%blocal)
       allocate(mesh%tri(it)%blocal(1,n,mesh%tri(it)%stencil_number_2nd))
       mesh%tri(it)%blocal(1,:,:) = b_local(:,1:mesh%tri(it)%stencil_number_2nd)
    END DO
    deallocate(b_local)
 
    allocate(mesh%edp_v_weight_dual_cell(16,2,mesh%ne_full))!v_weight_prime_cell
    allocate(mesh%tri_stencil_number_2nd(mesh%nt_full))
    allocate(mesh%tri_stencil_index_2nd(16,mesh%nt_full))
    do it = 1, mesh%nt_full
      mesh%tri_stencil_number_2nd(it) = mesh%tri(it)%stencil_number_2nd
      mesh%tri_stencil_index_2nd(1:mesh%tri(it)%stencil_number_2nd, it) = &
           mesh%tri(it)%stencil_index_2nd(1:mesh%tri(it)%stencil_number_2nd)
    end do
    do ie = 1, mesh%ne_full
      do k = 1, 2
        mesh%edp_v_weight_dual_cell(1:16, k, ie) = &
          mesh%edp(ie)%v_weight_dual_cell(k,1:16)
      end do
    end do
  end subroutine init_blocal_dual_cell

!================================================
!     B array for 4D polynomial on TR
!================================================

  subroutine init_blocal_dual_cell_4th_global(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp

    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny

    real(r8)                            :: xy_list(14)  !
    integer                             :: p1_index
    integer                             :: it        ! global index 
    integer                             :: ie        ! global index 
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin

! fixed 
    n       = 14
    m       = 3   ! 3 directions for LOB

    if(mpi_rank() == 0) print *,"global init_blocal_dual_cell_4th"
    DO it = 1, mesh%nt   ! traverse each dual cell

!================================================
!            common for this tr cell
!================================================
       call create_index_list_tr4th(mesh,it)
       num_point_of_stencil = mesh%tri(it)%stencil_number_4th

       op0     = mesh%tri(it)%c%p      ! OP0 vector (home cell)

       if(.not.allocated(b_local))                 allocate(b_local(n,num_point_of_stencil))
       if(.not.allocated(p_local))                 allocate(p_local(num_point_of_stencil,n))
       if(.not.allocated(mesh%tri(it)%blocal_4th)) allocate(mesh%tri(it)%blocal_4th(m,n,num_point_of_stencil))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)         ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_4th(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p4th(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

!================================================
!      calculate B using P for current TR
!================================================

          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal_4th(icell,:,:) = b_local

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing triangle 
 
    DO ie = 1, mesh%ne  ! store for each hx's edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_4th_dual_cell(mesh%tri(v1_index)%stencil_number_4th))
!       allocate(mesh%edp(ie)%v1_weight_4th_dual_cell(mesh%tri(v2_index)%stencil_number_4th))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_4th) &
                     = mesh%tri(v1_index)%blocal_4th(it_local,10,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v1, stop")

       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_4th) &
                     = mesh%tri(v2_index)%blocal_4th(it_local,10,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v2, stop")

    END DO

    DO it = 1, mesh%nt
      deallocate(mesh%tri(it)%blocal_4th)
    END DO
  end subroutine init_blocal_dual_cell_4th_global
  
  subroutine init_blocal_dual_cell_4th(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp

    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny

    real(r8)                            :: xy_list(14)  !
    integer                             :: p1_index
    integer                             :: it        ! global index 
    integer                             :: ie        ! global index 
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin
    type(scalar_2d_field)               :: tmp_data_exchange,tmp_data_exchange_bl
    integer(i4)                         :: max_num_4d,i,j,k,g_index, max_num_4d_global, err

! fixed 
    n = 14
    m = 3   ! 3 directions for LOB

    if(mpi_rank() == 0) print*,"local init_blocal_dual_cell_4th"
    mesh%nt = mesh%nt_compute
    DO it = 1, mesh%nt   ! traverse each dual cell

!================================================
!            common for this tr cell
!================================================
       call create_index_list_tr4th(mesh,it)
       num_point_of_stencil = mesh%tri(it)%stencil_number_4th

       op0     = mesh%tri(it)%c%p      ! OP0 vector (home cell)

       if(.not.allocated(b_local))                 allocate(b_local(n,num_point_of_stencil))
       if(.not.allocated(p_local))                 allocate(p_local(num_point_of_stencil,n))
       if(.not.allocated(mesh%tri(it)%blocal_4th)) allocate(mesh%tri(it)%blocal_4th(m,n,num_point_of_stencil))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)         ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_4th(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p4th(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

!================================================
!      calculate B using P for current TR
!================================================

          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal_4th(icell,:,:) = b_local

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing triangle 
    mesh%nt = mesh%nt_full

    !exchange halo data
    allocate(tmp_data_exchange%f(1,mesh%nt_full))
    do it=1,mesh%nt_compute
       tmp_data_exchange%f(1,it)=dble(mesh%tri(it)%stencil_number_4th)
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    max_num_4d = 0
    do it=1,mesh%nt_full
       if(it > mesh%nt_compute) mesh%tri(it)%stencil_number_4th = int(tmp_data_exchange%f(1,it))
       if(mesh%tri(it)%stencil_number_4th > max_num_4d) max_num_4d = mesh%tri(it)%stencil_number_4th
    end do
    deallocate(tmp_data_exchange%f)
    call MPI_ALLREDUCE(max_num_4d, max_num_4d_global, 1, MPI_INT, MPI_MAX, grist_comm(), err)

    do it=mesh%nt_compute+1,mesh%nt_full
       if(.not.allocated(mesh%tri(it)%stencil_index_4th)) &
         allocate(mesh%tri(it)%stencil_index_4th(mesh%tri(it)%stencil_number_4th))
       if(.not.allocated(mesh%tri(it)%blocal_4th)) &
         allocate(mesh%tri(it)%blocal_4th(m,n,mesh%tri(it)%stencil_number_4th))
    end do

    allocate(tmp_data_exchange%f(max_num_4d_global,mesh%nt_full), source=0._r8)!stencil_index_4th
    allocate(tmp_data_exchange_bl%f(m*n*max_num_4d_global,mesh%nt_full), source=0._r8)!blocal_4th
    do it=1,mesh%nt_compute
      do k=1,mesh%tri(it)%stencil_number_4th
        do j=1,n
          do i=1,m
            tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it) = mesh%tri(it)%blocal_4th(i,j,k)
          end do
        end do
        g_index = mesh%t_index(mesh%tri(it)%stencil_index_4th(k))
        tmp_data_exchange%f(k,it)=dble(g_index)
      end do
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_bl)
    call exchange_data_2d(mesh%local_block,field_head_2d)
 
    do it=mesh%nt_compute+1,mesh%nt_full
      do k=1,mesh%tri(it)%stencil_number_4th
        do j=1,n
          do i=1,m
            mesh%tri(it)%blocal_4th(i,j,k) = tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it)
          end do
        end do
        g_index = int(tmp_data_exchange%f(k,it))
        mesh%tri(it)%stencil_index_4th(k) = mesh%map_g2l_t%find(g_index)
      end do
    end do

    deallocate(tmp_data_exchange%f)!stencil_index_4th
    deallocate(tmp_data_exchange_bl%f)!blocal_4th

    mesh%ne = mesh%ne_compute
    DO ie = 1, mesh%ne  ! store for each hx's edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_4th_dual_cell(mesh%tri(v1_index)%stencil_number_4th))
!       allocate(mesh%edp(ie)%v1_weight_4th_dual_cell(mesh%tri(v2_index)%stencil_number_4th))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_4th) &
                     = mesh%tri(v1_index)%blocal_4th(it_local,10,:)
             mesh%edp(ie)%v_weight_4th_dual_cell(1,mesh%tri(v1_index)%stencil_number_4th+1:26) = 0._r8
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v1, stop")

       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_4th) &
                     = mesh%tri(v2_index)%blocal_4th(it_local,10,:)
             mesh%edp(ie)%v_weight_4th_dual_cell(2,mesh%tri(v2_index)%stencil_number_4th+1:26) = 0._r8
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v2, stop")

    END DO
    mesh%ne = mesh%ne_full

    !exchange halo data
    allocate(tmp_data_exchange%f(2*26,mesh%ne_full))
    do ie=1,mesh%ne_compute
      do k=1,26
        do j=1,2
          tmp_data_exchange%f(j+(k-1)*2,ie) = mesh%edp(ie)%v_weight_4th_dual_cell(j,k)
        end do
      end do
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do ie=mesh%ne_compute+1,mesh%ne_full
      do k=1,26
        do j=1,2
          mesh%edp(ie)%v_weight_4th_dual_cell(j,k) = tmp_data_exchange%f(j+(k-1)*2,ie)
        end do
      end do
    end do
    deallocate(tmp_data_exchange%f)
 
    DO it = 1, mesh%nt
      deallocate(mesh%tri(it)%blocal_4th)
    END DO
  end subroutine init_blocal_dual_cell_4th

!================================================
!     B array for 2D polynomial on TR, MBS
!================================================

  subroutine init_blocal_dual_cell_mbs_global(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp

    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny

    real(r8)                            :: xy_list(5) ! x, y, x2, xy, y2

    integer                             :: p1_index
    integer                             :: it         ! global index
    integer                             :: ie         ! global index
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index

!================================================
!          Arrays for Linear Algebra
!================================================

    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin
 
    n = 5   ! Xi, Yi, Xi^2, XiYi, Yi^2
 
    if(mpi_rank() == 0) print*,"global init_blocal_dual_cell_mbs"
    DO it = 1, mesh%nt   ! traverse each dual cell

!================================================
!            common for this tr cell
!================================================

       call create_index_list_tr_mbs(mesh,it)
       num_point_of_stencil = mesh%tri(it)%stencil_number_2nd

       op0     = mesh%tri(it)%c%p      ! OP0 vector (home cell)
       m       = 3                    ! 3 directions for LOB

       if(.not.allocated(mesh%tri(it)%blocal)) allocate(mesh%tri(it)%blocal(m,n,num_point_of_stencil))
       if(.not.allocated(b_local))             allocate(b_local(n,num_point_of_stencil))
       if(.not.allocated(p_local))             allocate(p_local(num_point_of_stencil,n))
       if(.not.allocated(mesh%tri(it)%lob_nx)) allocate(mesh%tri(it)%lob_nx(m,3))
       if(.not.allocated(mesh%tri(it)%lob_ny)) allocate(mesh%tri(it)%lob_ny(m,3))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)         ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          mesh%tri(it)%lob_nx(icell,:) = local_nx
          mesh%tri(it)%lob_ny(icell,:) = local_ny

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_2nd(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p2nd(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

!================================================
!      calculate B using P for current TR
!================================================

          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal(icell,:,:) = b_local

          ! calculate each neighbouring edge's midpoint location in LOB
          do it_local = 1,3

              p1_index       = mesh%tri(it)%ed(it_local)
              op1(1:3)       = mesh%edt(p1_index)%c%p(1:3)      ! ed1's midpoint vector
              call project_sphere2tplane(op1,op0,p0p1_tp)

            ! the 1st dim stores the direction of EACH LOB
            ! the 2nd dim stores each neighbouring edge of tr(it)
            ! the 3rd dim stores x and y locations in LOB
            ! e.g., (2,3,1) means the x location of 3rd neighbouring edge of
            ! tr(it) under the LOB formed by its 2th neighouring point
              mesh%tri(it)%mid_point_lob(icell,it_local,1) = dot_product(p0p1_tp,local_nx) !x
              mesh%tri(it)%mid_point_lob(icell,it_local,2) = dot_product(p0p1_tp,local_ny) !y
          end do

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing triangle 

 
    DO ie = 1, mesh%ne  ! store for each hx'S edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_dual_cell(mesh%tri(v1_index)%stencil_number_2nd))
!       allocate(mesh%edp(ie)%v1_weight_dual_cell(mesh%tri(v2_index)%stencil_number_2nd))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_2nd) &
                     = mesh%tri(v1_index)%blocal(it_local,3,:)
             mesh%edt(ie)%local_index_dual_cell(1)        = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(1,:)  = mesh%tri(v1_index)%mid_point_lob(1,it_local,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell_mbs, v1, stop")

       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_2nd) &
                     = mesh%tri(v2_index)%blocal(it_local,3,:)
             mesh%edt(ie)%local_index_dual_cell(2)        = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(2,:) = mesh%tri(v2_index)%mid_point_lob(1,it_local,:) 
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell_mbs, v2, stop")

    END DO
 

    DO it = 1, mesh%nt ! Only blocal(1,:,:) needs to be stored
       allocate(b_local(n,mesh%tri(it)%stencil_number_2nd))
       b_local(:,:) = mesh%tri(it)%blocal(1,:,:)
       deallocate(mesh%tri(it)%blocal)
       allocate(mesh%tri(it)%blocal(1,n,mesh%tri(it)%stencil_number_2nd))
       mesh%tri(it)%blocal(1,:,:) = b_local(:,:)
       deallocate(b_local)
    END DO

  end subroutine init_blocal_dual_cell_mbs_global

  subroutine init_blocal_dual_cell_mbs(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp

    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny

    real(r8)                            :: xy_list(5) ! x, y, x2, xy, y2

    integer                             :: p1_index
    integer                             :: it         ! global index
    integer                             :: ie         ! global index
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index

!================================================
!          Arrays for Linear Algebra
!================================================

    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin
    
    type(scalar_2d_field)               :: tmp_data_exchange,tmp_data_exchange_bl,&
                                           tmp_data_exchange_x,tmp_data_exchange_y, &
                                           tmp_data_exchange1
    integer(i4)                         :: max_num_2nd,i,j,k,g_index  !max_num_2nd = 16?
 
    n = 5   ! Xi, Yi, Xi^2, XiYi, Yi^2
    m = 3   ! 3 directions for LOB
    max_num_2nd = 16
    
    if(mpi_rank() == 0) print*,"local init_blocal_dual_cell_mbs"
    allocate(tmp_data_exchange%f(1,mesh%nt_full))
    allocate(tmp_data_exchange1%f(max_num_2nd,mesh%nt_full), source=0._r8)!stencil_index_2nd
    allocate(tmp_data_exchange_bl%f(max_num_2nd*m*n,mesh%nt_full), source=0._r8)!mesh%tri(it)%blocal
    allocate(tmp_data_exchange_x%f(3*m,mesh%nt_full))!mesh%tri(it)%lob_nx
    allocate(tmp_data_exchange_y%f(3*m,mesh%nt_full))!mesh%tri(it)%lob_ny

#ifdef LAM_DOMAIN
    mesh%nt = mesh%nt_halo(1)
#else
    mesh%nt = mesh%nt_compute
#endif
    DO it = 1, mesh%nt   ! traverse each dual cell

!================================================
!            common for this tr cell
!================================================

       call create_index_list_tr_mbs(mesh,it)
       num_point_of_stencil = mesh%tri(it)%stencil_number_2nd

       op0     = mesh%tri(it)%c%p      ! OP0 vector (home cell)

       if(.not.allocated(mesh%tri(it)%blocal)) allocate(mesh%tri(it)%blocal(m,n,num_point_of_stencil))
       if(.not.allocated(b_local))             allocate(b_local(n,num_point_of_stencil))
       if(.not.allocated(p_local))             allocate(p_local(num_point_of_stencil,n))
       if(.not.allocated(mesh%tri(it)%lob_nx)) allocate(mesh%tri(it)%lob_nx(m,3))
       if(.not.allocated(mesh%tri(it)%lob_ny)) allocate(mesh%tri(it)%lob_ny(m,3))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)         ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          mesh%tri(it)%lob_nx(icell,:) = local_nx
          mesh%tri(it)%lob_ny(icell,:) = local_ny

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_2nd(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p2nd(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

!================================================
!      calculate B using P for current TR
!================================================

          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal(icell,:,:) = b_local

          ! calculate each neighbouring edge's midpoint location in LOB
          do it_local = 1,3

              p1_index       = mesh%tri(it)%ed(it_local)
              op1(1:3)       = mesh%edt(p1_index)%c%p(1:3)      ! ed1's midpoint vector
              call project_sphere2tplane(op1,op0,p0p1_tp)

            ! the 1st dim stores the direction of EACH LOB
            ! the 2nd dim stores each neighbouring edge of tr(it)
            ! the 3rd dim stores x and y locations in LOB
            ! e.g., (2,3,1) means the x location of 3rd neighbouring edge of
            ! tr(it) under the LOB formed by its 2th neighouring point
              mesh%tri(it)%mid_point_lob(icell,it_local,1) = dot_product(p0p1_tp,local_nx) !x
              mesh%tri(it)%mid_point_lob(icell,it_local,2) = dot_product(p0p1_tp,local_ny) !y
          end do

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)

       ! for halo exchange
       tmp_data_exchange%f(1,it)=mesh%tri(it)%stencil_number_2nd
       do k=1,mesh%tri(it)%stencil_number_2nd
         do j=1,n
           do i=1,m
             tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it) = mesh%tri(it)%blocal(i,j,k)
           end do
         end do
         g_index = mesh%t_index(mesh%tri(it)%stencil_index_2nd(k))
         tmp_data_exchange1%f(k,it)=dble(g_index)
       end do
       do j=1,3
         do i=1,m
           tmp_data_exchange_x%f(i+(j-1)*m,it) = mesh%tri(it)%lob_nx(i,j)
           tmp_data_exchange_y%f(i+(j-1)*m,it) = mesh%tri(it)%lob_ny(i,j)
         end do
       end do 
    END DO ! end of traversing triangle 
    mesh%nt = mesh%nt_full

    !exchange halo data
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    !max_num_2nd = 0
    do it=mesh%nt_compute+1,mesh%nt_full
       mesh%tri(it)%stencil_number_2nd = int(tmp_data_exchange%f(1,it))
       !if(mesh%tri(it)%stencil_number_2nd > max_num_2nd) max_num_2nd = mesh%tri(it)%stencil_number_2nd
#ifdef LAM_DOMAIN
       if (mesh%tri(it)%stencil_number_2nd .le. 0 .or. mesh%tri(it)%stencil_number_2nd .gt. max_num_2nd) then
           mesh%tri(it)%stencil_number_2nd = 1
       end if
#endif
    end do
    deallocate(tmp_data_exchange%f)

    do it = mesh%nt_compute+1,mesh%nt_full
       if(.not.allocated(mesh%tri(it)%blocal))           allocate(mesh%tri(it)%blocal(m,n,mesh%tri(it)%stencil_number_2nd))
       if(.not.allocated(mesh%tri(it)%lob_nx))           allocate(mesh%tri(it)%lob_nx(m,3))
       if(.not.allocated(mesh%tri(it)%lob_ny))           allocate(mesh%tri(it)%lob_ny(m,3))
       if(.not.allocated(mesh%tri(it)%stencil_index_2nd)) &
          allocate(mesh%tri(it)%stencil_index_2nd(mesh%tri(it)%stencil_number_2nd))
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange1)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_bl)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_x)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_y)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do it=mesh%nt_compute+1,mesh%nt_full
      do k=1,mesh%tri(it)%stencil_number_2nd
        do j=1,n
          do i=1,m
             mesh%tri(it)%blocal(i,j,k) = tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it)
          end do
        end do 
        g_index = int(tmp_data_exchange1%f(k,it))
        mesh%tri(it)%stencil_index_2nd(k) = mesh%map_g2l_t%find(g_index)
      end do
      do j=1,3
        do i=1,m
          mesh%tri(it)%lob_nx(i,j) = tmp_data_exchange_x%f(i+(j-1)*m,it)
          mesh%tri(it)%lob_ny(i,j) = tmp_data_exchange_y%f(i+(j-1)*m,it)
        end do
      end do 
    end do

    deallocate(tmp_data_exchange1%f)
    deallocate(tmp_data_exchange_bl%f)
    deallocate(tmp_data_exchange_x%f)
    deallocate(tmp_data_exchange_y%f)
 
    allocate(tmp_data_exchange%f(2,mesh%ne_full))!local_index_dual_cell
    allocate(tmp_data_exchange_x%f(max_num_2nd*2,mesh%ne_full), source=0._r8)!v_weight_dual_cell
    allocate(tmp_data_exchange_y%f(2*2,mesh%ne_full))!v_mid_point_lob_dual_cell
#ifdef LAM_DOMAIN
    mesh%ne = mesh%ne_halo(1)
#else
    mesh%ne = mesh%ne_compute
#endif
    DO ie = 1, mesh%ne  ! store for each hx'S edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_dual_cell(mesh%tri(v1_index)%stencil_number_2nd))
!       allocate(mesh%edp(ie)%v1_weight_dual_cell(mesh%tri(v2_index)%stencil_number_2nd))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_2nd) &
                     = mesh%tri(v1_index)%blocal(it_local,3,:)
             mesh%edp(ie)%v_weight_dual_cell(1,mesh%tri(v1_index)%stencil_number_2nd+1:16) = 0._r8 
             mesh%edt(ie)%local_index_dual_cell(1)        = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(1,:)  = mesh%tri(v1_index)%mid_point_lob(1,it_local,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell_mbs, v1, stop")

       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_2nd) &
                     = mesh%tri(v2_index)%blocal(it_local,3,:)
             mesh%edp(ie)%v_weight_dual_cell(2,mesh%tri(v2_index)%stencil_number_2nd+1:16) = 0._r8
             mesh%edt(ie)%local_index_dual_cell(2)        = it_local
             mesh%edt(ie)%v_mid_point_lob_dual_cell(2,:) = mesh%tri(v2_index)%mid_point_lob(1,it_local,:) 
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell_mbs, v2, stop")

       do j=1,2
         do i=1,2
           tmp_data_exchange_y%f(i+(j-1)*2,ie) = mesh%edt(ie)%v_mid_point_lob_dual_cell(i,j)
         end do
         tmp_data_exchange%f(j,ie) = dble(mesh%edt(ie)%local_index_dual_cell(j))
       end do
       do j=1, max_num_2nd
         do i=1,2
           tmp_data_exchange_x%f(i+(j-1)*2,ie) = mesh%edp(ie)%v_weight_dual_cell(i,j)
         end do
       end do
    END DO
    mesh%ne = mesh%ne_full
    
    !exchange halo data

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_x)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_y)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    do ie=mesh%ne_compute+1,mesh%ne_full
      do j=1,2
        do i=1,2
          mesh%edt(ie)%v_mid_point_lob_dual_cell(i,j) = tmp_data_exchange_y%f(i+(j-1)*2,ie)
        end do
        mesh%edt(ie)%local_index_dual_cell(j) = int(tmp_data_exchange%f(j,ie))
      end do
      do j=1, max_num_2nd
        do i=1,2
          mesh%edp(ie)%v_weight_dual_cell(i,j) = tmp_data_exchange_x%f(i+(j-1)*2,ie)
        end do
      end do
    end do

    deallocate(tmp_data_exchange%f)
    deallocate(tmp_data_exchange_x%f)
    deallocate(tmp_data_exchange_y%f)
    
    allocate(b_local(n,max_num_2nd))
    DO it = 1, mesh%nt ! Only blocal(1,:,:) needs to be stored
       b_local(:,1:mesh%tri(it)%stencil_number_2nd) = mesh%tri(it)%blocal(1,:,:)
       deallocate(mesh%tri(it)%blocal)
       allocate(mesh%tri(it)%blocal(1,n,mesh%tri(it)%stencil_number_2nd))
       mesh%tri(it)%blocal(1,:,:) = b_local(:,1:mesh%tri(it)%stencil_number_2nd)
    END DO
    deallocate(b_local)

    allocate(mesh%edp_v_weight_dual_cell(16,2,mesh%ne_full))!v_weight_prime_cell
    allocate(mesh%tri_stencil_number_2nd(mesh%nt_full))
    allocate(mesh%tri_stencil_index_2nd(16,mesh%nt_full))
    do it = 1, mesh%nt_full
      mesh%tri_stencil_number_2nd(it) = mesh%tri(it)%stencil_number_2nd
      mesh%tri_stencil_index_2nd(1:mesh%tri(it)%stencil_number_2nd, it) = &
           mesh%tri(it)%stencil_index_2nd(1:mesh%tri(it)%stencil_number_2nd)
    end do
    do ie = 1, mesh%ne_full
      do k = 1, 2
        mesh%edp_v_weight_dual_cell(1:16, k, ie) = &
          mesh%edp(ie)%v_weight_dual_cell(k,1:16)
      end do
    end do
  end subroutine init_blocal_dual_cell_mbs
  
!================================================
!     B array for 4D polynomial on TR, MBS
!================================================

  subroutine init_blocal_dual_cell_4th_mbs_global(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp

    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny

    real(r8)                            :: xy_list(14)  !
    integer                             :: p1_index
    integer                             :: it        ! global index 
    integer                             :: ie        ! global index 
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin
!
! fixed
!

    n = 14
    m = 3   ! 3 directions for LOB
    if(mpi_rank() == 0) print*,"global init_blocal_dual_cell_4th_mbs"

    DO it = 1, mesh%nt   ! traverse each dual cell

!================================================
!            common for this tr cell
!================================================
       call create_index_list_tr4th_mbs(mesh,it)
       num_point_of_stencil = mesh%tri(it)%stencil_number_4th

       op0     = mesh%tri(it)%c%p      ! OP0 vector (home cell)

       if(.not.allocated(b_local))                 allocate(b_local(n,num_point_of_stencil))
       if(.not.allocated(p_local))                 allocate(p_local(num_point_of_stencil,n))
       if(.not.allocated(mesh%tri(it)%blocal_4th)) allocate(mesh%tri(it)%blocal_4th(m,n,num_point_of_stencil))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)    ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_4th(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p4th(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

!================================================
!      calculate B using P for current TR
!================================================

          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal_4th(icell,:,:) = b_local

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing triangle 
 
    DO ie = 1, mesh%ne  ! store for each hx's edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_4th_dual_cell(mesh%tri(v1_index)%stencil_number_4th))
!       allocate(mesh%edp(ie)%v1_weight_4th_dual_cell(mesh%tri(v2_index)%stencil_number_4th))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_2nd) &
                     = mesh%tri(v1_index)%blocal_4th(it_local,10,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v1, stop")


       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_2nd) &
                     = mesh%tri(v2_index)%blocal_4th(it_local,10,:)
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v2, stop")

    END DO

    DO it = 1, mesh%nt
      deallocate(mesh%tri(it)%blocal_4th)
    END DO
  end subroutine init_blocal_dual_cell_4th_mbs_global
  
  subroutine init_blocal_dual_cell_4th_mbs(mesh)
! io
    type(global_domain), intent(inout)  :: mesh
!================================================
!            Vars for local basis
!================================================
    real(r8), dimension(1:3)            :: op0
    real(r8), dimension(1:3)            :: op1
    real(r8), dimension(1:3)            :: p0p1_tp

    real(r8), dimension(1:3)            :: local_nx
    real(r8), dimension(1:3)            :: local_ny

    real(r8)                            :: xy_list(14)  !
    integer                             :: p1_index
    integer                             :: it        ! global index 
    integer                             :: ie        ! global index 
    integer                             :: icell
    integer                             :: it_local
    integer                             :: v1_index
    integer                             :: v2_index
!================================================
!          Arrays for Linear Algebra
!================================================
    real(r8), allocatable               :: p_local(:,:) ! P
    real(r8), allocatable               :: b_local(:,:) ! B
    integer(i4)                         :: m            ! number of directions
    integer(i4)                         :: n            ! 5 for 2d poly
    integer(i4)                         :: tmp
    integer(i4)                         :: num_point_of_stencil   ! how many points in the stencil except origin

    type(scalar_2d_field)               :: tmp_data_exchange,tmp_data_exchange_bl
    integer(i4)                         :: max_num_4d,i,j,k,g_index, max_num_4d_global, err
!
! fixed
!

    n = 14
    m = 3   ! 3 directions for LOB

    if(mpi_rank() == 0) print*,"local init_blocal_dual_cell_4th_mbs"
    mesh%nt = mesh%nt_compute
    DO it = 1, mesh%nt   ! traverse each dual cell

!================================================
!            common for this tr cell
!================================================
       call create_index_list_tr4th_mbs(mesh,it)
       num_point_of_stencil = mesh%tri(it)%stencil_number_4th

       op0     = mesh%tri(it)%c%p      ! OP0 vector (home cell)

       if(.not.allocated(b_local))                 allocate(b_local(n,num_point_of_stencil))
       if(.not.allocated(p_local))                 allocate(p_local(num_point_of_stencil,n))
       if(.not.allocated(mesh%tri(it)%blocal_4th)) allocate(mesh%tri(it)%blocal_4th(m,n,num_point_of_stencil))

!================================================
!           Build LOB using each nbrs
!================================================

       do icell = 1, m

          p1_index   = mesh%tri(it)%nb(icell)         ! using P1 to form a local basis (C0 in my paper)
          op1(1:3)   = mesh%tri(p1_index)%c%p(1:3)    ! OP1 vector (C1 in my paper)

          call project_sphere2tplane(op1,op0,p0p1_tp)
          call get_lob_xy(op0,p0p1_tp,local_nx,local_ny)

          ! for each cell in the list, create P array
          do it_local = 1, num_point_of_stencil
             p1_index            = mesh%tri(it)%stencil_index_4th(it_local)
             op1(1:3)            = mesh%tri(p1_index)%c%p(1:3)
             call project_sphere2tplane(op1,op0,p0p1_tp)
             call set_xy_p4th(p0p1_tp,local_nx,local_ny,xy_list)
             p_local(it_local,:) = xy_list
          end do

!================================================
!      calculate B using P for current TR
!================================================

          call calc_array_b(num_point_of_stencil,n,p_local,b_local)
          mesh%tri(it)%blocal_4th(icell,:,:) = b_local

       end do ! loop of each direction

       deallocate(p_local)
       deallocate(b_local)

    END DO ! end of traversing triangle 
    mesh%nt = mesh%nt_full
    
    !exchange halo data
    allocate(tmp_data_exchange%f(1,mesh%nt_full))
    do it=1,mesh%nt_compute
       tmp_data_exchange%f(1,it)=dble(mesh%tri(it)%stencil_number_4th)
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d(mesh%local_block,field_head_2d)

    max_num_4d = 0
    do it=1,mesh%nt_full
       if(it > mesh%nt_compute) mesh%tri(it)%stencil_number_4th = int(tmp_data_exchange%f(1,it))
       if(mesh%tri(it)%stencil_number_4th > max_num_4d) max_num_4d = mesh%tri(it)%stencil_number_4th
    end do
    deallocate(tmp_data_exchange%f)
    call MPI_ALLREDUCE(max_num_4d, max_num_4d_global, 1, MPI_INT, MPI_MAX, grist_comm(), err)
    
    do it=mesh%nt_compute+1,mesh%nt_full
      if(.not.allocated(mesh%tri(it)%stencil_index_4th)) &
        allocate(mesh%tri(it)%stencil_index_4th(mesh%tri(it)%stencil_number_4th))
      if(.not.allocated(mesh%tri(it)%blocal_4th)) &
        allocate(mesh%tri(it)%blocal_4th(m,n,mesh%tri(it)%stencil_number_4th))
    end do

    allocate(tmp_data_exchange%f(max_num_4d_global,mesh%nt_full), source=0._r8)!stencil_index_4th
    allocate(tmp_data_exchange_bl%f(m*n*max_num_4d_global,mesh%nt_full), source=0._r8)!blocal_4th
    do it=1,mesh%nt_compute
      do k=1,mesh%tri(it)%stencil_number_4th
        do j=1,n
          do i=1,m
            tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it) = mesh%tri(it)%blocal_4th(i,j,k)
          end do
        end do
        g_index = mesh%t_index(mesh%tri(it)%stencil_index_4th(k))
        tmp_data_exchange%f(k,it)=dble(g_index)
      end do
    end do
    
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange_bl)
    call exchange_data_2d(mesh%local_block,field_head_2d)
    
    do it=mesh%nt_compute+1,mesh%nt_full
      do k=1,mesh%tri(it)%stencil_number_4th
        do j=1,n
          do i=1,m
            mesh%tri(it)%blocal_4th(i,j,k) = tmp_data_exchange_bl%f(i+(j-1)*m+(k-1)*m*n,it)
          end do
      end do
      g_index = int(tmp_data_exchange%f(k,it))
      mesh%tri(it)%stencil_index_4th(k) = mesh%map_g2l_t%find(g_index)
      end do
    end do
    
    deallocate(tmp_data_exchange%f)!stencil_index_4th
    deallocate(tmp_data_exchange_bl%f)!blocal_4th
    
    mesh%ne = mesh%ne_compute
    DO ie = 1, mesh%ne  ! store for each hx's edge

       v1_index = mesh%edp(ie)%v(1)
       v2_index = mesh%edp(ie)%v(2)

!       allocate(mesh%edp(ie)%v0_weight_4th_dual_cell(mesh%tri(v1_index)%stencil_number_4th))
!       allocate(mesh%edp(ie)%v1_weight_4th_dual_cell(mesh%tri(v2_index)%stencil_number_4th))

       ! TRANSVERSE V1
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v1_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(1,1:mesh%tri(v1_index)%stencil_number_2nd) &
                     = mesh%tri(v1_index)%blocal_4th(it_local,10,:)
             mesh%edp(ie)%v_weight_4th_dual_cell(1,mesh%tri(v1_index)%stencil_number_2nd+1:26) = 0._r8
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v1, stop")


       ! Transverse V2
       tmp = -999
       do it_local = 1, 3
          if(mesh%tri(v2_index)%ed(it_local) .eq. ie)then
             mesh%edp(ie)%v_weight_4th_dual_cell(2,1:mesh%tri(v2_index)%stencil_number_2nd) &
                     = mesh%tri(v2_index)%blocal_4th(it_local,10,:)
             mesh%edp(ie)%v_weight_4th_dual_cell(2,mesh%tri(v2_index)%stencil_number_2nd+1:26) = 0._r8
             tmp = 0
          end if
       end do

       if(tmp.eq.-999) call endrun("we have problem in init_blocal_dual_cell, v2, stop")

    END DO
    mesh%ne = mesh%ne_full

    !exchange halo data
    allocate(tmp_data_exchange%f(2*26,mesh%ne_full))
    do ie=1,mesh%ne_compute
      do k=1,26
        do j=1,2
          tmp_data_exchange%f(j+(k-1)*2,ie) = mesh%edp(ie)%v_weight_4th_dual_cell(j,k)
        end do
      end do
    end do

    call exchange_data_2d_add(mesh,field_head_2d,tmp_data_exchange)
    call exchange_data_2d(mesh%local_block,field_head_2d)
    
    do ie=mesh%ne_compute+1,mesh%ne_full
      do k=1,26
        do j=1,2
          mesh%edp(ie)%v_weight_4th_dual_cell(j,k) = tmp_data_exchange%f(j+(k-1)*2,ie)
        end do
      end do
    end do
    deallocate(tmp_data_exchange%f)
 
    DO it = 1, mesh%nt
       deallocate(mesh%tri(it)%blocal_4th)
    END DO
  end subroutine init_blocal_dual_cell_4th_mbs

  subroutine init_plg_sub_triangle_info_global(mesh)
! io
    type(global_domain),  intent(inout)  :: mesh
! local
    integer(i4)                          :: iv, inb, v0, tr1, tr2
    real(r8), dimension(1:3)             :: p0p1, p0p2, p0p3
    real(r8), dimension(1:3)             :: tmparr1, tmparr2, tmparr3

    if(mpi_rank() == 0) print *, "global init_plg_sub_triangle_info"

    do iv = 1, mesh%nv  ! for each cell

      if(.not.allocated(mesh%plg(iv)%sub_triangle_midp_3d)) allocate( mesh%plg(iv)%sub_triangle_midp_3d(mesh%vtx(iv)%nnb,3,3))
      if(.not.allocated(mesh%plg(iv)%sub_triangle_midp_lob))allocate(mesh%plg(iv)%sub_triangle_midp_lob(mesh%vtx(iv)%nnb,3,2))
      if(.not.allocated(mesh%plg(iv)%sub_triangle_area))    allocate(mesh%plg(iv)%sub_triangle_area(mesh%vtx(iv)%nnb))

      do inb = 1, mesh%vtx(iv)%nnb
        tr1 = mesh%vtx(iv)%tr(inb)   ! tri index
        if(inb+1.gt.mesh%vtx(iv)%nnb)then
            tr2 = mesh%vtx(iv)%tr(1) ! tri index
        else
            tr2 = mesh%vtx(iv)%tr(inb+1)
        end if

        mesh%plg(iv)%sub_triangle_midp_3d(inb,1,1:3)=(mesh%vtx(iv)%p+mesh%tri(tr1)%c%p)/norm(mesh%vtx(iv)%p+mesh%tri(tr1)%c%p)
        mesh%plg(iv)%sub_triangle_midp_3d(inb,2,1:3)=(mesh%vtx(iv)%p+mesh%tri(tr2)%c%p)/norm(mesh%vtx(iv)%p+mesh%tri(tr2)%c%p)
        mesh%plg(iv)%sub_triangle_midp_3d(inb,3,1:3)=(mesh%tri(tr1)%c%p+mesh%tri(tr2)%c%p)/norm(mesh%tri(tr1)%c%p+mesh%tri(tr2)%c%p)

        tmparr1(1:3) = mesh%plg(iv)%sub_triangle_midp_3d(inb,1,1:3)
        tmparr2(1:3) = mesh%plg(iv)%sub_triangle_midp_3d(inb,2,1:3)
        tmparr3(1:3) = mesh%plg(iv)%sub_triangle_midp_3d(inb,3,1:3)

        call project_sphere2tplane(tmparr1, mesh%vtx(iv)%p, p0p1)
        call project_sphere2tplane(tmparr2, mesh%vtx(iv)%p, p0p2)
        call project_sphere2tplane(tmparr3, mesh%vtx(iv)%p, p0p3)
! x-component in lob
        mesh%plg(iv)%sub_triangle_midp_lob(inb,1,1)= dot_product(mesh%vtx(iv)%lob_nx(1,1:3), p0p1)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,2,1)= dot_product(mesh%vtx(iv)%lob_nx(1,1:3), p0p2)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,3,1)= dot_product(mesh%vtx(iv)%lob_nx(1,1:3), p0p3)
! y-component in lob
        mesh%plg(iv)%sub_triangle_midp_lob(inb,1,2)= dot_product(mesh%vtx(iv)%lob_ny(1,1:3), p0p1)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,2,2)= dot_product(mesh%vtx(iv)%lob_ny(1,1:3), p0p2)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,3,2)= dot_product(mesh%vtx(iv)%lob_ny(1,1:3), p0p3)

        mesh%plg(iv)%sub_triangle_area(inb)    = sphtriarea(mesh%vtx(iv)%p, mesh%tri(tr1)%c%p, mesh%tri(tr2)%c%p)
      end do
    end do

    return
  end subroutine init_plg_sub_triangle_info_global
  
  subroutine init_plg_sub_triangle_info(mesh)
! io
    type(global_domain),  intent(inout)  :: mesh
! local
    integer(i4)                          :: iv, inb, v0, tr1, tr2
    real(r8), dimension(1:3)             :: p0p1, p0p2, p0p3
    real(r8), dimension(1:3)             :: tmparr1, tmparr2, tmparr3

    if(mpi_rank() == 0) print *, "local init_plg_sub_triangle_info"

    do iv = 1, mesh%nv  ! for each cell

      if(.not.allocated(mesh%plg(iv)%sub_triangle_midp_3d)) allocate( mesh%plg(iv)%sub_triangle_midp_3d(mesh%vtx(iv)%nnb,3,3))
      if(.not.allocated(mesh%plg(iv)%sub_triangle_midp_lob))allocate(mesh%plg(iv)%sub_triangle_midp_lob(mesh%vtx(iv)%nnb,3,2))
      if(.not.allocated(mesh%plg(iv)%sub_triangle_area))    allocate(mesh%plg(iv)%sub_triangle_area(mesh%vtx(iv)%nnb))

      do inb = 1, mesh%vtx(iv)%nnb
        tr1 = mesh%vtx(iv)%tr(inb)   ! tri index
        if(inb+1.gt.mesh%vtx(iv)%nnb)then
            tr2 = mesh%vtx(iv)%tr(1) ! tri index
        else
            tr2 = mesh%vtx(iv)%tr(inb+1)
        end if

        mesh%plg(iv)%sub_triangle_midp_3d(inb,1,1:3)=(mesh%vtx(iv)%p+mesh%tri(tr1)%c%p)/norm(mesh%vtx(iv)%p+mesh%tri(tr1)%c%p)
        mesh%plg(iv)%sub_triangle_midp_3d(inb,2,1:3)=(mesh%vtx(iv)%p+mesh%tri(tr2)%c%p)/norm(mesh%vtx(iv)%p+mesh%tri(tr2)%c%p)
        mesh%plg(iv)%sub_triangle_midp_3d(inb,3,1:3)=(mesh%tri(tr1)%c%p+mesh%tri(tr2)%c%p)/norm(mesh%tri(tr1)%c%p+mesh%tri(tr2)%c%p)

        tmparr1(1:3) = mesh%plg(iv)%sub_triangle_midp_3d(inb,1,1:3)
        tmparr2(1:3) = mesh%plg(iv)%sub_triangle_midp_3d(inb,2,1:3)
        tmparr3(1:3) = mesh%plg(iv)%sub_triangle_midp_3d(inb,3,1:3)

        call project_sphere2tplane(tmparr1, mesh%vtx(iv)%p, p0p1)
        call project_sphere2tplane(tmparr2, mesh%vtx(iv)%p, p0p2)
        call project_sphere2tplane(tmparr3, mesh%vtx(iv)%p, p0p3)
! x-component in lob
        mesh%plg(iv)%sub_triangle_midp_lob(inb,1,1)= dot_product(mesh%vtx(iv)%lob_nx(1,1:3), p0p1)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,2,1)= dot_product(mesh%vtx(iv)%lob_nx(1,1:3), p0p2)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,3,1)= dot_product(mesh%vtx(iv)%lob_nx(1,1:3), p0p3)
! y-component in lob
        mesh%plg(iv)%sub_triangle_midp_lob(inb,1,2)= dot_product(mesh%vtx(iv)%lob_ny(1,1:3), p0p1)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,2,2)= dot_product(mesh%vtx(iv)%lob_ny(1,1:3), p0p2)
        mesh%plg(iv)%sub_triangle_midp_lob(inb,3,2)= dot_product(mesh%vtx(iv)%lob_ny(1,1:3), p0p3)

        mesh%plg(iv)%sub_triangle_area(inb)    = sphtriarea(mesh%vtx(iv)%p, mesh%tri(tr1)%c%p, mesh%tri(tr2)%c%p)
      end do
    end do

    return
  end subroutine init_plg_sub_triangle_info

!-------------------------------------------------------------
!  Private subroutines
!      create_index_list_xxx
!      calc_array_b 
!      project_xxx
!      get_lob_xxx
!      set_xy_xxx
!-------------------------------------------------------------

!================================================
!  Create index for 2nd-order polynomial of tr,
!  an old way
!================================================

  subroutine create_index_list_tr(mesh, it)
!
! io
!
    type(global_domain), intent(inout)  :: mesh
    integer(i4)        , intent(in)     :: it
! local
    type(index_list), pointer           :: tr_index_list,tmp_ptr
    integer(i4)                         :: ii, jj, kk
    integer(i4)                         :: idx 
    integer(i4)                         :: nt 

    nt = 0
    tr_index_list=>null()
    do ii = 1, 3
       idx = mesh%tri(it)%nb(ii)
       call check_and_add(tr_index_list,idx,nt)
       do jj = 1, 3
          if(mesh%tri(idx)%nb(jj).ne.it) then
             call check_and_add(tr_index_list, mesh%tri(idx)%nb(jj), nt)
#ifdef PV3SS
! when using a smaller stencil, exit once a unique tr is found, trivial sensitivity in jwbw
              exit
#endif
          end if
       end do
    end do

#ifdef PV3SS
    if(nt.ne.6) then
#else
    if(nt.ne.9) then
#endif
       print*,"error in create_index_list_tr, stop"
       stop
    end if

    mesh%tri(it)%stencil_number_2nd = nt
    if(.not.allocated(mesh%tri(it)%stencil_index_2nd))&
       allocate(mesh%tri(it)%stencil_index_2nd(nt))

    kk = 1
    do while(associated(tr_index_list))
       mesh%tri(it)%stencil_index_2nd(kk) = tr_index_list%gidx
       kk = kk+1
       !tr_index_list=>tr_index_list%next
       tmp_ptr=>tr_index_list%next
       deallocate(tr_index_list)
       tr_index_list=>tmp_ptr
    end do
    !deallocate(tr_index_list)

    return
  end subroutine create_index_list_tr

!================================================
! create index for 2nd-order polynomial of TR, 
! MBS, better way
!================================================

  subroutine create_index_list_tr_mbs(mesh,it)
! io
    type(global_domain), intent(inout)  :: mesh
    integer(i4)        , intent(in)     :: it
! local
    type(index_list), pointer           :: tr_index_list,tmp_ptr
    integer(i4)                         :: ii, jj, kk
    integer(i4)                         :: idx
    integer(i4)                         :: iv
    integer(i4)                         :: nt

    nt = 0
    tr_index_list => null()
    do iv = 1, 3 ! for three vertice of this triangle
      do ii = 1, mesh%vtx(mesh%tri(it)%v(iv))%nnb
         if(mesh%vtx(mesh%tri(it)%v(iv))%tr(ii).ne.it)then
            call check_and_add(tr_index_list,mesh%vtx(mesh%tri(it)%v(iv))%tr(ii),nt)
         end if
      end do
    end do
    mesh%tri(it)%stencil_number_2nd = nt

    if(.not.allocated(mesh%tri(it)%stencil_index_2nd))&
       allocate(mesh%tri(it)%stencil_index_2nd(nt))

    kk = 1
    do while(associated(tr_index_list))
       mesh%tri(it)%stencil_index_2nd(kk) = tr_index_list%gidx
       kk = kk +1
       !tr_index_list=>tr_index_list%next
       tmp_ptr=>tr_index_list%next
       deallocate(tr_index_list)
       tr_index_list=>tmp_ptr
    end do
    !deallocate(tr_index_list)

    return
  end subroutine create_index_list_tr_mbs

!================================================
!  create index for 4th-order polynomial of HX
!================================================

  subroutine create_index_list_hx4th(mesh,iv)
! io
    type(global_domain) , intent(inout) :: mesh
    integer(i4)         , intent(in)    :: iv
! local
    type(index_list), pointer           :: hx4th_index_list,tmp_ptr
    integer(i4)                         :: ii
    integer(i4)                         :: jj
    integer(i4)                         :: kk
    integer(i4)                         :: idx
    integer(i4)                         :: num_point_of_stencil

    num_point_of_stencil = 0
    hx4th_index_list => null()
! first ring
    do ii = 1, mesh%vtx(iv)%nnb  ! inner ring, all nbrs of cell iv
       call check_and_add(hx4th_index_list, mesh%vtx(iv)%nb(ii), num_point_of_stencil)
       idx = mesh%vtx(iv)%nb(ii)
! second ring
       do jj = 1, mesh%vtx(idx)%nnb
          if(mesh%vtx(idx)%nb(jj).ne.iv)&
            call check_and_add(hx4th_index_list, mesh%vtx(idx)%nb(jj), num_point_of_stencil)
       end do
    end do

    mesh%plg(iv)%stencil_number_4th = num_point_of_stencil

    !
    ! This condition only applies for quasi-uniform grid, not ccvt-refined grid
    ! if(num_point_of_stencil.gt.18) &
    !    call endrun("num_point_of_stencil >18, unrealistic, stop")

    if(.not.allocated(mesh%plg(iv)%stencil_index_4th)) &
       allocate(mesh%plg(iv)%stencil_index_4th(num_point_of_stencil))

    kk = 1
    do while(associated(hx4th_index_list))
        mesh%plg(iv)%stencil_index_4th(kk) = hx4th_index_list%gidx
        kk = kk +1
        !hx4th_index_list =>hx4th_index_list%next
        tmp_ptr=>hx4th_index_list%next
        deallocate(hx4th_index_list)
        hx4th_index_list=>tmp_ptr
    end do
    !deallocate(hx4th_index_list)

    return
  end subroutine create_index_list_hx4th

!================================================
!  create index for 4th-order polynomial of TR
!  an old way
!================================================

  subroutine create_index_list_tr4th(mesh,it)
! io
    type(global_domain) , intent(inout) :: mesh
    integer(i4)         , intent(in)    :: it
! local
    integer(i4)                         :: base_vertice(6) ! base vertice is six triangle points on a four-triangle stencil
    integer(i4)                         :: ii, jj, kk, tmp
    integer(i4)                         :: nt, iv
    type(index_list), pointer           :: tr_index_list,tmp_ptr

    tr_index_list => null()
! find three base vertice that does not belong to home cell
    tmp = 0
    do ii = 1, 3  ! for each nb of home cell
      do jj = 1,3 ! for each vertice of each nb of home cell
          if(mesh%tri(mesh%tri(it)%nb(ii))%v(jj).ne.mesh%tri(it)%v(1).and.&
             mesh%tri(mesh%tri(it)%nb(ii))%v(jj).ne.mesh%tri(it)%v(2).and.&
             mesh%tri(mesh%tri(it)%nb(ii))%v(jj).ne.mesh%tri(it)%v(3))then
             base_vertice(ii) = mesh%tri(mesh%tri(it)%nb(ii))%v(jj)
             tmp = tmp +1
          end if
      end do
    end do

    if(tmp.ne.3) call endrun("we have problem in create_index_list_tr4th, stop, tmp.ne.3 ")

! for each base vertice, get all related triangles

    nt = 0
    do iv = 1, 3
      do ii = 1, mesh%vtx(base_vertice(iv))%nnb
         if(mesh%vtx(base_vertice(iv))%tr(ii).ne.it)then
            call check_and_add(tr_index_list,mesh%vtx(base_vertice(iv))%tr(ii),nt)
         end if
      end do
    end do

    mesh%tri(it)%stencil_number_4th = nt

    if(.not.allocated(mesh%tri(it)%stencil_index_4th)) allocate(mesh%tri(it)%stencil_index_4th(nt))

    kk = 1
    do while(associated(tr_index_list))
       mesh%tri(it)%stencil_index_4th(kk) = tr_index_list%gidx
       kk = kk +1
       !tr_index_list=>tr_index_list%next
       tmp_ptr=>tr_index_list%next
       deallocate(tr_index_list)
       tr_index_list=>tmp_ptr
    end do
    !deallocate(tr_index_list)

    return
  end subroutine create_index_list_tr4th

!================================================
!  create index for 4th-order polynomial of TR
!  MBS 
!================================================

  subroutine create_index_list_tr4th_mbs(mesh,it)
! io
    type(global_domain) , intent(inout) :: mesh
    integer(i4)         , intent(in)    :: it
! local
    integer(i4)                         :: base_vertice(6) ! base vertice is six triangle points on a four-triangle stencil
    integer(i4)                         :: ii, jj, kk, tmp
    integer(i4)                         :: nt, iv
    type(index_list), pointer           :: tr_index_list,tmp_ptr

    tr_index_list => null()
!
! find three base vertice that does not belong to home cell
!
    tmp = 0
    do ii = 1, 3  ! for each nb of home cell
      do jj = 1,3 ! for each vertice of each nb of home cell
          if(mesh%tri(mesh%tri(it)%nb(ii))%v(jj).ne.mesh%tri(it)%v(1).and.&
             mesh%tri(mesh%tri(it)%nb(ii))%v(jj).ne.mesh%tri(it)%v(2).and.&
             mesh%tri(mesh%tri(it)%nb(ii))%v(jj).ne.mesh%tri(it)%v(3))then
             base_vertice(ii) = mesh%tri(mesh%tri(it)%nb(ii))%v(jj)
             tmp = tmp +1
          end if
      end do
    end do

    if(tmp.ne.3) call endrun("we have problem in create_index_list_tr4th, stop, tmp.ne.3 ")
    base_vertice(4:6) = mesh%tri(it)%v(1:3)

! for each base vertice, get all related triangles

    nt = 0
    do iv = 1, 6 ! for six vertice of this four-triangle
      do ii = 1, mesh%vtx(base_vertice(iv))%nnb
         if(mesh%vtx(base_vertice(iv))%tr(ii).ne.it)then
            call check_and_add(tr_index_list,mesh%vtx(base_vertice(iv))%tr(ii),nt)
         end if
      end do
    end do

    mesh%tri(it)%stencil_number_4th = nt

    if(.not.allocated(mesh%tri(it)%stencil_index_4th)) allocate(mesh%tri(it)%stencil_index_4th(nt))

    kk = 1
    do while(associated(tr_index_list))
       mesh%tri(it)%stencil_index_4th(kk) = tr_index_list%gidx
       kk = kk +1
       !tr_index_list=>tr_index_list%next
       tmp_ptr=>tr_index_list%next
       deallocate(tr_index_list)
       tr_index_list=>tmp_ptr
    end do
    !deallocate(tr_index_list)

    return
  end subroutine create_index_list_tr4th_mbs

!--------------------------------------
!
! END OF CREATE INDEX LIST ROUTINES
!
!--------------------------------------

!================================================
!         Calculate B array using LAPACK
!================================================

  subroutine calc_array_b0(m,n,p_local,b_local)
! io
    integer(i4),   intent(in)    :: m  ! cell number of stencil except home cell
    integer(i4),   intent(in)    :: n  ! 5 for 2nd, 14 for 4th
    real(r8)   ,   intent(in)    :: p_local(m,n)
    real(r8)   ,   intent(inout) :: b_local(n,m)
! local
    real(r8)                     :: pt_local(n,m)
    real(r8)                     :: ptp(m,n)
    real(r8)                     :: ptp_inv(n,m)
    real(r8)                     :: u_local(m,m)
    real(r8)                     :: s_local(m,n)
    real(r8)                     :: v_local(n,n)
   
    ! pt_local = transpose(p_local)
    ptp      = p_local
    call svd_lapack(m, n, ptp, u_local, s_local, v_local)
    call pseudo_inverse(m, n, u_local, s_local, v_local, ptp_inv)
    b_local  = ptp_inv ! matmul(ptp_inv,pt_local) ! b array with the icell th neighbor

    return
  end subroutine calc_array_b0

  subroutine calc_array_b(m,n,p_local,b_local)
! io
    integer(i4),   intent(in)    :: m  ! cell number of stencil except home cell
    integer(i4),   intent(in)    :: n  ! 5 for 2nd, 14 for 4th
    real(r8)   ,   intent(in)    :: p_local(m,n)
    real(r8)   ,   intent(inout) :: b_local(n,m)
! local
    real(r8)                     :: pt_local(n,m)
    real(r8)                     :: ptp(n,n)
    real(r8)                     :: ptp_inv(n,n)
    real(r8)                     :: u_local(n,n)
    real(r8)                     :: s_local(n,n)
    real(r8)                     :: v_local(n,n)

    pt_local = transpose(p_local)
    ptp      = matmul(pt_local,p_local)
    call svd_lapack(n, n, ptp, u_local, s_local, v_local)
    call pseudo_inverse(n, n, u_local, s_local, v_local, ptp_inv)
    b_local  = matmul(ptp_inv,pt_local) ! b array with the icell th neighbor

    return
  end subroutine calc_array_b

  subroutine calc_array_b1(m,n,p_local,b_local)
! io
    integer(i4),   intent(in)    :: m  ! cell number of stencil except home cell
    integer(i4),   intent(in)    :: n  ! 5 for 2nd, 14 for 4th
    real(r8)   ,   intent(in)    :: p_local(m,n)
    real(r8)   ,   intent(inout) :: b_local(n,m)
! local
    real(r8)                     :: pt_local(n,m)
    real(r8)                     :: ptp(n,n)
    real(r8)                     :: ptp_inv(n,n)
    integer                      :: indx(n)

    pt_local = transpose(p_local)
    ptp      = matmul(pt_local,p_local)
    call migs(ptp,n,ptp_inv,indx)
    b_local  = matmul(ptp_inv,pt_local) ! b array with the icell th neighbor

    return
  end subroutine calc_array_b1

!================================================
!   Project point from sphere to tangent plane
!================================================

  subroutine project_sphere2tplane_r8(op1,op0,p0p1_tp)
! io
    real(r8), dimension(1:3), intent(in)    :: op1      ! 2b projected
    real(r8), dimension(1:3), intent(in)    :: op0      ! home cell
    real(r8), dimension(1:3), intent(inout) :: p0p1_tp  ! projected
! local
    real(r8)                                :: cos_dlamda1
    real(r8)                                :: dlamda1
    real(r8), dimension(1:3)                :: op2
    real(r8), dimension(1:3)                :: p0p2

    cos_dlamda1 = dot_product(op1,op0)/(norm(op0)*norm(op1))
    dlamda1     = acos(cos_dlamda1)              ! radian between op0 and op1
!
! another way to compute dlamda1 is arclen(op1,op0), 
! same to machine precision
!
    op2         = op1/cos_dlamda1
    p0p2        = op2-op0
    p0p1_tp     = p0p2*dlamda1/tan(dlamda1)

    return
  end subroutine project_sphere2tplane_r8

  subroutine project_sphere2tplane_r4(op1,op0,p0p1_tp)
! io
    real(r4), dimension(1:3), intent(in)    :: op1      ! 2b projected
    real(r4), dimension(1:3), intent(in)    :: op0      ! home cell
    real(r4), dimension(1:3), intent(inout) :: p0p1_tp  ! projected
! local
    real(r4)                                :: cos_dlamda1
    real(r4)                                :: dlamda1
    real(r4), dimension(1:3)                :: op2
    real(r4), dimension(1:3)                :: p0p2

    cos_dlamda1 = dot_product(op1,op0)/(norm(op0)*norm(op1))
    dlamda1     = acos(cos_dlamda1)              ! radian between op0 and op1
!
! another way to compute dlamda1 is arclen(op1,op0), 
! same to machine precision
!
    op2         = op1/cos_dlamda1
    p0p2        = op2-op0
    p0p1_tp     = p0p2*dlamda1/tan(dlamda1)

    return
  end subroutine project_sphere2tplane_r4

!================================================
!          Set coordinate in LOB
!================================================

  subroutine get_lob_xy(op0,p0p1,local_nx,local_ny)
! io
    real(r8), dimension(1:3), intent(in)  :: op0
    real(r8), dimension(1:3), intent(in)  :: p0p1
    real(r8), dimension(1:3), intent(out) :: local_nx
    real(r8), dimension(1:3), intent(out) :: local_ny
! local
    real(r8), dimension(1:3)  :: local_nx_m

    local_nx   = p0p1/norm(p0p1)   ! what we want here
    local_nx_m = -1._r8*local_nx
    local_ny   = cross_product(local_nx_m,op0)
    local_ny   = local_ny/norm(local_ny) ! what we want here

    return
  end subroutine get_lob_xy

!================================================
!         Set X, Y for polynomial, 2nd
!================================================

  subroutine set_xy_p2nd(p0p1,local_nx,local_ny,xy_list)
! io
    real(r8), dimension(1:3), intent(in)  :: p0p1
    real(r8), dimension(1:3), intent(in)  :: local_nx
    real(r8), dimension(1:3), intent(in)  :: local_ny
    real(r8)                , intent(out) :: xy_list(5)
! local
    real(r8)                              :: x
    real(r8)                              :: y

    x = dot_product(p0p1,local_nx)
    y = dot_product(p0p1,local_ny)

    xy_list(1) = x 
    xy_list(2) = y 
    xy_list(3) = x**2
    xy_list(4) = x*y
    xy_list(5) = y**2
  
    return 
  end subroutine set_xy_p2nd

!================================================
!         Set X, Y for polynomial, 4th
!================================================

  subroutine set_xy_p4th(p0p1,local_nx,local_ny,xy_list)
! io
    real(r8), dimension(1:3), intent(in)  :: p0p1
    real(r8), dimension(1:3), intent(in)  :: local_nx
    real(r8), dimension(1:3), intent(in)  :: local_ny
    real(r8)                , intent(out) :: xy_list(14)
! local
    real(r8)                              :: x
    real(r8)                              :: y

    x = dot_product(p0p1,local_nx)
    y = dot_product(p0p1,local_ny)
! 1st
    xy_list(1)  =  x
    xy_list(2)  =  y
! 2nd
    xy_list(3)  =  x**2
    xy_list(4)  =  x*y
    xy_list(5)  =  y**2
! 3rd
    xy_list(6)  =  x**3
    xy_list(7)  = (x**2)*y
    xy_list(8)  = (y**2)*x
    xy_list(9)  =  y**3
! 4th
    xy_list(10) =  x**4
    xy_list(11) = (x**3)*y
    xy_list(12) = (x**2)*(y**2)
    xy_list(13) = (y**3)*x
    xy_list(14) =  y**4

    return 
  end subroutine set_xy_p4th

  end module grist_mesh_weight_icosh
