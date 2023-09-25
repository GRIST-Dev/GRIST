
 !----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: mesh data for holding grid file data
!
! Revision history:
!----------------------------------------------------------------------------

 module grist_grid_file_vars 

   use grist_constants,         only: i4, r8, i8
   use grist_domain_types,      only: global_domain
   use grist_data_types,        only: scalar_2d_field
   use grist_element_types_icosh,only: nDual

   implicit none

     public

!--------------------------
!      Global info (27)
!--------------------------

         integer(i8)           :: mesh_nt
         integer(i8)           :: mesh_ne
         integer(i8)           :: mesh_nv
         character(len=16)     :: mesh_kind
         character(len=16)     :: mesh_node
         character(len=16)     :: mesh_optm
         integer(i8)           :: mesh_glevel
         integer(i8)           :: mesh_maxvnb
         real(r8)              :: mesh_min_edt_dist,  mesh_max_edt_dist,  mesh_mean_edt_dist
         real(r8)              :: mesh_min_edp_dist,  mesh_max_edp_dist,  mesh_mean_edp_dist
         real(r8)              :: mesh_min_tri_area,  mesh_max_tri_area,  mesh_mean_tri_area
         real(r8)              :: mesh_min_plg_area,  mesh_max_plg_area,  mesh_mean_plg_area
         real(r8)              :: mesh_min_tri_angle, mesh_max_tri_angle, mesh_mean_tri_angle
         character(len=32)     :: mesh_name

!--------------------------
!     triangle point(8)
!--------------------------

         type(scalar_2d_field) :: tri_v         ! generating vertex of tr,          nt*nDual
         type(scalar_2d_field) :: tri_v_local   ! generating vertex of tr,          nt*nDual
         type(scalar_2d_field) :: tri_cc_p      ! vector of cc (circumcenter), nt*3
         type(scalar_2d_field) :: tri_cc_ltln   ! lat and lon of cc,           nt*2
         type(scalar_2d_field) :: tri_nb        ! edge-neighboring tr of this tr,   nt*nDual
         type(scalar_2d_field) :: tri_ed        ! generating edge of tr,            nt*nDual
         type(scalar_2d_field) :: tri_area      ! area geodesic,               nt*1
#ifdef CUBE
         type(scalar_2d_field) :: tri_nnb       !                                   nt*1
#endif

!--------------------------
!     edge point(12)
!--------------------------

         type(scalar_2d_field) :: edt_cc_p      ! tr's edge point's vector,         ne*3
         type(scalar_2d_field) :: edt_cc_ltln   ! tr's edge point's lat&lon,        ne*2
         type(scalar_2d_field) :: edt_nr        ! edt's normal vector,              ne*3
         type(scalar_2d_field) :: edt_tg        ! edt's tangent vector,             ne*3
         type(scalar_2d_field) :: edt_v         ! index of two ends of edt,         ne*2
         type(scalar_2d_field) :: edt_v_local   ! index of two ends of edt,         ne*2
         type(scalar_2d_field) :: edt_len       ! [1] lenp [2] leng,                ne*2

         type(scalar_2d_field) :: edp_cc_p      ! plg's edge point's vector,        ne*3
         type(scalar_2d_field) :: edp_cc_ltln   ! plg's edge point's lat&lon,       ne*2
         type(scalar_2d_field) :: edp_nr        ! edp's normal vector,              ne*3
         type(scalar_2d_field) :: edp_tg        ! edp's tangent vector,             ne*3
         type(scalar_2d_field) :: edp_v         ! index of two ends of edp,         ne*2
         type(scalar_2d_field) :: edp_len       ! [1] lenp [2] leng,                ne*2

!--------------------------
! vertex/hexagon point(10)
!--------------------------

         type(scalar_2d_field) :: vtx_p         ! vertex's vector,                  nv*3
         type(scalar_2d_field) :: vtx_ltln_nnb  ! vertex's [1] lat [2] lon [3] nnb, nv*3
         type(scalar_2d_field) :: vtx_ltln_nnb_local  ! vertex's [1] lat [2] lon [3] nnb,nv*3
         type(scalar_2d_field) :: vtx_nb        ! index of each nb,                 nv*maxvnb
         type(scalar_2d_field) :: vtx_nb_local  ! index of each nb,                 nv*maxvnb
         type(scalar_2d_field) :: vtx_ed        ! neighboring ed,                   nv*maxvnb
         type(scalar_2d_field) :: vtx_ed_local  ! neighboring ed,                   nv*maxvnb
         type(scalar_2d_field) :: vtx_tr        ! neighboring tr,                   nv*maxvnb
         type(scalar_2d_field) :: vtx_tr_local  ! neighboring tr,                   nv*maxvnb

         type(scalar_2d_field) :: plg_area      ! area geodesic                     nv*1
#ifdef LAM_DOMAIN
         type(scalar_2d_field) :: plg_mark      ! mark of bdry cells for LAM        nv*1
         type(scalar_2d_field) :: edt_mark      ! mark of bdry edgess for LAM       ne*1
         type(scalar_2d_field) :: tri_mark      ! mark of bdry tris for LAM         nt*1
#endif

   contains

   subroutine construct_grid_vars_for_write(mesh)
! io
     type(global_domain),  intent(in) :: mesh
! local
     integer(i4)     :: it
     integer(i4)     :: ie
     integer(i4)     :: iv

!--------------------------
!      Global info
!--------------------------

         mesh_nt            = mesh%nt
         mesh_ne            = mesh%ne
         mesh_nv            = mesh%nv
  
         mesh_kind          = mesh%kind
         mesh_node          = mesh%node
         mesh_optm          = mesh%optm
         mesh_glevel        = mesh%glevel
         mesh_maxvnb        = mesh%maxvnb

         mesh_min_edt_dist  = mesh%min_edt_dist
         mesh_max_edt_dist  = mesh%max_edt_dist
         mesh_mean_edt_dist = mesh%mean_edt_dist

         mesh_min_edp_dist  = mesh%min_edp_dist
         mesh_max_edp_dist  = mesh%max_edp_dist
         mesh_mean_edp_dist = mesh%mean_edp_dist

         mesh_min_tri_area  = mesh%min_tri_area
         mesh_max_tri_area  = mesh%max_tri_area
         mesh_mean_tri_area = mesh%mean_tri_area

         mesh_min_plg_area  = mesh%min_plg_area
         mesh_max_plg_area  = mesh%max_plg_area
         mesh_mean_plg_area = mesh%mean_plg_area

         mesh_min_tri_angle = mesh%min_tri_angle
         mesh_max_tri_angle = mesh%max_tri_angle
         mesh_mean_tri_angle= mesh%mean_tri_angle

         mesh_name          = mesh%gridFileName

!--------------------------
!   Triangle point
!--------------------------

     allocate(      tri_v%f(mesh%nt,nDual))
     allocate(   tri_cc_p%f(mesh%nt,3))
     allocate(tri_cc_ltln%f(mesh%nt,2))
     allocate(     tri_nb%f(mesh%nt,nDual))
     allocate(     tri_ed%f(mesh%nt,nDual))
     allocate(   tri_area%f(mesh%nt,1))

     tri_v%pos        = 1
     tri_cc_p%pos     = 1
     tri_cc_ltln%pos  = 1
     tri_nb%pos       = 1
     tri_ed%pos       = 1
     tri_area%pos     = 1

     do it = 1, mesh%nt

        tri_v%f(it,1:nDual) = mesh%tri(it)%v(1:nDual)
        tri_cc_p%f(it,1:3)  = mesh%tri(it)%c%p(1:3)
        tri_cc_ltln%f(it,1) = mesh%tri(it)%c%lat
        tri_cc_ltln%f(it,2) = mesh%tri(it)%c%lon
        tri_nb%f(it,1:nDual)= mesh%tri(it)%nb(1:nDual)
        tri_ed%f(it,1:nDual)= mesh%tri(it)%ed(1:nDual)
        tri_area%f(it,1)    = mesh%tri(it)%areag

     end do

!--------------------------
!   Edge point
!--------------------------

     allocate(   edt_cc_p%f(mesh%ne,3))
     allocate(edt_cc_ltln%f(mesh%ne,2))
     allocate(     edt_nr%f(mesh%ne,3))
     allocate(     edt_tg%f(mesh%ne,3))
     allocate(      edt_v%f(mesh%ne,2))
     allocate(    edt_len%f(mesh%ne,1))

     allocate(   edp_cc_p%f(mesh%ne,3))
     allocate(edp_cc_ltln%f(mesh%ne,2))
     allocate(     edp_nr%f(mesh%ne,3))
     allocate(     edp_tg%f(mesh%ne,3))
     allocate(      edp_v%f(mesh%ne,2))
     allocate(    edp_len%f(mesh%ne,1))

     edt_cc_p%pos     = 6
     edt_cc_ltln%pos  = 6
     edt_v%pos        = 6
     edt_tg%pos       = 6
     edt_nr%pos       = 6
     edt_len%pos      = 6

     edp_cc_p%pos     = 6
     edp_cc_ltln%pos  = 6
     edp_v%pos        = 6
     edp_tg%pos       = 6
     edp_nr%pos       = 6
     edp_len%pos      = 6

     do ie = 1, mesh%ne
         
        edt_cc_p%f(ie,1:3)   = mesh%edt(ie)%c%p(1:3)
        edt_cc_ltln%f(ie,1)  = mesh%edt(ie)%c%lat
        edt_cc_ltln%f(ie,2)  = mesh%edt(ie)%c%lon
        edt_nr%f(ie,1:3)     = mesh%edt(ie)%nr(1:3)
        edt_tg%f(ie,1:3)     = mesh%edt(ie)%tg(1:3)
        edt_v%f(ie,1:2)      = mesh%edt(ie)%v(1:2)
        edt_len%f(ie,1)      = mesh%edt(ie)%leng

        edp_cc_p%f(ie,1:3)   = mesh%edp(ie)%c%p(1:3)
        edp_cc_ltln%f(ie,1)  = mesh%edp(ie)%c%lat
        edp_cc_ltln%f(ie,2)  = mesh%edp(ie)%c%lon
        edp_nr%f(ie,1:3)     = mesh%edp(ie)%nr(1:3)
        edp_tg%f(ie,1:3)     = mesh%edp(ie)%tg(1:3)
        edp_v%f(ie,1:2)      = mesh%edp(ie)%v(1:2)
        edp_len%f(ie,1)      = mesh%edp(ie)%leng

     end do

!--------------------------
!  Vertex/Hexagon point
!--------------------------

     allocate(       vtx_p%f(mesh%nv,3))
     allocate(vtx_ltln_nnb%f(mesh%nv,3))
     allocate(      vtx_nb%f(mesh%nv,mesh_maxvnb))
     allocate(      vtx_ed%f(mesh%nv,mesh_maxvnb)) 
     allocate(      vtx_tr%f(mesh%nv,mesh_maxvnb)) 
          
     allocate(plg_area%f(mesh%nv,1))

     vtx_p%pos           = 0
     vtx_ltln_nnb%pos    = 0
     vtx_nb%pos          = 0
     vtx_ed%pos          = 0
     vtx_tr%pos          = 0

     plg_area%pos        = 0

     do iv = 1, mesh%nv

        vtx_p%f(iv,1:3)                      = mesh%vtx(iv)%p(1:3)
        vtx_ltln_nnb%f(iv,1)                 = mesh%vtx(iv)%lat
        vtx_ltln_nnb%f(iv,2)                 = mesh%vtx(iv)%lon
        vtx_ltln_nnb%f(iv,3)                 = mesh%vtx(iv)%nnb
        vtx_nb%f(iv,  1:mesh%vtx(iv)%nnb)    = mesh%vtx(iv)%nb(:)
        vtx_ed%f(iv,  1:mesh%vtx(iv)%nnb)    = mesh%vtx(iv)%ed(:)
        vtx_tr%f(iv,  1:mesh%vtx(iv)%nnb)    = mesh%vtx(iv)%tr(:)

        plg_area%f(iv,1)                     = mesh%plg(iv)%areag

        if(mesh%vtx(iv)%nnb.eq.5.and.mesh%glevel.gt.0)then
           vtx_nb%f(iv,6)  = -1e6_r8
           vtx_ed%f(iv,6)  = -1e6_r8
           vtx_tr%f(iv,6)  = -1e6_r8
        end if

     end do

     return
   end subroutine construct_grid_vars_for_write

   subroutine construct_grid_vars_for_read_local_t(max_cache)
     integer(i4),intent(in)     :: max_cache

!--------------------------
!   Triangle point
!--------------------------

     allocate(tri_v_local%f(max_cache,nDual))
     allocate(   tri_cc_p%f(max_cache,3))
     allocate(tri_cc_ltln%f(max_cache,2))
     allocate(     tri_nb%f(max_cache,nDual))
     allocate(     tri_ed%f(max_cache,nDual))
     allocate(   tri_area%f(max_cache,1))
#ifdef CUBE
     allocate(    tri_nnb%f(max_cache,1))
#endif
#ifdef LAM_DOMAIN
     allocate(   tri_mark%f(max_cache,1))
#endif

     tri_v_local%pos  = 1
     tri_cc_p%pos     = 1
     tri_cc_ltln%pos  = 1
     tri_nb%pos       = 1
     tri_ed%pos       = 1
     tri_area%pos     = 1
#ifdef CUBE
     tri_nnb%pos      = 1
#endif
#ifdef LAM_DOMAIN
     tri_mark%pos     = 1
#endif

     return
   end subroutine construct_grid_vars_for_read_local_t

   subroutine construct_grid_vars_for_read_local_e(max_cache)
     integer(i4)     :: max_cache

!--------------------------
!   Edge point
!--------------------------

     allocate(   edt_cc_p%f(max_cache,3))
     allocate(edt_cc_ltln%f(max_cache,2))
     allocate(     edt_nr%f(max_cache,3))
     allocate(     edt_tg%f(max_cache,3))
     allocate(edt_v_local%f(max_cache,2))  !repeat to global,so named edt_v_local
     allocate(    edt_len%f(max_cache,1))
#ifdef LAM_DOMAIN
     allocate(    edt_mark%f(max_cache,1))
#endif

     allocate(   edp_cc_p%f(max_cache,3))
     allocate(edp_cc_ltln%f(max_cache,2))
     allocate(     edp_nr%f(max_cache,3))
     allocate(     edp_tg%f(max_cache,3))
     allocate(      edp_v%f(max_cache,2))
     allocate(    edp_len%f(max_cache,1))

     edt_cc_p%pos     = 6
     edt_cc_ltln%pos  = 6
     edt_nr%pos       = 6
     edt_tg%pos       = 6
     edt_v_local%pos  = 6
     edt_len%pos      = 6

     edp_cc_p%pos     = 6
     edp_cc_ltln%pos  = 6
     edp_nr%pos       = 6
     edp_tg%pos       = 6
     edp_v%pos        = 6
     edp_len%pos      = 6
#ifdef LAM_DOMAIN
     edt_mark%pos     = 6
#endif

     return
   end subroutine construct_grid_vars_for_read_local_e

   subroutine construct_grid_vars_for_read_local_v(max_cache)
     integer(i4)     :: max_cache

!--------------------------
!  Vertex/Hexagon point
!--------------------------

     allocate(       vtx_p%f(max_cache,3))
     allocate(vtx_ltln_nnb_local%f(max_cache,3))
     allocate(      vtx_nb_local%f(max_cache,mesh_maxvnb))
     allocate(      vtx_ed_local%f(max_cache,mesh_maxvnb))
     allocate(      vtx_tr_local%f(max_cache,mesh_maxvnb))
          
     allocate(   plg_area%f(max_cache,1))
#ifdef LAM_DOMAIN
     allocate(    plg_mark%f(max_cache,1))
#endif

     vtx_p%pos           = 0
     vtx_ltln_nnb_local%pos    = 0
     vtx_nb_local%pos    = 0
     vtx_ed_local%pos    = 0
     vtx_tr_local%pos    = 0

     plg_area%pos        = 0
#ifdef LAM_DOMAIN
     plg_mark%pos        = 0
#endif

     return
   end subroutine construct_grid_vars_for_read_local_v

   subroutine construct_grid_vars_for_read

!--------------------------
!   Triangle point
!--------------------------


     allocate(      tri_v%f(mesh_nt,3))

!--------------------------
!   Edge point
!--------------------------

     allocate(      edt_v%f(mesh_ne,2))

!--------------------------
!  Vertex/Hexagon point
!--------------------------

     allocate(vtx_ltln_nnb%f(mesh_nv,3))
     allocate(      vtx_nb%f(mesh_nv,mesh_maxvnb))
     allocate(      vtx_ed%f(mesh_nv,mesh_maxvnb))
     allocate(      vtx_tr%f(mesh_nv,mesh_maxvnb))
          
     return
   end subroutine construct_grid_vars_for_read

   subroutine destruct_grid_file_vars_t

!--------------------------
!   Triangle point
!--------------------------

     deallocate(tri_v_local%f)
     deallocate(   tri_cc_p%f)
     deallocate(tri_cc_ltln%f)
     deallocate(     tri_nb%f)
     deallocate(     tri_ed%f)
     deallocate(   tri_area%f)
#ifdef CUBE
     deallocate(   tri_nnb%f)
#endif
#ifdef LAM_DOMAIN
     deallocate(   tri_mark%f)
#endif

   end subroutine destruct_grid_file_vars_t

   subroutine destruct_grid_file_vars_e

!--------------------------
!   Edge point
!--------------------------

     deallocate(   edt_cc_p%f)
     deallocate(edt_cc_ltln%f)
     deallocate(     edt_nr%f)
     deallocate(     edt_tg%f)
     deallocate(edt_v_local%f)
     deallocate(    edt_len%f)
#ifdef LAM_DOMAIN
     deallocate(   edt_mark%f)
#endif

     deallocate(   edp_cc_p%f)
     deallocate(edp_cc_ltln%f)
     deallocate(     edp_nr%f)
     deallocate(     edp_tg%f)
     deallocate(      edp_v%f)
     deallocate(    edp_len%f)

   end subroutine destruct_grid_file_vars_e

   subroutine destruct_grid_file_vars_v

!--------------------------
!  Vertex/Hexagon point
!--------------------------

     deallocate(       vtx_p%f)
     deallocate(vtx_ltln_nnb_local%f)
     deallocate(      vtx_nb_local%f)
     deallocate(      vtx_ed_local%f)
     deallocate(      vtx_tr_local%f)
          
     deallocate(   plg_area%f)
#ifdef LAM_DOMAIN
     deallocate(   plg_mark%f)
#endif

   end subroutine destruct_grid_file_vars_v

   subroutine destruct_grid_file_vars

!--------------------------
!   Triangle point
!--------------------------

     deallocate(      tri_v%f)
     deallocate(   tri_cc_p%f)
     deallocate(tri_cc_ltln%f)
     deallocate(     tri_nb%f)
     deallocate(     tri_ed%f)
     deallocate(   tri_area%f)

!--------------------------
!   Edge point
!--------------------------

     deallocate(   edt_cc_p%f)
     deallocate(edt_cc_ltln%f)
     deallocate(     edt_nr%f)
     deallocate(     edt_tg%f)
     deallocate(      edt_v%f)
     deallocate(    edt_len%f)

     deallocate(   edp_cc_p%f)
     deallocate(edp_cc_ltln%f)
     deallocate(     edp_nr%f)
     deallocate(     edp_tg%f)
     deallocate(      edp_v%f)
     deallocate(    edp_len%f)

!--------------------------
!  Vertex/Hexagon point
!--------------------------

     deallocate(       vtx_p%f)
     deallocate(vtx_ltln_nnb%f)
     deallocate(      vtx_nb%f)
     deallocate(      vtx_ed%f)
     deallocate(      vtx_tr%f)
          
     deallocate(   plg_area%f)

   end subroutine destruct_grid_file_vars

 end module grist_grid_file_vars
