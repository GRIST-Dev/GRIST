
!-----------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Major model domain for computing
! Revision history:
!-----------------------------------------------------------------------

module grist_domain_types

   use grist_lib
   use grist_constants,     only: i4, r8, i8
   !use grist_constants_dbl, only: i4, r8, i8
   use grist_element_types_icosh,only: node_structure  , &
                                  edge_structure  , &
                                  dual_structure  , &
                                  prim_structure  , &
                                  node_structure_global  , &
                                  edge_structure_global  , &
                                  dual_structure_global
   use grist_data_types,    only: scalar_2d_field

   implicit none

   public :: group_comm         ,&
             global_domain       &
#ifndef SEQ_GRIST
            ,block_structure     &
            ,global_domain_data  &
#endif
            ,grist_domain_check_runtime_allocated&
            ,grist_domain_clean_runtime_allocated


  type group_comm
     ! group communicator
     integer                  :: color, gcomm, grank, gsize, all_comm

     ! All to All
     integer(i4), allocatable :: displs_sendv(:), counts_sendv(:)
     integer(i4), allocatable :: displs_recvv(:), counts_recvv(:)
     integer(i4), allocatable :: displs_sende(:), counts_sende(:)
     integer(i4), allocatable :: displs_recve(:), counts_recve(:)
     integer(i4), allocatable :: displs_sendt(:), counts_sendt(:)
     integer(i4), allocatable :: displs_recvt(:), counts_recvt(:)

     integer(i4), allocatable :: idx_sendv(:), idx_sende(:), idx_sendt(:)
     integer(i4), allocatable :: idx_recvv(:), idx_recve(:), idx_recvt(:)
     integer(i4), allocatable :: idx_localv(:), idx_locale(:), idx_localt(:)
     integer(i4), allocatable :: mystart_v, mycount_v, mystart_e, mycount_e, mystart_t, mycount_t
  end type group_comm
  
  type global_domain

     type(node_structure), allocatable :: vtx(:)   ! list of vertice/node,  size 1:nv
     type(edge_structure), allocatable :: edt(:)   ! list of tr edge,       size 1:ne
     type(edge_structure), allocatable :: edp(:)   ! list of polygon edge,  size 1:ne
     type(dual_structure), allocatable :: tri(:)   ! list of triangle,      size 1:nt
     type(prim_structure), allocatable :: plg(:)   ! list of polygon,       size 1:nv
     type(group_comm)                  :: gcomm_read, gcomm_write !used for group read and write

     real(r8)    :: min_edt_dist, max_edt_dist,  mean_edt_dist   ! min/max/mean angular distance between triangle's vertices (edt_v) in radian
     real(r8)    :: min_edp_dist, max_edp_dist,  mean_edp_dist   ! min/max/mean angular distance between triangle's (circum)centers (edp_v) in radian
     real(r8)    :: min_tri_area, max_tri_area,  mean_tri_area   ! min/max/mean triangle geodesic area
     real(r8)    :: min_tri_angle,max_tri_angle, mean_tri_angle  ! min/max/mean triangle geodesic internal angle
     real(r8)    :: min_plg_area, max_plg_area,  mean_plg_area   ! min/max/mean voronoi cell geodesic areas

     integer(i4) :: nt              ! number of triangle
     integer(i4) :: ne              ! number of edges
     integer(i4) :: nv              ! number of vertice/polygon
     integer(i4) :: maxvnb          ! maximum number of neighbours for all nodes (triangle vertices)
     integer(i4) :: glevel          ! number of bisection (g-level)
     integer(i4) :: directly        ! used by gg, 0: directly glevel, 1: one by one
     integer(i4) :: num_cc_out_of_tr

     character(len=16)  :: kind          ! mesh kind: icos, rand
     character(len=16)  :: node          ! mesh position: eqs, pol, ran, ref
     character(len=16)  :: optm          ! optimization method: nopt, ccvt
     character(len=32)  :: gridFileName  ! grid name, used for file names as outputs
     character(len=128) :: meshFileName  ! file name if nodes are read from file

     !====add to par====
     !> globally unique index for tr-based values 
     integer(i4),              allocatable  :: t_index(:)

     !> globally unique index for ed-based values
     integer(i4),              allocatable  :: e_index(:)

     !> globally unique index for v-based values 
     integer(i4),              allocatable  :: v_index(:)   
#ifndef SEQ_GRIST
     !> linked list when searching
     type(set) :: t_index_list 
     type(set) :: e_index_list
     type(set) :: v_index_list

     !> given any sub-domain, hash table for connecting 
     !> key="global index" and val="local domain index"
     !> or use other alternative
     type(map) :: map_g2l_t   ! t index
     type(map) :: map_g2l_e   ! e index
     type(map) :: map_g2l_v   ! v index
#endif
     integer(i4), allocatable  :: index_in_full
 
     !> a pointer to local block
#ifndef SEQ_GRIST
     type(block_structure) ,pointer   :: local_block
#endif
 
     integer(i4) :: nt_compute          !> number of triangle/corner
     integer(i4) :: ne_compute          !> number of edges
     integer(i4) :: nv_compute          !> number of cell/polygon

     integer(i4) :: nt_full             !> number of triangle/corner
     integer(i4) :: ne_full             !> number of edges
     integer(i4) :: nv_full             !> number of cell/polygon
 
     integer(i4) :: nt_iner             !> number of triangle/corner
     integer(i4) :: ne_iner             !> number of edges
     integer(i4) :: nv_iner             !> number of cell/polygon
 
     integer(i4) :: nt_all              !> all number of triangle
     integer(i4) :: ne_all              !> all number of edges
     integer(i4) :: nv_all              !> all number of vertice/polygon
 
     integer(i4), allocatable :: nt_halo(:)              !> number of triangle
     integer(i4), allocatable :: ne_halo(:)              !> number of edges
     integer(i4), allocatable :: nv_halo(:)              !> number of vertice/polygon

     integer(i4), allocatable :: nt_bdry(:)              !> number of triangle
     integer(i4), allocatable :: ne_bdry(:)              !> number of edges
     integer(i4), allocatable :: nv_bdry(:)              !> number of vertice/polygon
 
    !> the decomposed partition flag
     integer(i4), allocatable :: parts(:)
     
!
! compute pattern 1, added for runtime
! frequently-used computational pattern (pattern1: 34 currently
! we should use all this for runtime! normally 34-4-2=28
!

     !> for calculating the coriolis term
     !real(r8),    allocatable :: edp_trsk_on_edge(:,:)
     !real(r8),    allocatable :: edp_edpl_on_edge(:,:)
     real(r8),    allocatable :: edp_trsk_edpl_on_edge(:,:) ! trsk weight*edpl
     integer(i4), allocatable :: edp_edge_on_edge(:,:)
     integer(i4), allocatable :: edp_nedge(:)
     !> for calculating hx 2nd order derivative
     real(r8),    allocatable :: edt_v_weight_prime_cell(:,:,:)
     integer(i4), allocatable :: plg_stencil_number_2nd(:)
     integer(i4), allocatable :: plg_stencil_index_2nd(:,:)
     !> for calculating hx 4th order derivative
     real(r8),    allocatable :: edt_v_weight_4th_prime_cell(:,:,:)      ! only for 4th-der, DUSE_HALO2 wlll ignore
     real(r8),    allocatable :: edt_v_weight_4th_2nd_prime_cell(:,:,:)  ! ...
     integer(i4), allocatable :: plg_stencil_number_4th(:)               ! ...
     integer(i4), allocatable :: plg_stencil_index_4th(:,:)              ! ...
     !> for calculating tr 2nd order derivative
     real(r8),    allocatable :: edp_v_weight_dual_cell(:,:,:) ! for 3d atmos, this is only for pv3, pv2 will ignore
     integer(i4), allocatable :: tri_stencil_number_2nd(:)     ! ...
     integer(i4), allocatable :: tri_stencil_index_2nd(:,:)    ! ...
     !> for calculating digvergence
     integer(i4), allocatable :: vtx_nnb(:)
     integer(i4), allocatable :: vtx_ed(:,:)
     integer(i4), allocatable :: plg_nr(:,:)
     real(r8),    allocatable :: edp_leng(:)
     real(r8),    allocatable :: plg_areag(:)
     !> for calculating prim-cell flux
     integer(i4), allocatable :: edt_v(:,:)
     integer(i4), allocatable :: edt_edpNr_edtTg(:)
     real(r8),    allocatable :: edt_leng(:)
     integer(i4), allocatable :: edt_my_edge_on_edge(:,:,:)  ! only for O58, DUSE_HALO2 will ignore
     integer(i4), allocatable :: edt_ur_cell_on_edge(:,:,:)  ! only for O58, DUSE_HALO2 will ignore
     !> for calculating dual-cell flux
     integer(i4), allocatable :: edp_v(:,:)
     integer(i4), allocatable :: edp12_edp_tg(:)
     integer(i4), allocatable :: edp12_edt_nr(:)
     !> for calculating vorticity(curl)
     integer(i4), allocatable :: tri_nnb(:)
     integer(i4), allocatable :: tri_ed(:,:)
     integer(i4), allocatable :: tri_nr(:,:)
     integer(i4), allocatable :: edt_edpTg_edtNr(:)
     real(r8),    allocatable :: tri_areag(:)
     real(r8),    allocatable :: tri_c_lat(:)
     !> for calculating kinetic_energy
     integer(i4), allocatable :: vtx_tr(:,:)
     real(r8),    allocatable :: plg_kite_area(:,:)
     !> added for VR
     real(r8),    allocatable :: vtxCellLeng(:)
     real(r8),    allocatable :: edt_scale_2nd(:)
     real(r8),    allocatable :: edt_scale_4th(:)
     real(r8),    allocatable :: edt_scale_6th(:)
     real(r8)                 :: minCellLeng
!
! pattern2 (20 in total, normal 15, following pattern1, to completely eliminate runtime mesh%${five_element}
! so as to improve efficiency
     real(r8),    allocatable :: vtx_lon(:)
     real(r8),    allocatable :: vtx_lat(:)
     real(r8),    allocatable :: vtx_p(:,:)
     real(r8),    allocatable :: vtx_nb(:,:)
     real(r8),    allocatable :: vtx_lob_nx(:,:,:)   ! used by ffsl tracer transport, non-ffsl will ignore
     real(r8),    allocatable :: vtx_lob_ny(:,:,:)   ! ffsl: see above

     real(r8),    allocatable :: plg_blocal(:,:,:)                   ! ffsl: m*n*stencil_number X nv
     real(r8),    allocatable :: plg_sub_triangle_midp_lob(:,:,:,:)  ! ffsl: nnb*3edges*2(lob_x,lob_y) X nplg
     real(r8),    allocatable :: plg_sub_triangle_area(:,:)          ! ffsl: nnb X nv

     real(r8),    allocatable :: tri_c_p(:,:)       ! 3, nv
     real(r8),    allocatable :: tri_c_lon(:)       ! nv
     real(r8),    allocatable :: tri_v(:,:)         ! nb=3, nv
     real(r8),    allocatable :: tri_kite_area(:,:) ! nb=3, nv

     real(r8),    allocatable :: edt_c_p(:,:)  ! 3*ne
     real(r8),    allocatable :: edt_c_lon(:)  ! ne
     real(r8),    allocatable :: edt_c_lat(:)  ! ne
     real(r8),    allocatable :: edt_nr(:,:)   ! 3*ne

     real(r8),    allocatable :: edp_c_p(:,:)  ! 3*ne
     real(r8),    allocatable :: edp_nr(:,:)   ! 3*ne
     real(r8),    allocatable :: edp_tg(:,:)   ! 3*ne

#ifdef LAM_DOMAIN
     integer(i4), allocatable :: plg_mark(:), edt_mark(:), tri_mark(:) ! comp: 199; relx: 188; spec: 177
     logical                  :: hasBdyDomain
#endif
  end type global_domain

#ifndef SEQ_GRIST  
  type comm_sites
     integer :: cpuid
     ! integer, allocatable :: sites(:) ! global index
     ! integer, allocatable :: local_sites(:) !local index
     integer(i4), allocatable :: v(:), e(:), t(:)
     integer(i4), allocatable :: local_v(:), local_e(:), local_t(:)
     integer(i4)              :: nv, nt, ne
  end type comm_sites

  type output_structure
     integer(i4), allocatable :: v(:)
     integer(i4), allocatable :: e(:)
     integer(i4), allocatable :: t(:)
  end type output_structure

  type block_structure

     !> communicator
     integer :: comm

     !> process id of this block
     integer(i4)               :: cpuid

     !> number of processes
     integer                   :: nprocs

     !> number of neighboring blocks
     integer(i4)               :: nnb

     !> cpuid  of neighboring blocks
     integer(i4), allocatable  :: nbcpuid(:)

     !> cpuid of neighboring blocks
     type(set)                 :: nb_list

     !> sub-domain info, no need to use them all
     !> halo1+halo2+compute, data array/mesh elements conformed 
     type(global_domain)          :: full_domain

     !> contains elements needed to be updated, belong to this block
     type(global_domain)          :: compute_domain

     !> iner domain
     type(global_domain)          :: iner_domain

     !> bdy1 + bdy2
     type(global_domain)          :: bdry_domain
     
     !> halo and boundary regions
     type(global_domain), allocatable  :: bdry(:)
     type(global_domain), allocatable  :: halo(:)

     type(comm_sites), allocatable  :: send_sites(:)
     type(comm_sites), allocatable  :: recv_sites(:)

     integer(i4)  :: max_nv_send, max_nv_recv
     integer(i4)  :: max_nt_send, max_nt_recv
     integer(i4)  :: max_ne_send, max_ne_recv

     type(output_structure) :: output
  end type block_structure

  type global_domain_data
! global data read on an entire mesh
     type(node_structure_global),  allocatable :: vtx(:)   ! list of vertice/node,  size 1:nv
     type(edge_structure_global),  allocatable :: edt(:)   ! list of tr edge,       size 1:ne
     type(dual_structure_global),  allocatable :: tri(:)   ! list of triangle,      size 1:nt
     
     real(r8)    :: min_edt_dist, max_edt_dist,  mean_edt_dist   ! min/max/mean angular distance between tri vertices in radian
     real(r8)    :: min_edp_dist, max_edp_dist,  mean_edp_dist   ! min/max/mean angular distance between tri circumcenters in radian
     real(r8)    :: min_tri_area, max_tri_area,  mean_tri_area   ! min/max/mean triangle geodesic area
     real(r8)    :: min_tri_angle,max_tri_angle, mean_tri_angle  ! min/max/mean triangle geodesic internal angle
     real(r8)    :: min_plg_area, max_plg_area,  mean_plg_area   ! min/max/mean voronoi cell geodesic areas
     integer(i4) :: nt              ! number of triangle
     integer(i4) :: ne              ! number of edges
     integer(i4) :: nv              ! number of vertice/polygon
     integer(i4) :: maxvnb          ! maximum number of neighbours for all nodes (triangle vertices)
     integer(i4) :: glevel          ! number of bisection (g-level)
     integer(i4) :: loadable        ! flag for loadable grid (1- yes, 0- no)
     integer(i4) :: directly        ! used by gg, 0: directly glevel, 1: one by one
     integer(i4) :: num_cc_out_of_tr

     character(len=16)  :: kind     ! mesh kind: icos, rand
     character(len=16)  :: node     ! mesh node position: eqs, pol, ran, ref
     character(len=16)  :: optm     ! optimization method: nopt, ccvt
     character(len=32)  :: gridFileName     ! grist grid name, used for file names as outputs
     character(len=128) :: meshFileName     ! read  file name, used by gg
     integer(i4), allocatable :: parts(:)

     type(scalar_2d_field)    :: tri_v, edt_v, vtx_nb, vtx_ed, vtx_tr, tri_nnb
#ifdef LAM_DOMAIN
     type(scalar_2d_field)    :: plg_mark, edt_mark, tri_mark
#endif
     integer(i4), allocatable :: vtx_nnb(:)
     real(r8),    allocatable :: vtx_ltln(:,:)

  end type global_domain_data
#endif

   contains

   subroutine grist_domain_check_runtime_allocated(mesh)
    use grist_lib
    use grist_nml_module, only: isolate_tracer_test
     type(global_domain), intent(inout)  :: mesh
     integer(i4) :: num_element

! remove unsued runtime based on check

     if(allocated(mesh%edt_c_lat).and. .not.isolate_tracer_test)  deallocate(mesh%edt_c_lat)
     if(allocated(mesh%edt_c_lon).and. .not.isolate_tracer_test)  deallocate(mesh%edt_c_lon)
     if(allocated(mesh%edt_nr))     deallocate(mesh%edt_nr)
     if(allocated(mesh%edp_c_p))    deallocate(mesh%edp_c_p)

     if(mpi_rank().eq.0)then
        num_element = 1
     !> for calculating the coriolis term
        !if(allocated(mesh%edp_trsk_on_edge)) call write_element_info("edp_trsk_on_edge is ON, num is ", num_element)
        !if(allocated(mesh%edp_edpl_on_edge)) call write_element_info("edp_edpl_on_edge is ON, num is ", num_element)
        if(allocated(mesh%edp_trsk_edpl_on_edge)) call write_element_info("edp_trsk_edpl_on_edge is ON, num is ", num_element)
        if(allocated(mesh%edp_edge_on_edge))      call write_element_info("edp_edge_on_edge is ON, num is ", num_element)
        if(allocated(mesh%edp_nedge))             call write_element_info("edp_nedge        is ON, num is ", num_element)

        !> for calculating hx 2nd order derivative
        if(allocated(mesh%edt_v_weight_prime_cell)) call write_element_info("edt_v_weight_prime_cell is ON, num is ", num_element)
        if(allocated(mesh%plg_stencil_number_2nd))  call write_element_info("plg_stencil_number_2nd  is ON, num is ", num_element)
        if(allocated(mesh%plg_stencil_index_2nd))   call write_element_info("plg_stencil_index_2nd   is ON, num is ", num_element)
        !> for calculating hx 4th order derivative
        if(allocated(mesh%edt_v_weight_4th_prime_cell))     call write_element_info("edt_v_weight_4th_prime_cell     is ON, num is ", num_element)
        if(allocated(mesh%edt_v_weight_4th_2nd_prime_cell)) call write_element_info("edt_v_weight_4th_2nd_prime_cell is ON, num is ", num_element)
        if(allocated(mesh%plg_stencil_number_4th))          call write_element_info("plg_stencil_number_4th is ON, num is ", num_element)
        if(allocated(mesh%plg_stencil_index_4th))           call write_element_info("plg_stencil_index_4th  is ON, num is ", num_element)
        !> for calculating tr 2nd order derivative
        if(allocated(mesh%edp_v_weight_dual_cell)) call write_element_info("edp_v_weight_dual_cell is ON, num is ", num_element)
        if(allocated(mesh%tri_stencil_number_2nd)) call write_element_info("tri_stencil_number_2nd is ON, num is ", num_element)
        if(allocated(mesh%tri_stencil_index_2nd))  call write_element_info("tri_stencil_index_2nd  is ON, num is ", num_element)
        !> for calculating digvergence
        if(allocated(mesh%vtx_nnb))   call write_element_info("vtx_nnb   is ON, num is ", num_element )
        if(allocated(mesh%vtx_ed))    call write_element_info("vtx_ed    is ON, num is ", num_element )
        if(allocated(mesh%plg_nr))    call write_element_info("plg_nr    is ON, num is ", num_element )
        if(allocated(mesh%edp_leng))  call write_element_info("edp_leng  is ON, num is ", num_element )
        if(allocated(mesh%plg_areag)) call write_element_info("plg_areag is ON, num is ", num_element )
        !> for calculating prim-cell flux
        if(allocated(mesh%edt_v))               call write_element_info("edt_v         is ON, num is ", num_element )
        if(allocated(mesh%edt_edpNr_edtTg))     call write_element_info("edt_edpNr_edtTg is ON, num is ", num_element )
        if(allocated(mesh%edt_leng))            call write_element_info("edt_leng      is ON, num is ", num_element )
        if(allocated(mesh%edt_my_edge_on_edge)) call write_element_info("edt_my_edge_on_edge is ON, num is ", num_element )
        if(allocated(mesh%edt_ur_cell_on_edge)) call write_element_info("edt_ur_cell_on_edge is ON, num is ", num_element )
        !> for calculating dual-cell flux
        if(allocated(mesh%edp_v))         call write_element_info("edp_v        is ON, num is ", num_element )
        if(allocated(mesh%edp12_edp_tg))  call write_element_info("edp12_edp_tg is ON, num is ", num_element )
        if(allocated(mesh%edp12_edt_nr))  call write_element_info("edp12_edt_nr is ON, num is ", num_element )
        !> for calculating vorticity(curl)
        if(allocated(mesh%tri_nnb))       call write_element_info("tri_nnb is ON, num is ",num_element )
        if(allocated(mesh%tri_ed))        call write_element_info("tri_ed is ON, num is ", num_element )
        if(allocated(mesh%tri_nr))        call write_element_info("tri_nr is ON, num is ", num_element )
        if(allocated(mesh%edt_edpTg_edtNr)) call write_element_info("edt_edpTg_edtNr is ON, num is ", num_element )
        if(allocated(mesh%tri_areag))     call write_element_info("tri_areag     is ON, num is ", num_element )
        if(allocated(mesh%tri_c_lat))     call write_element_info("tri_c_lat     is ON, num is ", num_element )
        !> for calculating kinetic_energy
        if(allocated(mesh%vtx_tr))        call write_element_info("vtx_tr        is ON, num is ", num_element)
        if(allocated(mesh%plg_kite_area)) call write_element_info("plg_kite_area is ON, num is ", num_element)
        !> added for VR
        if(allocated(mesh%vtxCellLeng))   call write_element_info("vtxCellLeng   is ON, num is ", num_element)
        if(allocated(mesh%edt_scale_2nd)) call write_element_info("edt_scale_2nd is ON, num is ", num_element)
        if(allocated(mesh%edt_scale_4th)) call write_element_info("edt_scale_4th is ON, num is ", num_element)
        if(allocated(mesh%edt_scale_6th)) call write_element_info("edt_scale_6th is ON, num is ", num_element)
        !real(r8)                 :: minCellLeng
!
! pattern2 (20 in total, normal 15, following pattern1, to completely eliminate runtime mesh%${five_element}
! so as to improve efficiency
!
        if(allocated(mesh%vtx_lon))    call write_element_info("vtx_lon    is ON, num is ", num_element )         
        if(allocated(mesh%vtx_lat))    call write_element_info("vtx_lat    is ON, num is ", num_element )      
        if(allocated(mesh%vtx_p))      call write_element_info("vtx_p      is ON, num is ", num_element )
        if(allocated(mesh%vtx_nb))     call write_element_info("vtx_nb     is ON, num is ", num_element )
        if(allocated(mesh%vtx_lob_nx)) call write_element_info("vtx_lob_nx is ON, num is ", num_element )
        if(allocated(mesh%vtx_lob_ny)) call write_element_info("vtx_lob_ny is ON, num is ", num_element )

        if(allocated(mesh%plg_blocal))                call write_element_info("plg_blocal is ON, num is ", num_element )
        if(allocated(mesh%plg_sub_triangle_midp_lob)) call write_element_info("plg_sub_triangle_midp_lob is ON, num is ", num_element )
        if(allocated(mesh%plg_sub_triangle_area))     call write_element_info("plg_sub_triangle_area is ON, num is ", num_element )

        if(allocated(mesh%tri_c_p))       call write_element_info("tri_c_p       is ON, num is ", num_element )
        if(allocated(mesh%tri_c_lon))     call write_element_info("tri_c_lon     is ON, num is ", num_element )
        if(allocated(mesh%tri_v))         call write_element_info("tri_v         is ON, num is ", num_element )
        if(allocated(mesh%tri_kite_area)) call write_element_info("tri_kite_area is ON, num is ", num_element )

        if(allocated(mesh%edt_c_p))       call write_element_info("edt_c_p   is ON, num is ", num_element )
        if(allocated(mesh%edt_c_lon))     call write_element_info("edt_c_lon is ON, num is ", num_element )
        if(allocated(mesh%edt_c_lat))     call write_element_info("edt_c_lat is ON, num is ", num_element )
        if(allocated(mesh%edt_nr))        call write_element_info("edt_nr    is ON, num is ", num_element )

        if(allocated(mesh%edp_c_p))       call write_element_info("edp_c_p is ON, num is ", num_element )
        if(allocated(mesh%edp_nr))        call write_element_info("edp_nr  is ON, num is ", num_element )
        if(allocated(mesh%edp_tg))        call write_element_info("edp_tg  is ON, num is ", num_element )

        if(allocated(mesh%vtx))then
           print*, "vtx is ON, We gonner deallocate it"
           deallocate(mesh%vtx)
        end if
        if(allocated(mesh%tri))then
           print*, "tri is ON, We gonner deallocate it"
           deallocate(mesh%tri)
        end if
        if(allocated(mesh%edt))then
           print*, "edt is ON, We gonner deallocate it"
           deallocate(mesh%edt)
        end if
        if(allocated(mesh%edp))then
           print*, "edp is ON, We gonner deallocate it"
           deallocate(mesh%edp)
        end if
        if(allocated(mesh%plg))then
           print*, "plg is ON, We gonner deallocate it"
           deallocate(mesh%plg)
        end if

      end if

   end subroutine grist_domain_check_runtime_allocated

   subroutine write_element_info(information,num_element)
! io
    character(len=*), intent(in)    :: information
    integer(i4),        intent(inout) :: num_element

    print*,trim(information), num_element
    num_element=num_element + 1
    
    return
   end subroutine write_element_info  

   subroutine grist_domain_clean_runtime_allocated(mesh)
    use grist_lib
     type(global_domain), intent(inout)  :: mesh
     integer(i4) :: num_element

! remove unsued runtime based on check

     !> for calculating the coriolis term
        !if(allocated(mesh%edp_trsk_on_edge))          deallocate(mesh%edp_trsk_on_edge)
        !if(allocated(mesh%edp_edpl_on_edge))          deallocate(mesh%edp_edpl_on_edge)
        if(allocated(mesh%edp_trsk_edpl_on_edge))           deallocate(mesh%edp_trsk_edpl_on_edge)
        if(allocated(mesh%edp_edge_on_edge))                deallocate(mesh%edp_edge_on_edge)
        if(allocated(mesh%edp_nedge))                       deallocate(mesh%edp_nedge)

        !> for calculating hx 2nd order derivative
        if(allocated(mesh%edt_v_weight_prime_cell))         deallocate(mesh%edt_v_weight_prime_cell)
        if(allocated(mesh%plg_stencil_number_2nd))          deallocate(mesh%plg_stencil_number_2nd)
        if(allocated(mesh%plg_stencil_index_2nd))           deallocate(mesh%plg_stencil_index_2nd)
        !> for calculating hx 4th order derivative
        if(allocated(mesh%edt_v_weight_4th_prime_cell))     deallocate(mesh%edt_v_weight_4th_prime_cell)
        if(allocated(mesh%edt_v_weight_4th_2nd_prime_cell)) deallocate(mesh%edt_v_weight_4th_2nd_prime_cell)
        if(allocated(mesh%plg_stencil_number_4th))          deallocate(mesh%plg_stencil_number_4th)
        if(allocated(mesh%plg_stencil_index_4th))           deallocate(mesh%plg_stencil_index_4th)
        !> for calculating tr 2nd order derivative
        if(allocated(mesh%edp_v_weight_dual_cell))          deallocate(mesh%edp_v_weight_dual_cell)
        if(allocated(mesh%tri_stencil_number_2nd))          deallocate(mesh%tri_stencil_number_2nd)
        if(allocated(mesh%tri_stencil_index_2nd))           deallocate(mesh%tri_stencil_index_2nd)
        !> for calculating digvergence
        if(allocated(mesh%vtx_nnb))                         deallocate(mesh%vtx_nnb)
        if(allocated(mesh%vtx_ed))                          deallocate(mesh%vtx_ed)
        if(allocated(mesh%plg_nr))                          deallocate(mesh%plg_nr)
        if(allocated(mesh%edp_leng))                        deallocate(mesh%edp_leng)
        if(allocated(mesh%plg_areag))                       deallocate(mesh%plg_areag)
        !> for calculating prim-cell flux
        if(allocated(mesh%edt_v))                           deallocate(mesh%edt_v)
        if(allocated(mesh%edt_edpNr_edtTg))                 deallocate(mesh%edt_edpNr_edtTg)
        if(allocated(mesh%edt_leng))                        deallocate(mesh%edt_leng)
        if(allocated(mesh%edt_my_edge_on_edge))             deallocate(mesh%edt_my_edge_on_edge)
        if(allocated(mesh%edt_ur_cell_on_edge))             deallocate(mesh%edt_ur_cell_on_edge)
        !> for calculating dual-cell flux
        if(allocated(mesh%edp_v))                           deallocate(mesh%edp_v)
        if(allocated(mesh%edp12_edp_tg))                    deallocate(mesh%edp12_edp_tg)
        if(allocated(mesh%edp12_edt_nr))                    deallocate(mesh%edp12_edt_nr)
        !> for calculating vorticity(curl)       
        if(allocated(mesh%tri_nnb))                         deallocate(mesh%tri_nnb)
        if(allocated(mesh%tri_ed))                          deallocate(mesh%tri_ed)
        if(allocated(mesh%tri_nr))                          deallocate(mesh%tri_nr)
        if(allocated(mesh%edt_edpTg_edtNr))                 deallocate(mesh%edt_edpTg_edtNr)
        if(allocated(mesh%tri_areag))                       deallocate(mesh%tri_areag)
        if(allocated(mesh%tri_c_lat))                       deallocate(mesh%tri_c_lat)
        !> for calculating kinetic_energy           
        if(allocated(mesh%vtx_tr))                          deallocate(mesh%vtx_tr)
        if(allocated(mesh%plg_kite_area))                   deallocate(mesh%plg_kite_area)
        !> added for VR
        if(allocated(mesh%vtxCellLeng))                     deallocate(mesh%vtxCellLeng)
        if(allocated(mesh%edt_scale_2nd))                   deallocate(mesh%edt_scale_2nd)
        if(allocated(mesh%edt_scale_4th))                   deallocate(mesh%edt_scale_4th)
        if(allocated(mesh%edt_scale_6th))                   deallocate(mesh%edt_scale_6th)
        !real(r8)                 :: minCellLeng
!
! pattern2 (20 in total, normal 15, following pattern1, to completely eliminate runtime mesh%${five_element}
! so as to improve efficiency
!
        if(allocated(mesh%vtx_lon))                         deallocate(mesh%vtx_lon)
        if(allocated(mesh%vtx_lat))                         deallocate(mesh%vtx_lat)
        if(allocated(mesh%vtx_p))                           deallocate(mesh%vtx_p)
        if(allocated(mesh%vtx_nb))                          deallocate(mesh%vtx_nb)
        if(allocated(mesh%vtx_lob_nx))                      deallocate(mesh%vtx_lob_nx)
        if(allocated(mesh%vtx_lob_ny))                      deallocate(mesh%vtx_lob_ny)
                                                                                                   
        if(allocated(mesh%plg_blocal))                      deallocate(mesh%plg_blocal)
        if(allocated(mesh%plg_sub_triangle_midp_lob))       deallocate(mesh%plg_sub_triangle_midp_lob)
        if(allocated(mesh%plg_sub_triangle_area))           deallocate(mesh%plg_sub_triangle_area)

        if(allocated(mesh%tri_c_p))                         deallocate(mesh%tri_c_p)
        if(allocated(mesh%tri_c_lon))                       deallocate(mesh%tri_c_lon)
        if(allocated(mesh%tri_v))                           deallocate(mesh%tri_v)
        if(allocated(mesh%tri_kite_area))                   deallocate(mesh%tri_kite_area)
                                                                         
        if(allocated(mesh%edt_c_p))                         deallocate(mesh%edt_c_p)
        if(allocated(mesh%edt_c_lon))                       deallocate(mesh%edt_c_lon)
        if(allocated(mesh%edt_c_lat))                       deallocate(mesh%edt_c_lat)
        if(allocated(mesh%edt_nr))                          deallocate(mesh%edt_nr)
                                                                         
        if(allocated(mesh%edp_c_p))                         deallocate(mesh%edp_c_p)
        if(allocated(mesh%edp_nr))                          deallocate(mesh%edp_nr)
        if(allocated(mesh%edp_tg))                          deallocate(mesh%edp_tg)

   end subroutine grist_domain_clean_runtime_allocated

 end module grist_domain_types
