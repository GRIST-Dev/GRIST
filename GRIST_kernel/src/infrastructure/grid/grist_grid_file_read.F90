
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: import mesh data from the grid file
!
! Revision history:
!----------------------------------------------------------------------------

module grist_grid_file_read
  
  use grist_lib
  use grist_wrap_pf
  ! core
  use grist_constants,             only: r8, i8, i4
  use grist_domain_types,          only: global_domain, global_domain_data
  use grist_nml_module,            only: gridFilePath, gridFileNameHead, index_flag, read_partition
  use grist_fileio_0d_module_gcm,      only: wrap_read_0d, wrap_bcast_0d
  use grist_fileio_list_2d_module_par, only: wrap_read_2d, wrap_bcast_2d, &
                                             wrap_read_2d_group
  use grist_util_module,           only: write_string
  use grist_data_types,            only: scalar_2d_field
  use grist_element_types_icosh,   only: nDual
  ! grid
  use grist_grid_file_vars

  implicit none

  private

  public  :: grist_read_grid_file, &
             grist_read_grid_file_par

  contains

  subroutine grist_read_grid_file_par(dm)
  
    include 'pnetcdf.inc'
    include 'mpif.h'
    type(global_domain), intent(inout), target :: dm
    ! local
    integer(i8)             :: dim_three
    integer(i8)             :: dim_two
    integer(i8)             :: dim_one
    integer(i8)             :: dim_Dual
    integer(i4)             :: startin
    integer(i4)             :: countin
    integer(i8)             :: starti
    integer(i8)             :: counti
!    character(128)          :: c_glevel
    character(128)          :: grid_file_0d_name
    character(128)          :: grid_file_2d_name
    integer                 :: i, j, k, option

    dim_three = 3
    dim_two   = 2
    dim_one   = 1
    dim_Dual  = nDual

    !================================================
    ! Create group
    !================================================

    !================================================
    ! [1]  construct filename
    !================================================
!    call write_string(mesh_glv, c_glevel)

    grid_file_0d_name  = trim(gridFileNameHead)//".0d.nc"
    grid_file_2d_name  = trim(gridFileNameHead)//".2d.nc"

    !================================================
    ! [2] read data
    !================================================
    !
    ! 2d
    !
    ! nt
    !
    call construct_grid_vars_for_read_local_t(size(dm%t_index))
    allocate(dm%tri(size(dm%t_index)))
    option = 1
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_v'      ,dim_Dual,  option, tri_v_local%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_cc_p'   ,dim_three, option, tri_cc_p%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_cc_ltln',dim_two,   option, tri_cc_ltln%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_ed'     ,dim_Dual,  option, tri_ed%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_nb'     ,dim_Dual,  option, tri_nb%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_area'   ,dim_one,   option, tri_area%f)
#ifdef CUBE
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_nnb'    ,dim_one,   option, tri_nnb%f)
#endif
#ifdef LAM_DOMAIN
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'tri_mark'   ,dim_one,   option, tri_mark%f)
    allocate(dm%tri_mark(size(dm%t_index)))
#endif

    do i = 1, size(dm%t_index)
       dm%tri(i)%c%p(1:3)    = tri_cc_p%f(i,1:3)
       dm%tri(i)%c%lat       = tri_cc_ltln%f(i,1)
       dm%tri(i)%c%lon       = tri_cc_ltln%f(i,2)
       dm%tri(i)%v(1:nDual)  = tri_v_local%f(i,1:nDual)
       dm%tri(i)%ed(1:nDual) = tri_ed%f(i,1:nDual)
       dm%tri(i)%nb(1:nDual) = tri_nb%f(i,1:nDual)
       dm%tri(i)%areag       = tri_area%f(i,1)
       dm%tri(i)%nnb         = 3      ! default for Voronoi-Delaunay grid 
#ifdef CUBE
       dm%tri(i)%nnb         = tri_nnb%f(i,1)
#endif
#ifdef LAM_DOMAIN
       dm%tri_mark(i)        = int(tri_mark%f(i,1))
#endif
   enddo

   call destruct_grid_file_vars_t

    !================================================
    ! [2] read data
    !================================================
    !
    ! 2d
    !
    ! ne
    !
    call construct_grid_vars_for_read_local_e(size(dm%e_index))
    allocate(dm%edt(size(dm%e_index)))
    allocate(dm%edp(size(dm%e_index)))
    option = 6
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edt_cc_p'   ,dim_three, option, edt_cc_p%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edt_cc_ltln',dim_two,   option, edt_cc_ltln%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edt_v'      ,dim_two,   option, edt_v_local%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edt_tg'     ,dim_three, option, edt_tg%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edt_nr'     ,dim_three, option, edt_nr%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edt_len'    ,dim_one,   option, edt_len%f)
#ifdef LAM_DOMAIN
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edt_mark'   ,dim_one,   option, edt_mark%f)
    allocate(dm%edt_mark(size(dm%e_index)))
#endif

    ! yizhang: we use intersecting point, so samee as edt's cc
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edp_cc_p'   ,dim_three, option, edp_cc_p%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edp_cc_ltln',dim_two,   option, edp_cc_ltln%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edp_v'      ,dim_two,   option, edp_v%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edp_tg'     ,dim_three, option, edp_tg%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edp_nr'     ,dim_three, option, edp_nr%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'edp_len'    ,dim_one,   option, edp_len%f)

    do i = 1, size(dm%e_index)
        dm%edt(i)%c%p(1:3) = edt_cc_p%f(i,1:3)
        dm%edt(i)%c%lat    = edt_cc_ltln%f(i,1)
        dm%edt(i)%c%lon    = edt_cc_ltln%f(i,2)
        dm%edt(i)%v(1:2)   = edt_v_local%f(i,1:2)
        dm%edt(i)%tg(1:3)  = edt_tg%f(i,1:3)
        dm%edt(i)%nr(1:3)  = edt_nr%f(i,1:3)
        dm%edt(i)%leng     = edt_len%f(i,1)

        dm%edp(i)%c%p(1:3) = edp_cc_p%f(i,1:3)
        dm%edp(i)%c%lat    = edp_cc_ltln%f(i,1)
        dm%edp(i)%c%lon    = edp_cc_ltln%f(i,2)
        dm%edp(i)%v(1:2)   = edp_v%f(i,1:2)
        dm%edp(i)%tg(1:3)  = edp_tg%f(i,1:3)
        dm%edp(i)%nr(1:3)  = edp_nr%f(i,1:3)
        dm%edp(i)%leng     = edp_len%f(i,1)
#ifdef LAM_DOMAIN
        dm%edt_mark(i)     = int(edt_mark%f(i,1))
#endif
    enddo
    call destruct_grid_file_vars_e

    !================================================
    ! [2] read data
    !================================================
    !
    ! 2d
    !
    ! nv
    !
    call construct_grid_vars_for_read_local_v(size(dm%v_index))
    allocate(dm%vtx(size(dm%v_index)))
    allocate(dm%plg(size(dm%v_index)))
    option = 0
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'vtx_p'       ,dim_three,   option, vtx_p%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'vtx_ltln_nnb',dim_three,   option, vtx_ltln_nnb_local%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'vtx_nb'      ,mesh_maxvnb, option, vtx_nb_local%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'vtx_ed'      ,mesh_maxvnb, option, vtx_ed_local%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'vtx_tr'      ,mesh_maxvnb, option, vtx_tr_local%f)
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'plg_area'    ,dim_one,     option, plg_area%f)
#ifdef LAM_DOMAIN
    call wrap_read_2d_group(dm%gcomm_read, gridFilePath,grid_file_2d_name,'plg_mark'    ,dim_one,     option, plg_mark%f)
    allocate(dm%plg_mark(size(dm%v_index)))
#endif

    do i = 1, size(dm%v_index)
        dm%vtx(i)%p(1:3)     = vtx_p%f(i,1:3)
        dm%vtx(i)%lat        = vtx_ltln_nnb_local%f(i,1)
        dm%vtx(i)%lon        = vtx_ltln_nnb_local%f(i,2)
        dm%vtx(i)%nnb        = vtx_ltln_nnb_local%f(i,3)

        allocate(  dm%vtx(i)%nb(1:dm%vtx(i)%nnb))
        allocate(  dm%vtx(i)%ed(1:dm%vtx(i)%nnb))
        allocate(  dm%vtx(i)%tr(1:dm%vtx(i)%nnb))
                               
        dm%vtx(i)%nb(:)      = vtx_nb_local%f(i,  1:dm%vtx(i)%nnb)
        dm%vtx(i)%ed(:)      = vtx_ed_local%f(i,  1:dm%vtx(i)%nnb)
        dm%vtx(i)%tr(:)      = vtx_tr_local%f(i,  1:dm%vtx(i)%nnb)
        dm%plg(i)%areag      = plg_area%f(i,1)
#ifdef LAM_DOMAIN
        dm%plg_mark(i)       = int(plg_mark%f(i,1))
#endif
   enddo
   call destruct_grid_file_vars_v

#ifdef LAM_DOMAIN
   !if(any(dm%plg_mark(:).eq.177).and.any(dm%edt_mark(:)/10.eq.177).and.any(dm%tri_mark(:).eq.177))then
   if(any(dm%plg_mark(:).eq.177))then
      dm%hasBdyDomain  = .true.
      print*,"this is a domain has assign-value region, i.e., mark index is 177"
   end if
#endif

   return
  end subroutine grist_read_grid_file_par

! read some global mesh info
  subroutine grist_read_grid_file(mesh, comm)
  
    include 'pnetcdf.inc'
    include 'mpif.h'

    type(global_domain_data), intent(inout) :: mesh
    integer,intent(in)                      :: comm
    ! local
    type(scalar_2d_field) :: vtx_ltln_nnb_l
    integer(i8)           :: dim_three
    integer(i8)           :: dim_two
    integer(i8)           :: dim_one
    integer(i8)           :: dim_Dual
!    character(128)        :: c_glevel
    character(128)        :: grid_file_0d_name
    character(128)        :: grid_file_2d_name

    dim_three = 3
    dim_two   = 2
    dim_one   = 1
    dim_Dual  = nDual

    !================================================
    ! Create group
    !================================================

    !================================================
    ! [1]  construct filename
    !================================================
!    call write_string(mesh_glv, c_glevel)

    grid_file_0d_name  = trim(gridFileNameHead)//".0d.nc"
    
    grid_file_2d_name  = trim(gridFileNameHead)//".2d.nc"

    if(mpi_rank() == 0)then
      print *,"Input 0d file:",trim(grid_file_0d_name)
      print *,"Input 2d file:",trim(grid_file_2d_name)
    end if
    !================================================
    ! [2] read data
    !================================================
    !
    ! 0d
    !
    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_nt'      , mesh_nt)
    end if
    call wrap_bcast_0d(comm, 0, mesh_nt)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_ne'      , mesh_ne)
    end if
    call wrap_bcast_0d(comm, 0, mesh_ne)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_nv'      , mesh_nv)
    end if
    call wrap_bcast_0d(comm, 0, mesh_nv)
 
    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_maxvnb'  , mesh_maxvnb)
    end if
    call wrap_bcast_0d(comm, 0, mesh_maxvnb)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_kind'          , mesh_kind)
    end if
    call wrap_bcast_0d(comm, 0, mesh_kind)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_optm'          , mesh_optm)
    end if
    call wrap_bcast_0d(comm, 0, mesh_optm)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_glevel'        , mesh_glevel)
    end if
    call wrap_bcast_0d(comm, 0, mesh_glevel)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_min_edt_dist'  , mesh_min_edt_dist)
    end if
    call wrap_bcast_0d(comm, 0, mesh_min_edt_dist)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_max_edt_dist'  , mesh_max_edt_dist)
    end if
    call wrap_bcast_0d(comm, 0, mesh_max_edt_dist)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_mean_edt_dist' , mesh_mean_edt_dist)
    end if
    call wrap_bcast_0d(comm, 0, mesh_mean_edt_dist)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_min_edp_dist'    , mesh_min_edp_dist)
    end if
    call wrap_bcast_0d(comm, 0, mesh_min_edp_dist)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_max_edp_dist'    , mesh_max_edp_dist)
    end if
    call wrap_bcast_0d(comm, 0, mesh_max_edp_dist)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_mean_edp_dist'   , mesh_mean_edp_dist)
    end if
    call wrap_bcast_0d(comm, 0, mesh_mean_edp_dist)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_min_tri_area'  , mesh_min_tri_area)
    end if
    call wrap_bcast_0d(comm, 0, mesh_min_tri_area)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_max_tri_area'  , mesh_max_tri_area)
    end if
    call wrap_bcast_0d(comm, 0, mesh_max_tri_area )

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_mean_tri_area' , mesh_mean_tri_area)
    end if
    call wrap_bcast_0d(comm, 0, mesh_mean_tri_area)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_min_plg_area'  , mesh_min_plg_area)
    end if
    call wrap_bcast_0d(comm, 0, mesh_min_plg_area)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_max_plg_area'  , mesh_max_plg_area)
    end if
    call wrap_bcast_0d(comm, 0, mesh_max_plg_area)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_mean_plg_area' , mesh_mean_plg_area)
    end if
    call wrap_bcast_0d(comm, 0, mesh_mean_plg_area)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_min_tri_angle' , mesh_min_tri_angle)
    end if
    call wrap_bcast_0d(comm, 0, mesh_min_tri_angle)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_max_tri_angle' , mesh_max_tri_angle)
    end if
    call wrap_bcast_0d(comm, 0, mesh_max_tri_angle)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_mean_tri_angle', mesh_mean_tri_angle)
    end if
    call wrap_bcast_0d(comm, 0, mesh_mean_tri_angle)

    if(mpi_rank()==0)then
      call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_name'          , mesh_name)
    end if
    call wrap_bcast_0d(comm, 0, mesh_name)

!    call construct_grid_vars_for_read

    if (.not. read_partition) then
       ! nt
       allocate(mesh%tri_v%f(mesh_nt,nDual))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'tri_v', mesh_nt, dim_Dual, mesh%tri_v)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nt, dim_Dual, mesh%tri_v)

       allocate(mesh%tri_nnb%f(mesh_nt,1))
#ifdef CUBE
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'tri_nnb', mesh_nt, dim_one, mesh%tri_nnb)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nt, dim_one, mesh%tri_nnb)
#else
       mesh%tri_nnb%f = 3
#endif

#ifdef LAM_DOMAIN
       allocate(mesh%tri_mark%f(mesh_nt,1))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'tri_mark', mesh_nt, dim_one, mesh%tri_mark)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nt, dim_one, mesh%tri_mark)
#endif

       ! ne
       allocate(mesh%edt_v%f(mesh_ne,2))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name,'edt_v', mesh_ne, dim_two, mesh%edt_v)
       end if
       call wrap_bcast_2d(comm, 0, mesh_ne, dim_two, mesh%edt_v)

#ifdef LAM_DOMAIN
       allocate(mesh%edt_mark%f(mesh_ne,1))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'edt_mark', mesh_ne, dim_one, mesh%edt_mark)
       end if
       call wrap_bcast_2d(comm, 0, mesh_ne, dim_one, mesh%edt_mark)
#endif

       ! nv
       allocate(vtx_ltln_nnb_l%f(mesh_nv,3))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_ltln_nnb', mesh_nv, dim_three, vtx_ltln_nnb_l)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nv, dim_three, vtx_ltln_nnb_l)
       allocate(mesh%vtx_nnb(mesh_nv))
       mesh%vtx_nnb(:) = int(vtx_ltln_nnb_l%f(:,3))
       if(index_flag .eq. "morton") then
         allocate(mesh%vtx_ltln(mesh_nv,2))
         mesh%vtx_ltln(:,:) = vtx_ltln_nnb_l%f(:,1:2)
       end if
       deallocate(vtx_ltln_nnb_l%f)
       allocate(mesh%vtx_nb%f(mesh_nv,mesh_maxvnb))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_nb', mesh_nv, mesh_maxvnb, mesh%vtx_nb)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nv, mesh_maxvnb, mesh%vtx_nb)

       allocate(mesh%vtx_ed%f(mesh_nv,mesh_maxvnb))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_ed', mesh_nv, mesh_maxvnb, mesh%vtx_ed)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nv, mesh_maxvnb, mesh%vtx_ed)    

       allocate(mesh%vtx_tr%f(mesh_nv,mesh_maxvnb))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_tr', mesh_nv, mesh_maxvnb, mesh%vtx_tr)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nv, mesh_maxvnb, mesh%vtx_tr)

#ifdef LAM_DOMAIN
       allocate(mesh%plg_mark%f(mesh_nv,1))
       if(mpi_rank()==0)then
         call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'plg_mark', mesh_nv, dim_one, mesh%plg_mark)
       end if
       call wrap_bcast_2d(comm, 0, mesh_nv, dim_one, mesh%plg_mark)
#endif
    end if

    return
  end subroutine grist_read_grid_file

end module grist_grid_file_read
