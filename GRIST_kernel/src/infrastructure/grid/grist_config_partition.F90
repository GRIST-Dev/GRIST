
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: domain decomposition and exchange_data functions
!
! Revision history:
!----------------------------------------------------------------------------

module grist_config_partition
  use grist_lib
  use grist_nml_module
  use grist_grid_file_read,      only: grist_read_grid_file_par
  use grist_constants,           only: i4, r8, i8, pi, r4
  use grist_domain_types,        only: global_domain, global_domain_data, block_structure, group_comm
  use grist_list,                only: index_list
  use grist_element_types_icosh, only: node_structure  , &
                                       edge_structure  , &
                                       dual_structure  , &
                                       prim_structure

  use grist_util_module,         only: write_string
  !use grist_domain_types
  use grist_data_types,          only: scalar_1d_field,&
                                       scalar_2d_field,&
                                       scalar_3d_field,&
                                       exchange_field_list_1d,&
                                       exchange_field_list_2d,&
                                       exchange_field_list_3d
  use grist_wrap_pf,             only: wrap_open, wrap_close
  use grist_grid_file_vars

  public  :: config_sub_domain , &
             domain_decompse ,&
             debug_data_1d ,&
             debug_data_2d ,&
             exchange_data_1d_add , &
             exchange_data_1d     , &
             exchange_data_2d_add , &
             exchange_data_2d     , &
             exchange_data_3d_add , &
             exchange_data_3d     , &
             exchange_data_1d_r4  , &
             exchange_data_2d_r4  , &
             exchange_data_3d_r4  , &
             show_basic_info

  private
 
!grid nv per core
  integer              :: nv_per_core = 20
!for LAM: record the number of cells in the compute and relax regions
  integer              :: region_cells_compute
  integer, parameter   :: maxi4 = huge(i4)
  integer, allocatable :: v_index_nomask(:)
  
contains

  !> partition the mesh, partition flag for each node
  !> is saved in mesh%parts
  subroutine domain_decompse(mesh, comm)
    implicit none
    type(global_domain_data), intent(inout) :: mesh
    integer                                 :: comm !MPI communicator
    integer(i4), allocatable                :: xadj(:), adjncy(:)
    integer(i4)                             :: i, j, x, act_nnb
    integer                                 :: rank, nprocs, err, c_idx, nb_idx, mark
    integer(i4),              allocatable   :: parts_metis(:)
    type(set)                               :: nomask
    type(map)                               :: map_nomask

    if (.not. read_partition) then
       allocate(mesh%parts(mesh%nv))
       if (mpi_rank() .eq. 0) then
          allocate(xadj(mesh%nv+1), adjncy(mesh%nv*mesh%maxvnb))
#ifdef LAM_DOMAIN
          region_cells_compute = 0
          ! initialize unmarked set and map
          do i = 1, mesh%nv
             mark = mesh%plg_mark%f(i, 1)
             if (mark .ne. 177) then
                region_cells_compute = region_cells_compute + 1
                call nomask%insert(i)
                call map_nomask%insert(i, region_cells_compute)
             end if
          end do
          call nomask%dump(v_index_nomask)
          ! set xadj and adjncy
          x = 1
          xadj(1) = 0
          do i = 2, region_cells_compute+1
             act_nnb = 0
             c_idx = v_index_nomask(i-1)
             do j = 1, mesh%vtx_nnb(c_idx)
                nb_idx = mesh%vtx_nb%f(c_idx, j)
                mark = mesh%plg_mark%f(nb_idx, 1)
                if (mark .ne. 177) then
                   adjncy(x) = map_nomask%find(nb_idx)-1
                   x = x + 1
                   act_nnb = act_nnb + 1
                end if
             end do
             xadj(i) = xadj(i-1) + act_nnb
          end do
          ! call metis
          nprocs = mpi_size(comm)
          allocate(parts_metis(region_cells_compute))
          call metis_decomp(region_cells_compute, adjncy, xadj, nprocs, parts_metis)
          mesh%parts = maxi4
          do i = 1, region_cells_compute
             mesh%parts(v_index_nomask(i)) = parts_metis(i)
          end do
          deallocate(parts_metis)
          call nomask%final
          call map_nomask%final
#else
          ! set xadj and adjncy
          x = 1
          xadj(1) = 0
          do i = 2, mesh%nv+1
             do j = 1, mesh%vtx_nnb(i-1)
                adjncy(x) = mesh%vtx_nb%f(i-1, j)-1
                x = x + 1
             end do
             xadj(i) = xadj(i-1) + mesh%vtx_nnb(i-1)
          end do
          ! call metis
          nprocs = mpi_size(comm)
          call metis_decomp(mesh%nv, adjncy, xadj, nprocs, mesh%parts)
          deallocate(xadj, adjncy)
#endif
       end if
#ifdef LAM_DOMAIN
       call MPI_Bcast(region_cells_compute, 1, MPI_INTEGER, 0, comm, err)
       if (mpi_rank() .ne. 0) allocate(v_index_nomask(region_cells_compute))
       call MPI_Bcast(v_index_nomask, region_cells_compute, MPI_INTEGER, 0, comm, err)
#endif
       call MPI_Bcast(mesh%parts, mesh%nv, MPI_INTEGER, 0, comm, err)
    end if

  end subroutine


  subroutine search_io(gmesh, local_block, comm)
    implicit none
    type(global_domain_data),      intent(in)    :: gmesh !global mesh
    type(block_structure), target, intent(inout) :: local_block
    integer                                      :: comm !MPI communicator
    integer                                      :: iv, ie, k, iv_global, i, ierr
    type(global_domain),   pointer               :: full
    type(set)                                    :: e_list, t_list
    logical                                      :: insert_ed, insert_tr
    integer(i4)                                  :: ed, tr, nv, ne, nt, v_tmp, act_iv
    integer                                      :: belong_cpu, ele3(3)
    integer,               allocatable           :: points(:)
    integer,               allocatable           :: counts(:), displs(:), ele3_all(:), &
                                                    buffer_all(:), buffer(:)
    integer                                      :: iblock, nprocs, pos, sums

    iblock = local_block%cpuid
    nprocs = local_block%nprocs
    if (.not. read_partition) then

       full => local_block%full_domain
       ! initialize output v_index
#ifdef LAM_DOMAIN
       allocate(points(full%nv))
       act_iv = 0
       do iv = 1, full%nv
          if (gmesh%parts(full%v_index(iv)) .eq. iblock) then
             act_iv = act_iv + 1
             points(act_iv) = full%v_index(iv)
          end if
       end do
       allocate(local_block%output%v(act_iv))
       local_block%output%v = points(1:act_iv)
       deallocate(points)
#else
       allocate(local_block%output%v(full%nv_compute))
       local_block%output%v  = full%v_index(1:full%nv_compute)
#endif

       do iv = 1, size(local_block%output%v)
          iv_global = local_block%output%v(iv)

          do ie = 1, gmesh%vtx_nnb(iv_global)
             ed = gmesh%vtx_ed%f(iv_global, ie)
             tr = gmesh%vtx_tr%f(iv_global, ie)

             insert_ed = .true.
             insert_tr = .true.

             !check trianglar point belong to which part
             do k = 1, int(gmesh%tri_nnb%f(tr,1))
#ifdef LAM_DOMAIN
                v_tmp = gmesh%tri_v%f(tr, k)
                if (v_tmp .gt. 0) then
#endif
                belong_cpu = gmesh%parts(int(gmesh%tri_v%f(tr, k)))
                if(belong_cpu > local_block%cpuid) then
                   insert_tr = .false.
                end if
#ifdef LAM_DOMAIN
                end if
#endif
             end do

             !check edge point belong to which part
             do k = 1, 2
#ifdef LAM_DOMAIN
                v_tmp = gmesh%edt_v%f(ed, k)
                if (v_tmp .gt. 0) then
#endif
                belong_cpu = gmesh%parts(int(gmesh%edt_v%f(ed, k)))
                if(belong_cpu > local_block%cpuid) then
                   insert_ed = .false.
                end if
#ifdef LAM_DOMAIN
                end if
#endif
             end do

             if(insert_ed) call e_list%insert(ed)

             if(insert_tr) call t_list%insert(tr)

          end do
       end do

       call e_list%dump(local_block%output%e)
       call t_list%dump(local_block%output%t)
       call e_list%final
       call t_list%final
       ! test
       !open(101, file="output-v-"//i2s(iblock+1)//".txt", form='formatted')
       !write(101, *) size(local_block%output%v)
       !do iv = 1, size(local_block%output%v)
       !   write(101, *) local_block%output%v(iv)
       !end do
       !close(101)
       !open(101, file="output-e-"//i2s(iblock+1)//".txt", form='formatted')
       !write(101, *) size(local_block%output%e)
       !do iv = 1, size(local_block%output%e)
       !   write(101, *) local_block%output%e(iv)
       !end do
       !close(101)
       !open(101, file="output-t-"//i2s(iblock+1)//".txt", form='formatted')
       !write(101, *) size(local_block%output%t)
       !do iv = 1, size(local_block%output%t)
       !   write(101, *) local_block%output%t(iv)
       !end do
       !close(101)
    else
      if(iblock .eq. 0) then
        open(101, file=trim(pardir)//"output-"//i2s(nprocs), form='unformatted', access="stream")
        allocate(ele3_all(3*nprocs))
        read(101) ele3_all
        allocate(counts(nprocs))
        pos = 1
        do i = 1, nprocs
          counts(i) = ele3_all(pos) + ele3_all(pos+1) + ele3_all(pos+2)
          pos = pos + 3
        end do
        sums = sum(counts)
        allocate(buffer_all(sums))
        read(101) buffer_all
        close(101)
        allocate(displs(nprocs))
        displs(1) = 0
        do i = 2, nprocs
          displs(i) = displs(i-1) + counts(i-1)
        end do
      else
        allocate(ele3_all(1))
        allocate(counts(1))
        allocate(displs(1))
        allocate(buffer_all(1))
      end if
      call MPI_Scatter(ele3_all, 3, MPI_INTEGER, ele3, &
                       3, MPI_INTEGER, 0, local_block%comm, ierr)
      sums = sum(ele3)
      allocate(buffer(sums))
      call MPI_Scatterv(buffer_all, counts, displs, MPI_INTEGER, buffer, &
                        sums, MPI_INTEGER, 0, local_block%comm, ierr)
      allocate(local_block%output%v(ele3(1)))
      allocate(local_block%output%t(ele3(2)))
      allocate(local_block%output%e(ele3(3)))
      local_block%output%v = buffer(1:ele3(1))
      local_block%output%t = buffer(ele3(1)+1:ele3(1)+ele3(2))
      local_block%output%e = buffer(ele3(1)+ele3(2)+1:ele3(1)+ele3(2)+ele3(3))
      deallocate(ele3_all)
      deallocate(counts)
      deallocate(displs)
      deallocate(buffer_all)
      deallocate(buffer)
    end if

  end subroutine


  subroutine write_output(local_block)
    implicit none
    type(block_structure), intent(in) :: local_block
    
    open(1, file="write_v."//i2s(local_block%cpuid)//".txt")
    write(1, "(I0)", advance='yes') local_block%output%v
    close(1)

    open(1, file="write_e."//i2s(local_block%cpuid)//".txt")
    write(1, "(I0)", advance='yes') local_block%output%e
    close(1)

    open(1, file="write_t."//i2s(local_block%cpuid)//".txt")
    write(1, "(I0)", advance='yes') local_block%output%t
    close(1)
  end subroutine

  
  subroutine write_domain(dm, tmp_string)
    implicit none
    type(global_domain), intent(in) :: dm
    character(len=*) :: tmp_string
    integer :: i

    open(1,file=trim(tmp_string)//trim(".t.txt"))
    do i = 1, dm%nt
       write(1,*) dm%t_index(i)
    end do
    close(1)

    open(1,file=trim(tmp_string)//trim(".e.txt"))
    do i = 1, dm%ne
       write(1,*) dm%e_index(i)
    end do
    close(1)

    open(1,file=trim(tmp_string)//trim(".v.txt"))
    do i = 1, dm%nv
       write(1,*) dm%v_index(i)
    end do
    close(1)

  end subroutine


  subroutine config_sub_domain(mesh, local_block, comm)
    !use grist_fileio_list_1d_module_par, only: wrap_output_init_1d    ,&
    !                                           wrap_add_field_1d      ,&
    !                                           wrap_output_1d_group   ,&
    !                                           wrap_output_1d_group_sp,&
    !                                           wrap_output_clean_1d
    implicit none
    !IO
    type(global_domain_data),      intent(inout) :: mesh
    type(block_structure), target, intent(inout) :: local_block ! partition
    integer,                       intent(in)    :: comm        ! MPI communicator

    !local
    integer(i4)                                  :: iblock
    ! mesh, domain, local nb
    logical                                      :: write_info
    integer                                      :: sw
    !type(scalar_1d_field)                        :: scalar_io_1d
    !type(exchange_field_list_1d),  pointer       :: field_head_1d

    iblock = mpi_rank(comm)

    local_block%comm              = comm
    local_block%cpuid             = iblock
    local_block%nprocs            = mpi_size(comm)
    local_block%full_domain%nv    = 0
    local_block%iner_domain%nv    = 0
    local_block%full_domain%glevel = mesh%glevel

    local_block%full_domain%nv_all = mesh%nv
    local_block%full_domain%nt_all = mesh%nt
    local_block%full_domain%ne_all = mesh%ne
    local_block%full_domain%maxvnb = mesh%maxvnb

    allocate(local_block%bdry(stencil_width))
    allocate(local_block%halo(stencil_width))

    allocate(local_block%full_domain%nv_halo(stencil_width))
    allocate(local_block%full_domain%ne_halo(stencil_width))
    allocate(local_block%full_domain%nt_halo(stencil_width))

    allocate(local_block%full_domain%nv_bdry(stencil_width))
    allocate(local_block%full_domain%ne_bdry(stencil_width))
    allocate(local_block%full_domain%nt_bdry(stencil_width))
    ! search and acquire all t, e, v based global to local domain's array

    call search_compute(mesh,    local_block)

    call search_bdry_halo(mesh,  local_block)

    call search_exchange(mesh, local_block)

    ! pull all element structs for necessary domain

    call pull_element_types_full(local_block)

    call search_io(mesh, local_block, comm)

! zhangyi comment on sgpi     
!    call write_output(local_block)

    call global2local(local_block)

    local_block%full_domain%local_block => local_block
    call grist_output_init(local_block%full_domain)
    ! write dm info
    write_info = .false.
    if(write_info)then

       call write_domain(local_block%full_domain, "full_"//i2s(iblock))
       call write_domain(local_block%compute_domain, "compute_"//i2s(iblock))
       call write_domain(local_block%iner_domain, "iner_"//i2s(iblock))

       !------------BDRY and HALO DOMAINS----------
       do sw = 1, stencil_width
          call write_domain(local_block%bdry(sw), "bdry"//i2s(sw)//"_"//i2s(iblock))
          call write_domain(local_block%halo(sw), "halo"//i2s(sw)//"_"//i2s(iblock))
       end do

    end if
    !! test write var
    !allocate(scalar_io_1d%f(local_block%full_domain%ne_full))
    !scalar_io_1d%f = mpi_rank()
    !scalar_io_1d%pos = 6
    !call wrap_output_init_1d(local_block%full_domain)
    !call wrap_add_field_1d(scalar_io_1d, "var" ,'','')
    !call wrap_output_1d_group(local_block%full_domain, './', "test.nc")
    !call wrap_output_clean_1d()

    !! test exchange var
    !scalar_io_1d%f = 0
    !do sw = 1, local_block%full_domain%ne_compute
    !   scalar_io_1d%f(sw) = local_block%full_domain%e_index(sw)
    !end do
    !field_head_1d=>null()
    !call exchange_data_1d_add(local_block%full_domain,field_head_1d,scalar_io_1d)
    !call exchange_data_1d(local_block,field_head_1d)

    !call wrap_output_init_1d(local_block%full_domain)
    !call wrap_add_field_1d(scalar_io_1d, "var" ,'','')
    !call wrap_output_1d_group(local_block%full_domain, './', "test2.nc")
    !call wrap_output_clean_1d()
    !deallocate(scalar_io_1d%f)

    return
  end subroutine config_sub_domain

  !====================================================
  !               PRIVATE ROUTINES
  !====================================================

  subroutine search_compute(mesh,local_block)
    implicit none
    ! io
    type(global_domain_data),      intent(in)    :: mesh
    type(block_structure), target, intent(inout) :: local_block

    ! local
    integer(i4)                                  :: iv_domain
    integer(i4)                                  :: iv_global
    integer(i4)                                  :: iv_local
    integer(i4)                                  :: idx, sums
    type(global_domain),   pointer               :: compute
    integer                                      :: iblock, nprocs
    integer(i4)                                  :: i, j, it_local

    type(set)                                    :: bfs_index, hilbert_index, morton_index
    integer(i4),           allocatable           :: hilbert_cartesian(:,:),queue(:),z_index(:)
    integer(i4)                                  :: head, fina, nnb, nb, ierr, f_head,  h_index, &
                                                    curr, minx, maxx, miny, maxy, maxi, maxj, &
                                                    maxtimes, times, square, starti, finallyi, &
                                                    startj, finallyj, list_n, hilbert_pos, hilbert_n,&
                                                    minsquare, temp1, temp2, cx, cy, diffindex, precs,&
                                                    ztmpid, x, y, z, start, diffstep
    logical                                      :: flag, flag1, flag2
    integer(i4),           allocatable           :: binary_code(:,:)    
    character(len = 5)                           :: ztmp 
    character(len = 5),    allocatable           :: z_list(:)
    character(len = 1),    allocatable           :: geo_hash(:)
    real(r8)                                     :: lon, lat, lond, latd
    integer,               allocatable           :: counts(:), displs(:), buffer_all(:)
    ! acquire t and e based index

    iblock = local_block%cpuid
    nprocs = local_block%nprocs
    compute => local_block%compute_domain

    if ( .not. read_partition) then

#ifdef LAM_DOMAIN
       do i = 1, region_cells_compute
          curr = v_index_nomask(i)
          if(mesh%parts(curr) == iblock) then
             call compute%v_index_list%insert(curr)
          end if
       end do
#else
       do i = 1, mesh%nv
          if(mesh%parts(i) == iblock) then
             call compute%v_index_list%insert(i)
          end if
       end do
#endif

       call compute%v_index_list%dump(compute%v_index)
       compute%nv = size(compute%v_index)

       !===================================
       ! reordering index of compute domain
       !===================================
       SELECT CASE(index_flag)
       CASE('bfs')
          if(iblock == 0) print*, "reordering index using method: [BFS]"
          allocate(queue(compute%nv))
          do i = 1, compute%nv
              queue(i) = -1
          enddo

          queue(1) = compute%v_index(1)
          head = 1
          fina = 2
          do while(head .ne. fina)
              call bfs_index%insert(queue(head))
              nnb = mesh%vtx_nnb(queue(head))
              do i = 1, nnb
                  nb = mesh%vtx_nb%f(queue(head), i)
                  flag1 = compute%v_index_list%find(nb)
                  flag2 = .true.
                  do j = 1, fina-1
                      if(nb == queue(j)) then
                          flag2 = .false.
                      endif
                  enddo
                  if(flag1.and.flag2) then
                      queue(fina) = nb
                      fina = fina + 1
                  endif
              enddo
              head = head +1
          enddo

          do i = 1, compute%nv
             ! Insert isolated point
             if (.not. bfs_index%find(compute%v_index(i))) then
                call bfs_index%insert(compute%v_index(i))
             end if
          end do
          call bfs_index%dump(compute%v_index)
          call bfs_index%final

          deallocate(queue)
          if(iblock == 0) print*,"index reordering by [BFS] end"

       CASE('hilbert')
          ! Currently, the algorithm used in the Hilbert curve only
          ! sort the vertices with 6 neighbors
          if(iblock == 0) print*, "reordering index using SFC: [Hilbert]"
          allocate(queue(compute%nv))
          allocate(hilbert_cartesian(5,compute%nv))
          head = 1
          fina = head + 1
          hilbert_cartesian = 0
          ! Find the first vtx with 6 neighbors
          do i = 1, compute%nv
             if (mesh%vtx_nnb(compute%v_index(i)) == 6) then
                queue(head) = compute%v_index(i)
                hilbert_cartesian(1,head) = queue(head)
                hilbert_cartesian(2,head) = 0
                hilbert_cartesian(3,head) = 0
                hilbert_cartesian(4,head) = -1
                hilbert_cartesian(5,head) = -1
                exit
             end if
          end do
          curr = 0
          ! For vertices with 6 neighbors, get the x, y coordinates
          ! (hilbert_cartesian(2,:) hilbert_cartesian(3,:)), for non
          ! 6-neighbor vertices, x = 0, y = 0 (default)
          do while(head .ne. fina)
             nnb = mesh%vtx_nnb(queue(head))
             x = hilbert_cartesian(2,head)
             y = hilbert_cartesian(3,head)
             f_head = hilbert_cartesian(4,head)
             h_index = hilbert_cartesian(5,head)
             if (nnb == 6) then
                if (head == 1) then
                  diffindex = 0
                else
                   do i = 1, nnb
                      nb = mesh%vtx_nb%f(queue(head), i)
                      if (nb == hilbert_cartesian(1,f_head)) then
                         diffindex = h_index + 3 - i
                         exit
                      endif
                   enddo
                endif

                do i = 1, nnb
                   nb = mesh%vtx_nb%f(queue(head), i)
                   flag1 = compute%v_index_list%find(nb)
                   flag2 = .true.
                   do j = 1, fina-1
                      if (nb == queue(j)) then
                         flag2 = .false.
                      endif
                   enddo

                   if (flag1.and.flag2) then
                      curr = i + diffindex
                      if (curr .gt. 6) then
                         curr = mod(curr,6)
                      endif
                      hilbert_cartesian(1,fina) = nb
                      if (curr == 1) then
                         hilbert_cartesian(2,fina) = x
                         hilbert_cartesian(3,fina) = y + 1
                         hilbert_cartesian(4,fina) = head
                         hilbert_cartesian(5,fina) = 1
                      else if(curr == 2) then
                         hilbert_cartesian(2,fina) = x - 1
                         hilbert_cartesian(3,fina) = y + 1
                         hilbert_cartesian(4,fina) = head
                         hilbert_cartesian(5,fina) = 2
                      else if(curr == 3) then
                         hilbert_cartesian(2,fina) = x - 1
                         hilbert_cartesian(3,fina) = y
                         hilbert_cartesian(4,fina) = head
                         hilbert_cartesian(5,fina) = 3
                      else if(curr == 4) then
                         hilbert_cartesian(2,fina) = x
                         hilbert_cartesian(3,fina) = y - 1
                         hilbert_cartesian(4,fina) = head
                         hilbert_cartesian(5,fina) = 4
                      else if(curr == 5) then
                         hilbert_cartesian(2,fina) = x + 1
                         hilbert_cartesian(3,fina) = y - 1
                         hilbert_cartesian(4,fina) = head
                         hilbert_cartesian(5,fina) = 5
                      else
                         hilbert_cartesian(2,fina) = x + 1
                         hilbert_cartesian(3,fina) = y
                         hilbert_cartesian(4,fina) = head
                         hilbert_cartesian(5,fina) = 6
                      endif
                      queue(fina) = nb
                      fina = fina + 1
                   endif
                enddo
             endif
             head = head + 1
          enddo

          ! Set the minimum x and y to 0
          minx = maxi4
          maxx = 0 - maxi4
          miny = maxi4
          maxy = 0 - maxi4
          do i = 1, compute%nv
             if(hilbert_cartesian(2,i).lt.minx) then
                minx = hilbert_cartesian(2,i)
             endif
             if(hilbert_cartesian(2,i).gt.maxx) then
                maxx = hilbert_cartesian(2,i)
             endif
             if(hilbert_cartesian(3,i).lt.miny) then
                miny = hilbert_cartesian(3,i)
             endif
             if(hilbert_cartesian(3,i).gt.maxy) then
                maxy = hilbert_cartesian(3,i)
             endif
          enddo
          do i = 1,compute%nv
             hilbert_cartesian(2,i) = hilbert_cartesian(2,i) -  minx
             hilbert_cartesian(3,i) = hilbert_cartesian(3,i) - miny
          enddo

          maxi = maxx - minx + 1
          maxj = maxy - miny + 1
          minsquare = 2
          do while((2**minsquare.lt.maxi).or.(2**minsquare.lt.maxj))
             minsquare = minsquare + 1
          enddo
          minsquare = minsquare + 1

          ! Get the Hilbert pos according to the x, y coordinates
          allocate(binary_code(2, compute%nv))
          square = 2**minsquare
          do i = 1, compute%nv
             cx = hilbert_cartesian(2,i)
             cy = hilbert_cartesian(3,i)
             hilbert_pos = 0
             call xy2d(square, cx, cy, hilbert_pos)
             binary_code(1,i) = hilbert_cartesian(1,i)
             binary_code(2,i) = hilbert_pos
          enddo
          ! Sort the vtx indices according to the Hilbert pos
          do i = 1, compute%nv - 1
             do j = 1, compute%nv - i
                if(binary_code(2,j) > binary_code(2,j+1)) then
                   temp1 = binary_code(2,j)
                   binary_code(2,j) = binary_code(2,j+1)
                   binary_code(2,j+1) = temp1
                   temp2 = binary_code(1,j)
                   binary_code(1,j) = binary_code(1,j+1)
                   binary_code(1,j+1) = temp2
                endif
             enddo
          enddo
          ! Assignment
          do i = 1, compute%nv
             ! Avoid 0 index (default for isolated compute vtx)
             if(binary_code(1,i) .ne. 0) call hilbert_index%insert(binary_code(1,i))
          enddo
          do i = 1, compute%nv
             ! Push back isolated compute vtx
             if (.not. hilbert_index%find(compute%v_index(i))) then 
                 call hilbert_index%insert(compute%v_index(i))
             end if
          end do
          call hilbert_index%dump(compute%v_index)
          call hilbert_index%final

          deallocate(queue)
          deallocate(hilbert_cartesian)
          deallocate(binary_code)
          if(iblock == 0) print*,"index reordering by [Hilbert] end"

       CASE('morton')
          if(iblock == 0) print*, "reordering index using SFC: [Morton]"
          precs = 5
          allocate(z_index(compute%nv))
          allocate(z_list(compute%nv))
          allocate(geo_hash(1:precs))
          ! Get hash code for each vtx (z_list)
          do i = 1, compute%nv
             iv_global = compute%v_index(i)
             lat = mesh%vtx_ltln(iv_global, 1)
             lon = mesh%vtx_ltln(iv_global, 2)
             latd = lat * 180 / pi
             lond = lon * 180 / pi
             call geohash_encode(latd, lond, precs, geo_hash)

             ztmp = geo_hash(1)
             do j = 2, precs
                ztmp = trim(ztmp)//geo_hash(j)
             end do
                
             z_list(i) = ztmp
             z_index(i) = iv_global
          enddo
          ! Sort the index according the hash code
          do i = 1, compute%nv - 1
             do j = 1, compute%nv - i
                if(z_list(j) > z_list(j+1)) then
                   ztmp = z_list(j)
                   z_list(j) = z_list(j+1)
                   z_list(j+1) = ztmp
                   ztmpid = z_index(j)
                   z_index(j) = z_index(j+1)
                   z_index(j+1) = ztmpid
                endif
             enddo
          enddo
          ! Assignment
          do i = 1, compute%nv
             call morton_index%insert(z_index(i))
          enddo
          call morton_index%dump(compute%v_index)
          call morton_index%final

          deallocate(z_index)
          deallocate(z_list)
          deallocate(geo_hash)
          if(iblock == 0) print*,"index reordering by [Morton] end"

       CASE DEFAULT
          if(iblock == 0) print*,"no selection of SFC: [random] "
       END SELECT

       do iv_domain = 1, compute%nv 
          iv_global = compute%v_index(iv_domain)

          do it_local = 1, mesh%vtx_nnb(iv_global)

             idx = mesh%vtx_tr%f(iv_global, it_local)
             call compute%t_index_list%insert(idx)

             idx = mesh%vtx_ed%f(iv_global, it_local)
             call compute%e_index_list%insert(idx)

          end do
       end do

       ! from linked list to array
       call compute%t_index_list%dump(compute%t_index)
       call compute%e_index_list%dump(compute%e_index)

       compute%nt = size(compute%t_index)
       compute%ne = size(compute%e_index)
       call compute%e_index_list%final
       call compute%t_index_list%final
#ifdef LAM_DOMAIN
       deallocate(v_index_nomask)
#endif
       !!test
       !open(101, file="compute-"//i2s(iblock+1)//".txt", form='formatted')
       !write(101, *) compute%nv
       !do i = 1, compute%nv
       !   write(101, *) compute%v_index(i)
       !end do
       !close(101)
    else
      if(iblock .eq. 0) then
        open(101, file=trim(pardir)//"compute-"//i2s(nprocs), form='unformatted', access="stream")
        allocate(counts(nprocs))
        read(101) counts 
        sums = sum(counts)
        allocate(buffer_all(sums))
        read(101) buffer_all 
        close(101)
        allocate(displs(nprocs))
        displs(1) = 0
        do i = 2, nprocs
          displs(i) = displs(i-1) + counts(i-1)
        end do
      else
        allocate(counts(1))
        allocate(buffer_all(1))
        allocate(displs(1))
      end if
      call MPI_Scatter(counts, 1, MPI_INTEGER, compute%nv, &
                       1, MPI_INTEGER, 0, local_block%comm, ierr)
      allocate(compute%v_index(compute%nv))
      call MPI_Scatterv(buffer_all, counts, displs, MPI_INTEGER, compute%v_index, &
                        compute%nv, MPI_INTEGER, 0, local_block%comm, ierr)
      deallocate(counts)
      deallocate(displs)
      deallocate(buffer_all)
    end if 

  end subroutine search_compute


  subroutine geohash_encode(zlat, zlon, precs, geo_hash)
    implicit none
    ! io
    real(r8),         intent(in)  :: zlat, zlon
    integer(i4),      intent(in)  :: precs
    character(len=1), intent(out) :: geo_hash(1:precs)
    ! local
    integer                       :: bit, i, num
    logical                       :: is_odd
    real(r8)                      :: lats(1:2),lons(1:2)
    real(r8)                      :: mid
    integer                       :: tmp(1:5)
    character(len=1), parameter   :: BASE_32(0:31) = (/'0','1','2','3','4','5','6','7','8','9', &
                                       'b','c','d','e','f','g','h','j','k','m','n',&
                                       'p','q','r','s','t','u','v','w','x','y','z'/)

    is_odd = .false.

    bit = 1
    tmp = 0
    i = 1

    lats(1) = -90.d0
    lats(2) = 90.d0
    lons(1) = -180.d0
    lons(2) = 180.d0
    do while (i .le. precs)

       if (is_odd) then
          mid = (lats(1) + lats(2))/ 2;
          if (zlat .gt. mid) then
             tmp(bit) = 1
             lats(1) = mid
          else
             lats(2) = mid
          endif
       else
          mid = (lons(1) + lons(2))/ 2
          if (zlon .gt. mid) then
             tmp(bit) = 1
             lons(1) = mid
          else
             lons(2) = mid
          endif
       endif

       is_odd = .not. is_odd
       if (bit .lt. 5) then
          bit = bit + 1
       else
          num = tmp(1) * 16 + tmp(2) * 8 + tmp(3) * 4 + tmp(4) * 2 + tmp(5)
          geo_hash(i)= BASE_32(num)
          i = i + 1
          bit = 1
          tmp = 0
       endif
    enddo
  end subroutine geohash_encode


  subroutine xy2d(square, i, j, hilbert_pos)
    integer(i4), intent(in)    :: square
    integer,     intent(in)    :: i, j
    integer(i4), intent(inout) :: hilbert_pos
    integer                    :: rx, ry, s, iv, jv

    iv = i
    jv = j
    rx = 0
    ry = 0
    s = square/2
    do while(s > 0)
        if(iand(iv, s).gt.0) then
            rx = 1
        else
            rx = 0
        endif
        if(iand(jv, s).gt.0) then
            ry = 1
        else
            ry = 0
        endif

        hilbert_pos = hilbert_pos + s * s * (ieor(3*rx, ry))
        call rot(s, iv, jv, rx, ry)
        s = s/2
    end do
  end subroutine xy2d

  subroutine rot(s, x, y, rx, ry)
    integer, intent(in)    :: s, rx, ry
    integer, intent(inout) :: x, y
    integer                :: t

    if(ry == 0) then
      if(rx == 1) then
        x = s - 1 - x
        y = s - 1 - y
      endif

      t = x
      x = y
      y = t
    endif
  end subroutine rot

  !
  ! search bdy1, halo2: 
  ! for each cell in the compute_domain, if this cell's nb lies
  ! outside this domain, this cell is bdy1, this nb is halo2
  !

  subroutine search_bdry_halo(mesh, local_block)
    implicit none
    ! io
    type(global_domain_data),      intent(inout) :: mesh
    type(block_structure), target, intent(inout) :: local_block

    ! local
    integer(i4)                                  :: iv_domain
    integer(i4)                                  :: iv_global
    integer(i4)                                  :: iv_local, it_local
    integer(i4)                                  :: cnt, idx, sw, sw2, mark
    type(global_domain),   pointer               :: compute, b1, b2, h1, h2, halo, bdry
    type(global_domain),   pointer               :: i_domain, dm, o_domain, full, iner
    type(global_domain),   pointer               :: bdry_full
    logical                                      :: flag
    integer                                      :: n_max, n_min, n_avg, ierr
    character(20)                                :: str_tmp
    integer,               allocatable           :: counts(:), displs(:), ele3(:), &
                                                    ele3_all(:), buffer_all(:), buffer(:)
    integer                                      :: iblock, nprocs, nn, pos, i, sums
    integer,               allocatable           :: reduce_parts(:)

    iblock = local_block%cpuid
    nprocs = local_block%nprocs
    if ( .not. read_partition) then
       compute => local_block%compute_domain

       ! search based on v index
       b1 => local_block%bdry(1)
       h1 => local_block%halo(1)

       !find halo 1 and bdry 1
       do iv_domain = 1, compute%nv
          iv_global  = compute%v_index(iv_domain)
          do iv_local = 1, mesh%vtx_nnb(iv_global)

             idx = mesh%vtx_nb%f(iv_global, iv_local)
#ifdef LAM_DOMAIN
             ! determine partitions for cells in the specific region (halo1 first)
             mark = mesh%plg_mark%f(idx, 1)
             if (mark .eq. 177) then
                mesh%parts(idx) = min(mesh%parts(idx), iblock)
             end if
#endif
             flag = compute%v_index_list%find(idx)

             ! if not in, then this cell belongs to body1, this nb belongs halo2          
             if(.not.flag)then 
                call b1%v_index_list%insert(iv_global)
                call h1%v_index_list%insert(idx)
             end if
          end do
       end do

       call b1%v_index_list%dump(b1%v_index)
       call h1%v_index_list%dump(h1%v_index)

#ifdef LAM_DOMAIN
       ! reduce parts to get the correct partition for cells in the halo1 specific region
       allocate(reduce_parts(mesh%nv))
       call MPI_Allreduce(mesh%parts, reduce_parts, mesh%nv, &
                          MPI_INTEGER, MPI_MIN, local_block%comm, ierr)
       mesh%parts = reduce_parts
#endif
       b1%nv = size(b1%v_index)
       h1%nv = size(h1%v_index)

       ! find bdry(sw) for sw > 1
       do sw = 2, stencil_width
          i_domain => local_block%bdry(sw)
          dm       => local_block%bdry(sw-1)

          if(sw == 2) then
             o_domain => h1
          else
             o_domain => local_block%bdry(sw-2)
          end if

          do iv_domain = 1, dm%nv
             iv_global = dm%v_index(iv_domain)
             do iv_local = 1, mesh%vtx_nnb(iv_global)
                idx = mesh%vtx_nb%f(iv_global, iv_local)
                flag = (.not. dm%v_index_list%find(idx)) .and. &
                     (.not. o_domain%v_index_list%find(idx))

                if(flag) then
                   call i_domain%v_index_list%insert(idx)
                end if
             end do
          end do
          call i_domain%v_index_list%dump(i_domain%v_index)
          i_domain%nv = size(i_domain%v_index)
       end do

       ! find halo(sw) for sw > 1    
       do sw = 2, stencil_width
          i_domain => local_block%halo(sw)
          dm   => local_block%halo(sw-1)

          if(sw == 2) then
             o_domain => b1
          else
             o_domain => local_block%halo(sw-2)
          end if

          do iv_domain = 1, dm%nv
             iv_global = dm%v_index(iv_domain)
             do iv_local = 1, mesh%vtx_nnb(iv_global)
                idx = mesh%vtx_nb%f(iv_global, iv_local)
#ifdef LAM_DOMAIN
             if (idx .gt. 0) then
                ! determine partitions for cells in the specific region (halo2+)
                mark = mesh%plg_mark%f(idx, 1)
                if (mark .eq. 177 .and. mesh%parts(idx) .eq. maxi4) then
                   mesh%parts(idx) = min(mesh%parts(idx), iblock)
                end if
#endif
                flag = (.not. dm%v_index_list%find(idx)) .and. &
                     (.not. o_domain%v_index_list%find(idx))

                if(flag) then
                   call i_domain%v_index_list%insert(idx)
                end if
#ifdef LAM_DOMAIN
             end if
#endif
             end do
          end do
          call i_domain%v_index_list%dump(i_domain%v_index)
          i_domain%nv = size(i_domain%v_index)
       end do

       !
       ! halo 
       ! acquire t and e based global index
       !
       do sw = 1, stencil_width

          bdry => local_block%bdry(sw)       
          do iv_domain = 1, bdry%nv
             iv_global = bdry%v_index(iv_domain)
             do it_local = 1, mesh%vtx_nnb(iv_global)
                call bdry%t_index_list%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
                call bdry%e_index_list%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
             end do
          end do
          call bdry%t_index_list%dump(bdry%t_index)
          call bdry%e_index_list%dump(bdry%e_index)
          bdry%nt = size(bdry%t_index)
          bdry%ne = size(bdry%e_index)

          halo => local_block%halo(sw)       
          do iv_domain = 1, halo%nv
             iv_global = halo%v_index(iv_domain)
             do it_local = 1, mesh%vtx_nnb(iv_global)
                call halo%t_index_list%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
                call halo%e_index_list%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
             end do
          end do
          call halo%t_index_list%dump(halo%t_index)
          call halo%e_index_list%dump(halo%e_index)
          halo%nt = size(halo%t_index)
          halo%ne = size(halo%e_index)
          
       end do

       !!=====================================
       !! add bdry(1), bdry(2) to bdry domain
       !!=====================================
       !bdry_full => local_block%bdry_domain
       !do sw = stencil_width, 1, -1
       !   bdry => local_block%bdry(sw)
       !   do iv_domain = 1, bdry%nv
       !      iv_global = bdry%v_index(iv_domain)
       !      call bdry_full%v_index_list%insert(iv_global)
       !      do it_local = 1, mesh%vtx_nnb(iv_global)
       !         call bdry_full%t_index_list%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
       !         call bdry_full%e_index_list%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
       !      end do
       !   end do
       !end do

       !call bdry_full%v_index_list%dump(bdry_full%v_index)
       !call bdry_full%t_index_list%dump(bdry_full%t_index)
       !call bdry_full%e_index_list%dump(bdry_full%e_index)

       !bdry_full%nv = size(bdry_full%v_index)
       !bdry_full%nt = size(bdry_full%t_index)
       !bdry_full%ne = size(bdry_full%e_index)
       !call local_block%bdry_domain%v_index_list%final
       !call local_block%bdry_domain%t_index_list%final
       !call local_block%bdry_domain%e_index_list%final

#ifdef LAM_DOMAIN
       ! reduce parts to get the correct partition for cells in the halo2+ specific region
       call MPI_Allreduce(mesh%parts, reduce_parts, mesh%nv, &
                          MPI_INTEGER, MPI_MIN, local_block%comm, ierr)
       mesh%parts = reduce_parts
       deallocate(reduce_parts)
#endif

       !=================================
       ! configure the inner domain
       !=================================
#ifdef USE_INER
       iner => local_block%iner_domain
       do iv_domain = 1, compute%nv
          iv_global = compute%v_index(iv_domain)

          flag = .false.
          do sw = 1, stencil_width
             flag = flag .or. local_block%bdry(sw)%v_index_list%find(iv_global)
          end do
          if(.not.flag) then
             call iner%v_index_list%insert(iv_global)
          end if

       end do

       call iner%v_index_list%dump(iner%v_index)
       iner%nv = size(iner%v_index)

       do iv_domain = 1, iner%nv
          iv_global = iner%v_index(iv_domain)
          do it_local = 1, mesh%vtx_nnb(iv_global)
             call iner%t_index_list%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
             call iner%e_index_list%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
          end do
       end do
       call iner%t_index_list%dump(iner%t_index)
       call iner%e_index_list%dump(iner%e_index)
       iner%nt = size(iner%t_index)
       iner%ne = size(iner%e_index)
       call iner%v_index_list%final
       call iner%t_index_list%final
       call iner%e_index_list%final
#endif
       do sw = 1, stencil_width
          call local_block%bdry(sw)%v_index_list%final
          !call local_block%bdry(sw)%t_index_list%final
          !call local_block%bdry(sw)%e_index_list%final
          call local_block%halo(sw)%v_index_list%final
          !call local_block%halo(sw)%t_index_list%final
          !call local_block%halo(sw)%e_index_list%final
       end do

       !=================================
       ! configure the full domain
       !=================================
       full => local_block%full_domain

#ifdef USE_INER
       full%nv_iner = iner%nv
       full%ne_iner = iner%ne
       full%nt_iner = iner%nt
#endif

       full%nv_compute = compute%nv
       full%nv_bdry(1) = full%nv_compute
       do sw = 2, stencil_width
          full%nv_bdry(sw) = full%nv_bdry(sw-1)-local_block%bdry(sw-1)%nv
       end do

       full%nv_halo(1) = full%nv_compute+local_block%halo(1)%nv
       do sw = 2, stencil_width
          full%nv_halo(sw) = full%nv_halo(sw-1)+local_block%halo(sw)%nv
       end do
       full%nv = full%nv_halo(stencil_width)
       allocate(full%v_index(full%nv))

       cnt = 1
#ifdef USE_INER
       ! add iner
       do iv_domain = 1, iner%nv
          full%v_index(cnt) = iner%v_index(iv_domain)
          cnt = cnt + 1
       end do
       ! add bdry
       do sw = stencil_width, 1, -1
          bdry => local_block%bdry(sw)
          do iv_domain = 1, bdry%nv
             ! call full%v_index_list%insert(bdry%v_index(iv_domain))
             full%v_index(cnt) = bdry%v_index(iv_domain)
             cnt = cnt + 1
          end do
       end do
#else
       do iv_domain = 1, compute%nv
          full%v_index(cnt) = compute%v_index(iv_domain)
          cnt = cnt + 1
       end do
#endif
       ! add halo
       do sw = 1, stencil_width
          halo => local_block%halo(sw)
          do iv_domain = 1, halo%nv
             ! call full%v_index_list%insert(halo%v_index(iv_domain))
             full%v_index(cnt) = halo%v_index(iv_domain)
             cnt = cnt + 1
          end do
       end do

       ! call full%v_index_list%dump(full%v_index)
       ! full%nv = size(full%v_index)

       ! add t and e
       sw = 1
       sw2=stencil_width
       do iv_domain = 1, full%nv
          iv_global = full%v_index(iv_domain)

          do it_local = 1, mesh%vtx_nnb(iv_global)
             call full%t_index_list%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
             call full%e_index_list%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
          end do

          if(iv_domain .eq. full%nv_bdry(sw2))then
            full%nt_bdry(sw2) = full%t_index_list%size()
            full%ne_bdry(sw2) = full%e_index_list%size()
            if(sw2 .gt. 1) sw2=sw2-1
          end if

          if(iv_domain == full%nv_compute) then
             full%nt_compute = full%t_index_list%size()
             full%ne_compute = full%e_index_list%size()          
          end if

          if(iv_domain .eq. full%nv_halo(sw))then
             full%nt_halo(sw) = full%t_index_list%size()
             full%ne_halo(sw) = full%e_index_list%size()
             sw=sw+1
          end if
       end do

       call full%t_index_list%dump(full%t_index)
       call full%e_index_list%dump(full%e_index)
       full%nt = size(full%t_index)
       full%ne = size(full%e_index)
       !call full%v_index_list%final
       call full%t_index_list%final
       call full%e_index_list%final

       call compute%v_index_list%final
       !! test
       !open(101, file="full-v-"//i2s(iblock+1)//".txt", form='formatted')
       !write(101, *) full%nv
       !write(101, *) full%nv_compute
       !do sw = 1, stencil_width
       !   write(101, *) full%nv_bdry(sw)
       !end do
       !do sw = 1, stencil_width
       !   write(101, *) full%nv_halo(sw)
       !end do
       !do i = 1, full%nv
       !   write(101, *) full%v_index(i)
       !end do
       !close(101)
       !open(101, file="full-e-"//i2s(iblock+1)//".txt", form='formatted')
       !write(101, *) full%ne
       !write(101, *) full%ne_compute
       !do sw = 1, stencil_width
       !   write(101, *) full%ne_bdry(sw)
       !end do
       !do sw = 1, stencil_width
       !   write(101, *) full%ne_halo(sw)
       !end do
       !do i = 1, full%ne
       !   write(101, *) full%e_index(i)
       !end do
       !close(101)
       !open(101, file="full-t-"//i2s(iblock+1)//".txt", form='formatted')
       !write(101, *) full%nt
       !write(101, *) full%nt_compute
       !do sw = 1, stencil_width
       !   write(101, *) full%nt_bdry(sw)
       !end do
       !do sw = 1, stencil_width
       !   write(101, *) full%nt_halo(sw)
       !end do
       !do i = 1, full%nt
       !   write(101, *) full%t_index(i)
       !end do
       !close(101)
       !!open(101, file="parts-"//i2s(iblock+1)//".txt", form='formatted')
       !if (iblock .eq. 0) then
       !open(101, file="parts.txt", form='formatted')
       !do i = 1, mesh%nv
       !   write(101, *) mesh%parts(i)
       !end do
       !close(101)
       !end if
    else
      full => local_block%full_domain
      ! read full (read by PE 0, then scatter to all PEs)
      if(iblock .eq. 0) then
        open(101, file=trim(pardir)//"full-"//i2s(nprocs), form='unformatted', access="stream")
        nn = (6 + 6*stencil_width)*nprocs
#ifdef USE_INER
        nn = nn + 3*nprocs
#endif
        allocate(ele3_all(nn))
        read(101) ele3_all
        allocate(counts(nprocs))
        pos = 1
        do i = 1, nprocs
          counts(i) = ele3_all(pos) + ele3_all(pos+2*stencil_width+2) + &
                      ele3_all(pos+4*stencil_width+4)
          pos = pos + 6*stencil_width + 6
#ifdef USE_INER
          pos = pos + 3
#endif
        end do
        sums = sum(counts)
        allocate(buffer_all(sums))
        read(101) buffer_all
        close(101)

        nn = nn/nprocs
        allocate(displs(nprocs))
        displs(1) = 0
        do i = 2, nprocs
          displs(i) = displs(i-1) + counts(i-1)
        end do
      else
        allocate(ele3_all(1))
        allocate(counts(1))
        allocate(buffer_all(1))
        allocate(displs(1))
        nn = 6 + 6*stencil_width
#ifdef USE_INER
        nn = nn + 3
#endif
      end if
      allocate(ele3(nn))
      call MPI_Scatter(ele3_all, nn, MPI_INTEGER, ele3, &
                       nn, MPI_INTEGER, 0, local_block%comm, ierr)
      sums = ele3(1) + ele3(1+2*stencil_width+2) + ele3(1+4*stencil_width+4)
      allocate(buffer(sums))
      call MPI_Scatterv(buffer_all, counts, displs, MPI_INTEGER, buffer, &
                        sums, MPI_INTEGER, 0, local_block%comm, ierr)
      pos = 1
      full%nv = ele3(pos)
      pos = pos + 1
      full%nv_compute = ele3(pos)
      pos = pos + 1
      do sw = 1, stencil_width
        full%nv_bdry(sw) = ele3(pos)
        pos = pos + 1
      end do
      do sw = 1, stencil_width
        full%nv_halo(sw) = ele3(pos)
        pos = pos + 1
      end do
      full%nt = ele3(pos)
      pos = pos + 1
      full%nt_compute = ele3(pos)
      pos = pos + 1
      do sw = 1, stencil_width
        full%nt_bdry(sw) = ele3(pos)
        pos = pos + 1
      end do
      do sw = 1, stencil_width
        full%nt_halo(sw) = ele3(pos)
        pos = pos + 1
      end do
      full%ne = ele3(pos)
      pos = pos + 1
      full%ne_compute = ele3(pos)
      pos = pos + 1
      do sw = 1, stencil_width
        full%ne_bdry(sw) = ele3(pos)
        pos = pos + 1
      end do
      do sw = 1, stencil_width
        full%ne_halo(sw) = ele3(pos)
        pos = pos + 1
      end do
#ifdef USE_INER
      full%nv_iner = ele3(pos)
      pos = pos + 1
      full%nt_iner = ele3(pos)
      pos = pos + 1
      full%ne_iner = ele3(pos)
#endif
      allocate(full%v_index(full%nv))
      allocate(full%t_index(full%nt))
      allocate(full%e_index(full%ne))
      full%v_index = buffer(1:full%nv)
      full%t_index = buffer(full%nv+1:full%nv+full%nt)
      full%e_index = buffer(full%nv+full%nt+1:full%nv+full%nt+full%ne)
      deallocate(ele3_all)
      deallocate(counts)
      deallocate(buffer_all)
      deallocate(displs)
      deallocate(ele3)
      deallocate(buffer)
    end if

    full%nv_full = full%nv    
    full%nt_full = full%nt
    full%ne_full = full%ne
!
! asign values
! only mean_edt_dist is needed by MODEL
!
    full%min_edt_dist  = mesh%min_edt_dist
    full%max_edt_dist  = mesh%max_edt_dist
    full%mean_edt_dist = mesh%mean_edt_dist

    full%min_edp_dist  = mesh%min_edp_dist
    full%max_edp_dist  = mesh%max_edp_dist
    full%mean_edp_dist = mesh%mean_edp_dist

    full%min_tri_area  = mesh%min_tri_area
    full%max_tri_area  = mesh%max_tri_area
    full%mean_tri_area = mesh%mean_tri_area

    full%min_plg_area  = mesh%min_plg_area
    full%max_plg_area  = mesh%max_plg_area
    full%mean_plg_area = mesh%mean_plg_area

    full%min_tri_angle = mesh%min_tri_angle
    full%max_tri_angle = mesh%max_tri_angle
    full%mean_tri_angle= mesh%mean_tri_angle

    call reduce(full%nv_compute, n_min, 'min', local_block%comm)
    call reduce(full%nv_compute, n_max, 'max', local_block%comm)
    call reduce(full%nv_compute, n_avg, 'sum', local_block%comm)
    if(iblock .eq. 0) then
        write(*,'(a,i6,a)') ' Partition information for ', nprocs, ' PEs...'
        write(str_tmp,'(a)') '#cells compute'
        write(*, '(20x,a)') '        min       max       avg'
        write(*,'(x,a20,3i10)') str_tmp, n_min, n_max, n_avg/nprocs
    end if
    call reduce(full%nv_full, n_min, 'min', local_block%comm)
    call reduce(full%nv_full, n_max, 'max', local_block%comm)
    call reduce(full%nv_full, n_avg, 'sum', local_block%comm)
    if(iblock .eq. 0) then
        write(str_tmp,'(a)') '#cells with halo'
        write(*,'(x,a20,3i10)') str_tmp, n_min, n_max, n_avg/nprocs
    end if
    call reduce(full%nt_compute, n_min, 'min', local_block%comm)
    call reduce(full%nt_compute, n_max, 'max', local_block%comm)
    call reduce(full%nt_compute, n_avg, 'sum', local_block%comm)
    if(iblock .eq. 0) then
        write(str_tmp,'(a)') '#verts compute'
        write(*,'(x,a20,3i10)') str_tmp, n_min, n_max, n_avg/nprocs
    end if
    call reduce(full%nt_full, n_min, 'min', local_block%comm)
    call reduce(full%nt_full, n_max, 'max', local_block%comm)
    call reduce(full%nt_full, n_avg, 'sum', local_block%comm)
    if(iblock .eq. 0) then
        write(str_tmp,'(a)') '#verts with halo'
        write(*,'(x,a20,3i10)') str_tmp, n_min, n_max, n_avg/nprocs
    end if
    call reduce(full%ne_compute, n_min, 'min', local_block%comm)
    call reduce(full%ne_compute, n_max, 'max', local_block%comm)
    call reduce(full%ne_compute, n_avg, 'sum', local_block%comm)
    if(iblock .eq. 0) then
        write(str_tmp,'(a)') '#edges compute'
        write(*,'(x,a20,3i10)') str_tmp, n_min, n_max, n_avg/nprocs
    end if
    call reduce(full%ne_full, n_min, 'min', local_block%comm)
    call reduce(full%ne_full, n_max, 'max', local_block%comm)
    call reduce(full%ne_full, n_avg, 'sum', local_block%comm)
    if(iblock .eq. 0) then
        write(str_tmp,'(a)') '#edges with halo'
        write(*,'(x,a20,3i10)') str_tmp, n_min, n_max, n_avg/nprocs
        print*,"===================================================="
    end if
  end subroutine search_bdry_halo


  subroutine pull_element_types_full(local_block)
    ! io
    type(block_structure), target, intent(inout) :: local_block

    ! local
    integer(i4)                                  :: sw

    call pull_element_types_domain(local_block%full_domain, local_block%comm)

    !call pull_element_types_domain(local_block%compute_domain)
    !
    !call pull_element_types_domain(local_block%iner_domain)

    !call pull_element_types_domain(local_block%bdry_domain)
 
    !do sw = 1, stencil_width
    !   call pull_element_types_domain(local_block%bdry(sw))
    !   call pull_element_types_domain(local_block%halo(sw))
    !   
    !   if(check)then
    !     open(1, file=i2s(local_block%cpuid)//".bdry"//i2s(sw)//".v.txt")
    !     write(1, "(I0)", advance='yes') local_block%bdry(sw)%v_index
    !     close(1)

    !     open(1, file=i2s(local_block%cpuid)//".halo"//i2s(sw)//".v.txt")
    !     write(1, "(I0)", advance='yes') local_block%halo(sw)%v_index
    !     close(1)
    !   end if
    !end do

  end subroutine pull_element_types_full


  subroutine pull_element_types_domain(dm, comm)
    ! io
    type(global_domain),   intent(inout), target :: dm
    integer,               intent(in)            :: comm
    ! local

    call grist_input_init(dm, comm)
    call grist_read_grid_file_par(dm)

    return
  end subroutine pull_element_types_domain


  subroutine grist_input_init(dm, comm)
    implicit none
    ! io
    type(global_domain), target, intent(inout) :: dm
    integer,                     intent(in)    :: comm
    ! local
    integer(i4)                                :: i, ig, np_input 
    integer(i4)                                :: nv_avg, nv_res, vpos, &
                                                  ne_avg, ne_res, epos, &
                                                  nt_avg, nt_res, tpos
    type(set)                                  :: sv_list, se_list, st_list
    integer, allocatable                       :: buff_svi(:), buff_sti(:), buff_sei(:)

    type(group_comm), pointer                  :: grcomm
    integer                                    :: ierr, nprocs
    type(map)                                  :: l2l_v, l2l_t, l2l_e
    integer, allocatable                       :: sendvi(:), sendti(:), sendei(:)
    integer(i4)                                :: svtot, sttot, setot

    grcomm => dm%gcomm_read
    grcomm%all_comm = comm
    grcomm%gsize = comm_group_size
    grcomm%color = mod(mpi_rank(comm),grcomm%gsize)
    call mpi_comm_split(comm,grcomm%color,mpi_rank(comm),grcomm%gcomm,ierr)
    call mpi_comm_rank(grcomm%gcomm,grcomm%grank,ierr)

    nprocs = mpi_size(comm)
    np_input = ceiling(nprocs/real(comm_group_size, 8))

    ! for vtx
    ! sort v_index & get localv
    do i = 1, size(dm%v_index)
        call sv_list%insert(dm%v_index(i))
    enddo
    allocate(buff_svi(sv_list%size()))
    call sv_list%sort(buff_svi)
    call sv_list%final
    do i = 1, size(buff_svi)
        call l2l_v%insert(buff_svi(i), i)
    enddo
    allocate(grcomm%idx_localv(size(buff_svi)))
    do i = 1, size(buff_svi)
       grcomm%idx_localv(i) = l2l_v%find(dm%v_index(i))
    end do
    call l2l_v%final
    ! recvv counts & displs
    allocate(grcomm%counts_recvv(1:nprocs), source=0)
    allocate(grcomm%displs_recvv(1:nprocs), source=0)
    nv_avg = dm%nv_all/np_input
    nv_res = mod(dm%nv_all, np_input)
    vpos = (nv_avg+1)*nv_res
    do i = 1, size(buff_svi)
        if(buff_svi(i) .gt. vpos) then
            ig = (buff_svi(i)-vpos-1)/nv_avg + nv_res
        else
            ig = (buff_svi(i)-1)/(nv_avg+1)
        end if
        grcomm%counts_recvv(ig*comm_group_size+1) = grcomm%counts_recvv(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        grcomm%displs_recvv(i) = grcomm%displs_recvv(i-1) + grcomm%counts_recvv(i-1)
    end do
    ! sendv counts & displs
    allocate(grcomm%counts_sendv(nprocs), source=0)
    allocate(grcomm%displs_sendv(nprocs), source=0)
    call MPI_Alltoall(grcomm%counts_recvv, 1, MPI_INTEGER, grcomm%counts_sendv, 1, MPI_INTEGER, comm, ierr)
    if(grcomm%color .eq. 0) then
        svtot = grcomm%counts_sendv(1)
        do i = 2, nprocs
            grcomm%displs_sendv(i) = grcomm%displs_sendv(i-1) + grcomm%counts_sendv(i-1)
            svtot = svtot + grcomm%counts_sendv(i)
        end do
    else
        svtot = 1
    end if
    allocate(sendvi(svtot))
    call MPI_Alltoallv(buff_svi, grcomm%counts_recvv, grcomm%displs_recvv, MPI_INTEGER, sendvi, &
                       grcomm%counts_sendv, grcomm%displs_sendv, MPI_INTEGER, comm, ierr)
    ! sendv local idx
    allocate(grcomm%idx_sendv(svtot))
    if(grcomm%color .eq. 0) then
        if(grcomm%grank .le. nv_res) then
            grcomm%mystart_v = grcomm%grank*(nv_avg+1) + 1
        else
            grcomm%mystart_v = vpos + (grcomm%grank-nv_res)*nv_avg + 1
        end if
        if(grcomm%grank .lt. nv_res) then
            grcomm%mycount_v = nv_avg + 1
        else
            grcomm%mycount_v = nv_avg
        end if
        do i = 1, svtot
            grcomm%idx_sendv(i) = sendvi(i)-grcomm%mystart_v+1
        end do
    end if
    deallocate(buff_svi)
    deallocate(sendvi)

    ! for edg
    ! sort e_index & get locale
    do i = 1, size(dm%e_index)
        call se_list%insert(dm%e_index(i))
    enddo
    allocate(buff_sei(se_list%size()))
    call se_list%sort(buff_sei)
    call se_list%final
    do i = 1, size(buff_sei)
        call l2l_e%insert(buff_sei(i), i)
    enddo
    allocate(grcomm%idx_locale(size(buff_sei)))
    do i = 1, size(buff_sei)
       grcomm%idx_locale(i) = l2l_e%find(dm%e_index(i))
    end do
    call l2l_e%final
    ! recve counts & displs
    allocate(grcomm%counts_recve(1:nprocs), source=0)
    allocate(grcomm%displs_recve(1:nprocs), source=0)
    ne_avg = dm%ne_all/np_input
    ne_res = mod(dm%ne_all, np_input)
    epos = (ne_avg+1)*ne_res
    do i = 1, size(buff_sei)
        if(buff_sei(i) .gt. epos) then
            ig = (buff_sei(i)-epos-1)/ne_avg + ne_res
        else
            ig = (buff_sei(i)-1)/(ne_avg+1)
        end if
        grcomm%counts_recve(ig*comm_group_size+1) = grcomm%counts_recve(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        grcomm%displs_recve(i) = grcomm%displs_recve(i-1) + grcomm%counts_recve(i-1)
    end do
    ! sende counts & displs
    allocate(grcomm%counts_sende(nprocs), source=0)
    allocate(grcomm%displs_sende(nprocs), source=0)
    call MPI_Alltoall(grcomm%counts_recve, 1, MPI_INTEGER, grcomm%counts_sende, 1, MPI_INTEGER, comm, ierr)
    if(grcomm%color .eq. 0) then
        setot = grcomm%counts_sende(1)
        do i = 2, nprocs
            grcomm%displs_sende(i) = grcomm%displs_sende(i-1) + grcomm%counts_sende(i-1)
            setot = setot + grcomm%counts_sende(i)
        end do
    else
        setot = 1
    end if
    allocate(sendei(setot))
    call MPI_Alltoallv(buff_sei, grcomm%counts_recve, grcomm%displs_recve, MPI_INTEGER, sendei, &
                       grcomm%counts_sende, grcomm%displs_sende, MPI_INTEGER, comm, ierr)
    ! sende local idx
    allocate(grcomm%idx_sende(setot))
    if(grcomm%color .eq. 0) then
        if(grcomm%grank .le. ne_res) then
            grcomm%mystart_e = grcomm%grank*(ne_avg+1) + 1
        else
            grcomm%mystart_e = epos + (grcomm%grank-ne_res)*ne_avg + 1
        end if
        if(grcomm%grank .lt. ne_res) then
            grcomm%mycount_e = ne_avg + 1
        else
            grcomm%mycount_e = ne_avg
        end if
        do i = 1, setot
            grcomm%idx_sende(i) = sendei(i)-grcomm%mystart_e+1
        end do
    end if
    deallocate(buff_sei)
    deallocate(sendei)

    ! for tri
    ! sort t_index & get localt
    do i = 1, size(dm%t_index)
        call st_list%insert(dm%t_index(i))
    enddo
    allocate(buff_sti(st_list%size()))
    call st_list%sort(buff_sti)
    call st_list%final
    do i = 1, size(buff_sti)
        call l2l_t%insert(buff_sti(i), i)
    enddo
    allocate(grcomm%idx_localt(size(buff_sti)))
    do i = 1, size(buff_sti)
       grcomm%idx_localt(i) = l2l_t%find(dm%t_index(i))
    end do
    call l2l_t%final
    ! recvt counts & displs
    allocate(grcomm%counts_recvt(1:nprocs), source=0)
    allocate(grcomm%displs_recvt(1:nprocs), source=0)
    nt_avg = dm%nt_all/np_input
    nt_res = mod(dm%nt_all, np_input)
    tpos = (nt_avg+1)*nt_res
    do i = 1, size(buff_sti)
        if(buff_sti(i) .gt. tpos) then
            ig = (buff_sti(i)-tpos-1)/nt_avg + nt_res
        else
            ig = (buff_sti(i)-1)/(nt_avg+1)
        end if
        grcomm%counts_recvt(ig*comm_group_size+1) = grcomm%counts_recvt(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        grcomm%displs_recvt(i) = grcomm%displs_recvt(i-1) + grcomm%counts_recvt(i-1)
    end do
    ! sendt counts & displs
    allocate(grcomm%counts_sendt(nprocs), source=0)
    allocate(grcomm%displs_sendt(nprocs), source=0)
    call MPI_Alltoall(grcomm%counts_recvt, 1, MPI_INTEGER, grcomm%counts_sendt, 1, MPI_INTEGER, comm, ierr)
    if(grcomm%color .eq. 0) then
        sttot = grcomm%counts_sendt(1)
        do i = 2, nprocs
            grcomm%displs_sendt(i) = grcomm%displs_sendt(i-1) + grcomm%counts_sendt(i-1)
            sttot = sttot + grcomm%counts_sendt(i)
        end do
    else
        sttot = 1
    end if
    allocate(sendti(sttot))
    call MPI_Alltoallv(buff_sti, grcomm%counts_recvt, grcomm%displs_recvt, MPI_INTEGER, sendti, &
                       grcomm%counts_sendt, grcomm%displs_sendt, MPI_INTEGER, comm, ierr)
    ! sendt local idx
    allocate(grcomm%idx_sendt(sttot))
    if(grcomm%color .eq. 0) then
        if(grcomm%grank .le. nt_res) then
            grcomm%mystart_t = grcomm%grank*(nt_avg+1) + 1
        else
            grcomm%mystart_t = tpos + (grcomm%grank-nt_res)*nt_avg + 1
        end if
        if(grcomm%grank .lt. nt_res) then
            grcomm%mycount_t = nt_avg + 1
        else
            grcomm%mycount_t = nt_avg
        end if
        do i = 1, sttot
            grcomm%idx_sendt(i) = sendti(i)-grcomm%mystart_t+1
        end do
    end if
    deallocate(buff_sti)
    deallocate(sendti)

    return
  end subroutine grist_input_init


  subroutine grist_input_final(dm)
    implicit none
    ! io
    type(global_domain), target, intent(inout) :: dm
    ! local
    type(group_comm),    pointer               :: comm

    comm => dm%gcomm_read
    deallocate(comm%idx_localv)
    deallocate(comm%idx_localt)
    deallocate(comm%idx_locale)
    deallocate(comm%counts_recvv)
    deallocate(comm%counts_recvt)
    deallocate(comm%counts_recve)
    deallocate(comm%displs_recvv)
    deallocate(comm%displs_recvt)
    deallocate(comm%displs_recve)
    deallocate(comm%counts_sendv)
    deallocate(comm%counts_sendt)
    deallocate(comm%counts_sende)
    deallocate(comm%displs_sendv)
    deallocate(comm%displs_sendt)
    deallocate(comm%displs_sende)
    deallocate(comm%idx_sendv)
    deallocate(comm%idx_sendt)
    deallocate(comm%idx_sende)
  end subroutine grist_input_final


  subroutine grist_output_init(mesh)
    ! io
    type(global_domain),   intent(inout), target :: mesh
    ! local
    integer(i4)                                  :: i, ig, np_output
    integer(i4)                                  :: nv_avg, nv_res, vpos, &
                                                    ne_avg, ne_res, epos, &
                                                    nt_avg, nt_res, tpos
    type(set)                                    :: sv_list, se_list, st_list
    integer,               allocatable           :: buff_svi(:), buff_sti(:), buff_sei(:)

    type(block_structure), pointer               :: local_block
    type(group_comm),      pointer               :: comm
    integer                                      :: nprocs, ierr
    integer, allocatable                         :: recvvi(:), recvti(:), recvei(:)
    integer(i4)                                  :: rvtot, rttot, retot

    local_block => mesh%local_block
    comm => mesh%gcomm_write
    !comm%all_comm = local_block%comm
    comm%all_comm = grist_comm()
    comm%gsize = comm_group_size
    comm%color = mod(mpi_rank(comm%all_comm),comm%gsize)
    call mpi_comm_split(comm%all_comm, comm%color, mpi_rank(comm%all_comm), comm%gcomm,ierr)
    call mpi_comm_rank(comm%gcomm,comm%grank,ierr)

    nprocs = mpi_size(comm%all_comm)
    np_output = ceiling(nprocs/real(comm_group_size, 8))
    ! for vtx
    ! sort v_index & get localv
    do i = 1, size(local_block%output%v)
        call sv_list%insert(local_block%output%v(i))
    enddo
    allocate(buff_svi(sv_list%size()))
    call sv_list%sort(buff_svi)
    call sv_list%final
    allocate(comm%idx_localv(size(buff_svi)))
    do i = 1, size(buff_svi)
       comm%idx_localv(i) = mesh%map_g2l_v%find(buff_svi(i))
    end do
    ! sendv counts & displs
    allocate(comm%counts_sendv(1:nprocs), source=0)
    allocate(comm%displs_sendv(1:nprocs), source=0)
    nv_avg = mesh%nv_all/np_output
    nv_res = mod(mesh%nv_all, np_output)
    vpos = (nv_avg+1)*nv_res
    do i = 1, size(buff_svi)
        if(buff_svi(i) .gt. vpos) then
            ig = (buff_svi(i)-vpos-1)/nv_avg + nv_res
        else
            ig = (buff_svi(i)-1)/(nv_avg+1)
        end if
        comm%counts_sendv(ig*comm_group_size+1) = comm%counts_sendv(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sendv(i) = comm%displs_sendv(i-1) + comm%counts_sendv(i-1)
    end do
    ! recvv counts & displs
    allocate(comm%counts_recvv(nprocs), source=0)
    allocate(comm%displs_recvv(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sendv, 1, MPI_INTEGER, comm%counts_recvv, 1, MPI_INTEGER, grist_comm(), ierr)
    if(comm%color .eq. 0) then
        rvtot = comm%counts_recvv(1)
        do i = 2, nprocs
            comm%displs_recvv(i) = comm%displs_recvv(i-1) + comm%counts_recvv(i-1)
            rvtot = rvtot + comm%counts_recvv(i)
        end do
    else
        rvtot = 1
    end if
    allocate(recvvi(rvtot))
    call MPI_Alltoallv(buff_svi, comm%counts_sendv, comm%displs_sendv, MPI_INTEGER, recvvi, &
                       comm%counts_recvv, comm%displs_recvv, MPI_INTEGER, grist_comm(), ierr)
    !open(101, file="recvv"//i2s(mpi_rank()))
    !do i = 1, size(recvvi)
    !  write(101, *) recvvi(i)
    !end do
    !close(101)
    ! recvv local idx
    allocate(comm%idx_recvv(rvtot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. nv_res) then
            comm%mystart_v = comm%grank*(nv_avg+1) + 1
        else
            comm%mystart_v = vpos + (comm%grank-nv_res)*nv_avg + 1
        end if
        if(comm%grank .lt. nv_res) then
            comm%mycount_v = nv_avg + 1
        else
            comm%mycount_v = nv_avg
        end if
        do i = 1, rvtot
            comm%idx_recvv(i) = recvvi(i)-comm%mystart_v+1
        end do
    end if
    deallocate(buff_svi)
    deallocate(recvvi)

    ! for edg
    ! sort e_index & get locale
    do i = 1, size(local_block%output%e)
        call se_list%insert(local_block%output%e(i))
    enddo
    allocate(buff_sei(se_list%size()))
    call se_list%sort(buff_sei)
    call se_list%final
    allocate(comm%idx_locale(size(buff_sei)))
    do i = 1, size(buff_sei)
       comm%idx_locale(i) = mesh%map_g2l_e%find(buff_sei(i))
    end do
    ! sende counts & displs
    allocate(comm%counts_sende(1:nprocs), source=0)
    allocate(comm%displs_sende(1:nprocs), source=0)
    ne_avg = mesh%ne_all/np_output
    ne_res = mod(mesh%ne_all, np_output)
    epos = (ne_avg+1)*ne_res
    do i = 1, size(buff_sei)
        if(buff_sei(i) .gt. epos) then
            ig = (buff_sei(i)-epos-1)/ne_avg + ne_res
        else
            ig = (buff_sei(i)-1)/(ne_avg+1)
        end if
        comm%counts_sende(ig*comm_group_size+1) = comm%counts_sende(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sende(i) = comm%displs_sende(i-1) + comm%counts_sende(i-1)
    end do
    ! recve counts & displs
    allocate(comm%counts_recve(nprocs), source=0)
    allocate(comm%displs_recve(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sende, 1, MPI_INTEGER, comm%counts_recve, 1, MPI_INTEGER, grist_comm(), ierr)
    if(comm%color .eq. 0) then
        retot = comm%counts_recve(1)
        do i = 2, nprocs
            comm%displs_recve(i) = comm%displs_recve(i-1) + comm%counts_recve(i-1)
            retot = retot + comm%counts_recve(i)
        end do
    else
        retot = 1
    end if
    allocate(recvei(retot))
    call MPI_Alltoallv(buff_sei, comm%counts_sende, comm%displs_sende, MPI_INTEGER, recvei, &
                       comm%counts_recve, comm%displs_recve, MPI_INTEGER, grist_comm(), ierr)
    ! recve local idx
    allocate(comm%idx_recve(retot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. ne_res) then
            comm%mystart_e = comm%grank*(ne_avg+1) + 1
        else
            comm%mystart_e = epos + (comm%grank-ne_res)*ne_avg + 1
        end if
        if(comm%grank .lt. ne_res) then
            comm%mycount_e = ne_avg + 1
        else
            comm%mycount_e = ne_avg
        end if
        do i = 1, retot
            comm%idx_recve(i) = recvei(i)-comm%mystart_e+1
        end do
    end if
    deallocate(buff_sei)
    deallocate(recvei)

    ! for tri
    ! sort t_index & get localt
    do i = 1, size(local_block%output%t)
        call st_list%insert(local_block%output%t(i))
    enddo
    allocate(buff_sti(st_list%size()))
    call st_list%sort(buff_sti)
    call st_list%final
    allocate(comm%idx_localt(size(buff_sti)))
    do i = 1, size(buff_sti)
       comm%idx_localt(i) = mesh%map_g2l_t%find(buff_sti(i))
    end do
    ! sendt counts & displs
    allocate(comm%counts_sendt(1:nprocs), source=0)
    allocate(comm%displs_sendt(1:nprocs), source=0)
    nt_avg = mesh%nt_all/np_output
    nt_res = mod(mesh%nt_all, np_output)
    tpos = (nt_avg+1)*nt_res
    do i = 1, size(buff_sti)
        if(buff_sti(i) .gt. tpos) then
            ig = (buff_sti(i)-tpos-1)/nt_avg + nt_res
        else
            ig = (buff_sti(i)-1)/(nt_avg+1)
        end if
        comm%counts_sendt(ig*comm_group_size+1) = comm%counts_sendt(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sendt(i) = comm%displs_sendt(i-1) + comm%counts_sendt(i-1)
    end do
    ! recvt counts & displs
    allocate(comm%counts_recvt(nprocs), source=0)
    allocate(comm%displs_recvt(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sendt, 1, MPI_INTEGER, comm%counts_recvt, 1, MPI_INTEGER, grist_comm(), ierr)
    if(comm%color .eq. 0) then
        rttot = comm%counts_recvt(1)
        do i = 2, nprocs
            comm%displs_recvt(i) = comm%displs_recvt(i-1) + comm%counts_recvt(i-1)
            rttot = rttot + comm%counts_recvt(i)
        end do
    else
        rttot = 1
    end if
    allocate(recvti(rttot))
    call MPI_Alltoallv(buff_sti, comm%counts_sendt, comm%displs_sendt, MPI_INTEGER, recvti, &
                       comm%counts_recvt, comm%displs_recvt, MPI_INTEGER, grist_comm(), ierr)
    ! recvt local idx
    allocate(comm%idx_recvt(rttot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. nt_res) then
            comm%mystart_t = comm%grank*(nt_avg+1) + 1
        else
            comm%mystart_t = tpos + (comm%grank-nt_res)*nt_avg + 1
        end if
        if(comm%grank .lt. nt_res) then
            comm%mycount_t = nt_avg + 1
        else
            comm%mycount_t = nt_avg
        end if
        do i = 1, rttot
            comm%idx_recvt(i) = recvti(i)-comm%mystart_t+1
        end do
    end if
    deallocate(buff_sti)
    deallocate(recvti)

    return
  end subroutine grist_output_init

  subroutine grist_output_final(mesh)
    ! io
    type(global_domain), intent(inout), target :: mesh
    ! local
    type(group_comm),      pointer             :: comm
    type(block_structure), pointer             :: local_block

    local_block => mesh%local_block
    comm => mesh%gcomm_write

    deallocate(comm%idx_localv)
    deallocate(comm%idx_localt)
    deallocate(comm%idx_locale)
    deallocate(comm%counts_sendv)
    deallocate(comm%counts_sendt)
    deallocate(comm%counts_sende)
    deallocate(comm%displs_sendv)
    deallocate(comm%displs_sendt)
    deallocate(comm%displs_sende)
    deallocate(comm%counts_recvv)
    deallocate(comm%counts_recvt)
    deallocate(comm%counts_recve)
    deallocate(comm%displs_recvv)
    deallocate(comm%displs_recvt)
    deallocate(comm%displs_recve)
    deallocate(comm%idx_recvv)
    deallocate(comm%idx_recvt)
    deallocate(comm%idx_recve)

    return
  end subroutine grist_output_final


  subroutine search_exchange(mesh,local_block)
    implicit none
    ! io
    type(global_domain_data),      intent(in)    :: mesh
    type(block_structure), target, intent(inout) :: local_block

    ! local
    integer(i4)                                  :: iv_domain, iv_global, iv_local
    type(global_domain),   pointer               :: compute, halo
    integer                                      :: iblock, c, cpuid, mark, marke, markt
    integer(i4)                                  :: i, sw, nid

    type(set),             allocatable           :: recv_v(:), recv_e(:), recv_t(:)
    
    type(set)                                    :: recv_cpus_list
    integer,               allocatable           :: recv_cpus(:)

    integer,               allocatable           :: reqs_send(:), reqs_recv(:)
    integer,               allocatable           :: num_send_sites(:,:), num_recv_sites(:,:)
    integer                                      :: dst, ierr
    INTEGER,               allocatable           :: status(:,:) !STATUS(MPI_STATUS_SIZE)
    integer                                      :: ie, it, n, ne, nt, nv, pos, poe
    integer                                      :: nprocs, sums
    integer,               allocatable           :: counts(:), displs(:), nbcpus(:), ele6(:), &
                                                    ele6_all(:), buffer_all(:), buffer(:), &
                                                    counts_idx(:), displs_idx(:)
    logical                                      :: is_comm_nb

    iblock = local_block%cpuid
    nprocs = local_block%nprocs
    if (.not. read_partition) then
       compute => local_block%compute_domain

       allocate(recv_v(nprocs))
       allocate(recv_e(nprocs))
       allocate(recv_t(nprocs))

       do sw = 1, stencil_width
          halo => local_block%halo(sw)

          do iv_domain = 1, halo%nv 
             iv_global = halo%v_index(iv_domain)
#ifdef LAM_DOMAIN
             ! discard cells/edges/triangles in the specific region
             if (mesh%parts(iv_global) .eq. iblock) cycle
             is_comm_nb = .false.
             mark = mesh%plg_mark%f(iv_global, 1)
             if (mark .ne. 177) then
                is_comm_nb = .true.
             else
                do n = 1, mesh%vtx_nnb(iv_global)
                   ie = mesh%vtx_ed%f(iv_global, n)
                   it = mesh%vtx_tr%f(iv_global, n)
                   marke = mesh%edt_mark%f(ie, 1)
                   markt = mesh%tri_mark%f(it, 1)
                   if (marke/10 .ne. 177 .and. .not. local_block%bdry(1)%e_index_list%find(ie) &
                       .or. markt .ne. 177 .and. .not. local_block%bdry(1)%t_index_list%find(it)) then
                      is_comm_nb = .true.
                      exit
                   end if
                end do
             end if
             if (.not. is_comm_nb) cycle

             nid = mesh%parts(iv_global)
             if (mark .ne. 177) then
                 call recv_cpus_list%insert(nid)
                 call recv_v(nid+1)%insert(iv_global)
             end if
             do n = 1, mesh%vtx_nnb(iv_global)
                ie = mesh%vtx_ed%f(iv_global, n)
                it = mesh%vtx_tr%f(iv_global, n)
                marke = mesh%edt_mark%f(ie, 1)
                markt = mesh%tri_mark%f(it, 1)
                nid = mesh%parts(iv_global)
                if (marke/10 .ne. 177 .and. .not. local_block%bdry(1)%e_index_list%find(ie)) then
                    if (mark .ne. 177) then
                        call recv_e(nid+1)%insert(ie)
                    else
                        ! transfer from compute id
                        do i = 1, 2
                            iv_local = mesh%edt_v%f(ie, i)
                            mark = mesh%plg_mark%f(iv_local, 1)
                            if (mark .ne. 177) then
                                nid = mesh%parts(iv_local)
                                call recv_cpus_list%insert(nid)
                                call recv_e(nid+1)%insert(ie)
                                exit
                            end if
                        end do
                    end if
                end if
                nid = mesh%parts(iv_global)
                if (markt .ne. 177 .and. .not. local_block%bdry(1)%t_index_list%find(it)) then
                    if (mark .ne. 177) then
                        call recv_t(nid+1)%insert(it)
                    else
                        ! transfer from compute id
                        do i = 1, 3
                            iv_local = mesh%tri_v%f(it, i)
                            mark = mesh%plg_mark%f(iv_local, 1)
                            if (mark .ne. 177) then
                                nid = mesh%parts(iv_local)
                                call recv_cpus_list%insert(nid)
                                call recv_t(nid+1)%insert(it)
                                exit
                            end if
                        end do
                    end if
                end if
             end do
#else
             nid = mesh%parts(iv_global)
             call recv_cpus_list%insert(nid)
             call recv_v(nid+1)%insert(iv_global)
             do n = 1, mesh%vtx_nnb(iv_global)
                ie = mesh%vtx_ed%f(iv_global, n)
                it = mesh%vtx_tr%f(iv_global, n)
                call recv_e(nid+1)%insert(ie)
                call recv_t(nid+1)%insert(it)
             end do
#endif
          end do
#ifdef LAM_DOMAIN
          if (sw .eq. stencil_width) then
             do i = 1, halo%ne
                ie = halo%e_index(i)
                marke = mesh%edt_mark%f(ie, 1)
                if (marke .eq. 1771) then
                   nid = mesh%parts(int(mesh%edt_v%f(ie, 1)))
                   if (nid .eq. iblock) nid = mesh%parts(int(mesh%edt_v%f(ie, 2)))
                   if (nid .ne. iblock .and. .not. local_block%halo(stencil_width-1)%e_index_list%find(ie)) then
                      call recv_cpus_list%insert(nid)
                      call recv_e(nid+1)%insert(ie)
                   end if
                end if
             end do
          end if
#endif
       end do
       do sw = 1, stencil_width
          call local_block%bdry(sw)%t_index_list%final
          call local_block%bdry(sw)%e_index_list%final
          call local_block%halo(sw)%t_index_list%final
          call local_block%halo(sw)%e_index_list%final
       end do

       call recv_cpus_list%dump(recv_cpus)
       call recv_cpus_list%final

       local_block%nbcpuid = recv_cpus
       
       local_block%nnb = size(recv_cpus)
       
       allocate(local_block%recv_sites(local_block%nnb))

       do c = 1, local_block%nnb
          cpuid = recv_cpus(c)

          local_block%recv_sites(c)%cpuid = cpuid

          call recv_v(cpuid+1)%dump(local_block%recv_sites(c)%v)
          call recv_e(cpuid+1)%dump(local_block%recv_sites(c)%e)
          call recv_t(cpuid+1)%dump(local_block%recv_sites(c)%t)

          call recv_v(cpuid+1)%final
          call recv_e(cpuid+1)%final
          call recv_t(cpuid+1)%final
          ! if(iblock == 1) then
          !    open(10, file="recv_"//i2s(cpuid)//".txt")
          !    write(10, "(I0)", advance='yes') local_block%recv_sites(c)%v
          ! end if
       end do
       deallocate(recv_v)
       deallocate(recv_e)
       deallocate(recv_t)

       allocate(num_recv_sites(3, local_block%nnb))
       allocate(num_send_sites(3, local_block%nnb))

       allocate(reqs_send(local_block%nnb * 3))
       allocate(reqs_recv(local_block%nnb * 3))

       num_recv_sites = -1
       num_send_sites = -1
       !---------------------------------
       do c = 1, local_block%nnb
          num_recv_sites(1, c) = size(local_block%recv_sites(c)%v)
          num_recv_sites(2, c) = size(local_block%recv_sites(c)%e)
          num_recv_sites(3, c) = size(local_block%recv_sites(c)%t)
          
          dst = local_block%nbcpuid(c)
          
          call MPI_Isend(num_recv_sites(:,c), 3, MPI_INTEGER, dst, 101, &
               local_block%comm, reqs_send(c), ierr)

          call MPI_Irecv(num_send_sites(:,c), 3, MPI_INTEGER, dst, 101, &
               local_block%comm, reqs_recv(c), ierr)

       end do

       allocate(status(MPI_STATUS_SIZE, local_block%nnb*3))
       
       call MPI_WaitAll(local_block%nnb, reqs_send, status, ierr)
       call MPI_WaitAll(local_block%nnb, reqs_recv, status, ierr)
       
       allocate(local_block%send_sites(local_block%nnb))
       
       !---------------------------------
       do c = 1, local_block%nnb
          dst = local_block%nbcpuid(c)
          
          call MPI_Isend(local_block%recv_sites(c)%v, &
               size(local_block%recv_sites(c)%v), MPI_INTEGER, &
               dst, 102, local_block%comm, reqs_send(3*c-2), ierr)

          call MPI_Isend(local_block%recv_sites(c)%e, &
               size(local_block%recv_sites(c)%e), MPI_INTEGER, &
               dst, 102, local_block%comm, reqs_send(3*c-1), ierr)

          call MPI_Isend(local_block%recv_sites(c)%t, &
               size(local_block%recv_sites(c)%t), MPI_INTEGER, &
               dst, 102, local_block%comm, reqs_send(3*c),   ierr)
          
          local_block%send_sites(c)%cpuid = dst
          
          nv = num_send_sites(1, c)
          ne = num_send_sites(2, c)
          nt = num_send_sites(3, c)
          
          allocate(local_block%send_sites(c)%v(nv))
          allocate(local_block%send_sites(c)%e(ne))
          allocate(local_block%send_sites(c)%t(nt))
          
          call MPI_Irecv(local_block%send_sites(c)%v, nv, MPI_INTEGER, &
               dst, 102, local_block%comm, reqs_recv(3*c-2), ierr)

          call MPI_Irecv(local_block%send_sites(c)%e, ne, MPI_INTEGER, &
               dst, 102, local_block%comm, reqs_recv(3*c-1), ierr)

          call MPI_Irecv(local_block%send_sites(c)%t, nt, MPI_INTEGER, &
               dst, 102, local_block%comm, reqs_recv(3*c), ierr)
          
       end do

       call MPI_WaitAll(3*local_block%nnb, reqs_send, status, ierr)
       call MPI_WaitAll(3*local_block%nnb, reqs_recv, status, ierr)
       deallocate(num_recv_sites)
       deallocate(num_send_sites)
       deallocate(reqs_send)
       deallocate(reqs_recv)
       deallocate(status)
       !! test
       !do c = 1, local_block%nnb
       !   open(101, file="send-v-"//i2s(iblock+1)//"-"//i2s(local_block%nbcpuid(c))//".txt", form='formatted')
       !   do i = 1, size(local_block%send_sites(c)%v)
       !      write(101, *) local_block%send_sites(c)%v(i)
       !   end do
       !   close(101)
       !   open(101, file="recv-v-"//i2s(iblock+1)//"-"//i2s(local_block%nbcpuid(c))//".txt", form='formatted')
       !   do i = 1, size(local_block%recv_sites(c)%v)
       !      write(101, *) local_block%recv_sites(c)%v(i)
       !   end do
       !   close(101)
       !end do
    else
      ! read comm (read by PE 0, then scatter to all PEs)
      if(iblock .eq. 0) then
        open(1, file=trim(pardir)//"comm-"//i2s(nprocs), form='unformatted', access="stream")
        ! number of neighbor cpus
        allocate(counts(nprocs))
        read(1) counts
        sums = sum(counts)
        ! id of neighbor cpus
        allocate(nbcpus(sums))
        read(1) nbcpus
        ! 6 numbers for each neighbor cpu
        allocate(ele6_all(6*sums))
        read(1) ele6_all
        ! all indices for each neighbor cpu
        sums = sum(ele6_all)
        allocate(buffer_all(sums))
        read(1) buffer_all
        close(1)
        allocate(displs(nprocs))
        displs(1) = 0
        do i = 2, nprocs
          displs(i) = displs(i-1) + counts(i-1)
        end do
        allocate(counts_idx(nprocs))
        allocate(displs_idx(nprocs))
        pos = 1
        do i = 1, nprocs
          counts_idx(i) = 0
          do c = 1, counts(i)
            do n = 1, 6
              counts_idx(i) = counts_idx(i) + ele6_all(pos)
              pos = pos + 1
            end do
          end do
        end do
        displs_idx(1) = 0
        do i = 2, nprocs
          displs_idx(i) = displs_idx(i-1) + counts_idx(i-1)
        end do
      else
        allocate(counts(1))
        allocate(displs(1))
        allocate(nbcpus(1))
        allocate(ele6_all(1))
        allocate(buffer_all(1))
        allocate(counts_idx(1))
        allocate(displs_idx(1))
      end if
      call MPI_Scatter(counts, 1, MPI_INTEGER, local_block%nnb, &
                       1, MPI_INTEGER, 0, local_block%comm, ierr)
      allocate(local_block%nbcpuid(local_block%nnb))
      call MPI_Scatterv(nbcpus, counts, displs, MPI_INTEGER, local_block%nbcpuid, &
                        local_block%nnb, MPI_INTEGER, 0, local_block%comm, ierr)
      allocate(ele6(6*local_block%nnb))
      call MPI_Scatterv(ele6_all, 6*counts, 6*displs, MPI_INTEGER, ele6, &
                        6*local_block%nnb, MPI_INTEGER, 0, local_block%comm, ierr)
      sums = sum(ele6)
      allocate(buffer(sums))
      call MPI_Scatterv(buffer_all, counts_idx, displs_idx, MPI_INTEGER, buffer, &
                        sums, MPI_INTEGER, 0, local_block%comm, ierr)

      allocate(local_block%recv_sites(local_block%nnb))
      allocate(local_block%send_sites(local_block%nnb))
      do c = 1, local_block%nnb
        local_block%send_sites(c)%cpuid = local_block%nbcpuid(c)
        local_block%recv_sites(c)%cpuid = local_block%nbcpuid(c)
        allocate(local_block%send_sites(c)%v(ele6((c-1)*6+1)))
        allocate(local_block%recv_sites(c)%v(ele6((c-1)*6+2)))
        allocate(local_block%send_sites(c)%e(ele6((c-1)*6+3)))
        allocate(local_block%recv_sites(c)%e(ele6((c-1)*6+4)))
        allocate(local_block%send_sites(c)%t(ele6((c-1)*6+5)))
        allocate(local_block%recv_sites(c)%t(ele6((c-1)*6+6)))
      end do
      poe = 0
      do c = 1, local_block%nnb
        pos = poe + 1
        poe = pos + ele6((c-1)*6+1)-1
        local_block%send_sites(c)%v = buffer(pos:poe)
        pos = poe + 1
        poe = pos + ele6((c-1)*6+2)-1
        local_block%recv_sites(c)%v = buffer(pos:poe)
        pos = poe + 1
        poe = pos + ele6((c-1)*6+3)-1
        local_block%send_sites(c)%e = buffer(pos:poe)
        pos = poe + 1
        poe = pos + ele6((c-1)*6+4)-1
        local_block%recv_sites(c)%e = buffer(pos:poe)
        pos = poe + 1
        poe = pos + ele6((c-1)*6+5)-1
        local_block%send_sites(c)%t = buffer(pos:poe)
        pos = poe + 1
        poe = pos + ele6((c-1)*6+6)-1
        local_block%recv_sites(c)%t = buffer(pos:poe)
      end do
      deallocate(counts)
      deallocate(displs)
      deallocate(nbcpus)
      deallocate(ele6_all)
      deallocate(ele6)
      deallocate(buffer_all)
      deallocate(counts_idx)
      deallocate(displs_idx)
      deallocate(buffer)
    end if
    return
  end subroutine


  subroutine global2local(local_block)
    implicit none
    type(block_structure), target, intent(inout) :: local_block
    integer(i4)                                  :: iv_domain
    type(global_domain),   pointer               :: full_domain, iner_domain, halo, bdry    
    integer(i4)                                  :: sw, iv_global
    type(map),             pointer               :: g2l_v, g2l_t, g2l_e

    full_domain => local_block%full_domain

    g2l_v => full_domain%map_g2l_v
    g2l_e => full_domain%map_g2l_e
    g2l_t => full_domain%map_g2l_t

    do iv_domain = 1, full_domain%nv
       iv_global = full_domain%v_index(iv_domain)
       call g2l_v%insert(iv_global, iv_domain)
    end do

    do iv_domain = 1, full_domain%nt
       iv_global = full_domain%t_index(iv_domain)
       call g2l_t%insert(iv_global, iv_domain)
    end do

    do iv_domain = 1, full_domain%ne
       iv_global = full_domain%e_index(iv_domain)
       call g2l_e%insert(iv_global, iv_domain)
    end do

    call global2local_domain(g2l_v, g2l_e, g2l_t, full_domain)
    !iner_domain => local_block%iner_domain
    !call global2local_domain(g2l_v, g2l_e, g2l_t, iner_domain)

    !do sw = 1, stencil_width
    !   halo => local_block%halo(sw)
    !   bdry => local_block%bdry(sw)
    !   
    !   call global2local_domain(g2l_v, g2l_e, g2l_t, halo)
    !   call global2local_domain(g2l_v, g2l_e, g2l_t, bdry)
    !end do
    call global2local_comm(local_block)    

  end subroutine


  subroutine global2local_domain(g2l_v, g2l_e, g2l_t, dm)
    implicit none
    type(global_domain),  target, intent(inout) :: dm
    type(map),                    intent(in)    :: g2l_v, g2l_e, g2l_t
    integer(i4)                                 :: iv_domain, n
    type(node_structure), pointer               :: vtx
    type(dual_structure), pointer               :: tri
    type(edge_structure), pointer               :: edt
    type(edge_structure), pointer               :: edp
    type(prim_structure), pointer               :: plg

    do iv_domain = 1, dm%nv
       vtx => dm%vtx(iv_domain)
       plg => dm%plg(iv_domain)
       do n = 1, vtx%nnb
          vtx%nb(n) = g2l_v%find(vtx%nb(n))
          vtx%ed(n) = g2l_e%find(vtx%ed(n))
          vtx%tr(n) = g2l_t%find(vtx%tr(n))
       end do
    end do

    do iv_domain = 1, dm%nt
       tri => dm%tri(iv_domain)
       do n = 1, tri%nnb
          tri%v(n)  = g2l_v%find(tri%v(n))
          tri%ed(n) = g2l_e%find(tri%ed(n))
          tri%nb(n) = g2l_t%find(tri%nb(n))
       end do
    end do

    do iv_domain = 1, dm%ne
       edt => dm%edt(iv_domain)
       edp => dm%edp(iv_domain)
       do n = 1, 2
          edt%v(n) = g2l_v%find(edt%v(n))
          edp%v(n) = g2l_t%find(edp%v(n))
          !edp%sh(n) = g2l_v%find(edp%sh(n))
       end do
    end do

  end subroutine


  subroutine global2local_comm(local_block)
    implicit none
    type(block_structure), target, intent(inout) :: local_block
    integer                                      :: c, n_v, n_t, n_e, iv, it, ie
    integer(i4)                                  :: max_vr, max_vs, max_tr, max_ts, max_er, max_es
    type(map),             pointer               :: g2l_v, g2l_t, g2l_e

    g2l_v => local_block%full_domain%map_g2l_v
    g2l_e => local_block%full_domain%map_g2l_e
    g2l_t => local_block%full_domain%map_g2l_t

    max_vr=0; max_tr=0; max_er=0
    max_vs=0; max_ts=0; max_es=0

    do c = 1, local_block%nnb
      n_v = size(local_block%recv_sites(c)%v)
      n_t = size(local_block%recv_sites(c)%t)
      n_e = size(local_block%recv_sites(c)%e)
      local_block%recv_sites(c)%nv = n_v
      local_block%recv_sites(c)%nt = n_t
      local_block%recv_sites(c)%ne = n_e
      if(n_v .gt. max_vr) max_vr = n_v
      if(n_t .gt. max_tr) max_tr = n_t
      if(n_e .gt. max_er) max_er = n_e
      allocate(local_block%recv_sites(c)%local_v(n_v))
      allocate(local_block%recv_sites(c)%local_t(n_t))
      allocate(local_block%recv_sites(c)%local_e(n_e))
      do iv=1, n_v
        local_block%recv_sites(c)%local_v(iv) = g2l_v%find(local_block%recv_sites(c)%v(iv))
      end do
      do it=1, n_t
        local_block%recv_sites(c)%local_t(it) = g2l_t%find(local_block%recv_sites(c)%t(it))
      end do
      do ie=1, n_e
        local_block%recv_sites(c)%local_e(ie) = g2l_e%find(local_block%recv_sites(c)%e(ie))
      end do
    enddo

    do c = 1, local_block%nnb
      n_v = size(local_block%send_sites(c)%v)
      n_t = size(local_block%send_sites(c)%t)
      n_e = size(local_block%send_sites(c)%e)
      local_block%send_sites(c)%nv = n_v
      local_block%send_sites(c)%nt = n_t
      local_block%send_sites(c)%ne = n_e
      if(n_v .gt. max_vs) max_vs = n_v
      if(n_t .gt. max_ts) max_ts = n_t
      if(n_e .gt. max_es) max_es = n_e
      allocate(local_block%send_sites(c)%local_v(n_v))
      allocate(local_block%send_sites(c)%local_t(n_t))
      allocate(local_block%send_sites(c)%local_e(n_e))
      do iv=1, n_v
        local_block%send_sites(c)%local_v(iv) = g2l_v%find(local_block%send_sites(c)%v(iv))
      end do
      do it=1, n_t
        local_block%send_sites(c)%local_t(it) = g2l_t%find(local_block%send_sites(c)%t(it))
      end do
      do ie=1, n_e
        local_block%send_sites(c)%local_e(ie) = g2l_e%find(local_block%send_sites(c)%e(ie))
      end do
    enddo

    local_block%max_nv_send = max_vs
    local_block%max_nv_recv = max_vr
    local_block%max_nt_send = max_ts
    local_block%max_nt_recv = max_tr
    local_block%max_ne_send = max_es
    local_block%max_ne_recv = max_er

  end subroutine
  
  
  subroutine exchange_data_1d_add(mesh,field_head,scalar_field_data,mix_flag)
    !IO
    type(global_domain)                                  :: mesh
    type(exchange_field_list_1d), pointer, intent(inout) :: field_head
    type(scalar_1d_field)       , target,  intent(inout) :: scalar_field_data
    character(len=*), optional, intent(in)               :: mix_flag

    !local
    type(exchange_field_list_1d), pointer, save          :: field_list

    if(.not.associated(field_head))then
      allocate(field_head)
      field_list => field_head
    else if(.not.associated(field_list%next))then
      allocate(field_list%next)
      field_list => field_list%next
    else
      print*,"field_list error"
    end if
    
    field_list%field_data => scalar_field_data
    if(present(mix_flag).and.trim(mix_flag).eq."r4")then
       field_list%dim1 = size(scalar_field_data%f_r4)
    else
       field_list%dim1 = size(scalar_field_data%f)
    end if
    !IF(scalar_field_data%pos .eq. -1)then
      if((field_list%dim1 .eq. mesh%nv) .or. (field_list%dim1 .eq. mesh%nv_full)) then ! v location
        scalar_field_data%pos = 0
      else if((field_list%dim1 .eq. mesh%nt) .or. (field_list%dim1 .eq. mesh%nt_full)) then ! t location
        scalar_field_data%pos = 1
      else if((field_list%dim1 .eq. mesh%ne) .or. (field_list%dim1 .eq. mesh%ne_full)) then ! e location
        scalar_field_data%pos = 6
      else
        print*,"Invalid scalar field data type"
      end if
    !END IF
  end subroutine exchange_data_1d_add

  subroutine exchange_data_1d_clean(field_head)
    !IO
    type(exchange_field_list_1d), pointer, intent(inout) :: field_head

    !local                              
    type(exchange_field_list_1d), pointer                :: field_list, tmpf
    
    field_list => field_head
    do while(associated(field_list))
      tmpf => field_list%next
      deallocate(field_list)
      field_list => tmpf
    end do
    field_head => null()

  end subroutine exchange_data_1d_clean

  subroutine exchange_data_1d(local_block,field_head)
    !io
    use omp_lib
    type(block_structure),        target,  intent(in)    :: local_block
    type(exchange_field_list_1d), pointer, intent(inout) :: field_head

    !local
    integer,  allocatable                                :: reqs_send(:), reqs_recv(:), reqs(:)
    integer                                              :: ierr, c, iv, ie, it, dst_s, dst_r
    integer                                              :: i_list, iblock, iv_list, it_list, ie_list, n, nr, ns
    integer                                              :: list_n, tmp_n, n_v, n_t, n_e, list_nv, list_nt, list_ne
    INTEGER,  allocatable                                :: status(:,:) !STATUS(MPI_STATUS_SIZE)

    real(r8), allocatable                                :: field_data_nv_send(:,:)!bdry
    real(r8), allocatable                                :: field_data_ne_send(:,:)
    real(r8), allocatable                                :: field_data_nt_send(:,:)
    
    real(r8), allocatable                                :: field_data_nv_recv(:,:)!halo
    real(r8), allocatable                                :: field_data_ne_recv(:,:)
    real(r8), allocatable                                :: field_data_nt_recv(:,:)

    type(exchange_field_list_1d), pointer                :: tmpf

    tmpf => field_head
    list_n = 0
    iblock = mpi_rank()
     
    list_nv=0; list_nt=0; list_ne=0
    do while(associated(tmpf))
      select case(tmpf%field_data%pos)
      case(0)
        list_nv = list_nv + 1
      case(1)
        list_nt = list_nt + 1
      case(6)
        list_ne = list_ne + 1
      case default
        stop
      end select
      list_n = list_n + 1   
      tmpf =>  tmpf%next
    end do
    n = 0
    if(list_nv .gt. 0) n = n + 1
    if(list_nt .gt. 0) n = n + 1
    if(list_ne .gt. 0) n = n + 1
    
    if(list_nv .gt. 0) allocate(field_data_nv_send(list_nv*local_block%max_nv_send,local_block%nnb), source=0._r8)
    if(list_nv .gt. 0) allocate(field_data_nv_recv(list_nv*local_block%max_nv_recv,local_block%nnb))
    if(list_nt .gt. 0) allocate(field_data_nt_send(list_nt*local_block%max_nt_send,local_block%nnb), source=0._r8)
    if(list_nt .gt. 0) allocate(field_data_nt_recv(list_nt*local_block%max_nt_recv,local_block%nnb))
    if(list_ne .gt. 0) allocate(field_data_ne_send(list_ne*local_block%max_ne_send,local_block%nnb), source=0._r8)
    if(list_ne .gt. 0) allocate(field_data_ne_recv(list_ne*local_block%max_ne_recv,local_block%nnb))

    !assign send_data value
    allocate(reqs_send(local_block%nnb * n))
    allocate(reqs_recv(local_block%nnb * n))
    allocate(status(MPI_STATUS_SIZE, local_block%nnb*n))
    ns=1; nr=1
    do c = 1, local_block%nnb
      dst_s = local_block%send_sites(c)%cpuid
      dst_r = local_block%recv_sites(c)%cpuid
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%send_sites(c)%nv
!$omp parallel  private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%send_sites(c)%local_v(iv)
            field_data_nv_send(iv+n_v*iv_list,c) = tmpf%field_data%f(tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%send_sites(c)%nt
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%send_sites(c)%local_t(it)
            field_data_nt_send(it+n_t*it_list,c) = tmpf%field_data%f(tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%send_sites(c)%ne
!$omp parallel  private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%send_sites(c)%local_e(ie)
            field_data_ne_send(ie+n_e*ie_list,c) = tmpf%field_data%f(tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select  
        tmpf => tmpf%next
      end do

      if (list_nv .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_nv_send(:,c), &
                         local_block%send_sites(c)%nv*list_nv, MPI_REAL8, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,c), &
                         local_block%recv_sites(c)%nv*list_nv, MPI_REAL8, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_nv_send(:,c), &
                         local_block%send_sites(c)%nv*list_nv, MPI_REAL, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,c), &
                         local_block%recv_sites(c)%nv*list_nv, MPI_REAL, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif

      if (list_nt .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_nt_send(:,c), &
                         local_block%send_sites(c)%nt*list_nt, MPI_REAL8, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,c), &
                         local_block%recv_sites(c)%nt*list_nt, MPI_REAL8, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_nt_send(:,c), &
                         local_block%send_sites(c)%nt*list_nt, MPI_REAL, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,c), &
                         local_block%recv_sites(c)%nt*list_nt, MPI_REAL, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif

      if (list_ne .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_ne_send(:,c), &
                         local_block%send_sites(c)%ne*list_ne, MPI_REAL8, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,c), &
                         local_block%recv_sites(c)%ne*list_ne, MPI_REAL8, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_ne_send(:,c), &
                         local_block%send_sites(c)%ne*list_ne, MPI_REAL, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,c), &
                         local_block%recv_sites(c)%ne*list_ne, MPI_REAL, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif

    end do

    call MPI_WaitAll(n*local_block%nnb, reqs_send, status, ierr)
    call MPI_WaitAll(n*local_block%nnb, reqs_recv, status, ierr)

    !get recv_data value
    do c = 1, local_block%nnb
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list=0,list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%recv_sites(c)%nv
!$omp parallel  private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%recv_sites(c)%local_v(iv)
            tmpf%field_data%f(tmp_n) = field_data_nv_recv(iv+n_v*iv_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%recv_sites(c)%nt
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%recv_sites(c)%local_t(it)
            tmpf%field_data%f(tmp_n) = field_data_nt_recv(it+n_t*it_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%recv_sites(c)%ne
!$omp parallel  private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%recv_sites(c)%local_e(ie)
            tmpf%field_data%f(tmp_n) = field_data_ne_recv(ie+n_e*ie_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select
        tmpf => tmpf%next
      end do
    end do

    call exchange_data_1d_clean(field_head)

    if(.false.)  then
        
      do c = 1, local_block%nnb
       open(101, file=i2s(iblock)//".send.cell."//i2s(local_block%send_sites(c)%cpuid)//".txt")
           n_v = size(local_block%send_sites(c)%v)
       write(101, "(G20.10)", advance='yes') field_data_nv_send(1:n_v,c)
       close(101)

       open(2, file=i2s(iblock)//".recv.cell."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
           n_v = size(local_block%recv_sites(c)%v)
       write(2, "(G20.10)", advance='yes') field_data_nv_recv(1:n_v,c)
       close(2)

       open(101, file=i2s(iblock)//".send.edge."//i2s(local_block%send_sites(c)%cpuid)//".txt")
           n_e = size(local_block%send_sites(c)%e)
       write(101, "(G20.10)", advance='yes') field_data_ne_send(1:n_e,c)
       close(101)

       open(2, file=i2s(iblock)//".recv.edge."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
           n_e = size(local_block%recv_sites(c)%e)
       write(2, "(G20.10)", advance='yes') field_data_ne_recv(1:n_e,c)
       close(2)

      end do
    end if

!clean  
    if(list_nv .gt. 0) deallocate(field_data_nv_send)
    if(list_nv .gt. 0) deallocate(field_data_nv_recv)    
    if(list_nt .gt. 0) deallocate(field_data_nt_send)
    if(list_nt .gt. 0) deallocate(field_data_nt_recv)    
    if(list_ne .gt. 0) deallocate(field_data_ne_send)
    if(list_ne .gt. 0) deallocate(field_data_ne_recv)
    deallocate(reqs_send)
    deallocate(reqs_recv)
    deallocate(status)

  end subroutine exchange_data_1d

  subroutine exchange_data_1d_r4(local_block,field_head)
    !io
    use omp_lib
    type(block_structure),        target,  intent(in)    :: local_block
    type(exchange_field_list_1d), pointer, intent(inout) :: field_head

    !local
    integer,  allocatable                                :: reqs_send(:), reqs_recv(:), reqs(:)
    integer                                              :: ierr, c, iv, ie, it, dst_s, dst_r
    integer                                              :: i_list, iblock, iv_list, it_list, ie_list, n, nr, ns
    integer                                              :: list_n, tmp_n, n_v, n_t, n_e, list_nv, list_nt, list_ne
    INTEGER,  allocatable                                :: status(:,:) !STATUS(MPI_STATUS_SIZE)

    real(r4), allocatable                                :: field_data_nv_send(:,:)!bdry
    real(r4), allocatable                                :: field_data_ne_send(:,:)
    real(r4), allocatable                                :: field_data_nt_send(:,:)
    
    real(r4), allocatable                                :: field_data_nv_recv(:,:)!halo
    real(r4), allocatable                                :: field_data_ne_recv(:,:)
    real(r4), allocatable                                :: field_data_nt_recv(:,:)

    type(exchange_field_list_1d), pointer                :: tmpf

    tmpf => field_head
    list_n = 0
    iblock = mpi_rank()
     
    list_nv=0; list_nt=0; list_ne=0
    do while(associated(tmpf))
      select case(tmpf%field_data%pos)
      case(0)
        list_nv = list_nv + 1
      case(1)
        list_nt = list_nt + 1
      case(6)
        list_ne = list_ne + 1
      case default
        stop
      end select
      list_n = list_n + 1   
      tmpf =>  tmpf%next
    end do
    n = 0
    if(list_nv .gt. 0) n = n + 1
    if(list_nt .gt. 0) n = n + 1
    if(list_ne .gt. 0) n = n + 1
    
    if(list_nv .gt. 0) allocate(field_data_nv_send(list_nv*local_block%max_nv_send,local_block%nnb), source=0._r4)
    if(list_nv .gt. 0) allocate(field_data_nv_recv(list_nv*local_block%max_nv_recv,local_block%nnb))
    if(list_nt .gt. 0) allocate(field_data_nt_send(list_nt*local_block%max_nt_send,local_block%nnb), source=0._r4)
    if(list_nt .gt. 0) allocate(field_data_nt_recv(list_nt*local_block%max_nt_recv,local_block%nnb))
    if(list_ne .gt. 0) allocate(field_data_ne_send(list_ne*local_block%max_ne_send,local_block%nnb), source=0._r4)
    if(list_ne .gt. 0) allocate(field_data_ne_recv(list_ne*local_block%max_ne_recv,local_block%nnb))

    !assign send_data value
    allocate(reqs_send(local_block%nnb * n))
    allocate(reqs_recv(local_block%nnb * n))
    allocate(status(MPI_STATUS_SIZE, local_block%nnb*n))
    ns=1; nr=1
    do c = 1, local_block%nnb
      dst_s = local_block%send_sites(c)%cpuid
      dst_r = local_block%recv_sites(c)%cpuid
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%send_sites(c)%nv
!$omp parallel  private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%send_sites(c)%local_v(iv)
            field_data_nv_send(iv+n_v*iv_list,c) = tmpf%field_data%f_r4(tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%send_sites(c)%nt
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%send_sites(c)%local_t(it)
            field_data_nt_send(it+n_t*it_list,c) = tmpf%field_data%f_r4(tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%send_sites(c)%ne
!$omp parallel  private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%send_sites(c)%local_e(ie)
            field_data_ne_send(ie+n_e*ie_list,c) = tmpf%field_data%f_r4(tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select  
        tmpf => tmpf%next
      end do

      if (list_nv .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_nv_send(:,c), &
                         local_block%send_sites(c)%nv*list_nv, MPI_REAL8, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,c), &
                         local_block%recv_sites(c)%nv*list_nv, MPI_REAL8, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_nv_send(:,c), &
                         local_block%send_sites(c)%nv*list_nv, MPI_REAL, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,c), &
                         local_block%recv_sites(c)%nv*list_nv, MPI_REAL, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif

      if (list_nt .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_nt_send(:,c), &
                         local_block%send_sites(c)%nt*list_nt, MPI_REAL8, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,c), &
                         local_block%recv_sites(c)%nt*list_nt, MPI_REAL8, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_nt_send(:,c), &
                         local_block%send_sites(c)%nt*list_nt, MPI_REAL, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,c), &
                         local_block%recv_sites(c)%nt*list_nt, MPI_REAL, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif

      if (list_ne .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_ne_send(:,c), &
                         local_block%send_sites(c)%ne*list_ne, MPI_REAL8, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,c), &
                         local_block%recv_sites(c)%ne*list_ne, MPI_REAL8, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_ne_send(:,c), &
                         local_block%send_sites(c)%ne*list_ne, MPI_REAL, &
                         dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,c), &
                         local_block%recv_sites(c)%ne*list_ne, MPI_REAL, &
                         dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif

    end do

    call MPI_WaitAll(n*local_block%nnb, reqs_send, status, ierr)
    call MPI_WaitAll(n*local_block%nnb, reqs_recv, status, ierr)

    !get recv_data value
    do c = 1, local_block%nnb
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list=0,list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%recv_sites(c)%nv
!$omp parallel  private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%recv_sites(c)%local_v(iv)
            tmpf%field_data%f_r4(tmp_n) = field_data_nv_recv(iv+n_v*iv_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%recv_sites(c)%nt
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%recv_sites(c)%local_t(it)
            tmpf%field_data%f_r4(tmp_n) = field_data_nt_recv(it+n_t*it_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%recv_sites(c)%ne
!$omp parallel  private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%recv_sites(c)%local_e(ie)
            tmpf%field_data%f_r4(tmp_n) = field_data_ne_recv(ie+n_e*ie_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select
        tmpf => tmpf%next
      end do
    end do

    call exchange_data_1d_clean(field_head)

    if(.false.)  then
        
      do c = 1, local_block%nnb
       open(101, file=i2s(iblock)//".send.cell."//i2s(local_block%send_sites(c)%cpuid)//".txt")
           n_v = size(local_block%send_sites(c)%v)
       write(101, "(G20.10)", advance='yes') field_data_nv_send(1:n_v,c)
       close(101)

       open(2, file=i2s(iblock)//".recv.cell."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
           n_v = size(local_block%recv_sites(c)%v)
       write(2, "(G20.10)", advance='yes') field_data_nv_recv(1:n_v,c)
       close(2)

       open(101, file=i2s(iblock)//".send.edge."//i2s(local_block%send_sites(c)%cpuid)//".txt")
           n_e = size(local_block%send_sites(c)%e)
       write(101, "(G20.10)", advance='yes') field_data_ne_send(1:n_e,c)
       close(101)

       open(2, file=i2s(iblock)//".recv.edge."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
           n_e = size(local_block%recv_sites(c)%e)
       write(2, "(G20.10)", advance='yes') field_data_ne_recv(1:n_e,c)
       close(2)

      end do
    end if

!clean  
    if(list_nv .gt. 0) deallocate(field_data_nv_send)
    if(list_nv .gt. 0) deallocate(field_data_nv_recv)    
    if(list_nt .gt. 0) deallocate(field_data_nt_send)
    if(list_nt .gt. 0) deallocate(field_data_nt_recv)    
    if(list_ne .gt. 0) deallocate(field_data_ne_send)
    if(list_ne .gt. 0) deallocate(field_data_ne_recv)
    deallocate(reqs_send)
    deallocate(reqs_recv)
    deallocate(status)

  end subroutine exchange_data_1d_r4
  
  subroutine exchange_data_2d_add(mesh,field_head,scalar_field_data,mix_flag)
    !IO
    type(global_domain)                                  :: mesh
    type(exchange_field_list_2d), intent(inout), pointer :: field_head
    type(scalar_2d_field),        intent(inout), target  :: scalar_field_data
    character(len=*), optional, intent(in)               :: mix_flag

    !local
    type(exchange_field_list_2d), save,          pointer :: field_list
    
    if(.not.associated(field_head))then
      allocate(field_head)
      field_list => field_head
    else if(.not.associated(field_list%next))then
      allocate(field_list%next)
      field_list => field_list%next
    else
      print*,"field_list error"
    end if

    field_list%field_data => scalar_field_data
    if(present(mix_flag).and.trim(mix_flag).eq."r4")then
       field_list%dim1 = size(scalar_field_data%f_r4, 1)
       field_list%dim2 = size(scalar_field_data%f_r4, 2)
    else
       field_list%dim1 = size(scalar_field_data%f, 1)
       field_list%dim2 = size(scalar_field_data%f, 2)
    end if
    
    !IF(scalar_field_data%pos .eq. -1)then
      if((field_list%dim2 .eq. mesh%nv) .or. (field_list%dim2 .eq. mesh%nv_full)) then ! v location
        scalar_field_data%pos = 0
      else if((field_list%dim2 .eq. mesh%nt) .or. (field_list%dim2 .eq. mesh%nt_full)) then ! t location
        scalar_field_data%pos = 1
      else if((field_list%dim2 .eq. mesh%ne) .or. (field_list%dim2 .eq. mesh%ne_full)) then ! e location
        scalar_field_data%pos = 6
      else
        print*,"Invalid scalar field data type"
      end if
    !END IF

  end subroutine exchange_data_2d_add

  subroutine exchange_data_2d_clean(field_head)
    !IO
    type(exchange_field_list_2d), pointer, intent(inout) :: field_head

    !local
    type(exchange_field_list_2d), pointer                :: field_list, tmpf

    field_list => field_head
    do while(associated(field_list))
      tmpf => field_list%next
      deallocate(field_list)
      field_list => tmpf
    end do
    field_head => null()

  end subroutine exchange_data_2d_clean

  subroutine exchange_data_2d(local_block,field_head)
    !io
    use omp_lib
    type(block_structure),        target,  intent(in)    :: local_block
    type(exchange_field_list_2d), pointer, intent(inout) :: field_head

    !local
    integer,  allocatable                                :: reqs_send(:), reqs_recv(:), reqs(:)
    integer                                              :: ierr, c, iv, ie, it, dst_s, dst_r
    integer                                              :: i_list, iblock, iv_list, it_list, ie_list, n, nr, ns
    integer                                              :: list_n, tmp_n, n_v, n_t, n_e, list_nv, list_nt, list_ne
    integer                                              :: maxdim1, sizedim1

    INTEGER,  allocatable                                :: status(:,:) !STATUS(MPI_STATUS_SIZE)

    real(r8), allocatable                                :: field_data_nv_send(:,:,:)!bdry
    real(r8), allocatable                                :: field_data_ne_send(:,:,:)
    real(r8), allocatable                                :: field_data_nt_send(:,:,:)

    real(r8), allocatable                                :: field_data_nv_recv(:,:,:)!halo
    real(r8), allocatable                                :: field_data_ne_recv(:,:,:)
    real(r8), allocatable                                :: field_data_nt_recv(:,:,:)
    type(exchange_field_list_2d), pointer                :: tmpf

    tmpf => field_head
    maxdim1 = 0
    do while(associated(tmpf))
      maxdim1 = max(maxdim1, tmpf%dim1)
      tmpf =>  tmpf%next
    end do
    list_n = 0
    iblock = mpi_rank()

    list_nv=0; list_nt=0; list_ne=0
    tmpf => field_head
    do while(associated(tmpf))
      select case(tmpf%field_data%pos)
      case(0)
        list_nv = list_nv + 1
      case(1)
        list_nt = list_nt + 1
      case(6)
        list_ne = list_ne + 1
      case default
        stop
      end select
      list_n = list_n + 1   
      tmpf =>  tmpf%next
    end do
    n=0
    if(list_nv.gt.0) n=n+1
    if(list_nt.gt.0) n=n+1
    if(list_ne.gt.0) n=n+1

    if(list_nv.gt.0) allocate(field_data_nv_send(maxdim1,list_nv*local_block%max_nv_send,local_block%nnb), source=0._r8)
    if(list_nv.gt.0) allocate(field_data_nv_recv(maxdim1,list_nv*local_block%max_nv_recv,local_block%nnb))
    if(list_nt.gt.0) allocate(field_data_nt_send(maxdim1,list_nt*local_block%max_nt_send,local_block%nnb), source=0._r8)
    if(list_nt.gt.0) allocate(field_data_nt_recv(maxdim1,list_nt*local_block%max_nt_recv,local_block%nnb))
    if(list_ne.gt.0) allocate(field_data_ne_send(maxdim1,list_ne*local_block%max_ne_send,local_block%nnb), source=0._r8)
    if(list_ne.gt.0) allocate(field_data_ne_recv(maxdim1,list_ne*local_block%max_ne_recv,local_block%nnb))

    !assign send_data value
    allocate(reqs_send(local_block%nnb * n))
    allocate(reqs_recv(local_block%nnb * n))
    allocate(status(MPI_STATUS_SIZE, local_block%nnb*n))
    ns=1; nr=1

    do c = 1, local_block%nnb
      dst_s = local_block%send_sites(c)%cpuid
      dst_r = local_block%recv_sites(c)%cpuid
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%send_sites(c)%nv
          sizedim1 = tmpf%dim1
!$omp parallel private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%send_sites(c)%local_v(iv)
            field_data_nv_send(1:sizedim1,iv+n_v*iv_list,c) = tmpf%field_data%f(1:sizedim1,tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%send_sites(c)%nt
          sizedim1 = tmpf%dim1
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%send_sites(c)%local_t(it)
            field_data_nt_send(1:sizedim1,it+n_t*it_list,c) = tmpf%field_data%f(1:sizedim1,tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%send_sites(c)%ne
          sizedim1 = tmpf%dim1
!$omp parallel private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%send_sites(c)%local_e(ie)
            field_data_ne_send(1:sizedim1,ie+n_e*ie_list,c) = tmpf%field_data%f(1:sizedim1,tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select
        tmpf => tmpf%next
      end do
      if (list_nv .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_nv_send(:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_nv_send(:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_nt .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_nt_send(:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_nt_send(:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_ne .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_ne_send(:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_ne_send(:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
    end do

    call MPI_WaitAll(n*local_block%nnb, reqs_send, status, ierr)
    call MPI_WaitAll(n*local_block%nnb, reqs_recv, status, ierr)

    !get recv_data value
    do c = 1, local_block%nnb
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0;
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%recv_sites(c)%nv
          sizedim1 = tmpf%dim1
!$omp parallel  private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%recv_sites(c)%local_v(iv)
            tmpf%field_data%f(1:sizedim1,tmp_n) = field_data_nv_recv(1:sizedim1,iv+n_v*iv_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%recv_sites(c)%nt
          sizedim1 = tmpf%dim1
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%recv_sites(c)%local_t(it)
            tmpf%field_data%f(1:sizedim1,tmp_n) = field_data_nt_recv(1:sizedim1,it+n_t*it_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%recv_sites(c)%ne
          sizedim1 = tmpf%dim1
!$omp parallel  private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%recv_sites(c)%local_e(ie)
            tmpf%field_data%f(1:sizedim1,tmp_n) = field_data_ne_recv(1:sizedim1,ie+n_e*ie_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select  
          tmpf => tmpf%next
      end do
    end do

    call exchange_data_2d_clean(field_head)

    if(.false.)  then
      do c = 1, local_block%nnb
         open(101, file=i2s(iblock)//".send.cell."//i2s(local_block%send_sites(c)%cpuid)//".txt")
         n_v = size(local_block%send_sites(c)%v)
         write(101, "(G20.10)", advance='yes') field_data_nv_send(:,1:n_v,c)
         close(101) 
         open(2, file=i2s(iblock)//".recv.cell."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
         n_v = size(local_block%recv_sites(c)%v)
         write(2, "(G20.10)", advance='yes') field_data_nv_recv(:,1:n_v,c)
         close(2)
         open(101, file=i2s(iblock)//".send.edge."//i2s(local_block%send_sites(c)%cpuid)//".txt")
         n_e = size(local_block%send_sites(c)%e)
         write(101, "(G20.10)", advance='yes') field_data_ne_send(:,1:n_e,c)
         close(101) 
         open(2, file=i2s(iblock)//".recv.edge."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
         n_e = size(local_block%recv_sites(c)%e)
         write(2, "(G20.10)", advance='yes') field_data_ne_recv(:,1:n_e,c)
         close(2)
      end do
    end if

    !clean  
    if(list_nv .gt. 0) deallocate(field_data_nv_send)
    if(list_nv .gt. 0) deallocate(field_data_nv_recv)    
    if(list_nt .gt. 0) deallocate(field_data_nt_send)
    if(list_nt .gt. 0) deallocate(field_data_nt_recv)    
    if(list_ne .gt. 0) deallocate(field_data_ne_send)
    if(list_ne .gt. 0) deallocate(field_data_ne_recv)
    deallocate(reqs_send)
    deallocate(reqs_recv)
    deallocate(status)

  end subroutine exchange_data_2d

 subroutine exchange_data_2d_r4(local_block,field_head)
    !io
    use omp_lib
    type(block_structure),        target,  intent(in)    :: local_block
    type(exchange_field_list_2d), pointer, intent(inout) :: field_head

    !local
    integer,  allocatable                                :: reqs_send(:), reqs_recv(:), reqs(:)
    integer                                              :: ierr, c, iv, ie, it, dst_s, dst_r
    integer                                              :: i_list, iblock, iv_list, it_list, ie_list, n, nr, ns
    integer                                              :: list_n, tmp_n, n_v, n_t, n_e, list_nv, list_nt, list_ne
    integer                                              :: maxdim1, sizedim1

    INTEGER,  allocatable                                :: status(:,:) !STATUS(MPI_STATUS_SIZE)

    real(r4), allocatable                                :: field_data_nv_send(:,:,:)!bdry
    real(r4), allocatable                                :: field_data_ne_send(:,:,:)
    real(r4), allocatable                                :: field_data_nt_send(:,:,:)

    real(r4), allocatable                                :: field_data_nv_recv(:,:,:)!halo
    real(r4), allocatable                                :: field_data_ne_recv(:,:,:)
    real(r4), allocatable                                :: field_data_nt_recv(:,:,:)
    type(exchange_field_list_2d), pointer                :: tmpf

    tmpf => field_head
    maxdim1 = 0
    do while(associated(tmpf))
      maxdim1 = max(maxdim1, tmpf%dim1)
      tmpf =>  tmpf%next
    end do
    list_n = 0
    iblock = mpi_rank()

    list_nv=0; list_nt=0; list_ne=0
    tmpf => field_head
    do while(associated(tmpf))
      select case(tmpf%field_data%pos)
      case(0)
        list_nv = list_nv + 1
      case(1)
        list_nt = list_nt + 1
      case(6)
        list_ne = list_ne + 1
      case default
        stop
      end select
      list_n = list_n + 1   
      tmpf =>  tmpf%next
    end do
    n=0
    if(list_nv.gt.0) n=n+1
    if(list_nt.gt.0) n=n+1
    if(list_ne.gt.0) n=n+1

    if(list_nv.gt.0) allocate(field_data_nv_send(maxdim1,list_nv*local_block%max_nv_send,local_block%nnb), source=0._r4)
    if(list_nv.gt.0) allocate(field_data_nv_recv(maxdim1,list_nv*local_block%max_nv_recv,local_block%nnb))
    if(list_nt.gt.0) allocate(field_data_nt_send(maxdim1,list_nt*local_block%max_nt_send,local_block%nnb), source=0._r4)
    if(list_nt.gt.0) allocate(field_data_nt_recv(maxdim1,list_nt*local_block%max_nt_recv,local_block%nnb))
    if(list_ne.gt.0) allocate(field_data_ne_send(maxdim1,list_ne*local_block%max_ne_send,local_block%nnb), source=0._r4)
    if(list_ne.gt.0) allocate(field_data_ne_recv(maxdim1,list_ne*local_block%max_ne_recv,local_block%nnb))

    !assign send_data value
    allocate(reqs_send(local_block%nnb * n))
    allocate(reqs_recv(local_block%nnb * n))
    allocate(status(MPI_STATUS_SIZE, local_block%nnb*n))
    ns=1; nr=1

    do c = 1, local_block%nnb
      dst_s = local_block%send_sites(c)%cpuid
      dst_r = local_block%recv_sites(c)%cpuid
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%send_sites(c)%nv
          sizedim1 = tmpf%dim1
!$omp parallel private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%send_sites(c)%local_v(iv)
            field_data_nv_send(1:sizedim1,iv+n_v*iv_list,c) = tmpf%field_data%f_r4(1:sizedim1,tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%send_sites(c)%nt
          sizedim1 = tmpf%dim1
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%send_sites(c)%local_t(it)
            field_data_nt_send(1:sizedim1,it+n_t*it_list,c) = tmpf%field_data%f_r4(1:sizedim1,tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%send_sites(c)%ne
          sizedim1 = tmpf%dim1
!$omp parallel private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%send_sites(c)%local_e(ie)
            field_data_ne_send(1:sizedim1,ie+n_e*ie_list,c) = tmpf%field_data%f_r4(1:sizedim1,tmp_n)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select
        tmpf => tmpf%next
      end do
      if (list_nv .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_nv_send(:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_nv_send(:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_nt .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_nt_send(:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_nt_send(:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_ne .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_ne_send(:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_ne_send(:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
    end do

    call MPI_WaitAll(n*local_block%nnb, reqs_send, status, ierr)
    call MPI_WaitAll(n*local_block%nnb, reqs_recv, status, ierr)

    !get recv_data value
    do c = 1, local_block%nnb
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0;
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%recv_sites(c)%nv
          sizedim1 = tmpf%dim1
!$omp parallel  private(iv,tmp_n) 
!$omp do schedule(dynamic,5)
          do iv = 1, n_v
            tmp_n = local_block%recv_sites(c)%local_v(iv)
            tmpf%field_data%f_r4(1:sizedim1,tmp_n) = field_data_nv_recv(1:sizedim1,iv+n_v*iv_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%recv_sites(c)%nt
          sizedim1 = tmpf%dim1
!$omp parallel  private(it,tmp_n) 
!$omp do schedule(dynamic,5)
          do it = 1, n_t
            tmp_n = local_block%recv_sites(c)%local_t(it)
            tmpf%field_data%f_r4(1:sizedim1,tmp_n) = field_data_nt_recv(1:sizedim1,it+n_t*it_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%recv_sites(c)%ne
          sizedim1 = tmpf%dim1
!$omp parallel  private(ie,tmp_n) 
!$omp do schedule(dynamic,5)
          do ie = 1, n_e
            tmp_n = local_block%recv_sites(c)%local_e(ie)
            tmpf%field_data%f_r4(1:sizedim1,tmp_n) = field_data_ne_recv(1:sizedim1,ie+n_e*ie_list,c)
          end do
!$omp end do nowait 
!$omp end parallel 
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select  
          tmpf => tmpf%next
      end do
    end do

    call exchange_data_2d_clean(field_head)

    if(.false.)  then
      do c = 1, local_block%nnb
         open(101, file=i2s(iblock)//".send.cell."//i2s(local_block%send_sites(c)%cpuid)//".txt")
         n_v = size(local_block%send_sites(c)%v)
         write(101, "(G20.10)", advance='yes') field_data_nv_send(:,1:n_v,c)
         close(101) 
         open(2, file=i2s(iblock)//".recv.cell."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
         n_v = size(local_block%recv_sites(c)%v)
         write(2, "(G20.10)", advance='yes') field_data_nv_recv(:,1:n_v,c)
         close(2)
         open(101, file=i2s(iblock)//".send.edge."//i2s(local_block%send_sites(c)%cpuid)//".txt")
         n_e = size(local_block%send_sites(c)%e)
         write(101, "(G20.10)", advance='yes') field_data_ne_send(:,1:n_e,c)
         close(101) 
         open(2, file=i2s(iblock)//".recv.edge."//i2s(local_block%recv_sites(c)%cpuid)//".txt")
         n_e = size(local_block%recv_sites(c)%e)
         write(2, "(G20.10)", advance='yes') field_data_ne_recv(:,1:n_e,c)
         close(2)
      end do
    end if

    !clean  
    if(list_nv .gt. 0) deallocate(field_data_nv_send)
    if(list_nv .gt. 0) deallocate(field_data_nv_recv)    
    if(list_nt .gt. 0) deallocate(field_data_nt_send)
    if(list_nt .gt. 0) deallocate(field_data_nt_recv)    
    if(list_ne .gt. 0) deallocate(field_data_ne_send)
    if(list_ne .gt. 0) deallocate(field_data_ne_recv)
    deallocate(reqs_send)
    deallocate(reqs_recv)
    deallocate(status)

  end subroutine exchange_data_2d_r4


  subroutine exchange_data_3d_add(mesh,field_head,scalar_field_data,mix_flag)
    !IO
    type(global_domain)                                  :: mesh
    type(exchange_field_list_3d), pointer, intent(inout) :: field_head
    type(scalar_3d_field)       , target,  intent(inout) :: scalar_field_data
    character(len=*), optional, intent(in)               :: mix_flag

    !local
    type(exchange_field_list_3d), pointer, save          :: field_list

    if(.not.associated(field_head))then
      allocate(field_head)
      field_list => field_head
    else if(.not.associated(field_list%next))then
      allocate(field_list%next)
      field_list => field_list%next
    else
      print*,"field_list error"
    end if

    field_list%field_data => scalar_field_data
    if(present(mix_flag).and.trim(mix_flag).eq."r4")then
       field_list%dim1 = size(scalar_field_data%f_r4, 1)
       field_list%dim2 = size(scalar_field_data%f_r4, 2)
       field_list%dim3 = size(scalar_field_data%f_r4, 3)
    else
       field_list%dim1 = size(scalar_field_data%f, 1)
       field_list%dim2 = size(scalar_field_data%f, 2)
       field_list%dim3 = size(scalar_field_data%f, 3)
    end if

    !IF(scalar_field_data%pos .eq. -1)then
      if((field_list%dim3 .eq. mesh%nv) .or. (field_list%dim3 .eq. mesh%nv_full)) then ! v location
        scalar_field_data%pos = 0
      else if((field_list%dim3 .eq. mesh%nt) .or. (field_list%dim3 .eq. mesh%nt_full)) then ! t location
        scalar_field_data%pos = 1
      else if((field_list%dim3 .eq. mesh%ne) .or. (field_list%dim3 .eq. mesh%ne_full)) then ! e location
        scalar_field_data%pos = 6
      else
        print*,"Invalid scalar field data type"
      end if
    !END IF

  end subroutine exchange_data_3d_add


  subroutine exchange_data_3d_clean(field_head)
    !IO
    type(exchange_field_list_3d), pointer, intent(inout) :: field_head

    !local
    type(exchange_field_list_3d), pointer                :: field_list,tmpf
    
    field_list => field_head
    do while(associated(field_list))
      tmpf => field_list%next
      deallocate(field_list)
      field_list => tmpf
    end do
    field_head => null()

  end subroutine exchange_data_3d_clean


  subroutine exchange_data_3d(local_block,field_head)
    !io
    type(block_structure),        target,  intent(in)    :: local_block
    type(exchange_field_list_3d), pointer, intent(inout) :: field_head

    !local
    integer, allocatable                                 :: reqs_send(:), reqs_recv(:), reqs(:)
    integer                                              :: ierr, c, iv, ie, it, dst_s, dst_r
    integer                                              :: i_list, iblock, iv_list, it_list, ie_list, n, nr, ns
    integer                                              :: list_n, tmp_n, n_v, n_t, n_e, list_nv, list_nt, list_ne
    integer                                              :: maxdim1, maxdim2, sizedim1, sizedim2

    INTEGER, allocatable                                 :: status(:,:) !STATUS(MPI_STATUS_SIZE)

    real(r8),allocatable                                 :: field_data_nv_send(:,:,:,:)!bdry
    real(r8),allocatable                                 :: field_data_ne_send(:,:,:,:)
    real(r8),allocatable                                 :: field_data_nt_send(:,:,:,:)

    real(r8),allocatable                                 :: field_data_nv_recv(:,:,:,:)!halo
    real(r8),allocatable                                 :: field_data_ne_recv(:,:,:,:)
    real(r8),allocatable                                 :: field_data_nt_recv(:,:,:,:)
    type(exchange_field_list_3d), pointer                :: tmpf

    tmpf => field_head
    maxdim1 = 0
    maxdim2 = 0
    do while(associated(tmpf))
      maxdim1 = max(maxdim1, tmpf%dim1)
      maxdim2 = max(maxdim2, tmpf%dim2)
      tmpf =>  tmpf%next
    end do

    list_n = 0
    iblock = mpi_rank()
    tmpf => field_head

    list_nv=0; list_nt=0; list_ne=0
    do while(associated(tmpf))
      select case(tmpf%field_data%pos)
      case(0)
        list_nv = list_nv + 1
      case(1)
        list_nt = list_nt + 1
      case(6)
        list_ne = list_ne + 1
      case default
        stop
      end select
      list_n = list_n + 1   
      tmpf =>  tmpf%next
    end do
    n=0
    if(list_nv.gt.0) n=n+1
    if(list_nt.gt.0) n=n+1
    if(list_ne.gt.0) n=n+1

    ! nlev  nlev+1??
    if(list_nv.gt.0) allocate(field_data_nv_send(maxdim1,maxdim2,list_nv*local_block%max_nv_send,local_block%nnb), source=0._r8)
    if(list_nv.gt.0) allocate(field_data_nv_recv(maxdim1,maxdim2,list_nv*local_block%max_nv_recv,local_block%nnb))
    if(list_nt.gt.0) allocate(field_data_nt_send(maxdim1,maxdim2,list_nt*local_block%max_nt_send,local_block%nnb), source=0._r8)
    if(list_nt.gt.0) allocate(field_data_nt_recv(maxdim1,maxdim2,list_nt*local_block%max_nt_recv,local_block%nnb))
    if(list_ne.gt.0) allocate(field_data_ne_send(maxdim1,maxdim2,list_ne*local_block%max_ne_send,local_block%nnb), source=0._r8)
    if(list_ne.gt.0) allocate(field_data_ne_recv(maxdim1,maxdim2,list_ne*local_block%max_ne_recv,local_block%nnb))

    !assign send_data value
    allocate(reqs_send(local_block%nnb * n))
    allocate(reqs_recv(local_block%nnb * n))
    allocate(status(MPI_STATUS_SIZE, local_block%nnb*n))
    ns=1; nr=1
    do c = 1, local_block%nnb
      dst_s = local_block%send_sites(c)%cpuid
      dst_r = local_block%recv_sites(c)%cpuid
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%send_sites(c)%nv
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do iv = 1, n_v
            tmp_n = local_block%send_sites(c)%local_v(iv)
            field_data_nv_send(1:sizedim1,1:sizedim2,iv+n_v*iv_list,c) = tmpf%field_data%f(1:sizedim1,1:sizedim2,tmp_n)
          end do
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%send_sites(c)%nt
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do it = 1, n_t
            tmp_n = local_block%send_sites(c)%local_t(it)
            field_data_nt_send(1:sizedim1,1:sizedim2,it+n_t*it_list,c) = tmpf%field_data%f(1:sizedim1,1:sizedim2,tmp_n)
          end do
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%send_sites(c)%ne
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do ie = 1, n_e
            tmp_n = local_block%send_sites(c)%local_e(ie)
            field_data_ne_send(1:sizedim1,1:sizedim2,ie+n_e*ie_list,c) = tmpf%field_data%f(1:sizedim1,1:sizedim2,tmp_n)
          end do
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select
        tmpf => tmpf%next
      end do
      if (list_nv .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_nv_send(:,:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_nv_send(:,:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_nt .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_nt_send(:,:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_nt_send(:,:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_ne .gt. 0) then
        if (r8 .eq. 8) then
          call MPI_Isend(field_data_ne_send(:,:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r8 .eq. 4) then
          call MPI_Isend(field_data_ne_send(:,:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
    end do

    call MPI_WaitAll(n*local_block%nnb, reqs_send, status, ierr)
    call MPI_WaitAll(n*local_block%nnb, reqs_recv, status, ierr)

    !get recv_data value
    do c = 1, local_block%nnb
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0;
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%recv_sites(c)%nv
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do iv = 1, n_v
            tmp_n = local_block%recv_sites(c)%local_v(iv)
            tmpf%field_data%f(1:sizedim1,1:sizedim2,tmp_n) = field_data_nv_recv(1:sizedim1,1:sizedim2,iv+n_v*iv_list,c)
          end do
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%recv_sites(c)%nt
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do it = 1, n_t
            tmp_n = local_block%recv_sites(c)%local_t(it)
            tmpf%field_data%f(1:sizedim1,1:sizedim2,tmp_n) = field_data_nt_recv(1:sizedim1,1:sizedim2,it+n_t*it_list,c)
          end do
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%recv_sites(c)%ne
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do ie = 1, n_e
            tmp_n = local_block%recv_sites(c)%local_e(ie)
            tmpf%field_data%f(1:sizedim1,1:sizedim2,tmp_n) = field_data_ne_recv(1:sizedim1,1:sizedim2,ie+n_e*ie_list,c)
          end do
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select  
          tmpf => tmpf%next
      end do
    end do

    call exchange_data_3d_clean(field_head)

    !clean  
    if(list_nv .gt. 0) deallocate(field_data_nv_send)
    if(list_nv .gt. 0) deallocate(field_data_nv_recv)    
    if(list_nt .gt. 0) deallocate(field_data_nt_send)
    if(list_nt .gt. 0) deallocate(field_data_nt_recv)    
    if(list_ne .gt. 0) deallocate(field_data_ne_send)
    if(list_ne .gt. 0) deallocate(field_data_ne_recv)
    deallocate(reqs_send)
    deallocate(reqs_recv)
    deallocate(status)

  end subroutine exchange_data_3d

  subroutine exchange_data_3d_r4(local_block,field_head)
    !io
    type(block_structure),        target,  intent(in)    :: local_block
    type(exchange_field_list_3d), pointer, intent(inout) :: field_head

    !local
    integer, allocatable                                 :: reqs_send(:), reqs_recv(:), reqs(:)
    integer                                              :: ierr, c, iv, ie, it, dst_s, dst_r
    integer                                              :: i_list, iblock, iv_list, it_list, ie_list, n, nr, ns
    integer                                              :: list_n, tmp_n, n_v, n_t, n_e, list_nv, list_nt, list_ne
    integer                                              :: maxdim1, maxdim2, sizedim1, sizedim2

    INTEGER, allocatable                                 :: status(:,:) !STATUS(MPI_STATUS_SIZE)

    real(r4),allocatable                                 :: field_data_nv_send(:,:,:,:)!bdry
    real(r4),allocatable                                 :: field_data_ne_send(:,:,:,:)
    real(r4),allocatable                                 :: field_data_nt_send(:,:,:,:)

    real(r4),allocatable                                 :: field_data_nv_recv(:,:,:,:)!halo
    real(r4),allocatable                                 :: field_data_ne_recv(:,:,:,:)
    real(r4),allocatable                                 :: field_data_nt_recv(:,:,:,:)
    type(exchange_field_list_3d), pointer                :: tmpf

    tmpf => field_head
    maxdim1 = 0
    maxdim2 = 0
    do while(associated(tmpf))
      maxdim1 = max(maxdim1, tmpf%dim1)
      maxdim2 = max(maxdim2, tmpf%dim2)
      tmpf =>  tmpf%next
    end do

    list_n = 0
    iblock = mpi_rank()
    tmpf => field_head

    list_nv=0; list_nt=0; list_ne=0
    do while(associated(tmpf))
      select case(tmpf%field_data%pos)
      case(0)
        list_nv = list_nv + 1
      case(1)
        list_nt = list_nt + 1
      case(6)
        list_ne = list_ne + 1
      case default
        stop
      end select
      list_n = list_n + 1   
      tmpf =>  tmpf%next
    end do
    n=0
    if(list_nv.gt.0) n=n+1
    if(list_nt.gt.0) n=n+1
    if(list_ne.gt.0) n=n+1

    ! nlev  nlev+1??
    if(list_nv.gt.0) allocate(field_data_nv_send(maxdim1,maxdim2,list_nv*local_block%max_nv_send,local_block%nnb), source=0._r4)
    if(list_nv.gt.0) allocate(field_data_nv_recv(maxdim1,maxdim2,list_nv*local_block%max_nv_recv,local_block%nnb))
    if(list_nt.gt.0) allocate(field_data_nt_send(maxdim1,maxdim2,list_nt*local_block%max_nt_send,local_block%nnb), source=0._r4)
    if(list_nt.gt.0) allocate(field_data_nt_recv(maxdim1,maxdim2,list_nt*local_block%max_nt_recv,local_block%nnb))
    if(list_ne.gt.0) allocate(field_data_ne_send(maxdim1,maxdim2,list_ne*local_block%max_ne_send,local_block%nnb), source=0._r4)
    if(list_ne.gt.0) allocate(field_data_ne_recv(maxdim1,maxdim2,list_ne*local_block%max_ne_recv,local_block%nnb))

    !assign send_data value
    allocate(reqs_send(local_block%nnb * n))
    allocate(reqs_recv(local_block%nnb * n))
    allocate(status(MPI_STATUS_SIZE, local_block%nnb*n))
    ns=1; nr=1
    do c = 1, local_block%nnb
      dst_s = local_block%send_sites(c)%cpuid
      dst_r = local_block%recv_sites(c)%cpuid
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%send_sites(c)%nv
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do iv = 1, n_v
            tmp_n = local_block%send_sites(c)%local_v(iv)
            field_data_nv_send(1:sizedim1,1:sizedim2,iv+n_v*iv_list,c) = tmpf%field_data%f_r4(1:sizedim1,1:sizedim2,tmp_n)
          end do
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%send_sites(c)%nt
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do it = 1, n_t
            tmp_n = local_block%send_sites(c)%local_t(it)
            field_data_nt_send(1:sizedim1,1:sizedim2,it+n_t*it_list,c) = tmpf%field_data%f_r4(1:sizedim1,1:sizedim2,tmp_n)
          end do
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%send_sites(c)%ne
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do ie = 1, n_e
            tmp_n = local_block%send_sites(c)%local_e(ie)
            field_data_ne_send(1:sizedim1,1:sizedim2,ie+n_e*ie_list,c) = tmpf%field_data%f_r4(1:sizedim1,1:sizedim2,tmp_n)
          end do
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select
        tmpf => tmpf%next
      end do
      if (list_nv .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_nv_send(:,:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_nv_send(:,:,:,c), &
               local_block%send_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nv_recv(:,:,:,c), &
               local_block%recv_sites(c)%nv*list_nv*maxdim2*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_nt .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_nt_send(:,:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_nt_send(:,:,:,c), &
               local_block%send_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_nt_recv(:,:,:,c), &
               local_block%recv_sites(c)%nt*list_nt*maxdim2*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
      if (list_ne .gt. 0) then
        if (r4 .eq. 8) then
          call MPI_Isend(field_data_ne_send(:,:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL8, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL8, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        elseif (r4 .eq. 4) then
          call MPI_Isend(field_data_ne_send(:,:,:,c), &
               local_block%send_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL, &
               dst_s, 102, local_block%comm, reqs_send(ns), ierr)
          call MPI_Irecv(field_data_ne_recv(:,:,:,c), &
               local_block%recv_sites(c)%ne*list_ne*maxdim2*maxdim1, MPI_REAL, &
               dst_r, 102, local_block%comm, reqs_recv(nr), ierr)
        end if
        ns=ns+1; nr=nr+1
      endif
    end do

    call MPI_WaitAll(n*local_block%nnb, reqs_send, status, ierr)
    call MPI_WaitAll(n*local_block%nnb, reqs_recv, status, ierr)

    !get recv_data value
    do c = 1, local_block%nnb
      tmpf => field_head
      iv_list=0; it_list=0; ie_list=0;
      do i_list = 0, list_n-1
        select case(tmpf%field_data%pos)
        !nv
        case(0)
          n_v = local_block%recv_sites(c)%nv
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do iv = 1, n_v
            tmp_n = local_block%recv_sites(c)%local_v(iv)
            tmpf%field_data%f_r4(1:sizedim1,1:sizedim2,tmp_n) = field_data_nv_recv(1:sizedim1,1:sizedim2,iv+n_v*iv_list,c)
          end do
          iv_list = iv_list + 1
        !nt
        case(1)
          n_t = local_block%recv_sites(c)%nt
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do it = 1, n_t
            tmp_n = local_block%recv_sites(c)%local_t(it)
            tmpf%field_data%f_r4(1:sizedim1,1:sizedim2,tmp_n) = field_data_nt_recv(1:sizedim1,1:sizedim2,it+n_t*it_list,c)
          end do
          it_list = it_list + 1
        !ne
        case(6)
          n_e = local_block%recv_sites(c)%ne
          sizedim1 = tmpf%dim1
          sizedim2 = tmpf%dim2
          do ie = 1, n_e
            tmp_n = local_block%recv_sites(c)%local_e(ie)
            tmpf%field_data%f_r4(1:sizedim1,1:sizedim2,tmp_n) = field_data_ne_recv(1:sizedim1,1:sizedim2,ie+n_e*ie_list,c)
          end do
          ie_list = ie_list + 1
        case default
          print*,"Invalid scalar field data type"
        end select  
          tmpf => tmpf%next
      end do
    end do

    call exchange_data_3d_clean(field_head)

    !clean  
    if(list_nv .gt. 0) deallocate(field_data_nv_send)
    if(list_nv .gt. 0) deallocate(field_data_nv_recv)    
    if(list_nt .gt. 0) deallocate(field_data_nt_send)
    if(list_nt .gt. 0) deallocate(field_data_nt_recv)    
    if(list_ne .gt. 0) deallocate(field_data_ne_send)
    if(list_ne .gt. 0) deallocate(field_data_ne_recv)
    deallocate(reqs_send)
    deallocate(reqs_recv)
    deallocate(status)

  end subroutine exchange_data_3d_r4

  subroutine debug_data_1d(timestep,g_index,varsize,varname,vardata)
    implicit none
    !IO
    integer(i4), intent(in)      :: timestep
    integer(i4), intent(in)      :: g_index(:)
    integer(i4), intent(in)      :: varsize
    character(len=*), intent(in) :: varname
    real(r8),intent(in)          :: vardata(:)

    !local
    integer            :: i,ilevl
    character(len=100) :: filepath

    filepath = "../varData/"

    open(1,file=trim(filepath)//"Data-"//i2s(mpi_rank())&
         //"-"//trim(varname)//"-"//i2s(timestep)//".txt")
      do i=1,varsize
        write(1,'(ES25.15)') vardata(i)
      end do
      close(1)

      open(1,file=trim(filepath)//"Index-"//i2s(mpi_rank())&
           //"-"//trim(varname)//"-"//i2s(timestep)//".txt")
      do i=1,varsize
        write(1,'(I0)') g_index(i)
      end do
    close(1)

  end subroutine debug_data_1d

  
  subroutine debug_data_2d(timestep,nlevl,g_index,varsize,varname,vardata)
    implicit none
    !IO
    integer(i4), intent(in)      :: timestep
    integer(i4), intent(in)      :: nlevl
    integer(i4), intent(in)      :: g_index(:)
    integer(i4), intent(in)      :: varsize
    character(len=*), intent(in) :: varname
    real(r8),intent(in)          :: vardata(:,:)

    !local
    integer            :: i,ilevl
    character(len=100) :: filepath

    filepath = "../varData/"

    do ilevl=1,nlevl
    open(1,file=trim(filepath)//"Data-"//i2s(mpi_rank())&
         //"-"//trim(varname)//"-levl-"//i2s(ilevl)//"-"//i2s(timestep)//".txt")
      do i=1,varsize
         write(1,'(ES25.15)') vardata(ilevl,i)
      end do
      close(1)

      open(1,file=trim(filepath)//"Index-"//i2s(mpi_rank())&
           //"-"//trim(varname)//"-levl-"//i2s(ilevl)//"-"//i2s(timestep)//".txt")
      do i=1,varsize
         write(1,'(I0)') g_index(i)
      end do
      close(1)
    end do

  end subroutine debug_data_2d

  
  subroutine show_basic_info(mesh,flag_rank)
    implicit none
    !IO
    type(global_domain_data), intent(in) :: mesh
    logical,                  intent(in) :: flag_rank

    !local
    logical                              :: enable

    call check_process(mesh%nv, mpi_size())

    !if(flag_rank) then
    !   if(mpi_rank() == 0) then
    !     enable = .true.
    !   else
    !     enable = .false.
    !   end if

    !  if(enable) then
    !    print*,"========================================"
    !    print*,"  Mpi process Num: ",mpi_size()
    !    print*,"========================================"
    !  end if
    !endif

  end subroutine show_basic_info


  subroutine check_process(nv, mpisize)
    integer(i4), intent(in)         :: nv
    integer,     intent(in)         :: mpisize

    if(nv/mpisize.lt.nv_per_core) then
      if(mpi_rank() == 0) then
          print*, "the number of process is too large! you can try small number"
      endif
      stop
    endif
    if(mpisize.lt.2) then
      print*, "the number of process must larger than 1"
      stop
    endif
  end subroutine

end module grist_config_partition
