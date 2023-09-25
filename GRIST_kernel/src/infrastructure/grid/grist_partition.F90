!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: assign mesh data
!
! Revision history:
!----------------------------------------------------------------------------

!================================================
! PAR_GRIST partition the global mesh
! Zhuang Liu, 2020, SEP
!================================================

module partition_routine
  use grist_constants_dbl,             only: i4, r8, i8, pi
  use grist_domain_types,              only: global_domain_data
  use grist_lib
  use grist_nml_module,                only: index_flag, stencil_width, gridFileNameHead, gridFilePath, nml_set_global_vars
  use grist_clocks,                    only: clock_begin, clock_end

  type index_list
      integer(i4), allocatable         :: v_index(:)
      integer(i4), allocatable         :: t_index(:)
      integer(i4), allocatable         :: e_index(:)
      integer(i4)                      :: nv, nt, ne, &
                                          nv_iner, nt_iner, ne_iner, &
                                          nv_compute, nt_compute, ne_compute, &
                                          nv_full, nt_full, ne_full
      integer(i4), allocatable         :: nv_bdry(:), nt_bdry(:), ne_bdry(:), &
                                          nv_halo(:), nt_halo(:), ne_halo(:)
  end type index_list

  type(index_list), allocatable        :: domain_compute(:)
  type(index_list), allocatable, target:: domain_halo(:,:)
  type(set),        allocatable        :: compute_v_list(:), bdry1_t_list(:), &
                                          bdry1_e_list(:), halo1_e_list(:)
  type(index_list), allocatable, target:: domain_full(:)

  !for LAM: record the number of cells in the compute and relax regions
  integer                              :: region_cells_compute
  integer,          parameter          :: maxi4 = huge(i4)
  integer(i4),      allocatable        :: v_index_nomask(:)
contains
  !> partition the mesh, partition flag for each node
  !> is saved in mesh%parts
  subroutine decompose(mesh, nprocs)
    implicit none
    type(global_domain_data), intent(inout) :: mesh
    integer,                  intent(in)    :: nprocs
    integer(i4),              allocatable   :: xadj(:), adjncy(:)
    integer(i4)                             :: i, j, x, act_nnb, c_idx, nb_idx, mark
    integer(i4),              allocatable   :: parts_metis(:)
    type(set)                               :: nomask
    type(map)                               :: map_nomask

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
    call metis_decomp(mesh%nv, adjncy, xadj, nprocs, mesh%parts)
#endif

    deallocate(xadj, adjncy)

  end subroutine decompose

  subroutine init_compute_index(mesh, nprocs)
    implicit none
    type(global_domain_data), intent(in) :: mesh
    integer,                  intent(in) :: nprocs

    ! local variables
    integer                              :: iproc
    integer(i4)                          :: iv, it_local, idx, iv_global
    type(set), allocatable               :: compute_t_list(:), &
                                            compute_e_list(:)

    integer(i4)                          :: head, fina, nnb, nb, f_head,  h_index, &
                                            i, j, curr, minx, maxx, miny, maxy, maxi, maxj, &
                                            minsquare, temp1, temp2, cx, cy, diffindex, precs,&
                                            square, ztmpid, x, y, z, hilbert_pos
    logical                              :: flag, flag1, flag2
    integer(i4),   allocatable           :: hilbert_cartesian(:,:),queue(:),z_index(:)
    integer(i4), allocatable             :: binary_code(:,:)    
    character(len = 5)                   :: ztmp 
    character(len = 5), allocatable      :: z_list(:)
    character(len = 1), allocatable      :: geo_hash(:)
    real(8)                              :: lon, lat, lond, latd
    type(set)                            :: bfs_index, hilbert_index, morton_index

    allocate(domain_compute(nprocs))

    allocate(compute_v_list(nprocs))
    allocate(compute_t_list(nprocs))
    allocate(compute_e_list(nprocs))

#ifdef LAM_DOMAIN
    do iv = 1, region_cells_compute
       curr = v_index_nomask(iv)
       call compute_v_list(mesh%parts(curr)+1)%insert(curr)
    end do
#else
    do iv = 1, mesh%nv
       call compute_v_list(mesh%parts(iv)+1)%insert(iv)
    end do
#endif

    do iproc = 1, nprocs
       call compute_v_list(iproc)%dump(domain_compute(iproc)%v_index)
       domain_compute(iproc)%nv = size(domain_compute(iproc)%v_index)
       !===================================
       ! reordering index of compute domain
       !===================================
       SELECT CASE(index_flag)
       CASE('bfs')
          if(iproc == 1) print*, "reordering index using method: [BFS]"
          allocate(queue(domain_compute(iproc)%nv))
          do i = 1, domain_compute(iproc)%nv
              queue(i) = -1
          enddo

          queue(1) = domain_compute(iproc)%v_index(1)
          head = 1
          fina = 2
          do while(head .ne. fina)
              call bfs_index%insert(queue(head))
              nnb = mesh%vtx_nnb(queue(head))
              do iv = 1, nnb
                  nb = mesh%vtx_nb%f(queue(head), iv)
                  flag1 = compute_v_list(iproc)%find(nb)
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

          do iv = 1, domain_compute(iproc)%nv
             ! Insert isolated compute point
             if (.not. bfs_index%find(domain_compute(iproc)%v_index(iv))) then
                call bfs_index%insert(domain_compute(iproc)%v_index(iv))
             end if
          end do
          call bfs_index%dump(domain_compute(iproc)%v_index)

          call bfs_index%final
          deallocate(queue)
          if(iproc == 1) print*,"index reordering by [BFS] end"

       CASE('hilbert')
          ! Currently, the algorithm used in the Hilbert curve only
          ! sort the vertices with 6 neighbors
          if(iproc == 1) print*, "reordering index using SFC: [Hilbert]"
          allocate(queue(domain_compute(iproc)%nv))
          allocate(hilbert_cartesian(5, domain_compute(iproc)%nv))
          head = 1
          fina = head + 1
          hilbert_cartesian = 0
          ! Find the first vtx with 6 neighbors
          do iv = 1, domain_compute(iproc)%nv
             if (mesh%vtx_nnb(domain_compute(iproc)%v_index(iv)) == 6) then
                queue(head) = domain_compute(iproc)%v_index(iv)
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
                   flag1 = compute_v_list(iproc)%find(nb)
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
          do i = 1, domain_compute(iproc)%nv
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
          do i = 1, domain_compute(iproc)%nv
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
          allocate(binary_code(2, domain_compute(iproc)%nv))
          square = 2**minsquare
          do i = 1, domain_compute(iproc)%nv
             cx = hilbert_cartesian(2,i)
             cy = hilbert_cartesian(3,i)
             hilbert_pos = 0
             call xy2d(square, cx, cy, hilbert_pos)
             binary_code(1,i) = hilbert_cartesian(1,i)
             binary_code(2,i) = hilbert_pos
          enddo
          ! Sort the vtx indices according to the Hilbert pos
          do i = 1, domain_compute(iproc)%nv - 1
             do j = 1, domain_compute(iproc)%nv - i
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
          do i = 1, domain_compute(iproc)%nv
             ! Avoid 0 index (default for isolated compute vtx)
             if(binary_code(1,i) .ne. 0) call hilbert_index%insert(binary_code(1,i))
          enddo
          do i = 1, domain_compute(iproc)%nv
             ! Push back isolated compute vtx
             if (.not. hilbert_index%find(domain_compute(iproc)%v_index(i))) then 
                call hilbert_index%insert(domain_compute(iproc)%v_index(i))
             end if
          end do
          call hilbert_index%dump(domain_compute(iproc)%v_index)

          call hilbert_index%final
          deallocate(queue)
          deallocate(hilbert_cartesian)
          deallocate(binary_code)
          if(iproc == 1) print*,"index reordering by [Hilbert] end"

       CASE('morton')
          if(iproc == 1) print*, "reordering index using SFC: [Morton]"
          precs = 5
          allocate(z_index(domain_compute(iproc)%nv))
          allocate(z_list(domain_compute(iproc)%nv))
          allocate(geo_hash(1:precs))
          ! Get hash code for each vtx (z_list)
          do i = 1, domain_compute(iproc)%nv
             iv_global = domain_compute(iproc)%v_index(i)
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
          do i = 1, domain_compute(iproc)%nv - 1
             do j = 1, domain_compute(iproc)%nv - i
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
          do i = 1, domain_compute(iproc)%nv
             call morton_index%insert(z_index(i))
          enddo
          call morton_index%dump(domain_compute(iproc)%v_index)

          call morton_index%final
          deallocate(z_index)
          deallocate(z_list)
          deallocate(geo_hash)
          if(iproc == 1) print*,"index reordering by [Morton] end"

       CASE DEFAULT
          if(iproc == 1) print*,"no selection of SFC: [random] "
       END SELECT
    end do

#ifdef LAM_DOMAIN
    do iv = 1, region_cells_compute
       curr = v_index_nomask(iv)
       do it_local = 1, mesh%vtx_nnb(curr)

          idx = mesh%vtx_tr%f(curr, it_local)
          call compute_t_list(mesh%parts(curr)+1)%insert(idx)

          idx = mesh%vtx_ed%f(curr, it_local)
          call compute_e_list(mesh%parts(curr)+1)%insert(idx)

       end do
    end do
#else
    do iv = 1, mesh%nv
       do it_local = 1, mesh%vtx_nnb(iv)

          idx = mesh%vtx_tr%f(iv, it_local)
          call compute_t_list(mesh%parts(iv)+1)%insert(idx)

          idx = mesh%vtx_ed%f(iv, it_local)
          call compute_e_list(mesh%parts(iv)+1)%insert(idx)

       end do
    end do
#endif

    ! from linked list to array
    do iproc = 1, nprocs
       call compute_t_list(iproc)%dump(domain_compute(iproc)%t_index)
       call compute_e_list(iproc)%dump(domain_compute(iproc)%e_index)
       domain_compute(iproc)%nt = size(domain_compute(iproc)%t_index)
       domain_compute(iproc)%ne = size(domain_compute(iproc)%e_index)

       call compute_t_list(iproc)%final
       call compute_e_list(iproc)%final
    end do

    call clock_begin("write compute")
    ! write index
    open(101, file="compute-"//i2s(nprocs), form='unformatted', access="stream")
    do iproc = 1, nprocs
       write(101) domain_compute(iproc)%nv
    end do
    do iproc = 1, nprocs
       write(101) domain_compute(iproc)%v_index
    end do
    !!test
    !do iproc = 1, nprocs
    !   open(101, file="compute-"//i2s(iproc)//".txt", form='formatted')
    !   write(101, *) domain_compute(iproc)%nv
    !   do iv = 1, domain_compute(iproc)%nv
    !      write(101, *) domain_compute(iproc)%v_index(iv)
    !   end do
    !   close(101)
    !end do
    call clock_end("write compute")

#ifdef LAM_DOMAIN
    deallocate(v_index_nomask)
#endif
    deallocate(compute_t_list)
    deallocate(compute_e_list)
    do iproc = 1, nprocs
       deallocate(domain_compute(iproc)%t_index, domain_compute(iproc)%e_index)
    end do

  end subroutine init_compute_index 

  subroutine init_bdry_halo_index(mesh, nprocs)
    implicit none
    type(global_domain_data), intent(inout) :: mesh
    integer,                  intent(in)    :: nprocs

    integer                                 :: iproc, sw, sw2
    integer(i4)                             :: iv, iv_local, iv_global, it_local, idx, cnt, mark
    logical                                 :: flag

    type(set), allocatable, target          :: bdry_v_list(:), bdry_t_list(:), bdry_e_list(:)
    type(set), allocatable, target          :: halo_v_list(:), halo_t_list(:), halo_e_list(:)
    type(set), target                       :: iner_v_list, iner_t_list, iner_e_list
    type(set), target                       :: full_t_list, full_e_list
    type(set), pointer                      :: c_list, i_list, o_list
    type(set), pointer                      :: b1, h1
    type(index_list), allocatable, target   :: domain_bdry(:,:)
    type(index_list), allocatable, target   :: domain_iner(:)
    type(index_list), pointer               :: c_domain, i_domain, o_domain

    allocate(bdry_v_list(stencil_width), bdry_t_list(stencil_width), bdry_e_list(stencil_width))
    allocate(halo_v_list(stencil_width), halo_t_list(stencil_width), halo_e_list(stencil_width))

    allocate(domain_bdry(stencil_width, nprocs))
    allocate(domain_halo(stencil_width, nprocs))
    allocate(domain_iner(nprocs))
    allocate(domain_full(nprocs))

#ifdef LAM_DOMAIN
    ! for init_send_recv
    allocate(bdry1_t_list(nprocs), bdry1_e_list(nprocs), halo1_e_list(nprocs))
    ! determine partitions for cells in the specific region (halo1 first)
    do iproc = 1, nprocs
       do iv = 1, domain_compute(iproc)%nv
          iv_global  = domain_compute(iproc)%v_index(iv)
          do iv_local = 1, mesh%vtx_nnb(iv_global)
             idx = mesh%vtx_nb%f(iv_global, iv_local)
             mark = mesh%plg_mark%f(idx, 1)
             if (mark .eq. 177) then
                mesh%parts(idx) = min(mesh%parts(idx), iproc-1)
             end if
          end do
       end do
    end do
#endif
    !find halo 1 and bdry 1
    do iproc = 1, nprocs
       b1 => bdry_v_list(1)
       h1 => halo_v_list(1)
       do iv = 1, domain_compute(iproc)%nv
          iv_global  = domain_compute(iproc)%v_index(iv)
          do iv_local = 1, mesh%vtx_nnb(iv_global)

             idx = mesh%vtx_nb%f(iv_global, iv_local)
             flag = compute_v_list(iproc)%find(idx)

             ! if not in, then this cell belongs to bdry1, this nb belongs halo1          
             if(.not.flag)then 
                call b1%insert(iv_global)
                call h1%insert(idx)
             end if
          end do
       end do
       call b1%dump(domain_bdry(1, iproc)%v_index)
       call h1%dump(domain_halo(1, iproc)%v_index)
       domain_bdry(1, iproc)%nv = size(domain_bdry(1, iproc)%v_index)
       domain_halo(1, iproc)%nv = size(domain_halo(1, iproc)%v_index)
       ! find bdry(sw) for sw > 1
       do sw = 2, stencil_width
          i_list => bdry_v_list(sw)
          c_list => bdry_v_list(sw-1)
          i_domain => domain_bdry(sw,   iproc)
          c_domain => domain_bdry(sw-1, iproc)

          if(sw == 2) then
             o_list => h1
          else
             o_list => bdry_v_list(sw-2)
          end if

          do iv = 1, c_domain%nv
             iv_global = c_domain%v_index(iv)
             do iv_local = 1, mesh%vtx_nnb(iv_global)
                idx = mesh%vtx_nb%f(iv_global, iv_local)
                flag = (.not. c_list%find(idx)) .and. (.not. o_list%find(idx))

                if(flag) then
                   call i_list%insert(idx)
                end if
             end do
          end do
          call i_list%dump(i_domain%v_index)
          i_domain%nv = size(i_domain%v_index)
       end do
       ! find halo(sw) for sw > 1    
       do sw = 2, stencil_width
          i_list => halo_v_list(sw)
          c_list => halo_v_list(sw-1)
          i_domain => domain_halo(sw,   iproc)
          c_domain => domain_halo(sw-1, iproc)

          if(sw == 2) then
             o_list => b1
          else
             o_list => halo_v_list(sw-2)
          end if

          do iv = 1, c_domain%nv
             iv_global = c_domain%v_index(iv)
             do iv_local = 1, mesh%vtx_nnb(iv_global)
                idx = mesh%vtx_nb%f(iv_global, iv_local)
#ifdef LAM_DOMAIN
             if (idx .gt. 0) then
                ! determine partitions for cells in the specific region (halo2+)
                mark = mesh%plg_mark%f(idx, 1)
                if (mark .eq. 177 .and. mesh%parts(idx) .eq. maxi4) then
                   mesh%parts(idx) = min(mesh%parts(idx), iproc-1)
                end if
#endif
                flag = (.not. c_list%find(idx)) .and. (.not. o_list%find(idx))

                if(flag) then
                   call i_list%insert(idx)
                end if
#ifdef LAM_DOMAIN
             end if
#endif
             end do
          end do
          call i_list%dump(i_domain%v_index)
          i_domain%nv = size(i_domain%v_index)
       end do

#ifdef LAM_DOMAIN
       c_domain => domain_bdry(1, iproc)
       do iv = 1, c_domain%nv
          iv_global = c_domain%v_index(iv)
          do it_local = 1, mesh%vtx_nnb(iv_global)
             call bdry1_t_list(iproc)%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
             call bdry1_e_list(iproc)%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
          end do
       end do
       c_domain => domain_halo(sw-1, iproc)
       do iv = 1, c_domain%nv
          iv_global = c_domain%v_index(iv)
          do it_local = 1, mesh%vtx_nnb(iv_global)
             call halo1_e_list(iproc)%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
          end do
       end do
#endif

#ifdef LAM_DOMAIN
       ! acquire t and e based global index
       do sw = 1, stencil_width

       !   c_domain => domain_bdry(sw, iproc)
       !   do iv = 1, c_domain%nv
       !      iv_global = c_domain%v_index(iv)
       !      do it_local = 1, mesh%vtx_nnb(iv_global)
       !         call bdry_t_list(sw)%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
       !         call bdry_e_list(sw)%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
       !      end do
       !   end do
       !   call bdry_t_list(sw)%dump(c_domain%t_index)
       !   call bdry_e_list(sw)%dump(c_domain%e_index)
       !   c_domain%nt = size(c_domain%t_index)
       !   c_domain%ne = size(c_domain%e_index)

          c_domain => domain_halo(sw, iproc)
          do iv = 1, c_domain%nv
             iv_global = c_domain%v_index(iv)
             do it_local = 1, mesh%vtx_nnb(iv_global)
                call halo_t_list(sw)%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
                call halo_e_list(sw)%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
             end do
          end do
          call halo_t_list(sw)%dump(c_domain%t_index)
          call halo_e_list(sw)%dump(c_domain%e_index)
          c_domain%nt = size(c_domain%t_index)
          c_domain%ne = size(c_domain%e_index)
       end do
#endif
       !=================================
       ! configure the inner domain
       !=================================
#ifdef USE_INER
       c_domain => domain_iner(iproc)
       do iv = 1, domain_compute(iproc)%nv
          iv_global = domain_compute(iproc)%v_index(iv)

          flag = .false.
          do sw = 1, stencil_width
             flag = flag .or. bdry_v_list(sw)%find(iv_global)
          end do

          if(.not.flag) then
             call iner_v_list%insert(iv_global)
          end if

       end do

       call iner_v_list%dump(c_domain%v_index)
       c_domain%nv = size(c_domain%v_index)
       !open(1, position='Append',file="z_index_v."//i2s(iproc-1)//".txt")
       !do i = 1, domain_iner(iproc)%nv
       !   write(1, "(I0)", advance='yes') domain_iner(iproc)%v_index(i)
       !enddo
       !close(1)

       ! t & e
       do iv = 1, domain_iner(iproc)%nv
          iv_global = domain_iner(iproc)%v_index(iv)
          do it_local = 1, mesh%vtx_nnb(iv_global)
             call iner_t_list%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
             call iner_e_list%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
          end do
       end do
       call iner_t_list%dump(domain_iner(iproc)%t_index)
       call iner_e_list%dump(domain_iner(iproc)%e_index)
       domain_iner(iproc)%nt = size(domain_iner(iproc)%t_index)
       domain_iner(iproc)%ne = size(domain_iner(iproc)%e_index)

       call iner_v_list%final
       call iner_t_list%final
       call iner_e_list%final
       deallocate(domain_iner(iproc)%v_index)
       deallocate(domain_iner(iproc)%t_index)
       deallocate(domain_iner(iproc)%e_index)
#endif
       ! deallocate
       do sw = 1, stencil_width
          call bdry_v_list(sw)%final
          call bdry_t_list(sw)%final
          call bdry_e_list(sw)%final
          call halo_v_list(sw)%final
          call halo_t_list(sw)%final
          call halo_e_list(sw)%final
       end do

       !=================================
       ! configure the full domain
       !=================================
#ifdef USE_INER
       domain_full(iproc)%nv_iner = domain_iner(iproc)%nv
       domain_full(iproc)%nt_iner = domain_iner(iproc)%nt
       domain_full(iproc)%ne_iner = domain_iner(iproc)%ne

       domain_full(iproc)%nv = domain_iner(iproc)%nv
       domain_full(iproc)%nt = domain_iner(iproc)%nt
       domain_full(iproc)%ne = domain_iner(iproc)%ne
#endif

       allocate(domain_full(iproc)%nv_bdry(sw))
       allocate(domain_full(iproc)%nv_halo(sw))
       allocate(domain_full(iproc)%nt_bdry(sw))
       allocate(domain_full(iproc)%nt_halo(sw))
       allocate(domain_full(iproc)%ne_bdry(sw))
       allocate(domain_full(iproc)%ne_halo(sw))

       domain_full(iproc)%nv_compute = domain_compute(iproc)%nv
       domain_full(iproc)%nv_bdry(1) = domain_full(iproc)%nv_compute
       do sw = 2, stencil_width
          domain_full(iproc)%nv_bdry(sw) = domain_full(iproc)%nv_bdry(sw-1)-domain_bdry(sw-1, iproc)%nv
       end do

       domain_full(iproc)%nv_halo(1) = domain_full(iproc)%nv_compute+domain_halo(1,iproc)%nv
       do sw = 2, stencil_width
          domain_full(iproc)%nv_halo(sw) = domain_full(iproc)%nv_halo(sw-1)+ domain_halo(sw, iproc)%nv
       end do
       domain_full(iproc)%nv = domain_full(iproc)%nv_halo(stencil_width)

       allocate(domain_full(iproc)%v_index(domain_full(iproc)%nv))

       cnt = 1
#ifdef USE_INER
       ! add iner
       do iv = 1, domain_iner(iproc)%nv
          domain_full(iproc)%v_index(cnt) = domain_iner(iproc)%v_index(iv)
          cnt = cnt + 1
       end do
       ! add bdry
       do sw = stencil_width, 1, -1
          do iv = 1, domain_bdry(sw, iproc)%nv
             domain_full(iproc)%v_index(cnt) = domain_bdry(sw, iproc)%v_index(iv)
             cnt = cnt + 1
          end do
       end do
#else
       do iv = 1, domain_compute(iproc)%nv
          domain_full(iproc)%v_index(cnt) = domain_compute(iproc)%v_index(iv)
          cnt = cnt + 1
       end do
#endif
       ! add halo
       do sw = 1, stencil_width
          do iv = 1, domain_halo(sw, iproc)%nv
             domain_full(iproc)%v_index(cnt) = domain_halo(sw, iproc)%v_index(iv)
             cnt = cnt + 1
          end do
       end do

       ! add t and e
       sw = 1
       sw2=stencil_width
       do iv = 1, domain_full(iproc)%nv
          iv_global = domain_full(iproc)%v_index(iv)

          do it_local = 1, mesh%vtx_nnb(iv_global)
             call full_t_list%insert(int(mesh%vtx_tr%f(iv_global, it_local)))
             call full_e_list%insert(int(mesh%vtx_ed%f(iv_global, it_local)))
          end do

          if(iv .eq. domain_full(iproc)%nv_bdry(sw2))then
            domain_full(iproc)%nt_bdry(sw2) = full_t_list%size()
            domain_full(iproc)%ne_bdry(sw2) = full_e_list%size()
            if(sw2 .gt. 1) sw2=sw2-1
          end if

          if(iv .eq. domain_full(iproc)%nv_compute) then
             domain_full(iproc)%nt_compute = full_t_list%size()
             domain_full(iproc)%ne_compute = full_e_list%size()          
          end if

          if(iv .eq. domain_full(iproc)%nv_halo(sw))then
             domain_full(iproc)%nt_halo(sw) = full_t_list%size()
             domain_full(iproc)%ne_halo(sw) = full_e_list%size()
             sw=sw+1
          end if
       end do

       call full_t_list%dump(domain_full(iproc)%t_index)
       call full_e_list%dump(domain_full(iproc)%e_index)
       domain_full(iproc)%nt = size(domain_full(iproc)%t_index)
       domain_full(iproc)%ne = size(domain_full(iproc)%e_index)

       domain_full(iproc)%nv_full = domain_full(iproc)%nv    
       domain_full(iproc)%nt_full = domain_full(iproc)%nt
       domain_full(iproc)%ne_full = domain_full(iproc)%ne

       call full_t_list%final
       call full_e_list%final
    end do

    call clock_begin("write bdry_halo")
    open(101, file="full-"//i2s(nprocs), form='unformatted', access="stream")
    do iproc = 1, nprocs
       write(101) domain_full(iproc)%nv
       write(101) domain_full(iproc)%nv_compute
       do sw = 1, stencil_width
          write(101) domain_full(iproc)%nv_bdry(sw)
       end do
       do sw = 1, stencil_width
          write(101) domain_full(iproc)%nv_halo(sw)
       end do
       write(101) domain_full(iproc)%nt
       write(101) domain_full(iproc)%nt_compute
       do sw = 1, stencil_width
          write(101) domain_full(iproc)%nt_bdry(sw)
       end do
       do sw = 1, stencil_width
          write(101) domain_full(iproc)%nt_halo(sw)
       end do
       write(101) domain_full(iproc)%ne
       write(101) domain_full(iproc)%ne_compute
       do sw = 1, stencil_width
          write(101) domain_full(iproc)%ne_bdry(sw)
       end do
       do sw = 1, stencil_width
          write(101) domain_full(iproc)%ne_halo(sw)
       end do
#ifdef USE_INER
       write(101) domain_full(iproc)%nv_iner
#endif
#ifdef USE_INER
       write(101) domain_full(iproc)%nt_iner
#endif
#ifdef USE_INER
       write(101) domain_full(iproc)%ne_iner
#endif
    end do

    do iproc = 1, nprocs
       write(101) domain_full(iproc)%v_index
       write(101) domain_full(iproc)%t_index
       write(101) domain_full(iproc)%e_index
    end do
    close(101)
    !! test
    !do iproc = 1, nprocs
    !   open(101, file="full-v-"//i2s(iproc)//".txt", form='formatted')
    !   write(101, *) domain_full(iproc)%nv
    !   write(101, *) domain_full(iproc)%nv_compute
    !   do sw = 1, stencil_width
    !      write(101, *) domain_full(iproc)%nv_bdry(sw)
    !   end do
    !   do sw = 1, stencil_width
    !      write(101, *) domain_full(iproc)%nv_halo(sw)
    !   end do
    !   do iv = 1, domain_full(iproc)%nv
    !      write(101, *) domain_full(iproc)%v_index(iv)
    !   end do
    !   close(101)
    !   open(101, file="full-e-"//i2s(iproc)//".txt", form='formatted')
    !   write(101, *) domain_full(iproc)%ne
    !   write(101, *) domain_full(iproc)%ne_compute
    !   do sw = 1, stencil_width
    !      write(101, *) domain_full(iproc)%ne_bdry(sw)
    !   end do
    !   do sw = 1, stencil_width
    !      write(101, *) domain_full(iproc)%ne_halo(sw)
    !   end do
    !   do iv = 1, domain_full(iproc)%ne
    !      write(101, *) domain_full(iproc)%e_index(iv)
    !   end do
    !   close(101)
    !   open(101, file="full-t-"//i2s(iproc)//".txt", form='formatted')
    !   write(101, *) domain_full(iproc)%nt
    !   write(101, *) domain_full(iproc)%nt_compute
    !   do sw = 1, stencil_width
    !      write(101, *) domain_full(iproc)%nt_bdry(sw)
    !   end do
    !   do sw = 1, stencil_width
    !      write(101, *) domain_full(iproc)%nt_halo(sw)
    !   end do
    !   do iv = 1, domain_full(iproc)%nt
    !      write(101, *) domain_full(iproc)%t_index(iv)
    !   end do
    !   close(101)
    !end do
    !open(101, file="parts.txt", form='formatted')
    !do iv = 1, mesh%nv
    !   write(101, *) mesh%parts(iv)
    !end do
    !close(101)
    call clock_end("write bdry_halo")

    do iproc = 1, nprocs
       call compute_v_list(iproc)%final
       do sw = 1, stencil_width
          deallocate(domain_bdry(sw, iproc)%v_index)
       end do
    end do

    deallocate(bdry_v_list, bdry_t_list, bdry_e_list)
    deallocate(halo_v_list, halo_t_list, halo_e_list)

    deallocate(compute_v_list)
    deallocate(domain_bdry)
    deallocate(domain_iner)
  end subroutine init_bdry_halo_index 

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
    character(len = 1), parameter :: BASE_32(0:31) = (/'0','1','2','3','4','5','6','7','8','9', &
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
    integer(i4), intent(in)                 :: square
    integer, intent(in)                     :: i, j
    integer(i4), intent(inout)              :: hilbert_pos
    integer                                 :: rx, ry, s, iv, jv

    iv = i
    jv = j
  ! do s = square/2, 1, s = s/2
  !     rx = (i & s) > 0
  !     ry = (j & s) > 0

  !     hilbert_pos = hilbert_pos + s * s * ((3 * rx) ^ ry)
  !     call rot(s, i, j, rx, ry)
  ! enddo
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
    enddo
  end subroutine xy2d

  subroutine rot(s, x, y, rx, ry)
    integer, intent(in)                     :: s, rx, ry
    integer, intent(inout)                  :: x, y
    integer                                 :: t

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

  subroutine init_send_recv_index(mesh, nprocs)
    implicit none
    ! io
    type(global_domain_data), intent(in) :: mesh
    integer,                  intent(in) :: nprocs

    ! local
    integer(i4)                          :: iv, iv_global
    integer                              :: c, cpuid, iproc, nnb, mark, marke, markt
    integer(i4)                          :: i, it_local, sw, nid
    type(set), allocatable               :: recv_v(:,:), recv_t(:,:), recv_e(:,:)

    type(set), allocatable               :: recv_cpus_list(:)
    type(index_list), allocatable        :: recv_cpus(:)

    type(index_list), allocatable        :: recv_sites(:,:), send_sites(:,:)

    integer                              :: ie, it, n, ne, nt, nv, kcpu, pos, max_cpu_nnb
    type(map), allocatable               :: cpumap(:)
    logical                              :: is_comm_nb

    max_cpu_nnb = 50
    allocate(recv_v(max_cpu_nnb, nprocs), recv_e(max_cpu_nnb, nprocs), recv_t(max_cpu_nnb, nprocs))
    allocate(recv_cpus_list(nprocs))
    allocate(recv_cpus(nprocs))
    allocate(recv_sites(max_cpu_nnb, nprocs), send_sites(max_cpu_nnb, nprocs))
    allocate(cpumap(nprocs))

    do iproc = 1, nprocs
       kcpu = 0
       call cpumap(iproc)%init
       do sw = 1, stencil_width
          do iv = 1, domain_halo(sw, iproc)%nv 

             iv_global = domain_halo(sw, iproc)%v_index(iv)
#ifdef LAM_DOMAIN
             ! discard cells/edges/triangles in the specific region
             if (mesh%parts(iv_global) .eq. iproc-1) cycle
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
                   if (marke/10 .ne. 177 .and. .not. bdry1_e_list(iproc)%find(ie) &
                       .or. markt .ne. 177 .and. .not. bdry1_t_list(iproc)%find(it)) then
                      is_comm_nb = .true.
                      exit
                   end if
                end do
             end if
             if (.not. is_comm_nb) cycle
#endif
             nid = mesh%parts(iv_global)

             pos = cpumap(iproc)%find(nid+1)
             if(pos .eq. -1) then
                 kcpu = kcpu + 1
                 call cpumap(iproc)%insert(nid+1, kcpu)
                 pos = kcpu
             end if

             call recv_cpus_list(iproc)%insert(nid)

#ifdef LAM_DOMAIN
             if (mark .ne. 177) call recv_v(pos, iproc)%insert(iv_global)
             do n = 1, mesh%vtx_nnb(iv_global)
                ie = mesh%vtx_ed%f(iv_global, n)
                it = mesh%vtx_tr%f(iv_global, n)
                marke = mesh%edt_mark%f(ie, 1)
                markt = mesh%tri_mark%f(it, 1)
                if (marke/10 .ne. 177 .and. .not. bdry1_e_list(iproc)%find(ie)) call recv_e(pos, iproc)%insert(ie)
                if (markt .ne. 177 .and. .not. bdry1_t_list(iproc)%find(it)) call recv_t(pos, iproc)%insert(it)
             end do
#else
             call recv_v(pos, iproc)%insert(iv_global)

             do n = 1, mesh%vtx_nnb(iv_global)
                ie = mesh%vtx_ed%f(iv_global, n)
                it = mesh%vtx_tr%f(iv_global, n)
                call recv_e(pos, iproc)%insert(ie)
                call recv_t(pos, iproc)%insert(it)
             end do
#endif
          end do
#ifdef LAM_DOMAIN
          if (sw .eq. stencil_width) then
             do i = 1, domain_halo(sw, iproc)%ne
                ie = domain_halo(sw, iproc)%e_index(i)
                marke = mesh%edt_mark%f(ie, 1)
                if (marke .eq. 1771) then
                   nid = mesh%parts(int(mesh%edt_v%f(ie, 1)))
                   if (nid .eq. iproc-1) nid = mesh%parts(int(mesh%edt_v%f(ie, 2)))
                   if (nid .ne. iproc-1 .and. .not. halo1_e_list(iproc)%find(ie)) then
                      pos = cpumap(iproc)%find(nid+1)
                      if(pos .eq. -1) then
                          kcpu = kcpu + 1
                          call cpumap(iproc)%insert(nid+1, kcpu)
                          pos = kcpu
                      end if
                      call recv_cpus_list(iproc)%insert(nid)
                      call recv_e(pos, iproc)%insert(ie)
                   end if
                end if
             end do
          end if
#endif
          deallocate(domain_halo(sw, iproc)%v_index)
#ifdef LAM_DOMAIN
          deallocate(domain_halo(sw, iproc)%e_index)
          deallocate(domain_halo(sw, iproc)%t_index)
#endif
       end do

       call recv_cpus_list(iproc)%dump(recv_cpus(iproc)%v_index)
       recv_cpus(iproc)%nv = size(recv_cpus(iproc)%v_index)

       do c = 1, recv_cpus(iproc)%nv
          call recv_v(c, iproc)%dump(recv_sites(c, iproc)%v_index)
          call recv_e(c, iproc)%dump(recv_sites(c, iproc)%e_index)
          call recv_t(c, iproc)%dump(recv_sites(c, iproc)%t_index)
       end do
#ifdef LAM_DOMAIN
       call bdry1_e_list(iproc)%final
       call bdry1_t_list(iproc)%final
       call halo1_e_list(iproc)%final
#endif
    end do
#ifdef LAM_DOMAIN
    deallocate(bdry1_e_list, bdry1_t_list, halo1_e_list)
#endif

    do iproc = 1, nprocs
       do c = 1, recv_cpus(iproc)%nv
          cpuid = recv_cpus(iproc)%v_index(c)
          pos = cpumap(cpuid+1)%find(iproc)

          call recv_v(pos, cpuid+1)%dump(send_sites(c, iproc)%v_index)
          call recv_e(pos, cpuid+1)%dump(send_sites(c, iproc)%e_index)
          call recv_t(pos, cpuid+1)%dump(send_sites(c, iproc)%t_index)
       end do
    end do

    call clock_begin("write comm")
    open(101, file="comm-"//i2s(nprocs), form='unformatted', access="stream")
    do iproc = 1, nprocs
        write(101) recv_cpus(iproc)%nv
    end do
    do iproc = 1, nprocs
        write(101) recv_cpus(iproc)%v_index
    end do
    do iproc = 1, nprocs
       do c = 1, recv_cpus(iproc)%nv
          write(101) size(send_sites(c, iproc)%v_index)
          write(101) size(recv_sites(c, iproc)%v_index)
          write(101) size(send_sites(c, iproc)%e_index)
          write(101) size(recv_sites(c, iproc)%e_index)
          write(101) size(send_sites(c, iproc)%t_index)
          write(101) size(recv_sites(c, iproc)%t_index)
       end do
    end do
    do iproc = 1, nprocs
       do c = 1, recv_cpus(iproc)%nv
          write(101) send_sites(c, iproc)%v_index
          write(101) recv_sites(c, iproc)%v_index
          write(101) send_sites(c, iproc)%e_index
          write(101) recv_sites(c, iproc)%e_index
          write(101) send_sites(c, iproc)%t_index
          write(101) recv_sites(c, iproc)%t_index
       end do
    end do
    close(101)
    !! test
    !do iproc = 1, nprocs
    !   do c = 1, recv_cpus(iproc)%nv
    !      open(101, file="send-v-"//i2s(iproc)//"-"//i2s(recv_cpus(iproc)%v_index(c))//".txt", form='formatted')
    !      do iv = 1, size(send_sites(c, iproc)%v_index)
    !         write(101, *) send_sites(c, iproc)%v_index(iv)
    !      end do
    !      close(101)
    !      open(101, file="recv-v-"//i2s(iproc)//"-"//i2s(recv_cpus(iproc)%v_index(c))//".txt", form='formatted')
    !      do iv = 1, size(recv_sites(c, iproc)%v_index)
    !         write(101, *) recv_sites(c, iproc)%v_index(iv)
    !      end do
    !      close(101)
    !   end do
    !end do
    call clock_end("write comm")

    do iproc = 1, nprocs
       do c = 1, max_cpu_nnb
          call recv_v(c, iproc)%final
          call recv_t(c, iproc)%final
          call recv_e(c, iproc)%final
          if(allocated(recv_sites(c, iproc)%v_index)) deallocate(recv_sites(c,iproc)%v_index)
          if(allocated(recv_sites(c, iproc)%t_index)) deallocate(recv_sites(c,iproc)%t_index)
          if(allocated(recv_sites(c, iproc)%e_index)) deallocate(recv_sites(c,iproc)%e_index)
          if(allocated(send_sites(c, iproc)%v_index)) deallocate(send_sites(c,iproc)%v_index)
          if(allocated(send_sites(c, iproc)%t_index)) deallocate(send_sites(c,iproc)%t_index)
          if(allocated(send_sites(c, iproc)%e_index)) deallocate(send_sites(c,iproc)%e_index)
       end do
       call recv_cpus_list(iproc)%final
       call cpumap(iproc)%final
       deallocate(recv_cpus(iproc)%v_index)
    end do

    deallocate(recv_v, recv_e, recv_t)
    deallocate(recv_cpus_list)
    deallocate(recv_cpus)
    deallocate(recv_sites, send_sites)
    deallocate(cpumap)
    return
  end subroutine init_send_recv_index

  subroutine init_output_index(mesh, nprocs)
    implicit none
    type(global_domain_data), intent(in) :: mesh
    integer,                  intent(in) :: nprocs

    integer                              :: iv, ie, it, k, sw, iv_global, iproc
    type(set)                            :: e_list, t_list
    logical                              :: insert_ed, insert_tr
    integer(i4)                          :: ed, tr, nv, ne, nt, v_tmp, act_iv
    integer                              :: belong_cpu
    integer,          allocatable        :: points(:)
    type(index_list), allocatable        :: domain_output(:)

    allocate(domain_output(nprocs))

    ! initialize output v_index
#ifdef LAM_DOMAIN
    do iproc = 1, nprocs
       allocate(points(domain_full(iproc)%nv))
       act_iv = 0
       do iv = 1, domain_full(iproc)%nv
          if (mesh%parts(domain_full(iproc)%v_index(iv)) .eq. iproc-1) then
             act_iv = act_iv + 1
             points(act_iv) = domain_full(iproc)%v_index(iv)
          end if
       end do
       allocate(domain_output(iproc)%v_index(act_iv))
       domain_output(iproc)%v_index = points(1:act_iv)
       domain_output(iproc)%nv = act_iv
       deallocate(points)
    end do
#else
    do iproc = 1, nprocs
       allocate(domain_output(iproc)%v_index(domain_compute(iproc)%nv))
       domain_output(iproc)%v_index = domain_compute(iproc)%v_index
       domain_output(iproc)%nv = domain_compute(iproc)%nv
    end do
#endif

    do iproc = 1, nprocs
       do iv = 1, domain_output(iproc)%nv
          iv_global = domain_output(iproc)%v_index(iv)

          do ie = 1, mesh%vtx_nnb(iv_global)
             ed = mesh%vtx_ed%f(iv_global, ie)
             tr = mesh%vtx_tr%f(iv_global, ie)
             insert_ed = .true.
             insert_tr = .true.

             !check trianglar point belong to which part
             do k = 1, 3
#ifdef LAM_DOMAIN
                v_tmp = mesh%tri_v%f(tr, k)
                if (v_tmp .gt. 0) then
#endif
                belong_cpu = mesh%parts(int(mesh%tri_v%f(tr, k)))
                if(belong_cpu > iproc-1) then
                   insert_tr = .false.
                end if
#ifdef LAM_DOMAIN
                end if
#endif
             end do
             !check edge point belong to which part
             do k = 1, 2
#ifdef LAM_DOMAIN
                v_tmp = mesh%edt_v%f(ed, k)
                if (v_tmp .gt. 0) then
#endif
                belong_cpu = mesh%parts(int(mesh%edt_v%f(ed, k)))
                if(belong_cpu > iproc-1) then
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

       call e_list%dump(domain_output(iproc)%e_index)
       call t_list%dump(domain_output(iproc)%t_index)

       call e_list%final
       call t_list%final
       domain_output(iproc)%nt = size(domain_output(iproc)%t_index)
       domain_output(iproc)%ne = size(domain_output(iproc)%e_index)
    end do

    call clock_begin("write output")
    open(101, file="output-"//i2s(nprocs), form='unformatted', access="stream")
    do iproc = 1, nprocs
       write(101) domain_output(iproc)%nv
       write(101) domain_output(iproc)%nt
       write(101) domain_output(iproc)%ne
    end do
    do iproc = 1, nprocs
       write(101) domain_output(iproc)%v_index
       write(101) domain_output(iproc)%t_index
       write(101) domain_output(iproc)%e_index
    end do
    close(101)
    call clock_end("write output")
    !! test
    !do iproc = 1, nprocs
    !   open(101, file="output-v-"//i2s(iproc)//".txt", form='formatted')
    !   write(101, *) domain_output(iproc)%nv
    !   do iv = 1, domain_output(iproc)%nv
    !      write(101, *) domain_output(iproc)%v_index(iv)
    !   end do
    !   close(101)
    !   open(101, file="output-e-"//i2s(iproc)//".txt", form='formatted')
    !   write(101, *) domain_output(iproc)%ne
    !   do iv = 1, domain_output(iproc)%ne
    !      write(101, *) domain_output(iproc)%e_index(iv)
    !   end do
    !   close(101)
    !   open(101, file="output-t-"//i2s(iproc)//".txt", form='formatted')
    !   write(101, *) domain_output(iproc)%nt
    !   do iv = 1, domain_output(iproc)%nt
    !      write(101, *) domain_output(iproc)%t_index(iv)
    !   end do
    !   close(101)
    !end do
    do iproc = 1, nprocs
       deallocate(domain_output(iproc)%v_index)
       deallocate(domain_output(iproc)%t_index)
       deallocate(domain_output(iproc)%e_index)
    end do
    deallocate(domain_output)
    do iproc = 1, nprocs
       deallocate(domain_compute(iproc)%v_index)
    end do
    deallocate(domain_compute)

    do iproc = 1, nprocs
       deallocate(domain_full(iproc)%nv_bdry)
       deallocate(domain_full(iproc)%nv_halo)
       deallocate(domain_full(iproc)%nt_bdry)
       deallocate(domain_full(iproc)%nt_halo)
       deallocate(domain_full(iproc)%ne_bdry)
       deallocate(domain_full(iproc)%ne_halo)
       deallocate(domain_full(iproc)%v_index)
       deallocate(domain_full(iproc)%t_index)
       deallocate(domain_full(iproc)%e_index)
    end do
    deallocate(domain_full)
  end subroutine init_output_index

end module partition_routine

PROGRAM GRIST_partition
  use grist_constants_dbl,             only: i4, r8, i8, pi
  use grist_data_types,                only: scalar_2d_field
  use grist_domain_types,              only: global_domain_data
  use grist_fileio_0d_module_gcm,      only: wrap_read_0d
  use grist_fileio_list_2d_module_par, only: wrap_read_2d
  use grist_lib
  use grist_nml_module,                only: index_flag
  use partition_routine
  use grist_clocks,                    only: clock_id, clock_begin, clock_end, clock_summary, clock0

    character(128)          :: grid_file_0d_name
    character(128)          :: grid_file_2d_name
    integer(i8)             :: nv_i8, nt_i8, ne_i8, three_i8, two_i8, one_i8, maxvnb_i8
    integer                 :: nprocs, ierr
    character*80            :: buff

    type(scalar_2d_field)   :: vtx_ltln_nnb_l
    type(global_domain_data):: mesh
    integer                 :: readclock, metisclock, initcomputeclock, initbhclock, initsrclock, initoutputclock

    call grist_comm_set(MPI_COMM_WORLD)
    comm = grist_comm()

    call mpi_init(ierr)

    clock0           = clock_id('Total')
    readclock        = clock_id('Read time')
    metisclock       = clock_id('METIS time')
    initcomputeclock = clock_id('Init compute time')
    initbhclock      = clock_id('Init bdry halo time')
    initsrclock      = clock_id('Init send recv time')
    initoutputclock  = clock_id('Init output time')

    call clock_begin(clock0)

    call getarg(1,buff)
    read(buff, *) nprocs
    ! read namelist
    call nml_set_global_vars()
    three_i8 = 3
    two_i8   = 2
    one_i8   = 1

    ! read 0d file
    call clock_begin(readclock)
    grid_file_0d_name  = trim(gridFileNameHead)//".0d.nc"
    call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_nv',     nv_i8)
    call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_nt',     nt_i8)
    call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_ne',     ne_i8)
    call wrap_read_0d(MPI_COMM_SELF, gridFilePath, grid_file_0d_name, 'mesh_maxvnb', maxvnb_i8)
    mesh%nv     = nv_i8
    mesh%nt     = nt_i8
    mesh%ne     = ne_i8
    mesh%maxvnb = maxvnb_i8

    allocate(mesh%tri_v%f(mesh%nt,3))
    allocate(mesh%edt_v%f(mesh%ne,2))
    allocate(vtx_ltln_nnb_l%f(mesh%nv,3))
    allocate(mesh%vtx_nb%f(mesh%nv,mesh%maxvnb))
    allocate(mesh%vtx_ed%f(mesh%nv,mesh%maxvnb))
    allocate(mesh%vtx_tr%f(mesh%nv,mesh%maxvnb))
#ifdef LAM_DOMAIN
    allocate(mesh%plg_mark%f(mesh%nv,1))
    allocate(mesh%edt_mark%f(mesh%ne,1))
    allocate(mesh%tri_mark%f(mesh%nt,1))
#endif
    ! read 2d file
    grid_file_2d_name  = trim(gridFileNameHead)//".2d.nc"
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'tri_v',  nt_i8, three_i8, mesh%tri_v)
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'edt_v',  ne_i8, two_i8,   mesh%edt_v)
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_ltln_nnb', nv_i8, three_i8, vtx_ltln_nnb_l)
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_nb', nv_i8, maxvnb_i8, mesh%vtx_nb)
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_ed', nv_i8, maxvnb_i8, mesh%vtx_ed)
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'vtx_tr', nv_i8, maxvnb_i8, mesh%vtx_tr)
#ifdef LAM_DOMAIN
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'plg_mark', nv_i8, one_i8, mesh%plg_mark)
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'edt_mark', nv_i8, one_i8, mesh%edt_mark)
    call wrap_read_2d(MPI_COMM_SELF, gridFilePath, grid_file_2d_name, 'tri_mark', nv_i8, one_i8, mesh%tri_mark)
#endif

    allocate(mesh%vtx_nnb(mesh%nv))
    mesh%vtx_nnb(:) = int(vtx_ltln_nnb_l%f(:,3))
    if(index_flag .eq. "morton") then
      allocate(mesh%vtx_ltln(mesh%nv,2))
      mesh%vtx_ltln(:,:) = vtx_ltln_nnb_l%f(:,1:2)
    end if
    deallocate(vtx_ltln_nnb_l%f)
    call clock_end(readclock)

    allocate(mesh%parts(mesh%nv))
    call clock_begin(metisclock)
    call decompose(mesh, nprocs)
    call clock_end(metisclock)
    call clock_begin(initcomputeclock)
    call init_compute_index(mesh, nprocs)
    call clock_end(initcomputeclock)
    call clock_begin(initbhclock)
    call init_bdry_halo_index(mesh, nprocs)
    call clock_end(initbhclock)
    call clock_begin(initsrclock)
    call init_send_recv_index(mesh, nprocs)
    call clock_end(initsrclock)
    call clock_begin(initoutputclock)
    call init_output_index(mesh, nprocs)
    call clock_end(initoutputclock)

    deallocate(domain_halo)

    deallocate(mesh%parts)
    deallocate(mesh%tri_v%f)
    deallocate(mesh%edt_v%f)
    deallocate(mesh%vtx_nb%f)
    deallocate(mesh%vtx_ed%f)
    deallocate(mesh%vtx_tr%f)
    deallocate(mesh%vtx_nnb)
    if(index_flag .eq. "morton") then
      deallocate(mesh%vtx_ltln)
    end if
#ifdef LAM_DOMAIN
    deallocate(mesh%plg_mark%f)
    deallocate(mesh%edt_mark%f)
    deallocate(mesh%tri_mark%f)
#endif
    call clock_end(clock0)
    call clock_summary

    call mpi_finalize(ierr)
end program

