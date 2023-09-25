
#include "config.h"

!> this module perform domain decomposition using metis
module grist_metis
  use grist_utils

contains
  ! ! !> this subroutine decomposes a domain
  ! subroutine metis_decomp(vtx_nb, nnb, nparts, parts)
  !   use iso_c_binding
  !   implicit none

  !   interface
  !      subroutine metis_partgraphkway(nvtx, ncon, xadj, adjncy, &
  !           vwght, vsize, adjwgt, &
  !           nparts, tpwgts, ubvec, &
  !           options, edgecut, part) bind(c)
  !        use iso_c_binding 
  !        integer :: nvtx, ncon, xadj(*), adjncy(*), &
  !             nparts, options(*), edgecut, part(*)
  !        type(c_ptr), value :: vwght, vsize, adjwgt, tpwgts, ubvec
  !      endsubroutine
  !   end interface
  !   real(8) :: start, finish

  !   !> describe the vertex-neigbhour relationship.
  !   !> shape(vtx_nb) = (max_num_neighbours, num_vertex).
  !   idx_t, dimension(:,:), intent(in)  :: vtx_nb

  !   !> 1D vector that describes the number of neighbours
  !   !> for each vertex. shape(nnb) = (num_vertex)
  !   idx_t, dimension(:),   intent(in)  :: nnb

  !   !> output parameter. subdomain flags (integer) for each
  !   !> vertex. shape(parts) = (num_vertex)
  !   idx_t, dimension(:),   intent(out) :: parts

  !   !> indicates how many sub domains should be decomposed.
  !   idx_t, intent(in) :: nparts

  !   integer :: vtx_nb_dim(2), max_nnb, vtx_num, i, j, cnt
  !   integer, allocatable :: xadj(:), adjncy(:)
  !   integer :: options(40), ncon, ierr, edgecut

  !   vtx_nb_dim = shape(vtx_nb)
  !   max_nnb = vtx_nb_dim(1)
  !   vtx_num = vtx_nb_dim(2)

  !   call assert(vtx_num > 2, &
  !        "number of vertices must be larger than 2.")

  !   call assert(nparts <= vtx_num, &
  !        "number of parts must be less than vtx_num.")

  !   allocate(xadj(vtx_num + 1), adjncy(size(vtx_nb)))

  !   cnt = 0
  !   do j = 1, vtx_num
  !      xadj(j) = cnt 
  !      do i = 1, max_nnb
  !         if(vtx_nb(i, j) > 0) then
  !            adjncy(cnt+1) = vtx_nb(i, j)
  !            cnt = cnt + 1
  !         end if
  !      end do
  !   end do
  !   xadj(vtx_num+1) = cnt - 1

  !   call METIS_SetDefaultOptions(options)
  !   options(18) = 1
  !   ncon = 1

  !   ! call cpu_time(start)

  !   call metis_partgraphkway(vtx_num, ncon, xadj, adjncy, &
  !        c_null_ptr, c_null_ptr, c_null_ptr,  &
  !        nparts, c_null_ptr, c_null_ptr, options, edgecut, parts)

  !   ! call cpu_time(finish)

  !   ! print*,"time elapsed for domain decomposition is ", &
  !   !      finish - start, " secs."

  !   parts = parts - 1

  !   deallocate(xadj, adjncy)
  ! end subroutine

  !> this subroutine decomposes a domain
  !> adjncy stores the nearest neighbors
  !> Note: indices in adjncy and xadj starts from 1.
  subroutine metis_decomp(num_vtx, adjncy, xadj, nparts, parts)
    use iso_c_binding
    implicit none

    interface
       subroutine metis_partgraphkway(nvtx, ncon, xadj, adjncy, &
            vwght, vsize, adjwgt, &
            nparts, tpwgts, ubvec, &
            options, edgecut, part) bind(c)
         use iso_c_binding 
         integer :: nvtx, ncon, xadj(*), adjncy(*), &
              nparts, options(*), edgecut, part(*)
         type(c_ptr), value :: vwght, vsize, adjwgt, tpwgts, ubvec
       endsubroutine
    end interface
    real(8) :: start, finish

    idx_t, intent(in) :: num_vtx

    !> describe the vertex-neigbhour relationship.
    !> shape(vtx_nb) = (max_num_neighbours, num_vertex).
    idx_t, dimension(:), intent(in)  :: adjncy

    !> 1D vector that describes the number of neighbours
    !> for each vertex. shape(nnb) = (num_vertex)
    idx_t, dimension(:),   intent(in)  :: xadj

    !> output parameter. subdomain flags (integer) for each
    !> vertex. shape(parts) = (num_vertex)
    idx_t, dimension(:),   intent(out), allocatable :: parts

    !> indicates how many sub domains should be decomposed.
    idx_t, intent(in) :: nparts

    integer :: vtx_nb_dim(2), max_nnb, i, j, cnt
    integer :: options(40), ncon, ierr, edgecut

    call assert(num_vtx > 2, &
         "number of vertices must be larger than 2.")

    call assert(nparts <= num_vtx, &
         "number of parts must be less than num_vtx.")

    if(.not. allocated(parts)) allocate(parts(num_vtx))
    
    call METIS_SetDefaultOptions(options)
    ncon = 1

    ! call cpu_time(start)
    ! print*, "num_vtx = ", num_vtx
    ! print*, "xadj = ", xadj
    ! print*, "adjncy = ", adjncy

    call metis_partgraphkway(num_vtx, ncon, xadj, adjncy, &
         c_null_ptr, c_null_ptr, c_null_ptr,  &
         nparts, c_null_ptr, c_null_ptr, options, edgecut, parts)

    ! call cpu_time(finish)
    ! print*,"time elapsed for domain decomposition is ", &
    !      finish - start, " secs."

    !parts = parts

  end subroutine

end module
