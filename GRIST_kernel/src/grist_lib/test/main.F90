
#include "../src/config.h"

module test
  use grist_lib

contains
  subroutine test_mpi()
    implicit none
    integer :: ierr
    call mpi_init(ierr)
    
    print*, "my rank = ", mpi_rank(), &
         " size = ", mpi_size()

    call mpi_finalize(ierr)
    
  end subroutine

  subroutine test_map()
    implicit none
    type(map) :: s
    idx_t, allocatable :: v(:)

    call s%init()
    call s%insert(3, 5)
    call s%insert(2, 4)
    call s%insert(3, 7)
    call s%insert(4, 6)
    
    print*, s%find(3)
    print*, s%find(2)
    print*, s%find(4)
  end subroutine

  
  subroutine test_set()
    implicit none
    type(set) :: s
    idx_t, allocatable :: v(:)

    call s%init()
    call s%insert(3)
    call s%insert(2)
    call s%insert(3)
    call s%insert(4)
    
    call s%dump(v)
    print*, "v = ", v

    print*, s%find(2)
    print*, s%find(5)
    print*, s%find(4)
  end subroutine

  subroutine test_metis()
    implicit none
    integer :: vtx_nb(0:100), xadj(0:100), i, j, NI, NJ, o(3), xi,xj,ii,jj
    integer :: cnt, cnt_xadj
    integer, allocatable :: parts(:)

    o = [-1, 0, 1]

    NJ = 4
    NI = 3
    cnt = 0

    xadj(0) = 0
    do j = 0, NJ - 1 
       do i = 0, NI - 1
          cnt_xadj = 0

          if(i - 1 >= 0) then
             vtx_nb(cnt) = (i-1) + j * NI
             cnt = cnt + 1
             cnt_xadj = cnt_xadj + 1           
          end if

          if(j - 1 >= 0) then
             vtx_nb(cnt) = i + (j-1) * NI
             cnt = cnt + 1
             cnt_xadj = cnt_xadj + 1           
          end if

          if(i + 1 < NI) then
             vtx_nb(cnt) = (i+1) + j * NI
             cnt = cnt + 1
             cnt_xadj = cnt_xadj + 1           
          end if

          if(j + 1 < NJ) then
             vtx_nb(cnt) = i + (j+1) * NI
             cnt = cnt + 1
             cnt_xadj = cnt_xadj + 1           
          end if

          xadj(i + j * NI + 1) = xadj(i + j * NI) + cnt_xadj
       end do
    end do

    ! print*, cnt
    ! print*, vtx_nb(0:cnt)
    ! print*, xadj(0:NI*NJ)
    parts = 0
    print*, "vtx_nb = ", vtx_nb

    call metis_decomp(NI*NJ, vtx_nb+1, xadj+1, 2, parts)

    print*, parts

  end subroutine
end module test

program main
  use grist_lib
  use test
  implicit none
  type(set) :: ss

  ! call test_set()

  ! call test_mpi()

  ! call test_metis()

  call test_map()
end program
