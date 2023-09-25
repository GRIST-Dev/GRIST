module cmparray_mod

  use grist_constants,      only: r8
 
  implicit none
  private
  save
  
  public expdaynite, cmpdaynite

  interface CmpDayNite
    module procedure CmpDayNite_1d_R
    module procedure CmpDayNite_2d_R
    module procedure CmpDayNite_3d_R
    module procedure CmpDayNite_1d_R_Copy
    module procedure CmpDayNite_2d_R_Copy
    module procedure CmpDayNite_3d_R_Copy
    module procedure CmpDayNite_1d_I
    module procedure CmpDayNite_2d_I
    module procedure CmpDayNite_3d_I
  end interface ! CmpDayNite

  interface ExpDayNite
    module procedure ExpDayNite_1d_R
    module procedure ExpDayNite_2d_R
    module procedure ExpDayNite_3d_R
    module procedure ExpDayNite_1d_I
    module procedure ExpDayNite_2d_I
    module procedure ExpDayNite_3d_I
  end interface ! ExpDayNite

  interface cmparray
    module procedure cmparray_1d_R
    module procedure cmparray_2d_R
    module procedure cmparray_3d_R
  end interface ! cmparray

  interface chksum
    module procedure chksum_1d_R
    module procedure chksum_2d_R
    module procedure chksum_3d_R
    module procedure chksum_1d_I
    module procedure chksum_2d_I
    module procedure chksum_3d_I
  end interface ! chksum

  contains

  subroutine CmpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1) :: Array

    call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, 1, 1, il1, iu1)

    return
  end subroutine CmpDayNite_1d_R

  subroutine CmpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il2, iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il2:iu2,il1:iu1) :: Array

    call CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, il2, iu2, il1, iu1)

    return
  end subroutine CmpDayNite_2d_R

  subroutine CmpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il3, iu3, il2,iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il3:iu3,il2:iu2,il1:iu1) :: Array

    real(r8), dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nnite) = Array(k,j,IdxNite(1:Nnite))
        Array(k,j,il1:il1+Nday-1) = Array(k,j,IdxDay(1:Nday))
        Array(k,j,il1+Nday:il1+Nday+Nnite-1) = tmp(1:Nnite)

      end do
    end do

    return
  end subroutine CmpDayNite_3d_R

  subroutine CmpDayNite_1d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il1:iu1) :: InArray
    real(r8), intent(out), dimension(il1:iu1) :: OutArray

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, 1, 1, 1, 1, il1, iu1)

    return
  end subroutine CmpDayNite_1d_R_Copy

  subroutine CmpDayNite_2d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il2, iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il2:iu2,il1:iu1) :: InArray
    real(r8), intent(out), dimension(il2:iu2,il1:iu1) :: OutArray

    call CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, 1, 1, il2, iu2, il1, iu1)

    return
  end subroutine CmpDayNite_2d_R_Copy

  subroutine CmpDayNite_3d_R_Copy(InArray, OutArray, Nday, IdxDay, Nnite, IdxNite, il3, iu3, il2,iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(in), dimension(il3:iu3,il2:iu2,il1:iu1) :: InArray
    real(r8), intent(out), dimension(il3:iu3,il2:iu2,il1:iu1) :: OutArray

    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

         do i=il1,il1+Nday-1
            OutArray(k,j,i) = InArray(k,j,IdxDay(i-il1+1))
         enddo
         do i=il1+Nday,il1+Nday+Nnite-1
            OutArray(k,j,i) = InArray(k,j,IdxNite(i-(il1+Nday)+1))
         enddo
        

      end do
    end do

    return
  end subroutine CmpDayNite_3d_R_Copy

  subroutine CmpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1) :: Array

    call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, 1, 1, il1, iu1)

    return
  end subroutine CmpDayNite_1d_I

  subroutine CmpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il2, iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il2:iu2,il1:iu1) :: Array

    call CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, il2, iu2, il1, iu1)

    return
  end subroutine CmpDayNite_2d_I

  subroutine CmpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il3, iu3, il2,iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il3:iu3,il2:iu2,il1:iu1) :: Array

    integer, dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nnite) = Array(k,j,IdxNite(1:Nnite))
        Array(k,j,il1:il1+Nday-1) = Array(k,j,IdxDay(1:Nday))
        Array(k,j,il1+Nday:il1+Nday+Nnite-1) = tmp(1:Nnite)

      end do
    end do

    return
  end subroutine CmpDayNite_3d_I

  subroutine ExpDayNite_1d_R(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il1:iu1) :: Array

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, 1, 1, il1, iu1)

    return
  end subroutine ExpDayNite_1d_R

  subroutine ExpDayNite_2d_R(Array, Nday, IdxDay, Nnite, IdxNite, il2, iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il2:iu2,il1:iu1) :: Array

    call ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, il2, iu2, il1, iu1)

    return
  end subroutine ExpDayNite_2d_R

  subroutine ExpDayNite_3d_R(Array, Nday, IdxDay, Nnite, IdxNite, il3, iu3, il2,iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    real(r8), intent(inout), dimension(il3:iu3,il2:iu2,il1:iu1) :: Array

    real(r8), dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nday) = Array(k,j,1:Nday)
        Array(k,j,IdxNite(1:Nnite)) = Array(k,j,il1+Nday:il1+Nday+Nnite-1)
        Array(k,j,IdxDay(1:Nday)) = tmp(1:Nday)

      end do
    end do

    return
  end subroutine ExpDayNite_3d_R

  subroutine ExpDayNite_1d_I(Array, Nday, IdxDay, Nnite, IdxNite, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il1:iu1) :: Array

    call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, 1, 1, il1, iu1)

    return
  end subroutine ExpDayNite_1d_I

  subroutine ExpDayNite_2d_I(Array, Nday, IdxDay, Nnite, IdxNite, il2, iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il2:iu2,il1:iu1) :: Array

    call ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, 1, 1, il2, iu2, il1, iu1)

    return
  end subroutine ExpDayNite_2d_I

  subroutine ExpDayNite_3d_I(Array, Nday, IdxDay, Nnite, IdxNite, il3, iu3, il2,iu2, il1, iu1)
    integer, intent(in) :: Nday, Nnite
    integer, intent(in) :: il1, iu1
    integer, intent(in) :: il2, iu2
    integer, intent(in) :: il3, iu3
    integer, intent(in), dimension(Nday) :: IdxDay
    integer, intent(in), dimension(Nnite) :: IdxNite
    integer, intent(inout), dimension(il3:iu3,il2:iu2,il1:iu1) :: Array

    integer, dimension(il1:iu1) :: tmp
    integer :: i, j, k


    do k = il3, iu3
      do j = il2, iu2

        tmp(1:Nday) = Array(k,j,1:Nday)
        Array(k,j,IdxNite(1:Nnite)) = Array(k,j,il1+Nday:il1+Nday+Nnite-1)
        Array(k,j,IdxDay(1:Nday)) = tmp(1:Nday)

      end do
    end do

    return
  end subroutine ExpDayNite_3d_I

!******************************************************************************!
!                                                                              !
!                                 DEBUG                                        !
!                                                                              !
!******************************************************************************!

  subroutine cmparray_1d_R(name, Ref, New, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    real(r8), intent(in), dimension(id1) :: Ref
    real(r8), intent(in), dimension(id1) :: New

    call cmparray_3d_R(name, Ref, New, 1, 1, 1, 1, 1, 1, id1, is1, ie1)
  end subroutine cmparray_1d_R

  subroutine cmparray_2d_R(name, Ref, New, id2, is2, ie2, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    real(r8), intent(in), dimension(id2, id1) :: Ref
    real(r8), intent(in), dimension(id2, id1) :: New

    call cmparray_3d_R(name, Ref, New, 1, 1, 1, id2, is2, ie2, id1, is1, ie1)
  end subroutine cmparray_2d_R

  subroutine cmparray_3d_R(name, Ref, New, id3, is3, ie3, id2, is2, ie2, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
    real(r8), intent(in), dimension(id3, id2, id1) :: Ref
    real(r8), intent(in), dimension(id3, id2, id1) :: New

    integer :: i, j, k
    integer :: nerr
    logical :: found
    real(r8):: rdiff
    real(r8), parameter :: rtol = 1.0e-13_r8

    nerr = 0

    do k = is3, ie3
      do j = is2, ie2

        found = .false.
        do i = is1, ie1
          rdiff = abs(New(k,j,i)-Ref(k,j,i))
          rdiff = rdiff / merge(abs(Ref(k,j,i)), 1.0_r8, Ref(k,j,i) /= 0.0_r8)
          if ( rdiff > rtol ) then
            found = .true.
            exit
          end if
        end do

        if ( found ) then
          do i = is1, ie1
            rdiff = abs(New(k,j,i)-Ref(k,j,i))
            rdiff = rdiff / merge(abs(Ref(k,j,i)), 1.0_r8, Ref(k,j,i) /= 0.0_r8)
            if ( rdiff > rtol ) then
              print 666, name, k, j, i, Ref(k,j,i), New(k,j,i), rdiff
              nerr = nerr + 1
              if ( nerr > 10 ) stop
            end if
          end do
        end if

      end do
    end do

    return
666 format('cmp3d: ', a10, 3(1x, i4), 3(1x, e20.14))

  end subroutine cmparray_3d_R

  subroutine chksum_1d_R(name, Ref, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    real(r8), intent(in), dimension(id1) :: Ref

    call chksum_3d_R(name, Ref, 1, 1, 1, 1, 1, 1, id1, is1, ie1)
  end subroutine chksum_1d_R

  subroutine chksum_1d_I(name, Ref, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in), dimension(id1) :: Ref

    call chksum_3d_I(name, Ref, 1, 1, 1, 1, 1, 1, id1, is1, ie1)
  end subroutine chksum_1d_I

  subroutine chksum_2d_R(name, Ref, id2, is2, ie2, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    real(r8), intent(in), dimension(id2, id1) :: Ref

    call chksum_3d_R(name, Ref, 1, 1, 1, id2, is2, ie2, id1, is1, ie1)
  end subroutine chksum_2d_R

  subroutine chksum_2d_I(name, Ref, id2, is2, ie2, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in), dimension(id2, id1) :: Ref

    call chksum_3d_I(name, Ref, 1, 1, 1, id2, is2, ie2, id1, is1, ie1)
  end subroutine chksum_2d_I

  subroutine chksum_3d_R(name, Ref, id3, is3, ie3, id2, is2, ie2, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
!orig    real(r8), intent(in), dimension(id1, id2, id3) :: Ref
    real(r8), intent(in), dimension(is3:ie3, is2:ie2, is1:ie1) :: Ref

    real(r8) :: chksum
    real(r8) :: rmin, rmax
    integer :: i, j, k
    integer :: imin, jmin, kmin
    integer :: imax, jmax, kmax

    imin = is1 ; jmin = is2 ; kmin = is3
    imax = is1 ; jmax = is2 ; kmax = is3
    rmin = Ref(is3, is2, is1) ; rmax = rmin

    chksum = 0.0_r8

    do k = is3, ie3
      do j = is2, ie2
        do i = is1, ie1
          chksum = chksum + abs(Ref(k,j,i))
          if ( Ref(k,j,i) < rmin ) then
            rmin = Ref(k,j,i)
            imin = i ; jmin = j ; kmin = k
          end if
          if ( Ref(k,j,i) > rmax ) then
            rmax = Ref(k,j,i)
            imax = i ; jmax = j ; kmax = k
          end if
        end do
      end do
    end do

    print 666, name, chksum, imin, jmin, kmin, imax, jmax, kmax
666 format('chksum: ', a8, 1x, e20.14, 6(1x, i4))

  end subroutine chksum_3d_R

  subroutine chksum_3d_I(name, Ref, id3, is3, ie3, id2, is2, ie2, id1, is1, ie1)
    character(*), intent(in) :: name
    integer,  intent(in) :: id1, is1, ie1
    integer,  intent(in) :: id2, is2, ie2
    integer,  intent(in) :: id3, is3, ie3
    integer,  intent(in), dimension(id3, id2, id1) :: Ref

    integer :: i, j, k
    integer :: chksum
    chksum = 0

    do k = is3, ie3
      do j = is2, ie2
        do i = is1, ie1
          chksum = chksum + abs(Ref(k,j,i))
        end do
      end do
    end do

    print 666, name, chksum
666 format('chksum: ', a8, 1x, i8)

  end subroutine chksum_3d_I

end module cmparray_mod
