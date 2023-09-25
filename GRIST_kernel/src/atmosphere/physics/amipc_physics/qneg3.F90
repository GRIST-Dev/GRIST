!===================================================================================================
!
!  Created by LiXiaohan on 19/07/10, adopted from CAM5
!
! Purpose: 
!  Check moisture and tracers for minimum value, reset any below
!  minimum value to minimum value and return information to allow
!  warning message to be printed. The global average is NOT preserved.
!  
!  Author: J. Rosinski
!===================================================================================================

subroutine qneg3 (subnam  ,ncol  ,lver  ,lconst_beg  ,lconst_end  ,qmin  ,q )

   use grist_constants,                    only: r8, i4 
   use grist_mpi
   implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
   character*(*), intent(in) :: subnam ! name of calling routine

   integer, intent(in) :: ncol         ! number of atmospheric columns
   integer, intent(in) :: lver         ! number of vertical levels in column
   integer, intent(in) :: lconst_beg   ! beginning constituent
   integer, intent(in) :: lconst_end   ! ending    constituent

   real(r8), intent(in) :: qmin(lconst_beg:lconst_end)      ! Global minimum constituent concentration

!
! Input/Output arguments
!
   real(r8), intent(inout) :: q(lconst_beg:lconst_end,lver,ncol) ! moisture/tracer field
!
!---------------------------Local workspace-----------------------------
!
   integer indx(lver,ncol)  ! array of indices of points < qmin
   integer nval(lver)       ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer nn
   integer iwtmp
   integer i,ii,k           ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator

   logical found            ! true => at least 1 minimum violator found

   real(r8) worst           ! biggest violator

!
!-----------------------------------------------------------------------
!

   do m=lconst_beg,lconst_end
      nvals = 0
      found = .false.
      worst = 1.e35_r8
      iw = -1
!
! Test all field values for being less than minimum value. Set q = qmin
! for all such points. Trace offenders and identify worst one.
!
      do k=1,lver
         nval(k) = 0
         nn = 0
         do i=1,ncol
            if (q(m,k,i) < qmin(m)) then
               nn = nn + 1
               indx(k,nn) = i
            end if
         end do
         nval(k) = nn
      end do

      do k=1,lver
         if (nval(k) > 0) then
            found = .true.
            nvals = nvals + nval(k)
            iwtmp = -1
!cdir nodep,altcode=loopcnt
            do ii=1,nval(k)
               i = indx(k,ii)
               if (q(m,k,i) < worst) then
                  worst = q(m,k,i)
                  iwtmp = ii
               end if
            end do
            if (iwtmp /= -1 ) kw = k
            if (iwtmp /= -1 ) iw = indx(k,iwtmp)
!cdir nodep,altcode=loopcnt
            do ii=1,nval(k)
               i = indx(k,ii)
               q(m,k,i) = qmin(m)
            end do
         end if
      end do
      if (found .and. abs(worst)>1.e-12_r8) then
     !    write(6,9000)subnam,m,nvals,qmin(m),worst,iw,kw
     !    print*,'rank=',mpi_rank()
      end if
   end do
!

   return
9000 format(' QNEG3 from ',a,':m=',i3, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)

end subroutine qneg3




