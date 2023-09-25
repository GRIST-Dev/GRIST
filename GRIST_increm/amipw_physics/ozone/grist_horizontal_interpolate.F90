module grist_horizontal_interpolate

    use grist_constants,                    only: r8, rad2deg
    use grist_nml_module,                   only: working_mode
    use grist_handle_error,                 only: endrun
    use grist_mpi

  implicit none
  private
  save 

  public :: xy_interp_init, xy_interp,  &
            interp_type, lininterp_init, lininterp, lininterp_finish


  type interp_type
     real(r8), pointer :: wgts(:)
     real(r8), pointer :: wgtn(:)
     integer, pointer  :: jjm(:)
     integer, pointer  :: jjp(:)
  end type interp_type
  interface lininterp
     module procedure lininterp_original
     module procedure lininterp_full1d
     module procedure lininterp1d
     module procedure lininterp2d1d
     module procedure lininterp3d2d
  end interface
 
contains
 
  subroutine xy_interp_init(ncol,grid_lon,grid_lat,im1,jm1,lon0,lat0,weight_x,weight_y,index_x,index_y)  
!------------------------------------------------------------------------------------------------------------
! This program computes weighting functions to map a variable of (im1,jm1) resolution to GRIST resolution
! weight_x(ncol) is the weighting function for zonal interpolation
! weight_y(ncol) is the weighting function for meridional interpolation
! 
! linear interpolation, LiXH
!------------------------------------------------------------------------------------------------------------
  !io
  integer,  intent(in)  :: ncol
  integer,  intent(in)  :: im1, jm1
  real(r8), intent(in)  :: lon0(im1), lat0(jm1)
  real(r8), intent(in)  :: grid_lon(ncol), grid_lat(ncol)
  real(r8), intent(out) :: weight_x(ncol,2), weight_y(ncol,2)
  integer,  intent(out) :: index_x(ncol,2), index_y(ncol,2)
  ! local  
  integer  :: icol, i, j
  logical  :: find_lon, find_Lat
  real(r8) :: lon1(im1), lat1(jm1)
  real(r8) :: lon2(ncol), lat2(ncol)


  weight_x = 0.0_r8
  weight_y = 0.0_r8

  lon1 = lon0*rad2deg
  lat1 = lat0*rad2deg

  lon2 = grid_lon*rad2deg
  lat2 = grid_lat*rad2deg


  do icol = 1, ncol
    ! lon
    if(lon2(icol) .lt. 0.)lon2(icol) = lon2(icol) + 360.
    find_lon = .false.
    do i = 1, im1-1
        if(lon1(i) .lt. lon2(icol) .and. lon1(i+1) .ge. lon2(icol))then
            find_lon = .true.
            index_x(icol,1) = i
            index_x(icol,2) = i+1
            weight_x(icol,1) = (lon2(icol)-lon1(i))/(lon1(i+1)-lon1(i))
            weight_x(icol,2) = 1.-weight_x(icol,1)
            exit
        end if
    end do

    !consider boundary:
    if(lon1(1) .ge. lon2(icol))then
        find_lon = .true.
        index_x(icol,1) = im1
        index_x(icol,2) = 1
        weight_x(icol,1) = (lon2(icol)-lon1(im1)+360.)/(lon1(1)-lon1(im1)+360.)
        weight_x(icol,2) = 1.-weight_x(icol,1)
    else if(lon1(im1) .lt. lon2(icol))then
        find_lon = .true.
        index_x(icol,1) = im1
        index_x(icol,2) = 1
        weight_x(icol,1) = (lon2(icol)-lon1(im1))/(lon1(1)-lon1(im1)+360.)
        weight_x(icol,2) = 1.-weight_x(icol,1)
    end if

    ! lat
    find_lat = .false.
    do j = 1, jm1
        if(lat1(j) .lt. lat2(icol) .and. lat1(j+1) .ge. lat2(icol))then
            find_lat = .true.
            index_y(icol,1) = j
            index_y(icol,2) = j+1
            weight_y(icol,1) = (lat2(icol)-lat1(j))/(lat1(j+1)-lon1(j))
            weight_y(icol,2) = 1.-weight_y(icol,1)
            exit
        end if
    end do

    !consider boundary:
    if(lat1(jm1) .lt. lat2(icol) )then
        find_lat = .true.
        index_y(icol,1) = jm1
        index_y(icol,2) = jm1
        weight_y(icol,1) = 1.
        weight_y(icol,2) = 0.
    else if(lat1(1) .ge. lat2(icol))then
        find_lat = .true.
        index_y(icol,1) = 1
        index_y(icol,2) = 1
        weight_y(icol,1) = 1.
        weight_y(icol,2) = 0.
    end if

  end do


  if( .not. find_lat .or. .not. find_lon)then
      print*,'Rank ',mpi_rank(),'can not complete hori interpolation!!'
      call endrun('check hori interpolation in chemistry_aero.')
  end if

  end subroutine xy_interp_init


  subroutine xy_interp(im1,jm1,km1,ncol,weight_x,weight_y,var_src,var_trg,index_x,index_y)
  !io
  integer,  intent(in)  :: im1   ! source number of longitudes
  integer,  intent(in)  :: jm1   ! source number of latitudes
  integer,  intent(in)  :: km1   ! source/target number of levels
  integer,  intent(in)  :: ncol
  real(r8), intent(in)  :: weight_x(ncol,2), weight_y(ncol,2)
  integer,  intent(in)  :: index_x(ncol,2), index_y(ncol,2)
  real(r8), intent(in)  :: var_src(im1,jm1,km1)
  real(r8), intent(out) :: var_trg(km1,ncol)
  ! local
  integer  :: k1, n
  real(r8) :: var_tmp1(km1), var_tmp2(km1)

  var_trg = 0.0_r8

  do n=1,ncol
    var_tmp1(:)  = var_src(index_x(n,1),index_y(n,1),:)*weight_x(n,1)+var_src(index_x(n,2),index_y(n,1),:)*weight_x(n,2)
    var_tmp2(:)  = var_src(index_x(n,1),index_y(n,2),:)*weight_x(n,1)+var_src(index_x(n,2),index_y(n,2),:)*weight_x(n,2)
    var_trg(:,n) = var_tmp1(:)*weight_y(n,1)+var_tmp2(:)*weight_y(n,2)
  end do

  end subroutine xy_interp



  subroutine lininterp_full1d (arrin, yin, nin, arrout, yout, nout)
    integer, intent(in) :: nin, nout
    real(r8), intent(in) :: arrin(nin), yin(nin), yout(nout)
    real(r8), intent(out) :: arrout(nout)
    type (interp_type) :: interp_wgts

    call lininterp_init(yin, nin, yout, nout, 1, interp_wgts)
    call lininterp1d(arrin, nin, arrout, nout, interp_wgts)
    call lininterp_finish(interp_wgts)

  end subroutine lininterp_full1d


  subroutine lininterp_init(yin, nin, yout, nout, extrap_method, interp_wgts, &
       cyclicmin, cyclicmax)
!
! Description:
!   Initialize a variable of type(interp_type) with weights for linear interpolation.
!       this variable can then be used in calls to lininterp1d and lininterp2d.
!   yin is a 1d array of length nin of locations to interpolate from - this array must 
!       be monotonic but can be increasing or decreasing
!   yout is a 1d array of length nout of locations to interpolate to, this array need
!       not be ordered
!   extrap_method determines how to handle yout points beyond the bounds of yin
!       if 0 set values outside output grid to 0 
!       if 1 set to boundary value
!       if 2 set to cyclic boundaries
!         optional values cyclicmin and cyclicmax can be used to set the bounds of the 
!         cyclic mapping - these default to 0 and 360.
!

    integer, intent(in) :: nin
    integer, intent(in) :: nout
    real(r8), intent(in) :: yin(:)           ! input mesh
    real(r8), intent(in) :: yout(:)         ! output mesh
    integer, intent(in) :: extrap_method       ! if 0 set values outside output grid to 0 
                                               ! if 1 set to boundary value
                                               ! if 2 set to cyclic boundaries
    real(r8), intent(in), optional :: cyclicmin, cyclicmax

    type (interp_type), intent(out) :: interp_wgts
    real(r8) :: cmin, cmax
    real(r8) :: extrap
    real(r8) :: dyinwrap
    real(r8) :: ratio
    real(r8) :: avgdyin
    integer :: i, j, icount
    integer :: jj
    real(r8), pointer :: wgts(:)
    real(r8), pointer :: wgtn(:)
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)
    logical :: increasing
    !
    ! Check validity of input coordinate arrays: must be monotonically increasing,
    ! and have a total of at least 2 elements
    !
    if (nin.lt.2) then
       call endrun('LININTERP: Must have at least 2 input points for interpolation')
    end if
    if(present(cyclicmin)) then
       cmin=cyclicmin
    else
       cmin=0._r8
    end if
    if(present(cyclicmax)) then
       cmax=cyclicmax
    else
       cmax=360._r8
    end if
    if(cmax<=cmin) then
       call endrun('LININTERP: cyclic min value must be < max value')
    end if
    increasing=.true.
    icount = 0
    do j=1,nin-1
       if (yin(j).gt.yin(j+1)) icount = icount + 1
    end do
    if(icount.eq.nin-1) then
       increasing = .false.
       icount=0
    endif
    if (icount.gt.0) then
       call endrun('LININTERP: Non-monotonic input coordinate array found')
    end if
    allocate(interp_wgts%jjm(nout), &
         interp_wgts%jjp(nout), &
         interp_wgts%wgts(nout), &
         interp_wgts%wgtn(nout))

    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn

    !
    ! Initialize index arrays for later checking
    !
    jjm = 0
    jjp = 0

    extrap = 0._r8
    if(extrap_method.eq.0) then
       !
       ! For values which extend beyond boundaries, set weights
       ! such that values will be 0.
       !
       do j=1,nout
          if(increasing) then
             if (yout(j).lt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).lt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 0._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          end if
       end do
    else if(extrap_method.eq.1) then
       !
       ! For values which extend beyond boundaries, set weights
       ! such that values will just be copied.
       !
       do j=1,nout
          if(increasing) then
             if (yout(j).le.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).gt.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          else
             if (yout(j).gt.yin(1)) then
                jjm(j) = 1
                jjp(j) = 1
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             else if (yout(j).le.yin(nin)) then
                jjm(j) = nin
                jjp(j) = nin
                wgts(j) = 1._r8
                wgtn(j) = 0._r8
                extrap = extrap + 1._r8
             end if
          end if
       end do
    else if(extrap_method.eq.2) then
       !
       ! For values which extend beyond boundaries, set weights
       ! for circular boundaries 
       !
       dyinwrap = yin(1) + (cmax-cmin) - yin(nin)
       avgdyin = abs(yin(nin)-yin(1))/(nin-1._r8)
       ratio = dyinwrap/avgdyin
       if (ratio < 0.9_r8 .or. ratio > 1.1_r8) then
          print*, 'rank=',mpi_rank(),'Lininterp: Bad dyinwrap value =',dyinwrap,&
               ' avg=', avgdyin, yin(1),yin(nin)
          call endrun('interpolate_data')
       end if

       do j=1,nout
          if(increasing) then
             if (yout(j) <= yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) > yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)-yin(nin))/dyinwrap
             end if
          else
             if (yout(j) > yin(1)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin) - yin(nin))/dyinwrap
             else if (yout(j) <= yin(nin)) then
                jjm(j) = nin
                jjp(j) = 1
                wgts(j) = (yin(1)+(cmax-cmin)-yout(j))/dyinwrap
                wgtn(j) = (yout(j)+(cmax-cmin)-yin(nin))/dyinwrap
             end if

          endif
       end do
    end if

    !
    ! Loop though output indices finding input indices and weights
    !
    if(increasing) then
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    else
       do j=1,nout
          do jj=1,nin-1
             if (yout(j).le.yin(jj) .and. yout(j).gt.yin(jj+1)) then
                jjm(j) = jj
                jjp(j) = jj + 1
                wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
                wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
                exit
             end if
          end do
       end do
    end if

    extrap = 100._r8*extrap/real(nout,r8)
    if (extrap.gt.50._r8 .and. working_mode .ne. 'scm') then
       print*, 'interpolate_data:','yout=',minval(yout),maxval(yout),increasing,nout
       print*, 'interpolate_data:','yin=',yin(1),yin(nin)
       print*, 'interpolate_data:',extrap,' % of output grid will have to be extrapolated'
       call endrun('interpolate_data: ')
    end if

    !
    ! Check that interp/extrap points have been found for all outputs
    !
    icount = 0
    do j=1,nout
       if (jjm(j).eq.0 .or. jjp(j).eq.0) icount = icount + 1
       ratio=wgts(j)+wgtn(j)
       if((ratio<0.9_r8.or.ratio>1.1_r8).and.extrap_method.ne.0) then
          print*, j, wgts(j),wgtn(j),jjm(j),jjp(j), increasing,extrap_method
          call endrun('Bad weight computed in LININTERP_init')
       end if
    end do
    if (icount.gt.0) then
       call endrun('LININTERP: Point found without interp indices')
    end if

  end subroutine lininterp_init


  subroutine lininterp1d (arrin, n1, arrout, m1, interp_wgts)
    !-----------------------------------------------------------------------
    !
    ! Purpose: Do a linear interpolation from input mesh to output
    !          mesh with weights as set in lininterp_init.
    !
    !
    ! Author: Jim Edwards
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: n1                 ! number of input latitudes
    integer, intent(in) :: m1                ! number of output latitudes
                                                                                                  
    real(r8), intent(in) :: arrin(n1)    ! input array of values to interpolate
    type(interp_type), intent(in) :: interp_wgts
    real(r8), intent(out) :: arrout(m1) ! interpolated array
                                                                                                  
    !
    ! Local workspace
    !
    integer j                ! latitude indices
    integer, pointer :: jjm(:)
    integer, pointer :: jjp(:)
                                                                                                  
    real(r8), pointer :: wgts(:)
    real(r8), pointer :: wgtn(:)
                                                                                                  
                                                                                                  
    jjm => interp_wgts%jjm
    jjp => interp_wgts%jjp
    wgts =>  interp_wgts%wgts
    wgtn =>  interp_wgts%wgtn
                                                                                                  
    !
    ! Do the interpolation
    !
    do j=1,m1
      arrout(j) = arrin(jjm(j))*wgts(j) + arrin(jjp(j))*wgtn(j)
    end do
                                                                                                  
    return
  end subroutine lininterp1d


  subroutine lininterp2d1d(arrin, n1, n2, arrout, m1, wgt1, wgt2, fldname)
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: n1, n2, m1
    real(r8), intent(in) :: arrin(n1,n2)    ! input array of values to interpolate
    type(interp_type), intent(in) :: wgt1, wgt2
    real(r8), intent(out) :: arrout(m1) ! interpolated array
    character(len=*), intent(in), optional :: fldname(:)
    !
    ! locals
    !
    integer i                ! indices
    integer, pointer :: iim(:), jjm(:)
    integer, pointer :: iip(:), jjp(:)
                                                                                                  
    real(r8), pointer :: wgts(:), wgte(:)
    real(r8), pointer :: wgtn(:), wgtw(:)

    jjm => wgt2%jjm
    jjp => wgt2%jjp
    wgts => wgt2%wgts
    wgtn => wgt2%wgtn
                                                                                                  
    iim => wgt1%jjm
    iip => wgt1%jjp
    wgtw => wgt1%wgts
    wgte => wgt1%wgtn

    do i=1,m1
       arrout(i) = arrin(iim(i),jjm(i))*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i))*wgte(i)*wgts(i) + &
                   arrin(iim(i),jjp(i))*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i))*wgte(i)*wgtn(i)
    end do

                                                                                                  
  end subroutine lininterp2d1d
 

  subroutine lininterp3d2d(arrin, n1, n2, n3, arrout, m1, wgt1, wgt2)
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: n1, n2, n3, m1  ! m1 is to len1 as ncols is to pcols 
    real(r8), intent(in) :: arrin(n1,n2,n3)    ! input array of values to interpolate
    type(interp_type), intent(in) :: wgt1, wgt2
    real(r8), intent(out) :: arrout(n3, m1) ! interpolated array

    !
    ! locals
    !
    integer i, k               ! indices
    integer, pointer :: iim(:), jjm(:)
    integer, pointer :: iip(:), jjp(:)
                                                                                                  
    real(r8), pointer :: wgts(:), wgte(:)
    real(r8), pointer :: wgtn(:), wgtw(:)

    jjm => wgt2%jjm
    jjp => wgt2%jjp
    wgts => wgt2%wgts
    wgtn => wgt2%wgtn
                                                                                                  
    iim => wgt1%jjm
    iip => wgt1%jjp
    wgtw => wgt1%wgts
    wgte => wgt1%wgtn

    do k=1,n3
       do i=1,m1
          arrout(k,i) = arrin(iim(i),jjm(i),k)*wgtw(i)*wgts(i)+arrin(iip(i),jjm(i),k)*wgte(i)*wgts(i) + &
               arrin(iim(i),jjp(i),k)*wgtw(i)*wgtn(i)+arrin(iip(i),jjp(i),k)*wgte(i)*wgtn(i)
       end do
    end do
                                                                                                  
  end subroutine lininterp3d2d
                                                                                                  

  subroutine lininterp_finish(interp_wgts)
    type(interp_type) :: interp_wgts

    deallocate(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)

    nullify(interp_wgts%jjm, &
         interp_wgts%jjp, &
         interp_wgts%wgts, &
         interp_wgts%wgtn)
  end subroutine lininterp_finish


  subroutine lininterp_original (arrin, yin, nlev, nlatin, arrout, &
       yout, nlatout)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Do a linear interpolation from input mesh defined by yin to output
    !          mesh defined by yout.  Where extrapolation is necessary, values will
    !          be copied from the extreme edge of the input grid.  Vectorization is over
    !          the vertical (nlev) dimension.
    ! 
    ! Method: Check validity of input, then determine weights, then do the N-S interpolation.
    ! 
    ! Author: Jim Rosinski
    ! Modified: Jim Edwards so that there is no requirement of monotonacity on the yout array
    !
    !-----------------------------------------------------------------------
    implicit none
    !-----------------------------------------------------------------------
    !
    ! Arguments
    !
    integer, intent(in) :: nlev                   ! number of vertical levels
    integer, intent(in) :: nlatin                 ! number of input latitudes
    integer, intent(in) :: nlatout                ! number of output latitudes

    real(r8), intent(in) :: arrin(nlev,nlatin)    ! input array of values to interpolate
    real(r8), intent(in) :: yin(nlatin)           ! input mesh
    real(r8), intent(in) :: yout(nlatout)         ! output mesh

    real(r8), intent(out) :: arrout(nlev,nlatout) ! interpolated array
    !
    ! Local workspace
    !
    integer j, jj              ! latitude indices
    integer jjprev             ! latitude indices
    integer k                  ! level index
    integer icount             ! number of values

    real(r8) extrap            ! percent grid non-overlap
    !
    ! Dynamic
    !
    integer :: jjm(nlatout)
    integer :: jjp(nlatout)

    real(r8) :: wgts(nlatout)
    real(r8) :: wgtn(nlatout)
    !
    ! Check validity of input coordinate arrays: must be monotonically increasing,
    ! and have a total of at least 2 elements
    !
    if (nlatin.lt.2) then
       call endrun('LININTERP: Must have at least 2 input points for interpolation')
    end if

    icount = 0
    do j=1,nlatin-1
       if (yin(j).gt.yin(j+1)) icount = icount + 1
    end do


    if (icount.gt.0) then
       call endrun('LININTERP: Non-monotonic coordinate array(s) found')
    end if
    !
    ! Initialize index arrays for later checking
    !
    do j=1,nlatout
       jjm(j) = 0
       jjp(j) = 0
    end do
    !
    ! For values which extend beyond N and S boundaries, set weights
    ! such that values will just be copied.
    !
    extrap = 0._r8

    do j=1,nlatout
       if (yout(j).le.yin(1)) then
          jjm(j) = 1
          jjp(j) = 1
          wgts(j) = 1._r8
          wgtn(j) = 0._r8
          extrap=extrap+1._r8
       else if (yout(j).gt.yin(nlatin)) then
          jjm(j) = nlatin
          jjp(j) = nlatin
          wgts(j) = 1._r8
          wgtn(j) = 0._r8
          extrap=extrap+1._r8
       endif
    end do

    !
    ! Loop though output indices finding input indices and weights
    !
    do j=1,nlatout
       do jj=1,nlatin-1
          if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
             jjm(j) = jj
             jjp(j) = jj + 1
             wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
             wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
             exit
          end if
       end do
    end do
    !
    ! Check that interp/extrap points have been found for all outputs
    !
    icount = 0
    do j=1,nlatout
       if (jjm(j).eq.0 .or. jjp(j).eq.0) then
          icount = icount + 1
       end if
    end do
    if (icount.gt.0) then
       call endrun('LININTERP: Point found without interp indices')
    end if
    !
    ! Do the interpolation
    !
    do j=1,nlatout
       do k=1,nlev
          arrout(k,j) = arrin(k,jjm(j))*wgts(j) + arrin(k,jjp(j))*wgtn(j)
       end do
    end do

    return
  end subroutine lininterp_original


 
end module grist_horizontal_interpolate
