!#include <misc.h>
!#include <params.h>

!================================================
!  Adapted from cam31 version for PhysW
!  2021-01-02, yizhang
!================================================

module cam3_cloud_fraction
  !
  ! Module for cloud fraction routines.
  !
  ! $Id: cloud_fraction.F90,v 1.1.4.9 2004/12/10 16:25:27 jeff Exp $
  !
  !use shr_kind_mod, only: r8 => shr_kind_r8
  use grist_constants, only: r8, i4,gravit=>gravity, cappa, rair=>rdry, zero,latvap,latice,cp,rvap,epsilo
  use wv_saturation,   only: aqsat, gestbl
  use grist_wv_saturation, only: qsat, wv_sat_init

  implicit none

  private
  save
  !
  ! Public interfaces
  !
  public cam3_cldfrc_init        ! Inititialization of cloud_fraction run-time parameters
  public cam3_cldfrc             ! Computation of cloud fraction
  public pcond                   ! Computation of cloud fraction
  !
  ! Private data
  !
  real(r8) rhminl                ! minimum rh for low stable clouds
  real(r8) rhminh                ! minimum rh for high stable clouds
  real(r8) sh1,sh2               ! parameters for shallow convection cloud fraction
  real(r8) dp1,dp2               ! parameters for deep convection cloud fraction
  real(r8) premit                ! top pressure bound for mid level cloud
  integer(i4)          :: k700   ! model level nearest 700 mb
  real(r8), public, parameter :: tmelt = 273.16_r8         ! Freezing point of water
  real(r8) :: rh2o 
contains  

  subroutine cam3_cldfrc_init(plev, hypm)
    integer(i4), intent(in) :: plev
    real(r8),    intent(in) :: hypm(plev)
    integer(i4) :: k
    real(r8) :: tmn, tmx, trice, epsil, cpair, tmeltx
    logical  :: ip

     tmn    = 173.16_r8
     tmx    = 375.16_r8
     trice  = 20.00_r8
     ip     =.true.
     epsil  = epsilo ! ratio of h2o to dry air molecular weights 
     rh2o   = rvap
     cpair  = cp
     tmeltx = tmelt 

     call gestbl(tmn     ,tmx     ,trice   ,ip      ,epsil   , &
                 latvap  ,latice  ,rh2o    ,cpair   ,tmeltx   )
     call wv_sat_init
    !
    ! Purpose:
    ! Initialize cloud fraction run-time parameters
    !
    ! Author: J. McCaa
    !    
!    use dycore, only: dycore_is, get_resolution
!    use ppgrid, only: pver          

!     if ( dycore_is ('LR') ) then
!        if ( get_resolution() == '1x1.25' ) then
          rhminl = .88_r8
          rhminh = .77_r8
#ifdef AMIPW_CLIMATE
          rhminl = .88_r8
          rhminh = .72_r8
#endif
          sh1 = 0.04_r8
          sh2 = 500.0_r8
          dp1 = 0.10_r8
#ifdef AMIPW_CLIMATE
          dp1 = 0.15_r8
#endif
          dp2 = 500.0_r8
          premit = 750.e2_r8  ! top of area defined to be mid-level cloud
!        elseif ( get_resolution() == '4x5' .and. pver == 66 ) then
!           rhminl = .90
!           rhminh = .90
!           sh1 = 0.04
!           sh2 = 500.0
!           dp1 = 0.10
!           dp2 = 500.0
!           premit = 750.e2  ! top of area defined to be mid-level cloud
!        else
!           rhminl = .90
!           rhminh = .80
!           sh1 = 0.04
!           sh2 = 500.0
!           dp1 = 0.10
!           dp2 = 500.0
!           premit = 750.e2  ! top of area defined to be mid-level cloud
!        endif
!     else
!        if ( get_resolution() == 'T85' ) then
!           rhminl = .91
!           rhminh = .70
!           sh1 = 0.07
!           sh2 = 500.0
!           dp1 = 0.14
!           dp2 = 500.0
!           premit = 250.e2  ! top of area defined to be mid-level cloud
!        elseif ( get_resolution() == 'T31' ) then
!           rhminl = .88
!           rhminh = .80
!           sh1 = 0.07
!           sh2 = 500.0
!           dp1 = 0.14
!           dp2 = 500.0
!           premit = 750.e2  ! top of area defined to be mid-level cloud
!        else
!           rhminl = .90
!           rhminh = .80
!           sh1 = 0.07
!           sh2 = 500.0
!           dp1 = 0.14
!           dp2 = 500.0
!           premit = 750.e2  ! top of area defined to be mid-level cloud
!        endif
!     endif

!
! Find vertical level nearest 700 mb
!
   k700 = 1
   do k=1,plev-1
      if (hypm(k) < 7.e4 .and. hypm(k+1) >= 7.e4) then
         if (7.e4-hypm(k) < hypm(k+1)-7.e4) then
            k700 = k
         else
            k700 = k + 1
         end if
         goto 20
      end if
   end do

   20 continue

   return
  end subroutine cam3_cldfrc_init

  subroutine cam3_cldfrc(ncol, pver, dindex)
    use grist_PhysW_data_structure,   only: pstate_wrf, psurf_wrf
    use grist_physics_data_structure, only: pstate
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Compute cloud fraction 
    ! 
    ! 
    ! Method: 
    ! This calculate cloud fraction using a relative humidity threshold
    ! The threshold depends upon pressure, and upon the presence or absence 
    ! of convection as defined by a reasonably large vertical mass flux 
    ! entering that layer from below.
    ! 
    ! Author: Many. Last modified by Jim McCaa
    ! 
    !-----------------------------------------------------------------------
!    use ppgrid
!    use physconst, only: cappa, gravit, rair, tmelt
!    use cldconst
!    use phys_grid,     only: get_rlat_all_p, get_rlon_all_p
!    use dycore,        only: dycore_is, get_resolution
    implicit none

    real(r8), parameter :: pnot = 1.e5       ! reference pressure
    real(r8), parameter :: lapse = 6.5e-3    ! U.S. Standard Atmsophere lapse rate
    real(r8), parameter :: premib = 750.e2   ! bottom pressure bound of middle cloud
    real(r8), parameter :: pretop = 1.0e2    ! pressure bounding high cloud
    !
    ! Arguments
    !
!    integer, intent(in) :: lchnk                 ! chunk identifier
    integer, intent(in)   :: ncol                 ! number of atmospheric columns
    integer, intent(in)   :: pver                 ! number of atmospheric columns, GRIST full level
    integer, intent(in)   :: dindex               ! 0 or 1 to perturb rh
!
    real(r8) :: pmid(ncol,pver)     ! midpoint pressures
    real(r8) :: pdel(ncol,pver)     ! pressure depth of layer
    real(r8) :: temp(ncol,pver)     ! temperature
    real(r8) :: q(ncol,pver)        ! specific humidity
    real(r8) :: phis(ncol)          ! surface geopotential
    real(r8) :: cmfmc(ncol,pver)    ! convective mass flux--m sub c
    real(r8) :: ts(ncol)            ! surface temperature
    real(r8) :: sst(ncol)           ! sea surface temperature
    real(r8) :: ps(ncol)            ! surface pressure
    real(r8) :: landfrac(ncol)      ! Land fraction
    real(r8) :: ocnfrac(ncol)       ! Ocean fraction
    real(r8) :: snowh(ncol)         ! snow depth (liquid water equivalent)
    real(r8) :: cloud(ncol,pver)    ! cloud fraction
    real(r8) :: concld(ncol,pver)   ! convective cloud cover
    real(r8) :: cldst(ncol,pver)    ! cloud fraction
!
    real(r8) :: rhu00(ncol,pver)     ! RH threshold for cloud
    real(r8) :: relhum(ncol,pver)    ! RH 
    real(r8) :: cmfmc2(ncol,pver)    ! shallow convective mass flux--m sub c, assume zero, all are deep
    !
    !---------------------------Local workspace-----------------------------
    !
    real(r8) cld                   ! intermediate scratch variable (low cld)
    real(r8) dthdpmn(ncol)         ! most stable lapse rate below 750 mb
    real(r8) dthdp                 ! lapse rate (intermediate variable)
    real(r8) es(ncol,pver)        ! saturation vapor pressure
    real(r8) qs(ncol,pver)        ! saturation specific humidity
    real(r8) rhwght                ! weighting function for rhlim transition
    real(r8) rh(ncol,pver)        ! relative humidity
    real(r8) rhdif                 ! intermediate scratch variable
    real(r8) strat                 ! intermediate scratch variable
    real(r8) theta(ncol,pver)     ! potential temperature
    real(r8) rhlim                 ! local rel. humidity threshold estimate
    real(r8) coef1                 ! coefficient to convert mass flux to mb/d
    real(r8) clrsky(ncol)         ! temporary used in random overlap calc
    real(r8) rpdeli(ncol,pver-1) ! 1./(pmid(k+1)-pmid(k))
    real(r8) rhpert                !the specified perturbation to rh
    real(r8) deepcu                ! deep convection cloud fraction
    real(r8) shallowcu             ! shallow convection cloud fraction

    logical cldbnd(ncol)          ! region below high cloud boundary

    integer i, ierror, k           ! column, level indices
    integer kp1
    integer kdthdp(ncol)
    integer numkcld                ! number of levels in which to allow clouds

    real(r8) thetas(ncol)                    ! ocean surface potential temperature
!    real(r8) :: clat(ncol)                   ! current latitudes(radians)
!    real(r8) :: clon(ncol)                   ! current longitudes(radians)
    !
    ! Statement functions
    !
!    logical land
    integer(i4) :: icell, ilev

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! CAM from top to bottom: 1->pver
! WRF from bottom to top
! Reverse data structure
!
    do icell = 1, ncol
       do ilev = 1, pver
          pmid (icell,ilev)  = pstate_wrf%p_phy(icell,pver+1-ilev,1)   ! midpoint pressures
          temp (icell,ilev)  = pstate_wrf%t_phy(icell,pver+1-ilev,1)   ! temperature
          q    (icell,ilev)  = pstate_wrf%moist(icell,pver+1-ilev,1,1) ! specific humidity
          cmfmc(icell,ilev)  = pstate_wrf%cmfmc(icell,pver+1-ilev,1)   ! convective mass flux--m sub c
       end do
       phis    (icell)       = psurf_wrf%ht (icell,1)  * gravit         ! 
       ts      (icell)       = psurf_wrf%tsk(icell,1)                    ! surface temperature
       sst     (icell)       = pstate%sst_at_pc_surface%f(icell)         ! sea surface temperature
       ps      (icell)       = psurf_wrf%psfc(icell,1)                   ! surface pressure
       landfrac(icell)       = pstate%landfrac_at_pc_surface%f(icell)    ! Land fraction
       ocnfrac (icell)       = pstate%ocnfrac_at_pc_surface%f(icell)     ! Ocean fraction
       snowh   (icell)       = zero                                      ! snow depth (liquid water equivalent)
    end do
    cmfmc2 = zero
    !cmfmc  = zero
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    land(i) = nint(landfrac(i)) == 1
!    call get_rlat_all_p(lchnk, ncol, clat)
!    call get_rlon_all_p(lchnk, ncol, clon)

    !==================================================================================
    ! PHILOSOPHY OF PRESENT IMPLEMENTATION
    !
    ! There are three co-existing cloud types: convective, inversion related low-level
    ! stratocumulus, and layered cloud (based on relative humidity).  Layered and 
    ! stratocumulus clouds do not compete with convective cloud for which one creates 
    ! the most cloud.  They contribute collectively to the total grid-box average cloud 
    ! amount.  This is reflected in the way in which the total cloud amount is evaluated 
    ! (a sum as opposed to a logical "or" operation)
    !
    !==================================================================================
    ! set defaults for rhu00
    rhu00(:,:) = 2.0
    ! define rh perturbation in order to estimate rhdfda
    rhpert = 0.01 

    !
    ! Evaluate potential temperature and relative humidity
    !
    call aqsat(temp    ,pmid    ,es      ,qs      ,ncol   , &
         ncol    ,pver    ,1       ,pver    )
    do k=1,pver
       do i=1,ncol
!================================================
! new added at 20210920 for checking
! call qsat(temp(i,k), pmid(i,k), es(i,k), qs(i,k))
! very-close results as CAM3 version
!================================================
          theta(i,k)  = temp(i,k)*(pnot/pmid(i,k))**cappa
          rh(i,k)     = q(i,k)/qs(i,k)*(1.0+float(dindex)*rhpert)
          !
          !  record relhum, rh itself will later be modified related with concld
          !
          relhum(i,k) = rh(i,k)
          cloud(i,k)  = 0.
          cldst(i,k)  = 0.
          concld(i,k) = 0.
       end do
    end do
    !
    ! Initialize other temporary variables
    !
    ierror = 0
    do i=1,ncol
       ! Adjust thetas(i) in the presence of non-zero ocean heights.
       ! This reduces the temperature for positive heights according to a standard lapse rate.
       if(ocnfrac(i).gt.0.01) thetas(i)  = &
            ( sst(i) - lapse * phis(i) / gravit) * (pnot/ps(i))**cappa
       if(ocnfrac(i).gt.0.01.and.sst(i).lt.260.) ierror = i
       !clc(i) = 0.0
    end do
    coef1 = gravit*864.0    ! conversion to millibars/day

    if (ierror > 0) then
       write(6,*) 'COLDSST: encountered in cldfrc:', ierror,ocnfrac(ierror),sst(ierror)
    endif

    do k=1,pver-1
       do i=1,ncol
          rpdeli(i,k) = 1./(pmid(i,k+1) - pmid(i,k))
       end do
    end do

    !
    ! Estimate of local convective cloud cover based on convective mass flux
    ! Modify local large-scale relative humidity to account for presence of 
    ! convective cloud when evaluating relative humidity based layered cloud amount
    !
    do k=1,pver
       do i=1,ncol
          concld(i,k) = 0.0
       end do
    end do
    !
    ! cloud mass flux in SI units of kg/m2/s; should produce typical numbers of 20%
    ! shallow and deep convective cloudiness are evaluated separately (since processes
    ! are evaluated separately) and summed
    !   
#ifndef PERGRO
    ! change from pverp to pver
! tdk does not seperate dp and sh, so just use a unfied equation here
    !do k=1,pver-1
    do k=1,pver
       do i=1,ncol
          !shallowcu = max(0.0,min(sh1*log(1.0+sh2*cmfmc2(i,k+1)),0.30))
          !deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k+1)-cmfmc2(i,k+1))),0.60))
          !-------------------------------------------
          !shallowcu = max(0.0,min(sh1*log(1.0+sh2*cmfmc2(i,k)),0.30))
          !deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k)-cmfmc2(i,k))),0.60))
          !concld(i,k) = min(shallowcu + deepcu,0.80)
          deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k)-cmfmc2(i,k))),0.80))
          concld(i,k) = min(deepcu,0.80)
#ifdef AMIPW_CLIMATE
          deepcu = max(0.0,min(dp1*log(1.0+dp2*(cmfmc(i,k)-cmfmc2(i,k))),0.80))
          concld(i,k) = min(deepcu,0.80)
#endif
          rh(i,k) = (rh(i,k) - concld(i,k))/(1.0 - concld(i,k))
       end do
    end do
#endif
    !==================================================================================
    !
    !          ****** Compute layer cloudiness ******
    !
    !====================================================================
    ! Begin the evaluation of layered cloud amount based on (modified) RH 
    !====================================================================
    !
    numkcld = pver
    do k=2,numkcld
       kp1 = min(k + 1,pver)
       do i=1,ncol
          !
          cldbnd(i) = pmid(i,k).ge.pretop
          !
          if ( pmid(i,k).ge.premib ) then
             !==============================================================
             ! This is the low cloud (below premib) block
             !==============================================================
             ! enhance low cloud activation over land with no snow cover
             if (nint(landfrac(i)).eq.1 .and. (snowh(i) <= 0.000001)) then
                rhlim = rhminl - 0.10
             else
                rhlim = rhminl
             endif
             !
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          else if ( pmid(i,k).lt.premit ) then
             !==============================================================
             ! This is the high cloud (above premit) block
             !==============================================================
             !
             rhlim = rhminh
             !
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          else
             !==============================================================
             ! This is the middle cloud block
             !==============================================================
             !
             !       linear rh threshold transition between thresholds for low & high cloud
             !
             rhwght = (premib-(max(pmid(i,k),premit)))/(premib-premit)
             
             if (nint(landfrac(i)).eq.1 .and. (snowh(i) <= 0.000001)) then
                rhlim = rhminh*rhwght + (rhminl - 0.10)*(1.0-rhwght)
             else
                rhlim = rhminh*rhwght + rhminl*(1.0-rhwght)
             endif
             rhdif = (rh(i,k) - rhlim)/(1.0_r8-rhlim)
             cloud(i,k) = min(0.999_r8,(max(rhdif,0.0_r8))**2)
          end if
          !==================================================================================
          ! WE NEED TO DOCUMENT THE PURPOSE OF THIS TYPE OF CODE (ASSOCIATED WITH 2ND CALL)
          !==================================================================================
          !      !
          !      ! save rhlim to rhu00, it handles well by itself for low/high cloud
          !      !
          rhu00(i,k)=rhlim
          !==================================================================================
       end do
       !
       ! Final evaluation of layered cloud fraction
       !
    end do
    !
    ! Add in the marine strat
    ! MARINE STRATUS SHOULD BE A SPECIAL CASE OF LAYERED CLOUD
    ! CLOUD CURRENTLY CONTAINS LAYERED CLOUD DETERMINED BY RH CRITERIA
    ! TAKE THE MAXIMUM OF THE DIAGNOSED LAYERED CLOUD OR STRATOCUMULUS
    !
    !===================================================================================
    !
    !  SOME OBSERVATIONS ABOUT THE FOLLOWING SECTION OF CODE (missed in earlier look)
    !  K700 IS SET AS A CONSTANT BASED ON HYBRID COORDINATE: IT DOES NOT DEPEND ON 
    !  LOCAL PRESSURE; THERE IS NO PRESSURE RAMP => LOOKS LEVEL DEPENDENT AND 
    !  DISCONTINUOUS IN SPACE (I.E., STRATUS WILL END SUDDENLY WITH NO TRANSITION)
    !
    !  IT APPEARS THAT STRAT IS EVALUATED ACCORDING TO KLEIN AND HARTMANN; HOWEVER,
    !  THE ACTUAL STRATUS AMOUNT (CLDST) APPEARS TO DEPEND DIRECTLY ON THE RH BELOW
    !  THE STRONGEST PART OF THE LOW LEVEL INVERSION.  
    !
    !==================================================================================
    !
    ! Find most stable level below 750 mb for evaluating stratus regimes
    !
    do i=1,ncol
       ! Nothing triggers unless a stability greater than this minimum threshold is found
       dthdpmn(i) = -0.125
       kdthdp(i) = 0
    end do
    !
    do k=2,pver
       do i=1,ncol
          if (pmid(i,k) >= premib .and. ocnfrac(i).gt. 0.01) then
             ! I think this is done so that dtheta/dp is in units of dg/mb (JJH)
             dthdp = 100.0*(theta(i,k) - theta(i,k-1))*rpdeli(i,k-1)
             if (dthdp < dthdpmn(i)) then
                dthdpmn(i) = dthdp
                kdthdp(i) = k     ! index of interface of max inversion
             end if
          end if
       end do
    end do

    ! Also check between the bottom layer and the surface
    ! Only perform this check if the criteria were not met above

    do i = 1,ncol
       if ( kdthdp(i) .eq. 0 .and. ocnfrac(i).gt.0.01) then
          dthdp = 100.0 * (thetas(i) - theta(i,pver)) / (ps(i)-pmid(i,pver))
          if (dthdp < dthdpmn(i)) then
             dthdpmn(i) = dthdp
             kdthdp(i) = pver     ! index of interface of max inversion
          endif
       endif
    enddo

    do i=1,ncol
       if (kdthdp(i) /= 0) then
          k = kdthdp(i)
          kp1 = min(k+1,pver)
          ! Note: strat will be zero unless ocnfrac > 0.01
          strat = min(1._r8,max(0._r8, ocnfrac(i) * ((theta(i,k700)-thetas(i))*.057-.5573) ) )
          !
          ! assign the stratus to the layer just below max inversion
          ! the relative humidity changes so rapidly across the inversion
          ! that it is not safe to just look immediately below the inversion
          ! so limit the stratus cloud by rh in both layers below the inversion
          !
          cldst(i,k) = min(strat,max(rh(i,k),rh(i,kp1)))
       end if
    end do
    !
    ! AGGREGATE CLOUD CONTRIBUTIONS (cldst should be zero everywhere except at level kdthdp(i))
    !
    do k=1,pver
       do i=1,ncol
          !
          !       which is greater; standard layered cloud amount or stratocumulus diagnosis
          !
          cloud(i,k) = max(cloud(i,k),cldst(i,k))
          !
          !       add in the contributions of convective cloud (determined separately and accounted
          !       for by modifications to the large-scale relative humidity.
          !
          cloud(i,k) = min(cloud(i,k)+concld(i,k), 1.0_r8)
          cldst(i,k) = min(max(cloud(i,k)-concld(i,k),0._r8),1._r8)
       end do
    end do

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! from top to bottom: 1->pver, reverse wrf data structure
!
    do ilev = 1, pver
       do icell = 1, ncol
          pstate_wrf%cldfra(icell,pver+1-ilev,1)  = cloud (icell,ilev)
          pstate_wrf%cldcum(icell,pver+1-ilev,1)  = concld(icell,ilev)
          pstate_wrf%cldstr(icell,pver+1-ilev,1)  = cldst (icell,ilev)
       end do
    end do
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    return
  end subroutine cam3_cldfrc

  subroutine pcond(plond, plon, plev, tdt     ,pmid    ,pdel    ,t       ,q       , &
                   qc      ,precl   )
!-----------------------------------------------------------------------
!
! Calculate large scale condensation
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      L. Buja, Jun 1992, Feb 1996
! Reviewed:          J. Hack, G. Taylor, Aug 1992
!                    J. Hack, Feb 1996 
!
!-----------------------------------------------------------------------
!
! $Id: cond.F,v 1.1 1998/04/01 07:21:30 ccm Exp $
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!
! Input arguments
!
      integer(i4), intent(in)::  plond, plon, plev
      real(r8), intent(in)   ::  tdt                    ! Physics time step (2 delta t)
      real(r8), intent(in)   ::  pmid(plond,plev)       ! Pressure at layer midpoints
      real(r8), intent(in)   ::  pdel(plond,plev)       ! Delta p at each model level
!
! Output arguments
!
      real(r8), intent(inout)   ::  t(plond,plev)       ! Temperature
      real(r8), intent(inout)   ::  q(plond,plev)       ! Specific humidity
      real(r8), intent(inout)   ::  qc(plond,plev)      ! Moisture tendency due to rainout
      real(r8), intent(out)     ::  precl(plond)        ! Large-scale precipitation rate (m/s)
!
!---------------------------Local variables-----------------------------
!
      real(r8)   ::  absqs                  ! Intermediate quantity
      real(r8)   ::  denom                  ! Intermediate quantity
      real(r8)   ::  dqsdt                  ! Change of qsat with respect to temp.
      real(r8)   ::  est(plond,plev)        ! Saturation vapor pressure
      real(r8)   ::  omeps                  ! 1 - 0.622
      real(r8)   ::  qsat(plond,plev)       ! Saturation specific humidity
      real(r8)   ::  rain(plond)            ! Rain (units of kg/m^2 s)
      real(r8)   ::  rga                    ! Reciprocal gravitatnl acceleration
      real(r8)   ::  rhm1                   ! RH - saturation RH
      real(r8)   ::  zqcd(plond)            ! Intermed quantity (several actually)
      real(r8)   ::  zqdt                   ! Reciprocal of tdt
      real(r8)   ::  cndwtr(plond,plev)     ! Water condensation rate (kg/m**2/s)
      real(r8)   ::  ke                     ! `disposable parameter' in evaporation
      real(r8)   ::  evap                   ! Water evaporation rate
      real(r8)   ::  relhum                 ! Relative humidity
      real(r8)   ::  dpovrg                 ! deltap/grav
      integer(i4)::  i                   ! Longitude index
      integer(i4)::  jiter               ! Iteration counter
      integer(i4)::  k                   ! Vertical index
      real(r8)   ::  cldcp, rhoh2o_comadj, clrh2o

! some constants
      cldcp         = latvap/cp
      rhoh2o_comadj = 1000._r8 ! density of liquid water, kg/m3
      clrh2o        = latvap/rh2o
!
!-----------------------------------------------------------------------
!
      rga   = 1./gravit
      zqdt  = 1./tdt
      omeps = 1. - epsilo
!
! First diagnose condensation rate due to stable processes
! Update column T and Q (evaporation process is `time-split')
! Condensation calculation is hard-wired for two iterations
!
      do k=1,plev
        do i=1,plon
          cndwtr(i,k) = 0.0
        end do
      end do
!
      do jiter=1,2
        call aqsat(t       ,pmid    ,est     ,qsat    ,plond   , &
                  plon    ,plev    ,1       ,plev    )
        do k=1,plev
!
! Calculate condensation-rate and new t- and q-values
!
          do i=1,plon
!
! Use of critical saturation vapor pressure requires coefficient on the
! term omeps*est(i,k) in the next statement (e.g. omeps*est(i,k)*escrit)
! Corresponding changes must also be incorporated into estabv.for (e.g.,
! terms est(i,k) in qsat evaluation become escrit*est(i,k))
!
            denom   = (pmid(i,k) - omeps*est(i,k))*t(i,k)**2
            dqsdt   = clrh2o*qsat(i,k)*pmid(i,k)/denom
            absqs   = abs(qsat(i,k))
            rhm1    = q(i,k)/qsat(i,k) - 1.
            zqcd(i) = max(absqs*rhm1/(1. + cldcp*dqsdt),0.)
            if (q(i,k) .lt. 0.0) zqcd(i) = 0.
            q(i,k)  = q(i,k) - zqcd(i)
            t(i,k)  = t(i,k) + zqcd(i)*cldcp
            cndwtr(i,k) = cndwtr(i,k) + zqcd(i)*pdel(i,k)*rga*zqdt
            qc    (i,k) = qc(i,k)     + zqcd(i)*zqdt
          end do
        end do
      end do
!
! Initialize rain vector (will be updated as rain falls through column)
!
      do i=1,plon
        rain(i) = max(cndwtr(i,1),0.0)
      end do
      call aqsat(t       ,pmid    ,est     ,qsat    ,plond   , &
                 plon    ,plev    ,1       ,plev    )
!
! Evaporate condensate on the way down (see Sundqvist, 1988: Physically
! Based Modelling ..., pp 433-461, Schlesinger, Ed., Kluwer Academic)
! variable evap has units of 1/s; variable rain has units of kg/m**2/s
! rain is used to accumuluate unevaporated rain water on the way down
!
      ke = 1.0e-5                     ! set in common block in final code
      do k=2,plev
        do i=1,plon
          dpovrg  = pdel(i,k)*rga
          relhum  = q(i,k)/qsat(i,k)
          evap    = max(ke*(1.0 - relhum)*sqrt(rain(i)), 0.0)
          evap    = min(evap, (qsat(i,k)-q(i,k))/tdt)
          evap    = min(rain(i)/dpovrg,evap)
          qc(i,k) = qc(i,k) - evap
          q(i,k)  = q(i,k) + evap*tdt
          t(i,k)  = t(i,k) - evap*tdt*cldcp
          rain(i) = max(rain(i) - evap*dpovrg + cndwtr(i,k),0.0)
        end do
      end do
      do i=1,plon
        precl(i) = rain(i)/rhoh2o_comadj
      end do
!
      return
    end subroutine pcond

end module cam3_cloud_fraction