Module double_plume_cape_zm

  use grist_constants,                    only: r8, epsilo, tmelt, &
                                                rdry, cpwv, cpliq, rvap
  use grist_nml_module,                   only: nlev, nlevp
  use cloud_fraction,                     only: cldfrc_fice
  use grist_handle_error,                 only: endrun

implicit none

   public :: double_plume_cape_dilute,  &
             double_plume_zm_init

   private

   real(r8) :: rl         ! wg latent heat of vaporization.
   real(r8) :: latice
   real(r8) :: cpres      ! specific heat at constant pressure in j/kg-degk.
   real(r8) :: tfreez
   real(r8) :: eps1

!moved from moistconvection.F90
   real(r8) :: rgrav       ! reciprocal of grav
   real(r8) :: rgas        ! gas constant for dry air
   real(r8) :: grav        ! = gravity
   real(r8) :: cp          ! = cpres = cp
   
   integer  limcnv       ! top interface level limit for convection

   real(r8),parameter ::  tiedke_add = 0.5_r8   

contains


subroutine double_plume_zm_init(limcnv_in, cp_in, grav_in, rl_in, latice_in)


   integer, intent(in)           :: limcnv_in       ! top interface level limit for convection
   real(r8), intent(in)          :: cp_in
   real(r8), intent(in)          :: grav_in
   real(r8), intent(in)          :: rl_in
   real(r8), intent(in)          :: latice_in

   ! local variables

   ! Initialization of ZM constants
   limcnv = limcnv_in
   tfreez = tmelt
   eps1   = epsilo
   rl     = rl_in
   latice = latice_in
   cpres  = cp_in
   rgrav  = 1.0_r8/grav_in
   rgas   = rdry
   grav   = grav_in
   cp     = cpres

end subroutine double_plume_zm_init

subroutine double_plume_cape_dilute(ncol    , &
                  q       ,t       ,p       ,z       ,pf      , &
                  cape    , pcape, pblt,   msg     , tpert   )

! 
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! input arguments
!
   integer, intent(in) :: ncol                  ! number of atmospheric columns

   real(r8), intent(in) :: q(nlev, ncol)        ! spec. humidity
   real(r8), intent(in) :: t(nlev, ncol)        ! temperature
   real(r8), intent(in) :: p(nlev, ncol)        ! pressure
   real(r8), intent(in) :: z(nlev, ncol)        ! height
   real(r8), intent(in) :: pf(nlevp, ncol)     ! pressure at interfaces
   real(r8), intent(in) :: pblt(ncol)          ! index of pbl depth
   real(r8), intent(in) :: tpert(ncol)         ! perturbation temperature by pbl processes
   integer , intent(in) :: msg
!
! output arguments
!
   
   real(r8), intent(out) :: cape(ncol)          ! convective aval. pot. energy.
   real(r8), intent(out) :: pcape(ncol)
!
!--------------------------Local Variables------------------------------
!
   integer lcl(ncol)        !
   integer lel(ncol)        !
   integer lon(ncol)        ! level of onset of deep convection
   integer mx(ncol)         ! level of max moist static energy
   real(r8) tp(nlev, ncol)       ! parcel temperature
   real(r8) qstp(nlev, ncol)     ! saturation mixing ratio of parcel (only above lcl, just q below).
   real(r8) tl(ncol)            ! parcel temperature at lcl
   real(r8) capeten(ncol, 5)     ! provisional value of cape
   real(r8) pcapeten(ncol, 5)
   real(r8) tv(nlev, ncol)       !
   real(r8) tpv(nlev, ncol)      !
   real(r8) buoy(nlev, ncol)

   real(r8) a1(ncol)
   real(r8) a2(ncol)
   real(r8) estp(ncol)
   real(r8) pl(ncol)
   real(r8) plexp(ncol)
   real(r8) hmax(ncol)
   real(r8) hmn(ncol)
   real(r8) y(ncol)

   logical plge600(ncol)
   integer knt(ncol)
   integer lelten(ncol, 5)

   real(r8) e

   integer i
   integer k
   integer n

!---------------LiXH does not use PERGRO------------->
!#ifdef PERGRO
!   real(r8) rhd
!#endif
!<--------------LiXH does not use PERGRO--------------

   do n = 1,5
      do i = 1,ncol
         lelten(i,n) = nlev
         capeten(i,n) = 0._r8
         pcapeten(i,n) = 0._r8
      end do
   end do
!
   do i = 1,ncol
      lon(i) = nlev
      knt(i) = 0
      lel(i) = nlev
      mx(i) = lon(i)
      cape(i) = 0._r8
      pcape(i) = 0._r8
      hmax(i) = 0._r8
   end do

   tp(:,:ncol) = t(:,:ncol)
   qstp(:,:ncol) = q(:,:ncol)

!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
   tv(:,:ncol) = t(:,:ncol) *(1._r8+1.608_r8*q(:,:ncol))/ (1._r8+q(:,:ncol))
   tpv(:,:ncol) = tv(:,:ncol)
   buoy(:,:ncol) = 0._r8

!
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
!
!---------------LiXH does not use PERGRO------------->
!#ifdef PERGRO
!   do k = nlev,msg + 1,-1
!      do i = 1,ncol
!         hmn(i) = cp*t(k,i) + grav*z(k,i) + rl*q(k,i)
!!
!! Reset max moist static energy level when relative difference exceeds 1.e-4
!!
!         rhd = (hmn(i) - hmax(i))/(hmn(i) + hmax(i))
!        ! if (k >= nint(pblt(i)) .and. k <= lon(i) .and. rhd > -1.e-4_r8) then
!          if (pf(k,i) >= 800 .and. k <= lon(i) .and. rhd > -1.e-4_r8) then
!            hmax(i) = hmn(i)
!            mx(i) = k
!         end if
!      end do
!   end do
!#else
   do i = 1,ncol
      do k = nlev,msg + 1,-1
         hmn(i) = cp*t(k,i) + grav*z(k,i) + rl*q(k,i)
       !  if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
       !  if (pf(k,i) >= 800._r8 .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
          if (pf(k,i) >= 600._r8 .and. k <= lon(i) .and. hmn(i) > hmax(i)) then     !LiXH: ULL in Xie et al. (2018)
            hmax(i) = hmn(i)
            mx(i) = k
         end if
      end do
   end do
!#endif
!<--------------LiXH does not use PERGRO--------------

! LCL dilute calculation - initialize to mx(i)
! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
! Original code actually sets LCL as level above wher condensate forms.
! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.

   do i = 1,ncol ! Initialise LCL variables.
      lcl(i) = mx(i)
      tl(i) = t(mx(i),i)
      pl(i) = p(mx(i),i)
   end do

!
! main buoyancy calculation.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!!   RBN 9/9/04   !!!
   !write(iulog,*) "mx=",mx,"p=",p,"t=",t

   call parcel_dilute(ncol, msg, mx, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)


! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
!
   do i = 1,ncol
      plge600(i) = pl(i).ge.600._r8 ! Just change to always allow buoy calculation.
   end do

!
! Main buoyancy calculation.
!
   do i=1,ncol
      do k = nlev,msg + 1,-1
         if (k <= mx(i) .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
            tv(k,i) = t(k,i)* (1._r8+1.608_r8*q(k,i))/ (1._r8+q(k,i))
            buoy(k,i) = tpv(k,i) - tv(k,i) + tiedke_add  ! +0.5K or not?
         else
            qstp(k,i) = q(k,i)
            tp(k,i)   = t(k,i)            
            tpv(k,i)  = tv(k,i)
         endif
      end do
   end do



!-------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
   do i = 1,ncol
      do k = msg + 2,nlev
         if (k < lcl(i) .and. plge600(i)) then
            if (buoy(k+1,i) > 0._r8 .and. buoy(k,i) <= 0._r8) then
               knt(i) = min(5,knt(i) + 1)
               lelten(i,knt(i)) = k
            end if
         end if
      end do
   end do
!
! calculate convective available potential energy (cape).
!
   do n = 1,5
      do i = 1,ncol
         do k = msg + 1,nlev
            if (plge600(i) .and. k <= mx(i) .and. k > lelten(i,n)) then
               capeten(i,n) = capeten(i,n) + rgas*buoy(k,i)*log(pf(k+1,i)/pf(k,i))
               pcapeten(i,n)= pcapeten(i,n)+ buoy(k,i)/tv(k,i)*(pf(k+1,i)-pf(k,i))*100._r8
            end if
         end do
      end do
   end do
!
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
!
   do n = 1,5
      do i = 1,ncol
         if (capeten(i,n) > cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
         end if
         if (pcapeten(i,n)>pcape(i)) then
           pcape(i)=pcapeten(i,n)
         endif
      end do
   end do
!
! put lower bound on cape for diagnostic purposes.
!
   do i = 1,ncol
      cape(i) = max(cape(i), 0._r8)
      pcape(i)= max(pcape(i),0._r8)
   end do
!
   return
end subroutine double_plume_cape_dilute


subroutine parcel_dilute (ncol, msg, klaunch, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)

! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

!--------------------
implicit none
!--------------------

integer, intent(in) :: ncol
integer, intent(in) :: msg

integer, intent(in), dimension(ncol) :: klaunch

real(r8), intent(in), dimension(nlev, ncol) :: p
real(r8), intent(in), dimension(nlev, ncol) :: t
real(r8), intent(in), dimension(nlev, ncol) :: q
real(r8), intent(in), dimension(ncol) :: tpert ! PBL temperature perturbation.

real(r8), intent(inout), dimension(nlev, ncol) :: tp    ! Parcel temp.
real(r8), intent(inout), dimension(nlev, ncol) :: qstp  ! Parcel water vapour (sat value above lcl).
real(r8), intent(inout), dimension(ncol) :: tl         ! Actual temp of LCL.
real(r8), intent(inout), dimension(ncol) :: pl          ! Actual pressure of LCL. 

integer, intent(inout), dimension(ncol) :: lcl ! Lifting condesation level (first model level with saturation).

real(r8), intent(out), dimension(nlev, ncol) :: tpv   ! Define tpv within this routine.

!--------------------

! Have to be careful as s is also dry static energy.


! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.


real(r8) tmix(nlev, ncol)        ! Tempertaure of the entraining parcel.
real(r8) qtmix(nlev, ncol)       ! Total water of the entraining parcel.
real(r8) qsmix(nlev, ncol)       ! Saturated mixing ratio at the tmix.
real(r8) smix(nlev, ncol)        ! Entropy of the entraining parcel.
real(r8) xsh2o(nlev, ncol)       ! Precipitate lost from parcel.
real(r8) ds_xsh2o(nlev, ncol)    ! Entropy change due to loss of condensate.
real(r8) ds_freeze(nlev, ncol)   ! Entropy change sue to freezing of precip.

real(r8) mp(ncol)    ! Parcel mass flux.
real(r8) qtp(ncol)   ! Parcel total water.
real(r8) sp(ncol)    ! Parcel entropy.

real(r8) sp0(ncol)    ! Parcel launch entropy.
real(r8) qtp0(ncol)   ! Parcel launch total water.
real(r8) mp0(ncol)    ! Parcel launch relative mass flux.

real(r8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
real(r8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
!real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
real(r8) dmpdz      ! Parcel fractional mass entrainment rate (/m)
real(r8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
real(r8) senv       ! Environmental entropy at each grid point.
real(r8) qtenv      ! Environmental total water "   "   ".
real(r8) penv       ! Environmental total pressure "   "   ".
real(r8) tenv       ! Environmental total temperature "   "   ".
real(r8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
real(r8) new_q      ! Hold value for total water after condensation/freezing adjustments.
real(r8) dp         ! Layer thickness (center to center)
real(r8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
real(r8) tscool     ! Super cooled temperature offset (in degC) (eg -35).

real(r8) qxsk, qxskp1        ! LCL excess water (k, k+1)
real(r8) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
real(r8) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.

integer rcall       ! Number of ientropy call for errors recording
integer nit_lheat     ! Number of iterations for condensation/freezing loop.
integer i,k,ii   ! Loop counters.

!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!

nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.
dmpdz=-0.5e-3_r8        ! Entrainment rate. (-ve for /m)
!dmpdpc = 3.e-2_r8   ! In cloud entrainment rate (/mb).
lwmax = 1.e-3_r8    ! Need to put formula in for this.
tscool = 0.0_r8   ! Temp at which water loading freezes in the cloud.

qtmix=0._r8
smix=0._r8

qtenv = 0._r8
senv = 0._r8
tenv = 0._r8
penv = 0._r8

qtp0 = 0._r8
sp0  = 0._r8
mp0 = 0._r8

qtp = 0._r8
sp = 0._r8
mp = 0._r8

new_q = 0._r8
new_s = 0._r8

! **** Begin loops ****

do i=1,ncol 
   do k = nlev, msg+1, -1

! Initialize parcel values at launch level.

      if (k == klaunch(i)) then 
         qtp0(i) = q(k,i)   ! Parcel launch total water (assuming subsaturated) - OK????.
         sp0(i)  = entropy(t(k,i),p(k,i),qtp0(i))  ! Parcel launch entropy.
         mp0(i)  = 1._r8       ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute). 
         smix(k,i)  = sp0(i)
         qtmix(k,i) = qtp0(i)
         tfguess = t(k,i)
         rcall = 1
         !write(iulog,*) "smix=",smix,"p=",p(k,i),"qtmix=",qtmix(k,i)
         call ientropy (rcall,i,smix(k,i),p(k,i),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)
         !write(iulog,*) "tmix=",tmix(k,i)
      end if

! Entraining levels
      
      if (k < klaunch(i)) then 

! Set environmental values for this level.                 
         
         dp = (p(k,i)-p(k+1,i)) ! In -ve mb as p decreasing with height - difference between center of layers.
         qtenv = 0.5_r8*(q(k,i)+q(k+1,i))         ! Total water of environment.
         tenv  = 0.5_r8*(t(k,i)+t(k+1,i)) 
         penv  = 0.5_r8*(p(k,i)+p(k+1,i))

         senv  = entropy(tenv,penv,qtenv)  ! Entropy of environment.   

! Determine fractional entrainment rate /pa given value /m.

         dpdz = -(penv*grav)/(rgas*tenv) ! in mb/m since  p in mb.
         dzdp = 1._r8/dpdz                  ! in m/mb
         dmpdp = dmpdz*dzdp              ! /mb Fractional entrainment

! Sum entrainment to current level
! entrains q,s out of intervening dp layers, in which linear variation is assumed
! so really it entrains the mean of the 2 stored values.

         sp(i)  = sp(i)  - dmpdp*dp*senv 
         qtp(i) = qtp(i) - dmpdp*dp*qtenv 
         mp(i)  = mp(i)  - dmpdp*dp
            
! Entrain s and qt to next level.

         smix(k,i)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))
         qtmix(k,i) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))

! Invert entropy from s and q to determine T and saturation-capped q of mixture.
! t(k,i) used as a first guess so that it converges faster.

         tfguess = tmix(k+1,i)
         rcall = 2
        ! write(iulog,*) "tfguess=",tfguess,"p=",p(k,i)
         call ientropy(rcall,i,smix(k,i),p(k,i),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)   
         !write(iulog,*) "tmix=",tmix(k,i)

!
! Determine if this is lcl of this column if qsmix <= qtmix.
! FIRST LEVEL where this happens on ascending.

         if (qsmix(k,i) <= qtmix(k,i) .and. qsmix(k+1,i) > qtmix(k+1,i)) then
            lcl(i) = k
            qxsk   = qtmix(k,i) - qsmix(k,i)
            qxskp1 = qtmix(k+1,i) - qsmix(k+1,i)
            dqxsdp = (qxsk - qxskp1)/dp
            pl(i)  = p(k+1,i) - qxskp1/dqxsdp    ! pressure level of actual lcl.
            dsdp   = (smix(k,i)  - smix(k+1,i))/dp
            dqtdp  = (qtmix(k,i) - qtmix(k+1,i))/dp
            slcl   = smix(k+1,i)  +  dsdp* (pl(i)-p(k+1,i))  
            qtlcl  = qtmix(k+1,i) +  dqtdp*(pl(i)-p(k+1,i))

            tfguess = tmix(k,i)
            rcall = 3
            call ientropy (rcall,i,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess)

!            write(iulog,*)' '
!            write(iulog,*)' p',p(k+1,i),pl(i),p(lcl(i),i)
!            write(iulog,*)' t',tmix(k+1,i),tl(i),tmix(lcl(i),i)
!            write(iulog,*)' s',smix(k+1,i),slcl,smix(lcl(i),i)
!            write(iulog,*)'qt',qtmix(k+1,i),qtlcl,qtmix(lcl(i),i)
!            write(iulog,*)'qs',qsmix(k+1,i),qslcl,qsmix(lcl(i),i)

         endif
!         
      end if !  k < klaunch

 
   end do ! Levels loop
end do ! Columns loop

!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the 
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or 
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
!! provides latent heating to the mixed parcel and so this has to be added back 
!! to it. But does this also increase qsmix as well? Also freezing processes
 

xsh2o = 0._r8
ds_xsh2o = 0._r8
ds_freeze = 0._r8

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy



do i=1,ncol    
   do k = nlev, msg+1, -1
      
! Initialize variables at k=klaunch
      
      if (k == klaunch(i)) then

! Set parcel values at launch level assume no liquid water.            

         tp(k,i)    = tmix(k,i)
         qstp(k,i)  = q(k,i) 
         tpv(k,i)   =  (tp(k,i) + tpert(i)) * (1._r8+1.608_r8*qstp(k,i)) / (1._r8+qstp(k,i))
         
      end if

      if (k < klaunch(i)) then
            
! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.

! Iterate nit_lheat times for s,qt changes.

         do ii=0,nit_lheat-1            

! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).

            xsh2o(k,i) = max (0._r8, qtmix(k,i) - qsmix(k,i) - lwmax)

! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)                     
                     
            ds_xsh2o(k,i) = ds_xsh2o(k+1,i) - cpliq * log (tmix(k,i)/tfreez) * max(0._r8,(xsh2o(k,i)-xsh2o(k+1,i)))
!
! Entropy of freezing: latice times amount of water involved divided by T.
!
 
            if (tmix(k,i) <= tfreez+tscool .and. ds_freeze(k+1,i) == 0._r8) then ! One off freezing of condensate. 
               ds_freeze(k,i) = (latice/tmix(k,i)) * max(0._r8,qtmix(k,i)-qsmix(k,i)-xsh2o(k,i)) ! Gain of LH
            end if
            
            if (tmix(k,i) <= tfreez+tscool .and. ds_freeze(k+1,i) /= 0._r8) then ! Continual freezing of additional condensate.
               ds_freeze(k,i) = ds_freeze(k+1,i)+(latice/tmix(k,i)) * max(0._r8,(qsmix(k+1,i)-qsmix(k,i)))
            end if
            
! Adjust entropy and accordingly to sum of ds (be careful of signs).

            new_s = smix(k,i) + ds_xsh2o(k,i) + ds_freeze(k,i) 

! Adjust liquid water and accordingly to xsh2o.

            new_q = qtmix(k,i) - xsh2o(k,i)

! Invert entropy to get updated Tmix and qsmix of parcel.

            tfguess = tmix(k,i)
            rcall =4
            call ientropy (rcall,i,new_s, p(k,i), new_q, tmix(k,i), qsmix(k,i), tfguess)
            
         end do  ! Iteration loop for freezing processes.

! tp  - Parcel temp is temp of mixture.
! tpv - Parcel v. temp should be density temp with new_q total water. 

         tp(k,i)    = tmix(k,i)

! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)

         if (new_q > qsmix(k,i)) then  ! Super-saturated so condensate present - reduces buoyancy.
            qstp(k,i) = qsmix(k,i)
         else                          ! Just saturated/sub-saturated - no condensate virtual effects.
            qstp(k,i) = new_q
         end if

         tpv(k,i) = (tp(k,i)+tpert(i))* (1._r8+1.608_r8*qstp(k,i)) / (1._r8+ new_q) 

      end if ! k < klaunch
      
   end do  ! Loop for vertical levels.
end do ! Loop for columns
   


return
end subroutine parcel_dilute

!-----------------------------------------------------------------------------------------
real(r8) function entropy(TK,p,qtot)
!-----------------------------------------------------------------------------------------
!
! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
!
     real(r8), intent(in) :: p,qtot,TK
     real(r8) :: qv,qst,e,est,L,eref,pref

pref = 1000.0_r8           ! mb
eref = 6.106_r8            ! sat p at tfreez (mb)

L = rl - (cpliq - cpwv)*(TK-tfreez)         ! T IN CENTIGRADE

! Replace call to satmixutils.

call qmmr_hPa(TK, p, est, qst)

qv = min(qtot,qst)                         ! Partition qtot into vapor part only.
e = qv*p / (eps1 +qv)

entropy = (cpres + qtot*cpliq)*log( TK/tfreez) - rgas*log( (p-e)/pref ) + &
        L*qv/TK - qv*rvap*log(qv/qst)
! 
return
end FUNCTION entropy

!
!-----------------------------------------------------------------------------------------
   SUBROUTINE ientropy (rcall,icol,s,p,qt,T,qst,Tfg)
!-----------------------------------------------------------------------------------------
!
! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
! 

     integer, intent(in) :: icol, rcall
     real(r8), intent(in)  :: s, p, Tfg, qt
     real(r8), intent(out) :: qst, T
     real(r8) :: qv,Ts,dTs,fs1,fs2,est
     real(r8) :: pref,eref,L,e
     real(r8) :: this_lat,this_lon
     integer :: LOOPMAX,i

LOOPMAX = 100                   !* max number of iteration loops 

! Values for entropy
pref = 1000.0_r8           ! mb ref pressure.
eref = 6.106_r8           ! sat p at tfreez (mb)

! Invert the entropy equation -- use Newton's method

Ts = Tfg                  ! Better first guess based on Tprofile from conv.

converge: do i=0, LOOPMAX

   L = rl - (cpliq - cpwv)*(Ts-tfreez) 

   call qmmr_hPa(Ts, p, est, qst)
   qv = min(qt,qst) 
   e = qv*p / (eps1 +qv)  ! Bolton (eq. 16)
   fs1 = (cpres + qt*cpliq)*log( Ts/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/Ts - qv*rvap*log(qv/qst) - s
   
   L = rl - (cpliq - cpwv)*(Ts-1._r8-tfreez)         

   call qmmr_hPa(Ts-1._r8, p, est, qst)
   qv = min(qt,qst) 
   e = qv*p / (eps1 +qv)
   fs2 = (cpres + qt*cpliq)*log( (Ts-1._r8)/tfreez ) - rgas*log( (p-e)/pref ) + &
        L*qv/(Ts-1._r8) - qv*rvap*log(qv/qst) - s 
   
   dTs = fs1/(fs2 - fs1)
   Ts  = Ts+dTs
   if (abs(dTs).lt.0.001_r8) exit converge
   if (i .eq. LOOPMAX - 1) then
    !  this_lat = get_rlat_p(lchnk, icol)*57.296_r8
    !  this_lon = get_rlon_p(lchnk, icol)*57.296_r8
      write(*,*) '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
      write(*,100) 'ZM_CONV: IENTROPY. Details: call#,icol= ',rcall,icol, &
      ! ' lat: ',this_lat,' lon: ',this_lon, &
       ' P(mb)= ', p, ' Tfg(K)= ', Tfg, ' qt(g/kg) = ', 1000._r8*qt, &
       ' qst(g/kg) = ', 1000._r8*qst,', s(J/kg) = ',s
      call endrun('**** ZM_CONV IENTROPY: Tmix did not converge ****')
   end if
enddo converge

! Replace call to satmixutils.

call qmmr_hPa(Ts, p, est, qst)

qv = min(qt,qst)                             !       /* check for saturation */
T = Ts 

 100    format (A,I1,I4,5(A,F6.2))

return
end SUBROUTINE ientropy

! Wrapper for qmmr that does translation between Pa and hPa
! qmmr uses Pa internally, so get qmmr right, need to pass in Pa.
! Afterward, set es back to hPa.
elemental subroutine qmmr_hPa(t, p, es, qm)
  use grist_wv_saturation, only: qmmr

  ! Inputs
  real(r8), intent(in) :: t    ! Temperature (K)
  real(r8), intent(in) :: p    ! Pressure (hPa)
  ! Outputs
  real(r8), intent(out) :: es  ! Saturation vapor pressure (hPa)
  real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                               ! (vapor mass over dry mass, kg/kg)

  call qmmr(t, p*100._r8, es, qm)

  es = es*0.01_r8

end subroutine qmmr_hPa

END Module double_plume_cape_zm
