module grist_modal_aero_wateruptake

    use grist_constants,                    only: r8, pi, rhoh2o
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_nml_module,                   only: nlev
    use grist_rad_constituents,             only: rad_cnst_get_info,        &
                                                  rad_cnst_get_mode_props,  &
                                                  rad_cnst_get_aer_props,   &
                                                  rad_cnst_get_aer_mmr
    use grist_wv_saturation,                only: qsat_water
    use cloud_fraction,                     only: top_lev => clim_modal_aero_top_lev
    use grist_physics_update,               only: old_time_level
    use grist_handle_error,                 only: endrun
    use grist_mpi

implicit none
private
save


public ::   modal_aero_wateruptake_dr

    real(r8), parameter :: third = 1._r8/3._r8
    real(r8), parameter :: pi43  = pi*4.0_r8/3.0_r8


contains

subroutine modal_aero_wateruptake_dr(ncol, list_idx_in, dgnumdry_m, dgnumwet_m, &
                                     qaerwat_m, wetdens_m)
!-----------------------------------------------------------------------
!
! CAM specific driver for modal aerosol water uptake code.
!
! *** N.B. *** The calculation has been enabled for diagnostic mode lists
!              via optional arguments.  If the list_idx arg is present then
!              all the optional args must be present.
!
!-----------------------------------------------------------------------

   ! Arguments
   integer,                     intent(in)         :: ncol
   integer,  optional,          intent(in)         :: list_idx_in
   real(r8), optional,          intent(in)         :: dgnumdry_m(:,:,:)
   real(r8), optional, allocatable, intent(out)    :: dgnumwet_m(:,:,:)
   real(r8), optional, allocatable, intent(out)    :: qaerwat_m(:,:,:)
   real(r8), optional, allocatable, intent(out)    :: wetdens_m(:,:,:)

   ! local variables

   integer  :: list_idx           ! radiative constituents list index

   integer  :: i, k, l, m
   integer  :: nmodes
   integer  :: nspec

   real(r8) :: h2ommr(nlev,ncol) ! specific humidity
   real(r8) :: t(nlev,ncol)      ! temperatures (K)
   real(r8) :: pmid(nlev,ncol)   ! layer pressure (Pa)
   real(r8) :: cldn(nlev,ncol)   ! layer cloud fraction (0-1)
   real(r8) :: raer(nlev,ncol)   ! aerosol species MRs (kg/kg and #/kg)

   real(r8), allocatable :: dgncur_a(:,:,:)
   real(r8), allocatable :: dgncur_awet(:,:,:)
   real(r8), allocatable :: wetdens(:,:,:)
   real(r8), allocatable :: qaerwat(:,:,:)

   real(r8), allocatable :: maer(:,:,:)      ! aerosol wet mass MR (including water) (kg/kg-air)
   real(r8), allocatable :: hygro(:,:,:)     ! volume-weighted mean hygroscopicity (--)
   real(r8), allocatable :: naer(:,:,:)      ! aerosol number MR (bounded!) (#/kg-air)
   real(r8), allocatable :: dryvol(:,:,:)    ! single-particle-mean dry volume (m3)
   real(r8), allocatable :: drymass(:,:,:)   ! single-particle-mean dry mass  (kg)
   real(r8), allocatable :: dryrad(:,:,:)    ! dry volume mean radius of aerosol (m)

   real(r8), allocatable :: wetrad(:,:,:)    ! wet radius of aerosol (m)
   real(r8), allocatable :: wetvol(:,:,:)    ! single-particle-mean wet volume (m3)
   real(r8), allocatable :: wtrvol(:,:,:)    ! single-particle-mean water volume in wet aerosol (m3)

   real(r8), allocatable :: rhcrystal(:)
   real(r8), allocatable :: rhdeliques(:)
   real(r8), allocatable :: specdens_1(:)

   real(r8) :: dryvolmr(nlev,ncol)          ! volume MR for aerosol mode (m3/kg)
   real(r8) :: specdens
   real(r8) :: spechygro, spechygro_1
   real(r8) :: duma, dumb
   real(r8) :: sigmag
   real(r8) :: alnsg
   real(r8) :: v2ncur_a
   real(r8) :: drydens               ! dry particle density  (kg/m^3)
   real(r8) :: rh(nlev,ncol)         ! relative humidity (0-1)

   real(r8) :: es(nlev,ncol)         ! saturation vapor pressure
   real(r8) :: qs(nlev,ncol)         ! saturation specific humidity

   list_idx = 0
   if (present(list_idx_in)) list_idx = list_idx_in

   call rad_cnst_get_info(list_idx, nmodes=nmodes)

   if (list_idx /= 0) then
      ! check that all optional args are present
      if (.not. present(dgnumdry_m)   .or. .not. present(dgnumwet_m) .or. &
          .not. present(qaerwat_m) .or. .not. present(wetdens_m)) then
         call endrun('modal_aero_wateruptake_dr called for'// &
                     'diagnostic list but required args not present')
      end if
      allocate(dgnumwet_m(nmodes, nlev, ncol))
      allocate(qaerwat_m (nmodes, nlev, ncol))
      allocate(wetdens_m (nmodes, nlev, ncol))
   end if

   ! loop over all aerosol modes

   allocate( &
      maer(nmodes,nlev,ncol),     &
      hygro(nmodes,nlev,ncol),    &
      naer(nmodes,nlev,ncol),     &
      dryvol(nmodes,nlev,ncol),   &
      drymass(nmodes,nlev,ncol),  &
      dryrad(nmodes,nlev,ncol),   &
      wetrad(nmodes,nlev,ncol),   &
      wetvol(nmodes,nlev,ncol),   &
      wtrvol(nmodes,nlev,ncol),   &
      rhcrystal(nmodes),           &
      rhdeliques(nmodes),          &
      specdens_1(nmodes)           )

   maer(:,:,:)     = 0._r8
   hygro(:,:,:)    = 0._r8


   allocate(dgncur_a(nmodes,nlev,ncol)) 
   if (list_idx == 0) then
      dgncur_a = pstate_cam%aerosol_dgnum%f(:,:,1:ncol)
   else
      dgncur_a = dgnumdry_m(:,:,1:ncol)
   end if

   do m = 1, nmodes

      dryvolmr(:,:) = 0._r8

      ! get mode properties
      call rad_cnst_get_mode_props(list_idx, m, sigmag=sigmag,  &
         rhcrystal=rhcrystal(m), rhdeliques=rhdeliques(m))

      ! get mode info
      call rad_cnst_get_info(list_idx, m, nspec=nspec)

      do l = 1, nspec

         ! get species interstitial mixing ratio ('a')
         call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', raer)
         call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, hygro_aer=spechygro)

         if (l == 1) then
            ! save off these values to be used as defaults
            specdens_1(m)  = specdens
            spechygro_1    = spechygro
         end if

         do i = 1, ncol
            do k = top_lev, nlev
               duma          = raer(k,i)
               maer(m,k,i)   = maer(m,k,i) + duma
               dumb          = duma/specdens
               dryvolmr(k,i) = dryvolmr(k,i) + dumb
               hygro(m,k,i)  = hygro(m,k,i) + dumb*spechygro
            end do
         end do
      end do

      alnsg = log(sigmag)

      do i = 1, ncol
         do k = top_lev, nlev

            if (dryvolmr(k,i) > 1.0e-30_r8) then
               hygro(m,k,i) = hygro(m,k,i)/dryvolmr(k,i)
            else
               hygro(m,k,i) = spechygro_1
            end if

            ! dry aerosol properties

            v2ncur_a = 1._r8 / ( (pi/6._r8)*(dgncur_a(m,k,i)**3._r8)*exp(4.5_r8*alnsg**2._r8) )
            ! naer = aerosol number (#/kg)
            naer(m,k,i) = dryvolmr(k,i)*v2ncur_a

            ! compute mean (1 particle) dry volume and mass for each mode
            ! old coding is replaced because the new (1/v2ncur_a) is equal to
            ! the mean particle volume
            ! also moletomass forces maer >= 1.0e-30, so (maer/dryvolmr)
            ! should never cause problems (but check for maer < 1.0e-31 anyway)
            if (maer(m,k,i) .gt. 1.0e-31_r8) then
               drydens = maer(m,k,i)/dryvolmr(k,i)
            else
               drydens = 1.0_r8
            end if
            dryvol(m,k,i)   = 1.0_r8/v2ncur_a
            drymass(m,k,i)  = drydens*dryvol(m,k,i)
            dryrad(m,k,i)   = (dryvol(m,k,i)/pi43)**third

         end do
      end do

   end do    ! modes

   ! relative humidity calc

   h2ommr = pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol)
   t      = pstate%temp_at_pc_full_level%f(:,1:ncol)
   pmid   = pstate%pressure_at_pc_full_level%f(:,1:ncol)

   cldn   = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)

   call qsat_water(t(:,1:ncol), pmid(:,1:ncol), es(:,1:ncol), qs(:,1:ncol))
   do i = 1, ncol
      do k = top_lev, nlev
         rh(k,i) = h2ommr(k,i)/qs(k,i)
         rh(k,i) = max(rh(k,i), 0.0_r8)
         rh(k,i) = min(rh(k,i), 0.98_r8)
         if (cldn(k,i) .lt. 1.0_r8) then
            rh(k,i) = (rh(k,i) - cldn(k,i)) / (1.0_r8 - cldn(k,i))  ! clear portion
         end if
         rh(k,i) = max(rh(k,i), 0.0_r8)
      end do
   end do

   call modal_aero_wateruptake_sub( &
      ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
      hygro, rh, dryvol, wetrad, wetvol,           &
      wtrvol)

   allocate(dgncur_awet(nmodes,nlev,ncol), &
            wetdens(nmodes,nlev,ncol),     &
            qaerwat(nmodes,nlev,ncol) )

   do i = 1, ncol
      do k = top_lev, nlev
         do m = 1, nmodes

            dgncur_awet(m,k,i) = dgncur_a(m,k,i) * (wetrad(m,k,i)/dryrad(m,k,i))
            qaerwat(m,k,i)     = rhoh2o*naer(m,k,i)*wtrvol(m,k,i)

            ! compute aerosol wet density (kg/m3)
            if (wetvol(m,k,i) > 1.0e-30_r8) then
               wetdens(m,k,i) = (drymass(m,k,i) + rhoh2o*wtrvol(m,k,i))/wetvol(m,k,i)
            else
               wetdens(m,k,i) = specdens_1(m)
            end if

         end do
      end do
   end do

   if (list_idx == 0) then
      pstate_cam%aerosol_dgnumwet%f(:,:,1:ncol) = dgncur_awet
      pstate_cam%aerosol_qaerwat%f(:,:,1:ncol)  = qaerwat
      pstate_cam%aerosol_wetdens_ap%f(:,:,1:ncol)  = wetdens
   else
      ! for diagnostic calcs just return results
      dgnumwet_m = dgncur_awet
      qaerwat_m  = qaerwat
      wetdens_m  = wetdens
   end if

   deallocate( &
      maer,   hygro,  naer,   dryvol,    drymass,    dryrad, &
      wetrad, wetvol, wtrvol, rhcrystal, rhdeliques, specdens_1)

   deallocate(dgncur_a, dgncur_awet, qaerwat, wetdens)

end subroutine modal_aero_wateruptake_dr



subroutine modal_aero_wateruptake_sub( &
   ncol, nmodes, rhcrystal, rhdeliques, dryrad, &
   hygro, rh, dryvol, wetrad, wetvol,           &
   wtrvol)

!-----------------------------------------------------------------------
!
! Purpose: Compute aerosol wet radius
!
! Method:  Kohler theory
!
! Author:  S. Ghan
!
!-----------------------------------------------------------------------

   ! Arguments
   integer, intent(in)  :: ncol                    ! number of columns
   integer, intent(in)  :: nmodes

   real(r8), intent(in) :: rhcrystal(:)
   real(r8), intent(in) :: rhdeliques(:)
   real(r8), intent(in) :: dryrad(:,:,:)         ! dry volume mean radius of aerosol (m)
   real(r8), intent(in) :: hygro(:,:,:)          ! volume-weighted mean hygroscopicity (--)
   real(r8), intent(in) :: rh(:,:)               ! relative humidity (0-1)
   real(r8), intent(in) :: dryvol(:,:,:)

   real(r8), intent(out) :: wetrad(:,:,:)        ! wet radius of aerosol (m)
   real(r8), intent(out) :: wetvol(:,:,:)        ! single-particle-mean wet volume (m3)
   real(r8), intent(out) :: wtrvol(:,:,:)        ! single-particle-mean water volume in wet aerosol (m3)

   ! local variables

   integer :: i, k, m

   real(r8) :: hystfac                ! working variable for hysteresis
   !-----------------------------------------------------------------------


   ! loop over all aerosol modes
   do m = 1, nmodes

      hystfac = 1.0_r8 / max(1.0e-5_r8, (rhdeliques(m) - rhcrystal(m)))

      do i = 1, ncol
         do k = top_lev, nlev

            ! compute wet radius for each mode
            call modal_aero_kohler(dryrad(m,k,i:i), hygro(m,k,i:i), rh(k,i:i), wetrad(m,k,i:i), 1, 1)

            wetrad(m,k,i) = max(wetrad(m,k,i), dryrad(m,k,i))
            wetvol(m,k,i) = pi43*wetrad(m,k,i)**3
            wetvol(m,k,i) = max(wetvol(m,k,i), dryvol(m,k,i))
            wtrvol(m,k,i) = wetvol(m,k,i) - dryvol(m,k,i)
            wtrvol(m,k,i) = max(wtrvol(m,k,i), 0.0_r8)

            ! apply simple treatment of deliquesence/crystallization hysteresis
            ! for rhcrystal < rh < rhdeliques, aerosol water is a fraction of
            ! the "upper curve" value, and the fraction is a linear function of rh
            if (rh(k,i) < rhcrystal(m)) then
               wetrad(m,k,i) = dryrad(m,k,i)
               wetvol(m,k,i) = dryvol(m,k,i)
               wtrvol(m,k,i) = 0.0_r8
            else if (rh(k,i) < rhdeliques(m)) then
               wtrvol(m,k,i) = wtrvol(m,k,i)*hystfac*(rh(k,i) - rhcrystal(m))
               wtrvol(m,k,i) = max(wtrvol(m,k,i), 0.0_r8)
               wetvol(m,k,i) = dryvol(m,k,i) + wtrvol(m,k,i)
               wetrad(m,k,i) = (wetvol(m,k,i)/pi43)**third
            end if

         end do  ! columns
      end do     ! levels

   end do ! modes

end subroutine modal_aero_wateruptake_sub



subroutine modal_aero_kohler(   &
           rdry_in, hygro, s, rwet_out, im, imx )

! calculates equlibrium radius r of haze droplets as function of
! dry particle mass and relative humidity s using kohler solution
! given in pruppacher and klett (eqn 6-35)

! for multiple aerosol types, assumes an internal mixture of aerosols

      implicit none

! arguments
      integer :: im         ! number of grid points to be processed
      integer :: imx        ! dimensioned number of grid points
      real(r8) :: rdry_in(imx)    ! aerosol dry radius (m)
      real(r8) :: hygro(imx)      ! aerosol volume-mean hygroscopicity (--)
      real(r8) :: s(imx)          ! relative humidity (1 = saturated)
      real(r8) :: rwet_out(imx)   ! aerosol wet radius (m)

! local variables
      integer, parameter :: imax=200
      integer :: i, n, nsol

      real(r8) :: a, b
      real(r8) :: p40(imax),p41(imax),p42(imax),p43(imax) ! coefficients of polynomial
      real(r8) :: p30(imax),p31(imax),p32(imax) ! coefficients of polynomial
      real(r8) :: p
      real(r8) :: r3, r4
      real(r8) :: r(imx)        ! wet radius (microns)
      real(r8) :: rdry(imax)    ! radius of dry particle (microns)
      real(r8) :: ss            ! relative humidity (1 = saturated)
      real(r8) :: slog(imax)    ! log relative humidity
      real(r8) :: vol(imax)     ! total volume of particle (microns**3)
      real(r8) :: xi, xr

      complex(r8) :: cx4(4,imax),cx3(3,imax)

      real(r8), parameter :: eps = 1.e-4_r8
      real(r8), parameter :: mw = 18._r8
      !real(r8), parameter :: pi = 3.14159_r8
      real(r8), parameter :: rhow = 1._r8
      real(r8), parameter :: surften = 76._r8
      real(r8), parameter :: tair = 273._r8
      real(r8), parameter :: third = 1._r8/3._r8
      real(r8), parameter :: ugascon = 8.3e7_r8


!     effect of organics on surface tension is neglected
      a=2.e4_r8*mw*surften/(ugascon*tair*rhow)

      do i=1,im
           rdry(i) = rdry_in(i)*1.0e6_r8   ! convert (m) to (microns)
           vol(i) = rdry(i)**3          ! vol is r**3, not volume
           b = vol(i)*hygro(i)

!          quartic
           ss=min(s(i),1._r8-eps)
           ss=max(ss,1.e-10_r8)
           slog(i)=log(ss)
           p43(i)=-a/slog(i)
           p42(i)=0._r8
           p41(i)=b/slog(i)-vol(i)
           p40(i)=a*vol(i)/slog(i)
!          cubic for rh=1
           p32(i)=0._r8
           p31(i)=-b/a
           p30(i)=-vol(i)
      end do


       do 100 i=1,im

!       if(vol(i).le.1.e-20)then
        if(vol(i).le.1.e-12_r8)then
           r(i)=rdry(i)
           go to 100
        endif

        p=abs(p31(i))/(rdry(i)*rdry(i))
        if(p.lt.eps)then
!          approximate solution for small particles
           r(i)=rdry(i)*(1._r8+p*third/(1._r8-slog(i)*rdry(i)/a))
        else
           call makoh_quartic(cx4(1,i),p43(i),p42(i),p41(i),p40(i),1)
!          find smallest real(r8) solution
           r(i)=1000._r8*rdry(i)
           nsol=0
           do n=1,4
              xr=real(cx4(n,i))
              xi=aimag(cx4(n,i))
              if(abs(xi).gt.abs(xr)*eps) cycle  
              if(xr.gt.r(i)) cycle  
              if(xr.lt.rdry(i)*(1._r8-eps)) cycle  
              if(xr.ne.xr) cycle  
              r(i)=xr
              nsol=n
           end do  
           if(nsol.eq.0)then
              print*,'rank =',mpi_rank(),'ccm kohlerc - no real(r8) solution found (quartic)'
              print*,'roots =', (cx4(n,i),n=1,4)
              print*,'p0-p3 =', p40(i), p41(i), p42(i), p43(i)
              print*,'rh=',s(i)
              print*,'setting radius to dry radius=',rdry(i)
              r(i)=rdry(i)
!             stop
           endif
        endif

        if(s(i).gt.1._r8-eps)then
!          save quartic solution at s=1-eps
           r4=r(i)
!          cubic for rh=1
           p=abs(p31(i))/(rdry(i)*rdry(i))
           if(p.lt.eps)then
              r(i)=rdry(i)*(1._r8+p*third)
           else
              call makoh_cubic(cx3,p32,p31,p30,im)
!             find smallest real(r8) solution
              r(i)=1000._r8*rdry(i)
              nsol=0
              do n=1,3
                 xr=real(cx3(n,i))
                 xi=aimag(cx3(n,i))
                 if(abs(xi).gt.abs(xr)*eps) cycle  
                 if(xr.gt.r(i)) cycle  
                 if(xr.lt.rdry(i)*(1._r8-eps)) cycle  
                 if(xr.ne.xr) cycle  
                 r(i)=xr
                 nsol=n
              end do  
              if(nsol.eq.0)then
                  print*,'rank =',mpi_rank(),'ccm kohlerc - no real(r8) solution found (cubic)'
                 print*,'roots =', (cx3(n,i),n=1,3)
                 print*,'p0-p2 =', p30(i), p31(i), p32(i)
                 print*,'rh=',s(i)
                 print*,'setting radius to dry radius=',rdry(i)
                 r(i)=rdry(i)
!                stop
              endif
           endif
           r3=r(i)
!          now interpolate between quartic, cubic solutions
           r(i)=(r4*(1._r8-s(i))+r3*(s(i)-1._r8+eps))/eps
        endif

  100 continue

! bound and convert from microns to m
      do i=1,im
         r(i) = min(r(i),30._r8) ! upper bound based on 1 day lifetime
         rwet_out(i) = r(i)*1.e-6_r8
      end do

      return
end subroutine modal_aero_kohler

subroutine makoh_cubic( cx, p2, p1, p0, im )
!
!     solves  x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx)
      complex(r8) :: cx(3,imx)

      integer :: i
      real(r8) :: eps, q(imx), r(imx), sqrt3, third
      complex(r8) :: ci, cq, crad(imx), cw, cwsq, cy(imx), cz(imx)

      save eps
      data eps/1.e-20_r8/

      third=1._r8/3._r8
      ci=cmplx(0._r8,1._r8,r8)
      sqrt3=sqrt(3._r8)
      cw=0.5_r8*(-1+ci*sqrt3)
      cwsq=0.5_r8*(-1-ci*sqrt3)

      do i=1,im
      if(p1(i).eq.0._r8)then
!        completely insoluble particle
         cx(1,i)=(-p0(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
      else
         q(i)=p1(i)/3._r8
         r(i)=p0(i)/2._r8
         crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
         crad(i)=sqrt(crad(i))

         cy(i)=r(i)-crad(i)
         if (abs(cy(i)).gt.eps) cy(i)=cy(i)**third
         cq=q(i)
         cz(i)=-cq/cy(i)

         cx(1,i)=-cy(i)-cz(i)
         cx(2,i)=-cw*cy(i)-cwsq*cz(i)
         cx(3,i)=-cwsq*cy(i)-cw*cz(i)
      endif
      enddo

      return
end subroutine makoh_cubic

subroutine makoh_quartic( cx, p3, p2, p1, p0, im )

!     solves x**4 + p3 x**3 + p2 x**2 + p1 x + p0 = 0
!     where p0, p1, p2, p3 are real
!
      integer, parameter :: imx=200
      integer :: im
      real(r8) :: p0(imx), p1(imx), p2(imx), p3(imx)
      complex(r8) :: cx(4,imx)

      integer :: i
      real(r8) :: third, q(imx), r(imx)
      complex(r8) :: cb(imx), cb0(imx), cb1(imx),   &
                     crad(imx), cy(imx), czero


      czero=cmplx(0.0_r8,0.0_r8)
      third=1._r8/3._r8

      do 10 i=1,im

      q(i)=-p2(i)*p2(i)/36._r8+(p3(i)*p1(i)-4*p0(i))/12._r8
      r(i)=-(p2(i)/6)**3+p2(i)*(p3(i)*p1(i)-4*p0(i))/48._r8   &
       +(4*p0(i)*p2(i)-p0(i)*p3(i)*p3(i)-p1(i)*p1(i))/16

      crad(i)=r(i)*r(i)+q(i)*q(i)*q(i)
      crad(i)=sqrt(crad(i))

      cb(i)=r(i)-crad(i)
      if(cb(i).eq.czero)then
!        insoluble particle
         cx(1,i)=(-p1(i))**third
         cx(2,i)=cx(1,i)
         cx(3,i)=cx(1,i)
         cx(4,i)=cx(1,i)
      else
         cb(i)=cb(i)**third

         cy(i)=-cb(i)+q(i)/cb(i)+p2(i)/6

         cb0(i)=sqrt(cy(i)*cy(i)-p0(i))
         cb1(i)=(p3(i)*cy(i)-p1(i))/(2*cb0(i))

         cb(i)=p3(i)/2+cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)+cb0(i))
         crad(i)=sqrt(crad(i))
         cx(1,i)=(-cb(i)+crad(i))/2._r8
         cx(2,i)=(-cb(i)-crad(i))/2._r8

         cb(i)=p3(i)/2-cb1(i)
         crad(i)=cb(i)*cb(i)-4*(cy(i)-cb0(i))
         crad(i)=sqrt(crad(i))
         cx(3,i)=(-cb(i)+crad(i))/2._r8
         cx(4,i)=(-cb(i)-crad(i))/2._r8
      endif
   10 continue

      return
end subroutine makoh_quartic


end module grist_modal_aero_wateruptake
