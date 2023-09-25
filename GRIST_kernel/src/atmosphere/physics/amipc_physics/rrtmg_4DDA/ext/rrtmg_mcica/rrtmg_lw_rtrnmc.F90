!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrnmc.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.3 $
!     created:   $Date: 2008/04/24 16:17:28 $
!
      module rrtmg_lw_rtrnmc

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! --------- Modules ----------
      use grist_constants,            only: r8

!      use parkind, only : jpim, jprb 
      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_con, only: fluxfac, heatfac
      use rrlw_wvn, only: delwave, ngb, ngs
      use rrlw_tbl, only: tblint, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrrtc, hnamrtc
      use rrlw_wvn, only: wavenum1, wavenum2
      !use grist_mpi

      implicit none

      contains

!-----------------------------------------------------------------------------
      subroutine rtrnmc(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                        cldfmc, taucmc, planklay, planklev, plankbnd, &
                        pwvcm, fracs, taut, &
                        totuflux, totdflux, fnet, htr, &
                        totuclfl, totdclfl, fnetc, htrc, totufluxs, totdfluxs ) 
!-----------------------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with the McICA stochastic approach and maximum-random
!  cloud overlap. 
!***************************************************************************

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: nlayers                    ! total number of layers
      integer, intent(in) :: istart                     ! beginning band of calculation
      integer, intent(in) :: iend                       ! ending band of calculation
      integer, intent(in) :: iout                       ! output option flag

! Atmosphere
      real(kind=r8), intent(in) :: pz(0:)               ! level (interface) pressures (hPa, mb)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: pwvcm                ! precipitable water vapor (cm)
      real(kind=r8), intent(in) :: semiss(:)            ! lw surface emissivity
                                                        !    Dimensions: (nbndlw)
      real(kind=r8), intent(in) :: planklay(:,:)        ! 
                                                        !    Dimensions: (nbndlw,nlayers)
      real(kind=r8), intent(in) :: planklev(:,0:)       ! 
                                                        !    Dimensions: (nbndlw,0:nlayers)
      real(kind=r8), intent(in) :: plankbnd(:)          ! 
                                                        !    Dimensions: (nbndlw)
      real(kind=r8), intent(in) :: fracs(:,:)           ! 
                                                        !    Dimensions: (ngptw,nlayers)
      real(kind=r8), intent(in) :: taut(:,:)            ! gaseous + aerosol optical depths
                                                        !    Dimensions: (ngptlw,nlayers)

! Clouds
      integer, intent(in) :: ncbands                    ! number of cloud spectral bands
      real(kind=r8), intent(in) :: cldfmc(:,:)          ! layer cloud fraction [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: taucmc(:,:)          ! layer cloud optical depth [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)

! ----- Output -----
      real(kind=r8), intent(out) :: totuflux(0:)        ! upward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdflux(0:)        ! downward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnet(0:)            ! net longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htr(0:)             ! longwave heating rate (k/day)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totuclfl(0:)        ! clear sky upward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdclfl(0:)        ! clear sky downward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnetc(0:)           ! clear sky net longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htrc(0:)            ! clear sky longwave heating rate (k/day)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totufluxs(:,0:)     ! upward longwave flux spectral (w/m2) 
                                                        !    Dimensions: (nbndlw, 0:nlayers)
      real(kind=r8), intent(out) :: totdfluxs(:,0:)     ! downward longwave flux spectral (w/m2)
                                                        !    Dimensions: (nbndlw, 0:nlayers)

! ----- Local -----
! Declarations for radiative transfer
      real(kind=r8) :: abscld(ngptlw,nlayers)
      real(kind=r8) :: atot(nlayers)
      real(kind=r8) :: atrans(nlayers)
      real(kind=r8) :: bbugas(nlayers)
      real(kind=r8) :: bbutot(nlayers)
      real(kind=r8) :: clrurad(0:nlayers)
      real(kind=r8) :: clrdrad(0:nlayers)
      real(kind=r8) :: efclfrac(ngptlw,nlayers)
      real(kind=r8) :: uflux(0:nlayers)
      real(kind=r8) :: dflux(0:nlayers)
      real(kind=r8) :: urad(0:nlayers)
      real(kind=r8) :: drad(0:nlayers)
      real(kind=r8) :: uclfl(0:nlayers)
      real(kind=r8) :: dclfl(0:nlayers)
      real(kind=r8) :: odcld(ngptlw,nlayers)


      real(kind=r8) :: secdiff(nbndlw)                   ! secant of diffusivity angle
      real(kind=r8) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)  ! diffusivity angle adjustment coefficients
      real(kind=r8) :: wtdiff, rec_6
      real(kind=r8) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=r8) :: odepth, odtot, odepth_rec, odtot_rec, gassrc
      real(kind=r8) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real(kind=r8) :: rad0, reflect, radlu, radclru

      integer :: icldlyr(nlayers)                        ! flag for cloud in layer
      integer :: ibnd, ib, iband, lay, lev, l, ig        ! loop indices
      integer :: igc                                     ! g-point interval counter
      integer :: iclddn                                  ! flag for cloud in down path
      integer :: ittot, itgas, itr                       ! lookup table indices


! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdgas                       ! gas-only planck function for downward rt
!    bbugas                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance 
!    radclru                      ! spectrally summed clear sky upward radiance 
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance 
!    radclrd                      ! spectrally summed clear sky downward radiance 
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)


! This secant and weight corresponds to the standard diffusivity 
! angle.  This initial value is redefined below for some bands.
      data wtdiff /0.5_r8/
      data rec_6 /0.166667_r8/

! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
! and 1.80) as a function of total column water vapor.  The function
! has been defined to minimize flux and cooling rate errors in these bands
! over a wide range of precipitable water values.
      data a0 / 1.66_r8,  1.55_r8,  1.58_r8,  1.66_r8, &
                1.54_r8, 1.454_r8,  1.89_r8,  1.33_r8, &
               1.668_r8,  1.66_r8,  1.66_r8,  1.66_r8, &
                1.66_r8,  1.66_r8,  1.66_r8,  1.66_r8 /
      data a1 / 0.00_r8,  0.25_r8,  0.22_r8,  0.00_r8, &
                0.13_r8, 0.446_r8, -0.10_r8,  0.40_r8, &
              -0.006_r8,  0.00_r8,  0.00_r8,  0.00_r8, &
                0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8 /
      data a2 / 0.00_r8, -12.0_r8, -11.7_r8,  0.00_r8, &
               -0.72_r8,-0.243_r8,  0.19_r8,-0.062_r8, &
               0.414_r8,  0.00_r8,  0.00_r8,  0.00_r8, &
                0.00_r8,  0.00_r8,  0.00_r8,  0.00_r8 /

      hvrrtc = '$Revision: 1.3 $'

      do ibnd = 1,nbndlw
         if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
           secdiff(ibnd) = 1.66_r8
         else
           secdiff(ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm)
         endif
      enddo
      if (pwvcm.lt.1.0) secdiff(6) = 1.80_r8
      if (pwvcm.gt.7.1) secdiff(7) = 1.50_r8

      urad(0) = 0.0_r8
      drad(0) = 0.0_r8
      totuflux(0) = 0.0_r8
      totdflux(0) = 0.0_r8
      clrurad(0) = 0.0_r8
      clrdrad(0) = 0.0_r8
      totuclfl(0) = 0.0_r8
      totdclfl(0) = 0.0_r8

      do lay = 1, nlayers
         urad(lay) = 0.0_r8
         drad(lay) = 0.0_r8
         totuflux(lay) = 0.0_r8
         totdflux(lay) = 0.0_r8
         clrurad(lay) = 0.0_r8
         clrdrad(lay) = 0.0_r8
         totuclfl(lay) = 0.0_r8
         totdclfl(lay) = 0.0_r8
         icldlyr(lay) = 0

! Change to band loop?
         do ig = 1, ngptlw
            if (cldfmc(ig,lay) .eq. 1._r8) then
               ib = ngb(ig)
               odcld(ig,lay) = secdiff(ib) * taucmc(ig,lay)
               transcld = exp(-odcld(ig,lay))
               abscld(ig,lay) = 1._r8 - transcld
               efclfrac(ig,lay) = abscld(ig,lay) * cldfmc(ig,lay)
               icldlyr(lay) = 1
            else
               odcld(ig,lay) = 0.0_r8
               abscld(ig,lay) = 0.0_r8
               efclfrac(ig,lay) = 0.0_r8
            endif
         enddo

      enddo

      igc = 1
! Loop over frequency bands.
      do iband = istart, iend

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1

! Loop over g-channels.
 1000    continue

! Radiative transfer starts here.
         radld = 0._r8
         radclrd = 0._r8
         iclddn = 0

! Downward radiative transfer loop.  

         do lev = nlayers, 1, -1
               plfrac = fracs(igc,lev)
               blay = planklay(iband,lev)
               dplankup = planklev(iband,lev) - blay
               dplankdn = planklev(iband,lev-1) - blay
               odepth = secdiff(iband) * taut(igc,lev)
               if (odepth .lt. 0.0_r8) odepth = 0.0_r8
!  Cloudy layer
               if (icldlyr(lev).eq.1) then
                  iclddn = 1
                  odtot = odepth + odcld(igc,lev)
                  if (odtot .lt. 0.06_r8) then
                     atrans(lev) = odepth - 0.5_r8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     atot(lev) =  odtot - 0.5_r8*odtot*odtot
                     odtot_rec = rec_6*odtot
                     bbdtot =  plfrac * (blay+dplankdn*odtot_rec)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(igc,lev) * (1. - atrans(lev))) + &
                         gassrc + cldfmc(igc,lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld
                  
                     bbugas(lev) =  plfrac * (blay+dplankup*odepth_rec)
                     bbutot(lev) =  plfrac * (blay+dplankup*odtot_rec)

                  elseif (odepth .le. 0.06_r8) then
                     atrans(lev) = odepth - 0.5_r8*odepth*odepth
                     odepth_rec = rec_6*odepth
                     gassrc = plfrac*(blay+dplankdn*odepth_rec)*atrans(lev)

                     odtot = odepth + odcld(igc,lev)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_r8
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+dplankdn*odepth_rec)
                     atot(lev) = 1. - exp_tbl(ittot)

                     radld = radld - radld * (atrans(lev) + &
                         efclfrac(igc,lev) * (1._r8 - atrans(lev))) + &
                         gassrc + cldfmc(igc,lev) * &
                         (bbdtot * atot(lev) - gassrc)
                     drad(lev-1) = drad(lev-1) + radld

                     bbugas(lev) = plfrac * (blay + dplankup*odepth_rec)
                     bbutot(lev) = plfrac * (blay + tfactot * dplankup)

                  else

                     tblind = odepth/(bpade+odepth)
                     itgas = tblint*tblind+0.5_r8
                     odepth = tau_tbl(itgas)
                     atrans(lev) = 1._r8 - exp_tbl(itgas)
                     tfacgas = tfn_tbl(itgas)
                     gassrc = atrans(lev) * plfrac * (blay + tfacgas*dplankdn)

                     odtot = odepth + odcld(igc,lev)
                     tblind = odtot/(bpade+odtot)
                     ittot = tblint*tblind + 0.5_r8
                     tfactot = tfn_tbl(ittot)
                     bbdtot = plfrac * (blay + tfactot*dplankdn)
                     bbd = plfrac*(blay+tfacgas*dplankdn)
                     atot(lev) = 1._r8 - exp_tbl(ittot)

                  radld = radld - radld * (atrans(lev) + &
                    efclfrac(igc,lev) * (1._r8 - atrans(lev))) + &
                    gassrc + cldfmc(igc,lev) * &
                    (bbdtot * atot(lev) - gassrc)
                  drad(lev-1) = drad(lev-1) + radld
                  bbugas(lev) = plfrac * (blay + tfacgas * dplankup)
                  bbutot(lev) = plfrac * (blay + tfactot * dplankup)
                  endif
!  Clear layer
               else
                  if (odepth .le. 0.06_r8) then
                     atrans(lev) = odepth-0.5_r8*odepth*odepth
                     odepth = rec_6*odepth
                     bbd = plfrac*(blay+dplankdn*odepth)
                     bbugas(lev) = plfrac*(blay+dplankup*odepth)
                  else
                     tblind = odepth/(bpade+odepth)
                     itr = tblint*tblind+0.5_r8
                     transc = exp_tbl(itr)
                     atrans(lev) = 1._r8-transc
                     tausfac = tfn_tbl(itr)
                     bbd = plfrac*(blay+tausfac*dplankdn)
                     bbugas(lev) = plfrac * (blay + tausfac * dplankup)
                  endif   
                  radld = radld + (bbd-radld)*atrans(lev)
                  drad(lev-1) = drad(lev-1) + radld
               endif
!  Set clear sky stream to total sky stream as long as layers
!  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
!  and clear sky stream must be computed separately from that point.
                  if (iclddn.eq.1) then
                     radclrd = radclrd + (bbd-radclrd) * atrans(lev) 
                     clrdrad(lev-1) = clrdrad(lev-1) + radclrd
                  else
                     radclrd = radld
                     clrdrad(lev-1) = drad(lev-1)
                  endif
            enddo

! Spectral emissivity & reflectance
!  Include the contribution of spectrally varying longwave emissivity
!  and reflection from the surface to the upward radiative transfer.
!  Note: Spectral and Lambertian reflection are identical for the
!  diffusivity angle flux integration used here.

         rad0 = fracs(igc,1) * plankbnd(iband)
!  Add in specular reflection of surface downward radiance.
         reflect = 1._r8 - semiss(iband)
         radlu = rad0 + reflect * radld
         radclru = rad0 + reflect * radclrd


! Upward radiative transfer loop.
         urad(0) = urad(0) + radlu
         clrurad(0) = clrurad(0) + radclru

         do lev = 1, nlayers
!  Cloudy layer
            if (icldlyr(lev) .eq. 1) then
               gassrc = bbugas(lev) * atrans(lev)
               radlu = radlu - radlu * (atrans(lev) + &
                   efclfrac(igc,lev) * (1._r8 - atrans(lev))) + &
                   gassrc + cldfmc(igc,lev) * &
                   (bbutot(lev) * atot(lev) - gassrc)
               urad(lev) = urad(lev) + radlu
!  Clear layer
            else
               radlu = radlu + (bbugas(lev)-radlu)*atrans(lev)
               urad(lev) = urad(lev) + radlu
            endif
!  Set clear sky stream to total sky stream as long as all layers
!  are clear (iclddn=0).  Streams must be calculated separately at 
!  all layers when a cloud is present (ICLDDN=1), because surface 
!  reflectance is different for each stream.
               if (iclddn.eq.1) then
                  radclru = radclru + (bbugas(lev)-radclru)*atrans(lev) 
                  clrurad(lev) = clrurad(lev) + radclru
               else
                  radclru = radlu
                  clrurad(lev) = urad(lev)
               endif
         enddo

! Increment g-point counter
         igc = igc + 1
! Return to continue radiative transfer for all g-channels in present band
         if (igc .le. ngs(iband)) go to 1000

! Process longwave output from band for total and clear streams.
! Calculate upward, downward, and net flux.
         do lev = nlayers, 0, -1
            uflux(lev) = urad(lev)*wtdiff
            dflux(lev) = drad(lev)*wtdiff
            urad(lev) = 0.0_r8
            drad(lev) = 0.0_r8
            totuflux(lev) = totuflux(lev) + uflux(lev) * delwave(iband)
            totdflux(lev) = totdflux(lev) + dflux(lev) * delwave(iband)
            uclfl(lev) = clrurad(lev)*wtdiff
            dclfl(lev) = clrdrad(lev)*wtdiff
            clrurad(lev) = 0.0_r8
            clrdrad(lev) = 0.0_r8
            totuclfl(lev) = totuclfl(lev) + uclfl(lev) * delwave(iband)
            totdclfl(lev) = totdclfl(lev) + dclfl(lev) * delwave(iband)
            totufluxs(iband,lev) = uflux(lev) * delwave(iband)
            totdfluxs(iband,lev) = dflux(lev) * delwave(iband)
         enddo

! End spectral band loop
      enddo

! Calculate fluxes at surface
      totuflux(0) = totuflux(0) * fluxfac
      totdflux(0) = totdflux(0) * fluxfac
      totufluxs(:,0) = totufluxs(:,0) * fluxfac
      totdfluxs(:,0) = totdfluxs(:,0) * fluxfac
      fnet(0) = totuflux(0) - totdflux(0)
      totuclfl(0) = totuclfl(0) * fluxfac
      totdclfl(0) = totdclfl(0) * fluxfac
      fnetc(0) = totuclfl(0) - totdclfl(0)

! Calculate fluxes at model levels
      do lev = 1, nlayers
         totuflux(lev) = totuflux(lev) * fluxfac
         totdflux(lev) = totdflux(lev) * fluxfac
         totufluxs(:,lev) = totufluxs(:,lev) * fluxfac
         totdfluxs(:,lev) = totdfluxs(:,lev) * fluxfac
         fnet(lev) = totuflux(lev) - totdflux(lev)
         totuclfl(lev) = totuclfl(lev) * fluxfac
         totdclfl(lev) = totdclfl(lev) * fluxfac
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
         l = lev - 1

! Calculate heating rates at model layers
         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev)) 
      enddo

! Set heating rate to zero in top layer
      htr(nlayers) = 0.0_r8
      htrc(nlayers) = 0.0_r8

      end subroutine rtrnmc

!--------------------------linhan--------------------------------
      subroutine rtrnmc_4DDA(nlayers, istart, iend, iout, pz, semiss, ncbands, &
                        cldfmc, taucmc, planklay, planklev, plankbnd, &
                        pwvcm, fracs, taut, &
                        ssacmc, asmcmc, tz, amu0, & !linhan
                        totuflux, totdflux, fnet, htr, &
                        totuclfl, totdclfl, fnetc, htrc, totufluxs, totdfluxs ) 
!-----------------------------------------------------------------------------
!
!  Original version:   E. J. Mlawer, et al. RRTM_V3.0
!  Revision for GCMs:  Michael J. Iacono; October, 2002
!  Revision for F90:  Michael J. Iacono; June, 2006
!
!  This program calculates the upward fluxes, downward fluxes, and
!  heating rates for an arbitrary clear or cloudy atmosphere.  The input
!  to this program is the atmospheric profile, all Planck function
!  information, and the cloud fraction by layer.  A variable diffusivity 
!  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
!  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
!  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
!  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
!  use of the emissivity angle for the flux integration can cause errors of 
!  1 to 4 W/m2 within cloudy layers.  
!  Clouds are treated with the McICA stochastic approach and maximum-random
!  cloud overlap. 
!***************************************************************************

! ------- Declarations -------

! ----- Input -----
      integer, intent(in) :: nlayers                    ! total number of layers
      integer, intent(in) :: istart                     ! beginning band of calculation
      integer, intent(in) :: iend                       ! ending band of calculation
      integer, intent(in) :: iout                       ! output option flag

! Atmosphere
      real(kind=r8), intent(in) :: pz(0:)               ! level (interface) pressures (hPa, mb)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: tz(0:)             ! level (interface) temperatures (K) linhan
                                                      !    Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: pwvcm                ! precipitable water vapor (cm)
      real(kind=r8), intent(in) :: semiss(:)            ! lw surface emissivity
                                                        !    Dimensions: (nbndlw)
      real(kind=r8), intent(in) :: planklay(:,:)        ! 
                                                        !    Dimensions: (nbndlw,nlayers)
      real(kind=r8), intent(in) :: planklev(:,0:)       ! 
                                                        !    Dimensions: (nbndlw,0:nlayers)
      real(kind=r8), intent(in) :: plankbnd(:)          ! 
                                                        !    Dimensions: (nbndlw)
      real(kind=r8), intent(in) :: fracs(:,:)           ! 
                                                        !    Dimensions: (ngptw,nlayers)
      real(kind=r8), intent(in) :: taut(:,:)            ! gaseous + aerosol optical depths
                                                        !    Dimensions: (ngptlw,nlayers)

! Clouds
      integer, intent(in) :: ncbands                    ! number of cloud spectral bands
      real(kind=r8), intent(in) :: cldfmc(:,:)          ! layer cloud fraction [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: taucmc(:,:)          ! layer cloud optical depth [mcica]
                                                        !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: ssacmc(:,:)        ! cloud single scattering albedo [mcica] linhan
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=r8), intent(in) :: asmcmc(:,:)        ! cloud asymmetry parameter [mcica] linhan
                                                      !    Dimensions: (ngptlw,nlayers) 
      real(kind=r8), intent(in) :: amu0        ! cosine of solar zenith angle linhan

! ----- Output -----
      real(kind=r8), intent(out) :: totuflux(0:)        ! upward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdflux(0:)        ! downward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnet(0:)            ! net longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htr(0:)             ! longwave heating rate (k/day)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totuclfl(0:)        ! clear sky upward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: totdclfl(0:)        ! clear sky downward longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: fnetc(0:)           ! clear sky net longwave flux (w/m2)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(out) :: htrc(0:)            ! clear sky longwave heating rate (k/day)
                                                        !    Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: totufluxs(:,0:)     ! upward longwave flux spectral (w/m2)!linhan change out to in useless
                                                        !    Dimensions: (nbndlw, 0:nlayers)
      real(kind=r8), intent(in) :: totdfluxs(:,0:)     ! downward longwave flux spectral (w/m2)!linhan change out to in useless
                                                        !    Dimensions: (nbndlw, 0:nlayers)

! ----- Local -----
! Declarations for radiative transfer
      real(kind=r8) :: abscld(ngptlw,nlayers)
      real(kind=r8) :: atot(nlayers)
      real(kind=r8) :: atrans(nlayers)
      real(kind=r8) :: bbugas(nlayers)
      real(kind=r8) :: bbutot(nlayers)
      real(kind=r8) :: clrurad(0:nlayers)
      real(kind=r8) :: clrdrad(0:nlayers)
      real(kind=r8) :: efclfrac(ngptlw,nlayers)
      real(kind=r8) :: uflux(0:nlayers)
      real(kind=r8) :: dflux(0:nlayers)
      real(kind=r8) :: urad(0:nlayers)
      real(kind=r8) :: drad(0:nlayers)
      real(kind=r8) :: uclfl(0:nlayers)
      real(kind=r8) :: dclfl(0:nlayers)
      real(kind=r8) :: odcld(ngptlw,nlayers)


      real(kind=r8) :: secdiff(nbndlw)                   ! secant of diffusivity angle
      real(kind=r8) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)  ! diffusivity angle adjustment coefficients
      real(kind=r8) :: wtdiff, rec_6
      real(kind=r8) :: transcld, radld, radclrd, plfrac, blay, dplankup, dplankdn
      real(kind=r8) :: odepth, odtot, odepth_rec, odtot_rec, gassrc
      real(kind=r8) :: tblind, tfactot, bbd, bbdtot, tfacgas, transc, tausfac
      real(kind=r8) :: rad0, reflect, radlu, radclru

      integer :: icldlyr(nlayers)                        ! flag for cloud in layer
      integer :: ibnd, ib, iband, lay, lev, l, ig        ! loop indices
      integer :: igc                                     ! g-point interval counter
      integer :: iclddn                                  ! flag for cloud in down path
      integer :: ittot, itgas, itr                       ! lookup table indices

      !------------------linhan------------------------
      real(kind=r8) :: uflux_sol(0:nlayers)
      real(kind=r8) :: dflux_sol(0:nlayers)
      real(kind=r8) :: uflux_ir(0:nlayers)
      real(kind=r8) :: dflux_ir(0:nlayers)
	  integer :: K,KM1,LP1,iclr_cld,ibegin !linhan 20190927

	  real(kind=r8) :: tau(nlayers), ssa(nlayers), asm(nlayers),bf(nlayers+1),bs,PKAG(nlayers+1)
	  real(kind=r8) :: a(4),u(4)
      real(kind=r8) :: Ru0(nlayers,2),Tu0(nlayers,2),IRu0(nlayers,2),ITu0(nlayers,2), RBAR(nlayers,2,2),TBAR(nlayers,2,2)
      real(kind=r8) :: IRu01(nlayers+1,2),ITu01(nlayers+1,2), IRu0N(nlayers+1,2),ITu0N(nlayers+1,2),&
                       RBARS(nlayers+1,2,2),TBARS(nlayers+1,2,2),RBARN(nlayers+1,2,2),TBARN(nlayers+1,2,2),&
                       RU01(nlayers+1,2),TU01(nlayers+1,2),RU0N(nlayers+1,2),TU0N(nlayers+1,2),TDIR(nlayers+1),DTR(nlayers+1)

	  real(kind=r8) :: C11,C12,C21,C22,PMOD,TMM11,TMM12,TMM21,TMM22,TMR11,TMR12,TMR21,TMR22
	  real(kind=r8) :: RT1,RT2,UU1,UU2,DU1,DU2,BRT11,BRT12,BRT21,BRT22,BMR11,BMR12,BMR21,BMR22 
	  real(kind=r8) :: rt1ir,rt2ir,uu1ir,uu2ir,du1ir,du2ir
	  real(kind=r8) :: fsol
	  real(kind=r8) :: pi
	  real(kind=r8) :: fracs_lh(ngptlw,nlayers), amu0_loc
      logical:: clear_flg=.false.  !true for calculate the clearsky
	  data pi /3.1415927_r8/
      real(kind=r8) :: wt_lh(140) = (/ &
          0.1527534276, 0.1491729617, 0.2737848013, 0.2201246098, 0.1459487156, 0.0471194894, 0.0068539977, &
          0.0036339760, 0.0005330000, 0.0000750000, 0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544, &
          0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116, 0.0471194894, 0.0068539977, 0.0036339760, &
          0.0006080000, 0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544, 0.1181945205, 0.1019300893, &
          0.0832767040, 0.0626720116, 0.0424925000, 0.0046269894, 0.0038279891, 0.0030260086, 0.0022199750, &
          0.0014140010, 0.0005330000, 0.0000750000, 0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544, &
          0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116, 0.0424925000, 0.0046269894, 0.0038279891, &
          0.0030260086, 0.0022199750, 0.0020220010, 0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544, &
          0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116, 0.0424925000, 0.0046269894, 0.0038279891, &
          0.0030260086, 0.0022199750, 0.0014140010, 0.0005330000, 0.0000750000, 0.3019263893, 0.2737848013, &
          0.2201246098, 0.1459487156, 0.0471194894, 0.0068539977, 0.0036339760, 0.0006080000, 0.3019263893, &
          0.2737848013, 0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116, 0.0424925000, 0.0046269894, &
          0.0038279891, 0.0030260086, 0.0036339760, 0.0006080000, 0.3019263893, 0.2737848013, 0.2201246098, &
          0.1459487156, 0.0471194894, 0.0068539977, 0.0036339760, 0.0006080000, 0.1527534276, 0.1491729617, &
          0.1420961469, 0.1316886544, 0.1181945205, 0.1019300893, 0.0832767040, 0.0626720116, 0.0471194894, &
          0.0068539977, 0.0036339760, 0.0006080000, 0.3019263893, 0.2737848013, 0.2201246098, 0.1459487156, &
          0.0539734871, 0.0042419760, 0.1527534276, 0.1491729617, 0.2737848013, 0.2201246098, 0.1459487156, &
          0.0471194894, 0.0090739727, 0.0020220010, 0.1527534276, 0.1491729617, 0.1420961469, 0.1316886544, &
          0.2201246098, 0.1459487156, 0.0539734871, 0.0042419760, 0.4440225362, 0.3518132642, 0.1930682050, &
          0.0110959737, 0.9417845160, 0.0582154631, 0.9417845160, 0.0582154631, 0.5757111906, 0.4242887885/)
! This secant and weight corresponds to the standard diffusivity 
! angle for 2GQ.  This initial value is redefined below.	          
 	data a / 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8 /
	data u / -0.7886752_r8, -0.2113247_r8, 0.2113247_r8, 0.7886752_r8/    
      !------------------linhan end--------------------

! ------- Definitions -------
! input
!    nlayers                      ! number of model layers
!    ngptlw                       ! total number of g-point subintervals
!    nbndlw                       ! number of longwave spectral bands
!    ncbands                      ! number of spectral bands for clouds
!    secdiff                      ! diffusivity angle
!    wtdiff                       ! weight for radiance to flux conversion
!    pavel                        ! layer pressures (mb)
!    pz                           ! level (interface) pressures (mb)
!    tavel                        ! layer temperatures (k)
!    tz                           ! level (interface) temperatures(mb)
!    tbound                       ! surface temperature (k)
!    cldfrac                      ! layer cloud fraction
!    taucloud                     ! layer cloud optical depth
!    itr                          ! integer look-up table index
!    icldlyr                      ! flag for cloudy layers
!    iclddn                       ! flag for cloud in column at any layer
!    semiss                       ! surface emissivities for each band
!    reflect                      ! surface reflectance
!    bpade                        ! 1/(pade constant)
!    tau_tbl                      ! clear sky optical depth look-up table
!    exp_tbl                      ! exponential look-up table for transmittance
!    tfn_tbl                      ! tau transition function look-up table

! local
!    atrans                       ! gaseous absorptivity
!    abscld                       ! cloud absorptivity
!    atot                         ! combined gaseous and cloud absorptivity
!    odclr                        ! clear sky (gaseous) optical depth
!    odcld                        ! cloud optical depth
!    odtot                        ! optical depth of gas and cloud
!    tfacgas                      ! gas-only pade factor, used for planck fn
!    tfactot                      ! gas and cloud pade factor, used for planck fn
!    bbdgas                       ! gas-only planck function for downward rt
!    bbugas                       ! gas-only planck function for upward rt
!    bbdtot                       ! gas and cloud planck function for downward rt
!    bbutot                       ! gas and cloud planck function for upward calc.
!    gassrc                       ! source radiance due to gas only
!    efclfrac                     ! effective cloud fraction
!    radlu                        ! spectrally summed upward radiance 
!    radclru                      ! spectrally summed clear sky upward radiance 
!    urad                         ! upward radiance by layer
!    clrurad                      ! clear sky upward radiance by layer
!    radld                        ! spectrally summed downward radiance 
!    radclrd                      ! spectrally summed clear sky downward radiance 
!    drad                         ! downward radiance by layer
!    clrdrad                      ! clear sky downward radiance by layer

! output
!    totuflux                     ! upward longwave flux (w/m2)
!    totdflux                     ! downward longwave flux (w/m2)
!    fnet                         ! net longwave flux (w/m2)
!    htr                          ! longwave heating rate (k/day)
!    totuclfl                     ! clear sky upward longwave flux (w/m2)
!    totdclfl                     ! clear sky downward longwave flux (w/m2)
!    fnetc                        ! clear sky net longwave flux (w/m2)
!    htrc                         ! clear sky longwave heating rate (k/day)



      fracs_lh=fracs
      hvrrtc = '$Revision: 1.7 $'

      totuflux(0) = 0.0_r8
      totdflux(0) = 0.0_r8
      totuclfl(0) = 0.0_r8
      totdclfl(0) = 0.0_r8

      do lay = 1, nlayers
         totuflux(lay) = 0.0_r8
         totdflux(lay) = 0.0_r8
         totuclfl(lay) = 0.0_r8
         totdclfl(lay) = 0.0_r8
         icldlyr(lay) = 0
      end do

      igc = 1
! Loop over frequency bands.
      do iband = istart, iend
 

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1

! Loop over g-channels.
 1000    continue

!-------------------for fsol in long wave----------------------
         if (amu0 .le. 1.E-10_r8) then
             amu0_loc=1.E-10_r8
         else
             amu0_loc=amu0
         endif

         !linhan changes solar constant 
         if(iband.eq.1) then
         fsol = 0.033726151480236055_r8 * amu0_loc
         elseif(iband.eq.2)  then
         fsol = 6.8623923076864615E-002_r8 * amu0_loc 
         elseif(iband.eq.3)  then
         fsol = 0.1038364918079918_r8 * amu0_loc
         elseif(iband.eq.4)  then
         fsol = 7.6823076651195171E-002_r8 * amu0_loc
         elseif(iband.eq.5)  then
         fsol = 0.1713105360113988_r8 * amu0_loc
         elseif(iband.eq.6)  then
         fsol = 0.3180066751843222_r8 * amu0_loc
         elseif(iband.eq.7)  then
         fsol = 0.2578194707607223_r8 * amu0_loc
         elseif(iband.eq.8)  then
         fsol = 0.3083760699177852_r8 * amu0_loc
         elseif(iband.eq.9)  then
         fsol = 0.8307485663607102_r8 * amu0_loc
         elseif(iband.eq.10)  then
         fsol = 0.4392074384391506_r8 * amu0_loc
         elseif(iband.eq.11)  then
         fsol = 2.003914372237450_r8 * amu0_loc
         elseif(iband.eq.12)  then
         fsol = 2.353256949281163_r8 * amu0_loc
         elseif(iband.eq.13)  then
         fsol = 1.717422236461998_r8 * amu0_loc
         elseif(iband.eq.14)  then
         fsol = 1.527722952002332_r8 * amu0_loc
         elseif(iband.eq.15)  then
         fsol = 3.015639541935780_r8 * amu0_loc
         elseif(iband.eq.16)  then
         !fsol = 11.99388477137011 * amu0_loc
         fsol = 0._r8
         endif
        
        fsol = fsol * wt_lh(igc) !linhan
!fsol = 0.0
!-----------------------------------------------------------------

! Radiative transfer starts here.
!open(8,file='ssa.txt',form='formatted')
        uflux_sol(:) = 0.0_r8
        dflux_sol(:) = 0.0_r8
        uflux_ir(:) = 0.0_r8
        dflux_ir(:) = 0.0_r8
        uflux(:) = 0.0_r8
        dflux(:) = 0.0_r8
    
    
        do  lay=1,nlayers+1
        PKAG(lay) = PLKAVG(WAVENUM1(IBAND),WAVENUM2(IBAND),tz(nlayers-lay+1))
        enddo
        
        do  lay=1,nlayers
        if(fracs_lh(igc, nlayers-lay+1) .eq. 0.0) fracs_lh(igc, nlayers-lay+1) = 10E-16
        bf(lay) = PKAG(lay) * fracs_lh(igc, nlayers-lay+1)
        !linhan debug
        if (bf(lay)<=0)then
            print*,'IR',PKAG(lay),fracs_lh(igc, nlayers-lay+1)
        endif
        enddo
        
        if(fracs_lh(igc, 1) .eq. 0.0) fracs_lh(igc, 1) = 10E-16
        bf(nlayers+1) = PKAG(nlayers+1) * fracs_lh(igc, 1)
        bs = bf(nlayers+1)
        if (.not.clear_flg) then !linhan
            uclfl=0._r8
            dclfl=0._r8
        endif
        ibegin=2 !linhan
        if (clear_flg) ibegin=1 !linhan
        do iclr_cld=ibegin,2  !linhan 1 for clear, 2 for cloud radiation
        do lay = 1, nlayers
           tau(lay) = taut(igc,nlayers-lay+1) + taucmc(igc,nlayers-lay+1)
           if(tau(lay).lt.1.E-15_r8) tau(lay) = 1.E-15_r8
           ssa(lay) = ssacmc(igc,nlayers-lay+1)*taucmc(igc,nlayers-lay+1)/tau(lay)
           ssa(lay) = min(ssa(lay),1._r8-1.E-8_r8)
           ssa(lay) = max(ssa(lay),1.E-8_r8)
           asm(lay) = asmcmc(igc,nlayers-lay+1) 
           asm(lay) = min(asm(lay),1._r8-1.E-8_r8)   
           asm(lay) = max(asm(lay),0.0_r8)   
           if (iclr_cld==1) then 
               tau(lay) = taut(igc,nlayers-lay+1) 
               if(tau(lay).lt.1.E-15_r8) tau(lay) = 1.E-15_r8
               ssa(lay) = 1.E-12_r8
               asm(lay) = 1.E-12_r8
           endif
        end do

! Radiative transfer loop.  

        CALL delfour(tau,ssa,asm,amu0_loc,bf,RBAR,TBAR,RU0,TU0,DTR,IRU0,ITU0,nlayers)

! ... initializaton for TOA.      

        ITu01(1,1)     = 0.0_r8 
        ITu01(1,2)     = 0.0_r8
        Rbars(1,1:2,1) = 0.0_r8
        Rbars(1,1:2,2) = 0.0_r8
        Tu01(1,1)      = 0.0_r8
        Tu01(1,2)      = 0.0_r8
        Tdir(1)        = 1.0_r8

! ... initializaton for the lev layer (surface).                       

        Ru0N(nlayers+1,1:2)    =  1._r8-semiss(iband)
        IRu0N(nlayers+1,1:2)   =  semiss(iband) * bs
        RbarN(nlayers+1,1:2,1) =  (1._r8-semiss(iband))*u(3)
        RbarN(nlayers+1,1:2,2) =  (1._r8-semiss(iband))*u(4)

! ...  ADD THE LAYERS DOWNWARD FROM ONE LAYER BELOW LEV1 TO THE SURFACE.   

      DO K = 2, nlayers+1
        KM1 = K - 1
!-----------------------TMM------------------------------------------
        C11           =  1._r8 - RBAR(KM1,1,1) * RBARS(KM1,1,1) -RBAR(KM1,1,2) * RBARS(KM1,2,1)
        C12           =     - RBAR(KM1,1,1) * RBARS(KM1,1,2) -RBAR(KM1,1,2) * RBARS(KM1,2,2)
        C21           =     - RBAR(KM1,2,1) * RBARS(KM1,1,1) -RBAR(KM1,2,2) * RBARS(KM1,2,1)
        C22           =  1._r8 - RBAR(KM1,2,1) * RBARS(KM1,1,2) -RBAR(KM1,2,2) * RBARS(KM1,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD
!------------------------------BMR----------------------------------------
        BMR11         =  TMM11 * TBAR(KM1,1,1) + TMM12 * TBAR(KM1,2,1)
        BMR12         =  TMM11 * TBAR(KM1,1,2) + TMM12 * TBAR(KM1,2,2)
        BMR21         =  TMM21 * TBAR(KM1,1,1) + TMM22 * TBAR(KM1,2,1)
        BMR22         =  TMM21 * TBAR(KM1,1,2) + TMM22 * TBAR(KM1,2,2)
!------------------------------TDIR---------------------------------------
        TDIR(K)       =  TDIR(KM1) * DTR(KM1)
!------------------------------SOL TU01---------------------------------------
        RT1           =  RU0(KM1,1) * TDIR(KM1) + RBAR(KM1,1,1) *TU01(KM1,1) + RBAR(KM1,1,2) * TU01(KM1,2)
        RT2           =  RU0(KM1,2) * TDIR(KM1) + RBAR(KM1,2,1) *TU01(KM1,1) + RBAR(KM1,2,2) * TU01(KM1,2)
        UU1           =  TMM11 * RT1 + TMM12 * RT2
        UU2           =  TMM21 * RT1 + TMM22 * RT2
        DU1           =  TU01(KM1,1) + RBARS(KM1,1,1) * UU1 +RBARS(KM1,1,2) * UU2
        DU2           =  TU01(KM1,2) + RBARS(KM1,2,1) * UU1 +RBARS(KM1,2,2) * UU2
        TU01(K,1)     =  TU0(KM1,1) * TDIR(KM1) + TBAR(KM1,1,1) *DU1 + TBAR(KM1,1,2) * DU2
        TU01(K,2)     =  TU0(KM1,2) * TDIR(KM1) + TBAR(KM1,2,1) *DU1 + TBAR(KM1,2,2) * DU2
!------------------------------IR ITU01---------------------------------------
        RT1IR         =  IRU0(KM1,1)  + RBAR(KM1,1,1) *ITU01(KM1,1) + RBAR(KM1,1,2) * ITU01(KM1,2)
        RT2IR         =  IRU0(KM1,2)  + RBAR(KM1,2,1) *ITU01(KM1,1) + RBAR(KM1,2,2) * ITU01(KM1,2)
        UU1IR         =  TMM11 * RT1IR + TMM12 * RT2IR
        UU2IR         =  TMM21 * RT1IR + TMM22 * RT2IR
        DU1IR         =  ITU01(KM1,1) + RBARS(KM1,1,1) * UU1IR +RBARS(KM1,1,2) * UU2IR
        DU2IR         =  ITU01(KM1,2) + RBARS(KM1,2,1) * UU1IR +RBARS(KM1,2,2) * UU2IR
        ITU01(K,1)    =  ITU0(KM1,1)  + TBAR(KM1,1,1) * DU1IR + TBAR(KM1,1,2) * DU2IR
        ITU01(K,2)    =  ITU0(KM1,2)  + TBAR(KM1,2,1) * DU1IR + TBAR(KM1,2,2) * DU2IR
!------------------------------- RBARS-------------------------------------
        BRT11         =  TBAR(KM1,1,1) * RBARS(KM1,1,1) +TBAR(KM1,1,2) * RBARS(KM1,2,1)
        BRT12         =  TBAR(KM1,1,1) * RBARS(KM1,1,2) +TBAR(KM1,1,2) * RBARS(KM1,2,2)
        BRT21         =  TBAR(KM1,2,1) * RBARS(KM1,1,1) +TBAR(KM1,2,2) * RBARS(KM1,2,1)
        BRT22         =  TBAR(KM1,2,1) * RBARS(KM1,1,2) +TBAR(KM1,2,2) * RBARS(KM1,2,2)
        RBARS(K,1,1)  =  RBAR(KM1,1,1) + BRT11 * BMR11 + BRT12 * BMR21
        RBARS(K,1,2)  =  RBAR(KM1,1,2) + BRT11 * BMR12 + BRT12 * BMR22
        RBARS(K,2,1)  =  RBAR(KM1,2,1) + BRT21 * BMR11 + BRT22 * BMR21
        RBARS(K,2,2)  =  RBAR(KM1,2,2) + BRT21 * BMR12 + BRT22 * BMR22


! ... ADD THE LAYERS UPWARD FROM LAYER ABOVE SURFACE TO THE LEV1.

        L = nlayers+1 - K + 1
        LP1 = L + 1
!-----------------------TMM--------------------------------------
        C11           =  1._r8 - RBARN(LP1,1,1) * RBAR(L,1,1) -  RBARN(LP1,1,2) * RBAR(L,2,1)
        C12           =     - RBARN(LP1,1,1) * RBAR(L,1,2) -  RBARN(LP1,1,2) * RBAR(L,2,2)
        C21           =     - RBARN(LP1,2,1) * RBAR(L,1,1) -  RBARN(LP1,2,2) * RBAR(L,2,1)
        C22           =  1._r8 - RBARN(LP1,2,1) * RBAR(L,1,2) -  RBARN(LP1,2,2) * RBAR(L,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD
!-----------------------SOL RU0N-----------------------------------
        RT1           =  RU0N(LP1,1) * DTR(L) + RBARN(LP1,1,1) *  TU0(L,1) + RBARN(LP1,1,2) * TU0(L,2)
        RT2           =  RU0N(LP1,2) * DTR(L) + RBARN(LP1,2,1) *  TU0(L,1) + RBARN(LP1,2,2) * TU0(L,2)
        UU1           =  TMM11 * RT1 + TMM12 * RT2
        UU2           =  TMM21 * RT1 + TMM22 * RT2
        DU1           =  TU0(L,1) + RBAR(L,1,1) * UU1 +  RBAR(L,1,2) * UU2
        DU2           =  TU0(L,2) + RBAR(L,2,1) * UU1 +  RBAR(L,2,2) * UU2
        RU0N(L,1)     =  RU0(L,1) + TBAR(L,1,1) * UU1 +  TBAR(L,1,2) * UU2
        RU0N(L,2)     =  RU0(L,2) + TBAR(L,2,1) * UU1 +  TBAR(L,2,2) * UU2
!-----------------------IR RU0N-----------------------------------
        RT1IR         =  IRU0N(LP1,1) + RBARN(LP1,1,1) *  ITU0(L,1) + RBARN(LP1,1,2) * ITU0(L,2)
        RT2IR         =  IRU0N(LP1,2) + RBARN(LP1,2,1) *  ITU0(L,1) + RBARN(LP1,2,2) * ITU0(L,2)
        UU1IR         =  TMM11 * RT1IR + TMM12 * RT2IR
        UU2IR         =  TMM21 * RT1IR + TMM22 * RT2IR
        DU1IR         =  ITU0(L,1) + RBAR(L,1,1) * UU1IR +  RBAR(L,1,2) * UU2IR
        DU2IR         =  ITU0(L,2) + RBAR(L,2,1) * UU1IR +  RBAR(L,2,2) * UU2IR
        IRU0N(L,1)    =  IRU0(L,1) + TBAR(L,1,1) * UU1IR +  TBAR(L,1,2) * UU2IR
        IRU0N(L,2)    =  IRU0(L,2) + TBAR(L,2,1) * UU1IR +  TBAR(L,2,2) * UU2IR
!----------------------------TMM----------------------------------
        C11           =   1._r8 - RBAR(L,1,1) * RBARN(LP1,1,1) -  RBAR(L,1,2) * RBARN(LP1,2,1)
        C12           =      - RBAR(L,1,1) * RBARN(LP1,1,2) -  RBAR(L,1,2) * RBARN(LP1,2,2)
        C21           =      - RBAR(L,2,1) * RBARN(LP1,1,1) -  RBAR(L,2,2) * RBARN(LP1,2,1)
        C22           =   1._r8 - RBAR(L,2,1) * RBARN(LP1,1,2) -  RBAR(L,2,2) * RBARN(LP1,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD
!----------------------------TMR------------------------------------
        TMR11         =  TBAR(L,1,1) * (RBARN(LP1,1,1) * TMM11 +  RBARN(LP1,1,2) * TMM21) + &
                         TBAR(L,1,2) * (RBARN(LP1,2,1) * TMM11 +  RBARN(LP1,2,2) * TMM21)
        TMR12         =  TBAR(L,1,1) * (RBARN(LP1,1,1) * TMM12 +  RBARN(LP1,1,2) * TMM22) + &
                         TBAR(L,1,2) * (RBARN(LP1,2,1) * TMM12 +  RBARN(LP1,2,2) * TMM22)
        TMR21         =  TBAR(L,2,1) * (RBARN(LP1,1,1) * TMM11 +  RBARN(LP1,1,2) * TMM21) + &
                         TBAR(L,2,2) * (RBARN(LP1,2,1) * TMM11 +  RBARN(LP1,2,2) * TMM21)
        TMR22         =  TBAR(L,2,1) * (RBARN(LP1,1,1) * TMM12 +  RBARN(LP1,1,2) * TMM22) + &
                         TBAR(L,2,2) * (RBARN(LP1,2,1) * TMM12 +  RBARN(LP1,2,2) * TMM22)
!-----------------------------RBARN(---------------------------------
        RBARN(L,1,1)  =  RBAR(L,1,1) + TMR11 * TBAR(L,1,1) + TMR12 * TBAR(L,2,1)
        RBARN(L,1,2)  =  RBAR(L,1,2) + TMR11 * TBAR(L,1,2) + TMR12 * TBAR(L,2,2)
        RBARN(L,2,1)  =  RBAR(L,2,1) + TMR21 * TBAR(L,1,1) + TMR22 * TBAR(L,2,1)
        RBARN(L,2,2)  =  RBAR(L,2,2) + TMR21 * TBAR(L,1,2) + TMR22 * TBAR(L,2,2)

      end do

!----------------------------------------------------------------------C
!     ADD DOWNWARD TO CALCULATE THE RESULTANT REFLECTANCES AND         C
!     TRANSMITTANCE AT FLUX LEVELS.                                    C
!----------------------------------------------------------------------C
! open(10,file='flux.txt')

      DO K = 1, nlayers+1
!---------------------------TMM----------------------------------
        C11           =  1._r8 - RBARN(K,1,1) * RBARS(K,1,1) -  RBARN(K,1,2) * RBARS(K,2,1)
        C12           =     - RBARN(K,1,1) * RBARS(K,1,2) -  RBARN(K,1,2) * RBARS(K,2,2)
        C21           =     - RBARN(K,2,1) * RBARS(K,1,1) -  RBARN(K,2,2) * RBARS(K,2,1)
        C22           =  1._r8 - RBARN(K,2,1) * RBARS(K,1,2) -  RBARN(K,2,2) * RBARS(K,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD
!--------------------------SOL FLUX--------------------------------
        RT1           =  RU0N(K,1) * TDIR(K) + RBARN(K,1,1) * TU01(K,1) + RBARN(K,1,2) * TU01(K,2)
        RT2           =  RU0N(K,2) * TDIR(K) + RBARN(K,2,1) * TU01(K,1) + RBARN(K,2,2) * TU01(K,2)
        UU1           =  TMM11 * RT1 + TMM12 * RT2
        UU2           =  TMM21 * RT1 + TMM22 * RT2
        DU1           =  TU01(K,1) + RBARS(K,1,1) * UU1 + RBARS(K,1,2) * UU2
        DU2           =  TU01(K,2) + RBARS(K,2,1) * UU1 + RBARS(K,2,2) * UU2

        uflux_sol(nlayers+1-K)   =  (u(3)*UU1+u(4)*UU2) *  FSOL
        dflux_sol(nlayers+1-K)   =  (u(3)*DU1+u(4)*DU2) *  FSOL + TDIR(K) *  FSOL
!--------------------------IR FLUX--------------------------------
        RT1IR         =  IRU0N(K,1) + RBARN(K,1,1) * ITU01(K,1) + RBARN(K,1,2) * ITU01(K,2)
        RT2IR         =  IRU0N(K,2) + RBARN(K,2,1) * ITU01(K,1) + RBARN(K,2,2) * ITU01(K,2)
        UU1IR         =  TMM11 * RT1IR + TMM12 * RT2IR
        UU2IR         =  TMM21 * RT1IR + TMM22 * RT2IR
        DU1IR         =  ITU01(K,1) + RBARS(K,1,1) * UU1IR + RBARS(K,1,2) * UU2IR
        DU2IR         =  ITU01(K,2) + RBARS(K,2,1) * UU1IR + RBARS(K,2,2) * UU2IR

        uflux_ir(nlayers+1-K)     =  (u(3)*UU1IR+u(4)*UU2IR) * pi
        dflux_ir(nlayers+1-K)     =  (u(3)*DU1IR+u(4)*DU2IR) * pi

        uflux(nlayers+1-K)     =  uflux_sol(nlayers+1-K) + uflux_ir(nlayers+1-K)
        dflux(nlayers+1-K)     =  dflux_sol(nlayers+1-K) + dflux_ir(nlayers+1-K)
      end do

      if (iclr_cld==1)then  !linhan
          uclfl=uflux
          dclfl=dflux
      endif

      enddo ! linhan cloud and clear cycle


! Process longwave output from band for total and clear streams.
! Calculate upward, downward, and net flux.
         do lev = nlayers, 0, -1
            totuflux(lev) = totuflux(lev) + uflux(lev) 
            totdflux(lev) = totdflux(lev) + dflux(lev)
            totuclfl(lev) = totuclfl(lev) + uclfl(lev)
            totdclfl(lev) = totdclfl(lev) + dclfl(lev)
         enddo


! Increment g-point counter
         igc = igc + 1

! Return to continue radiative transfer for all g-channels in present band
         if (igc .le. ngs(iband)) go to 1000

! End spectral band loop
      enddo

! Calculate fluxes at surface
      totuflux(0) = totuflux(0) 
      totdflux(0) = totdflux(0) 
      fnet(0) = totuflux(0) - totdflux(0)
      totuclfl(0) = totuclfl(0) 
      totdclfl(0) = totdclfl(0) 
      fnetc(0) = totuclfl(0) - totdclfl(0)

! Calculate fluxes at model levels
      do lev = 1, nlayers
         totuflux(lev) = totuflux(lev) 
         totdflux(lev) = totdflux(lev) 
         fnet(lev) = totuflux(lev) - totdflux(lev)
         totuclfl(lev) = totuclfl(lev) 
         totdclfl(lev) = totdclfl(lev) 
         fnetc(lev) = totuclfl(lev) - totdclfl(lev)
         l = lev - 1

! Calculate heating rates at model layers
         htr(l)=heatfac*(fnet(l)-fnet(lev))/(pz(l)-pz(lev)) 
         htrc(l)=heatfac*(fnetc(l)-fnetc(lev))/(pz(l)-pz(lev)) 
      enddo

! Set heating rate to zero in top layer
      htr(nlayers) = 0.0_r8
      htrc(nlayers) = 0.0_r8

      end subroutine rtrnmc_4DDA

      Subroutine  delfour(OTAU,OOM,gg,AMU0,bf,RBAR,TBAR,RU0,Tu0,Tdir,IRU0,ITu0,nlayers)

! ---------------------------------------------------------------------------
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      gg     = assymetry factor
!      bf     = planck function
!      OTAU    = optical thickness
!      OOM     = single scattering albedo
!
! outputs
! -------
!      IRU0    : collimated beam reflectivity
!      IRBAR   : diffuse beam reflectivity 
!      ITu0    : collimated beam transmissivity
!      ITBAR   : diffuse beam transmissivity
!
!
! Method:
! -------
! LW 4DDA
!---------------------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer, intent(in) :: nlayers

      real(kind=r8), intent(in) :: gg(:)                      ! asymmetry parameter
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: OTAU(:)                     ! optical depth
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: OOM(:)                      ! single scattering albedo 
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: bf(:)
! ------- Output -------
      real(kind=r8), intent(inout) :: RU0(:,:)                 ! direct beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: RBAR(:,:,:)              ! diffuse beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: TU0(:,:)                 ! direct beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: TBAR(:,:,:)              ! diffuse beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: IRu0(:,:)                 ! direct beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: ITu0(:,:)                 ! direct beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: Tdir(:)             
! ------- Local -------

      integer :: K,L

      real(kind=r8) :: ex1,ex2,qq,q1,q2
      real(kind=r8) :: F,FW,OM,TAU,OMEGA1,OMEGA2,OMEGA3
      real(kind=r8) :: AMU0,RMU0,RMU1,RMU2,eft0,eft1,eft2
      real(kind=r8) :: xx,zz,xx1,xx2,yy1,yy2
      real(kind=r8) :: of0,p2u1,p2u2,p2u3,p2u4,p2f0,p3u1,p3u2,p3u3,p3u4,p3f0
      real(kind=r8) :: bu11,bu12,bu13,bu14,bu21,bu22,bu23,bu24,bu01,bu02,bu03,bu04
      real(kind=r8) :: oa1,oa2,c11,c22,c33,c44,c12p,c13p,c14p,c23p,c24p,c34p
      real(kind=r8) :: c12,c13,c14,c21,c23,c24,c31,c32,c34,c41,c42,c43
      real(kind=r8) :: b22p,b22n,b21p,b21n,b12p,b12n,b11p,b11n
      real(kind=r8) :: a22,a22c,a21,a21c,a12,a12c,a11,a11c
      real(kind=r8) :: bsi,csi,an,rk1,rk2,AA1, AA2,p1p,p1n,p2p,p2n,v1p,v1n,v2p, v2n,psi1,psi2    
      real(kind=r8) :: b2u1p,b2u1n,b1u1p,b1u1n,b2u2p,b2u2n,b1u2p,b1u2n,b2u0p,b2u0n,b1u0p,b1u0n 
      real(kind=r8) :: d11,d12,d13,d14,d21,d22,d23,d24,d31,d32,d33,d34,fa1,fa2,fa3
      real(kind=r8) :: eta11,eta12,eta13,eta14,eta21,eta22,eta23,eta24,eta31,eta32,eta33,eta34
      real(kind=r8) :: H101,H102,H103,H104,H201,H202,H203,H204
      real(kind=r8) :: W11,W12,W13,W14,W21,W22,W23,W24
      real(kind=r8) :: WA,WB,WC,WD,WE,WF,DET,YA,YB,YC,YD,YE,YF,YG,YH
      real(kind=r8) :: V11,V12,V13,V14,V21,V22,V23,V24
      real(kind=r8) :: C1,C2,D1,D2
      real(kind=r8) :: H111,H112,H113,H114,H211,H212,H213,H214
      real(kind=r8) :: H121,H122,H123,H124,H221,H222,H223,H224
      real(kind=r8) :: rmu0ir,eft0ir,bu01ir,bu02ir,bu03ir,bu04ir
      real(kind=r8) :: b2u0nir,b1u0nir,d33ir,d31ir,d34ir,d32ir,fa3ir
      real(kind=r8) :: eta31ir,eta32ir,eta33ir,eta34ir,h101ir,h102ir,h103ir,h104ir
      real(kind=r8) :: h201ir,h202ir,h203ir,h204ir,c1ir,d1ir,c2ir,d2ir

      real(kind=r8), parameter :: eps = 1.e-08_r8

      real(kind=r8) :: a(4),u(4)
      data a / 0.5_r8, 0.5_r8, 0.5_r8, 0.5_r8 /
      data u / -0.7886752_r8, -0.2113247_r8, 0.2113247_r8, 0.7886752_r8/    
!-------------------------------------------------------------------

  DO K = 1, nlayers

!--------------- adjust -----------------------------------------
      F            =  gg(K)**4
      FW           =  1.0_r8 - OOM(K) * F
      OM           =  OOM(K) * (1.0_r8 - F) / FW
      TAU          =  OTAU(K) * FW

      OMEGA1       =  (3.0_r8*gg(K) - 3.0_r8 * F) / (1.0_r8 - F)
      OMEGA2       =  (5.0_r8*gg(K)**2 - 5.0_r8 * F) / (1.0_r8 - F)
      OMEGA3       =  (7.0_r8*gg(K)**3 - 7.0_r8 * F) / (1.0_r8 - F)
!----------------------------------------------------------------
      q1          =  dlog ( bf(K+1) / bf(K) )
      q2          =  1.0_r8 /  tau
      RMU0IR       =  - q1 * q2   ! f0

      RMU0         =  1.0_r8 / AMU0   ! f0
      RMU1         =  1.0_r8 / u(3)   ! f1
      RMU2         =  1.0_r8 / u(4)   ! f2
      eft0         =  exp(-RMU0*tau)
      eft0IR       =  exp(-RMU0IR*tau)
      eft1         =  exp(-RMU1*tau)
      eft2         =  exp(-RMU2*tau)
     
!----------------------------SOL--------------------------------------
      of0          =  OM/4.0_r8
        
      p2u1         =  1.5_r8 * u(1)**2 - 0.5_r8
      p2u2         =  1.5_r8 * u(2)**2 - 0.5_r8
      p2u3         =  1.5_r8 * u(3)**2 - 0.5_r8
      p2u4         =  1.5_r8 * u(4)**2 - 0.5_r8
      p2f0         =  1.5_r8 * (-AMU0)**2 - 0.5_r8

      p3u1         =  ( 2.5_r8 * u(1)**2 - 1.5_r8 ) * u(1)
      p3u2         =  ( 2.5_r8 * u(2)**2 - 1.5_r8 ) * u(2)
      p3u3         =  ( 2.5_r8 * u(3)**2 - 1.5_r8 ) * u(3)
      p3u4         =  ( 2.5_r8 * u(4)**2 - 1.5_r8 ) * u(4)
      p3f0         =  ( 2.5_r8 * (-AMU0)**2 - 1.5_r8 ) * (-AMU0)
    
        bu11         =  of0/u(1)*( 1.0+OMEGA1*u(1)*(-u(3)) +OMEGA2*p2u1*p2u2+OMEGA3*p3u1*p3u2 ) 
        bu12         =  of0/u(2)*( 1.0+OMEGA1*u(2)*(-u(3))+OMEGA2*p2u2*p2u2+OMEGA3*p3u2*p3u2 )  
        bu13         =  of0/u(3)*( 1.0+OMEGA1*u(3)*(-u(3))+OMEGA2*p2u3*p2u2+OMEGA3*p3u3*p3u2 )  
        bu14         =  of0/u(4)*( 1.0+OMEGA1*u(4)*(-u(3))+OMEGA2*p2u4*p2u2+OMEGA3*p3u4*p3u2 ) 

        bu21         =  of0/u(1)*( 1.0+OMEGA1*u(1)*(-u(4))+OMEGA2*p2u1*p2u1+OMEGA3*p3u1*p3u1 )
        bu22         =  of0/u(2)*( 1.0+OMEGA1*u(2)*(-u(4))+OMEGA2*p2u2*p2u1+OMEGA3*p3u2*p3u1 )
        bu23         =  of0/u(3)*( 1.0+OMEGA1*u(3)*(-u(4))+OMEGA2*p2u3*p2u1+OMEGA3*p3u3*p3u1 )
        bu24         =  of0/u(4)*( 1.0+OMEGA1*u(4)*(-u(4))+OMEGA2*p2u4*p2u1+OMEGA3*p3u4*p3u1 )

        bu01         =  of0/u(1)*( 1.0+OMEGA1*u(1)*(-AMU0)+OMEGA2*p2u1*p2f0+OMEGA3*p3u1*p3f0 )
        bu02         =  of0/u(2)*( 1.0+OMEGA1*u(2)*(-AMU0)+OMEGA2*p2u2*p2f0+OMEGA3*p3u2*p3f0 )
        bu03         =  of0/u(3)*( 1.0+OMEGA1*u(3)*(-AMU0)+OMEGA2*p2u3*p2f0+OMEGA3*p3u3*p3f0 )
        bu04         =  of0/u(4)*( 1.0+OMEGA1*u(4)*(-AMU0)+OMEGA2*p2u4*p2f0+OMEGA3*p3u4*p3f0 )
!--------------------------IR------------------------------------------       
        xx             =  (1.0_r8-OM)*bf(K)
        bu01IR         =  xx/u(1)  
        bu02IR         =  xx/u(2)
        bu03IR         =  -bu02IR
        bu04IR         =  -bu01IR
!--------------------------- c(4,4) -----------------------------------
        oa1          =  OM/2.0_r8*a(1)   ! oa1=oa4
        oa2          =  OM/2.0_r8*a(2)   ! oa2=oa3

        c11          =  (oa1*(1.0_r8+OMEGA1*u(1)**2+OMEGA2*p2u1**2+OMEGA3*p3u1**2)-1.0_r8)/u(1)
        c22          =  (oa2*(1.0_r8+OMEGA1*u(2)**2+OMEGA2*p2u2**2+OMEGA3*p3u2**2)-1.0_r8)/u(2)
        c33          =  (oa2*(1.0_r8+OMEGA1*u(3)**2+OMEGA2*p2u3**2+OMEGA3*p3u3**2)-1.0_r8)/u(3)
        c44          =  (oa1*(1.0_r8+OMEGA1*u(4)**2+OMEGA2*p2u4**2+OMEGA3*p3u4**2)-1.0_r8)/u(4)

        c12p         =  1.0_r8+OMEGA1*u(1)*u(2)+OMEGA2*p2u1*p2u2+OMEGA3*p3u1*p3u2
        c13p         =  1.0_r8+OMEGA1*u(1)*u(3)+OMEGA2*p2u1*p2u3+OMEGA3*p3u1*p3u3
        c14p         =  1.0_r8+OMEGA1*u(1)*u(4)+OMEGA2*p2u1*p2u4+OMEGA3*p3u1*p3u4
        c23p         =  1.0_r8+OMEGA1*u(2)*u(3)+OMEGA2*p2u2*p2u3+OMEGA3*p3u2*p3u3
        c24p         =  1.0_r8+OMEGA1*u(2)*u(4)+OMEGA2*p2u2*p2u4+OMEGA3*p3u2*p3u4
        c34p         =  1.0_r8+OMEGA1*u(3)*u(4)+OMEGA2*p2u3*p2u4+OMEGA3*p3u3*p3u4

        c12          =  c12p*oa2/u(1)
        c13          =  c13p*oa2/u(1)
        c14          =  c14p*oa1/u(1)
        c21          =  c12p*oa1/u(2)  ! c21p=c12p
        c23          =  c23p*oa2/u(2)  
        c24          =  c24p*oa1/u(2)
        c31          =  c13p*oa1/u(3)  ! c31p=c13p
        c32          =  c23p*oa2/u(3)  ! c32p=c23p
        c34          =  c34p*oa1/u(3)  
        c41          =  c14p*oa1/u(4)  ! c41p=c14p
        c42          =  c24p*oa2/u(4)  ! c42p=c24p
        c43          =  c34p*oa2/u(4)  ! c43p=c34p
!-----------------------------------------------------------------------
      b22p         =  c44+c41
      b22n         =  c44-c41
      b21p         =  c43+c42
      b21n         =  c43-c42
      b12p         =  c34+c31
      b12n         =  c34-c31
      b11p         =  c33+c32
      b11n         =  c33-c32
!--------------------a22 a22C----------------------------------------
      a22          =  b22p*b22n+b12p*b21n
      a22c         =  b22n*b22p+b12n*b21p
      a21          =  b22n*b21p+b21n*b11p
      a21c         =  b22p*b21n+b21p*b11n
      a12          =  b12n*b22p+b11n*b12p
      a12c         =  b12p*b22n+b11p*b12n
      a11          =  b12n*b21p+b11n*b11p
      a11c         =  b12p*b21n+b11p*b11n

!-------------------b,an--a^-,ap---a^+,rk(1)---k(1),rk(2)----k(2)----
      bsi          =  a22+a11
      csi          =  a21*a12-a11*a22
      an           =  b22n*b11n-b12n*b21n
      
      if ((bsi**2+4_r8*csi)<0) then
        zz         = 1.E-16_r8 ! =  -(-(bsi**2+4*csi))**0.5/2.0 !linhan 
      else
        zz           =  (bsi**2+4_r8*csi)**0.5_r8/2.0_r8
      endif
      rk1          =  (bsi/2.0_r8+zz)**0.5_r8
      rk2          =  (bsi/2.0_r8-zz)**0.5_r8
      AA1          =  (rk1**2-a22)/a21
      AA2          =  (rk2**2-a22)/a21
!-------------------------------- phi ----------------------------------
      xx1          =  (AA1*b22n-b12n)/an*rk1
      xx2          =  (AA2*b22n-b12n)/an*rk2
      p1p          =  0.5_r8*(AA1+xx1)
      p1n          =  0.5_r8*(AA1-xx1)
      p2p          =  0.5_r8*(AA2+xx2)
      p2n          =  0.5_r8*(AA2-xx2)
!---------------------- varphi ---- psi --------------------------------
      yy1          =  (b11n-AA1*b21n)/an*rk1
      yy2          =  (b11n-AA2*b21n)/an*rk2
      v1p          =  0.5_r8*(1+yy1)
      v1n          =  0.5_r8*(1-yy1)
      v2p          =  0.5_r8*(1+yy2)
      v2n          =  0.5_r8*(1-yy2)
      psi1         =  exp(-rk1*tau)
      psi2         =  exp(-rk2*tau)
!---------------dependent on u1,u2,u0---------------
      b2u1p        =  bu14+bu11 !SOL
      b2u1n        =  bu14-bu11
      b1u1p        =  bu13+bu12
      b1u1n        =  bu13-bu12
    
      b2u2p        =  bu24+bu21
      b2u2n        =  bu24-bu21
      b1u2p        =  bu23+bu22
      b1u2n        =  bu23-bu22
    
      b2u0p        =  bu04+bu01
      b2u0n        =  bu04-bu01
      b1u0p        =  bu03+bu02
      b1u0n        =  bu03-bu02
    
      b2u0nIR      =  bu04IR-bu01IR !IR
      b1u0nIR      =  bu03IR-bu02IR

!----------d(:,1,:) normal; d(:,2,:)=d` acute;d(:,:,1) u1;d(:,:,3) u2; d(:,:,3) u0;------------------------c
!             fa=f`
        d13        =  b22n*b2u1n+b21n*b1u1n+b2u1p/u(3)   !SOL
        d11        =  b12n*b2u1n+b11n*b1u1n+b1u1p/u(3)
        d14        =  b22p*b2u1p+b21p*b1u1p+b2u1n/u(3)
        d12        =  b12p*b2u1p+b11p*b1u1p+b1u1n/u(3)
    
        d23        =  b22n*b2u2n+b21n*b1u2n+b2u2p/u(4)
        d21        =  b12n*b2u2n+b11n*b1u2n+b1u2p/u(4)
        d24        =  b22p*b2u2p+b21p*b1u2p+b2u2n/u(4)
        d22        =  b12p*b2u2p+b11p*b1u2p+b1u2n/u(4)
    
        d33        =  b22n*b2u0n+b21n*b1u0n+b2u0p*RMU0
        d31        =  b12n*b2u0n+b11n*b1u0n+b1u0p*RMU0
        d34        =  b22p*b2u0p+b21p*b1u0p+b2u0n*RMU0
        d32        =  b12p*b2u0p+b11p*b1u0p+b1u0n*RMU0
    
        d33IR      =  b22n*b2u0nIR+b21n*b1u0nIR !IR
        d31IR      =  b12n*b2u0nIR+b11n*b1u0nIR
        d34IR      =  b2u0nIR*RMU0IR
        d32IR      =  b1u0nIR*RMU0IR
    
        fa1        =  RMU1**4-bsi*RMU1**2-csi 
        fa2        =  RMU2**4-bsi*RMU2**2-csi
        fa3        =  RMU0**4-bsi*RMU0**2-csi
        fa3IR      =  RMU0IR**4-bsi*RMU0IR**2-csi  !IR

!-----------------eta(1,:,:) eta1;eta(:,1,:) normal; eta(:,2,:) acute;d(:,:,1) u1;d(:,:,3) u2; d(:,:,3) u0;
      eta11        =  (d11*RMU1**2+a12*d13-a22*d11)/fa1 !SOL
      eta13        =  (d13*RMU1**2+a21*d11-a11*d13)/fa1
      eta12        =  (d12*RMU1**2+a12c*d14-a22c*d12)/fa1
      eta14        =  (d14*RMU1**2+a21c*d12-a11c*d14)/fa1
    
      eta21        =  (d21*RMU2**2+a12*d23-a22*d21)/fa2
      eta23        =  (d23*RMU2**2+a21*d21-a11*d23)/fa2
      eta22        =  (d22*RMU2**2+a12c*d24-a22c*d22)/fa2
      eta24        =  (d24*RMU2**2+a21c*d22-a11c*d24)/fa2
    
      eta31        =  (d31*RMU0**2+a12*d33-a22*d31)/fa3
      eta33        =  (d33*RMU0**2+a21*d31-a11*d33)/fa3
      eta32        =  (d32*RMU0**2+a12c*d34-a22c*d32)/fa3
      eta34        =  (d34*RMU0**2+a21c*d32-a11c*d34)/fa3
    
      eta31IR      =  (d31IR*RMU0IR**2+a12*d33IR-a22*d31IR)/fa3IR   !IR
      eta33IR      =  (d33IR*RMU0IR**2+a21*d31IR-a11*d33IR)/fa3IR
      eta32IR      =  (d32IR*RMU0IR**2+a12c*d34IR-a22c*d32IR)/fa3IR
      eta34IR      =  (d34IR*RMU0IR**2+a21c*d32IR-a11c*d34IR)/fa3IR

! ...  DIRECT SOLUTION

!---------------------- H1, H2 -----------------------------------
      H101          =  -0.5_r8*(eta33+eta34)*eft0   !SOL
      H102          =  -0.5_r8*(eta31+eta32)*eft0
      H103          =  -0.5_r8*(eta31-eta32)
      H104          =  -0.5_r8*(eta33-eta34)
    
      H201          =  0.5_r8*(eta33+eta34)
      H202          =  0.5_r8*(eta31+eta32)
      H203          =  0.5_r8*(eta31-eta32)*eft0
      H204          =  0.5_r8*(eta33-eta34)*eft0
    
      H101IR        =  -0.5_r8*(eta33IR+eta34IR)*eft0IR   !IR
      H102IR        =  -0.5_r8*(eta31IR+eta32IR)*eft0IR
      H103IR        =  -0.5_r8*(eta31IR-eta32IR)
      H104IR        =  -0.5_r8*(eta33IR-eta34IR)
    
      H201IR        =  0.5_r8*(eta33IR+eta34IR)
      H202IR        =  0.5_r8*(eta31IR+eta32IR)
      H203IR        =  0.5_r8*(eta31IR-eta32IR)*eft0IR
      H204IR        =  0.5_r8*(eta33IR-eta34IR)*eft0IR
!------------------- A1, A1I-------------------------------------- 
!   W31=W24  W41=W14          WI31=WI24  WI41=WI14
!   W32=W23  W42=W13          WI32=WI23  WI42=WI13
!   W33=W22  W43=W12          WI33=WI22  WI43=WI12
!   W34=W21  W44=W11          WI34=WI21  WI44=WI11
!-----------------------------------------------------------------
        W11          =  v2p*psi2  
        W12          =  v1p*psi1 
        W13          =  v1n
        W14          =  v2n
        W21          =  p2p*psi2
        W22          =  p1p*psi1
        W23          =  p1n
        W24          =  p2n

        WA           =  W11 * W22 - W21 * W12  
        WB           =  W14 * W23 - W24 * W13
        WC           =  W11 * W23 - W21 * W13 
        WD           =  W11 * W24 - W21 * W14
        WE           =  W12 * W23 - W22 * W13
        WF           =  W12 * W24 - W22 * W14

        DET          =  1.0_r8 / (2.0_r8 * WD * WE + WA * WA - WC * WC + WB * WB - WF * WF)
        YA           = ( W22 * WA - W23 * WC + W24 * WE) * DET
        YB           = (-W12 * WA + W13 * WC - W14 * WE) * DET
        YC           = ( W11 * WE - W12 * WF - W13 * WB) * DET
        YD           = (-W21 * WE + W22 * WF + W23 * WB) * DET
        YE           = (-W21 * WA + W23 * WD - W24 * WF) * DET
        YF           = ( W11 * WA - W13 * WD + W14 * WF) * DET
        YG           = (-W11 * WC + W12 * WD + W14 * WB) * DET
        YH           = ( W21 * WC - W22 * WD - W24 * WB) * DET
!------------------------SOL G1 -----------------------------------------
        C1           =  YA * H101 + YB * H102 + YC * H103 + YD * H104
        D1           =  YE * H101 + YF * H102 + YG * H103 + YH * H104
        C2           =  YH * H101 + YG * H102 + YF * H103 + YE * H104
        D2           =  YD * H101 + YC * H102 + YB * H103 + YA * H104
!------------------------IRL G1 -----------------------------------------
        C1IR         =  YA * H101IR + YB * H102IR + YC * H103IR + YD * H104IR
        D1IR         =  YE * H101IR + YF * H102IR + YG * H103IR + YH * H104IR
        C2IR         =  YH * H101IR + YG * H102IR + YF * H103IR + YE * H104IR
        D2IR         =  YD * H101IR + YC * H102IR + YB * H103IR + YA * H104IR
!------------------------ A2 -----------------------------------------
        V11          =  v2p 
        V12          =  v1p 
        V13          =  v1n*psi1
        V14          =  v2n*psi2 
        V21          =  p2p
        V22          =  p1p
        V23          =  p1n*psi1
        V24          =  p2n*psi2
!------------------------ Ru0, Tu0 ----------------------------------
        RU0(K,1)     = (V21*C1+V22*D1+V23*C2+V24*D2+H202) * RMU0
        RU0(K,2)     = (V11*C1+V12*D1+V13*C2+V14*D2+H201) * RMU0
        TU0(K,1)     = (V24*C1+V23*D1+V22*C2+V21*D2+H203) * RMU0
        TU0(K,2)     = (V14*C1+V13*D1+V12*C2+V11*D2+H204) * RMU0
        TDIR(K)      =  eft0
!------------------------ IRu0, ITu0 ----------------------------------
        IRU0(K,1)    = V21*C1IR+V22*D1IR+V23*C2IR+V24*D2IR+H202IR
        IRU0(K,2)    = V11*C1IR+V12*D1IR+V13*C2IR+V14*D2IR+H201IR
        ITU0(K,1)    = V24*C1IR+V23*D1IR+V22*C2IR+V21*D2IR+H203IR
        ITU0(K,2)    = V14*C1IR+V13*D1IR+V12*C2IR+V11*D2IR+H204IR

! ...   DIFFUSE SOLUTION

!----------------------- H1 -------------------------------------
       H111          =  0._r8
       H112          =  0._r8
       H113          =  1._r8
       H114          =  0._r8
!------------------------ G1 -----------------------------------------
        C1           =  YA * H111 + YB * H112 + YC * H113 + YD * H114
        D1           =  YE * H111 + YF * H112 + YG * H113 + YH * H114
        C2           =  YH * H111 + YG * H112 + YF * H113 + YE * H114
        D2           =  YD * H111 + YC * H112 + YB * H113 + YA * H114

        RBAR(K,1,1)  =  V21*C1+V22*D1+V23*C2+V24*D2
        RBAR(K,2,1)  =  V11*C1+V12*D1+V13*C2+V14*D2
        TBAR(K,1,1)  =  V24*C1+V23*D1+V22*C2+V21*D2
        TBAR(K,2,1)  =  V14*C1+V13*D1+V12*C2+V11*D2
!------------------------ H1, H2 -------------------------------------
       H121          =  0._r8
       H122          =  0._r8
       H123          =  0._r8
       H124          =  1._r8
!------------------------ G1 -----------------------------------------
        C1           =  YA * H121 + YB * H122 + YC * H123 + YD * H124
        D1           =  YE * H121 + YF * H122 + YG * H123 + YH * H124
        C2           =  YH * H121 + YG * H122 + YF * H123 + YE * H124
        D2           =  YD * H121 + YC * H122 + YB * H123 + YA * H124

        RBAR(K,1,2)  =  V21*C1+V22*D1+V23*C2+V24*D2
        RBAR(K,2,2)  =  V11*C1+V12*D1+V13*C2+V14*D2
        TBAR(K,1,2)  =  V24*C1+V23*D1+V22*C2+V21*D2
        TBAR(K,2,2)  =  V14*C1+V13*D1+V12*C2+V11*D2
  end do
      END SUBROUTINE delfour

REAL FUNCTION PLKAVG( WNUMLO, WNUMHI, T )

!        Computes Planck function integrated between two wavenumbers
!
!  INPUT :  WNUMLO : Lower wavenumber (inv cm) of spectral interval
!
!           WNUMHI : Upper wavenumber
!
!           T      : Temperature (K)
!
!  OUTPUT : PLKAVG : Integrated Planck function ( Watts/sq m )
!                      = Integral (WNUMLO to WNUMHI) of
!                        2h c**2  nu**3 / ( EXP(hc nu/kT) - 1)
!                        (where h=Plancks constant, c=speed of
!                         light, nu=wavenumber, T=temperature,
!                         and k = Boltzmann constant)
!
!  Reference : Specifications of the Physical World: New Value
!                 of the Fundamental Constants, Dimensions/N.B.S.,
!                 Jan. 1974
!   Calls- D1MACH, ERRMSG
! ----------------------------------------------------------------------

!     .. Parameters ..

      REAL      A1, A2, A3, A4, A5, A6
      PARAMETER ( A1 = 1. / 3., A2 = -1. / 8., A3 = 1. / 60., &
               A4 = -1. / 5040., A5 = 1. / 272160., A6 = -1. / 13305600. )

!     .. Scalar Arguments ..

      REAL(kind=r8)     T, WNUMHI, WNUMLO

!     .. Local Scalars ..

      INTEGER   I, K, M, MMAX, N, SMALLV
      REAL      C2, CONC, DEL, EPSIL, EX, EXM, HH, MV, OLDVAL, PI, &
               SIGDPI, SIGMA, VAL, VAL0, VCUT, VMAX, VSQ, X	, B

!     .. Local Arrays ..

      REAL      D( 2 ), P( 2 ), V( 2 ), VCP( 7 )

!     .. Intrinsic Functions ..

      INTRINSIC ABS, ASIN, EXP, LOG, MOD

!     .. Statement Functions ..

      REAL      PLKF

      SAVE      PI, CONC, VMAX, EPSIL, SIGDPI

      DATA      C2 / 1.438786 / , SIGMA / 5.67032E-8 / , VCUT / 1.5 / ,	 &
                VCP / 10.25, 5.7, 3.9, 2.9, 2.3, 1.9, 0.0 /
      DATA      PI / 0.0 /

!     .. Statement Function definitions ..

      PLKF( X ) = X**3 / ( EXP( X ) - 1 )
      B = RADIX(X)
!     ..

      IF( PI .EQ. 0.0 ) THEN

         PI     = 2.*ASIN( 1.0 )
         VMAX   = LOG( HUGE(X) )
         EPSIL  = B**(1-DIGITS(X))
         SIGDPI = SIGMA / PI
         CONC   = 15. / PI**4

      END IF

      IF( T .LT. 1.E-4 ) THEN

         PLKAVG = 0.0
         RETURN

      END IF


      V( 1 ) = C2*WNUMLO / T
      V( 2 ) = C2*WNUMHI / T

      IF( V( 1 ).GT.EPSIL .AND. V( 2 ).LT.VMAX .AND.( WNUMHI - WNUMLO ) / WNUMHI .LT. 1.E-2 ) THEN

!                          ** Wavenumbers are very close.  Get integral
!                          ** by iterating Simpson rule to convergence.

         HH     = V( 2 ) - V( 1 )
         OLDVAL = 0.0
         VAL0   = PLKF( V( 1 ) ) + PLKF( V( 2 ) )

         DO N = 1, 10

            DEL  = HH / ( 2*N )
            VAL  = VAL0

            DO K = 1, 2*N - 1
               VAL  = VAL + 2*( 1 + MOD( K,2 ) )*PLKF( V( 1 ) + K*DEL )
            enddo

            VAL  = DEL / 3.*VAL
            IF( ABS( ( VAL - OLDVAL ) / VAL ).LE.1.E-6 ) GO TO  30
            OLDVAL = VAL

         enddo

   30    CONTINUE

         PLKAVG = SIGDPI * T**4 * CONC * VAL

         RETURN

      END IF

!                          *** General case ***
      SMALLV = 0

      DO 60 I = 1, 2

         IF( V( I ).LT.VCUT ) THEN
!                                   ** Use power series
            SMALLV = SMALLV + 1
            VSQ    = V( I )**2
            P( I ) = CONC*VSQ*V( I )*( A1 +V( I )*( A2 + V( I )*( A3 + VSQ*( A4 + VSQ*( A5 +VSQ*A6 ) ) ) ) )

         ELSE
!                      ** Use exponential series
            MMAX  = 0
!                                ** Find upper limit of series
   40       CONTINUE
            MMAX  = MMAX + 1

            IF( V(I) .LT. VCP( MMAX ) ) GO TO  40

            EX     = EXP( - V(I) )
            EXM    = 1.0
            D( I ) = 0.0

            DO 50 M = 1, MMAX
               MV     = M*V( I )
               EXM    = EX*EXM
               D( I ) = D( I ) + EXM*( 6.+ MV*( 6.+ MV*( 3.+ MV ) ) )/ M**4
   50       CONTINUE

            D( I ) = CONC*D( I )

         END IF

   60 CONTINUE

!                              ** Handle ill-conditioning
      IF( SMALLV.EQ.2 ) THEN
!                                    ** WNUMLO and WNUMHI both small
         PLKAVG = P( 2 ) - P( 1 )

      ELSE IF( SMALLV.EQ.1 ) THEN
!                                    ** WNUMLO small, WNUMHI large
         PLKAVG = 1.- P( 1 ) - D( 2 )

      ELSE
!                                    ** WNUMLO and WNUMHI both large
         PLKAVG = D( 1 ) - D( 2 )

      END IF

      PLKAVG = SIGDPI * T**4 * PLKAVG

  end function PLKAVG
!--------------------------linhan end----------------------------

      end module rrtmg_lw_rtrnmc

