!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_spcvrt.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:14 $

      module rrtmg_sw_spcvrt

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------
      use grist_constants,            only: r8

!      use parkind, only : jpim, jprb
      use parrrsw, only : nbndsw, ngptsw, mxmol, jpband
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrspv, hnamspv
      use rrsw_wvn, only : ngc, ngs
      use rrtmg_sw_reftra, only: reftra_sw
      use rrtmg_sw_taumol, only: taumol_sw
      use rrtmg_sw_vrtqdr, only: vrtqdr_sw

      implicit none

      contains

! ---------------------------------------------------------------------------
      subroutine spcvrt_sw &
            (ncol, iplon, nlayers, istart, iend, icpr, idelm, iout, dotau, &
             pavel, tavel, pz, tz, tbound, palbd, palbp, &
             pclfr, ptauc, pasyc, pomgc, ptaucorig, &
             ptaua, pasya, pomga, prmu0, coldry, wkl, adjflux, &
             laytrop, layswtch, laylow, jp, jt, jt1, &
             co2mult, colch4, colco2, colh2o, colmol, coln2o, colo2, colo3, &
             fac00, fac01, fac10, fac11, &
             selffac, selffrac, indself, forfac, forfrac, indfor, &
             pbbfd, pbbfu, pbbcd, pbbcu, puvfd, puvcd, pnifd, pnicd, pnifu, pnicu, &
             pbbfddir, pbbcddir, puvfddir, puvcddir, pnifddir, pnicddir)
! ---------------------------------------------------------------------------
!
! Purpose: Contains spectral loop to compute the shortwave radiative fluxes, 
!          using the two-stream method of H. Barker. 
!
! Interface:  *spcvrt_sw* is called from *rrtmg_sw.F90* or rrtmg_sw.1col.F90*
!
! Method:
!    Adapted from two-stream model of H. Barker;
!    Two-stream model options (selected with kmodts in rrtmg_sw_reftra.F90):
!        1: Eddington, 2: PIFM, Zdunkowski et al., 3: discret ordinates
!
! Modifications:
!
! Original: H. Barker
! Revision: Merge with RRTMG_SW: J.-J.Morcrette, ECMWF, Feb 2003
! Revision: Add adjustment for Earth/Sun distance : MJIacono, AER, Oct 2003
! Revision: Bug fix for use of PALBP and PALBD: MJIacono, AER, Nov 2003
! Revision: Bug fix to apply delta scaling to clear sky: AER, Dec 2004
! Revision: Code modified so that delta scaling is not done in cloudy profiles
!           if routine cldprop is used; delta scaling can be applied by swithcing
!           code below if cldprop is not used to get cloud properties. 
!           AER, Jan 2005
! Revision: Uniform formatting for RRTMG: MJIacono, AER, Jul 2006 
! Revision: Use exponential lookup table for transmittance: MJIacono, AER, 
!           Aug 2007 
!
! ------------------------------------------------------------------

! ------- Declarations ------

! -------- Input -------

      integer, intent(in) :: ncol
      integer, intent(in) :: iplon                      ! column loop index
      integer, intent(in) :: nlayers
      integer, intent(in) :: istart
      integer, intent(in) :: iend
      integer, intent(in) :: icpr
      integer, intent(in) :: idelm     ! delta-m scaling flag
                                       ! [0 = direct and diffuse fluxes are unscaled]
                                       ! [1 = direct and diffuse fluxes are scaled]
      integer, intent(in) :: iout
      logical, intent(in) :: dotau                      ! True -> do tau calculation

      integer, intent(in) :: laytrop
      integer, intent(in) :: layswtch
      integer, intent(in) :: laylow

      integer, intent(in) :: indfor(:)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: indself(:)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: jp(:)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: jt(:)
                                                                 !   Dimensions: (nlayers)
      integer, intent(in) :: jt1(:)
                                                                 !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: pavel(:)                    ! layer pressure (hPa, mb) 
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: tavel(:)                    ! layer temperature (K)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: pz(0:)                      ! level (interface) pressure (hPa, mb)
                                                                 !   Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: tz(0:)                      ! level temperatures (hPa, mb)
                                                                 !   Dimensions: (0:nlayers)
      real(kind=r8), intent(in) :: tbound                      ! surface temperature (K)
      real(kind=r8), intent(in) :: wkl(:,:)                    ! molecular amounts (mol/cm2) 
                                                                 !   Dimensions: (mxmol,nlayers)
      real(kind=r8), intent(in) :: coldry(:)                   ! dry air column density (mol/cm2)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colmol(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: adjflux(:)                  ! Earth/Sun distance adjustment
                                                                 !   Dimensions: (jpband)

      real(kind=r8), intent(in) :: palbd(:)                    ! surface albedo (diffuse)
                                                                 !   Dimensions: (nbndsw)
      real(kind=r8), intent(in) :: palbp(:)                    ! surface albedo (direct)
                                                                 !   Dimensions: (nbndsw)
      real(kind=r8), intent(in) :: prmu0                       ! cosine of solar zenith angle
      real(kind=r8), intent(in) :: pclfr(:)                    ! cloud fraction
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: ptauc(:,:)                  ! cloud optical depth
                                                                 !   Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: pasyc(:,:)                  ! cloud asymmetry parameter
                                                                 !   Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: pomgc(:,:)                  ! cloud single scattering albedo
                                                                 !   Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: ptaucorig(:,:)              ! cloud optical depth, non-delta scaled
                                                                 !   Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: ptaua(:,:)                  ! aerosol optical depth
                                                                 !   Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: pasya(:,:)                  ! aerosol asymmetry parameter
                                                                 !   Dimensions: (nbndsw,nlayers)
      real(kind=r8), intent(in) :: pomga(:,:)                  ! aerosol single scattering albedo
                                                                 !   Dimensions: (nbndsw,nlayers)

      real(kind=r8), intent(in) :: colh2o(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colco2(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colch4(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: co2mult(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colo3(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: colo2(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: coln2o(:)
                                                                 !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: forfac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: forfrac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: selffac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: selffrac(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac00(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac01(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac10(:)
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: fac11(:)
                                                                 !   Dimensions: (nlayers)

! ------- Output -------
                                                                 !   All Dimensions: (nlayers+1)
      real(kind=r8), intent(out) :: pbbcd(:)
      real(kind=r8), intent(out) :: pbbcu(:)
      real(kind=r8), intent(out) :: pbbfd(:)
      real(kind=r8), intent(out) :: pbbfu(:)
      real(kind=r8), intent(out) :: pbbfddir(:)
      real(kind=r8), intent(out) :: pbbcddir(:)

      real(kind=r8), intent(out) :: puvcd(:)
      real(kind=r8), intent(out) :: puvfd(:)
      real(kind=r8), intent(out) :: puvcddir(:)
      real(kind=r8), intent(out) :: puvfddir(:)

      real(kind=r8), intent(out) :: pnicd(:)
      real(kind=r8), intent(out) :: pnifd(:)
      real(kind=r8), intent(out) :: pnicddir(:)
      real(kind=r8), intent(out) :: pnifddir(:)

! Added for net near-IR flux diagnostic
      real(kind=r8), intent(out) :: pnicu(:)
      real(kind=r8), intent(out) :: pnifu(:)

! Output - inactive                                              !   All Dimensions: (nlayers+1)
!      real(kind=r8), intent(out) :: puvcu(:)
!      real(kind=r8), intent(out) :: puvfu(:)
!      real(kind=r8), intent(out) :: pnicu(:)
!      real(kind=r8), intent(out) :: pnifu(:)
!      real(kind=r8), intent(out) :: pvscd(:)
!      real(kind=r8), intent(out) :: pvscu(:)
!      real(kind=r8), intent(out) :: pvsfd(:)
!      real(kind=r8), intent(out) :: pvsfu(:)


! ------- Local -------

      logical :: lrtchkclr(nlayers),lrtchkcld(nlayers)

      integer  :: klev
      integer :: ib1, ib2, ibm, igt, ikl, ikp, ikx
      integer :: iw, jb, jg, jl, jk
!      integer, parameter :: nuv = ?? 
!      integer, parameter :: nvs = ?? 
      integer :: itind

      real(kind=r8) :: tblind, ze1
      real(kind=r8) :: zclear, zcloud
      real(kind=r8) :: zdbt(nlayers+1), zdbt_nodel(nlayers+1)
      real(kind=r8) :: zgc(nlayers), zgcc(nlayers), zgco(nlayers)
      real(kind=r8) :: zomc(nlayers), zomcc(nlayers), zomco(nlayers)
      real(kind=r8) :: zrdnd(nlayers+1), zrdndc(nlayers+1)
      real(kind=r8) :: zref(nlayers+1), zrefc(nlayers+1), zrefo(nlayers+1)
      real(kind=r8) :: zrefd(nlayers+1), zrefdc(nlayers+1), zrefdo(nlayers+1)
      real(kind=r8) :: zrup(nlayers+1), zrupd(nlayers+1)
      real(kind=r8) :: zrupc(nlayers+1), zrupdc(nlayers+1)
      real(kind=r8) :: zs1(nlayers+1)
      real(kind=r8) :: ztauc(nlayers), ztauo(nlayers)
      real(kind=r8) :: ztdn(nlayers+1), ztdnd(nlayers+1), ztdbt(nlayers+1)
      real(kind=r8) :: ztoc(nlayers), ztor(nlayers)
      real(kind=r8) :: ztra(nlayers+1), ztrac(nlayers+1), ztrao(nlayers+1)
      real(kind=r8) :: ztrad(nlayers+1), ztradc(nlayers+1), ztrado(nlayers+1)
      real(kind=r8) :: zdbtc(nlayers+1), ztdbtc(nlayers+1)
      real(kind=r8) :: zincflx(ngptsw), zdbtc_nodel(nlayers+1) 
      real(kind=r8) :: ztdbt_nodel(nlayers+1), ztdbtc_nodel(nlayers+1)

      real(kind=r8) :: zdbtmc, zdbtmo, zf, zgw, zreflect
      real(kind=r8) :: zwf, tauorig, repclc
!     real(kind=r8) :: zincflux                                   ! inactive

! Arrays from rrtmg_sw_taumoln routines

!      real(kind=r8) :: ztaug(nlayers,16), ztaur(nlayers,16)
!      real(kind=r8) :: zsflxzen(16)
      real(kind=r8) :: ztaug(ngptsw,nlayers), ztaur(ngptsw,nlayers)
      real(kind=r8) :: zsflxzen(ngptsw)

! Arrays from rrtmg_sw_vrtqdr routine

      real(kind=r8) :: zcd(ngptsw,nlayers+1), zcu(ngptsw,nlayers+1)
      real(kind=r8) :: zfd(ngptsw,nlayers+1), zfu(ngptsw,nlayers+1)

! Inactive arrays
!     real(kind=r8) :: zbbcd(nlayers+1), zbbcu(nlayers+1)
!     real(kind=r8) :: zbbfd(nlayers+1), zbbfu(nlayers+1)
!     real(kind=r8) :: zbbfddir(nlayers+1), zbbcddir(nlayers+1)

! Local storage arrays for zsflxzen, ztaug, and ztaur 

      real(kind=r8), allocatable, save :: ztaugst(:,:,:)      ! Gaseous optical depth storage array
      real(kind=r8), allocatable, save :: ztaurst(:,:,:)      ! Rayleigh optical depth storage array
      real(kind=r8), allocatable, save :: zsflxzenst(:,:)     ! Solar source function storage array

! Allocate storage arrays
      if (.not.allocated(ztaugst)) allocate (ztaugst(ngptsw,nlayers,ncol))
      if (.not.allocated(ztaurst)) allocate (ztaurst(ngptsw,nlayers,ncol))
      if (.not.allocated(zsflxzenst)) allocate (zsflxzenst(ngptsw,ncol))
       
! ------------------------------------------------------------------

! Initializations

      ib1 = istart
      ib2 = iend
      klev = nlayers
      iw = 0
      repclc = 1.e-12_r8
!      zincflux = 0.0_r8

      do jk=1,klev+1
         pbbcd(jk)=0._r8
         pbbcu(jk)=0._r8
         pbbfd(jk)=0._r8
         pbbfu(jk)=0._r8
         pbbcddir(jk)=0._r8
         pbbfddir(jk)=0._r8
         puvcd(jk)=0._r8
         puvfd(jk)=0._r8
         puvcddir(jk)=0._r8
         puvfddir(jk)=0._r8
         pnicd(jk)=0._r8
         pnifd(jk)=0._r8
         pnicddir(jk)=0._r8
         pnifddir(jk)=0._r8
         pnicu(jk)=0._r8
         pnifu(jk)=0._r8
      enddo


! Perform calculations of optical depth and solar source function at interval
! specified by dotau and store output arrays zsflxzen, ztaug, and ztaur in memory.
! In intervening time steps, obtain zsflxzen, ztaug, and ztaur from stored arrays. 
! Memory storage should be replaced with I/O to netCDF file.

         if (dotau) then
 
! Calculate the optical depths for gaseous absorption and Rayleigh scattering

            call taumol_sw(klev, &
                     colh2o, colco2, colch4, colo2, colo3, colmol, &
                     laytrop, jp, jt, jt1, &
                     fac00, fac01, fac10, fac11, &
                     selffac, selffrac, indself, forfac, forfrac, indfor, &
                     zsflxzen, ztaug, ztaur)

! Store ztaug, ztaur, and zsflxzen for use when dotau is false
            ztaugst(:ngptsw,:nlayers,iplon) = ztaug(:ngptsw,:nlayers)
            ztaurst(:ngptsw,:nlayers,iplon) = ztaur(:ngptsw,:nlayers)
            zsflxzenst(:ngptsw,iplon) = zsflxzen(:ngptsw)

         else

! Restore ztaug, ztaur, and zsflxzen from storage when dotau is false
            ztaug(:ngptsw,:nlayers) = ztaugst(:ngptsw,:nlayers,iplon)
            ztaur(:ngptsw,:nlayers) = ztaurst(:ngptsw,:nlayers,iplon)
            zsflxzen(:ngptsw) = zsflxzenst(:ngptsw,iplon)

         endif

! Top of shortwave spectral band loop, jb = 16 -> 29; ibm = 1 -> 14

      jb = ib1-1                  ! ???
      do jb = ib1, ib2
         ibm = jb-15
         igt = ngc(ibm)

! Reinitialize g-point counter for each band if output for each band is requested.
         if (iout.gt.0.and.ibm.ge.2) iw = ngs(ibm-1)

!        do jk=1,klev+1
!           zbbcd(jk)=0.0_r8
!           zbbcu(jk)=0.0_r8
!           zbbfd(jk)=0.0_r8
!           zbbfu(jk)=0.0_r8
!        enddo

! Top of g-point interval loop within each band (iw is cumulative counter) 
         do jg = 1,igt
            iw = iw+1

! Apply adjustments for correct Earth/Sun distance and zenith angle to incoming solar flux
            zincflx(iw) = adjflux(jb) * zsflxzen(iw) * prmu0
!             zincflux = zincflux + adjflux(jb) * zsflxzen(iw) * prmu0           ! inactive

! Compute layer reflectances and transmittances for direct and diffuse sources, 
! first clear then cloudy

! zrefc(jk)  direct albedo for clear
! zrefo(jk)  direct albedo for cloud
! zrefdc(jk) diffuse albedo for clear
! zrefdo(jk) diffuse albedo for cloud
! ztrac(jk)  direct transmittance for clear
! ztrao(jk)  direct transmittance for cloudy
! ztradc(jk) diffuse transmittance for clear
! ztrado(jk) diffuse transmittance for cloudy
!  
! zref(jk)   direct reflectance
! zrefd(jk)  diffuse reflectance
! ztra(jk)   direct transmittance
! ztrad(jk)  diffuse transmittance
!
! zdbtc(jk)  clear direct beam transmittance
! zdbto(jk)  cloudy direct beam transmittance
! zdbt(jk)   layer mean direct beam transmittance
! ztdbt(jk)  total direct beam transmittance at levels

! Clear-sky    
!   TOA direct beam    
            ztdbtc(1)=1.0_r8
            ztdbtc_nodel(1)=1.0_r8
!   Surface values
            zdbtc(klev+1) =0.0_r8
            ztrac(klev+1) =0.0_r8
            ztradc(klev+1)=0.0_r8
            zrefc(klev+1) =palbp(ibm)
            zrefdc(klev+1)=palbd(ibm)
            zrupc(klev+1) =palbp(ibm)
            zrupdc(klev+1)=palbd(ibm)
           
! Total sky    
!   TOA direct beam    
            ztdbt(1)=1.0_r8
            ztdbt_nodel(1)=1.0_r8
!   Surface values
            zdbt(klev+1) =0.0_r8
            ztra(klev+1) =0.0_r8
            ztrad(klev+1)=0.0_r8
            zref(klev+1) =palbp(ibm)
            zrefd(klev+1)=palbd(ibm)
            zrup(klev+1) =palbp(ibm)
            zrupd(klev+1)=palbd(ibm)
    
    
! Top of layer loop
            do jk=1,klev

! Note: two-stream calculations proceed from top to bottom; 
!   RRTMG_SW quantities are given bottom to top and are reversed here

               ikl=klev+1-jk

! Set logical flag to do REFTRA calculation
!   Do REFTRA for all clear layers
               lrtchkclr(jk)=.true.

!   Do REFTRA only for cloudy layers in profile, since already done for clear layers
               lrtchkcld(jk)=.false.
               lrtchkcld(jk)=(pclfr(ikl) > repclc)

! Clear-sky optical parameters - this section inactive     
!   Original
!               ztauc(jk) = ztaur(iw,ikl) + ztaug(iw,ikl)
!               zomcc(jk) = ztaur(iw,ikl) / ztauc(jk)
!               zgcc(jk) = 0.0001_r8
!   Total sky optical parameters        
!               ztauo(jk) = ztaur(iw,ikl) + ztaug(iw,ikl) + ptauc(ibm,ikl)
!               zomco(jk) = ptauc(ibm,ikl) * pomgc(ibm,ikl) + ztaur(iw,ikl)
!               zgco (jk) = (ptauc(ibm,ikl) * pomgc(ibm,ikl) * pasyc(ibm,ikl) + &
!                           ztaur(iw,ikl) * 0.0001_r8) / zomco(jk)
!               zomco(jk) = zomco(jk) / ztauo(jk)

! Clear-sky optical parameters including aerosols
               ztauc(jk) = ztaur(iw,ikl) + ztaug(iw,ikl) + ptaua(ibm,ikl)
               zomcc(jk) = ztaur(iw,ikl) * 1.0_r8 + ptaua(ibm,ikl) * pomga(ibm,ikl)
               zgcc(jk) = pasya(ibm,ikl) * pomga(ibm,ikl) * ptaua(ibm,ikl) / zomcc(jk)
               zomcc(jk) = zomcc(jk) / ztauc(jk)

! Pre-delta-scaling clear and cloudy direct beam transmittance (must use 'orig', unscaled cloud OD)       
!   \/\/\/ This block of code is only needed for unscaled direct beam calculation
               if (idelm .eq. 0) then
!     
                  zclear = 1.0_r8 - pclfr(ikl)
                  zcloud = pclfr(ikl)

! Clear
!                   zdbtmc = exp(-ztauc(jk) / prmu0)
 
! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ze1 = ztauc(jk) / prmu0
                  if (ze1 .le. od_lo) then
                     zdbtmc = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
                  else 
                     tblind = ze1 / (bpade + ze1)
                     itind = tblint * tblind + 0.5_r8
                     zdbtmc = exp_tbl(itind)
                  endif

                  zdbtc_nodel(jk) = zdbtmc
                  ztdbtc_nodel(jk+1) = zdbtc_nodel(jk) * ztdbtc_nodel(jk)

! Clear + Cloud
                  tauorig = ztauc(jk) + ptaucorig(ibm,ikl)
!                   zdbtmo = exp(-tauorig / prmu0)

! Use exponential lookup table for transmittance, or expansion of exponential for low tau
                  ze1 = tauorig / prmu0
                  if (ze1 .le. od_lo) then
                     zdbtmo = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
                  else
                     tblind = ze1 / (bpade + ze1)
                     itind = tblint * tblind + 0.5_r8
                     zdbtmo = exp_tbl(itind)
                  endif

                  zdbt_nodel(jk) = zclear * zdbtmc + zcloud * zdbtmo
                  ztdbt_nodel(jk+1) = zdbt_nodel(jk) * ztdbt_nodel(jk)

               endif
!   /\/\/\ Above code only needed for unscaled direct beam calculation


! Delta scaling - clear   
               zf = zgcc(jk) * zgcc(jk)
               zwf = zomcc(jk) * zf
               ztauc(jk) = (1.0_r8 - zwf) * ztauc(jk)
               zomcc(jk) = (zomcc(jk) - zwf) / (1.0_r8 - zwf)
               zgcc (jk) = (zgcc(jk) - zf) / (1.0_r8 - zf)

! Total sky optical parameters (cloud properties already delta-scaled)
!   Use this code if cloud properties are derived in rrtmg_sw_cldprop       
               if (icpr .ge. 1) then
                  ztauo(jk) = ztauc(jk) + ptauc(ibm,ikl)
                  zomco(jk) = ztauc(jk) * zomcc(jk) + ptauc(ibm,ikl) * pomgc(ibm,ikl) 
                  zgco (jk) = (ptauc(ibm,ikl) * pomgc(ibm,ikl) * pasyc(ibm,ikl) + &
                              ztauc(jk) * zomcc(jk) * zgcc(jk)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Total sky optical parameters (if cloud properties not delta scaled)
!   Use this code if cloud properties are not derived in rrtmg_sw_cldprop       
               elseif (icpr .eq. 0) then
                  ztauo(jk) = ztaur(iw,ikl) + ztaug(iw,ikl) + ptaua(ibm,ikl) + ptauc(ibm,ikl)
                  zomco(jk) = ptaua(ibm,ikl) * pomga(ibm,ikl) + ptauc(ibm,ikl) * pomgc(ibm,ikl) + &
                              ztaur(iw,ikl) * 1.0_r8
                  zgco (jk) = (ptauc(ibm,ikl) * pomgc(ibm,ikl) * pasyc(ibm,ikl) + &
                              ptaua(ibm,ikl)*pomga(ibm,ikl)*pasya(ibm,ikl)) / zomco(jk)
                  zomco(jk) = zomco(jk) / ztauo(jk)

! Delta scaling - clouds 
!   Use only if subroutine rrtmg_sw_cldprop is not used to get cloud properties and to apply delta scaling
                  zf = zgco(jk) * zgco(jk)
                  zwf = zomco(jk) * zf
                  ztauo(jk) = (1._r8 - zwf) * ztauo(jk)
                  zomco(jk) = (zomco(jk) - zwf) / (1.0_r8 - zwf)
                  zgco (jk) = (zgco(jk) - zf) / (1.0_r8 - zf)
               endif 

! End of layer loop
            enddo    


! Clear sky reflectivities
            call reftra_sw (klev, &
                            lrtchkclr, zgcc, prmu0, ztauc, zomcc, &
                            zrefc, zrefdc, ztrac, ztradc)

! Total sky reflectivities      
            call reftra_sw (klev, &
                            lrtchkcld, zgco, prmu0, ztauo, zomco, &
                            zrefo, zrefdo, ztrao, ztrado)


            do jk=1,klev

! Combine clear and cloudy contributions for total sky
               ikl = klev+1-jk 
               zclear = 1.0_r8 - pclfr(ikl)
               zcloud = pclfr(ikl)

               zref(jk) = zclear*zrefc(jk) + zcloud*zrefo(jk)
               zrefd(jk)= zclear*zrefdc(jk) + zcloud*zrefdo(jk)
               ztra(jk) = zclear*ztrac(jk) + zcloud*ztrao(jk)
               ztrad(jk)= zclear*ztradc(jk) + zcloud*ztrado(jk)

! Direct beam transmittance        

! Clear
!                zdbtmc = exp(-ztauc(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauc(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmc = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  zdbtmc = exp_tbl(itind)
               endif

               zdbtc(jk) = zdbtmc
               ztdbtc(jk+1) = zdbtc(jk)*ztdbtc(jk)

! Clear + Cloud
!                zdbtmo = exp(-ztauo(jk) / prmu0)

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               ze1 = ztauo(jk) / prmu0
               if (ze1 .le. od_lo) then
                  zdbtmo = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  zdbtmo = exp_tbl(itind)
               endif

               zdbt(jk) = zclear*zdbtmc + zcloud*zdbtmo
               ztdbt(jk+1) = zdbt(jk)*ztdbt(jk)
        
            enddo           
                 
! Vertical quadrature for clear-sky fluxes

            call vrtqdr_sw (klev, iw, &
                            zrefc, zrefdc, ztrac, ztradc, &
                            zdbtc, zrdndc, zrupc, zrupdc, ztdbtc, &
                            zcd, zcu)
      
! Vertical quadrature for cloudy fluxes

            call vrtqdr_sw (klev, iw, &
                            zref, zrefd, ztra, ztrad, &
                            zdbt, zrdnd, zrup, zrupd, ztdbt, &
                            zfd, zfu)

! Upwelling and downwelling fluxes at levels
!   Two-stream calculations go from top to bottom; 
!   layer indexing is reversed to go bottom to top for output arrays

            do jk=1,klev+1
               ikl=klev+2-jk

! Accumulate spectral fluxes over bands - inactive
!               zbbfu(ikl) = zbbfu(ikl) + zincflx(iw)*zfu(iw,jk)  
!               zbbfd(ikl) = zbbfd(ikl) + zincflx(iw)*zfd(iw,jk)
!               zbbcu(ikl) = zbbcu(ikl) + zincflx(iw)*zcu(iw,jk)
!               zbbcd(ikl) = zbbcd(ikl) + zincflx(iw)*zcd(iw,jk)
!               zbbfddir(ikl) = zbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
!               zbbcddir(ikl) = zbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)

! Accumulate spectral fluxes over whole spectrum  
               pbbfu(ikl) = pbbfu(ikl) + zincflx(iw)*zfu(iw,jk)
               pbbfd(ikl) = pbbfd(ikl) + zincflx(iw)*zfd(iw,jk)
               pbbcu(ikl) = pbbcu(ikl) + zincflx(iw)*zcu(iw,jk)
               pbbcd(ikl) = pbbcd(ikl) + zincflx(iw)*zcd(iw,jk)
               if (idelm .eq. 0) then
                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
               elseif (idelm .eq. 1) then
                  pbbfddir(ikl) = pbbfddir(ikl) + zincflx(iw)*ztdbt(jk)
                  pbbcddir(ikl) = pbbcddir(ikl) + zincflx(iw)*ztdbtc(jk)
               endif

! Accumulate direct fluxes for UV/visible bands
               if (ibm >= 10 .and. ibm <= 13) then
                  puvcd(ikl) = puvcd(ikl) + zincflx(iw)*zcd(iw,jk)
                  puvfd(ikl) = puvfd(ikl) + zincflx(iw)*zfd(iw,jk)
                  if (idelm .eq. 0) then
                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  elseif (idelm .eq. 1) then
                     puvfddir(ikl) = puvfddir(ikl) + zincflx(iw)*ztdbt(jk)
                     puvcddir(ikl) = puvcddir(ikl) + zincflx(iw)*ztdbtc(jk)
                  endif
! Accumulate direct fluxes for near-IR bands
               else if (ibm == 14 .or. ibm <= 9) then  
                  pnicd(ikl) = pnicd(ikl) + zincflx(iw)*zcd(iw,jk)
                  pnifd(ikl) = pnifd(ikl) + zincflx(iw)*zfd(iw,jk)
                  if (idelm .eq. 0) then
                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt_nodel(jk)
                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc_nodel(jk)
                  elseif (idelm .eq. 1) then
                     pnifddir(ikl) = pnifddir(ikl) + zincflx(iw)*ztdbt(jk)
                     pnicddir(ikl) = pnicddir(ikl) + zincflx(iw)*ztdbtc(jk)
                  endif
! Added for net near-IR flux diagnostic 
                  pnicu(ikl) = pnicu(ikl) + zincflx(iw)*zcu(iw,jk)
                  pnifu(ikl) = pnifu(ikl) + zincflx(iw)*zfu(iw,jk)
               endif

            enddo

! End loop on jg, g-point interval
         enddo             

! End loop on jb, spectral band
      enddo                    

      end subroutine spcvrt_sw

      end module rrtmg_sw_spcvrt


