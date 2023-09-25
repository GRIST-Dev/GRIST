!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_init.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/22 19:20:03 $
!
      module rrtmg_lw_init

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2007, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

      use grist_constants,  only: r8

!      use parkind, only : jpim, jprb 
      use rrlw_wvn
      use rrtmg_lw_setcoef, only: lwatmref, lwavplank

      implicit none

      contains

! **************************************************************************
      subroutine rrtmg_lw_ini
! **************************************************************************
!
!  Original version:       Michael J. Iacono; July, 1998
!  First revision for NCAR CCM:   September, 1998
!  Second revision for RRTM_V3.0:  September, 2002
!
!  This subroutine performs calculations necessary for the initialization
!  of the longwave model.  Lookup tables are computed for use in the LW
!  radiative transfer, and input absorption coefficient data for each
!  spectral band are reduced from 256 g-point intervals to 140.
! **************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw
      use rrlw_tbl, only: ntbl, tblint, pade, bpade, tau_tbl, exp_tbl, tfn_tbl
      use rrlw_vsn, only: hvrini, hnamini

! ------- Local -------

      integer :: itr, ibnd, igc, ig, ind, ipr 
      integer :: igcsm, iprsm

      real(kind=r8) :: wtsum, wtsm(mg)        !
      real(kind=r8) :: tfn                    !

! ------- Definitions -------
!     Arrays for 10000-point look-up tables:
!     TAU_TBL Clear-sky optical depth (used in cloudy radiative transfer)
!     EXP_TBL Exponential lookup table for ransmittance
!     TFN_TBL Tau transition function; i.e. the transition of the Planck
!             function from that for the mean layer temperature to that for
!             the layer boundary temperature as a function of optical depth.
!             The "linear in tau" method is used to make the table.
!     PADE    Pade approximation constant (= 0.278)
!     BPADE   Inverse of the Pade approximation constant
!

      hvrini = '$Revision: 1.2 $'

! Initialize model data
      call lwdatinit
      call lwcmbdat               ! g-point interval reduction data
      call lwcldpr                ! cloud optical properties
      call lwatmref               ! reference MLS profile
      call lwavplank              ! Planck function 
      call lw_kgb01               ! molecular absorption coefficients
      call lw_kgb02
      call lw_kgb03
      call lw_kgb04
      call lw_kgb05
      call lw_kgb06
      call lw_kgb07
      call lw_kgb08
      call lw_kgb09
      call lw_kgb10
      call lw_kgb11
      call lw_kgb12
      call lw_kgb13
      call lw_kgb14
      call lw_kgb15
      call lw_kgb16

! Compute lookup tables for transmittance, tau transition function,
! and clear sky tau (for the cloudy sky radiative transfer).  Tau is 
! computed as a function of the tau transition function, transmittance 
! is calculated as a function of tau, and the tau transition function 
! is calculated using the linear in tau formulation at values of tau 
! above 0.01.  TF is approximated as tau/6 for tau < 0.01.  All tables 
! are computed at intervals of 0.001.  The inverse of the constant used
! in the Pade approximation to the tau transition function is set to b.

      tau_tbl(0) = 0.0_r8
      tau_tbl(ntbl) = 1.e10_r8
      exp_tbl(0) = 1.0_r8
      exp_tbl(ntbl) = 0.0_r8
      tfn_tbl(0) = 0.0_r8
      tfn_tbl(ntbl) = 1.0_r8
      bpade = 1.0_r8 / pade
      do itr = 1, ntbl-1
         tfn = float(itr) / float(ntbl)
         tau_tbl(itr) = bpade * tfn / (1._r8 - tfn)
         exp_tbl(itr) = exp(-tau_tbl(itr))
         if (tau_tbl(itr) .lt. 0.06_r8) then
            tfn_tbl(itr) = tau_tbl(itr)/6._r8
         else
            tfn_tbl(itr) = 1._r8-2._r8*((1._r8/tau_tbl(itr))-(exp_tbl(itr)/(1.-exp_tbl(itr))))
         endif
      enddo

! Perform g-point reduction from 16 per band (256 total points) to
! a band dependant number (140 total points) for all absorption
! coefficient input data and Planck fraction input data.
! Compute relative weighting for new g-point combinations.

      igcsm = 0
      do ibnd = 1,nbndlw
         iprsm = 0
         if (ngc(ibnd).lt.mg) then
            do igc = 1,ngc(ibnd) 
               igcsm = igcsm + 1
               wtsum = 0._r8
               do ipr = 1, ngn(igcsm)
                  iprsm = iprsm + 1
                  wtsum = wtsum + wt(iprsm)
               enddo
               wtsm(igc) = wtsum
            enddo
            do ig = 1, ng(ibnd)
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = wt(ig)/wtsm(ngm(ind))
            enddo
         else
            do ig = 1, ng(ibnd)
               igcsm = igcsm + 1
               ind = (ibnd-1)*mg + ig
               rwgt(ind) = 1.0_r8
            enddo
         endif
      enddo

! Reduce g-points for absorption coefficient data in each LW spectral band.

      call cmbgb1
      call cmbgb2
      call cmbgb3
      call cmbgb4
      call cmbgb5
      call cmbgb6
      call cmbgb7
      call cmbgb8
      call cmbgb9
      call cmbgb10
      call cmbgb11
      call cmbgb12
      call cmbgb13
      call cmbgb14
      call cmbgb15
      call cmbgb16

      end subroutine rrtmg_lw_ini

!***************************************************************************
      subroutine lwdatinit
!***************************************************************************

! --------- Modules ----------

      use parrrtm, only : maxxsec, maxinpx
      use rrlw_con, only: heatfac, grav, planck, boltz, &
                          clight, avogad, alosmt, gascon, radcn1, radcn2 
      use grist_constants, only: shr_const_avogad=>avogad, cday, gravit=>gravity, cpair=>cp

      use rrlw_vsn

      save 
 
! Longwave spectral band limits (wavenumbers)
      wavenum1(:) = (/ 10._r8, 350._r8, 500._r8, 630._r8, 700._r8, 820._r8, &
                      980._r8,1080._r8,1180._r8,1390._r8,1480._r8,1800._r8, &
                     2080._r8,2250._r8,2390._r8,2600._r8/)
      wavenum2(:) = (/350._r8, 500._r8, 630._r8, 700._r8, 820._r8, 980._r8, &
                     1080._r8,1180._r8,1390._r8,1480._r8,1800._r8,2080._r8, &
                     2250._r8,2390._r8,2600._r8,3250._r8/)
      delwave(:) =  (/340._r8, 150._r8, 130._r8,  70._r8, 120._r8, 160._r8, &
                      100._r8, 100._r8, 210._r8,  90._r8, 320._r8, 280._r8, &
                      170._r8, 130._r8, 220._r8, 650._r8/)

! Spectral band information
      ng(:) = (/16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16/)
      nspa(:) = (/1,1,9,9,9,1,9,1,9,1,1,9,9,1,9,9/)
      nspb(:) = (/1,1,5,5,5,0,1,1,1,1,1,0,0,1,0,0/)

! Use constants set in CAM for consistency
      grav = gravit
      avogad = shr_const_avogad * 1.e-3_r8

!     Heatfac is the factor by which one must multiply delta-flux/ 
!     delta-pressure, with flux in w/m-2 and pressure in mbar, to get 
!     the heating rate in units of degrees/day.  It is equal to 
!           (g)x(#sec/day)x(1e-5)/(specific heat of air at const. p)
!        =  (9.8066)(86400)(1e-5)/(1.004)
!      heatfac = 8.4391_r8

!     Modified values for consistency with CAM:
!        =  (9.80616)(86400)(1e-5)/(1.00464)
!      heatfac = 8.43339130434_r8

!     Calculate heatfac directly from CAM constants:
      heatfac = grav * cday * 1.e-5_r8 / (cpair * 1.e-3_r8)

!     nxmol     - number of cross-sections input by user
!     ixindx(i) - index of cross-section molecule corresponding to Ith
!                 cross-section specified by user
!                 = 0 -- not allowed in rrtm
!                 = 1 -- ccl4
!                 = 2 -- cfc11
!                 = 3 -- cfc12
!                 = 4 -- cfc22
      nxmol = 4
      ixindx(1) = 1
      ixindx(2) = 2
      ixindx(3) = 3
      ixindx(4) = 4
      ixindx(5:maxinpx) = 0

!    Constants from NIST 01/11/2002

!      grav = 9.8066_r8
      planck = 6.62606876e-27_r8
      boltz = 1.3806503e-16_r8
      clight = 2.99792458e+10_r8
!      avogad = 6.02214199e+23_r8
      alosmt = 2.6867775e+19_r8
      gascon = 8.31447200e+07_r8
      radcn1 = 1.191042722e-12_r8
      radcn2 = 1.4387752_r8

!
!     units are generally cgs
!
!     The first and second radiation constants are taken from NIST.
!     They were previously obtained from the relations:
!          radcn1 = 2.*planck*clight*clight*1.e-07
!          radcn2 = planck*clight/boltz

      end subroutine lwdatinit

!***************************************************************************
      subroutine lwcmbdat
!***************************************************************************

      save
 
! ------- Definitions -------
!     Arrays for the g-point reduction from 256 to 140 for the 16 LW bands:
!     This mapping from 256 to 140 points has been carefully selected to 
!     minimize the effect on the resulting fluxes and cooling rates, and
!     caution should be used if the mapping is modified.  The full 256
!     g-point set can be restored with ngptlw=256, ngc=16*16, ngn=256*1., etc.
!     ngptlw  The total number of new g-points
!     ngc     The number of new g-points in each band
!     ngs     The cumulative sum of new g-points for each band
!     ngm     The index of each new g-point relative to the original
!             16 g-points for each band.  
!     ngn     The number of original g-points that are combined to make
!             each new g-point in each band.
!     ngb     The band index for each new g-point.
!     wt      RRTM weights for 16 g-points.

! ------- Data statements -------
      ngc(:) = (/10,12,16,14,16,8,12,8,12,6,8,8,4,2,2,2/)
      ngs(:) = (/10,22,38,52,68,76,88,96,108,114,122,130,134,136,138,140/)
      ngm(:) = (/1,2,3,3,4,4,5,5,6,6,7,7,8,8,9,10, &          ! band 1
                 1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 2
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 3
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,14,14, &    ! band 4
                 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, &    ! band 5
                 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 6
                 1,1,2,2,3,4,5,6,7,8,9,10,11,11,12,12, &      ! band 7
                 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8, &           ! band 8
                 1,2,3,4,5,6,7,8,9,9,10,10,11,11,12,12, &     ! band 9
                 1,1,2,2,3,3,4,4,5,5,5,5,6,6,6,6, &           ! band 10
                 1,2,3,3,4,4,5,5,6,6,7,7,7,8,8,8, &           ! band 11
                 1,2,3,4,5,5,6,6,7,7,7,7,8,8,8,8, &           ! band 12
                 1,1,1,2,2,2,3,3,3,3,4,4,4,4,4,4, &           ! band 13
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 14
                 1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2, &           ! band 15
                 1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2/)            ! band 16
      ngn(:) = (/1,1,2,2,2,2,2,2,1,1, &                       ! band 1
                 1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 2
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 3
                 1,1,1,1,1,1,1,1,1,1,1,1,1,3, &               ! band 4
                 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1, &           ! band 5
                 2,2,2,2,2,2,2,2, &                           ! band 6
                 2,2,1,1,1,1,1,1,1,1,2,2, &                   ! band 7
                 2,2,2,2,2,2,2,2, &                           ! band 8
                 1,1,1,1,1,1,1,1,2,2,2,2, &                   ! band 9
                 2,2,2,2,4,4, &                               ! band 10
                 1,1,2,2,2,2,3,3, &                           ! band 11
                 1,1,1,1,2,2,4,4, &                           ! band 12
                 3,3,4,6, &                                   ! band 13
                 8,8, &                                       ! band 14
                 8,8, &                                       ! band 15
                 4,12/)                                       ! band 16
      ngb(:) = (/1,1,1,1,1,1,1,1,1,1, &                       ! band 1
                 2,2,2,2,2,2,2,2,2,2,2,2, &                   ! band 2
                 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3, &           ! band 3
                 4,4,4,4,4,4,4,4,4,4,4,4,4,4, &               ! band 4
                 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5, &           ! band 5
                 6,6,6,6,6,6,6,6, &                           ! band 6
                 7,7,7,7,7,7,7,7,7,7,7,7, &                   ! band 7
                 8,8,8,8,8,8,8,8, &                           ! band 8
                 9,9,9,9,9,9,9,9,9,9,9,9, &                   ! band 9
                 10,10,10,10,10,10, &                         ! band 10
                 11,11,11,11,11,11,11,11, &                   ! band 11
                 12,12,12,12,12,12,12,12, &                   ! band 12
                 13,13,13,13, &                               ! band 13
                 14,14, &                                     ! band 14
                 15,15, &                                     ! band 15
                 16,16/)                                      ! band 16
      wt(:) = (/ 0.1527534276_r8, 0.1491729617_r8, 0.1420961469_r8, &
                 0.1316886544_r8, 0.1181945205_r8, 0.1019300893_r8, &
                 0.0832767040_r8, 0.0626720116_r8, 0.0424925000_r8, &
                 0.0046269894_r8, 0.0038279891_r8, 0.0030260086_r8, &
                 0.0022199750_r8, 0.0014140010_r8, 0.0005330000_r8, &
                 0.0000750000_r8/)

      end subroutine lwcmbdat

!***************************************************************************
      subroutine cmbgb1
!***************************************************************************
!
!  Original version:    MJIacono; July 1998
!  Revision for GCMs:   MJIacono; September 1998
!  Revision for RRTMG:  MJIacono, September 2002
!  Revision for F90 reformatting:  MJIacono, June 2006
!
!  The subroutines CMBGB1->CMBGB16 input the absorption coefficient
!  data for each band, which are defined for 16 g-points and 16 spectral
!  bands. The data are combined with appropriate weighting following the
!  g-point mapping arrays specified in RRTMINIT.  Plank fraction data
!  in arrays FRACREFA and FRACREFB are combined without weighting.  All
!  g-point reduced data are put into new arrays for use in RRTM.
!
!  band 1:  10-350 cm-1 (low key - h2o; low minor - n2)
!                       (high key - h2o; high minor - n2)
!  note: previous versions of rrtm band 1: 
!        10-250 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng1
      use rrlw_kg01, only: fracrefao, fracrefbo, kao, kbo, kao_mn2, kbo_mn2, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mn2, kb_mn2, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(1)
               sumk = 0.
               do ipr = 1, ngn(igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(1)
               sumk = 0.
               do ipr = 1, ngn(igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(1)
            sumk = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(1)
            sumk = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(1)
            sumk1 = 0.
            sumk2 = 0.
            do ipr = 1, ngn(igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mn2(jt,iprsm)*rwgt(iprsm)
               sumk2 = sumk2 + kbo_mn2(jt,iprsm)*rwgt(iprsm)
            enddo
            ka_mn2(jt,igc) = sumk1
            kb_mn2(jt,igc) = sumk2
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(1)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb1

!***************************************************************************
      subroutine cmbgb2
!***************************************************************************
!
!     band 2:  350-500 cm-1 (low key - h2o; high key - h2o)
!
!     note: previous version of rrtm band 2: 
!           250 - 500 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng2
      use rrlw_kg02, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(2)
               sumk = 0.
               do ipr = 1, ngn(ngs(1)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+16)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(2)
               sumk = 0.
               do ipr = 1, ngn(ngs(1)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+16)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(2)
            sumk = 0.
            do ipr = 1, ngn(ngs(1)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+16)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(2)
            sumk = 0.
            do ipr = 1, ngn(ngs(1)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+16)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(2)
         sumf1 = 0.
         sumf2 = 0.
         do ipr = 1, ngn(ngs(1)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb2

!***************************************************************************
      subroutine cmbgb3
!***************************************************************************
!
!     band 3:  500-630 cm-1 (low key - h2o,co2; low minor - n2o)
!                           (high key - h2o,co2; high minor - n2o)
!
! old band 3:  500-630 cm-1 (low - h2o,co2; high - h2o,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng3
      use rrlw_kg03, only: fracrefao, fracrefbo, kao, kbo, kao_mn2o, kbo_mn2o, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mn2o, kb_mn2o, &
                           selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(3)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(2)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+32)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(3)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(2)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+32)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(3)
              sumk = 0.
               do ipr = 1, ngn(ngs(2)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
               enddo
               ka_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jn = 1,5
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(3)
              sumk = 0.
               do ipr = 1, ngn(ngs(2)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo_mn2o(jn,jt,iprsm)*rwgt(iprsm+32)
               enddo
               kb_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(3)
            sumk = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+32)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(3)
            sumk = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+32)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(3)
            sumf = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(3)
            sumf = 0.
            do ipr = 1, ngn(ngs(2)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb3

!***************************************************************************
      subroutine cmbgb4
!***************************************************************************
!
!     band 4:  630-700 cm-1 (low key - h2o,co2; high key - o3,co2)
!
! old band 4:  630-700 cm-1 (low - h2o,co2; high - o3,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng4
      use rrlw_kg04, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(4)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(3)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+48)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(4)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(3)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+48)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(4)
            sumk = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+48)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(4)
            sumk = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+48)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(4)
            sumf = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(4)
            sumf = 0.
            do ipr = 1, ngn(ngs(3)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb4

!***************************************************************************
      subroutine cmbgb5
!***************************************************************************
!
!     band 5:  700-820 cm-1 (low key - h2o,co2; low minor - o3, ccl4)
!                           (high key - o3,co2)
!
! old band 5:  700-820 cm-1 (low - h2o,co2; high - o3,co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng5
      use rrlw_kg05, only: fracrefao, fracrefbo, kao, kbo, kao_mo3, ccl4o, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mo3, ccl4, &
                           selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(5)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(4)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+64)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jn = 1,5
         do jt = 1,5
            do jp = 13,59
               iprsm = 0
               do igc = 1,ngc(5)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(4)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kbo(jn,jt,jp,iprsm)*rwgt(iprsm+64)
                  enddo
                  kb(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(5)
              sumk = 0.
               do ipr = 1, ngn(ngs(4)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mo3(jn,jt,iprsm)*rwgt(iprsm+64)
               enddo
               ka_mo3(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(5)
            sumk = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+64)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(5)
            sumk = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+64)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(5)
            sumf = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      do jp = 1,5
         iprsm = 0
         do igc = 1,ngc(5)
            sumf = 0.
            do ipr = 1, ngn(ngs(4)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefbo(iprsm,jp)
            enddo
            fracrefb(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(5)
         sumk = 0.
         do ipr = 1, ngn(ngs(4)+igc)
            iprsm = iprsm + 1
            sumk = sumk + ccl4o(iprsm)*rwgt(iprsm+64)
         enddo
         ccl4(igc) = sumk
      enddo

      end subroutine cmbgb5

!***************************************************************************
      subroutine cmbgb6
!***************************************************************************
!
!     band 6:  820-980 cm-1 (low key - h2o; low minor - co2)
!                           (high key - nothing; high minor - cfc11, cfc12)
!
! old band 6:  820-980 cm-1 (low - h2o; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng6
      use rrlw_kg06, only: fracrefao, kao, kao_mco2, cfc11adjo, cfc12o, &
                           selfrefo, forrefo, &
                           fracrefa, ka, ka_mco2, cfc11adj, cfc12, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf, sumk1, sumk2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(6)
               sumk = 0.
               do ipr = 1, ngn(ngs(5)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+80)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kao_mco2(jt,iprsm)*rwgt(iprsm+80)
            enddo
            ka_mco2(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+80)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(6)
            sumk = 0.
            do ipr = 1, ngn(ngs(5)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+80)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(6)
         sumf = 0.
         sumk1= 0.
         sumk2= 0.
         do ipr = 1, ngn(ngs(5)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefao(iprsm)
            sumk1= sumk1+ cfc11adjo(iprsm)*rwgt(iprsm+80)
            sumk2= sumk2+ cfc12o(iprsm)*rwgt(iprsm+80)
         enddo
         fracrefa(igc) = sumf
         cfc11adj(igc) = sumk1
         cfc12(igc) = sumk2
      enddo

      end subroutine cmbgb6

!***************************************************************************
      subroutine cmbgb7
!***************************************************************************
!
!     band 7:  980-1080 cm-1 (low key - h2o,o3; low minor - co2)
!                            (high key - o3; high minor - co2)
!
! old band 7:  980-1080 cm-1 (low - h2o,o3; high - o3)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng7
      use rrlw_kg07, only: fracrefao, fracrefbo, kao, kbo, kao_mco2, kbo_mco2, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, ka_mco2, kb_mco2, &
                           selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(7)
                 sumk = 0.
                  do ipr = 1, ngn(ngs(6)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+96)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(7)
               sumk = 0.
               do ipr = 1, ngn(ngs(6)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+96)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(7)
              sumk = 0.
               do ipr = 1, ngn(ngs(6)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+96)
               enddo
               ka_mco2(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mco2(jt,iprsm)*rwgt(iprsm+96)
            enddo
            kb_mco2(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+96)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(7)
            sumk = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+96)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(7)
            sumf = 0.
            do ipr = 1, ngn(ngs(6)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(7)
         sumf = 0.
         do ipr = 1, ngn(ngs(6)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      end subroutine cmbgb7

!***************************************************************************
      subroutine cmbgb8
!***************************************************************************
!
!     band 8:  1080-1180 cm-1 (low key - h2o; low minor - co2,o3,n2o)
!                             (high key - o3; high minor - co2, n2o)
!
! old band 8:  1080-1180 cm-1 (low (i.e.>~300mb) - h2o; high - o3)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng8
      use rrlw_kg08, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mn2o, &
                           kao_mo3, kbo, kbo_mco2, kbo_mn2o, selfrefo, forrefo, &
                           cfc12o, cfc22adjo, &
                           fracrefa, fracrefb, ka, ka_mco2, ka_mn2o, &
                           ka_mo3, kb, kb_mco2, kb_mn2o, selfref, forref, &
                           cfc12, cfc22adj

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumk3, sumk4, sumk5, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(8)
              sumk = 0.
               do ipr = 1, ngn(ngs(7)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+112)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(8)
               sumk = 0.
               do ipr = 1, ngn(ngs(7)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+112)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(8)
            sumk = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+112)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(8)
            sumk = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+112)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(8)
            sumk1 = 0.
            sumk2 = 0.
            sumk3 = 0.
            sumk4 = 0.
            sumk5 = 0.
            do ipr = 1, ngn(ngs(7)+igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mco2(jt,iprsm)*rwgt(iprsm+112)
               sumk2 = sumk2 + kbo_mco2(jt,iprsm)*rwgt(iprsm+112)
               sumk3 = sumk3 + kao_mo3(jt,iprsm)*rwgt(iprsm+112)
               sumk4 = sumk4 + kao_mn2o(jt,iprsm)*rwgt(iprsm+112)
               sumk5 = sumk5 + kbo_mn2o(jt,iprsm)*rwgt(iprsm+112)
            enddo
            ka_mco2(jt,igc) = sumk1
            kb_mco2(jt,igc) = sumk2
            ka_mo3(jt,igc) = sumk3
            ka_mn2o(jt,igc) = sumk4
            kb_mn2o(jt,igc) = sumk5
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(8)
         sumf1= 0.
         sumf2= 0.
         sumk1= 0.
         sumk2= 0.
         do ipr = 1, ngn(ngs(7)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
            sumk1= sumk1+ cfc12o(iprsm)*rwgt(iprsm+112)
            sumk2= sumk2+ cfc22adjo(iprsm)*rwgt(iprsm+112)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
         cfc12(igc) = sumk1
         cfc22adj(igc) = sumk2
      enddo

      end subroutine cmbgb8

!***************************************************************************
      subroutine cmbgb9
!***************************************************************************
!
!     band 9:  1180-1390 cm-1 (low key - h2o,ch4; low minor - n2o)
!                             (high key - ch4; high minor - n2o)!

! old band 9:  1180-1390 cm-1 (low - h2o,ch4; high - ch4)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng9
      use rrlw_kg09, only: fracrefao, fracrefbo, kao, kao_mn2o, &
                           kbo, kbo_mn2o, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, ka_mn2o, &
                           kb, kb_mn2o, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(9)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(8)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+128)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(9)
               sumk = 0.
               do ipr = 1, ngn(ngs(8)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+128)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(9)
              sumk = 0.
               do ipr = 1, ngn(ngs(8)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2o(jn,jt,iprsm)*rwgt(iprsm+128)
               enddo
               ka_mn2o(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mn2o(jt,iprsm)*rwgt(iprsm+128)
            enddo
            kb_mn2o(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+128)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(9)
            sumk = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+128)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(9)
            sumf = 0.
            do ipr = 1, ngn(ngs(8)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(9)
         sumf = 0.
         do ipr = 1, ngn(ngs(8)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      end subroutine cmbgb9

!***************************************************************************
      subroutine cmbgb10
!***************************************************************************
!
!     band 10:  1390-1480 cm-1 (low key - h2o; high key - h2o)
!
! old band 10:  1390-1480 cm-1 (low - h2o; high - h2o)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng10
      use rrlw_kg10, only: fracrefao, fracrefbo, kao, kbo, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(10)
               sumk = 0.
               do ipr = 1, ngn(ngs(9)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+144)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(10)
               sumk = 0.
               do ipr = 1, ngn(ngs(9)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+144)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(10)
            sumk = 0.
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(10)
            sumk = 0.
            do ipr = 1, ngn(ngs(9)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+144)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(10)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(9)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb10

!***************************************************************************
      subroutine cmbgb11
!***************************************************************************
!
!     band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!
! old band 11:  1480-1800 cm-1 (low - h2o; low minor - o2)
!                              (high key - h2o; high minor - o2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng11
      use rrlw_kg11, only: fracrefao, fracrefbo, kao, kao_mo2, &
                           kbo, kbo_mo2, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, ka_mo2, &
                           kb, kb_mo2, selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(11)
               sumk = 0.
               do ipr = 1, ngn(ngs(10)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+160)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo
      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(11)
               sumk = 0.
               do ipr = 1, ngn(ngs(10)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+160)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(11)
            sumk1 = 0.
            sumk2 = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk1 = sumk1 + kao_mo2(jt,iprsm)*rwgt(iprsm+160)
               sumk2 = sumk2 + kbo_mo2(jt,iprsm)*rwgt(iprsm+160)
            enddo
            ka_mo2(jt,igc) = sumk1
            kb_mo2(jt,igc) = sumk2
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(11)
            sumk = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+160)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(11)
            sumk = 0.
            do ipr = 1, ngn(ngs(10)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+160)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(11)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(10)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb11

!***************************************************************************
      subroutine cmbgb12
!***************************************************************************
!
!     band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!
! old band 12:  1800-2080 cm-1 (low - h2o,co2; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng12
      use rrlw_kg12, only: fracrefao, kao, selfrefo, forrefo, &
                           fracrefa, ka, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(12)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(11)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+176)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(12)
            sumk = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+176)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(12)
            sumk = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+176)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(12)
            sumf = 0.
            do ipr = 1, ngn(ngs(11)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb12

!***************************************************************************
      subroutine cmbgb13
!***************************************************************************
!
!     band 13:  2080-2250 cm-1 (low key - h2o,n2o; high minor - o3 minor)
!
! old band 13:  2080-2250 cm-1 (low - h2o,n2o; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng13
      use rrlw_kg13, only: fracrefao, fracrefbo, kao, kao_mco2, kao_mco, &
                           kbo_mo3, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, ka_mco2, ka_mco, &
                           kb_mo3, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumk1, sumk2, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(13)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(12)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+192)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(13)
              sumk1 = 0.
              sumk2 = 0.
               do ipr = 1, ngn(ngs(12)+igc)
                  iprsm = iprsm + 1
                  sumk1 = sumk1 + kao_mco2(jn,jt,iprsm)*rwgt(iprsm+192)
                  sumk2 = sumk2 + kao_mco(jn,jt,iprsm)*rwgt(iprsm+192)
               enddo
               ka_mco2(jn,jt,igc) = sumk1
               ka_mco(jn,jt,igc) = sumk2
            enddo
         enddo
      enddo

      do jt = 1,19
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + kbo_mo3(jt,iprsm)*rwgt(iprsm+192)
            enddo
            kb_mo3(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+192)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(13)
            sumk = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+192)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(13)
         sumf = 0.
         do ipr = 1, ngn(ngs(12)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(13)
            sumf = 0.
            do ipr = 1, ngn(ngs(12)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb13

!***************************************************************************
      subroutine cmbgb14
!***************************************************************************
!
!     band 14:  2250-2380 cm-1 (low - co2; high - co2)
!
! old band 14:  2250-2380 cm-1 (low - co2; high - co2)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng14
      use rrlw_kg14, only: fracrefao, fracrefbo, kao, kbo, &
                           selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, &
                           selfref, forref

! ------- Local -------
      integer :: jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf1, sumf2


      do jt = 1,5
         do jp = 1,13
            iprsm = 0
            do igc = 1,ngc(14)
               sumk = 0.
               do ipr = 1, ngn(ngs(13)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao(jt,jp,iprsm)*rwgt(iprsm+208)
               enddo
               ka(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(14)
               sumk = 0.
               do ipr = 1, ngn(ngs(13)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+208)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(14)
            sumk = 0.
            do ipr = 1, ngn(ngs(13)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+208)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(14)
            sumk = 0.
            do ipr = 1, ngn(ngs(13)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+208)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(14)
         sumf1= 0.
         sumf2= 0.
         do ipr = 1, ngn(ngs(13)+igc)
            iprsm = iprsm + 1
            sumf1= sumf1+ fracrefao(iprsm)
            sumf2= sumf2+ fracrefbo(iprsm)
         enddo
         fracrefa(igc) = sumf1
         fracrefb(igc) = sumf2
      enddo

      end subroutine cmbgb14

!***************************************************************************
      subroutine cmbgb15
!***************************************************************************
!
!     band 15:  2380-2600 cm-1 (low - n2o,co2; low minor - n2)
!                              (high - nothing)
!
! old band 15:  2380-2600 cm-1 (low - n2o,co2; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng15
      use rrlw_kg15, only: fracrefao, kao, kao_mn2, selfrefo, forrefo, &
                           fracrefa, ka, ka_mn2, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(15)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(14)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+224)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jn = 1,9
         do jt = 1,19
            iprsm = 0
            do igc = 1,ngc(15)
              sumk = 0.
               do ipr = 1, ngn(ngs(14)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kao_mn2(jn,jt,iprsm)*rwgt(iprsm+224)
               enddo
               ka_mn2(jn,jt,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(15)
            sumk = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+224)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(15)
            sumk = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+224)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(15)
            sumf = 0.
            do ipr = 1, ngn(ngs(14)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb15

!***************************************************************************
      subroutine cmbgb16
!***************************************************************************
!
!     band 16:  2600-3250 cm-1 (low key- h2o,ch4; high key - ch4)
!
! old band 16:  2600-3000 cm-1 (low - h2o,ch4; high - nothing)
!***************************************************************************

      use parrrtm, only : mg, nbndlw, ngptlw, ng16
      use rrlw_kg16, only: fracrefao, fracrefbo, kao, kbo, selfrefo, forrefo, &
                           fracrefa, fracrefb, ka, kb, selfref, forref

! ------- Local -------
      integer :: jn, jt, jp, igc, ipr, iprsm 
      real(kind=r8) :: sumk, sumf


      do jn = 1,9
         do jt = 1,5
            do jp = 1,13
               iprsm = 0
               do igc = 1,ngc(16)
                  sumk = 0.
                  do ipr = 1, ngn(ngs(15)+igc)
                     iprsm = iprsm + 1
                     sumk = sumk + kao(jn,jt,jp,iprsm)*rwgt(iprsm+240)
                  enddo
                  ka(jn,jt,jp,igc) = sumk
               enddo
            enddo
         enddo
      enddo

      do jt = 1,5
         do jp = 13,59
            iprsm = 0
            do igc = 1,ngc(16)
               sumk = 0.
               do ipr = 1, ngn(ngs(15)+igc)
                  iprsm = iprsm + 1
                  sumk = sumk + kbo(jt,jp,iprsm)*rwgt(iprsm+240)
               enddo
               kb(jt,jp,igc) = sumk
            enddo
         enddo
      enddo

      do jt = 1,10
         iprsm = 0
         do igc = 1,ngc(16)
            sumk = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumk = sumk + selfrefo(jt,iprsm)*rwgt(iprsm+240)
            enddo
            selfref(jt,igc) = sumk
         enddo
      enddo

      do jt = 1,4
         iprsm = 0
         do igc = 1,ngc(16)
            sumk = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumk = sumk + forrefo(jt,iprsm)*rwgt(iprsm+240)
            enddo
            forref(jt,igc) = sumk
         enddo
      enddo

      iprsm = 0
      do igc = 1,ngc(16)
         sumf = 0.
         do ipr = 1, ngn(ngs(15)+igc)
            iprsm = iprsm + 1
            sumf = sumf + fracrefbo(iprsm)
         enddo
         fracrefb(igc) = sumf
      enddo

      do jp = 1,9
         iprsm = 0
         do igc = 1,ngc(16)
            sumf = 0.
            do ipr = 1, ngn(ngs(15)+igc)
               iprsm = iprsm + 1
               sumf = sumf + fracrefao(iprsm,jp)
            enddo
            fracrefa(igc,jp) = sumf
         enddo
      enddo

      end subroutine cmbgb16

!***************************************************************************
      subroutine lwcldpr
!***************************************************************************

! --------- Modules ----------

      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3, &
                          asyice2, asyice3, ssaice2, ssaice3, & !linhan
                          asyliq, ssaliq !linhan

      save

! ABSCLDn is the liquid water absorption coefficient (m2/g). 
! For INFLAG = 1.
      abscld1 = 0.0602410_r8
!  
! Everything below is for INFLAG = 2.

! ABSICEn(J,IB) are the parameters needed to compute the liquid water 
! absorption coefficient in spectral region IB for ICEFLAG=n.  The units
! of ABSICEn(1,IB) are m2/g and ABSICEn(2,IB) has units (microns (m2/g)).
! For ICEFLAG = 0.

      absice0(:)= (/0.005_r8,  1.0_r8/)

! For ICEFLAG = 1.
      absice1(1,:) = (/0.0036_r8, 0.0068_r8, 0.0003_r8, 0.0016_r8, 0.0020_r8/)
      absice1(2,:) = (/1.136_r8 , 0.600_r8 , 1.338_r8 , 1.166_r8 , 1.118_r8 /)

! For ICEFLAG = 2.  In each band, the absorption
! coefficients are listed for a range of effective radii from 5.0
! to 131.0 microns in increments of 3.0 microns.
! Spherical Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice2(:,1) = (/ &
! band 1
       7.798999e-02_r8,6.340479e-02_r8,5.417973e-02_r8,4.766245e-02_r8,4.272663e-02_r8, &
       3.880939e-02_r8,3.559544e-02_r8,3.289241e-02_r8,3.057511e-02_r8,2.855800e-02_r8, &
       2.678022e-02_r8,2.519712e-02_r8,2.377505e-02_r8,2.248806e-02_r8,2.131578e-02_r8, &
       2.024194e-02_r8,1.925337e-02_r8,1.833926e-02_r8,1.749067e-02_r8,1.670007e-02_r8, &
       1.596113e-02_r8,1.526845e-02_r8,1.461739e-02_r8,1.400394e-02_r8,1.342462e-02_r8, &
       1.287639e-02_r8,1.235656e-02_r8,1.186279e-02_r8,1.139297e-02_r8,1.094524e-02_r8, &
       1.051794e-02_r8,1.010956e-02_r8,9.718755e-03_r8,9.344316e-03_r8,8.985139e-03_r8, &
       8.640223e-03_r8,8.308656e-03_r8,7.989606e-03_r8,7.682312e-03_r8,7.386076e-03_r8, &
       7.100255e-03_r8,6.824258e-03_r8,6.557540e-03_r8/)
      absice2(:,2) = (/ &
! band 2
       2.784879e-02_r8,2.709863e-02_r8,2.619165e-02_r8,2.529230e-02_r8,2.443225e-02_r8, &
       2.361575e-02_r8,2.284021e-02_r8,2.210150e-02_r8,2.139548e-02_r8,2.071840e-02_r8, &
       2.006702e-02_r8,1.943856e-02_r8,1.883064e-02_r8,1.824120e-02_r8,1.766849e-02_r8, &
       1.711099e-02_r8,1.656737e-02_r8,1.603647e-02_r8,1.551727e-02_r8,1.500886e-02_r8, &
       1.451045e-02_r8,1.402132e-02_r8,1.354084e-02_r8,1.306842e-02_r8,1.260355e-02_r8, &
       1.214575e-02_r8,1.169460e-02_r8,1.124971e-02_r8,1.081072e-02_r8,1.037731e-02_r8, &
       9.949167e-03_r8,9.526021e-03_r8,9.107615e-03_r8,8.693714e-03_r8,8.284096e-03_r8, &
       7.878558e-03_r8,7.476910e-03_r8,7.078974e-03_r8,6.684586e-03_r8,6.293589e-03_r8, &
       5.905839e-03_r8,5.521200e-03_r8,5.139543e-03_r8/)
      absice2(:,3) = (/ &
! band 3
       1.065397e-01_r8,8.005726e-02_r8,6.546428e-02_r8,5.589131e-02_r8,4.898681e-02_r8, &
       4.369932e-02_r8,3.947901e-02_r8,3.600676e-02_r8,3.308299e-02_r8,3.057561e-02_r8, &
       2.839325e-02_r8,2.647040e-02_r8,2.475872e-02_r8,2.322164e-02_r8,2.183091e-02_r8, &
       2.056430e-02_r8,1.940407e-02_r8,1.833586e-02_r8,1.734787e-02_r8,1.643034e-02_r8, &
       1.557512e-02_r8,1.477530e-02_r8,1.402501e-02_r8,1.331924e-02_r8,1.265364e-02_r8, &
       1.202445e-02_r8,1.142838e-02_r8,1.086257e-02_r8,1.032445e-02_r8,9.811791e-03_r8, &
       9.322587e-03_r8,8.855053e-03_r8,8.407591e-03_r8,7.978763e-03_r8,7.567273e-03_r8, &
       7.171949e-03_r8,6.791728e-03_r8,6.425642e-03_r8,6.072809e-03_r8,5.732424e-03_r8, &
       5.403748e-03_r8,5.086103e-03_r8,4.778865e-03_r8/)
      absice2(:,4) = (/ &
! band 4
       1.804566e-01_r8,1.168987e-01_r8,8.680442e-02_r8,6.910060e-02_r8,5.738174e-02_r8, &
       4.902332e-02_r8,4.274585e-02_r8,3.784923e-02_r8,3.391734e-02_r8,3.068690e-02_r8, &
       2.798301e-02_r8,2.568480e-02_r8,2.370600e-02_r8,2.198337e-02_r8,2.046940e-02_r8, &
       1.912777e-02_r8,1.793016e-02_r8,1.685420e-02_r8,1.588193e-02_r8,1.499882e-02_r8, &
       1.419293e-02_r8,1.345440e-02_r8,1.277496e-02_r8,1.214769e-02_r8,1.156669e-02_r8, &
       1.102694e-02_r8,1.052412e-02_r8,1.005451e-02_r8,9.614854e-03_r8,9.202335e-03_r8, &
       8.814470e-03_r8,8.449077e-03_r8,8.104223e-03_r8,7.778195e-03_r8,7.469466e-03_r8, &
       7.176671e-03_r8,6.898588e-03_r8,6.634117e-03_r8,6.382264e-03_r8,6.142134e-03_r8, &
       5.912913e-03_r8,5.693862e-03_r8,5.484308e-03_r8/)
      absice2(:,5) = (/ &
! band 5
       2.131806e-01_r8,1.311372e-01_r8,9.407171e-02_r8,7.299442e-02_r8,5.941273e-02_r8, &
       4.994043e-02_r8,4.296242e-02_r8,3.761113e-02_r8,3.337910e-02_r8,2.994978e-02_r8, &
       2.711556e-02_r8,2.473461e-02_r8,2.270681e-02_r8,2.095943e-02_r8,1.943839e-02_r8, &
       1.810267e-02_r8,1.692057e-02_r8,1.586719e-02_r8,1.492275e-02_r8,1.407132e-02_r8, &
       1.329989e-02_r8,1.259780e-02_r8,1.195618e-02_r8,1.136761e-02_r8,1.082583e-02_r8, &
       1.032552e-02_r8,9.862158e-03_r8,9.431827e-03_r8,9.031157e-03_r8,8.657217e-03_r8, &
       8.307449e-03_r8,7.979609e-03_r8,7.671724e-03_r8,7.382048e-03_r8,7.109032e-03_r8, &
       6.851298e-03_r8,6.607615e-03_r8,6.376881e-03_r8,6.158105e-03_r8,5.950394e-03_r8, &
       5.752942e-03_r8,5.565019e-03_r8,5.385963e-03_r8/)
      absice2(:,6) = (/ &
! band 6
       1.546177e-01_r8,1.039251e-01_r8,7.910347e-02_r8,6.412429e-02_r8,5.399997e-02_r8, &
       4.664937e-02_r8,4.104237e-02_r8,3.660781e-02_r8,3.300218e-02_r8,3.000586e-02_r8, &
       2.747148e-02_r8,2.529633e-02_r8,2.340647e-02_r8,2.174723e-02_r8,2.027731e-02_r8, &
       1.896487e-02_r8,1.778492e-02_r8,1.671761e-02_r8,1.574692e-02_r8,1.485978e-02_r8, &
       1.404543e-02_r8,1.329489e-02_r8,1.260066e-02_r8,1.195636e-02_r8,1.135657e-02_r8, &
       1.079664e-02_r8,1.027257e-02_r8,9.780871e-03_r8,9.318505e-03_r8,8.882815e-03_r8, &
       8.471458e-03_r8,8.082364e-03_r8,7.713696e-03_r8,7.363817e-03_r8,7.031264e-03_r8, &
       6.714725e-03_r8,6.413021e-03_r8,6.125086e-03_r8,5.849958e-03_r8,5.586764e-03_r8, &
       5.334707e-03_r8,5.093066e-03_r8,4.861179e-03_r8/)
      absice2(:,7) = (/ &
! band 7
       7.583404e-02_r8,6.181558e-02_r8,5.312027e-02_r8,4.696039e-02_r8,4.225986e-02_r8, &
       3.849735e-02_r8,3.538340e-02_r8,3.274182e-02_r8,3.045798e-02_r8,2.845343e-02_r8, &
       2.667231e-02_r8,2.507353e-02_r8,2.362606e-02_r8,2.230595e-02_r8,2.109435e-02_r8, &
       1.997617e-02_r8,1.893916e-02_r8,1.797328e-02_r8,1.707016e-02_r8,1.622279e-02_r8, &
       1.542523e-02_r8,1.467241e-02_r8,1.395997e-02_r8,1.328414e-02_r8,1.264164e-02_r8, &
       1.202958e-02_r8,1.144544e-02_r8,1.088697e-02_r8,1.035218e-02_r8,9.839297e-03_r8, &
       9.346733e-03_r8,8.873057e-03_r8,8.416980e-03_r8,7.977335e-03_r8,7.553066e-03_r8, &
       7.143210e-03_r8,6.746888e-03_r8,6.363297e-03_r8,5.991700e-03_r8,5.631422e-03_r8, &
       5.281840e-03_r8,4.942378e-03_r8,4.612505e-03_r8/)
      absice2(:,8) = (/ &
! band 8
       9.022185e-02_r8,6.922700e-02_r8,5.710674e-02_r8,4.898377e-02_r8,4.305946e-02_r8, &
       3.849553e-02_r8,3.484183e-02_r8,3.183220e-02_r8,2.929794e-02_r8,2.712627e-02_r8, &
       2.523856e-02_r8,2.357810e-02_r8,2.210286e-02_r8,2.078089e-02_r8,1.958747e-02_r8, &
       1.850310e-02_r8,1.751218e-02_r8,1.660205e-02_r8,1.576232e-02_r8,1.498440e-02_r8, &
       1.426107e-02_r8,1.358624e-02_r8,1.295474e-02_r8,1.236212e-02_r8,1.180456e-02_r8, &
       1.127874e-02_r8,1.078175e-02_r8,1.031106e-02_r8,9.864433e-03_r8,9.439878e-03_r8, &
       9.035637e-03_r8,8.650140e-03_r8,8.281981e-03_r8,7.929895e-03_r8,7.592746e-03_r8, &
       7.269505e-03_r8,6.959238e-03_r8,6.661100e-03_r8,6.374317e-03_r8,6.098185e-03_r8, &
       5.832059e-03_r8,5.575347e-03_r8,5.327504e-03_r8/)
      absice2(:,9) = (/ &
! band 9
       1.294087e-01_r8,8.788217e-02_r8,6.728288e-02_r8,5.479720e-02_r8,4.635049e-02_r8, &
       4.022253e-02_r8,3.555576e-02_r8,3.187259e-02_r8,2.888498e-02_r8,2.640843e-02_r8, &
       2.431904e-02_r8,2.253038e-02_r8,2.098024e-02_r8,1.962267e-02_r8,1.842293e-02_r8, &
       1.735426e-02_r8,1.639571e-02_r8,1.553060e-02_r8,1.474552e-02_r8,1.402953e-02_r8, &
       1.337363e-02_r8,1.277033e-02_r8,1.221336e-02_r8,1.169741e-02_r8,1.121797e-02_r8, &
       1.077117e-02_r8,1.035369e-02_r8,9.962643e-03_r8,9.595509e-03_r8,9.250088e-03_r8, &
       8.924447e-03_r8,8.616876e-03_r8,8.325862e-03_r8,8.050057e-03_r8,7.788258e-03_r8, &
       7.539388e-03_r8,7.302478e-03_r8,7.076656e-03_r8,6.861134e-03_r8,6.655197e-03_r8, &
       6.458197e-03_r8,6.269543e-03_r8,6.088697e-03_r8/)
      absice2(:,10) = (/ &
! band 10
       1.593628e-01_r8,1.014552e-01_r8,7.458955e-02_r8,5.903571e-02_r8,4.887582e-02_r8, &
       4.171159e-02_r8,3.638480e-02_r8,3.226692e-02_r8,2.898717e-02_r8,2.631256e-02_r8, &
       2.408925e-02_r8,2.221156e-02_r8,2.060448e-02_r8,1.921325e-02_r8,1.799699e-02_r8, &
       1.692456e-02_r8,1.597177e-02_r8,1.511961e-02_r8,1.435289e-02_r8,1.365933e-02_r8, &
       1.302890e-02_r8,1.245334e-02_r8,1.192576e-02_r8,1.144037e-02_r8,1.099230e-02_r8, &
       1.057739e-02_r8,1.019208e-02_r8,9.833302e-03_r8,9.498395e-03_r8,9.185047e-03_r8, &
       8.891237e-03_r8,8.615185e-03_r8,8.355325e-03_r8,8.110267e-03_r8,7.878778e-03_r8, &
       7.659759e-03_r8,7.452224e-03_r8,7.255291e-03_r8,7.068166e-03_r8,6.890130e-03_r8, &
       6.720536e-03_r8,6.558794e-03_r8,6.404371e-03_r8/)
      absice2(:,11) = (/ &
! band 11
       1.656227e-01_r8,1.032129e-01_r8,7.487359e-02_r8,5.871431e-02_r8,4.828355e-02_r8, &
       4.099989e-02_r8,3.562924e-02_r8,3.150755e-02_r8,2.824593e-02_r8,2.560156e-02_r8, &
       2.341503e-02_r8,2.157740e-02_r8,2.001169e-02_r8,1.866199e-02_r8,1.748669e-02_r8, &
       1.645421e-02_r8,1.554015e-02_r8,1.472535e-02_r8,1.399457e-02_r8,1.333553e-02_r8, &
       1.273821e-02_r8,1.219440e-02_r8,1.169725e-02_r8,1.124104e-02_r8,1.082096e-02_r8, &
       1.043290e-02_r8,1.007336e-02_r8,9.739338e-03_r8,9.428223e-03_r8,9.137756e-03_r8, &
       8.865964e-03_r8,8.611115e-03_r8,8.371686e-03_r8,8.146330e-03_r8,7.933852e-03_r8, &
       7.733187e-03_r8,7.543386e-03_r8,7.363597e-03_r8,7.193056e-03_r8,7.031072e-03_r8, &
       6.877024e-03_r8,6.730348e-03_r8,6.590531e-03_r8/)
      absice2(:,12) = (/ &
! band 12
       9.194591e-02_r8,6.446867e-02_r8,4.962034e-02_r8,4.042061e-02_r8,3.418456e-02_r8, &
       2.968856e-02_r8,2.629900e-02_r8,2.365572e-02_r8,2.153915e-02_r8,1.980791e-02_r8, &
       1.836689e-02_r8,1.714979e-02_r8,1.610900e-02_r8,1.520946e-02_r8,1.442476e-02_r8, &
       1.373468e-02_r8,1.312345e-02_r8,1.257858e-02_r8,1.209010e-02_r8,1.164990e-02_r8, &
       1.125136e-02_r8,1.088901e-02_r8,1.055827e-02_r8,1.025531e-02_r8,9.976896e-03_r8, &
       9.720255e-03_r8,9.483022e-03_r8,9.263160e-03_r8,9.058902e-03_r8,8.868710e-03_r8, &
       8.691240e-03_r8,8.525312e-03_r8,8.369886e-03_r8,8.224042e-03_r8,8.086961e-03_r8, &
       7.957917e-03_r8,7.836258e-03_r8,7.721400e-03_r8,7.612821e-03_r8,7.510045e-03_r8, &
       7.412648e-03_r8,7.320242e-03_r8,7.232476e-03_r8/)
      absice2(:,13) = (/ &
! band 13
       1.437021e-01_r8,8.872535e-02_r8,6.392420e-02_r8,4.991833e-02_r8,4.096790e-02_r8, &
       3.477881e-02_r8,3.025782e-02_r8,2.681909e-02_r8,2.412102e-02_r8,2.195132e-02_r8, &
       2.017124e-02_r8,1.868641e-02_r8,1.743044e-02_r8,1.635529e-02_r8,1.542540e-02_r8, &
       1.461388e-02_r8,1.390003e-02_r8,1.326766e-02_r8,1.270395e-02_r8,1.219860e-02_r8, &
       1.174326e-02_r8,1.133107e-02_r8,1.095637e-02_r8,1.061442e-02_r8,1.030126e-02_r8, &
       1.001352e-02_r8,9.748340e-03_r8,9.503256e-03_r8,9.276155e-03_r8,9.065205e-03_r8, &
       8.868808e-03_r8,8.685571e-03_r8,8.514268e-03_r8,8.353820e-03_r8,8.203272e-03_r8, &
       8.061776e-03_r8,7.928578e-03_r8,7.803001e-03_r8,7.684443e-03_r8,7.572358e-03_r8, &
       7.466258e-03_r8,7.365701e-03_r8,7.270286e-03_r8/)
      absice2(:,14) = (/ &
! band 14
       1.288870e-01_r8,8.160295e-02_r8,5.964745e-02_r8,4.703790e-02_r8,3.888637e-02_r8, &
       3.320115e-02_r8,2.902017e-02_r8,2.582259e-02_r8,2.330224e-02_r8,2.126754e-02_r8, &
       1.959258e-02_r8,1.819130e-02_r8,1.700289e-02_r8,1.598320e-02_r8,1.509942e-02_r8, &
       1.432666e-02_r8,1.364572e-02_r8,1.304156e-02_r8,1.250220e-02_r8,1.201803e-02_r8, &
       1.158123e-02_r8,1.118537e-02_r8,1.082513e-02_r8,1.049605e-02_r8,1.019440e-02_r8, &
       9.916989e-03_r8,9.661116e-03_r8,9.424457e-03_r8,9.205005e-03_r8,9.001022e-03_r8, &
       8.810992e-03_r8,8.633588e-03_r8,8.467646e-03_r8,8.312137e-03_r8,8.166151e-03_r8, &
       8.028878e-03_r8,7.899597e-03_r8,7.777663e-03_r8,7.662498e-03_r8,7.553581e-03_r8, &
       7.450444e-03_r8,7.352662e-03_r8,7.259851e-03_r8/)
      absice2(:,15) = (/ &
! band 15
       8.254229e-02_r8,5.808787e-02_r8,4.492166e-02_r8,3.675028e-02_r8,3.119623e-02_r8, &
       2.718045e-02_r8,2.414450e-02_r8,2.177073e-02_r8,1.986526e-02_r8,1.830306e-02_r8, &
       1.699991e-02_r8,1.589698e-02_r8,1.495199e-02_r8,1.413374e-02_r8,1.341870e-02_r8, &
       1.278883e-02_r8,1.223002e-02_r8,1.173114e-02_r8,1.128322e-02_r8,1.087900e-02_r8, &
       1.051254e-02_r8,1.017890e-02_r8,9.873991e-03_r8,9.594347e-03_r8,9.337044e-03_r8, &
       9.099589e-03_r8,8.879842e-03_r8,8.675960e-03_r8,8.486341e-03_r8,8.309594e-03_r8, &
       8.144500e-03_r8,7.989986e-03_r8,7.845109e-03_r8,7.709031e-03_r8,7.581007e-03_r8, &
       7.460376e-03_r8,7.346544e-03_r8,7.238978e-03_r8,7.137201e-03_r8,7.040780e-03_r8, &
       6.949325e-03_r8,6.862483e-03_r8,6.779931e-03_r8/)
      absice2(:,16) = (/ &
! band 16
       1.382062e-01_r8,8.643227e-02_r8,6.282935e-02_r8,4.934783e-02_r8,4.063891e-02_r8, &
       3.455591e-02_r8,3.007059e-02_r8,2.662897e-02_r8,2.390631e-02_r8,2.169972e-02_r8, &
       1.987596e-02_r8,1.834393e-02_r8,1.703924e-02_r8,1.591513e-02_r8,1.493679e-02_r8, &
       1.407780e-02_r8,1.331775e-02_r8,1.264061e-02_r8,1.203364e-02_r8,1.148655e-02_r8, &
       1.099099e-02_r8,1.054006e-02_r8,1.012807e-02_r8,9.750215e-03_r8,9.402477e-03_r8, &
       9.081428e-03_r8,8.784143e-03_r8,8.508107e-03_r8,8.251146e-03_r8,8.011373e-03_r8, &
       7.787140e-03_r8,7.577002e-03_r8,7.379687e-03_r8,7.194071e-03_r8,7.019158e-03_r8, &
       6.854061e-03_r8,6.697986e-03_r8,6.550224e-03_r8,6.410138e-03_r8,6.277153e-03_r8, &
       6.150751e-03_r8,6.030462e-03_r8,5.915860e-03_r8/)

! ICEFLAG = 3; Fu parameterization. Particle size 5 - 140 micron in 
! increments of 3 microns.
! units = m2/g
! Hexagonal Ice Particle Parameterization
! absorption units (abs coef/iwc): [(m^-1)/(g m^-3)]
      absice3(:,1) = (/ &
! band 1
       3.110649e-03_r8,4.666352e-02_r8,6.606447e-02_r8,6.531678e-02_r8,6.012598e-02_r8, &
       5.437494e-02_r8,4.906411e-02_r8,4.441146e-02_r8,4.040585e-02_r8,3.697334e-02_r8, &
       3.403027e-02_r8,3.149979e-02_r8,2.931596e-02_r8,2.742365e-02_r8,2.577721e-02_r8, &
       2.433888e-02_r8,2.307732e-02_r8,2.196644e-02_r8,2.098437e-02_r8,2.011264e-02_r8, &
       1.933561e-02_r8,1.863992e-02_r8,1.801407e-02_r8,1.744812e-02_r8,1.693346e-02_r8, &
       1.646252e-02_r8,1.602866e-02_r8,1.562600e-02_r8,1.524933e-02_r8,1.489399e-02_r8, &
       1.455580e-02_r8,1.423098e-02_r8,1.391612e-02_r8,1.360812e-02_r8,1.330413e-02_r8, &
       1.300156e-02_r8,1.269801e-02_r8,1.239127e-02_r8,1.207928e-02_r8,1.176014e-02_r8, &
       1.143204e-02_r8,1.109334e-02_r8,1.074243e-02_r8,1.037786e-02_r8,9.998198e-03_r8, &
       9.602126e-03_r8/)
      absice3(:,2) = (/ &
! band 2
       3.984966e-04_r8,1.681097e-02_r8,2.627680e-02_r8,2.767465e-02_r8,2.700722e-02_r8, &
       2.579180e-02_r8,2.448677e-02_r8,2.323890e-02_r8,2.209096e-02_r8,2.104882e-02_r8, &
       2.010547e-02_r8,1.925003e-02_r8,1.847128e-02_r8,1.775883e-02_r8,1.710358e-02_r8, &
       1.649769e-02_r8,1.593449e-02_r8,1.540829e-02_r8,1.491429e-02_r8,1.444837e-02_r8, &
       1.400704e-02_r8,1.358729e-02_r8,1.318654e-02_r8,1.280258e-02_r8,1.243346e-02_r8, &
       1.207750e-02_r8,1.173325e-02_r8,1.139941e-02_r8,1.107487e-02_r8,1.075861e-02_r8, &
       1.044975e-02_r8,1.014753e-02_r8,9.851229e-03_r8,9.560240e-03_r8,9.274003e-03_r8, &
       8.992020e-03_r8,8.713845e-03_r8,8.439074e-03_r8,8.167346e-03_r8,7.898331e-03_r8, &
       7.631734e-03_r8,7.367286e-03_r8,7.104742e-03_r8,6.843882e-03_r8,6.584504e-03_r8, &
       6.326424e-03_r8/)
      absice3(:,3) = (/ &
! band 3
       6.933163e-02_r8,8.540475e-02_r8,7.701816e-02_r8,6.771158e-02_r8,5.986953e-02_r8, &
       5.348120e-02_r8,4.824962e-02_r8,4.390563e-02_r8,4.024411e-02_r8,3.711404e-02_r8, &
       3.440426e-02_r8,3.203200e-02_r8,2.993478e-02_r8,2.806474e-02_r8,2.638464e-02_r8, &
       2.486516e-02_r8,2.348288e-02_r8,2.221890e-02_r8,2.105780e-02_r8,1.998687e-02_r8, &
       1.899552e-02_r8,1.807490e-02_r8,1.721750e-02_r8,1.641693e-02_r8,1.566773e-02_r8, &
       1.496515e-02_r8,1.430509e-02_r8,1.368398e-02_r8,1.309865e-02_r8,1.254634e-02_r8, &
       1.202456e-02_r8,1.153114e-02_r8,1.106409e-02_r8,1.062166e-02_r8,1.020224e-02_r8, &
       9.804381e-03_r8,9.426771e-03_r8,9.068205e-03_r8,8.727578e-03_r8,8.403876e-03_r8, &
       8.096160e-03_r8,7.803564e-03_r8,7.525281e-03_r8,7.260560e-03_r8,7.008697e-03_r8, &
       6.769036e-03_r8/)
      absice3(:,4) = (/ &
! band 4
       1.765735e-01_r8,1.382700e-01_r8,1.095129e-01_r8,8.987475e-02_r8,7.591185e-02_r8, &
       6.554169e-02_r8,5.755500e-02_r8,5.122083e-02_r8,4.607610e-02_r8,4.181475e-02_r8, &
       3.822697e-02_r8,3.516432e-02_r8,3.251897e-02_r8,3.021073e-02_r8,2.817876e-02_r8, &
       2.637607e-02_r8,2.476582e-02_r8,2.331871e-02_r8,2.201113e-02_r8,2.082388e-02_r8, &
       1.974115e-02_r8,1.874983e-02_r8,1.783894e-02_r8,1.699922e-02_r8,1.622280e-02_r8, &
       1.550296e-02_r8,1.483390e-02_r8,1.421064e-02_r8,1.362880e-02_r8,1.308460e-02_r8, &
       1.257468e-02_r8,1.209611e-02_r8,1.164628e-02_r8,1.122287e-02_r8,1.082381e-02_r8, &
       1.044725e-02_r8,1.009154e-02_r8,9.755166e-03_r8,9.436783e-03_r8,9.135163e-03_r8, &
       8.849193e-03_r8,8.577856e-03_r8,8.320225e-03_r8,8.075451e-03_r8,7.842755e-03_r8, &
       7.621418e-03_r8/)
      absice3(:,5) = (/ &
! band 5
       2.339673e-01_r8,1.692124e-01_r8,1.291656e-01_r8,1.033837e-01_r8,8.562949e-02_r8, &
       7.273526e-02_r8,6.298262e-02_r8,5.537015e-02_r8,4.927787e-02_r8,4.430246e-02_r8, &
       4.017061e-02_r8,3.669072e-02_r8,3.372455e-02_r8,3.116995e-02_r8,2.894977e-02_r8, &
       2.700471e-02_r8,2.528842e-02_r8,2.376420e-02_r8,2.240256e-02_r8,2.117959e-02_r8, &
       2.007567e-02_r8,1.907456e-02_r8,1.816271e-02_r8,1.732874e-02_r8,1.656300e-02_r8, &
       1.585725e-02_r8,1.520445e-02_r8,1.459852e-02_r8,1.403419e-02_r8,1.350689e-02_r8, &
       1.301260e-02_r8,1.254781e-02_r8,1.210941e-02_r8,1.169468e-02_r8,1.130118e-02_r8, &
       1.092675e-02_r8,1.056945e-02_r8,1.022757e-02_r8,9.899560e-03_r8,9.584021e-03_r8, &
       9.279705e-03_r8,8.985479e-03_r8,8.700322e-03_r8,8.423306e-03_r8,8.153590e-03_r8, &
       7.890412e-03_r8/)
      absice3(:,6) = (/ &
! band 6
       1.145369e-01_r8,1.174566e-01_r8,9.917866e-02_r8,8.332990e-02_r8,7.104263e-02_r8, &
       6.153370e-02_r8,5.405472e-02_r8,4.806281e-02_r8,4.317918e-02_r8,3.913795e-02_r8, &
       3.574916e-02_r8,3.287437e-02_r8,3.041067e-02_r8,2.828017e-02_r8,2.642292e-02_r8, &
       2.479206e-02_r8,2.335051e-02_r8,2.206851e-02_r8,2.092195e-02_r8,1.989108e-02_r8, &
       1.895958e-02_r8,1.811385e-02_r8,1.734245e-02_r8,1.663573e-02_r8,1.598545e-02_r8, &
       1.538456e-02_r8,1.482700e-02_r8,1.430750e-02_r8,1.382150e-02_r8,1.336499e-02_r8, &
       1.293447e-02_r8,1.252685e-02_r8,1.213939e-02_r8,1.176968e-02_r8,1.141555e-02_r8, &
       1.107508e-02_r8,1.074655e-02_r8,1.042839e-02_r8,1.011923e-02_r8,9.817799e-03_r8, &
       9.522962e-03_r8,9.233688e-03_r8,8.949041e-03_r8,8.668171e-03_r8,8.390301e-03_r8, &
       8.114723e-03_r8/)
      absice3(:,7) = (/ &
! band 7
       1.222345e-02_r8,5.344230e-02_r8,5.523465e-02_r8,5.128759e-02_r8,4.676925e-02_r8, &
       4.266150e-02_r8,3.910561e-02_r8,3.605479e-02_r8,3.342843e-02_r8,3.115052e-02_r8, &
       2.915776e-02_r8,2.739935e-02_r8,2.583499e-02_r8,2.443266e-02_r8,2.316681e-02_r8, &
       2.201687e-02_r8,2.096619e-02_r8,2.000112e-02_r8,1.911044e-02_r8,1.828481e-02_r8, &
       1.751641e-02_r8,1.679866e-02_r8,1.612598e-02_r8,1.549360e-02_r8,1.489742e-02_r8, &
       1.433392e-02_r8,1.380002e-02_r8,1.329305e-02_r8,1.281068e-02_r8,1.235084e-02_r8, &
       1.191172e-02_r8,1.149171e-02_r8,1.108936e-02_r8,1.070341e-02_r8,1.033271e-02_r8, &
       9.976220e-03_r8,9.633021e-03_r8,9.302273e-03_r8,8.983216e-03_r8,8.675161e-03_r8, &
       8.377478e-03_r8,8.089595e-03_r8,7.810986e-03_r8,7.541170e-03_r8,7.279706e-03_r8, &
       7.026186e-03_r8/)
      absice3(:,8) = (/ &
! band 8
       6.711058e-02_r8,6.918198e-02_r8,6.127484e-02_r8,5.411944e-02_r8,4.836902e-02_r8, &
       4.375293e-02_r8,3.998077e-02_r8,3.683587e-02_r8,3.416508e-02_r8,3.186003e-02_r8, &
       2.984290e-02_r8,2.805671e-02_r8,2.645895e-02_r8,2.501733e-02_r8,2.370689e-02_r8, &
       2.250808e-02_r8,2.140532e-02_r8,2.038609e-02_r8,1.944018e-02_r8,1.855918e-02_r8, &
       1.773609e-02_r8,1.696504e-02_r8,1.624106e-02_r8,1.555990e-02_r8,1.491793e-02_r8, &
       1.431197e-02_r8,1.373928e-02_r8,1.319743e-02_r8,1.268430e-02_r8,1.219799e-02_r8, &
       1.173682e-02_r8,1.129925e-02_r8,1.088393e-02_r8,1.048961e-02_r8,1.011516e-02_r8, &
       9.759543e-03_r8,9.421813e-03_r8,9.101089e-03_r8,8.796559e-03_r8,8.507464e-03_r8, &
       8.233098e-03_r8,7.972798e-03_r8,7.725942e-03_r8,7.491940e-03_r8,7.270238e-03_r8, &
       7.060305e-03_r8/)
      absice3(:,9) = (/ &
! band 9
       1.236780e-01_r8,9.222386e-02_r8,7.383997e-02_r8,6.204072e-02_r8,5.381029e-02_r8, &
       4.770678e-02_r8,4.296928e-02_r8,3.916131e-02_r8,3.601540e-02_r8,3.335878e-02_r8, &
       3.107493e-02_r8,2.908247e-02_r8,2.732282e-02_r8,2.575276e-02_r8,2.433968e-02_r8, &
       2.305852e-02_r8,2.188966e-02_r8,2.081757e-02_r8,1.982974e-02_r8,1.891599e-02_r8, &
       1.806794e-02_r8,1.727865e-02_r8,1.654227e-02_r8,1.585387e-02_r8,1.520924e-02_r8, &
       1.460476e-02_r8,1.403730e-02_r8,1.350416e-02_r8,1.300293e-02_r8,1.253153e-02_r8, &
       1.208808e-02_r8,1.167094e-02_r8,1.127862e-02_r8,1.090979e-02_r8,1.056323e-02_r8, &
       1.023786e-02_r8,9.932665e-03_r8,9.646744e-03_r8,9.379250e-03_r8,9.129409e-03_r8, &
       8.896500e-03_r8,8.679856e-03_r8,8.478852e-03_r8,8.292904e-03_r8,8.121463e-03_r8, &
       7.964013e-03_r8/)
      absice3(:,10) = (/ &
! band 10
       1.655966e-01_r8,1.134205e-01_r8,8.714344e-02_r8,7.129241e-02_r8,6.063739e-02_r8, &
       5.294203e-02_r8,4.709309e-02_r8,4.247476e-02_r8,3.871892e-02_r8,3.559206e-02_r8, &
       3.293893e-02_r8,3.065226e-02_r8,2.865558e-02_r8,2.689288e-02_r8,2.532221e-02_r8, &
       2.391150e-02_r8,2.263582e-02_r8,2.147549e-02_r8,2.041476e-02_r8,1.944089e-02_r8, &
       1.854342e-02_r8,1.771371e-02_r8,1.694456e-02_r8,1.622989e-02_r8,1.556456e-02_r8, &
       1.494415e-02_r8,1.436491e-02_r8,1.382354e-02_r8,1.331719e-02_r8,1.284339e-02_r8, &
       1.239992e-02_r8,1.198486e-02_r8,1.159647e-02_r8,1.123323e-02_r8,1.089375e-02_r8, &
       1.057679e-02_r8,1.028124e-02_r8,1.000607e-02_r8,9.750376e-03_r8,9.513303e-03_r8, &
       9.294082e-03_r8,9.092003e-03_r8,8.906412e-03_r8,8.736702e-03_r8,8.582314e-03_r8, &
       8.442725e-03_r8/)
      absice3(:,11) = (/ &
! band 11
       1.775615e-01_r8,1.180046e-01_r8,8.929607e-02_r8,7.233500e-02_r8,6.108333e-02_r8, &
       5.303642e-02_r8,4.696927e-02_r8,4.221206e-02_r8,3.836768e-02_r8,3.518576e-02_r8, &
       3.250063e-02_r8,3.019825e-02_r8,2.819758e-02_r8,2.643943e-02_r8,2.487953e-02_r8, &
       2.348414e-02_r8,2.222705e-02_r8,2.108762e-02_r8,2.004936e-02_r8,1.909892e-02_r8, &
       1.822539e-02_r8,1.741975e-02_r8,1.667449e-02_r8,1.598330e-02_r8,1.534084e-02_r8, &
       1.474253e-02_r8,1.418446e-02_r8,1.366325e-02_r8,1.317597e-02_r8,1.272004e-02_r8, &
       1.229321e-02_r8,1.189350e-02_r8,1.151915e-02_r8,1.116859e-02_r8,1.084042e-02_r8, &
       1.053338e-02_r8,1.024636e-02_r8,9.978326e-03_r8,9.728357e-03_r8,9.495613e-03_r8, &
       9.279327e-03_r8,9.078798e-03_r8,8.893383e-03_r8,8.722488e-03_r8,8.565568e-03_r8, &
       8.422115e-03_r8/)
      absice3(:,12) = (/ &
! band 12
       9.465447e-02_r8,6.432047e-02_r8,5.060973e-02_r8,4.267283e-02_r8,3.741843e-02_r8, &
       3.363096e-02_r8,3.073531e-02_r8,2.842405e-02_r8,2.651789e-02_r8,2.490518e-02_r8, &
       2.351273e-02_r8,2.229056e-02_r8,2.120335e-02_r8,2.022541e-02_r8,1.933763e-02_r8, &
       1.852546e-02_r8,1.777763e-02_r8,1.708528e-02_r8,1.644134e-02_r8,1.584009e-02_r8, &
       1.527684e-02_r8,1.474774e-02_r8,1.424955e-02_r8,1.377957e-02_r8,1.333549e-02_r8, &
       1.291534e-02_r8,1.251743e-02_r8,1.214029e-02_r8,1.178265e-02_r8,1.144337e-02_r8, &
       1.112148e-02_r8,1.081609e-02_r8,1.052642e-02_r8,1.025178e-02_r8,9.991540e-03_r8, &
       9.745130e-03_r8,9.512038e-03_r8,9.291797e-03_r8,9.083980e-03_r8,8.888195e-03_r8, &
       8.704081e-03_r8,8.531306e-03_r8,8.369560e-03_r8,8.218558e-03_r8,8.078032e-03_r8, &
       7.947730e-03_r8/)
      absice3(:,13) = (/ &
! band 13
       1.560311e-01_r8,9.961097e-02_r8,7.502949e-02_r8,6.115022e-02_r8,5.214952e-02_r8, &
       4.578149e-02_r8,4.099731e-02_r8,3.724174e-02_r8,3.419343e-02_r8,3.165356e-02_r8, &
       2.949251e-02_r8,2.762222e-02_r8,2.598073e-02_r8,2.452322e-02_r8,2.321642e-02_r8, &
       2.203516e-02_r8,2.096002e-02_r8,1.997579e-02_r8,1.907036e-02_r8,1.823401e-02_r8, &
       1.745879e-02_r8,1.673819e-02_r8,1.606678e-02_r8,1.544003e-02_r8,1.485411e-02_r8, &
       1.430574e-02_r8,1.379215e-02_r8,1.331092e-02_r8,1.285996e-02_r8,1.243746e-02_r8, &
       1.204183e-02_r8,1.167164e-02_r8,1.132567e-02_r8,1.100281e-02_r8,1.070207e-02_r8, &
       1.042258e-02_r8,1.016352e-02_r8,9.924197e-03_r8,9.703953e-03_r8,9.502199e-03_r8, &
       9.318400e-03_r8,9.152066e-03_r8,9.002749e-03_r8,8.870038e-03_r8,8.753555e-03_r8, &
       8.652951e-03_r8/)
      absice3(:,14) = (/ &
! band 14
       1.559547e-01_r8,9.896700e-02_r8,7.441231e-02_r8,6.061469e-02_r8,5.168730e-02_r8, &
       4.537821e-02_r8,4.064106e-02_r8,3.692367e-02_r8,3.390714e-02_r8,3.139438e-02_r8, &
       2.925702e-02_r8,2.740783e-02_r8,2.578547e-02_r8,2.434552e-02_r8,2.305506e-02_r8, &
       2.188910e-02_r8,2.082842e-02_r8,1.985789e-02_r8,1.896553e-02_r8,1.814165e-02_r8, &
       1.737839e-02_r8,1.666927e-02_r8,1.600891e-02_r8,1.539279e-02_r8,1.481712e-02_r8, &
       1.427865e-02_r8,1.377463e-02_r8,1.330266e-02_r8,1.286068e-02_r8,1.244689e-02_r8, &
       1.205973e-02_r8,1.169780e-02_r8,1.135989e-02_r8,1.104492e-02_r8,1.075192e-02_r8, &
       1.048004e-02_r8,1.022850e-02_r8,9.996611e-03_r8,9.783753e-03_r8,9.589361e-03_r8, &
       9.412924e-03_r8,9.253977e-03_r8,9.112098e-03_r8,8.986903e-03_r8,8.878039e-03_r8, &
       8.785184e-03_r8/)
      absice3(:,15) = (/ &
! band 15
       1.102926e-01_r8,7.176622e-02_r8,5.530316e-02_r8,4.606056e-02_r8,4.006116e-02_r8, &
       3.579628e-02_r8,3.256909e-02_r8,3.001360e-02_r8,2.791920e-02_r8,2.615617e-02_r8, &
       2.464023e-02_r8,2.331426e-02_r8,2.213817e-02_r8,2.108301e-02_r8,2.012733e-02_r8, &
       1.925493e-02_r8,1.845331e-02_r8,1.771269e-02_r8,1.702531e-02_r8,1.638493e-02_r8, &
       1.578648e-02_r8,1.522579e-02_r8,1.469940e-02_r8,1.420442e-02_r8,1.373841e-02_r8, &
       1.329931e-02_r8,1.288535e-02_r8,1.249502e-02_r8,1.212700e-02_r8,1.178015e-02_r8, &
       1.145348e-02_r8,1.114612e-02_r8,1.085730e-02_r8,1.058633e-02_r8,1.033263e-02_r8, &
       1.009564e-02_r8,9.874895e-03_r8,9.669960e-03_r8,9.480449e-03_r8,9.306014e-03_r8, &
       9.146339e-03_r8,9.001138e-03_r8,8.870154e-03_r8,8.753148e-03_r8,8.649907e-03_r8, &
       8.560232e-03_r8/)
      absice3(:,16) = (/ &
! band 16
       1.688344e-01_r8,1.077072e-01_r8,7.994467e-02_r8,6.403862e-02_r8,5.369850e-02_r8, &
       4.641582e-02_r8,4.099331e-02_r8,3.678724e-02_r8,3.342069e-02_r8,3.065831e-02_r8, &
       2.834557e-02_r8,2.637680e-02_r8,2.467733e-02_r8,2.319286e-02_r8,2.188299e-02_r8, &
       2.071701e-02_r8,1.967121e-02_r8,1.872692e-02_r8,1.786931e-02_r8,1.708641e-02_r8, &
       1.636846e-02_r8,1.570743e-02_r8,1.509665e-02_r8,1.453052e-02_r8,1.400433e-02_r8, &
       1.351407e-02_r8,1.305631e-02_r8,1.262810e-02_r8,1.222688e-02_r8,1.185044e-02_r8, &
       1.149683e-02_r8,1.116436e-02_r8,1.085153e-02_r8,1.055701e-02_r8,1.027961e-02_r8, &
       1.001831e-02_r8,9.772141e-03_r8,9.540280e-03_r8,9.321966e-03_r8,9.116517e-03_r8, &
       8.923315e-03_r8,8.741803e-03_r8,8.571472e-03_r8,8.411860e-03_r8,8.262543e-03_r8, &
       8.123136e-03_r8/)

! For LIQFLAG = 0.
      absliq0 = 0.0903614_r8

! For LIQFLAG = 1.  In each band, the absorption
! coefficients are listed for a range of effective radii from 2.5
! to 59.5 microns in increments of 1.0 micron.
      absliq1(:, 1) = (/ &
! band  1
       1.64047e-03_r8, 6.90533e-02_r8, 7.72017e-02_r8, 7.78054e-02_r8, 7.69523e-02_r8, &
       7.58058e-02_r8, 7.46400e-02_r8, 7.35123e-02_r8, 7.24162e-02_r8, 7.13225e-02_r8, &
       6.99145e-02_r8, 6.66409e-02_r8, 6.36582e-02_r8, 6.09425e-02_r8, 5.84593e-02_r8, &
       5.61743e-02_r8, 5.40571e-02_r8, 5.20812e-02_r8, 5.02245e-02_r8, 4.84680e-02_r8, &
       4.67959e-02_r8, 4.51944e-02_r8, 4.36516e-02_r8, 4.21570e-02_r8, 4.07015e-02_r8, &
       3.92766e-02_r8, 3.78747e-02_r8, 3.64886e-02_r8, 3.53632e-02_r8, 3.41992e-02_r8, &
       3.31016e-02_r8, 3.20643e-02_r8, 3.10817e-02_r8, 3.01490e-02_r8, 2.92620e-02_r8, &
       2.84171e-02_r8, 2.76108e-02_r8, 2.68404e-02_r8, 2.61031e-02_r8, 2.53966e-02_r8, &
       2.47189e-02_r8, 2.40678e-02_r8, 2.34418e-02_r8, 2.28392e-02_r8, 2.22586e-02_r8, &
       2.16986e-02_r8, 2.11580e-02_r8, 2.06356e-02_r8, 2.01305e-02_r8, 1.96417e-02_r8, &
       1.91682e-02_r8, 1.87094e-02_r8, 1.82643e-02_r8, 1.78324e-02_r8, 1.74129e-02_r8, &
       1.70052e-02_r8, 1.66088e-02_r8, 1.62231e-02_r8/)
      absliq1(:, 2) = (/ &
! band  2
       2.19486e-01_r8, 1.80687e-01_r8, 1.59150e-01_r8, 1.44731e-01_r8, 1.33703e-01_r8, &
       1.24355e-01_r8, 1.15756e-01_r8, 1.07318e-01_r8, 9.86119e-02_r8, 8.92739e-02_r8, &
       8.34911e-02_r8, 7.70773e-02_r8, 7.15240e-02_r8, 6.66615e-02_r8, 6.23641e-02_r8, &
       5.85359e-02_r8, 5.51020e-02_r8, 5.20032e-02_r8, 4.91916e-02_r8, 4.66283e-02_r8, &
       4.42813e-02_r8, 4.21236e-02_r8, 4.01330e-02_r8, 3.82905e-02_r8, 3.65797e-02_r8, &
       3.49869e-02_r8, 3.35002e-02_r8, 3.21090e-02_r8, 3.08957e-02_r8, 2.97601e-02_r8, &
       2.86966e-02_r8, 2.76984e-02_r8, 2.67599e-02_r8, 2.58758e-02_r8, 2.50416e-02_r8, &
       2.42532e-02_r8, 2.35070e-02_r8, 2.27997e-02_r8, 2.21284e-02_r8, 2.14904e-02_r8, &
       2.08834e-02_r8, 2.03051e-02_r8, 1.97536e-02_r8, 1.92271e-02_r8, 1.87239e-02_r8, &
       1.82425e-02_r8, 1.77816e-02_r8, 1.73399e-02_r8, 1.69162e-02_r8, 1.65094e-02_r8, &
       1.61187e-02_r8, 1.57430e-02_r8, 1.53815e-02_r8, 1.50334e-02_r8, 1.46981e-02_r8, &
       1.43748e-02_r8, 1.40628e-02_r8, 1.37617e-02_r8/)
      absliq1(:, 3) = (/ &
! band  3
       2.95174e-01_r8, 2.34765e-01_r8, 1.98038e-01_r8, 1.72114e-01_r8, 1.52083e-01_r8, &
       1.35654e-01_r8, 1.21613e-01_r8, 1.09252e-01_r8, 9.81263e-02_r8, 8.79448e-02_r8, &
       8.12566e-02_r8, 7.44563e-02_r8, 6.86374e-02_r8, 6.36042e-02_r8, 5.92094e-02_r8, &
       5.53402e-02_r8, 5.19087e-02_r8, 4.88455e-02_r8, 4.60951e-02_r8, 4.36124e-02_r8, &
       4.13607e-02_r8, 3.93096e-02_r8, 3.74338e-02_r8, 3.57119e-02_r8, 3.41261e-02_r8, &
       3.26610e-02_r8, 3.13036e-02_r8, 3.00425e-02_r8, 2.88497e-02_r8, 2.78077e-02_r8, &
       2.68317e-02_r8, 2.59158e-02_r8, 2.50545e-02_r8, 2.42430e-02_r8, 2.34772e-02_r8, &
       2.27533e-02_r8, 2.20679e-02_r8, 2.14181e-02_r8, 2.08011e-02_r8, 2.02145e-02_r8, &
       1.96561e-02_r8, 1.91239e-02_r8, 1.86161e-02_r8, 1.81311e-02_r8, 1.76673e-02_r8, &
       1.72234e-02_r8, 1.67981e-02_r8, 1.63903e-02_r8, 1.59989e-02_r8, 1.56230e-02_r8, &
       1.52615e-02_r8, 1.49138e-02_r8, 1.45791e-02_r8, 1.42565e-02_r8, 1.39455e-02_r8, &
       1.36455e-02_r8, 1.33559e-02_r8, 1.30761e-02_r8/)
      absliq1(:, 4) = (/ &
! band  4
       3.00925e-01_r8, 2.36949e-01_r8, 1.96947e-01_r8, 1.68692e-01_r8, 1.47190e-01_r8, &
       1.29986e-01_r8, 1.15719e-01_r8, 1.03568e-01_r8, 9.30028e-02_r8, 8.36658e-02_r8, &
       7.71075e-02_r8, 7.07002e-02_r8, 6.52284e-02_r8, 6.05024e-02_r8, 5.63801e-02_r8, &
       5.27534e-02_r8, 4.95384e-02_r8, 4.66690e-02_r8, 4.40925e-02_r8, 4.17664e-02_r8, &
       3.96559e-02_r8, 3.77326e-02_r8, 3.59727e-02_r8, 3.43561e-02_r8, 3.28662e-02_r8, &
       3.14885e-02_r8, 3.02110e-02_r8, 2.90231e-02_r8, 2.78948e-02_r8, 2.69109e-02_r8, &
       2.59884e-02_r8, 2.51217e-02_r8, 2.43058e-02_r8, 2.35364e-02_r8, 2.28096e-02_r8, &
       2.21218e-02_r8, 2.14700e-02_r8, 2.08515e-02_r8, 2.02636e-02_r8, 1.97041e-02_r8, &
       1.91711e-02_r8, 1.86625e-02_r8, 1.81769e-02_r8, 1.77126e-02_r8, 1.72683e-02_r8, &
       1.68426e-02_r8, 1.64344e-02_r8, 1.60427e-02_r8, 1.56664e-02_r8, 1.53046e-02_r8, &
       1.49565e-02_r8, 1.46214e-02_r8, 1.42985e-02_r8, 1.39871e-02_r8, 1.36866e-02_r8, &
       1.33965e-02_r8, 1.31162e-02_r8, 1.28453e-02_r8/)
      absliq1(:, 5) = (/ &
! band  5
       2.64691e-01_r8, 2.12018e-01_r8, 1.78009e-01_r8, 1.53539e-01_r8, 1.34721e-01_r8, &
       1.19580e-01_r8, 1.06996e-01_r8, 9.62772e-02_r8, 8.69710e-02_r8, 7.87670e-02_r8, &
       7.29272e-02_r8, 6.70920e-02_r8, 6.20977e-02_r8, 5.77732e-02_r8, 5.39910e-02_r8, &
       5.06538e-02_r8, 4.76866e-02_r8, 4.50301e-02_r8, 4.26374e-02_r8, 4.04704e-02_r8, &
       3.84981e-02_r8, 3.66948e-02_r8, 3.50394e-02_r8, 3.35141e-02_r8, 3.21038e-02_r8, &
       3.07957e-02_r8, 2.95788e-02_r8, 2.84438e-02_r8, 2.73790e-02_r8, 2.64390e-02_r8, &
       2.55565e-02_r8, 2.47263e-02_r8, 2.39437e-02_r8, 2.32047e-02_r8, 2.25056e-02_r8, &
       2.18433e-02_r8, 2.12149e-02_r8, 2.06177e-02_r8, 2.00495e-02_r8, 1.95081e-02_r8, &
       1.89917e-02_r8, 1.84984e-02_r8, 1.80269e-02_r8, 1.75755e-02_r8, 1.71431e-02_r8, &
       1.67283e-02_r8, 1.63303e-02_r8, 1.59478e-02_r8, 1.55801e-02_r8, 1.52262e-02_r8, &
       1.48853e-02_r8, 1.45568e-02_r8, 1.42400e-02_r8, 1.39342e-02_r8, 1.36388e-02_r8, &
       1.33533e-02_r8, 1.30773e-02_r8, 1.28102e-02_r8/)
      absliq1(:, 6) = (/ &
! band  6
       8.81182e-02_r8, 1.06745e-01_r8, 9.79753e-02_r8, 8.99625e-02_r8, 8.35200e-02_r8, &
       7.81899e-02_r8, 7.35939e-02_r8, 6.94696e-02_r8, 6.56266e-02_r8, 6.19148e-02_r8, &
       5.83355e-02_r8, 5.49306e-02_r8, 5.19642e-02_r8, 4.93325e-02_r8, 4.69659e-02_r8, &
       4.48148e-02_r8, 4.28431e-02_r8, 4.10231e-02_r8, 3.93332e-02_r8, 3.77563e-02_r8, &
       3.62785e-02_r8, 3.48882e-02_r8, 3.35758e-02_r8, 3.23333e-02_r8, 3.11536e-02_r8, &
       3.00310e-02_r8, 2.89601e-02_r8, 2.79365e-02_r8, 2.70502e-02_r8, 2.62618e-02_r8, &
       2.55025e-02_r8, 2.47728e-02_r8, 2.40726e-02_r8, 2.34013e-02_r8, 2.27583e-02_r8, &
       2.21422e-02_r8, 2.15522e-02_r8, 2.09869e-02_r8, 2.04453e-02_r8, 1.99260e-02_r8, &
       1.94280e-02_r8, 1.89501e-02_r8, 1.84913e-02_r8, 1.80506e-02_r8, 1.76270e-02_r8, &
       1.72196e-02_r8, 1.68276e-02_r8, 1.64500e-02_r8, 1.60863e-02_r8, 1.57357e-02_r8, &
       1.53975e-02_r8, 1.50710e-02_r8, 1.47558e-02_r8, 1.44511e-02_r8, 1.41566e-02_r8, &
       1.38717e-02_r8, 1.35960e-02_r8, 1.33290e-02_r8/)
      absliq1(:, 7) = (/ &
! band  7
       4.32174e-02_r8, 7.36078e-02_r8, 6.98340e-02_r8, 6.65231e-02_r8, 6.41948e-02_r8, &
       6.23551e-02_r8, 6.06638e-02_r8, 5.88680e-02_r8, 5.67124e-02_r8, 5.38629e-02_r8, &
       4.99579e-02_r8, 4.86289e-02_r8, 4.70120e-02_r8, 4.52854e-02_r8, 4.35466e-02_r8, &
       4.18480e-02_r8, 4.02169e-02_r8, 3.86658e-02_r8, 3.71992e-02_r8, 3.58168e-02_r8, &
       3.45155e-02_r8, 3.32912e-02_r8, 3.21390e-02_r8, 3.10538e-02_r8, 3.00307e-02_r8, &
       2.90651e-02_r8, 2.81524e-02_r8, 2.72885e-02_r8, 2.62821e-02_r8, 2.55744e-02_r8, &
       2.48799e-02_r8, 2.42029e-02_r8, 2.35460e-02_r8, 2.29108e-02_r8, 2.22981e-02_r8, &
       2.17079e-02_r8, 2.11402e-02_r8, 2.05945e-02_r8, 2.00701e-02_r8, 1.95663e-02_r8, &
       1.90824e-02_r8, 1.86174e-02_r8, 1.81706e-02_r8, 1.77411e-02_r8, 1.73281e-02_r8, &
       1.69307e-02_r8, 1.65483e-02_r8, 1.61801e-02_r8, 1.58254e-02_r8, 1.54835e-02_r8, &
       1.51538e-02_r8, 1.48358e-02_r8, 1.45288e-02_r8, 1.42322e-02_r8, 1.39457e-02_r8, &
       1.36687e-02_r8, 1.34008e-02_r8, 1.31416e-02_r8/)
      absliq1(:, 8) = (/ &
! band  8
       1.41881e-01_r8, 7.15419e-02_r8, 6.30335e-02_r8, 6.11132e-02_r8, 6.01931e-02_r8, &
       5.92420e-02_r8, 5.78968e-02_r8, 5.58876e-02_r8, 5.28923e-02_r8, 4.84462e-02_r8, &
       4.60839e-02_r8, 4.56013e-02_r8, 4.45410e-02_r8, 4.31866e-02_r8, 4.17026e-02_r8, &
       4.01850e-02_r8, 3.86892e-02_r8, 3.72461e-02_r8, 3.58722e-02_r8, 3.45749e-02_r8, &
       3.33564e-02_r8, 3.22155e-02_r8, 3.11494e-02_r8, 3.01541e-02_r8, 2.92253e-02_r8, &
       2.83584e-02_r8, 2.75488e-02_r8, 2.67925e-02_r8, 2.57692e-02_r8, 2.50704e-02_r8, &
       2.43918e-02_r8, 2.37350e-02_r8, 2.31005e-02_r8, 2.24888e-02_r8, 2.18996e-02_r8, &
       2.13325e-02_r8, 2.07870e-02_r8, 2.02623e-02_r8, 1.97577e-02_r8, 1.92724e-02_r8, &
       1.88056e-02_r8, 1.83564e-02_r8, 1.79241e-02_r8, 1.75079e-02_r8, 1.71070e-02_r8, &
       1.67207e-02_r8, 1.63482e-02_r8, 1.59890e-02_r8, 1.56424e-02_r8, 1.53077e-02_r8, &
       1.49845e-02_r8, 1.46722e-02_r8, 1.43702e-02_r8, 1.40782e-02_r8, 1.37955e-02_r8, &
       1.35219e-02_r8, 1.32569e-02_r8, 1.30000e-02_r8/)
      absliq1(:, 9) = (/ &
! band  9
       6.72726e-02_r8, 6.61013e-02_r8, 6.47866e-02_r8, 6.33780e-02_r8, 6.18985e-02_r8, &
       6.03335e-02_r8, 5.86136e-02_r8, 5.65876e-02_r8, 5.39839e-02_r8, 5.03536e-02_r8, &
       4.71608e-02_r8, 4.63630e-02_r8, 4.50313e-02_r8, 4.34526e-02_r8, 4.17876e-02_r8, &
       4.01261e-02_r8, 3.85171e-02_r8, 3.69860e-02_r8, 3.55442e-02_r8, 3.41954e-02_r8, &
       3.29384e-02_r8, 3.17693e-02_r8, 3.06832e-02_r8, 2.96745e-02_r8, 2.87374e-02_r8, &
       2.78662e-02_r8, 2.70557e-02_r8, 2.63008e-02_r8, 2.52450e-02_r8, 2.45424e-02_r8, &
       2.38656e-02_r8, 2.32144e-02_r8, 2.25885e-02_r8, 2.19873e-02_r8, 2.14099e-02_r8, &
       2.08554e-02_r8, 2.03230e-02_r8, 1.98116e-02_r8, 1.93203e-02_r8, 1.88482e-02_r8, &
       1.83944e-02_r8, 1.79578e-02_r8, 1.75378e-02_r8, 1.71335e-02_r8, 1.67440e-02_r8, &
       1.63687e-02_r8, 1.60069e-02_r8, 1.56579e-02_r8, 1.53210e-02_r8, 1.49958e-02_r8, &
       1.46815e-02_r8, 1.43778e-02_r8, 1.40841e-02_r8, 1.37999e-02_r8, 1.35249e-02_r8, &
       1.32585e-02_r8, 1.30004e-02_r8, 1.27502e-02_r8/)
      absliq1(:,10) = (/ &
! band 10
       7.97040e-02_r8, 7.63844e-02_r8, 7.36499e-02_r8, 7.13525e-02_r8, 6.93043e-02_r8, &
       6.72807e-02_r8, 6.50227e-02_r8, 6.22395e-02_r8, 5.86093e-02_r8, 5.37815e-02_r8, &
       5.14682e-02_r8, 4.97214e-02_r8, 4.77392e-02_r8, 4.56961e-02_r8, 4.36858e-02_r8, &
       4.17569e-02_r8, 3.99328e-02_r8, 3.82224e-02_r8, 3.66265e-02_r8, 3.51416e-02_r8, &
       3.37617e-02_r8, 3.24798e-02_r8, 3.12887e-02_r8, 3.01812e-02_r8, 2.91505e-02_r8, &
       2.81900e-02_r8, 2.72939e-02_r8, 2.64568e-02_r8, 2.54165e-02_r8, 2.46832e-02_r8, &
       2.39783e-02_r8, 2.33017e-02_r8, 2.26531e-02_r8, 2.20314e-02_r8, 2.14359e-02_r8, &
       2.08653e-02_r8, 2.03187e-02_r8, 1.97947e-02_r8, 1.92924e-02_r8, 1.88106e-02_r8, &
       1.83483e-02_r8, 1.79043e-02_r8, 1.74778e-02_r8, 1.70678e-02_r8, 1.66735e-02_r8, &
       1.62941e-02_r8, 1.59286e-02_r8, 1.55766e-02_r8, 1.52371e-02_r8, 1.49097e-02_r8, &
       1.45937e-02_r8, 1.42885e-02_r8, 1.39936e-02_r8, 1.37085e-02_r8, 1.34327e-02_r8, &
       1.31659e-02_r8, 1.29075e-02_r8, 1.26571e-02_r8/)
      absliq1(:,11) = (/ &
! band 11
       1.49438e-01_r8, 1.33535e-01_r8, 1.21542e-01_r8, 1.11743e-01_r8, 1.03263e-01_r8, &
       9.55774e-02_r8, 8.83382e-02_r8, 8.12943e-02_r8, 7.42533e-02_r8, 6.70609e-02_r8, &
       6.38761e-02_r8, 5.97788e-02_r8, 5.59841e-02_r8, 5.25318e-02_r8, 4.94132e-02_r8, &
       4.66014e-02_r8, 4.40644e-02_r8, 4.17706e-02_r8, 3.96910e-02_r8, 3.77998e-02_r8, &
       3.60742e-02_r8, 3.44947e-02_r8, 3.30442e-02_r8, 3.17079e-02_r8, 3.04730e-02_r8, &
       2.93283e-02_r8, 2.82642e-02_r8, 2.72720e-02_r8, 2.61789e-02_r8, 2.53277e-02_r8, &
       2.45237e-02_r8, 2.37635e-02_r8, 2.30438e-02_r8, 2.23615e-02_r8, 2.17140e-02_r8, &
       2.10987e-02_r8, 2.05133e-02_r8, 1.99557e-02_r8, 1.94241e-02_r8, 1.89166e-02_r8, &
       1.84317e-02_r8, 1.79679e-02_r8, 1.75238e-02_r8, 1.70983e-02_r8, 1.66901e-02_r8, &
       1.62983e-02_r8, 1.59219e-02_r8, 1.55599e-02_r8, 1.52115e-02_r8, 1.48761e-02_r8, &
       1.45528e-02_r8, 1.42411e-02_r8, 1.39402e-02_r8, 1.36497e-02_r8, 1.33690e-02_r8, &
       1.30976e-02_r8, 1.28351e-02_r8, 1.25810e-02_r8/)
      absliq1(:,12) = (/ &
! band 12
       3.71985e-02_r8, 3.88586e-02_r8, 3.99070e-02_r8, 4.04351e-02_r8, 4.04610e-02_r8, &
       3.99834e-02_r8, 3.89953e-02_r8, 3.74886e-02_r8, 3.54551e-02_r8, 3.28870e-02_r8, &
       3.32576e-02_r8, 3.22444e-02_r8, 3.12384e-02_r8, 3.02584e-02_r8, 2.93146e-02_r8, &
       2.84120e-02_r8, 2.75525e-02_r8, 2.67361e-02_r8, 2.59618e-02_r8, 2.52280e-02_r8, &
       2.45327e-02_r8, 2.38736e-02_r8, 2.32487e-02_r8, 2.26558e-02_r8, 2.20929e-02_r8, &
       2.15579e-02_r8, 2.10491e-02_r8, 2.05648e-02_r8, 1.99749e-02_r8, 1.95704e-02_r8, &
       1.91731e-02_r8, 1.87839e-02_r8, 1.84032e-02_r8, 1.80315e-02_r8, 1.76689e-02_r8, &
       1.73155e-02_r8, 1.69712e-02_r8, 1.66362e-02_r8, 1.63101e-02_r8, 1.59928e-02_r8, &
       1.56842e-02_r8, 1.53840e-02_r8, 1.50920e-02_r8, 1.48080e-02_r8, 1.45318e-02_r8, &
       1.42631e-02_r8, 1.40016e-02_r8, 1.37472e-02_r8, 1.34996e-02_r8, 1.32586e-02_r8, &
       1.30239e-02_r8, 1.27954e-02_r8, 1.25728e-02_r8, 1.23559e-02_r8, 1.21445e-02_r8, &
       1.19385e-02_r8, 1.17376e-02_r8, 1.15417e-02_r8/)
      absliq1(:,13) = (/ &
! band 13
       3.11868e-02_r8, 4.48357e-02_r8, 4.90224e-02_r8, 4.96406e-02_r8, 4.86806e-02_r8, &
       4.69610e-02_r8, 4.48630e-02_r8, 4.25795e-02_r8, 4.02138e-02_r8, 3.78236e-02_r8, &
       3.74266e-02_r8, 3.60384e-02_r8, 3.47074e-02_r8, 3.34434e-02_r8, 3.22499e-02_r8, &
       3.11264e-02_r8, 3.00704e-02_r8, 2.90784e-02_r8, 2.81463e-02_r8, 2.72702e-02_r8, &
       2.64460e-02_r8, 2.56698e-02_r8, 2.49381e-02_r8, 2.42475e-02_r8, 2.35948e-02_r8, &
       2.29774e-02_r8, 2.23925e-02_r8, 2.18379e-02_r8, 2.11793e-02_r8, 2.07076e-02_r8, &
       2.02470e-02_r8, 1.97981e-02_r8, 1.93613e-02_r8, 1.89367e-02_r8, 1.85243e-02_r8, &
       1.81240e-02_r8, 1.77356e-02_r8, 1.73588e-02_r8, 1.69935e-02_r8, 1.66392e-02_r8, &
       1.62956e-02_r8, 1.59624e-02_r8, 1.56393e-02_r8, 1.53259e-02_r8, 1.50219e-02_r8, &
       1.47268e-02_r8, 1.44404e-02_r8, 1.41624e-02_r8, 1.38925e-02_r8, 1.36302e-02_r8, &
       1.33755e-02_r8, 1.31278e-02_r8, 1.28871e-02_r8, 1.26530e-02_r8, 1.24253e-02_r8, &
       1.22038e-02_r8, 1.19881e-02_r8, 1.17782e-02_r8/)
      absliq1(:,14) = (/ &
! band 14
       1.58988e-02_r8, 3.50652e-02_r8, 4.00851e-02_r8, 4.07270e-02_r8, 3.98101e-02_r8, &
       3.83306e-02_r8, 3.66829e-02_r8, 3.50327e-02_r8, 3.34497e-02_r8, 3.19609e-02_r8, &
       3.13712e-02_r8, 3.03348e-02_r8, 2.93415e-02_r8, 2.83973e-02_r8, 2.75037e-02_r8, &
       2.66604e-02_r8, 2.58654e-02_r8, 2.51161e-02_r8, 2.44100e-02_r8, 2.37440e-02_r8, &
       2.31154e-02_r8, 2.25215e-02_r8, 2.19599e-02_r8, 2.14282e-02_r8, 2.09242e-02_r8, &
       2.04459e-02_r8, 1.99915e-02_r8, 1.95594e-02_r8, 1.90254e-02_r8, 1.86598e-02_r8, &
       1.82996e-02_r8, 1.79455e-02_r8, 1.75983e-02_r8, 1.72584e-02_r8, 1.69260e-02_r8, &
       1.66013e-02_r8, 1.62843e-02_r8, 1.59752e-02_r8, 1.56737e-02_r8, 1.53799e-02_r8, &
       1.50936e-02_r8, 1.48146e-02_r8, 1.45429e-02_r8, 1.42782e-02_r8, 1.40203e-02_r8, &
       1.37691e-02_r8, 1.35243e-02_r8, 1.32858e-02_r8, 1.30534e-02_r8, 1.28270e-02_r8, &
       1.26062e-02_r8, 1.23909e-02_r8, 1.21810e-02_r8, 1.19763e-02_r8, 1.17766e-02_r8, &
       1.15817e-02_r8, 1.13915e-02_r8, 1.12058e-02_r8/)
      absliq1(:,15) = (/ &
! band 15
       5.02079e-03_r8, 2.17615e-02_r8, 2.55449e-02_r8, 2.59484e-02_r8, 2.53650e-02_r8, &
       2.45281e-02_r8, 2.36843e-02_r8, 2.29159e-02_r8, 2.22451e-02_r8, 2.16716e-02_r8, &
       2.11451e-02_r8, 2.05817e-02_r8, 2.00454e-02_r8, 1.95372e-02_r8, 1.90567e-02_r8, &
       1.86028e-02_r8, 1.81742e-02_r8, 1.77693e-02_r8, 1.73866e-02_r8, 1.70244e-02_r8, &
       1.66815e-02_r8, 1.63563e-02_r8, 1.60477e-02_r8, 1.57544e-02_r8, 1.54755e-02_r8, &
       1.52097e-02_r8, 1.49564e-02_r8, 1.47146e-02_r8, 1.43684e-02_r8, 1.41728e-02_r8, &
       1.39762e-02_r8, 1.37797e-02_r8, 1.35838e-02_r8, 1.33891e-02_r8, 1.31961e-02_r8, &
       1.30051e-02_r8, 1.28164e-02_r8, 1.26302e-02_r8, 1.24466e-02_r8, 1.22659e-02_r8, &
       1.20881e-02_r8, 1.19131e-02_r8, 1.17412e-02_r8, 1.15723e-02_r8, 1.14063e-02_r8, &
       1.12434e-02_r8, 1.10834e-02_r8, 1.09264e-02_r8, 1.07722e-02_r8, 1.06210e-02_r8, &
       1.04725e-02_r8, 1.03269e-02_r8, 1.01839e-02_r8, 1.00436e-02_r8, 9.90593e-03_r8, &
       9.77080e-03_r8, 9.63818e-03_r8, 9.50800e-03_r8/)
      absliq1(:,16) = (/ &
! band 16
       5.64971e-02_r8, 9.04736e-02_r8, 8.11726e-02_r8, 7.05450e-02_r8, 6.20052e-02_r8, &
       5.54286e-02_r8, 5.03503e-02_r8, 4.63791e-02_r8, 4.32290e-02_r8, 4.06959e-02_r8, &
       3.74690e-02_r8, 3.52964e-02_r8, 3.33799e-02_r8, 3.16774e-02_r8, 3.01550e-02_r8, &
       2.87856e-02_r8, 2.75474e-02_r8, 2.64223e-02_r8, 2.53953e-02_r8, 2.44542e-02_r8, &
       2.35885e-02_r8, 2.27894e-02_r8, 2.20494e-02_r8, 2.13622e-02_r8, 2.07222e-02_r8, &
       2.01246e-02_r8, 1.95654e-02_r8, 1.90408e-02_r8, 1.84398e-02_r8, 1.80021e-02_r8, &
       1.75816e-02_r8, 1.71775e-02_r8, 1.67889e-02_r8, 1.64152e-02_r8, 1.60554e-02_r8, &
       1.57089e-02_r8, 1.53751e-02_r8, 1.50531e-02_r8, 1.47426e-02_r8, 1.44428e-02_r8, &
       1.41532e-02_r8, 1.38734e-02_r8, 1.36028e-02_r8, 1.33410e-02_r8, 1.30875e-02_r8, &
       1.28420e-02_r8, 1.26041e-02_r8, 1.23735e-02_r8, 1.21497e-02_r8, 1.19325e-02_r8, &
       1.17216e-02_r8, 1.15168e-02_r8, 1.13177e-02_r8, 1.11241e-02_r8, 1.09358e-02_r8, &
       1.07525e-02_r8, 1.05741e-02_r8, 1.04003e-02_r8/)

!linhan add
      !     SINGLE-SCATTERING ALBEDO: Unitless
!BAND 1
    SSALIQ(:,1) = (/ &
        4.180200e-02_r8,1.617930e-01_r8,2.609520e-01_r8,3.264170e-01_r8,3.691660e-01_r8,&
        3.976570e-01_r8,4.171290e-01_r8,4.307740e-01_r8,4.405960e-01_r8,4.478810e-01_r8,&
        4.534720e-01_r8,4.579280e-01_r8,4.616190e-01_r8,4.647910e-01_r8,4.676080e-01_r8,&
        4.701740e-01_r8,4.725590e-01_r8,4.748050e-01_r8,4.769410e-01_r8,4.789850e-01_r8,&
        4.809470e-01_r8,4.828360e-01_r8,4.846570e-01_r8,4.864160e-01_r8,4.881160e-01_r8,&
        4.897600e-01_r8,4.913530e-01_r8,4.928960e-01_r8,4.943920e-01_r8,4.958440e-01_r8/)
!BAND 2
    SSALIQ(:,2) = (/ &
        1.167210e-01_r8,2.846730e-01_r8,3.692340e-01_r8,4.123190e-01_r8,4.358070e-01_r8,&
        4.493970e-01_r8,4.578450e-01_r8,4.635970e-01_r8,4.679350e-01_r8,4.715250e-01_r8,&
        4.747120e-01_r8,4.776700e-01_r8,4.804830e-01_r8,4.831890e-01_r8,4.858040e-01_r8,&
        4.883300e-01_r8,4.907660e-01_r8,4.931070e-01_r8,4.953500e-01_r8,4.974930e-01_r8,&
        4.995350e-01_r8,5.014780e-01_r8,5.033250e-01_r8,5.050790e-01_r8,5.067460e-01_r8,&
        5.083310e-01_r8,5.098390e-01_r8,5.112740e-01_r8,5.126430e-01_r8,5.139490e-01_r8/)
!BAND 3
    SSALIQ(:,3) = (/ &
        1.432150e-01_r8,2.919690e-01_r8,3.646840e-01_r8,4.037180e-01_r8,4.268360e-01_r8,&
        4.418900e-01_r8,4.526070e-01_r8,4.608320e-01_r8,4.675200e-01_r8,4.731850e-01_r8,&
        4.781180e-01_r8,4.824950e-01_r8,4.864330e-01_r8,4.900110e-01_r8,4.932840e-01_r8,&
        4.962950e-01_r8,4.990740e-01_r8,5.016460e-01_r8,5.040300e-01_r8,5.062440e-01_r8,&
        5.083030e-01_r8,5.102190e-01_r8,5.120070e-01_r8,5.136780e-01_r8,5.152430e-01_r8,&
        5.167100e-01_r8,5.180910e-01_r8,5.193910e-01_r8,5.206190e-01_r8,5.217820e-01_r8/)
!BAND 4
    SSALIQ(:,4) = (/ &
        1.410170e-01_r8,2.781060e-01_r8,3.495930e-01_r8,3.908550e-01_r8,4.167900e-01_r8,&
        4.344190e-01_r8,4.472450e-01_r8,4.571030e-01_r8,4.650100e-01_r8,4.715580e-01_r8,&
        4.771140e-01_r8,4.819190e-01_r8,4.861380e-01_r8,4.898880e-01_r8,4.932530e-01_r8,&
        4.962980e-01_r8,4.990690e-01_r8,5.016020e-01_r8,5.039260e-01_r8,5.060650e-01_r8,&
        5.080380e-01_r8,5.098630e-01_r8,5.115550e-01_r8,5.131280e-01_r8,5.145940e-01_r8,&
        5.159630e-01_r8,5.172460e-01_r8,5.184510e-01_r8,5.195860e-01_r8,5.206560e-01_r8/)
!BAND 5
    SSALIQ(:,5) = (/ &
        1.235130e-01_r8,2.492240e-01_r8,3.220120e-01_r8,3.673180e-01_r8,3.972600e-01_r8,&
        4.181930e-01_r8,4.335650e-01_r8,4.453200e-01_r8,4.546120e-01_r8,4.621570e-01_r8,&
        4.684220e-01_r8,4.737240e-01_r8,4.782860e-01_r8,4.822670e-01_r8,4.857830e-01_r8,&
        4.889180e-01_r8,4.917370e-01_r8,4.942870e-01_r8,4.966060e-01_r8,4.987240e-01_r8,&
        5.006640e-01_r8,5.024490e-01_r8,5.040960e-01_r8,5.056210e-01_r8,5.070370e-01_r8,&
        5.083550e-01_r8,5.095870e-01_r8,5.107420e-01_r8,5.118260e-01_r8,5.128470e-01_r8/)
!BAND 6
    SSALIQ(:,6) = (/ &
        1.588680e-01_r8,3.049070e-01_r8,3.870520e-01_r8,4.354490e-01_r8,4.647830e-01_r8,&
        4.826800e-01_r8,4.934340e-01_r8,4.996270e-01_r8,5.028930e-01_r8,5.043060e-01_r8,&
        5.045830e-01_r8,5.042010e-01_r8,5.034740e-01_r8,5.026070e-01_r8,5.017240e-01_r8,&
        5.009030e-01_r8,5.001820e-01_r8,4.995810e-01_r8,4.991030e-01_r8,4.987460e-01_r8,&
        4.984990e-01_r8,4.983520e-01_r8,4.982930e-01_r8,4.983100e-01_r8,4.983930e-01_r8,&
        4.985320e-01_r8,4.987180e-01_r8,4.989440e-01_r8,4.992020e-01_r8,4.994870e-01_r8/)
!BAND 7
    SSALIQ(:,7) = (/ &
        4.130460e-01_r8,6.041380e-01_r8,6.683640e-01_r8,6.893730e-01_r8,6.905330e-01_r8,&
        6.811520e-01_r8,6.661220e-01_r8,6.485150e-01_r8,6.303660e-01_r8,6.129690e-01_r8,&
        5.970560e-01_r8,5.829550e-01_r8,5.707270e-01_r8,5.602760e-01_r8,5.514340e-01_r8,&
        5.440020e-01_r8,5.377840e-01_r8,5.325980e-01_r8,5.282840e-01_r8,5.247000e-01_r8,&
        5.217290e-01_r8,5.192690e-01_r8,5.172360e-01_r8,5.155600e-01_r8,5.141820e-01_r8,&
        5.130540e-01_r8,5.121360e-01_r8,5.113940e-01_r8,5.108000e-01_r8,5.103320e-01_r8/)
!BAND 8
    SSALIQ(:,8) = (/ &
        5.620170e-01_r8,7.177250e-01_r8,7.536060e-01_r8,7.535810e-01_r8,7.377370e-01_r8,&
        7.143160e-01_r8,6.880930e-01_r8,6.621510e-01_r8,6.383080e-01_r8,6.174040e-01_r8,&
        5.996170e-01_r8,5.847530e-01_r8,5.724530e-01_r8,5.623210e-01_r8,5.539860e-01_r8,&
        5.471260e-01_r8,5.414730e-01_r8,5.368070e-01_r8,5.329480e-01_r8,5.297530e-01_r8,&
        5.271030e-01_r8,5.249040e-01_r8,5.230780e-01_r8,5.215620e-01_r8,5.203040e-01_r8,&
        5.192630e-01_r8,5.184030e-01_r8,5.176960e-01_r8,5.171170e-01_r8,5.166490e-01_r8/)
!BAND 9
    SSALIQ(:,9) = (/ &
        6.724270e-01_r8,7.795980e-01_r8,7.896900e-01_r8,7.709470e-01_r8,7.403710e-01_r8,&
        7.065150e-01_r8,6.743390e-01_r8,6.462250e-01_r8,6.228000e-01_r8,6.037570e-01_r8,&
        5.884350e-01_r8,5.761270e-01_r8,5.662070e-01_r8,5.581670e-01_r8,5.516060e-01_r8,&
        5.462180e-01_r8,5.417680e-01_r8,5.380750e-01_r8,5.349970e-01_r8,5.324240e-01_r8,&
        5.302690e-01_r8,5.284600e-01_r8,5.269400e-01_r8,5.256640e-01_r8,5.245930e-01_r8,&
        5.236950e-01_r8,5.229440e-01_r8,5.223190e-01_r8,5.218010e-01_r8,5.213760e-01_r8/)
!BAND 10
    SSALIQ(:,10) = (/ &
        7.270150e-01_r8,7.993900e-01_r8,7.899010e-01_r8,7.560280e-01_r8,7.150160e-01_r8,&
        6.761530e-01_r8,6.434330e-01_r8,6.174490e-01_r8,5.973030e-01_r8,5.817600e-01_r8,&
        5.697130e-01_r8,5.602950e-01_r8,5.528590e-01_r8,5.469300e-01_r8,5.421620e-01_r8,&
        5.383010e-01_r8,5.351560e-01_r8,5.325840e-01_r8,5.304750e-01_r8,5.287430e-01_r8,&
        5.273200e-01_r8,5.261530e-01_r8,5.251970e-01_r8,5.244180e-01_r8,5.237860e-01_r8,&
        5.232780e-01_r8,5.228740e-01_r8,5.225570e-01_r8,5.223150e-01_r8,5.221350e-01_r8/)
!BAND 11
    SSALIQ(:,11) = (/ &
        6.629450e-01_r8,7.311740e-01_r8,7.171340e-01_r8,6.826510e-01_r8,6.464450e-01_r8,&
        6.158400e-01_r8,5.922890e-01_r8,5.748430e-01_r8,5.620420e-01_r8,5.526090e-01_r8,&
        5.455860e-01_r8,5.402950e-01_r8,5.362610e-01_r8,5.331540e-01_r8,5.307410e-01_r8,&
        5.288530e-01_r8,5.273710e-01_r8,5.262050e-01_r8,5.252880e-01_r8,5.245690e-01_r8,&
        5.240070e-01_r8,5.235730e-01_r8,5.232430e-01_r8,5.229970e-01_r8,5.228190e-01_r8,&
        5.226990e-01_r8,5.226250e-01_r8,5.225910e-01_r8,5.225880e-01_r8,5.226120e-01_r8/)
!BAND 12
    SSALIQ(:,12) = (/ &
        9.098160e-01_r8,9.211850e-01_r8,8.975640e-01_r8,8.626740e-01_r8,8.271890e-01_r8,&
        7.958910e-01_r8,7.693990e-01_r8,7.469350e-01_r8,7.276110e-01_r8,7.107420e-01_r8,&
        6.958380e-01_r8,6.825420e-01_r8,6.705880e-01_r8,6.597770e-01_r8,6.499510e-01_r8,&
        6.409920e-01_r8,6.328050e-01_r8,6.253120e-01_r8,6.184500e-01_r8,6.121610e-01_r8,&
        6.063950e-01_r8,6.011050e-01_r8,5.962470e-01_r8,5.917820e-01_r8,5.876730e-01_r8,&
        5.838880e-01_r8,5.803950e-01_r8,5.771680e-01_r8,5.741810e-01_r8,5.714140e-01_r8/)
!BAND 13
    SSALIQ(:,13) = (/ &
        9.009350e-01_r8,8.980340e-01_r8,8.565590e-01_r8,8.080760e-01_r8,7.660510e-01_r8,&
        7.326240e-01_r8,7.059420e-01_r8,6.840930e-01_r8,6.657860e-01_r8,6.502020e-01_r8,&
        6.367950e-01_r8,6.251780e-01_r8,6.150500e-01_r8,6.061750e-01_r8,5.983630e-01_r8,&
        5.914600e-01_r8,5.853410e-01_r8,5.799050e-01_r8,5.750650e-01_r8,5.707500e-01_r8,&
        5.668970e-01_r8,5.634530e-01_r8,5.603690e-01_r8,5.576040e-01_r8,5.551220e-01_r8,&
        5.528920e-01_r8,5.508830e-01_r8,5.490730e-01_r8,5.474390e-01_r8,5.459620e-01_r8/)
!BAND 14
    SSALIQ(:,14) = (/ &
        9.391850e-01_r8,9.301700e-01_r8,8.947360e-01_r8,8.554280e-01_r8,8.215400e-01_r8,&
        7.936740e-01_r8,7.703040e-01_r8,7.501750e-01_r8,7.324800e-01_r8,7.167060e-01_r8,&
        7.025120e-01_r8,6.896550e-01_r8,6.779510e-01_r8,6.672580e-01_r8,6.574630e-01_r8,&
        6.484760e-01_r8,6.402220e-01_r8,6.326380e-01_r8,6.256700e-01_r8,6.192660e-01_r8,&
        6.133800e-01_r8,6.079670e-01_r8,6.029870e-01_r8,5.984010e-01_r8,5.941720e-01_r8,&
        5.902680e-01_r8,5.866600e-01_r8,5.833200e-01_r8,5.802230e-01_r8,5.773480e-01_r8/)
!BAND 15
    SSALIQ(:,15) = (/ &
        9.663240e-01_r8,9.554500e-01_r8,9.284420e-01_r8,9.004680e-01_r8,8.764910e-01_r8,&
        8.560280e-01_r8,8.379490e-01_r8,8.215420e-01_r8,8.064340e-01_r8,7.924280e-01_r8,&
        7.793970e-01_r8,7.672410e-01_r8,7.558680e-01_r8,7.451960e-01_r8,7.351570e-01_r8,&
        7.256970e-01_r8,7.167740e-01_r8,7.083550e-01_r8,7.004120e-01_r8,6.929190e-01_r8,&
        6.858510e-01_r8,6.791830e-01_r8,6.728900e-01_r8,6.669470e-01_r8,6.613280e-01_r8,&
        6.560110e-01_r8,6.509730e-01_r8,6.461940e-01_r8,6.416560e-01_r8,6.373400e-01_r8/)
!BAND 16
    SSALIQ(:,16) = (/ &
        9.265390e-01_r8,8.912150e-01_r8,8.530260e-01_r8,8.248650e-01_r8,8.036150e-01_r8,&
        7.862980e-01_r8,7.714110e-01_r8,7.582230e-01_r8,7.463480e-01_r8,7.355510e-01_r8,&
        7.256670e-01_r8,7.165640e-01_r8,7.081340e-01_r8,7.002840e-01_r8,6.929420e-01_r8,&
        6.860540e-01_r8,6.795770e-01_r8,6.734780e-01_r8,6.677310e-01_r8,6.623130e-01_r8,&
        6.572020e-01_r8,6.523780e-01_r8,6.478220e-01_r8,6.435140e-01_r8,6.394360e-01_r8,&
        6.355700e-01_r8,6.319000e-01_r8,6.284120e-01_r8,6.250910e-01_r8,6.219260e-01_r8/)
!BAND 1
    ASYLIQ(:,1) = (/ &
        6.361100e-02_r8,2.164160e-01_r8,3.638400e-01_r8,4.741720e-01_r8,5.533910e-01_r8,&
        6.113580e-01_r8,6.550910e-01_r8,6.890470e-01_r8,7.160740e-01_r8,7.380450e-01_r8,&
        7.562290e-01_r8,7.715130e-01_r8,7.845330e-01_r8,7.957530e-01_r8,8.055170e-01_r8,&
        8.140820e-01_r8,8.216450e-01_r8,8.283590e-01_r8,8.343470e-01_r8,8.397080e-01_r8,&
        8.445250e-01_r8,8.488670e-01_r8,8.527930e-01_r8,8.563560e-01_r8,8.596000e-01_r8,&
        8.625610e-01_r8,8.652740e-01_r8,8.677670e-01_r8,8.700640e-01_r8,8.721880e-01_r8/)
!BAND 2
    ASYLIQ(:,2) = (/ &
        1.630810e-01_r8,4.336220e-01_r8,5.996230e-01_r8,6.916060e-01_r8,7.469470e-01_r8,&
        7.832600e-01_r8,8.088230e-01_r8,8.277810e-01_r8,8.423760e-01_r8,8.539200e-01_r8,&
        8.632340e-01_r8,8.708690e-01_r8,8.772110e-01_r8,8.825400e-01_r8,8.870620e-01_r8,&
        8.909320e-01_r8,8.942700e-01_r8,8.971680e-01_r8,8.996990e-01_r8,9.019210e-01_r8,&
        9.038820e-01_r8,9.056200e-01_r8,9.071680e-01_r8,9.085520e-01_r8,9.097970e-01_r8,&
        9.109200e-01_r8,9.119380e-01_r8,9.128640e-01_r8,9.137100e-01_r8,9.144860e-01_r8/)
!BAND 3
    ASYLIQ(:,3) = (/ &
        2.452440e-01_r8,5.450590e-01_r8,6.937090e-01_r8,7.676020e-01_r8,8.098990e-01_r8,&
        8.368200e-01_r8,8.552990e-01_r8,8.686760e-01_r8,8.787380e-01_r8,8.865260e-01_r8,&
        8.926910e-01_r8,8.976650e-01_r8,9.017440e-01_r8,9.051370e-01_r8,9.079970e-01_r8,&
        9.104320e-01_r8,9.125280e-01_r8,9.143450e-01_r8,9.159320e-01_r8,9.173270e-01_r8,&
        9.185600e-01_r8,9.196560e-01_r8,9.206350e-01_r8,9.215130e-01_r8,9.223050e-01_r8,&
        9.230230e-01_r8,9.236750e-01_r8,9.242710e-01_r8,9.248170e-01_r8,9.253200e-01_r8/)
!BAND 4
    ASYLIQ(:,4) = (/ &
        3.062050e-01_r8,6.130050e-01_r8,7.482310e-01_r8,8.115950e-01_r8,8.467500e-01_r8,&
        8.686650e-01_r8,8.834710e-01_r8,8.940570e-01_r8,9.019460e-01_r8,9.080110e-01_r8,&
        9.127920e-01_r8,9.166390e-01_r8,9.197910e-01_r8,9.224130e-01_r8,9.246250e-01_r8,&
        9.265120e-01_r8,9.281390e-01_r8,9.295530e-01_r8,9.307900e-01_r8,9.318810e-01_r8,&
        9.328470e-01_r8,9.337080e-01_r8,9.344780e-01_r8,9.351710e-01_r8,9.357970e-01_r8,&
        9.363650e-01_r8,9.368820e-01_r8,9.373560e-01_r8,9.377910e-01_r8,9.381920e-01_r8/)
!BAND 5
    ASYLIQ(:,5) = (/ &
        3.648540e-01_r8,6.718040e-01_r8,7.941090e-01_r8,8.491930e-01_r8,8.791490e-01_r8,&
        8.975880e-01_r8,9.099340e-01_r8,9.187050e-01_r8,9.252120e-01_r8,9.302010e-01_r8,&
        9.341270e-01_r8,9.372850e-01_r8,9.398720e-01_r8,9.420280e-01_r8,9.438480e-01_r8,&
        9.454030e-01_r8,9.467460e-01_r8,9.479150e-01_r8,9.489400e-01_r8,9.498440e-01_r8,&
        9.506470e-01_r8,9.513640e-01_r8,9.520060e-01_r8,9.525850e-01_r8,9.531080e-01_r8,&
        9.535840e-01_r8,9.540180e-01_r8,9.544150e-01_r8,9.547810e-01_r8,9.551180e-01_r8/)
!BAND 6
    ASYLIQ(:,6) = (/ &
        4.644690e-01_r8,7.541690e-01_r8,8.525530e-01_r8,8.951070e-01_r8,9.177020e-01_r8,&
        9.313290e-01_r8,9.403070e-01_r8,9.466240e-01_r8,9.513040e-01_r8,9.549220e-01_r8,&
        9.578200e-01_r8,9.602130e-01_r8,9.622380e-01_r8,9.639890e-01_r8,9.655260e-01_r8,&
        9.668910e-01_r8,9.681150e-01_r8,9.692180e-01_r8,9.702160e-01_r8,9.711220e-01_r8,&
        9.719450e-01_r8,9.726950e-01_r8,9.733790e-01_r8,9.740030e-01_r8,9.745740e-01_r8,&
        9.750970e-01_r8,9.755770e-01_r8,9.760190e-01_r8,9.764250e-01_r8,9.768000e-01_r8/)
!BAND 7
    ASYLIQ(:,7) = (/ &
        5.573020e-01_r8,7.945420e-01_r8,8.676470e-01_r8,8.977900e-01_r8,9.125030e-01_r8,&
        9.204880e-01_r8,9.253020e-01_r8,9.286440e-01_r8,9.313840e-01_r8,9.339590e-01_r8,&
        9.365600e-01_r8,9.392390e-01_r8,9.419790e-01_r8,9.447310e-01_r8,9.474450e-01_r8,&
        9.500730e-01_r8,9.525790e-01_r8,9.549370e-01_r8,9.571340e-01_r8,9.591630e-01_r8,&
        9.610250e-01_r8,9.627250e-01_r8,9.642710e-01_r8,9.656750e-01_r8,9.669470e-01_r8,&
        9.680980e-01_r8,9.691410e-01_r8,9.700860e-01_r8,9.709430e-01_r8,9.717200e-01_r8/)
!BAND 8
    ASYLIQ(:,8) = (/ &
        6.081420e-01_r8,8.098950e-01_r8,8.677740e-01_r8,8.898190e-01_r8,8.991640e-01_r8,&
        9.035520e-01_r8,9.062960e-01_r8,9.089190e-01_r8,9.120270e-01_r8,9.157270e-01_r8,&
        9.198830e-01_r8,9.242840e-01_r8,9.287300e-01_r8,9.330650e-01_r8,9.371820e-01_r8,&
        9.410200e-01_r8,9.445470e-01_r8,9.477550e-01_r8,9.506520e-01_r8,9.532540e-01_r8,&
        9.555820e-01_r8,9.576600e-01_r8,9.595150e-01_r8,9.611680e-01_r8,9.626430e-01_r8,&
        9.639610e-01_r8,9.651390e-01_r8,9.661950e-01_r8,9.671430e-01_r8,9.679960e-01_r8/)
!BAND 9
    ASYLIQ(:,9) = (/ &
        6.495040e-01_r8,8.221570e-01_r8,8.654680e-01_r8,8.788560e-01_r8,8.832110e-01_r8,&
        8.857420e-01_r8,8.891590e-01_r8,8.940770e-01_r8,9.001690e-01_r8,9.068490e-01_r8,&
        9.135990e-01_r8,9.200640e-01_r8,9.260440e-01_r8,9.314560e-01_r8,9.362840e-01_r8,&
        9.405530e-01_r8,9.443090e-01_r8,9.476030e-01_r8,9.504880e-01_r8,9.530140e-01_r8,&
        9.552280e-01_r8,9.571690e-01_r8,9.588750e-01_r8,9.603780e-01_r8,9.617050e-01_r8,&
        9.628810e-01_r8,9.639240e-01_r8,9.648540e-01_r8,9.656840e-01_r8,9.664280e-01_r8/)
!BAND 10
    ASYLIQ(:,10) = (/ &
        6.777400e-01_r8,8.294410e-01_r8,8.623000e-01_r8,8.700110e-01_r8,8.728460e-01_r8,&
        8.770910e-01_r8,8.839030e-01_r8,8.924460e-01_r8,9.015510e-01_r8,9.103580e-01_r8,&
        9.184010e-01_r8,9.255050e-01_r8,9.316610e-01_r8,9.369430e-01_r8,9.414560e-01_r8,&
        9.453080e-01_r8,9.485990e-01_r8,9.514150e-01_r8,9.538330e-01_r8,9.559140e-01_r8,&
        9.577120e-01_r8,9.592700e-01_r8,9.606260e-01_r8,9.618090e-01_r8,9.628470e-01_r8,&
        9.637590e-01_r8,9.645650e-01_r8,9.652790e-01_r8,9.659140e-01_r8,9.664810e-01_r8/)
!BAND 11
    ASYLIQ(:,11) = (/ &
        7.161050e-01_r8,8.481250e-01_r8,8.765770e-01_r8,8.858100e-01_r8,8.921160e-01_r8,&
        8.993790e-01_r8,9.076200e-01_r8,9.159760e-01_r8,9.237640e-01_r8,9.306500e-01_r8,&
        9.365590e-01_r8,9.415480e-01_r8,9.457330e-01_r8,9.492400e-01_r8,9.521830e-01_r8,&
        9.546640e-01_r8,9.567640e-01_r8,9.585510e-01_r8,9.600790e-01_r8,9.613930e-01_r8,&
        9.625270e-01_r8,9.635110e-01_r8,9.643700e-01_r8,9.651220e-01_r8,9.657840e-01_r8,&
        9.663690e-01_r8,9.668880e-01_r8,9.673510e-01_r8,9.677650e-01_r8,9.681380e-01_r8/)
!BAND 12
    ASYLIQ(:,12) = (/ &
        7.669540e-01_r8,8.461300e-01_r8,8.475170e-01_r8,8.390320e-01_r8,8.369280e-01_r8,&
        8.423430e-01_r8,8.519440e-01_r8,8.628480e-01_r8,8.734790e-01_r8,8.831800e-01_r8,&
        8.917750e-01_r8,8.993060e-01_r8,9.058940e-01_r8,9.116750e-01_r8,9.167760e-01_r8,&
        9.213020e-01_r8,9.253420e-01_r8,9.289660e-01_r8,9.322300e-01_r8,9.351810e-01_r8,&
        9.378570e-01_r8,9.402910e-01_r8,9.425110e-01_r8,9.445410e-01_r8,9.464010e-01_r8,&
        9.481100e-01_r8,9.496830e-01_r8,9.511350e-01_r8,9.524780e-01_r8,9.537220e-01_r8/)
!BAND 13
    ASYLIQ(:,13) = (/ &
        7.842630e-01_r8,8.421260e-01_r8,8.359550e-01_r8,8.319060e-01_r8,8.393510e-01_r8,&
        8.529060e-01_r8,8.674980e-01_r8,8.808610e-01_r8,8.923860e-01_r8,9.021260e-01_r8,&
        9.103320e-01_r8,9.172760e-01_r8,9.231950e-01_r8,9.282840e-01_r8,9.326920e-01_r8,&
        9.365370e-01_r8,9.399100e-01_r8,9.428840e-01_r8,9.455170e-01_r8,9.478550e-01_r8,&
        9.499380e-01_r8,9.517990e-01_r8,9.534660e-01_r8,9.549620e-01_r8,9.563090e-01_r8,&
        9.575240e-01_r8,9.586240e-01_r8,9.596210e-01_r8,9.605270e-01_r8,9.613520e-01_r8/)
!BAND 14
    ASYLIQ(:,14) = (/ &
        7.854430e-01_r8,8.309500e-01_r8,8.164570e-01_r8,8.122940e-01_r8,8.221140e-01_r8,&
        8.371240e-01_r8,8.520970e-01_r8,8.653110e-01_r8,8.765320e-01_r8,8.859890e-01_r8,&
        8.939990e-01_r8,9.008520e-01_r8,9.067840e-01_r8,9.119780e-01_r8,9.165730e-01_r8,&
        9.206740e-01_r8,9.243600e-01_r8,9.276910e-01_r8,9.307140e-01_r8,9.334670e-01_r8,&
        9.359800e-01_r8,9.382790e-01_r8,9.403880e-01_r8,9.423250e-01_r8,9.441090e-01_r8,&
        9.457550e-01_r8,9.472750e-01_r8,9.486830e-01_r8,9.499900e-01_r8,9.512030e-01_r8/)
!BAND 15
    ASYLIQ(:,15) = (/ &
        7.834020e-01_r8,8.117940e-01_r8,7.910930e-01_r8,7.910280e-01_r8,8.049980e-01_r8,&
        8.214430e-01_r8,8.360980e-01_r8,8.482930e-01_r8,8.583740e-01_r8,8.668210e-01_r8,&
        8.740330e-01_r8,8.803000e-01_r8,8.858310e-01_r8,8.907770e-01_r8,8.952470e-01_r8,&
        8.993180e-01_r8,9.030510e-01_r8,9.064900e-01_r8,9.096680e-01_r8,9.126140e-01_r8,&
        9.153520e-01_r8,9.179010e-01_r8,9.202800e-01_r8,9.225030e-01_r8,9.245860e-01_r8,&
        9.265410e-01_r8,9.283800e-01_r8,9.301120e-01_r8,9.317480e-01_r8,9.332950e-01_r8/)
!BAND 16
    ASYLIQ(:,16) = (/ &
        7.825820e-01_r8,7.913460e-01_r8,7.884810e-01_r8,8.065030e-01_r8,8.271350e-01_r8,&
        8.441100e-01_r8,8.573200e-01_r8,8.677240e-01_r8,8.761400e-01_r8,8.831260e-01_r8,&
        8.890490e-01_r8,8.941590e-01_r8,8.986300e-01_r8,9.025900e-01_r8,9.061320e-01_r8,&
        9.093280e-01_r8,9.122300e-01_r8,9.148810e-01_r8,9.173130e-01_r8,9.195520e-01_r8,&
        9.216210e-01_r8,9.235380e-01_r8,9.253200e-01_r8,9.269800e-01_r8,9.285320e-01_r8,&
        9.299860e-01_r8,9.313520e-01_r8,9.326390e-01_r8,9.338540e-01_r8,9.350030e-01_r8/)
! --------------linhan add-------------------------------------
      SSAICE2(:,1) = (/ &
!    BAND 1
       2.384497e-01_r8,2.909285e-01_r8,3.267788e-01_r8,3.540130e-01_r8,3.759747e-01_r8, &
       3.943837e-01_r8,4.102388e-01_r8,4.241709e-01_r8,4.366037e-01_r8,4.478358e-01_r8, &
       4.580852e-01_r8,4.675163e-01_r8,4.762559e-01_r8,4.844041e-01_r8,4.920411e-01_r8, &
       4.992324e-01_r8,5.060323e-01_r8,5.124859e-01_r8,5.186316e-01_r8,5.245023e-01_r8, &
       5.301260e-01_r8,5.355275e-01_r8,5.407281e-01_r8,5.457469e-01_r8,5.506008e-01_r8, &
       5.553048e-01_r8,5.598726e-01_r8,5.643165e-01_r8,5.686477e-01_r8,5.728765e-01_r8, &
       5.770124e-01_r8,5.810643e-01_r8,5.850404e-01_r8,5.889484e-01_r8,5.927958e-01_r8, &
       5.965895e-01_r8,6.003362e-01_r8,6.040423e-01_r8,6.077142e-01_r8,6.113580e-01_r8, &
       6.149798e-01_r8,6.185855e-01_r8,6.221813e-01_r8/)
      SSAICE2(:,2) = (/ &
!    BAND 2
       8.221222e-01_r8,7.943076e-01_r8,7.738458e-01_r8,7.572275e-01_r8,7.430030e-01_r8, &
       7.304254e-01_r8,7.190571e-01_r8,7.086179e-01_r8,6.989172e-01_r8,6.898189e-01_r8, &
       6.812222e-01_r8,6.730501e-01_r8,6.652427e-01_r8,6.577517e-01_r8,6.505382e-01_r8, &
       6.435701e-01_r8,6.368202e-01_r8,6.302659e-01_r8,6.238876e-01_r8,6.176684e-01_r8, &
       6.115934e-01_r8,6.056497e-01_r8,5.998256e-01_r8,5.941108e-01_r8,5.884958e-01_r8, &
       5.829720e-01_r8,5.775313e-01_r8,5.721662e-01_r8,5.668698e-01_r8,5.616352e-01_r8, &
       5.564557e-01_r8,5.513249e-01_r8,5.462362e-01_r8,5.411828e-01_r8,5.361576e-01_r8, &
       5.311530e-01_r8,5.261609e-01_r8,5.211719e-01_r8,5.161756e-01_r8,5.111597e-01_r8, &
       5.061093e-01_r8,5.010065e-01_r8,4.958288e-01_r8/)
      SSAICE2(:,3) = (/ &
!    BAND 3
       6.411479e-01_r8,6.217354e-01_r8,6.080982e-01_r8,5.975189e-01_r8,5.888491e-01_r8, &
       5.814910e-01_r8,5.750920e-01_r8,5.694265e-01_r8,5.643407e-01_r8,5.597253e-01_r8, &
       5.554996e-01_r8,5.516022e-01_r8,5.479855e-01_r8,5.446115e-01_r8,5.414499e-01_r8, &
       5.384755e-01_r8,5.356676e-01_r8,5.330090e-01_r8,5.304849e-01_r8,5.280827e-01_r8, &
       5.257918e-01_r8,5.236027e-01_r8,5.215074e-01_r8,5.194988e-01_r8,5.175707e-01_r8, &
       5.157175e-01_r8,5.139344e-01_r8,5.122171e-01_r8,5.105619e-01_r8,5.089654e-01_r8, & 
       5.074245e-01_r8,5.059368e-01_r8,5.044998e-01_r8,5.031117e-01_r8,5.017706e-01_r8, &
       5.004753e-01_r8,4.992246e-01_r8,4.980175e-01_r8,4.968537e-01_r8,4.957327e-01_r8, &
       4.946547e-01_r8,4.936202e-01_r8,4.926301e-01_r8/)
      SSAICE2(:,4) = (/ &
!    BAND 4
       5.308658e-01_r8,5.286308e-01_r8,5.270808e-01_r8,5.259007e-01_r8,5.249552e-01_r8, &
       5.241732e-01_r8,5.235120e-01_r8,5.229445e-01_r8,5.224518e-01_r8,5.220204e-01_r8, & 
       5.216405e-01_r8,5.213044e-01_r8,5.210062e-01_r8,5.207411e-01_r8,5.205054e-01_r8, &
       5.202958e-01_r8,5.201099e-01_r8,5.199454e-01_r8,5.198004e-01_r8,5.196735e-01_r8, &
       5.195632e-01_r8,5.194684e-01_r8,5.193881e-01_r8,5.193214e-01_r8,5.192675e-01_r8, &
       5.192259e-01_r8,5.191958e-01_r8,5.191767e-01_r8,5.191683e-01_r8,5.191702e-01_r8, &
       5.191819e-01_r8,5.192032e-01_r8,5.192338e-01_r8,5.192736e-01_r8,5.193222e-01_r8, &
       5.193797e-01_r8,5.194458e-01_r8,5.195204e-01_r8,5.196036e-01_r8,5.196951e-01_r8, &
       5.197952e-01_r8,5.199036e-01_r8,5.200205e-01_r8/)
      SSAICE2(:,5) = (/ &
!    BAND 5
       4.443291e-01_r8,4.591129e-01_r8,4.693524e-01_r8,4.772410e-01_r8,4.836841e-01_r8, &
       4.891456e-01_r8,4.938956e-01_r8,4.981055e-01_r8,5.018908e-01_r8,5.053336e-01_r8, &
       5.084939e-01_r8,5.114170e-01_r8,5.141381e-01_r8,5.166852e-01_r8,5.190805e-01_r8, &
       5.213424e-01_r8,5.234860e-01_r8,5.255239e-01_r8,5.274670e-01_r8,5.293243e-01_r8, &
       5.311038e-01_r8,5.328121e-01_r8,5.344554e-01_r8,5.360388e-01_r8,5.375669e-01_r8, &
       5.390438e-01_r8,5.404732e-01_r8,5.418583e-01_r8,5.432020e-01_r8,5.445071e-01_r8, &
       5.457758e-01_r8,5.470105e-01_r8,5.482129e-01_r8,5.493851e-01_r8,5.505286e-01_r8, &
       5.516449e-01_r8,5.527355e-01_r8,5.538017e-01_r8,5.548446e-01_r8,5.558653e-01_r8, &
       5.568650e-01_r8,5.578445e-01_r8,5.588047e-01_r8/)
      SSAICE2(:,6) = (/ &
!    BAND 6
       3.723265e-01_r8,3.996998e-01_r8,4.181439e-01_r8,4.320349e-01_r8,4.431625e-01_r8, &
       4.524352e-01_r8,4.603781e-01_r8,4.673218e-01_r8,4.734879e-01_r8,4.790323e-01_r8, &
       4.840687e-01_r8,4.886825e-01_r8,4.929397e-01_r8,4.968922e-01_r8,5.005815e-01_r8, &
       5.040415e-01_r8,5.073000e-01_r8,5.103805e-01_r8,5.133026e-01_r8,5.160831e-01_r8, &
       5.187363e-01_r8,5.212749e-01_r8,5.237097e-01_r8,5.260504e-01_r8,5.283055e-01_r8, &
       5.304825e-01_r8,5.325884e-01_r8,5.346293e-01_r8,5.366107e-01_r8,5.385378e-01_r8, &
       5.404153e-01_r8,5.422477e-01_r8,5.440388e-01_r8,5.457926e-01_r8,5.475127e-01_r8, &
       5.492024e-01_r8,5.508652e-01_r8,5.525042e-01_r8,5.541224e-01_r8,5.557230e-01_r8, &
       5.573090e-01_r8,5.588835e-01_r8,5.604495e-01_r8/)
      SSAICE2(:,7) = (/ &
!    BAND 7
       6.749288e-01_r8,6.475680e-01_r8,6.291381e-01_r8,6.152112e-01_r8,6.040010e-01_r8, &
       5.946084e-01_r8,5.865167e-01_r8,5.794023e-01_r8,5.730486e-01_r8,5.673037e-01_r8, &
       5.620572e-01_r8,5.572259e-01_r8,5.527461e-01_r8,5.485674e-01_r8,5.446498e-01_r8, &
       5.409604e-01_r8,5.374723e-01_r8,5.341632e-01_r8,5.310141e-01_r8,5.280089e-01_r8, &
       5.251338e-01_r8,5.223770e-01_r8,5.197280e-01_r8,5.171778e-01_r8,5.147184e-01_r8, &
       5.123427e-01_r8,5.100445e-01_r8,5.078180e-01_r8,5.056582e-01_r8,5.035605e-01_r8, &
       5.015208e-01_r8,4.995353e-01_r8,4.976005e-01_r8,4.957133e-01_r8,4.938707e-01_r8, &
       4.920701e-01_r8,4.903088e-01_r8,4.885844e-01_r8,4.868948e-01_r8,4.852377e-01_r8, &
       4.836110e-01_r8,4.820127e-01_r8,4.804408e-01_r8/)
      SSAICE2(:,8) = (/ &
!    BAND 8
       7.148414e-01_r8,6.801247e-01_r8,6.565983e-01_r8,6.387794e-01_r8,6.244318e-01_r8, &
       6.124207e-01_r8,6.020902e-01_r8,5.930271e-01_r8,5.849538e-01_r8,5.776752e-01_r8, &
       5.710486e-01_r8,5.649666e-01_r8,5.593463e-01_r8,5.541225e-01_r8,5.492428e-01_r8, &
       5.446646e-01_r8,5.403528e-01_r8,5.362779e-01_r8,5.324151e-01_r8,5.287436e-01_r8, &
       5.252450e-01_r8,5.219039e-01_r8,5.187066e-01_r8,5.156411e-01_r8,5.126970e-01_r8, &
       5.098650e-01_r8,5.071367e-01_r8,5.045048e-01_r8,5.019627e-01_r8,4.995044e-01_r8, &
       4.971244e-01_r8,4.948179e-01_r8,4.925804e-01_r8,4.904078e-01_r8,4.882964e-01_r8, &
       4.862427e-01_r8,4.842437e-01_r8,4.822964e-01_r8,4.803980e-01_r8,4.785462e-01_r8, &
       4.767386e-01_r8,4.749730e-01_r8,4.732474e-01_r8/)
      SSAICE2(:,9) = (/ &
!    BAND 9
       6.712279e-01_r8,6.442293e-01_r8,6.257659e-01_r8,6.116928e-01_r8,6.003067e-01_r8, &
       5.907386e-01_r8,5.824839e-01_r8,5.752233e-01_r8,5.687419e-01_r8,5.628879e-01_r8, &
       5.575502e-01_r8,5.526449e-01_r8,5.481072e-01_r8,5.438859e-01_r8,5.399399e-01_r8, &
       5.362355e-01_r8,5.327451e-01_r8,5.294455e-01_r8,5.263172e-01_r8,5.233434e-01_r8, &
       5.205098e-01_r8,5.178041e-01_r8,5.152154e-01_r8,5.127343e-01_r8,5.103524e-01_r8, &
       5.080623e-01_r8,5.058575e-01_r8,5.037321e-01_r8,5.016807e-01_r8,4.996987e-01_r8, &
       4.977816e-01_r8,4.959258e-01_r8,4.941275e-01_r8,4.923836e-01_r8,4.906910e-01_r8, &
       4.890472e-01_r8,4.874496e-01_r8,4.858958e-01_r8,4.843839e-01_r8,4.829118e-01_r8, &
       4.814778e-01_r8,4.800801e-01_r8,4.787172e-01_r8/)
      SSAICE2(:,10) = (/ &
!    BAND 10
       6.254590e-01_r8,6.055970e-01_r8,5.921137e-01_r8,5.818892e-01_r8,5.736492e-01_r8, &
       5.667460e-01_r8,5.608054e-01_r8,5.555910e-01_r8,5.509442e-01_r8,5.467533e-01_r8, &
       5.429368e-01_r8,5.394330e-01_r8,5.361945e-01_r8,5.331840e-01_r8,5.303713e-01_r8, &
       5.277321e-01_r8,5.252462e-01_r8,5.228967e-01_r8,5.206694e-01_r8,5.185522e-01_r8, &
       5.165348e-01_r8,5.146082e-01_r8,5.127644e-01_r8,5.109968e-01_r8,5.092993e-01_r8, &
       5.076665e-01_r8,5.060936e-01_r8,5.045765e-01_r8,5.031113e-01_r8,5.016946e-01_r8, &
       5.003232e-01_r8,4.989945e-01_r8,4.977057e-01_r8,4.964546e-01_r8,4.952390e-01_r8, & 
       4.940570e-01_r8,4.929068e-01_r8,4.917867e-01_r8,4.906952e-01_r8,4.896308e-01_r8, &
       4.885923e-01_r8,4.875785e-01_r8,4.865881e-01_r8/)
      SSAICE2(:,11) = (/ &
!    BAND 11
       6.232263e-01_r8,6.037961e-01_r8,5.906828e-01_r8,5.807779e-01_r8,5.728182e-01_r8, &
       5.661645e-01_r8,5.604483e-01_r8,5.554375e-01_r8,5.509768e-01_r8,5.469570e-01_r8, &
       5.432984e-01_r8,5.399411e-01_r8,5.368390e-01_r8,5.339556e-01_r8,5.312619e-01_r8, &
       5.287341e-01_r8,5.263529e-01_r8,5.241018e-01_r8,5.219671e-01_r8,5.199373e-01_r8, &
       5.180022e-01_r8,5.161534e-01_r8,5.143831e-01_r8,5.126849e-01_r8,5.110531e-01_r8, &
       5.094823e-01_r8,5.079682e-01_r8,5.065066e-01_r8,5.050939e-01_r8,5.037269e-01_r8, &
       5.024025e-01_r8,5.011181e-01_r8,4.998713e-01_r8,4.986598e-01_r8,4.974815e-01_r8, &
       4.963348e-01_r8,4.952177e-01_r8,4.941288e-01_r8,4.930666e-01_r8,4.920297e-01_r8, &
       4.910170e-01_r8,4.900272e-01_r8,4.890593e-01_r8/)
      SSAICE2(:,12) = (/ &
!    BAND 12
       8.165189e-01_r8,7.690707e-01_r8,7.369135e-01_r8,7.125590e-01_r8,6.929511e-01_r8, &
       6.765386e-01_r8,6.624248e-01_r8,6.500445e-01_r8,6.390180e-01_r8,6.290783e-01_r8, &
       6.200302e-01_r8,6.117269e-01_r8,6.040549e-01_r8,5.969249e-01_r8,5.902653e-01_r8, &
       5.840177e-01_r8,5.781341e-01_r8,5.725743e-01_r8,5.673045e-01_r8,5.622957e-01_r8, &
       5.575233e-01_r8,5.529660e-01_r8,5.486050e-01_r8,5.444242e-01_r8,5.404091e-01_r8, &
       5.365471e-01_r8,5.328269e-01_r8,5.292383e-01_r8,5.257724e-01_r8,5.224210e-01_r8, &
       5.191768e-01_r8,5.160329e-01_r8,5.129835e-01_r8,5.100229e-01_r8,5.071461e-01_r8, &
       5.043485e-01_r8,5.016258e-01_r8,4.989740e-01_r8,4.963895e-01_r8,4.938691e-01_r8, &
       4.914094e-01_r8,4.890078e-01_r8,4.866614e-01_r8/)
      SSAICE2(:,13) = (/ &
!    BAND 13
       7.081550e-01_r8,6.764607e-01_r8,6.549873e-01_r8,6.387269e-01_r8,6.256366e-01_r8, &
       6.146798e-01_r8,6.052573e-01_r8,5.969915e-01_r8,5.896290e-01_r8,5.829913e-01_r8, &
       5.769484e-01_r8,5.714020e-01_r8,5.662765e-01_r8,5.615123e-01_r8,5.570617e-01_r8, &
       5.528856e-01_r8,5.489522e-01_r8,5.452345e-01_r8,5.417100e-01_r8,5.383595e-01_r8, &
       5.351665e-01_r8,5.321168e-01_r8,5.291979e-01_r8,5.263991e-01_r8,5.237107e-01_r8, &
       5.211244e-01_r8,5.186325e-01_r8,5.162284e-01_r8,5.139061e-01_r8,5.116601e-01_r8, &
       5.094855e-01_r8,5.073778e-01_r8,5.053332e-01_r8,5.033477e-01_r8,5.014182e-01_r8, &
       4.995414e-01_r8,4.977147e-01_r8,4.959352e-01_r8,4.942007e-01_r8,4.925089e-01_r8, &
       4.908577e-01_r8,4.892452e-01_r8,4.876697e-01_r8/)
      SSAICE2(:,14) = (/ &
!    BAND 14
       7.353033e-01_r8,6.997772e-01_r8,6.756819e-01_r8,6.574216e-01_r8,6.427122e-01_r8, &
       6.303940e-01_r8,6.197965e-01_r8,6.104970e-01_r8,6.022115e-01_r8,5.947404e-01_r8, &
       5.879375e-01_r8,5.816930e-01_r8,5.759219e-01_r8,5.705575e-01_r8,5.655461e-01_r8, &
       5.608441e-01_r8,5.564153e-01_r8,5.522298e-01_r8,5.482621e-01_r8,5.444907e-01_r8, &
       5.408970e-01_r8,5.374650e-01_r8,5.341808e-01_r8,5.310321e-01_r8,5.280082e-01_r8, &
       5.250995e-01_r8,5.222976e-01_r8,5.195949e-01_r8,5.169845e-01_r8,5.144605e-01_r8, &
       5.120172e-01_r8,5.096496e-01_r8,5.073532e-01_r8,5.051238e-01_r8,5.029576e-01_r8, &
       5.008511e-01_r8,4.988011e-01_r8,4.968047e-01_r8,4.948591e-01_r8,4.929617e-01_r8, &
       4.911103e-01_r8,4.893027e-01_r8,4.875368e-01_r8/)
      SSAICE2(:,15) = (/ &
!    BAND 15
       8.248475e-01_r8,7.814674e-01_r8,7.518947e-01_r8,7.294008e-01_r8,7.112284e-01_r8, &
       6.959732e-01_r8,6.828218e-01_r8,6.712600e-01_r8,6.609423e-01_r8,6.516250e-01_r8, &
       6.431297e-01_r8,6.353220e-01_r8,6.280981e-01_r8,6.213760e-01_r8,6.150901e-01_r8, &
       6.091866e-01_r8,6.036215e-01_r8,5.983576e-01_r8,5.933637e-01_r8,5.886134e-01_r8, &
       5.840837e-01_r8,5.797549e-01_r8,5.756098e-01_r8,5.716333e-01_r8,5.678121e-01_r8, &
       5.641345e-01_r8,5.605899e-01_r8,5.571691e-01_r8,5.538635e-01_r8,5.506656e-01_r8, &
       5.475687e-01_r8,5.445663e-01_r8,5.416530e-01_r8,5.388234e-01_r8,5.360730e-01_r8, &
       5.333974e-01_r8,5.307925e-01_r8,5.282548e-01_r8,5.257807e-01_r8,5.233673e-01_r8, &
       5.210114e-01_r8,5.187106e-01_r8,5.164621e-01_r8/)
      SSAICE2(:,16) = (/ &
!    BAND 16
       6.630615e-01_r8,6.451169e-01_r8,6.333696e-01_r8,6.246927e-01_r8,6.178420e-01_r8, &
       6.121976e-01_r8,6.074069e-01_r8,6.032505e-01_r8,5.995830e-01_r8,5.963030e-01_r8, &
       5.933372e-01_r8,5.906311e-01_r8,5.881427e-01_r8,5.858395e-01_r8,5.836955e-01_r8, &
       5.816896e-01_r8,5.798046e-01_r8,5.780264e-01_r8,5.763429e-01_r8,5.747441e-01_r8, &
       5.732213e-01_r8,5.717672e-01_r8,5.703754e-01_r8,5.690403e-01_r8,5.677571e-01_r8, &
       5.665215e-01_r8,5.653297e-01_r8,5.641782e-01_r8,5.630643e-01_r8,5.619850e-01_r8, &
       5.609381e-01_r8,5.599214e-01_r8,5.589328e-01_r8,5.579707e-01_r8,5.570333e-01_r8, &
       5.561193e-01_r8,5.552272e-01_r8,5.543558e-01_r8,5.535041e-01_r8,5.526708e-01_r8, &
       5.518551e-01_r8,5.510561e-01_r8,5.502729e-01_r8/)

!     ASYMMETRY FACTOR: Unitless
      ASYICE2(:,1) = (/ &
!    BAND 1
       2.255639e-01_r8,4.645627e-01_r8,5.756219e-01_r8,6.411863e-01_r8,6.850450e-01_r8, &
       7.167515e-01_r8,7.409211e-01_r8,7.600705e-01_r8,7.756950e-01_r8,7.887423e-01_r8, &
       7.998436e-01_r8,8.094368e-01_r8,8.178357e-01_r8,8.252715e-01_r8,8.319187e-01_r8, &
       8.379115e-01_r8,8.433551e-01_r8,8.483330e-01_r8,8.529127e-01_r8,8.571493e-01_r8, &
       8.610881e-01_r8,8.647670e-01_r8,8.682179e-01_r8,8.714677e-01_r8,8.745396e-01_r8, &
       8.774535e-01_r8,8.802267e-01_r8,8.828741e-01_r8,8.854091e-01_r8,8.878434e-01_r8, &
       8.901874e-01_r8,8.924504e-01_r8,8.946407e-01_r8,8.967659e-01_r8,8.988330e-01_r8, & 
       9.008482e-01_r8,9.028173e-01_r8,9.047457e-01_r8,9.066384e-01_r8,9.085002e-01_r8, &
       9.103353e-01_r8,9.121481e-01_r8,9.139426e-01_r8/)
      ASYICE2(:,2) = (/ &
!    BAND 2
       5.393286e-01_r8,6.558766e-01_r8,7.164199e-01_r8,7.545486e-01_r8,7.811948e-01_r8, &
       8.010778e-01_r8,8.165998e-01_r8,8.291243e-01_r8,8.394884e-01_r8,8.482370e-01_r8, &
       8.557418e-01_r8,8.622656e-01_r8,8.680002e-01_r8,8.730891e-01_r8,8.776420e-01_r8, &
       8.817442e-01_r8,8.854634e-01_r8,8.888538e-01_r8,8.919594e-01_r8,8.948166e-01_r8, &
       8.974552e-01_r8,8.999005e-01_r8,9.021735e-01_r8,9.042922e-01_r8,9.062719e-01_r8, &
       9.081256e-01_r8,9.098647e-01_r8,9.114988e-01_r8,9.130362e-01_r8,9.144842e-01_r8, &
       9.158490e-01_r8,9.171357e-01_r8,9.183488e-01_r8,9.194918e-01_r8,9.205677e-01_r8, &
       9.215785e-01_r8,9.225256e-01_r8,9.234093e-01_r8,9.242292e-01_r8,9.249837e-01_r8, &
       9.256698e-01_r8,9.262828e-01_r8,9.268157e-01_r8/)
      ASYICE2(:,3) = (/ &
!    BAND 3
       6.402550e-01_r8,7.366100e-01_r8,7.861283e-01_r8,8.170823e-01_r8,8.385946e-01_r8, &
       8.545773e-01_r8,8.670110e-01_r8,8.770146e-01_r8,8.852724e-01_r8,8.922285e-01_r8, &
       8.981848e-01_r8,9.033542e-01_r8,9.078918e-01_r8,9.119133e-01_r8,9.155069e-01_r8, &
       9.187415e-01_r8,9.216713e-01_r8,9.243398e-01_r8,9.267824e-01_r8,9.290280e-01_r8, &
       9.311009e-01_r8,9.330210e-01_r8,9.348055e-01_r8,9.364687e-01_r8,9.380231e-01_r8, &
       9.394792e-01_r8,9.408463e-01_r8,9.421324e-01_r8,9.433444e-01_r8,9.444886e-01_r8, &
       9.455702e-01_r8,9.465940e-01_r8,9.475642e-01_r8,9.484844e-01_r8,9.493581e-01_r8, &
       9.501880e-01_r8,9.509766e-01_r8,9.517263e-01_r8,9.524388e-01_r8,9.531158e-01_r8, &
       9.537587e-01_r8,9.543686e-01_r8,9.549462e-01_r8/)
      ASYICE2(:,4) = (/ &
!    BAND 4
       6.868425e-01_r8,7.885874e-01_r8,8.342997e-01_r8,8.602518e-01_r8,8.769749e-01_r8, &
       8.886475e-01_r8,8.972569e-01_r8,9.038686e-01_r8,9.091054e-01_r8,9.133553e-01_r8, &
       9.168731e-01_r8,9.198324e-01_r8,9.223563e-01_r8,9.245338e-01_r8,9.264314e-01_r8, &
       9.280994e-01_r8,9.295769e-01_r8,9.308944e-01_r8,9.320762e-01_r8,9.331421e-01_r8, &
       9.341079e-01_r8,9.349870e-01_r8,9.357901e-01_r8,9.365266e-01_r8,9.372040e-01_r8, &
       9.378290e-01_r8,9.384073e-01_r8,9.389435e-01_r8,9.394420e-01_r8,9.399062e-01_r8, &
       9.403395e-01_r8,9.407445e-01_r8,9.411237e-01_r8,9.414794e-01_r8,9.418132e-01_r8, &
       9.421271e-01_r8,9.424225e-01_r8,9.427007e-01_r8,9.429630e-01_r8,9.432104e-01_r8, &
       9.434440e-01_r8,9.436646e-01_r8,9.438730e-01_r8/)
      ASYICE2(:,5) = (/ &
!    BAND 5
       7.273309e-01_r8,8.266992e-01_r8,8.655715e-01_r8,8.855658e-01_r8,8.974914e-01_r8, &
       9.053017e-01_r8,9.107583e-01_r8,9.147555e-01_r8,9.177919e-01_r8,9.201655e-01_r8, &
       9.220647e-01_r8,9.236137e-01_r8,9.248978e-01_r8,9.259770e-01_r8,9.268948e-01_r8, &
       9.276835e-01_r8,9.283676e-01_r8,9.289656e-01_r8,9.294923e-01_r8,9.299591e-01_r8, &
       9.303752e-01_r8,9.307482e-01_r8,9.310841e-01_r8,9.313880e-01_r8,9.316640e-01_r8, &
       9.319156e-01_r8,9.321458e-01_r8,9.323571e-01_r8,9.325516e-01_r8,9.327311e-01_r8, &
       9.328973e-01_r8,9.330515e-01_r8,9.331948e-01_r8,9.333284e-01_r8,9.334532e-01_r8, &
       9.335698e-01_r8,9.336792e-01_r8,9.337819e-01_r8,9.338784e-01_r8,9.339694e-01_r8, &
       9.340551e-01_r8,9.341361e-01_r8,9.342127e-01_r8/)
      ASYICE2(:,6) = (/ &
!    BAND 6
       7.887503e-01_r8,8.770704e-01_r8,9.104538e-01_r8,9.272855e-01_r8,9.371982e-01_r8, &
       9.436352e-01_r8,9.481052e-01_r8,9.513645e-01_r8,9.538305e-01_r8,9.557506e-01_r8, &
       9.572803e-01_r8,9.585216e-01_r8,9.595441e-01_r8,9.603968e-01_r8,9.611150e-01_r8, &
       9.617249e-01_r8,9.622461e-01_r8,9.626938e-01_r8,9.630796e-01_r8,9.634129e-01_r8, &
       9.637009e-01_r8,9.639497e-01_r8,9.641639e-01_r8,9.643477e-01_r8,9.645041e-01_r8, &
       9.646359e-01_r8,9.647453e-01_r8,9.648341e-01_r8,9.649039e-01_r8,9.649559e-01_r8, &
       9.649912e-01_r8,9.650105e-01_r8,9.650147e-01_r8,9.650041e-01_r8,9.649792e-01_r8, &
       9.649403e-01_r8,9.648876e-01_r8,9.648210e-01_r8,9.647406e-01_r8,9.646461e-01_r8, &
       9.645374e-01_r8,9.644141e-01_r8,9.642758e-01_r8/)
      ASYICE2(:,7) = (/ &
!    BAND 7
       8.378017e-01_r8,8.852928e-01_r8,9.090067e-01_r8,9.234678e-01_r8,9.333026e-01_r8, &
       9.404700e-01_r8,9.459500e-01_r8,9.502900e-01_r8,9.538212e-01_r8,9.567564e-01_r8, &
       9.592387e-01_r8,9.613686e-01_r8,9.632181e-01_r8,9.648409e-01_r8,9.662775e-01_r8, &
       9.675592e-01_r8,9.687106e-01_r8,9.697512e-01_r8,9.706969e-01_r8,9.715604e-01_r8, &
       9.723525e-01_r8,9.730820e-01_r8,9.737564e-01_r8,9.743818e-01_r8,9.749637e-01_r8, &
       9.755068e-01_r8,9.760148e-01_r8,9.764914e-01_r8,9.769395e-01_r8,9.773618e-01_r8, &
       9.777606e-01_r8,9.781380e-01_r8,9.784958e-01_r8,9.788357e-01_r8,9.791591e-01_r8, &
       9.794674e-01_r8,9.797619e-01_r8,9.800437e-01_r8,9.803137e-01_r8,9.805730e-01_r8, &
       9.808224e-01_r8,9.810629e-01_r8,9.812952e-01_r8/)
      ASYICE2(:,8) = (/ &
!    BAND 8
       8.410085e-01_r8,8.742948e-01_r8,8.935566e-01_r8,9.065956e-01_r8,9.162142e-01_r8, &
       9.237079e-01_r8,9.297715e-01_r8,9.348164e-01_r8,9.391043e-01_r8,9.428109e-01_r8, &
       9.460592e-01_r8,9.489383e-01_r8,9.515147e-01_r8,9.538391e-01_r8,9.559510e-01_r8, &
       9.578816e-01_r8,9.596561e-01_r8,9.612949e-01_r8,9.628148e-01_r8,9.642301e-01_r8, &
       9.655524e-01_r8,9.667918e-01_r8,9.679568e-01_r8,9.690549e-01_r8,9.700923e-01_r8, &
       9.710747e-01_r8,9.720070e-01_r8,9.728933e-01_r8,9.737376e-01_r8,9.745431e-01_r8, &
       9.753129e-01_r8,9.760497e-01_r8,9.767558e-01_r8,9.774336e-01_r8,9.780849e-01_r8, &
       9.787115e-01_r8,9.793151e-01_r8,9.798972e-01_r8,9.804592e-01_r8,9.810022e-01_r8, &
       9.815276e-01_r8,9.820363e-01_r8,9.825294e-01_r8/)
      ASYICE2(:,9) = (/ &
!    BAND 9
       8.447005e-01_r8,8.780132e-01_r8,8.968896e-01_r8,9.095104e-01_r8,9.187442e-01_r8, &
       9.258958e-01_r8,9.316569e-01_r8,9.364336e-01_r8,9.404822e-01_r8,9.439738e-01_r8, &
       9.470275e-01_r8,9.497295e-01_r8,9.521436e-01_r8,9.543184e-01_r8,9.562917e-01_r8, &
       9.580932e-01_r8,9.597468e-01_r8,9.612720e-01_r8,9.626849e-01_r8,9.639986e-01_r8, &
       9.652243e-01_r8,9.663716e-01_r8,9.674484e-01_r8,9.684618e-01_r8,9.694176e-01_r8, &
       9.703212e-01_r8,9.711770e-01_r8,9.719892e-01_r8,9.727611e-01_r8,9.734961e-01_r8, &
       9.741968e-01_r8,9.748658e-01_r8,9.755053e-01_r8,9.761174e-01_r8,9.767038e-01_r8, &
       9.772662e-01_r8,9.778062e-01_r8,9.783250e-01_r8,9.788239e-01_r8,9.793041e-01_r8, &
       9.797666e-01_r8,9.802124e-01_r8,9.806423e-01_r8/)
      ASYICE2(:,10) = (/ &
!    BAND 10
       8.564947e-01_r8,8.915851e-01_r8,9.102676e-01_r8,9.221993e-01_r8,9.306152e-01_r8, &
       9.369369e-01_r8,9.418971e-01_r8,9.459156e-01_r8,9.492519e-01_r8,9.520761e-01_r8, &
       9.545046e-01_r8,9.566202e-01_r8,9.584835e-01_r8,9.601400e-01_r8,9.616244e-01_r8, &
       9.629641e-01_r8,9.641805e-01_r8,9.652912e-01_r8,9.663102e-01_r8,9.672493e-01_r8, &
       9.681181e-01_r8,9.689247e-01_r8,9.696762e-01_r8,9.703783e-01_r8,9.710361e-01_r8, &
       9.716540e-01_r8,9.722357e-01_r8,9.727846e-01_r8,9.733035e-01_r8,9.737950e-01_r8, &
       9.742615e-01_r8,9.747048e-01_r8,9.751268e-01_r8,9.755291e-01_r8,9.759132e-01_r8, &
       9.762803e-01_r8,9.766317e-01_r8,9.769684e-01_r8,9.772913e-01_r8,9.776014e-01_r8, &
       9.778994e-01_r8,9.781862e-01_r8,9.784623e-01_r8/)
      ASYICE2(:,11) = (/ &
!    BAND 11
       8.697900e-01_r8,9.013810e-01_r8,9.181518e-01_r8,9.288524e-01_r8,9.363985e-01_r8, &
       9.420679e-01_r8,9.465182e-01_r8,9.501255e-01_r8,9.531224e-01_r8,9.556611e-01_r8, &
       9.578459e-01_r8,9.597507e-01_r8,9.614299e-01_r8,9.629240e-01_r8,9.642643e-01_r8, &
       9.654751e-01_r8,9.665758e-01_r8,9.675818e-01_r8,9.685059e-01_r8,9.693585e-01_r8, &
       9.701484e-01_r8,9.708827e-01_r8,9.715676e-01_r8,9.722085e-01_r8,9.728097e-01_r8, &
       9.733753e-01_r8,9.739087e-01_r8,9.744127e-01_r8,9.748899e-01_r8,9.753428e-01_r8, &
       9.757732e-01_r8,9.761830e-01_r8,9.765738e-01_r8,9.769470e-01_r8,9.773040e-01_r8, &
       9.776459e-01_r8,9.779737e-01_r8,9.782883e-01_r8,9.785908e-01_r8,9.788818e-01_r8, &
       9.791621e-01_r8,9.794323e-01_r8,9.796931e-01_r8/)
      ASYICE2(:,12) = (/ &
!    BAND 12
       8.283310e-01_r8,8.543591e-01_r8,8.712236e-01_r8,8.835970e-01_r8,8.933168e-01_r8, &
       9.012909e-01_r8,9.080327e-01_r8,9.138601e-01_r8,9.189835e-01_r8,9.235486e-01_r8, &
       9.276609e-01_r8,9.313988e-01_r8,9.348222e-01_r8,9.379780e-01_r8,9.409034e-01_r8, &
       9.436284e-01_r8,9.461776e-01_r8,9.485715e-01_r8,9.508272e-01_r8,9.529591e-01_r8, &
       9.549796e-01_r8,9.568992e-01_r8,9.587273e-01_r8,9.604716e-01_r8,9.621394e-01_r8, &
       9.637367e-01_r8,9.652691e-01_r8,9.667413e-01_r8,9.681578e-01_r8,9.695225e-01_r8, &
       9.708388e-01_r8,9.721100e-01_r8,9.733388e-01_r8,9.745280e-01_r8,9.756799e-01_r8, &
       9.767967e-01_r8,9.778803e-01_r8,9.789326e-01_r8,9.799553e-01_r8,9.809499e-01_r8, &
       9.819179e-01_r8,9.828606e-01_r8,9.837793e-01_r8/)
      ASYICE2(:,13) = (/ &
!    BAND 13
       8.306222e-01_r8,8.696458e-01_r8,8.911483e-01_r8,9.052375e-01_r8,9.153823e-01_r8, &
       9.231361e-01_r8,9.293121e-01_r8,9.343824e-01_r8,9.386425e-01_r8,9.422880e-01_r8, &
       9.454540e-01_r8,9.482376e-01_r8,9.507103e-01_r8,9.529263e-01_r8,9.549272e-01_r8, &
       9.567459e-01_r8,9.584087e-01_r8,9.599368e-01_r8,9.613475e-01_r8,9.626553e-01_r8, &
       9.638721e-01_r8,9.650082e-01_r8,9.660721e-01_r8,9.670713e-01_r8,9.680121e-01_r8, &
       9.689001e-01_r8,9.697400e-01_r8,9.705361e-01_r8,9.712922e-01_r8,9.720114e-01_r8, &
       9.726968e-01_r8,9.733509e-01_r8,9.739760e-01_r8,9.745744e-01_r8,9.751477e-01_r8, &
       9.756979e-01_r8,9.762263e-01_r8,9.767345e-01_r8,9.772236e-01_r8,9.776949e-01_r8, &
       9.781495e-01_r8,9.785882e-01_r8,9.790121e-01_r8/)
      ASYICE2(:,14) = (/ &
!    BAND 14
       8.240566e-01_r8,8.644102e-01_r8,8.868070e-01_r8,9.015382e-01_r8,9.121690e-01_r8, &
       9.203056e-01_r8,9.267921e-01_r8,9.321201e-01_r8,9.365980e-01_r8,9.404302e-01_r8, &
       9.437583e-01_r8,9.466840e-01_r8,9.492823e-01_r8,9.516101e-01_r8,9.537112e-01_r8, &
       9.556203e-01_r8,9.573649e-01_r8,9.589674e-01_r8,9.604460e-01_r8,9.618159e-01_r8, &
       9.630899e-01_r8,9.642786e-01_r8,9.653911e-01_r8,9.664352e-01_r8,9.674178e-01_r8, &
       9.683445e-01_r8,9.692205e-01_r8,9.700502e-01_r8,9.708376e-01_r8,9.715862e-01_r8, &
       9.722990e-01_r8,9.729789e-01_r8,9.736282e-01_r8,9.742491e-01_r8,9.748438e-01_r8, &
       9.754140e-01_r8,9.759613e-01_r8,9.764872e-01_r8,9.769931e-01_r8,9.774802e-01_r8, &
       9.779496e-01_r8,9.784025e-01_r8,9.788396e-01_r8/)
      ASYICE2(:,15) = (/ &
!    BAND 15
       7.868651e-01_r8,8.322305e-01_r8,8.577581e-01_r8,8.747072e-01_r8,8.870264e-01_r8, &
       8.965099e-01_r8,9.041069e-01_r8,9.103734e-01_r8,9.156595e-01_r8,9.201985e-01_r8, &
       9.241524e-01_r8,9.276378e-01_r8,9.307412e-01_r8,9.335281e-01_r8,9.360494e-01_r8, &
       9.383451e-01_r8,9.404472e-01_r8,9.423817e-01_r8,9.441700e-01_r8,9.458298e-01_r8, &
       9.473759e-01_r8,9.488208e-01_r8,9.501753e-01_r8,9.514485e-01_r8,9.526482e-01_r8, &
       9.537815e-01_r8,9.548541e-01_r8,9.558715e-01_r8,9.568383e-01_r8,9.577585e-01_r8, &
       9.586359e-01_r8,9.594736e-01_r8,9.602747e-01_r8,9.610417e-01_r8,9.617771e-01_r8, &
       9.624829e-01_r8,9.631611e-01_r8,9.638136e-01_r8,9.644418e-01_r8,9.650473e-01_r8, &
       9.656314e-01_r8,9.661955e-01_r8,9.667405e-01_r8/)
      ASYICE2(:,16) = (/ &
!    BAND 16
       7.946655e-01_r8,8.547685e-01_r8,8.806016e-01_r8,8.949880e-01_r8,9.041676e-01_r8, &
       9.105399e-01_r8,9.152249e-01_r8,9.188160e-01_r8,9.216573e-01_r8,9.239620e-01_r8, &
       9.258695e-01_r8,9.274745e-01_r8,9.288441e-01_r8,9.300267e-01_r8,9.310584e-01_r8, & 
       9.319665e-01_r8,9.327721e-01_r8,9.334918e-01_r8,9.341387e-01_r8,9.347236e-01_r8, &
       9.352551e-01_r8,9.357402e-01_r8,9.361850e-01_r8,9.365942e-01_r8,9.369722e-01_r8, &
       9.373225e-01_r8,9.376481e-01_r8,9.379516e-01_r8,9.382352e-01_r8,9.385010e-01_r8, &
       9.387505e-01_r8,9.389854e-01_r8,9.392070e-01_r8,9.394163e-01_r8,9.396145e-01_r8, &
       9.398024e-01_r8,9.399809e-01_r8,9.401508e-01_r8,9.403126e-01_r8,9.404670e-01_r8, &
       9.406144e-01_r8,9.407555e-01_r8,9.408906e-01_r8/)

!      SINGLE-SCATTERING ALBEDO: Unitless
      SSAICE3(:,1) = (/ &
!    BAND 1
       8.733202e-02_r8,1.733042e-01_r8,2.494597e-01_r8,2.853637e-01_r8,3.129915e-01_r8, &
       3.366261e-01_r8,3.574925e-01_r8,3.761110e-01_r8,3.927678e-01_r8,4.076523e-01_r8, &
       4.209083e-01_r8,4.326557e-01_r8,4.430008e-01_r8,4.520419e-01_r8,4.598721e-01_r8, &
       4.665813e-01_r8,4.722567e-01_r8,4.769845e-01_r8,4.808494e-01_r8,4.839354e-01_r8, &
       4.863257e-01_r8,4.881034e-01_r8,4.893508e-01_r8,4.901503e-01_r8,4.905838e-01_r8, &
       4.907333e-01_r8,4.906804e-01_r8,4.905066e-01_r8,4.902935e-01_r8,4.901225e-01_r8, &
       4.900749e-01_r8,4.902321e-01_r8,4.906751e-01_r8,4.914853e-01_r8,4.927439e-01_r8, &
       4.945318e-01_r8,4.969304e-01_r8,5.000205e-01_r8,5.038834e-01_r8,5.086000e-01_r8, &
       5.142513e-01_r8,5.209185e-01_r8,5.286824e-01_r8,5.376241e-01_r8,5.478246e-01_r8, &
       5.593649e-01_r8/)
      SSAICE3(:,2) = (/ &
!    BAND 2
       7.617260e-01_r8,8.272990e-01_r8,8.134738e-01_r8,8.040737e-01_r8,7.953976e-01_r8, &
       7.869914e-01_r8,7.787330e-01_r8,7.705796e-01_r8,7.625153e-01_r8,7.545359e-01_r8, &
       7.466425e-01_r8,7.388394e-01_r8,7.311327e-01_r8,7.235295e-01_r8,7.160381e-01_r8, &
       7.086672e-01_r8,7.014261e-01_r8,6.943246e-01_r8,6.873731e-01_r8,6.805821e-01_r8, &
       6.739626e-01_r8,6.675263e-01_r8,6.612848e-01_r8,6.552505e-01_r8,6.494361e-01_r8, &
       6.438547e-01_r8,6.385199e-01_r8,6.334458e-01_r8,6.286469e-01_r8,6.241384e-01_r8, &
       6.199357e-01_r8,6.160553e-01_r8,6.125138e-01_r8,6.093286e-01_r8,6.065180e-01_r8, &
       6.041006e-01_r8,6.020961e-01_r8,6.005248e-01_r8,5.994078e-01_r8,5.987672e-01_r8, &
       5.986259e-01_r8,5.990079e-01_r8,5.999382e-01_r8,6.014428e-01_r8,6.035491e-01_r8, &
       6.062855e-01_r8/)
      SSAICE3(:,3) = (/ &
!    BAND 3
       6.485070e-01_r8,6.545008e-01_r8,6.483874e-01_r8,6.411445e-01_r8,6.338616e-01_r8, &
       6.268166e-01_r8,6.201055e-01_r8,6.137647e-01_r8,6.078060e-01_r8,6.022297e-01_r8, &
       5.970301e-01_r8,5.921982e-01_r8,5.877229e-01_r8,5.835920e-01_r8,5.797922e-01_r8, &
       5.763099e-01_r8,5.731308e-01_r8,5.702407e-01_r8,5.676249e-01_r8,5.652685e-01_r8, &
       5.631568e-01_r8,5.612747e-01_r8,5.596072e-01_r8,5.581391e-01_r8,5.568553e-01_r8, &
       5.557405e-01_r8,5.547796e-01_r8,5.539571e-01_r8,5.532579e-01_r8,5.526665e-01_r8, &
       5.521676e-01_r8,5.517458e-01_r8,5.513856e-01_r8,5.510717e-01_r8,5.507886e-01_r8, &
       5.505208e-01_r8,5.502528e-01_r8,5.499691e-01_r8,5.496542e-01_r8,5.492923e-01_r8, &
       5.488681e-01_r8,5.483658e-01_r8,5.477697e-01_r8,5.470643e-01_r8,5.462337e-01_r8, &
       5.452623e-01_r8/)
      SSAICE3(:,4) = (/ &
!    BAND 4
       5.367012e-01_r8,5.372179e-01_r8,5.366731e-01_r8,5.359685e-01_r8,5.352768e-01_r8, &
       5.346515e-01_r8,5.341121e-01_r8,5.336653e-01_r8,5.333122e-01_r8,5.330513e-01_r8, &
       5.328793e-01_r8,5.327924e-01_r8,5.327859e-01_r8,5.328549e-01_r8,5.329943e-01_r8, &
       5.331989e-01_r8,5.334632e-01_r8,5.337819e-01_r8,5.341494e-01_r8,5.345602e-01_r8, &
       5.350088e-01_r8,5.354897e-01_r8,5.359972e-01_r8,5.365259e-01_r8,5.370701e-01_r8, &
       5.376244e-01_r8,5.381831e-01_r8,5.387408e-01_r8,5.392920e-01_r8,5.398311e-01_r8, &
       5.403526e-01_r8,5.408512e-01_r8,5.413213e-01_r8,5.417575e-01_r8,5.421545e-01_r8, &
       5.425067e-01_r8,5.428090e-01_r8,5.430558e-01_r8,5.432420e-01_r8,5.433623e-01_r8, &
       5.434113e-01_r8,5.433838e-01_r8,5.432748e-01_r8,5.430789e-01_r8,5.427911e-01_r8, &
       5.424063e-01_r8/)
      SSAICE3(:,5) = (/ &
!    BAND 5
       4.144473e-01_r8,4.233085e-01_r8,4.317951e-01_r8,4.397881e-01_r8,4.472768e-01_r8, &
       4.542704e-01_r8,4.607843e-01_r8,4.668359e-01_r8,4.724439e-01_r8,4.776272e-01_r8, &
       4.824055e-01_r8,4.867983e-01_r8,4.908252e-01_r8,4.945063e-01_r8,4.978612e-01_r8, &
       5.009100e-01_r8,5.036726e-01_r8,5.061689e-01_r8,5.084192e-01_r8,5.104433e-01_r8, &
       5.122613e-01_r8,5.138934e-01_r8,5.153598e-01_r8,5.166804e-01_r8,5.178756e-01_r8, &
       5.189655e-01_r8,5.199702e-01_r8,5.209102e-01_r8,5.218055e-01_r8,5.226765e-01_r8, &
       5.235435e-01_r8,5.244268e-01_r8,5.253468e-01_r8,5.263238e-01_r8,5.273783e-01_r8, &
       5.285306e-01_r8,5.298012e-01_r8,5.312106e-01_r8,5.327792e-01_r8,5.345277e-01_r8, &
       5.364765e-01_r8,5.386462e-01_r8,5.410574e-01_r8,5.437308e-01_r8,5.466871e-01_r8, &
       5.499469e-01_r8/)
      SSAICE3(:,6) = (/ &
!    BAND 6
       3.685509e-01_r8,4.062125e-01_r8,4.219575e-01_r8,4.336348e-01_r8,4.434522e-01_r8, &
       4.520497e-01_r8,4.596957e-01_r8,4.665324e-01_r8,4.726498e-01_r8,4.781130e-01_r8, &
       4.829739e-01_r8,4.872771e-01_r8,4.910627e-01_r8,4.943681e-01_r8,4.972286e-01_r8, &
       4.996785e-01_r8,5.017514e-01_r8,5.034799e-01_r8,5.048966e-01_r8,5.060334e-01_r8, &
       5.069224e-01_r8,5.075951e-01_r8,5.080833e-01_r8,5.084185e-01_r8,5.086320e-01_r8, &
       5.087554e-01_r8,5.088200e-01_r8,5.088572e-01_r8,5.088985e-01_r8,5.089752e-01_r8, &
       5.091187e-01_r8,5.093605e-01_r8,5.097319e-01_r8,5.102646e-01_r8,5.109899e-01_r8, &
       5.119394e-01_r8,5.131447e-01_r8,5.146373e-01_r8,5.164489e-01_r8,5.186112e-01_r8, &
       5.211558e-01_r8,5.241145e-01_r8,5.275191e-01_r8,5.314014e-01_r8,5.357933e-01_r8, &
       5.407268e-01_r8/)
      SSAICE3(:,7) = (/ &
!    BAND 7
       7.347648e-01_r8,6.945725e-01_r8,6.844409e-01_r8,6.757818e-01_r8,6.676970e-01_r8, &
       6.599913e-01_r8,6.525964e-01_r8,6.454817e-01_r8,6.386312e-01_r8,6.320353e-01_r8, &
       6.256876e-01_r8,6.195836e-01_r8,6.137200e-01_r8,6.080940e-01_r8,6.027033e-01_r8, &
       5.975460e-01_r8,5.926203e-01_r8,5.879246e-01_r8,5.834574e-01_r8,5.792174e-01_r8, &
       5.752032e-01_r8,5.714137e-01_r8,5.678476e-01_r8,5.645040e-01_r8,5.613817e-01_r8, &
       5.584797e-01_r8,5.557971e-01_r8,5.533330e-01_r8,5.510865e-01_r8,5.490569e-01_r8, &
       5.472434e-01_r8,5.456452e-01_r8,5.442618e-01_r8,5.430926e-01_r8,5.421369e-01_r8, &
       5.413943e-01_r8,5.408645e-01_r8,5.405470e-01_r8,5.404415e-01_r8,5.405479e-01_r8, &
       5.408659e-01_r8,5.413956e-01_r8,5.421369e-01_r8,5.430900e-01_r8,5.442550e-01_r8, &
       5.456321e-01_r8/)
      SSAICE3(:,8) = (/ &
!    BAND 8
       7.565911e-01_r8,7.410307e-01_r8,7.267244e-01_r8,7.132696e-01_r8,7.005889e-01_r8, &
       6.886431e-01_r8,6.774024e-01_r8,6.668406e-01_r8,6.569327e-01_r8,6.476543e-01_r8, &
       6.389815e-01_r8,6.308907e-01_r8,6.233581e-01_r8,6.163602e-01_r8,6.098736e-01_r8, &
       6.038748e-01_r8,5.983403e-01_r8,5.932468e-01_r8,5.885708e-01_r8,5.842888e-01_r8, &
       5.803776e-01_r8,5.768135e-01_r8,5.735734e-01_r8,5.706336e-01_r8,5.679708e-01_r8, &
       5.655616e-01_r8,5.633825e-01_r8,5.614101e-01_r8,5.596210e-01_r8,5.579915e-01_r8, &
       5.564984e-01_r8,5.551181e-01_r8,5.538271e-01_r8,5.526020e-01_r8,5.514191e-01_r8, &
       5.502550e-01_r8,5.490862e-01_r8,5.478890e-01_r8,5.466400e-01_r8,5.453154e-01_r8, &
       5.438918e-01_r8,5.423455e-01_r8,5.406528e-01_r8,5.387902e-01_r8,5.367338e-01_r8, &
       5.344600e-01_r8/)
      SSAICE3(:,9) = (/ &
!    BAND 9
       7.253131e-01_r8,7.112494e-01_r8,6.974209e-01_r8,6.843066e-01_r8,6.719929e-01_r8, &
       6.604896e-01_r8,6.497830e-01_r8,6.398502e-01_r8,6.306641e-01_r8,6.221956e-01_r8, &
       6.144141e-01_r8,6.072887e-01_r8,6.007877e-01_r8,5.948792e-01_r8,5.895312e-01_r8, &
       5.847116e-01_r8,5.803881e-01_r8,5.765282e-01_r8,5.730996e-01_r8,5.700698e-01_r8, &
       5.674064e-01_r8,5.650768e-01_r8,5.630484e-01_r8,5.612888e-01_r8,5.597653e-01_r8, &
       5.584452e-01_r8,5.572961e-01_r8,5.562853e-01_r8,5.553800e-01_r8,5.545478e-01_r8, &
       5.537558e-01_r8,5.529714e-01_r8,5.521620e-01_r8,5.512948e-01_r8,5.503371e-01_r8, &
       5.492561e-01_r8,5.480192e-01_r8,5.465935e-01_r8,5.449463e-01_r8,5.430449e-01_r8, &
       5.408564e-01_r8,5.383480e-01_r8,5.354869e-01_r8,5.322403e-01_r8,5.285753e-01_r8, &
       5.244591e-01_r8/)
      SSAICE3(:,10) = (/ &
!    BAND 10
       6.692312e-01_r8,6.569887e-01_r8,6.455728e-01_r8,6.349744e-01_r8,6.251679e-01_r8, &
       6.161241e-01_r8,6.078124e-01_r8,6.002017e-01_r8,5.932608e-01_r8,5.869582e-01_r8, &
       5.812625e-01_r8,5.761422e-01_r8,5.715658e-01_r8,5.675017e-01_r8,5.639183e-01_r8, &
       5.607841e-01_r8,5.580676e-01_r8,5.557370e-01_r8,5.537609e-01_r8,5.521078e-01_r8, &
       5.507459e-01_r8,5.496437e-01_r8,5.487696e-01_r8,5.480921e-01_r8,5.475796e-01_r8, &
       5.472004e-01_r8,5.469230e-01_r8,5.467159e-01_r8,5.465474e-01_r8,5.463859e-01_r8, &
       5.461999e-01_r8,5.459577e-01_r8,5.456279e-01_r8,5.451788e-01_r8,5.445788e-01_r8, &
       5.437964e-01_r8,5.427999e-01_r8,5.415578e-01_r8,5.400386e-01_r8,5.382105e-01_r8, &
       5.360422e-01_r8,5.335019e-01_r8,5.305581e-01_r8,5.271792e-01_r8,5.233336e-01_r8, &
       5.189898e-01_r8/)
      SSAICE3(:,11) = (/ &
!    BAND 11
       6.597440e-01_r8,6.486312e-01_r8,6.385754e-01_r8,6.293389e-01_r8,6.208334e-01_r8, &
       6.130062e-01_r8,6.058168e-01_r8,5.992306e-01_r8,5.932156e-01_r8,5.877415e-01_r8, &
       5.827791e-01_r8,5.782999e-01_r8,5.742760e-01_r8,5.706800e-01_r8,5.674845e-01_r8, &
       5.646626e-01_r8,5.621876e-01_r8,5.600328e-01_r8,5.581718e-01_r8,5.565782e-01_r8, &
       5.552258e-01_r8,5.540883e-01_r8,5.531398e-01_r8,5.523542e-01_r8,5.517054e-01_r8, &
       5.511675e-01_r8,5.507147e-01_r8,5.503210e-01_r8,5.499606e-01_r8,5.496077e-01_r8, &
       5.492364e-01_r8,5.488211e-01_r8,5.483359e-01_r8,5.477550e-01_r8,5.470529e-01_r8, &
       5.462036e-01_r8,5.451816e-01_r8,5.439610e-01_r8,5.425163e-01_r8,5.408217e-01_r8, &
       5.388515e-01_r8,5.365800e-01_r8,5.339817e-01_r8,5.310307e-01_r8,5.277016e-01_r8, &
       5.239685e-01_r8/)
      SSAICE3(:,12) = (/ &
!    BAND 12
       8.415691e-01_r8,8.237634e-01_r8,8.070239e-01_r8,7.912304e-01_r8,7.763216e-01_r8, &
       7.622533e-01_r8,7.489874e-01_r8,7.364895e-01_r8,7.247271e-01_r8,7.136691e-01_r8, & 
       7.032852e-01_r8,6.935461e-01_r8,6.844229e-01_r8,6.758871e-01_r8,6.679108e-01_r8, &
       6.604662e-01_r8,6.535258e-01_r8,6.470623e-01_r8,6.410487e-01_r8,6.354581e-01_r8, &
       6.302637e-01_r8,6.254387e-01_r8,6.209567e-01_r8,6.167911e-01_r8,6.129154e-01_r8, &
       6.093033e-01_r8,6.059286e-01_r8,6.027648e-01_r8,5.997858e-01_r8,5.969653e-01_r8, &
       5.942772e-01_r8,5.916953e-01_r8,5.891935e-01_r8,5.867455e-01_r8,5.843253e-01_r8, &
       5.819067e-01_r8,5.794637e-01_r8,5.769700e-01_r8,5.743997e-01_r8,5.717265e-01_r8, &
       5.689243e-01_r8,5.659670e-01_r8,5.628284e-01_r8,5.594825e-01_r8,5.559030e-01_r8, &
       5.520638e-01_r8/)
      SSAICE3(:,13) = (/ &
!    BAND 13
       7.418384e-01_r8,7.252704e-01_r8,7.097171e-01_r8,6.951785e-01_r8,6.816320e-01_r8, &
       6.690487e-01_r8,6.573967e-01_r8,6.466425e-01_r8,6.367511e-01_r8,6.276869e-01_r8, & 
       6.194133e-01_r8,6.118931e-01_r8,6.050886e-01_r8,5.989616e-01_r8,5.934736e-01_r8, &
       5.885856e-01_r8,5.842585e-01_r8,5.804528e-01_r8,5.771289e-01_r8,5.742470e-01_r8, &
       5.717672e-01_r8,5.696493e-01_r8,5.678531e-01_r8,5.663385e-01_r8,5.650650e-01_r8, &
       5.639922e-01_r8,5.630797e-01_r8,5.622870e-01_r8,5.615735e-01_r8,5.608986e-01_r8, &
       5.602218e-01_r8,5.595025e-01_r8,5.586999e-01_r8,5.577736e-01_r8,5.566826e-01_r8, &
       5.553865e-01_r8,5.538444e-01_r8,5.520156e-01_r8,5.498592e-01_r8,5.473345e-01_r8, &
       5.444005e-01_r8,5.410162e-01_r8,5.371407e-01_r8,5.327328e-01_r8,5.277513e-01_r8, &
       5.221549e-01_r8/)
      SSAICE3(:,14) = (/ &
!    BAND 14
       7.431625e-01_r8,7.253592e-01_r8,7.089350e-01_r8,6.937690e-01_r8,6.797727e-01_r8, &
       6.668739e-01_r8,6.550099e-01_r8,6.441244e-01_r8,6.341649e-01_r8,6.250822e-01_r8, &
       6.168287e-01_r8,6.093588e-01_r8,6.026277e-01_r8,5.965913e-01_r8,5.912065e-01_r8, &
       5.864305e-01_r8,5.822207e-01_r8,5.785351e-01_r8,5.753318e-01_r8,5.725690e-01_r8, &
       5.702053e-01_r8,5.681992e-01_r8,5.665094e-01_r8,5.650948e-01_r8,5.639141e-01_r8, &
       5.629263e-01_r8,5.620905e-01_r8,5.613656e-01_r8,5.607109e-01_r8,5.600853e-01_r8, &
       5.594482e-01_r8,5.587586e-01_r8,5.579758e-01_r8,5.570591e-01_r8,5.559676e-01_r8, &
       5.546606e-01_r8,5.530972e-01_r8,5.512367e-01_r8,5.490383e-01_r8,5.464611e-01_r8, &
       5.434641e-01_r8,5.400063e-01_r8,5.360468e-01_r8,5.315444e-01_r8,5.264578e-01_r8, &
       5.207459e-01_r8/)
      SSAICE3(:,15) = (/ &
!    BAND 15
       8.214523e-01_r8,8.021560e-01_r8,7.840431e-01_r8,7.670541e-01_r8,7.511382e-01_r8, &
       7.362494e-01_r8,7.223443e-01_r8,7.093816e-01_r8,6.973208e-01_r8,6.861227e-01_r8, &
       6.757482e-01_r8,6.661588e-01_r8,6.573164e-01_r8,6.491829e-01_r8,6.417204e-01_r8, &
       6.348911e-01_r8,6.286574e-01_r8,6.229816e-01_r8,6.178260e-01_r8,6.131531e-01_r8, &
       6.089254e-01_r8,6.051053e-01_r8,6.016552e-01_r8,5.985377e-01_r8,5.957152e-01_r8, &
       5.931503e-01_r8,5.908054e-01_r8,5.886430e-01_r8,5.866256e-01_r8,5.847156e-01_r8, &
       5.828757e-01_r8,5.810683e-01_r8,5.792558e-01_r8,5.774008e-01_r8,5.754657e-01_r8, &
       5.734130e-01_r8,5.712052e-01_r8,5.688048e-01_r8,5.661741e-01_r8,5.632757e-01_r8, &
       5.600721e-01_r8,5.565255e-01_r8,5.525985e-01_r8,5.482534e-01_r8,5.434526e-01_r8, &
       5.381586e-01_r8/)
      SSAICE3(:,16) = (/ &
!    BAND 16
       6.749442e-01_r8,6.649947e-01_r8,6.565828e-01_r8,6.489928e-01_r8,6.420046e-01_r8, &
       6.355231e-01_r8,6.294964e-01_r8,6.238901e-01_r8,6.186783e-01_r8,6.138395e-01_r8, &
       6.093543e-01_r8,6.052049e-01_r8,6.013742e-01_r8,5.978457e-01_r8,5.946030e-01_r8, &
       5.916302e-01_r8,5.889115e-01_r8,5.864310e-01_r8,5.841731e-01_r8,5.821221e-01_r8, &
       5.802624e-01_r8,5.785785e-01_r8,5.770549e-01_r8,5.756759e-01_r8,5.744262e-01_r8, &
       5.732901e-01_r8,5.722524e-01_r8,5.712974e-01_r8,5.704097e-01_r8,5.695739e-01_r8, &
       5.687747e-01_r8,5.679964e-01_r8,5.672238e-01_r8,5.664415e-01_r8,5.656340e-01_r8, &
       5.647860e-01_r8,5.638821e-01_r8,5.629070e-01_r8,5.618452e-01_r8,5.606815e-01_r8, &
       5.594006e-01_r8,5.579870e-01_r8,5.564255e-01_r8,5.547008e-01_r8,5.527976e-01_r8, &
       5.507005e-01_r8/)
!     ASYMMETRY FACTOR: Unitless
      ASYICE3(:,1) = (/ &
!    BAND 1
       6.596879e-01_r8,6.405129e-01_r8,6.397929e-01_r8,6.577098e-01_r8,6.788112e-01_r8, &
       6.996961e-01_r8,7.195365e-01_r8,7.380555e-01_r8,7.551580e-01_r8,7.708248e-01_r8, &
       7.850743e-01_r8,7.979456e-01_r8,8.094902e-01_r8,8.197678e-01_r8,8.288438e-01_r8, &
       8.367874e-01_r8,8.436711e-01_r8,8.495695e-01_r8,8.545592e-01_r8,8.587185e-01_r8, &
       8.621270e-01_r8,8.648654e-01_r8,8.670153e-01_r8,8.686596e-01_r8,8.698817e-01_r8, &
       8.707659e-01_r8,8.713971e-01_r8,8.718612e-01_r8,8.722444e-01_r8,8.726337e-01_r8, &
       8.731167e-01_r8,8.737814e-01_r8,8.747166e-01_r8,8.760115e-01_r8,8.777560e-01_r8, &
       8.800403e-01_r8,8.829553e-01_r8,8.865923e-01_r8,8.910430e-01_r8,8.963998e-01_r8, &
       9.027552e-01_r8,9.102023e-01_r8,9.188345e-01_r8,9.287454e-01_r8,9.400292e-01_r8, &
       9.523461e-01_r8/)
      ASYICE3(:,2) = (/ &
!    BAND 2
       7.389769e-01_r8,7.411960e-01_r8,7.490019e-01_r8,7.564984e-01_r8,7.637359e-01_r8, &
       7.707210e-01_r8,7.774570e-01_r8,7.839465e-01_r8,7.901930e-01_r8,7.962001e-01_r8, &
       8.019724e-01_r8,8.075150e-01_r8,8.128335e-01_r8,8.179342e-01_r8,8.228237e-01_r8, &
       8.275093e-01_r8,8.319986e-01_r8,8.362997e-01_r8,8.404208e-01_r8,8.443707e-01_r8, &
       8.481583e-01_r8,8.517929e-01_r8,8.552840e-01_r8,8.586412e-01_r8,8.618744e-01_r8, &
       8.649936e-01_r8,8.680092e-01_r8,8.709315e-01_r8,8.737713e-01_r8,8.765394e-01_r8, &
       8.792470e-01_r8,8.819057e-01_r8,8.845275e-01_r8,8.871247e-01_r8,8.897104e-01_r8, &
       8.922982e-01_r8,8.949025e-01_r8,8.975389e-01_r8,9.002236e-01_r8,9.029745e-01_r8, &
       9.058105e-01_r8,9.087524e-01_r8,9.118223e-01_r8,9.150447e-01_r8,9.184456e-01_r8, &
       9.220537e-01_r8/)
      ASYICE3(:,3) = (/ &
!    BAND 3
       7.533538e-01_r8,7.631828e-01_r8,7.742293e-01_r8,7.848919e-01_r8,7.950315e-01_r8, &
       8.046253e-01_r8,8.136762e-01_r8,8.221958e-01_r8,8.301991e-01_r8,8.377028e-01_r8, &
       8.447247e-01_r8,8.512829e-01_r8,8.573959e-01_r8,8.630826e-01_r8,8.683620e-01_r8, &
       8.732532e-01_r8,8.777756e-01_r8,8.819487e-01_r8,8.857920e-01_r8,8.893253e-01_r8, &
       8.925684e-01_r8,8.955411e-01_r8,8.982635e-01_r8,9.007556e-01_r8,9.030377e-01_r8, &
       9.051300e-01_r8,9.070528e-01_r8,9.088266e-01_r8,9.104718e-01_r8,9.120090e-01_r8, &
       9.134588e-01_r8,9.148418e-01_r8,9.161788e-01_r8,9.174904e-01_r8,9.187976e-01_r8, &
       9.201211e-01_r8,9.214818e-01_r8,9.229007e-01_r8,9.243986e-01_r8,9.259967e-01_r8, &
       9.277160e-01_r8,9.295774e-01_r8,9.316023e-01_r8,9.338116e-01_r8,9.362266e-01_r8, &
       9.388685e-01_r8/)
      ASYICE3(:,4) = (/ &
!    BAND 4
       7.840105e-01_r8,7.959983e-01_r8,8.074838e-01_r8,8.182465e-01_r8,8.282675e-01_r8, &
       8.375608e-01_r8,8.461497e-01_r8,8.540609e-01_r8,8.613228e-01_r8,8.679643e-01_r8, &
       8.740148e-01_r8,8.795041e-01_r8,8.844620e-01_r8,8.889182e-01_r8,8.929028e-01_r8, &
       8.964456e-01_r8,8.995766e-01_r8,9.023257e-01_r8,9.047230e-01_r8,9.067984e-01_r8, &
       9.085817e-01_r8,9.101031e-01_r8,9.113922e-01_r8,9.124791e-01_r8,9.133937e-01_r8, & 
       9.141659e-01_r8,9.148254e-01_r8,9.154023e-01_r8,9.159263e-01_r8,9.164274e-01_r8, &
       9.169353e-01_r8,9.174799e-01_r8,9.180910e-01_r8,9.187986e-01_r8,9.196324e-01_r8, &
       9.206223e-01_r8,9.217981e-01_r8,9.231897e-01_r8,9.248270e-01_r8,9.267399e-01_r8, &
       9.289582e-01_r8,9.315119e-01_r8,9.344308e-01_r8,9.377451e-01_r8,9.414845e-01_r8, &
       9.456792e-01_r8/)
      ASYICE3(:,5) = (/ &
!    BAND 5
       8.304283e-01_r8,8.399040e-01_r8,8.485749e-01_r8,8.565538e-01_r8,8.638881e-01_r8, &
       8.706126e-01_r8,8.767578e-01_r8,8.823528e-01_r8,8.874253e-01_r8,8.920026e-01_r8, &
       8.961113e-01_r8,8.997780e-01_r8,9.030287e-01_r8,9.058894e-01_r8,9.083861e-01_r8, &
       9.105442e-01_r8,9.123895e-01_r8,9.139473e-01_r8,9.152430e-01_r8,9.163020e-01_r8, &
       9.171497e-01_r8,9.178111e-01_r8,9.183116e-01_r8,9.186762e-01_r8,9.189303e-01_r8, &
       9.190989e-01_r8,9.192071e-01_r8,9.192802e-01_r8,9.193432e-01_r8,9.194211e-01_r8, &
       9.195392e-01_r8,9.197225e-01_r8,9.199961e-01_r8,9.203850e-01_r8,9.209143e-01_r8, &
       9.216090e-01_r8,9.224942e-01_r8,9.235947e-01_r8,9.249355e-01_r8,9.265416e-01_r8, &
       9.284376e-01_r8,9.306484e-01_r8,9.331987e-01_r8,9.361131e-01_r8,9.394162e-01_r8, &
       9.431324e-01_r8/)
      ASYICE3(:,6) = (/ &
!    BAND 6
       8.792712e-01_r8,8.938493e-01_r8,9.015383e-01_r8,9.078755e-01_r8,9.134852e-01_r8, &
       9.185471e-01_r8,9.231387e-01_r8,9.273046e-01_r8,9.310755e-01_r8,9.344765e-01_r8, &
       9.375293e-01_r8,9.402543e-01_r8,9.426709e-01_r8,9.447982e-01_r8,9.466549e-01_r8, &
       9.482596e-01_r8,9.496311e-01_r8,9.507879e-01_r8,9.517486e-01_r8,9.525319e-01_r8, &
       9.531565e-01_r8,9.536411e-01_r8,9.540046e-01_r8,9.542656e-01_r8,9.544431e-01_r8, &
       9.545559e-01_r8,9.546229e-01_r8,9.546631e-01_r8,9.546954e-01_r8,9.547387e-01_r8, &
       9.548121e-01_r8,9.549345e-01_r8,9.551251e-01_r8,9.554028e-01_r8,9.557866e-01_r8, &
       9.562958e-01_r8,9.569494e-01_r8,9.577665e-01_r8,9.587663e-01_r8,9.599680e-01_r8, & 
       9.613908e-01_r8,9.630539e-01_r8,9.649767e-01_r8,9.671785e-01_r8,9.696786e-01_r8, &
       9.724966e-01_r8/)
      ASYICE3(:,7) = (/ &
!   BAND 7
       8.992856e-01_r8,9.093587e-01_r8,9.140642e-01_r8,9.183077e-01_r8,9.222541e-01_r8, &
       9.259455e-01_r8,9.294018e-01_r8,9.326363e-01_r8,9.356600e-01_r8,9.384824e-01_r8, &
       9.411129e-01_r8,9.435602e-01_r8,9.458332e-01_r8,9.479406e-01_r8,9.498907e-01_r8, &
       9.516924e-01_r8,9.533540e-01_r8,9.548842e-01_r8,9.562914e-01_r8,9.575842e-01_r8, &
       9.587710e-01_r8,9.598605e-01_r8,9.608610e-01_r8,9.617811e-01_r8,9.626294e-01_r8, &
       9.634142e-01_r8,9.641441e-01_r8,9.648277e-01_r8,9.654732e-01_r8,9.660894e-01_r8, &
       9.666845e-01_r8,9.672673e-01_r8,9.678460e-01_r8,9.684293e-01_r8,9.690257e-01_r8, &
       9.696436e-01_r8,9.702917e-01_r8,9.709786e-01_r8,9.717127e-01_r8,9.725029e-01_r8, &
       9.733577e-01_r8,9.742858e-01_r8,9.752962e-01_r8,9.763975e-01_r8,9.775987e-01_r8, &
       9.789088e-01_r8/)
      ASYICE3(:,8) = (/ &
!    BAND 8
       8.768357e-01_r8,8.832796e-01_r8,8.887092e-01_r8,8.937467e-01_r8,8.984873e-01_r8, &
       9.029645e-01_r8,9.071961e-01_r8,9.111946e-01_r8,9.149701e-01_r8,9.185316e-01_r8, &
       9.218877e-01_r8,9.250464e-01_r8,9.280157e-01_r8,9.308033e-01_r8,9.334170e-01_r8, &
       9.358643e-01_r8,9.381529e-01_r8,9.402904e-01_r8,9.422842e-01_r8,9.441418e-01_r8, &
       9.458708e-01_r8,9.474786e-01_r8,9.489727e-01_r8,9.503604e-01_r8,9.516494e-01_r8, &
       9.528469e-01_r8,9.539604e-01_r8,9.549973e-01_r8,9.559651e-01_r8,9.568711e-01_r8, &
       9.577229e-01_r8,9.585276e-01_r8,9.592929e-01_r8,9.600261e-01_r8,9.607346e-01_r8, &
       9.614258e-01_r8,9.621071e-01_r8,9.627859e-01_r8,9.634698e-01_r8,9.641660e-01_r8, &
       9.648820e-01_r8,9.656253e-01_r8,9.664032e-01_r8,9.672233e-01_r8,9.680930e-01_r8, &
       9.690198e-01_r8/)
      ASYICE3(:,9) = (/ &
!    BAND 9
       8.609711e-01_r8,8.676429e-01_r8,8.740194e-01_r8,8.800768e-01_r8,8.858172e-01_r8, &
       8.912476e-01_r8,8.963767e-01_r8,9.012137e-01_r8,9.057682e-01_r8,9.100498e-01_r8, &
       9.140683e-01_r8,9.178335e-01_r8,9.213552e-01_r8,9.246434e-01_r8,9.277079e-01_r8, &
       9.305587e-01_r8,9.332058e-01_r8,9.356590e-01_r8,9.379284e-01_r8,9.400240e-01_r8, &
       9.419558e-01_r8,9.437340e-01_r8,9.453685e-01_r8,9.468694e-01_r8,9.482469e-01_r8, &
       9.495111e-01_r8,9.506721e-01_r8,9.517401e-01_r8,9.527252e-01_r8,9.536376e-01_r8, &
       9.544876e-01_r8,9.552853e-01_r8,9.560409e-01_r8,9.567647e-01_r8,9.574670e-01_r8, &
       9.581578e-01_r8,9.588476e-01_r8,9.595466e-01_r8,9.602650e-01_r8,9.610130e-01_r8, &
       9.618010e-01_r8,9.626392e-01_r8,9.635379e-01_r8,9.645074e-01_r8,9.655578e-01_r8, &
       9.666994e-01_r8/)
      ASYICE3(:,10) = (/ &
!    BAND 10
       8.805232e-01_r8,8.872100e-01_r8,8.934945e-01_r8,8.993817e-01_r8,9.048833e-01_r8, &
       9.100128e-01_r8,9.147839e-01_r8,9.192106e-01_r8,9.233070e-01_r8,9.270873e-01_r8, &
       9.305656e-01_r8,9.337561e-01_r8,9.366732e-01_r8,9.393308e-01_r8,9.417434e-01_r8, &
       9.439252e-01_r8,9.458903e-01_r8,9.476530e-01_r8,9.492276e-01_r8,9.506283e-01_r8, &
       9.518694e-01_r8,9.529650e-01_r8,9.539295e-01_r8,9.547771e-01_r8,9.555221e-01_r8, &
       9.561786e-01_r8,9.567610e-01_r8,9.572836e-01_r8,9.577605e-01_r8,9.582059e-01_r8, &
       9.586343e-01_r8,9.590598e-01_r8,9.594966e-01_r8,9.599591e-01_r8,9.604615e-01_r8, &
       9.610179e-01_r8,9.616428e-01_r8,9.623503e-01_r8,9.631546e-01_r8,9.640701e-01_r8, &
       9.651110e-01_r8,9.662916e-01_r8,9.676260e-01_r8,9.691286e-01_r8,9.708136e-01_r8, &
       9.726952e-01_r8/)
      ASYICE3(:,11) = (/ &
!    BAND 11
       8.860103e-01_r8,8.920722e-01_r8,8.976788e-01_r8,9.029125e-01_r8,9.078046e-01_r8, &
       9.123744e-01_r8,9.166373e-01_r8,9.206068e-01_r8,9.242956e-01_r8,9.277160e-01_r8, &
       9.308801e-01_r8,9.337998e-01_r8,9.364867e-01_r8,9.389526e-01_r8,9.412092e-01_r8, &
       9.432681e-01_r8,9.451410e-01_r8,9.468395e-01_r8,9.483754e-01_r8,9.497603e-01_r8, &
       9.510059e-01_r8,9.521239e-01_r8,9.531262e-01_r8,9.540244e-01_r8,9.548305e-01_r8, &
       9.555563e-01_r8,9.562136e-01_r8,9.568142e-01_r8,9.573703e-01_r8,9.578935e-01_r8, &
       9.583959e-01_r8,9.588895e-01_r8,9.593861e-01_r8,9.598976e-01_r8,9.604362e-01_r8, &
       9.610135e-01_r8,9.616417e-01_r8,9.623325e-01_r8,9.630979e-01_r8,9.639496e-01_r8, &
       9.648996e-01_r8,9.659594e-01_r8,9.671409e-01_r8,9.684557e-01_r8,9.699153e-01_r8, &
       9.715312e-01_r8/)
      ASYICE3(:,12) = (/ &
!    BAND 12
       8.097123e-01_r8,8.172929e-01_r8,8.245083e-01_r8,8.314251e-01_r8,8.380674e-01_r8, &
       8.444481e-01_r8,8.505763e-01_r8,8.564596e-01_r8,8.621044e-01_r8,8.675169e-01_r8, &
       8.727030e-01_r8,8.776684e-01_r8,8.824186e-01_r8,8.869590e-01_r8,8.912951e-01_r8, &
       8.954321e-01_r8,8.993753e-01_r8,9.031299e-01_r8,9.067009e-01_r8,9.100936e-01_r8, &
       9.133130e-01_r8,9.163640e-01_r8,9.192519e-01_r8,9.219815e-01_r8,9.245578e-01_r8, &
       9.269860e-01_r8,9.292709e-01_r8,9.314175e-01_r8,9.334309e-01_r8,9.353161e-01_r8, &
       9.370781e-01_r8,9.387218e-01_r8,9.402523e-01_r8,9.416746e-01_r8,9.429938e-01_r8, &
       9.442149e-01_r8,9.453428e-01_r8,9.463827e-01_r8,9.473394e-01_r8,9.482180e-01_r8, &
       9.490235e-01_r8,9.497606e-01_r8,9.504344e-01_r8,9.510495e-01_r8,9.516106e-01_r8, &
       9.521225e-01_r8/)
      ASYICE3(:,13) = (/ &
!    BAND 13
       8.377231e-01_r8,8.469574e-01_r8,8.556844e-01_r8,8.638984e-01_r8,8.716085e-01_r8, &
       8.788282e-01_r8,8.855727e-01_r8,8.918582e-01_r8,8.977012e-01_r8,9.031189e-01_r8, &
       9.081287e-01_r8,9.127481e-01_r8,9.169950e-01_r8,9.208873e-01_r8,9.244432e-01_r8, &
       9.276807e-01_r8,9.306183e-01_r8,9.332744e-01_r8,9.356674e-01_r8,9.378159e-01_r8, &
       9.397385e-01_r8,9.414539e-01_r8,9.429808e-01_r8,9.443380e-01_r8,9.455443e-01_r8, &
       9.466185e-01_r8,9.475794e-01_r8,9.484460e-01_r8,9.492372e-01_r8,9.499718e-01_r8, &
       9.506689e-01_r8,9.513473e-01_r8,9.520260e-01_r8,9.527240e-01_r8,9.534602e-01_r8, &
       9.542535e-01_r8,9.551230e-01_r8,9.560875e-01_r8,9.571659e-01_r8,9.583773e-01_r8, &
       9.597404e-01_r8,9.612742e-01_r8,9.629975e-01_r8,9.649292e-01_r8,9.670880e-01_r8, &
       9.694928e-01_r8/)
      ASYICE3(:,14) = (/ &
!    BAND 14
       8.344994e-01_r8,8.443682e-01_r8,8.536201e-01_r8,8.622890e-01_r8,8.703985e-01_r8, &
       8.779698e-01_r8,8.850230e-01_r8,8.915780e-01_r8,8.976545e-01_r8,9.032725e-01_r8, &
       9.084518e-01_r8,9.132124e-01_r8,9.175742e-01_r8,9.215575e-01_r8,9.251823e-01_r8, &
       9.284689e-01_r8,9.314374e-01_r8,9.341082e-01_r8,9.365017e-01_r8,9.386382e-01_r8, &
       9.405382e-01_r8,9.422220e-01_r8,9.437101e-01_r8,9.450231e-01_r8,9.461814e-01_r8, &
       9.472056e-01_r8,9.481162e-01_r8,9.489338e-01_r8,9.496790e-01_r8,9.503724e-01_r8, &
       9.510345e-01_r8,9.516860e-01_r8,9.523475e-01_r8,9.530396e-01_r8,9.537830e-01_r8, &
       9.545982e-01_r8,9.555059e-01_r8,9.565267e-01_r8,9.576813e-01_r8,9.589902e-01_r8, &
       9.604741e-01_r8,9.621535e-01_r8,9.640490e-01_r8,9.661813e-01_r8,9.685709e-01_r8, &
       9.712382e-01_r8/)
      ASYICE3(:,15) = (/ &
!    BAND 15
       7.853957e-01_r8,7.968767e-01_r8,8.077401e-01_r8,8.180018e-01_r8,8.276797e-01_r8, &
       8.367923e-01_r8,8.453586e-01_r8,8.533977e-01_r8,8.609288e-01_r8,8.679711e-01_r8, &
       8.745440e-01_r8,8.806667e-01_r8,8.863587e-01_r8,8.916394e-01_r8,8.965282e-01_r8, &
       9.010445e-01_r8,9.052079e-01_r8,9.090379e-01_r8,9.125540e-01_r8,9.157758e-01_r8, &
       9.187228e-01_r8,9.214147e-01_r8,9.238711e-01_r8,9.261117e-01_r8,9.281562e-01_r8, &
       9.300242e-01_r8,9.317354e-01_r8,9.333098e-01_r8,9.347669e-01_r8,9.361267e-01_r8, &
       9.374088e-01_r8,9.386331e-01_r8,9.398195e-01_r8,9.409877e-01_r8,9.421576e-01_r8, &
       9.433489e-01_r8,9.445817e-01_r8,9.458755e-01_r8,9.472504e-01_r8,9.487260e-01_r8, &
       9.503221e-01_r8,9.520585e-01_r8,9.539550e-01_r8,9.560313e-01_r8,9.583069e-01_r8, &
       9.608016e-01_r8/)
      ASYICE3(:,16) = (/ &
!    BAND 16
       8.340752e-01_r8,8.435170e-01_r8,8.517487e-01_r8,8.592064e-01_r8,8.660387e-01_r8, &
       8.723204e-01_r8,8.780997e-01_r8,8.834137e-01_r8,8.882934e-01_r8,8.927662e-01_r8, &
       8.968577e-01_r8,9.005914e-01_r8,9.039899e-01_r8,9.070745e-01_r8,9.098659e-01_r8, &
       9.123836e-01_r8,9.146466e-01_r8,9.166734e-01_r8,9.184817e-01_r8,9.200886e-01_r8, &
       9.215109e-01_r8,9.227648e-01_r8,9.238661e-01_r8,9.248304e-01_r8,9.256727e-01_r8, &
       9.264078e-01_r8,9.270505e-01_r8,9.276150e-01_r8,9.281156e-01_r8,9.285662e-01_r8, &
       9.289806e-01_r8,9.293726e-01_r8,9.297557e-01_r8,9.301435e-01_r8,9.305491e-01_r8, &
       9.309859e-01_r8,9.314671e-01_r8,9.320055e-01_r8,9.326140e-01_r8,9.333053e-01_r8, &
       9.340919e-01_r8,9.349861e-01_r8,9.360000e-01_r8,9.371451e-01_r8,9.384329e-01_r8, &
       9.398744e-01_r8/)
! --------------linhan add end---------------------------------

      end subroutine lwcldpr

      end module rrtmg_lw_init

