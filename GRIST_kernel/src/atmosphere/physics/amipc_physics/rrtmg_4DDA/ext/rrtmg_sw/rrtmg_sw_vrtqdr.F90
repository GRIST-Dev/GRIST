!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_vrtqdr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:15 $
!
      module rrtmg_sw_vrtqdr

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

!      use parkind, only: jpim, jprb
!      use parrrsw, only: ngptsw

      implicit none

      contains

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw(klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      integer, intent (in) :: klev                   ! number of model layers
      integer, intent (in) :: kw                     ! g-point index

      real(kind=r8), intent(in) :: pref(:)                    ! direct beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: prefd(:)                   ! diffuse beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptra(:)                    ! direct beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptrad(:)                   ! diffuse beam transmissivity
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(in) :: pdbt(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptdbt(:)
                                                                 !   Dimensions: (nlayers+1)

      real(kind=r8), intent(inout) :: prdnd(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prup(:)
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prupd(:)
                                                                 !   Dimensions: (nlayers+1)

! Output
      real(kind=r8), intent(out) :: pfd(:,:)                   ! downwelling flux (W/m2)
                                                                 !   Dimensions: (ngptsw,nlayers+1)
                                                                 ! unadjusted for earth/sun distance or zenith angle
      real(kind=r8), intent(out) :: pfu(:,:)                   ! upwelling flux (W/m2)
                                                                 !   Dimensions: (ngptsw,nlayers+1)
                                                                 ! unadjusted for earth/sun distance or zenith angle

! Local

      integer :: ikp, ikx, jk

      real(kind=r8) :: zreflect
      real(kind=r8) :: ztdn(klev+1)  

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------
                   
! Link lowest layer with surface
             
      zreflect = 1._r8 / (1._r8 - prefd(klev+1) * prefd(klev))
      prup(klev) = pref(klev) + (ptrad(klev) * &
                 ((ptra(klev) - pdbt(klev)) * prefd(klev+1) + &
                   pdbt(klev) * pref(klev+1))) * zreflect
      prupd(klev) = prefd(klev) + ptrad(klev) * ptrad(klev) * &
                    prefd(klev+1) * zreflect

! Pass from bottom to top 

      do jk = 1,klev-1
         ikp = klev+1-jk                       
         ikx = ikp-1
         zreflect = 1._r8 / (1._r8 -prupd(ikp) * prefd(ikx))
         prup(ikx) = pref(ikx) + (ptrad(ikx) * &
                   ((ptra(ikx) - pdbt(ikx)) * prupd(ikp) + &
                     pdbt(ikx) * prup(ikp))) * zreflect
         prupd(ikx) = prefd(ikx) + ptrad(ikx) * ptrad(ikx) * &
                      prupd(ikp) * zreflect
      enddo
    
! Upper boundary conditions

      ztdn(1) = 1._r8
      prdnd(1) = 0._r8
      ztdn(2) = ptra(1)
      prdnd(2) = prefd(1)

! Pass from top to bottom

      do jk = 2,klev
         ikp = jk+1
         zreflect = 1._r8 / (1._r8 - prefd(jk) * prdnd(jk))
         ztdn(ikp) = ptdbt(jk) * ptra(jk) + &
                    (ptrad(jk) * ((ztdn(jk) - ptdbt(jk)) + &
                     ptdbt(jk) * pref(jk) * prdnd(jk))) * zreflect
         prdnd(ikp) = prefd(jk) + ptrad(jk) * ptrad(jk) * &
                      prdnd(jk) * zreflect
      enddo
    
! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         zreflect = 1._r8 / (1._r8 - prdnd(jk) * prupd(jk))
         pfu(kw,jk) = (ptdbt(jk) * prup(jk) + &
                      (ztdn(jk) - ptdbt(jk)) * prupd(jk)) * zreflect
         pfd(kw,jk) = ptdbt(jk) + (ztdn(jk) - ptdbt(jk)+ &
                      ptdbt(jk) * prup(jk) * prdnd(jk)) * zreflect
      enddo

      end subroutine vrtqdr_sw

      subroutine vrtqdr_4DDA_sw(klev, kw, &
                           RU0, RBAR, Tu0, TBAR, DTR, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
						   albedop, albedod, &
                           FLXD, FLXU)
! --------------------------------------------------------------------------
 
! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
! 
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input
      integer, intent (in) :: klev                   ! number of model layers
      integer, intent (in) :: kw                     ! g-point index

      real(kind=r8), intent(inout) :: RU0(:,:)                 ! direct beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: RBAR(:,:,:)              ! diffuse beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: Tu0(:,:)                 ! direct beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: TBAR(:,:,:)              ! diffuse beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: DTR(:)
	   
	  real(kind=r8), intent(in) :: pdbt(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: ptdbt(:)
                                                              !   Dimensions: (nlayers+1)

      real(kind=r8), intent(inout) :: prdnd(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prup(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prupd(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=r8), intent(in) :: albedod                    ! surface albedo (diffuse)
                                                               !   Dimensions: (nbndsw)
      real(kind=r8), intent(in) :: albedop                   ! surface albedo (direct)
                                                               !   Dimensions: (nbndsw)

! Output
      real(kind=r8), intent(out) :: FLXD(:,:)                  ! downwelling flux (W/m2)
                                                              !   Dimensions: (ngptsw,nlayers+1)
                                                              ! unadjusted for earth/sun distance or zenith angle
      real(kind=r8), intent(out) :: FLXU(:,:)                  ! upwelling flux (W/m2)
                                                              !   Dimensions: (ngptsw,nlayers+1)
                                                              ! unadjusted for earth/sun distance or zenith angle

! Local

      integer :: ikp, ikx, jk, K, KM1, L, LP1
      real(kind=r8) :: RBAR1(klev+1,2,2),TBAR1(klev+1,2,2),RBARS(klev+1,2,2),TBARS1(klev+1,2,2),RBARN(klev+1,2,2),TBARN(klev+1,2,2)
	  real(kind=r8) :: RU01(klev+1,2),TU01(klev+1,2),RU0N(klev+1,2),TU0N(klev+1,2),TDIR(klev+1)
	  real(kind=r8) :: C11,C12,C21,C22,PMOD,TMM11,TMM12,TMM21,TMM22,TMR11,TMR12,TMR21,TMR22
	  real(kind=r8) :: RT1,RT2,UU1,UU2,DU1,DU2,BRT11,BRT12,BRT21,BRT22,BMR11,BMR12,BMR21,BMR22 

	real(kind=r8) :: a(4),u(4)  	          
 	data a / 0.5, 0.5, 0.5, 0.5 /
	data u / -0.7886752, -0.2113247, 0.2113247, 0.7886752/    
! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------

! ... initializaton for the TOA.   
        RBARS(1,1:2,1)= 0.0
        RBARS(1,1:2,2)= 0.0
        Tu01(1,1)= 0.0
        Tu01(1,2)= 0.0
        TDIR(1)=1.0

! ... initializaton for the lev layer (surface).                       

        Ru0N(klev+1,1:2)  = albedop
        RbarN(klev+1,1:2,1) = albedod *u(3)*2*a(3)
        RbarN(klev+1,1:2,2) = albedod *u(4)*2*a(4)                 

! ...  ADD THE LAYERS UPWARD FROM ONE LAYER ABOVE SURFACE TO TOA. 

      DO K = 2, klev+1
        KM1 = K - 1
        C11           =  1. - RBAR(KM1,1,1) * RBARS(KM1,1,1) -RBAR(KM1,1,2) * RBARS(KM1,2,1)
        C12           =     - RBAR(KM1,1,1) * RBARS(KM1,1,2) -RBAR(KM1,1,2) * RBARS(KM1,2,2)
        C21           =     - RBAR(KM1,2,1) * RBARS(KM1,1,1) -RBAR(KM1,2,2) * RBARS(KM1,2,1)
        C22           =  1. - RBAR(KM1,2,1) * RBARS(KM1,1,2) -RBAR(KM1,2,2) * RBARS(KM1,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD

        BMR11         =  TMM11 * TBAR(KM1,1,1) + TMM12 * TBAR(KM1,2,1)
        BMR12         =  TMM11 * TBAR(KM1,1,2) + TMM12 * TBAR(KM1,2,2)
        BMR21         =  TMM21 * TBAR(KM1,1,1) + TMM22 * TBAR(KM1,2,1)
        BMR22         =  TMM21 * TBAR(KM1,1,2) + TMM22 * TBAR(KM1,2,2)

        TDIR(K)       =  TDIR(KM1) * DTR(KM1)
        RT1           =  RU0(KM1,1) * TDIR(KM1) + RBAR(KM1,1,1) *TU01(KM1,1) + RBAR(KM1,1,2) * TU01(KM1,2)
        RT2           =  RU0(KM1,2) * TDIR(KM1) + RBAR(KM1,2,1) *TU01(KM1,1) + RBAR(KM1,2,2) * TU01(KM1,2)
        UU1           =  TMM11 * RT1 + TMM12 * RT2
        UU2           =  TMM21 * RT1 + TMM22 * RT2
        DU1           =  TU01(KM1,1) + RBARS(KM1,1,1) * UU1 +RBARS(KM1,1,2) * UU2
        DU2           =  TU01(KM1,2) + RBARS(KM1,2,1) * UU1 +RBARS(KM1,2,2) * UU2
        TU01(K,1)     =  TU0(KM1,1) * TDIR(KM1) + TBAR(KM1,1,1) *DU1 + TBAR(KM1,1,2) * DU2
        TU01(K,2)     =  TU0(KM1,2) * TDIR(KM1) + TBAR(KM1,2,1) *DU1 + TBAR(KM1,2,2) * DU2

        BRT11         =  TBAR(KM1,1,1) * RBARS(KM1,1,1) +TBAR(KM1,1,2) * RBARS(KM1,2,1)
        BRT12         =  TBAR(KM1,1,1) * RBARS(KM1,1,2) +TBAR(KM1,1,2) * RBARS(KM1,2,2)
        BRT21         =  TBAR(KM1,2,1) * RBARS(KM1,1,1) +TBAR(KM1,2,2) * RBARS(KM1,2,1)
        BRT22         =  TBAR(KM1,2,1) * RBARS(KM1,1,2) +TBAR(KM1,2,2) * RBARS(KM1,2,2)
        RBARS(K,1,1)  =  RBAR(KM1,1,1) + BRT11 * BMR11 + BRT12 * BMR21
        RBARS(K,1,2)  =  RBAR(KM1,1,2) + BRT11 * BMR12 + BRT12 * BMR22
        RBARS(K,2,1)  =  RBAR(KM1,2,1) + BRT21 * BMR11 + BRT22 * BMR21
        RBARS(K,2,2)  =  RBAR(KM1,2,2) + BRT21 * BMR12 + BRT22 * BMR22

! ... ADD THE LAYERS UPWARD FROM LAYER ABOVE SURFACE TO THE LEV1.

        L = klev+1 - K + 1
        LP1 = L + 1
        C11           =  1. - RBARN(LP1,1,1) * RBAR(L,1,1) - RBARN(LP1,1,2) * RBAR(L,2,1)
        C12           =     - RBARN(LP1,1,1) * RBAR(L,1,2) - RBARN(LP1,1,2) * RBAR(L,2,2)
        C21           =     - RBARN(LP1,2,1) * RBAR(L,1,1) - RBARN(LP1,2,2) * RBAR(L,2,1)
        C22           =  1. - RBARN(LP1,2,1) * RBAR(L,1,2) - RBARN(LP1,2,2) * RBAR(L,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD
        RT1           =  RU0N(LP1,1) * DTR(L) + RBARN(LP1,1,1) * TU0(L,1) + RBARN(LP1,1,2) * TU0(L,2)
        RT2           =  RU0N(LP1,2) * DTR(L) + RBARN(LP1,2,1) * TU0(L,1) + RBARN(LP1,2,2) * TU0(L,2)
        UU1           =  TMM11 * RT1 + TMM12 * RT2
        UU2           =  TMM21 * RT1 + TMM22 * RT2
        DU1           =  TU0(L,1) + RBAR(L,1,1) * UU1 + RBAR(L,1,2) * UU2
        DU2           =  TU0(L,2) + RBAR(L,2,1) * UU1 + RBAR(L,2,2) * UU2
        RU0N(L,1)     =  RU0(L,1) + TBAR(L,1,1) * UU1 + TBAR(L,1,2) * UU2
        RU0N(L,2)     =  RU0(L,2) + TBAR(L,2,1) * UU1 + TBAR(L,2,2) * UU2

        C11           =   1. - RBAR(L,1,1) * RBARN(LP1,1,1) - RBAR(L,1,2) * RBARN(LP1,2,1)
        C12           =      - RBAR(L,1,1) * RBARN(LP1,1,2) - RBAR(L,1,2) * RBARN(LP1,2,2)
        C21           =      - RBAR(L,2,1) * RBARN(LP1,1,1) - RBAR(L,2,2) * RBARN(LP1,2,1)
        C22           =   1. - RBAR(L,2,1) * RBARN(LP1,1,2) - RBAR(L,2,2) * RBARN(LP1,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD

        TMR11         =  TBAR(L,1,1) * (RBARN(LP1,1,1) * TMM11 + RBARN(LP1,1,2) * TMM21) + &
                         TBAR(L,1,2) * (RBARN(LP1,2,1) * TMM11 + RBARN(LP1,2,2) * TMM21)
        TMR12         =  TBAR(L,1,1) * (RBARN(LP1,1,1) * TMM12 + RBARN(LP1,1,2) * TMM22) + &
                         TBAR(L,1,2) * (RBARN(LP1,2,1) * TMM12 + RBARN(LP1,2,2) * TMM22)
        TMR21         =  TBAR(L,2,1) * (RBARN(LP1,1,1) * TMM11 + RBARN(LP1,1,2) * TMM21) + &
                         TBAR(L,2,2) * (RBARN(LP1,2,1) * TMM11 + RBARN(LP1,2,2) * TMM21)
        TMR22         =  TBAR(L,2,1) * (RBARN(LP1,1,1) * TMM12 + RBARN(LP1,1,2) * TMM22) + &
                         TBAR(L,2,2) * (RBARN(LP1,2,1) * TMM12 + RBARN(LP1,2,2) * TMM22)
 
        RBARN(L,1,1)  =  RBAR(L,1,1) + TMR11 * TBAR(L,1,1) + TMR12 * TBAR(L,2,1)
        RBARN(L,1,2)  =  RBAR(L,1,2) + TMR11 * TBAR(L,1,2) + TMR12 * TBAR(L,2,2)
        RBARN(L,2,1)  =  RBAR(L,2,1) + TMR21 * TBAR(L,1,1) + TMR22 * TBAR(L,2,1)
        RBARN(L,2,2)  =  RBAR(L,2,2) + TMR21 * TBAR(L,1,2) + TMR22 * TBAR(L,2,2)
  end do

!----------------------------------------------------------------------C
!     ADD DOWNWARD TO CALCULATE THE RESULTANT REFLECTANCES AND         C
!     TRANSMITTANCE AT FLUX LEVELS.                                    C
!----------------------------------------------------------------------C

      DO K = 1, klev+1
        C11           =  1. - RBARN(K,1,1) * RBARS(K,1,1) - RBARN(K,1,2) * RBARS(K,2,1)
        C12           =     - RBARN(K,1,1) * RBARS(K,1,2) - RBARN(K,1,2) * RBARS(K,2,2)
        C21           =     - RBARN(K,2,1) * RBARS(K,1,1) - RBARN(K,2,2) * RBARS(K,2,1)
        C22           =  1. - RBARN(K,2,1) * RBARS(K,1,2) - RBARN(K,2,2) * RBARS(K,2,2)
        PMOD          =    C11 * C22 - C12 * C21
        TMM11         =    C22 / PMOD
        TMM12         =  - C12 / PMOD
        TMM21         =  - C21 / PMOD
        TMM22         =    C11 / PMOD

        RT1           =  RU0N(K,1) * TDIR(K) + RBARN(K,1,1) *TU01(K,1) + RBARN(K,1,2) * TU01(K,2)
        RT2           =  RU0N(K,2) * TDIR(K) + RBARN(K,2,1) * TU01(K,1) + RBARN(K,2,2) * TU01(K,2)
        UU1           =  TMM11 * RT1 + TMM12 * RT2
        UU2           =  TMM21 * RT1 + TMM22 * RT2
        DU1           =  TU01(K,1) + RBARS(K,1,1) * UU1 + RBARS(K,1,2) * UU2
        DU2           =  TU01(K,2) + RBARS(K,2,1) * UU1 + RBARS(K,2,2) * UU2

        FLXU(kw,K)       =  2*(a(3)*u(3)*UU1+a(4)*u(4)*UU2)
        FLXD(kw,K)       =  2*(a(3)*u(3)*DU1+a(4)*u(4)*DU2)+TDIR(K)

!open(10,file='FLUD.txt',form='formatted')
!write(10,*) "Ru0",K,RU0(K,1),RU0(K,2),TU0(K,1),TU0(K,2)
!write(10,*) "Rbar",K,RBAR(K,1,1),RBAR(K,1,2),RBAR(K,2,1),RBAR(K,2,2)
!write(10,*) "Flux",K,FLXD(K,kw),FLXU(K,kw)
!write(10,*) "UU",K,UU1,UU2,DU1,DU2

  end do


      end subroutine vrtqdr_4DDA_sw
      end module rrtmg_sw_vrtqdr
