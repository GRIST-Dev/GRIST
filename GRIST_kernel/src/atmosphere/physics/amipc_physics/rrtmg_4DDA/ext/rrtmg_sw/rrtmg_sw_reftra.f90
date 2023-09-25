!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_reftra.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.2 $
!     created:   $Date: 2007/08/23 20:40:14 $

      module rrtmg_sw_reftra

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
      use rrsw_tbl, only : tblint, bpade, od_lo, exp_tbl
      use rrsw_vsn, only : hvrrft, hnamrft

      implicit none

      contains

! --------------------------------------------------------------------
      subroutine reftra_sw(nlayers, lrtchk, pgg, prmuz, ptau, pw, &
                           pref, prefd, ptra, ptrad)
! --------------------------------------------------------------------
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      lrtchk  = .t. for all layers in clear profile
!      lrtchk  = .t. for cloudy layers in cloud profile 
!              = .f. for clear layers in cloud profile
!      pgg     = assymetry factor
!      prmuz   = cosine solar zenith angle
!      ptau    = optical thickness
!      pw      = single scattering albedo
!
! outputs
! -------
!      pref    : collimated beam reflectivity
!      prefd   : diffuse beam reflectivity 
!      ptra    : collimated beam transmissivity
!      ptrad   : diffuse beam transmissivity
!
!
! Method:
! -------
!      standard delta-eddington, p.i.f.m., or d.o.m. layer calculations.
!      kmodts  = 1 eddington (joseph et al., 1976)
!              = 2 pifm (zdunkowski et al., 1980)
!              = 3 discrete ordinates (liou, 1973)
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer, intent(in) :: nlayers

      logical, intent(in) :: lrtchk(:)                           ! Logical flag for reflectivity and
                                                                 ! and transmissivity calculation; 
                                                                 !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: pgg(:)                      ! asymmetry parameter
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: ptau(:)                     ! optical depth
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: pw(:)                       ! single scattering albedo 
                                                                 !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: prmuz                       ! cosine of solar zenith angle

! ------- Output -------

      real(kind=r8), intent(inout) :: pref(:)                    ! direct beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: prefd(:)                   ! diffuse beam reflectivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: ptra(:)                    ! direct beam transmissivity
                                                                 !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: ptrad(:)                   ! diffuse beam transmissivity
                                                                 !   Dimensions: (nlayers+1)

! ------- Local -------

      integer :: jk, jl, kmodts
      integer :: itind

      real(kind=r8) :: tblind
      real(kind=r8) :: za, za1, za2
      real(kind=r8) :: zbeta, zdend, zdenr, zdent
      real(kind=r8) :: ze1, ze2, zem1, zem2, zemm, zep1, zep2
      real(kind=r8) :: zg, zg3, zgamma1, zgamma2, zgamma3, zgamma4, zgt
      real(kind=r8) :: zr1, zr2, zr3, zr4, zr5
      real(kind=r8) :: zrk, zrk2, zrkg, zrm1, zrp, zrp1, zrpp
      real(kind=r8) :: zsr3, zt1, zt2, zt3, zt4, zt5, zto1
      real(kind=r8) :: zw, zwcrit, zwo

      real(kind=r8), parameter :: eps = 1.e-08_r8

!     ------------------------------------------------------------------

! Initialize

      hvrrft = '$Revision: 1.2 $'

      zsr3=sqrt(3._r8)
      zwcrit=0.9999995_r8
      kmodts=2

      do jk=1, nlayers
         if (.not.lrtchk(jk)) then
            pref(jk) =0._r8
            ptra(jk) =1._r8
            prefd(jk)=0._r8
            ptrad(jk)=1._r8
         else
            zto1=ptau(jk)
            zw  =pw(jk)
            zg  =pgg(jk)  

! General two-stream expressions

            zg3= 3._r8 * zg
            if (kmodts == 1) then
               zgamma1= (7._r8 - zw * (4._r8 + zg3)) * 0.25_r8
               zgamma2=-(1._r8 - zw * (4._r8 - zg3)) * 0.25_r8
               zgamma3= (2._r8 - zg3 * prmuz ) * 0.25_r8
            else if (kmodts == 2) then  
               zgamma1= (8._r8 - zw * (5._r8 + zg3)) * 0.25_r8
               zgamma2=  3._r8 *(zw * (1._r8 - zg )) * 0.25_r8
               zgamma3= (2._r8 - zg3 * prmuz ) * 0.25_r8
            else if (kmodts == 3) then  
               zgamma1= zsr3 * (2._r8 - zw * (1._r8 + zg)) * 0.5_r8
               zgamma2= zsr3 * zw * (1._r8 - zg ) * 0.5_r8
               zgamma3= (1._r8 - zsr3 * zg * prmuz ) * 0.5_r8
            end if
            zgamma4= 1._r8 - zgamma3
    
! Recompute original s.s.a. to test for conservative solution

            zwo= zw / (1._r8 - (1._r8 - zw) * (zg / (1._r8 - zg))**2)
    
            if (zwo >= zwcrit) then
! Conservative scattering

               za  = zgamma1 * prmuz 
               za1 = za - zgamma3
               zgt = zgamma1 * zto1
        
! Homogeneous reflectance and transmittance,
! collimated beam

               ze1 = min ( zto1 / prmuz , 500._r8)
!               ze2 = exp( -ze1 )

! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  ze2 = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  ze2 = exp_tbl(itind)
               endif
!

               pref(jk) = (zgt - za1 * (1._r8 - ze2)) / (1._r8 + zgt)
               ptra(jk) = 1._r8 - pref(jk)

! isotropic incidence

               prefd(jk) = zgt / (1._r8 + zgt)
               ptrad(jk) = 1._r8 - prefd(jk)        

! This is applied for consistency between total (delta-scaled) and direct (unscaled) 
! calculations at very low optical depths (tau < 1.e-4) when the exponential lookup
! table returns a transmittance of 1.0.
               if (ze2 .eq. 1.0_r8) then 
                  pref(jk) = 0.0_r8
                  ptra(jk) = 1.0_r8
                  prefd(jk) = 0.0_r8
                  ptrad(jk) = 1.0_r8
               endif

            else
! Non-conservative scattering

               za1 = zgamma1 * zgamma4 + zgamma2 * zgamma3
               za2 = zgamma1 * zgamma3 + zgamma2 * zgamma4
               zrk = sqrt ( zgamma1**2 - zgamma2**2)
               zrp = zrk * prmuz               
               zrp1 = 1._r8 + zrp
               zrm1 = 1._r8 - zrp
               zrk2 = 2._r8 * zrk
               zrpp = 1._r8 - zrp*zrp
               zrkg = zrk + zgamma1
               zr1  = zrm1 * (za2 + zrk * zgamma3)
               zr2  = zrp1 * (za2 - zrk * zgamma3)
               zr3  = zrk2 * (zgamma3 - za2 * prmuz )
               zr4  = zrpp * zrkg
               zr5  = zrpp * (zrk - zgamma1)
               zt1  = zrp1 * (za1 + zrk * zgamma4)
               zt2  = zrm1 * (za1 - zrk * zgamma4)
               zt3  = zrk2 * (zgamma4 + za1 * prmuz )
               zt4  = zr4
               zt5  = zr5
               zbeta = (zgamma1 - zrk) / zrkg !- zr5 / zr4
        
! Homogeneous reflectance and transmittance

               ze1 = min ( zrk * zto1, 500._r8)
               ze2 = min ( zto1 / prmuz , 500._r8)
!
! Original
!              zep1 = exp( ze1 )
!              zem1 = exp(-ze1 )
!              zep2 = exp( ze2 )
!              zem2 = exp(-ze2 )
!
! Revised original, to reduce exponentials
!              zep1 = exp( ze1 )
!              zem1 = 1._r8 / zep1
!              zep2 = exp( ze2 )
!              zem2 = 1._r8 / zep2
!
! Use exponential lookup table for transmittance, or expansion of 
! exponential for low tau
               if (ze1 .le. od_lo) then 
                  zem1 = 1._r8 - ze1 + 0.5_r8 * ze1 * ze1
                  zep1 = 1._r8 / zem1
               else
                  tblind = ze1 / (bpade + ze1)
                  itind = tblint * tblind + 0.5_r8
                  zem1 = exp_tbl(itind)
                  zep1 = 1._r8 / zem1
               endif

               if (ze2 .le. od_lo) then 
                  zem2 = 1._r8 - ze2 + 0.5_r8 * ze2 * ze2
                  zep2 = 1._r8 / zem2
               else
                  tblind = ze2 / (bpade + ze2)
                  itind = tblint * tblind + 0.5_r8
                  zem2 = exp_tbl(itind)
                  zep2 = 1._r8 / zem2
               endif

! collimated beam

               zdenr = zr4*zep1 + zr5*zem1
               zdent = zt4*zep1 + zt5*zem1
               if (zdenr .ge. -eps .and. zdenr .le. eps) then
                  pref(jk) = eps
                  ptra(jk) = zem2
               else
                  pref(jk) = zw * (zr1*zep1 - zr2*zem1 - zr3*zem2) / zdenr
                  ptra(jk) = zem2 - zem2 * zw * (zt1*zep1 - zt2*zem1 - zt3*zep2) / zdent
               endif 

! diffuse beam

               zemm = zem1*zem1
               zdend = 1._r8 / ( (1._r8 - zbeta*zemm ) * zrkg)
               prefd(jk) =  zgamma2 * (1._r8 - zemm) * zdend
               ptrad(jk) =  zrk2*zem1*zdend

            endif

         endif         

      enddo    

      end subroutine reftra_sw
! --------------------------4DDA linhan--------------------------------------
      subroutine reftra_4DDA_sw(nlayers, lrtchk, pgg,pw2,pw3, AMU0, OTAU, OOM, &
                           RU0,RBAR,Tu0,TBAR,TDIR)
! --------------------------------------------------------------------
  
! Purpose: computes the reflectivity and transmissivity of a clear or 
!   cloudy layer using a choice of various approximations.
!
! Interface:  *rrtmg_sw_reftra* is called by *rrtmg_sw_spcvrt*
!
! Description:
! explicit arguments :
! --------------------
! inputs
! ------ 
!      lrtchk  = .t. for all layers in clear profile
!      lrtchk  = .t. for cloudy layers in cloud profile 
!              = .f. for clear layers in cloud profile
!      pgg     = assymetry factor
!      AMU0    = cosine solar zenith angle
!      OTAU    = optical thickness
!      OOM     = single scattering albedo
!
! outputs
! -------
!      RU0    : collimated beam reflectivity
!      RBAR   : diffuse beam reflectivity 
!      Tu0    : collimated beam transmissivity
!      TBAR   : diffuse beam transmissivity
!
!
! Method:
! -------
! SW 4DDA
!
!
! Modifications:
! --------------
! Original: J-JMorcrette, ECMWF, Feb 2003
! Revised for F90 reformatting: MJIacono, AER, Jul 2006
! Revised to add exponential lookup table: MJIacono, AER, Aug 2007
!
! ------------------------------------------------------------------

! ------- Declarations ------

! ------- Input -------

      integer, intent(in) :: nlayers

      logical, intent(in) :: lrtchk(:)                         ! Logical flag for reflectivity and
                                                               ! and transmissivity calculation; 
                                                               !   Dimensions: (nlayers)

      real(kind=r8), intent(in) :: pgg(:)                      ! asymmetry parameter
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: pw2(:)                      ! asymmetry parameter
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: pw3(:)                      ! asymmetry parameter
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: OTAU(:)                     ! optical depth
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: OOM(:)                      ! single scattering albedo 
                                                               !   Dimensions: (nlayers)
      real(kind=r8), intent(in) :: AMU0                        ! cosine of solar zenith angle

! ------- Output -------
													 
      real(kind=r8), intent(inout) :: RU0(:,:)                 ! direct beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: RBAR(:,:,:)              ! diffuse beam reflectivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: Tu0(:,:)                 ! direct beam transmissivity
                                                               !   Dimensions: (nlayers+1)
      real(kind=r8), intent(inout) :: TBAR(:,:,:)              ! diffuse beam transmissivity
                                                               !   Dimensions: (nlayers+1)
	  real(kind=r8), intent(inout) :: TDIR(:) 
! ------- Local -------

      integer :: jk

      real(kind=r8) :: D9,F,FW,OM,TAU,OMEGA1,OMEGA2,OMEGA3
      real(kind=r8) :: RMU0,RMU1,RMU2,eft0,eft1,eft2
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

      real(kind=r8), parameter :: eps = 1.e-08_r8

	real(kind=r8) :: a(4),u(4)  	          
 	data a / 0.5, 0.5, 0.5, 0.5 /
	data u / -0.7886752, -0.2113247, 0.2113247, 0.7886752/    
!     ------------------------------------------------------------------
! Initialize

      hvrrft = '$Revision: 29812 $'

      do jk=1,nlayers
         if (.not.lrtchk(jk)) then
            RU0(jk,:) =0._r8
            Tu0(jk,:) =1._r8
            RBAR(jk,:,:) =0._r8
            TBAR(jk,:,1) =1._r8 *2*u(3)*a(3)
			TBAR(jk,:,2) =1._r8 *2*u(4)*a(4)
         else
!----------------- adjust -----------------------------------------
!        F            =  0._r8!pgg(jk)**4
!	    FW           =  1.0_r8 - OOM(jk) * F
        OM           =  OOM(jk) !* (1.0_r8 - F) / FW
		OM           =  max(OM,1.E-8_r8)
		OM           =  min(OM,1._r8-1.E-8_r8)
        TAU          =  OTAU(jk)! * FW
        if (TAU<1.E-15_r8) TAU=1.E-15_r8

!        OMEGA1       =  (3._r8* pgg(jk) - 3.0_r8 * F) / (1.0_r8 - F)
!        OMEGA2       =  (5._r8* pw2(jk) - 5.0_r8 * F) / (1.0_r8 - F)
!        OMEGA3       =  (7._r8* pw3(jk) - 7.0_r8 * F) / (1.0_r8 - F)
        OMEGA1       =  3._r8* pgg(jk) !- 3.0_r8 * F) / (1.0_r8 - F)
        OMEGA2       =  5._r8* pw2(jk) !- 5.0_r8 * F) / (1.0_r8 - F)
        OMEGA3       =  7._r8* pw3(jk) !- 7.0_r8 * F) / (1.0_r8 - F)
!-----------------------------------------------------------------
        RMU0         =  1.0_r8 / AMU0   ! f0
        RMU1         =  1.0_r8 / u(3)   ! f1
        RMU2         =  1.0_r8 / u(4)   ! f2
	    eft0         =  exp(-RMU0*tau)
	    eft1         =  exp(-RMU1*tau)
	    eft2         =  exp(-RMU2*tau)
!-----------------------------------------------------------------
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

	    bu11         =  of0/u(1)*( 1.0_r8+OMEGA1*u(1)*(-u(3))+OMEGA2*p2u1*p2u2+OMEGA3*p3u1*p3u2 ) 
        bu12         =  of0/u(2)*( 1.0_r8+OMEGA1*u(2)*(-u(3))+OMEGA2*p2u2*p2u2+OMEGA3*p3u2*p3u2 )  
        bu13         =  of0/u(3)*( 1.0_r8+OMEGA1*u(3)*(-u(3))+OMEGA2*p2u3*p2u2+OMEGA3*p3u3*p3u2 )  
        bu14         =  of0/u(4)*( 1.0_r8+OMEGA1*u(4)*(-u(3))+OMEGA2*p2u4*p2u2+OMEGA3*p3u4*p3u2 ) 

        bu21         =  of0/u(1)*( 1.0_r8+OMEGA1*u(1)*(-u(4))+OMEGA2*p2u1*p2u1+OMEGA3*p3u1*p3u1 )
        bu22         =  of0/u(2)*( 1.0_r8+OMEGA1*u(2)*(-u(4))+OMEGA2*p2u2*p2u1+OMEGA3*p3u2*p3u1 )
        bu23         =  of0/u(3)*( 1.0_r8+OMEGA1*u(3)*(-u(4))+OMEGA2*p2u3*p2u1+OMEGA3*p3u3*p3u1 )
        bu24         =  of0/u(4)*( 1.0_r8+OMEGA1*u(4)*(-u(4))+OMEGA2*p2u4*p2u1+OMEGA3*p3u4*p3u1 )

        bu01         =  of0/u(1)*( 1.0_r8+OMEGA1*u(1)*(-AMU0)+OMEGA2*p2u1*p2f0+OMEGA3*p3u1*p3f0 )
        bu02         =  of0/u(2)*( 1.0_r8+OMEGA1*u(2)*(-AMU0)+OMEGA2*p2u2*p2f0+OMEGA3*p3u2*p3f0 )
        bu03         =  of0/u(3)*( 1.0_r8+OMEGA1*u(3)*(-AMU0)+OMEGA2*p2u3*p2f0+OMEGA3*p3u3*p3f0 )
        bu04         =  of0/u(4)*( 1.0_r8+OMEGA1*u(4)*(-AMU0)+OMEGA2*p2u4*p2f0+OMEGA3*p3u4*p3f0 )
!--------------------------- c(4,4) -----------------------------------
        oa1          =  OM/2.0_r8*a(1)   ! oa1=oa4
        oa2          =  OM/2.0_r8*a(2)   ! oa2=oa3

        c11          =  (oa1*(1.0_r8+OMEGA1*u(1)**2+OMEGA2*p2u1**2+OMEGA3*p3u1**2)-1.0)/u(1)
        c22          =  (oa2*(1.0_r8+OMEGA1*u(2)**2+OMEGA2*p2u2**2+OMEGA3*p3u2**2)-1.0)/u(2)
        c33          =  (oa2*(1.0_r8+OMEGA1*u(3)**2+OMEGA2*p2u3**2+OMEGA3*p3u3**2)-1.0)/u(3)
        c44          =  (oa1*(1.0_r8+OMEGA1*u(4)**2+OMEGA2*p2u4**2+OMEGA3*p3u4**2)-1.0)/u(4)

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
	   rk1          =  (bsi/2.0_r8+(bsi**2+4_r8*csi)**0.5/2.0_r8)**0.5
       rk2          =  bsi/2.0-(bsi**2+4_r8*csi)**0.5/2.0_r8
	   rk2          =  rk2**0.5
       AA1          =  (rk1**2-a22)/a21
	   AA2          =  (rk2**2-a22)/a21
!-------------------------------- phi ----------------------------------
       p1p          =  0.5_r8*(AA1+(AA1*b22n-b12n)/an*rk1)
	   p1n          =  0.5_r8*(AA1-(AA1*b22n-b12n)/an*rk1)
	   p2p          =  0.5_r8*(AA2+(AA2*b22n-b12n)/an*rk2)
	   p2n          =  0.5_r8*(AA2-(AA2*b22n-b12n)/an*rk2)
!---------------------- varphi ---- psi --------------------------------
       v1p          =  0.5_r8*(1._r8+(b11n-AA1*b21n)/an*rk1)
       v1n          =  0.5_r8*(1._r8-(b11n-AA1*b21n)/an*rk1)
	   v2p          =  0.5_r8*(1._r8+(b11n-AA2*b21n)/an*rk2)
       v2n          =  0.5_r8*(1._r8-(b11n-AA2*b21n)/an*rk2)
       psi1         =  exp(-rk1*tau)
       psi2         =  exp(-rk2*tau)
!---------------dependent on u1,u2,u0---------------
       b2u1p        =  bu14+bu11
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

!----------d(:,1,:) normal; d(:,2,:)=d` acute;d(:,:,1) u1;d(:,:,3) u2; d(:,:,3) u0;------------------------c
!             fa=f`
        d13        =  b22n*b2u1n+b21n*b1u1n+b2u1p/u(3)
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

	    fa1        =  RMU1**4-bsi*RMU1**2-csi 
	    fa2        =  RMU2**4-bsi*RMU2**2-csi
	    fa3        =  RMU0**4-bsi*RMU0**2-csi

!-----------------eta(1,:,:) eta1;eta(:,1,:) normal; eta(:,2,:) acute;d(:,:,1) u1;d(:,:,3) u2; d(:,:,3) u0;
	    eta11        =  (d11*RMU1**2+a12*d13-a22*d11)/fa1
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

! ...  DIRECT SOLUTION

!---------------------- H1, H2 ----------------------------------- 
        H101          =  -0.5_r8*(eta33+eta34)*eft0
	    H102          =  -0.5_r8*(eta31+eta32)*eft0
        H103          =  -0.5_r8*(eta31-eta32)
	    H104          =  -0.5_r8*(eta33-eta34)
 
        H201          =  0.5_r8*(eta33+eta34)
	    H202          =  0.5_r8*(eta31+eta32)
        H203          =  0.5_r8*(eta31-eta32)*eft0
	    H204          =  0.5_r8*(eta33-eta34)*eft0
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
!------------------------ G1 -----------------------------------------
        C1           =  YA * H101 + YB * H102 + YC * H103 + YD * H104
        D1           =  YE * H101 + YF * H102 + YG * H103 + YH * H104
        C2           =  YH * H101 + YG * H102 + YF * H103 + YE * H104
        D2           =  YD * H101 + YC * H102 + YB * H103 + YA * H104
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
        RU0(jk,1)     = (V21*C1+V22*D1+V23*C2+V24*D2+H202) * RMU0
        RU0(jk,2)     = (V11*C1+V12*D1+V13*C2+V14*D2+H201) * RMU0
        TU0(jk,1)     = (V24*C1+V23*D1+V22*C2+V21*D2+H203) * RMU0
        TU0(jk,2)     = (V14*C1+V13*D1+V12*C2+V11*D2+H204) * RMU0
        TDIR(jk)      =  eft0
! ...   DIFFUSE SOLUTION

!----------------------- H1, H2 -------------------------------------
       H111          =  -0.5_r8*(eta13+eta14)*eft1
	   H112          =  -0.5_r8*(eta11+eta12)*eft1
       H113          =  -0.5_r8*(eta11-eta12)
	   H114          =  -0.5_r8*(eta13-eta14)

       H211          =  0.5_r8*(eta13+eta14)
	   H212          =  0.5_r8*(eta11+eta12)
       H213          =  0.5_r8*(eta11-eta12)*eft1
	   H214          =  0.5_r8*(eta13-eta14)*eft1
!------------------------ G1 -----------------------------------------
        C1           =  YA * H111 + YB * H112 + YC * H113 + YD * H114
        D1           =  YE * H111 + YF * H112 + YG * H113 + YH * H114
        C2           =  YH * H111 + YG * H112 + YF * H113 + YE * H114
        D2           =  YD * H111 + YC * H112 + YB * H113 + YA * H114

        RBAR(jk,1,1)  =  (V21*C1+V22*D1+V23*C2+V24*D2+H212)*2._r8*a(3)
        RBAR(jk,2,1)  =  (V11*C1+V12*D1+V13*C2+V14*D2+H211)*2._r8*a(4) 
        TBAR(jk,1,1)  =  (V24*C1+V23*D1+V22*C2+V21*D2+H213)*2._r8*a(3)+eft1
        TBAR(jk,2,1)  =  (V14*C1+V13*D1+V12*C2+V11*D2+H214)*2._r8*a(4)
!------------------------ H1, H2 -------------------------------------
       H121          =  -0.5_r8*(eta23+eta24)*eft2
	   H122          =  -0.5_r8*(eta21+eta22)*eft2
       H123          =  -0.5_r8*(eta21-eta22)
	   H124          =  -0.5_r8*(eta23-eta24)

       H221          =  0.5_r8*(eta23+eta24)
	   H222          =  0.5_r8*(eta21+eta22)
       H223          =  0.5_r8*(eta21-eta22)*eft2
	   H224          =  0.5_r8*(eta23-eta24)*eft2
!------------------------ G1 -----------------------------------------
        C1           =  YA * H121 + YB * H122 + YC * H123 + YD * H124
        D1           =  YE * H121 + YF * H122 + YG * H123 + YH * H124
        C2           =  YH * H121 + YG * H122 + YF * H123 + YE * H124
        D2           =  YD * H121 + YC * H122 + YB * H123 + YA * H124

        RBAR(jk,1,2)  =  (V21*C1+V22*D1+V23*C2+V24*D2+H222)*2*a(3)
        RBAR(jk,2,2)  =  (V11*C1+V12*D1+V13*C2+V14*D2+H221)*2*a(4) 
        TBAR(jk,1,2)  =  (V24*C1+V23*D1+V22*C2+V21*D2+H223)*2*a(3)
        TBAR(jk,2,2)  =  (V14*C1+V13*D1+V12*C2+V11*D2+H224)*2*a(4)+eft2

!open(10,file='FLUD',form='formatted')
!write(10,*) "Ru0Tu0",jk,RU0(jk,1),RU0(jk,2),TU0(jk,1),TU0(jk,2)

	 end if
   end do

      end subroutine reftra_4DDA_sw

      end module rrtmg_sw_reftra

