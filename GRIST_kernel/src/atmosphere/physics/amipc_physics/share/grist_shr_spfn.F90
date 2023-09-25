!======================================================
!
!  Created by LiXiaohan on 19/7/5.
!  Shared ERF_INTRINSICS for GRIST physics package
!======================================================


 module grist_shr_spfn
    use grist_constants,                    only: i4, r8, pi

    implicit none
    private
    save

    public ::   shr_spfn_erf,              & 
                shr_spfn_erfc,             &
                shr_spfn_erfc_scaled,      &
                shr_spfn_gamma,            &
                shr_spfn_igamma

! Mathematical constants
real(r8), parameter :: sqrtpi = 1.77245385090551602729_r8

! Define machine-specific constants needed in this module
! These were used by the original gamma and calerf functions
! to guarantee safety against overflow, and precision, on
! many different machines.

! By defining the constants in this way, we assume
! that 1/xmin is representable (i.e. does not overflow
! the real type). This assumption was not in the original
! code, but is valid for IEEE single and double precision.

! Double precision
!---------------------------------------------------------------------
! Machine epsilon
real(r8), parameter :: epsr8 = epsilon(1._r8)
! "Huge" value is returned when actual value would be infinite.
real(r8), parameter :: xinfr8 = huge(1._r8)
! Smallest normal value.
real(r8), parameter :: xminr8 = tiny(1._r8)
! Largest number that, when added to 1., yields 1.
real(r8), parameter :: xsmallr8 = epsr8/2._r8
! Largest argument for which erfcx > 0.
real(r8), parameter :: xmaxr8 = 1/(sqrtpi*xminr8)

! For gamma/igamma
! Approximate value of largest acceptable argument to gamma,
! for IEEE double-precision.
real(r8), parameter :: xbig_gamma = 171.624_r8


contains

! Purpose: This subprogram computes approximate values for erf(x).
!          (see comments heading CALERF)
!          Author/date: W. J. Cody, January 8, 1985
    function shr_spfn_erf(x)
    real(r8), intent(in) :: x
    real(r8)             :: shr_spfn_erf
    integer(i4)          :: jint = 0

    call calerf_r8(x, shr_spfn_erf, jint)

    end function shr_spfn_erf


! Purpose: This subprogram computes approximate values for erfc(x).
!          (see comments heading CALERF)
!          Author/date: W. J. Cody, January 8, 1985
    function shr_spfn_erfc(x)
    real(r8), intent(in) :: x
    real(r8)             :: shr_spfn_erfc
    integer(i4)          :: jint = 1

    call calerf_r8(x, shr_spfn_erfc, jint)

    end function shr_spfn_erfc


! Purpose: This subprogram computes approximate values for exp(x*x)*erfc(x).
!          (see comments heading CALERF)
!          Author/date: W. J. Cody, January 8, 1985
    function shr_spfn_erfc_scaled(x)
    real(r8), intent(in) :: x
    real(r8)             :: shr_spfn_erfc_scaled
    integer(i4)          :: jint = 2

    call calerf_r8(x, shr_spfn_erfc_scaled, jint)

    end function shr_spfn_erfc_scaled


! Purpose: Gamma functions
!----------------------------------------------------------------------
!
! THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
!   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
!   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
!   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
!   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
!   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
!   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
!   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
!   MACHINE-DEPENDENT CONSTANTS.
!
!
! EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
!
! BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
! MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
! XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
!          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
!                  GAMMA(XBIG) = BETA**MAXEXP
! XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
!          APPROXIMATELY BETA**MAXEXP
! EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1.0+EPS .GT. 1.0
! XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
!          1/XMININ IS MACHINE REPRESENTABLE
!
!     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
!
!                            BETA       MAXEXP        XBIG
!
! CRAY-1         (S.P.)        2         8191        966.961
! CYBER 180/855
!   UNDER NOS    (S.P.)        2         1070        177.803
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)        2          128        35.040
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)        2         1024        171.624
! IBM 3033       (D.P.)       16           63        57.574
! VAX D-FORMAT   (D.P.)        2          127        34.844
! VAX G-FORMAT   (D.P.)        2         1023        171.489
!
!                            XINF         EPS        XMININ
!
! CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
! CYBER 180/855
!   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
! IEEE (IBM/XT,
!   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
! IEEE (IBM/XT,
!   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
! IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
! VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
! VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
!
!*******************************************************************
!*******************************************************************
!
! ERROR RETURNS
!
!  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
!     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
!     TO BE FREE OF UNDERFLOW AND OVERFLOW.
!
!
!  INTRINSIC FUNCTIONS REQUIRED ARE:
!
!     INT, DBLE, EXP, LOG, REAL, SIN
!
!
! REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
!              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
!              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
!              (ED.), SPRINGER VERLAG, BERLIN, 1976.
!
!              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
!              SONS, NEW YORK, 1968.
!
!  LATEST MODIFICATION: OCTOBER 12, 1989
!
!  AUTHORS: W. J. CODY AND L. STOLTZ
!           APPLIED MATHEMATICS DIVISION
!           ARGONNE NATIONAL LABORATORY
!           ARGONNE, IL 60439
!
    function shr_spfn_gamma(x) result(gamma)
    real(r8), intent(in) :: x
    real(r8)             :: gamma
    real(r8) :: fact, res, sum, xden, xnum, y, y1, ysq, z

    integer :: i, n
    logical :: negative_odd

    ! log(2*pi)/2
    real(r8), parameter :: logsqrt2pi = 0.9189385332046727417803297E0_r8

!----------------------------------------------------------------------
!  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
!     APPROXIMATION OVER (1,2).
!----------------------------------------------------------------------
    real(r8), parameter :: P(8) = &
       (/-1.71618513886549492533811E+0_r8, 2.47656508055759199108314E+1_r8, &
         -3.79804256470945635097577E+2_r8, 6.29331155312818442661052E+2_r8, &
          8.66966202790413211295064E+2_r8,-3.14512729688483675254357E+4_r8, &
         -3.61444134186911729807069E+4_r8, 6.64561438202405440627855E+4_r8 /)
    real(r8), parameter :: Q(8) = &
       (/-3.08402300119738975254353E+1_r8, 3.15350626979604161529144E+2_r8, &
         -1.01515636749021914166146E+3_r8,-3.10777167157231109440444E+3_r8, &
          2.25381184209801510330112E+4_r8, 4.75584627752788110767815E+3_r8, &
         -1.34659959864969306392456E+5_r8,-1.15132259675553483497211E+5_r8 /)
!----------------------------------------------------------------------
!  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
!----------------------------------------------------------------------
    real(r8), parameter :: C(7) = &
       (/-1.910444077728E-03_r8,          8.4171387781295E-04_r8, &
         -5.952379913043012E-04_r8,       7.93650793500350248E-04_r8, &
         -2.777777777777681622553E-03_r8, 8.333333333333333331554247E-02_r8, &
          5.7083835261E-03_r8 /)

  negative_odd = .false.
  fact = 1._r8
  n = 0
  y = x
  if (y <= 0._r8) then
!----------------------------------------------------------------------
!  ARGUMENT IS NEGATIVE
!----------------------------------------------------------------------
     y = -x
     y1 = aint(y)
     res = y - y1
     if (res /= 0._r8) then
        negative_odd = (y1 /= aint(y1*0.5_r8)*2._r8)
        fact = -pi/sin(pi*res)
        y = y + 1._r8
     else
        gamma = xinfr8
        return
     end if
  end if
!----------------------------------------------------------------------
!  ARGUMENT IS POSITIVE
!----------------------------------------------------------------------
  if (y < epsr8) then
!----------------------------------------------------------------------
!  ARGUMENT .LT. EPS
!----------------------------------------------------------------------
     if (y >= xminr8) then
        res = 1._r8/y
     else
        gamma = xinfr8
        return
     end if
  elseif (y < 12._r8) then
     y1 = y
     if (y < 1._r8) then
!----------------------------------------------------------------------
!  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
        z = y
        y = y + 1._r8
     else
!----------------------------------------------------------------------
!  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
!----------------------------------------------------------------------
        n = int(y) - 1
        y = y - real(n, r8)
        z = y - 1._r8
     end if
!----------------------------------------------------------------------
!  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
!----------------------------------------------------------------------
     xnum = 0._r8
     xden = 1._r8
     do i=1,8
        xnum = (xnum+P(i))*z
        xden = xden*z + Q(i)
     end do
     res = xnum/xden + 1._r8
     if (y1 < y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
!----------------------------------------------------------------------
        res = res/y1
     elseif (y1 > y) then
!----------------------------------------------------------------------
!  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
!----------------------------------------------------------------------
        do i = 1,n
           res = res*y
           y = y + 1._r8
        end do
     end if
  else
!----------------------------------------------------------------------
!  EVALUATE FOR ARGUMENT .GE. 12.0,
!----------------------------------------------------------------------
     if (y <= xbig_gamma) then
        ysq = y*y
        sum = C(7)
        do i=1,6
           sum = sum/ysq + C(i)
        end do
        sum = sum/y - y + logsqrt2pi
        sum = sum + (y-0.5_r8)*log(y)
        res = exp(sum)
     else
        gamma = xinfr8
        return
     end if
  end if
!----------------------------------------------------------------------
!  FINAL ADJUSTMENTS AND RETURN
!----------------------------------------------------------------------
  if (negative_odd)  res = -res
  if (fact /= 1._r8) res = fact/res
  gamma = res
! ---------- LAST LINE OF GAMMA ----------
end function shr_spfn_gamma


!Purpose: Incomplete Gamma function
!         Upper incomplete gamma function.
!         Modified for inclusion in this module and made pure elemental
real(r8) function shr_spfn_igamma(a, x)

  real(r8), intent(in) ::      a
  real(r8), intent(in) ::      x

  ! local variable
  real(r8) :: xam, gin, s, r, t0
  integer  :: k


  if (x == 0.0_r8) then
     shr_spfn_igamma = shr_spfn_gamma(a)
     return
  end if

  xam = -x + a * log(x)
  
  if ((xam > 700.0_r8) .or. (a > xbig_gamma)) then
     ! Out of bounds
     ! Return "huge" value.
     shr_spfn_igamma = xinfr8
     return

  else if (x <= (1.0_r8 + a)) then
     s = 1.0_r8 / a
     r = s

     do  k = 1,60
        r = r * x / (a+k)
        s = s + r

        if (abs(r/s) < 1.0e-15_r8) exit
     end do
        
     gin = exp(xam) * s           
     shr_spfn_igamma = shr_spfn_gamma(a) - gin
        
  else
     t0 = 0.0_r8

     do k = 60,1,-1
        t0 = (k - a) / (1.0_r8 + k / (x + t0))
     end do
     shr_spfn_igamma = exp(xam) / (x + t0)
  endif

end function shr_spfn_igamma



SUBROUTINE CALERF_r8(ARG, RESULT, JINT)

   !------------------------------------------------------------------
   !  This version uses 8-byte reals
   !------------------------------------------------------------------

   ! arguments
   real(r8),     intent(in)  :: arg
   integer(i4),  intent(in)  :: jint
   real(r8),     intent(out) :: result

   ! local variables
   INTEGER(i4) :: I

   real(r8) :: X, Y, YSQ, XNUM, XDEN, DEL

   !------------------------------------------------------------------
   !  Mathematical constants
   !------------------------------------------------------------------
   real(r8), parameter :: ZERO   = 0.0E0_r8
   real(r8), parameter :: FOUR   = 4.0E0_r8
   real(r8), parameter :: ONE    = 1.0E0_r8
   real(r8), parameter :: HALF   = 0.5E0_r8
   real(r8), parameter :: TWO    = 2.0E0_r8
   ! 1/sqrt(pi)
   real(r8), parameter :: SQRPI  = 5.6418958354775628695E-1_r8
   real(r8), parameter :: THRESH = 0.46875E0_r8
   real(r8), parameter :: SIXTEN = 16.0E0_r8

   !------------------------------------------------------------------
   !  Machine-dependent constants: IEEE double precision values
   !------------------------------------------------------------------
   real(r8), parameter :: XNEG   = -26.628E0_r8
   real(r8), parameter :: XBIG   =  26.543E0_r8
   real(r8), parameter :: XHUGE  =   6.71E7_r8

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erf  in first interval
   !------------------------------------------------------------------
   real(r8), parameter :: A(5) = (/ 3.16112374387056560E00_r8, 1.13864154151050156E02_r8, &
                                    3.77485237685302021E02_r8, 3.20937758913846947E03_r8, &
                                    1.85777706184603153E-1_r8 /)
   real(r8), parameter :: B(4) = (/ 2.36012909523441209E01_r8, 2.44024637934444173E02_r8, &
                                    1.28261652607737228E03_r8, 2.84423683343917062E03_r8 /)

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erfc  in second interval
   !------------------------------------------------------------------
   real(r8), parameter :: C(9) = (/ 5.64188496988670089E-1_r8, 8.88314979438837594E00_r8, &
                                    6.61191906371416295E01_r8, 2.98635138197400131E02_r8, &
                                    8.81952221241769090E02_r8, 1.71204761263407058E03_r8, &
                                    2.05107837782607147E03_r8, 1.23033935479799725E03_r8, &
                                    2.15311535474403846E-8_r8 /)
   real(r8), parameter :: D(8) = (/ 1.57449261107098347E01_r8, 1.17693950891312499E02_r8, &
                                    5.37181101862009858E02_r8, 1.62138957456669019E03_r8, &
                                    3.29079923573345963E03_r8, 4.36261909014324716E03_r8, &
                                    3.43936767414372164E03_r8, 1.23033935480374942E03_r8 /)

   !------------------------------------------------------------------
   !  Coefficients for approximation to  erfc  in third interval
   !------------------------------------------------------------------
   real(r8), parameter :: P(6) = (/ 3.05326634961232344E-1_r8, 3.60344899949804439E-1_r8, &
                                    1.25781726111229246E-1_r8, 1.60837851487422766E-2_r8, &
                                    6.58749161529837803E-4_r8, 1.63153871373020978E-2_r8 /)
   real(r8), parameter :: Q(5) = (/ 2.56852019228982242E00_r8, 1.87295284992346047E00_r8, &
                                    5.27905102951428412E-1_r8, 6.05183413124413191E-2_r8, &
                                    2.33520497626869185E-3_r8 /)

   !------------------------------------------------------------------
   X = ARG
   Y = ABS(X)
   IF (Y .LE. THRESH) THEN
      !------------------------------------------------------------------
      !  Evaluate  erf  for  |X| <= 0.46875
      !------------------------------------------------------------------
      YSQ = ZERO
      IF (Y .GT. XSMALLR8) YSQ = Y * Y
      XNUM = A(5)*YSQ
      XDEN = YSQ
      DO I = 1, 3
         XNUM = (XNUM + A(I)) * YSQ
         XDEN = (XDEN + B(I)) * YSQ
      end do
      RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
      IF (JINT .NE. 0) RESULT = ONE - RESULT
      IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
      GO TO 80
   ELSE IF (Y .LE. FOUR) THEN
      !------------------------------------------------------------------
      !  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
      !------------------------------------------------------------------
      XNUM = C(9)*Y
      XDEN = Y
      DO I = 1, 7
         XNUM = (XNUM + C(I)) * Y
         XDEN = (XDEN + D(I)) * Y
      end do
      RESULT = (XNUM + C(8)) / (XDEN + D(8))
      IF (JINT .NE. 2) THEN
         YSQ = AINT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      END IF
   ELSE
      !------------------------------------------------------------------
      !  Evaluate  erfc  for |X| > 4.0
      !------------------------------------------------------------------
      RESULT = ZERO
      IF (Y .GE. XBIG) THEN
         IF ((JINT .NE. 2) .OR. (Y .GE. XMAXR8)) GO TO 30
         IF (Y .GE. XHUGE) THEN
            RESULT = SQRPI / Y
            GO TO 30
         END IF
      END IF
      YSQ = ONE / (Y * Y)
      XNUM = P(6)*YSQ
      XDEN = YSQ
      DO I = 1, 4
         XNUM = (XNUM + P(I)) * YSQ
         XDEN = (XDEN + Q(I)) * YSQ
      end do
      RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
      RESULT = (SQRPI -  RESULT) / Y
      IF (JINT .NE. 2) THEN
         YSQ = AINT(Y*SIXTEN)/SIXTEN
         DEL = (Y-YSQ)*(Y+YSQ)
         RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
      END IF
   END IF
30 continue
   !------------------------------------------------------------------
   !  Fix up for negative argument, erf, etc.
   !------------------------------------------------------------------
   IF (JINT .EQ. 0) THEN
      RESULT = (HALF - RESULT) + HALF
      IF (X .LT. ZERO) RESULT = -RESULT
   ELSE IF (JINT .EQ. 1) THEN
      IF (X .LT. ZERO) RESULT = TWO - RESULT
   ELSE
      IF (X .LT. ZERO) THEN
         IF (X .LT. XNEG) THEN
            RESULT = XINFR8
         ELSE
            YSQ = AINT(X*SIXTEN)/SIXTEN
            DEL = (X-YSQ)*(X+YSQ)
            Y = EXP(YSQ*YSQ) * EXP(DEL)
            RESULT = (Y+Y) - RESULT
         END IF
      END IF
   END IF
80 continue
end SUBROUTINE CALERF_r8


 end module grist_shr_spfn
