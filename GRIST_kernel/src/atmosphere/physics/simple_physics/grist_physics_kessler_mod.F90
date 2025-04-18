!
! Code impoted from CESM2.0
!
module grist_physics_kessler_mod

  use grist_constants, only: r8, rdry, cp, p00

  implicit none
  private
  save

  public :: kessler_cam ! Main routine
!  public :: kessler_set_const ! Initialize physical constants

  ! Private module data (constants set at initialization)
  real(r8) :: rd    = rdry   ! gas constant for dry air, J/(kgK)
!  real(r8) :: cp    = cp   ! heat capacity at constant pressure, J/(kgK)
  real(r8) :: lv    = 2.501e6_R8    ! latent heat of vaporization, J/kg
  real(r8) :: psl   = 1000._r8  ! reference pressure at sea level, mb
  real(r8) :: rhoqr = 1e3_r8! density of liquid water, kg/m^3

CONTAINS

!  subroutine kessler_set_const(rd_in, cp_in, lv_in, psl_in, rhoqr_in)
    ! Set physical constants to be consistent with calling model
!    real(r8), intent(in) :: rd_in    ! gas constant for dry air, J/(kgK)
!    real(r8), intent(in) :: cp_in    ! heat capacity at constant pres., J/(kgK)
!    real(r8), intent(in) :: lv_in    ! latent heat of vaporization, J/kg
!    real(r8), intent(in) :: psl_in   ! reference pressure at sea level, mb
!    real(r8), intent(in) :: rhoqr_in ! density of liquid water, kg/m^3

!    rd    = rd_in
!    cp    = cp_in
!    lv    = lv_in
!    psl   = psl_in
!    rhoqr = rhoqr_in
    
!  end subroutine kessler_set_const

  !-----------------------------------------------------------------------
  !
  !  Version:  2.0
  !
  !  Date:  January 22nd, 2015
  !
  !  Change log:
  !  v2 - Added sub-cycling of rain sedimentation so as not to violate
  !       CFL condition.
  !
  !  The KESSLER subroutine implements the Kessler (1969) microphysics
  !  parameterization as described by Soong and Ogura (1973) and Klemp
  !  and Wilhelmson (1978, KW). KESSLER is called at the end of each
  !  time step and makes the final adjustments to the potential
  !  temperature and moisture variables due to microphysical processes
  !  occurring during that time step. KESSLER is called once for each
  !  vertical column of grid cells. Increments are computed and added
  !  into the respective variables. The Kessler scheme contains three
  !  moisture categories: water vapor, cloud water (liquid water that
  !  moves with the flow), and rain water (liquid water that falls
  !  relative to the surrounding air). There  are no ice categories.
  !  Variables in the column are ordered from the surface to the top.
  !
  !  SUBROUTINE KESSLER(theta, qv, qc, qr, rho, pk, dt, z, nz, rainnc)
  !
  !  Input variables:
  !     temp   - temperature (K)
  !     qv     - water vapor mixing ratio (gm/gm)
  !     qc     - cloud water mixing ratio (gm/gm)
  !     qr     - rain  water mixing ratio (gm/gm)
  !     rho    - dry air density (not mean state as in KW) (kg/m^3)
  !     pk     - Exner function  (not mean state as in KW) (p/p0)**(R/cp)
  !     dt     - time step (s)
  !     z      - heights of thermodynamic levels in the grid column (m)
  !     nz     - number of thermodynamic levels in the column
  !     precl  - Precipitation rate (m_water/s)
  !
  ! Output variables:
  !     Increments are added into t, qv, qc, qr, and rainnc which are
  !     returned to the routine from which KESSLER was called. To obtain
  !     the total precip qt, after calling the KESSLER routine, compute:
  !
  !       qt = sum over surface grid cells of (rainnc * cell area)  (kg)
  !       [here, the conversion to kg uses (10^3 kg/m^3)*(10^-3 m/mm) = 1]
  !
  !
  !  Authors: Paul Ullrich
  !           University of California, Davis
  !           Email: paullrich@ucdavis.edu
  !
  !           Based on a code by Joseph Klemp
  !           (National Center for Atmospheric Research)
  !
  !  Reference:
  !
  !    Klemp, J. B., W. C. Skamarock, W. C., and S.-H. Park, 2015:
  !    Idealized Global Nonhydrostatic Atmospheric Test Cases on a Reduced
  !    Radius Sphere. Journal of Advances in Modeling Earth Systems. 
  !    doi:10.1002/2015MS000435
  !
  !=======================================================================

  SUBROUTINE KESSLER_CAM(nz, dt, rho, z, pk, theta, qv, qc, qr, precl,&
                         dtheta_dt,dqv_dt,dqc_dt,dqr_dt)
    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------

!    integer,          intent(in)    :: ncol      ! Number of columns
    integer,          intent(in)    :: nz        ! Number of vertical levels
    real(r8),         intent(in)    :: dt        ! Time step (s)
    real(r8),         intent(in)    :: rho(:)  ! Dry air density (kg/m^3)
    real(r8),         intent(in)    :: z(:)    ! Heights of thermo. levels (m)
    real(r8),         intent(in)    :: pk(:)   ! Exner function (p/p0)**(R/cp)

    real(r8),         intent(inout) :: theta(:) ! temperature (K)
    real(r8),         intent(inout) :: qv(:)   ! Water vapor mixing ratio (gm/gm)
    real(r8),         intent(inout) :: qc(:)   ! Cloud water mixing ratio (gm/gm)
    real(r8),         intent(inout) :: qr(:)   ! Rain  water mixing ratio (gm/gm)

    real(r8),         intent(out)   :: precl  ! Precipitation rate (m_water / s)

    real(r8),         intent(out) :: dtheta_dt(:) ! temperature (K)
    real(r8),         intent(out) :: dqv_dt(:)   ! Water vapor mixing ratio (gm/gm)
    real(r8),         intent(out) :: dqc_dt(:)   ! Cloud water mixing ratio (gm/gm)
    real(r8),         intent(out) :: dqr_dt(:)   ! Rain  water mixing ratio (gm/gm)
    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    real(r8) :: r(nz), rhalf(nz), velqr(nz), sed(nz), pc(nz)
    real(r8) :: f5, f2x, xk, ern, qrprod, prod, qvs, dt_max, dt0
    integer  ::  k, rainsplit, nt
    real(r8) :: theta_ini(nz), qv_ini(nz), qc_ini(nz), qr_ini(nz)

    theta_ini = theta
    qv_ini    = qv
    qc_ini    = qc
    qr_ini    = qr

    ! Initialize output variables
    precl = 0._r8

    ! Check inputs
    if (dt <= 0._r8) then
      write(6,*) 'KESSLER called with nonpositive dt'
      return
    end if

    !------------------------------------------------
    !   Begin calculation
    !------------------------------------------------
    f2x   = 17.27_r8
    f5    = 237.3_r8 * f2x * lv / cp
    xk    = .2875_r8  !  kappa (r/cp)
    ! Loop through columns
!    do col = 1, ncol
      do k = 1, nz
        r(k)     = 0.001_r8 * rho( k)
        rhalf(k) = sqrt(rho( 1) / rho( k))
        pc(k)    = 3.8_r8 / (pk( k)**(1._r8/xk)*psl)
        !
        ! if qr is (round-off) negative then the computation of
        ! velqr triggers floating point exception error when running
        ! in debugging mode with NAG
        !
        qr(k) = MAX(qr(k),0.0_r8)
        !qr(k) = MAX(qr(k),1e-80.0_r8)
        !
        ! Liquid water terminal velocity (m/s) following KW eq. 2.15
        velqr(k)  = 36.34_r8 * rhalf(k) * (qr( k) * r(k))**0.1364_r8
      end do

      ! Maximum time step size in accordance with CFL condition
      dt_max = dt
      do k = 1, nz - 1
!        if (velqr(k) /= 0._r8) then !this causes rainsplit to become NaN on Hobart
        if (abs(velqr(k)) > 1.0E-12_r8) then
          dt_max = min(dt_max, 0.8_r8*(z( k+1) - z( k)) / velqr(k))
        end if
      end do

      ! Number of subcycles
      rainsplit = ceiling(dt / dt_max)
      if (rainsplit < 1) then
        write(6, *) 'KESSLER: bad rainsplit ',dt,dt_max,rainsplit
        return
      end if
      dt0 = dt / real(rainsplit, r8)

      ! Subcycle through rain process
      nt = 1
      do while (nt.le.rainsplit)

        ! Precipitation rate (m/s)
        precl = precl + rho( 1) * qr( 1) * velqr(1) / rhoqr

        ! Sedimentation term using upstream differencing
        do k = 1, nz-1
          sed(k) = dt0 * ((r(k+1) * qr( k+1) * velqr(k+1)) -              &
                          (r(k) *   qr( k)   * velqr(k))) /               &
                          (r(k) * (z( k+1) - z( k)))
        end do
        sed(nz) = -dt0 * qr( nz) * velqr(nz) / (0.5_r8 * (z( nz)-z( nz-1)))

        ! Adjustment terms
        do k = 1, nz

          ! Autoconversion and accretion rates following KW eq. 2.13a,b
          qrprod = qc( k) - (qc( k) - dt0 * max(.001_r8 * (qc( k)-.001_r8), 0._r8)) / &
               (1._r8 + dt0 * 2.2_r8 * qr( k)**.875_r8)
          qc( k) = max(qc( k) - qrprod, 0._r8)
          qr( k) = max(qr( k) + qrprod + sed(k), 0._r8)

          ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
          qvs = pc(k) * exp(f2x*(pk( k)*theta( k) - 273._r8) / (pk( k)*theta( k) - 36._r8))
          prod = (qv( k) - qvs) / (1._r8 + qvs*f5 / (pk( k)*theta( k) - 36._r8)**2)


          ! Evaporation rate following KW eq. 2.14a,b
          ern = min(dt0 * (((1.6_r8 + 124.9_r8*(r(k)*qr( k))**.2046_r8) * &
               (r(k) * qr( k))**.525_r8) /                                &
               (2550000._r8 * pc(k) / (3.8_r8*qvs) + 540000._r8)) *           &
               (dim(qvs,qv( k)) / (r(k)*qvs)),                            &
               max(-prod-qc( k),0._r8),qr( k))

          ! Saturation adjustment following KW eq. 3.10
          theta( k)= theta( k) + (lv / (cp * pk( k)) * (max(prod,-qc( k)) - ern))
          qv( k) = max(qv( k) - max(prod, -qc( k)) + ern, 0._r8)
          qc( k) = qc( k) + max(prod, -qc( k))
          qr( k) = qr( k) - ern
        end do

        ! Recalculate liquid water terminal velocity
        if (nt /= rainsplit) then
          do k = 1, nz
            velqr(k)  = 36.34_r8 * rhalf(k) * (qr( k)*r(k))**0.1364_r8
          end do
          !
          ! recompute rainsplit since velqr has changed
          !
          do k = 1, nz - 1
            if (abs(velqr(k)) > 1.0E-12_r8) then
              dt_max = min(dt_max, 0.8_r8*(z( k+1) - z( k)) / velqr(k))
            end if
          end do
          ! Number of subcycles
          rainsplit = ceiling(dt / dt_max)
          if (rainsplit < 1) then
            write(6, *) 'KESSLER: bad rainsplit ',dt,dt_max,rainsplit
            return
          end if
          dt0 = dt / real(rainsplit, r8)
        end if
        nt=nt+1
      end do

      precl = precl / real(rainsplit, r8)

!    end do ! column loop
     dtheta_dt = (theta-theta_ini)/dt
     dqv_dt = (qv-qv_ini)/dt
     dqc_dt = (qc-qc_ini)/dt
     dqr_dt = (qr-qr_ini)/dt
     return
  END SUBROUTINE KESSLER_CAM

!=======================================================================
end module grist_physics_kessler_mod
