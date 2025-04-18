
MODULE grist_dtp_dcmip2016_sc

!=======================================================================
!
!  Date:  April 22, 2016
!
!  Functions for setting up idealized initial conditions for the
!  Klemp et al. supercell test.  Before sampling the result,
!  supercell_test_init() must be called.
!
!  SUBROUTINE supercell_test(
!    lon,lat,p,z,zcoords,u,v,t,thetav,ps,rho,q,pert)
!
!  Given a point specified by: 
!      lon    longitude (radians) 
!      lat    latitude (radians) 
!      p/z    pressure (Pa) / height (m)
!  zcoords    1 if z is specified, 0 if p is specified
!     pert    1 if thermal perturbation included, 0 if not
!
!  the functions will return:
!        p    pressure if z is specified (Pa)
!        z    geopotential height if p is specified (m)
!        u    zonal wind (m s^-1)
!        v    meridional wind (m s^-1)
!        t    temperature (K)
!   thetav    virtual potential temperature (K)
!       ps    surface pressure (Pa)
!      rho    density (kj m^-3)
!        q    water vapor mixing ratio (kg/kg)
!
!  Author: Paul Ullrich
!          University of California, Davis
!          Email: paullrich@ucdavis.edu
!
!          Based on a code by Joseph Klemp
!          (National Center for Atmospheric Research)
!  
!  Yi Zhang modified to include dry mass as an option
!
!=======================================================================

   use grist_constants,  only: i4, r8, zero, gravity, rdry, one

  IMPLICIT NONE
  private

  public :: supercell_init, supercell_test, supercell_z
!=======================================================================
!    Physical constants
!=======================================================================

  real(r8), PARAMETER ::               &
       a     = 6371220.0d0,           & ! Reference Earth's Radius (m)
       Rd    = 287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
       g     = 9.80616d0,             & ! Gravity (m s^2)
       cp    = 1004.5d0,              & ! Specific heat capacity (J kg^-1 K^1)
       Lvap  = 2.5d6,                 & ! Latent heat of vaporization of water
       Rvap  = 461.5d0,               & ! Ideal gas constnat for water vapor
       Mvap  = 0.608d0,               & ! Ratio of molar mass of dry air/water
       pi    = 3.14159265358979d0,    & ! pi
       p0    = 100000.0d0,            & ! surface pressure (Pa)
       kappa = 2.d0/7.d0,             & ! Ratio of Rd to cp
       omega = 7.29212d-5,            & ! Reference rotation rate of the Earth (s^-1)
       deg2rad  = pi/180.d0             ! Conversion factor of degrees to radians

!=======================================================================
!    Test case parameters
!=======================================================================
  INTEGER(4), PARAMETER ::            &
       nz         = 30         ,      & ! number of vertical levels in init
       nphi       = 16                  ! number of meridional points in init

  real(r8), PARAMETER ::               &
       z1         = 0.0d0      ,      & ! lower sample altitude
       z2         = 50000.0d0           ! upper sample altitude

  real(r8), PARAMETER ::               &
       X          = 120.d0     ,      & ! Earth reduction factor
       theta0     = 300.d0     ,      & ! theta at the equatorial surface
       theta_tr   = 343.d0     ,      & ! theta at the tropopause
       z_tr       = 12000.d0   ,      & ! altitude at the tropopause
       T_tr       = 213.d0     ,      & ! temperature at the tropopause
       pseq       = 100000.0d0          ! surface pressure at equator (Pa)
       !pseq       = 95690.0d0           ! surface pressure at equator (Pa)

  real(r8), PARAMETER ::               &
       us         = 30.d0      ,      & ! maximum zonal wind velocity
       uc         = 15.d0      ,      & ! coordinate reference velocity
       zs         = 5000.d0    ,      & ! lower altitude of maximum velocity
       zt         = 1000.d0             ! transition distance of velocity
 
  real(r8), PARAMETER ::               &
       pert_dtheta = 3.d0         ,   & ! perturbation magnitude
       pert_lonc   = 0.d0         ,   & ! perturbation longitude
       pert_latc   = 0.d0         ,   & ! perturbation latitude
       pert_rh     = 10000.d0 * X ,   & ! perturbation horiz. halfwidth
       pert_zc     = 1500.d0      ,   & ! perturbation center altitude
       pert_rz     = 1500.d0            ! perturbation vert. halfwidth

!-----------------------------------------------------------------------
!    Coefficients computed from initialization
!-----------------------------------------------------------------------

  INTEGER(4)                  :: initialized = 0

  real(r8), DIMENSION(nphi)    :: phicoord
  real(r8), DIMENSION(nz)      :: zcoord
  real(r8), DIMENSION(nphi,nz) :: thetavyz
  real(r8), DIMENSION(nphi,nz) :: exneryz
  real(r8), DIMENSION(nz)      :: qveq

!
! Weights and arguments for Gaussian Qudrature
!
#ifndef GAUSS_30
  integer , parameter  :: nGauss = 20
  real(r8), parameter, dimension(nGauss), private :: gaussx = (/-0.0765265211334973, 0.0765265211334973, &
                                                                -0.2277858511416451, 0.2277858511416451, &
                                                                -0.3737060887154195, 0.3737060887154195, &
                                                                -0.5108670019508271, 0.5108670019508271, &
                                                                -0.6360536807265150, 0.6360536807265150, &
                                                                -0.7463319064601508, 0.7463319064601508, &
                                                                -0.8391169718222188, 0.8391169718222188, &
                                                                -0.9122344282513259, 0.9122344282513259, &
                                                                -0.9639719272779138, 0.9639719272779138, &
                                                                -0.9931285991850949, 0.9931285991850949/)

  real(r8), parameter, dimension(nGauss), private :: gaussw = (/0.1527533871307258 , 0.1527533871307258, &
                                                                0.1491729864726037 , 0.1491729864726037, &
                                                                0.1420961093183820 , 0.1420961093183820, &
                                                                0.1316886384491766 , 0.1316886384491766, &
                                                                0.1181945319615184 , 0.1181945319615184, &
                                                                0.1019301198172404 , 0.1019301198172404, &
                                                                0.0832767415767048 , 0.0832767415767048, &
                                                                0.0626720483341091 , 0.0626720483341091, &
                                                                0.0406014298003869 , 0.0406014298003869, &
                                                                0.0176140071391521 , 0.0176140071391521/) 
#endif
#ifdef GAUSS_30
  integer , parameter  :: nGauss = 30
  real(r8), parameter, dimension(nGauss), private :: gaussx = (/-0.0514718425553177, 0.0514718425553177, &
                                                                -0.1538699136085835, 0.1538699136085835, &
                                                                -0.2546369261678899, 0.2546369261678899, &
                                                                -0.3527047255308781, 0.3527047255308781, &
                                                                -0.4470337695380892, 0.4470337695380892, &
                                                                -0.5366241481420199, 0.5366241481420199, &
                                                                -0.6205261829892429, 0.6205261829892429, &
                                                                -0.6978504947933158, 0.6978504947933158, &
                                                                -0.7677774321048262, 0.7677774321048262, &
                                                                -0.8295657623827684, 0.8295657623827684, &
                                                                -0.8825605357920527, 0.8825605357920527, &
                                                                -0.9262000474292743, 0.9262000474292743, &
                                                                -0.9600218649683075, 0.9600218649683075, &
                                                                -0.9836681232797472, 0.9836681232797472, &
                                                                -0.9968934840746495, 0.9968934840746495/)
                                                                                                          
  real(r8), parameter, dimension(nGauss), private :: gaussw = (/0.1028526528935588, 0.1028526528935588, &
                                                                0.1017623897484055, 0.1017623897484055, &
                                                                0.0995934205867953, 0.0995934205867953, &
                                                                0.0963687371746443, 0.0963687371746443, &
                                                                0.0921225222377861, 0.0921225222377861, &
                                                                0.0868997872010830, 0.0868997872010830, &
                                                                0.0807558952294202, 0.0807558952294202, &
                                                                0.0737559747377052, 0.0737559747377052, &
                                                                0.0659742298821805, 0.0659742298821805, &
                                                                0.0574931562176191, 0.0574931562176191, &
                                                                0.0484026728305941, 0.0484026728305941, &
                                                                0.0387991925696271, 0.0387991925696271, &
                                                                0.0287847078833234, 0.0287847078833234, &
                                                                0.0184664683110910, 0.0184664683110910, &
                                                                0.0079681924961666, 0.0079681924961666/)
#endif   
                                                                                                         
! temporaily value                                                                                       
   real(r8), parameter  ::  eps = 1e-12_r8

CONTAINS

!=======================================================================
!    Generate the supercell initial conditions
!=======================================================================
  SUBROUTINE supercell_init() &
    BIND(c, name = "supercell_init")

    IMPLICIT NONE

    ! d/dphi and int(dphi) operators
    real(r8), DIMENSION(nphi,nphi) :: ddphi, intphi

    ! d/dz and int(dz) operators
    real(r8), DIMENSION(nz, nz) :: ddz, intz

    ! Buffer matrices for computing SVD of d/dphi operator
    real(r8), DIMENSION(nphi,nphi) :: ddphibak
    real(r8), DIMENSION(nphi,nphi) :: svdpu, svdpvt
    real(r8), DIMENSION(nphi)      :: svdps
    real(r8), DIMENSION(5*nphi)    :: pwork

    ! Buffer matrices for computing SVD of d/dz operator
    real(r8), DIMENSION(nz, nz) :: ddzbak
    real(r8), DIMENSION(nz, nz) :: svdzu, svdzvt
    real(r8), DIMENSION(nz)     :: svdzs
    real(r8), DIMENSION(5*nz)   :: zwork

    ! Buffer data for calculation of SVD
    INTEGER(4) :: lwork, info

    ! Sampled values of ueq**2 and d/dz(ueq**2)
    real(r8), DIMENSION(nphi, nz) :: ueq2, dueq2

    ! Buffer matrices for iteration
    real(r8), DIMENSION(nphi, nz) :: phicoordmat, dztheta, rhs, irhs
  
    ! Buffer for sampled potential temperature at equator
    real(r8), DIMENSION(nz) :: thetaeq

    ! Buffer for computed equatorial Exner pressure and relative humidity
    real(r8), DIMENSION(nz) :: exnereq, H

    ! Variables for calculation of equatorial profile
    real(r8) :: exnereqs, p, T, qvs, qv

    ! Error metric
    real(r8) :: err

    ! Loop indices
    INTEGER(4) :: i, k, iter

    ! Chebyshev nodes in the phi direction
    do i = 1, nphi
      phicoord(i) = - cos(dble(i-1) * pi / dble(nphi-1))
      phicoord(i) = 0.25d0 * pi * (phicoord(i) + 1.0d0)
    end do

    ! Matrix of phis
    do k = 1, nz
      phicoordmat(:,k) = phicoord
    end do

    ! Chebyshev nodes in the z direction
    do k = 1, nz
      zcoord(k) = - cos(dble(k-1) * pi / dble(nz-1))
      zcoord(k) = z1 + 0.5d0*(z2-z1)*(zcoord(k)+1.0d0)
    end do

    ! Compute the d/dphi operator
    do i = 1, nphi
      call diff_lagrangian_polynomial_coeffs( &
        nphi, phicoord, ddphi(:,i), phicoord(i))
    end do

    ! Zero derivative at pole
    ddphi(:,nphi) = 0.0d0

    ! Compute the d/dz operator
    do k = 1, nz
      call diff_lagrangian_polynomial_coeffs( &
        nz, zcoord, ddz(:,k), zcoord(k))
    end do

    ! Compute the int(dphi) operator via pseudoinverse
    lwork = 5*nphi

    ddphibak = ddphi
    call DGESVD('A', 'A', &
       nphi, nphi, ddphibak, nphi, &
       svdps, svdpu, nphi, svdpvt, nphi, &
       pwork, lwork, info)

    if (info .ne. 0) then
      write(*,*) 'Unable to compute SVD of d/dphi matrix'
      stop
    end if

    do i = 1, nphi
      if (abs(svdps(i)) .le. 1.0d-12) then
        call DSCAL(nphi, 0.0d0, svdpu(1,i), 1)
      else
        call DSCAL(nphi, 1.0d0 / svdps(i), svdpu(1,i), 1)
      end if
    end do
    call DGEMM('T', 'T', &
      nphi, nphi, nphi, 1.0d0, svdpvt, nphi, svdpu, nphi, 0.0d0, &
      intphi, nphi)

    ! Compute the int(dz) operator via pseudoinverse
    lwork = 5*nz

    ddzbak = ddz
    call DGESVD('A', 'A', &
       nz, nz, ddzbak, nz, &
       svdzs, svdzu, nz, svdzvt, nz, &
       zwork, lwork, info)

    if (info .ne. 0) then
      write(*,*) 'Unable to compute SVD of d/dz matrix'
      stop
    end if

    do i = 1, nz
      if (abs(svdzs(i)) .le. 1.0d-12) then
        call DSCAL(nz, 0.0d0, svdzu(1,i), 1)
      else
        call DSCAL(nz, 1.0d0 / svdzs(i), svdzu(1,i), 1)
      end if
    end do
    call DGEMM('T', 'T', &
      nz, nz, nz, 1.0d0, svdzvt, nz, svdzu, nz, 0.0d0, &
      intz, nz)

    ! Sample the equatorial velocity field and its derivative
    do k = 1, nz
      ueq2(1,k) = zonal_velocity(zcoord(k), zero)
      ueq2(1,k) = ueq2(1,k)**2
    end do
    do k = 1, nz
      dueq2(1,k) = dot_product(ddz(:,k), ueq2(1,:))
    end do
    do i = 2, nphi
      ueq2(i,:) = ueq2(1,:)
      dueq2(i,:) = dueq2(1,:)
    end do

    ! Initialize potential temperature at equator
    do k = 1, nz
      thetaeq(k) = equator_theta(zcoord(k))
      H(k) = equator_relative_humidity(zcoord(k))
    end do
    thetavyz(1,:) = thetaeq

    ! Exner pressure at the equatorial surface
    exnereqs = (pseq / p0)**(Rd/cp)

    ! Iterate on equatorial profile
    do iter = 1, 12

      ! Calculate Exner pressure in equatorial column (p0 at surface)
      rhs(1,:) = - g / cp / thetavyz(1,:)
      do k = 1, nz
        exnereq(k) = dot_product(intz(:,k), rhs(1,:))
      end do
      do k = 2, nz
        exnereq(k) = exnereq(k) + (exnereqs - exnereq(1))
      end do
      exnereq(1) = exnereqs

      ! Calculate new pressure and temperature
      do k = 1, nz
        p = p0 * exnereq(k)**(cp/Rd)
        T = thetaeq(k) * exnereq(k)

        qvs = saturation_mixing_ratio(p, T)
        qveq(k) = qvs * H(k)

        thetavyz(1,k) = thetaeq(k) * (1.d0 + 0.61d0 * qveq(k))
      end do
    end do

    !do k = 1, nz
    !  write(*,*) exnereq(k) * thetaeq(k)
    !end do

    ! Iterate on remainder of domain
    do iter = 1, 12

      ! Compute d/dz(theta)
      do i = 1, nphi
        do k = 1, nz
          dztheta(i,k) = dot_product(ddz(:,k), thetavyz(i,:))
        end do
      end do

      ! Compute rhs
      rhs = sin(2.0d0*phicoordmat)/(2.0d0*g) &
            * (ueq2 * dztheta - thetavyz * dueq2)

      ! Integrate
      do k = 1, nz
        do i = 1, nphi
          irhs(i,k) = dot_product(intphi(:,i), rhs(:,k))
        end do
      end do

      ! Apply boundary conditions (fixed Dirichlet condition at equator)
      do i = 2, nphi
        irhs(i,:) = irhs(i,:) + (thetavyz(1,:) - irhs(1,:))
      end do
      irhs(1,:) = thetavyz(1,:)

      ! Compute difference after iteration
      !err = sum(irhs - thetavyz)
      !write(*,*) iter, err

      ! Update iteration
      thetavyz = irhs
    end do

    ! Calculate pressure through remainder of domain
    rhs = - ueq2 * sin(phicoordmat) * cos(phicoordmat) / cp / thetavyz

    do k = 1, nz
      do i = 1, nphi
        exneryz(i,k) = dot_product(intphi(:,i), rhs(:,k))
      end do
      do i = 2, nphi
        exneryz(i,k) = exneryz(i,k) + (exnereq(k) - exneryz(1,k))
      end do

      exneryz(1,k) = exnereq(k)
    end do

    ! Initialization successful
    initialized = 1

  END SUBROUTINE supercell_init

!-----------------------------------------------------------------------
!    Evaluate the supercell initial conditions
!-----------------------------------------------------------------------
  SUBROUTINE supercell_test(lon,lat,p,z,ptop,ztop,zcoords,u,v,t,thetav,thetap, ps,rho,q,pert) &
    BIND(c, name = "supercell_test")
 
    IMPLICIT NONE

    !------------------------------------------------
    !   Input / output parameters
    !------------------------------------------------
    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat           ! Latitude (radians)

    real(r8), INTENT(INOUT) :: &
                p,            & ! Pressure (Pa) or Dry Mass
                z,            & ! Altitude (m)
                ptop, ztop      ! p/m and z at model top, only used by supercell_m

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified
                                       ! 2 if dry mass coordinates are specified
    real(r8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                t,          & ! Temperature (K)
                thetav,     & ! Virtual potential Temperature (K)
                thetap,     & ! perturbed potential Temperature (K)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q             ! water vapor mixing ratio (kg/kg)

    INTEGER, INTENT(IN) :: pert  ! 1 if perturbation should be included
                                 ! 0 if no perturbation should be included

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------

    ! Absolute latitude
    real(r8) :: nh_lat

    ! Check that we are initialized
    if (initialized .ne. 1) then
      write(*,*) 'supercell_init() has not been called'
      stop
    end if

    !------------------------------------------------
    !   Begin sampling
    !------------------------------------------------

    ! Northern hemisphere latitude
    if (lat .le. 0.0d0) then
      nh_lat = -lat
    else
      nh_lat = lat
    end if

    ! Sample surface pressure
    CALL supercell_z(lon, lat, zero, ps, thetav, thetap, rho, q, pert)

    ! Calculate dependent variables
    if (zcoords .eq. 1) then
      CALL supercell_z(lon, lat, z, p, thetav, thetap, rho, q, pert)
    else if (zcoords .eq. 0) then
      CALL supercell_p(lon, lat, p, z, thetav, thetap, rho, q, pert)
    else if (zcoords .eq. 2) then
      CALL supercell_m(lon, lat, p, ptop, ztop, z, thetav, thetap, rho, q, pert)
    else
      print*,"supercell_test: zcoords number is wrong"
    end if

    ! Sample the zonal velocity
    u = zonal_velocity(z, lat)

    ! Zero meridional velocity
    v = 0.d0

    ! Temperature
    t = thetav / (1.d0 + 0.61d0 * q) * (p / p0)**(Rd/cp)

  END SUBROUTINE supercell_test

!-----------------------------------------------------------------------
!    Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
  SUBROUTINE supercell_z(lon, lat, z, p, thetav, thetap, rho, q, pert)

    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    INTEGER, INTENT(IN) :: pert  ! 1 if perturbation should be included
                                 ! 0 if no perturbation should be included

    ! Evaluated variables
    real(r8), INTENT(OUT) :: p, thetav, thetap, rho, q

    ! Northern hemisphere latitude
    real(r8) :: nh_lat

    ! Pointwise Exner pressure
    real(r8) :: exner

    ! Assembled variable values in a column
    real(r8), DIMENSION(nz) :: varcol

    ! Coefficients for computing a polynomial fit in each coordinate
    real(r8), DIMENSION(nphi) :: fitphi
    real(r8), DIMENSION(nz)   :: fitz

    ! Loop indices
    INTEGER(4) :: k

    ! Northern hemisphere latitude
    if (lat .le. 0.0d0) then
      nh_lat = -lat
    else
      nh_lat = lat
    end if

    ! Perform fit
    CALL lagrangian_polynomial_coeffs(nz, zcoord, fitz, z)
    CALL lagrangian_polynomial_coeffs(nphi, phicoord, fitphi, nh_lat)

    ! Obtain exner pressure of background state
    do k = 1, nz
      varcol(k) = dot_product(fitphi, exneryz(:,k))
    end do
    exner = dot_product(fitz, varcol)
    p = p0 * exner**(cp/Rd)

    ! Sample the initialized fit at this point for theta_v
    do k = 1, nz
      varcol(k) = dot_product(fitphi, thetavyz(:,k))
    end do
    thetav = dot_product(fitz, varcol)

    ! Sample water vapor mixing ratio
    q = dot_product(fitz, qveq)

    ! Fixed density
    rho = p / (Rd * exner * thetav)

    ! Modified virtual potential temperature
    if (pert .ne. 0) then
        thetav = thetav &
           + thermal_perturbation(lon, lat, z) * (1.d0 + 0.61d0 * q)
        thetap = thermal_perturbation(lon, lat, z) !* (1.d0 + 0.61d0 * q)
    end if

    ! Updated pressure
    p = p0 * (rho * Rd * thetav / p0)**(cp/(cp-Rd))

  END SUBROUTINE supercell_z

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  SUBROUTINE supercell_p(lon, lat, p, z, thetav, thetap, rho, q, pert)

    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p             ! Pressure (Pa)

    INTEGER, INTENT(IN) :: pert  ! 1 if perturbation should be included
                                 ! 0 if no perturbation should be included

    ! Evaluated variables
    real(r8), INTENT(OUT) :: z, thetav, thetap, rho, q

    ! Bounding interval and sampled values
    real(r8) :: za, zb, zc, pa, pb, pc

    ! Iterate
    INTEGER(4) :: iter

    za = z1
    zb = z2

    CALL supercell_z(lon, lat, za, pa, thetav, thetap, rho, q, pert)
    CALL supercell_z(lon, lat, zb, pb, thetav, thetap, rho, q, pert)

    if (pa .lt. p) then
      write(*,*) 'Requested pressure out of range on bottom, adjust sample interval'
      write(*,*) pa, p
      stop
    end if
    if (pb .gt. p) then
      write(*,*) 'Requested pressure out of range on top, adjust sample interval'
      write(*,*) pb, p
      stop
    end if

    ! Iterate using fixed point method
    !do iter = 1, 20
    zc = (za * (pb - p) - zb * (pa - p)) / (pb - pa)
    CALL supercell_z(lon, lat, zc, pc, thetav, thetap, rho, q, pert)

    do while(abs((pc - p) / p) .gt. 1.e-12)

      zc = (za * (pb - p) - zb * (pa - p)) / (pb - pa)

      CALL supercell_z(lon, lat, zc, pc, thetav, thetap, rho, q, pert)

      !write(*,*) pc

      !if (abs((pc - p) / p) .lt. 1.d-12) then
      !  exit
      !end if

      if (pc .gt. p) then
        za = zc
        pa = pc
      else
        zb = zc
        pb = pc
      end if
    end do

    !if (iter .eq. 21) then
    !  write(*,*) 'Iteration failed to converge->iter=',iter
    !  stop
    !end if

    z = zc

  END SUBROUTINE supercell_p

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  real(r8) FUNCTION thermal_perturbation(lon, lat, z)

    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    ! Great circle radius from the perturbation centerpoint
    real(r8) :: gr

    ! Approximately spherical radius from the perturbation centerpoint
    real(r8) :: Rtheta

    gr = a*acos(sin(pert_latc*deg2rad)*sin(lat) + &
         (cos(pert_latc*deg2rad)*cos(lat)*cos(lon-pert_lonc*deg2rad)))

    Rtheta = sqrt((gr/pert_rh)**2 + ((z - pert_zc) / pert_rz)**2)

    if (Rtheta .le. 1.d0) then
      thermal_perturbation = pert_dtheta * (cos(0.5d0 * pi * Rtheta))**2
    else
      thermal_perturbation = 0.0d0
    end if

  END FUNCTION thermal_perturbation

!-----------------------------------------------------------------------
!    Calculate the reference zonal velocity
!-----------------------------------------------------------------------
  real(r8) FUNCTION zonal_velocity(z, lat)

    IMPLICIT NONE

    real(r8), INTENT(IN) :: z, lat

    if (z .le. zs - zt) then
      zonal_velocity = us * (z / zs) - uc
    elseif (abs(z - zs) .le. zt) then
      zonal_velocity = &
        (-4.0d0/5.0d0 + 3.0d0*z/zs - 5.0d0/4.0d0*(z**2)/(zs**2)) * us - uc
    else
      zonal_velocity = us - uc
    end if

    zonal_velocity = zonal_velocity * cos(lat)
  
  END FUNCTION zonal_velocity

!-----------------------------------------------------------------------
!    Calculate pointwise theta at the equator at the given altitude
!-----------------------------------------------------------------------
  real(r8) FUNCTION equator_theta(z)

    IMPLICIT NONE

    real(r8), INTENT(IN) :: z

    if (z .le. z_tr) then
      equator_theta = &
        theta0 + (theta_tr - theta0) * (z / z_tr)**(1.25d0)
    else
      equator_theta = &
        theta_tr * exp(g/cp/T_tr * (z - z_tr))
    end if

  END FUNCTION equator_theta

!-----------------------------------------------------------------------
!    Calculate pointwise relative humidity (in %) at the equator at the
!    given altitude
!-----------------------------------------------------------------------
  real(r8) FUNCTION equator_relative_humidity(z)

    IMPLICIT NONE

    real(r8), INTENT(IN) :: z

    if (z .le. z_tr) then
      equator_relative_humidity = 1.0d0 - 0.75d0 * (z / z_tr)**(1.25d0)
    else
      equator_relative_humidity = 0.25d0
    end if

  END FUNCTION equator_relative_humidity

!-----------------------------------------------------------------------
!    Calculate saturation mixing ratio (in kg/kg) in terms of pressure
!    (in Pa) and temperature (in K)
!-----------------------------------------------------------------------
  real(r8) FUNCTION saturation_mixing_ratio(p, T)

    IMPLICIT NONE

    real(r8), INTENT(IN)  :: &
                p,        & ! Pressure in Pa
                T           ! Temperature

    saturation_mixing_ratio = &
      380.d0 / p * exp(17.27d0 * (T - 273.d0) / (T - 36.d0))

    if (saturation_mixing_ratio > 0.014) then
      saturation_mixing_ratio = 0.014
    end if

  END FUNCTION saturation_mixing_ratio

!-----------------------------------------------------------------------
!    Calculate coefficients for a Lagrangian polynomial
!-----------------------------------------------------------------------
  SUBROUTINE lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

    IMPLICIT NONE

    ! Number of points to fit
    INTEGER(4), INTENT(IN) :: npts

    ! Sample points to fit
    real(r8), DIMENSION(npts), INTENT(IN) :: x

    ! Computed coefficients
    real(r8), DIMENSION(npts), INTENT(OUT) :: coeffs

    ! Point at which sample is taken
    real(r8), INTENT(IN) :: xs

    ! Loop indices
    INTEGER(4) :: i, j
    
    ! Compute the Lagrangian polynomial coefficients
    do i = 1, npts
      coeffs(i) = 1.0d0
      do j = 1, npts
        if (i .eq. j) then
          cycle
        end if
        coeffs(i) = coeffs(i) * (xs - x(j)) / (x(i) - x(j))
      end do
    end do

  END SUBROUTINE lagrangian_polynomial_coeffs

!-----------------------------------------------------------------------
!    Calculate coefficients of the derivative of a Lagrangian polynomial
!-----------------------------------------------------------------------
  SUBROUTINE diff_lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

    IMPLICIT NONE

    ! Number of points to fit
    INTEGER(4), INTENT(IN) :: npts

    ! Sample points to fit
    real(r8), DIMENSION(npts), INTENT(IN) :: x

    ! Computed coefficients
    real(r8), DIMENSION(npts), INTENT(OUT) :: coeffs

    ! Point at which sample is taken
    real(r8), INTENT(IN) :: xs

    ! Loop indices
    INTEGER(4) :: i, j, imatch

    ! Buffer sum
    real(r8) :: coeffsum, differential

    ! Check if xs is equivalent to one of the values of x
    imatch = (-1)
    do i = 1, npts
      if (abs(xs - x(i)) < 1.0d-14) then
        imatch = i
        exit
      end if
    end do

    ! Equivalence detected; special treatment required
    if (imatch .ne. (-1)) then
      do i = 1, npts
        coeffs(i) = 1.0d0
        coeffsum = 0.0d0

        do j = 1, npts
          if ((j .eq. i) .or. (j .eq. imatch)) then
            cycle
          end if

          coeffs(i) = coeffs(i) * (xs - x(j)) / (x(i) - x(j))
          coeffsum = coeffsum + 1.0 / (xs - x(j))
        end do

        if (i .ne. imatch) then
          coeffs(i) = coeffs(i)                   &
            * (1.0 + (xs - x(imatch)) * coeffsum) &
            / (x(i) - x(imatch))
        else
          coeffs(i) = coeffs(i) * coeffsum
        end if
      end do

    ! No equivalence; simply differentiate Lagrangian fit
    else
      call lagrangian_polynomial_coeffs(npts, x, coeffs, xs)

      do i = 1, npts
        differential = 0.0d0
        do j = 1, npts
          if (i .eq. j) then
            cycle
          end if
          differential = differential + 1.0 / (xs - x(j))
        end do
        coeffs(i) = coeffs(i) * differential
      end do
    end if

  END SUBROUTINE diff_lagrangian_polynomial_coeffs

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given Dry mass (still p)
!-----------------------------------------------------------------------

  SUBROUTINE supercell_m(lon, lat, p, ptop, ztop, z, thetav, thetap, rho, q, pert)

    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p  ,        & ! Dry Mass (Pa)
                ptop,       & ! Dry mass at top
                ztop


    ! Evaluated variables
    real(r8), INTENT(OUT) :: z, thetav, thetap, rho, q
    INTEGER, INTENT(IN)  :: pert  ! 1 if perturbation should be included
                                  ! 0 if no perturbation should be included
    ! Bounding interval and sampled values
    real(r8) :: za, zb, zc, pa, pb, pc

    ! Iterate
    INTEGER(4) :: iter

    za = z1
    zb = z2

    !CALL supercell_z(lon, lat, za, pa, thetav, thetap, rho, q, pert)
    !CALL supercell_z(lon, lat, zb, pb, thetav, thetap, rho, q, pert)
    pa = get_dryAirMass_from_z(lon, lat, za, ptop, ztop, pert)
    pb = get_dryAirMass_from_z(lon, lat, zb, ptop, ztop, pert)

    if (pa .lt. p) then
      write(*,*) 'Requested pressure out of range on bottom, adjust sample interval'
      write(*,*) pa, p
      stop
    end if
    if (pb .gt. p) then
      write(*,*) 'Requested pressure out of range on top, adjust sample interval'
      write(*,*) pb, p
      stop
    end if

    ! Iterate using fixed point method
    !do iter = 1, 20
    zc = (za * (pb - p) - zb * (pa - p)) / (pb - pa)
    !CALL supercell_z(lon, lat, zc, pc, thetav, thetap, rho, q, pert)
    pc = get_dryAirMass_from_z(lon, lat, zc, ptop, ztop, pert)

    do while(abs((pc - p) / p) .gt. eps)

      zc = (za * (pb - p) - zb * (pa - p)) / (pb - pa)

      !CALL supercell_z(lon, lat, zc, pc, thetav, thetap, rho, q, pert)
      pc = get_dryAirMass_from_z(lon, lat, zc, ptop, ztop, pert)

      !write(*,*) pc

      !if (abs((pc - p) / p) .lt. 1.d-12) then
      !  exit
      !end if

      if (pc .gt. p) then
        za = zc
        pa = pc
      else
        zb = zc
        pb = pc
      end if
    end do

    !if (iter .eq. 21) then
    !  write(*,*) 'Iteration failed to converge->iter=',iter
    !  stop
    !end if

    z = zc
! use supercell_z to retrieve final
    CALL supercell_z(lon, lat, z, pc, thetav, thetap, rho, q, pert)

  END SUBROUTINE supercell_m

!------------------------------------------------------------------
!  Vertical Gaussian quadrature for integrating dry air mass,
!  based on full pressure, Rd, Tv (imply alpham using gas law),
!  specific humidity (not dry mixing ratio)
!------------------------------------------------------------------

  real(r8) function get_dryAirMass_from_z(lon, lat, z, ptop,  ztop, pert)
! io
    real(r8),    intent(in)  :: lon
    real(r8),    intent(in)  :: lat
    real(r8),    intent(in)  :: z
    real(r8),    intent(in)  :: ptop
    real(r8),    intent(in)  :: ztop
    integer(i4), intent(in)  :: pert 
! local
    real(r8)                 :: xm, xr, integral
    real(r8)                 :: qv, z1, z2, thetav, thetap, rho, Tv, pfull, ztmp
    integer                  :: jgw

      z1       = z
      z2       = ztop
      xm       = 0.5*(z1+z2)
      xr       = 0.5*(z2-z1)
      integral = zero

      do jgw = 1, nGauss
         ztmp = xm+gaussx(jgw)*xr
         call supercell_z(lon, lat, ztmp, pfull, thetav, thetap, rho, qv, pert)
         Tv       = thetav*(pfull/p0)**(rdry/cp)
         integral = integral + gaussw(jgw)*gravity*pfull*(one-qv)/(rdry*Tv)
      enddo
      integral             = 0.5_r8*(z2-z1)*integral
      get_dryAirMass_from_z = integral+ptop
  end function get_dryAirMass_from_z

END MODULE grist_dtp_dcmip2016_sc
