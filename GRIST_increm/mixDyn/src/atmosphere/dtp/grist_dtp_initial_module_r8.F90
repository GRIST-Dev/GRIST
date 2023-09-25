module grist_dtp_initial_module_r8

    use grist_constants,             only: i4, rearth, gravity, omega, rdry, pi, rvap, cp
    use grist_domain_types,          only: global_domain
! hpe
    use grist_hpe_constants,         only: eta_face_a, eta_face_b, eta_full_a, eta_full_b, p0

    implicit none

    ! private

    integer, parameter      ::  r8 = 8
    real(r8), parameter     ::  a     = rearth , & !6371220.0d0,           & ! Reference Earth's Radius (m)
                                Rd    = rdry   , & !287.0d0,               & ! Ideal gas const dry air (J kg^-1 K^1)
                                g     = gravity, & !9.80616d0,             & ! Gravity (m s^2)
                                Mvap  = 0.608d0,               & ! Ratio of molar mass of dry air/water
                                T0E        = 310.d0     ,      & ! temperature at equatorial surface (K)
                                T0P        = 240.d0     ,      & ! temperature at polar surface (K)
                                B          = 2.d0       ,      & ! jet half-width parameter
                                K          = 3.d0       ,      & ! jet width parameter
                                pertu0     = 0.5d0      ,      & ! SF Perturbation wind velocity (m/s)
                                pertr      = 1.d0/6.d0  ,      & ! SF Perturbation radius (Earth radii)
                                pertup     = 1.0d0      ,      & ! Exp. perturbation wind velocity (m/s)
                                pertexpr   = 0.1d0      ,      & ! Exp. perturbation radius (Earth radii)
#ifdef IDEAL_CASE
                                pertlon    = pi/3.d0    ,      & ! Perturbation longitude
#else
                                pertlon    = pi/9.d0    ,      & ! Perturbation longitude
#endif
                                pertlat    = 2.d0*pi/9.d0,     & ! Perturbation latitude
                                pertz      = 15000.d0   ,      & ! Perturbation height cap
                                dxepsilon  = 1.d-5      ,      & ! Small value for numerical derivatives
                                moistqp    = 34000.d0,         & ! Humidity vertical pressure width
                                moisttr    = 0.1d0,            & ! Vertical cut-off pressure for humidity
                                moistqs    = 1.d-12,           & ! Humidity above cut-off
                                moistq0    = 0.018d0,          & ! Maximum specific humidity
                                moistqlat  = 2.d0*pi/9.d0,     & ! Humidity latitudinal width
                                lapse      = 0.005d0             ! lapse rate parameter

    real(r8)  :: one  = 1._r8
    real(r8)  :: zero = 0._r8
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

contains

!------------------------------------------------------------------
!  evaluate dry mass and geometric height at each level
!------------------------------------------------------------------

  subroutine evaluate_dry_mass_height_r8(mesh, ncell, nlev, testcase         , &
                                         scalar_hpressure_at_pc_surface      , &
                                         scalar_hpressure_at_pc_face_level   , &
                                         scalar_hpressure_at_pc_full_level   , &
                                         scalar_zzz_at_pc_face_level         , &
                                         scalar_zzz_at_pc_full_level)

    use grist_constants, only: r4 => r8

    implicit none
!io
   type(global_domain),  intent(in)   :: mesh
   integer(i4),          intent(in)   :: ncell
   integer(i4),          intent(in)   :: nlev
   character(len=*),     intent(in)   :: testcase
   real(r4),             intent(out)  :: scalar_hpressure_at_pc_surface(:)
   real(r4),             intent(out)  :: scalar_hpressure_at_pc_face_level(:,:)
   real(r4),             intent(out)  :: scalar_hpressure_at_pc_full_level(:,:)
   real(r4),             intent(out)  :: scalar_zzz_at_pc_face_level(:,:)
   real(r4),             intent(out)  :: scalar_zzz_at_pc_full_level(:,:)
! local
   real(r8), allocatable              :: scalar_hpressure_at_pc_surface_r8(:)
   real(r8), allocatable              :: scalar_hpressure_at_pc_face_level_r8(:,:)
   real(r8), allocatable              :: scalar_hpressure_at_pc_full_level_r8(:,:)
   real(r8), allocatable              :: scalar_zzz_at_pc_face_level_r8(:,:)
   real(r8), allocatable              :: scalar_zzz_at_pc_full_level_r8(:,:)

   allocate(scalar_hpressure_at_pc_surface_r8,    source = real(scalar_hpressure_at_pc_surface, r8)   )
   allocate(scalar_hpressure_at_pc_face_level_r8, source = real(scalar_hpressure_at_pc_face_level, r8))
   allocate(scalar_hpressure_at_pc_full_level_r8, source = real(scalar_hpressure_at_pc_full_level, r8))
   allocate(scalar_zzz_at_pc_face_level_r8,       source = real(scalar_zzz_at_pc_face_level, r8)      )
   allocate(scalar_zzz_at_pc_full_level_r8,       source = real(scalar_zzz_at_pc_full_level, r8)      )

   call evaluate_dry_mass_height(mesh, ncell, nlev, testcase            , &
                                 scalar_hpressure_at_pc_surface_r8      , &
                                 scalar_hpressure_at_pc_face_level_r8   , &
                                 scalar_hpressure_at_pc_full_level_r8   , &
                                 scalar_zzz_at_pc_face_level_r8         , &
                                 scalar_zzz_at_pc_full_level_r8)

    scalar_hpressure_at_pc_surface    = scalar_hpressure_at_pc_surface_r8
    scalar_hpressure_at_pc_face_level = scalar_hpressure_at_pc_face_level_r8
    scalar_hpressure_at_pc_full_level = scalar_hpressure_at_pc_full_level_r8
    scalar_zzz_at_pc_face_level       = scalar_zzz_at_pc_face_level_r8
    scalar_zzz_at_pc_full_level       = scalar_zzz_at_pc_full_level_r8

  end subroutine evaluate_dry_mass_height_r8

  subroutine evaluate_dry_mass_height(mesh, ncell, nlev, testcase         , &
                                      scalar_hpressure_at_pc_surface      , &
                                      scalar_hpressure_at_pc_face_level   , &
                                      scalar_hpressure_at_pc_full_level   , &
                                      scalar_zzz_at_pc_face_level         , &
                                      scalar_zzz_at_pc_full_level)
!io
   type(global_domain),  intent(in)   :: mesh
   integer(i4),          intent(in)   :: ncell
   integer(i4),          intent(in)   :: nlev
   character(len=*),     intent(in)   :: testcase
   real(r8),             intent(out)  :: scalar_hpressure_at_pc_surface(:)
   real(r8),             intent(out)  :: scalar_hpressure_at_pc_face_level(:,:)
   real(r8),             intent(out)  :: scalar_hpressure_at_pc_full_level(:,:)
   real(r8),             intent(out)  :: scalar_zzz_at_pc_face_level(:,:)
   real(r8),             intent(out)  :: scalar_zzz_at_pc_full_level(:,:)
! local
   integer(i4)         :: iv, ilev
   real(r8)            :: pres, dr_mass, wv_mass, ptop, ztop
!
! obtain top level height
!
    ptop = eta_face_a(1)*p0

    DO iv = 1, ncell
! for each profile
       ztop  = get_ztop(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), ptop, testcase) ! analytic
       !ztop  = get_z_from_pressure(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat, ptop, ptop, 0.0_r8, testcase,.false.) ! iteration

       !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
       !ztop1 = ztop(iv)
       !call  supercell_test(mesh%vtx(iv)%lon, mesh%vtx(iv)%lat, ptop, ztop(iv), 0, tmp1, tmp2, tmp3, tmp4, &
       !                       tmp5, tmp6, tmp7, tmp8, 1)
       !ztop2 = ztop(iv)
       !ztop2 = ztop1-ztop2
       !print*,"ztop2=",ztop2, such difference is 1e-9 for supercell due to different configuration for interation
       !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!
! Gaussian quadrature for surface dry and wv mass, their sum should be close to analytic (height=0) produced pressure;
! this has been confirmed;
!
       dr_mass  = get_dryAirMass_from_z(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8),zero,ptop,ztop,testcase)
       !wv_mass  = get_waterVaporMass_from_z(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,zero,ztop,testcase)
       !pres     =   get_dcmipPressure_from_z(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,zero     ,testcase)
       !scalar_hpressure_at_pc_surface(iv) = pres-wv_mass
       scalar_hpressure_at_pc_surface(iv) = dr_mass
!
!  get surface hpressure based on analytical function
!
#ifdef DTPHPS
       scalar_hpressure_at_pc_surface(iv) = get_surface_hps(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8),testcase)
#endif
!
! below are same
!
       do ilev = 1, nlev+1
          scalar_hpressure_at_pc_face_level(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*scalar_hpressure_at_pc_surface(iv)
       end do

       do ilev = 1, nlev
          scalar_hpressure_at_pc_full_level(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)*scalar_hpressure_at_pc_surface(iv)
       end do
!
! iterate to obtain z-face
! 
       do ilev = 1, nlev+1
          scalar_zzz_at_pc_face_level(ilev,iv)  = get_z_from_pressure(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), &
                                                  scalar_hpressure_at_pc_face_level(ilev,iv), ptop, ztop, testcase,.true.)
       end do
       if(scalar_zzz_at_pc_face_level(nlev+1,iv).ne.0._r8)then
          !print*,"reset to zero", scalar_zzz_at_pc_face_level(nlev+1,iv)
          scalar_zzz_at_pc_face_level(nlev+1,iv) = zero
       end if

       do ilev = 1, nlev
          scalar_zzz_at_pc_full_level(ilev,iv)  = 0.5_r8*(scalar_zzz_at_pc_face_level(ilev,iv)+scalar_zzz_at_pc_face_level(ilev+1,iv))
       end do

    END DO

   return
  end subroutine evaluate_dry_mass_height

!------------------------------------------------------------------
! Evaluate model level height based on pressure or dry mass (dry 
! hydrostatic pressure) using fixed-point iteration; if dry mass 
! is used, the function is an vertical Gauss-integral function; 
! if (full) pressure is used, use the DCMIP analytic function now
!------------------------------------------------------------------

  real(r8) function get_z_from_pressure(lon, lat, pres, ptop, ztop, testcase, is_dry_mass)
! io
    real(r8), intent(in)  ::  lon
    real(r8), intent(in)  ::  lat
    real(r8), intent(in)  ::  pres     ! Pressure (Pa)
    real(r8), intent(in)  ::  ptop     ! Pressure (Pa)
    real(r8), intent(in)  ::  ztop
    character(len=*), intent(in) :: testcase
    logical,  intent(in)  :: is_dry_mass
! local
    integer    :: ix
    real(r8)   :: z0, z1, z2, za, zb, zc
    real(r8)   :: p0, p1, p2, pa, pb, pc


    select case(trim(testcase))

    case('DCMIP2016-BW','DCMIP2016-TC')
!
! some initial
!
    z0 = 0.0_r8
    z1 = 10000.0_r8

    if (is_dry_mass) then
       p0 = get_dryAirMass_from_z(lon,lat,z0,ptop,ztop,testcase)
       p1 = get_dryAirMass_from_z(lon,lat,z1,ptop,ztop,testcase)
    else
       p0 = get_dcmipPressure_from_z(lon,lat,z0,testcase)
       p1 = get_dcmipPressure_from_z(lon,lat,z1,testcase)
    endif
!
! fixed point iteration
!
    DO ix = 1, 1000
       z2 = z1 - (p1 - pres) * (z1 - z0) / (p1 - p0)
       if (is_dry_mass) then
          p2 = get_dryAirMass_from_z(lon,lat,z2,ptop,ztop,testcase)
       else
          p2 = get_dcmipPressure_from_z(lon,lat,z2,testcase)
       end if

       IF (ABS(p2 - pres)/pres < eps.or.ABS(z1-z2)<eps.or.ABS(p1-p2)<eps) THEN
          EXIT
       END IF

       z0 = z1
       p0 = p1

       z1 = z2
       p1 = p2
    END DO

    if (ix==1001) then
      write(*,*) "pres,p1,z1",pres,p1,z1
      print*, 'iteration did not converge in get_z_from_pressure'
      call mpi_abort
    end if
    get_z_from_pressure = z2

    case default

    end select

  end function get_z_from_pressure

!------------------------------------------------------------------
!  Vertical Gaussian quadrature for integrating dry air mass,
!  based on full pressure, Rd, Tv (imply alpham using gas law),
!  specific humidity (not dry mixing ratio)
!------------------------------------------------------------------

  real(r8) function get_dryAirMass_from_z(lon, lat, z, ptop,  ztop, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    real(r8), intent(in)  :: z
    real(r8), intent(in)  :: ptop
    real(r8), intent(in)  :: ztop
    character(len=*), intent(in) :: testcase
! local
    real(r8)              :: xm, xr, integral
    real(r8)              :: qv, z1, z2, Tv, pfull, ztmp
    integer               :: jgw
    real(r8)              :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

      z1       = z
      z2       = ztop
      xm       = 0.5*(z1+z2)
      xr       = 0.5*(z2-z1)
      integral = zero

      select case (trim(testcase))
      case('DCMIP2016-BW')
         do jgw = 1, nGauss
            ztmp = xm + gaussx(jgw)*xr
            call baroclinic_wave_test(0, 1, 0, one, lon, lat, pfull, ztmp, 1, &
                                      tmp1, tmp2, tmp3, Tv, tmp5, tmp6, tmp7, qv)
            !this is actually ptv, transform to Tv
             Tv       = Tv*(pfull/p0)**(rdry/cp)
             integral = integral + gaussw(jgw)*gravity*pfull*(one-qv)/(rdry*Tv)
         end do

      case default
         print*,"get_dryAirMass_from_z: this cannot happen"
      end select

      integral             = 0.5_r8*(z2-z1)*integral
      get_dryAirMass_from_z = integral+ptop
  end function get_dryAirMass_from_z

!=======================================================================
!    Generate the baroclinic instability initial conditions
!=======================================================================
  SUBROUTINE baroclinic_wave_test(deep,moist,pertt,X,lon,lat,p,z,zcoords,u,v,t,thetav,phis,ps,rho,q) &
    BIND(c, name = "baroclinic_wave_test_r8")
 
    IMPLICIT NONE

!-----------------------------------------------------------------------
!     input/output params parameters at given location
!-----------------------------------------------------------------------
    INTEGER, INTENT(IN)  :: &
                deep,       & ! Deep (1) or Shallow (0) test case
                moist,      & ! Moist (1) or Dry (0) test case
                pertt         ! Perturbation type

    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                X             ! Earth scaling parameter

    real(r8), INTENT(INOUT) :: &
                p,            & ! Pressure (Pa)
                z               ! Altitude (m)

    INTEGER, INTENT(IN) :: zcoords     ! 1 if z coordinates are specified
                                       ! 0 if p coordinates are specified

    real(r8), INTENT(OUT) :: &
                u,          & ! Zonal wind (m s^-1)
                v,          & ! Meridional wind (m s^-1)
                t,          & ! Temperature (K)
                thetav,     & ! Virtual potential temperature (K)
                phis,       & ! Surface Geopotential (m^2 s^-2)
                ps,         & ! Surface Pressure (Pa)
                rho,        & ! density (kg m^-3)
                q             ! water vapor mixing ratio (kg/kg)

    !------------------------------------------------
    !   Local variables
    !------------------------------------------------
    real(r8) :: aref, omegaref
    real(r8) :: T0, constH, constC, scaledZ, inttau2, rratio
    real(r8) :: inttermU, bigU, rcoslat, omegarcoslat
    real(r8) :: eta, qratio, qnum, qden

    !------------------------------------------------
    !   Pressure and temperature
    !------------------------------------------------
    if (zcoords .eq. 1) then
      CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)
    else
      CALL evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    end if

    !------------------------------------------------
    !   Compute test case constants
    !------------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)

    constH = Rd * T0 / g

    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)

    scaledZ = z / (B * constH)

    inttau2 = constC * z * exp(- scaledZ**2)

    ! radius ratio
    if (deep .eq. 0) then
      rratio = 1.d0
    else
      rratio = (z + aref) / aref;
    end if

    !-----------------------------------------------------
    !   Initialize surface pressure
    !-----------------------------------------------------
    ps = p0

    !-----------------------------------------------------
    !   Initialize velocity field
    !-----------------------------------------------------
    inttermU = (rratio * cos(lat))**(K - 1.d0) - (rratio * cos(lat))**(K + 1.d0)
    bigU = g / aref * K * inttau2 * inttermU * t
    if (deep .eq. 0) then
      rcoslat = aref * cos(lat)
    else
      rcoslat = (z + aref) * cos(lat)
    end if

    omegarcoslat = omegaref * rcoslat
    
    u = - omegarcoslat + sqrt(omegarcoslat**2 + rcoslat * bigU)
    v = 0.d0

    !-----------------------------------------------------
    !   Add perturbation to the velocity field
    !-----------------------------------------------------

    ! Exponential type
    if (pertt .eq. 0) then
      u = u + evaluate_exponential(lon, lat, z)

    ! Stream function type
    elseif (pertt .eq. 1) then
      u = u - 1.d0 / (2.d0 * dxepsilon) *                       &
          ( evaluate_streamfunction(lon, lat + dxepsilon, z)    &
          - evaluate_streamfunction(lon, lat - dxepsilon, z))

      v = v + 1.d0 / (2.d0 * dxepsilon * cos(lat)) *            &
          ( evaluate_streamfunction(lon + dxepsilon, lat, z)    &
          - evaluate_streamfunction(lon - dxepsilon, lat, z))
    end if

    !-----------------------------------------------------
    !   Initialize surface geopotential
    !-----------------------------------------------------
    phis = 0.d0

    !-----------------------------------------------------
    !   Initialize density
    !-----------------------------------------------------
    rho = p / (Rd * t)

    !-----------------------------------------------------
    !   Initialize specific humidity
    !-----------------------------------------------------
    if (moist .eq. 1) then
      eta = p/p0

      if (eta .gt. moisttr) then
        q = moistq0 * exp(- (lat/moistqlat)**4)          &
                    * exp(- ((eta-1.d0)*p0/moistqp)**2)
      else
        q = moistqs
      end if

      ! Convert virtual temperature to temperature
      t = t / (1.d0 + Mvap * q)

    else
      q = 0.d0
    end if

    !-----------------------------------------------------
    !   Initialize virtual potential temperature
    !-----------------------------------------------------
    thetav = t * (1.d0 + 0.61d0 * q) * (p0 / p)**(Rd / cp)

  END SUBROUTINE baroclinic_wave_test
  
!-----------------------------------------------------------------------
!    Calculate pointwise pressure and temperature
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_pressure_temperature(deep, X, lon, lat, z, p, t)

    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    real(r8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (m)

    real(r8), INTENT(OUT) :: &
                p,          & ! Pressure (Pa)
                t             ! Temperature (K)

    real(r8) :: aref, omegaref
    real(r8) :: T0, constA, constB, constC, constH, scaledZ
    real(r8) :: tau1, tau2, inttau1, inttau2
    real(r8) :: rratio, inttermT

    !--------------------------------------------
    ! Constants
    !--------------------------------------------
    aref = a / X
    omegaref = omega * X

    T0 = 0.5d0 * (T0E + T0P)
    constA = 1.d0 / lapse
    constB = (T0 - T0P) / (T0 * T0P)
    constC = 0.5d0 * (K + 2.d0) * (T0E - T0P) / (T0E * T0P)
    constH = Rd * T0 / g

    scaledZ = z / (B * constH)

    !--------------------------------------------
    !    tau values
    !--------------------------------------------
    tau1 = constA * lapse / T0 * exp(lapse * z / T0) &
         + constB * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)
    tau2 = constC * (1.d0 - 2.d0 * scaledZ**2) * exp(- scaledZ**2)

    inttau1 = constA * (exp(lapse * z / T0) - 1.d0) &
            + constB * z * exp(- scaledZ**2)
    inttau2 = constC * z * exp(- scaledZ**2)

    !--------------------------------------------
    !    radius ratio
    !--------------------------------------------
    if (deep .eq. 0) then
      rratio = 1.d0
    else
      rratio = (z + aref) / aref;
    end if

    !--------------------------------------------
    !    interior term on temperature expression
    !--------------------------------------------
    inttermT = (rratio * cos(lat))**K &
             - K / (K + 2.d0) * (rratio * cos(lat))**(K + 2.d0)

    !--------------------------------------------
    !    temperature
    !--------------------------------------------
    t = 1.d0 / (rratio**2 * (tau1 - tau2 * inttermT))

    !--------------------------------------------
    !    hydrostatic pressure
    !--------------------------------------------
    p = p0 * exp(- g / Rd * (inttau1 - inttau2 * inttermT))

  END SUBROUTINE evaluate_pressure_temperature

!-----------------------------------------------------------------------
!    Calculate pointwise z and temperature given pressure
!-----------------------------------------------------------------------
  SUBROUTINE evaluate_z_temperature(deep, X, lon, lat, p, z, t)
    
    INTEGER, INTENT(IN)  :: deep ! Deep (1) or Shallow (0) test case

    real(r8), INTENT(IN)  :: &
                X,          & ! Earth scaling ratio
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                p             ! Pressure (Pa)

    real(r8), INTENT(OUT) :: &
                z,          & ! Altitude (m)
                t             ! Temperature (K)

    INTEGER :: ix

    real(r8) :: z0, z1, z2
    real(r8) :: p0, p1, p2

    z0 = 0.d0
    z1 = 10000.d0

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z0, p0, t)
    CALL evaluate_pressure_temperature(deep, X, lon, lat, z1, p1, t)

    DO ix = 1, 100
      z2 = z1 - (p1 - p) * (z1 - z0) / (p1 - p0)

      CALL evaluate_pressure_temperature(deep, X, lon, lat, z2, p2, t)

      IF (ABS((p2 - p)/p) .lt. 1.0d-13) THEN
        EXIT
      END IF

      z0 = z1
      p0 = p1

      z1 = z2
      p1 = p2
    END DO

    z = z2

    CALL evaluate_pressure_temperature(deep, X, lon, lat, z, p0, t)

  END SUBROUTINE evaluate_z_temperature

!-----------------------------------------------------------------------
!    Exponential perturbation function
!-----------------------------------------------------------------------
  real(r8) FUNCTION evaluate_exponential(lon, lat, z)

    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    real(r8) :: greatcircler, perttaper

    ! Great circle distance
    greatcircler = 1.d0 / pertexpr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
      perttaper = 0.d0
    end if

    ! Zonal velocity perturbation
    if (greatcircler < 1.d0) then
      evaluate_exponential = pertup * perttaper * exp(- greatcircler**2)
    else
      evaluate_exponential = 0.d0
    end if

  END FUNCTION evaluate_exponential
  
!-----------------------------------------------------------------------
!    Stream function perturbation function
!-----------------------------------------------------------------------
  real(r8) FUNCTION evaluate_streamfunction(lon, lat, z)

    real(r8), INTENT(IN)  :: &
                lon,        & ! Longitude (radians)
                lat,        & ! Latitude (radians)
                z             ! Altitude (meters)

    real(r8) :: greatcircler, perttaper, cospert

    ! Great circle distance
    greatcircler = 1.d0 / pertr &
      * acos(sin(pertlat) * sin(lat) + cos(pertlat) * cos(lat) * cos(lon - pertlon))

    ! Vertical tapering of stream function
    if (z < pertz) then
      perttaper = 1.d0 - 3.d0 * z**2 / pertz**2 + 2.d0 * z**3 / pertz**3
    else
      perttaper = 0.d0
    end if

    ! Horizontal tapering of stream function
    if (greatcircler .lt. 1.d0) then
      cospert = cos(0.5d0 * pi * greatcircler)
    else
      cospert = 0.d0
    end if

    evaluate_streamfunction = &
        (- pertu0 * pertr * perttaper * cospert**4)

  END FUNCTION evaluate_streamfunction

  real(r8) function get_ztop(lon, lat, ptop, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    real(r8), intent(in)  :: ptop
    character(len=*), intent(in) :: testcase
! local
    real(r8)  :: pfull, zlocal
    real(r8)  :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

    pfull = ptop

    select case (trim(testcase))
    case('DCMIP2016-BW')
        call baroclinic_wave_test(0, 1, 0, one, lon, lat , pfull, zlocal, 0, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
    end select
    get_ztop = zlocal
 
  end function get_ztop


!------------------------------------------------------------------
!  produce pressure based on DCMIP analytic functions with height
!------------------------------------------------------------------

  real(r8) function get_dcmipPressure_from_z(lon, lat, z, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    real(r8), intent(in)  :: z
    character(len=*), intent(in) :: testcase
! local
    real(r8)  :: pfull, zlocal
    real(r8)  :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
    zlocal = z

    select case (trim(testcase))
    case('DCMIP2016-BW')
        call baroclinic_wave_test(0, 1, 0, one, lon, lat , pfull, zlocal, 1, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
    case default
       print*,"get_dcmipPressure_from_z: this cannot happen"
    end select
    get_dcmipPressure_from_z = pfull
 
  end function get_dcmipPressure_from_z

!------------------------------------------------------------------
!  produce pressure based on DCMIP analytic functions with height
!------------------------------------------------------------------

  real(r8) function get_surface_hps(lon, lat, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    character(len=*), intent(in) :: testcase
! local
    real(r8)  :: pfull, zlocal
    real(r8)  :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
    zlocal = zero

    select case (trim(testcase))
    case('DCMIP2016-BW')
        call baroclinic_wave_test(0, 1, 0, one, lon, lat , pfull, zlocal, 1, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
        pfull = pfull/(one+rvap/rdry*tmp8)
    case default
       print*,"get_surface_hps: this cannot happen"
    end select
    get_surface_hps = pfull
 
  end function get_surface_hps

end module grist_dtp_initial_module_r8