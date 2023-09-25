module grist_constants

  implicit none
  save 
  public

  !---------------------------------------------------
  !                 Kind attributions
  !---------------------------------------------------

  integer,  parameter  :: i2      = selected_real_kind(2,20)
  integer,  parameter  :: i4      = selected_real_kind(6,20)
  integer,  parameter  :: i8      = selected_real_kind(14,40)
  integer,  parameter  :: r4      = selected_real_kind(6,37)
#ifdef SPCODE
  integer,  parameter  :: r8      = selected_real_kind(6,37)
#else
  integer, parameter   :: r8      = selected_real_kind(12,100)
#endif
#ifdef MIXCODE

#ifdef MIXCHECK
  integer, parameter   :: ns      = r8
#else
  integer, parameter   :: ns      = r4
#endif

#elif SPCODE
  integer, parameter   :: ns      = r4
#else
  integer, parameter   :: ns      = r8
#endif
  integer,  parameter  :: r16     = max(r8,selected_real_kind(27,2400))

  real(r8), parameter  :: zero    = 0._r8
  real(ns), parameter  :: zero_ns = 0._ns
  real(r8), parameter  :: half    = 0.5_r8
  real(ns), parameter  :: half_ns = 0.5_ns
  real(r8), parameter  :: one     = 1._r8
  real(ns), parameter  :: one_ns  = 1._ns
  real(r8), parameter  :: two     = 2._r8
  real(ns), parameter  :: two_ns  = 2._ns
  real(r8), parameter  :: pi      = 3.14159265358979323846_r8
  real(ns), parameter  :: pi_ns   = pi
  real(r8), parameter  :: pi2     = 2._r8*pi
  real(r8), parameter  :: pio2    = pi/2._r8
  real(r8), parameter  :: piby2   = pi*0.5_r8
  real(r8), parameter  :: pio4    = pi/4._r8
  real(r8), parameter  :: deg2rad = pi / 180._r8
  real(r8), parameter  :: rad2deg = 1._r8/deg2rad
  real(r8), parameter  :: eps     = epsilon(one)
  !real(r8), parameter  :: eps     = epsilon(1.)
  real(r8), parameter  :: eps2    = epsilon(pi)
#ifdef SMALL_EARTH
#ifdef DCMIP21
  real(r8), parameter  :: rearth  = 6.37122e6_r8/166.7_r8
#endif
#ifdef DCMIP31
  real(r8), parameter  :: rearth  = 6.37122e6_r8/125._r8
#endif
#ifdef DCMIP_SUPERCELL
  real(r8), parameter  :: rearth  = 6.37122e6_r8/120._r8
#endif
#else
  real(r8), parameter  :: rearth  = 6.37122e6_r8
#endif
  real(ns), parameter  :: rearth_ns = rearth
  real(r8), parameter  :: eradi   = 1._r8/6.37122e6_r8
  real(r8), parameter  :: gravity = 9.80616_r8
  real(r8), parameter  :: gravi   = 1._r8/9.80616_r8

  real(ns), parameter  :: gravity_ns = gravity
! Angular velocity of the Earth (rot/s)
#ifdef NONROT
  real(r8), parameter  :: omega   = 0._r8
#else
  real(r8), parameter  :: omega   = 7.2921e-5_r8
#endif
  real(ns), parameter  :: omega_ns = omega
  real(r8), parameter  :: day2sec  = 86400_r8
  real(r8), parameter  :: sec2day  = 1._r8/86400_r8
#ifndef OLDSIX
  real(r8), parameter  :: rdry     = 287._r8 
  real(r8), parameter  :: cp       = 1004.5_r8 ! dry air
  real(r8), parameter  :: cv       = 717.5_r8  ! dry air
  real(r8), parameter  :: rvap     = 461.5_r8
  real(r8), parameter  :: p00      = 1.e5_r8         
  real(r8), parameter  :: t00      = 273.15_r8
#endif
  real(r8), parameter  :: cpv      = 1810._r8
  real(r8), parameter  :: cvv      = 1348.5_r8
#ifdef OLDSIX
  real(r8), parameter  :: rdry     = 287.
  real(r8), parameter  :: cp       = 1004.
  real(r8), parameter  :: cv       = 717.
  real(r8), parameter  :: rvap     = 461.
  real(r8), parameter  :: p00      = 1.e5
  real(r8), parameter  :: t00      = 273.15
#endif
#ifdef MHC
  real(r8), parameter  :: ptfactor = (rvap-rdry)/rdry ! for potential temp
  real(r8), parameter  :: prfactor = zero             ! for pressure
#else
  real(r8), parameter  :: ptfactor = rvap/rdry
  real(r8), parameter  :: prfactor = rvap/rdry
#endif
  real(r8), parameter  :: fillvalue= -1e20_r8
  real(r8), parameter  :: dbleps   = EPSILON(1._r8)

  !---------------------------------------------------
  !           Constants for CAM physics
  !---------------------------------------------------

  real(r8), parameter  :: cday     = 86400.0_r8         ! sec in calendar day (sec) 
  real(r8), parameter  :: cpliq    = 4.188e3_r8         ! specific heat of fresh h2o
  real(r8), parameter  :: cpwv     = 1.810e3_r8         ! specific heat of water vap
  real(r8), parameter  :: cappa    = rdry/cp            ! Rd/Cp
  real(r8), parameter  :: zvir     = rvap/rdry-1.       ! Rv/Rd-1
  real(r8), parameter  :: mwdry    = 28.966_r8          ! Molecular weight dry air
  real(r8), parameter  :: mwh2o    = 18.016_r8          ! Molecular weight h2o
  real(r8), parameter  :: avogad   = 6.02214e26         ! Avogadro's number (molecules/kmole)
  real(r8), parameter  :: boltz    = 1.38065e-23        ! Boltzman's constant (J/K/molecule)
  real(r8), parameter  :: latvap   = 2.501e6            ! Latent heat of vaporization (J/kg)
  real(r8), parameter  :: latice   = 3.337e5            ! Latent heat of fusion (J/kg)
  real(r8), parameter  :: karman   = 0.4_r8             ! Von Karman constant
  real(r8), parameter  :: epsilo   = mwh2o/mwdry        ! ratio of h2o to dry air molecular weights
  real(r8), parameter  :: tmelt    = 273.15_r8          ! Freezing point of water (K)
  real(r8), parameter  :: h2otrip  = 273.16_r8          ! Triple point temperature of water (K)
  real(r8), parameter  :: rhoh2o   = 1.000e3_r8         ! density of fresh water (kg/m^3)
  real(r8), parameter  :: r_universal = avogad*boltz    ! Universal gas constant (J/K/kmole)
  real(r8), parameter  :: stebol   = 5.67e-8_r8         ! Stefan-Boltzmann constant (W/m^2/K^4)

  real(r8), parameter  :: c0       = 2.99792458e8_r8    ! Speed of light in a vacuum (m/s)
  real(r8), parameter  :: planck   = 6.6260755e-34_r8   ! Planck's constant (J.s)

  ! Molecular weights
  real(r8), parameter  :: mwco2    = 44._r8             ! molecular weight co2
  real(r8), parameter  :: mwn2o    = 44._r8             ! molecular weight n2o
  real(r8), parameter  :: mwch4    = 16._r8             ! molecular weight ch4
  real(r8), parameter  :: mwf11    = 136._r8            ! molecular weight cfc11
  real(r8), parameter  :: mwf12    = 120._r8            ! molecular weight cfc12
  real(r8), parameter  :: mwo3     = 48._r8             ! molecular weight O3

end module grist_constants
