!=======================================================================
!---! This module defines a variety of physical and numerical constants
!---! used throughout the ice model
!---!
!---! author C.M. Bitz
!---!
!---! code originally based on constants.F in POP and ice_constant.F 
!---!  in CSIM
!=======================================================================

  module ice_constants

      use shr_const_mod
      use grist_constants, only: dbl_kind=>r8
!      use ice_kinds_mod
      !use pmgrid
!      use comsrf, only: plevmx

      implicit none
      save

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

!      integer, parameter :: ni = plevmx      ! number of ice levels
      real (kind=dbl_kind), parameter ::    gravit    = SHR_CONST_G  ,& ! gravitational acceleration (m/s^2)
    rhoi      = SHR_CONST_RHOICE        ,&  ! density of ice (kg/m^3)
    rhow      = SHR_CONST_RHOSW         ,&  ! density of seawater (kg/m^3)
    rhos      = 330.0_dbl_kind          ,&  ! density of snow (kg/m^3)
    rhofresh  = SHR_CONST_RHOFW         ,&  ! density of fresh water (kg/m^3)
    cp_air    = SHR_CONST_CPDAIR        ,&  ! heat capacity of air (J/kg/K)
    cpwv   =  SHR_CONST_CPWV            ,&  ! Specific heat of water vapor
    stefan_boltzmann = SHR_CONST_STEBOL ,&  !  W/m^2/K^4
    Tffresh = SHR_CONST_TKFRZ           ,&  ! Freezing temp of fresh ice (K)
    depressT = 0.054_dbl_kind           ,&  ! Tf:brine salinity ratio (C/ppt)
    cp_ice = SHR_CONST_CPICE            ,&  ! heat capacity of fresh ice (J/kg/K)
    cp_sno = 0.0_dbl_kind               ,&  ! heat capacity of snow     (J/kg/K)
    cp_ocn = SHR_CONST_CPSW             ,&  ! heat capacity of ocn    (J/kg/K)
    Lsub = SHR_CONST_LATSUB             ,&  ! latent heat, sublimation freshwater (J/kg)
    Lvap = SHR_CONST_LATVAP             ,&  ! latent heat, vaporization freshwater (J/kg)
    Lfus = SHR_CONST_LATICE             ,&  ! latent heat of fusion freshwater          (J/kg)
    Timelt = 0.0_dbl_kind               ,&  ! melting temperature of ice top surface  (C)
    Tsmelt = 0.0_dbl_kind               ,&  ! melting temperature of snow top surface (C)
      ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
    kappav = 1.4_dbl_kind  ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
    kappan = 17.6_dbl_kind ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
      ! (Briegleb JGR 97 11475-11485  July 1992)
    emissivity_ice = 0.95_dbl_kind ,&    ! emissivity of snow and ice

    kice  = 2.0340_dbl_kind  ,&          ! thermal conductivity of fresh ice(W/m/deg)
    ksno  = 0.3100_dbl_kind  ,&          ! thermal conductivity of snow  (W/m/deg)
!!   &,  ksno  = 2.0340_dbl_kind   ! thermal conductivity of snow  (W/m/deg)

    vonkar  = SHR_CONST_KARMAN  ,& ! von Karman constant
    zref    = 10._dbl_kind      ,& ! reference height for stability (m)
    iceruf  = 0.0005_dbl_kind    ! ice surface roughness (m)

   real (kind=dbl_kind), parameter ::  &
   qqqice       = 11637800._dbl_kind  ,&     ! for qsat over ice
   TTTice       = 5897.8_dbl_kind     ,&     ! for qsat over ice
   qqqocn       = 627572.4_dbl_kind   ,&     ! for qsat over ocn
   TTTocn       = 5107.4_dbl_kind          ! for qsat over ocn

   !-----------------------------------------------------------------
   ! derived and miscellaneous parameters
   !-----------------------------------------------------------------

   real (kind=dbl_kind), parameter :: rLfi = Lfus*rhoi   ,&! latent heat of fusion ice              (J/m**3)
    rLfs = Lfus*rhos   ,&! specific latent heat of fusion snow    (J/m**3)
    rLvi = Lvap*rhoi   ,&! specific latent heat of vapor*rhoice    (J/m**3)
    rLvs = Lvap*rhos   ,&! specific latent heat of vapor*rhosno    (J/m**3)
    rcpi = cp_ice*rhoi ,&! specific heat capacity of fresh ice     (J/m**3)
    rcps = cp_sno*rhos ,&! specific heat capacity of snow          (J/m**3)
    rcpidepressT = rcpi*depressT, & ! param for finding T(z) from q (J/m**3)
    rLfidepressT = rLfi*depressT ! param for heat capacity   (J deg/m**3)

   real (kind=dbl_kind), parameter :: hi_min = 0.1 ! minimum ice thickness m

   !-----------------------------------------------------------------
   ! numbers
   !-----------------------------------------------------------------

   real (kind=dbl_kind), parameter ::  c0   = 0.0_dbl_kind, &
    c1   = 1.0_dbl_kind,       & 
    c2   = 2.0_dbl_kind,       &
    c3   = 3.0_dbl_kind,       &
    c4   = 4.0_dbl_kind,       &
    c5   = 5.0_dbl_kind,       &
    c6   = 6.0_dbl_kind,       &
    c7   = 7.0_dbl_kind,       &
    c8   = 8.0_dbl_kind,       &
    c10  = 10.0_dbl_kind,      &
    c15  = 15.0_dbl_kind,      &
    c25  = 25.0_dbl_kind,      &
    c16  = 16.0_dbl_kind,      &
    c20  = 20.0_dbl_kind,      &
    c90  = 90.0_dbl_kind,      &
    c100 = 100.0_dbl_kind,     &
    c360 = 360.0_dbl_kind,     &
    c1000= 1000.0_dbl_kind,    &
    p01  = 0.01_dbl_kind,      &
    p1   = 0.1_dbl_kind,       &
    p2   = 0.2_dbl_kind,       &
    p33  = c1/c3,              &
    p66  = c2/c3,              &
    p5   = 0.5_dbl_kind,       &
    p25  = 0.25_dbl_kind,      &
    eps04  = 1.0e-4_dbl_kind,  &
    eps11  = 1.0e-11_dbl_kind, &
    eps12  = 1.0e-12_dbl_kind, &
    eps13  = 1.0e-13_dbl_kind, &
    eps15  = 1.0e-15_dbl_kind, &
    puny = eps13

   !-----------------------------------------------------------------
   ! conversion factors
   !-----------------------------------------------------------------

   real (kind=dbl_kind), parameter :: cm_to_m   = 0.01_dbl_kind, &        ! cm to meters
   m_to_cm       = 100._dbl_kind  ,&      ! meters to cm
   m2_to_km2     = 1.e-6_dbl_kind ,&      ! m^2 to km^2
   mps_to_cmpdy  = 8.64e6_dbl_kind,&      ! m per s to cm per day
   mps_to_cmpyr  = mps_to_cmpdy*365._dbl_kind ! m per s to cm per yr

    real (kind=dbl_kind), parameter :: & 
   saltmax = 3.2_dbl_kind, &   ! max ice salinity
   salnew = 4.0_dbl_kind   ! new ice salinity

    real (kind=dbl_kind), parameter :: snowpatch = 0.02_dbl_kind 

   real (kind=dbl_kind) ::  pi
!  &,  saltz(plevmx+1) ! salinity of each layer for each cat        (ppm)
!  &,  tmelz(plevmx  ) ! melting temper. of each layer for each cat   (C)

   real (kind=dbl_kind), parameter  ::    Tfrez  = -1.8_dbl_kind,&    ! freezing temperature of sea water (C)
                                          TfrezK = Tfrez+Tffresh      ! freezing temperature of sea water (K)



#ifdef USE_ICE_INIT
!c=======================================================================

      contains

!c=======================================================================
      subroutine init_constants

!---!-------------------------------------------------------------------
!---! Initializes constants that are best defined at run time 
!---!-------------------------------------------------------------------

      integer :: layer
      real(kind=dbl_kind) ::  zn

      pi  = c4*atan(c1)

      ! salinity and melting temperature profile
      do layer=1,plevmx
        zn=(layer-p5)/plevmx
        saltz(layer)=(saltmax/c2)* 
     &       (c1-cos(pi*zn**(0.407_dbl_kind/(0.573_dbl_kind+zn))))
!c     saltz(layer,nc)=saltmax ! for isosaline ice
        tmelz(layer)=-saltz(layer)*depressT
      enddo
      saltz(plevmx+1)=saltmax
      
!c     if (my_task.eq.master_task) then
!c     write (6,*) '   salt profile, cat ',nc
!c     write (6,*) (saltz(layer,nc),layer=1,nilay(nc)+1)
!c     write (6,*) '   melt temp, cat ',nc
!c     write (6,*) (tmelz(layer,nc),layer=1,nilay(nc))
!c     endif

      end subroutine init_constants
#endif

   end module ice_constants

