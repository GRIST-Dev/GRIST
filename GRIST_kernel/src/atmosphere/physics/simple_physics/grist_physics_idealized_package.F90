
!==========================================================================
! This module contains various idealized physics packages, dry or moist
! currently we use them for mimicing the behavior of full physics such that
! dtpcpl can be tested.
! Yi Zhang, 2019
!==========================================================================

 module grist_physics_idealized_package

   use grist_constants,     only: r8, i4, rdry, cp, cv, p00, one, zero
   use grist_domain_types,  only: global_domain
 
   implicit none
   private
   
   public :: grist_idealized_physics_hsdry_exp, &
             grist_idealized_physics_hsdry_imp, &
             grist_idealized_physics_dcmip2016, &
             grist_idealized_physics_mitc     , &
             grist_idealized_physics_kessler  , &
             grist_idealized_physics_kessler_klemp15  , &
             phys_tend_hsdry_exp, &    ! obsolete
             phys_tend_hsdry_imp       ! obsolete

   contains

!================================================
! dry physics: Held-Suarez, explicit
!================================================
 
  subroutine grist_idealized_physics_hsdry_exp(mesh, nlev                              , &
                                             scalar_pressure_at_pc_full_level        , &
                                             scalar_pressure_at_pc_face_level        , &
                                             scalar_U_wind_at_pc_full_level          , &
                                             scalar_V_wind_at_pc_full_level          , &
                                             scalar_potential_temp_at_pc_full_level  , &
                                             tend_U_wind_at_pc_full_level            , &
                                             tend_V_wind_at_pc_full_level            , &
                                             tend_potential_temp_at_pc_full_level)
! io   
      type(global_domain), intent(in) :: mesh
      integer(i4)        , intent(in) :: nlev
      real(r8),  intent(in)           :: scalar_pressure_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_pressure_at_pc_face_level(:,:)
      real(r8),  intent(in)           :: scalar_U_wind_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_V_wind_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_potential_temp_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_U_wind_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_V_wind_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_potential_temp_at_pc_full_level(:,:)
! local
      real(r8),  parameter            :: eta_b        = 0.7_r8
      real(r8),  parameter            :: kf           = 1._r8/86400._r8
      real(r8),  parameter            :: ka           = 1._r8/(40*86400._r8)
      real(r8),  parameter            :: ks           = 1._r8/(4*86400._r8)
      real(r8),  parameter            :: delta_t      = 60._r8
      real(r8),  parameter            :: delta_pt     = 10._r8
      real(r8),  parameter            :: zero         = 0._r8
      real(r8)                        :: kappa
      real(r8)                        :: ptemp_eq
      real(r8)                        :: kt
      real(r8)                        :: kv
      real(r8)                        :: tmpa
      real(r8)                        :: tmpb
      real(r8)                        :: eta, lat
      integer(i4)                     :: iv, ie, ilev, icell1, icell2, nlevp

          kappa     = rdry/cp
          nlevp     = nlev+1

          do iv     = 1, mesh%nv
             do ilev   = 1, nlev
! horizontal momentum
                eta    = scalar_pressure_at_pc_full_level(ilev,iv)/scalar_pressure_at_pc_face_level(nlevp,iv)
                kv     = kf*max(zero,(eta-eta_b)/(1._r8-eta_b))
                tend_U_wind_at_pc_full_level(ilev,iv) = -kv*scalar_U_wind_at_pc_full_level(ilev,iv)
                tend_V_wind_at_pc_full_level(ilev,iv) = -kv*scalar_V_wind_at_pc_full_level(ilev,iv)
! potential temperature
                lat       = mesh%vtx_lat(iv)
                kt        = ka+(ks-ka)*max(zero,(eta-eta_b)/(1._r8-eta_b))*(cos(lat)**4)
                tmpa      = 200._r8/((scalar_pressure_at_pc_full_level(ilev,iv)/p00)**kappa)
                tmpb      = 315._r8-delta_t*(sin(lat)**2)-delta_pt*log(scalar_pressure_at_pc_full_level(ilev,iv)/p00)*(cos(lat)**2)
                ptemp_eq  = max(tmpa,tmpb)
                tend_potential_temp_at_pc_full_level(ilev,iv) = -kt*(scalar_potential_temp_at_pc_full_level(ilev,iv)-ptemp_eq)
             end do
          end do

       return
  end subroutine grist_idealized_physics_hsdry_exp

!================================================
! dry physics: Held-Suarez, explicit
!================================================
 
  subroutine grist_idealized_physics_hsdry_imp(mesh, nlev, dtime                       , &
                                               scalar_pressure_at_pc_full_level        , &
                                               scalar_pressure_at_pc_face_level        , &
                                               scalar_U_wind_at_pc_full_level          , &
                                               scalar_V_wind_at_pc_full_level          , &
                                               scalar_mass_pt_at_pc_full_level         , &
                                               scalar_potential_temp_at_pc_full_level  , &
                                               tend_U_wind_at_pc_full_level            , &
                                               tend_V_wind_at_pc_full_level            , &
                                               tend_mass_pt_at_pc_full_level)
! io   
      type(global_domain), intent(in) :: mesh
      integer(i4)        , intent(in) :: nlev
      real(r8),  intent(in)           :: dtime
      real(r8),  intent(in)           :: scalar_pressure_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_pressure_at_pc_face_level(:,:)
      real(r8),  intent(in)           :: scalar_U_wind_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_V_wind_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_mass_pt_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_potential_temp_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_U_wind_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_V_wind_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_mass_pt_at_pc_full_level(:,:)
! local
      real(r8),  parameter            :: eta_b        = 0.7_r8
      real(r8),  parameter            :: kf           = 1._r8/86400._r8
      real(r8),  parameter            :: ka           = 1._r8/(40*86400._r8)
      real(r8),  parameter            :: ks           = 1._r8/(4*86400._r8)
      real(r8),  parameter            :: delta_t      = 60._r8
      real(r8),  parameter            :: delta_pt     = 10._r8
      real(r8),  parameter            :: zero         = 0._r8
      real(r8)                        :: kappa
      real(r8)                        :: temp_eq
      real(r8)                        :: kt
      real(r8)                        :: kv
      real(r8)                        :: tmpa
      real(r8)                        :: tmpb
      real(r8)                        :: eta, lat
      real(r8)                        :: han, parhan, temp
      integer(i4)                     :: iv, ie, ilev, icell1, icell2, nlevp

          kappa     = rdry/cp
          nlevp     = nlev+1

          do iv     = 1, mesh%nv
             do ilev   = 1, nlev
! horizontal momentum
                eta    = scalar_pressure_at_pc_full_level(ilev,iv)/scalar_pressure_at_pc_face_level(nlevp,iv)
                kv     = kf*max(zero,(eta-eta_b)/(1._r8-eta_b))
                tend_U_wind_at_pc_full_level(ilev,iv) = -kv*scalar_U_wind_at_pc_full_level(ilev,iv)/(one+kv*dtime)
                tend_V_wind_at_pc_full_level(ilev,iv) = -kv*scalar_V_wind_at_pc_full_level(ilev,iv)/(one+kv*dtime)
! temperature
                lat       = mesh%vtx(iv)%lat
                kt        = ka+(ks-ka)*max(zero,(eta-eta_b)/(1._r8-eta_b))*(cos(lat)**4)
                tmpa      = 200._r8
                tmpb      = (315._r8-delta_t*(sin(lat)**2)-delta_pt*log(scalar_pressure_at_pc_full_level(ilev,iv)/p00)*(cos(lat)**2))*&
                            ((scalar_pressure_at_pc_full_level(ilev,iv)/p00)**kappa)
                temp_eq   = max(tmpa,tmpb)
! UJ2012
                temp      = scalar_potential_temp_at_pc_full_level(ilev,iv)*((scalar_pressure_at_pc_full_level(ilev,iv)/p00)**kappa)
                han       = -kt/(cp/cv)*(one-temp_eq/temp)*scalar_mass_pt_at_pc_full_level(ilev,iv)
                parhan    = -kt/(cp/cv)*(one+(cp/cv-one)*temp_eq/temp)
                tend_mass_pt_at_pc_full_level(ilev,iv) = han/(one-dtime*parhan)
             end do
          end do

       return
    end subroutine grist_idealized_physics_hsdry_imp

    subroutine grist_idealized_physics_hsdry_impbad(mesh, nlev, dtime                       , &
                                               scalar_pressure_at_pc_full_level        , &
                                               scalar_pressure_at_pc_face_level        , &
                                               scalar_U_wind_at_pc_full_level          , &
                                               scalar_V_wind_at_pc_full_level          , &
                                               scalar_potential_temp_at_pc_full_level  , &
                                               tend_U_wind_at_pc_full_level            , &
                                               tend_V_wind_at_pc_full_level            , &
                                               tend_potential_temp_at_pc_full_level)
! io   
      type(global_domain), intent(in) :: mesh
      integer(i4)        , intent(in) :: nlev
      real(r8),  intent(in)           :: dtime
      real(r8),  intent(in)           :: scalar_pressure_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_pressure_at_pc_face_level(:,:)
      real(r8),  intent(in)           :: scalar_U_wind_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_V_wind_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_potential_temp_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_U_wind_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_V_wind_at_pc_full_level(:,:)
      real(r8),  intent(inout)        :: tend_potential_temp_at_pc_full_level(:,:)
! local
      real(r8),  parameter            :: eta_b        = 0.7_r8
      real(r8),  parameter            :: kf           = 1._r8/86400._r8
      real(r8),  parameter            :: ka           = 1._r8/(40*86400._r8)
      real(r8),  parameter            :: ks           = 1._r8/(4*86400._r8)
      real(r8),  parameter            :: delta_t      = 60._r8
      real(r8),  parameter            :: delta_pt     = 10._r8
      real(r8),  parameter            :: zero         = 0._r8
      real(r8)                        :: kappa
      real(r8)                        :: ptemp_eq
      real(r8)                        :: kt
      real(r8)                        :: kv
      real(r8)                        :: tmpa
      real(r8)                        :: tmpb
      real(r8)                        :: eta, lat
      integer(i4)                     :: iv, ie, ilev, icell1, icell2, nlevp

          kappa     = rdry/cp
          nlevp     = nlev+1

          do iv     = 1, mesh%nv
             do ilev   = 1, nlev
! horizontal momentum
                eta    = scalar_pressure_at_pc_full_level(ilev,iv)/scalar_pressure_at_pc_face_level(nlevp,iv)
                kv     = kf*max(zero,(eta-eta_b)/(1._r8-eta_b))
                tend_U_wind_at_pc_full_level(ilev,iv) = -kv*scalar_U_wind_at_pc_full_level(ilev,iv)/(one+kv*dtime)
                tend_V_wind_at_pc_full_level(ilev,iv) = -kv*scalar_V_wind_at_pc_full_level(ilev,iv)/(one+kv*dtime)
! potential temperature
                lat       = mesh%vtx(iv)%lat
                kt        = ka+(ks-ka)*max(zero,(eta-eta_b)/(1._r8-eta_b))*(cos(lat)**4)
                tmpa      = 200._r8/((scalar_pressure_at_pc_full_level(ilev,iv)/p00)**kappa)
                tmpb      = 315._r8-delta_t*(sin(lat)**2)-delta_pt*log(scalar_pressure_at_pc_full_level(ilev,iv)/p00)*(cos(lat)**2)
                ptemp_eq  = max(tmpa,tmpb)
                tend_potential_temp_at_pc_full_level(ilev,iv) = kt*(ptemp_eq-scalar_potential_temp_at_pc_full_level(ilev,iv))/(one+kt*dtime)
             end do
          end do

       return
   end subroutine grist_idealized_physics_hsdry_impbad

!-----------------------------------------------------------------------
!  
!  Date:  April 26, 2016 (Version 6)
!
!  Simple Physics Package ( Obtained From DCMIP2016 Experimental Protocols)
!
!  SIMPLE_PHYSICS includes large-scale precipitation, surface fluxes and
!  boundary-leyer mixing. The processes are time-split in that order.
!  A partially implicit formulation is used to foster numerical
!  stability. The routine assumes that the model levels are ordered
!  in a top-down approach, e.g. level 1 denotes the uppermost full model
!  level.
!
!  This routine is based on an implementation which was developed for
!  the NCAR Community Atmosphere Model (CAM). Adjustments for other
!  models may be necessary.
!
!  The routine provides both updates of the state variables u, v, T, q
!  (these are local copies of u,v,T,q within this physics routine) and
!  also collects their time tendencies. The latter might be used to
!  couple the physics and dynamics in a process-split way. For a
!  time-split coupling, the final state should be given to the
!  dynamical core for the next time step.
!
! Test:      0 = Reed and Jablonowski (2011) tropical cyclone test
!            1 = Moist baroclinic instability test
! RJ2012_precip:
!         true  = Turn on Reed and Jablonowski (2012) precip scheme 
!         false = Turn off Reed and Jablonowski (2012) precip scheme
! TC_PBL_mod:
!         true  = Turn on George Bryan PBL mod for tropical cyclone test
!         false = Turn off George Bryan PBL mod (i.e., run as in Reed and Jablonowski (2012))
!
!  SUBROUTINE SIMPLE_PHYSICS(pcols, pver, dtime, lat, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, test, RJ2012_precip, TC_PBL_mod)
!
!  Input variables:
!     pcols  - number of atmospheric columns (#)
!     pver   - number of model levels (#)
!     dtime  - time step (s)
!     lat    - latitude (radians)
!     t      - temperature at model levels (K)
!     q      - specific humidity at model levels (gm/gm)
!     u      - zonal wind at model levels (m/s)
!     v      - meridional wind at model levels (m/s)
!     pmid   - pressure at model levels (Pa)
!     pint   - pressure at interfaces (Pa)
!     pdel   - layer thickness (Pa)
!     rpdel  - reciprocal of layer thickness (1/Pa)
!     ps     - surface pressure (Pa)
!     test   - test case to use for sea-surface temperatures
!     RJ2012_precip - RJ2012 precip flag
!     TC_PBL_mod    - PCL modification for TC test 
!
!  Output variables:
!     Increments are added into t, q, u, v, pmid, pint, pdel, rpdel and ps
!     which are returned to the routine from which SIMPLE_PHYSICS was
!     called.  Precpitation is returned via precl.
!
!  Change log:
!  v2: removal of some NCAR CAM-specific 'use' associations
!  v3: corrected precl(i) computation, the precipitation rate is now
!      computed via a vertical integral, the previous single-level
!      computation in v2 was a bug
!  v3: corrected dtdt(i,1) computation, the term '-(i,1)' was missing
!      the temperature variable: '-t(i,1)'
!  v4: modified and enhanced parameter list to make the routine truly
!      standalone, the number of columns and vertical levels have been
!      added: pcols, pver
!  v4: 'ncol' has been removed, 'pcols' is used instead
!  v5: the sea surface temperature (SST) field Tsurf is now an array,
!      the SST now depends on the latitude
!  v5: addition of the latitude array 'lat' and the flag 'test' in the
!      parameter list
!      if test = 0: constant SST is used, correct setting for the
!                   tropical cyclone test
!      if test = 1: newly added latitude-dependent SST is used,
!                   correct setting for the moist baroclinic wave test
!                   with simple-physics
!  v6: addition of flags for a modified PBL for the TC test and
!      to turn off large-scale condensation scheme when using Kessler physics
!      Included virtual temperature in density calculation in PBL scheme
!      Also, included the virtual temperature, instead of temperature, for
!      the calculation of rho in the PBL scheme
!      (v6_1) Minor specification and generalization fixes.
!      
! Reference: Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone
!            simulations of intermediate complexity: A test case for AGCMs,
!            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
!-----------------------------------------------------------------------

  subroutine grist_idealized_physics_dcmip2016(pcols, pver, dtime, lat, t, q, u, v, pmid, pint, pdel, rpdel, ps, precl, test, &
                                               RJ2012_precip, TC_PBL_mod,dudt,dvdt,dtdt,dqdt)

  ! use physics_types     , only: physics_dme_adjust   ! This is for CESM/CAM
  ! use cam_diagnostics,    only: diag_phys_writeout   ! This is for CESM/CAM

   implicit none
!
! Input arguments - MODEL DEPENDENT
!
   integer, intent(in)  :: pcols        ! Set number of atmospheric columns
   integer, intent(in)  :: pver         ! Set number of model levels
   real(r8), intent(in) :: dtime        ! Set model physics timestep
   real(r8), intent(in) :: lat(pcols)   ! Latitude
   integer, intent(in)  :: test         ! Test number
   logical, intent(in)  :: RJ2012_precip
   logical, intent(in)  :: TC_PBL_mod

!
! Input/Output arguments
!
!  pcols is the maximum number of vertical columns per 'chunk' of atmosphere
!
   real(r8), intent(inout) :: t(pcols,pver)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(pcols,pver)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(pcols,pver)      ! Meridional wind at full-model level (m/s)
   real(r8), intent(in) :: pmid(pcols,pver)   ! Pressure is full-model level (Pa)
   real(r8), intent(in) :: pint(pcols,pver+1) ! Pressure at model interfaces (Pa)
   real(r8), intent(in) :: pdel(pcols,pver)   ! Layer thickness (Pa)
   real(r8), intent(in) :: rpdel(pcols,pver)  ! Reciprocal of layer thickness (1/Pa)
   real(r8), intent(in) :: ps(pcols)          ! Surface Pressue (Pa)

!
! Output arguments
!
   real(r8), intent(out) :: precl(pcols)         ! Precipitation rate (m_water / s)
! Physics Tendency Arrays
   real(r8), intent(out) :: dudt(pcols,pver)     ! Zonal wind tendency
   real(r8), intent(out) :: dvdt(pcols,pver)     ! Meridional wind tendency
   real(r8), intent(out) :: dtdt(pcols,pver)     ! Temperature tendency
   real(r8), intent(out) :: dqdt(pcols,pver)     ! Specific humidity tendency

!
!---------------------------Local workspace-----------------------------
!

! Integers for loops

   integer  i,k                         ! Longitude, level indices

! Physical Constants - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
   real(r8) a                           ! Reference Earth's Radius (m)
   real(r8) omega                       ! Reference rotation rate of the Earth (s^-1)
   real(r8) pi                          ! pi

! Simple Physics Specific Constants 

!++++++++
   real(r8) Tsurf(pcols)                ! Sea Surface Temperature (constant for tropical cyclone)
!++++++++                                 Tsurf needs to be dependent on latitude for the
                                        ! moist baroclinic wave test, adjust

   real(r8) SST_TC                      ! Sea Surface Temperature for tropical cyclone test
   real(r8) T0                          ! Control temp for calculation of qsat
   real(r8) e0                          ! Saturation vapor pressure at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation
   real(r8) T00                         ! Horizontal mean T at surface for moist baro test
   real(r8) u0                          ! Zonal wind constant for moist baro test
   real(r8) latw                        ! halfwidth for  for baro test
   real(r8) eta0                        ! Center of jets (hybrid) for baro test
   real(r8) etav                        ! Auxiliary variable for baro test
   real(r8) q0                          ! Maximum specific humidity for baro test
   real(r8) kappa                       ! von Karman constant


! Temporary variables for tendency calculations

   real(r8) tmp                         ! Temporary
   real(r8) qsat                        ! Saturation vapor pressure
   real(r8) qsats                       ! Saturation vapor pressure of SST

! Variables for Boundary Layer Calculation

   real(r8) wind(pcols)                 ! Magnitude of Wind
   real(r8) Cd(pcols)                   ! Drag coefficient for momentum
   real(r8) Km(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations 
   real(r8) Ke(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations
   real(r8) rho                         ! Density at lower/upper interface
   real(r8) za(pcols)                   ! Heights at midpoints of first model level
   real(r8) zi(pcols,pver+1)            ! Heights at model interfaces
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Top of boundary layer
   real(r8) zpbltop                     ! Top of boundary layer for George Bryan Modifcation
   real(r8) pblconst                    ! Constant for the calculation of the decay of diffusivity 
   real(r8) CA(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pcols,pver+1)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme


! Variable for Dry Mass Adjustment, this dry air adjustment is necessary to
! conserve the mass of the dry air

   real(r8) qini(pcols,pver)            ! Initial specific humidity

!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT
!
!===============================================================================
   gravit = 9.80616_r8                   ! Gravity (9.80616 m/s^2)
   rair   = 287.0_r8                     ! Gas constant for dry air: 287 J/(kg K)
   cpair  = 1.0045e3_r8                  ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   latvap = 2.5e6_r8                     ! Latent heat of vaporization (J/kg)
   rh2o   = 461.5_r8                     ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = 6371220.0_r8                 ! Reference Earth's Radius (m)
   omega  = 7.29212d-5                   ! Reference rotation rate of the Earth (s^-1)
   pi     = 4._r8*atan(1._r8)            ! pi

!===============================================================================
!
! Local Constants for Simple Physics
!
!===============================================================================
      C        = 0.0011_r8      ! From Simth and Vogl 2008
      SST_TC   = 302.15_r8      ! Constant Value for SST
      T0       = 273.16_r8      ! control temp for calculation of qsat
      e0       = 610.78_r8      ! saturation vapor pressure at T0 for calculation of qsat
      rhow     = 1000.0_r8      ! Density of Liquid Water 
      Cd0      = 0.0007_r8      ! Constant for Cd calc. Simth and Vogl 2008
      Cd1      = 0.000065_r8    ! Constant for Cd calc. Simth and Vogl 2008
      Cm       = 0.002_r8       ! Constant for Cd calc. Simth and Vogl 2008
      v20      = 20.0_r8        ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
      p0       = 100000.0_r8    ! Constant for potential temp calculation
      pbltop   = 85000._r8      ! Top of boundary layer in p
      zpbltop  = 1000._r8       ! Top of boundary layer in z
      pblconst = 10000._r8      ! Constant for the calculation of the decay of diffusivity
      T00      = 288.0_r8         ! Horizontal mean T at surface for moist baro test
      u0       = 35.0_r8          ! Zonal wind constant for moist baro test
      latw     = 2.0_r8*pi/9.0_r8 ! Halfwidth for  for baro test
      eta0     = 0.252_r8         ! Center of jets (hybrid) for baro test
      etav     = (1._r8-eta0)*0.5_r8*pi ! Auxiliary variable for baro test
      q0       = 0.021_r8         ! Maximum specific humidity for baro test
      kappa    = 0.4_r8         ! von Karman constant

!===============================================================================
!
! Definition of local arrays
!
!===============================================================================
!
! Calculate hydrostatic height
!
     do i=1,pcols
        dlnpint = log(ps(i)) - log(pint(i,pver))  ! ps(i) is identical to pint(i,pver+1), note: this is the correct sign (corrects typo in JAMES paper) 
        za(i) = rair/gravit*t(i,pver)*(1._r8+zvir*q(i,pver))*0.5_r8*dlnpint
        zi(i,pver+1) = 0.0_r8
     end do
!
! Set Initial Specific Humidity
!
     qini(:pcols,:pver) = q(:pcols,:pver)
!
! Set Sea Surface Temperature (constant for tropical cyclone)
! Tsurf needs to be dependent on latitude for moist baroclinic wave test
! Tsurf needs to be constant for tropical cyclone test
!
     if (test .eq. 1) then ! Moist Baroclinic Wave Test
        do i=1,pcols
           Tsurf(i) = (T00 + pi*u0/rair * 1.5_r8 * sin(etav) * (cos(etav))**0.5_r8 *                 &
                     ((-2._r8*(sin(lat(i)))**6 * ((cos(lat(i)))**2 + 1._r8/3._r8) + 10._r8/63._r8)* &
                     u0 * (cos(etav))**1.5_r8  +                                                    &
                     (8._r8/5._r8*(cos(lat(i)))**3 * ((sin(lat(i)))**2 + 2._r8/3._r8) - pi/4._r8)*a*omega*0.5_r8 ))/ &
                     (1._r8+zvir*q0*exp(-(lat(i)/latw)**4))

        end do
     else if (test .eq. 0) then ! Tropical Cyclone Test
        do i=1,pcols
           Tsurf(i) = SST_TC
        end do
     end if

!===============================================================================
!
! Set initial physics time tendencies and precipitation field to zero
!
!===============================================================================
     dtdt(:pcols,:pver)  = 0._r8            ! initialize temperature tendency with zero
     dqdt(:pcols,:pver)  = 0._r8            ! initialize specific humidity tendency with zero
     dudt(:pcols,:pver)  = 0._r8            ! initialize zonal wind tendency with zero
     dvdt(:pcols,:pver)  = 0._r8            ! initialize meridional wind tendency with zero
     precl(:pcols) = 0._r8                  ! initialize precipitation rate with zero

!===============================================================================
!
! Large-Scale Condensation and Precipitation
!
!===============================================================================

      if (RJ2012_precip) then
!
! Calculate Tendencies
!
      do k=1,pver
         do i=1,pcols
            qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/t(i,k))-1._r8/T0))
            if (q(i,k) > qsat) then
               tmp  = 1._r8/dtime*(q(i,k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(i,k)**2)))
               dtdt(i,k) = dtdt(i,k)+latvap/cpair*tmp
               dqdt(i,k) = dqdt(i,k)-tmp
               precl(i)  = precl(i)+tmp*pdel(i,k)/(gravit*rhow)
            end if
         end do
      end do
!
! Update moisture and temperature fields from Larger-Scale Precipitation Scheme
!
      do k=1,pver
         do i=1,pcols
            t(i,k) =  t(i,k) + dtdt(i,k)*dtime
            q(i,k) =  q(i,k) + dqdt(i,k)*dtime
         end do
      end do

!===============================================================================
! Send variables to history file - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: The variables, as done in many GCMs, are written to the history file
!       after the moist physics process.  This ensures that the moisture fields
!       are somewhat in equilibrium.
!===============================================================================
  !  call diag_phys_writeout(state)   ! This is for CESM/CAM
    
      end if
     
!===============================================================================
!
! Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
! sensible heat and latent heat
!
! We are using Simplified Ekman theory to compute the diffusion coefficients
! Kx for the boundary-layer mixing. The Kx values are calculated at each time step
! and in each column.
!
!===============================================================================!
! Compute magnitude of the wind and drag coeffcients for turbulence scheme
!
     do i=1,pcols
        wind(i) = sqrt(u(i,pver)**2+v(i,pver)**2)
     end do
     do i=1,pcols
        if( wind(i) .lt. v20) then
           Cd(i) = Cd0+Cd1*wind(i) 
        else
           Cd(i) = Cm
        endif
     end do

     if (TC_PBL_mod) then !Bryan TC PBL Modification 
     do k=pver,1,-1
        do i=1,pcols
           dlnpint = log(pint(i,k+1)) - log(pint(i,k))
           zi(i,k) = zi(i,k+1)+rair/gravit*t(i,k)*(1._r8+zvir*q(i,k))*dlnpint
           if( zi(i,k) .le. zpbltop) then
              Km(i,k) = kappa*sqrt(Cd(i))*wind(i)*zi(i,k)*(1._r8-zi(i,k)/zpbltop)*(1._r8-zi(i,k)/zpbltop)
              Ke(i,k) = kappa*sqrt(C)*wind(i)*zi(i,k)*(1._r8-zi(i,k)/zpbltop)*(1._r8-zi(i,k)/zpbltop) 
           else
              Km(i,k) = 0.0_r8
              Ke(i,k) = 0.0_r8
           end if 
        end do
     end do     
     else ! Reed and Jablonowski (2012) Configuration
     do k=1,pver
        do i=1,pcols
           if( pint(i,k) .ge. pbltop) then
              Km(i,k) = Cd(i)*wind(i)*za(i) 
              Ke(i,k) = C*wind(i)*za(i)
           else
              Km(i,k) = Cd(i)*wind(i)*za(i)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)
              Ke(i,k) = C*wind(i)*za(i)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)
           end if
        end do
     end do
     end if
!===============================================================================
! Update the state variables u, v, t, q with the surface fluxes at the
! lowest model level, this is done with an implicit approach
! see Reed and Jablonowski (JAMES, 2012)
!
! Sea Surface Temperature Tsurf is constant for tropical cyclone test 5-1
! Tsurf needs to be dependent on latitude for the
! moist baroclinic wave test 
!===============================================================================
#ifndef NO_TSURF
     do i=1,pcols
        qsats = epsilo*e0/ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))
        dudt(i,pver) = dudt(i,pver) + (u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))-u(i,pver))/dtime
        dvdt(i,pver) = dvdt(i,pver) + (v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))-v(i,pver))/dtime
        u(i,pver) = u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        v(i,pver) = v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))
        dtdt(i,pver) = dtdt(i,pver) +((t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))-t(i,pver))/dtime 
        t(i,pver) = (t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))  
        dqdt(i,pver) = dqdt(i,pver) +((q(i,pver)+C*wind(i)*qsats*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))-q(i,pver))/dtime
        q(i,pver) = (q(i,pver)+C*wind(i)*qsats*dtime/za(i))/(1._r8+C*wind(i)*dtime/za(i))
     end do
#endif

!===============================================================================
! Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
!===============================================================================
! Calculate Diagonal Variables for Implicit PBL Scheme
!
      do k=1,pver-1
         do i=1,pcols
            rho = (pint(i,k+1)/(rair*(t(i,k+1)*(1._r8+zvir*q(i,k+1))+t(i,k)*(1._r8+zvir*q(i,k)))/2.0_r8)) 
            CAm(i,k) = rpdel(i,k)*dtime*gravit*gravit*Km(i,k+1)*rho*rho/(pmid(i,k+1)-pmid(i,k))    
            CCm(i,k+1) = rpdel(i,k+1)*dtime*gravit*gravit*Km(i,k+1)*rho*rho/(pmid(i,k+1)-pmid(i,k))
            CA(i,k) = rpdel(i,k)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho/(pmid(i,k+1)-pmid(i,k))
            CC(i,k+1) = rpdel(i,k+1)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho/(pmid(i,k+1)-pmid(i,k))
         end do
      end do
      do i=1,pcols
         CAm(i,pver) = 0._r8
         CCm(i,1) = 0._r8
         CEm(i,pver+1) = 0._r8
         CA(i,pver) = 0._r8
         CC(i,1) = 0._r8
         CE(i,pver+1) = 0._r8
         CFu(i,pver+1) = 0._r8
         CFv(i,pver+1) = 0._r8
         CFt(i,pver+1) = 0._r8
         CFq(i,pver+1) = 0._r8 
      end do
      do i=1,pcols
         do k=pver,1,-1
            CE(i,k) = CC(i,k)/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFu(i,k) = (u(i,k)+CAm(i,k)*CFu(i,k+1))/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFv(i,k) = (v(i,k)+CAm(i,k)*CFv(i,k+1))/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFt(i,k) = ((p0/pmid(i,k))**(rair/cpair)*t(i,k)+CA(i,k)*CFt(i,k+1))/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CFq(i,k) = (q(i,k)+CA(i,k)*CFq(i,k+1))/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
       end do
      end do
!
! Calculate the updated temperaure and specific humidity and wind tendencies
!
! First we need to calculate the tendencies at the top model level
!
      do i=1,pcols
            dudt(i,1)  = dudt(i,1)+(CFu(i,1)-u(i,1))/dtime
            dvdt(i,1)  = dvdt(i,1)+(CFv(i,1)-v(i,1))/dtime
            u(i,1)    = CFu(i,1)
            v(i,1)    = CFv(i,1)
            dtdt(i,1)  = dtdt(i,1)+(CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)-t(i,1))/dtime
            t(i,1)    = CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)
            dqdt(i,1)  = dqdt(i,1)+(CFq(i,1)-q(i,1))/dtime
            q(i,1)  = CFq(i,1)
      end do

      do i=1,pcols
         do k=2,pver
            dudt(i,k)  = dudt(i,k)+(CEm(i,k)*u(i,k-1)+CFu(i,k)-u(i,k))/dtime
            dvdt(i,k)  = dvdt(i,k)+(CEm(i,k)*v(i,k-1)+CFv(i,k)-v(i,k))/dtime
            u(i,k)    = CEm(i,k)*u(i,k-1)+CFu(i,k) 
            v(i,k)    = CEm(i,k)*v(i,k-1)+CFv(i,k)
            dtdt(i,k)  = dtdt(i,k)+((CE(i,k)*t(i,k-1)*(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k))*(pmid(i,k)/p0)**(rair/cpair)-t(i,k))/dtime 
            t(i,k)    = (CE(i,k)*t(i,k-1)*(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k))*(pmid(i,k)/p0)**(rair/cpair)
            dqdt(i,k)  = dqdt(i,k)+(CE(i,k)*q(i,k-1)+CFq(i,k)-q(i,k))/dtime
            q(i,k)  = CE(i,k)*q(i,k-1)+CFq(i,k)
         end do
      end do

!===============================================================================
!
! Dry Mass Adjustment - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: Care needs to be taken to ensure that the model conserves the total
!       dry air mass. Add your own routine here.
!===============================================================================
  !  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM

   return
end subroutine grist_idealized_physics_dcmip2016

subroutine grist_idealized_physics_mitc(pcols, pver, dtime, lat, t, q, u, v, pmid, pint, etamid, ps, precl, &
                                        dtdt, dqdt,dudt,dvdt)
!----------------------------------------------------------------------- 
! 
! Purpose: Moist Idealized Test Case (MITC) Physics Package
!          This is a moist variant of the Held-Suarez test
!
! Author: D. R. Thatcher (University of Michigan, dtatch@umich.edu)
!         July/15/2015
!
! Based on the Simple-Physics Package by K. A. Reed
!              (University of Michigan, now Stony Brook University,
!               email: kevin.a.reed@stonybrook.edu)
!               and the Held-Suarez test for dry dynamical cores
!
! 
! Description: Includes large-scale precipitation, surface fluxes & boundary-leyer mixing
!              of heat, moisture and momentum, and radiation based on the Simple-Physics
!              Package (Reed and Jablonowski 2012) and the Held-Suarez test.
!              The processes are time-split in that order. A partially-implicit
!              formulation is used to foster numerical stability.
!              The routine assumes that the model levels are ordered
!              in a top-down approach, e.g. level 1 denotes the uppermost
!              full model level, the level pver is the lowest model level.
!
!              Sea surface temperature is defined as a Gaussian profile
!              dependent on latitude. Idealized radiation and dissipation of
!              low level wind is based on the Held-Suarez test for dry 
!              dynamical cores.
!
!              This routine is based on an implementation which was
!              developed for the NCAR Community Atmosphere Model (CAM) which uses
!              a hybrid pressure-based vertical coordinate eta.
!              Adjustments for other models will be necessary.
!
!              The routine provides both updates of the state variables
!              u, v, T, q (these are local copies of u,v,T,q within this physics
!              routine) and also collects their time tendencies.
!              The latter might be used to couple the physics and dynamics
!              in a process-split way. They can also be used for diagnostic
!              purposes to monitor the physics forcing. For a time-split coupling, the
!              final state should be given to the dynamical core for the next time step.
!
!
! Reference: Thatcher, D. R. and C. Jablonowski (2015), A moist aquaplanet variant 
!            of the Held-Suarez test for atmospheric model dynamical cores, 
!            Geoscientic Model Development Discussions
!
!            Reed, K. A. and C. Jablonowski (2012), Idealized tropical cyclone
!            simulations of intermediate complexity: A test case for AGCMs, 
!            J. Adv. Model. Earth Syst., Vol. 4, M04001, doi:10.1029/2011MS000099
!
!            Held, I. M., and M. J. Suarez (1994), A proposal for the
!            intercomparison of the dynamical cores of atmospheric general
!            circulation models, Bulletin of the Amer. Meteor. Soc., Vol. 75,
!            pp. 1825-1830, doi:10.1175/1520-0477(1994)075<1825:APFTIO>2.0.CO;2
!-----------------------------------------------------------------------

! This is for NCAR's CESM/CAM model ------------------------------------
!   use shr_kind_mod, only: r8 => shr_kind_r8
!   use pmgrid            , only: plev,plat,plevp
!   use ppgrid
!   use phys_grid         , only: get_lat_all_p, get_rlat_all_p
!   use physics_types     , only: physics_state, physics_tend, physics_dme_adjust, set_dry_to_wet
!   use geopotential      , only: geopotential_t
!   use cam_history,        only: outfld
!   use physconst,          only: gravit, rair, cpair, latvap, rh2o, epsilo, karman, zvir, pi
!   use camsrfexch,         only: cam_out_t
!   use cam_diagnostics,    only: diag_phys_writeout
!   use hycoef,             only: ps0, etamid
!   use dycore,             only: dycore_is
!   use perf_mod
! ----------------------------------------------------------------------

   implicit none
!
! Input arguments - MODEL DEPENDENT
!
   integer, intent(in)  :: pcols        ! Number of atmospheric columns
   integer, intent(in)  :: pver         ! Number of full model levels
   real(r8), intent(in) :: dtime        ! Physics timestep in seconds, if a leapfrong time-stepping scheme is
                                        ! used dtime needs to be initialized with the doubled physics time step
   real(r8), intent(in) :: lat(pcols)   ! Latitude of all atmospheric columns in radians
   real(r8), intent(in) :: etamid(pver) ! Hybrid vertical coordinate midpoint (lies between 0-1), 
                                        ! can be replaced by sigma=pmid/ps in models that do not
                                        ! use the hybrid pressure-based coordinate system in the vertical

!
! Input/Output arguments 
!
!  pcols is the maximum number of vertical columns per 'chunk' of atmosphere
!
   real(r8), intent(inout) :: t(pcols,pver)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(pcols,pver)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(pcols,pver)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(pcols,pver)      ! Meridional wind at full-model level (m/s)
   real(r8), intent(in)    :: pmid(pcols,pver)   ! Pressure at full-model level (Pa)
   real(r8), intent(in)    :: pint(pcols,pver+1) ! Pressure at model interfaces (Pa), position pver+1 is the surface
   real(r8), intent(in)    :: ps(pcols)          ! Surface Pressue (Pa)

!
! Output arguments 
!
   real(r8), intent(out) :: precl(pcols)         ! Large-scale precipitation rate (m_water / s)
   real(r8), intent(out) :: dtdt(pcols,pver)     ! Temperature at full-model level (K)
   real(r8), intent(out) :: dqdt(pcols,pver)     ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(out) :: dudt(pcols,pver)     ! Zonal wind at full-model level (m/s)
   real(r8), intent(out) :: dvdt(pcols,pver)     ! Meridional wind at full-model level (m/s)

!
!---------------------------Local workspace-----------------------------
!

! Integers for loops

   integer  i,k                         ! Longitude and level indices

! Physical Constants, all in SI units - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air
   real(r8) latvap                      ! Latent heat of vaporization
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
   real(r8) a                           ! Reference Earth's Radius (m)
   real(r8) omega                       ! Reference rotation rate of the Earth (s^-1)
   real(r8) pi                          ! pi

! Moist Idealized Physics Specific Constants

   real(r8) Tsurf(pcols)                ! Sea Surface Temperature (K)
   real(r8) T_min                       ! Minimum sea surface temperature (K)
   real(r8) del_T                       ! difference in sea surface temperature from equator to pole (K)
   real(r8) T_width                     ! Width parameter for sea surface temperature profile
   real(r8) T0                          ! Control temp for calculation of qsat (K),  triple point of water
   real(r8) e0                          ! Saturation vapor pressure (Pa) at T0 for calculation of qsat
   real(r8) rhow                        ! Density of Liquid Water (kg/m3)
   real(r8) p0                          ! Constant for calculation of potential temperature
   real(r8) Cd0                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cd1                         ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) Cm                          ! Constant for calculating Cd from Smith and Vogl 2008
   real(r8) v20                         ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   real(r8) C                           ! Drag coefficient for sensible heat and evaporation

! Constants and variables for the modified Held Suarez forcing

   real(r8) sec_per_day                 ! Number of seconds per day
   real(r8) kf                          ! 1./efolding_time for wind dissipation
   real(r8) ka                          ! 1./efolding_time for temperature diss.
   real(r8) ks                          ! 1./efolding_time for temperature diss.
   real(r8) sigmab                      ! threshold sigma level (PBL level)
   real(r8) onemsig                     ! 1. - sigma_reference
   real(r8) ps0                         ! Base state surface pressure (Pa)
   real(r8) t00                         ! minimum reference temperature (K)
   real(r8) t_max                       ! modified maximum HS equilibrium temperature (original is 315 K)
   real(r8) delta_T                     ! difference in eq-polar HS equilibrium temperature (K)
   real(r8) delta_theta                 ! parameter for vertical temperature gradient (K)
   real(r8) kv                          ! 1./efolding_time (normalized) for wind
   real(r8) kt                          ! 1./efolding_time for temperature diss.
   real(r8) trefa                       ! "radiative equilibrium" T
   real(r8) trefc                       ! used in calc of "radiative equilibrium" T

! Physics Tendency Arrays
! move to output
!   real(r8) dtdt(pcols,pver)            ! Temperature tendency
!   real(r8) dqdt(pcols,pver)            ! Specific humidity tendency
!   real(r8) dudt(pcols,pver)            ! Zonal wind tendency
!   real(r8) dvdt(pcols,pver)            ! Meridional wind tendency

! Surface fluxes for diagnostic purposes

   real(r8) taux(pcols)                 ! surface momentum flux in the longitudinal direction
   real(r8) tauy(pcols)                 ! surface momentum flux in the latitudinal direction
   real(r8) shflx(pcols)                ! sensible heat flux at the surface
   real(r8) lhflx(pcols)                ! latent heat flux at the surface

! Temporary variables and arrays

   real(r8) tmp                         ! holds coeffiecients
   real(r8) qsat                        ! Saturation vapor pressure (Pa)
   real(r8) qsats                       ! Saturation vapor pressure (Pa) at the surface with prescribed SST
   real(r8) pdel(pcols,pver)            ! Pressure thickness (Pa) of the layer, pdel is the positive
                                        ! difference between the surrounding model interface pressures
   real(r8) rpdel(pcols,pver)           ! Reciprocal of the pressure layer thickness (1/Pa)


! Variables for Boundary Layer Calculation

   real(r8) wind(pcols)                 ! Magnitude of Wind
   real(r8) Cd(pcols)                   ! Drag coefficient for momentum
   real(r8) Km(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations: u and v
   real(r8) Ke(pcols,pver+1)            ! Eddy diffusivity for boundary layer calculations: T and q
   real(r8) rho                         ! Moist density
   real(r8) za(pcols)                   ! Heights at midpoints of first model level (m)
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Pressure position (Pa) of the top of boundary layer
   real(r8) pblconst                    ! Pressure constant (Pa) for the calculation of the decay of diffusivity
   real(r8) CA(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pcols,pver)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pcols,pver+1)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pcols,pver)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pcols,pver+1)           ! Matrix Coefficents for PBL Scheme


! Save the initial specific humidity contents to be able to compute the moist pressure
! adjustments at the end of the routine, whether this is needed depends on the model

   real(r8) qini(pcols,pver)            ! Initial specific humidity

!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT 
!
!===============================================================================
   gravit = 9.80616_r8                   ! Gravity (9.80616 m/s^2)
   rair   = 287.0_r8                     ! Gas constant for dry air: 287 J/(kg K)
   cpair  = 1.0045e3_r8                  ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   latvap = 2.5e6_r8                     ! Latent heat of vaporization (J/kg)
   rh2o   = 461.5_r8                     ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = 6371220.0_r8                 ! Reference Earth's Radius (m)
   omega  = 7.29212d-5                   ! Reference rotation rate of the Earth (s^-1)
   pi     = 4._r8*atan(1._r8)            ! pi

!===============================================================================
!
! Local Constants for Moist Idealized Physics
!
!===============================================================================
   T_min    = 271._r8             ! Minimum sea surface temperature (K)
   del_T    = 29._r8              ! SST difference between equator and poles (K)
   T_width  = 26.0_r8*pi/180.0_r8 ! width parameter for sea surface temperature
   C        = 0.0044_r8           ! Increased by a factor of 4 in comaprison to Smith and Vogl 2008
   T0       = 273.16_r8           ! control temp (K) for calculation of qsat, triple point of water
   e0       = 610.78_r8           ! saturation vapor pressure (Pa) at T0 for calculation of qsat
   rhow     = 1000.0_r8           ! Density of Liquid Water (kg/m3)
   Cd0      = 0.0007_r8           ! Constant for Cd calc. Smith and Vogl 2008
   Cd1      = 0.000065_r8         ! Constant for Cd calc. Smith and Vogl 2008
   Cm       = 0.002_r8            ! Constant for Cd calc. Smith and Vogl 2008
   v20      = 20.0_r8             ! Threshold wind speed for calculating Cd from Smith and Vogl 2008
   p0       = 100000.0_r8         ! Pressure constant (Pa) for potential temp calculation
   pbltop   = 85000._r8           ! Pressure position (Pa) of the top of the boundary layer
   pblconst = 10000._r8           ! Pressure constant (Pa) for the calculation of the decay of diffusivity

   sec_per_day = 86400._r8                  ! Number of seconds per day
   kf          = 1._r8/( 1._r8*sec_per_day) ! 1./efolding_time for wind dissipation
   ka          = 1._r8/(40._r8*sec_per_day) ! 1./efolding_time for temperature diss.
   ks          = 1._r8/( 4._r8*sec_per_day) ! 1./efolding_time for temperature diss.
   sigmab      = 0.7_r8                     ! threshold sigma level (PBL level)
   onemsig     = 1._r8-sigmab               ! 1. - sigma_reference
   ps0         = 100000._r8                 ! Base state surface pressure (Pa)
   t00         = 200._r8                    ! minimum reference temperature (K)
#ifdef MITCHS94
   t_max       = 315._r8                    ! modified maximum HS equilibrium temperature (original is 315 K)
   delta_T     = 60._r8                     ! difference in eq-polar HS equilibrium temperature (K)
#else
   t_max       = 294._r8                    ! modified maximum HS equilibrium temperature (original is 315 K)
   delta_T     = 65._r8                     ! difference in eq-polar HS equilibrium temperature (K)
#endif
   delta_theta = 10._r8                     ! parameter for vertical temperature gradient (K)

!===============================================================================
!
! Definition of local arrays
!
!===============================================================================
!
! Calculate hydrostatic height za of the lowest model level
!
     do i=1,pcols 
        dlnpint = log(ps(i)) - log(pint(i,pver))  ! ps(i) is identical to pint(i,pver+1), 
                                                  ! note: this gives the correct positive sign and corrects
                                                  ! the typo in the Reed and Jablonowski (2012) paper
        za(i) = rair/gravit*t(i,pver)*(1._r8+zvir*q(i,pver))*0.5_r8*dlnpint
     end do
!
! Set Initial Specific Humidity - For dry mass adjustment at the end
!
     qini(:pcols,:pver) = q(:pcols,:pver)

!--------------------------------------------------------------
! Set Sea Surface Temperature
!--------------------------------------------------------------
     do i=1,pcols
        Tsurf(i) = del_T*exp(-(((lat(i))**2.0_r8)/(2.0_r8*(T_width**2.0_r8)))) + T_min
     end do

!===============================================================================
!
! Set initial physics time tendencies and precipitation field to zero
!
!===============================================================================
     dtdt(:pcols,:pver)  = 0._r8            ! initialize temperature tendency with zero
     dqdt(:pcols,:pver)  = 0._r8            ! initialize specific humidity tendency with zero
     dudt(:pcols,:pver)  = 0._r8            ! initialize zonal wind tendency with zero
     dvdt(:pcols,:pver)  = 0._r8            ! initialize meridional wind tendency with zero
     taux(:pcols)        = 0._r8            ! surface momentum flux in the longitudinal direction
     tauy(:pcols)        = 0._r8            ! surface momentum flux in the latitudinal direction
     shflx(:pcols)       = 0._r8            ! sensible heat flux at the surface
     lhflx(:pcols)       = 0._r8            ! latent heat flux at the surface
     precl(:pcols)       = 0._r8            ! initialize precipitation rate with zero

!===============================================================================
!
! Large-Scale Condensation and Large-Scale Precipitation Rate
!
!===============================================================================
!
! Calculate Tendencies
!
      do k=1,pver
         do i=1,pcols
!------------
!           Initialize the pressure thickness of the layer and its reciprocal
!------------
            pdel(i,k)  = abs(pint(i,k+1)-pint(i,k))                                 ! Pressure thickness (Pa)
            rpdel(i,k) = 1._r8 / pdel(i,k)                                          ! reciprocal of the pressure thickness
!------------
!           Calculate the moisture removal
!------------
            qsat = epsilo*e0/pmid(i,k)*exp(-latvap/rh2o*((1._r8/t(i,k))-1._r8/T0))  ! saturation specific humidity
            if (q(i,k) > qsat) then                                                 ! saturated?
               tmp  = 1._r8/dtime*(q(i,k)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(i,k)**2))) ! condensation rate
               dtdt(i,k) = dtdt(i,k)+latvap/cpair*tmp                               ! temperature tendency
               dqdt(i,k) = dqdt(i,k)-tmp                                            ! specific humidity tendency
               precl(i) = precl(i) + tmp*pdel(i,k)/(gravit*rhow)                    ! precipitation rate
            end if
         end do
      end do
!
! Update moisture and temperature fields from Large-Scale Precipitation Scheme
!
      do k=1,pver
         do i=1,pcols
            t(i,k) =  t(i,k) + dtdt(i,k)*dtime    ! update the state variable T
            q(i,k) =  q(i,k) + dqdt(i,k)*dtime    ! update the state variable q
         end do
      end do

!===============================================================================
! Send the state variables u, v, t, q, ps to an output file - THIS PROCESS WILL BE MODEL SPECIFIC
!
! note: The variables, as done in many GCMs, are written to the output file
!       after the moist physics adjustments.  This ensures that the moisture fields
!       are somewhat in equilibrium (e.g. no supersaturation in archived data).
!===============================================================================
  !  call diag_phys_writeout(state)   ! This is for CESM/CAM

!===============================================================================
!
! Turbulent mixing coefficients for the PBL mixing of horizontal momentum,
! sensible heat and latent heat
!
! The turbulent mixing of horizontal momentum as defined in Red and Jablonowski (2012)
! is not applied in the moist Held-Suarez approach by Thatcher and Jablonowski (2015),
! and replaced by the Held-Suarez Rayleigh friction. Nevertheless, we leave the code
! for the Km coefficients in here since it could easily be used for other applications.
!
! We are using Simplified Ekman theory to compute the diffusion coefficients 
! Kx for the boundary-layer mixing. The Kx values are calculated at each time step
! and in each column.
!
!===============================================================================
!
! Compute magnitude of the wind and drag coeffcients for turbulence scheme:
! they depend on the conditions at the lowest model level and stay constant
! up to the 850 hPa level. Above this level the coefficients are decreased
! and tapered to zero. At the 700 hPa level the strength of the K coefficients
! is about 10% of the maximum strength. 
!
     do i=1,pcols
        wind(i) = sqrt(u(i,pver)**2+v(i,pver)**2)    ! wind magnitude at the lowest level
     end do
     do i=1,pcols
        Ke(i,pver+1) = C*wind(i)*za(i)               ! Ke: latent and sensible heat mixing coefficient
        if( wind(i) .lt. v20) then                   ! Km: momentum mixing coefficient
           Cd(i) = Cd0+Cd1*wind(i) 
           Km(i,pver+1) = Cd(i)*wind(i)*za(i)
        else
           Cd(i) = Cm
           Km(i,pver+1) = Cm*wind(i)*za(i)
        endif
     end do

      do k=1,pver                                     ! adjustments in the vertical direction
         do i=1,pcols
            if( pint(i,k) .ge. pbltop) then
               Km(i,k) = Km(i,pver+1)                 ! constant Km below 850 hPa level
               Ke(i,k) = Ke(i,pver+1)                 ! constant Ke below 850 hPa level
            else
               Km(i,k) = Km(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Km tapered to 0
               Ke(i,k) = Ke(i,pver+1)*exp(-(pbltop-pint(i,k))**2/(pblconst)**2)  ! Ke tapered to 0
            end if
         end do
      end do     


!===============================================================================
! Update the state variables t, q with the surface fluxes at the
! lowest model level, this is done with an implicit approach
! see Reed and Jablonowski (2012)
!
! Note that Reed and Jablonowski (2012) also update u and v with the momentum fluxes
! at the surface which is commented out here. The code is left since it could be
! useful for other applications.
!===============================================================================

#ifndef NO_TSURF
     do i=1,pcols
        rho   = pmid(i,pver)/(rair * t(i,pver) * (1._r8+zvir*q(i,pver)))      ! moist density at lowest model level rho = p/(Rd Tv)
        qsats = epsilo*e0/ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))      ! saturation specific humidity at the surface
!---------
!       Surface momentum fluxes are replaced by Rayleigh friction, commented out
!
!       dudt(i,pver) = dudt(i,pver) + (u(i,pver) &                                ! replaced by Rayleigh friction below
!                           /(1._r8+Cd(i)*wind(i)*dtime/za(i))-u(i,pver))/dtime
!       dvdt(i,pver) = dvdt(i,pver) + (v(i,pver) &                                ! replaced by Rayleigh friction below
!                          /(1._r8+Cd(i)*wind(i)*dtime/za(i))-v(i,pver))/dtime
!       taux(i)     = -rho * Cd(i) * wind(i) * u(i,pver)                          ! surface momentum flux x in N/m2
!       tauy(i)     = -rho * Cd(i) * wind(i) * v(i,pver)                          ! surface momentum flux y in N/m2
!       u(i,pver)   = u(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))                 ! update of the state variable u
!       v(i,pver)   = v(i,pver)/(1._r8+Cd(i)*wind(i)*dtime/za(i))                 ! update of the state variable v
!--------

!---------------------------
!       Sensible heat fluxes
!---------------------------
        dtdt(i,pver) = dtdt(i,pver) +((t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) & ! tendency due to sensible heat flux
                            /(1._r8+C*wind(i)*dtime/za(i))-t(i,pver))/dtime
        shflx(i)     = rho * cpair  * C*wind(i)*(Tsurf(i)-t(i,pver))               ! sensible heat flux in W/m2
        t(i,pver)    = (t(i,pver)+C*wind(i)*Tsurf(i)*dtime/za(i)) &                ! update of the state variable t
                            /(1._r8+C*wind(i)*dtime/za(i))  
!-------------------------
!       Latent heat fluxes
!-------------------------
        dqdt(i,pver) = dqdt(i,pver) +((q(i,pver)+C*wind(i)*qsats*dtime/za(i)) &    ! tendency due to latent heat flux
                            /(1._r8+C*wind(i)*dtime/za(i))-q(i,pver))/dtime
        lhflx(i)     = rho * latvap * C*wind(i)*(qsats-q(i,pver))                  ! latent heat flux in W/m2
        q(i,pver)    = (q(i,pver)+C*wind(i)*qsats*dtime/za(i)) &                   ! update of the state variable q
                            /(1._r8+C*wind(i)*dtime/za(i))
     end do
!===============================================================================
#endif


!===============================================================================
! Boundary layer mixing, see Reed and Jablonowski (JAMES, 2012)
!
! Note that Reed and Jablonowski (2012) also update u and v with the boundary-layer
! scheme which is commented out. The code is left here since it could be
! useful for other applications.
!===============================================================================
! Calculate Diagonal Variables for Implicit PBL Scheme
!
      do k=1,pver-1
         do i=1,pcols
            rho = (pint(i,k+1)/(rair*(t(i,k+1)*(1._r8+zvir*q(i,k+1))+t(i,k)*(1._r8+zvir*q(i,k)))/2.0_r8))! density of moist air
            CAm(i,k)   = rpdel(i,k)*dtime*gravit*gravit*Km(i,k+1)*rho*rho   &
                         /(pmid(i,k+1)-pmid(i,k))    
            CCm(i,k+1) = rpdel(i,k+1)*dtime*gravit*gravit*Km(i,k+1)*rho*rho &
                         /(pmid(i,k+1)-pmid(i,k))
            CA(i,k)    = rpdel(i,k)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho   &
                         /(pmid(i,k+1)-pmid(i,k))
            CC(i,k+1)  = rpdel(i,k+1)*dtime*gravit*gravit*Ke(i,k+1)*rho*rho &
                         /(pmid(i,k+1)-pmid(i,k))
         end do
      end do
      do i=1,pcols
         CAm(i,pver)   = 0._r8
         CCm(i,1)      = 0._r8
         CEm(i,pver+1) = 0._r8
         CA(i,pver)    = 0._r8
         CC(i,1)       = 0._r8
         CE(i,pver+1)  = 0._r8
         CFu(i,pver+1) = 0._r8
         CFv(i,pver+1) = 0._r8
         CFt(i,pver+1) = 0._r8
         CFq(i,pver+1) = 0._r8 
      end do
      do i=1,pcols
         do k=pver,1,-1
            CE(i,k)  = CC(i,k)/(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CEm(i,k) = CCm(i,k)/(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFu(i,k) = (u(i,k)+CAm(i,k)*CFu(i,k+1)) &
                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFv(i,k) = (v(i,k)+CAm(i,k)*CFv(i,k+1)) &
                       /(1._r8+CAm(i,k)+CCm(i,k)-CAm(i,k)*CEm(i,k+1))
            CFt(i,k) = ((p0/pmid(i,k))**(rair/cpair)*t(i,k)+CA(i,k)*CFt(i,k+1)) &
                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1)) 
            CFq(i,k) = (q(i,k)+CA(i,k)*CFq(i,k+1)) &
                       /(1._r8+CA(i,k)+CC(i,k)-CA(i,k)*CE(i,k+1))
        end do
      end do

!
! Calculate the effects of the boundary-layer mixing and update temperature and specific humidity
!
! First we need to calculate the updates at the top model level
!
      do i=1,pcols
!
!            Mixing of momentum is done via the Rayleigh friction
!
!            dudt(i,1)  = dudt(i,1)+(CFu(i,1)-u(i,1))/dtime
!            dvdt(i,1)  = dvdt(i,1)+(CFv(i,1)-v(i,1))/dtime
!            u(i,1)    = CFu(i,1)
!            v(i,1)    = CFv(i,1)
!
            dtdt(i,1)  = dtdt(i,1)+(CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)-t(i,1))/dtime
            t(i,1)    = CFt(i,1)*(pmid(i,1)/p0)**(rair/cpair)
            dqdt(i,1)  = dqdt(i,1)+(CFq(i,1)-q(i,1))/dtime
            q(i,1)  = CFq(i,1)
      end do
!
! Loop over the remaining level
!
      do i=1,pcols
         do k=2,pver
!
!            Mixing of momentum is done via the Rayleigh friction
!
!            dudt(i,k)  = dudt(i,k)+(CEm(i,k)*u(i,k-1)+CFu(i,k)-u(i,k))/dtime
!            dvdt(i,k)  = dvdt(i,k)+(CEm(i,k)*v(i,k-1)+CFv(i,k)-v(i,k))/dtime
!            u(i,k)    = CEm(i,k)*u(i,k-1)+CFu(i,k)
!            v(i,k)    = CEm(i,k)*v(i,k-1)+CFv(i,k)
!
            dtdt(i,k)  = dtdt(i,k)+((CE(i,k)*t(i,k-1) &
                              *(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
                              *(pmid(i,k)/p0)**(rair/cpair)-t(i,k))/dtime 
            t(i,k)    = (CE(i,k)*t(i,k-1)*(p0/pmid(i,k-1))**(rair/cpair)+CFt(i,k)) &
                              *(pmid(i,k)/p0)**(rair/cpair)
            dqdt(i,k)  = dqdt(i,k)+(CE(i,k)*q(i,k-1)+CFq(i,k)-q(i,k))/dtime
            q(i,k)  = CE(i,k)*q(i,k-1)+CFq(i,k)
         end do
      end do

!===============================================================================
! HS forcing
!------------------------------------------------------------------
! Held/Suarez IDEALIZED physics algorithm:
!
!   Held, I. M., and M. J. Suarez, 1994: A proposal for the
!   intercomparison of the dynamical cores of atmospheric general
!   circulation models.
!   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
!===============================================================================

!---------------------------------------------------------------------------------------
! Apply HS Rayleigh friction near the surface (below eta=0.7), replaces surface
! stresses and PBL diffusion which are not applied in the simple-physics package above
!---------------------------------------------------------------------------------------
      do k=1,pver
         if (etamid(k) > sigmab) then                      ! below approximately 700 hPa
            kv  = kf*(etamid(k) - sigmab)/onemsig          !
            tmp = -kv                                      ! Rayleigh friction coefficient
            do i=1,pcols
               dudt(i,k) = dudt(i,k) + tmp*u(i,k)          ! compute u tendency
               dvdt(i,k) = dvdt(i,k) + tmp*v(i,k)          ! compute v tendency
               u(i,k)    = u(i,k) +  tmp*u(i,k)*dtime      ! update u
               v(i,k)    = v(i,k) +  tmp*v(i,k)*dtime      ! update v
            end do
         endif
      end do

!-----------------------------------------------------------------------
! Compute idealized radiative heating rates (with modified HS equilibrium temperature)
!-----------------------------------------------------------------------
      do k=1,pver
         if (etamid(k) > sigmab) then                      ! below approximately 700 hPa, relaxation coefficient varies
            do i=1,pcols
               kt = ka + (ks - ka)*cos(lat(i))**4*(etamid(k) - sigmab)/onemsig
               tmp = kt
               trefc = T_max - delta_T*sin(lat(i))**2
               trefa = (trefc - delta_theta*cos(lat(i))**2*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**(rair/cpair)
               trefa = max(t00,trefa)
               dtdt(i,k) = dtdt(i,k) + (trefa - t(i,k))*tmp             ! compute T tendency
               t(i,k)    = t(i,k) + (trefa - t(i,k))*tmp*dtime          ! update T
            end do
         else
            tmp = ka
            do i=1,pcols
               trefc = T_max - delta_T*sin(lat(i))**2
               trefa = (trefc - delta_theta*cos(lat(i))**2*log((pmid(i,k)/ps0)))*(pmid(i,k)/ps0)**(rair/cpair)
               trefa = max(t00,trefa)
               dtdt(i,k) = dtdt(i,k) + (trefa - t(i,k))*tmp             ! compute T tendency
               t(i,k)    = t(i,k) + (trefa - t(i,k))*tmp*dtime          ! update T
            end do
         endif
      end do
!===============================================================================

!===============================================================================
! If desired send the diagnostic variables dudt, dvdt, dtdt, dqdt, shflx and lhflx
! to an output file - THIS PROCESS WILL BE MODEL SPECIFIC
!===============================================================================

!===============================================================================
! Adjustment of the moist pressure values after moisture has been removed or added to the
! vertical columns - THIS PROCESS WILL BE MODEL SPECIFIC 
!
! Insert or call the model's own pressure adjustment mechanism if needed
!
! note: Care needs to be taken to ensure that the model conserves the total
!       dry air mass. Add your own routine here.
!===============================================================================
!  call physics_dme_adjust(state, tend, qini, dtime)   ! This is for CESM/CAM
      return
   end subroutine grist_idealized_physics_mitc


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
!     theta  - potential temperature (K)
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

subroutine grist_idealized_physics_kessler(theta, qv, qc, qr, rho, pk, dt, z, nz, precl,&
                                           dtheta_dt, dqv_dt, dqc_dt, dqr_dt)
  IMPLICIT NONE

  !------------------------------------------------
  !   Input / output parameters
  !------------------------------------------------

  REAL(8), DIMENSION(nz), INTENT(INOUT) :: &
            theta   ,     & ! Potential temperature (K)
            qv      ,     & ! Water vapor mixing ratio (gm/gm)
            qc      ,     & ! Cloud water mixing ratio (gm/gm)
            qr              ! Rain  water mixing ratio (gm/gm)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            rho             ! Dry air density (not mean state as in KW) (kg/m^3)

  REAL(8), INTENT(OUT) :: &
            precl          ! Precipitation rate (m_water / s)

  REAL(8), DIMENSION(nz), INTENT(IN) :: &
            z       ,     & ! Heights of thermo. levels in the grid column (m)
            pk              ! Exner function (p/p0)**(R/cp)

  REAL(8), INTENT(IN) :: & 
            dt              ! Time step (s)

  INTEGER, INTENT(IN) :: nz ! Number of thermodynamic levels in the column

  REAL(8), DIMENSION(nz), INTENT(OUT) :: &
            dtheta_dt   ,     & ! Time tendency of Potential temperature (K)
            dqv_dt      ,     & ! Time tendency of Water vapor mixing ratio (gm/gm)
            dqc_dt      ,     & ! Time tendency of Cloud water mixing ratio (gm/gm)
            dqr_dt              ! Time tendency of Rain  water mixing ratio (gm/gm)

  !------------------------------------------------
  !   Local variables
  !------------------------------------------------
  REAL, DIMENSION(nz) :: r, rhalf, velqr, sed, pc

  REAL(8) :: f5, f2x, xk, ern, qrprod, prod, qvs, psl, rhoqr, dt_max, dt0

  INTEGER :: k, rainsplit, nt

  REAL(8), DIMENSION(nz)  :: theta_ini, qv_ini, qc_ini, qr_ini
!
!copy initial inputs
!
  theta_ini = theta
  qv_ini    = qv
  qc_ini    = qc
  qr_ini    = qr

  !------------------------------------------------
  !   Begin calculation
  !------------------------------------------------
  f2x = 17.27d0
  f5 = 237.3d0 * f2x * 2500000.d0 / 1003.d0
  xk = .2875d0      !  kappa (r/cp)
  psl    = 1000.d0  !  pressure at sea level (mb)
  rhoqr  = 1000.d0  !  density of liquid water (kg/m^3)

  do k=1,nz
    r(k)     = 0.001d0*rho(k)
    rhalf(k) = sqrt(rho(1)/rho(k))
    pc(k)    = 3.8d0/(pk(k)**(1./xk)*psl)

    ! Liquid water terminal velocity (m/s) following KW eq. 2.15
    velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364_r8*rhalf(k)

  end do

  ! Maximum time step size in accordance with CFL condition
  if (dt .le. 0.d0) then
    write(*,*) 'kessler.f90 called with nonpositive dt'
    stop
  end if

  dt_max = dt
  do k=1,nz-1
    if (velqr(k) .ne. 0.d0) then
      dt_max = min(dt_max, 0.8d0*(z(k+1)-z(k))/velqr(k))
    end if
  end do

  ! Number of subcycles
  rainsplit = ceiling(dt / dt_max)
  dt0 = dt / real(rainsplit,8)

  ! Subcycle through rain process
  precl = 0.d0

  do nt=1,rainsplit

    ! Precipitation rate (m/s)
    precl = precl + rho(1) * qr(1) * velqr(1) / rhoqr

    ! Sedimentation term using upstream differencing
    do k=1,nz-1
      sed(k) = dt0*(r(k+1)*qr(k+1)*velqr(k+1)-r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)))
    end do
    sed(nz)  = -dt0*qr(nz)*velqr(nz)/(.5*(z(nz)-z(nz-1)))

    ! Adjustment terms
    do k=1,nz

      ! Autoconversion and accretion rates following KW eq. 2.13a,b
      qrprod = qc(k) - (qc(k)-dt0*max(.001*(qc(k)-.001d0),0.))/(1.d0+dt0*2.2d0*qr(k)**.875)
      qc(k) = max(qc(k)-qrprod,0.)
      qr(k) = max(qr(k)+qrprod+sed(k),0.)

      ! Saturation vapor mixing ratio (gm/gm) following KW eq. 2.11
      qvs = pc(k)*exp(f2x*(pk(k)*theta(k)-273.d0)   &
             /(pk(k)*theta(k)- 36.d0))
      prod = (qv(k)-qvs)/(1.d0+qvs*f5/(pk(k)*theta(k)-36.d0)**2)

      ! Evaporation rate following KW eq. 2.14a,b
      ern = min(dt0*(((1.6d0+124.9d0*(r(k)*qr(k))**.2046_r8)  &
            *(r(k)*qr(k))**.525)/(2550000d0*pc(k)            &
            /(3.8d0 *qvs)+540000d0))*(dim(qvs,qv(k))         &
            /(r(k)*qvs)),max(-prod-qc(k),0.),qr(k))

      ! Saturation adjustment following KW eq. 3.10
      theta(k)= theta(k) + 2500000d0/(1003.d0*pk(k))*(max( prod,-qc(k))-ern)
      qv(k) = max(qv(k)-max(prod,-qc(k))+ern,0.)
      qc(k) = qc(k)+max(prod,-qc(k))
      qr(k) = qr(k)-ern
    end do

    ! Recalculate liquid water terminal velocity
    if (nt .ne. rainsplit) then
      do k=1,nz
        velqr(k)  = 36.34d0*(qr(k)*r(k))**0.1364_r8*rhalf(k)
      end do
    end if
  end do

  precl = precl / dble(rainsplit)
!
! convert changes to time tendency and do update outside
!
  dtheta_dt = (theta-theta_ini)/dt
  dqv_dt    = (qv-qv_ini)/dt
  dqc_dt    = (qc-qc_ini)/dt
  dqr_dt    = (qr-qr_ini)/dt

  return
  end subroutine grist_idealized_physics_kessler

!
! Implementation exactly follows that in Klemp et al. (2015), JAMES
!

 subroutine grist_idealized_physics_kessler_klemp15(pt,qv,qc,qr,rho,pk,dt,z,nz,rainnc,dpt_dt,dqv_dt,dqc_dt,dqr_dt)

    real(r8),    intent(inout) :: pt(nz), qv(nz), qc(nz), qr(nz)
    real(r8),    intent(in)    :: rho(nz), pk(nz), dt, z(nz)
    integer(i4), intent(in)    :: nz
    real(r8),    intent(out)   :: rainnc
    real(r8),    intent(out)   :: dpt_dt(nz), dqv_dt(nz), dqc_dt(nz), dqr_dt(nz)
! local
    real(r8)    :: r(nz), rhalf(nz), velqr(nz), sed(nz), pc(nz)
    real(r8)    :: rd, lv, ern, qrprod, prod, qvs, psl, rhoqr
    integer(i4) :: k
    real(r8)    :: pt_ini(nz), qv_ini(nz), qc_ini(nz), qr_ini(nz)

! copy initials
    pt_ini = pt
    qv_ini = qv
    qc_ini = qc
    qr_ini = qr

! constants
    rd    = rdry
!    cp    = cp
    lv    = 2.5e6_r8
    psl   = 1000._r8
    rhoqr = 1000._r8

    do k = 1, nz
       r(k)     = 0.001_r8*rho(k)
       rhalf(k) = sqrt(rho(1)/rho(k))
       velqr(k) = 36.34_r8*(qr(k)*r(k))**0.1364_r8*rhalf(k)
    end do

    rainnc = zero
    !rainnc = rainnc+1000._r8*rho(1)*qr(1)*velqr(1)*dt/rhoqr
    rainnc = rainnc+1000._r8*rho(1)*qr(1)*velqr(1)/rhoqr ! use rain rate (mm/s)

    do k = 1, nz-1
       sed(k) = dt*(r(k+1)*qr(k+1)*velqr(k+1) &
              -r(k)*qr(k)*velqr(k))/(r(k)*(z(k+1)-z(k)))
    end do
    sed(nz)   = -dt*qr(nz)*velqr(nz)/(0.5_r8*(z(nz)-z(nz-1)))

    do k = 1, nz
       qrprod = qc(k)-(qc(k)-dt*max(0.001_r8*(qc(k)-0.001_r8),zero))&
               /(one+dt*2.2_r8*qr(k)**0.875_r8)
       qc(k)  = max(qc(k)-qrprod,zero)
       qr(k)  = max(qr(k)+qrprod+sed(k),zero)

       pc(k)  = 3.8_r8/(pk(k)**(cp/rd)*psl)
       qvs    = pc(k)*exp(17.27*(pk(k)*pt(k)-273._r8) &
                       /(pk(k)*pt(k)-36._r8))

! water vapor adjustment to reach saturation
       prod   = (qv(k)-qvs)/(one+qvs*(4093._r8*lv/cp)/(pk(k)*pt(k)-36._r8)**2)
! evaporation
       ern    = min(dt*(((1.6_r8+124.9_r8*(r(k)*qr(k))**0.2046_r8)&
                      *(r(k)*qr(k))**0.525_r8)/(2.55e6*pc(k)  &
                      /(3.8_r8*qvs)+5.4e5_r8))*(dim(qvs,qv(k))   &
                       /(r(k)*qvs)),max(-prod-qc(k),zero),qr(k))
! saturation adjustment
       pt(k)  = pt(k)+lv/(cp*pk(k))*(max(prod,-qc(k))-ern)
       qv(k)  = max(qv(k)-max(prod,-qc(k))+ern,zero)
       qc(k)  = qc(k)+max(prod,-qc(k))
       qr(k)  = qr(k)-ern
    end do

! evaluate tendency
    dpt_dt = (pt-pt_ini)/dt
    dqv_dt = (qv-qv_ini)/dt
    dqc_dt = (qc-qc_ini)/dt
    dqr_dt = (qr-qr_ini)/dt

    return
 end subroutine grist_idealized_physics_kessler_klemp15
!=======================================================================

   subroutine phys_tend_hsdry_exp(mesh, nlev, dtime, &
                                         scalar_pressure_at_pc_full_level         , &
                                         scalar_pressure_at_pc_face_level         , &
                                         scalar_potential_temp_at_pc_full_level   , &
                                         scalar_normal_velocity_at_edge_full_level, &
                                         scalar_www_at_pc_face_level              , &
                                         tend_normal_velocity_at_edge_full_level  , &
                                         tend_www_at_pc_face_level                , &
                                         tend_potential_temp_at_pc_full_level)
! io   
      type(global_domain), intent(in) :: mesh
      integer(i4)        , intent(in) :: nlev
      real(r8),  intent(in)           :: dtime
      real(r8),  intent(in)           :: scalar_pressure_at_pc_full_level(:,:)       
      real(r8),  intent(in)           :: scalar_pressure_at_pc_face_level(:,:)         
      real(r8),  intent(in)           :: scalar_potential_temp_at_pc_full_level(:,:)
      real(r8),  intent(in)           :: scalar_normal_velocity_at_edge_full_level(:,:)
      real(r8),  intent(in)           :: scalar_www_at_pc_face_level(:,:)
      real(r8),  intent(inout)        :: tend_normal_velocity_at_edge_full_level(:,:)
      real(r8),  intent(inout)        :: tend_www_at_pc_face_level(:,:)    
      real(r8),  intent(inout)        :: tend_potential_temp_at_pc_full_level(:,:) 
! local
      real(r8),  parameter            :: eta_b        = 0.7_r8
      real(r8),  parameter            :: kf           = 1._r8/86400._r8
      real(r8),  parameter            :: ka           = 1._r8/(40*86400._r8)
      real(r8),  parameter            :: ks           = 1._r8/(4*86400._r8)
      real(r8),  parameter            :: delta_t      = 60_r8 
      real(r8),  parameter            :: delta_pt     = 10_r8
      real(r8),  parameter            :: zero         = 0_r8
      real(r8)                        :: kappa
      real(r8)                        :: ptemp_eq
      real(r8)                        :: kt
      real(r8)                        :: kv
      real(r8)                        :: tmpa
      real(r8)                        :: tmpb
      real(r8)                        :: eta1,eta2,eta, lat
      integer(i4)                     :: iv, ie, ilev, icell1, icell2, nlevp

          kappa     = rdry/cp
          nlevp     = nlev+1

! vertically vary
          do ie     = 1, mesh%ne
             do ilev   = 1, nlev
                icell1 = mesh%edt(ie)%v(1)
                icell2 = mesh%edt(ie)%v(2)
                eta1   = 0.5_r8*(scalar_pressure_at_pc_full_level(ilev,icell1)+scalar_pressure_at_pc_full_level(ilev,icell2))
                eta2   = 0.5_r8*(scalar_pressure_at_pc_face_level(nlevp,icell1)+scalar_pressure_at_pc_face_level(nlevp,icell2))
                eta    = eta1/eta2
                kv     = kf*max(zero,(eta-eta_b)/(1._r8-eta_b))
                tend_normal_velocity_at_edge_full_level(ilev,ie) = -kv*scalar_normal_velocity_at_edge_full_level(ilev,ie)
             end do
          end do

          do iv     = 1, mesh%nv
             do ilev   = 1, nlev+1
                eta    = scalar_pressure_at_pc_face_level(ilev,iv)/scalar_pressure_at_pc_face_level(nlevp,iv)
                kv     = kf*max(zero,(eta-eta_b)/(1._r8-eta_b))
                !tend_www_at_pc_face_level(ilev,iv) = -kv*scalar_www_at_pc_face_level(ilev,iv)
                tend_www_at_pc_face_level(ilev,iv) = 0._r8 !-kv*scalar_www_at_pc_face_level(ilev,iv)
             end do
          end do

          do iv     = 1, mesh%nv
             do ilev = 1, nlev
                lat       = mesh%vtx(iv)%lat
                eta       = scalar_pressure_at_pc_full_level(ilev,iv)/scalar_pressure_at_pc_face_level(nlevp,iv)
                kt        = ka+(ks-ka)*max(zero,(eta-eta_b)/(1._r8-eta_b))*(cos(lat)**4)
                tmpa      = 200._r8/((scalar_pressure_at_pc_full_level(ilev,iv)/p00)**kappa)
                tmpb      = 315._r8-delta_t*(sin(lat)**2)-delta_pt*log(scalar_pressure_at_pc_full_level(ilev,iv)/p00)*(cos(lat)**2)
                ptemp_eq  = max(tmpa,tmpb)
                tend_potential_temp_at_pc_full_level(ilev,iv) = -kt*(scalar_potential_temp_at_pc_full_level(ilev,iv)-ptemp_eq)
             end do
          end do

       return
    end subroutine phys_tend_hsdry_exp

!------------------------------------------------------
! implicit approach directly update the state
! we may also diagnose the tendency, send back
! and do an update in the driver part,
! doing so is same for primitive prognostic variables,
! while differ for compund prognostic variable pt
!------------------------------------------------------

    subroutine phys_tend_hsdry_imp(mesh, nlev, dtime, &
                                         scalar_pressure_at_pc_full_level         , &
                                         scalar_pressure_at_pc_face_level         , &
                                         scalar_potential_temp_at_pc_full_level   , &
                                         scalar_normal_velocity_at_edge_full_level, &
                                         scalar_www_at_pc_face_level              , &
                                         tend_normal_velocity_at_edge_full_level  , &
                                         tend_www_at_pc_face_level                , &
                                         tend_potential_temp_at_pc_full_level)
! io
      type(global_domain), intent(in) :: mesh
      integer(i4),intent(in)          :: nlev
      real(r8),   intent(in)          :: dtime
      real(r8),   intent(in)          :: scalar_pressure_at_pc_full_level(:,:)          ! full pressure (n+1)
      real(r8),   intent(in)          :: scalar_pressure_at_pc_face_level(:,:)          ! face pressure (n+1)
      real(r8),   intent(in)          :: scalar_potential_temp_at_pc_full_level(:,:)    ! state of potential temperature
      real(r8),   intent(in)          :: scalar_normal_velocity_at_edge_full_level(:,:) ! state of normal velocity (in:n, out: n+1)
      real(r8),   intent(in)          :: scalar_www_at_pc_face_level(:,:)               ! state of vertical velocity (in:n, out: n+1)
      real(r8),   intent(out)         :: tend_normal_velocity_at_edge_full_level(:,:)   ! implicit tendency
      real(r8),   intent(out)         :: tend_www_at_pc_face_level(:,:)                 ! implicit tendency 
      real(r8),   intent(out)         :: tend_potential_temp_at_pc_full_level(:,:)      ! implicit tendency 
! local
      real(r8),  parameter            :: eta_b        = 0.7_r8
      real(r8),  parameter            :: kf           = 1._r8/86400._r8
      real(r8),  parameter            :: ka           = 1._r8/(40*86400._r8)
      real(r8),  parameter            :: ks           = 1._r8/(4*86400._r8)
      real(r8),  parameter            :: delta_t      = 60_r8 
      real(r8),  parameter            :: delta_pt     = 10_r8
      real(r8),  parameter            :: zero         = 0_r8
      real(r8)                        :: kappa
      real(r8)                        :: ptemp_eq
      real(r8)                        :: kt
      real(r8)                        :: kv
      real(r8)                        :: tmpa
      real(r8)                        :: tmpb
      real(r8)                        :: eta1,eta2,eta, lat
      integer(i4)                     :: iv, ie, ilev, icell1, icell2, nlevp


          nlevp     = nlev+1
          kappa     = rdry/cp

          do ie     = 1, mesh%ne
             do ilev   = 1, nlev
                icell1 = mesh%edt(ie)%v(1)
                icell2 = mesh%edt(ie)%v(2)
                eta1   = 0.5_r8*(scalar_pressure_at_pc_full_level(ilev,icell1)+scalar_pressure_at_pc_full_level(ilev,icell2))
                eta2   = 0.5_r8*(scalar_pressure_at_pc_face_level(nlevp,icell1)+scalar_pressure_at_pc_face_level(nlevp,icell2))
                eta    = eta1/eta2
                kv     = kf*max(zero,(eta-eta_b)/(1._r8-eta_b))
                !scalar_normal_velocity_at_edge_full_level(ilev,ie) = scalar_normal_velocity_at_edge_full_level(ilev,ie)/(1._r8+kv*dtime)
                tend_normal_velocity_at_edge_full_level(ilev,ie) = -kv*scalar_normal_velocity_at_edge_full_level(ilev,ie)/(1._r8+dtime*kv)
             end do
          end do

          do iv     = 1, mesh%nv
             do ilev   = 1, nlev+1
                eta    = scalar_pressure_at_pc_face_level(ilev,iv)/scalar_pressure_at_pc_face_level(nlevp,iv)
                kv     = kf*max(zero,(eta-eta_b)/(1._r8-eta_b))
                !scalar_www_at_pc_face_level(ilev,iv) = scalar_www_at_pc_face_level(ilev,iv)/(1._r8+kv*dtime)
                tend_www_at_pc_face_level(ilev,iv) = -kv*scalar_www_at_pc_face_level(ilev,iv)/(1._r8+dtime*kv)
             end do
          end do

          do iv     = 1, mesh%nv
             do ilev = 1, nlev
                lat       = mesh%vtx(iv)%lat
                eta       = scalar_pressure_at_pc_full_level(ilev,iv)/scalar_pressure_at_pc_face_level(nlevp,iv)
                kt        = ka+(ks-ka)*max(zero,(eta-eta_b)/(1._r8-eta_b))*(cos(lat)**4)
                tmpa      = 200._r8/((scalar_pressure_at_pc_full_level(ilev,iv)/p00)**kappa)
                tmpb      = 315._r8-delta_t*(sin(lat)**2)-delta_pt*log(scalar_pressure_at_pc_full_level(ilev,iv)/p00)*(cos(lat)**2)
                ptemp_eq  = max(tmpa,tmpb)
                tend_potential_temp_at_pc_full_level(ilev,iv) = kt*(ptemp_eq-scalar_potential_temp_at_pc_full_level(ilev,iv))/(1._r8+kt*dtime)
             end do
          end do
      return
    end subroutine phys_tend_hsdry_imp

 end module grist_physics_idealized_package
