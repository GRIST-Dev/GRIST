!===================================================================================================
!
!  Created by LiXiaohan on 20/06/20, adopted from CAM5
!
!   MG microphysics version 1.5 - Update of MG microphysics and jumping-off
!       point for the development of MG2
!
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Peter Caldwell, Xiaohong Liu and Steve Ghan
! Version 2: Development begun: September 2011
! invoked in CAM by specifying -microphys=mg1.5
! 
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!
! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
! microphysics in cooperation with the MG liquid microphysics. This is
! controlled by the do_cldice variable.
!
! NOTE: If do_cldice is false, then MG microphysics should not update CLDICE
! or NUMICE; however, it is assumed that the other microphysics scheme will have
! updated CLDICE and NUMICE. The other microphysics should handle the following
! processes that would have been done by MG:
!   - Detrainment (liquid and ice)
!   - Homogeneous ice nucleation
!   - Heterogeneous ice nucleation
!   - Bergeron process
!   - Melting of ice
!   - Freezing of cloud drops
!   - Autoconversion (ice -> snow)
!   - Growth/Sublimation of ice
!   - Sedimentation of ice
!---------------------------------------------------------------------------------
! Based on micro_mg (restructuring of former cldwat2m_micro)
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
! Code comments added by HM, 093011
! General code structure:
!
! Code is divided into two main subroutines:
!   subroutine micro_mg_init --> initializes microphysics routine, should be called
!                                  once at start of simulation
!   subroutine micro_mg_tend --> main microphysics routine to be called each time step
!
! List of external functions:
!   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
!   qsat_ice --> for calculating saturation vapor pressure with respect to ice
!   gamma   --> standard mathematical gamma function
! .........................................................................
! List of inputs through use statement in fortran90:
! Variable Name                      Description                Units
! .........................................................................
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                  J kg-1 K-1
! tmelt           temperature of melting point for water          K
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! rh2o            gas constant for water vapor                  J kg-1 K-1
! latvap          latent heat of vaporization                   J kg-1
! latice          latent heat of fusion                         J kg-1
! qsat_water      external function for calculating liquid water
!                 saturation vapor pressure/humidity              -
! qsat_ice        external function for calculating ice
!                 saturation vapor pressure/humidity              pa
! rhmini          relative humidity threshold parameter for
!                 nucleating ice                                  -
! .........................................................................
! NOTE: List of all inputs/outputs passed through the call/subroutine statement
!       for micro_mg_tend is given below at the start of subroutine micro_mg_tend.
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure and specific humidity over water
! 3) svp over ice
!
!
!===================================================================================================

module micro_mg1_5

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8, i4
    use wv_sat_methods,                     only: qsat_water => wv_sat_qsat_water, &
                                                  qsat_ice   => wv_sat_qsat_ice
    use grist_mpi

implicit none
private
save

public :: micro_mg1_5_init ,  &
          micro_mg1_5_tend ,  & 
          micro_mg1_5_get_cols

! switch for specification rather than prediction of droplet and crystal number
! note: number will be adjusted as needed to keep mean size within bounds,
! even when specified droplet or ice number is used

! ***note: Even if constant cloud ice number is set, ice number is allowed
! to evolve based on process rates. This is needed in order to calculate
! the change in mass due to ice nucleation. All other ice microphysical
! processes are consistent with the specified constant ice number if
! this switch is turned on.

! nccons = .true. to specify constant cloud droplet number
! cicons = .true. to specify constant cloud ice number

logical, parameter, public :: nccons = .false.
logical, parameter, public :: nicons = .false.

!=========================================================
! Private module parameters
!=========================================================

real(r8), parameter :: pi = 3.14159265358979323846_r8 ! To 20 digits; more than enough
                                                          ! to reach the limit of double precision

real(r8), parameter :: omsm   = 0.99999_r8    ! number near unity for round-off issues

! parameters for specified ice and droplet number concentration
! note: these are local in-cloud values, not grid-mean
real(r8), parameter :: ncnst = 100.e6_r8    ! droplet num concentration when nccons=.true. (m-3)
real(r8), parameter :: ninst = 0.1e6_r8     ! ice num concentration when nicons=.true. (m-3)

! rhow used to be passed in in init, but then overwritten.
! For now, just setting it as a parameter here
real(r8), parameter :: rhosn = 250._r8  ! bulk density snow
real(r8), parameter :: rhoi = 500._r8   ! bulk density ice
real(r8), parameter :: rhow = 1000._r8  ! bulk density liquid

! fall speed parameters, V = aD^b (V is in m/s)
! droplets
real(r8), parameter :: ac = 3.e7_r8
real(r8), parameter :: bc = 2._r8
! snow
real(r8), parameter :: as = 11.72_r8    !Fall speed parameter for snow, 5.86-23.44, default:11.72
real(r8), parameter :: bs = 0.41_r8
! cloud ice
real(r8), parameter :: ai = 700._r8     !Fall speed parameter for cloud ice, 350-1400, default:700
real(r8), parameter :: bi = 1._r8
! rain
real(r8), parameter :: ar = 841.99667_r8
real(r8), parameter :: br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d
! exponent
real(r8), parameter :: dsph = 3._r8

! ventilation parameters
! for snow
real(r8), parameter :: f1s = 0.86_r8
real(r8), parameter :: f2s = 0.28_r8
! for rain
real(r8), parameter :: f1r = 0.78_r8
real(r8), parameter :: f2r = 0.308_r8

! collection efficiencies
! aggregation of cloud ice and snow, 0.001 - 1
real(r8), parameter :: eii = 0.1_r8

! autoconversion size threshold for cloud ice to snow (m), 100 - 500
real(r8), parameter :: dcs = 150.e-6_r8 !LiXH Test. default: 250.e-6_r8, DP:150.e-6_r8

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8  

! alternate threshold used for some in-cloud mmr
real(r8), parameter :: icsmall = 1.e-8_r8

! immersion freezing parameters, bigg 1953
real(r8), parameter :: bimm = 100._r8
real(r8), parameter :: aimm = 0.66_r8

! mass of new crystal due to aerosol freezing and growth (kg)
real(r8), parameter :: mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)

!Range of cloudsat reflectivities (dBz) for analytic simulator
real(r8), parameter :: csmin = -30._r8
real(r8), parameter :: csmax = 26._r8
real(r8), parameter :: mindbz = -99._r8
real(r8), parameter :: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: g           ! gravity
real(r8) :: r           ! dry air gas constant
real(r8) :: rv          ! water vapor gas constant
real(r8) :: cpp         ! specific heat of dry air
real(r8) :: tmelt       ! freezing point of water (K)

! latent heats of:
real(r8) :: xxlv        ! vaporization
real(r8) :: xlf         ! freezing
real(r8) :: xxls        ! sublimation

real(r8) :: rhmini      ! Minimum rh for ice cloud fraction > 0.

! flags
logical :: microp_uniform
logical :: do_cldice

real(r8) :: rhosu       ! typical 850mn air density

real(r8) :: icenuct     ! ice nucleation temperature: currently -5 degrees C

real(r8) :: snowmelt    ! what temp to melt all snow: currently 2 degrees C
real(r8) :: rainfrze    ! what temp to freeze all rain: currently -5 degrees C

! additional constants to help speed up code
real(r8) :: cons1
real(r8) :: cons4
real(r8) :: cons5
real(r8) :: cons7
real(r8) :: cons8
real(r8) :: cons11
real(r8) :: cons13
real(r8) :: cons14
real(r8) :: cons16
real(r8) :: cons17
real(r8) :: cons22
real(r8) :: cons23
real(r8) :: cons24
real(r8) :: cons25
real(r8) :: cons27
real(r8) :: cons28

! Generic interface for packing routines
interface pack_array
   module procedure pack_array_1Dr8
   module procedure pack_array_2Dr8
   module procedure pack_array_3Dr8
end interface

interface unpack_array
   module procedure unpack_array_1Dr8
   module procedure unpack_array_1Dr8_arrayfill
   module procedure unpack_array_2Dr8
   module procedure unpack_array_2Dr8_arrayfill
end interface

!===============================================================================
contains
!===============================================================================

subroutine micro_mg1_5_init( &
     gravit, rair, rh2o, cpair,    &
     tmelt_in, latvap, latice,           &
     rhmini_in, microp_uniform_in, do_cldice_in )
 
  !----------------------------------------------------------------------- 
  ! 
  ! Purpose: 
  ! initialize constants for MG microphysics
  ! 
  ! Author: Andrew Gettelman Dec 2005
  ! 
  !-----------------------------------------------------------------------
  
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: rhmini_in    ! Minimum rh for ice cloud fraction > 0.

  logical,  intent(in)  :: microp_uniform_in    ! .true. = configure uniform for sub-columns 
                                            ! .false. = use w/o sub-columns (standard)
  logical,  intent(in)  :: do_cldice_in     ! .true. = do all processes (standard)
                                            ! .false. = skip all processes affecting
                                            !           cloud ice

  !-----------------------------------------------------------------------

  ! declarations for MG code (transforms variable names)

  g= gravit                 ! gravity
  r= rair                   ! dry air gas constant: note units(phys_constants are in J/K/kmol)
  rv= rh2o                  ! water vapor gas constant
  cpp = cpair               ! specific heat of dry air
  tmelt = tmelt_in
  rhmini = rhmini_in

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! flags
  microp_uniform = microp_uniform_in
  do_cldice  = do_cldice_in

  ! typical air density at 850 mb

  rhosu = 85000._r8/(rair * tmelt)

  ! Maximum temperature at which snow is allowed to exist
  snowmelt = tmelt + 2._r8
  ! Minimum temperature at which rain is allowed to exist
  rainfrze = tmelt - 5._r8

  ! Ice nucleation temperature
  icenuct  = tmelt - 5._r8

  ! Define constants to help speed up code (this limits calls to gamma function)
  ! Unused names: cons6, cons15, cons21, cons26
  cons1=gamma(1._r8+dsph)
  cons4=gamma(1._r8+br)
  cons5=gamma(4._r8+br)
  cons7=gamma(1._r8+bs)     
  cons8=gamma(4._r8+bs)     
  cons11=gamma(3._r8+bs)
  cons13=gamma(5._r8/2._r8+br/2._r8)
  cons14=gamma(5._r8/2._r8+bs/2._r8)
  cons16=gamma(1._r8+bi)
  cons17=gamma(4._r8+bi)
  cons22=(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
  cons23=dcs**3
  cons24=dcs**2
  cons25=dcs**bs
  cons27=xxlv**2
  cons28=xxls**2

end subroutine micro_mg1_5_init

!===============================================================================
!microphysics routine for each timestep goes here...

subroutine micro_mg1_5_tend( mgncol,   mgcols,   nlev,     top_lev,  deltatin,           &
                             tn,                 qn,                                     &
                             qcn,                          qin,                          &
                             ncn,                          nin,                          &
                             relvarn,            accre_enhann,                           &     
                             pn,                 pdeln,              pint,               &
                             cldn,               liqcldf,            icecldf,            &
                             rate1ord_cw2pr_st,  naain,    npccnin,  rndstn,   naconin,  &
                             tlato,    qvlato,   qctendo,  qitendo,  nctendo,  nitendo,  &
                             effco,    effco_fn, effio,              precto,   precio,   &
                             nevapro, evapsnowo, praino,   prodsnowo,cmeouto,  deffio,   &
                             pgamrado, lamcrado, qsouto,   dsouto,   rsouto,             &
                             rflxo,    sflxo,    qrouto,   reff_raino,reff_snowo,        &
                             qcsevapo, qisevapo, qvreso,   cmeiout,  vtrmco,   vtrmio,   &
                             qcsedteno,qisedteno,prao,     prco,     mnuccco,  mnuccto,  &
                             msacwio,  psacwso,  bergso,   bergo,    melto,    homoo,    &
                             qcreso,             prcio,    praio,    qireso,             &
                             mnuccro,  pracso,   meltsdto, frzrdto,  mnuccdo,            &
                             nrouto,   nsouto,   reflo,    areflo,   areflzo,  freflo,   &
                             csrflo,   acsrflo,  fcsrflo,            rercldo,            &
                             ncaio,    ncalo,    qrouto2,  qsouto2,  nrouto2,  nsouto2,  &
                             drouto2,  dsouto2,  freqso,   freqro,   nficeo,             &
                             tnd_qsnown,         tnd_nsnown,         re_icen  )

  ! input arguments
  integer,  intent(in) :: mgncol                ! number of microphysics columns
  integer,  intent(in) :: mgcols(:)             ! list of microphysics columns
  integer,  intent(in) :: nlev                  ! number of layers
  integer,  intent(in) :: top_lev               ! top level to do microphysics
  real(r8), intent(in) :: deltatin              ! time step (s)
  real(r8), intent(in) :: tn(:,:)               ! input temperature (K)
  real(r8), intent(in) :: qn(:,:)               ! input h20 vapor mixing ratio (kg/kg)
  real(r8), intent(in) :: relvarn(:,:)          ! relative variance of cloud water (-)
  real(r8), intent(in) :: accre_enhann(:,:)     ! optional accretion enhancement factor (-)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(in) :: qcn(:,:)       ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:)       ! cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: ncn(:,:)       ! cloud water number conc (1/kg)
  real(r8), intent(in) :: nin(:,:)       ! cloud ice number conc (1/kg)
  real(r8), intent(in) :: pn(:,:)         ! air pressure (pa)
  real(r8), intent(in) :: pdeln(:,:)      ! pressure difference across level (pa)
  ! hm add 11-16-11, interface pressure
  real(r8), intent(in) :: pint(:,:)    ! level interface pressure (pa)
  real(r8), intent(in) :: cldn(:,:)      ! cloud fraction (no units)
  real(r8), intent(in) :: liqcldf(:,:)   ! liquid cloud fraction (no units)
  real(r8), intent(in) :: icecldf(:,:)   ! ice cloud fraction (no units)
  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8), intent(in) :: naain(:,:)     ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8), intent(in) :: npccnin(:,:)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8), intent(in) :: rndstn(:,:,:)  ! radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8), intent(in) :: naconin(:,:,:) ! number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)

  ! Used with CARMA cirrus microphysics
  ! (or similar external microphysics model)
  real(r8), intent(in) :: tnd_qsnown(:,:) ! snow mass tendency (kg/kg/s)
  real(r8), intent(in) :: tnd_nsnown(:,:) ! snow number tendency (#/kg/s)
  real(r8), intent(in) :: re_icen(:,:)    ! ice effective radius (m)

  ! output arguments

  real(r8), intent(out) :: rate1ord_cw2pr_st(:,:)    ! 1st order rate for
  ! direct cw to precip conversion
  real(r8), intent(out) :: tlato(:,:)         ! latent heating rate       (W/kg)
  real(r8), intent(out) :: qvlato(:,:)        ! microphysical tendency qv (1/s)
  real(r8), intent(out) :: qctendo(:,:)       ! microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitendo(:,:)       ! microphysical tendency qi (1/s)
  real(r8), intent(out) :: nctendo(:,:)       ! microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitendo(:,:)       ! microphysical tendency ni (1/(kg*s))
  real(r8), intent(out) :: effco(:,:)         ! droplet effective radius (micron)
  real(r8), intent(out) :: effco_fn(:,:)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8), intent(out) :: effio(:,:)         ! cloud ice effective radius (micron)
  real(r8), intent(out) :: precto(:)          ! surface precip rate (m/s)
  real(r8), intent(out) :: precio(:)          ! cloud ice/snow precip rate (m/s)
  real(r8), intent(out) :: nevapro(:,:)       ! evaporation rate of rain + snow (1/s)
  real(r8), intent(out) :: evapsnowo(:,:)     ! sublimation rate of snow (1/s)
  real(r8), intent(out) :: praino(:,:)        ! production of rain + snow (1/s)
  real(r8), intent(out) :: prodsnowo(:,:)     ! production of snow (1/s)
  real(r8), intent(out) :: cmeouto(:,:)       ! evap/sub of cloud (1/s)
  real(r8), intent(out) :: deffio(:,:)        ! ice effective diameter for optics (radiation) (micron)
  real(r8), intent(out) :: pgamrado(:,:)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8), intent(out) :: lamcrado(:,:)      ! slope of droplet distribution for optics (radiation) (1/m)
  real(r8), intent(out) :: qsouto(:,:)        ! snow mixing ratio (kg/kg)
  real(r8), intent(out) :: dsouto(:,:)        ! snow diameter (m)
  real(r8), intent(out) :: rsouto(:,:)        ! snow radius (m)     LiXH add for RRTMG 4DDA
  real(r8), intent(out) :: rflxo(:,:)         ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8), intent(out) :: sflxo(:,:)         ! grid-box average snow flux (kg m^-2 s^-1)
  real(r8), intent(out) :: qrouto(:,:)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8), intent(out) :: reff_raino(:,:)    ! rain effective radius (micron)
  real(r8), intent(out) :: reff_snowo(:,:)    ! snow effective radius (micron)
  real(r8), intent(out) :: qcsevapo(:,:)      ! cloud water evaporation due to sedimentation (1/s)
  real(r8), intent(out) :: qisevapo(:,:)      ! cloud ice sublimation due to sublimation (1/s)
  real(r8), intent(out) :: qvreso(:,:)        ! residual condensation term to ensure RH < 100% (1/s)
  real(r8), intent(out) :: cmeiout(:,:)       ! grid-mean cloud ice sub/dep (1/s)
  real(r8), intent(out) :: vtrmco(:,:)        ! mass-weighted cloud water fallspeed (m/s)
  real(r8), intent(out) :: vtrmio(:,:)        ! mass-weighted cloud ice fallspeed (m/s)
  real(r8), intent(out) :: qcsedteno(:,:)     ! qc sedimentation tendency (1/s)
  real(r8), intent(out) :: qisedteno(:,:)     ! qi sedimentation tendency (1/s)

  ! microphysical process rates for output (mixing ratio tendencies) (all have units of 1/s)
  real(r8), intent(out) :: prao(:,:)         ! accretion of cloud by rain 
  real(r8), intent(out) :: prco(:,:)         ! autoconversion of cloud to rain
  real(r8), intent(out) :: mnuccco(:,:)      ! mixing ratio tend due to immersion freezing
  real(r8), intent(out) :: mnuccto(:,:)      ! mixing ratio tend due to contact freezing
  real(r8), intent(out) :: msacwio(:,:)      ! mixing ratio tend due to H-M splintering
  real(r8), intent(out) :: psacwso(:,:)      ! collection of cloud water by snow
  real(r8), intent(out) :: bergso(:,:)       ! bergeron process on snow
  real(r8), intent(out) :: bergo(:,:)        ! bergeron process on cloud ice
  real(r8), intent(out) :: melto(:,:)        ! melting of cloud ice
  real(r8), intent(out) :: homoo(:,:)        ! homogeneous freezing cloud water
  real(r8), intent(out) :: qcreso(:,:)       ! residual cloud condensation due to removal of excess supersat
  real(r8), intent(out) :: prcio(:,:)        ! autoconversion of cloud ice to snow
  real(r8), intent(out) :: praio(:,:)        ! accretion of cloud ice by snow
  real(r8), intent(out) :: qireso(:,:)       ! residual ice deposition due to removal of excess supersat
  real(r8), intent(out) :: mnuccro(:,:)      ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8), intent(out) :: pracso(:,:)       ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8), intent(out) :: meltsdto(:,:)     ! latent heating rate due to melting of snow  (W/kg)
  real(r8), intent(out) :: frzrdto(:,:)      ! latent heating rate due to homogeneous freezing of rain (W/kg)
  real(r8), intent(out) :: mnuccdo(:,:)      ! mass tendency from ice nucleation
  real(r8), intent(out) :: nrouto(:,:)        ! rain number concentration (1/m3)
  real(r8), intent(out) :: nsouto(:,:)        ! snow number concentration (1/m3)
  real(r8), intent(out) :: reflo(:,:)         ! analytic radar reflectivity        
  real(r8), intent(out) :: areflo(:,:)        ! average reflectivity will zero points outside valid range
  real(r8), intent(out) :: areflzo(:,:)       ! average reflectivity in z.
  real(r8), intent(out) :: freflo(:,:)        ! fractional occurrence of radar reflectivity
  real(r8), intent(out) :: csrflo(:,:)        ! cloudsat reflectivity 
  real(r8), intent(out) :: acsrflo(:,:)       ! cloudsat average
  real(r8), intent(out) :: fcsrflo(:,:)       ! cloudsat fractional occurrence of radar reflectivity
  real(r8), intent(out) :: rercldo(:,:)       ! effective radius calculation for rain + cloud
  real(r8), intent(out) :: ncaio(:,:)        ! output number conc of ice nuclei available (1/m3)
  real(r8), intent(out) :: ncalo(:,:)        ! output number conc of CCN (1/m3)
  real(r8), intent(out) :: qrouto2(:,:)       ! copy of qrout as used to compute drout2
  real(r8), intent(out) :: qsouto2(:,:)       ! copy of qsout as used to compute dsout2
  real(r8), intent(out) :: nrouto2(:,:)       ! copy of nrout as used to compute drout2
  real(r8), intent(out) :: nsouto2(:,:)       ! copy of nsout as used to compute dsout2
  real(r8), intent(out) :: drouto2(:,:)       ! mean rain particle diameter (m)
  real(r8), intent(out) :: dsouto2(:,:)       ! mean snow particle diameter (m)
  real(r8), intent(out) :: freqso(:,:)        ! fractional occurrence of snow
  real(r8), intent(out) :: freqro(:,:)        ! fractional occurrence of rain
  real(r8), intent(out) :: nficeo(:,:)        ! fractional occurrence of ice

  ! local workspace
  ! all units mks unless otherwise stated

  ! parameters
  real(r8), parameter :: mincld = 0.0001_r8     ! minimum allowed cloud fraction
  real(r8), parameter :: cdnl   = 0.e6_r8       ! Lower bound on droplet number, 0 - 1e+7, default: 0
  ! assign number of sub-steps to iter
  ! use 2 sub-steps, following tests described in MG2008
  integer, parameter :: iter = 2

  ! local copies of input variables
  real(r8) :: q(nlev,mgncol)           ! water vapor mixing ratio (kg/kg)
  real(r8) :: t(nlev,mgncol)           ! temperature (K)
  real(r8) :: qc(nlev,mgncol)          ! cloud liquid mixing ratio (kg/kg)
  real(r8) :: qi(nlev,mgncol)          ! cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(nlev,mgncol)          ! cloud liquid number concentration (1/kg)
  real(r8) :: ni(nlev,mgncol)          ! cloud liquid number concentration (1/kg)
  real(r8) :: p(nlev,mgncol)           ! pressure (Pa)
  real(r8) :: pdel(nlev,mgncol)        ! pressure difference across level (Pa)
  real(r8) :: relvar(nlev,mgncol)      ! relative variance of cloud water (-)
  real(r8) :: accre_enhan(nlev,mgncol) ! optional accretion enhancement factor (-)

  real(r8) :: naai(nlev,mgncol)    ! ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8) :: npccn(nlev,mgncol)   ! ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  real(r8), allocatable :: rndst(:,:,:)
  real(r8), allocatable :: nacon(:,:,:)

  real(r8) :: tnd_qsnow(nlev,mgncol) ! snow mass tendency (kg/kg/s)
  real(r8) :: tnd_nsnow(nlev,mgncol) ! snow number tendency (#/kg/s)
  real(r8) :: re_ice(nlev,mgncol)    ! ice effective radius (m)

  ! Packed copies of output variables
  real(r8) :: tlat(nlev,mgncol)         ! latent heating rate       (W/kg)
  real(r8) :: qvlat(nlev,mgncol)        ! microphysical tendency qv (1/s)
  real(r8) :: qctend(nlev,mgncol)       ! microphysical tendency qc (1/s)
  real(r8) :: qitend(nlev,mgncol)       ! microphysical tendency qi (1/s)
  real(r8) :: nctend(nlev,mgncol)       ! microphysical tendency nc (1/(kg*s))
  real(r8) :: nitend(nlev,mgncol)       ! microphysical tendency ni (1/(kg*s))

  real(r8) :: effc(nlev,mgncol)         ! droplet effective radius (micron)
  real(r8) :: effc_fn(nlev,mgncol)      ! droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8) :: effi(nlev,mgncol)         ! cloud ice effective radius (micron)

  real(r8) :: prect(mgncol)          ! surface precip rate (m/s)
  real(r8) :: preci(mgncol)          ! cloud ice/snow precip rate (m/s)

  real(r8) :: nevapr(nlev,mgncol)       ! evaporation rate of rain + snow (1/s)
  real(r8) :: evapsnow(nlev,mgncol)     ! sublimation rate of snow (1/s)
  real(r8) :: prain(nlev,mgncol)        ! production of rain + snow (1/s)
  real(r8) :: prodsnow(nlev,mgncol)     ! production of snow (1/s)
  real(r8) :: cmeout(nlev,mgncol)       ! evap/sub of cloud (1/s)
  real(r8) :: deffi(nlev,mgncol)        ! ice effective diameter for optics (radiation) (micron)
  real(r8) :: pgamrad(nlev,mgncol)      ! ice gamma parameter for optics (radiation) (no units)
  real(r8) :: lamcrad(nlev,mgncol)      ! slope of droplet distribution for optics (radiation) (1/m)


  real(r8) :: qsout(nlev,mgncol)        ! snow mixing ratio (kg/kg)
  real(r8) :: qsout2(nlev,mgncol)       ! copy of qsout as used to compute dsout2
  real(r8) :: nsout(nlev,mgncol)        ! snow number concentration (1/m3)
  real(r8) :: nsout2(nlev,mgncol)       ! copy of nsout as used to compute dsout2
  real(r8) :: dsout(nlev,mgncol)        ! snow diameter (m)
  real(r8) :: dsout2(nlev,mgncol)       ! mean snow particle diameter (m)
  real(r8) :: rsout(nlev,mgncol)        ! snow radius (m)

  real(r8) :: qrout(nlev,mgncol)        ! grid-box average rain mixing ratio (kg/kg)
  real(r8) :: qrout2(nlev,mgncol)       ! copy of qrout as used to compute drout2
  real(r8) :: nrout(nlev,mgncol)        ! rain number concentration (1/m3)
  real(r8) :: nrout2(nlev,mgncol)       ! copy of nrout as used to compute drout2
  real(r8) :: drout2(nlev,mgncol)       ! mean rain particle diameter (m)

  real(r8) :: reff_rain(nlev,mgncol)    ! rain effective radius (micron)
  real(r8) :: reff_snow(nlev,mgncol)    ! snow effective radius (micron)

  real(r8) :: freqs(nlev,mgncol)        ! fractional occurrence of snow
  real(r8) :: freqr(nlev,mgncol)        ! fractional occurrence of rain

  real(r8) :: rflx(nlev+1,mgncol)       ! grid-box average rain flux (kg m^-2 s^-1)
  real(r8) :: sflx(nlev+1,mgncol)       ! grid-box average snow flux (kg m^-2 s^-1)

  real(r8) :: qcsevap(nlev,mgncol)      ! cloud water evaporation due to sedimentation (1/s)
  real(r8) :: qisevap(nlev,mgncol)      ! cloud ice sublimation due to sublimation (1/s)
  real(r8) :: qvres(nlev,mgncol)        ! residual condensation term to ensure RH < 100% (1/s)
  real(r8) :: cmeitot(nlev,mgncol)      ! grid-mean cloud ice sub/dep (1/s)
  real(r8) :: vtrmc(nlev,mgncol)        ! mass-weighted cloud water fallspeed (m/s)
  real(r8) :: vtrmi(nlev,mgncol)        ! mass-weighted cloud ice fallspeed (m/s)
  real(r8) :: qcsedten(nlev,mgncol)     ! qc sedimentation tendency (1/s)
  real(r8) :: qisedten(nlev,mgncol)     ! qi sedimentation tendency (1/s)

  real(r8) :: pratot(nlev,mgncol)         ! accretion of cloud by rain 
  real(r8) :: prctot(nlev,mgncol)         ! autoconversion of cloud to rain
  real(r8) :: mnuccctot(nlev,mgncol)      ! mixing ratio tend due to immersion freezing
  real(r8) :: mnuccttot(nlev,mgncol)      ! mixing ratio tend due to contact freezing
  real(r8) :: msacwitot(nlev,mgncol)      ! mixing ratio tend due to H-M splintering
  real(r8) :: psacwstot(nlev,mgncol)      ! collection of cloud water by snow
  real(r8) :: bergstot(nlev,mgncol)       ! bergeron process on snow
  real(r8) :: bergtot(nlev,mgncol)        ! bergeron process on cloud ice
  real(r8) :: melttot(nlev,mgncol)        ! melting of cloud ice
  real(r8) :: homotot(nlev,mgncol)        ! homogeneous freezing cloud water
  real(r8) :: qcrestot(nlev,mgncol)       ! residual cloud condensation due to removal of excess supersat
  real(r8) :: prcitot(nlev,mgncol)        ! autoconversion of cloud ice to snow
  real(r8) :: praitot(nlev,mgncol)        ! accretion of cloud ice by snow
  real(r8) :: qirestot(nlev,mgncol)       ! residual ice deposition due to removal of excess supersat
  real(r8) :: mnuccrtot(nlev,mgncol)      ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8) :: pracstot(nlev,mgncol)      ! mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8) :: mnuccdtot(nlev,mgncol)      ! mass tendency from ice nucleation
  real(r8) :: meltsdttot(nlev,mgncol)      ! latent heating rate due to melting of snow  (W/kg)
  real(r8) :: frzrdttot(nlev,mgncol)      ! latent heating rate due to homogeneous freezing of rain (W/kg)

  real(r8) :: refl(nlev,mgncol)         ! analytic radar reflectivity        
  real(r8) :: arefl(nlev,mgncol)        ! average reflectivity will zero points outside valid range
  real(r8) :: areflz(nlev,mgncol)       ! average reflectivity in z.
  real(r8) :: frefl(nlev,mgncol)        ! fractional occurrence of radar reflectivity
  real(r8) :: csrfl(nlev,mgncol)        ! cloudsat reflectivity 
  real(r8) :: acsrfl(nlev,mgncol)       ! cloudsat average
  real(r8) :: fcsrfl(nlev,mgncol)       ! cloudsat fractional occurrence of radar reflectivity

  real(r8) :: rercld(nlev,mgncol)       ! effective radius calculation for rain + cloud

  real(r8) :: nfice(nlev,mgncol)        ! fractional occurrence of ice

  real(r8) :: ncai(nlev,mgncol)    ! output number conc of ice nuclei available (1/m3)
  real(r8) :: ncal(nlev,mgncol)    ! output number conc of CCN (1/m3)

  ! general purpose variables
  real(r8) :: deltat            ! sub-time step (s)
  real(r8) :: mtime             ! the assumed ice nucleation timescale

  real(r8) :: dz(nlev,mgncol)         ! height difference across model vertical level
  real(r8) :: int_to_mid(nlev,mgncol) ! Coefficients for linear interpolation from
                                      ! interface to mid-level

  ! temporary variables for sub-stepping
  real(r8) :: tlat1(nlev,mgncol)
  real(r8) :: qvlat1(nlev,mgncol)
  real(r8) :: qctend1(nlev,mgncol)
  real(r8) :: qitend1(nlev,mgncol)
  real(r8) :: nctend1(nlev,mgncol)
  real(r8) :: nitend1(nlev,mgncol)
  real(r8) :: prect1(mgncol)
  real(r8) :: preci1(mgncol)

  ! physical properties of the air at a given point
  real(r8) :: rho(nlev,mgncol)    ! density (kg m-3)
  real(r8) :: dv(nlev,mgncol)     ! diffusivity of water vapor
  real(r8) :: mu(nlev,mgncol)     ! viscosity
  real(r8) :: sc(nlev,mgncol)     ! schmidt number
  real(r8) :: rhof(nlev,mgncol)   ! density correction factor for fallspeed

  ! cloud fractions
  real(r8) :: cldmax(nlev,mgncol) ! precip fraction assuming maximum overlap
  real(r8) :: cldm(nlev,mgncol)   ! cloud fraction
  real(r8) :: icldm(nlev,mgncol)  ! ice cloud fraction
  real(r8) :: lcldm(nlev,mgncol)  ! liq cloud fraction

  ! mass mixing ratios
  real(r8) :: qcic(nlev,mgncol)   ! in-cloud cloud liquid
  real(r8) :: qiic(nlev,mgncol)   ! in-cloud cloud ice
  real(r8) :: qsic(nlev,mgncol)   ! in-precip snow
  real(r8) :: qric(nlev,mgncol)   ! in-precip rain

  ! number concentrations
  real(r8) :: ncic(nlev,mgncol)   ! in-cloud droplet
  real(r8) :: niic(nlev,mgncol)   ! in-cloud cloud ice
  real(r8) :: nsic(nlev,mgncol)   ! in-precip snow
  real(r8) :: nric(nlev,mgncol)   ! in-precip rain
  ! maximum allowed ni value
  real(r8) :: nimax(nlev,mgncol)

  ! Size distribution parameters for:
  ! cloud ice
  real(r8) :: lami(nlev,mgncol)   ! slope
  real(r8) :: n0i(nlev,mgncol)    ! intercept
  ! cloud liquid
  real(r8) :: lamc(nlev,mgncol)   ! slope
  real(r8) :: pgam(nlev,mgncol)   ! spectral width parameter
  real(r8) :: cdist1(nlev,mgncol) ! droplet freezing calculation
  ! snow
  real(r8) :: lams(nlev,mgncol)   ! slope
  real(r8) :: n0s(nlev,mgncol)    ! intercept
  ! rain
  real(r8) :: lamr(nlev,mgncol)   ! slope
  real(r8) :: n0r(nlev,mgncol)    ! intercept

  ! combined size of precip & cloud drops
  integer :: arcld(nlev,mgncol)  ! averaging control flag

  ! Rates/tendencies due to:
  ! deposition/sublimation of cloud ice
  real(r8) :: cmei(nlev,mgncol)
  ! ice nucleation
  real(r8) :: nnuccd(nlev,mgncol) ! number rate from deposition/cond.-freezing
  real(r8) :: mnuccd(nlev,mgncol) ! mass mixing ratio
  ! freezing of cloud water
  real(r8) :: mnuccc(nlev,mgncol) ! mass mixing ratio
  real(r8) :: nnuccc(nlev,mgncol) ! number concentration
  ! contact freezing of cloud water
  real(r8) :: mnucct(nlev,mgncol) ! mass mixing ratio
  real(r8) :: nnucct(nlev,mgncol) ! number concentration
  ! HM ice multiplication
  real(r8) :: msacwi(nlev,mgncol) ! mass mixing ratio
  real(r8) :: nsacwi(nlev,mgncol) ! number conc
  ! autoconversion of cloud droplets
  real(r8) :: prc(nlev,mgncol)    ! mass mixing ratio
  real(r8) :: nprc(nlev,mgncol)   ! number concentration (rain)
  real(r8) :: nprc1(nlev,mgncol)  ! number concentration (cloud droplets)
  ! self-aggregation of snow
  real(r8) :: nsagg(nlev,mgncol)  ! number concentration
  ! self-collection of rain
  real(r8) :: nragg(nlev,mgncol)  ! number concentration
  ! collection of droplets by snow
  real(r8) :: psacws(nlev,mgncol)     ! mass mixing ratio
  real(r8) :: npsacws(nlev,mgncol)    ! number concentration
  ! collection of rain by snow
  real(r8) :: pracs(nlev,mgncol)  ! mass mixing ratio
  real(r8) :: npracs(nlev,mgncol) ! number concentration
  ! freezing of rain
  real(r8) :: mnuccr(nlev,mgncol) ! mass mixing ratio
  real(r8) :: nnuccr(nlev,mgncol) ! number concentration
  ! accretion of droplets by rain
  real(r8) :: pra(nlev,mgncol)    ! mass mixing ratio
  real(r8) :: npra(nlev,mgncol)   ! number concentration
  ! autoconversion of cloud ice to snow
  real(r8) :: prci(nlev,mgncol)   ! mass mixing ratio
  real(r8) :: nprci(nlev,mgncol)  ! number concentration
  ! accretion of cloud ice by snow
  real(r8) :: prai(nlev,mgncol)   ! mass mixing ratio
  real(r8) :: nprai(nlev,mgncol)  ! number concentration
  ! evaporation of rain
  real(r8) :: pre(nlev,mgncol)    ! mass mixing ratio
  ! sublimation of snow
  real(r8) :: prds(nlev,mgncol)   ! mass mixing ratio
  ! number evaporation
  real(r8) :: nsubi(nlev,mgncol)  ! cloud ice
  real(r8) :: nsubc(nlev,mgncol)  ! droplet
  real(r8) :: nsubs(nlev,mgncol)  ! snow
  real(r8) :: nsubr(nlev,mgncol)  ! rain
  ! bergeron process
  real(r8) :: berg(nlev,mgncol)   ! mass mixing ratio (cloud ice)
  real(r8) :: bergs(nlev,mgncol)  ! mass mixing ratio (snow)

  ! fallspeeds
  ! number-weighted
  real(r8) :: uns(nlev,mgncol)    ! snow
  real(r8) :: unr(nlev,mgncol)    ! rain
  ! mass-weighted
  real(r8) :: ums(nlev,mgncol)    ! snow
  real(r8) :: umr(nlev,mgncol)    ! rain
  ! air density corrected fallspeed parameters
  real(r8) :: arn(nlev,mgncol)    ! rain
  real(r8) :: asn(nlev,mgncol)    ! snow
  real(r8) :: acn(nlev,mgncol)    ! cloud droplet
  real(r8) :: ain(nlev,mgncol)    ! cloud ice

  ! saturation vapor pressures
  real(r8) :: esl(nlev,mgncol)    ! liquid
  real(r8) :: esi(nlev,mgncol)    ! ice
  real(r8) :: esn               ! checking for RH after rain evap

  ! saturation vapor mixing ratios
  real(r8) :: qvl(nlev,mgncol)        ! liquid
  real(r8) :: qvi(nlev,mgncol)        ! ice
  real(r8) :: qsn                   ! checking for RH after rain evap

  ! relative humidity
  real(r8) :: relhum(nlev,mgncol)

  ! parameters for cloud water and cloud ice sedimentation calculations
  real(r8) :: fc(nlev)
  real(r8) :: fnc(nlev)
  real(r8) :: fi(nlev)
  real(r8) :: fni(nlev)

  real(r8) :: faloutc(nlev)
  real(r8) :: faloutnc(nlev)
  real(r8) :: falouti(nlev)
  real(r8) :: faloutni(nlev)

  real(r8) :: faltndc
  real(r8) :: faltndnc
  real(r8) :: faltndi
  real(r8) :: faltndni
  real(r8) :: faltndqie
  real(r8) :: faltndqce

  ! sum of source/sink terms for diagnostic precip
  real(r8) :: qstend(nlev,mgncol)    ! snow mixing ratio
  real(r8) :: nstend(nlev,mgncol)     ! snow number concentration
  real(r8) :: qrtend(nlev,mgncol)     ! rain mixing ratio
  real(r8) :: nrtend(nlev,mgncol)     ! rain number concentration
  ! vertically integrated source/sink terms
  real(r8) :: qrtot(mgncol)           ! rain mixing ratio
  real(r8) :: nrtot(mgncol)           ! rain number concentration
  real(r8) :: qstot(mgncol)           ! snow mixing ratio
  real(r8) :: nstot(mgncol)           ! snow number concentration

  ! for calculation of rate1ord_cw2pr_st
  real(r8) :: qcsinksum_rate1ord(nlev,mgncol) ! sum over iterations of cw to precip sink
  real(r8) :: qcsum_rate1ord(nlev,mgncol)     ! sum over iterations of cloud water

  real(r8) :: rainrt(nlev,mgncol)    ! rain rate for reflectivity calculation

  ! dummy variables
  real(r8) :: dum
  real(r8) :: dum1
  ! dummies for checking RH
  real(r8) :: qtmp
  real(r8) :: ttmp
  ! dummies for conservation check
  real(r8) :: qce                   ! qc
  real(r8) :: qie                   ! qi
  real(r8) :: nce                   ! nc
  real(r8) :: nie                   ! ni
  real(r8) :: ratio
  ! dummies for in-cloud variables
  real(r8) :: dumc(nlev,mgncol)   ! qc
  real(r8) :: dumnc(nlev,mgncol)  ! nc
  real(r8) :: dumi(nlev,mgncol)   ! qi
  real(r8) :: dumni(nlev,mgncol)  ! ni
  real(r8) :: dumr(nlev,mgncol)   ! rain mixing ratio
  real(r8) :: dumnr(nlev,mgncol)  ! rain number concentration
  ! Array dummy variable
  real(r8) :: dum_2D(nlev,mgncol)

  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  ! "ii" and "kk" are used for indices into input/output buffers
  ! "it" is substepping variable
  ! "n" is used for other iterations (currently just sedimentation)
  integer i, ii, k, kk, it, n

  ! number of iterations for loops over "n"
  integer nstep

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! Process inputs

  ! assign variable deltat for sub-stepping...
  deltat = deltatin

  call pack_array(          qn, mgcols, top_lev,           q)
  call pack_array(          tn, mgcols, top_lev,           t)
  call pack_array(         qcn, mgcols, top_lev,          qc)
  call pack_array(         qin, mgcols, top_lev,          qi)
  call pack_array(         ncn, mgcols, top_lev,          nc)
  call pack_array(         nin, mgcols, top_lev,          ni)
  call pack_array(          pn, mgcols, top_lev,           p)
  call pack_array(       pdeln, mgcols, top_lev,        pdel)
  call pack_array(     relvarn, mgcols, top_lev,      relvar)
  call pack_array(accre_enhann, mgcols, top_lev, accre_enhan)

  call pack_array(  naain, mgcols, top_lev,  naai)
  call pack_array(npccnin, mgcols, top_lev, npccn)

  ! These are allocated instead of used as automatic arrays
  ! purely to work around a PGI bug.
  allocate(rndst(size(rndstn,1),nlev,mgncol))
  allocate(nacon(size(rndstn,1),nlev,mgncol))
  call pack_array( rndstn, mgcols, top_lev, rndst)
  call pack_array(naconin, mgcols, top_lev, nacon)

  if (.not. do_cldice) then
     call pack_array(tnd_qsnown, mgcols, top_lev, tnd_qsnow)
     call pack_array(tnd_nsnown, mgcols, top_lev, tnd_nsnow)
     call pack_array(   re_icen, mgcols, top_lev,    re_ice)
  end if

  ! Some inputs are only used once, to create a reference array that is
  ! used repeatedly later. Rather than bothering to pack these, just
  ! set the local reference directly from the inputs.

  ! pint: used to set int_to_mid
  ! interface to mid-level linear interpolation
  do k = 1,nlev
     do i = 1, mgncol
        ! Set ii and kk to values that correspond to i and k.
        ii = mgcols(i)
        kk = k + top_lev - 1

        int_to_mid(k,i) = (p(k,i) - pint(kk,ii))/ &
             (pint(kk+1,ii) - pint(kk,ii))
     end do
  end do

  ! cldn: used to set cldm, unused for subcolumns
  ! liqcldf: used to set lcldm, unused for subcolumns
  ! icecldf: used to set icldm, unused for subcolumns

  if (microp_uniform) then
     ! subcolumns, set cloud fraction variables to one
     ! if cloud water or ice is present, if not present
     ! set to mincld (mincld used instead of zero, to prevent
     ! possible division by zero errors).

     where (qc >= qsmall)
        lcldm = 1._r8
     elsewhere
        lcldm = mincld
     end where

     where (qi >= qsmall)
        icldm = 1._r8
     elsewhere
        icldm = mincld
     end where

     cldm = max(icldm, lcldm)

  else
     ! get cloud fraction, check for minimum

     do k = 1,nlev
        do i = 1, mgncol
           ii = mgcols(i)
           kk = k + top_lev - 1

           cldm(k,i) = max(cldn(kk,ii),mincld)
           lcldm(k,i) = max(liqcldf(kk,ii),mincld)
           icldm(k,i) = max(icecldf(kk,ii),mincld)
        end do
     end do
  end if

  ! Initialize local variables

  ! local physical properties
  rho = p/(r*t)
  dv = 8.794E-5_r8 * t**1.81_r8 / p
  mu = 1.496E-6_r8 * t**1.5_r8 / (t + 120._r8)
  sc = mu/(rho*dv)

  ! get dz from dp and hydrostatic approx
  ! keep dz positive (define as layer k-1 - layer k)
  dz = pdel/(rho*g)

  ! air density adjustment for fallspeed parameters
  ! includes air density correction factor to the
  ! power of 0.54 following Heymsfield and Bansemer 2007

  rhof=(rhosu/rho)**0.54_r8

  arn=ar*rhof
  asn=as*rhof
  acn=g*rhow/(18._r8*mu)
  ain=ai*(rhosu/rho)**0.35_r8

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Get humidity and saturation vapor pressures

  do k=1,nlev
     do i=1,mgncol

        call qsat_water(t(k,i), p(k,i), esl(k,i), qvl(k,i))

        ! hm fix, make sure when above freezing that esi=esl, not active yet
        if (t(k,i) >= tmelt) then
           esi(k,i)=esl(k,i)
           qvi(k,i)=qvl(k,i)
        else
           call qsat_ice(t(k,i), p(k,i), esi(k,i), qvi(k,i))
        end if

     end do
  end do

  where (qvl <= 0.0_r8)
     relhum = q
  elsewhere
     relhum = q / min(1.0_r8,qvl)
  end where

  !===============================================
  ! Processes done before substepping
  !===============================================

  ! Initial deposition/sublimation of ice
  !===========================================

  if (do_cldice) then

     call ice_deposition_sublimation_init(deltat, t, q, qc, qi, ni,     &
                                          lcldm, icldm, naai, rho, dv,  &
                                          esl, esi, qvl, qvi, relhum,   &
                                          berg, cmei)

  else
     berg = 0._r8
     cmei = 0._r8
  end if  ! end do_cldice

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! droplet activation
  ! hm, modify 5/12/11 
  ! get provisional droplet number after activation. This is used for
  ! all microphysical process calculations, for consistency with update of
  ! droplet mass before microphysics

  ! calculate potential for droplet activation if cloud water is present
  ! tendency from activation (npccn) is read in from companion routine
  
  ! output activated liquid and ice (convert from #/kg -> #/m3)
  !--------------------------------------------------
  where (qc >= qsmall)
     nc = nc + npccn*deltat
     ncal = max(nc/lcldm,cdnl/rho)*rho ! sghan minimum in #/cm3
  elsewhere
     ncal = 0._r8
  end where

  where (t < icenuct)
     ncai = naai*rho
  elsewhere
     ncai = 0._r8
  end where

  !INITIALIZE STUFF FOR SUBSTEPPING
  !===============================================

  ! get sub-step time step
  deltat=deltat/real(iter)

  ! hm, set mtime here to avoid answer-changing
  mtime=deltat

  ! initialize tendencies to zero
  tlat1 = 0._r8
  qvlat1 = 0._r8
  qctend1 = 0._r8
  qitend1 = 0._r8
  nctend1 = 0._r8
  nitend1 = 0._r8

  ! initialize microphysics output
  qcsevap=0._r8
  qisevap=0._r8
  qvres  =0._r8
  cmeitot =0._r8
  vtrmc =0._r8
  vtrmi =0._r8
  qcsedten =0._r8
  qisedten =0._r8

  pratot=0._r8
  prctot=0._r8
  mnuccctot=0._r8
  mnuccttot=0._r8
  msacwitot=0._r8
  psacwstot=0._r8
  bergstot=0._r8
  bergtot=0._r8
  melttot=0._r8
  homotot=0._r8
  qcrestot=0._r8
  prcitot=0._r8
  praitot=0._r8
  qirestot=0._r8
  mnuccrtot=0._r8
  pracstot=0._r8
  meltsdttot=0._r8
  frzrdttot=0._r8
  mnuccdtot=0._r8

  rflx=0._r8
  sflx=0._r8

  ! initialize precip output
  
  qrout=0._r8
  qsout=0._r8
  nrout=0._r8
  nsout=0._r8

  ! for refl calc
  rainrt = 0._r8

  ! initialize rain size
  rercld=0._r8
  arcld = 0

  qcsinksum_rate1ord = 0._r8 
  qcsum_rate1ord     = 0._r8 

  ! initialize variables for trop_mozart
  nevapr = 0._r8
  evapsnow = 0._r8
  prain = 0._r8
  prodsnow = 0._r8
  cmeout = 0._r8

  prect1 = 0._r8
  preci1 = 0._r8

  cldmax = mincld

  lamc=0._r8


  !*********DO SUBSTEPPING!***************
  !============================================
  substepping: do it=1,iter

     ! initialize sub-step microphysical tendencies
     tlat=0._r8
     qvlat=0._r8
     qctend=0._r8
     qitend=0._r8
     qstend = 0._r8
     qrtend = 0._r8
     nctend=0._r8
     nitend=0._r8
     nrtend = 0._r8
     nstend = 0._r8

     ! initialize diagnostic precipitation to zero
     qcic  = 0._r8
     qiic  = 0._r8
     qsic  = 0._r8
     qric  = 0._r8

     ncic  = 0._r8
     niic  = 0._r8
     nsic  = 0._r8
     nric  = 0._r8

     ! initialize precip at surface
     prect = 0._r8
     preci = 0._r8

     ! initialize vertically-integrated rain and snow tendencies
     qrtot = 0._r8
     nrtot = 0._r8
     qstot = 0._r8
     nstot = 0._r8

     ! recalculate saturation vapor pressure for liquid and ice
     do k = 1, nlev
        do i = 1, mgncol

           call qsat_water(t(k,i), p(k,i), esl(k,i), qvl(k,i))

           ! hm fix, make sure when above freezing that esi=esl, not active yet
           if (t(k,i) >= tmelt) then
              esi(k,i)=esl(k,i)
              qvi(k,i)=qvl(k,i)
           else
              call qsat_ice(t(k,i), p(k,i), esi(k,i), qvi(k,i))
           end if

        end do
     end do

     where (qvl <= 0.0_r8)
        relhum = q
     elsewhere
        relhum = q / min(1.0_r8,qvl)
     end where

     ! decrease in number concentration due to sublimation/evap
     !-------------------------------------------------------
     ! divide by cloud fraction to get in-cloud decrease
     ! don't reduce Nc due to bergeron process

     where (cmei < 0._r8 .and. qi > qsmall .and. icldm > mincld)
        nsubi = cmei / qi * ni / icldm
     elsewhere
        nsubi = 0._r8
     end where

     nsubc = 0._r8

     ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
     !-------------------------------------------------------

     if (do_cldice) then
        where (naai > 0._r8 .and. t < icenuct .and. &
          relhum*esl/esi > rhmini+0.05_r8)

           !if NAAI > 0. then set numice = naai (as before)
           !note: this is gridbox averaged
           ! hm, modify to use mtime
           nnuccd = (naai-ni/icldm)/mtime*icldm
           nnuccd = max(nnuccd,0._r8)
           nimax = naai*icldm

           !Calc mass of new particles using new crystal mass...
           !also this will be multiplied by mtime as nnuccd is...

           mnuccd = nnuccd * mi0

           !  add mnuccd to cmei....
           cmei = cmei + mnuccd
        
           !  limit cmei
           !-------------------------------------------------------
           cmei = min(cmei,(q-qvi)/calc_ab(t, qvi, xxls)/deltat)

           ! limit for roundoff error
           cmei = cmei * omsm

        elsewhere
           nnuccd = 0._r8
           nimax = 0._r8
           mnuccd = 0._r8
        end where

     end if

     pre_vert_loop: do k=1,nlev

        pre_col_loop: do i=1,mgncol

           ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
           !-------------------------------------------------------
           ! for microphysical process calculations
           ! units are kg/kg for mixing ratio, 1/kg for number conc

           if (qc(k,i) - berg(k,i)*deltat.ge.qsmall) then
              ! limit in-cloud values to 0.005 kg/kg
              qcic(k,i)=min(qc(k,i)/lcldm(k,i),5.e-3_r8)
              ncic(k,i)=max(nc(k,i)/lcldm(k,i),0._r8)

              ! hm add 6/2/11 specify droplet concentration
              if (nccons) then
                 ncic(k,i)=ncnst/rho(k,i)
              end if
           else
              qcic(k,i)=0._r8
              ncic(k,i)=0._r8

              berg(k,i)=qc(k,i)/deltat*omsm
           end if

           if (qi(k,i)+(cmei(k,i)+berg(k,i))*deltat.ge.qsmall) then
              ! limit in-cloud values to 0.005 kg/kg
              qiic(k,i)=min(qi(k,i)/icldm(k,i),5.e-3_r8)
              niic(k,i)=max(ni(k,i)/icldm(k,i),0._r8)

              ! hm add 6/2/11 switch for specification of cloud ice number
              if (nicons) then
                 niic(k,i)=ninst/rho(k,i)
              end if
           else
              qiic(k,i)=0._r8
              niic(k,i)=0._r8

              if (do_cldice) then
                 cmei(k,i)=(-qi(k,i)/deltat-berg(k,i))*omsm
              end if
           end if

        end do pre_col_loop
     end do pre_vert_loop

     ! add to cme output
     cmeout = cmeout + cmei

     !=========================================================
     ! Main microphysical loop
     !=========================================================

     ! initialize precip fallspeeds to zero
     ums = 0._r8
     uns = 0._r8
     umr = 0._r8
     unr = 0._r8

     ! for sub-columns cldm has already been set to 1 if cloud
     ! water or ice is present, so cldmax will be correctly set below
     ! and nothing extra needs to be done here

     cldmax = cldm

     micro_vert_loop: do k=1,nlev

        ! calculate precip fraction based on maximum overlap assumption

        ! if rain or snow mix ratios are smaller than threshold, 
        ! then leave cldmax as cloud fraction at current level
        if (k /= 1) then
           where (qric(k-1,:).ge.qsmall .or. qsic(k-1,:).ge.qsmall)
              cldmax(k,:)=max(cldmax(k-1,:),cldmax(k,:))
           end where
        end if
        
        do i = 1, mgncol

           !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
           ! get size distribution parameters based on in-cloud cloud water
           ! these calculations also ensure consistency between number and mixing ratio
           !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

           ! cloud liquid
           !-------------------------------------------

           ! ".true." below turns on adjustment of ncic for consistency.
           call size_dist_param_liq(qcic(k,i), ncic(k,i), cdnl, rho(k,i), .true., &
                                    pgam(k,i), lamc(k,i))

           if (lamc(k,i) > 0._r8) then

              ! parameter to calculate droplet freezing
              cdist1(k,i) = ncic(k,i)/gamma(pgam(k,i)+1._r8)

           else
              cdist1(k,i) = 0._r8
           end if

        end do

        !========================================================================
        ! autoconversion of cloud liquid water to rain
        ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
        ! minimum qc of 1 x 10^-8 prevents floating point error

        call kk2000_liq_autoconversion(qcic(k,:), ncic(k,:), rho(k,:), relvar(k,:), &
                                       prc(k,:), nprc(k,:), nprc1(k,:))

        ! add autoconversion to precip from above to get provisional rain mixing ratio
        ! and number concentration (qric and nric)

        ! hm 11-16-11, modify, divide dz by 2 to use modified mid-point method
        ! This estimates rain and snow mass and number mixing ratios at
        ! mid-point to calculate process rates at mid-point, with final
        ! values of rain and snow mass and number mixing ratios calculated
        ! on interfaces

        if (k .eq. 1) then
           dum=0.45_r8

           qric(k,:)= prc(k,:)*lcldm(k,:)*dz(k,:)/2._r8/cldmax(k,:)/dum
           nric(k,:)=nprc(k,:)*lcldm(k,:)*dz(k,:)/2._r8/cldmax(k,:)/dum
        else

           ! no autoconversion of rain number if rain/snow falling from above
           ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
           ! by the existing rain/snow particles from above

           where (qric(k-1,:).ge.1.e-9_r8 .or. qsic(k-1,:).ge.1.e-9_r8)
              nprc(k,:) = 0._r8
           end where

           do i = 1,mgncol
              if (qric(k-1,i).ge.qsmall) then
                 dum=umr(k-1,i)
                 dum1=unr(k-1,i)
              else
                 ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)

                 dum=0.45_r8
                 dum1=0.45_r8
              end if

              qric(k,i) = (rho(k-1,i)*umr(k-1,i)*qric(k-1,i)*cldmax(k-1,i)+ &
                   (rho(k,i)*dz(k,i)/2._r8*((pra(k-1,i)+prc(k,i))*lcldm(k,i)+ &
                   (pre(k-1,i)-pracs(k-1,i)-mnuccr(k-1,i))*cldmax(k,i)))) &
                   /(dum*rho(k,i)*cldmax(k,i))
              nric(k,i) = (rho(k-1,i)*unr(k-1,i)*nric(k-1,i)*cldmax(k-1,i)+ &
                   (rho(k,i)*dz(k,i)/2._r8*(nprc(k,i)*lcldm(k,i)+ &
                   (nsubr(k-1,i)-npracs(k-1,i)-nnuccr(k-1,i)+nragg(k-1,i))*cldmax(k,i))))&
                   /(dum1*rho(k,i)*cldmax(k,i))
           end do

        end if

        ! if precip mix ratio is zero so should number concentration

        where (qric(k,:).lt.qsmall)
           qric(k,:)=0._r8
           nric(k,:)=0._r8
        end where

        ! make sure number concentration is a positive number to avoid 
        ! taking root of negative later

        nric(k,:)=max(nric(k,:),0._r8)

        ! Get size distribution parameters for cloud ice

        call size_dist_param_ice(qiic(k,:), niic(k,:), lami(k,:), n0i(k,:))

        !.......................................................................
        ! Autoconversion of cloud ice to snow
        ! similar to Ferrier (1994)

        if (do_cldice) then
           call ice_autoconversion(t(k,:), qiic(k,:), lami(k,:), n0i(k,:), &
                                   prci(k,:), nprci(k,:))
        else
           ! Add in the particles that we have already converted to snow, and
           ! don't do any further autoconversion of ice.
           prci(k,:)  = tnd_qsnow(k,:) / cldm(k,:)
           nprci(k,:) = tnd_nsnow(k,:) / cldm(k,:)
        end if

        do i=1,mgncol

           ! add autoconversion to flux from level above to get provisional snow mixing ratio
           ! and number concentration (qsic and nsic)

           ! hm 11-16-11 modify for mid-point method, see comments above

           if (k == 1) then
              dum=(asn(k,i)*cons25)
              qsic(k,i)=prci(k,i)*icldm(k,i)*dz(k,i)/2._r8/cldmax(k,i)/dum
              nsic(k,i)=nprci(k,i)*icldm(k,i)*dz(k,i)/2._r8/cldmax(k,i)/dum
           else
              if (qsic(k-1,i) >= qsmall) then
                 dum=ums(k-1,i)
                 dum1=uns(k-1,i)
              else
                 dum = asn(k,i)*cons25
                 dum1 = dum
              end if

              qsic(k,i) = (rho(k-1,i)*ums(k-1,i)*qsic(k-1,i)*cldmax(k-1,i)+ &
                   (rho(k,i)*dz(k,i)/2._r8*((prci(k,i)+prai(k-1,i)+psacws(k-1,i)+bergs(k-1,i))*icldm(k,i)+ &
                   (prds(k-1,i)+pracs(k-1,i)+mnuccr(k-1,i))*cldmax(k,i))))&
                   /(dum*rho(k,i)*cldmax(k,i))

              nsic(k,i) = (rho(k-1,i)*uns(k-1,i)*nsic(k-1,i)*cldmax(k-1,i)+ &
                   (rho(k,i)*dz(k,i)/2._r8*(nprci(k,i)*icldm(k,i)+ &
                   (nsubs(k-1,i)+nsagg(k-1,i)+nnuccr(k-1,i))*cldmax(k,i)))) &
                   /(dum1*rho(k,i)*cldmax(k,i))

           end if

        end do

        ! if precip mix ratio is zero so should number concentration
        where (qsic(k,:) < qsmall)
           qsic(k,:)=0._r8
           nsic(k,:)=0._r8
        end where

        ! make sure number concentration is a positive number to avoid 
        ! taking root of negative later

        nsic(k,:)=max(nsic(k,:),0._r8)

        !.......................................................................
        ! get size distribution parameters for precip
        !......................................................................
        ! rain

        call size_dist_param_rain(qric(k,:), nric(k,:), lamr(k,:), n0r(k,:))

        where (lamr(k,:) >= qsmall)

           ! provisional rain number and mass weighted mean fallspeed (m/s)

           unr(k,:) = min(arn(k,:)*cons4/lamr(k,:)**br,9.1_r8*rhof(k,:))
           umr(k,:) = min(arn(k,:)*cons5/(6._r8*lamr(k,:)**br),9.1_r8*rhof(k,:))

        elsewhere
           umr(k,:) = 0._r8
           unr(k,:) = 0._r8
        end where

        !......................................................................
        ! snow

        call size_dist_param_snow(qsic(k,:), nsic(k,:), lams(k,:), n0s(k,:))

        where (lams(k,:) > 0._r8)

           ! provisional snow number and mass weighted mean fallspeed (m/s)

           ums(k,:) = min(asn(k,:)*cons8/(6._r8*lams(k,:)**bs),1.2_r8*rhof(k,:))
           uns(k,:) = min(asn(k,:)*cons7/lams(k,:)**bs,1.2_r8*rhof(k,:))

        elsewhere
           ums(k,:) = 0._r8
           uns(k,:) = 0._r8
        end where

        if (do_cldice) then

           ! heterogeneous freezing of cloud water
           !----------------------------------------------

           call immersion_freezing(t(k,:), pgam(k,:), lamc(k,:), cdist1(k,:), qcic(k,:), relvar(k,:),  &
                                   mnuccc(k,:), nnuccc(k,:))

           ! make sure number of droplets frozen does not exceed available ice nuclei concentration
           ! this prevents 'runaway' droplet freezing


           where (qcic(k,:).ge.qsmall .and. t(k,:).lt.269.15_r8)
              where (nnuccc(k,:)*lcldm(k,:).gt.nnuccd(k,:))
                 ! scale mixing ratio of droplet freezing with limit
                 mnuccc(k,:)=mnuccc(k,:)*(nnuccd(k,:)/(nnuccc(k,:)*lcldm(k,:)))
                 nnuccc(k,:)=nnuccd(k,:)/lcldm(k,:)
              end where
           end where

           call contact_freezing(t(k,:), p(k,:), rndst(:,k,:), nacon(:,k,:), pgam(k,:),  &
                                 lamc(k,:), cdist1(k,:), qcic(k,:), relvar(k,:),         &
                                 mnucct(k,:), nnucct(k,:))

        else
           mnuccc(k,:)=0._r8
           nnuccc(k,:)=0._r8
           mnucct(k,:)=0._r8
           nnucct(k,:)=0._r8
        end if

        call snow_self_aggregation(t(k,:), rho(k,:), asn(k,:), qsic(k,:), nsic(k,:), &
                                   nsagg(k,:))

        call accrete_cloud_water_snow(t(k,:), rho(k,:), asn(k,:), uns(k,:), mu(k,:), &
                                      qcic(k,:), ncic(k,:), qsic(k,:), pgam(k,:),    &
                                      lamc(k,:), lams(k,:), n0s(k,:),                &
                                      psacws(k,:), npsacws(k,:))

        if (do_cldice) then
           call secondary_ice_production(t(k,:), psacws(k,:), msacwi(k,:), nsacwi(k,:))
        else
           nsacwi(k,:) = 0.0_r8
           msacwi(k,:) = 0.0_r8
        end if

        call accrete_rain_snow(t(k,:), rho(k,:), umr(k,:), ums(k,:), unr(k,:), uns(k,:),        &
                               qric(k,:), qsic(k,:), lamr(k,:), n0r(k,:), lams(k,:), n0s(k,:),  &
                               pracs(k,:), npracs(k,:))

        call heterogeneous_rain_freezing(t(k,:), qric(k,:), nric(k,:), lamr(k,:), &
                                         mnuccr(k,:), nnuccr(k,:))

        call accrete_cloud_water_rain(qric(k,:), qcic(k,:), ncic(k,:), &
                                      relvar(k,:), accre_enhan(k,:),   &
                                      pra(k,:), npra(k,:))

        call self_collection_rain(rho(k,:), qric(k,:), nric(k,:), nragg(k,:))

        if (do_cldice) then
           call accrete_cloud_ice_snow(t(k,:), rho(k,:), asn(k,:), qiic(k,:), niic(k,:), qsic(k,:), &
                                       lams(k,:), n0s(k,:), prai(k,:), nprai(k,:))
        else
           prai(k,:) = 0._r8
           nprai(k,:) = 0._r8
        end if
              
        call evaporate_sublimate_precip(deltat, t(k,:), p(k,:), rho(k,:), dv(k,:), mu(k,:), sc(k,:),        &
                                        q(k,:), qvl(k,:), qvi(k,:), lcldm(k,:), cldmax(k,:),                &
                                        arn(k,:), asn(k,:), qcic(k,:), qiic(k,:), qric(k,:), qsic(k,:),     &
                                        lamr(k,:), n0r(k,:), lams(k,:), n0s(k,:), cmei(k,:),                &
                                        pre(k,:), prds(k,:))

        call bergeron_process(t(k,:), rho(k,:), dv(k,:), mu(k,:), sc(k,:), qvl(k,:), qvi(k,:), asn(k,:),    &
                              qcic(k,:), qsic(k,:), lams(k,:), n0s(k,:), bergs(k,:))

        ! Big "administration" loop enforces conservation, updates variables
        ! that accumulate over substeps, and sets output variables.
        do i=1,mgncol

           ! conservation to ensure no negative values of cloud water/precipitation
           ! in case microphysical process rates are large
           !===================================================================

           ! make sure to use end-of-time step values for cloud water, ice, due
           ! condensation/deposition

           ! note: for check on conservation, processes are multiplied by omsm
           ! to prevent problems due to round off error

           qce=(qc(k,i) - berg(k,i)*deltat)
           nce=nc(k,i)
           qie=(qi(k,i)+(cmei(k,i)+berg(k,i))*deltat)
           nie=(ni(k,i)+nnuccd(k,i)*deltat)

           ! conservation of qc
           !-------------------------------------------------------------------

           dum = (prc(k,i)+pra(k,i)+mnuccc(k,i)+mnucct(k,i)+msacwi(k,i)+ &
                psacws(k,i)+bergs(k,i))*lcldm(k,i)*deltat

           if (dum.gt.qce) then
              ratio = qce/deltat/lcldm(k,i)/(prc(k,i)+pra(k,i)+mnuccc(k,i)+mnucct(k,i)+msacwi(k,i)+psacws(k,i)+bergs(k,i))*omsm 

              prc(k,i) = prc(k,i)*ratio
              pra(k,i) = pra(k,i)*ratio
              mnuccc(k,i) = mnuccc(k,i)*ratio
              mnucct(k,i) = mnucct(k,i)*ratio  
              msacwi(k,i) = msacwi(k,i)*ratio  
              psacws(k,i) = psacws(k,i)*ratio
              bergs(k,i) = bergs(k,i)*ratio
           end if

           ! conservation of nc
           !-------------------------------------------------------------------
           dum = (nprc1(k,i)+npra(k,i)+nnuccc(k,i)+nnucct(k,i)+ &
                npsacws(k,i)-nsubc(k,i))*lcldm(k,i)*deltat

           if (dum.gt.nce) then
              ratio = nce/deltat/((nprc1(k,i)+npra(k,i)+nnuccc(k,i)+nnucct(k,i)+&
                   npsacws(k,i)-nsubc(k,i))*lcldm(k,i))*omsm

              nprc1(k,i) = nprc1(k,i)*ratio
              npra(k,i) = npra(k,i)*ratio
              nnuccc(k,i) = nnuccc(k,i)*ratio
              nnucct(k,i) = nnucct(k,i)*ratio
              npsacws(k,i) = npsacws(k,i)*ratio
              nsubc(k,i)=nsubc(k,i)*ratio
           end if

           if (do_cldice) then

              ! conservation of qi
              !-------------------------------------------------------------------
              dum = ((-mnuccc(k,i)-mnucct(k,i)-msacwi(k,i))*lcldm(k,i)+(prci(k,i)+ &
                   prai(k,i))*icldm(k,i))*deltat

              if (dum.gt.qie) then

                 ratio = (qie/deltat+(mnuccc(k,i)+mnucct(k,i)+msacwi(k,i))*lcldm(k,i))/ &
                      ((prci(k,i)+prai(k,i))*icldm(k,i))*omsm
                 prci(k,i) = prci(k,i)*ratio
                 prai(k,i) = prai(k,i)*ratio
              end if

              ! conservation of ni
              !-------------------------------------------------------------------
              dum = ((-nnucct(k,i)-nsacwi(k,i))*lcldm(k,i)+(nprci(k,i)+ &
                   nprai(k,i)-nsubi(k,i))*icldm(k,i))*deltat

              if (dum.gt.nie) then

                 ratio = (nie/deltat+(nnucct(k,i)+nsacwi(k,i))*lcldm(k,i))/ &  
                      ((nprci(k,i)+nprai(k,i)-nsubi(k,i))*icldm(k,i))*omsm
                 nprci(k,i) = nprci(k,i)*ratio
                 nprai(k,i) = nprai(k,i)*ratio
                 nsubi(k,i) = nsubi(k,i)*ratio
              end if
           end if

           ! for precipitation conservation, use logic that vertical integral 
           ! of tendency from current level to top of model (i.e., qrtot) cannot be negative

           ! conservation of rain mixing rat
           !-------------------------------------------------------------------
           if (((prc(k,i)+pra(k,i))*lcldm(k,i)+(-mnuccr(k,i)+pre(k,i)-pracs(k,i))*&
                cldmax(k,i))*dz(k,i)*rho(k,i)+qrtot(i).lt.0._r8) then

              if (-pre(k,i)+pracs(k,i)+mnuccr(k,i).ge.qsmall) then

                 ratio = (qrtot(i)/(dz(k,i)*rho(k,i))+(prc(k,i)+pra(k,i))*lcldm(k,i))/&
                      ((-pre(k,i)+pracs(k,i)+mnuccr(k,i))*cldmax(k,i))*omsm 

                 pre(k,i) = pre(k,i)*ratio
                 pracs(k,i) = pracs(k,i)*ratio
                 mnuccr(k,i) = mnuccr(k,i)*ratio
              end if
           end if

           ! conservation of nr
           !-------------------------------------------------------------------
           ! for now neglect evaporation of nr
           nsubr(k,i)=0._r8

           if ((nprc(k,i)*lcldm(k,i)+(-nnuccr(k,i)+nsubr(k,i)-npracs(k,i)&
                +nragg(k,i))*cldmax(k,i))*dz(k,i)*rho(k,i)+nrtot(i).lt.0._r8) then

              if (-nsubr(k,i)-nragg(k,i)+npracs(k,i)+nnuccr(k,i).ge.qsmall) then

                 ratio = (nrtot(i)/(dz(k,i)*rho(k,i))+nprc(k,i)*lcldm(k,i))/&
                      ((-nsubr(k,i)-nragg(k,i)+npracs(k,i)+nnuccr(k,i))*cldmax(k,i))*omsm

                 nsubr(k,i) = nsubr(k,i)*ratio
                 npracs(k,i) = npracs(k,i)*ratio
                 nnuccr(k,i) = nnuccr(k,i)*ratio
                 nragg(k,i) = nragg(k,i)*ratio
              end if
           end if

           ! conservation of snow mix ratio
           !-------------------------------------------------------------------
           if (((bergs(k,i)+psacws(k,i))*lcldm(k,i)+(prai(k,i)+prci(k,i))*icldm(k,i)+(pracs(k,i)+&
                mnuccr(k,i)+prds(k,i))*cldmax(k,i))*dz(k,i)*rho(k,i)+qstot(i).lt.0._r8) then

              if (-prds(k,i).ge.qsmall) then

                 ratio = (qstot(i)/(dz(k,i)*rho(k,i))+(bergs(k,i)+psacws(k,i))*lcldm(k,i)+(prai(k,i)+prci(k,i))*icldm(k,i)+&
                      (pracs(k,i)+mnuccr(k,i))*cldmax(k,i))/(-prds(k,i)*cldmax(k,i))*omsm

                 prds(k,i) = prds(k,i)*ratio
              else
                 prds(k,i) = 0._r8
              end if
           end if

           ! conservation of ns
           !-------------------------------------------------------------------
           ! calculate loss of number due to sublimation
           ! for now neglect sublimation of ns
           nsubs(k,i)=0._r8

           if ((nprci(k,i)*icldm(k,i)+(nnuccr(k,i)+nsubs(k,i)+nsagg(k,i))*cldmax(k,i))*&
                dz(k,i)*rho(k,i)+nstot(i).lt.0._r8) then

              if (-nsubs(k,i)-nsagg(k,i).ge.qsmall) then

                 ratio = (nstot(i)/(dz(k,i)*rho(k,i))+nprci(k,i)*icldm(k,i)+&
                      nnuccr(k,i)*cldmax(k,i))/((-nsubs(k,i)-nsagg(k,i))*cldmax(k,i))*omsm

                 nsubs(k,i) = nsubs(k,i)*ratio
                 nsagg(k,i) = nsagg(k,i)*ratio
              end if
           end if

           ! get tendencies due to microphysical conversion processes
           !==========================================================
           ! note: tendencies are multiplied by appropriate cloud/precip 
           ! fraction to get grid-scale values
           ! note: cmei is already grid-average values

           qvlat(k,i) = qvlat(k,i)-(pre(k,i)+prds(k,i))*cldmax(k,i)-cmei(k,i) 

           tlat(k,i) = tlat(k,i)+((pre(k,i)*cldmax(k,i)) &
                *xxlv+(prds(k,i)*cldmax(k,i)+cmei(k,i))*xxls+ &
                ((bergs(k,i)+psacws(k,i)+mnuccc(k,i)+mnucct(k,i)+msacwi(k,i))*lcldm(k,i)+(mnuccr(k,i)+ &
                pracs(k,i))*cldmax(k,i)+berg(k,i))*xlf)

           qctend(k,i) = qctend(k,i)+ &
                (-pra(k,i)-prc(k,i)-mnuccc(k,i)-mnucct(k,i)-msacwi(k,i)- & 
                psacws(k,i)-bergs(k,i))*lcldm(k,i)-berg(k,i)

           if (do_cldice) then
              qitend(k,i) = qitend(k,i)+ &
                   (mnuccc(k,i)+mnucct(k,i)+msacwi(k,i))*lcldm(k,i)+(-prci(k,i)- & 
                   prai(k,i))*icldm(k,i)+cmei(k,i)+berg(k,i)
           end if

           qrtend(k,i) = qrtend(k,i)+ &
                (pra(k,i)+prc(k,i))*lcldm(k,i)+(pre(k,i)-pracs(k,i)- &
                mnuccr(k,i))*cldmax(k,i)

           qstend(k,i) = qstend(k,i)+ &
                (prai(k,i)+prci(k,i))*icldm(k,i)+(psacws(k,i)+bergs(k,i))*lcldm(k,i)+(prds(k,i)+ &
                pracs(k,i)+mnuccr(k,i))*cldmax(k,i)

           ! add output for cmei (accumulate)
           cmeitot(k,i) = cmeitot(k,i) + cmei(k,i)

           ! assign variables for trop_mozart, these are grid-average
           !-------------------------------------------------------------------
           ! evaporation/sublimation is stored here as positive term

           evapsnow(k,i) = evapsnow(k,i)-prds(k,i)*cldmax(k,i)
           nevapr(k,i) = nevapr(k,i)-pre(k,i)*cldmax(k,i)

           ! change to make sure prain is positive: do not remove snow from
           ! prain used for wet deposition
           prain(k,i) = prain(k,i)+(pra(k,i)+prc(k,i))*lcldm(k,i)+(-pracs(k,i)- &
                mnuccr(k,i))*cldmax(k,i)
           prodsnow(k,i) = prodsnow(k,i)+(prai(k,i)+prci(k,i))*icldm(k,i)+(psacws(k,i)+bergs(k,i))*lcldm(k,i)+(&
                pracs(k,i)+mnuccr(k,i))*cldmax(k,i)

           ! following are used to calculate 1st order conversion rate of cloud water
           !    to rain and snow (1/s), for later use in aerosol wet removal routine
           ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
           !    used to calculate pra, prc, ... in this routine
           ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
           !                      (no cloud ice or bergeron terms)
           ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }

           qcsinksum_rate1ord(k,i) = qcsinksum_rate1ord(k,i) + (pra(k,i)+prc(k,i)+psacws(k,i))*lcldm(k,i) 
           qcsum_rate1ord(k,i) = qcsum_rate1ord(k,i) + qc(k,i) 

           ! microphysics output, note this is grid-averaged
           pratot(k,i)=pratot(k,i)+pra(k,i)*lcldm(k,i)
           prctot(k,i)=prctot(k,i)+prc(k,i)*lcldm(k,i)
           mnuccctot(k,i)=mnuccctot(k,i)+mnuccc(k,i)*lcldm(k,i)
           mnuccttot(k,i)=mnuccttot(k,i)+mnucct(k,i)*lcldm(k,i)
           mnuccdtot(k,i)=mnuccdtot(k,i)+mnuccd(k,i)*lcldm(k,i)
           msacwitot(k,i)=msacwitot(k,i)+msacwi(k,i)*lcldm(k,i)
           psacwstot(k,i)=psacwstot(k,i)+psacws(k,i)*lcldm(k,i)
           bergstot(k,i)=bergstot(k,i)+bergs(k,i)*lcldm(k,i)
           bergtot(k,i)=bergtot(k,i)+berg(k,i)
           prcitot(k,i)=prcitot(k,i)+prci(k,i)*icldm(k,i)
           praitot(k,i)=praitot(k,i)+prai(k,i)*icldm(k,i)
           mnuccrtot(k,i)=mnuccrtot(k,i)+mnuccr(k,i)*cldmax(k,i)
           pracstot(k,i)=pracstot(k,i)+pracs(k,i)*cldmax(k,i)

           nctend(k,i) = nctend(k,i)+&
                (-nnuccc(k,i)-nnucct(k,i)-npsacws(k,i)+nsubc(k,i) & 
                -npra(k,i)-nprc1(k,i))*lcldm(k,i)

           if (do_cldice) then
              nitend(k,i) = nitend(k,i)+ nnuccd(k,i)+ &
                   (nnucct(k,i)+nsacwi(k,i))*lcldm(k,i)+(nsubi(k,i)-nprci(k,i)- &
                   nprai(k,i))*icldm(k,i)
           end if

           nstend(k,i) = nstend(k,i)+(nsubs(k,i)+ &
                nsagg(k,i)+nnuccr(k,i))*cldmax(k,i)+nprci(k,i)*icldm(k,i)

           nrtend(k,i) = nrtend(k,i)+ &
                nprc(k,i)*lcldm(k,i)+(nsubr(k,i)-npracs(k,i)-nnuccr(k,i) &
                +nragg(k,i))*cldmax(k,i)

           ! make sure that ni at advanced time step does not exceed
           ! maximum (existing N + source terms*dt), which is possible if mtime < deltat
           ! note that currently mtime = deltat
           !================================================================

           if (do_cldice .and. nitend(k,i).gt.0._r8.and.ni(k,i)+nitend(k,i)*deltat.gt.nimax(k,i)) then
              nitend(k,i)=max(0._r8,(nimax(k,i)-ni(k,i))/deltat)
           end if

        end do
           
        ! End of "administration" loop

        ! get final values for precipitation q and N, based on
        ! flux of precip from above, source/sink term, and terminal fallspeed
        ! see eq. 15-16 in MG2008

        do i = 1, mgncol

           ! rain

           if (qric(k,i).ge.qsmall) then
              if (k .eq. 1) then
                 qric(k,i)=qrtend(k,i)*dz(k,i)/cldmax(k,i)/umr(k,i)
                 nric(k,i)=nrtend(k,i)*dz(k,i)/cldmax(k,i)/unr(k,i)
              else
                 qric(k,i) = (rho(k-1,i)*umr(k-1,i)*qric(k-1,i)*cldmax(k-1,i)+ &
                      (rho(k,i)*dz(k,i)*qrtend(k,i)))/(umr(k,i)*rho(k,i)*cldmax(k,i))
                 nric(k,i) = (rho(k-1,i)*unr(k-1,i)*nric(k-1,i)*cldmax(k-1,i)+ &
                      (rho(k,i)*dz(k,i)*nrtend(k,i)))/(unr(k,i)*rho(k,i)*cldmax(k,i))

              end if
           else
              qric(k,i)=0._r8
              nric(k,i)=0._r8
           end if

           ! snow

           if (qsic(k,i).ge.qsmall) then
              if (k .eq. 1) then
                 qsic(k,i)=qstend(k,i)*dz(k,i)/cldmax(k,i)/ums(k,i)
                 nsic(k,i)=nstend(k,i)*dz(k,i)/cldmax(k,i)/uns(k,i)
              else
                 qsic(k,i) = (rho(k-1,i)*ums(k-1,i)*qsic(k-1,i)*cldmax(k-1,i)+ &
                      (rho(k,i)*dz(k,i)*qstend(k,i)))/(ums(k,i)*rho(k,i)*cldmax(k,i))
                 nsic(k,i) = (rho(k-1,i)*uns(k-1,i)*nsic(k-1,i)*cldmax(k-1,i)+ &
                      (rho(k,i)*dz(k,i)*nstend(k,i)))/(uns(k,i)*rho(k,i)*cldmax(k,i))
              end if
           else
              qsic(k,i)=0._r8
              nsic(k,i)=0._r8
           end if

           ! calculate precipitation flux at surface
           !=========================================================
           ! divide by density of water to get units of m/s

           prect(i) = prect(i)+(qrtend(k,i)*dz(k,i)*rho(k,i)+&
                qstend(k,i)*dz(k,i)*rho(k,i))/rhow
           preci(i) = preci(i)+qstend(k,i)*dz(k,i)*rho(k,i)/rhow

           ! convert rain rate from m/s to mm/hr

           rainrt(k,i)=rainrt(k,i) + (qric(k,i)*rho(k,i)*umr(k,i)/rhow*3600._r8*1000._r8)

           ! vertically-integrated precip source/sink terms (note: grid-averaged)

           qrtot(i) = max(qrtot(i)+qrtend(k,i)*dz(k,i)*rho(k,i),0._r8)
           qstot(i) = max(qstot(i)+qstend(k,i)*dz(k,i)*rho(k,i),0._r8)
           nrtot(i) = max(nrtot(i)+nrtend(k,i)*dz(k,i)*rho(k,i),0._r8)
           nstot(i) = max(nstot(i)+nstend(k,i)*dz(k,i)*rho(k,i),0._r8)

           ! calculate melting and freezing of precip
           !=========================================================

           ! melt snow at +2 C

           ! Difference in amount of heat between temperature at which
           ! all snow melts, and the current state.
           dum = (snowmelt - t(k,i) - tlat(k,i)/cpp*deltat) * cpp

           ! Test if temperature is above threshold and snow is present.
           if (dum < 0._r8 .and. qstot(i) > 0._r8) then

              ! (negative) heat produced if all snow is melted
              dum1 = -xlf * qstot(i)/(dz(k,i)*rho(k,i))*deltat

              ! ratio of heating needed to get to the threshold, to
              ! the total heating from melting everything, is equal
              ! to the proportion of snow that actually melts
              dum = min(1._r8, dum/dum1)
              dum = max(0._r8, dum)

              ! Melt snow
              qric(k,i)=qric(k,i)+dum*qsic(k,i)
              nric(k,i)=nric(k,i)+dum*nsic(k,i)
              qsic(k,i)=(1._r8-dum)*qsic(k,i)
              nsic(k,i)=(1._r8-dum)*nsic(k,i)

              qrtot(i)=qrtot(i)+dum*qstot(i)
              nrtot(i)=nrtot(i)+dum*nstot(i)
              qstot(i)=(1._r8-dum)*qstot(i)
              nstot(i)=(1._r8-dum)*nstot(i)

              preci(i)=(1._r8-dum)*preci(i)

              ! Get heating tendency based on proportion of snow that
              ! actually melts.
              dum1 = dum * dum1/deltat

              meltsdttot(k,i)=meltsdttot(k,i) + dum1
              tlat(k,i)=tlat(k,i)+dum1

           end if

           ! freeze all rain at -5C for Arctic
           !=========================================================

           ! Difference in amount of heat between temperature at which
           ! all rain freezes, and the current state.
           dum = (rainfrze - t(k,i) - tlat(k,i)/cpp*deltat) * cpp

           ! Test if temperature is below threshold and snow is present.
           if (dum > 0._r8 .and. qrtot(i) > 0._r8) then

              ! heat produced if all rain freezes
              dum1 = xlf * qrtot(i)/(dz(k,i)*rho(k,i))*deltat

              ! ratio of heating needed to get to the threshold, to
              ! the total heating from freezing everything, is equal
              ! to the proportion of rain that actually freezes
              dum = min(1._r8, dum/dum1)
              dum = max(0._r8, dum)

              ! Freeze rain
              qsic(k,i)=qsic(k,i)+dum*qric(k,i)
              nsic(k,i)=nsic(k,i)+dum*nric(k,i)
              qric(k,i)=(1._r8-dum)*qric(k,i)
              nric(k,i)=(1._r8-dum)*nric(k,i)

              qstot(i)=qstot(i)+dum*qrtot(i)
              nstot(i)=nstot(i)+dum*nrtot(i)
              qrtot(i)=(1._r8-dum)*qrtot(i)
              nrtot(i)=(1._r8-dum)*nrtot(i)

              preci(i)=preci(i)+dum*(prect(i)-preci(i))

              ! Get heating tendency based on proportion of rain that
              ! actually freezes.
              dum1 = dum * dum1/deltat

              frzrdttot(k,i)=frzrdttot(k,i) + dum1
              tlat(k,i)=tlat(k,i)+dum1

           end if

           ! if rain/snow mix ratio is zero so should number concentration
           !=========================================================
           
           if (qsic(k,i) < qsmall) then
              qsic(k,i)=0._r8
              nsic(k,i)=0._r8
           end if

           if (qric(k,i) < qsmall) then
              qric(k,i)=0._r8
              nric(k,i)=0._r8
           end if

           ! make sure number concentration is a positive number to avoid 
           ! taking root of negative

           nric(k,i)=max(nric(k,i),0._r8)
           nsic(k,i)=max(nsic(k,i),0._r8)

           ! get size distribution parameters for fallspeed calculations
           !=========================================================

           ! rain

           call size_dist_param_rain(qric(k,i), nric(k,i), lamr(k,i), n0r(k,i))

           if (lamr(k,i).ge.qsmall) then

              ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

              unr(k,i) = min(arn(k,i)*cons4/lamr(k,i)**br,9.1_r8*rhof(k,i))
              umr(k,i) = min(arn(k,i)*cons5/(6._r8*lamr(k,i)**br),9.1_r8*rhof(k,i))

           else
              umr(k,i)=0._r8
              unr(k,i)=0._r8
           end if

           !......................................................................
           ! snow

           call size_dist_param_snow(qsic(k,i), nsic(k,i), lams(k,i), n0s(k,i))

           if (lams(k,i) > 0._r8) then

              ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

              ums(k,i) = min(asn(k,i)*cons8/(6._r8*lams(k,i)**bs),1.2_r8*rhof(k,i))
              uns(k,i) = min(asn(k,i)*cons7/lams(k,i)**bs,1.2_r8*rhof(k,i))

           else
              ums(k,i) = 0._r8
              uns(k,i) = 0._r8
           end if

        end do

        ! Done with vertical dependencies from precipitation.

     end do micro_vert_loop ! end k loop

     ! sum over sub-step for average process rates
     !-----------------------------------------------------
     ! convert rain/snow q and N for output to history, note, 
     ! output is for gridbox average

     ! calculate precip fluxes and adding them to summing sub-stepping variables
     ! calculate the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
     ! ---------------------------------------------------------------------

     rflx(2:,:)  = rflx(2:,:) + (qric*rho*umr*cldmax)

     dumr(1,:)   = int_to_mid(1,:) * qric(1,:)
     qrout(1,:)  = qrout(1,:)  + int_to_mid(1,:)*qric(1,:)*cldmax(1,:)

     dumr(2:,:)  = interp_to_mid(qric,int_to_mid(2:,:))
     qrout(2:,:) = qrout(2:,:) + interp_to_mid(qric * cldmax, int_to_mid(2:,:))

     dumnr(1,:)  = int_to_mid(1,:) * nric(1,:)
     nrout(1,:)  = nrout(1,:)  + &
          (int_to_mid(1,:)*nric(1,:)*cldmax(1,:)*rho(1,:))

     dumnr(2:,:) = interp_to_mid(nric,int_to_mid(2:,:))
     nrout(2:,:) = nrout(2:,:) + interp_to_mid(nric * cldmax * rho, int_to_mid(2:,:))
     ! Calculate rercld

     ! calculate mean size of combined rain and cloud water
 
     ! hm 11-22-11 modify to interpolate rain from interface to mid-point
     ! logic is to interpolate rain mass and number, then recalculate PSD
     ! parameters to get relevant parameters for mean size

     ! interpolate rain mass and number, store in dummy variables
     
     ! calculate n0r and lamr from interpolated mid-point rain mass and number
     ! divide by precip fraction to get in-precip (local) values of
     ! rain mass and number, divide by rhow to get rain number in kg^-1

     call size_dist_param_rain(dumr, dumnr, lamr, n0r)

     call calc_rercld(lamr, n0r, lamc, cdist1, pgam, dumr, qcic, &
                      arcld, rercld)

     nsout(1,:)  = nsout(1,:)  + &
          (int_to_mid(1,:)*nsic(1,:)*cldmax(1,:)*rho(1,:))
     nsout(2:,:) = nsout(2:,:) + interp_to_mid(nsic * cldmax * rho, int_to_mid(2:,:))

     qsout(1,:)  = qsout(1,:)  + int_to_mid(1,:)*qsic(1,:)*cldmax(1,:)
     qsout(2:,:) = qsout(2:,:) + interp_to_mid(qsic * cldmax, int_to_mid(2:,:))

     sflx(2:,:) = sflx(2:,:) + (qsic*rho*ums*cldmax)

     ! Sum into other variables that accumulate over substeps.
     tlat1 = tlat1 + tlat
     t = t + tlat*deltat/cpp

     qvlat1 = qvlat1 + qvlat
     q = q + qvlat*deltat

     qctend1 = qctend1 + qctend
     qc = qc + qctend*deltat

     qitend1 = qitend1 + qitend
     qi = qi + qitend*deltat

     nctend1 = nctend1 + nctend
     nc = nc + nctend*deltat

     nitend1 = nitend1 + nitend
     ni = ni + nitend*deltat

     prect1 = prect1 + prect
     preci1 = preci1 + preci

  end do substepping ! it loop, sub-step

  ! divide rain radius over substeps for average
  where (arcld > 0) rercld = rercld/arcld

  ! convert dt from sub-step back to full time step
  !-------------------------------------------------------------------
  deltat = deltatin

  ! assign variables back to start-of-timestep values before updating after sub-steps 
  !================================================================================
  
  call pack_array(qn, mgcols, top_lev, q)
  call pack_array(tn, mgcols, top_lev, t)
  call pack_array(qcn, mgcols, top_lev, qc)
  call pack_array(ncn, mgcols, top_lev, nc)
  call pack_array(qin, mgcols, top_lev, qi)
  call pack_array(nin, mgcols, top_lev, ni)

  !.............................................................................

  ! divide precip rate by number of sub-steps to get average over time step
  
  prect = prect1/real(iter)
  preci = preci1/real(iter)

  ! divide microphysical tendencies by number of sub-steps to get average over time step
  !================================================================================

  tlat = tlat1/real(iter)
  qvlat = qvlat1/real(iter)
  qctend = qctend1/real(iter)
  qitend = qitend1/real(iter)
  nctend = nctend1/real(iter)
  nitend = nitend1/real(iter)

  ! Re-apply droplet activation tendency
  nctend = nctend + npccn

  rainrt = rainrt/real(iter)

  ! divide by number of sub-steps to find final values
  rflx = rflx/real(iter)
  sflx = sflx/real(iter)

  ! divide output precip q and N by number of sub-steps to get average over time step
  !================================================================================

  qrout = qrout/real(iter)
  qsout = qsout/real(iter)
  nrout = nrout/real(iter)
  nsout = nsout/real(iter)

  ! divide trop_mozart variables by number of sub-steps to get average over time step 
  !================================================================================

  nevapr = nevapr/real(iter)
  evapsnow = evapsnow/real(iter)
  prain = prain/real(iter)
  prodsnow = prodsnow/real(iter)

  ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
  nevapr = nevapr + evapsnow
  prain = prain + prodsnow

  cmeout = cmeout/real(iter)

  cmeitot = cmeitot/real(iter)
  meltsdttot = meltsdttot/real(iter)
  frzrdttot  = frzrdttot /real(iter)

  ! microphysics output
  pratot=pratot/real(iter)
  prctot=prctot/real(iter)
  mnuccctot=mnuccctot/real(iter)
  mnuccttot=mnuccttot/real(iter)
  msacwitot=msacwitot/real(iter)
  psacwstot=psacwstot/real(iter)
  bergstot=bergstot/real(iter)
  bergtot=bergtot/real(iter)
  prcitot=prcitot/real(iter)
  praitot=praitot/real(iter)

  mnuccrtot=mnuccrtot/real(iter)
  pracstot =pracstot /real(iter)

  mnuccdtot=mnuccdtot/real(iter)

  sed_col_loop: do i=1,mgncol

     do k=1,nlev

        ! calculate sedimentation for cloud water and ice
        !================================================================================

        ! update in-cloud cloud mixing ratio and number concentration 
        ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
        ! note: these are in-cloud values***, hence we divide by cloud fraction

        dumc(k,i) = (qc(k,i)+qctend(k,i)*deltat)/lcldm(k,i)
        dumi(k,i) = (qi(k,i)+qitend(k,i)*deltat)/icldm(k,i)
        dumnc(k,i) = max((nc(k,i)+nctend(k,i)*deltat)/lcldm(k,i),0._r8)
        dumni(k,i) = max((ni(k,i)+nitend(k,i)*deltat)/icldm(k,i),0._r8)

        ! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(k,i)=ncnst/rho(k,i)
        end if

        ! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
           dumni(k,i)=ninst/rho(k,i)
        end if

        ! obtain new slope parameter to avoid possible singularity

        call size_dist_param_ice(dumi(k,i), dumni(k,i), lami(k,i), n0i(k,i))

        call size_dist_param_liq(dumc(k,i), dumnc(k,i), cdnl, rho(k,i), .true., &
                                 pgam(k,i), lamc(k,i))

        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------


        if (dumc(k,i).ge.qsmall) then

           vtrmc(k,i)=acn(k,i)*gamma(4._r8+bc+pgam(k,i))/ &
                (lamc(k,i)**bc*gamma(pgam(k,i)+4._r8))

           fc(k) = g*rho(k,i)*vtrmc(k,i)

           fnc(k) = g*rho(k,i)* &
                acn(k,i)*gamma(1._r8+bc+pgam(k,i))/ &
                (lamc(k,i)**bc*gamma(pgam(k,i)+1._r8))
        else
           fc(k) = 0._r8
           fnc(k)= 0._r8
        end if

        ! calculate number and mass weighted fall velocity for cloud ice

        if (dumi(k,i).ge.qsmall) then

           vtrmi(k,i)=min(ain(k,i)*cons17/(6._r8*lami(k,i)**bi), &
                1.2_r8*rhof(k,i))

           fi(k) = g*rho(k,i)*vtrmi(k,i)
           fni(k) = g*rho(k,i)* &
                min(ain(k,i)*cons16/lami(k,i)**bi,1.2_r8*rhof(k,i))
        else
           fi(k) = 0._r8
           fni(k)= 0._r8
        end if

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation

        dumc(k,i) = (qc(k,i)+qctend(k,i)*deltat)
        dumi(k,i) = (qi(k,i)+qitend(k,i)*deltat)
        dumnc(k,i) = max((nc(k,i)+nctend(k,i)*deltat),0._r8)
        dumni(k,i) = max((ni(k,i)+nitend(k,i)*deltat),0._r8)

        if (dumc(k,i).lt.qsmall) dumnc(k,i)=0._r8
        if (dumi(k,i).lt.qsmall) dumni(k,i)=0._r8

     end do       !!! vertical loop

     ! initialize nstep for sedimentation sub-steps

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + int(max( &
          maxval( fi/pdel(:,i)), &
          maxval( fc/pdel(:,i)), &
          maxval(fni/pdel(:,i)), &
          maxval(fnc/pdel(:,i))) &
          * deltat)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     do n = 1,nstep

        if (do_cldice) then
           falouti  = fi  * dumi(:,i)
           faloutni = fni * dumni(:,i)
        else
           falouti  = 0._r8
           faloutni = 0._r8
        end if
        faloutc  = fc  * dumc(:,i)
        faloutnc = fnc * dumnc(:,i)

        ! top of model

        k = 1
        faltndi = falouti(k)/pdel(k,i)
        faltndni = faloutni(k)/pdel(k,i)
        faltndc = faloutc(k)/pdel(k,i)
        faltndnc = faloutnc(k)/pdel(k,i)

        ! add fallout terms to microphysical tendencies

        qitend(k,i) = qitend(k,i)-faltndi/nstep
        nitend(k,i) = nitend(k,i)-faltndni/nstep
        qctend(k,i) = qctend(k,i)-faltndc/nstep
        nctend(k,i) = nctend(k,i)-faltndnc/nstep

        ! sedimentation tendencies for output
        qcsedten(k,i)=qcsedten(k,i)-faltndc/nstep
        qisedten(k,i)=qisedten(k,i)-faltndi/nstep

        dumi(k,i) = dumi(k,i)-faltndi*deltat/nstep
        dumni(k,i) = dumni(k,i)-faltndni*deltat/nstep
        dumc(k,i) = dumc(k,i)-faltndc*deltat/nstep
        dumnc(k,i) = dumnc(k,i)-faltndnc*deltat/nstep

        do k = 2,nlev

           ! for cloud liquid and ice, if cloud fraction increases with height
           ! then add flux from above to both vapor and cloud water of current level
           ! this means that flux entering clear portion of cell from above evaporates
           ! instantly

           dum=lcldm(k,i)/lcldm(k-1,i)
           dum=min(dum,1._r8)
           dum1=icldm(k,i)/icldm(k-1,i)
           dum1=min(dum1,1._r8)

           faltndqie=(falouti(k)-falouti(k-1))/pdel(k,i)
           faltndi=(falouti(k)-dum1*falouti(k-1))/pdel(k,i)
           faltndni=(faloutni(k)-dum1*faloutni(k-1))/pdel(k,i)
           faltndqce=(faloutc(k)-faloutc(k-1))/pdel(k,i)
           faltndc=(faloutc(k)-dum*faloutc(k-1))/pdel(k,i)
           faltndnc=(faloutnc(k)-dum*faloutnc(k-1))/pdel(k,i)

           ! add fallout terms to eulerian tendencies

           qitend(k,i) = qitend(k,i)-faltndi/nstep
           nitend(k,i) = nitend(k,i)-faltndni/nstep
           qctend(k,i) = qctend(k,i)-faltndc/nstep
           nctend(k,i) = nctend(k,i)-faltndnc/nstep

           ! sedimentation tendencies for output
           qcsedten(k,i)=qcsedten(k,i)-faltndc/nstep
           qisedten(k,i)=qisedten(k,i)-faltndi/nstep

           ! add terms to to evap/sub of cloud water

           qvlat(k,i)=qvlat(k,i)-(faltndqie-faltndi)/nstep
           ! for output
           qisevap(k,i)=qisevap(k,i)-(faltndqie-faltndi)/nstep
           qvlat(k,i)=qvlat(k,i)-(faltndqce-faltndc)/nstep
           ! for output
           qcsevap(k,i)=qcsevap(k,i)-(faltndqce-faltndc)/nstep

           tlat(k,i)=tlat(k,i)+(faltndqie-faltndi)*xxls/nstep
           tlat(k,i)=tlat(k,i)+(faltndqce-faltndc)*xxlv/nstep

           dumi(k,i) = dumi(k,i)-faltndi*deltat/nstep
           dumni(k,i) = dumni(k,i)-faltndni*deltat/nstep
           dumc(k,i) = dumc(k,i)-faltndc*deltat/nstep
           dumnc(k,i) = dumnc(k,i)-faltndnc*deltat/nstep

        end do   !! k loop

        ! units below are m/s
        ! cloud water/ice sedimentation flux at surface 
        ! is added to precip flux at surface to get total precip (cloud + precip water)
        ! rate

        prect(i) = prect(i)+(faloutc(nlev)+falouti(nlev))/g/nstep/1000._r8
        preci(i) = preci(i)+(              falouti(nlev))/g/nstep/1000._r8

     end do   !! nstep loop

     ! end sedimentation
     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     ! get new update for variables that includes sedimentation tendency
     ! note : here dum variables are grid-average, NOT in-cloud

     do k=1,nlev

        dumc(k,i) = max(qc(k,i)+qctend(k,i)*deltat,0._r8)
        dumi(k,i) = max(qi(k,i)+qitend(k,i)*deltat,0._r8)
        dumnc(k,i) = max(nc(k,i)+nctend(k,i)*deltat,0._r8)
        dumni(k,i) = max(ni(k,i)+nitend(k,i)*deltat,0._r8)

        ! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(k,i)=ncnst/rho(k,i)*lcldm(k,i)
        end if

        ! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
           dumni(k,i)=ninst/rho(k,i)*icldm(k,i)
        end if

        if (dumc(k,i).lt.qsmall) dumnc(k,i)=0._r8
        if (dumi(k,i).lt.qsmall) dumni(k,i)=0._r8

        ! calculate instantaneous processes (melting, homogeneous freezing)
        !====================================================================

        if (do_cldice) then
           if (t(k,i)+tlat(k,i)/cpp*deltat > tmelt) then
              if (dumi(k,i) > 0._r8) then

                 ! limit so that melting does not push temperature below freezing
                 !-----------------------------------------------------------------
                 dum = -dumi(k,i)*xlf/cpp
                 if (t(k,i)+tlat(k,i)/cpp*deltat+dum.lt.tmelt) then
                    dum = (t(k,i)+tlat(k,i)/cpp*deltat-tmelt)*cpp/xlf
                    dum = dum/dumi(k,i)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qctend(k,i)=qctend(k,i)+dum*dumi(k,i)/deltat

                 ! for output
                 melttot(k,i)=dum*dumi(k,i)/deltat

                 ! assume melting ice produces droplet
                 ! mean volume radius of 8 micron

                 nctend(k,i)=nctend(k,i)+3._r8*dum*dumi(k,i)/deltat/ &
                      (4._r8*pi*5.12e-16_r8*rhow)

                 qitend(k,i)=((1._r8-dum)*dumi(k,i)-qi(k,i))/deltat
                 nitend(k,i)=((1._r8-dum)*dumni(k,i)-ni(k,i))/deltat
                 tlat(k,i)=tlat(k,i)-xlf*dum*dumi(k,i)/deltat
              end if
           end if

           ! homogeneously freeze droplets at -40 C
           !-----------------------------------------------------------------

           if (t(k,i)+tlat(k,i)/cpp*deltat < 233.15_r8) then
              if (dumc(k,i) > 0._r8) then

                 ! limit so that freezing does not push temperature above threshold
                 dum = dumc(k,i)*xlf/cpp
                 if (t(k,i)+tlat(k,i)/cpp*deltat+dum.gt.233.15_r8) then
                    dum = -(t(k,i)+tlat(k,i)/cpp*deltat-233.15_r8)*cpp/xlf
                    dum = dum/dumc(k,i)
                    dum = max(0._r8,dum)
                    dum = min(1._r8,dum)
                 else
                    dum = 1._r8
                 end if

                 qitend(k,i)=qitend(k,i)+dum*dumc(k,i)/deltat
                 ! for output
                 homotot(k,i)=dum*dumc(k,i)/deltat

                 ! assume 25 micron mean volume radius of homogeneously frozen droplets
                 ! consistent with size of detrained ice in stratiform.F90
                 nitend(k,i)=nitend(k,i)+dum*3._r8*dumc(k,i)/(4._r8*3.14_r8*1.563e-14_r8* &
                      500._r8)/deltat
                 qctend(k,i)=((1._r8-dum)*dumc(k,i)-qc(k,i))/deltat
                 nctend(k,i)=((1._r8-dum)*dumnc(k,i)-nc(k,i))/deltat
                 tlat(k,i)=tlat(k,i)+xlf*dum*dumc(k,i)/deltat
              end if
           end if

           ! remove any excess over-saturation, which is possible due to non-linearity when adding 
           ! together all microphysical processes
           !-----------------------------------------------------------------
           ! follow code similar to old CAM scheme

           qtmp=q(k,i)+qvlat(k,i)*deltat
           ttmp=t(k,i)+tlat(k,i)/cpp*deltat

           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p(k,i), esn, qsn)
           qsn = min(qsn,1._r8)

           if (qtmp > qsn .and. qsn > 0) then
              ! expression below is approximate since there may be ice deposition
              dum = (qtmp-qsn)/(1._r8+cons27*qsn/(cpp*rv*ttmp**2))/deltat
              ! add to output cme
              cmeout(k,i) = cmeout(k,i)+dum
              ! now add to tendencies, partition between liquid and ice based on temperature
              if (ttmp > 268.15_r8) then
                 dum1=0.0_r8
                 ! now add to tendencies, partition between liquid and ice based on te
                 !-------------------------------------------------------
              else if (ttmp < 238.15_r8) then
                 dum1=1.0_r8
              else
                 dum1=(268.15_r8-ttmp)/30._r8
              end if

              dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                   *qsn/(cpp*rv*ttmp**2))/deltat
              qctend(k,i)=qctend(k,i)+dum*(1._r8-dum1)
              ! for output
              qcrestot(k,i)=dum*(1._r8-dum1)
              qitend(k,i)=qitend(k,i)+dum*dum1
              qirestot(k,i)=dum*dum1
              qvlat(k,i)=qvlat(k,i)-dum
              ! for output
              qvres(k,i)=-dum
              tlat(k,i)=tlat(k,i)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
           end if
        end if

        ! calculate effective radius for pass to radiation code
        !=========================================================
        ! if no cloud water, default value is 10 micron for droplets,
        ! 25 micron for cloud ice

        ! update cloud variables after instantaneous processes to get effective radius
        ! variables are in-cloud to calculate size dist parameters

        dumc(k,i) = max(qc(k,i)+qctend(k,i)*deltat,0._r8)/lcldm(k,i)
        dumi(k,i) = max(qi(k,i)+qitend(k,i)*deltat,0._r8)/icldm(k,i)
        dumnc(k,i) = max(nc(k,i)+nctend(k,i)*deltat,0._r8)/lcldm(k,i)
        dumni(k,i) = max(ni(k,i)+nitend(k,i)*deltat,0._r8)/icldm(k,i)

        ! hm add 6/2/11 switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(k,i)=ncnst/rho(k,i)
        end if

        ! hm add 6/2/11 switch for specification of cloud ice number
        if (nicons) then
           dumni(k,i)=ninst/rho(k,i)
        end if

        ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1

        dumc(k,i)=min(dumc(k,i),5.e-3_r8)
        dumi(k,i)=min(dumi(k,i),5.e-3_r8)

        ! cloud ice effective radius
        !-----------------------------------------------------------------

        if (do_cldice) then
           if (dumi(k,i).ge.qsmall) then

              dum_2D(k,i) = dumni(k,i)
              call size_dist_param_ice(dumi(k,i), dumni(k,i), lami(k,i), n0i(k,i))

              if (dumni(k,i) /=dum_2D(k,i)) then
                 ! adjust number conc if needed to keep mean size in reasonable range
                 nitend(k,i)=(dumni(k,i)*icldm(k,i)-ni(k,i))/deltat
              end if

              effi(k,i) = 1.5_r8/lami(k,i)*1.e6_r8

           else
              effi(k,i) = 25._r8
           end if

           ! ice effective diameter for david mitchell's optics
           deffi(k,i)=effi(k,i)*rhoi/917._r8*2._r8
        else
           ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
           ! radius has already been determined from the size distribution.
           effi(k,i) = re_ice(k,i) * 1.e6_r8      ! m -> um
           deffi(k,i)=effi(k,i) * 2._r8
        end if

        ! cloud droplet effective radius
        !-----------------------------------------------------------------
        if (dumc(k,i).ge.qsmall) then


           ! hm add 6/2/11 switch for specification of droplet and crystal number
           if (nccons) then
              ! make sure nc is consistence with the constant N by adjusting tendency, need
              ! to multiply by cloud fraction
              ! note that nctend may be further adjusted below if mean droplet size is
              ! out of bounds

              nctend(k,i)=(ncnst/rho(k,i)*lcldm(k,i)-nc(k,i))/deltat

           end if           

           dum = dumnc(k,i)

           call size_dist_param_liq(dumc(k,i), dumnc(k,i), cdnl, rho(k,i), .true., &
                                    pgam(k,i), lamc(k,i))

           if (dum /= dumnc(k,i)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nctend(k,i)=(dumnc(k,i)*lcldm(k,i)-nc(k,i))/deltat
           end if

           effc(k,i) = gamma(pgam(k,i)+4._r8)/ &
                gamma(pgam(k,i)+3._r8)/lamc(k,i)/2._r8*1.e6_r8
           !assign output fields for shape here
           lamcrad(k,i)=lamc(k,i)
           pgamrad(k,i)=pgam(k,i)


           ! recalculate effective radius for constant number, in order to separate
           ! first and second indirect effects
           !======================================
           ! assume constant number of 10^8 kg-1

           dumnc(k,i)=1.e8_r8

           ! Pass in "false" adjust flag to prevent number from being changed within
           ! size distribution subroutine.
           call size_dist_param_liq(dumc(k,i), dumnc(k,i), cdnl, rho(k,i), .false., &
                                    pgam(k,i), lamc(k,i))

           effc_fn(k,i) = gamma(pgam(k,i)+4._r8)/ &
                gamma(pgam(k,i)+3._r8)/lamc(k,i)/2._r8*1.e6_r8

        else
           effc(k,i) = 10._r8
           lamcrad(k,i)=0._r8
           pgamrad(k,i)=0._r8
           effc_fn(k,i) = 10._r8
        end if

     end do ! vertical k loop

     do k=1,nlev
        ! if updated q (after microphysics) is zero, then ensure updated n is also zero
        !=================================================================================
        if (qc(k,i)+qctend(k,i)*deltat.lt.qsmall) nctend(k,i)=-nc(k,i)/deltat
        if (do_cldice .and. qi(k,i)+qitend(k,i)*deltat.lt.qsmall) nitend(k,i)=-ni(k,i)/deltat
     end do

  end do sed_col_loop! i loop

  ! DO STUFF FOR OUTPUT:
  !==================================================

  ! qc and qi are only used for output calculations past here,
  ! so add qctend and qitend back in one more time
  qc = qc + qctend*deltat
  qi = qi + qitend*deltat

  ! averaging for snow and rain number and diameter
  !--------------------------------------------------

  ! drout2/dsout2:
  ! diameter of rain and snow
  ! dsout:
  ! scaled diameter of snow (passed to radiation in CAM)
  ! reff_rain/reff_snow:
  ! calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual

  where (qrout .gt. 1.e-7_r8 &
       .and. nrout.gt.0._r8)
     qrout2 = qrout * cldmax
     nrout2 = nrout * cldmax
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just drout2 times constants.
     drout2 = avg_diameter(qrout, nrout, rho, rhow) * cldmax
     freqr = cldmax

     reff_rain=1.5_r8*drout2*1.e6_r8 
  elsewhere
     qrout2 = 0._r8
     nrout2 = 0._r8
     drout2 = 0._r8
     freqr = 0._r8
     reff_rain = 0._r8
  end where

  where (qsout .gt. 1.e-7_r8 &
       .and. nsout.gt.0._r8)
     qsout2 = qsout * cldmax
     nsout2 = nsout * cldmax
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just dsout2 times constants.
     dsout2 = avg_diameter(qsout, nsout, rho, rhosn)
     freqs = cldmax

     dsout=3._r8*rhosn/917._r8*dsout2

     !----------------LiXH add for RRTMG_4DDA--------------------
     rsout=dsout*917._r8/(2._r8*rhosn)
     !----------------LiXH add for RRTMG_4DDA--------------------

     dsout2 = dsout2 * cldmax

     reff_snow=1.5_r8*dsout2*1.e6_r8
  elsewhere
     dsout  = 0._r8
     qsout2 = 0._r8
     nsout2 = 0._r8
     dsout2 = 0._r8
     freqs  = 0._r8
     reff_snow=0._r8
     !----------------LiXH add for RRTMG_4DDA--------------------
     rsout  = 0._r8 
     !----------------LiXH add for RRTMG_4DDA--------------------
  end where

  ! analytic radar reflectivity
  !--------------------------------------------------
  ! formulas from Matthew Shupe, NOAA/CERES
  ! *****note: radar reflectivity is local (in-precip average)
  ! units of mm^6/m^3

  do i = 1,mgncol
     do k=1,nlev
        if (qc(k,i).ge.qsmall) then
           dum=(qc(k,i)/lcldm(k,i)*rho(k,i)*1000._r8)**2 &
                /(0.109_r8*(nc(k,i)+nctend(k,i)*deltat)/lcldm(k,i)*rho(k,i)/1.e6_r8)*lcldm(k,i)/cldmax(k,i)
        else
           dum=0._r8
        end if
        if (qi(k,i).ge.qsmall) then
           dum1=(qi(k,i)*rho(k,i)/icldm(k,i)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(k,i)/cldmax(k,i)
        else 
           dum1=0._r8
        end if

        if (qsout(k,i).ge.qsmall) then
           dum1=dum1+(qsout(k,i)*rho(k,i)*1000._r8/0.1_r8)**(1._r8/0.63_r8)
        end if

        refl(k,i)=dum+dum1

        ! add rain rate, but for 37 GHz formulation instead of 94 GHz
        ! formula approximated from data of Matrasov (2007)
        ! rainrt is the rain rate in mm/hr
        ! reflectivity (dum) is in DBz

        if (rainrt(k,i).ge.0.001_r8) then
           dum=log10(rainrt(k,i)**6._r8)+16._r8

           ! convert from DBz to mm^6/m^3

           dum = 10._r8**(dum/10._r8)
        else
           ! don't include rain rate in R calculation for values less than 0.001 mm/hr
           dum=0._r8
        end if

        ! add to refl

        refl(k,i)=refl(k,i)+dum

        !output reflectivity in Z.
        areflz(k,i)=refl(k,i) * cldmax(k,i)

        ! convert back to DBz 

        if (refl(k,i).gt.minrefl) then 
           refl(k,i)=10._r8*log10(refl(k,i))
        else
           refl(k,i)=-9999._r8
        end if

        !set averaging flag
        if (refl(k,i).gt.mindbz) then 
           arefl(k,i)=refl(k,i) * cldmax(k,i)
           frefl(k,i)=cldmax(k,i)
        else
           arefl(k,i)=0._r8
           areflz(k,i)=0._r8
           frefl(k,i)=0._r8
        end if

        ! bound cloudsat reflectivity

        csrfl(k,i)=min(csmax,refl(k,i))

        !set averaging flag
        if (csrfl(k,i).gt.csmin) then 
           acsrfl(k,i)=refl(k,i) * cldmax(k,i)
           fcsrfl(k,i)=cldmax(k,i)
        else
           acsrfl(k,i)=0._r8
           fcsrfl(k,i)=0._r8
        end if

     end do
  end do

  !redefine fice here....
  dum_2D = qsout + qrout + qc + qi
  dumi = qsout + qi
  where (dumi .gt. qsmall .and. dum_2D .gt. qsmall)
     nfice=min(dumi/dum_2D,1._r8)
  elsewhere
     nfice=0._r8
  end where

  ! Unpack all outputs

  ! Avoid zero/near-zero division.
  qcsinksum_rate1ord = qcsinksum_rate1ord/max(qcsum_rate1ord,1.0e-30_r8)

  call unpack_array(qcsinksum_rate1ord, mgcols, top_lev, 0._r8, rate1ord_cw2pr_st)

  call unpack_array(tlat, mgcols, top_lev, 0._r8, tlato)
  call unpack_array(qvlat, mgcols, top_lev, 0._r8, qvlato)

  call unpack_array(qctend, mgcols, top_lev, 0._r8, qctendo)
  call unpack_array(qitend, mgcols, top_lev, 0._r8, qitendo)

  ! Note that where there is no water, we set nctend and nitend to remove number
  ! concentration as well.
  call unpack_array(nctend, mgcols, top_lev, -ncn/deltat, nctendo)
  if (do_cldice) then
     call unpack_array(nitend, mgcols, top_lev, -nin/deltat, nitendo)
  else
     call unpack_array(nitend, mgcols, top_lev, 0._r8, nitendo)
  end if

  call unpack_array(effc, mgcols, top_lev, 10._r8, effco)
  call unpack_array(effc_fn, mgcols, top_lev, 10._r8, effco_fn)
  call unpack_array(effi, mgcols, top_lev, 25._r8, effio)

  call unpack_array(prect, mgcols, 0._r8, precto)
  call unpack_array(preci, mgcols, 0._r8, precio)

  call unpack_array(nevapr, mgcols, top_lev, 0._r8, nevapro)
  call unpack_array(evapsnow, mgcols, top_lev, 0._r8, evapsnowo)
  call unpack_array(prain, mgcols, top_lev, 0._r8, praino)
  call unpack_array(prodsnow, mgcols, top_lev, 0._r8, prodsnowo)
  call unpack_array(cmeout, mgcols, top_lev, 0._r8, cmeouto)

  call unpack_array(lamcrad, mgcols, top_lev, 0._r8, lamcrado)
  call unpack_array(pgamrad, mgcols, top_lev, 0._r8, pgamrado)
  call unpack_array(deffi, mgcols, top_lev, 0._r8, deffio)

  call unpack_array(qsout, mgcols, top_lev, 0._r8, qsouto)
  call unpack_array(qsout2, mgcols, top_lev, 0._r8, qsouto2)
  call unpack_array(nsout, mgcols, top_lev, 0._r8, nsouto)
  call unpack_array(nsout2, mgcols, top_lev, 0._r8, nsouto2)
  call unpack_array(dsout, mgcols, top_lev, 0._r8, dsouto)
  call unpack_array(dsout2, mgcols, top_lev, 0._r8, dsouto2)
  call unpack_array(rsout, mgcols, top_lev, 0._r8, rsouto)

  call unpack_array(qrout, mgcols, top_lev, 0._r8, qrouto)
  call unpack_array(qrout2, mgcols, top_lev, 0._r8, qrouto2)
  call unpack_array(nrout, mgcols, top_lev, 0._r8, nrouto)
  call unpack_array(nrout2, mgcols, top_lev, 0._r8, nrouto2)
  call unpack_array(drout2, mgcols, top_lev, 0._r8, drouto2)

  call unpack_array(reff_rain, mgcols, top_lev, 0._r8, reff_raino)
  call unpack_array(reff_snow, mgcols, top_lev, 0._r8, reff_snowo)

  call unpack_array(freqs, mgcols, top_lev, 0._r8, freqso)
  call unpack_array(freqr, mgcols, top_lev, 0._r8, freqro)

  call unpack_array(rflx, mgcols, top_lev, 0._r8, rflxo)
  call unpack_array(sflx, mgcols, top_lev, 0._r8, sflxo)

  call unpack_array(qcsevap, mgcols, top_lev, 0._r8, qcsevapo)
  call unpack_array(qisevap, mgcols, top_lev, 0._r8, qisevapo)
  call unpack_array(qvres, mgcols, top_lev, 0._r8, qvreso)
  call unpack_array(cmeitot, mgcols, top_lev, 0._r8, cmeiout)
  call unpack_array(vtrmc, mgcols, top_lev, 0._r8, vtrmco)
  call unpack_array(vtrmi, mgcols, top_lev, 0._r8, vtrmio)
  call unpack_array(qcsedten, mgcols, top_lev, 0._r8, qcsedteno)
  call unpack_array(qisedten, mgcols, top_lev, 0._r8, qisedteno)

  call unpack_array(pratot, mgcols,top_lev, 0._r8, prao)
  call unpack_array(prctot, mgcols,top_lev, 0._r8, prco)
  call unpack_array(mnuccctot, mgcols,top_lev, 0._r8, mnuccco)
  call unpack_array(mnuccttot, mgcols,top_lev, 0._r8, mnuccto)
  call unpack_array(msacwitot, mgcols,top_lev, 0._r8, msacwio)
  call unpack_array(psacwstot, mgcols,top_lev, 0._r8, psacwso)
  call unpack_array(bergstot, mgcols,top_lev, 0._r8, bergso)
  call unpack_array(bergtot, mgcols,top_lev, 0._r8, bergo)
  call unpack_array(melttot, mgcols,top_lev, 0._r8, melto)
  call unpack_array(homotot, mgcols,top_lev, 0._r8, homoo)
  call unpack_array(qcrestot, mgcols,top_lev, 0._r8, qcreso)
  call unpack_array(prcitot, mgcols,top_lev, 0._r8, prcio)
  call unpack_array(praitot, mgcols,top_lev, 0._r8, praio)
  call unpack_array(qirestot, mgcols,top_lev, 0._r8, qireso)
  call unpack_array(mnuccrtot, mgcols,top_lev, 0._r8, mnuccro)
  call unpack_array(pracstot, mgcols,top_lev, 0._r8, pracso)
  call unpack_array(mnuccdtot, mgcols,top_lev, 0._r8, mnuccdo)
  call unpack_array(meltsdttot, mgcols,top_lev, 0._r8, meltsdto)
  call unpack_array(frzrdttot, mgcols,top_lev, 0._r8, frzrdto)

  call unpack_array(refl, mgcols, top_lev, -9999._r8, reflo)
  call unpack_array(arefl, mgcols, top_lev, 0._r8, areflo)
  call unpack_array(areflz, mgcols, top_lev, 0._r8, areflzo)
  call unpack_array(frefl, mgcols, top_lev, 0._r8, freflo)
  call unpack_array(csrfl, mgcols, top_lev, -9999._r8, csrflo)
  call unpack_array(acsrfl, mgcols, top_lev, 0._r8, acsrflo)
  call unpack_array(fcsrfl, mgcols, top_lev, 0._r8, fcsrflo)

  call unpack_array(rercld, mgcols, top_lev, 0._r8, rercldo)

  call unpack_array(nfice, mgcols, top_lev, 0._r8, nficeo)

  call unpack_array(ncai, mgcols, top_lev, 0._r8, ncaio)
  call unpack_array(ncal, mgcols, top_lev, 0._r8, ncalo)

end subroutine micro_mg1_5_tend


!========================================================================
!FORMULAS
!========================================================================

! Calculate correction due to latent heat for evaporation/sublimation
elemental function calc_ab(t, qv, xxl) result(ab)
  real(r8), intent(in) :: t     ! Temperature
  real(r8), intent(in) :: qv    ! Saturation vapor pressure
  real(r8), intent(in) :: xxl   ! Latent heat

  real(r8) :: ab

  real(r8) :: dqsdt

  dqsdt = xxl*qv / (rv * t**2)
  ab = 1._r8 + dqsdt*xxl/cpp

end function calc_ab

! get cloud droplet size distribution parameters
elemental subroutine size_dist_param_liq(qcic, ncic, cdnl, rho, nadjflag,     &
                                         pgam, lamc)

  real(r8), intent(in) :: qcic
  real(r8), intent(inout) :: ncic
  real(r8), intent(in) :: cdnl
  real(r8), intent(in) :: rho
  logical,  intent(in) :: nadjflag ! Whether to adjust number concentration to fall
                                   ! within certain bounds

  real(r8), intent(out) :: pgam
  real(r8), intent(out) :: lamc

  real(r8) :: dumgam1
  real(r8) :: dumgam2
  real(r8) :: lammin
  real(r8) :: lammax

  if (qcic > qsmall) then

     if (nadjflag) then
        ! add upper limit to in-cloud number concentration to prevent numerical error
        ncic=min(ncic,qcic*1.e20_r8)
        ! add lower limit to in-cloud number concentration
        ncic=max(ncic,cdnl/rho) ! sghan minimum in #/cm
     end if

     ! get pgam from fit to observations of martin et al. 1994

     pgam=0.0005714_r8*(ncic/1.e6_r8*rho)+0.2714_r8
     pgam=1._r8/(pgam**2)-1._r8
     pgam=max(pgam,2._r8)
     pgam=min(pgam,15._r8)

     ! calculate lamc
     dumgam1 = gamma(pgam+1._r8)
     dumgam2 = gamma(pgam+4._r8)

     lamc = (pi/6._r8*rhow*ncic*dumgam2/ &
          (qcic*dumgam1))**(1._r8/3._r8)

     ! lammin, 50 micron diameter max mean size

     lammin = (pgam+1._r8)/50.e-6_r8
     lammax = (pgam+1._r8)/2.e-6_r8

     if (lamc < lammin) then
        lamc = lammin

        if (nadjflag) then
           ncic = 6._r8 * lamc**3 * qcic * dumgam1/ &
                (pi * rhow * dumgam2)
        end if
     else if (lamc > lammax) then
        lamc = lammax

        if (nadjflag) then
           ncic = 6._r8 * lamc**3 * qcic * dumgam1/ &
                (pi * rhow * dumgam2)
        end if
     end if

  else
     ! pgam not calculated in this case, so set it to a value likely to cause an error
     ! if it is accidentally used
     ! (gamma function undefined for negative integers)
     pgam = -100._r8
     lamc = 0._r8
  end if

end subroutine size_dist_param_liq

! get ice size distribution parameters
elemental subroutine size_dist_param_ice(qiic, niic, lami, n0i)
  real(r8), intent(in) :: qiic
  real(r8), intent(inout) :: niic

  real(r8), intent(out) :: lami
  real(r8), intent(out) :: n0i

  ! particle mass-diameter relationship
  ! currently we assume spherical particles for cloud ice/snow
  ! m = cD^d
  ! cloud ice mass-diameter relationship
  real(r8), parameter :: ci = rhoi*pi/6._r8

  ! local parameters
  real(r8), parameter :: lammaxi = 1._r8/10.e-6_r8
  real(r8), parameter :: lammini = 1._r8/(2._r8*dcs)

  if (qiic > qsmall) then

     ! add upper limit to in-cloud number concentration to prevent numerical error
     niic = min(niic, qiic * 1.e20_r8)

     lami = (cons1*ci*niic/qiic)**(1._r8/dsph)
     n0i = niic * lami
  
     ! check for slope
     ! adjust vars
     if (lami < lammini) then
        lami = lammini
        n0i = lami**(dsph+1._r8) * qiic/(ci*cons1)
        niic = n0i/lami
     else if (lami > lammaxi) then
        lami = lammaxi
        n0i = lami**(dsph+1._r8) * qiic/(ci*cons1)
        niic = n0i/lami
     end if
  else
     lami = 0._r8
     n0i  = 0._r8
  end if

end subroutine size_dist_param_ice

! get rain size distribution parameters
elemental subroutine size_dist_param_rain(qric, nric, lamr, n0r)
  real(r8), intent(in) :: qric
  real(r8), intent(inout) :: nric

  real(r8), intent(out) :: lamr
  real(r8), intent(out) :: n0r

  ! particle mass-diameter relationship
  ! currently we assume spherical particles for cloud ice/snow
  ! m = cD^d
  ! rain mass-diameter relationship
  ! Rain is hard-coded for spherical drops
  ! Therefore cr is rhow*pi, rather than
  ! using rhow*pi/6 and later multiplying
  ! by gamma(1+dsph) == 6.
  real(r8), parameter :: cr = rhow*pi

  ! local parameters
  real(r8), parameter :: lammaxr = 1._r8/20.e-6_r8
  real(r8), parameter :: lamminr = 1._r8/500.e-6_r8

  if (qric > qsmall) then

     lamr = (cr*nric/qric)**(1._r8/3._r8)
     n0r = nric * lamr
  
     ! check for slope
     ! adjust vars

     if (lamr < lamminr) then
        lamr = lamminr
        n0r = lamr**4 * qric/cr
        nric = n0r/lamr
     else if (lamr > lammaxr) then
        lamr = lammaxr
        n0r = lamr**4 * qric/cr
        nric = n0r/lamr
     end if

  else
     lamr = 0._r8
     n0r  = 0._r8
  end if

end subroutine size_dist_param_rain

! get snow size distribution parameters
elemental subroutine size_dist_param_snow(qsic, nsic, lams, n0s)
  real(r8), intent(in) :: qsic
  real(r8), intent(inout) :: nsic

  real(r8), intent(out) :: lams
  real(r8), intent(out) :: n0s

  ! particle mass-diameter relationship
  ! currently we assume spherical particles for cloud ice/snow
  ! m = cD^d
  ! cloud ice mass-diameter relationship
  real(r8), parameter :: cs = rhosn*pi/6._r8

  ! local parameters
  real(r8), parameter :: lammaxs = 1._r8/10.e-6_r8
  real(r8), parameter :: lammins = 1._r8/2000.e-6_r8

  if (qsic > qsmall) then

     lams = (cons1*cs*nsic/qsic)**(1._r8/dsph)
     n0s = nsic * lams
  
     ! check for slope
     ! adjust vars
     if (lams < lammins) then
        lams = lammins
        n0s = lams**(dsph+1._r8) * qsic/(cs*cons1)
        nsic = n0s/lams
     else if (lams > lammaxs) then
        lams = lammaxs
        n0s = lams**(dsph+1._r8) * qsic/(cs*cons1)
        nsic = n0s/lams
     end if
  
  else
     lams = 0._r8
     n0s  = 0._r8
  end if

end subroutine size_dist_param_snow

real(r8) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  real(r8), intent(in) :: q         ! mass mixing ratio
  real(r8), intent(in) :: n         ! number concentration
  real(r8), intent(in) :: rho_air   ! local density of the air
  real(r8), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi * rho_sub * n/(q*rho_air))**(-1._r8/3._r8)

end function avg_diameter

real(r8) elemental function var_coef(relvar, a)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in) :: relvar
  real(r8), intent(in) :: a

  var_coef = gamma(relvar + a) / (gamma(relvar) * relvar**a)

end function var_coef

!========================================================================
!MICROPHYSICAL PROCESS CALCULATIONS
!========================================================================

!========================================================================
! Initial ice deposition and sublimation loop.
! Run before the main loop

elemental subroutine ice_deposition_sublimation_init(deltat, t, q, qc, qi, ni,    &
                                                     lcldm, icldm, naai, rho, dv, &
                                                     esl, esi, qvl, qvi, relhum,  &
                                                     berg, cmei)

  ! Inputs
  real(r8), intent(in) :: deltat
  real(r8), intent(in) :: t
  real(r8), intent(in) :: q

  real(r8), intent(in) :: qc
  real(r8), intent(in) :: qi
  real(r8), intent(in) :: ni

  real(r8), intent(in) :: lcldm
  real(r8), intent(in) :: icldm

  real(r8), intent(in) :: naai
  real(r8), intent(in) :: rho
  real(r8), intent(in) :: dv

  real(r8), intent(in) :: esl
  real(r8), intent(in) :: esi
  real(r8), intent(in) :: qvl
  real(r8), intent(in) :: qvi
  real(r8), intent(in) :: relhum

  ! Outputs
  real(r8), intent(out) :: berg
  real(r8), intent(out) :: cmei

  ! Internal variables
  real(r8) :: qiic
  real(r8) :: niic

  real(r8) :: prd               ! provisional deposition rate of cloud ice at water sat 
  real(r8) :: bergtsf           ! bergeron timescale to remove all liquid
  real(r8) :: rhin              ! modified RH for vapor deposition

  real(r8) :: epsi              ! 1/sat relaxation timescale for ice

  real(r8) :: ab
  real(r8) :: lami
  real(r8) :: n0i
  real(r8) :: dum

  ! initialize bergeron process to zero
  berg = 0._r8

  ! Initialize CME components
  cmei = 0._r8
  
  if (t < icenuct) then
     ! provisional nucleation rate
     dum = max((naai - ni/icldm)/deltat*icldm,0._r8)
  else
     dum = 0._r8
  end if

  ! get in-cloud qi and ni after nucleation
  qiic = (qi + dum*deltat*mi0)/icldm
  niic = (ni + dum*deltat)/icldm

  ! hm add 6/2/11 switch for specification of cloud ice number
  if (nicons) niic = ninst/rho

  !ICE DEPOSITION:
  !=============================================
  !if ice exists
  if (t < tmelt .and. qi >= qsmall) then

     ab = calc_ab(t, qvi, xxls)

     ! get ice size distribution parameters
     call size_dist_param_ice(qiic, niic, lami, n0i)

     epsi = 2._r8*pi*n0i*rho*Dv/(lami*lami)
              
     !if liquid exists
     if (qc >= qsmall) then
                 
        ! calculate Bergeron process
        berg = epsi*(qvl-qvi)/ab
                 
        ! multiply by cloud fraction
        berg = berg*min(icldm,lcldm)
                 
        ! Must be positive
        if (berg <= 0._r8) then
           berg = 0._r8
        else
            !BERGERON LIMITING WHEN ALL LIQUID DEPLETED IN 1 TIMESTEP
           !-------------------------------------------------------------
                 
           bergtsf = (qc/berg) / deltat ! bergeron time scale (fraction of timestep)
                    
           if (bergtsf < 1._r8) then
              berg = qc/deltat
                 
              rhin  = (1.0_r8 + relhum) / 2._r8
              !assume RH for frac of step w/ no liq is 1/2 way btwn cldy & cell-ave RH.

              if ((rhin*esl/esi) > 1._r8) then !if ice saturated (but all liquid evap'd)
                 prd = epsi*(rhin*qvl-qvi)/ab

                 ! multiply by cloud fraction assuming liquid/ice maximum overlap
                 prd = prd*min(icldm,lcldm)
   
                 ! add to cmei
                 cmei = cmei + (prd * (1._r8- bergtsf))
              end if ! rhin 
                    
           end if
                    
        end if
     end if

     !Ice deposition in frac of cell with no liquid
     !-------------------------------------------------------------
     ! store liquid cloud fraction in 'dum'
     if (qc >= qsmall) then
        dum = lcldm
     else
        ! for case of no liquid, need to set liquid cloud fraction to zero
        dum = 0._r8
     end if

     if (icldm > dum) then

        ! set RH to grid-mean value for pure ice cloud
        rhin = relhum
                 
        if ((rhin*esl/esi) > 1._r8) then !if rh over ice>ice saturation.
                    
           prd = epsi*(rhin*qvl-qvi)/ab
                    
           ! multiply by relevant cloud fraction for pure ice cloud
           ! assuming maximum overlap of liquid/ice
           prd = prd*(icldm-dum) !apply to ice-only part of cld.
           cmei = cmei + prd
                    
        end if ! rhin
     end if ! qc or icldm > lcldm


     !if grid-mean is ice saturated & qi formed in non-liq cld part, 
     !limit ice formation to avoid mean becoming undersaturated.
     !-------------------------------------------------------------
     if(cmei > 0.0_r8 .and. (relhum*esl/esi) > 1._r8 ) &
          ! max berg is val which removes all ice supersaturation from vapor phase. 
          cmei=min(cmei,(q-qvl*esi/esl)/ab/deltat)

  end if ! end ice exists and t < tmelt

  !ICE SUBLIMATION:
  !=========================================
  !If ice-subsaturated and ice exists:

  if ((relhum*esl/esi) < 1._r8 .and. qiic >= qsmall ) then

     ab = calc_ab(t, qvi, xxls)

     ! get ice size distribution parameters
     call size_dist_param_ice(qiic, niic, lami, n0i)

     epsi = 2._r8*pi*n0i*rho*Dv/(lami*lami)

     ! modify for ice fraction below
     prd = epsi*(relhum*qvl-qvi)/ab * icldm
     cmei=min(prd,0._r8)
           
  endif !subsaturated and ice exists

  ! sublimation should not exceed available ice
  cmei = max(cmei, -qi/deltat)
        
  ! sublimation should not increase grid mean rhi above 1.0 
  if(cmei < 0.0_r8 .and. (relhum*esl/esi) < 1._r8 ) &
       cmei=min(0._r8,max(cmei,(q-qvl*esi/esl)/ab/deltat))

  ! limit cmei due for roundoff error
  cmei = cmei*omsm

end subroutine ice_deposition_sublimation_init

!========================================================================
! autoconversion of cloud liquid water to rain
! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
! minimum qc of 1 x 10^-8 prevents floating point error
elemental subroutine kk2000_liq_autoconversion(qcic, ncic, rho, relvar, &
                                               prc, nprc, nprc1)
  
  real(r8), intent(in) :: qcic
  real(r8), intent(in) :: ncic
  real(r8), intent(in) :: rho
  real(r8), intent(in) :: relvar

  real(r8), intent(out) :: prc
  real(r8), intent(out) :: nprc
  real(r8), intent(out) :: nprc1

  real(r8) :: prc_coef
  real(r8) :: nprc_denom

  ! subcolumn modifications change coefficients
  if (microp_uniform) then
     prc_coef = 1._r8
     nprc_denom = (4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
  else
     prc_coef = var_coef(relvar, 2.47_r8)
     nprc_denom = cons22
  end if

  if (qcic .ge. icsmall) then

     ! nprc is increase in rain number conc due to autoconversion
     ! nprc1 is decrease in cloud droplet conc due to autoconversion
     
     ! assume exponential sub-grid distribution of qc, resulting in additional
     ! factor related to qcvar below
     ! hm switch for sub-columns, don't include sub-grid qc
     
     prc = prc_coef * &
          1350._r8 * qcic**2.47_r8 * (ncic/1.e6_r8*rho)**(-1.79_r8)
     nprc = prc/nprc_denom
     nprc1 = prc/(qcic/ncic)

  else
     prc=0._r8
     nprc=0._r8
     nprc1=0._r8
  end if

end subroutine kk2000_liq_autoconversion

!========================================================================
! Autoconversion of cloud ice to snow
! similar to Ferrier (1994)

elemental subroutine ice_autoconversion(t, qiic, lami, n0i,   &
                                        prci, nprci)

  real(r8), intent(in) :: t
  real(r8), intent(in) :: qiic
  real(r8), intent(in) :: lami
  real(r8), intent(in) :: n0i
  
  real(r8), intent(out) :: prci
  real(r8), intent(out) :: nprci

  if (t .le. tmelt .and.qiic.ge.qsmall) then

     ! note: assumes autoconversion timescale of 180 sec

     nprci = n0i/(lami*180._r8)*exp(-lami*dcs)

     prci = pi*rhoi*n0i/(6._r8*180._r8)* &
          (cons23/lami+3._r8*cons24/lami**2+ &
          6._r8*dcs/lami**3+6._r8/lami**4)*exp(-lami*dcs)          

  else
     prci=0._r8
     nprci=0._r8
  end if

end subroutine ice_autoconversion

! immersion freezing (Bigg, 1953)
!===================================
elemental subroutine immersion_freezing(t, pgam, lamc, cdist1, qcic, relvar,  &
                                        mnuccc, nnuccc)

  ! Temperature
  real(r8), intent(in) :: t

  ! Cloud droplet size distribution parameters
  real(r8), intent(in) :: pgam
  real(r8), intent(in) :: lamc
  real(r8), intent(in) :: cdist1

  ! MMR of in-cloud liquid water
  real(r8), intent(in) :: qcic
  
  ! Relative variance of cloud water
  real(r8), intent(in) :: relvar

  ! Output tendencies
  real(r8), intent(out) :: mnuccc ! MMR
  real(r8), intent(out) :: nnuccc ! Number
  
  ! Coefficients that will be omitted for sub-columns
  real(r8) :: dum, dum1


  if (microp_uniform) then
     dum = 1._r8
     dum1 = 1._r8
  else
     dum = var_coef(relvar, 2._r8)
     dum1 = var_coef(relvar, 1._r8)
  end if

  if (qcic >= qsmall .and. t < 269.15_r8) then

     mnuccc = dum * &
          pi*pi/36._r8*rhow* &
          cdist1*gamma(7._r8+pgam)* &
          bimm*(exp(aimm*(tmelt - t))-1._r8)/lamc**3/lamc**3

     nnuccc = dum1 * &
          pi/6._r8*cdist1*gamma(pgam+4._r8) &
          *bimm*(exp(aimm*(tmelt - t))-1._r8)/lamc**3

  else
     mnuccc = 0._r8
     nnuccc = 0._r8
  end if ! qcic > qsmall and t < 4 deg C

end subroutine immersion_freezing

! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
!===================================================================
! dust size and number in multiple bins are read in from companion routine
pure subroutine contact_freezing (t, p, rndst, nacon, pgam, lamc, cdist1, qcic, relvar,  &
                                  mnucct, nnucct)
  
  real(r8), intent(in) :: t(:)            ! Temperature
  real(r8), intent(in) :: p(:)            ! Pressure
  real(r8), intent(in) :: rndst(:,:)      ! Radius (for multiple dust bins)
  real(r8), intent(in) :: nacon(:,:)      ! Number (for multiple dust bins)

  ! Size distribution parameters for cloud droplets
  real(r8), intent(in) :: pgam(:)
  real(r8), intent(in) :: lamc(:)
  real(r8), intent(in) :: cdist1(:)

  ! MMR of in-cloud liquid water
  real(r8), intent(in) :: qcic(:)
  
  ! Relative cloud water variance
  real(r8), intent(in) :: relvar(:)

  ! Output tendencies
  real(r8), intent(out) :: mnucct(:) ! MMR
  real(r8), intent(out) :: nnucct(:) ! Number

  real(r8) :: tcnt                  ! scaled relative temperature
  real(r8) :: viscosity             ! temperature-specific viscosity (kg/m/s)
  real(r8) :: mfp                   ! temperature-specific mean free path (m)

  ! Dimension these according to number of dust bins, inferred from rndst size
  real(r8) :: nslip(size(rndst,1))  ! slip correction factors
  real(r8) :: ndfaer(size(rndst,1)) ! aerosol diffusivities (m^2/sec)

  ! Coefficients not used for subcolumns
  real(r8) :: dum, dum1

  integer  :: i

  ! subcolumns

  do i = 1,size(t)

     if (qcic(i) >= qsmall .and. t(i) < 269.15_r8) then

        if (microp_uniform) then
           dum = 1._r8
           dum1 = 1._r8
        else
           dum = var_coef(relvar(i), 4._r8/3._r8)
           dum1 = var_coef(relvar(i), 1._r8/3._r8)
        endif

        tcnt=(270.16_r8-t(i))**1.3_r8
        viscosity = 1.8e-5_r8*(t(i)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
        mfp = 2.0_r8*viscosity/ &                         ! Mean free path (m)
                     (p(i)*sqrt( 8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(i)) ))

        ! Note that these two are vectors.
        nslip = 1.0_r8+(mfp/rndst(:,i))*(1.257_r8+(0.4_r8*exp(-(1.1_r8*rndst(:,i)/mfp))))! Slip correction factor
        ndfaer = 1.381e-23_r8*t(i)*nslip/(6._r8*pi*viscosity*rndst(:,i))  ! aerosol diffusivity (m2/s)

        mnucct(i) = dum *  &
             dot_product(ndfaer,nacon(:,i)*tcnt)*pi*pi/3._r8*rhow* &
             cdist1(i)*gamma(pgam(i)+5._r8)/lamc(i)**4

        nnucct(i) =  dum1 *  &
             dot_product(ndfaer,nacon(:,i)*tcnt)*2._r8*pi*  &
             cdist1(i)*gamma(pgam(i)+2._r8)/lamc(i)

     else

        mnucct(i)=0._r8
        nnucct(i)=0._r8

     end if ! qcic > qsmall and t < 4 deg C
  end do

end subroutine contact_freezing

! snow self-aggregation from passarelli, 1978, used by reisner, 1998
!===================================================================
! this is hard-wired for bs = 0.4 for now
! ignore self-collection of cloud ice

elemental subroutine snow_self_aggregation(t, rho, asn, qsic, nsic, nsagg)

  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: rho  ! Density
  real(r8), intent(in) :: asn  ! fall speed parameter for snow

  ! In-cloud snow
  real(r8), intent(in) :: qsic ! MMR
  real(r8), intent(in) :: nsic ! Number

  ! Output number tendency
  real(r8), intent(out) :: nsagg

  if (qsic >= qsmall .and. t <= tmelt) then
     nsagg = -1108._r8*asn*eii* &
          pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)* &
          rho**((2._r8+bs)/3._r8)*qsic**((2._r8+bs)/3._r8)* &
          (nsic*rho)**((4._r8-bs)/3._r8) /(4._r8*720._r8*rho)
  else
     nsagg=0._r8
  end if

end subroutine snow_self_aggregation

! accretion of cloud droplets onto snow/graupel
!===================================================================
! here use continuous collection equation with
! simple gravitational collection kernel
! ignore collisions between droplets/cloud ice
! since minimum size ice particle for accretion is 50 - 150 micron
elemental subroutine accrete_cloud_water_snow(t, rho, asn, uns, mu,             &
                                              qcic, ncic, qsic, pgam,           &
                                              lamc, lams, n0s,                  & 
                                              psacws, npsacws)

  real(r8), intent(in) :: t   ! Temperature
  real(r8), intent(in) :: rho ! Density
  real(r8), intent(in) :: asn ! Fallspeed parameter (snow)
  real(r8), intent(in) :: uns ! Current fallspeed   (snow)
  real(r8), intent(in) :: mu  ! Viscosity

  ! In-cloud liquid water
  real(r8), intent(in) :: qcic ! MMR
  real(r8), intent(in) :: ncic ! Number

  ! In-cloud snow
  real(r8), intent(in) :: qsic ! MMR

  ! Cloud droplet size parameters
  real(r8), intent(in) :: pgam
  real(r8), intent(in) :: lamc

  ! Snow size parameters
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: psacws  ! Mass mixing ratio
  real(r8), intent(out) :: npsacws ! Number concentration

  real(r8) :: dc0 ! Provisional mean droplet size
  real(r8) :: dum
  real(r8) :: eci ! collection efficiency for riming of snow by droplets

  ! ignore collision of snow with droplets above freezing

  if (qsic >= qsmall .and. t <= tmelt .and. qcic >= qsmall) then

     ! put in size dependent collection efficiency
     ! mean diameter of snow is area-weighted, since
     ! accretion is function of crystal geometric area
     ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

     dc0 = (pgam+1._r8)/lamc
     dum = dc0*dc0*uns*rhow/(9._r8*mu*(1._r8/lams))
     eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

     eci = max(eci,0._r8)
     eci = min(eci,1._r8)

     ! no impact of sub-grid distribution of qc since psacws
     ! is linear in qc

     psacws = pi/4._r8*asn*qcic*rho*n0s*eci*cons11 / lams**(bs+3._r8)
     npsacws = pi/4._r8*asn*ncic*rho*n0s*eci*cons11 / lams**(bs+3._r8)
  else
     psacws = 0._r8
     npsacws = 0._r8
  end if

end subroutine accrete_cloud_water_snow

! add secondary ice production due to accretion of droplets by snow 
!===================================================================
! (Hallet-Mossop process) (from Cotton et al., 1986)
elemental subroutine secondary_ice_production(t, psacws, msacwi, nsacwi)
  real(r8), intent(in) :: t ! Temperature

  ! Accretion of cloud water to snow tendencies
  real(r8), intent(inout) :: psacws ! MMR

  ! Output (ice) tendencies
  real(r8), intent(out) :: msacwi ! MMR
  real(r8), intent(out) :: nsacwi ! Number

  if((t < 270.16_r8) .and. (t >= 268.16_r8)) then
     nsacwi = 3.5e8_r8*(270.16_r8-t)/2.0_r8*psacws
     msacwi = min(nsacwi*mi0, psacws)
  else if((t < 268.16_r8) .and. (t >= 265.16_r8)) then
     nsacwi = 3.5e8_r8*(t-265.16_r8)/3.0_r8*psacws
     msacwi = min(nsacwi*mi0, psacws)
  else
     nsacwi = 0.0_r8
     msacwi = 0.0_r8
  endif

  psacws = max(0.0_r8,psacws - nsacwi*mi0)

end subroutine secondary_ice_production

! accretion of rain water by snow
!===================================================================
! formula from ikawa and saito, 1991, used by reisner et al., 1998
elemental subroutine accrete_rain_snow(t, rho, umr, ums, unr, uns,        &
                                       qric, qsic, lamr, n0r, lams, n0s,  &
                                       pracs, npracs )

  real(r8), intent(in) :: t   ! Temperature
  real(r8), intent(in) :: rho ! Density

   ! Fallspeeds
  ! mass-weighted
  real(r8), intent(in) :: umr ! rain
  real(r8), intent(in) :: ums ! snow
  ! number-weighted
  real(r8), intent(in) :: unr ! rain
  real(r8), intent(in) :: uns ! snow

  ! In cloud MMRs
  real(r8), intent(in) :: qric ! rain
  real(r8), intent(in) :: qsic ! snow

  ! Size distribution parameters
  ! rain
  real(r8), intent(in) :: lamr
  real(r8), intent(in) :: n0r
  ! snow
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: pracs  ! MMR
  real(r8), intent(out) :: npracs ! Number

  ! Collection efficiency for accretion of rain by snow
  real(r8), parameter :: ecr = 1.0_r8

  if (qric >= icsmall .and. qsic >= icsmall .and. t <= tmelt) then

     pracs = pi*pi*ecr*(((1.2_r8*umr-0.95_r8*ums)**2 + &
          0.08_r8*ums*umr)**0.5_r8 *  &
          rhow * rho * n0r * n0s * &
          (5._r8/(lamr**6 * lams)+ &
          2._r8/(lamr**5 * lams**2)+ &
          0.5_r8/(lamr**4 * lams**3)))

     npracs = pi/2._r8*rho*ecr* (1.7_r8*(unr-uns)**2 + &
          0.3_r8*unr*uns)**0.5_r8 * &
          n0r*n0s* &
          (1._r8/(lamr**3 * lams)+ &
          1._r8/(lamr**2 * lams**2)+ &
          1._r8/(lamr * lams**3))

  else
     pracs = 0._r8
     npracs = 0._r8
  end if

end subroutine accrete_rain_snow

! heterogeneous freezing of rain drops
!===================================================================
! follows from Bigg (1953)
elemental subroutine heterogeneous_rain_freezing(t, qric, nric, lamr,     &
                                                 mnuccr, nnuccr)

  real(r8), intent(in) :: t    ! Temperature

  ! In-cloud rain
  real(r8), intent(in) :: qric ! MMR
  real(r8), intent(in) :: nric ! Number
  real(r8), intent(in) :: lamr ! size parameter

  ! Output tendencies
  real(r8), intent(out) :: mnuccr ! MMR
  real(r8), intent(out) :: nnuccr ! Number

  if (t < 269.15_r8 .and. qric >= qsmall) then

     ! Division by lamr**3 twice is old workaround to avoid overflow.
     ! Probably no longer necessary
     mnuccr = 20._r8*pi*pi*rhow*nric*bimm* &
          (exp(aimm*(tmelt - t))-1._r8)/lamr**3 &
          /lamr**3

     nnuccr = pi*nric*bimm* &
          (exp(aimm*(tmelt - t))-1._r8)/lamr**3
  else
     mnuccr = 0._r8
     nnuccr = 0._r8
  end if
end subroutine heterogeneous_rain_freezing

! accretion of cloud liquid water by rain
!===================================================================
! formula from Khrouditnov and Kogan (2000)
! gravitational collection kernel, droplet fall speed neglected
elemental subroutine accrete_cloud_water_rain(qric, qcic, ncic, relvar, accre_enhan,  &
                                              pra, npra)

  ! In-cloud rain
  real(r8), intent(in) :: qric ! MMR

  ! Cloud droplets
  real(r8), intent(in) :: qcic ! MMR
  real(r8), intent(in) :: ncic ! Number
  
  ! SGS variability
  real(r8), intent(in) :: relvar
  real(r8), intent(in) :: accre_enhan

  ! Output tendencies
  real(r8), intent(out) :: pra  ! MMR
  real(r8), intent(out) :: npra ! Number

  ! Coefficient that varies for subcolumns
  real(r8) :: pra_coef

  if (microp_uniform) then
     pra_coef = 1._r8
  else
     pra_coef = accre_enhan * var_coef(relvar, 1.15_r8)
  end if

  if (qric >= qsmall .and. qcic >= qsmall) then

     ! include sub-grid distribution of cloud water
     pra = pra_coef * 67._r8*(qcic*qric)**1.15_r8

     npra = pra/(qcic/ncic)

  else
     pra = 0._r8
     npra = 0._r8
  end if
end subroutine accrete_cloud_water_rain

! Self-collection of rain drops
!===================================================================
! from Beheng(1994)
elemental subroutine self_collection_rain(rho, qric, nric, nragg)

  real(r8), intent(in) :: rho  ! Air density

  ! Rain
  real(r8), intent(in) :: qric ! MMR
  real(r8), intent(in) :: nric ! Number

  ! Output number tendency
  real(r8), intent(out) :: nragg

  if (qric >= qsmall) then
     nragg = -8._r8*nric*qric*rho
  else
     nragg = 0._r8
  end if

end subroutine self_collection_rain

! Accretion of cloud ice by snow
!===================================================================
! For this calculation, it is assumed that the Vs >> Vi
! and Ds >> Di for continuous collection
elemental subroutine accrete_cloud_ice_snow(t, rho, asn, qiic, niic, qsic, &
                                            lams, n0s, prai, nprai)

  real(r8), intent(in) :: t    ! Temperature
  real(r8), intent(in) :: rho  ! Density

  real(r8), intent(in) :: asn  ! Snow fallspeed parameter

  ! Cloud ice
  real(r8), intent(in) :: qiic ! MMR
  real(r8), intent(in) :: niic ! Number

  real(r8), intent(in) :: qsic ! Snow MMR

  ! Snow size parameters
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s  

  ! Output tendencies
  real(r8), intent(out) :: prai  ! MMR
  real(r8), intent(out) :: nprai ! Number

  if (qsic >= qsmall .and. qiic >= qsmall .and. t <= tmelt) then

     prai = pi/4._r8 * asn * qiic * rho * n0s * eii * cons11/ &
          lams**(bs+3._r8)

     nprai = pi/4._r8 * asn * niic * rho * n0s * eii * cons11/ &
          lams**(bs+3._r8)
  else
     prai = 0._r8
     nprai = 0._r8
  end if

end subroutine accrete_cloud_ice_snow

! calculate evaporation/sublimation of rain and snow
!===================================================================
! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
! in-cloud condensation/deposition of rain and snow is neglected
! except for transfer of cloud water to snow through bergeron process
elemental subroutine evaporate_sublimate_precip(deltat, t, p, rho, dv, mu, sc, q, qvl, qvi, lcldm, cldmax,    &
                                                arn, asn, qcic, qiic, qric, qsic, lamr, n0r, lams, n0s, cmei, &
                                                pre, prds)

  real(r8), intent(in) :: deltat ! timestep

  real(r8), intent(in) :: t    ! temperature
  real(r8), intent(in) :: p    ! pressure
  real(r8), intent(in) :: rho  ! air density
  real(r8), intent(in) :: dv   ! water vapor diffusivity
  real(r8), intent(in) :: mu   ! viscosity
  real(r8), intent(in) :: sc   ! schmidt number
  real(r8), intent(in) :: q    ! humidity
  real(r8), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), intent(in) :: qvi  ! saturation humidity (ice)
  real(r8), intent(in) :: lcldm  ! liquid cloud fraction
  real(r8), intent(in) :: cldmax ! precipitation fraction (maximum overlap)

  ! fallspeed parameters
  real(r8), intent(in) :: arn  ! rain
  real(r8), intent(in) :: asn  ! snow

  ! In-cloud MMRs
  real(r8), intent(in) :: qcic ! cloud liquid
  real(r8), intent(in) :: qiic ! cloud ice
  real(r8), intent(in) :: qric ! rain
  real(r8), intent(in) :: qsic ! snow

  ! Size parameters
  ! rain
  real(r8), intent(in) :: lamr
  real(r8), intent(in) :: n0r
  ! snow
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! cloud ice sublimation/deposition tendency
  real(r8), intent(in) :: cmei

  ! Output tendencies
  real(r8), intent(out) :: pre
  real(r8), intent(out) :: prds

  ! checking for RH after rain evap
  real(r8) :: esn    ! saturation pressure
  real(r8) :: qsn    ! saturation humidity

  real(r8) :: qclr   ! water vapor mixing ratio in clear air
  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  ! Temps/dummies
  real(r8) :: qtmp
  real(r8) :: ttmp

  real(r8) :: dum, dum1

  ! set temporary cloud fraction to zero if cloud water + ice is very small
  ! this will ensure that evaporation/sublimation of precip occurs over
  ! entire grid cell, since min cloud fraction is specified otherwise
  if (qcic+qiic < 1.e-6_r8) then
     dum = 0._r8
  else
     dum = lcldm
  end if

  ! only calculate if there is some precip fraction > cloud fraction

  if (cldmax > dum) then

     ! calculate q for out-of-cloud region
     qsn = min(qvl,1._r8)
     qclr=(q-dum*qsn)/(1._r8-dum)

     ! evaporation of rain
     if (qric.ge.qsmall) then

        ab = calc_ab(t, qvl, xxlv)
        eps = 2._r8*pi*n0r*rho*Dv* &
             (f1r/(lamr*lamr)+ &
             f2r*(arn*rho/mu)**0.5_r8* &
             sc**(1._r8/3._r8)*cons13/ &
             (lamr**(5._r8/2._r8+br/2._r8)))

        pre = eps*(qclr-qvl)/ab

        ! only evaporate in out-of-cloud region
        ! and distribute across cldmax
        pre=min(pre*(cldmax-dum),0._r8)
        pre=pre/cldmax
     else
        pre = 0._r8
     end if

     ! sublimation of snow
     if (qsic.ge.qsmall) then
        ab = calc_ab(t, qvi, xxls)
        eps = 2._r8*pi*n0s*rho*Dv* &
             (f1s/(lams*lams)+ &
             f2s*(asn*rho/mu)**0.5_r8* &
             sc**(1._r8/3._r8)*cons14/ &
             (lams**(5._r8/2._r8+bs/2._r8)))
        prds = eps*(qclr-qvi)/ab

        ! only sublimate in out-of-cloud region and distribute over cldmax
        prds=min(prds*(cldmax-dum),0._r8)
        prds=prds/cldmax
     else
        prds = 0._r8
     end if

     ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
     ! get updated RH at end of time step based on cloud water/ice condensation/evap

     qtmp=q-(cmei+(pre+prds)*cldmax)*deltat
     ttmp=t+((pre*cldmax)*xxlv+ &
          (cmei+prds*cldmax)*xxls)*deltat/cpp

     !limit range of temperatures!
     ttmp=max(180._r8,min(ttmp,323._r8))

     ! use rhw to allow ice supersaturation
     call qsat_water(ttmp, p, esn, qsn)
     qsn=min(qsn,1._r8)

     ! modify precip evaporation rate if q > qsat
     if (qtmp.gt.qsn) then
        if (pre+prds.lt.-1.e-20_r8) then
           dum1=pre/(pre+prds)
           ! recalculate q and t after cloud water cond but without precip evap
           qtmp=q-(cmei)*deltat
           ttmp=t+(cmei*xxls)*deltat/cpp
           ! use rhw to allow ice supersaturation
           call qsat_water(ttmp, p, esn, qsn)
           qsn=min(qsn,1._r8)

           dum=(qtmp-qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))
           dum=min(dum,0._r8)

           ! modify rates if needed, divide by cldmax to get local (in-precip) value
           pre=dum*dum1/deltat/cldmax

           ! do separately using RHI for prds....
           ! use rhi to allow ice supersaturation
           call qsat_ice(ttmp, p, esn, qsn)
           qsn=min(qsn,1._r8)

           dum=(qtmp-qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp**2))
           dum=min(dum,0._r8)

           ! modify rates if needed, divide by cldmax to get local (in-precip) value
           prds=dum*(1._r8-dum1)/deltat/cldmax
        end if
     end if

  else
     prds = 0._r8
     pre = 0._r8
  end if

end subroutine evaporate_sublimate_precip

! bergeron process - evaporation of droplets and deposition onto snow
!===================================================================
elemental subroutine bergeron_process(t, rho, dv, mu, sc, qvl, qvi, asn, &
                                      qcic, qsic, lams, n0s, bergs)

  real(r8), intent(in) :: t    ! temperature
  real(r8), intent(in) :: rho  ! air density
  real(r8), intent(in) :: dv   ! water vapor diffusivity
  real(r8), intent(in) :: mu   ! viscosity
  real(r8), intent(in) :: sc   ! schmidt number
  real(r8), intent(in) :: qvl  ! saturation humidity (water)
  real(r8), intent(in) :: qvi  ! saturation humidity (ice)

  ! fallspeed parameter for snow
  real(r8), intent(in) :: asn

  ! In-cloud MMRs
  real(r8), intent(in) :: qcic ! cloud liquid
  real(r8), intent(in) :: qsic ! snow

  ! Size parameters for snow
  real(r8), intent(in) :: lams
  real(r8), intent(in) :: n0s

  ! Output tendencies
  real(r8), intent(out) :: bergs

  real(r8) :: ab     ! correction to account for latent heat
  real(r8) :: eps    ! 1/ sat relaxation timescale

  if (qsic >= qsmall.and. qcic >= qsmall .and. t < tmelt) then
     ab = calc_ab(t, qvi, xxls)
     eps = 2._r8*pi*n0s*rho*Dv* &
          (f1s/(lams*lams)+ &
          f2s*(asn*rho/mu)**0.5_r8* &
          sc**(1._r8/3._r8)*cons14/ &
          (lams**(5._r8/2._r8+bs/2._r8)))
     bergs = eps*(qvl-qvi)/ab
  else
     bergs = 0._r8
  end if

end subroutine bergeron_process

!========================================================================
!OUTPUT CALCULATIONS
!========================================================================
elemental subroutine calc_rercld(lamr, n0r, lamc, cdist1, pgam, dumr, qcic, &
                                 arcld, rercld)
  real(r8), intent(in) :: lamr          ! rain size parameter (slope)
  real(r8), intent(in) :: n0r           ! rain size parameter (intercept)
  real(r8), intent(in) :: lamc          ! size distribution parameter (slope)
  real(r8), intent(in) :: cdist1        ! for droplet freezing
  real(r8), intent(in) :: pgam          ! droplet size parameter
  real(r8), intent(in) :: dumr          ! in-cloud rain mass mixing ratio
  real(r8), intent(in) :: qcic          ! in-cloud cloud liquid

  integer,  intent(inout) :: arcld      ! number of substeps rercld has been through
  real(r8), intent(inout) :: rercld     ! effective radius calculation for rain + cloud

  ! combined size of precip & cloud drops
  real(r8) :: Atmp

  ! Rain drops
  if (lamr > 0._r8) then
     Atmp = n0r * pi / (2._r8 * lamr**3._r8)
  else
     Atmp = 0._r8
  end if

  ! Add cloud drops
  if (lamc > 0._r8) then
     Atmp = Atmp + cdist1 * pi * gamma(pgam+3._r8)/(4._r8 * lamc**2._r8)
  end if

  if (Atmp > 0._r8) then
     rercld = rercld + 3._r8 *(dumr + qcic) / (4._r8 * rhow * Atmp)
     arcld = arcld+1
  end if

end subroutine calc_rercld


!========================================================================
!UTILITIES
!========================================================================

pure subroutine micro_mg1_5_get_cols(ncol, nlev, top_lev, qcn, qin, &
                                     mgncol, mgcols)

  ! Determines which columns microphysics should operate over by
  ! checking for non-zero cloud water/ice.
  ! io  
  integer, intent(in) :: ncol      ! Number of columns with meaningful data
  integer, intent(in) :: nlev      ! Number of levels to use
  integer, intent(in) :: top_lev   ! Top level for microphysics

  real(r8), intent(in) :: qcn(:,:) ! cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(:,:) ! cloud ice mixing ratio (kg/kg)

  integer, intent(out) :: mgncol   ! Number of columns MG will use
  integer, allocatable, intent(out) :: mgcols(:) ! column indices

  ! local
  integer :: lev_offset  ! top_lev - 1 (defined here for consistency)
  logical :: ltrue(ncol) ! store tests for each column

  integer :: i, ii ! column indices

  if (allocated(mgcols)) deallocate(mgcols)

  lev_offset = top_lev - 1

  ! Using "any" along dimension 2 collapses across levels, but
  ! not columns, so we know if water is present at any level
  ! in each column.

  ltrue = any(qcn(top_lev:(nlev+lev_offset),:ncol) >= qsmall, 1)
  ltrue = ltrue .or. any(qin(top_lev:(nlev+lev_offset),:ncol) >= qsmall, 1)

  ! Scan for true values to get a usable list of indices.
  
  mgncol = count(ltrue)
  allocate(mgcols(mgncol))
  i = 0
  do ii = 1,ncol
     if (ltrue(ii)) then
        i = i + 1
        mgcols(i) = ii
     end if
  end do

end subroutine micro_mg1_5_get_cols

pure function interp_to_mid(orig_val, weights) result(new_val)
  ! Linear interpolation, here used to move from interfaces to midlevel
  real(r8), intent(in) :: orig_val(:,:)
  real(r8), intent(in) :: weights(:,:)

  ! New value is slightly smaller than the old value
  real(r8) :: new_val(size(orig_val,1)-1,size(orig_val,2))

  new_val = orig_val(:size(new_val,1),:)
  new_val = new_val + weights * (orig_val(2:,:) - new_val)

end function interp_to_mid

! Subroutines to pack arrays into smaller, contiguous pieces
!========================================================================
! Rank 1 array of reals, columns only
pure subroutine pack_array_1Dr8(old_array, cols, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:)   ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include

  ! Output
  real(r8), intent(out) :: new_array(:)
  
  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = old_array(cols)
  end if

end subroutine pack_array_1Dr8

! Rank 2 array of reals, columns and levels
pure subroutine pack_array_2Dr8(old_array, cols, top_lev, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  integer,  intent(in)  :: top_lev        ! First level to use

  ! Output
  real(r8), intent(out) :: new_array(:,:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = old_array(top_lev:, cols)
  end if

end subroutine pack_array_2Dr8

! Rank 3 array of reals, assume last index is extra
pure subroutine pack_array_3Dr8(old_array, cols, top_lev, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)          ! List of columns to include
  integer,  intent(in)  :: top_lev          ! First level to use

  ! Output
  real(r8), intent(out) :: new_array(:,:,:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = old_array(:, top_lev:, cols)
  end if

end subroutine pack_array_3Dr8

! Subroutines to unpack arrays for output
!========================================================================
! Rank 1 array of reals, columns only
pure subroutine unpack_array_1Dr8(old_array, cols, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:)   ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  real(r8), intent(in)  :: fill           ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:)
  
  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(cols) = old_array
  end if

end subroutine unpack_array_1Dr8

! Rank 1 array of reals, columns only, "fill" value is an array
pure subroutine unpack_array_1Dr8_arrayfill(old_array, cols, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:)   ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  real(r8), intent(in)  :: fill(:)        ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(cols) = old_array
  end if
  
end subroutine unpack_array_1Dr8_arrayfill

! Rank 2 array of reals, columns and levels
pure subroutine unpack_array_2Dr8(old_array, cols, top_lev, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  integer,  intent(in)  :: top_lev        ! First level to use
  real(r8), intent(in)  :: fill           ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:,:)

  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(top_lev:, cols) = old_array
  end if
  
end subroutine unpack_array_2Dr8

! Rank 2 array of reals, columns and levels, "fill" value is an array
pure subroutine unpack_array_2Dr8_arrayfill(old_array, cols, top_lev, fill, new_array)
  ! Inputs
  real(r8), intent(in)  :: old_array(:,:) ! Array to be packed
  integer,  intent(in)  :: cols(:)        ! List of columns to include
  integer,  intent(in)  :: top_lev        ! First level to use
  real(r8), intent(in)  :: fill(:,:)      ! Value with which to fill unused
                                          ! sections of new_array.

  ! Output
  real(r8), intent(out) :: new_array(:,:)
  
  ! Attempt to speed up packing if it is unnecessary.
  if (size(new_array) == size(old_array)) then
     new_array = old_array
  else
     new_array = fill

     new_array(top_lev:, cols) = old_array
  end if
  
end subroutine unpack_array_2Dr8_arrayfill

end module micro_mg1_5
