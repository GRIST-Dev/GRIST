!===================================================================================================
!
!  Created by LiXiaohan on 19/08/1, adopted from CAM5
!
!     MG microphysics version 1.0
!
!     NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
!     microphysics in cooperation with the MG liquid microphysics. This is
!     controlled by the do_cldice variable.
!    
!     NOTE: If do_cldice is false, then MG microphysics should not update CLDICE
!     or NUMICE; however, it is assumed that the other microphysics scheme will have
!     updated CLDICE and NUMICE. The other microphysics should handle the following
!     processes that would have been done by MG:
!       - Detrainment (liquid and ice)
!       - Homogeneous ice nucleation
!       - Heterogeneous ice nucleation
!       - Bergeron process
!       - Melting of ice
!       - Freezing of cloud drops
!       - Autoconversion (ice -> snow)
!       - Growth/Sublimation of ice
!       - Sedimentation of ice
!
!     Procedures required:
!     1) An implementation of the gamma function (if not intrinsic).
!     2) saturation vapor pressure to specific humidity formula
!     3) svp over water
!     4) svp over ice
!
!     Author: Andrew Gettelman, Hugh Morrison.
!     Contributions from: Xiaohong Liu and Steve Ghan
!     December 2005-May 2010
!     Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                     Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
!
!===================================================================================================

module micro_mg1_0

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8, i4
!    use grist_shr_spfn,                     only: gamma => shr_spfn_gamma
    use grist_mpi

implicit none
private
save

public :: micro_mg1_0_init ,  &
          micro_mg1_0_tend

real(r8) :: g              !gravity
real(r8) :: r              !Dry air Gas constant
real(r8) :: rv             !water vapor gas contstant
real(r8) :: cpp            !specific heat of dry air
real(r8) :: rhow           !density of liquid water
real(r8) :: tmelt          ! Freezing point of water (K)
real(r8) :: xxlv           ! latent heat of vaporization
real(r8) :: xlf            !latent heat of freezing
real(r8) :: xxls           !latent heat of sublimation
real(r8) :: rhmini     ! Minimum rh for ice cloud fraction > 0.

real(r8) :: rhosn  ! bulk density snow
real(r8) :: rhoi   ! bulk density ice

real(r8) :: ac,bc,as,bs,ai,bi,ar,br  ! fall speed parameters 
real(r8) :: ci,di     ! ice mass-diameter relation parameters
real(r8) :: cs,ds     ! snow mass-diameter relation parameters
real(r8) :: cr,dr     ! drop mass-diameter relation parameters
real(r8) :: f1s,f2s   ! ventilation param for snow
real(r8) :: Eii       ! collection efficiency aggregation of ice
real(r8) :: Ecr       ! collection efficiency cloud droplets/rain
real(r8) :: f1r,f2r   ! ventilation param for rain
real(r8) :: DCS       ! autoconversion size threshold
real(r8) :: qsmall    ! min mixing ratio 
real(r8) :: bimm,aimm ! immersion freezing
real(r8) :: rhosu     ! typical 850mn air density
real(r8) :: mi0       ! new crystal mass
real(r8) :: rin       ! radius of contact nuclei
real(r8) :: pi        ! pi

! Additional constants to help speed up code
real(r8) :: cons1
real(r8) :: cons2
real(r8) :: cons3
real(r8) :: cons4
real(r8) :: cons5
real(r8) :: cons6
real(r8) :: cons7
real(r8) :: cons8
real(r8) :: cons9
real(r8) :: cons10
real(r8) :: cons11
real(r8) :: cons12
real(r8) :: cons13
real(r8) :: cons14
real(r8) :: cons15
real(r8) :: cons16
real(r8) :: cons17
real(r8) :: cons18
real(r8) :: cons19
real(r8) :: cons20
real(r8) :: cons22
real(r8) :: cons23
real(r8) :: cons24
real(r8) :: cons25
real(r8) :: cons27
real(r8) :: cons28

real(r8) :: lammini
real(r8) :: lammaxi
real(r8) :: lamminr
real(r8) :: lammaxr
real(r8) :: lammins
real(r8) :: lammaxs


! parameters for snow/rain fraction for convective clouds
real(r8) :: tmax_fsnow ! max temperature for transition to convective snow
real(r8) :: tmin_fsnow ! min temperature for transition to convective snow

!needed for findsp
real(r8) :: tt0       ! Freezing temperature

real(r8) :: csmin,csmax,minrefl,mindbz

contains

! Purpose: initialize constants for the morrison microphysics
subroutine micro_mg1_0_init( gravit_in, rair_in, rh2o_in, cpair_in, rhoh2o_in,  &
                          tmelt_in, latvap_in, latice_in, rhmini_in, pi_in)
! io
real(r8),         intent(in)  :: pi_in
real(r8),         intent(in)  :: gravit_in
real(r8),         intent(in)  :: rair_in
real(r8),         intent(in)  :: rh2o_in
real(r8),         intent(in)  :: cpair_in
real(r8),         intent(in)  :: rhoh2o_in
real(r8),         intent(in)  :: tmelt_in        ! Freezing point of water (K)
real(r8),         intent(in)  :: latvap_in
real(r8),         intent(in)  :: latice_in
real(r8),         intent(in)  :: rhmini_in       ! Minimum rh for ice cloud fraction > 0.
! local
integer k
integer l,m, iaer
real(r8) surften       ! surface tension of water w/respect to air (N/m)
real(r8) arg
!-----------------------------------------------------------------------

!declarations for morrison codes (transforms variable names)

g= gravit_in                  ! gravity
r= rair_in                    ! Dry air Gas constant: note units(phys_constants are in J/K/kmol)
rv= rh2o_in                   ! water vapor gas contstant
cpp = cpair_in                ! specific heat of dry air
rhow = rhoh2o_in              ! density of liquid water
tmelt = tmelt_in
rhmini = rhmini_in

! latent heats
xxlv = latvap_in              ! latent heat vaporization
xlf = latice_in               ! latent heat freezing
xxls = xxlv + xlf             ! latent heat of sublimation

! parameters for snow/rain fraction for convective clouds
tmax_fsnow = tmelt_in
tmin_fsnow = tmelt_in-5._r8

! parameters below from Reisner et al. (1998)
! density parameters (kg/m3)
rhosn = 250._r8    ! bulk density snow  (++ ceh)
rhoi = 500._r8     ! bulk density ice
rhow = 1000._r8    ! bulk density liquid

! fall speed parameters, V = aD^b
! V is in m/s

! droplets
ac = 3.e7_r8
bc = 2._r8

! snow
as = 11.72_r8
bs = 0.41_r8

! cloud ice
ai = 700._r8
bi = 1._r8

! rain
ar = 841.99667_r8
br = 0.8_r8

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d

pi= pi_in


! cloud ice mass-diameter relationship
ci = rhoi*pi/6._r8
di = 3._r8

! snow mass-diameter relationship
cs = rhosn*pi/6._r8
ds = 3._r8

! drop mass-diameter relationship
cr = rhow*pi/6._r8
dr = 3._r8

! ventilation parameters for snow
! hall and prupacher
f1s = 0.86_r8
f2s = 0.28_r8

! collection efficiency, aggregation of cloud ice and snow
Eii = 0.1_r8

! collection efficiency, accretion of cloud water by rain
Ecr = 1.0_r8

! ventilation constants for rain
f1r = 0.78_r8
f2r = 0.32_r8

! autoconversion size threshold for cloud ice to snow (m)
Dcs = 400.e-6_r8

! smallest mixing ratio considered in microphysics
qsmall = 1.e-18_r8  

! immersion freezing parameters, bigg 1953
bimm = 100._r8
aimm = 0.66_r8

! typical air density at 850 mb
rhosu = 85000._r8/(r * tmelt)

! mass of new crystal due to aerosol freezing and growth (kg)
mi0 = 4._r8/3._r8*pi*rhoi*(10.e-6_r8)*(10.e-6_r8)*(10.e-6_r8)

! radius of contact nuclei aerosol (m)
rin = 0.1e-6_r8

! freezing temperature
tt0=273.15_r8

!Range of cloudsat reflectivities (dBz) for analytic simulator
csmin= -30._r8
csmax= 26._r8
mindbz = -99._r8
!      minrefl = 10._r8**(mindbz/10._r8)
minrefl = 1.26e-10_r8

! Define constants to help speed up code (limit calls to gamma function)
cons1=gamma(1._r8+di)
cons4=gamma(1._r8+br)
cons5=gamma(4._r8+br)
cons6=gamma(1._r8+ds)
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

lammaxi = 1._r8/10.e-6_r8
lammini = 1._r8/(2._r8*dcs)
lammaxr = 1._r8/20.e-6_r8
lamminr = 1._r8/500.e-6_r8
lammaxs = 1._r8/10.e-6_r8
lammins = 1._r8/2000.e-6_r8


end subroutine micro_mg1_0_init


subroutine micro_mg1_0_tend( microp_uniform, pver, ncol, top_lev, deltatin,   &
                             tn, qn, qc, qi, nc,                              &
                             ni, p, pdel, cldn, liqcldf,                      &
                             relvar, accre_enhan,                             &
                             icecldf, cldo, rate1ord_cw2pr_st, naai, npccnin, &
                             rndst, nacon, tlat, qvlat, qctend,               &
                             qitend, nctend, nitend, effc, effc_fn,           &
                             effi, prect, preci, nevapr, evapsnow,            &
                             prain, prodsnow, cmeout, deffi, pgamrad,         &
                             lamcrad, qsout, dsout, rsout, rflx, sflx,        &
                             qrout, reff_rain, reff_snow, qcsevap, qisevap,   &
                             qvres, cmeiout, vtrmc, vtrmi, qcsedten,          &
                             qisedten, prao, prco, mnuccco, mnuccto,          &
                             msacwio, psacwso, bergso, bergo, melto,          &
                             homoo, qcreso, prcio, praio, qireso,             &
                             mnuccro, pracso, meltsdt, frzrdt, mnuccdo,       &
                             nrout, nsout, refl, arefl, areflz,               &
                             frefl, csrfl, acsrfl, fcsrfl, rercld,            &
                             ncai, ncal, qrout2, qsout2, nrout2,              &
                             nsout2, drout2, dsout2, freqs, freqr,            &
                             nfice, do_cldice, tnd_qsnow,                     &
                             tnd_nsnow, re_ice )

use wv_sat_methods, only:    svp_water => wv_sat_svp_water, &
                             svp_ice => wv_sat_svp_ice,     &
                             svp_to_qsat => wv_sat_svp_to_qsat

! io
logical,  intent(in) :: microp_uniform  ! True = configure uniform for sub-columns  False = use w/o sub-columns (standard)
integer,  intent(in) :: pver                 ! number of layers in columns
integer,  intent(in) :: ncol                 ! number of columns
integer,  intent(in) :: top_lev              ! top level microphys is applied
real(r8), intent(in) :: deltatin             ! time step (s)
real(r8), intent(in) :: tn(pver, ncol)       ! input temperature (K)
real(r8), intent(in) :: qn(pver, ncol)       ! input h20 vapor mixing ratio (kg/kg)
real(r8), intent(in) :: relvar(pver, ncol)   ! relative variance of cloud water (-)
real(r8), intent(in) :: accre_enhan(pver, ncol) ! optional accretion enhancement factor (-)

! note: all input cloud variables are grid-averaged
real(r8), intent(inout) :: qc(pver, ncol)    ! cloud water mixing ratio (kg/kg)
real(r8), intent(inout) :: qi(pver, ncol)    ! cloud ice mixing ratio (kg/kg)
real(r8), intent(inout) :: nc(pver, ncol)    ! cloud water number conc (1/kg)
real(r8), intent(inout) :: ni(pver, ncol)    ! cloud ice number conc (1/kg)
real(r8), intent(in) :: p(pver, ncol)        ! air pressure (pa)
real(r8), intent(in) :: pdel(pver, ncol)     ! pressure difference across level (pa)
real(r8), intent(in) :: cldn(pver, ncol)     ! cloud fraction
real(r8), intent(in) :: icecldf(pver, ncol)  ! ice cloud fraction   
real(r8), intent(in) :: liqcldf(pver, ncol)  ! liquid cloud fraction
real(r8), intent(inout) :: cldo(pver, ncol)  ! old cloud fraction
real(r8), intent(out) :: rate1ord_cw2pr_st(pver, ncol) ! 1st order rate for direct cw to precip conversion 
! used for scavenging
! Inputs for aerosol activation
real(r8), intent(in) :: naai(pver, ncol)      ! ice nulceation number (from microp_aero_ts) 
real(r8), intent(in) :: npccnin(pver, ncol)   ! ccn activated number tendency (from microp_aero_ts)
real(r8), intent(in) :: rndst(4,pver, ncol)   ! radius of 4 dust bins for contact freezing (from microp_aero_ts)
real(r8), intent(in) :: nacon(4,pver, ncol)   ! number in 4 dust bins for contact freezing  (from microp_aero_ts)

! Used with CARMA cirrus microphysics
! (or similar external microphysics model)
logical,  intent(in) :: do_cldice             ! Prognosing cldice
real(r8), intent(in) :: tnd_qsnow(pver, ncol) ! snow mass tendency (kg/kg/s)
real(r8), intent(in) :: tnd_nsnow(pver, ncol) ! snow number tendency (#/kg/s)
real(r8), intent(in) :: re_ice(pver, ncol)    ! ice effective radius (m)

! output arguments
real(r8), intent(out) :: tlat(pver, ncol)    ! latent heating rate       (W/kg)
real(r8), intent(out) :: qvlat(pver, ncol)   ! microphysical tendency qv (1/s)
real(r8), intent(out) :: qctend(pver, ncol)  ! microphysical tendency qc (1/s) 
real(r8), intent(out) :: qitend(pver, ncol)  ! microphysical tendency qi (1/s)
real(r8), intent(out) :: nctend(pver, ncol)  ! microphysical tendency nc (1/(kg*s))
real(r8), intent(out) :: nitend(pver, ncol)  ! microphysical tendency ni (1/(kg*s))
real(r8), intent(out) :: effc(pver, ncol)    ! droplet effective radius (micron)
real(r8), intent(out) :: effc_fn(pver, ncol) ! droplet effective radius, assuming nc = 1.e8 kg-1
real(r8), intent(out) :: effi(pver, ncol)    ! cloud ice effective radius (micron)
real(r8), intent(out) :: prect(ncol)         ! surface precip rate (m/s)
real(r8), intent(out) :: preci(ncol)         ! cloud ice/snow precip rate (m/s)
real(r8), intent(out) :: nevapr(pver, ncol)  ! evaporation rate of rain + snow
real(r8), intent(out) :: evapsnow(pver, ncol)! sublimation rate of snow
real(r8), intent(out) :: prain(pver, ncol)   ! production of rain + snow
real(r8), intent(out) :: prodsnow(pver, ncol)! production of snow
real(r8), intent(out) :: cmeout(pver, ncol)  ! evap/sub of cloud
real(r8), intent(out) :: deffi(pver, ncol)   ! ice effective diameter for optics (radiation)
real(r8), intent(out) :: pgamrad(pver, ncol) ! ice gamma parameter for optics (radiation)
real(r8), intent(out) :: lamcrad(pver, ncol) ! slope of droplet distribution for optics (radiation)
real(r8), intent(out) :: qsout(pver, ncol)   ! snow mixing ratio (kg/kg)
real(r8), intent(out) :: dsout(pver, ncol)   ! snow diameter (m)
real(r8), intent(out) :: rsout(pver, ncol)   ! snow radius (m)     Lihan add for RRTMG 4DDA
real(r8), intent(out) :: rflx(pver+1, ncol)  ! grid-box average rain flux (kg m^-2 s^-1)
real(r8), intent(out) :: sflx(pver+1, ncol)  ! grid-box average snow flux (kg m^-2 s^-1)
real(r8), intent(out) :: qrout(pver, ncol)   ! grid-box average rain mixing ratio (kg/kg)
real(r8), intent(inout) :: reff_rain(pver, ncol) ! rain effective radius (micron)
real(r8), intent(inout) :: reff_snow(pver, ncol) ! snow effective radius (micron)
real(r8), intent(out) :: qcsevap(pver, ncol) ! cloud water evaporation due to sedimentation
real(r8), intent(out) :: qisevap(pver, ncol) ! cloud ice sublimation due to sublimation
real(r8), intent(out) :: qvres(pver, ncol)   ! residual condensation term to ensure RH < 100%
real(r8), intent(out) :: cmeiout(pver, ncol) ! grid-mean cloud ice sub/dep
real(r8), intent(out) :: vtrmc(pver, ncol)   ! mass-weighted cloud water fallspeed
real(r8), intent(out) :: vtrmi(pver, ncol)   ! mass-weighted cloud ice fallspeed
real(r8), intent(out) :: qcsedten(pver, ncol)! qc sedimentation tendency
real(r8), intent(out) :: qisedten(pver, ncol)! qi sedimentation tendency
! microphysical process rates for output (mixing ratio tendencies)
real(r8), intent(out) :: prao(pver, ncol) ! accretion of cloud by rain 
real(r8), intent(out) :: prco(pver, ncol) ! autoconversion of cloud to rain
real(r8), intent(out) :: mnuccco(pver, ncol) ! mixing rat tend due to immersion freezing
real(r8), intent(out) :: mnuccto(pver, ncol) ! mixing ratio tend due to contact freezing
real(r8), intent(out) :: msacwio(pver, ncol) ! mixing ratio tend due to H-M splintering
real(r8), intent(out) :: psacwso(pver, ncol) ! collection of cloud water by snow
real(r8), intent(out) :: bergso(pver, ncol) ! bergeron process on snow
real(r8), intent(out) :: bergo(pver, ncol) ! bergeron process on cloud ice
real(r8), intent(out) :: melto(pver, ncol) ! melting of cloud ice
real(r8), intent(out) :: homoo(pver, ncol) ! homogeneos freezign cloud water
real(r8), intent(out) :: qcreso(pver, ncol) ! residual cloud condensation due to removal of excess supersat
real(r8), intent(out) :: prcio(pver, ncol) ! autoconversion of cloud ice to snow
real(r8), intent(out) :: praio(pver, ncol) ! accretion of cloud ice by snow
real(r8), intent(out) :: qireso(pver, ncol) ! residual ice deposition due to removal of excess supersat
real(r8), intent(out) :: mnuccro(pver, ncol) ! mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
real(r8), intent(out) :: pracso (pver, ncol) ! mixing ratio tendency due to accretion of rain by snow (1/s)
real(r8), intent(out) :: meltsdt(pver, ncol) ! latent heating rate due to melting of snow  (W/kg)
real(r8), intent(out) :: frzrdt (pver, ncol) ! latent heating rate due to homogeneous freezing of rain (W/kg)
real(r8), intent(out) :: mnuccdo(pver, ncol) ! mass tendency from ice nucleation
real(r8), intent(out) :: nrout(pver, ncol) ! rain number concentration (1/m3)
real(r8), intent(out) :: nsout(pver, ncol) ! snow number concentration (1/m3)
real(r8), intent(out) :: refl(pver, ncol)    ! analytic radar reflectivity        
real(r8), intent(out) :: arefl(pver, ncol)   ! average reflectivity will zero points outside valid range
real(r8), intent(out) :: areflz(pver, ncol)  ! average reflectivity in z.
real(r8), intent(out) :: frefl(pver, ncol)
real(r8), intent(out) :: csrfl(pver, ncol)  !cloudsat reflectivity 
real(r8), intent(out) :: acsrfl(pver, ncol) !cloudsat average
real(r8), intent(out) :: fcsrfl(pver, ncol)
real(r8), intent(out) :: rercld(pver, ncol) ! effective radius calculation for rain + cloud
real(r8), intent(out) :: ncai(pver, ncol) ! output number conc of ice nuclei available (1/m3)
real(r8), intent(out) :: ncal(pver, ncol) ! output number conc of CCN (1/m3)
real(r8), intent(out) :: qrout2(pver, ncol)
real(r8), intent(out) :: qsout2(pver, ncol)
real(r8), intent(out) :: nrout2(pver, ncol)
real(r8), intent(out) :: nsout2(pver, ncol)
real(r8), intent(out) :: drout2(pver, ncol) ! mean rain particle diameter (m)
real(r8), intent(out) :: dsout2(pver, ncol) ! mean snow particle diameter (m)
real(r8), intent(out) :: freqs(pver, ncol)
real(r8), intent(out) :: freqr(pver, ncol)
real(r8), intent(out) :: nfice(pver, ncol)
! local workspace
! all units mks unless otherwise stated

! temporary variables for sub-stepping 
real(r8) :: t1(pver, ncol)
real(r8) :: q1(pver, ncol)
real(r8) :: qc1(pver, ncol)
real(r8) :: qi1(pver, ncol)
real(r8) :: nc1(pver, ncol)
real(r8) :: ni1(pver, ncol)
real(r8) :: tlat1(pver, ncol)
real(r8) :: qvlat1(pver, ncol)
real(r8) :: qctend1(pver, ncol)
real(r8) :: qitend1(pver, ncol)
real(r8) :: nctend1(pver, ncol)
real(r8) :: nitend1(pver, ncol)
real(r8) :: prect1(ncol)
real(r8) :: preci1(ncol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(r8) :: deltat        ! sub-time step (s)
real(r8) :: omsm    ! number near unity for round-off issues
real(r8) :: dto2    ! dt/2 (s)
real(r8) :: mincld  ! minimum allowed cloud fraction
real(r8) :: q(pver, ncol) ! water vapor mixing ratio (kg/kg)
real(r8) :: t(pver, ncol) ! temperature (K)
real(r8) :: rho(pver, ncol) ! air density (kg m-3)
real(r8) :: dv(pver, ncol)  ! diffusivity of water vapor in air
real(r8) :: mu(pver, ncol)  ! viscocity of air
real(r8) :: sc(pver, ncol)  ! schmidt number
real(r8) :: kap(pver, ncol) ! thermal conductivity of air
real(r8) :: rhof(pver, ncol) ! air density correction factor for fallspeed
real(r8) :: cldmax(pver, ncol) ! precip fraction assuming maximum overlap
real(r8) :: cldm(pver, ncol)   ! cloud fraction
real(r8) :: icldm(pver, ncol)   ! ice cloud fraction
real(r8) :: lcldm(pver, ncol)   ! liq cloud fraction
real(r8) :: icwc(ncol)    ! in cloud water content (liquid+ice)
real(r8) :: calpha(ncol)  ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbeta(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cbetah(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamma(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cgamah(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: rcgama(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec1(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec2(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec3(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: cmec4(ncol) ! parameter for cond/evap (Zhang et al. 2003)
real(r8) :: qtmp ! dummy qv 
real(r8) :: dum  ! temporary dummy variable
real(r8) :: cme(pver, ncol)  ! total (liquid+ice) cond/evap rate of cloud
real(r8) :: cmei(pver, ncol) ! dep/sublimation rate of cloud ice
real(r8) :: cwml(pver, ncol) ! cloud water mixing ratio
real(r8) :: cwmi(pver, ncol) ! cloud ice mixing ratio
real(r8) :: nnuccd(pver)   ! ice nucleation rate from deposition/cond.-freezing
real(r8) :: mnuccd(pver)   ! mass tendency from ice nucleation
real(r8) :: qcld              ! total cloud water
real(r8) :: lcldn(pver, ncol) ! fractional coverage of new liquid cloud
real(r8) :: lcldo(pver, ncol) ! fractional coverage of old liquid cloud
real(r8) :: nctend_mixnuc(pver, ncol)
real(r8) :: arg ! argument of erfc

! for calculation of rate1ord_cw2pr_st
real(r8) :: qcsinksum_rate1ord(pver)   ! sum over iterations of cw to precip sink
real(r8) :: qcsum_rate1ord(pver)    ! sum over iterations of cloud water       

real(r8) :: alpha

real(r8) :: dum1,dum2   !general dummy variables

real(r8) :: npccn(pver)     ! droplet activation rate
real(r8) :: qcic(pver, ncol) ! in-cloud cloud liquid mixing ratio
real(r8) :: qiic(pver, ncol) ! in-cloud cloud ice mixing ratio
real(r8) :: qniic(pver, ncol) ! in-precip snow mixing ratio
real(r8) :: qric(pver, ncol) ! in-precip rain mixing ratio
real(r8) :: ncic(pver, ncol) ! in-cloud droplet number conc
real(r8) :: niic(pver, ncol) ! in-cloud cloud ice number conc
real(r8) :: nsic(pver, ncol) ! in-precip snow number conc
real(r8) :: nric(pver, ncol) ! in-precip rain number conc
real(r8) :: lami(pver) ! slope of cloud ice size distr
real(r8) :: n0i(pver) ! intercept of cloud ice size distr
real(r8) :: lamc(pver) ! slope of cloud liquid size distr
real(r8) :: n0c(pver) ! intercept of cloud liquid size distr
real(r8) :: lams(pver) ! slope of snow size distr
real(r8) :: n0s(pver) ! intercept of snow size distr
real(r8) :: lamr(pver) ! slope of rain size distr
real(r8) :: n0r(pver) ! intercept of rain size distr
real(r8) :: cdist1(pver) ! size distr parameter to calculate droplet freezing
! combined size of precip & cloud drops
real(r8) :: arcld(pver, ncol) ! averaging control flag
real(r8) :: Actmp  !area cross section of drops
real(r8) :: Artmp  !area cross section of rain

real(r8) :: pgam(pver) ! spectral width parameter of droplet size distr
real(r8) :: lammax  ! maximum allowed slope of size distr
real(r8) :: lammin  ! minimum allowed slope of size distr
real(r8) :: nacnt   ! number conc of contact ice nuclei
real(r8) :: mnuccc(pver) ! mixing ratio tendency due to freezing of cloud water
real(r8) :: nnuccc(pver) ! number conc tendency due to freezing of cloud water

real(r8) :: mnucct(pver) ! mixing ratio tendency due to contact freezing of cloud water
real(r8) :: nnucct(pver) ! number conc tendency due to contact freezing of cloud water
real(r8) :: msacwi(pver) ! mixing ratio tendency due to HM ice multiplication
real(r8) :: nsacwi(pver) ! number conc tendency due to HM ice multiplication

real(r8) :: prc(pver) ! qc tendency due to autoconversion of cloud droplets
real(r8) :: nprc(pver) ! number conc tendency due to autoconversion of cloud droplets
real(r8) :: nprc1(pver) ! qr tendency due to autoconversion of cloud droplets
real(r8) :: nsagg(pver) ! ns tendency due to self-aggregation of snow
real(r8) :: dc0  ! mean size droplet size distr
real(r8) :: ds0  ! mean size snow size distr (area weighted)
real(r8) :: eci  ! collection efficiency for riming of snow by droplets
real(r8) :: psacws(pver) ! mixing rat tendency due to collection of droplets by snow
real(r8) :: npsacws(pver) ! number conc tendency due to collection of droplets by snow
real(r8) :: uni ! number-weighted cloud ice fallspeed
real(r8) :: umi ! mass-weighted cloud ice fallspeed
real(r8) :: uns(pver) ! number-weighted snow fallspeed
real(r8) :: ums(pver) ! mass-weighted snow fallspeed
real(r8) :: unr(pver) ! number-weighted rain fallspeed
real(r8) :: umr(pver) ! mass-weighted rain fallspeed
real(r8) :: unc ! number-weighted cloud droplet fallspeed
real(r8) :: umc ! mass-weighted cloud droplet fallspeed
real(r8) :: pracs(pver) ! mixing rat tendency due to collection of rain	by snow
real(r8) :: npracs(pver) ! number conc tendency due to collection of rain by snow
real(r8) :: mnuccr(pver) ! mixing rat tendency due to freezing of rain
real(r8) :: nnuccr(pver) ! number conc tendency due to freezing of rain
real(r8) :: pra(pver) ! mixing rat tendnency due to accretion of droplets by rain
real(r8) :: npra(pver) ! nc tendnency due to accretion of droplets by rain
real(r8) :: nragg(pver) ! nr tendency due to self-collection of rain
real(r8) :: prci(pver) ! mixing rat tendency due to autoconversion of cloud ice to snow
real(r8) :: nprci(pver) ! number conc tendency due to autoconversion of cloud ice to snow
real(r8) :: prai(pver) ! mixing rat tendency due to accretion of cloud ice by snow
real(r8) :: nprai(pver) ! number conc tendency due to accretion of cloud ice by snow
real(r8) :: qvs ! liquid saturation vapor mixing ratio
real(r8) :: qvi ! ice saturation vapor mixing ratio
real(r8) :: dqsdt ! change of sat vapor mixing ratio with temperature
real(r8) :: dqsidt ! change of ice sat vapor mixing ratio with temperature
real(r8) :: ab ! correction factor for rain evap to account for latent heat
real(r8) :: qclr ! water vapor mixing ratio in clear air
real(r8) :: abi ! correction factor for snow sublimation to account for latent heat
real(r8) :: epss ! 1/ sat relaxation timescale for snow
real(r8) :: epsr ! 1/ sat relaxation timescale for rain
real(r8) :: pre(pver) ! rain mixing rat tendency due to evaporation
real(r8) :: prds(pver) ! snow mixing rat tendency due to sublimation
real(r8) :: qce ! dummy qc for conservation check
real(r8) :: qie ! dummy qi for conservation check
real(r8) :: nce ! dummy nc for conservation check
real(r8) :: nie ! dummy ni for conservation check
real(r8) :: ratio ! parameter for conservation check
real(r8) :: dumc(pver, ncol) ! dummy in-cloud qc
real(r8) :: dumnc(pver, ncol) ! dummy in-cloud nc
real(r8) :: dumi(pver, ncol) ! dummy in-cloud qi
real(r8) :: dumni(pver, ncol) ! dummy in-cloud ni
real(r8) :: dums(pver, ncol) ! dummy in-cloud snow mixing rat
real(r8) :: dumns(pver, ncol) ! dummy in-cloud snow number conc
real(r8) :: dumr(pver, ncol) ! dummy in-cloud rain mixing rat
real(r8) :: dumnr(pver, ncol) ! dummy in-cloud rain number conc
! below are parameters for cloud water and cloud ice sedimentation calculations
real(r8) :: fr(pver)
real(r8) :: fnr(pver)
real(r8) :: fc(pver)
real(r8) :: fnc(pver)
real(r8) :: fi(pver)
real(r8) :: fni(pver)
real(r8) :: fs(pver)
real(r8) :: fns(pver)
real(r8) :: faloutr(pver)
real(r8) :: faloutnr(pver)
real(r8) :: faloutc(pver)
real(r8) :: faloutnc(pver)
real(r8) :: falouti(pver)
real(r8) :: faloutni(pver)
real(r8) :: falouts(pver)
real(r8) :: faloutns(pver)
real(r8) :: faltndr
real(r8) :: faltndnr
real(r8) :: faltndc
real(r8) :: faltndnc
real(r8) :: faltndi
real(r8) :: faltndni
real(r8) :: faltnds
real(r8) :: faltndns
real(r8) :: faltndqie
real(r8) :: faltndqce
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(r8) :: relhum(pver, ncol) ! relative humidity
real(r8) :: csigma(ncol) ! parameter for cond/evap of cloud water/ice
real(r8) :: rgvm ! max fallspeed for all species
real(r8) :: arn(pver, ncol) ! air density corrected rain fallspeed parameter
real(r8) :: asn(pver, ncol) ! air density corrected snow fallspeed parameter
real(r8) :: acn(pver, ncol) ! air density corrected cloud droplet fallspeed parameter
real(r8) :: ain(pver, ncol) ! air density corrected cloud ice fallspeed parameter
real(r8) :: nsubi(pver) ! evaporation of cloud ice number
real(r8) :: nsubc(pver) ! evaporation of droplet number
real(r8) :: nsubs(pver) ! evaporation of snow number
real(r8) :: nsubr(pver) ! evaporation of rain number
real(r8) :: mtime ! factor to account for droplet activation timescale
real(r8) :: dz(pver, ncol) ! height difference across model vertical level

!! add precip flux variables for sub-stepping
real(r8) :: rflx1(pver+1, ncol)
real(r8) :: sflx1(pver+1, ncol)

! returns from function/subroutine calls
real(r8) :: tsp(pver, ncol)      ! saturation temp (K)
real(r8) :: qsp(pver, ncol)      ! saturation mixing ratio (kg/kg)
real(r8) :: qsphy(pver, ncol)    ! saturation mixing ratio (kg/kg): hybrid rh
real(r8) :: qs(ncol)             ! liquid-ice weighted sat mixing rat (kg/kg)
real(r8) :: es(ncol)             ! liquid-ice weighted sat vapor press (pa)
real(r8) :: esl(pver, ncol)      ! liquid sat vapor pressure (pa)
real(r8) :: esi(pver, ncol)      ! ice sat vapor pressure (pa)

! sum of source/sink terms for diagnostic precip

real(r8) :: qnitend(pver, ncol) ! snow mixing ratio source/sink term
real(r8) :: nstend(pver, ncol)  ! snow number concentration source/sink term
real(r8) :: qrtend(pver, ncol) ! rain mixing ratio source/sink term
real(r8) :: nrtend(pver, ncol)  ! rain number concentration source/sink term
real(r8) :: qrtot ! vertically-integrated rain mixing rat source/sink term
real(r8) :: nrtot ! vertically-integrated rain number conc source/sink term
real(r8) :: qstot ! vertically-integrated snow mixing rat source/sink term
real(r8) :: nstot ! vertically-integrated snow number conc source/sink term

! new terms for Bergeron process

real(r8) :: dumnnuc ! provisional ice nucleation rate (for calculating bergeron)
real(r8) :: ninew  ! provisional cloud ice number conc (for calculating bergeron)
real(r8) :: qinew ! provisional cloud ice mixing ratio (for calculating bergeron)
real(r8) :: qvl  ! liquid sat mixing ratio   
real(r8) :: epsi ! 1/ sat relaxation timecale for cloud ice
real(r8) :: prd ! provisional deposition rate of cloud ice at water sat 
real(r8) :: berg(pver, ncol) ! mixing rat tendency due to bergeron process for cloud ice
real(r8) :: bergs(pver) ! mixing rat tendency due to bergeron process for snow

!bergeron terms
real(r8) :: bergtsf   !bergeron timescale to remove all liquid
real(r8) :: rhin      !modified RH for vapor deposition

! diagnostic rain/snow for output to history
! values are in-precip (local) !!!!
real(r8) :: drout(pver, ncol)     ! rain diameter (m)

!averageed rain/snow for history
real(r8) :: dumfice

!ice nucleation, droplet activation
real(r8) :: dum2i(pver, ncol) ! number conc of ice nuclei available (1/kg)
real(r8) :: dum2l(pver, ncol) ! number conc of CCN (1/kg)
real(r8) :: ncmax
real(r8) :: nimax

real(r8) :: qcvar     ! 1/relative variance of sub-grid qc

! loop array variables
integer i,k,nstep,n, l
integer ii,kk, m

! loop variables for sub-step solution
integer iter,it,ltrue(ncol)

! used in contact freezing via dust particles
real(r8)  tcnt, viscosity, mfp
real(r8)  slip1, slip2, slip3, slip4
!        real(r8)  dfaer1, dfaer2, dfaer3, dfaer4
!        real(r8)  nacon1,nacon2,nacon3,nacon4
real(r8)  ndfaer1, ndfaer2, ndfaer3, ndfaer4
real(r8)  nslip1, nslip2, nslip3, nslip4

! used in ice effective radius
real(r8)  bbi, cci, ak, iciwc, rvi

! used in Bergeron processe and water vapor deposition
real(r8)  Tk, deles, Aprpr, Bprpr, Cice, qi0, Crate, qidep

! mean cloud fraction over the time step
real(r8)  cldmw(pver, ncol)

! used in secondary ice production
real(r8) ni_secp

! variabels to check for RH after rain evap

real(r8) :: esn
real(r8) :: qsn
real(r8) :: ttmp

real(r8) :: rainrt(pver, ncol)  ! rain rate for reflectivity calculation
real(r8) :: rainrt1(pver, ncol)
real(r8) :: tmp

real(r8) dmc,ssmc,dstrn  ! variables for modal scheme.
real(r8), parameter :: cdnl    = 0.e6_r8    ! cloud droplet number limiter

! initialize  output fields for number conc qand ice nucleation
ncai(1:pver,1:ncol)=0._r8 
ncal(1:pver,1:ncol)=0._r8  

!Initialize rain size
rercld(1:pver,1:ncol)=0._r8
arcld(1:pver,1:ncol)=0._r8

!initialize radiation output variables
pgamrad(1:pver,1:ncol)=0._r8 ! liquid gamma parameter for optics (radiation)
lamcrad(1:pver,1:ncol)=0._r8 ! slope of droplet distribution for optics (radiation)
deffi  (1:pver,1:ncol)=0._r8 ! slope of droplet distribution for optics (radiation)
!initialize radiation output variables
!initialize water vapor tendency term output
qcsevap(1:pver,1:ncol)=0._r8 
qisevap(1:pver,1:ncol)=0._r8 
qvres  (1:pver,1:ncol)=0._r8 
cmeiout (1:pver,1:ncol)=0._r8
vtrmc (1:pver,1:ncol)=0._r8
vtrmi (1:pver,1:ncol)=0._r8
qcsedten (1:pver,1:ncol)=0._r8
qisedten (1:pver,1:ncol)=0._r8    

prao(1:pver,1:ncol)=0._r8 
prco(1:pver,1:ncol)=0._r8 
mnuccco(1:pver,1:ncol)=0._r8 
mnuccto(1:pver,1:ncol)=0._r8 
msacwio(1:pver,1:ncol)=0._r8 
psacwso(1:pver,1:ncol)=0._r8 
bergso(1:pver,1:ncol)=0._r8 
bergo(1:pver,1:ncol)=0._r8 
melto(1:pver,1:ncol)=0._r8 
homoo(1:pver,1:ncol)=0._r8 
qcreso(1:pver,1:ncol)=0._r8 
prcio(1:pver,1:ncol)=0._r8 
praio(1:pver,1:ncol)=0._r8 
qireso(1:pver,1:ncol)=0._r8 
mnuccro(1:pver,1:ncol)=0._r8 
pracso (1:pver,1:ncol)=0._r8 
meltsdt(1:pver,1:ncol)=0._r8
frzrdt (1:pver,1:ncol)=0._r8
mnuccdo(1:pver,1:ncol)=0._r8

rflx(:,:)=0._r8
sflx(:,:)=0._r8
effc(:,:)=0._r8
effc_fn(:,:)=0._r8
effi(:,:)=0._r8

! assign variable deltat for sub-stepping...
deltat=deltatin

! parameters for scheme
omsm=0.99999_r8
dto2=0.5_r8*deltat
mincld=0.0001_r8

! initialize multi-level fields
q(1:pver,1:ncol)=qn(1:pver,1:ncol)
t(1:pver,1:ncol)=tn(1:pver,1:ncol)

! initialize time-varying parameters

do k=1,pver
   do i=1,ncol
      rho(k,i)=p(k,i)/(r*t(k,i))
      dv(k,i) = 8.794E-5_r8*t(k,i)**1.81_r8/p(k,i)
      mu(k,i) = 1.496E-6_r8*t(k,i)**1.5_r8/(t(k,i)+120._r8) 
      sc(k,i) = mu(k,i)/(rho(k,i)*dv(k,i))
      kap(k,i) = 1.414e3_r8*1.496e-6_r8*t(k,i)**1.5_r8/(t(k,i)+120._r8) 

      ! air density adjustment for fallspeed parameters
      ! includes air density correction factor to the
      ! power of 0.54 following Heymsfield and Bansemer 2007

      rhof(k,i)=(rhosu/rho(k,i))**0.54_r8

      arn(k,i)=ar*rhof(k,i)
      asn(k,i)=as*rhof(k,i)
      acn(k,i)=ac*rhof(k,i)
      ain(k,i)=ai*rhof(k,i)

      ! get dz from dp and hydrostatic approx
      ! keep dz positive (define as layer k-1 - layer k)

      dz(k,i)= pdel(k,i)/(rho(k,i)*g)

   end do
end do

! initialization
qc(1:top_lev-1,1:ncol) = 0._r8
qi(1:top_lev-1,1:ncol) = 0._r8
nc(1:top_lev-1,1:ncol) = 0._r8
ni(1:top_lev-1,1:ncol) = 0._r8
t1(1:pver,1:ncol) = t(1:pver,1:ncol)
q1(1:pver,1:ncol) = q(1:pver,1:ncol)
qc1(1:pver,1:ncol) = qc(1:pver,1:ncol)
qi1(1:pver,1:ncol) = qi(1:pver,1:ncol)
nc1(1:pver,1:ncol) = nc(1:pver,1:ncol)
ni1(1:pver,1:ncol) = ni(1:pver,1:ncol)

! initialize tendencies to zero
tlat1(1:pver,1:ncol)=0._r8
qvlat1(1:pver,1:ncol)=0._r8
qctend1(1:pver,1:ncol)=0._r8
qitend1(1:pver,1:ncol)=0._r8
nctend1(1:pver,1:ncol)=0._r8
nitend1(1:pver,1:ncol)=0._r8

! initialize precip output
qrout(1:pver,1:ncol)=0._r8
qsout(1:pver,1:ncol)=0._r8
nrout(1:pver,1:ncol)=0._r8
nsout(1:pver,1:ncol)=0._r8
rsout(1:pver,1:ncol)=0._r8
dsout(1:pver,1:ncol)=0._r8

drout(1:pver,1:ncol)=0._r8

reff_rain(1:pver,1:ncol)=0._r8
reff_snow(1:pver,1:ncol)=0._r8

! initialize variables for trop_mozart
nevapr(1:pver,1:ncol) = 0._r8
evapsnow(1:pver,1:ncol) = 0._r8
prain(1:pver,1:ncol) = 0._r8
prodsnow(1:pver,1:ncol) = 0._r8
cmeout(1:pver,1:ncol) = 0._r8

! for refl calc
rainrt1(1:pver,1:ncol) = 0._r8

! initialize precip fraction and output tendencies
cldmax(1:pver,1:ncol)=mincld

!initialize aerosol number
!        naer2(1:ncol,1:pver,:)=0._r8
dum2l(1:pver,1:ncol)=0._r8
dum2i(1:pver,1:ncol)=0._r8

! initialize avg precip rate
prect1(1:ncol)=0._r8
preci1(1:ncol)=0._r8

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Get humidity and saturation vapor pressures

do k=top_lev,pver
   do i=1,ncol

      ! find wet bulk temperature and saturation value for provisional t and q without
      ! condensation
      
      es(i) = svp_water(t(k,i))
      qs(i) = svp_to_qsat(es(i), p(k,i))

      ! Prevents negative values.
      if (qs(i) < 0.0_r8) then
         qs(i) = 1.0_r8
         es(i) = p(k,i)
      end if

      esl(k,i)=svp_water(t(k,i))
      esi(k,i)=svp_ice(t(k,i))

      ! hm fix, make sure when above freezing that esi=esl, not active yet
      if (t(k,i).gt.tmelt)esi(k,i)=esl(k,i)

      relhum(k,i)=q(k,i)/qs(i)

      ! get cloud fraction, check for minimum

      cldm(k,i)=max(cldn(k,i),mincld)
      cldmw(k,i)=max(cldn(k,i),mincld)

      icldm(k,i)=max(icecldf(k,i),mincld)
      lcldm(k,i)=max(liqcldf(k,i),mincld)

      ! subcolumns, set cloud fraction variables to one
      ! if cloud water or ice is present, if not present
      ! set to mincld (mincld used instead of zero, to prevent
      ! possible division by zero errors

      if (microp_uniform) then

         cldm(k,i)=mincld
         cldmw(k,i)=mincld
         icldm(k,i)=mincld
         lcldm(k,i)=mincld

         if (qc(k,i).ge.qsmall) then
            lcldm(k,i)=1._r8           
            cldm(k,i)=1._r8
            cldmw(k,i)=1._r8
         end if

         if (qi(k,i).ge.qsmall) then             
            cldm(k,i)=1._r8
            icldm(k,i)=1._r8
         end if

      end if               ! sub-columns

      ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)

      nfice(k,i)=0._r8
      dumfice=qc(k,i)+qi(k,i)
      if (dumfice.gt.qsmall .and. qi(k,i).gt.qsmall) then
         nfice(k,i)=qi(k,i)/dumfice
      endif

      if (do_cldice .and. (t(k,i).lt.tmelt - 5._r8)) then

         ! if aerosols interact with ice set number of activated ice nuclei
         dum2=naai(k,i)

         dumnnuc=(dum2-ni(k,i)/icldm(k,i))/deltat*icldm(k,i)
         dumnnuc=max(dumnnuc,0._r8)
         ! get provisional ni and qi after nucleation in order to calculate
         ! Bergeron process below
         ninew=ni(k,i)+dumnnuc*deltat
         qinew=qi(k,i)+dumnnuc*deltat*mi0

         !T>268
      else
         ninew=ni(k,i)
         qinew=qi(k,i)
      end if

      ! Initialize CME components

      cme(k,i) = 0._r8
      cmei(k,i)=0._r8

      !-------------------------------------------------------------------
      !Bergeron process

      ! make sure to initialize bergeron process to zero
      berg(k,i)=0._r8
      prd = 0._r8

      !condensation loop.

      ! get in-cloud qi and ni after nucleation
      if (icldm(k,i) .gt. 0._r8) then 
         qiic(k,i)=qinew/icldm(k,i)
         niic(k,i)=ninew/icldm(k,i)
      else
         qiic(k,i)=0._r8
         niic(k,i)=0._r8
      endif

      !if T < 0 C then bergeron.
      if (do_cldice .and. (t(k,i).lt.273.15_r8)) then

         !if ice exists
         if (qi(k,i).gt.qsmall) then

            bergtsf = 0._r8 ! bergeron time scale (fraction of timestep)

            qvi = svp_to_qsat(esi(k,i), p(k,i))
            qvl = svp_to_qsat(esl(k,i), p(k,i))

            dqsidt =  xxls*qvi/(rv*t(k,i)**2)
            abi = 1._r8+dqsidt*xxls/cpp

            ! get ice size distribution parameters

            if (qiic(k,i).ge.qsmall) then
               lami(k) = (cons1*ci* &
                    niic(k,i)/qiic(k,i))**(1._r8/di)
               n0i(k) = niic(k,i)*lami(k)

               ! check for slope
               ! adjust vars
               if (lami(k).lt.lammini) then
                  lami(k) = lammini
                  n0i(k) = lami(k)**(di+1._r8)*qiic(k,i)/(ci*cons1)
               else if (lami(k).gt.lammaxi) then
                  lami(k) = lammaxi
                  n0i(k) = lami(k)**(di+1._r8)*qiic(k,i)/(ci*cons1)
               end if

               epsi = 2._r8*pi*n0i(k)*rho(k,i)*Dv(k,i)/(lami(k)*lami(k))

               !if liquid exists  
               if (qc(k,i).gt. qsmall) then 

                  !begin bergeron process
                  !     do bergeron (vapor deposition with RHw=1)
                  !     code to find berg (a rate) goes here

                  ! calculate Bergeron process

                  prd = epsi*(qvl-qvi)/abi

               else
                  prd = 0._r8
               end if

               ! multiply by cloud fraction

               prd = prd*min(icldm(k,i),lcldm(k,i))

               !     transfer of existing cloud liquid to ice

               berg(k,i)=max(0._r8,prd)

            end if  !end liquid exists bergeron

            if (berg(k,i).gt.0._r8) then
               bergtsf=max(0._r8,(qc(k,i)/berg(k,i))/deltat) 

               if(bergtsf.lt.1._r8) berg(k,i) = max(0._r8,qc(k,i)/deltat)

            endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (bergtsf.lt.1._r8.or.icldm(k,i).gt.lcldm(k,i)) then

               if (qiic(k,i).ge.qsmall) then

                  ! first case is for case when liquid water is present, but is completely depleted in time step, i.e., bergrsf > 0 but < 1

                  if (qc(k,i).ge.qsmall) then
                     rhin  = (1.0_r8 + relhum(k,i)) / 2._r8
                     if ((rhin*esl(k,i)/esi(k,i)) > 1._r8) then
                        prd = epsi*(rhin*qvl-qvi)/abi

                        ! multiply by cloud fraction assuming liquid/ice maximum overlap
                        prd = prd*min(icldm(k,i),lcldm(k,i))

                        ! add to cmei
                        cmei(k,i) = cmei(k,i) + (prd * (1._r8- bergtsf))

                     end if ! rhin 
                  end if ! qc > qsmall

                  ! second case is for pure ice cloud, either no liquid, or icldm > lcldm

                  if (qc(k,i).lt.qsmall.or.icldm(k,i).gt.lcldm(k,i)) then

                     ! note: for case of no liquid, need to set liquid cloud fraction to zero
                     ! store liquid cloud fraction in 'dum'

                     if (qc(k,i).lt.qsmall) then 
                        dum=0._r8 
                     else
                        dum=lcldm(k,i)
                     end if

                     ! set RH to grid-mean value for pure ice cloud
                     rhin = relhum(k,i)

                     if ((rhin*esl(k,i)/esi(k,i)) > 1._r8) then

                        prd = epsi*(rhin*qvl-qvi)/abi

                        ! multiply by relevant cloud fraction for pure ice cloud
                        ! assuming maximum overlap of liquid/ice
                        prd = prd*max((icldm(k,i)-dum),0._r8)
                        cmei(k,i) = cmei(k,i) + prd

                     end if ! rhin
                  end if ! qc or icldm > lcldm
               end if ! qiic
            end if ! bergtsf or icldm > lcldm

            !     if deposition, it should not reduce grid mean rhi below 1.0
            if(cmei(k,i) > 0.0_r8 .and. (relhum(k,i)*esl(k,i)/esi(k,i)) > 1._r8 ) &
                 cmei(k,i)=min(cmei(k,i),(q(k,i)-qs(i)*esi(k,i)/esl(k,i))/abi/deltat)

         end if            !end ice exists loop
         !this ends temperature < 0. loop

         !-------------------------------------------------------------------
      end if  ! 
      !..............................................................

      ! evaporation should not exceed available water

      if ((-berg(k,i)).lt.-qc(k,i)/deltat) berg(k,i) = max(qc(k,i)/deltat,0._r8)

      !sublimation process...
      if (do_cldice .and. ((relhum(k,i)*esl(k,i)/esi(k,i)).lt.1._r8 .and. qiic(k,i).ge.qsmall )) then

         qvi = svp_to_qsat(esi(k,i), p(k,i))
         qvl = svp_to_qsat(esl(k,i), p(k,i))
         dqsidt =  xxls*qvi/(rv*t(k,i)**2)
         abi = 1._r8+dqsidt*xxls/cpp

         ! get ice size distribution parameters

         lami(k) = (cons1*ci* niic(k,i)/qiic(k,i))**(1._r8/di)

         n0i(k) = niic(k,i)*lami(k)

         ! check for slope
         ! adjust vars
         if (lami(k).lt.lammini) then
            lami(k) = lammini
            n0i(k) = lami(k)**(di+1._r8)*qiic(k,i)/(ci*cons1)
         else if (lami(k).gt.lammaxi) then
            lami(k) = lammaxi
            n0i(k) = lami(k)**(di+1._r8)*qiic(k,i)/(ci*cons1)
         end if

         epsi = 2._r8*pi*n0i(k)*rho(k,i)*Dv(k,i)/(lami(k)*lami(k))

         ! modify for ice fraction below
         prd = epsi*(relhum(k,i)*qvl-qvi)/abi * icldm(k,i)
         cmei(k,i)=min(prd,0._r8)

      endif

      ! sublimation should not exceed available ice
      if (cmei(k,i).lt.-qi(k,i)/deltat) cmei(k,i)=-qi(k,i)/deltat

      ! sublimation should not increase grid mean rhi above 1.0 
      if(cmei(k,i) < 0.0_r8 .and. (relhum(k,i)*esl(k,i)/esi(k,i)) < 1._r8 ) &
           cmei(k,i)=min(0._r8,max(cmei(k,i),(q(k,i)-qs(i)*esi(k,i)/esl(k,i))/abi/deltat))

      ! limit cmei due for roundoff error

      cmei(k,i)=cmei(k,i)*omsm

      ! conditional for ice nucleation 
      if (do_cldice .and. (t(k,i).lt.(tmelt - 5._r8))) then 

         ! using Liu et al. (2007) ice nucleation with hooks into simulated aerosol
         ! ice nucleation rate (dum2) has already been calculated and read in (naai)

         dum2i(k,i)=naai(k,i)
      else
         dum2i(k,i)=0._r8
      end if

   end do ! i loop
end do ! k loop

cldo(:,:ncol)=cldn(:,:ncol)
!! initialize sub-step precip flux variables
do i=1,ncol
   !! flux is zero at top interface, so these should stay as 0.
   rflx1(1,i)=0._r8
   sflx1(1,i)=0._r8
   do k=top_lev,pver

      ! initialize normal and sub-step precip flux variables
      rflx1(k+1,i)=0._r8
      sflx1(k+1,i)=0._r8
   end do ! i loop
end do ! k loop
!! initialize final precip flux variables.
do i=1,ncol
   !! flux is zero at top interface, so these should stay as 0.
   rflx(1,i)=0._r8
   sflx(1,i)=0._r8
   do k=top_lev,pver
      ! initialize normal and sub-step precip flux variables
      rflx(k+1,i)=0._r8
      sflx(k+1,i)=0._r8
   end do ! i loop
end do ! k loop	

do i=1,ncol
   ltrue(i)=0
   do k=top_lev,pver
      ! skip microphysical calculations if no cloud water

      if (qc(k,i).ge.qsmall.or.qi(k,i).ge.qsmall.or.cmei(k,i).ge.qsmall) ltrue(i)=1
   end do
end do

! assign number of sub-steps to iter
! use 2 sub-steps, following tests described in MG2008
iter = 2

! get sub-step time step
deltat=deltat/real(iter)

! since activation/nucleation processes are fast, need to take into account
! factor mtime = mixing timescale in cloud / model time step
! mixing time can be interpreted as cloud depth divided by sub-grid vertical velocity
! for now mixing timescale is assumed to be 1 timestep for modal aerosols, 20 min bulk

!        note: mtime for bulk aerosols was set to: mtime=deltat/1200._r8

mtime=1._r8
rate1ord_cw2pr_st(:,:)=0._r8 ! rce 2010/05/01

!!!! skip calculations if no cloud water
do i=1,ncol
   if (ltrue(i).eq.0) then
      tlat(1:pver,i)=0._r8
      qvlat(1:pver,i)=0._r8
      qctend(1:pver,i)=0._r8
      qitend(1:pver,i)=0._r8
      qnitend(1:pver,i)=0._r8
      qrtend(1:pver,i)=0._r8
      nctend(1:pver,i)=0._r8
      nitend(1:pver,i)=0._r8
      nrtend(1:pver,i)=0._r8
      nstend(1:pver,i)=0._r8
      prect(i)=0._r8
      preci(i)=0._r8
      qniic(1:pver,i)=0._r8
      qric(1:pver,i)=0._r8
      nsic(1:pver,i)=0._r8
      nric(1:pver,i)=0._r8
      rainrt(1:pver,i)=0._r8
      goto 300
   end if

   qcsinksum_rate1ord(1:pver)=0._r8 
   qcsum_rate1ord(1:pver)=0._r8 


!!!!!!!!! begin sub-step!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !.....................................................................................................
   do it=1,iter

      ! initialize sub-step microphysical tendencies

      tlat(1:pver,i)=0._r8
      qvlat(1:pver,i)=0._r8
      qctend(1:pver,i)=0._r8
      qitend(1:pver,i)=0._r8
      qnitend(1:pver,i)=0._r8
      qrtend(1:pver,i)=0._r8
      nctend(1:pver,i)=0._r8
      nitend(1:pver,i)=0._r8
      nrtend(1:pver,i)=0._r8
      nstend(1:pver,i)=0._r8

      ! initialize diagnostic precipitation to zero

      qniic(1:pver,i)=0._r8
      qric(1:pver,i)=0._r8
      nsic(1:pver,i)=0._r8
      nric(1:pver,i)=0._r8

      rainrt(1:pver,i)=0._r8

      ! begin new i,k loop, calculate new cldmax after adjustment to cldm above

      ! initialize vertically-integrated rain and snow tendencies

      qrtot = 0._r8
      nrtot = 0._r8
      qstot = 0._r8
      nstot = 0._r8

      ! initialize precip at surface

      prect(i)=0._r8
      preci(i)=0._r8

      do k=top_lev,pver
      
         qcvar=relvar(k,i)
         cons2=gamma(qcvar+2.47_r8)
         cons3=gamma(qcvar)
         cons9=gamma(qcvar+2._r8)
         cons10=gamma(qcvar+1._r8)
         cons12=gamma(qcvar+1.15_r8) 
         cons15=gamma(qcvar+bc/3._r8)
         cons18=qcvar**2.47_r8
         cons19=qcvar**2
         cons20=qcvar**1.15_r8

         ! set cwml and cwmi to current qc and qi

         cwml(k,i)=qc(k,i)
         cwmi(k,i)=qi(k,i)

         ! initialize precip fallspeeds to zero

         ums(k)=0._r8 
         uns(k)=0._r8 
         umr(k)=0._r8 
         unr(k)=0._r8

         ! calculate precip fraction based on maximum overlap assumption

         ! for sub-columns cldm has already been set to 1 if cloud
         ! water or ice is present, so cldmax will be correctly set below
         ! and nothing extra needs to be done here

         if (k.eq.top_lev) then
            cldmax(k,i)=cldm(k,i)
         else
            ! if rain or snow mix ratio is smaller than
            ! threshold, then set cldmax to cloud fraction at current level
            if (qric(k-1,i).ge.qsmall.or.qniic(k-1,i).ge.qsmall) then
               cldmax(k,i)=max(cldmax(k-1,i),cldm(k,i))
            else
               cldmax(k,i)=cldm(k,i)
            end if
         end if

         ! decrease in number concentration due to sublimation/evap
         ! divide by cloud fraction to get in-cloud decrease
         ! don't reduce Nc due to bergeron process

         if (cmei(k,i) < 0._r8 .and. qi(k,i) > qsmall .and. cldm(k,i) > mincld) then
            nsubi(k)=cmei(k,i)/qi(k,i)*ni(k,i)/cldm(k,i)
         else
            nsubi(k)=0._r8
         end if
         nsubc(k)=0._r8

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%

         if (do_cldice .and. dum2i(k,i).gt.0._r8.and.t(k,i).lt.(tmelt - 5._r8).and. &
              relhum(k,i)*esl(k,i)/esi(k,i).gt. rhmini+0.05_r8) then

            !if NCAI > 0. then set numice = ncai (as before)
            !note: this is gridbox averaged

            nnuccd(k)=(dum2i(k,i)-ni(k,i)/icldm(k,i))/deltat*icldm(k,i)
            nnuccd(k)=max(nnuccd(k),0._r8)
            nimax = dum2i(k,i)*icldm(k,i)

            !Calc mass of new particles using new crystal mass...
            !also this will be multiplied by mtime as nnuccd is...

            mnuccd(k) = nnuccd(k) * mi0

            !  add mnuccd to cmei....
            cmei(k,i)= cmei(k,i) + mnuccd(k) * mtime

            !  limit cmei
            qvi = svp_to_qsat(esi(k,i), p(k,i))
            dqsidt =  xxls*qvi/(rv*t(k,i)**2)
            abi = 1._r8+dqsidt*xxls/cpp
            cmei(k,i)=min(cmei(k,i),(q(k,i)-qvi)/abi/deltat)

            ! limit for roundoff error
            cmei(k,i)=cmei(k,i)*omsm

         else
            nnuccd(k)=0._r8
            nimax = 0._r8
            mnuccd(k) = 0._r8
         end if

         !c............................................................................
         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
         ! for microphysical process calculations
         ! units are kg/kg for mixing ratio, 1/kg for number conc

         ! limit in-cloud values to 0.005 kg/kg

         qcic(k,i)=min(cwml(k,i)/lcldm(k,i),5.e-3_r8)
         qiic(k,i)=min(cwmi(k,i)/icldm(k,i),5.e-3_r8)
         ncic(k,i)=max(nc(k,i)/lcldm(k,i),0._r8)
         niic(k,i)=max(ni(k,i)/icldm(k,i),0._r8)

         if (qc(k,i) - berg(k,i)*deltat.lt.qsmall) then
            qcic(k,i)=0._r8
            ncic(k,i)=0._r8
            if (qc(k,i)-berg(k,i)*deltat.lt.0._r8) then
               berg(k,i)=qc(k,i)/deltat*omsm
            end if
         end if

         if (do_cldice .and. qi(k,i)+(cmei(k,i)+berg(k,i))*deltat.lt.qsmall) then
            qiic(k,i)=0._r8
            niic(k,i)=0._r8
            if (qi(k,i)+(cmei(k,i)+berg(k,i))*deltat.lt.0._r8) then
               cmei(k,i)=(-qi(k,i)/deltat-berg(k,i))*omsm
            end if
         end if

         ! add to cme output

         cmeout(k,i) = cmeout(k,i)+cmei(k,i)

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! droplet activation
         ! calculate potential for droplet activation if cloud water is present
         ! formulation from Abdul-Razzak and Ghan (2000) and Abdul-Razzak et al. (1998), AR98
         ! number tendency (npccnin) is read in from companion routine

         ! assume aerosols already activated are equal to number of existing droplets for simplicity
         ! multiply by cloud fraction to obtain grid-average tendency

         if (qcic(k,i).ge.qsmall) then   
            npccn(k) = max(0._r8,npccnin(k,i))  
            dum2l(k,i)=(nc(k,i)+npccn(k)*deltat)/lcldm(k,i)
            dum2l(k,i)=max(dum2l(k,i),cdnl/rho(k,i)) ! sghan minimum in #/cm3  
            ncmax = dum2l(k,i)*lcldm(k,i)
         else
            npccn(k)=0._r8
            dum2l(k,i)=0._r8
            ncmax = 0._r8
         end if

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! get size distribution parameters based on in-cloud cloud water/ice 
         ! these calculations also ensure consistency between number and mixing ratio
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         !......................................................................
         ! cloud ice

         if (qiic(k,i).ge.qsmall) then

            ! add upper limit to in-cloud number concentration to prevent numerical error
            niic(k,i)=min(niic(k,i),qiic(k,i)*1.e20_r8)

            lami(k) = (cons1*ci*niic(k,i)/qiic(k,i))**(1._r8/di)
            n0i(k) = niic(k,i)*lami(k)

            ! check for slope
            ! adjust vars

            if (lami(k).lt.lammini) then

               lami(k) = lammini
               n0i(k) = lami(k)**(di+1._r8)*qiic(k,i)/(ci*cons1)
               niic(k,i) = n0i(k)/lami(k)
            else if (lami(k).gt.lammaxi) then
               lami(k) = lammaxi
               n0i(k) = lami(k)**(di+1._r8)*qiic(k,i)/(ci*cons1)
               niic(k,i) = n0i(k)/lami(k)
            end if

         else
            lami(k) = 0._r8
            n0i(k) = 0._r8
         end if

         if (qcic(k,i).ge.qsmall) then

            ! add upper limit to in-cloud number concentration to prevent numerical error
            ncic(k,i)=min(ncic(k,i),qcic(k,i)*1.e20_r8)
            ncic(k,i)=max(ncic(k,i),cdnl/rho(k,i)) ! sghan minimum in #/cm  

            ! get pgam from fit to observations of martin et al. 1994

            pgam(k)=0.0005714_r8*(ncic(k,i)/1.e6_r8*rho(k,i))+0.2714_r8
            pgam(k)=1._r8/(pgam(k)**2)-1._r8
            pgam(k)=max(pgam(k),2._r8)
            pgam(k)=min(pgam(k),15._r8)

            ! calculate lamc

            lamc(k) = (pi/6._r8*rhow*ncic(k,i)*gamma(pgam(k)+4._r8)/ &
                 (qcic(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)

            ! lammin, 50 micron diameter max mean size

            lammin = (pgam(k)+1._r8)/50.e-6_r8
            lammax = (pgam(k)+1._r8)/2.e-6_r8

            if (lamc(k).lt.lammin) then
               lamc(k) = lammin
               ncic(k,i) = 6._r8*lamc(k)**3*qcic(k,i)* &
                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
            else if (lamc(k).gt.lammax) then
               lamc(k) = lammax
               ncic(k,i) = 6._r8*lamc(k)**3*qcic(k,i)* &
                    gamma(pgam(k)+1._r8)/(pi*rhow*gamma(pgam(k)+4._r8))
            end if

            ! parameter to calculate droplet freezing

            cdist1(k) = ncic(k,i)/gamma(pgam(k)+1._r8) 

         else
            lamc(k) = 0._r8
            cdist1(k) = 0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! begin micropysical process calculations 
         !.................................................................
         ! autoconversion of cloud liquid water to rain
         ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
         ! minimum qc of 1 x 10^-8 prevents floating point error

         if (qcic(k,i).ge.1.e-8_r8) then

            ! nprc is increase in rain number conc due to autoconversion
            ! nprc1 is decrease in cloud droplet conc due to autoconversion

            ! assume exponential sub-grid distribution of qc, resulting in additional
            ! factor related to qcvar below

            ! hm switch for sub-columns, don't include sub-grid qc
            if (microp_uniform) then

               prc(k) = 1350._r8*qcic(k,i)**2.47_r8* &
                    (ncic(k,i)/1.e6_r8*rho(k,i))**(-1.79_r8)
               nprc(k) = prc(k)/(4._r8/3._r8*pi*rhow*(25.e-6_r8)**3)
               nprc1(k) = prc(k)/(qcic(k,i)/ncic(k,i))

            else

               prc(k) = cons2/(cons3*cons18)*1350._r8*qcic(k,i)**2.47_r8* &
                    (ncic(k,i)/1.e6_r8*rho(k,i))**(-1.79_r8)
               nprc(k) = prc(k)/cons22
               nprc1(k) = prc(k)/(qcic(k,i)/ncic(k,i))

            end if               ! sub-column switch

         else
            prc(k)=0._r8
            nprc(k)=0._r8
            nprc1(k)=0._r8
         end if

         ! add autoconversion to precip from above to get provisional rain mixing ratio
         ! and number concentration (qric and nric)

         ! 0.45 m/s is fallspeed of new rain drop (80 micron diameter)

         dum=0.45_r8
         dum1=0.45_r8

         if (k.eq.top_lev) then
            qric(k,i)=prc(k)*lcldm(k,i)*dz(k,i)/cldmax(k,i)/dum
            nric(k,i)=nprc(k)*lcldm(k,i)*dz(k,i)/cldmax(k,i)/dum
         else
            if (qric(k-1,i).ge.qsmall) then
               dum=umr(k-1)
               dum1=unr(k-1)
            end if

            ! no autoconversion of rain number if rain/snow falling from above
            ! this assumes that new drizzle drops formed by autoconversion are rapidly collected
            ! by the existing rain/snow particles from above

            if (qric(k-1,i).ge.1.e-9_r8.or.qniic(k-1,i).ge.1.e-9_r8) then
               nprc(k)=0._r8
            end if

            qric(k,i) = (rho(k-1,i)*umr(k-1)*qric(k-1,i)*cldmax(k-1,i)+ &
                 (rho(k,i)*dz(k,i)*((pra(k-1)+prc(k))*lcldm(k,i)+(pre(k-1)-pracs(k-1)-mnuccr(k-1))*cldmax(k,i))))&
                 /(dum*rho(k,i)*cldmax(k,i))
            nric(k,i) = (rho(k-1,i)*unr(k-1)*nric(k-1,i)*cldmax(k-1,i)+ &
                 (rho(k,i)*dz(k,i)*(nprc(k)*lcldm(k,i)+(nsubr(k-1)-npracs(k-1)-nnuccr(k-1)+nragg(k-1))*cldmax(k,i))))&
                 /(dum1*rho(k,i)*cldmax(k,i))

         end if

         !.......................................................................
         ! Autoconversion of cloud ice to snow
         ! similar to Ferrier (1994)

         if (do_cldice) then
            if (t(k,i).le.273.15_r8.and.qiic(k,i).ge.qsmall) then

               ! note: assumes autoconversion timescale of 180 sec
               
               nprci(k) = n0i(k)/(lami(k)*180._r8)*exp(-lami(k)*dcs)

               prci(k) = pi*rhoi*n0i(k)/(6._r8*180._r8)* &
                    (cons23/lami(k)+3._r8*cons24/lami(k)**2+ &
                    6._r8*dcs/lami(k)**3+6._r8/lami(k)**4)*exp(-lami(k)*dcs)
            else
               prci(k)=0._r8
               nprci(k)=0._r8
            end if
         else
            ! Add in the particles that we have already converted to snow, and
            ! don't do any further autoconversion of ice.
            prci(k)  = tnd_qsnow(k,i) / cldm(k,i)
            nprci(k) = tnd_nsnow(k,i) / cldm(k,i)
         end if

         ! add autoconversion to flux from level above to get provisional snow mixing ratio
         ! and number concentration (qniic and nsic)

         dum=(asn(k,i)*cons25)
         dum1=(asn(k,i)*cons25)

         if (k.eq.top_lev) then
            qniic(k,i)=prci(k)*icldm(k,i)*dz(k,i)/cldmax(k,i)/dum
            nsic(k,i)=nprci(k)*icldm(k,i)*dz(k,i)/cldmax(k,i)/dum
         else
            if (qniic(k-1,i).ge.qsmall) then
               dum=ums(k-1)
               dum1=uns(k-1)
            end if

            qniic(k,i) = (rho(k-1,i)*ums(k-1)*qniic(k-1,i)*cldmax(k-1,i)+ &
                 (rho(k,i)*dz(k,i)*((prci(k)+prai(k-1)+psacws(k-1)+bergs(k-1))*icldm(k,i)+(prds(k-1)+ &
                 pracs(k-1)+mnuccr(k-1))*cldmax(k,i))))&
                 /(dum*rho(k,i)*cldmax(k,i))

            nsic(k,i) = (rho(k-1,i)*uns(k-1)*nsic(k-1,i)*cldmax(k-1,i)+ &
                 (rho(k,i)*dz(k,i)*(nprci(k)*icldm(k,i)+(nsubs(k-1)+nsagg(k-1)+nnuccr(k-1))*cldmax(k,i))))&
                 /(dum1*rho(k,i)*cldmax(k,i))

         end if

         ! if precip mix ratio is zero so should number concentration

         if (qniic(k,i).lt.qsmall) then
            qniic(k,i)=0._r8
            nsic(k,i)=0._r8
         end if

         if (qric(k,i).lt.qsmall) then
            qric(k,i)=0._r8
            nric(k,i)=0._r8
         end if

         ! make sure number concentration is a positive number to avoid 
         ! taking root of negative later

         nric(k,i)=max(nric(k,i),0._r8)
         nsic(k,i)=max(nsic(k,i),0._r8)

         !.......................................................................
         ! get size distribution parameters for precip
         !......................................................................
         ! rain

         if (qric(k,i).ge.qsmall) then
            lamr(k) = (pi*rhow*nric(k,i)/qric(k,i))**(1._r8/3._r8)
            n0r(k) = nric(k,i)*lamr(k)

            ! check for slope
            ! adjust vars

            if (lamr(k).lt.lamminr) then

               lamr(k) = lamminr

               n0r(k) = lamr(k)**4*qric(k,i)/(pi*rhow)
               nric(k,i) = n0r(k)/lamr(k)
            else if (lamr(k).gt.lammaxr) then
               lamr(k) = lammaxr
               n0r(k) = lamr(k)**4*qric(k,i)/(pi*rhow)
               nric(k,i) = n0r(k)/lamr(k)
            end if

            ! provisional rain number and mass weighted mean fallspeed (m/s)

            unr(k) = min(arn(k,i)*cons4/lamr(k)**br,9.1_r8*rhof(k,i))
            umr(k) = min(arn(k,i)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(k,i))

         else
            lamr(k) = 0._r8
            n0r(k) = 0._r8
            umr(k) = 0._r8
            unr(k) = 0._r8
         end if

         !......................................................................
         ! snow

         if (qniic(k,i).ge.qsmall) then
            lams(k) = (cons6*cs*nsic(k,i)/qniic(k,i))**(1._r8/ds)
            n0s(k) = nsic(k,i)*lams(k)

            ! check for slope
            ! adjust vars

            if (lams(k).lt.lammins) then
               lams(k) = lammins
               n0s(k) = lams(k)**(ds+1._r8)*qniic(k,i)/(cs*cons6)
               nsic(k,i) = n0s(k)/lams(k)

            else if (lams(k).gt.lammaxs) then
               lams(k) = lammaxs
               n0s(k) = lams(k)**(ds+1._r8)*qniic(k,i)/(cs*cons6)
               nsic(k,i) = n0s(k)/lams(k)
            end if

            ! provisional snow number and mass weighted mean fallspeed (m/s)

            ums(k) = min(asn(k,i)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(k,i))
            uns(k) = min(asn(k,i)*cons7/lams(k)**bs,1.2_r8*rhof(k,i))

         else
            lams(k) = 0._r8
            n0s(k) = 0._r8
            ums(k) = 0._r8
            uns(k) = 0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         ! heterogeneous freezing of cloud water

         if (do_cldice .and. qcic(k,i).ge.qsmall .and. t(k,i).lt.269.15_r8) then
 
            ! immersion freezing (Bigg, 1953)


            ! subcolumns

            if (microp_uniform) then

               mnuccc(k) = &
                    pi*pi/36._r8*rhow* &
                    cdist1(k)*gamma(7._r8+pgam(k))* &
                    bimm*(exp(aimm*(273.15_r8-t(k,i)))-1._r8)/ &
                    lamc(k)**3/lamc(k)**3

               nnuccc(k) = &
                    pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                    *bimm* &
                    (exp(aimm*(273.15_r8-t(k,i)))-1._r8)/lamc(k)**3

            else

               mnuccc(k) = cons9/(cons3*cons19)* &
                    pi*pi/36._r8*rhow* &
                    cdist1(k)*gamma(7._r8+pgam(k))* &
                    bimm*(exp(aimm*(273.15_r8-t(k,i)))-1._r8)/ &
                    lamc(k)**3/lamc(k)**3

               nnuccc(k) = cons10/(cons3*qcvar)* &
                    pi/6._r8*cdist1(k)*gamma(pgam(k)+4._r8) &
                    *bimm* &
                    (exp(aimm*(273.15_r8-t(k,i)))-1._r8)/lamc(k)**3
            end if           ! sub-columns

            ! contact freezing (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
            ! dust size and number in 4 bins are read in from companion routine

            tcnt=(270.16_r8-t(k,i))**1.3_r8
            viscosity=1.8e-5_r8*(t(k,i)/298.0_r8)**0.85_r8    ! Viscosity (kg/m/s)
            mfp=2.0_r8*viscosity/(p(k,i)  &                   ! Mean free path (m)
                 *sqrt(8.0_r8*28.96e-3_r8/(pi*8.314409_r8*t(k,i))))           

            nslip1=1.0_r8+(mfp/rndst(1,k,i))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(1,k,i)/mfp))))! Slip correction factor
            nslip2=1.0_r8+(mfp/rndst(2,k,i))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(2,k,i)/mfp))))
            nslip3=1.0_r8+(mfp/rndst(3,k,i))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(3,k,i)/mfp))))
            nslip4=1.0_r8+(mfp/rndst(4,k,i))*(1.257_r8+(0.4_r8*Exp(-(1.1_r8*rndst(4,k,i)/mfp))))

            ndfaer1=1.381e-23_r8*t(k,i)*nslip1/(6._r8*pi*viscosity*rndst(1,k,i))  ! aerosol diffusivity (m2/s)
            ndfaer2=1.381e-23_r8*t(k,i)*nslip2/(6._r8*pi*viscosity*rndst(2,k,i))
            ndfaer3=1.381e-23_r8*t(k,i)*nslip3/(6._r8*pi*viscosity*rndst(3,k,i))
            ndfaer4=1.381e-23_r8*t(k,i)*nslip4/(6._r8*pi*viscosity*rndst(4,k,i))

            if (microp_uniform) then

               mnucct(k) = &
                    (ndfaer1*(nacon(1,k,i)*tcnt)+ndfaer2*(nacon(2,k,i)*tcnt)+ &
                    ndfaer3*(nacon(3,k,i)*tcnt)+ndfaer4*(nacon(4,k,i)*tcnt))*pi*pi/3._r8*rhow* &
                    cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4

               nnucct(k) = (ndfaer1*(nacon(1,k,i)*tcnt)+ndfaer2*(nacon(2,k,i)*tcnt)+ &
                    ndfaer3*(nacon(3,k,i)*tcnt)+ndfaer4*(nacon(4,k,i)*tcnt))*2._r8*pi*  &
                    cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)

            else

               mnucct(k) = gamma(qcvar+4._r8/3._r8)/(cons3*qcvar**(4._r8/3._r8))*  &
                    (ndfaer1*(nacon(1,k,i)*tcnt)+ndfaer2*(nacon(2,k,i)*tcnt)+ &
                    ndfaer3*(nacon(3,k,i)*tcnt)+ndfaer4*(nacon(4,k,i)*tcnt))*pi*pi/3._r8*rhow* &
                    cdist1(k)*gamma(pgam(k)+5._r8)/lamc(k)**4

               nnucct(k) =  gamma(qcvar+1._r8/3._r8)/(cons3*qcvar**(1._r8/3._r8))*  &
                    (ndfaer1*(nacon(1,k,i)*tcnt)+ndfaer2*(nacon(2,k,i)*tcnt)+ &
                    ndfaer3*(nacon(3,k,i)*tcnt)+ndfaer4*(nacon(4,k,i)*tcnt))*2._r8*pi*  &
                    cdist1(k)*gamma(pgam(k)+2._r8)/lamc(k)

            end if      ! sub-column switch
 
            ! make sure number of droplets frozen does not exceed available ice nuclei concentration
            ! this prevents 'runaway' droplet freezing

            if (nnuccc(k)*lcldm(k,i).gt.nnuccd(k)) then
               dum=(nnuccd(k)/(nnuccc(k)*lcldm(k,i)))
               ! scale mixing ratio of droplet freezing with limit
               mnuccc(k)=mnuccc(k)*dum
               nnuccc(k)=nnuccd(k)/lcldm(k,i)
            end if

         else
            mnuccc(k)=0._r8
            nnuccc(k)=0._r8
            mnucct(k)=0._r8
            nnucct(k)=0._r8
         end if

         !.......................................................................
         ! snow self-aggregation from passarelli, 1978, used by reisner, 1998
         ! this is hard-wired for bs = 0.4 for now
         ! ignore self-collection of cloud ice

         if (qniic(k,i).ge.qsmall .and. t(k,i).le.273.15_r8) then
            nsagg(k) = -1108._r8*asn(k,i)*Eii* &
                 pi**((1._r8-bs)/3._r8)*rhosn**((-2._r8-bs)/3._r8)*rho(k,i)** &
                 ((2._r8+bs)/3._r8)*qniic(k,i)**((2._r8+bs)/3._r8)* &
                 (nsic(k,i)*rho(k,i))**((4._r8-bs)/3._r8)/ &
                 (4._r8*720._r8*rho(k,i))
         else
            nsagg(k)=0._r8
         end if

         !.......................................................................
         ! accretion of cloud droplets onto snow/graupel
         ! here use continuous collection equation with
         ! simple gravitational collection kernel
         ! ignore collisions between droplets/cloud ice
         ! since minimum size ice particle for accretion is 50 - 150 micron

         ! ignore collision of snow with droplets above freezing

         if (qniic(k,i).ge.qsmall .and. t(k,i).le.tmelt .and. &
              qcic(k,i).ge.qsmall) then

            ! put in size dependent collection efficiency
            ! mean diameter of snow is area-weighted, since
            ! accretion is function of crystal geometric area
            ! collection efficiency is approximation based on stoke's law (Thompson et al. 2004)

            dc0 = (pgam(k)+1._r8)/lamc(k)
            ds0 = 1._r8/lams(k)
            dum = dc0*dc0*uns(k)*rhow/(9._r8*mu(k,i)*ds0)
            eci = dum*dum/((dum+0.4_r8)*(dum+0.4_r8))

            eci = max(eci,0._r8)
            eci = min(eci,1._r8)


            ! no impact of sub-grid distribution of qc since psacws
            ! is linear in qc

            psacws(k) = pi/4._r8*asn(k,i)*qcic(k,i)*rho(k,i)* &
                 n0s(k)*Eci*cons11/ &
                 lams(k)**(bs+3._r8)
            npsacws(k) = pi/4._r8*asn(k,i)*ncic(k,i)*rho(k,i)* &
                 n0s(k)*Eci*cons11/ &
                 lams(k)**(bs+3._r8)
         else
            psacws(k)=0._r8
            npsacws(k)=0._r8
         end if

         ! add secondary ice production due to accretion of droplets by snow 
         ! (Hallet-Mossop process) (from Cotton et al., 1986)

         if (.not. do_cldice) then
            ni_secp   = 0.0_r8
            nsacwi(k) = 0.0_r8
            msacwi(k) = 0.0_r8
         else if((t(k,i).lt.270.16_r8) .and. (t(k,i).ge.268.16_r8)) then
            ni_secp   = 3.5e8_r8*(270.16_r8-t(k,i))/2.0_r8*psacws(k)
            nsacwi(k) = ni_secp
            msacwi(k) = min(ni_secp*mi0,psacws(k))
         else if((t(k,i).lt.268.16_r8) .and. (t(k,i).ge.265.16_r8)) then
            ni_secp   = 3.5e8_r8*(t(k,i)-265.16_r8)/3.0_r8*psacws(k)
            nsacwi(k) = ni_secp
            msacwi(k) = min(ni_secp*mi0,psacws(k))
         else
            ni_secp   = 0.0_r8
            nsacwi(k) = 0.0_r8
            msacwi(k) = 0.0_r8
         endif
         psacws(k) = max(0.0_r8,psacws(k)-ni_secp*mi0)

         !.......................................................................
         ! accretion of rain water by snow
         ! formula from ikawa and saito, 1991, used by reisner et al., 1998

         if (qric(k,i).ge.1.e-8_r8 .and. qniic(k,i).ge.1.e-8_r8 .and. & 
              t(k,i).le.273.15_r8) then

            pracs(k) = pi*pi*ecr*(((1.2_r8*umr(k)-0.95_r8*ums(k))**2+ &
                 0.08_r8*ums(k)*umr(k))**0.5_r8*rhow*rho(k,i)* &
                 n0r(k)*n0s(k)* &
                 (5._r8/(lamr(k)**6*lams(k))+ &
                 2._r8/(lamr(k)**5*lams(k)**2)+ &
                 0.5_r8/(lamr(k)**4*lams(k)**3)))

            npracs(k) = pi/2._r8*rho(k,i)*ecr*(1.7_r8*(unr(k)-uns(k))**2+ &
                 0.3_r8*unr(k)*uns(k))**0.5_r8*n0r(k)*n0s(k)* &
                 (1._r8/(lamr(k)**3*lams(k))+ &
                 1._r8/(lamr(k)**2*lams(k)**2)+ &
                 1._r8/(lamr(k)*lams(k)**3))

         else
            pracs(k)=0._r8
            npracs(k)=0._r8
         end if

         !.......................................................................
         ! heterogeneous freezing of rain drops
         ! follows from Bigg (1953)

         if (t(k,i).lt.269.15_r8 .and. qric(k,i).ge.qsmall) then

            mnuccr(k) = 20._r8*pi*pi*rhow*nric(k,i)*bimm* &
                 (exp(aimm*(273.15_r8-t(k,i)))-1._r8)/lamr(k)**3 &
                 /lamr(k)**3

            nnuccr(k) = pi*nric(k,i)*bimm* &
                 (exp(aimm*(273.15_r8-t(k,i)))-1._r8)/lamr(k)**3
         else
            mnuccr(k)=0._r8
            nnuccr(k)=0._r8
         end if

         !.......................................................................
         ! accretion of cloud liquid water by rain
         ! formula from Khrouditnov and Kogan (2000)
         ! gravitational collection kernel, droplet fall speed neglected

         if (qric(k,i).ge.qsmall .and. qcic(k,i).ge.qsmall) then

            ! include sub-grid distribution of cloud water

            ! add sub-column switch

            if (microp_uniform) then

               pra(k) = 67._r8*(qcic(k,i)*qric(k,i))**1.15_r8
               npra(k) = pra(k)/(qcic(k,i)/ncic(k,i))

            else

               pra(k) = accre_enhan(k,i)*(cons12/(cons3*cons20)*67._r8*(qcic(k,i)*qric(k,i))**1.15_r8)
               npra(k) = pra(k)/(qcic(k,i)/ncic(k,i))

            end if               ! sub-column switch

         else
            pra(k)=0._r8
            npra(k)=0._r8
         end if

         !.......................................................................
         ! Self-collection of rain drops
         ! from Beheng(1994)

         if (qric(k,i).ge.qsmall) then
            nragg(k) = -8._r8*nric(k,i)*qric(k,i)*rho(k,i)
         else
            nragg(k)=0._r8
         end if

         !.......................................................................
         ! Accretion of cloud ice by snow
         ! For this calculation, it is assumed that the Vs >> Vi
         ! and Ds >> Di for continuous collection

         if (do_cldice .and. qniic(k,i).ge.qsmall.and.qiic(k,i).ge.qsmall &
              .and.t(k,i).le.273.15_r8) then

            prai(k) = pi/4._r8*asn(k,i)*qiic(k,i)*rho(k,i)* &
                 n0s(k)*Eii*cons11/ &
                 lams(k)**(bs+3._r8)
            nprai(k) = pi/4._r8*asn(k,i)*niic(k,i)* &
                 rho(k,i)*n0s(k)*Eii*cons11/ &
                 lams(k)**(bs+3._r8)
         else
            prai(k)=0._r8
            nprai(k)=0._r8
         end if

         !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! calculate evaporation/sublimation of rain and snow
         ! note: evaporation/sublimation occurs only in cloud-free portion of grid cell
         ! in-cloud condensation/deposition of rain and snow is neglected
         ! except for transfer of cloud water to snow through bergeron process

         ! initialize evap/sub tendncies
         pre(k)=0._r8
         prds(k)=0._r8

         ! evaporation of rain
         ! only calculate if there is some precip fraction > cloud fraction

         if (qcic(k,i)+qiic(k,i).lt.1.e-6_r8.or.cldmax(k,i).gt.lcldm(k,i)) then

            ! set temporary cloud fraction to zero if cloud water + ice is very small
            ! this will ensure that evaporation/sublimation of precip occurs over
            ! entire grid cell, since min cloud fraction is specified otherwise
            if (qcic(k,i)+qiic(k,i).lt.1.e-6_r8) then
               dum=0._r8
            else
               dum=lcldm(k,i)
            end if

            ! saturation vapor pressure
            esn=svp_water(t(k,i))
            qsn=svp_to_qsat(esn, p(k,i))

            ! recalculate saturation vapor pressure for liquid and ice
            esl(k,i)=esn
            esi(k,i)=svp_ice(t(k,i))
            ! hm fix, make sure when above freezing that esi=esl, not active yet
            if (t(k,i).gt.tmelt)esi(k,i)=esl(k,i)

            ! calculate q for out-of-cloud region
            qclr=(q(k,i)-dum*qsn)/(1._r8-dum)

            if (qric(k,i).ge.qsmall) then

               qvs=svp_to_qsat(esl(k,i), p(k,i))
               dqsdt = xxlv*qvs/(rv*t(k,i)**2)
               ab = 1._r8+dqsdt*xxlv/cpp
               epsr = 2._r8*pi*n0r(k)*rho(k,i)*Dv(k,i)* &
                    (f1r/(lamr(k)*lamr(k))+ &
                    f2r*(arn(k,i)*rho(k,i)/mu(k,i))**0.5_r8* &
                    sc(k,i)**(1._r8/3._r8)*cons13/ &
                    (lamr(k)**(5._r8/2._r8+br/2._r8)))

               pre(k) = epsr*(qclr-qvs)/ab

               ! only evaporate in out-of-cloud region
               ! and distribute across cldmax
               pre(k)=min(pre(k)*(cldmax(k,i)-dum),0._r8)
               pre(k)=pre(k)/cldmax(k,i)
            end if

            ! sublimation of snow
            if (qniic(k,i).ge.qsmall) then
               qvi=svp_to_qsat(esi(k,i), p(k,i))
               dqsidt =  xxls*qvi/(rv*t(k,i)**2)
               abi = 1._r8+dqsidt*xxls/cpp
               epss = 2._r8*pi*n0s(k)*rho(k,i)*Dv(k,i)* &
                    (f1s/(lams(k)*lams(k))+ &
                    f2s*(asn(k,i)*rho(k,i)/mu(k,i))**0.5_r8* &
                    sc(k,i)**(1._r8/3._r8)*cons14/ &
                    (lams(k)**(5._r8/2._r8+bs/2._r8)))
               prds(k) = epss*(qclr-qvi)/abi	

               ! only sublimate in out-of-cloud region and distribute over cldmax
               prds(k)=min(prds(k)*(cldmax(k,i)-dum),0._r8)
               prds(k)=prds(k)/cldmax(k,i)
            end if

            ! make sure RH not pushed above 100% due to rain evaporation/snow sublimation
            ! get updated RH at end of time step based on cloud water/ice condensation/evap

            qtmp=q(k,i)-(cmei(k,i)+(pre(k)+prds(k))*cldmax(k,i))*deltat
            ttmp=t(k,i)+((pre(k)*cldmax(k,i))*xxlv+ &
                 (cmei(k,i)+prds(k)*cldmax(k,i))*xxls)*deltat/cpp

            !limit range of temperatures!
            ttmp=max(180._r8,min(ttmp,323._r8))

            esn=svp_water(ttmp)  ! use rhw to allow ice supersaturation
            qsn=svp_to_qsat(esn, p(k,i))

            ! modify precip evaporation rate if q > qsat
            if (qtmp.gt.qsn) then
               if (pre(k)+prds(k).lt.-1.e-20_r8) then
                  dum1=pre(k)/(pre(k)+prds(k))
                  ! recalculate q and t after cloud water cond but without precip evap
                  qtmp=q(k,i)-(cmei(k,i))*deltat
                  ttmp=t(k,i)+(cmei(k,i)*xxls)*deltat/cpp
                  esn=svp_water(ttmp) ! use rhw to allow ice supersaturation
                  qsn=svp_to_qsat(esn, p(k,i))
                  dum=(qtmp-qsn)/(1._r8 + cons27*qsn/(cpp*rv*ttmp**2))
                  dum=min(dum,0._r8)

                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
                  pre(k)=dum*dum1/deltat/cldmax(k,i)

                  ! do separately using RHI for prds....
                  esn=svp_ice(ttmp) ! use rhi to allow ice supersaturation
                  qsn=svp_to_qsat(esn, p(k,i))
                  dum=(qtmp-qsn)/(1._r8 + cons28*qsn/(cpp*rv*ttmp**2))
                  dum=min(dum,0._r8)

                  ! modify rates if needed, divide by cldmax to get local (in-precip) value
                  prds(k)=dum*(1._r8-dum1)/deltat/cldmax(k,i)
               end if
            end if
         end if

         ! bergeron process - evaporation of droplets and deposition onto snow

         if (qniic(k,i).ge.qsmall.and.qcic(k,i).ge.qsmall.and.t(k,i).lt.tmelt) then
            qvi=svp_to_qsat(esi(k,i), p(k,i))
            qvs=svp_to_qsat(esl(k,i), p(k,i))
            dqsidt =  xxls*qvi/(rv*t(k,i)**2)
            abi = 1._r8+dqsidt*xxls/cpp
            epss = 2._r8*pi*n0s(k)*rho(k,i)*Dv(k,i)* &
                 (f1s/(lams(k)*lams(k))+ &
                 f2s*(asn(k,i)*rho(k,i)/mu(k,i))**0.5_r8* &
                 sc(k,i)**(1._r8/3._r8)*cons14/ &
                 (lams(k)**(5._r8/2._r8+bs/2._r8)))
            bergs(k)=epss*(qvs-qvi)/abi
         else
            bergs(k)=0._r8
         end if

         !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         ! conservation to ensure no negative values of cloud water/precipitation
         ! in case microphysical process rates are large

         ! make sure and use end-of-time step values for cloud water, ice, due
         ! condensation/deposition

         ! note: for check on conservation, processes are multiplied by omsm
         ! to prevent problems due to round off error

         ! include mixing timescale  (mtime)

         qce=(qc(k,i) - berg(k,i)*deltat)
         nce=(nc(k,i)+npccn(k)*deltat*mtime)
         qie=(qi(k,i)+(cmei(k,i)+berg(k,i))*deltat)
         nie=(ni(k,i)+nnuccd(k)*deltat*mtime)

         ! conservation of qc

         dum = (prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+ &
              psacws(k)+bergs(k))*lcldm(k,i)*deltat

         if (dum.gt.qce) then
            ratio = qce/deltat/lcldm(k,i)/(prc(k)+pra(k)+mnuccc(k)+mnucct(k)+msacwi(k)+psacws(k)+bergs(k))*omsm 

            prc(k) = prc(k)*ratio
            pra(k) = pra(k)*ratio
            mnuccc(k) = mnuccc(k)*ratio
            mnucct(k) = mnucct(k)*ratio  
            msacwi(k) = msacwi(k)*ratio  
            psacws(k) = psacws(k)*ratio
            bergs(k) = bergs(k)*ratio
         end if
 
         ! conservation of nc
         dum = (nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+ &
              npsacws(k)-nsubc(k))*lcldm(k,i)*deltat

         if (dum.gt.nce) then
            ratio = nce/deltat/((nprc1(k)+npra(k)+nnuccc(k)+nnucct(k)+&
                 npsacws(k)-nsubc(k))*lcldm(k,i))*omsm

            nprc1(k) = nprc1(k)*ratio
            npra(k) = npra(k)*ratio
            nnuccc(k) = nnuccc(k)*ratio
            nnucct(k) = nnucct(k)*ratio  
            npsacws(k) = npsacws(k)*ratio
            nsubc(k)=nsubc(k)*ratio
         end if

         ! conservation of qi

         if (do_cldice) then
            dum = ((-mnuccc(k)-mnucct(k)-msacwi(k))*lcldm(k,i)+(prci(k)+ &
                 prai(k))*icldm(k,i))*deltat

            if (dum.gt.qie) then

               ratio = (qie/deltat+(mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(k,i))/((prci(k)+prai(k))*icldm(k,i))*omsm 
               prci(k) = prci(k)*ratio
               prai(k) = prai(k)*ratio
            end if

            ! conservation of ni

            dum = ((-nnucct(k)-nsacwi(k))*lcldm(k,i)+(nprci(k)+ &
                 nprai(k)-nsubi(k))*icldm(k,i))*deltat

            if (dum.gt.nie) then

               ratio = (nie/deltat+(nnucct(k)+nsacwi(k))*lcldm(k,i))/ &  
                    ((nprci(k)+nprai(k)-nsubi(k))*icldm(k,i))*omsm
               nprci(k) = nprci(k)*ratio
               nprai(k) = nprai(k)*ratio
               nsubi(k) = nsubi(k)*ratio
            end if
         end if

         ! for precipitation conservation, use logic that vertical integral 
         ! of tendency from current level to top of model (i.e., qrtot) cannot be negative

         ! conservation of rain mixing rat

         if (((prc(k)+pra(k))*lcldm(k,i)+(-mnuccr(k)+pre(k)-pracs(k))*&
              cldmax(k,i))*dz(k,i)*rho(k,i)+qrtot.lt.0._r8) then

            if (-pre(k)+pracs(k)+mnuccr(k).ge.qsmall) then

               ratio = (qrtot/(dz(k,i)*rho(k,i))+(prc(k)+pra(k))*lcldm(k,i))/&
                    ((-pre(k)+pracs(k)+mnuccr(k))*cldmax(k,i))*omsm 

               pre(k) = pre(k)*ratio
               pracs(k) = pracs(k)*ratio
               mnuccr(k) = mnuccr(k)*ratio
            end if
         end if
 
         ! conservation of nr
         ! for now neglect evaporation of nr
         nsubr(k)=0._r8

         if ((nprc(k)*lcldm(k,i)+(-nnuccr(k)+nsubr(k)-npracs(k)&
              +nragg(k))*cldmax(k,i))*dz(k,i)*rho(k,i)+nrtot.lt.0._r8) then

            if (-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k).ge.qsmall) then

               ratio = (nrtot/(dz(k,i)*rho(k,i))+nprc(k)*lcldm(k,i))/&
                    ((-nsubr(k)-nragg(k)+npracs(k)+nnuccr(k))*cldmax(k,i))*omsm

               nsubr(k) = nsubr(k)*ratio
               npracs(k) = npracs(k)*ratio
               nnuccr(k) = nnuccr(k)*ratio
               nragg(k) = nragg(k)*ratio
            end if
         end if

         ! conservation of snow mix ratio

         if (((bergs(k)+psacws(k))*lcldm(k,i)+(prai(k)+prci(k))*icldm(k,i)+(pracs(k)+&
              mnuccr(k)+prds(k))*cldmax(k,i))*dz(k,i)*rho(k,i)+qstot.lt.0._r8) then

            if (-prds(k).ge.qsmall) then

               ratio = (qstot/(dz(k,i)*rho(k,i))+(bergs(k)+psacws(k))*lcldm(k,i)+(prai(k)+prci(k))*icldm(k,i)+&
                    (pracs(k)+mnuccr(k))*cldmax(k,i))/(-prds(k)*cldmax(k,i))*omsm

               prds(k) = prds(k)*ratio
            end if
         end if

         ! conservation of ns

         ! calculate loss of number due to sublimation
         ! for now neglect sublimation of ns
         nsubs(k)=0._r8

         if ((nprci(k)*icldm(k,i)+(nnuccr(k)+nsubs(k)+nsagg(k))*cldmax(k,i))*&
              dz(k,i)*rho(k,i)+nstot.lt.0._r8) then

            if (-nsubs(k)-nsagg(k).ge.qsmall) then

               ratio = (nstot/(dz(k,i)*rho(k,i))+nprci(k)*icldm(k,i)+&
                    nnuccr(k)*cldmax(k,i))/((-nsubs(k)-nsagg(k))*cldmax(k,i))*omsm

               nsubs(k) = nsubs(k)*ratio
               nsagg(k) = nsagg(k)*ratio
            end if
         end if

         ! get tendencies due to microphysical conversion processes
         ! note: tendencies are multiplied by appropaiate cloud/precip 
         ! fraction to get grid-scale values
         ! note: cmei is already grid-average values

         qvlat(k,i) = qvlat(k,i)-(pre(k)+prds(k))*cldmax(k,i)-cmei(k,i) 

         tlat(k,i) = tlat(k,i)+                                                      &
                    ((pre(k)*cldmax(k,i))*xxlv+(prds(k)*cldmax(k,i)+cmei(k,i))*xxls+ &
                    ((bergs(k)+psacws(k)+mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(k,i)+  &
                     (mnuccr(k)+ pracs(k))*cldmax(k,i)+berg(k,i))*xlf)
 
         qctend(k,i) = qctend(k,i)+ &
              (-pra(k)-prc(k)-mnuccc(k)-mnucct(k)-msacwi(k)- & 
              psacws(k)-bergs(k))*lcldm(k,i)-berg(k,i)

         if (do_cldice) then
            qitend(k,i) = qitend(k,i)+ &
                 (mnuccc(k)+mnucct(k)+msacwi(k))*lcldm(k,i)+(-prci(k)- & 
                 prai(k))*icldm(k,i)+cmei(k,i)+berg(k,i)
         end if

         qrtend(k,i) = qrtend(k,i)+ &
              (pra(k)+prc(k))*lcldm(k,i)+(pre(k)-pracs(k)- &
              mnuccr(k))*cldmax(k,i)

         qnitend(k,i) = qnitend(k,i)+ &
              (prai(k)+prci(k))*icldm(k,i)+(psacws(k)+bergs(k))*lcldm(k,i)+(prds(k)+ &
              pracs(k)+mnuccr(k))*cldmax(k,i)

         ! add output for cmei (accumulate)
         cmeiout(k,i) = cmeiout(k,i) + cmei(k,i)

         ! assign variables for trop_mozart, these are grid-average
         ! evaporation/sublimation is stored here as positive term

         evapsnow(k,i) = evapsnow(k,i)-prds(k)*cldmax(k,i)
         nevapr(k,i) = nevapr(k,i)-pre(k)*cldmax(k,i)

         ! change to make sure prain is positive: do not remove snow from
         ! prain used for wet deposition
         prain(k,i) = prain(k,i)+(pra(k)+prc(k))*lcldm(k,i)+(-pracs(k)- &
              mnuccr(k))*cldmax(k,i)
         prodsnow(k,i) = prodsnow(k,i)+(prai(k)+prci(k))*icldm(k,i)+(psacws(k)+bergs(k))*lcldm(k,i)+(&
              pracs(k)+mnuccr(k))*cldmax(k,i)

         ! following are used to calculate 1st order conversion rate of cloud water
         !    to rain and snow (1/s), for later use in aerosol wet removal routine
         ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
         !    used to calculate pra, prc, ... in this routine
         ! qcsinksum_rate1ord = sum over iterations{ rate of direct transfer of cloud water to rain & snow }
         !                      (no cloud ice or bergeron terms)
         ! qcsum_rate1ord     = sum over iterations{ qc used in calculation of the transfer terms }

         qcsinksum_rate1ord(k) = qcsinksum_rate1ord(k) + (pra(k)+prc(k)+psacws(k))*lcldm(k,i) 
         qcsum_rate1ord(k) = qcsum_rate1ord(k) + qc(k,i) 

         ! microphysics output, note this is grid-averaged
         prao(k,i)=prao(k,i)+pra(k)*lcldm(k,i)
         prco(k,i)=prco(k,i)+prc(k)*lcldm(k,i)
         mnuccco(k,i)=mnuccco(k,i)+mnuccc(k)*lcldm(k,i)
         mnuccto(k,i)=mnuccto(k,i)+mnucct(k)*lcldm(k,i)
         mnuccdo(k,i)=mnuccdo(k,i)+mnuccd(k)*lcldm(k,i)
         msacwio(k,i)=msacwio(k,i)+msacwi(k)*lcldm(k,i)
         psacwso(k,i)=psacwso(k,i)+psacws(k)*lcldm(k,i)
         bergso(k,i)=bergso(k,i)+bergs(k)*lcldm(k,i)
         bergo(k,i)=bergo(k,i)+berg(k,i)
         prcio(k,i)=prcio(k,i)+prci(k)*icldm(k,i)
         praio(k,i)=praio(k,i)+prai(k)*icldm(k,i)
         mnuccro(k,i)=mnuccro(k,i)+mnuccr(k)*cldmax(k,i)
         pracso (k,i)=pracso (k,i)+pracs (k)*cldmax(k,i)

         ! multiply activation/nucleation by mtime to account for fast timescale

         nctend(k,i) = nctend(k,i)+ npccn(k)*mtime+&
              (-nnuccc(k)-nnucct(k)-npsacws(k)+nsubc(k) & 
              -npra(k)-nprc1(k))*lcldm(k,i)      

         if (do_cldice) then
            nitend(k,i) = nitend(k,i)+ nnuccd(k)*mtime+ & 
                 (nnucct(k)+nsacwi(k))*lcldm(k,i)+(nsubi(k)-nprci(k)- &   
                 nprai(k))*icldm(k,i)
         end if

         nstend(k,i) = nstend(k,i)+(nsubs(k)+ &
              nsagg(k)+nnuccr(k))*cldmax(k,i)+nprci(k)*icldm(k,i)

         nrtend(k,i) = nrtend(k,i)+ &
              nprc(k)*lcldm(k,i)+(nsubr(k)-npracs(k)-nnuccr(k) &
              +nragg(k))*cldmax(k,i)

         ! make sure that nc and ni at advanced time step do not exceed
         ! maximum (existing N + source terms*dt), which is possible due to
         ! fast nucleation timescale

         if (nctend(k,i).gt.0._r8.and.nc(k,i)+nctend(k,i)*deltat.gt.ncmax) then
            nctend(k,i)=max(0._r8,(ncmax-nc(k,i))/deltat)
         end if

         if (do_cldice .and. nitend(k,i).gt.0._r8.and.ni(k,i)+nitend(k,i)*deltat.gt.nimax) then
            nitend(k,i)=max(0._r8,(nimax-ni(k,i))/deltat)
         end if

         ! get final values for precipitation q and N, based on
         ! flux of precip from above, source/sink term, and terminal fallspeed
         ! see eq. 15-16 in MG2008

         ! rain

         if (qric(k,i).ge.qsmall) then
            if (k.eq.top_lev) then
               qric(k,i)=qrtend(k,i)*dz(k,i)/cldmax(k,i)/umr(k)
               nric(k,i)=nrtend(k,i)*dz(k,i)/cldmax(k,i)/unr(k)
            else
               qric(k,i) = (rho(k-1,i)*umr(k-1)*qric(k-1,i)*cldmax(k-1,i)+ &
                    (rho(k,i)*dz(k,i)*qrtend(k,i)))/(umr(k)*rho(k,i)*cldmax(k,i))
               nric(k,i) = (rho(k-1,i)*unr(k-1)*nric(k-1,i)*cldmax(k-1,i)+ &
                    (rho(k,i)*dz(k,i)*nrtend(k,i)))/(unr(k)*rho(k,i)*cldmax(k,i))

            end if
         else
            qric(k,i)=0._r8
            nric(k,i)=0._r8
         end if

         ! snow

         if (qniic(k,i).ge.qsmall) then
            if (k.eq.top_lev) then
               qniic(k,i)=qnitend(k,i)*dz(k,i)/cldmax(k,i)/ums(k)
               nsic(k,i)=nstend(k,i)*dz(k,i)/cldmax(k,i)/uns(k)
            else
               qniic(k,i) = (rho(k-1,i)*ums(k-1)*qniic(k-1,i)*cldmax(k-1,i)+ &
                    (rho(k,i)*dz(k,i)*qnitend(k,i)))/(ums(k)*rho(k,i)*cldmax(k,i))
               nsic(k,i) = (rho(k-1,i)*uns(k-1)*nsic(k-1,i)*cldmax(k-1,i)+ &
                    (rho(k,i)*dz(k,i)*nstend(k,i)))/(uns(k)*rho(k,i)*cldmax(k,i))
            end if
         else
            qniic(k,i)=0._r8
            nsic(k,i)=0._r8
         end if

         ! calculate precipitation flux at surface
         ! divide by density of water to get units of m/s

         prect(i) = prect(i)+(qrtend(k,i)*dz(k,i)*rho(k,i)+&
                    qnitend(k,i)*dz(k,i)*rho(k,i))/rhow
         preci(i) = preci(i)+qnitend(k,i)*dz(k,i)*rho(k,i)/rhow

         ! convert rain rate from m/s to mm/hr

         rainrt(k,i)=qric(k,i)*rho(k,i)*umr(k)/rhow*3600._r8*1000._r8

         ! vertically-integrated precip source/sink terms (note: grid-averaged)

         qrtot = max(qrtot+qrtend(k,i)*dz(k,i)*rho(k,i),0._r8)
         qstot = max(qstot+qnitend(k,i)*dz(k,i)*rho(k,i),0._r8)
         nrtot = max(nrtot+nrtend(k,i)*dz(k,i)*rho(k,i),0._r8)
         nstot = max(nstot+nstend(k,i)*dz(k,i)*rho(k,i),0._r8)

         ! calculate melting and freezing of precip

         ! melt snow at +2 C

         if (t(k,i)+tlat(k,i)/cpp*deltat > 275.15_r8) then
            if (qstot > 0._r8) then

               ! make sure melting snow doesn't reduce temperature below threshold
               dum = -xlf/cpp*qstot/(dz(k,i)*rho(k,i))
               if (t(k,i)+tlat(k,i)/cpp*deltat+dum.lt.275.15_r8) then
                  dum = (t(k,i)+tlat(k,i)/cpp*deltat-275.15_r8)*cpp/xlf
                  dum = dum/(xlf/cpp*qstot/(dz(k,i)*rho(k,i)))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qric(k,i)=qric(k,i)+dum*qniic(k,i)
               nric(k,i)=nric(k,i)+dum*nsic(k,i)
               qniic(k,i)=(1._r8-dum)*qniic(k,i)
               nsic(k,i)=(1._r8-dum)*nsic(k,i)
               ! heating tendency 
               tmp=-xlf*dum*qstot/(dz(k,i)*rho(k,i))
               meltsdt(k,i)=meltsdt(k,i) + tmp

               tlat(k,i)=tlat(k,i)+tmp
               qrtot=qrtot+dum*qstot
               nrtot=nrtot+dum*nstot
               qstot=(1._r8-dum)*qstot
               nstot=(1._r8-dum)*nstot
               preci(i)=(1._r8-dum)*preci(i)
            end if
         end if
 
         ! freeze all rain at -5C for Arctic

         if (t(k,i)+tlat(k,i)/cpp*deltat < (tmelt - 5._r8)) then

            if (qrtot > 0._r8) then

               ! make sure freezing rain doesn't increase temperature above threshold
               dum = xlf/cpp*qrtot/(dz(k,i)*rho(k,i))
               if (t(k,i)+tlat(k,i)/cpp*deltat+dum.gt.(tmelt - 5._r8)) then
                  dum = -(t(k,i)+tlat(k,i)/cpp*deltat-(tmelt-5._r8))*cpp/xlf
                  dum = dum/(xlf/cpp*qrtot/(dz(k,i)*rho(k,i)))
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qniic(k,i)=qniic(k,i)+dum*qric(k,i)
               nsic(k,i)=nsic(k,i)+dum*nric(k,i)
               qric(k,i)=(1._r8-dum)*qric(k,i)
               nric(k,i)=(1._r8-dum)*nric(k,i)
               ! heating tendency 
               tmp = xlf*dum*qrtot/(dz(k,i)*rho(k,i))
               frzrdt(k,i)=frzrdt(k,i) + tmp

               tlat(k,i)=tlat(k,i)+tmp
               qstot=qstot+dum*qrtot
               qrtot=(1._r8-dum)*qrtot
               nstot=nstot+dum*nrtot
               nrtot=(1._r8-dum)*nrtot
               preci(i)=preci(i)+dum*(prect(i)-preci(i))
            end if
         end if
 
         ! if rain/snow mix ratio is zero so should number concentration

         if (qniic(k,i).lt.qsmall) then
            qniic(k,i)=0._r8
            nsic(k,i)=0._r8
         end if

         if (qric(k,i).lt.qsmall) then
            qric(k,i)=0._r8
            nric(k,i)=0._r8
         end if

         ! make sure number concentration is a positive number to avoid 
         ! taking root of negative

         nric(k,i)=max(nric(k,i),0._r8)
         nsic(k,i)=max(nsic(k,i),0._r8)

         !.......................................................................
         ! get size distribution parameters for fallspeed calculations
         !......................................................................
         ! rain

         if (qric(k,i).ge.qsmall) then
            lamr(k) = (pi*rhow*nric(k,i)/qric(k,i))**(1._r8/3._r8)
            n0r(k) = nric(k,i)*lamr(k)

            ! check for slope
            ! change lammax and lammin for rain and snow
            ! adjust vars

            if (lamr(k).lt.lamminr) then

               lamr(k) = lamminr

               n0r(k) = lamr(k)**4*qric(k,i)/(pi*rhow)
               nric(k,i) = n0r(k)/lamr(k)
            else if (lamr(k).gt.lammaxr) then
               lamr(k) = lammaxr
               n0r(k) = lamr(k)**4*qric(k,i)/(pi*rhow)
               nric(k,i) = n0r(k)/lamr(k)
            end if

            ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

            unr(k) = min(arn(k,i)*cons4/lamr(k)**br,9.1_r8*rhof(k,i))
            umr(k) = min(arn(k,i)*cons5/(6._r8*lamr(k)**br),9.1_r8*rhof(k,i))

         else
            lamr(k) = 0._r8
            n0r(k) = 0._r8
            umr(k)=0._r8
            unr(k)=0._r8
         end if

         !calculate mean size of combined rain and snow

         if (lamr(k).gt.0._r8) then
            Artmp = n0r(k) * pi / (2._r8 * lamr(k)**3._r8)
         else 
            Artmp = 0._r8
         endif

         if (lamc(k).gt.0._r8) then
            Actmp = cdist1(k) * pi * gamma(pgam(k)+3._r8)/(4._r8 * lamc(k)**2._r8)
         else 
            Actmp = 0._r8
         endif

         if (Actmp.gt.0_r8.or.Artmp.gt.0) then
            rercld(k,i)=rercld(k,i) + 3._r8 *(qric(k,i) + qcic(k,i)) / (4._r8 * rhow * (Actmp + Artmp))
            arcld(k,i)=arcld(k,i)+1._r8
         endif

         !......................................................................
         ! snow

         if (qniic(k,i).ge.qsmall) then
            lams(k) = (cons6*cs*nsic(k,i)/ &
                 qniic(k,i))**(1._r8/ds)
            n0s(k) = nsic(k,i)*lams(k)

            ! check for slope
            ! adjust vars

            if (lams(k).lt.lammins) then
               lams(k) = lammins
               n0s(k) = lams(k)**(ds+1._r8)*qniic(k,i)/(cs*cons6)
               nsic(k,i) = n0s(k)/lams(k)

            else if (lams(k).gt.lammaxs) then
               lams(k) = lammaxs
               n0s(k) = lams(k)**(ds+1._r8)*qniic(k,i)/(cs*cons6)
               nsic(k,i) = n0s(k)/lams(k)
            end if

            ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)

            ums(k) = min(asn(k,i)*cons8/(6._r8*lams(k)**bs),1.2_r8*rhof(k,i))
            uns(k) = min(asn(k,i)*cons7/lams(k)**bs,1.2_r8*rhof(k,i))

         else
            lams(k) = 0._r8
            n0s(k) = 0._r8
            ums(k) = 0._r8
            uns(k) = 0._r8
         end if

         !c........................................................................
         ! sum over sub-step for average process rates

         ! convert rain/snow q and N for output to history, note, 
         ! output is for gridbox average

         qrout(k,i)=qrout(k,i)+qric(k,i)*cldmax(k,i)
         qsout(k,i)=qsout(k,i)+qniic(k,i)*cldmax(k,i)
         nrout(k,i)=nrout(k,i)+nric(k,i)*rho(k,i)*cldmax(k,i)
         nsout(k,i)=nsout(k,i)+nsic(k,i)*rho(k,i)*cldmax(k,i)

         tlat1(k,i)=tlat1(k,i)+tlat(k,i)
         qvlat1(k,i)=qvlat1(k,i)+qvlat(k,i)
         qctend1(k,i)=qctend1(k,i)+qctend(k,i)
         qitend1(k,i)=qitend1(k,i)+qitend(k,i)
         nctend1(k,i)=nctend1(k,i)+nctend(k,i)
         nitend1(k,i)=nitend1(k,i)+nitend(k,i)

         t(k,i)=t(k,i)+tlat(k,i)*deltat/cpp
         q(k,i)=q(k,i)+qvlat(k,i)*deltat
         qc(k,i)=qc(k,i)+qctend(k,i)*deltat
         qi(k,i)=qi(k,i)+qitend(k,i)*deltat
         nc(k,i)=nc(k,i)+nctend(k,i)*deltat
         ni(k,i)=ni(k,i)+nitend(k,i)*deltat

         rainrt1(k,i)=rainrt1(k,i)+rainrt(k,i)

         !divide rain radius over substeps for average
         if (arcld(k,i) .gt. 0._r8) then
            rercld(k,i)=rercld(k,i)/arcld(k,i)
         end if

         !calculate precip fluxes and adding them to summing sub-stepping variables
         !! flux is zero at top interface
         rflx(1,i)=0.0_r8
         sflx(1,i)=0.0_r8

         !! calculating the precip flux (kg/m2/s) as mixingratio(kg/kg)*airdensity(kg/m3)*massweightedfallspeed(m/s)
         rflx(k+1,i)=qrout(k,i)*rho(k,i)*umr(k)
         sflx(k+1,i)=qsout(k,i)*rho(k,i)*ums(k)

         !! add to summing sub-stepping variable
         rflx1(k+1,i)=rflx1(k+1,i)+rflx(k+1,i)
         sflx1(k+1,i)=sflx1(k+1,i)+sflx(k+1,i)

         !c........................................................................

      end do ! k loop

      prect1(i)=prect1(i)+prect(i)
      preci1(i)=preci1(i)+preci(i)

   end do ! it loop, sub-step

   do k = top_lev, pver
      rate1ord_cw2pr_st(k,i) = qcsinksum_rate1ord(k)/max(qcsum_rate1ord(k),1.0e-30_r8) 
   end do

300 continue  ! continue if no cloud water
end do ! i loop

! convert dt from sub-step back to full time step
deltat=deltat*real(iter)

!c.............................................................................

do i=1,ncol

   ! skip all calculations if no cloud water
   if (ltrue(i).eq.0) then

      do k=1,top_lev-1
         ! assign zero values for effective radius above 1 mbar
         effc(k,i)=0._r8
         effi(k,i)=0._r8
         effc_fn(k,i)=0._r8
         lamcrad(k,i)=0._r8
         pgamrad(k,i)=0._r8
         deffi(k,i)=0._r8
      end do

      do k=top_lev,pver
         ! assign default values for effective radius
         effc(k,i)=10._r8
         effi(k,i)=25._r8
         effc_fn(k,i)=10._r8
         lamcrad(k,i)=0._r8
         pgamrad(k,i)=0._r8
         deffi(k,i)=0._r8
      end do
      goto 500
   end if

   ! initialize nstep for sedimentation sub-steps
   nstep = 1

   ! divide precip rate by number of sub-steps to get average over time step

   prect(i)=prect1(i)/real(iter)
   preci(i)=preci1(i)/real(iter)

   do k=top_lev,pver

      ! assign variables back to start-of-timestep values before updating after sub-steps 

      t(k,i)=t1(k,i)
      q(k,i)=q1(k,i)
      qc(k,i)=qc1(k,i)
      qi(k,i)=qi1(k,i)
      nc(k,i)=nc1(k,i)
      ni(k,i)=ni1(k,i)

      ! divide microphysical tendencies by number of sub-steps to get average over time step

      tlat(k,i)=tlat1(k,i)/real(iter)
      qvlat(k,i)=qvlat1(k,i)/real(iter)
      qctend(k,i)=qctend1(k,i)/real(iter)
      qitend(k,i)=qitend1(k,i)/real(iter)
      nctend(k,i)=nctend1(k,i)/real(iter)
      nitend(k,i)=nitend1(k,i)/real(iter)

      rainrt(k,i)=rainrt1(k,i)/real(iter)

      ! divide by number of sub-steps to find final values
      rflx(k+1,i)=rflx1(k+1,i)/real(iter)
      sflx(k+1,i)=sflx1(k+1,i)/real(iter)

      ! divide output precip q and N by number of sub-steps to get average over time step

      qrout(k,i)=qrout(k,i)/real(iter)
      qsout(k,i)=qsout(k,i)/real(iter)
      nrout(k,i)=nrout(k,i)/real(iter)
      nsout(k,i)=nsout(k,i)/real(iter)

      ! divide trop_mozart variables by number of sub-steps to get average over time step 

      nevapr(k,i) = nevapr(k,i)/real(iter)
      evapsnow(k,i) = evapsnow(k,i)/real(iter)
      prain(k,i) = prain(k,i)/real(iter)
      prodsnow(k,i) = prodsnow(k,i)/real(iter)
      cmeout(k,i) = cmeout(k,i)/real(iter)

      cmeiout(k,i) = cmeiout(k,i)/real(iter)
      meltsdt(k,i) = meltsdt(k,i)/real(iter)
      frzrdt (k,i) = frzrdt (k,i)/real(iter)

      ! microphysics output
      prao(k,i)=prao(k,i)/real(iter)
      prco(k,i)=prco(k,i)/real(iter)
      mnuccco(k,i)=mnuccco(k,i)/real(iter)
      mnuccto(k,i)=mnuccto(k,i)/real(iter)
      msacwio(k,i)=msacwio(k,i)/real(iter)
      psacwso(k,i)=psacwso(k,i)/real(iter)
      bergso(k,i)=bergso(k,i)/real(iter)
      bergo(k,i)=bergo(k,i)/real(iter)
      prcio(k,i)=prcio(k,i)/real(iter)
      praio(k,i)=praio(k,i)/real(iter)

      mnuccro(k,i)=mnuccro(k,i)/real(iter)
      pracso (k,i)=pracso (k,i)/real(iter)

      mnuccdo(k,i)=mnuccdo(k,i)/real(iter)

      ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
      nevapr(k,i) = nevapr(k,i) + evapsnow(k,i)
      prain(k,i) = prain(k,i) + prodsnow(k,i)

      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! calculate sedimentation for cloud water and ice
      !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ! update in-cloud cloud mixing ratio and number concentration 
      ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
      ! note: these are in-cloud values***, hence we divide by cloud fraction

      dumc(k,i) = (qc(k,i)+qctend(k,i)*deltat)/lcldm(k,i)
      dumi(k,i) = (qi(k,i)+qitend(k,i)*deltat)/icldm(k,i)
      dumnc(k,i) = max((nc(k,i)+nctend(k,i)*deltat)/lcldm(k,i),0._r8)
      dumni(k,i) = max((ni(k,i)+nitend(k,i)*deltat)/icldm(k,i),0._r8)

      ! obtain new slope parameter to avoid possible singularity

      if (dumi(k,i).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumni(k,i)=min(dumni(k,i),dumi(k,i)*1.e20_r8)

         lami(k) = (cons1*ci* &
              dumni(k,i)/dumi(k,i))**(1._r8/di)
         lami(k)=max(lami(k),lammini)
         lami(k)=min(lami(k),lammaxi)
      else
         lami(k)=0._r8
      end if

      if (dumc(k,i).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumnc(k,i)=min(dumnc(k,i),dumc(k,i)*1.e20_r8)
         ! add lower limit to in-cloud number concentration
         dumnc(k,i)=max(dumnc(k,i),cdnl/rho(k,i)) ! sghan minimum in #/cm3 
         pgam(k)=0.0005714_r8*(ncic(k,i)/1.e6_r8*rho(k,i))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(k,i)*gamma(pgam(k)+4._r8)/ &
              (dumc(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         lammax = (pgam(k)+1._r8)/2.e-6_r8
         lamc(k)=max(lamc(k),lammin)
         lamc(k)=min(lamc(k),lammax)
      else
         lamc(k)=0._r8
      end if

      ! calculate number and mass weighted fall velocity for droplets
      ! include effects of sub-grid distribution of cloud water


      if (dumc(k,i).ge.qsmall) then
         unc = acn(k,i)*gamma(1._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+1._r8))
         umc = acn(k,i)*gamma(4._r8+bc+pgam(k))/(lamc(k)**bc*gamma(pgam(k)+4._r8))
         ! fallspeed for output
         vtrmc(k,i)=umc
      else
         umc = 0._r8
         unc = 0._r8
      end if

      ! calculate number and mass weighted fall velocity for cloud ice

      if (dumi(k,i).ge.qsmall) then
         uni =  ain(k,i)*cons16/lami(k)**bi
         umi = ain(k,i)*cons17/(6._r8*lami(k)**bi)
         uni=min(uni,1.2_r8*rhof(k,i))
         umi=min(umi,1.2_r8*rhof(k,i))

         ! fallspeed
         vtrmi(k,i)=umi
      else
         umi = 0._r8
         uni = 0._r8
      end if

      fi(k) = g*rho(k,i)*umi
      fni(k) = g*rho(k,i)*uni
      fc(k) = g*rho(k,i)*umc
      fnc(k) = g*rho(k,i)*unc

      ! calculate number of split time steps to ensure courant stability criteria
      ! for sedimentation calculations

      rgvm = max(fi(k),fc(k),fni(k),fnc(k))
      nstep = max(int(rgvm*deltat/pdel(k,i)+1._r8),nstep)

      ! redefine dummy variables - sedimentation is calculated over grid-scale
      ! quantities to ensure conservation

      dumc(k,i) = (qc(k,i)+qctend(k,i)*deltat)
      dumi(k,i) = (qi(k,i)+qitend(k,i)*deltat)
      dumnc(k,i) = max((nc(k,i)+nctend(k,i)*deltat),0._r8)
      dumni(k,i) = max((ni(k,i)+nitend(k,i)*deltat),0._r8)

      if (dumc(k,i).lt.qsmall) dumnc(k,i)=0._r8
      if (dumi(k,i).lt.qsmall) dumni(k,i)=0._r8

   end do       !!! vertical loop


   do n = 1,nstep  !! loop over sub-time step to ensure stability	

      do k = top_lev,pver
         if (do_cldice) then
            falouti(k) = fi(k)*dumi(k,i)
            faloutni(k) = fni(k)*dumni(k,i)
         else
            falouti(k)  = 0._r8
            faloutni(k) = 0._r8
         end if

         faloutc(k) = fc(k)*dumc(k,i)
         faloutnc(k) = fnc(k)*dumnc(k,i)
      end do

      ! top of model

      k = top_lev
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

      do k = top_lev+1,pver

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

         Fni(K)=MAX(Fni(K)/pdel(k,i),Fni(K-1)/pdel(k-1,i))*pdel(k,i)
         FI(K)=MAX(FI(K)/pdel(k,i),FI(K-1)/pdel(k-1,i))*pdel(k,i)
         fnc(k)=max(fnc(k)/pdel(k,i),fnc(k-1)/pdel(k-1,i))*pdel(k,i)
         Fc(K)=MAX(Fc(K)/pdel(k,i),Fc(K-1)/pdel(k-1,i))*pdel(k,i)

      end do   !! k loop

      ! units below are m/s
      ! cloud water/ice sedimentation flux at surface 
      ! is added to precip flux at surface to get total precip (cloud + precip water)
      ! rate

      prect(i) = prect(i)+(faloutc(pver)+falouti(pver))/g/nstep/1000._r8  
      preci(i) = preci(i)+(falouti(pver))/g/nstep/1000._r8

   end do   !! nstep loop

   ! end sedimentation
   !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   ! get new update for variables that includes sedimentation tendency
   ! note : here dum variables are grid-average, NOT in-cloud

   do k=top_lev,pver

      dumc(k,i) = max(qc(k,i)+qctend(k,i)*deltat,0._r8)
      dumi(k,i) = max(qi(k,i)+qitend(k,i)*deltat,0._r8)
      dumnc(k,i) = max(nc(k,i)+nctend(k,i)*deltat,0._r8)
      dumni(k,i) = max(ni(k,i)+nitend(k,i)*deltat,0._r8)

      if (dumc(k,i).lt.qsmall) dumnc(k,i)=0._r8
      if (dumi(k,i).lt.qsmall) dumni(k,i)=0._r8

      ! calculate instantaneous processes (melting, homogeneous freezing)
      if (do_cldice) then

         if (t(k,i)+tlat(k,i)/cpp*deltat > tmelt) then
            if (dumi(k,i) > 0._r8) then

               ! limit so that melting does not push temperature below freezing
               dum = -dumi(k,i)*xlf/cpp
               if (t(k,i)+tlat(k,i)/cpp*deltat+dum.lt.tmelt) then
                  dum = (t(k,i)+tlat(k,i)/cpp*deltat-tmelt)*cpp/xlf
                  dum = dum/dumi(k,i)*xlf/cpp 
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qctend(k,i)=qctend(k,i)+dum*dumi(k,i)/deltat

               ! for output
               melto(k,i)=dum*dumi(k,i)/deltat

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

         if (t(k,i)+tlat(k,i)/cpp*deltat < 233.15_r8) then
            if (dumc(k,i) > 0._r8) then

               ! limit so that freezing does not push temperature above threshold
               dum = dumc(k,i)*xlf/cpp
               if (t(k,i)+tlat(k,i)/cpp*deltat+dum.gt.233.15_r8) then
                  dum = -(t(k,i)+tlat(k,i)/cpp*deltat-233.15_r8)*cpp/xlf
                  dum = dum/dumc(k,i)*xlf/cpp
                  dum = max(0._r8,dum)
                  dum = min(1._r8,dum)
               else
                  dum = 1._r8
               end if

               qitend(k,i)=qitend(k,i)+dum*dumc(k,i)/deltat
               ! for output
               homoo(k,i)=dum*dumc(k,i)/deltat

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
         ! follow code similar to old CAM scheme

         qtmp=q(k,i)+qvlat(k,i)*deltat
         ttmp=t(k,i)+tlat(k,i)/cpp*deltat

         esn = svp_water(ttmp)  ! use rhw to allow ice supersaturation
         qsn = svp_to_qsat(esn, p(k,i))

         if (qtmp > qsn .and. qsn > 0) then
            ! expression below is approximate since there may be ice deposition
            dum = (qtmp-qsn)/(1._r8+cons27*qsn/(cpp*rv*ttmp**2))/deltat
            ! add to output cme
            cmeout(k,i) = cmeout(k,i)+dum
            ! now add to tendencies, partition between liquid and ice based on temperature
            if (ttmp > 268.15_r8) then
               dum1=0.0_r8
               ! now add to tendencies, partition between liquid and ice based on te
            else if (ttmp < 238.15_r8) then
               dum1=1.0_r8
            else
               dum1=(268.15_r8-ttmp)/30._r8
            end if

            dum = (qtmp-qsn)/(1._r8+(xxls*dum1+xxlv*(1._r8-dum1))**2 &
                 *qsn/(cpp*rv*ttmp**2))/deltat
            qctend(k,i)=qctend(k,i)+dum*(1._r8-dum1)
            ! for output
            qcreso(k,i)=dum*(1._r8-dum1)
            qitend(k,i)=qitend(k,i)+dum*dum1
            qireso(k,i)=dum*dum1
            qvlat(k,i)=qvlat(k,i)-dum
            ! for output
            qvres(k,i)=-dum
            tlat(k,i)=tlat(k,i)+dum*(1._r8-dum1)*xxlv+dum*dum1*xxls
         end if
      end if

      !...............................................................................
      ! calculate effective radius for pass to radiation code
      ! if no cloud water, default value is 10 micron for droplets,
      ! 25 micron for cloud ice

      ! update cloud variables after instantaneous processes to get effective radius
      ! variables are in-cloud to calculate size dist parameters

      dumc(k,i) = max(qc(k,i)+qctend(k,i)*deltat,0._r8)/lcldm(k,i)
      dumi(k,i) = max(qi(k,i)+qitend(k,i)*deltat,0._r8)/icldm(k,i)
      dumnc(k,i) = max(nc(k,i)+nctend(k,i)*deltat,0._r8)/lcldm(k,i)
      dumni(k,i) = max(ni(k,i)+nitend(k,i)*deltat,0._r8)/icldm(k,i)

      ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
      dumc(k,i)=min(dumc(k,i),5.e-3_r8)
      dumi(k,i)=min(dumi(k,i),5.e-3_r8)

      !...................
      ! cloud ice effective radius

      if (dumi(k,i).ge.qsmall) then
         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumni(k,i)=min(dumni(k,i),dumi(k,i)*1.e20_r8)
         lami(k) = (cons1*ci*dumni(k,i)/dumi(k,i))**(1._r8/di)

         if (lami(k).lt.lammini) then
            lami(k) = lammini
            n0i(k) = lami(k)**(di+1._r8)*dumi(k,i)/(ci*cons1)
            niic(k,i) = n0i(k)/lami(k)
            ! adjust number conc if needed to keep mean size in reasonable range
            if (do_cldice) nitend(k,i)=(niic(k,i)*icldm(k,i)-ni(k,i))/deltat

         else if (lami(k).gt.lammaxi) then
            lami(k) = lammaxi
            n0i(k) = lami(k)**(di+1._r8)*dumi(k,i)/(ci*cons1)
            niic(k,i) = n0i(k)/lami(k)
            ! adjust number conc if needed to keep mean size in reasonable range
            if (do_cldice) nitend(k,i)=(niic(k,i)*icldm(k,i)-ni(k,i))/deltat
         end if
         effi(k,i) = 1.5_r8/lami(k)*1.e6_r8

      else
         effi(k,i) = 25._r8
      end if

      ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
      ! radius has already been determined from the size distribution.
      if (.not. do_cldice) then
         effi(k,i) = re_ice(k,i) * 1e6_r8      ! m -> um
      end if

      !...................
      ! cloud droplet effective radius

      if (dumc(k,i).ge.qsmall) then

         ! add upper limit to in-cloud number concentration to prevent numerical error
         dumnc(k,i)=min(dumnc(k,i),dumc(k,i)*1.e20_r8)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! set tendency to ensure minimum droplet concentration
         ! after update by microphysics, except when lambda exceeds bounds on mean drop
         ! size or if there is no cloud water
         if (dumnc(k,i).lt.cdnl/rho(k,i)) then   
            nctend(k,i)=(cdnl/rho(k,i)*lcldm(k,i)-nc(k,i))/deltat   
         end if
         dumnc(k,i)=max(dumnc(k,i),cdnl/rho(k,i)) ! sghan minimum in #/cm3 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         pgam(k)=0.0005714_r8*(ncic(k,i)/1.e6_r8*rho(k,i))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(k,i)*gamma(pgam(k)+4._r8)/ &
              (dumc(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         lammax = (pgam(k)+1._r8)/2.e-6_r8
         if (lamc(k).lt.lammin) then
            lamc(k) = lammin
            ncic(k,i) = 6._r8*lamc(k)**3*dumc(k,i)* &
                 gamma(pgam(k)+1._r8)/ &
                 (pi*rhow*gamma(pgam(k)+4._r8))
            ! adjust number conc if needed to keep mean size in reasonable range
            nctend(k,i)=(ncic(k,i)*lcldm(k,i)-nc(k,i))/deltat

         else if (lamc(k).gt.lammax) then
            lamc(k) = lammax
            ncic(k,i) = 6._r8*lamc(k)**3*dumc(k,i)* &
                 gamma(pgam(k)+1._r8)/ &
                 (pi*rhow*gamma(pgam(k)+4._r8))
            ! adjust number conc if needed to keep mean size in reasonable range
            nctend(k,i)=(ncic(k,i)*lcldm(k,i)-nc(k,i))/deltat
         end if

         effc(k,i) = &
              gamma(pgam(k)+4._r8)/ &
              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8
         !assign output fields for shape here
         lamcrad(k,i)=lamc(k)
         pgamrad(k,i)=pgam(k)

      else
         effc(k,i) = 10._r8
         lamcrad(k,i)=0._r8
         pgamrad(k,i)=0._r8
      end if

      ! ice effective diameter for david mitchell's optics
      if (do_cldice) then
         deffi(k,i)=effi(k,i)*rhoi/917._r8*2._r8
      else
         deffi(k,i)=effi(k,i) * 2._r8
      end if

!!! recalculate effective radius for constant number, in order to separate
      ! first and second indirect effects
      ! assume constant number of 10^8 kg-1

      dumnc(k,i)=1.e8_r8

      if (dumc(k,i).ge.qsmall) then
         pgam(k)=0.0005714_r8*(ncic(k,i)/1.e6_r8*rho(k,i))+0.2714_r8
         pgam(k)=1._r8/(pgam(k)**2)-1._r8
         pgam(k)=max(pgam(k),2._r8)
         pgam(k)=min(pgam(k),15._r8)

         lamc(k) = (pi/6._r8*rhow*dumnc(k,i)*gamma(pgam(k)+4._r8)/ &
              (dumc(k,i)*gamma(pgam(k)+1._r8)))**(1._r8/3._r8)
         lammin = (pgam(k)+1._r8)/50.e-6_r8
         lammax = (pgam(k)+1._r8)/2.e-6_r8
         if (lamc(k).lt.lammin) then
            lamc(k) = lammin
         else if (lamc(k).gt.lammax) then
            lamc(k) = lammax
         end if
         effc_fn(k,i) = &
              gamma(pgam(k)+4._r8)/ &
              gamma(pgam(k)+3._r8)/lamc(k)/2._r8*1.e6_r8

      else
         effc_fn(k,i) = 10._r8
      end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

   end do ! vertical k loop

500 continue

   do k=top_lev,pver
      ! if updated q (after microphysics) is zero, then ensure updated n is also zero

      if (qc(k,i)+qctend(k,i)*deltat.lt.qsmall) nctend(k,i)=-nc(k,i)/deltat
      if (do_cldice .and. qi(k,i)+qitend(k,i)*deltat.lt.qsmall) nitend(k,i)=-ni(k,i)/deltat
   end do

end do ! i loop

! add snow ouptut
do i = 1,ncol
   do k=top_lev,pver
      if (qsout(k,i).gt.1.e-7_r8.and.nsout(k,i).gt.0._r8) then
         !---------------Lin Han------------------> 
         !dsout(k,i)=3._r8*rhosn/917._r8*(pi * rhosn * nsout(k,i)/qsout(k,i))**(-1._r8/3._r8)

         rsout(k,i)=1.5_r8*(pi * rhosn * nsout(k,i)/qsout(k,i))**(-1._r8/3._r8)
         dsout(k,i)=2._r8*rsout(k,i)*rhosn/917._r8
         !<--------------Lin Han------------------- 
      endif
   end do
end do

!calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual
do i = 1,ncol
   do k=top_lev,pver
      !! RAIN
      if (qrout(k,i).gt.1.e-7_r8.and.nrout(k,i).gt.0._r8) then
         reff_rain(k,i)=1.5_r8*(pi * rhow * nrout(k,i)/qrout(k,i))**(-1._r8/3._r8)*1.e6_r8
      endif
      !! SNOW
      if (qsout(k,i).gt.1.e-7_r8.and.nsout(k,i).gt.0._r8) then
         reff_snow(k,i)=1.5_r8*(pi * rhosn * nsout(k,i)/qsout(k,i))**(-1._r8/3._r8)*1.e6_r8
      end if
   end do
end do

! analytic radar reflectivity
! formulas from Matthew Shupe, NOAA/CERES
! *****note: radar reflectivity is local (in-precip average)
! units of mm^6/m^3

do i = 1,ncol
   do k=top_lev,pver
      if (qc(k,i)+qctend(k,i)*deltat.ge.qsmall) then
         dum=((qc(k,i)+qctend(k,i)*deltat)/lcldm(k,i)*rho(k,i)*1000._r8)**2 &
              /(0.109_r8*(nc(k,i)+nctend(k,i)*deltat)/lcldm(k,i)*rho(k,i)/1.e6_r8)*lcldm(k,i)/cldmax(k,i)
      else
         dum=0._r8
      end if
      if (qi(k,i)+qitend(k,i)*deltat.ge.qsmall) then
         dum1=((qi(k,i)+qitend(k,i)*deltat)*rho(k,i)/icldm(k,i)*1000._r8/0.1_r8)**(1._r8/0.63_r8)*icldm(k,i)/cldmax(k,i)
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
      areflz(k,i)=refl(k,i)

      ! convert back to DBz 

      if (refl(k,i).gt.minrefl) then 
         refl(k,i)=10._r8*log10(refl(k,i))
      else
         refl(k,i)=-9999._r8
      end if

      !set averaging flag
      if (refl(k,i).gt.mindbz) then 
         arefl(k,i)=refl(k,i)
         frefl(k,i)=1.0_r8  
      else
         arefl(k,i)=0._r8
         areflz(k,i)=0._r8
         frefl(k,i)=0._r8
      end if

      ! bound cloudsat reflectivity

      csrfl(k,i)=min(csmax,refl(k,i))

      !set averaging flag
      if (csrfl(k,i).gt.csmin) then 
         acsrfl(k,i)=refl(k,i)
         fcsrfl(k,i)=1.0_r8  
      else
         acsrfl(k,i)=0._r8
         fcsrfl(k,i)=0._r8
      end if

   end do
end do

! averaging for snow and rain number and diameter

qrout2(:,:)=0._r8
qsout2(:,:)=0._r8
nrout2(:,:)=0._r8
nsout2(:,:)=0._r8
drout2(:,:)=0._r8
dsout2(:,:)=0._r8
freqs(:,:)=0._r8
freqr(:,:)=0._r8
do i = 1,ncol
   do k=top_lev,pver
      if (qrout(k,i).gt.1.e-7_r8.and.nrout(k,i).gt.0._r8) then
         qrout2(k,i)=qrout(k,i)
         nrout2(k,i)=nrout(k,i)
         drout2(k,i)=(pi * rhow * nrout(k,i)/qrout(k,i))**(-1._r8/3._r8)
         freqr(k,i)=1._r8
      endif
      if (qsout(k,i).gt.1.e-7_r8.and.nsout(k,i).gt.0._r8) then
         qsout2(k,i)=qsout(k,i)
         nsout2(k,i)=nsout(k,i)
         dsout2(k,i)=(pi * rhosn * nsout(k,i)/qsout(k,i))**(-1._r8/3._r8)
         freqs(k,i)=1._r8
      endif
   end do
end do

! output activated liquid and ice (convert from #/kg -> #/m3)
do i = 1,ncol
   do k=top_lev,pver
      ncai(k,i)=dum2i(k,i)*rho(k,i)
      ncal(k,i)=dum2l(k,i)*rho(k,i)
   end do
end do


!redefine fice here....
nfice(:,:)=0._r8
do k=top_lev,pver
   do i=1,ncol
      dumc(k,i) = (qc(k,i)+qctend(k,i)*deltat)
      dumi(k,i) = (qi(k,i)+qitend(k,i)*deltat)
      dumfice=qsout(k,i) + qrout(k,i) + dumc(k,i) + dumi(k,i)  

      if (dumfice.gt.qsmall.and.(qsout(k,i)+dumi(k,i).gt.qsmall)) then
         nfice(k,i)=(qsout(k,i) + dumi(k,i))/dumfice
      endif

      if (nfice(k,i).gt.1._r8) then
         nfice(k,i)=1._r8
      endif

   enddo
enddo

end subroutine micro_mg1_0_tend


end module micro_mg1_0
