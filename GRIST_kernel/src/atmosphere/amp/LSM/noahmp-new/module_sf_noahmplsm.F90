MODULE MODULE_SF_NOAHMPLSM
!-> xinyf
 USE grist_constants, only: r8
!<-
  IMPLICIT NONE

  public  :: noahmp_options
  public  :: NOAHMP_SFLX

  private :: ATM
  private :: PHENOLOGY
  private :: PRECIP_HEAT
  private :: ENERGY
  private ::       THERMOPROP
  private ::               CSNOW
  private ::               TDFCND
  private ::       RADIATION
  private ::               ALBEDO
  private ::                         SNOW_AGE
  private ::                         SNOWALB_BATS  
  private ::                         SNOWALB_CLASS
  private ::                         GROUNDALB
  private ::                         TWOSTREAM
  private ::               SURRAD
  private ::       VEGE_FLUX
  private ::               SFCDIF1                  
  private ::               SFCDIF2                
  private ::               STOMATA                  
  private ::               CANRES                  
  private ::               ESAT
  private ::               RAGRB
  private ::       BARE_FLUX
  private ::       TSNOSOI
  private ::               HRT
  private ::               HSTEP   
  private ::                         ROSR12
  private ::       PHASECHANGE
  private ::               FRH2O           

  private :: WATER
  private ::       CANWATER
  private ::       SNOWWATER
  private ::               SNOWFALL
  private ::               COMBINE
  private ::               DIVIDE
  private ::                         COMBO
  private ::               COMPACT
  private ::               SNOWH2O
  private ::       SOILWATER
  private ::               ZWTEQ
  private ::               INFIL
  private ::               SRT
  private ::                         WDFCND1        
  private ::                         WDFCND2       
  private ::               SSTEP
  private ::       GROUNDWATER
  private ::       SHALLOWWATERTABLE

  private :: CARBON
  private ::       CO2FLUX
!  private ::       BVOCFLUX
!  private ::       CH4FLUX

  private :: ERROR

! =====================================options for different schemes================================
! **recommended

  INTEGER :: DVEG     ! options for dynamic vegetation: 
                      !   1 -> off (use table LAI; use FVEG = SHDFAC from input)
                      !   2 -> on  (together with OPT_CRS = 1)
                      !   3 -> off (use table LAI; calculate FVEG)
                      ! **4 -> off (use table LAI; use maximum vegetation fraction)
                      ! **5 -> on  (use maximum vegetation fraction)
                      !   6 -> on  (use FVEG = SHDFAC from input)
                      !   7 -> off (use input LAI; use FVEG = SHDFAC from input)
                      !   8 -> off (use input LAI; calculate FVEG)
                      !   9 -> off (use input LAI; use maximum vegetation fraction)
                      !  10 -> crop model on (use maximum vegetation fraction)

  INTEGER :: OPT_CRS  ! options for canopy stomatal resistance
                      ! **1 -> Ball-Berry
		      !   2 -> Jarvis

  INTEGER :: OPT_BTR  ! options for soil moisture factor for stomatal resistance
                      ! **1 -> Noah (soil moisture) 
                      !   2 -> CLM  (matric potential)
                      !   3 -> SSiB (matric potential)

  INTEGER :: OPT_RUN  ! options for runoff and groundwater
                      ! **1 -> TOPMODEL with groundwater (Niu et al. 2007 JGR) ;
                      !   2 -> TOPMODEL with an equilibrium water table (Niu et al. 2005 JGR) ;
                      !   3 -> original surface and subsurface runoff (free drainage)
                      !   4 -> BATS surface and subsurface runoff (free drainage)
                      !   5 -> Miguez-Macho&Fan groundwater scheme (Miguez-Macho et al. 2007 JGR; Fan et al. 2007 JGR)
		      !          (needs further testing for public use)

  INTEGER :: OPT_SFC  ! options for surface layer drag coeff (CH & CM)
                      ! **1 -> M-O
		      ! **2 -> original Noah (Chen97)
		      ! **3 -> MYJ consistent; 4->YSU consistent. MB: removed in v3.7 for further testing

  INTEGER :: OPT_FRZ  ! options for supercooled liquid water (or ice fraction)
                      ! **1 -> no iteration (Niu and Yang, 2006 JHM)
		      !   2 -> Koren's iteration 

  INTEGER :: OPT_INF  ! options for frozen soil permeability
                      ! **1 -> linear effects, more permeable (Niu and Yang, 2006, JHM)
                      !   2 -> nonlinear effects, less permeable (old)

  INTEGER :: OPT_RAD  ! options for radiation transfer
                      !   1 -> modified two-stream (gap = F(solar angle, 3D structure ...)<1-FVEG)
                      !   2 -> two-stream applied to grid-cell (gap = 0)
                      ! **3 -> two-stream applied to vegetated fraction (gap=1-FVEG)

  INTEGER :: OPT_ALB  ! options for ground snow surface albedo
                      !   1 -> BATS
		      ! **2 -> CLASS

  INTEGER :: OPT_SNF  ! options for partitioning  precipitation into rainfall & snowfall
                      ! **1 -> Jordan (1991)
		      !   2 -> BATS: when SFCTMP<TFRZ+2.2 
		      !   3 -> SFCTMP < TFRZ
		      !   4 -> Use WRF microphysics output

  INTEGER :: OPT_TBOT ! options for lower boundary condition of soil temperature
                      !   1 -> zero heat flux from bottom (ZBOT and TBOT not used)
                      ! **2 -> TBOT at ZBOT (8m) read from a file (original Noah)

  INTEGER :: OPT_STC  ! options for snow/soil temperature time scheme (only layer 1)
                      ! **1 -> semi-implicit; flux top boundary condition
		      !   2 -> full implicit (original Noah); temperature top boundary condition
                      !   3 -> same as 1, but FSNO for TS calculation (generally improves snow; v3.7)

  INTEGER :: OPT_RSF  ! options for surface resistent to evaporation/sublimation
                      ! **1 -> Sakaguchi and Zeng, 2009
		      !   2 -> Sellers (1992)
                      !   3 -> adjusted Sellers to decrease RSURF for wet soil
		      !   4 -> option 1 for non-snow; rsurf = rsurf_snow for snow (set in MPTABLE); AD v3.8

!------------------------------------------------------------------------------------------!
! Physical Constants:                                                                      !
!------------------------------------------------------------------------------------------!

  REAL(r8), PARAMETER :: GRAV   = 9.80616   !acceleration due to gravity (m/s2)
  REAL(r8), PARAMETER :: SB     = 5.67E-08  !Stefan-Boltzmann constant (w/m2/k4)
  REAL(r8), PARAMETER :: VKC    = 0.40      !von Karman constant
  REAL(r8), PARAMETER :: TFRZ   = 273.16    !freezing/melting point (k)
  REAL(r8), PARAMETER :: HSUB   = 2.8440E06 !latent heat of sublimation (j/kg)
  REAL(r8), PARAMETER :: HVAP   = 2.5104E06 !latent heat of vaporization (j/kg)
  REAL(r8), PARAMETER :: HFUS   = 0.3336E06 !latent heat of fusion (j/kg)
  REAL(r8), PARAMETER :: CWAT   = 4.188E06  !specific heat capacity of water (j/m3/k)
  REAL(r8), PARAMETER :: CICE   = 2.094E06  !specific heat capacity of ice (j/m3/k)
  REAL(r8), PARAMETER :: CPAIR  = 1004.64   !heat capacity dry air at const pres (j/kg/k)
  REAL(r8), PARAMETER :: TKWAT  = 0.6       !thermal conductivity of water (w/m/k)
  REAL(r8), PARAMETER :: TKICE  = 2.2       !thermal conductivity of ice (w/m/k)
  REAL(r8), PARAMETER :: TKAIR  = 0.023     !thermal conductivity of air (w/m/k) (not used MB: 20140718)
  REAL(r8), PARAMETER :: RAIR   = 287.04    !gas constant for dry air (j/kg/k)
  REAL(r8), PARAMETER :: RW     = 461.269   !gas constant for  water vapor (j/kg/k)
  REAL(r8), PARAMETER :: DENH2O = 1000.     !density of water (kg/m3)
  REAL(r8), PARAMETER :: DENICE = 917.      !density of ice (kg/m3)

  INTEGER, PRIVATE, PARAMETER :: MBAND = 2
  INTEGER, PRIVATE, PARAMETER :: NSOIL = 4
  INTEGER, PRIVATE, PARAMETER :: NSTAGE = 8

  TYPE noahmp_parameters ! define a NoahMP parameters type

!------------------------------------------------------------------------------------------!
! From the veg section of MPTABLE.TBL
!------------------------------------------------------------------------------------------!

    LOGICAL :: URBAN_FLAG
    INTEGER :: ISWATER
    INTEGER :: ISBARREN
    INTEGER :: ISICE
    INTEGER :: EBLFOREST

    REAL(r8) :: CH2OP              !maximum intercepted h2o per unit lai+sai (mm)
    REAL(r8) :: DLEAF              !characteristic leaf dimension (m)
    REAL(r8) :: Z0MVT              !momentum roughness length (m)
    REAL(r8) :: HVT                !top of canopy (m)
    REAL(r8) :: HVB                !bottom of canopy (m)
    REAL(r8) :: DEN                !tree density (no. of trunks per m2)
    REAL(r8) :: RC                 !tree crown radius (m)
    REAL(r8) :: MFSNO              !snowmelt m parameter ()
    REAL(r8) :: SAIM(12)           !monthly stem area index, one-sided
    REAL(r8) :: LAIM(12)           !monthly leaf area index, one-sided
    REAL(r8) :: SLA                !single-side leaf area per Kg [m2/kg]
    REAL(r8) :: DILEFC             !coeficient for leaf stress death [1/s]
    REAL(r8) :: DILEFW             !coeficient for leaf stress death [1/s]
    REAL(r8) :: FRAGR              !fraction of growth respiration  !original was 0.3 
    REAL(r8) :: LTOVRC             !leaf turnover [1/s]

    REAL(r8) :: C3PSN              !photosynthetic pathway: 0. = c4, 1. = c3
    REAL(r8) :: KC25               !co2 michaelis-menten constant at 25c (pa)
    REAL(r8) :: AKC                !q10 for kc25
    REAL(r8) :: KO25               !o2 michaelis-menten constant at 25c (pa)
    REAL(r8) :: AKO                !q10 for ko25
    REAL(r8) :: VCMX25             !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    REAL(r8) :: AVCMX              !q10 for vcmx25
    REAL(r8) :: BP                 !minimum leaf conductance (umol/m**2/s)
    REAL(r8) :: MP                 !slope of conductance-to-photosynthesis relationship
    REAL(r8) :: QE25               !quantum efficiency at 25c (umol co2 / umol photon)
    REAL(r8) :: AQE                !q10 for qe25
    REAL(r8) :: RMF25              !leaf maintenance respiration at 25c (umol co2/m**2/s)
    REAL(r8) :: RMS25              !stem maintenance respiration at 25c (umol co2/kg bio/s)
    REAL(r8) :: RMR25              !root maintenance respiration at 25c (umol co2/kg bio/s)
    REAL(r8) :: ARM                !q10 for maintenance respiration
    REAL(r8) :: FOLNMX             !foliage nitrogen concentration when f(n)=1 (%)
    REAL(r8) :: TMIN               !minimum temperature for photosynthesis (k)
       
    REAL(r8) :: XL                 !leaf/stem orientation index
    REAL(r8) :: RHOL(MBAND)        !leaf reflectance: 1=vis, 2=nir
    REAL(r8) :: RHOS(MBAND)        !stem reflectance: 1=vis, 2=nir
    REAL(r8) :: TAUL(MBAND)        !leaf transmittance: 1=vis, 2=nir
    REAL(r8) :: TAUS(MBAND)        !stem transmittance: 1=vis, 2=nir

    REAL(r8) :: MRP                !microbial respiration parameter (umol co2 /kg c/ s)
    REAL(r8) :: CWPVT              !empirical canopy wind parameter

    REAL(r8) :: WRRAT              !wood to non-wood ratio
    REAL(r8) :: WDPOOL             !wood pool (switch 1 or 0) depending on woody or not [-]
    REAL(r8) :: TDLEF              !characteristic T for leaf freezing [K]

  INTEGER :: NROOT              !number of soil layers with root present
     REAL(r8) :: RGL                !Parameter used in radiation stress function
     REAL(r8) :: RSMIN              !Minimum stomatal resistance [s m-1]
     REAL(r8) :: HS                 !Parameter used in vapor pressure deficit function
     REAL(r8) :: TOPT               !Optimum transpiration air temperature [K]
     REAL(r8) :: RSMAX              !Maximal stomatal resistance [s m-1]

     REAL(r8) :: SLAREA
     REAL(r8) :: EPS(5)

!------------------------------------------------------------------------------------------!
! From the rad section of MPTABLE.TBL
!------------------------------------------------------------------------------------------!

     REAL(r8) :: ALBSAT(MBAND)       !saturated soil albedos: 1=vis, 2=nir
     REAL(r8) :: ALBDRY(MBAND)       !dry soil albedos: 1=vis, 2=nir
     REAL(r8) :: ALBICE(MBAND)       !albedo land ice: 1=vis, 2=nir
     REAL(r8) :: ALBLAK(MBAND)       !albedo frozen lakes: 1=vis, 2=nir
     REAL(r8) :: OMEGAS(MBAND)       !two-stream parameter omega for snow
     REAL(r8) :: BETADS              !two-stream parameter betad for snow
     REAL(r8) :: BETAIS              !two-stream parameter betad for snow
     REAL(r8) :: EG(2)               !emissivity

!------------------------------------------------------------------------------------------!
! From the globals section of MPTABLE.TBL
!------------------------------------------------------------------------------------------!
 
     REAL(r8) :: CO2          !co2 partial pressure
     REAL(r8) :: O2           !o2 partial pressure
     REAL(r8) :: TIMEAN       !gridcell mean topgraphic index (global mean)
     REAL(r8) :: FSATMX       !maximum surface saturated fraction (global mean)
     REAL(r8) :: Z0SNO        !snow surface roughness length (m) (0.002)
     REAL(r8) :: SSI          !liquid water holding capacity for snowpack (m3/m3)
     REAL(r8) :: SWEMX        !new snow mass to fully cover old snow (mm)
     REAL(r8) :: RSURF_SNOW   !surface resistance for snow(s/m)

!------------------------------------------------------------------------------------------!
! From the crop section of MPTABLE.TBL
!------------------------------------------------------------------------------------------!
 
  INTEGER :: PLTDAY           ! Planting date
  INTEGER :: HSDAY            ! Harvest date
     REAL(r8) :: PLANTPOP         ! Plant density [per ha] - used?
     REAL(r8) :: IRRI             ! Irrigation strategy 0= non-irrigation 1=irrigation (no water-stress)
     REAL(r8) :: GDDTBASE         ! Base temperature for GDD accumulation [C]
     REAL(r8) :: GDDTCUT          ! Upper temperature for GDD accumulation [C]
     REAL(r8) :: GDDS1            ! GDD from seeding to emergence
     REAL(r8) :: GDDS2            ! GDD from seeding to initial vegetative 
     REAL(r8) :: GDDS3            ! GDD from seeding to post vegetative 
     REAL(r8) :: GDDS4            ! GDD from seeding to intial reproductive
     REAL(r8) :: GDDS5            ! GDD from seeding to pysical maturity 
  INTEGER :: C3C4             ! photosynthetic pathway:  1 = c3 2 = c4
     REAL(r8) :: AREF             ! reference maximum CO2 assimulation rate 
     REAL(r8) :: PSNRF            ! CO2 assimulation reduction factor(0-1) (caused by non-modeling part,e.g.pest,weeds)
     REAL(r8) :: I2PAR            ! Fraction of incoming solar radiation to photosynthetically active radiation
     REAL(r8) :: TASSIM0          ! Minimum temperature for CO2 assimulation [C]
     REAL(r8) :: TASSIM1          ! CO2 assimulation linearly increasing until temperature reaches T1 [C]
     REAL(r8) :: TASSIM2          ! CO2 assmilation rate remain at Aref until temperature reaches T2 [C]
     REAL(r8) :: K                ! light extinction coefficient
     REAL(r8) :: EPSI             ! initial light use efficiency
     REAL(r8) :: Q10MR            ! q10 for maintainance respiration
     REAL(r8) :: FOLN_MX          ! foliage nitrogen concentration when f(n)=1 (%)
     REAL(r8) :: LEFREEZ          ! characteristic T for leaf freezing [K]
     REAL(r8) :: DILE_FC(NSTAGE)  ! coeficient for temperature leaf stress death [1/s]
     REAL(r8) :: DILE_FW(NSTAGE)  ! coeficient for water leaf stress death [1/s]
     REAL(r8) :: FRA_GR           ! fraction of growth respiration 
     REAL(r8) :: LF_OVRC(NSTAGE)  ! fraction of leaf turnover  [1/s]
     REAL(r8) :: ST_OVRC(NSTAGE)  ! fraction of stem turnover  [1/s]
     REAL(r8) :: RT_OVRC(NSTAGE)  ! fraction of root tunrover  [1/s]
     REAL(r8) :: LFMR25           ! leaf maintenance respiration at 25C [umol CO2/m**2  /s]
     REAL(r8) :: STMR25           ! stem maintenance respiration at 25C [umol CO2/kg bio/s]
     REAL(r8) :: RTMR25           ! root maintenance respiration at 25C [umol CO2/kg bio/s]
     REAL(r8) :: GRAINMR25        ! grain maintenance respiration at 25C [umol CO2/kg bio/s]
     REAL(r8) :: LFPT(NSTAGE)     ! fraction of carbohydrate flux to leaf
     REAL(r8) :: STPT(NSTAGE)     ! fraction of carbohydrate flux to stem
     REAL(r8) :: RTPT(NSTAGE)     ! fraction of carbohydrate flux to root
     REAL(r8) :: GRAINPT(NSTAGE)  ! fraction of carbohydrate flux to grain
     REAL(r8) :: BIO2LAI          ! leaf are per living leaf biomass [m^2/kg]

!------------------------------------------------------------------------------------------!
! From the SOILPARM.TBL tables, as functions of soil category.
!------------------------------------------------------------------------------------------!
     REAL(r8) :: BEXP(NSOIL)   !B parameter
     REAL(r8) :: SMCDRY(NSOIL) !dry soil moisture threshold where direct evap from top
                           !layer ends (volumetric) (not used MB: 20140718)
     REAL(r8) :: SMCWLT(NSOIL) !wilting point soil moisture (volumetric)
     REAL(r8) :: SMCREF(NSOIL) !reference soil moisture (field capacity) (volumetric)
     REAL(r8) :: SMCMAX(NSOIL) !porosity, saturated value of soil moisture (volumetric)
     REAL(r8) :: PSISAT(NSOIL) !saturated soil matric potential
     REAL(r8) :: DKSAT(NSOIL)  !saturated soil hydraulic conductivity
     REAL(r8) :: DWSAT(NSOIL)  !saturated soil hydraulic diffusivity
     REAL(r8) :: QUARTZ(NSOIL) !soil quartz content
     REAL(r8) :: F1            !soil thermal diffusivity/conductivity coef (not used MB: 20140718)
!------------------------------------------------------------------------------------------!
! From the GENPARM.TBL file
!------------------------------------------------------------------------------------------!
     REAL(r8) :: SLOPE       !slope index (0 - 1)
     REAL(r8) :: CSOIL       !vol. soil heat capacity [j/m3/K]
     REAL(r8) :: ZBOT        !Depth (m) of lower boundary soil temperature
     REAL(r8) :: CZIL        !Calculate roughness length of heat
     REAL(r8) :: REFDK
     REAL(r8) :: REFKDT

     REAL(r8) :: KDT         !used in compute maximum infiltration rate (in INFIL)
     REAL(r8) :: FRZX        !used in compute maximum infiltration rate (in INFIL)

  END TYPE noahmp_parameters

contains
!
!== begin noahmp_sflx ==============================================================================

  SUBROUTINE NOAHMP_SFLX (parameters, &
                   ILOC    , JLOC    , LAT     , YEARLEN , JULIAN  , COSZ    , & ! IN : Time/Space-related
                   DT      , DX      , DZ8W    , NSOIL   , ZSOIL   , NSNOW   , & ! IN : Model configuration 
                   SHDFAC  , SHDMAX  , VEGTYP  , ICE     , IST     ,           & ! IN : Vegetation/Soil characteristics
                   SMCEQ   ,                                                   & ! IN : Vegetation/Soil characteristics
                   SFCTMP  , SFCPRS  , PSFC    , UU      , VV      , Q2      , & ! IN : Forcing
                   QC      , SOLDN   , LWDN    ,                               & ! IN : Forcing
	           PRCPCONV, PRCPNONC, PRCPSHCV, PRCPSNOW, PRCPGRPL, PRCPHAIL, & ! IN : Forcing
                   TBOT    , CO2AIR  , O2AIR   , FOLN    , FICEOLD , ZLVL    , & ! IN : Forcing
                   ALBOLD  , SNEQVO  ,                                         & ! IN/OUT : 
                   STC     , SH2O    , SMC     , TAH     , EAH     , FWET    , & ! IN/OUT : 
                   CANLIQ  , CANICE  , TV      , TG      , QSFC    , QSNOW   , & ! IN/OUT : 
                   ISNOW   , ZSNSO   , SNOWH   , SNEQV   , SNICE   , SNLIQ   , & ! IN/OUT : 
                   ZWT     , WA      , WT      , WSLAKE  , LFMASS  , RTMASS  , & ! IN/OUT : 
                   STMASS  , WOOD    , STBLCP  , FASTCP  , LAI     , SAI     , & ! IN/OUT : 
                   CM      , CH      , TAUSS   ,                               & ! IN/OUT : 
                   GRAIN   , GDD     ,                                         & ! IN/OUT 
                   SMCWTD  ,DEEPRECH , RECH    ,                               & ! IN/OUT :
		   Z0WRF   , &
                   FSA     , FSR     , FIRA    , FSH     , SSOIL   , FCEV    , & ! OUT : 
                   FGEV    , FCTR    , ECAN    , ETRAN   , EDIR    , TRAD    , & ! OUT :
                   TGB     , TGV     , T2MV    , T2MB    , Q2V     , Q2B     , & ! OUT :
                   RUNSRF  , RUNSUB  , APAR    , PSN     , SAV     , SAG     , & ! OUT :
                   FSNO    , NEE     , GPP     , NPP     , FVEG    , ALBEDO  , & ! OUT :
                   QSNBOT  , PONDING , PONDING1, PONDING2, RSSUN   , RSSHA   , & ! OUT :
                   BGAP    , WGAP    , CHV     , CHB     , EMISSI  ,           & ! OUT :
		   SHG     , SHC     , SHB     , EVG     , EVB     , GHV     , & ! OUT :
		   GHB     , IRG     , IRC     , IRB     , TR      , EVC     , & ! OUT :
		   CHLEAF  , CHUC    , CHV2    , CHB2    , FPICE   , PAHV    , &
                   PAHG    , PAHB    , PAH,                             &
#ifdef WRF_HYDRO
                   ,SFCHEADRT,                                         & ! IN/OUT :
#endif
           TAUX, TAUY, FIRE, ALBD_OUT, ALBI_OUT        )
!<-
! --------------------------------------------------------------------------------------------------
! Initial code: Guo-Yue Niu, Oct. 2007
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), INTENT(IN) :: parameters

  INTEGER                        , INTENT(IN)    :: ICE    !ice (ice = 1)
  INTEGER                        , INTENT(IN)    :: IST    !surface type 1->soil; 2->lake
  INTEGER                        , INTENT(IN)    :: VEGTYP !vegetation type 
  INTEGER                        , INTENT(IN)    :: NSNOW  !maximum no. of snow layers        
  INTEGER                        , INTENT(IN)    :: NSOIL  !no. of soil layers        
  INTEGER                        , INTENT(IN)    :: ILOC   !grid index
  INTEGER                        , INTENT(IN)    :: JLOC   !grid index
  REAL(r8)                           , INTENT(IN)    :: DT     !time step [sec]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL  !layer-bottom depth from soil surf (m)
  REAL(r8)                           , INTENT(IN)    :: Q2     !mixing ratio (kg/kg) lowest model layer
  REAL(r8)                           , INTENT(IN)    :: SFCTMP !surface air temperature [K]
  REAL(r8)                           , INTENT(IN)    :: UU     !wind speed in eastward dir (m/s)
  REAL(r8)                           , INTENT(IN)    :: VV     !wind speed in northward dir (m/s)
  REAL(r8)                           , INTENT(IN)    :: SOLDN  !downward shortwave radiation (w/m2)
  REAL(r8)                           , INTENT(IN)    :: LWDN   !downward longwave radiation (w/m2)
  REAL(r8)                           , INTENT(IN)    :: SFCPRS !pressure (pa)
  REAL(r8)                           , INTENT(INOUT) :: ZLVL   !reference height (m)
  REAL(r8)                           , INTENT(IN)    :: COSZ   !cosine solar zenith angle [0-1]
  REAL(r8)                           , INTENT(IN)    :: TBOT   !bottom condition for soil temp. [K]
  REAL(r8)                           , INTENT(IN)    :: FOLN   !foliage nitrogen (%) [1-saturated]
  REAL(r8)                           , INTENT(IN)    :: SHDFAC !green vegetation fraction [0.0-1.0]
  INTEGER                        , INTENT(IN)    :: YEARLEN!Number of days in the particular year.
  REAL(r8)                           , INTENT(IN)    :: JULIAN !Julian day of year (floating point)
  REAL(r8)                           , INTENT(IN)    :: LAT    !latitude (radians)
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: FICEOLD!ice fraction at last timestep
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)    :: SMCEQ  !equilibrium soil water  content [m3/m3]
  REAL(r8)                           , INTENT(IN)    :: PRCPCONV ! convective precipitation entering  [mm/s]    ! MB/AN : v3.7
  REAL(r8)                           , INTENT(IN)    :: PRCPNONC ! non-convective precipitation entering [mm/s] ! MB/AN : v3.7
  REAL(r8)                           , INTENT(IN)    :: PRCPSHCV ! shallow convective precip entering  [mm/s]   ! MB/AN : v3.7
  REAL(r8)                           , INTENT(IN)    :: PRCPSNOW ! snow entering land model [mm/s]              ! MB/AN : v3.7
  REAL(r8)                           , INTENT(IN)    :: PRCPGRPL ! graupel entering land model [mm/s]           ! MB/AN : v3.7
  REAL(r8)                           , INTENT(IN)    :: PRCPHAIL ! hail entering land model [mm/s]              ! MB/AN : v3.7

!jref:start; in 
  REAL(r8)                           , INTENT(IN)    :: QC     !cloud water mixing ratio
  REAL(r8)                           , INTENT(INOUT)    :: QSFC   !mixing ratio at lowest model layer
  REAL(r8)                           , INTENT(IN)    :: PSFC   !pressure at lowest model layer
  REAL(r8)                           , INTENT(IN)    :: DZ8W   !thickness of lowest layer
  REAL(r8)                           , INTENT(IN)    :: DX
  REAL(r8)                           , INTENT(IN)    :: SHDMAX  !yearly max vegetation fraction
!jref:end

#ifdef WRF_HYDRO
  REAL(r8)                           , INTENT(INOUT)    :: sfcheadrt
#endif

! input/output : need arbitary intial values
  REAL(r8)                           , INTENT(INOUT) :: QSNOW  !snowfall [mm/s]
  REAL(r8)                           , INTENT(INOUT) :: FWET   !wetted or snowed fraction of canopy (-)
  REAL(r8)                           , INTENT(INOUT) :: SNEQVO !snow mass at last time step (mm)
  REAL(r8)                           , INTENT(INOUT) :: EAH    !canopy air vapor pressure (pa)
  REAL(r8)                           , INTENT(INOUT) :: TAH    !canopy air tmeperature (k)
  REAL(r8)                           , INTENT(INOUT) :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(r8)                           , INTENT(INOUT) :: CM     !momentum drag coefficient
  REAL(r8)                           , INTENT(INOUT) :: CH     !sensible heat exchange coefficient
  REAL(r8)                           , INTENT(INOUT) :: TAUSS  !non-dimensional snow age

! prognostic variables
  INTEGER                        , INTENT(INOUT) :: ISNOW  !actual no. of snow layers [-]
  REAL(r8)                           , INTENT(INOUT) :: CANLIQ !intercepted liquid water (mm)
  REAL(r8)                           , INTENT(INOUT) :: CANICE !intercepted ice mass (mm)
  REAL(r8)                           , INTENT(INOUT) :: SNEQV  !snow water eqv. [mm]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SMC    !soil moisture (ice + liq.) [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO  !layer-bottom depth from snow surf [m]
  REAL(r8)                           , INTENT(INOUT) :: SNOWH  !snow height [m]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  !snow layer ice [mm]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  !snow layer liquid water [mm]
  REAL(r8)                           , INTENT(INOUT) :: TV     !vegetation temperature (k)
  REAL(r8)                           , INTENT(INOUT) :: TG     !ground temperature (k)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    !snow/soil temperature [k]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   !liquid soil moisture [m3/m3]
  REAL(r8)                           , INTENT(INOUT) :: ZWT    !depth to water table [m]
  REAL(r8)                           , INTENT(INOUT) :: WA     !water storage in aquifer [mm]
  REAL(r8)                           , INTENT(INOUT) :: WT     !water in aquifer&saturated soil [mm]
  REAL(r8)                           , INTENT(INOUT) :: WSLAKE !lake water storage (can be neg.) (mm)
  REAL(r8),                            INTENT(INOUT) :: SMCWTD !soil water content between bottom of the soil and water table [m3/m3]
  REAL(r8),                            INTENT(INOUT) :: DEEPRECH !recharge to or from the water table when deep [m]
  REAL(r8),                            INTENT(INOUT) :: RECH !recharge to or from the water table when shallow [m] (diagnostic)

! output
  REAL(r8)                           , INTENT(OUT)   :: Z0WRF  !combined z0 sent to coupled model
  REAL(r8)                           , INTENT(OUT)   :: FSA    !total absorbed solar radiation (w/m2)
  REAL(r8)                           , INTENT(OUT)   :: FSR    !total reflected solar radiation (w/m2)
  REAL(r8)                           , INTENT(OUT)   :: FIRA   !total net LW rad (w/m2)  [+ to atm]
  REAL(r8)                           , INTENT(OUT)   :: FSH    !total sensible heat (w/m2) [+ to atm]
  REAL(r8)                           , INTENT(OUT)   :: FCEV   !canopy evap heat (w/m2) [+ to atm]
  REAL(r8)                           , INTENT(OUT)   :: FGEV   !ground evap heat (w/m2) [+ to atm]
  REAL(r8)                           , INTENT(OUT)   :: FCTR   !transpiration heat (w/m2) [+ to atm]
  REAL(r8)                           , INTENT(OUT)   :: SSOIL  !ground heat flux (w/m2)   [+ to soil]
  REAL(r8)                           , INTENT(OUT)   :: TRAD   !surface radiative temperature (k)
  REAL(r8)                                           :: TS     !surface temperature (k)
  REAL(r8)                           , INTENT(OUT)   :: ECAN   !evaporation of intercepted water (mm/s)
  REAL(r8)                           , INTENT(OUT)   :: ETRAN  !transpiration rate (mm/s)
  REAL(r8)                           , INTENT(OUT)   :: EDIR   !soil surface evaporation rate (mm/s]
  REAL(r8)                           , INTENT(OUT)   :: RUNSRF !surface runoff [mm/s] 
  REAL(r8)                           , INTENT(OUT)   :: RUNSUB !baseflow (saturation excess) [mm/s]
  REAL(r8)                           , INTENT(OUT)   :: PSN    !total photosynthesis (umol co2/m2/s) [+]
  REAL(r8)                           , INTENT(OUT)   :: APAR   !photosyn active energy by canopy (w/m2)
  REAL(r8)                           , INTENT(OUT)   :: SAV    !solar rad absorbed by veg. (w/m2)
  REAL(r8)                           , INTENT(OUT)   :: SAG    !solar rad absorbed by ground (w/m2)
  REAL(r8)                           , INTENT(OUT)   :: FSNO   !snow cover fraction on the ground (-)
  REAL(r8)                           , INTENT(OUT)   :: FVEG   !green vegetation fraction [0.0-1.0]
  REAL(r8)                           , INTENT(OUT)   :: ALBEDO !surface albedo [-]
  REAL(r8)                                           :: ERRWAT !water error [kg m{-2}]
  REAL(r8)                           , INTENT(OUT)   :: QSNBOT !snowmelt out bottom of pack [mm/s]
  REAL(r8)                           , INTENT(OUT)   :: PONDING!surface ponding [mm]
  REAL(r8)                           , INTENT(OUT)   :: PONDING1!surface ponding [mm]
  REAL(r8)                           , INTENT(OUT)   :: PONDING2!surface ponding [mm]

!jref:start; output
  REAL(r8)                           , INTENT(OUT)     :: T2MV   !2-m air temperature over vegetated part [k]
  REAL(r8)                           , INTENT(OUT)     :: T2MB   !2-m air temperature over bare ground part [k]
  REAL(r8), INTENT(OUT) :: RSSUN        !sunlit leaf stomatal resistance (s/m)
  REAL(r8), INTENT(OUT) :: RSSHA        !shaded leaf stomatal resistance (s/m)
  REAL(r8), INTENT(OUT) :: BGAP
  REAL(r8), INTENT(OUT) :: WGAP
  REAL(r8), INTENT(OUT) :: TGV
  REAL(r8), INTENT(OUT) :: TGB
  REAL(r8)              :: Q1
  REAL(r8), INTENT(OUT) :: EMISSI
!jref:end

! local
  INTEGER                                        :: IZ     !do-loop index
  INTEGER, DIMENSION(-NSNOW+1:NSOIL)             :: IMELT  !phase change index [1-melt; 2-freeze]
  REAL(r8)                                           :: CMC    !intercepted water (CANICE+CANLIQ) (mm)
!  REAL(r8)                                           :: TAUX   !wind stress: e-w (n/m2)
!  REAL(r8)                                           :: TAUY   !wind stress: n-s (n/m2)
  REAL(r8),intent(out)                                :: TAUX   !wind stress: e-w (n/m2)
  REAL(r8),intent(out)                                :: TAUY   !wind stress: n-s (n/m2)
  REAL(r8),intent(out)                                :: FIRE  ! cheyz for coupling  
  REAL(r8), DIMENSION(       1:    2)                :: ALBD_OUT,ALBI_OUT  !cheyz for coupling              
  REAL(r8)                                           :: RHOAIR !density air (kg/m3)
!  REAL(r8), DIMENSION(       1:    5)                :: VOCFLX !voc fluxes [ug C m-2 h-1]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL)                :: DZSNSO !snow/soil layer thickness [m]
  REAL(r8)                                           :: THAIR  !potential temperature (k)
  REAL(r8)                                           :: QAIR   !specific humidity (kg/kg) (q2/(1+q2))
  REAL(r8)                                           :: EAIR   !vapor pressure air (pa)
  REAL(r8), DIMENSION(       1:    2)                :: SOLAD  !incoming direct solar rad (w/m2)
  REAL(r8), DIMENSION(       1:    2)                :: SOLAI  !incoming diffuse solar rad (w/m2)
  REAL(r8)                                           :: QPRECC !convective precipitation (mm/s)
  REAL(r8)                                           :: QPRECL !large-scale precipitation (mm/s)
  REAL(r8)                                           :: IGS    !growing season index (0=off, 1=on)
  REAL(r8)                                           :: ELAI   !leaf area index, after burying by snow
  REAL(r8)                                           :: ESAI   !stem area index, after burying by snow
  REAL(r8)                                           :: BEVAP  !soil water evaporation factor (0 - 1)
  REAL(r8), DIMENSION(       1:NSOIL)                :: BTRANI !Soil water transpiration factor (0 - 1)
  REAL(r8)                                           :: BTRAN  !soil water transpiration factor (0 - 1)
  REAL(r8)                                           :: QIN    !groundwater recharge [mm/s]
  REAL(r8)                                           :: QDIS   !groundwater discharge [mm/s]
  REAL(r8), DIMENSION(       1:NSOIL)                :: SICE   !soil ice content (m3/m3)
  REAL(r8), DIMENSION(-NSNOW+1:    0)                :: SNICEV !partial volume ice of snow [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:    0)                :: SNLIQV !partial volume liq of snow [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:    0)                :: EPORE  !effective porosity [m3/m3]
  REAL(r8)                                           :: TOTSC  !total soil carbon (g/m2)
  REAL(r8)                                           :: TOTLB  !total living carbon (g/m2)
  REAL(r8)                                           :: T2M    !2-meter air temperature (k)
  REAL(r8)                                           :: QDEW   !ground surface dew rate [mm/s]
  REAL(r8)                                           :: QVAP   !ground surface evap. rate [mm/s]
  REAL(r8)                                           :: LATHEA !latent heat [j/kg]
  REAL(r8)                                           :: SWDOWN !downward solar [w/m2]
  REAL(r8)                                           :: QMELT  !snowmelt [mm/s]
  REAL(r8)                                           :: BEG_WB !water storage at begin of a step [mm]
  REAL(r8),INTENT(OUT)                                              :: IRC    !canopy net LW rad. [w/m2] [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: IRG    !ground net LW rad. [w/m2] [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: SHC    !canopy sen. heat [w/m2]   [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: SHG    !ground sen. heat [w/m2]   [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: EVG    !ground evap. heat [w/m2]  [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: GHV    !ground heat flux [w/m2]  [+ to soil]
  REAL(r8),INTENT(OUT)                                              :: IRB    !net longwave rad. [w/m2] [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: SHB    !sensible heat [w/m2]     [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: EVB    !evaporation heat [w/m2]  [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: GHB    !ground heat flux [w/m2] [+ to soil]
  REAL(r8),INTENT(OUT)                                              :: EVC    !canopy evap. heat [w/m2]  [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: TR     !transpiration heat [w/m2] [+ to atm]
  REAL(r8), INTENT(OUT)   :: FPICE   !snow fraction in precipitation
  REAL(r8), INTENT(OUT)   :: PAHV    !precipitation advected heat - vegetation net (W/m2)
  REAL(r8), INTENT(OUT)   :: PAHG    !precipitation advected heat - under canopy net (W/m2)
  REAL(r8), INTENT(OUT)   :: PAHB    !precipitation advected heat - bare ground net (W/m2)
  REAL(r8), INTENT(OUT)                                           :: PAH     !precipitation advected heat - total (W/m2)

!jref:start 
  REAL(r8)                                           :: FSRV
  REAL(r8)                                           :: FSRG
  REAL(r8),INTENT(OUT)                               :: Q2V
  REAL(r8),INTENT(OUT)                               :: Q2B
  REAL(r8) :: Q2E
  REAL(r8) :: QFX
  REAL(r8),INTENT(OUT)                               :: CHV    !sensible heat exchange coefficient over vegetated fraction
  REAL(r8),INTENT(OUT)                               :: CHB    !sensible heat exchange coefficient over bare-ground
  REAL(r8),INTENT(OUT)                               :: CHLEAF !leaf exchange coefficient
  REAL(r8),INTENT(OUT)                               :: CHUC   !under canopy exchange coefficient
  REAL(r8),INTENT(OUT)                               :: CHV2    !sensible heat exchange coefficient over vegetated fraction
  REAL(r8),INTENT(OUT)                               :: CHB2    !sensible heat exchange coefficient over bare-ground
!jref:end  

! carbon
! inputs
  REAL(r8)                           , INTENT(IN)    :: CO2AIR !atmospheric co2 concentration (pa)
  REAL(r8)                           , INTENT(IN)    :: O2AIR  !atmospheric o2 concentration (pa)

! inputs and outputs : prognostic variables
  REAL(r8)                        , INTENT(INOUT)    :: LFMASS !leaf mass [g/m2]
  REAL(r8)                        , INTENT(INOUT)    :: RTMASS !mass of fine roots [g/m2]
  REAL(r8)                        , INTENT(INOUT)    :: STMASS !stem mass [g/m2]
  REAL(r8)                        , INTENT(INOUT)    :: WOOD   !mass of wood (incl. woody roots) [g/m2]
  REAL(r8)                        , INTENT(INOUT)    :: STBLCP !stable carbon in deep soil [g/m2]
  REAL(r8)                        , INTENT(INOUT)    :: FASTCP !short-lived carbon, shallow soil [g/m2]
  REAL(r8)                        , INTENT(INOUT)    :: LAI    !leaf area index [-]
  REAL(r8)                        , INTENT(INOUT)    :: SAI    !stem area index [-]
  REAL(r8)                        , INTENT(INOUT)    :: GRAIN  !grain mass [g/m2]
  REAL(r8)                        , INTENT(INOUT)    :: GDD    !growing degree days

! outputs
  REAL(r8)                          , INTENT(OUT)    :: NEE    !net ecosys exchange (g/m2/s CO2)
  REAL(r8)                          , INTENT(OUT)    :: GPP    !net instantaneous assimilation [g/m2/s C]
  REAL(r8)                          , INTENT(OUT)    :: NPP    !net primary productivity [g/m2/s C]
  REAL(r8)                                           :: AUTORS !net ecosystem respiration (g/m2/s C)
  REAL(r8)                                           :: HETERS !organic respiration (g/m2/s C)
  REAL(r8)                                           :: TROOT  !root-zone averaged temperature (k)
  REAL(r8)                                           :: BDFALL   !bulk density of new snow (kg/m3)    ! MB/AN: v3.7
  REAL(r8)                                           :: RAIN     !rain rate                   (mm/s)  ! MB/AN: v3.7
  REAL(r8)                                           :: SNOW     !liquid equivalent snow rate (mm/s)  ! MB/AN: v3.7
  REAL(r8)                                           :: FP                                            ! MB/AN: v3.7
  REAL(r8)                                           :: PRCP                                          ! MB/AN: v3.7
!more local variables for precip heat MB
  REAL(r8)                                           :: QINTR   !interception rate for rain (mm/s)
  REAL(r8)                                           :: QDRIPR  !drip rate for rain (mm/s)
  REAL(r8)                                           :: QTHROR  !throughfall for rain (mm/s)
  REAL(r8)                                           :: QINTS   !interception (loading) rate for snowfall (mm/s)
  REAL(r8)                                           :: QDRIPS  !drip (unloading) rate for intercepted snow (mm/s)
  REAL(r8)                                           :: QTHROS  !throughfall of snowfall (mm/s)
  REAL(r8)                                           :: QRAIN   !rain at ground srf (mm/s) [+]
  REAL(r8)                                           :: SNOWHIN !snow depth increasing rate (m/s)
  REAL(r8)                                 :: LATHEAV !latent heat vap./sublimation (j/kg)
  REAL(r8)                                 :: LATHEAG !latent heat vap./sublimation (j/kg)
  LOGICAL                             :: FROZEN_GROUND ! used to define latent heat pathway
  LOGICAL                             :: FROZEN_CANOPY ! used to define latent heat pathway
  
  ! INTENT (OUT) variables need to be assigned a value.  These normally get assigned values
  ! only if DVEG == 2.
  nee = 0.0
  npp = 0.0
  gpp = 0.0
      PAHV  = 0.
      PAHG  = 0.
      PAHB  = 0.
      PAH  = 0.

! --------------------------------------------------------------------------------------------------
! re-process atmospheric forcing

   CALL ATM (parameters,SFCPRS  ,SFCTMP   ,Q2      ,                            &
             PRCPCONV, PRCPNONC,PRCPSHCV,PRCPSNOW,PRCPGRPL,PRCPHAIL, &
             SOLDN   ,COSZ     ,THAIR   ,QAIR    ,                   & 
             EAIR    ,RHOAIR   ,QPRECC  ,QPRECL  ,SOLAD   ,SOLAI   , &
             SWDOWN  ,BDFALL   ,RAIN    ,SNOW    ,FP      ,FPICE   , PRCP )     

! snow/soil layer thickness (m)

     DO IZ = ISNOW+1, NSOIL
         IF(IZ == ISNOW+1) THEN
           DZSNSO(IZ) = - ZSNSO(IZ)
         ELSE
           DZSNSO(IZ) = ZSNSO(IZ-1) - ZSNSO(IZ)
         END IF
     END DO

! root-zone temperature

     TROOT  = 0.
     DO IZ=1,parameters%NROOT
        TROOT = TROOT + STC(IZ)*DZSNSO(IZ)/(-ZSOIL(parameters%NROOT))
     ENDDO

! total water storage for water balance check
    
     IF(IST == 1) THEN
     BEG_WB = CANLIQ + CANICE + SNEQV + WA
     DO IZ = 1,NSOIL
        BEG_WB = BEG_WB + SMC(IZ) * DZSNSO(IZ) * 1000.
     END DO
     END IF

! vegetation phenology

     CALL PHENOLOGY (parameters,VEGTYP , SNOWH  , TV     , LAT   , YEARLEN , JULIAN , & !in
                     LAI    , SAI    , TROOT  , ELAI    , ESAI   ,IGS)

!input GVF should be consistent with LAI
     IF(DVEG == 1 .or. DVEG == 6 .or. DVEG == 7) THEN
        FVEG = SHDFAC
        IF(FVEG <= 0.05) FVEG = 0.05
     ELSE IF (DVEG == 2 .or. DVEG == 3 .or. DVEG == 8) THEN
        FVEG = 1.-EXP(-0.52*(LAI+SAI))
        IF(FVEG <= 0.05) FVEG = 0.05
     ELSE IF (DVEG == 4 .or. DVEG == 5 .or. DVEG == 9 .or. DVEG == 10) THEN
        FVEG = SHDMAX
        IF(FVEG <= 0.05) FVEG = 0.05
     ELSE
        WRITE(*,*) "-------- FATAL CALLED IN SFLX -----------"
        print*, ("Namelist parameter DVEG unknown") 
     ENDIF
     IF(parameters%urban_flag .OR. VEGTYP == parameters%ISBARREN) FVEG = 0.0
     IF(ELAI+ESAI == 0.0) FVEG = 0.0

    CALL PRECIP_HEAT(parameters,ILOC   ,JLOC   ,VEGTYP ,DT     ,UU     ,VV     , & !in
                     ELAI   ,ESAI   ,FVEG   ,IST    ,                 & !in
                     BDFALL ,RAIN   ,SNOW   ,FP     ,                 & !in
                     CANLIQ ,CANICE ,TV     ,SFCTMP ,TG     ,         & !in
                     QINTR  ,QDRIPR ,QTHROR ,QINTS  ,QDRIPS ,QTHROS , & !out
                     PAHV   ,PAHG   ,PAHB   ,QRAIN  ,QSNOW  ,SNOWHIN, & !out
	             FWET   ,CMC                                    )   !out

! compute energy budget (momentum & energy fluxes and phase changes) 

    CALL ENERGY (parameters,ICE    ,VEGTYP ,IST    ,NSNOW  ,NSOIL  , & !in
                 ISNOW  ,DT     ,RHOAIR ,SFCPRS ,QAIR   , & !in
                 SFCTMP ,THAIR  ,LWDN   ,UU     ,VV     ,ZLVL   , & !in
                 CO2AIR ,O2AIR  ,SOLAD  ,SOLAI  ,COSZ   ,IGS    , & !in
                 EAIR   ,TBOT   ,ZSNSO  ,ZSOIL  , & !in
                 ELAI   ,ESAI   ,FWET   ,FOLN   ,         & !in
                 FVEG   ,PAHV   ,PAHG   ,PAHB   ,                 & !in
                 QSNOW  ,DZSNSO ,LAT    ,CANLIQ ,CANICE ,iloc, jloc , & !in
		 Z0WRF  ,                                         &
                 IMELT  ,SNICEV ,SNLIQV ,EPORE  ,T2M    ,FSNO   , & !out
                 SAV    ,SAG    ,QMELT  ,FSA    ,FSR    ,TAUX   , & !out
                 TAUY   ,FIRA   ,FSH    ,FCEV   ,FGEV   ,FCTR   , & !out
                 TRAD   ,PSN    ,APAR   ,SSOIL  ,BTRANI ,BTRAN  , & !out
                 PONDING,TS     ,LATHEAV , LATHEAG , frozen_canopy,frozen_ground,                         & !out
                 TV     ,TG     ,STC    ,SNOWH  ,EAH    ,TAH    , & !inout
                 SNEQVO ,SNEQV  ,SH2O   ,SMC    ,SNICE  ,SNLIQ  , & !inout
                 ALBOLD ,CM     ,CH     ,DX     ,DZ8W   ,Q2     , & !inout
                 TAUSS  ,                                         & !inout
!jref:start
                 QC     ,QSFC   ,PSFC   , & !in 
                 T2MV   ,T2MB  ,FSRV   , &
                 FSRG   ,RSSUN   ,RSSHA ,BGAP   ,WGAP, TGV,TGB,&
                 Q1     ,Q2V    ,Q2B    ,Q2E    ,CHV   ,CHB     , & !out
                 EMISSI ,PAH    ,                                 &
                 SHG,SHC,SHB,EVG,EVB,GHV,GHB,IRG,IRC,IRB,TR,EVC,CHLEAF,CHUC,CHV2,CHB2, &
!-> xinyf and cheyz                 
                 FIRE, ALBD_OUT, ALBI_OUT)                                            !out
!<-
!jref:end

    SICE(:) = MAX(0.0, SMC(:) - SH2O(:))   
    SNEQVO  = SNEQV

    QVAP = MAX( FGEV/LATHEAG, 0.)       ! positive part of fgev; Barlage change to ground v3.6
    QDEW = ABS( MIN(FGEV/LATHEAG, 0.))  ! negative part of fgev
    EDIR = QVAP - QDEW

! compute water budgets (water storages, ET components, and runoff)

     CALL WATER (parameters,VEGTYP ,NSNOW  ,NSOIL  ,IMELT  ,DT     ,UU     , & !in
                 VV     ,FCEV   ,FCTR   ,QPRECC ,QPRECL ,ELAI   , & !in
                 ESAI   ,SFCTMP ,QVAP   ,QDEW   ,ZSOIL  ,BTRANI , & !in
                 FICEOLD,PONDING,TG     ,IST    ,FVEG   ,iloc,jloc , SMCEQ , & !in
                 BDFALL ,FP     ,RAIN   ,SNOW   ,                 & !in  MB/AN: v3.7
		 QSNOW  ,QRAIN  ,SNOWHIN,LATHEAV,LATHEAG,frozen_canopy,frozen_ground,  & !in  MB
                 ISNOW  ,CANLIQ ,CANICE ,TV     ,SNOWH  ,SNEQV  , & !inout
                 SNICE  ,SNLIQ  ,STC    ,ZSNSO  ,SH2O   ,SMC    , & !inout
                 SICE   ,ZWT    ,WA     ,WT     ,DZSNSO ,WSLAKE , & !inout
                 SMCWTD ,DEEPRECH,RECH                          , & !inout
                 CMC    ,ECAN   ,ETRAN  ,FWET   ,RUNSRF ,RUNSUB , & !out
                 QIN    ,QDIS   ,PONDING1       ,PONDING2,&
                 QSNBOT                             &
#ifdef WRF_HYDRO
                        ,sfcheadrt                     &
#endif
                 )  !out

!     write(*,'(a20,10F15.5)') 'SFLX:RUNOFF=',RUNSRF*DT,RUNSUB*DT,EDIR*DT

! compute carbon budgets (carbon storages and co2 & bvoc fluxes)

   IF (DVEG == 2 .OR. DVEG == 5 .OR. DVEG == 6) THEN
    CALL CARBON (parameters,NSNOW  ,NSOIL  ,VEGTYP ,DT     ,ZSOIL  , & !in
                 DZSNSO ,STC    ,SMC    ,TV     ,TG     ,PSN    , & !in
                 FOLN   ,BTRAN  ,APAR   ,FVEG   ,IGS    , & !in
                 TROOT  ,IST    ,LAT    ,iloc   ,jloc   , & !in
                 LFMASS ,RTMASS ,STMASS ,WOOD   ,STBLCP ,FASTCP , & !inout
                 GPP    ,NPP    ,NEE    ,AUTORS ,HETERS ,TOTSC  , & !out
                 TOTLB  ,LAI    ,SAI    )                   !out
   END IF

   IF (DVEG == 10) THEN !XING
    CALL CARBON_CROP (parameters,NSNOW  ,NSOIL  ,VEGTYP ,DT     ,ZSOIL  ,JULIAN , & !in 
                         DZSNSO ,STC    ,SMC    ,TV     ,PSN    ,FOLN   ,BTRAN  , & !in
			 SOLDN  ,T2M    ,                                         & !in
                         LFMASS ,RTMASS ,STMASS ,WOOD   ,STBLCP ,FASTCP ,GRAIN  , & !inout
			 LAI    ,SAI    ,GDD    ,                                 & !inout
                         GPP    ,NPP    ,NEE    ,AUTORS ,HETERS ,TOTSC  ,TOTLB    ) !out
   END IF
   

! water and energy balance check

     CALL ERROR (parameters,SWDOWN ,FSA    ,FSR    ,FIRA   ,FSH    ,FCEV   , & !in
                 FGEV   ,FCTR   ,SSOIL  ,BEG_WB ,CANLIQ ,CANICE , & !in
                 SNEQV  ,WA     ,SMC    ,DZSNSO ,PRCP   ,ECAN   , & !in
                 ETRAN  ,EDIR   ,RUNSRF ,RUNSUB ,DT     ,NSOIL  , & !in
                 NSNOW  ,IST    ,ERRWAT ,ILOC   , JLOC  ,FVEG   , &
                 SAV    ,SAG    ,FSRV   ,FSRG   ,ZWT    ,PAH    , &
                 PAHV   ,PAHG   ,PAHB   )   !in ( Except ERRWAT, which is out )

! urban - jref
    QFX = ETRAN + ECAN + EDIR
    IF ( parameters%urban_flag ) THEN
       QSFC = QFX/(RHOAIR*CH) + QAIR
       Q2B = QSFC
    END IF

    IF(SNOWH <= 1.E-6 .OR. SNEQV <= 1.E-3) THEN
     SNOWH = 0.0
     SNEQV = 0.0
    END IF

    IF(SWDOWN.NE.0.) THEN
      ALBEDO = FSR / SWDOWN
    ELSE
      ALBEDO = -999.9
    END IF
    

  END SUBROUTINE NOAHMP_SFLX

!== begin atm ======================================================================================

  SUBROUTINE ATM (parameters,SFCPRS  ,SFCTMP   ,Q2      ,                             &
                  PRCPCONV,PRCPNONC ,PRCPSHCV,PRCPSNOW,PRCPGRPL,PRCPHAIL , &
                  SOLDN   ,COSZ     ,THAIR   ,QAIR    ,                    & 
                  EAIR    ,RHOAIR   ,QPRECC  ,QPRECL  ,SOLAD   , SOLAI   , &
		  SWDOWN  ,BDFALL   ,RAIN    ,SNOW    ,FP      , FPICE   ,PRCP )     
! --------------------------------------------------------------------------------------------------
! re-process atmospheric forcing
! ----------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
  REAL(r8)                          , INTENT(IN)  :: SFCPRS !pressure (pa)
  REAL(r8)                          , INTENT(IN)  :: SFCTMP !surface air temperature [k]
  REAL(r8)                          , INTENT(IN)  :: Q2     !mixing ratio (kg/kg)
  REAL(r8)                          , INTENT(IN)  :: PRCPCONV ! convective precipitation entering  [mm/s]    ! MB/AN : v3.7
  REAL(r8)                          , INTENT(IN)  :: PRCPNONC ! non-convective precipitation entering [mm/s] ! MB/AN : v3.7
  REAL(r8)                          , INTENT(IN)  :: PRCPSHCV ! shallow convective precip entering  [mm/s]   ! MB/AN : v3.7
  REAL(r8)                          , INTENT(IN)  :: PRCPSNOW ! snow entering land model [mm/s]              ! MB/AN : v3.7
  REAL(r8)                          , INTENT(IN)  :: PRCPGRPL ! graupel entering land model [mm/s]           ! MB/AN : v3.7
  REAL(r8)                          , INTENT(IN)  :: PRCPHAIL ! hail entering land model [mm/s]              ! MB/AN : v3.7
  REAL(r8)                          , INTENT(IN)  :: SOLDN  !downward shortwave radiation (w/m2)
  REAL(r8)                          , INTENT(IN)  :: COSZ   !cosine solar zenith angle [0-1]

! outputs

  REAL(r8)                          , INTENT(OUT) :: THAIR  !potential temperature (k)
  REAL(r8)                          , INTENT(OUT) :: QAIR   !specific humidity (kg/kg) (q2/(1+q2))
  REAL(r8)                          , INTENT(OUT) :: EAIR   !vapor pressure air (pa)
  REAL(r8)                          , INTENT(OUT) :: RHOAIR !density air (kg/m3)
  REAL(r8)                          , INTENT(OUT) :: QPRECC !convective precipitation (mm/s)
  REAL(r8)                          , INTENT(OUT) :: QPRECL !large-scale precipitation (mm/s)
  REAL(r8), DIMENSION(       1:   2), INTENT(OUT) :: SOLAD  !incoming direct solar radiation (w/m2)
  REAL(r8), DIMENSION(       1:   2), INTENT(OUT) :: SOLAI  !incoming diffuse solar radiation (w/m2)
  REAL(r8)                          , INTENT(OUT) :: SWDOWN !downward solar filtered by sun angle [w/m2]
  REAL(r8)                          , INTENT(OUT) :: BDFALL  !!bulk density of snowfall (kg/m3) AJN
  REAL(r8)                          , INTENT(OUT) :: RAIN    !rainfall (mm/s) AJN
  REAL(r8)                          , INTENT(OUT) :: SNOW    !liquid equivalent snowfall (mm/s) AJN
  REAL(r8)                          , INTENT(OUT) :: FP      !fraction of area receiving precipitation  AJN
  REAL(r8)                          , INTENT(OUT) :: FPICE   !fraction of ice                AJN
  REAL(r8)                          , INTENT(OUT) :: PRCP    !total precipitation [mm/s]     ! MB/AN : v3.7

!locals

  REAL(r8)                                        :: PAIR   !atm bottom level pressure (pa)
  REAL(r8)                                        :: PRCP_FROZEN   !total frozen precipitation [mm/s] ! MB/AN : v3.7
  REAL(r8), PARAMETER                             :: RHO_GRPL = 500.0  ! graupel bulk density [kg/m3] ! MB/AN : v3.7
  REAL(r8), PARAMETER                             :: RHO_HAIL = 917.0  ! hail bulk density [kg/m3]    ! MB/AN : v3.7
! --------------------------------------------------------------------------------------------------

!jref: seems like PAIR should be P1000mb??
       PAIR   = SFCPRS                   ! atm bottom level pressure (pa)
       THAIR  = SFCTMP * (SFCPRS/PAIR)**(RAIR/CPAIR) 

       QAIR   = Q2                       ! In WRF, driver converts to specific humidity

       EAIR   = QAIR*SFCPRS / (0.622+0.378*QAIR)
       RHOAIR = (SFCPRS-0.378*EAIR) / (RAIR*SFCTMP)

       IF(COSZ <= 0.) THEN 
          SWDOWN = 0.
       ELSE
          SWDOWN = SOLDN
       END IF 

       SOLAD(1) = SWDOWN*0.7*0.5     ! direct  vis
       SOLAD(2) = SWDOWN*0.7*0.5     ! direct  nir
       SOLAI(1) = SWDOWN*0.3*0.5     ! diffuse vis
       SOLAI(2) = SWDOWN*0.3*0.5     ! diffuse nir

       PRCP = PRCPCONV + PRCPNONC + PRCPSHCV

       IF(OPT_SNF == 4) THEN
         QPRECC = PRCPCONV + PRCPSHCV
	 QPRECL = PRCPNONC
       ELSE
         QPRECC = 0.10 * PRCP          ! should be from the atmospheric model
         QPRECL = 0.90 * PRCP          ! should be from the atmospheric model
       END IF

! fractional area that receives precipitation (see, Niu et al. 2005)
   
    FP = 0.0
    IF(QPRECC + QPRECL > 0.) & 
       FP = (QPRECC + QPRECL) / (10.*QPRECC + QPRECL)

! partition precipitation into rain and snow. Moved from CANWAT MB/AN: v3.7

! Jordan (1991)

     IF(OPT_SNF == 1) THEN
       IF(SFCTMP > TFRZ+2.5)THEN
           FPICE = 0.
       ELSE
         IF(SFCTMP <= TFRZ+0.5)THEN
           FPICE = 1.0
         ELSE IF(SFCTMP <= TFRZ+2.)THEN
           FPICE = 1.-(-54.632 + 0.2*SFCTMP)
         ELSE
           FPICE = 0.6
         ENDIF
       ENDIF
     ENDIF

     IF(OPT_SNF == 2) THEN
       IF(SFCTMP >= TFRZ+2.2) THEN
           FPICE = 0.
       ELSE
           FPICE = 1.0
       ENDIF
     ENDIF

     IF(OPT_SNF == 3) THEN
       IF(SFCTMP >= TFRZ) THEN
           FPICE = 0.
       ELSE
           FPICE = 1.0
       ENDIF
     ENDIF

! Hedstrom NR and JW Pomeroy (1998), Hydrol. Processes, 12, 1611-1625
! fresh snow density

     BDFALL = MIN(120.,67.92+51.25*EXP((SFCTMP-TFRZ)/2.59))       !MB/AN: change to MIN  
     IF(OPT_SNF == 4) THEN
        PRCP_FROZEN = PRCPSNOW + PRCPGRPL + PRCPHAIL
        IF(PRCPNONC > 0. .and. PRCP_FROZEN > 0.) THEN
	  FPICE = MIN(1.0,PRCP_FROZEN/PRCPNONC)
	  FPICE = MAX(0.0,FPICE)
	  BDFALL = BDFALL*(PRCPSNOW/PRCP_FROZEN) + RHO_GRPL*(PRCPGRPL/PRCP_FROZEN) + &
	             RHO_HAIL*(PRCPHAIL/PRCP_FROZEN)
	ELSE
	  FPICE = 0.0
        ENDIF
	
     ENDIF

     RAIN   = PRCP * (1.-FPICE)
     SNOW   = PRCP * FPICE


  END SUBROUTINE ATM

!== begin phenology ================================================================================

  SUBROUTINE PHENOLOGY (parameters,VEGTYP , SNOWH  , TV     , LAT   , YEARLEN , JULIAN , & !in
                        LAI    , SAI    , TROOT  , ELAI    , ESAI   , IGS)

! --------------------------------------------------------------------------------------------------
! vegetation phenology considering vegeation canopy being buries by snow and evolution in time
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER                , INTENT(IN   ) :: VEGTYP !vegetation type 
  REAL(r8)                   , INTENT(IN   ) :: SNOWH  !snow height [m]
  REAL(r8)                   , INTENT(IN   ) :: TV     !vegetation temperature (k)
  REAL(r8)                   , INTENT(IN   ) :: LAT    !latitude (radians)
  INTEGER                , INTENT(IN   ) :: YEARLEN!Number of days in the particular year
  REAL(r8)                   , INTENT(IN   ) :: JULIAN !Julian day of year (fractional) ( 0 <= JULIAN < YEARLEN )
  real                   , INTENT(IN   ) :: TROOT  !root-zone averaged temperature (k)
  REAL(r8)                   , INTENT(INOUT) :: LAI    !LAI, unadjusted for burying by snow
  REAL(r8)                   , INTENT(INOUT) :: SAI    !SAI, unadjusted for burying by snow

! outputs
  REAL(r8)                   , INTENT(OUT  ) :: ELAI   !leaf area index, after burying by snow
  REAL(r8)                   , INTENT(OUT  ) :: ESAI   !stem area index, after burying by snow
  REAL(r8)                   , INTENT(OUT  ) :: IGS    !growing season index (0=off, 1=on)

! locals

  REAL(r8)                                   :: DB     !thickness of canopy buried by snow (m)
  REAL(r8)                                   :: FB     !fraction of canopy buried by snow
  REAL(r8)                                   :: SNOWHC !critical snow depth at which short vege
                                                   !is fully covered by snow

  INTEGER                                :: K       !index
  INTEGER                                :: IT1,IT2 !interpolation months
  REAL(r8)                                   :: DAY     !current day of year ( 0 <= DAY < YEARLEN )
  REAL(r8)                                   :: WT1,WT2 !interpolation weights
  REAL(r8)                                   :: T       !current month (1.00, ..., 12.00)
! --------------------------------------------------------------------------------------------------

  IF ( DVEG == 1 .or. DVEG == 3 .or. DVEG == 4 ) THEN

     IF (LAT >= 0.) THEN
        ! Northern Hemisphere
        DAY = JULIAN
     ELSE
        ! Southern Hemisphere.  DAY is shifted by 1/2 year.
        DAY = MOD ( JULIAN + ( 0.5 * YEARLEN ) , REAL(YEARLEN) )
     ENDIF

     T = 12. * DAY / REAL(YEARLEN)
     IT1 = T + 0.5
     IT2 = IT1 + 1
     WT1 = (IT1+0.5) - T
     WT2 = 1.-WT1
     IF (IT1 .LT.  1) IT1 = 12
     IF (IT2 .GT. 12) IT2 = 1

     LAI = WT1*parameters%LAIM(IT1) + WT2*parameters%LAIM(IT2)
     SAI = WT1*parameters%SAIM(IT1) + WT2*parameters%SAIM(IT2)
  ENDIF

  IF(DVEG == 7 .or. DVEG == 8 .or. DVEG == 9) THEN
    SAI = MAX(0.05,0.1 * LAI)  ! when reading LAI, set SAI to 10% LAI, but not below 0.05 MB: v3.8
    IF (LAI < 0.05) SAI = 0.0  ! if LAI below minimum, make sure SAI = 0
  ENDIF

  IF (SAI < 0.05 .and. DVEG /= 10) SAI = 0.0                  ! MB: SAI CHECK, change to 0.05 v3.6
  IF (LAI < 0.05 .OR. SAI == 0.0 .and. DVEG /= 10) LAI = 0.0  ! MB: LAI CHECK

  IF ( ( VEGTYP == parameters%iswater ) .OR. ( VEGTYP == parameters%ISBARREN ) .OR. &
       ( VEGTYP == parameters%ISICE   ) .or. ( parameters%urban_flag ) ) THEN
     LAI  = 0.
     SAI  = 0.
  ENDIF

!buried by snow

     DB = MIN( MAX(SNOWH - parameters%HVB,0.), parameters%HVT-parameters%HVB )
     FB = DB / MAX(1.E-06,parameters%HVT-parameters%HVB)

     IF(parameters%HVT> 0. .AND. parameters%HVT <= 1.0) THEN          !MB: change to 1.0 and 0.2 to reflect
       SNOWHC = parameters%HVT*EXP(-SNOWH/0.2)             !      changes to HVT in MPTABLE
       FB     = MIN(SNOWH,SNOWHC)/SNOWHC
     ENDIF

     ELAI =  LAI*(1.-FB)
     ESAI =  SAI*(1.-FB)
     IF (ESAI < 0.05 .and. DVEG /= 10) ESAI = 0.0                   ! MB: ESAI CHECK, change to 0.05 v3.6
     IF (ELAI < 0.05 .OR. ESAI == 0.0 .and. DVEG /= 10) ELAI = 0.0  ! MB: LAI CHECK

     IF (TV .GT. parameters%TMIN) THEN
         IGS = 1.
     ELSE
         IGS = 0.
     ENDIF

  END SUBROUTINE PHENOLOGY

!== begin precip_heat ==============================================================================

  SUBROUTINE PRECIP_HEAT (parameters,ILOC   ,JLOC   ,VEGTYP ,DT     ,UU     ,VV     , & !in
                          ELAI   ,ESAI   ,FVEG   ,IST    ,                 & !in
                          BDFALL ,RAIN   ,SNOW   ,FP     ,                 & !in
                          CANLIQ ,CANICE ,TV     ,SFCTMP ,TG     ,         & !in
                          QINTR  ,QDRIPR ,QTHROR ,QINTS  ,QDRIPS ,QTHROS , & !out
			  PAHV   ,PAHG   ,PAHB   ,QRAIN  ,QSNOW  ,SNOWHIN, & !out
			  FWET   ,CMC                                    )   !out

! ------------------------ code history ------------------------------
! Michael Barlage: Oct 2013 - split CANWATER to calculate precip movement for 
!                             tracking of advected heat
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! ------------------------ input/output variables --------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,INTENT(IN)  :: ILOC    !grid index
  INTEGER,INTENT(IN)  :: JLOC    !grid index
  INTEGER,INTENT(IN)  :: VEGTYP  !vegetation type
  INTEGER,INTENT(IN)  :: IST     !surface type 1-soil; 2-lake
  REAL(r8),   INTENT(IN)  :: DT      !main time step (s)
  REAL(r8),   INTENT(IN)  :: UU      !u-direction wind speed [m/s]
  REAL(r8),   INTENT(IN)  :: VV      !v-direction wind speed [m/s]
  REAL(r8),   INTENT(IN)  :: ELAI    !leaf area index, after burying by snow
  REAL(r8),   INTENT(IN)  :: ESAI    !stem area index, after burying by snow
  REAL(r8),   INTENT(IN)  :: FVEG    !greeness vegetation fraction (-)
  REAL(r8),   INTENT(IN)  :: BDFALL  !bulk density of snowfall (kg/m3)
  REAL(r8),   INTENT(IN)  :: RAIN    !rainfall (mm/s)
  REAL(r8),   INTENT(IN)  :: SNOW    !snowfall (mm/s)
  REAL(r8),   INTENT(IN)  :: FP      !fraction of the gridcell that receives precipitation
  REAL(r8),   INTENT(IN)  :: TV      !vegetation temperature (k)
  REAL(r8),   INTENT(IN)  :: SFCTMP  !model-level temperature (k)
  REAL(r8),   INTENT(IN)  :: TG      !ground temperature (k)

! input & output
  REAL(r8), INTENT(INOUT) :: CANLIQ  !intercepted liquid water (mm)
  REAL(r8), INTENT(INOUT) :: CANICE  !intercepted ice mass (mm)

! output
  REAL(r8), INTENT(OUT)   :: QINTR   !interception rate for rain (mm/s)
  REAL(r8), INTENT(OUT)   :: QDRIPR  !drip rate for rain (mm/s)
  REAL(r8), INTENT(OUT)   :: QTHROR  !throughfall for rain (mm/s)
  REAL(r8), INTENT(OUT)   :: QINTS   !interception (loading) rate for snowfall (mm/s)
  REAL(r8), INTENT(OUT)   :: QDRIPS  !drip (unloading) rate for intercepted snow (mm/s)
  REAL(r8), INTENT(OUT)   :: QTHROS  !throughfall of snowfall (mm/s)
  REAL(r8), INTENT(OUT)   :: PAHV    !precipitation advected heat - vegetation net (W/m2)
  REAL(r8), INTENT(OUT)   :: PAHG    !precipitation advected heat - under canopy net (W/m2)
  REAL(r8), INTENT(OUT)   :: PAHB    !precipitation advected heat - bare ground net (W/m2)
  REAL(r8), INTENT(OUT)   :: QRAIN   !rain at ground srf (mm/s) [+]
  REAL(r8), INTENT(OUT)   :: QSNOW   !snow at ground srf (mm/s) [+]
  REAL(r8), INTENT(OUT)   :: SNOWHIN !snow depth increasing rate (m/s)
  REAL(r8), INTENT(OUT)   :: FWET    !wetted or snowed fraction of the canopy (-)
  REAL(r8), INTENT(OUT)   :: CMC     !intercepted water (mm)
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  REAL(r8)                :: MAXSNO  !canopy capacity for snow interception (mm)
  REAL(r8)                :: MAXLIQ  !canopy capacity for rain interception (mm)
  REAL(r8)                :: FT      !temperature factor for unloading rate
  REAL(r8)                :: FV      !wind factor for unloading rate
  REAL(r8)                :: PAH_AC  !precipitation advected heat - air to canopy (W/m2)
  REAL(r8)                :: PAH_CG  !precipitation advected heat - canopy to ground (W/m2)
  REAL(r8)                :: PAH_AG  !precipitation advected heat - air to ground (W/m2)
  REAL(r8)                :: ICEDRIP !canice unloading
! --------------------------------------------------------------------
! initialization

      QINTR   = 0.
      QDRIPR  = 0.
      QTHROR  = 0.
      QINTR   = 0.
      QINTS   = 0.
      QDRIPS  = 0.
      QTHROS  = 0.
      PAH_AC  = 0.
      PAH_CG  = 0.
      PAH_AG  = 0.
      PAHV    = 0.
      PAHG    = 0.
      PAHB    = 0.
      QRAIN   = 0.0
      QSNOW   = 0.0
      SNOWHIN = 0.0
      ICEDRIP = 0.0
!      print*, "precip_heat begin canopy balance:",canliq+canice+(rain+snow)*dt
!      print*,  "precip_heat snow*3600.0:",snow*3600.0
!      print*,  "precip_heat rain*3600.0:",rain*3600.0
!      print*,  "precip_heat canice:",canice
!      print*,  "precip_heat canliq:",canliq

! --------------------------- liquid water ------------------------------
! maximum canopy water

      MAXLIQ =  parameters%CH2OP * (ELAI+ ESAI)

! average interception and throughfall

      IF((ELAI+ ESAI).GT.0.) THEN
         QINTR  = FVEG * RAIN * FP  ! interception capability
         QINTR  = MIN(QINTR, (MAXLIQ - CANLIQ)/DT * (1.-EXP(-RAIN*DT/MAXLIQ)) )
         QINTR  = MAX(QINTR, 0.)
         QDRIPR = FVEG * RAIN - QINTR
         QTHROR = (1.-FVEG) * RAIN
         CANLIQ=MAX(0.,CANLIQ+QINTR*DT)
      ELSE
         QINTR  = 0.
         QDRIPR = 0.
         QTHROR = RAIN
	 IF(CANLIQ > 0.) THEN             ! FOR CASE OF CANOPY GETTING BURIED
	   QDRIPR = QDRIPR + CANLIQ/DT
	   CANLIQ = 0.0
	 END IF
      END IF
      
! heat transported by liquid water

      PAH_AC = FVEG * RAIN * (CWAT/1000.0) * (SFCTMP - TV)
      PAH_CG = QDRIPR * (CWAT/1000.0) * (TV - TG)
      PAH_AG = QTHROR * (CWAT/1000.0) * (SFCTMP - TG)
!      print*, "precip_heat PAH_AC:",PAH_AC
!      print*, "precip_heat PAH_CG:",PAH_CG
!      print*, "precip_heat PAH_AG:",PAH_AG

! --------------------------- canopy ice ------------------------------
! for canopy ice

      MAXSNO = 6.6*(0.27+46./BDFALL) * (ELAI+ ESAI)

      IF((ELAI+ ESAI).GT.0.) THEN
         QINTS = FVEG * SNOW * FP
         QINTS = MIN(QINTS, (MAXSNO - CANICE)/DT * (1.-EXP(-SNOW*DT/MAXSNO)) )
         QINTS = MAX(QINTS, 0.)
         FT = MAX(0.0,(TV - 270.15) / 1.87E5)
         FV = SQRT(UU*UU + VV*VV) / 1.56E5
	 ! MB: changed below to reflect the rain assumption that all precip gets intercepted 
	 ICEDRIP = MAX(0.,CANICE) * (FV+FT)    !MB: removed /DT
         QDRIPS = (FVEG * SNOW - QINTS) + ICEDRIP
         QTHROS = (1.0-FVEG) * SNOW
         CANICE= MAX(0.,CANICE + (QINTS - ICEDRIP)*DT)
      ELSE
         QINTS  = 0.
         QDRIPS = 0.
         QTHROS = SNOW
	 IF(CANICE > 0.) THEN             ! FOR CASE OF CANOPY GETTING BURIED
	   QDRIPS = QDRIPS + CANICE/DT
	   CANICE = 0.0
	 END IF
      ENDIF
!      print*, "precip_heat canopy through:",3600.0*(FVEG * SNOW - QINTS)
!      print*, "precip_heat canopy drip:",3600.0*MAX(0.,CANICE) * (FV+FT)

! wetted fraction of canopy

      IF(CANICE.GT.0.) THEN
           FWET = MAX(0.,CANICE) / MAX(MAXSNO,1.E-06)
      ELSE
           FWET = MAX(0.,CANLIQ) / MAX(MAXLIQ,1.E-06)
      ENDIF
      FWET = MIN(FWET, 1.) ** 0.667

! total canopy water

      CMC = CANLIQ + CANICE

! heat transported by snow/ice

      PAH_AC = PAH_AC +  FVEG * SNOW * (CICE/1000.0) * (SFCTMP - TV)
      PAH_CG = PAH_CG + QDRIPS * (CICE/1000.0) * (TV - TG)
      PAH_AG = PAH_AG + QTHROS * (CICE/1000.0) * (SFCTMP - TG)
      
      PAHV = PAH_AC - PAH_CG
      PAHG = PAH_CG
      PAHB = PAH_AG
      
      IF (FVEG > 0.0 .AND. FVEG < 1.0) THEN
        PAHG = PAHG / FVEG         ! these will be multiplied by fraction later
	PAHB = PAHB / (1.0-FVEG)
      ELSEIF (FVEG <= 0.0) THEN
        PAHB = PAHG + PAHB         ! for case of canopy getting buried
        PAHG = 0.0
	PAHV = 0.0
      ELSEIF (FVEG >= 1.0) THEN
	PAHB = 0.0
      END IF
      
      PAHV = MAX(PAHV,-20.0)       ! Put some artificial limits here for stability
      PAHV = MIN(PAHV,20.0)
      PAHG = MAX(PAHG,-20.0)
      PAHG = MIN(PAHG,20.0)
      PAHB = MAX(PAHB,-20.0)
      PAHB = MIN(PAHB,20.0)
      
!      print*, 'precip_heat sfctmp,tv,tg:',sfctmp,tv,tg
!      print*, 'precip_heat 3600.0*qints+qdrips+qthros:',3600.0*(qints+qdrips+qthros)
!      print*, "precip_heat maxsno:",maxsno
!      print*, "precip_heat PAH_AC:",PAH_AC
!      print*, "precip_heat PAH_CG:",PAH_CG
!      print*, "precip_heat PAH_AG:",PAH_AG
      
!      print*, "precip_heat PAHV:",PAHV
!      print*, "precip_heat PAHG:",PAHG
!      print*, "precip_heat PAHB:",PAHB
!      print*, "precip_heat fveg:",fveg
!      print*,  "precip_heat qints*3600.0:",qints*3600.0
!      print*,  "precip_heat qdrips*3600.0:",qdrips*3600.0
!      print*,  "precip_heat qthros*3600.0:",qthros*3600.0
      
! rain or snow on the ground

      QRAIN   = QDRIPR + QTHROR
      QSNOW   = QDRIPS + QTHROS
      SNOWHIN = QSNOW/BDFALL

      IF (IST == 2 .AND. TG > TFRZ) THEN
         QSNOW   = 0.
         SNOWHIN = 0.
      END IF
!      print*,  "precip_heat qsnow*3600.0:",qsnow*3600.0
!      print*,  "precip_heat qrain*3600.0:",qrain*3600.0
!      print*,  "precip_heat SNOWHIN:",SNOWHIN
!      print*,  "precip_heat canice:",canice
!      print*,  "precip_heat canliq:",canliq
!      print*, "precip_heat end canopy balance:",canliq+canice+(qrain+qsnow)*dt
      

  END SUBROUTINE PRECIP_HEAT

!== begin error ====================================================================================

  SUBROUTINE ERROR (parameters,SWDOWN ,FSA    ,FSR    ,FIRA   ,FSH    ,FCEV   , &
                    FGEV   ,FCTR   ,SSOIL  ,BEG_WB ,CANLIQ ,CANICE , &
                    SNEQV  ,WA     ,SMC    ,DZSNSO ,PRCP   ,ECAN   , &
                    ETRAN  ,EDIR   ,RUNSRF ,RUNSUB ,DT     ,NSOIL  , &
                    NSNOW  ,IST    ,ERRWAT, ILOC   ,JLOC   ,FVEG   , &
                    SAV    ,SAG    ,FSRV   ,FSRG   ,ZWT    ,PAH    , &
                    PAHV   ,PAHG   ,PAHB   )
! --------------------------------------------------------------------------------------------------
! check surface energy balance and water balance
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER                        , INTENT(IN) :: NSNOW  !maximum no. of snow layers        
  INTEGER                        , INTENT(IN) :: NSOIL  !number of soil layers
  INTEGER                        , INTENT(IN) :: IST    !surface type 1->soil; 2->lake
  INTEGER                        , INTENT(IN) :: ILOC   !grid index
  INTEGER                        , INTENT(IN) :: JLOC   !grid index
  REAL(r8)                           , INTENT(IN) :: SWDOWN !downward solar filtered by sun angle [w/m2]
  REAL(r8)                           , INTENT(IN) :: FSA    !total absorbed solar radiation (w/m2)
  REAL(r8)                           , INTENT(IN) :: FSR    !total reflected solar radiation (w/m2)
  REAL(r8)                           , INTENT(IN) :: FIRA   !total net longwave rad (w/m2)  [+ to atm]
  REAL(r8)                           , INTENT(IN) :: FSH    !total sensible heat (w/m2)     [+ to atm]
  REAL(r8)                           , INTENT(IN) :: FCEV   !canopy evaporation heat (w/m2) [+ to atm]
  REAL(r8)                           , INTENT(IN) :: FGEV   !ground evaporation heat (w/m2) [+ to atm]
  REAL(r8)                           , INTENT(IN) :: FCTR   !transpiration heat flux (w/m2) [+ to atm]
  REAL(r8)                           , INTENT(IN) :: SSOIL  !ground heat flux (w/m2)        [+ to soil]
  REAL(r8)                           , INTENT(IN) :: FVEG
  REAL(r8)                           , INTENT(IN) :: SAV
  REAL(r8)                           , INTENT(IN) :: SAG
  REAL(r8)                           , INTENT(IN) :: FSRV
  REAL(r8)                           , INTENT(IN) :: FSRG
  REAL(r8)                           , INTENT(IN) :: ZWT

  REAL(r8)                           , INTENT(IN) :: PRCP   !precipitation rate (kg m-2 s-1)
  REAL(r8)                           , INTENT(IN) :: ECAN   !evaporation of intercepted water (mm/s)
  REAL(r8)                           , INTENT(IN) :: ETRAN  !transpiration rate (mm/s)
  REAL(r8)                           , INTENT(IN) :: EDIR   !soil surface evaporation rate[mm/s]
  REAL(r8)                           , INTENT(IN) :: RUNSRF !surface runoff [mm/s] 
  REAL(r8)                           , INTENT(IN) :: RUNSUB !baseflow (saturation excess) [mm/s]
  REAL(r8)                           , INTENT(IN) :: CANLIQ !intercepted liquid water (mm)
  REAL(r8)                           , INTENT(IN) :: CANICE !intercepted ice mass (mm)
  REAL(r8)                           , INTENT(IN) :: SNEQV  !snow water eqv. [mm]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: SMC    !soil moisture (ice + liq.) [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !snow/soil layer thickness [m]
  REAL(r8)                           , INTENT(IN) :: WA     !water storage in aquifer [mm]
  REAL(r8)                           , INTENT(IN) :: DT     !time step [sec]
  REAL(r8)                           , INTENT(IN) :: BEG_WB !water storage at begin of a timesetp [mm]
  REAL(r8)                           , INTENT(OUT) :: ERRWAT !error in water balance [mm/timestep]
  REAL(r8), INTENT(IN)   :: PAH     !precipitation advected heat - total (W/m2)
  REAL(r8), INTENT(IN)   :: PAHV    !precipitation advected heat - total (W/m2)
  REAL(r8), INTENT(IN)   :: PAHG    !precipitation advected heat - total (W/m2)
  REAL(r8), INTENT(IN)   :: PAHB    !precipitation advected heat - total (W/m2)

  INTEGER                                     :: IZ     !do-loop index
  REAL(r8)                                        :: END_WB !water storage at end of a timestep [mm]
  !KWM REAL(r8)                                        :: ERRWAT !error in water balance [mm/timestep]
  REAL(r8)                                        :: ERRENG !error in surface energy balance [w/m2]
  REAL(r8)                                        :: ERRSW  !error in shortwave radiation balance [w/m2]
  REAL(r8)                                        :: FSRVG
  CHARACTER(len=256)                          :: message
! --------------------------------------------------------------------------------------------------
!jref:start
   ERRSW   = SWDOWN - (FSA + FSR)
!   ERRSW   = SWDOWN - (SAV+SAG + FSRV+FSRG)
!   WRITE(*,*) "ERRSW =",ERRSW
   IF (ABS(ERRSW) > 0.01) THEN            ! w/m2
   WRITE(*,*) "VEGETATION!"
   WRITE(*,*) "SWDOWN*FVEG =",SWDOWN*FVEG
   WRITE(*,*) "FVEG*(SAV+SAG) =",FVEG*SAV + SAG
   WRITE(*,*) "FVEG*(FSRV +FSRG)=",FVEG*FSRV + FSRG
   WRITE(*,*) "GROUND!"
   WRITE(*,*) "(1-.FVEG)*SWDOWN =",(1.-FVEG)*SWDOWN
   WRITE(*,*) "(1.-FVEG)*SAG =",(1.-FVEG)*SAG
   WRITE(*,*) "(1.-FVEG)*FSRG=",(1.-FVEG)*FSRG
   WRITE(*,*) "FSRV   =",FSRV
   WRITE(*,*) "FSRG   =",FSRG
   WRITE(*,*) "FSR    =",FSR
   WRITE(*,*) "SAV    =",SAV
   WRITE(*,*) "SAG    =",SAG
   WRITE(*,*) "FSA    =",FSA
!jref:end   
      WRITE(message,*) 'ERRSW =',ERRSW
!-> xinyf 
!      call wrf_message(trim(message))
!      call wrf_error_fatal("Stop in Noah-MP")
       print*, (trim(message))
       print*, ("stop in noah-MP")
!<-
   END IF

   ERRENG = SAV+SAG-(FIRA+FSH+FCEV+FGEV+FCTR+SSOIL) +PAH
!   ERRENG = FVEG*SAV+SAG-(FIRA+FSH+FCEV+FGEV+FCTR+SSOIL)
!   WRITE(*,*) "ERRENG =",ERRENG
   IF(ABS(ERRENG) > 0.01) THEN
      write(message,*) 'ERRENG =',ERRENG,' at i,j: ',ILOC,JLOC
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Net solar:       ",FSA
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Net longwave:    ",FIRA
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Total sensible:  ",FSH
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Canopy evap:     ",FCEV
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Ground evap:     ",FGEV
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Transpiration:   ",FCTR
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Total ground:    ",SSOIL
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,4F10.4)') "Precip advected: ",PAH,PAHV,PAHG,PAHB
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Precip: ",PRCP
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
      WRITE(message,'(a17,F10.4)') "Veg fraction: ",FVEG
!xinyf      call wrf_message(trim(message))
      print*, (trim(message))
!      call wrf_error_fatal("Energy budget problem in NOAHMP LSM")
      print*, ("Energy budget problem in NoahMP LSM")
   END IF

   IF (IST == 1) THEN                                       !soil
        END_WB = CANLIQ + CANICE + SNEQV + WA
        DO IZ = 1,NSOIL
          END_WB = END_WB + SMC(IZ) * DZSNSO(IZ) * 1000.
        END DO
        ERRWAT = END_WB-BEG_WB-(PRCP-ECAN-ETRAN-EDIR-RUNSRF-RUNSUB)*DT

#ifndef WRF_HYDRO
        IF(ABS(ERRWAT) > 0.1) THEN
           if (ERRWAT > 0) then
              print*, ('The model is gaining water (ERRWAT is positive)')
           else
              print*,('The model is losing water (ERRWAT is negative)')
           endif
           write(message, *) 'ERRWAT =',ERRWAT, "kg m{-2} timestep{-1}"
           print*,(trim(message))
           WRITE(message, &
           '("    I      J     END_WB     BEG_WB       PRCP       ECAN       EDIR      ETRAN      RUNSRF     RUNSUB")')
           print*,(trim(message))
           WRITE(message,'(i6,1x,i6,1x,2f15.3,9f11.5)')ILOC,JLOC,END_WB,BEG_WB,PRCP*DT,ECAN*DT,&
                EDIR*DT,ETRAN*DT,RUNSRF*DT,RUNSUB*DT,ZWT
           print*,(trim(message))
           print*,("Water budget problem in NOAHMP LSM")
        END IF
#endif
   ELSE                 !KWM
      ERRWAT = 0.0      !KWM
   ENDIF

 END SUBROUTINE ERROR

!== begin energy ===================================================================================

  SUBROUTINE ENERGY (parameters,ICE    ,VEGTYP ,IST    ,NSNOW  ,NSOIL  , & !in
                     ISNOW  ,DT     ,RHOAIR ,SFCPRS ,QAIR   , & !in
                     SFCTMP ,THAIR  ,LWDN   ,UU     ,VV     ,ZREF   , & !in
                     CO2AIR ,O2AIR  ,SOLAD  ,SOLAI  ,COSZ   ,IGS    , & !in
                     EAIR   ,TBOT   ,ZSNSO  ,ZSOIL  , & !in
                     ELAI   ,ESAI   ,FWET   ,FOLN   ,         & !in
                     FVEG   ,PAHV   ,PAHG   ,PAHB   ,                 & !in
                     QSNOW  ,DZSNSO ,LAT    ,CANLIQ ,CANICE ,ILOC   , JLOC, & !in
		     Z0WRF  ,                                         &
                     IMELT  ,SNICEV ,SNLIQV ,EPORE  ,T2M    ,FSNO   , & !out
                     SAV    ,SAG    ,QMELT  ,FSA    ,FSR    ,TAUX   , & !out
                     TAUY   ,FIRA   ,FSH    ,FCEV   ,FGEV   ,FCTR   , & !out
                     TRAD   ,PSN    ,APAR   ,SSOIL  ,BTRANI ,BTRAN  , & !out
                     PONDING,TS     ,LATHEAV , LATHEAG , frozen_canopy,frozen_ground,                       & !out
                     TV     ,TG     ,STC    ,SNOWH  ,EAH    ,TAH    , & !inout
                     SNEQVO ,SNEQV  ,SH2O   ,SMC    ,SNICE  ,SNLIQ  , & !inout
                     ALBOLD ,CM     ,CH     ,DX     ,DZ8W   ,Q2     , &   !inout
                     TAUSS  ,                                         & !inout
!jref:start
                     QC     ,QSFC   ,PSFC   , & !in 
                     T2MV   ,T2MB   ,FSRV   , &
                     FSRG   ,RSSUN  ,RSSHA  ,BGAP   ,WGAP,TGV,TGB,&
                     Q1     ,Q2V    ,Q2B    ,Q2E    ,CHV  ,CHB, EMISSI,PAH  ,&
		     SHG,SHC,SHB,EVG,EVB,GHV,GHB,IRG,IRC,IRB,TR,EVC,CHLEAF,CHUC,CHV2,CHB2, &
!-> xinyf and cheyz
             FIRE, ALBD_OUT, ALBI_OUT)   !out 
!<-
!jref:end                            

! --------------------------------------------------------------------------------------------------
! we use different approaches to deal with subgrid features of radiation transfer and turbulent
! transfer. We use 'tile' approach to compute turbulent fluxes, while we use modified two-
! stream to compute radiation transfer. Tile approach, assemblying vegetation canopies together,
! may expose too much ground surfaces (either covered by snow or grass) to solar radiation. The
! modified two-stream assumes vegetation covers fully the gridcell but with gaps between tree
! crowns.
! --------------------------------------------------------------------------------------------------
! turbulence transfer : 'tile' approach to compute energy fluxes in vegetated fraction and
!                         bare fraction separately and then sum them up weighted by fraction
!                     --------------------------------------
!                    / O  O  O  O  O  O  O  O  /          / 
!                   /  |  |  |  |  |  |  |  | /          /
!                  / O  O  O  O  O  O  O  O  /          /
!                 /  |  |  |tile1|  |  |  | /  tile2   /
!                / O  O  O  O  O  O  O  O  /  bare    /
!               /  |  |  | vegetated |  | /          /
!              / O  O  O  O  O  O  O  O  /          /
!             /  |  |  |  |  |  |  |  | /          /
!            --------------------------------------
! --------------------------------------------------------------------------------------------------
! radiation transfer : modified two-stream (Yang and Friedl, 2003, JGR; Niu ang Yang, 2004, JGR)
!                     --------------------------------------  two-stream treats leaves as
!                    /   O   O   O   O   O   O   O   O    /  cloud over the entire grid-cell,
!                   /    |   |   |   |   |   |   |   |   / while the modified two-stream 
!                  /   O   O   O   O   O   O   O   O    / aggregates cloudy leaves into  
!                 /    |   |   |   |   |   |   |   |   / tree crowns with gaps (as shown in
!                /   O   O   O   O   O   O   O   O    / the left figure). We assume these
!               /    |   |   |   |   |   |   |   |   / tree crowns are evenly distributed
!              /   O   O   O   O   O   O   O   O    / within the gridcell with 100% veg
!             /    |   |   |   |   |   |   |   |   / fraction, but with gaps. The 'tile'
!            -------------------------------------- approach overlaps too much shadows.
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  integer                           , INTENT(IN)    :: ILOC
  integer                           , INTENT(IN)    :: JLOC
  INTEGER                           , INTENT(IN)    :: ICE    !ice (ice = 1)
  INTEGER                           , INTENT(IN)    :: VEGTYP !vegetation physiology type
  INTEGER                           , INTENT(IN)    :: IST    !surface type: 1->soil; 2->lake
  INTEGER                           , INTENT(IN)    :: NSNOW  !maximum no. of snow layers        
  INTEGER                           , INTENT(IN)    :: NSOIL  !number of soil layers
  INTEGER                           , INTENT(IN)    :: ISNOW  !actual no. of snow layers
  REAL(r8)                              , INTENT(IN)    :: DT     !time step [sec]
  REAL(r8)                              , INTENT(IN)    :: QSNOW  !snowfall on the ground (mm/s)
  REAL(r8)                              , INTENT(IN)    :: RHOAIR !density air (kg/m3)
  REAL(r8)                              , INTENT(IN)    :: EAIR   !vapor pressure air (pa)
  REAL(r8)                              , INTENT(IN)    :: SFCPRS !pressure (pa)
  REAL(r8)                              , INTENT(IN)    :: QAIR   !specific humidity (kg/kg)
  REAL(r8)                              , INTENT(IN)    :: SFCTMP !air temperature (k)
  REAL(r8)                              , INTENT(IN)    :: THAIR  !potential temperature (k)
  REAL(r8)                              , INTENT(IN)    :: LWDN   !downward longwave radiation (w/m2)
  REAL(r8)                              , INTENT(IN)    :: UU     !wind speed in e-w dir (m/s)
  REAL(r8)                              , INTENT(IN)    :: VV     !wind speed in n-s dir (m/s)
  REAL(r8)   , DIMENSION(       1:    2), INTENT(IN)    :: SOLAD  !incoming direct solar rad. (w/m2)
  REAL(r8)   , DIMENSION(       1:    2), INTENT(IN)    :: SOLAI  !incoming diffuse solar rad. (w/m2)
  REAL(r8)   , DIMENSION(       1:    2), INTENT(out)    :: ALBD_OUT  !cheyz for coupling                        
  REAL(r8)   , DIMENSION(       1:    2), INTENT(out)    :: ALBI_OUT  !
  REAL(r8)                              , INTENT(IN)    :: COSZ   !cosine solar zenith angle (0-1)
  REAL(r8)                              , INTENT(IN)    :: ELAI   !LAI adjusted for burying by snow
  REAL(r8)                              , INTENT(IN)    :: ESAI   !LAI adjusted for burying by snow
  REAL(r8)                              , INTENT(IN)    :: FWET   !fraction of canopy that is wet [-]
  REAL(r8)                              , INTENT(IN)    :: FVEG   !greeness vegetation fraction (-)
  REAL(r8)                              , INTENT(IN)    :: LAT    !latitude (radians)
  REAL(r8)                              , INTENT(IN)    :: CANLIQ !canopy-intercepted liquid water (mm)
  REAL(r8)                              , INTENT(IN)    :: CANICE !canopy-intercepted ice mass (mm)
  REAL(r8)                              , INTENT(IN)    :: FOLN   !foliage nitrogen (%)
  REAL(r8)                              , INTENT(IN)    :: CO2AIR !atmospheric co2 concentration (pa)
  REAL(r8)                              , INTENT(IN)    :: O2AIR  !atmospheric o2 concentration (pa)
  REAL(r8)                              , INTENT(IN)    :: IGS    !growing season index (0=off, 1=on)

  REAL(r8)                              , INTENT(IN)    :: ZREF   !reference height (m)
  REAL(r8)                              , INTENT(IN)    :: TBOT   !bottom condition for soil temp. (k) 
  REAL(r8)   , DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)    :: ZSNSO  !layer-bottom depth from snow surf [m]
  REAL(r8)   , DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL  !layer-bottom depth from soil surf [m]
  REAL(r8)   , DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)    :: DZSNSO !depth of snow & soil layer-bottom [m]
  REAL(r8), INTENT(IN)   :: PAHV    !precipitation advected heat - vegetation net (W/m2)
  REAL(r8), INTENT(IN)   :: PAHG    !precipitation advected heat - under canopy net (W/m2)
  REAL(r8), INTENT(IN)   :: PAHB    !precipitation advected heat - bare ground net (W/m2)

!jref:start; in 
  REAL(r8)                              , INTENT(IN)    :: QC     !cloud water mixing ratio
  REAL(r8)                              , INTENT(INOUT) :: QSFC   !mixing ratio at lowest model layer
  REAL(r8)                              , INTENT(IN)    :: PSFC   !pressure at lowest model layer
  REAL(r8)                              , INTENT(IN)    :: DX     !horisontal resolution
  REAL(r8)                              , INTENT(IN)    :: DZ8W   !thickness of lowest layer
  REAL(r8)                              , INTENT(IN)    :: Q2     !mixing ratio (kg/kg)
!jref:end

! outputs
  REAL(r8)                              , INTENT(OUT)   :: Z0WRF  !combined z0 sent to coupled model
  INTEGER, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT)   :: IMELT  !phase change index [1-melt; 2-freeze]
  REAL(r8)   , DIMENSION(-NSNOW+1:    0), INTENT(OUT)   :: SNICEV !partial volume ice [m3/m3]
  REAL(r8)   , DIMENSION(-NSNOW+1:    0), INTENT(OUT)   :: SNLIQV !partial volume liq. water [m3/m3]
  REAL(r8)   , DIMENSION(-NSNOW+1:    0), INTENT(OUT)   :: EPORE  !effective porosity [m3/m3]
  REAL(r8)                              , INTENT(OUT)   :: FSNO   !snow cover fraction (-)
  REAL(r8)                              , INTENT(OUT)   :: QMELT  !snowmelt [mm/s]
  REAL(r8)                              , INTENT(OUT)   :: PONDING!pounding at ground [mm]
  REAL(r8)                              , INTENT(OUT)   :: SAV    !solar rad. absorbed by veg. (w/m2)
  REAL(r8)                              , INTENT(OUT)   :: SAG    !solar rad. absorbed by ground (w/m2)
  REAL(r8)                              , INTENT(OUT)   :: FSA    !tot. absorbed solar radiation (w/m2)
  REAL(r8)                              , INTENT(OUT)   :: FSR    !tot. reflected solar radiation (w/m2)
  REAL(r8)                              , INTENT(OUT)   :: TAUX   !wind stress: e-w (n/m2)
  REAL(r8)                              , INTENT(OUT)   :: TAUY   !wind stress: n-s (n/m2)
  REAL(r8)                              , INTENT(OUT)   :: FIRA   !total net LW. rad (w/m2)   [+ to atm]
  REAL(r8)                              , INTENT(OUT)   :: FSH    !total sensible heat (w/m2) [+ to atm]
  REAL(r8)                              , INTENT(OUT)   :: FCEV   !canopy evaporation (w/m2)  [+ to atm]
  REAL(r8)                              , INTENT(OUT)   :: FGEV   !ground evaporation (w/m2)  [+ to atm]
  REAL(r8)                              , INTENT(OUT)   :: FCTR   !transpiration (w/m2)       [+ to atm]
  REAL(r8)                              , INTENT(OUT)   :: TRAD   !radiative temperature (k)
  REAL(r8)                              , INTENT(OUT)   :: T2M    !2 m height air temperature (k)
  REAL(r8)                              , INTENT(OUT)   :: PSN    !total photosyn. (umolco2/m2/s) [+]
  REAL(r8)                              , INTENT(OUT)   :: APAR   !total photosyn. active energy (w/m2)
  REAL(r8)                              , INTENT(OUT)   :: SSOIL  !ground heat flux (w/m2)   [+ to soil]
  REAL(r8)   , DIMENSION(       1:NSOIL), INTENT(OUT)   :: BTRANI !soil water transpiration factor (0-1)
  REAL(r8)                              , INTENT(OUT)   :: BTRAN  !soil water transpiration factor (0-1)
!  REAL(r8)                              , INTENT(OUT)   :: LATHEA !latent heat vap./sublimation (j/kg)
  REAL(r8)                              , INTENT(OUT)   :: LATHEAV !latent heat vap./sublimation (j/kg)
  REAL(r8)                              , INTENT(OUT)   :: LATHEAG !latent heat vap./sublimation (j/kg)
  LOGICAL                           , INTENT(OUT)   :: FROZEN_GROUND ! used to define latent heat pathway
  LOGICAL                           , INTENT(OUT)   :: FROZEN_CANOPY ! used to define latent heat pathway

!jref:start  
  REAL(r8)                              , INTENT(OUT)   :: FSRV    !veg. reflected solar radiation (w/m2)
  REAL(r8)                              , INTENT(OUT)   :: FSRG    !ground reflected solar radiation (w/m2)
  REAL(r8), INTENT(OUT) :: RSSUN        !sunlit leaf stomatal resistance (s/m)
  REAL(r8), INTENT(OUT) :: RSSHA        !shaded leaf stomatal resistance (s/m)
!jref:end - out for debug  

!jref:start; output
  REAL(r8)                              , INTENT(OUT)   :: T2MV   !2-m air temperature over vegetated part [k]
  REAL(r8)                              , INTENT(OUT)   :: T2MB   !2-m air temperature over bare ground part [k]
  REAL(r8)                              , INTENT(OUT)   :: BGAP
  REAL(r8)                              , INTENT(OUT)   :: WGAP
!jref:end

! input & output
  REAL(r8)                              , INTENT(INOUT) :: TS     !surface temperature (k)
  REAL(r8)                              , INTENT(INOUT) :: TV     !vegetation temperature (k)
  REAL(r8)                              , INTENT(INOUT) :: TG     !ground temperature (k)
  REAL(r8)   , DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    !snow/soil temperature [k]
  REAL(r8)                              , INTENT(INOUT) :: SNOWH  !snow height [m]
  REAL(r8)                              , INTENT(INOUT) :: SNEQV  !snow mass (mm)
  REAL(r8)                              , INTENT(INOUT) :: SNEQVO !snow mass at last time step (mm)
  REAL(r8)   , DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   !liquid soil moisture [m3/m3]
  REAL(r8)   , DIMENSION(       1:NSOIL), INTENT(INOUT) :: SMC    !soil moisture (ice + liq.) [m3/m3]
  REAL(r8)   , DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  !snow ice mass (kg/m2)
  REAL(r8)   , DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  !snow liq mass (kg/m2)
  REAL(r8)                              , INTENT(INOUT) :: EAH    !canopy air vapor pressure (pa)
  REAL(r8)                              , INTENT(INOUT) :: TAH    !canopy air temperature (k)
  REAL(r8)                              , INTENT(INOUT) :: ALBOLD !snow albedo at last time step(CLASS type)
  REAL(r8)                              , INTENT(INOUT) :: TAUSS  !non-dimensional snow age
  REAL(r8)                              , INTENT(INOUT) :: CM     !momentum drag coefficient
  REAL(r8)                              , INTENT(INOUT) :: CH     !sensible heat exchange coefficient
  REAL(r8)                              , INTENT(INOUT) :: Q1
!  REAL(r8)                                              :: Q2E
  REAL(r8),                               INTENT(OUT)   :: EMISSI
  REAL(r8),                               INTENT(OUT)   :: PAH    !precipitation advected heat - total (W/m2)

! local
  INTEGER                                           :: IZ     !do-loop index
  LOGICAL                                           :: VEG    !true if vegetated surface
  REAL(r8)                                              :: UR     !wind speed at height ZLVL (m/s)
  REAL(r8)                                              :: ZLVL   !reference height (m)
  REAL(r8)                                              :: FSUN   !sunlit fraction of canopy [-]
  REAL(r8)                                              :: RB     !leaf boundary layer resistance (s/m)
  REAL(r8)                                              :: RSURF  !ground surface resistance (s/m)
  REAL(r8)                                              :: L_RSURF!Dry-layer thickness for computing RSURF (Sakaguchi and Zeng, 2009)
  REAL(r8)                                              :: D_RSURF!Reduced vapor diffusivity in soil for computing RSURF (SZ09)
  REAL(r8)                                              :: BEVAP  !soil water evaporation factor (0- 1)
  REAL(r8)                                              :: MOL    !Monin-Obukhov length (m)
  REAL(r8)                                              :: VAI    !sum of LAI  + stem area index [m2/m2]
  REAL(r8)                                              :: CWP    !canopy wind extinction parameter
  REAL(r8)                                              :: ZPD    !zero plane displacement (m)
  REAL(r8)                                              :: Z0M    !z0 momentum (m)
  REAL(r8)                                              :: ZPDG   !zero plane displacement (m)
  REAL(r8)                                              :: Z0MG   !z0 momentum, ground (m)
  REAL(r8)                                              :: EMV    !vegetation emissivity
  REAL(r8)                                              :: EMG    !ground emissivity
  REAL(r8), intent(out)                                 :: FIRE   !emitted IR (w/m2)
  REAL(r8)                                              :: LAISUN !sunlit leaf area index (m2/m2)
  REAL(r8)                                              :: LAISHA !shaded leaf area index (m2/m2)
  REAL(r8)                                              :: PSNSUN !sunlit photosynthesis (umolco2/m2/s)
  REAL(r8)                                              :: PSNSHA !shaded photosynthesis (umolco2/m2/s)
!jref:start - for debug  
!  REAL(r8)                                              :: RSSUN  !sunlit stomatal resistance (s/m)
!  REAL(r8)                                              :: RSSHA  !shaded stomatal resistance (s/m)
!jref:end - for debug
  REAL(r8)                                              :: PARSUN !par absorbed per sunlit LAI (w/m2)
  REAL(r8)                                              :: PARSHA !par absorbed per shaded LAI (w/m2)

  REAL(r8), DIMENSION(-NSNOW+1:NSOIL)                   :: FACT   !temporary used in phase change
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL)                   :: DF     !thermal conductivity [w/m/k]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL)                   :: HCPCT  !heat capacity [j/m3/k]
  REAL(r8)                                              :: BDSNO  !bulk density of snow (kg/m3)
  REAL(r8)                                              :: FMELT  !melting factor for snow cover frac
  REAL(r8)                                              :: GX     !temporary variable
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL)                   :: PHI    !light through water (w/m2)
!  REAL(r8)                                              :: GAMMA  !psychrometric constant (pa/k)
  REAL(r8)                                              :: GAMMAV  !psychrometric constant (pa/k)
  REAL(r8)                                              :: GAMMAG  !psychrometric constant (pa/k)
  REAL(r8)                                              :: PSI    !surface layer soil matrix potential (m)
  REAL(r8)                                              :: RHSUR  !raltive humidity in surface soil/snow air space (-)

! temperature and fluxes over vegetated fraction

  REAL(r8)                                              :: TAUXV  !wind stress: e-w dir [n/m2]
  REAL(r8)                                              :: TAUYV  !wind stress: n-s dir [n/m2]
  REAL(r8),INTENT(OUT)                                              :: IRC    !canopy net LW rad. [w/m2] [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: IRG    !ground net LW rad. [w/m2] [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: SHC    !canopy sen. heat [w/m2]   [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: SHG    !ground sen. heat [w/m2]   [+ to atm]
!jref:start  
  REAL(r8),INTENT(OUT)                                  :: Q2V
  REAL(r8),INTENT(OUT)                                  :: Q2B
  REAL(r8),INTENT(OUT)                                  :: Q2E
!jref:end  
  REAL(r8),INTENT(OUT)                                              :: EVC    !canopy evap. heat [w/m2]  [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: EVG    !ground evap. heat [w/m2]  [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: TR     !transpiration heat [w/m2] [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: GHV    !ground heat flux [w/m2]  [+ to soil]
  REAL(r8),INTENT(OUT)                                  :: TGV    !ground surface temp. [k]
  REAL(r8)                                              :: CMV    !momentum drag coefficient
  REAL(r8),INTENT(OUT)                                  :: CHV    !sensible heat exchange coefficient

! temperature and fluxes over bare soil fraction

  REAL(r8)                                              :: TAUXB  !wind stress: e-w dir [n/m2]
  REAL(r8)                                              :: TAUYB  !wind stress: n-s dir [n/m2]
  REAL(r8),INTENT(OUT)                                              :: IRB    !net longwave rad. [w/m2] [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: SHB    !sensible heat [w/m2]     [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: EVB    !evaporation heat [w/m2]  [+ to atm]
  REAL(r8),INTENT(OUT)                                              :: GHB    !ground heat flux [w/m2] [+ to soil]
  REAL(r8),INTENT(OUT)                                  :: TGB    !ground surface temp. [k]
  REAL(r8)                                              :: CMB    !momentum drag coefficient
  REAL(r8),INTENT(OUT)                                  :: CHB    !sensible heat exchange coefficient
  REAL(r8),INTENT(OUT)                                  :: CHLEAF !leaf exchange coefficient
  REAL(r8),INTENT(OUT)                                  :: CHUC   !under canopy exchange coefficient
!jref:start  
  REAL(r8),INTENT(OUT)                                  :: CHV2    !sensible heat conductance, canopy air to ZLVL air (m/s)
  REAL(r8),INTENT(OUT)                                  :: CHB2    !sensible heat conductance, canopy air to ZLVL air (m/s)
  REAL(r8)                                  :: noahmpres

!jref:end  

  REAL(r8), PARAMETER                   :: MPE    = 1.E-6
  REAL(r8), PARAMETER                   :: PSIWLT = -150.  !metric potential for wilting point (m)
  REAL(r8), PARAMETER                   :: Z0     = 0.01   ! Bare-soil roughness length (m) (i.e., under the canopy)

! ---------------------------------------------------------------------------------------------------
! initialize fluxes from veg. fraction

    TAUXV     = 0.    
    TAUYV     = 0.
    IRC       = 0.
    SHC       = 0.
    IRG       = 0.
    SHG       = 0.
    EVG       = 0.       
    EVC       = 0.
    TR        = 0.
    GHV       = 0.       
    PSNSUN    = 0.
    PSNSHA    = 0.
    T2MV      = 0.
    Q2V       = 0.
    CHV       = 0.
    CHLEAF    = 0.
    CHUC      = 0.
    CHV2      = 0.

! wind speed at reference height: ur >= 1

    UR = MAX( SQRT(UU**2.+VV**2.), 1. )

! vegetated or non-vegetated

    VAI = ELAI + ESAI
    VEG = .FALSE.
    IF(VAI > 0.) VEG = .TRUE.

! ground snow cover fraction [Niu and Yang, 2007, JGR]

     FSNO = 0.
! zhangyi xinyf    IF(SNOWH.GT.0.)  THEN
     IF(SNOWH.GT. 0.)  THEN
         BDSNO    = SNEQV / SNOWH
         FMELT    = (BDSNO/100.)**parameters%MFSNO
         FSNO     = TANH( SNOWH /(2.5* Z0 * FMELT))
#ifndef REGLSM 
         FSNO    = SNOWH/(0.1+SNOWH)  ! from CLM3 
#endif
     ENDIF

! ground roughness length

     IF(IST == 2) THEN
       IF(TG .LE. TFRZ) THEN
         Z0MG = 0.01 * (1.0-FSNO) + FSNO * parameters%Z0SNO
       ELSE
         Z0MG = 0.01  
       END IF
     ELSE
       Z0MG = Z0 * (1.0-FSNO) + FSNO * parameters%Z0SNO
     END IF

! roughness length and displacement height

     ZPDG  = SNOWH
     IF(VEG) THEN
        Z0M  = parameters%Z0MVT
        ZPD  = 0.65 * parameters%HVT
        IF(SNOWH.GT.ZPD) ZPD  = SNOWH
     ELSE
        Z0M  = Z0MG
        ZPD  = ZPDG
     END IF

     ZLVL = MAX(ZPD,parameters%HVT) + ZREF
     IF(ZPDG >= ZLVL) ZLVL = ZPDG + ZREF
!     UR   = UR*LOG(ZLVL/Z0M)/LOG(10./Z0M)       !input UR is at 10m

! canopy wind absorption coeffcient

     CWP = parameters%CWPVT

! Thermal properties of soil, snow, lake, and frozen soil

  CALL THERMOPROP (parameters,NSOIL   ,NSNOW   ,ISNOW   ,IST     ,DZSNSO  , & !in
                   DT      ,SNOWH   ,SNICE   ,SNLIQ   , & !in
                   SMC     ,SH2O    ,TG      ,STC     ,UR      , & !in
                   LAT     ,Z0M     ,ZLVL    ,VEGTYP  , & !in
                   DF      ,HCPCT   ,SNICEV  ,SNLIQV  ,EPORE   , & !out
                   FACT    )                              !out

! Solar radiation: absorbed & reflected by the ground and canopy

  CALL  RADIATION (parameters,VEGTYP  ,IST     ,ICE     ,NSOIL   , & !in 
                   SNEQVO  ,SNEQV   ,DT      ,COSZ    ,SNOWH   , & !in
                   TG      ,TV      ,FSNO    ,QSNOW   ,FWET    , & !in
                   ELAI    ,ESAI    ,SMC     ,SOLAD   ,SOLAI   , & !in
                   FVEG    ,ILOC    ,JLOC    ,                   & !in
                   ALBOLD  ,TAUSS   ,                            & !inout
                   FSUN    ,LAISUN  ,LAISHA  ,PARSUN  ,PARSHA  , & !out
                   SAV     ,SAG     ,FSR     ,FSA     ,FSRV    , & 
                   FSRG    ,BGAP    ,WGAP,                       &
                   ALBD_OUT, ALBI_OUT)            !out cheyz

! vegetation and ground emissivity

     EMV = 1. - EXP(-(ELAI+ESAI)/1.0)
     IF (ICE == 1) THEN
       EMG = 0.98*(1.-FSNO) + 1.0*FSNO
     ELSE
       EMG = parameters%EG(IST)*(1.-FSNO) + 1.0*FSNO
     END IF

! soil moisture factor controlling stomatal resistance
   
     BTRAN = 0.

     IF(IST ==1 ) THEN
       DO IZ = 1, parameters%NROOT
          IF(OPT_BTR == 1) then                  ! Noah
            GX    = (SH2O(IZ)-parameters%SMCWLT(IZ)) / (parameters%SMCREF(IZ)-parameters%SMCWLT(IZ))
          END IF
          IF(OPT_BTR == 2) then                  ! CLM
            PSI   = MAX(PSIWLT,-parameters%PSISAT(IZ)*(MAX(0.01,SH2O(IZ))/parameters%SMCMAX(IZ))**(-parameters%BEXP(IZ)) )
            GX    = (1.-PSI/PSIWLT)/(1.+parameters%PSISAT(IZ)/PSIWLT)
          END IF
          IF(OPT_BTR == 3) then                  ! SSiB
            PSI   = MAX(PSIWLT,-parameters%PSISAT(IZ)*(MAX(0.01,SH2O(IZ))/parameters%SMCMAX(IZ))**(-parameters%BEXP(IZ)) )
            GX    = 1.-EXP(-5.8*(LOG(PSIWLT/PSI))) 
          END IF
       
          GX = MIN(1.,MAX(0.,GX))
          BTRANI(IZ) = MAX(MPE,DZSNSO(IZ) / (-ZSOIL(parameters%NROOT)) * GX)
          BTRAN      = BTRAN + BTRANI(IZ)
       END DO
       BTRAN = MAX(MPE,BTRAN)

       BTRANI(1:parameters%NROOT) = BTRANI(1:parameters%NROOT)/BTRAN
     END IF

! soil surface resistance for ground evap.

     BEVAP = MAX(0.0,SH2O(1)/parameters%SMCMAX(1))
     IF(IST == 2) THEN
       RSURF = 1.          ! avoid being divided by 0
       RHSUR = 1.0
     ELSE

       IF(OPT_RSF == 1 .OR. OPT_RSF == 4) THEN
         ! RSURF based on Sakaguchi and Zeng, 2009
         ! taking the "residual water content" to be the wilting point, 
         ! and correcting the exponent on the D term (typo in SZ09 ?)
         L_RSURF = (-ZSOIL(1)) * ( exp ( (1.0 - MIN(1.0,SH2O(1)/parameters%SMCMAX(1))) ** 5 ) - 1.0 ) / ( 2.71828 - 1.0 ) 
         D_RSURF = 2.2E-5 * parameters%SMCMAX(1) * parameters%SMCMAX(1) * ( 1.0 - parameters%SMCWLT(1) / parameters%SMCMAX(1) ) ** (2.0+3.0/parameters%BEXP(1))
         RSURF = L_RSURF / D_RSURF
       ELSEIF(OPT_RSF == 2) THEN
         RSURF = FSNO * 1. + (1.-FSNO)* EXP(8.25-4.225*BEVAP) !Sellers (1992) ! Older RSURF computations
       ELSEIF(OPT_RSF == 3) THEN
         RSURF = FSNO * 1. + (1.-FSNO)* EXP(8.25-6.0  *BEVAP) !adjusted to decrease RSURF for wet soil
       ENDIF

       IF(OPT_RSF == 4) THEN  ! AD: FSNO weighted; snow RSURF set in MPTABLE v3.8
         RSURF = 1. / (FSNO * (1./parameters%RSURF_SNOW) + (1.-FSNO) * (1./max(RSURF, 0.001)))
       ENDIF

       IF(SH2O(1) < 0.01 .and. SNOWH == 0.) RSURF = 1.E6
       PSI   = -parameters%PSISAT(1)*(MAX(0.01,SH2O(1))/parameters%SMCMAX(1))**(-parameters%BEXP(1))   
       RHSUR = FSNO + (1.-FSNO) * EXP(PSI*GRAV/(RW*TG)) 
     END IF

! urban - jref 
     IF (parameters%urban_flag .and. SNOWH == 0. ) THEN
        RSURF = 1.E6
     ENDIF

! set psychrometric constant

     IF (TV .GT. TFRZ) THEN           ! Barlage: add distinction between ground and 
        LATHEAV = HVAP                ! vegetation in v3.6
	frozen_canopy = .false.
     ELSE
        LATHEAV = HSUB
	frozen_canopy = .true.
     END IF
     GAMMAV = CPAIR*SFCPRS/(0.622*LATHEAV)

     IF (TG .GT. TFRZ) THEN
        LATHEAG = HVAP
	frozen_ground = .false.
     ELSE
        LATHEAG = HSUB
	frozen_ground = .true.
     END IF
     GAMMAG = CPAIR*SFCPRS/(0.622*LATHEAG)

!     IF (SFCTMP .GT. TFRZ) THEN
!        LATHEA = HVAP
!     ELSE
!        LATHEA = HSUB
!     END IF
!     GAMMA = CPAIR*SFCPRS/(0.622*LATHEA)

! Surface temperatures of the ground and canopy and energy fluxes

    IF (VEG .AND. FVEG > 0) THEN 
    TGV = TG
    CMV = CM
    CHV = CH
    CALL VEGE_FLUX (parameters,NSNOW   ,NSOIL   ,ISNOW   ,VEGTYP  ,VEG     , & !in
                    DT      ,SAV     ,SAG     ,LWDN    ,UR      , & !in
                    UU      ,VV      ,SFCTMP  ,THAIR   ,QAIR    , & !in
                    EAIR    ,RHOAIR  ,SNOWH   ,VAI     ,GAMMAV   ,GAMMAG   , & !in
                    FWET    ,LAISUN  ,LAISHA  ,CWP     ,DZSNSO  , & !in
                    ZLVL    ,ZPD     ,Z0M     ,FVEG    , & !in
                    Z0MG    ,EMV     ,EMG     ,CANLIQ  ,FSNO, & !in
                    CANICE  ,STC     ,DF      ,RSSUN   ,RSSHA   , & !in
                    RSURF   ,LATHEAV ,LATHEAG ,PARSUN  ,PARSHA  ,IGS     , & !in
                    FOLN    ,CO2AIR  ,O2AIR   ,BTRAN   ,SFCPRS  , & !in
                    RHSUR   ,ILOC    ,JLOC    ,Q2      ,PAHV  ,PAHG  , & !in
                    EAH     ,TAH     ,TV      ,TGV     ,CMV     , & !inout
                    CHV     ,DX      ,DZ8W    ,                   & !inout
                    TAUXV   ,TAUYV   ,IRG     ,IRC     ,SHG     , & !out
                    SHC     ,EVG     ,EVC     ,TR      ,GHV     , & !out
                    T2MV    ,PSNSUN  ,PSNSHA  ,                   & !out
!jref:start
                    QC      ,QSFC    ,PSFC    , & !in
                    Q2V     ,CHV2, CHLEAF, CHUC)               !inout 
!jref:end                            
    END IF

    TGB = TG
    CMB = CM
    CHB = CH
    CALL BARE_FLUX (parameters,NSNOW   ,NSOIL   ,ISNOW   ,DT      ,SAG     , & !in
                    LWDN    ,UR      ,UU      ,VV      ,SFCTMP  , & !in
                    THAIR   ,QAIR    ,EAIR    ,RHOAIR  ,SNOWH   , & !in
                    DZSNSO  ,ZLVL    ,ZPDG    ,Z0MG    ,FSNO,          & !in
                    EMG     ,STC     ,DF      ,RSURF   ,LATHEAG  , & !in
                    GAMMAG   ,RHSUR   ,ILOC    ,JLOC    ,Q2      ,PAHB  , & !in
                    TGB     ,CMB     ,CHB     ,                   & !inout
                    TAUXB   ,TAUYB   ,IRB     ,SHB     ,EVB     , & !out
                    GHB     ,T2MB    ,DX      ,DZ8W    ,VEGTYP  , & !out
!jref:start
                    QC      ,QSFC    ,PSFC    , & !in
                    SFCPRS  ,Q2B,   CHB2)                          !in 
!jref:end                            

!energy balance at vege canopy: SAV          =(IRC+SHC+EVC+TR)     *FVEG  at   FVEG 
!energy balance at vege ground: SAG*    FVEG =(IRG+SHG+EVG+GHV)    *FVEG  at   FVEG
!energy balance at bare ground: SAG*(1.-FVEG)=(IRB+SHB+EVB+GHB)*(1.-FVEG) at 1-FVEG

    IF (VEG .AND. FVEG > 0) THEN 
        TAUX  = FVEG * TAUXV     + (1.0 - FVEG) * TAUXB
        TAUY  = FVEG * TAUYV     + (1.0 - FVEG) * TAUYB
        FIRA  = FVEG * IRG       + (1.0 - FVEG) * IRB       + IRC
        FSH   = FVEG * SHG       + (1.0 - FVEG) * SHB       + SHC
        FGEV  = FVEG * EVG       + (1.0 - FVEG) * EVB
        SSOIL = FVEG * GHV       + (1.0 - FVEG) * GHB
        FCEV  = EVC
        FCTR  = TR
	PAH   = FVEG * PAHG      + (1.0 - FVEG) * PAHB   + PAHV
        TG    = FVEG * TGV       + (1.0 - FVEG) * TGB
        T2M   = FVEG * T2MV      + (1.0 - FVEG) * T2MB
        TS    = FVEG * TV        + (1.0 - FVEG) * TGB
        CM    = FVEG * CMV       + (1.0 - FVEG) * CMB      ! better way to average?
        CH    = FVEG * CHV       + (1.0 - FVEG) * CHB
        Q1    = FVEG * (EAH*0.622/(SFCPRS - 0.378*EAH)) + (1.0 - FVEG)*QSFC
        Q2E   = FVEG * Q2V       + (1.0 - FVEG) * Q2B
	Z0WRF = Z0M
    ELSE
        TAUX  = TAUXB
        TAUY  = TAUYB
        FIRA  = IRB
        FSH   = SHB
        FGEV  = EVB
        SSOIL = GHB
        TG    = TGB
        T2M   = T2MB
        FCEV  = 0.
        FCTR  = 0.
	PAH   = PAHB
        TS    = TG
        CM    = CMB
        CH    = CHB
        Q1    = QSFC
        Q2E   = Q2B
        RSSUN = 0.0
        RSSHA = 0.0
        TGV   = TGB
        CHV   = CHB
	Z0WRF = Z0MG
    END IF

    FIRE = LWDN + FIRA

    IF(FIRE <=0.) THEN
       WRITE(6,*) 'emitted longwave <0; skin T may be wrong due to inconsistent'
       WRITE(6,*) 'input of SHDFAC with LAI'
       WRITE(6,*) ILOC, JLOC, 'SHDFAC=',FVEG,'VAI=',VAI,'TV=',TV,'TG=',TG
       WRITE(6,*) 'LWDN=',LWDN,'FIRA=',FIRA,'SNOWH=',SNOWH
       print*, ("STOP in Noah-MP")
    END IF

    ! Compute a net emissivity
    EMISSI = FVEG * ( EMG*(1-EMV) + EMV + EMV*(1-EMV)*(1-EMG) ) + &
         (1-FVEG) * EMG

    ! When we're computing a TRAD, subtract from the emitted IR the
    ! reflected portion of the incoming LWDN, so we're just
    ! considering the IR originating in the canopy/ground system.
    
    TRAD = ( ( FIRE - (1-EMISSI)*LWDN ) / (EMISSI*SB) ) ** 0.25

    ! Old TRAD calculation not taking into account Emissivity:
    ! TRAD = (FIRE/SB)**0.25

    APAR = PARSUN*LAISUN + PARSHA*LAISHA
    PSN  = PSNSUN*LAISUN + PSNSHA*LAISHA

! 3L snow & 4L soil temperatures

    CALL TSNOSOI (parameters,ICE     ,NSOIL   ,NSNOW   ,ISNOW   ,IST     , & !in
                  TBOT    ,ZSNSO   ,SSOIL   ,DF      ,HCPCT   , & !in
                  SAG     ,DT      ,SNOWH   ,DZSNSO  , & !in
                  TG      ,ILOC    ,JLOC    ,                   & !in
                  STC     )                                       !inout

! adjusting snow surface temperature
     IF(OPT_STC == 2) THEN
      IF (SNOWH > 0.05 .AND. TG > TFRZ) THEN
        TGV = TFRZ
        TGB = TFRZ
          IF (VEG .AND. FVEG > 0) THEN
             TG    = FVEG * TGV       + (1.0 - FVEG) * TGB
             TS    = FVEG * TV        + (1.0 - FVEG) * TGB
          ELSE
             TG    = TGB
             TS    = TGB
          END IF
      END IF
     END IF

! Energy released or consumed by snow & frozen soil

 CALL PHASECHANGE (parameters,NSNOW   ,NSOIL   ,ISNOW   ,DT      ,FACT    , & !in
                   DZSNSO  ,HCPCT   ,IST     ,ILOC    ,JLOC    , & !in
                   STC     ,SNICE   ,SNLIQ   ,SNEQV   ,SNOWH   , & !inout
                   SMC     ,SH2O    ,                            & !inout
                   QMELT   ,IMELT   ,PONDING )                     !out


  END SUBROUTINE ENERGY

!== begin thermoprop ===============================================================================

  SUBROUTINE THERMOPROP (parameters,NSOIL   ,NSNOW   ,ISNOW   ,IST     ,DZSNSO  , & !in
                         DT      ,SNOWH   ,SNICE   ,SNLIQ   , & !in
                         SMC     ,SH2O    ,TG      ,STC     ,UR      , & !in
                         LAT     ,Z0M     ,ZLVL    ,VEGTYP  , & !in
                         DF      ,HCPCT   ,SNICEV  ,SNLIQV  ,EPORE   , & !out
                         FACT    )                                       !out
! ------------------------------------------------------------------------------------------------- 
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER                        , INTENT(IN)  :: NSOIL   !number of soil layers
  INTEGER                        , INTENT(IN)  :: NSNOW   !maximum no. of snow layers        
  INTEGER                        , INTENT(IN)  :: ISNOW   !actual no. of snow layers
  INTEGER                        , INTENT(IN)  :: IST     !surface type
  REAL(r8)                           , INTENT(IN)  :: DT      !time step [s]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(IN)  :: SNICE   !snow ice mass (kg/m2)
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(IN)  :: SNLIQ   !snow liq mass (kg/m2)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DZSNSO  !thickness of snow/soil layers [m]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)  :: SMC     !soil moisture (ice + liq.) [m3/m3]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)  :: SH2O    !liquid soil moisture [m3/m3]
  REAL(r8)                           , INTENT(IN)  :: SNOWH   !snow height [m]
  REAL(r8),                            INTENT(IN)  :: TG      !surface temperature (k)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: STC     !snow/soil/lake temp. (k)
  REAL(r8),                            INTENT(IN)  :: UR      !wind speed at ZLVL (m/s)
  REAL(r8),                            INTENT(IN)  :: LAT     !latitude (radians)
  REAL(r8),                            INTENT(IN)  :: Z0M     !roughness length (m)
  REAL(r8),                            INTENT(IN)  :: ZLVL    !reference height (m)
  INTEGER                        , INTENT(IN)  :: VEGTYP  !vegtyp type

! outputs
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: DF      !thermal conductivity [w/m/k]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: HCPCT   !heat capacity [j/m3/k]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNICEV  !partial volume of ice [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNLIQV  !partial volume of liquid water [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: EPORE   !effective porosity [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: FACT    !computing energy for phase change
! --------------------------------------------------------------------------------------------------
! locals

  INTEGER :: IZ
  REAL(r8), DIMENSION(-NSNOW+1:    0)              :: CVSNO   !volumetric specific heat (j/m3/k)
  REAL(r8), DIMENSION(-NSNOW+1:    0)              :: TKSNO   !snow thermal conductivity (j/m3/k)
  REAL(r8), DIMENSION(       1:NSOIL)              :: SICE    !soil ice content
! --------------------------------------------------------------------------------------------------

! compute snow thermal conductivity and heat capacity

    CALL CSNOW (parameters,ISNOW   ,NSNOW   ,NSOIL   ,SNICE   ,SNLIQ   ,DZSNSO  , & !in
                TKSNO   ,CVSNO   ,SNICEV  ,SNLIQV  ,EPORE   )   !out

    DO IZ = ISNOW+1, 0
      DF   (IZ) = TKSNO(IZ)
      HCPCT(IZ) = CVSNO(IZ)
    END DO

! compute soil thermal properties

    DO  IZ = 1, NSOIL
       SICE(IZ)  = SMC(IZ) - SH2O(IZ)
       HCPCT(IZ) = SH2O(IZ)*CWAT + (1.0-parameters%SMCMAX(IZ))*parameters%CSOIL &
                + (parameters%SMCMAX(IZ)-SMC(IZ))*CPAIR + SICE(IZ)*CICE
       CALL TDFCND (parameters,IZ,DF(IZ), SMC(IZ), SH2O(IZ))
    END DO
       
    IF ( parameters%urban_flag ) THEN
       DO IZ = 1,NSOIL
         DF(IZ) = 3.24
       END DO
    ENDIF

! heat flux reduction effect from the overlying green canopy, adapted from 
! section 2.1.2 of Peters-Lidard et al. (1997, JGR, VOL 102(D4)).
! not in use because of the separation of the canopy layer from the ground.
! but this may represent the effects of leaf litter (Niu comments)
!       DF1 = DF1 * EXP (SBETA * SHDFAC)

! compute lake thermal properties 
! (no consideration of turbulent mixing for this version)

    IF(IST == 2) THEN
       DO IZ = 1, NSOIL 
         IF(STC(IZ) > TFRZ) THEN
            HCPCT(IZ) = CWAT
            DF(IZ)    = TKWAT  !+ KEDDY * CWAT 
         ELSE
            HCPCT(IZ) = CICE
            DF(IZ)    = TKICE 
         END IF
       END DO
    END IF

! combine a temporary variable used for melting/freezing of snow and frozen soil

    DO IZ = ISNOW+1,NSOIL
     FACT(IZ) = DT/(HCPCT(IZ)*DZSNSO(IZ))
    END DO

! snow/soil interface

    IF(ISNOW == 0) THEN
       DF(1) = (DF(1)*DZSNSO(1)+0.35*SNOWH)      / (SNOWH    +DZSNSO(1)) 
    ELSE
       DF(1) = (DF(1)*DZSNSO(1)+DF(0)*DZSNSO(0)) / (DZSNSO(0)+DZSNSO(1))
    END IF


  END SUBROUTINE THERMOPROP

!== begin csnow ====================================================================================

  SUBROUTINE CSNOW (parameters,ISNOW   ,NSNOW   ,NSOIL   ,SNICE   ,SNLIQ   ,DZSNSO  , & !in
                    TKSNO   ,CVSNO   ,SNICEV  ,SNLIQV  ,EPORE   )   !out
! --------------------------------------------------------------------------------------------------
! Snow bulk density,volumetric capacity, and thermal conductivity
!---------------------------------------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                          INTENT(IN) :: ISNOW  !number of snow layers (-)            
  INTEGER                        ,  INTENT(IN) :: NSNOW  !maximum no. of snow layers        
  INTEGER                        ,  INTENT(IN) :: NSOIL  !number of soil layers
  REAL(r8), DIMENSION(-NSNOW+1:    0),  INTENT(IN) :: SNICE  !snow ice mass (kg/m2)
  REAL(r8), DIMENSION(-NSNOW+1:    0),  INTENT(IN) :: SNLIQ  !snow liq mass (kg/m2) 
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL),  INTENT(IN) :: DZSNSO !snow/soil layer thickness [m]

! outputs

  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: CVSNO  !volumetric specific heat (j/m3/k)
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: TKSNO  !thermal conductivity (w/m/k)
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNICEV !partial volume of ice [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: SNLIQV !partial volume of liquid water [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(OUT) :: EPORE  !effective porosity [m3/m3]

! locals

  INTEGER :: IZ
  REAL(r8), DIMENSION(-NSNOW+1:    0) :: BDSNOI  !bulk density of snow(kg/m3)

!---------------------------------------------------------------------------------------------------
! thermal capacity of snow

  DO IZ = ISNOW+1, 0
      SNICEV(IZ)   = MIN(1., SNICE(IZ)/(DZSNSO(IZ)*DENICE) )
      EPORE(IZ)    = 1. - SNICEV(IZ)
      SNLIQV(IZ)   = MIN(EPORE(IZ),SNLIQ(IZ)/(DZSNSO(IZ)*DENH2O))
  ENDDO

  DO IZ = ISNOW+1, 0
      BDSNOI(IZ) = (SNICE(IZ)+SNLIQ(IZ))/DZSNSO(IZ)
      CVSNO(IZ) = CICE*SNICEV(IZ)+CWAT*SNLIQV(IZ)
!      CVSNO(IZ) = 0.525E06                          ! constant
  enddo

! thermal conductivity of snow

  DO IZ = ISNOW+1, 0
     TKSNO(IZ) = 3.2217E-6*BDSNOI(IZ)**2.           ! Stieglitz(yen,1965)
!    TKSNO(IZ) = 2E-2+2.5E-6*BDSNOI(IZ)*BDSNOI(IZ)   ! Anderson, 1976
!    TKSNO(IZ) = 0.35                                ! constant
!    TKSNO(IZ) = 2.576E-6*BDSNOI(IZ)**2. + 0.074    ! Verseghy (1991)
!    TKSNO(IZ) = 2.22*(BDSNOI(IZ)/1000.)**1.88      ! Douvill(Yen, 1981)
  ENDDO

  END SUBROUTINE CSNOW

!== begin tdfcnd ===================================================================================

  SUBROUTINE TDFCND (parameters, ISOIL, DF, SMC, SH2O)
! --------------------------------------------------------------------------------------------------
! Calculate thermal diffusivity and conductivity of the soil.
! Peters-Lidard approach (Peters-Lidard et al., 1998)
! --------------------------------------------------------------------------------------------------
! Code history:
! June 2001 changes: frozen soil condition.
! --------------------------------------------------------------------------------------------------
    IMPLICIT NONE
  type (noahmp_parameters), intent(in) :: parameters
    INTEGER, INTENT(IN)    :: ISOIL  ! soil layer
    REAL(r8), INTENT(IN)       :: SMC    ! total soil water
    REAL(r8), INTENT(IN)       :: SH2O   ! liq. soil water
    REAL(r8), INTENT(OUT)      :: DF     ! thermal diffusivity

! local variables
    REAL(r8)  :: AKE
    REAL(r8)  :: GAMMD
    REAL(r8)  :: THKDRY
    REAL(r8)  :: THKO     ! thermal conductivity for other soil components         
    REAL(r8)  :: THKQTZ   ! thermal conductivity for quartz
    REAL(r8)  :: THKSAT   ! 
    REAL(r8)  :: THKS     ! thermal conductivity for the solids
    REAL(r8)  :: THKW     ! water thermal conductivity
    REAL(r8)  :: SATRATIO
    REAL(r8)  :: XU
    REAL(r8)  :: XUNFROZ
! --------------------------------------------------------------------------------------------------
! We now get quartz as an input argument (set in routine redprm):
!      DATA QUARTZ /0.82, 0.10, 0.25, 0.60, 0.52,
!     &             0.35, 0.60, 0.40, 0.82/
! --------------------------------------------------------------------------------------------------
! If the soil has any moisture content compute a partial sum/product
! otherwise use a constant value which works well with most soils
! --------------------------------------------------------------------------------------------------
!  QUARTZ ....QUARTZ CONTENT (SOIL TYPE DEPENDENT)
! --------------------------------------------------------------------------------------------------
! USE AS IN PETERS-LIDARD, 1998 (MODIF. FROM JOHANSEN, 1975).

!                                  PABLO GRUNMANN, 08/17/98
! Refs.:
!      Farouki, O.T.,1986: Thermal properties of soils. Series on Rock
!              and Soil Mechanics, Vol. 11, Trans Tech, 136 pp.
!      Johansen, O., 1975: Thermal conductivity of soils. PH.D. Thesis,
!              University of Trondheim,
!      Peters-Lidard, C. D., et al., 1998: The effect of soil thermal
!              conductivity parameterization on surface energy fluxes
!              and temperatures. Journal of The Atmospheric Sciences,
!              Vol. 55, pp. 1209-1224.
! --------------------------------------------------------------------------------------------------
! NEEDS PARAMETERS
! POROSITY(SOIL TYPE):
!      POROS = SMCMAX
! SATURATION RATIO:
! PARAMETERS  W/(M.K)
    SATRATIO = SMC / parameters%SMCMAX(ISOIL)
    THKW = 0.57
!      IF (QUARTZ .LE. 0.2) THKO = 3.0
    THKO = 2.0
! SOLIDS' CONDUCTIVITY
! QUARTZ' CONDUCTIVITY
    THKQTZ = 7.7

! UNFROZEN FRACTION (FROM 1., i.e., 100%LIQUID, TO 0. (100% FROZEN))
    THKS = (THKQTZ ** parameters%QUARTZ(ISOIL))* (THKO ** (1. - parameters%QUARTZ(ISOIL)))

! UNFROZEN VOLUME FOR SATURATION (POROSITY*XUNFROZ)
    XUNFROZ = 1.0                       ! Prevent divide by zero (suggested by D. Mocko)
    IF(SMC > 0.) XUNFROZ = SH2O / SMC
! SATURATED THERMAL CONDUCTIVITY
    XU = XUNFROZ * parameters%SMCMAX(ISOIL)

! DRY DENSITY IN KG/M3
    THKSAT = THKS ** (1. - parameters%SMCMAX(ISOIL))* TKICE ** (parameters%SMCMAX(ISOIL) - XU)* THKW **   &
         (XU)

! DRY THERMAL CONDUCTIVITY IN W.M-1.K-1
    GAMMD = (1. - parameters%SMCMAX(ISOIL))*2700.

    THKDRY = (0.135* GAMMD+ 64.7)/ (2700. - 0.947* GAMMD)
! FROZEN
    IF ( (SH2O + 0.0005) <  SMC ) THEN
       AKE = SATRATIO
! UNFROZEN
! RANGE OF VALIDITY FOR THE KERSTEN NUMBER (AKE)
    ELSE

! KERSTEN NUMBER (USING "FINE" FORMULA, VALID FOR SOILS CONTAINING AT
! LEAST 5% OF PARTICLES WITH DIAMETER LESS THAN 2.E-6 METERS.)
! (FOR "COARSE" FORMULA, SEE PETERS-LIDARD ET AL., 1998).

       IF ( SATRATIO >  0.1 ) THEN

          AKE = LOG10 (SATRATIO) + 1.0

! USE K = KDRY
       ELSE

          AKE = 0.0
       END IF
!  THERMAL CONDUCTIVITY

    END IF

    DF = AKE * (THKSAT - THKDRY) + THKDRY


  end subroutine TDFCND

!== begin radiation ================================================================================

  SUBROUTINE RADIATION (parameters,VEGTYP  ,IST     ,ICE     ,NSOIL   , & !in
                        SNEQVO  ,SNEQV   ,DT      ,COSZ    ,SNOWH   , & !in
                        TG      ,TV      ,FSNO    ,QSNOW   ,FWET    , & !in
                        ELAI    ,ESAI    ,SMC     ,SOLAD   ,SOLAI   , & !in
                        FVEG    ,ILOC    ,JLOC    ,                   & !in
                        ALBOLD  ,TAUSS   ,                            & !inout
                        FSUN    ,LAISUN  ,LAISHA  ,PARSUN  ,PARSHA  , & !out
                        SAV     ,SAG     ,FSR     ,FSA     ,FSRV    , &
                        FSRG    ,BGAP    ,WGAP,                       &
                        ALBD_OUT, ALBI_OUT)            !out
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER, INTENT(IN)                  :: ILOC
  INTEGER, INTENT(IN)                  :: JLOC
  INTEGER, INTENT(IN)                  :: VEGTYP !vegetation type
  INTEGER, INTENT(IN)                  :: IST    !surface type
  INTEGER, INTENT(IN)                  :: ICE    !ice (ice = 1)
  INTEGER, INTENT(IN)                  :: NSOIL  !number of soil layers

  REAL(r8), INTENT(IN)                     :: DT     !time step [s]
  REAL(r8), INTENT(IN)                     :: QSNOW  !snowfall (mm/s)
  REAL(r8), INTENT(IN)                     :: SNEQVO !snow mass at last time step(mm)
  REAL(r8), INTENT(IN)                     :: SNEQV  !snow mass (mm)
  REAL(r8), INTENT(IN)                     :: SNOWH  !snow height (mm)
  REAL(r8), INTENT(IN)                     :: COSZ   !cosine solar zenith angle (0-1)
  REAL(r8), INTENT(IN)                     :: TG     !ground temperature (k)
  REAL(r8), INTENT(IN)                     :: TV     !vegetation temperature (k)
  REAL(r8), INTENT(IN)                     :: ELAI   !LAI, one-sided, adjusted for burying by snow
  REAL(r8), INTENT(IN)                     :: ESAI   !SAI, one-sided, adjusted for burying by snow
  REAL(r8), INTENT(IN)                     :: FWET   !fraction of canopy that is wet
  REAL(r8), DIMENSION(1:NSOIL), INTENT(IN) :: SMC    !volumetric soil water [m3/m3]
  REAL(r8), DIMENSION(1:2)    , INTENT(IN) :: SOLAD  !incoming direct solar radiation (w/m2)
  REAL(r8), DIMENSION(1:2)    , INTENT(IN) :: SOLAI  !incoming diffuse solar radiation (w/m2)
  REAL(r8), INTENT(IN)                     :: FSNO   !snow cover fraction (-)
  REAL(r8), INTENT(IN)                     :: FVEG   !green vegetation fraction [0.0-1.0]

! inout
  REAL(r8),                  INTENT(INOUT) :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(r8),                  INTENT(INOUT) :: TAUSS  !non-dimensional snow age.

! output
  REAL(r8), INTENT(OUT)                    :: FSUN   !sunlit fraction of canopy (-)
  REAL(r8), INTENT(OUT)                    :: LAISUN !sunlit leaf area (-)
  REAL(r8), INTENT(OUT)                    :: LAISHA !shaded leaf area (-)
  REAL(r8), INTENT(OUT)                    :: PARSUN !average absorbed par for sunlit leaves (w/m2)
  REAL(r8), INTENT(OUT)                    :: PARSHA !average absorbed par for shaded leaves (w/m2)
  REAL(r8), INTENT(OUT)                    :: SAV    !solar radiation absorbed by vegetation (w/m2)
  REAL(r8), INTENT(OUT)                    :: SAG    !solar radiation absorbed by ground (w/m2)
  REAL(r8), INTENT(OUT)                    :: FSA    !total absorbed solar radiation (w/m2)
  REAL(r8), INTENT(OUT)                    :: FSR    !total reflected solar radiation (w/m2)

!jref:start  
  REAL(r8), INTENT(OUT)                    :: FSRV    !veg. reflected solar radiation (w/m2)
  REAL(r8), INTENT(OUT)                    :: FSRG    !ground reflected solar radiation (w/m2)
  REAL(r8), INTENT(OUT)                    :: BGAP
  REAL(r8), INTENT(OUT)                    :: WGAP
!jref:end  

! local
  REAL(r8)                                 :: FAGE   !snow age function (0 - new snow)
  REAL(r8), DIMENSION(1:2)                 :: ALBGRD !ground albedo (direct)
  REAL(r8), DIMENSION(1:2)                 :: ALBGRI !ground albedo (diffuse)
  REAL(r8), DIMENSION(1:2)                 :: ALBD   !surface albedo (direct)
  REAL(r8), DIMENSION(1:2)                 :: ALBI   !surface albedo (diffuse)
  REAL(r8), DIMENSION(1:2),intent(out)     :: ALBD_OUT   !cheyz for coupling
  REAL(r8), DIMENSION(1:2),intent(out)     :: ALBI_OUT   !cheyz for coupling
  REAL(r8), DIMENSION(1:2)                 :: FABD   !flux abs by veg (per unit direct flux)
  REAL(r8), DIMENSION(1:2)                 :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(r8), DIMENSION(1:2)                 :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(r8), DIMENSION(1:2)                 :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(r8), DIMENSION(1:2)                 :: FTII   !down diffuse flux below veg (per unit dif flux)
!jref:start  
  REAL(r8), DIMENSION(1:2)                 :: FREVI
  REAL(r8), DIMENSION(1:2)                 :: FREVD
  REAL(r8), DIMENSION(1:2)                 :: FREGI
  REAL(r8), DIMENSION(1:2)                 :: FREGD
!jref:end

  REAL(r8)                                 :: FSHA   !shaded fraction of canopy
  REAL(r8)                                 :: VAI    !total LAI + stem area index, one sided

  REAL(r8),PARAMETER :: MPE = 1.E-6
  LOGICAL VEG  !true: vegetated for surface temperature calculation

! --------------------------------------------------------------------------------------------------

! surface abeldo

   CALL ALBEDO (parameters,VEGTYP ,IST    ,ICE    ,NSOIL  , & !in
                DT     ,COSZ   ,FAGE   ,ELAI   ,ESAI   , & !in
                TG     ,TV     ,SNOWH  ,FSNO   ,FWET   , & !in
                SMC    ,SNEQVO ,SNEQV  ,QSNOW  ,FVEG   , & !in
                ILOC   ,JLOC   ,                         & !in
                ALBOLD ,TAUSS                          , & !inout
                ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
                FABI   ,FTDD   ,FTID   ,FTII   ,FSUN   , & !)   !out
                FREVI  ,FREVD   ,FREGD ,FREGI  ,BGAP   , & !inout
                WGAP)

!-> cheyz and xinyf 
   albd_out = albd
   albi_out = albi
!<-
! surface radiation

     FSHA = 1.-FSUN
     LAISUN = ELAI*FSUN
     LAISHA = ELAI*FSHA
     VAI = ELAI+ ESAI
     IF (VAI .GT. 0.) THEN
        VEG = .TRUE.
     ELSE
        VEG = .FALSE.
     END IF

   CALL SURRAD (parameters,MPE    ,FSUN   ,FSHA   ,ELAI   ,VAI    , & !in
                LAISUN ,LAISHA ,SOLAD  ,SOLAI  ,FABD   , & !in
                FABI   ,FTDD   ,FTID   ,FTII   ,ALBGRD , & !in
                ALBGRI ,ALBD   ,ALBI   ,ILOC   ,JLOC   , & !in
                PARSUN ,PARSHA ,SAV    ,SAG    ,FSA    , & !out
                FSR    ,                                 & !out
                FREVI  ,FREVD  ,FREGD  ,FREGI  ,FSRV   , & !inout
                FSRG)

  END SUBROUTINE RADIATION

!== begin albedo ===================================================================================

  SUBROUTINE ALBEDO (parameters,VEGTYP ,IST    ,ICE    ,NSOIL  , & !in
                     DT     ,COSZ   ,FAGE   ,ELAI   ,ESAI   , & !in
                     TG     ,TV     ,SNOWH  ,FSNO   ,FWET   , & !in
                     SMC    ,SNEQVO ,SNEQV  ,QSNOW  ,FVEG   , & !in
                     ILOC   ,JLOC   ,                         & !in
                     ALBOLD ,TAUSS                          , & !inout
                     ALBGRD ,ALBGRI ,ALBD   ,ALBI   ,FABD   , & !out
                     FABI   ,FTDD   ,FTID   ,FTII   ,FSUN   , & !out
                     FREVI  ,FREVD  ,FREGD  ,FREGI  ,BGAP   , & !out
                     WGAP)

! --------------------------------------------------------------------------------------------------
! surface albedos. also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation.
! also sunlit fraction of the canopy.
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                  INTENT(IN)  :: ILOC
  INTEGER,                  INTENT(IN)  :: JLOC
  INTEGER,                  INTENT(IN)  :: NSOIL  !number of soil layers
  INTEGER,                  INTENT(IN)  :: VEGTYP !vegetation type
  INTEGER,                  INTENT(IN)  :: IST    !surface type
  INTEGER,                  INTENT(IN)  :: ICE    !ice (ice = 1)

  REAL(r8),                     INTENT(IN)  :: DT     !time step [sec]
  REAL(r8),                     INTENT(IN)  :: QSNOW  !snowfall
  REAL(r8),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle for next time step
  REAL(r8),                     INTENT(IN)  :: SNOWH  !snow height (mm)
  REAL(r8),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(r8),                     INTENT(IN)  :: TV     !vegetation temperature (k)
  REAL(r8),                     INTENT(IN)  :: ELAI   !LAI, one-sided, adjusted for burying by snow
  REAL(r8),                     INTENT(IN)  :: ESAI   !SAI, one-sided, adjusted for burying by snow
  REAL(r8),                     INTENT(IN)  :: FSNO   !fraction of grid covered by snow
  REAL(r8),                     INTENT(IN)  :: FWET   !fraction of canopy that is wet
  REAL(r8),                     INTENT(IN)  :: SNEQVO !snow mass at last time step(mm)
  REAL(r8),                     INTENT(IN)  :: SNEQV  !snow mass (mm)
  REAL(r8),                     INTENT(IN)  :: FVEG   !green vegetation fraction [0.0-1.0]
  REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water (m3/m3)

! inout
  REAL(r8),                  INTENT(INOUT)  :: ALBOLD !snow albedo at last time step (CLASS type)
  REAL(r8),                  INTENT(INOUT)  :: TAUSS  !non-dimensional snow age

! output
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: ALBD   !surface albedo (direct)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: ALBI   !surface albedo (diffuse)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FABD   !flux abs by veg (per unit direct flux)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FABI   !flux abs by veg (per unit diffuse flux)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FTDD   !down direct flux below veg (per unit dir flux)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FTID   !down diffuse flux below veg (per unit dir flux)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FTII   !down diffuse flux below veg (per unit dif flux)
  REAL(r8),                     INTENT(OUT) :: FSUN   !sunlit fraction of canopy (-)
!jref:start
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FREVD
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FREVI
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FREGD
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: FREGI
  REAL(r8), INTENT(OUT) :: BGAP
  REAL(r8), INTENT(OUT) :: WGAP
!jref:end

! ------------------------------------------------------------------------
! ------------------------ local variables -------------------------------
! local
  REAL(r8)                 :: FAGE     !snow age function
  REAL(r8)                 :: ALB
  INTEGER              :: IB       !indices
  INTEGER              :: NBAND    !number of solar radiation wave bands
  INTEGER              :: IC       !direct beam: ic=0; diffuse: ic=1

  REAL(r8)                 :: WL       !fraction of LAI+SAI that is LAI
  REAL(r8)                 :: WS       !fraction of LAI+SAI that is SAI
  REAL(r8)                 :: MPE      !prevents overflow for division by zero

  REAL(r8), DIMENSION(1:2) :: RHO      !leaf/stem reflectance weighted by fraction LAI and SAI
  REAL(r8), DIMENSION(1:2) :: TAU      !leaf/stem transmittance weighted by fraction LAI and SAI
  REAL(r8), DIMENSION(1:2) :: FTDI     !down direct flux below veg per unit dif flux = 0
  REAL(r8), DIMENSION(1:2) :: ALBSND   !snow albedo (direct)
  REAL(r8), DIMENSION(1:2) :: ALBSNI   !snow albedo (diffuse)

  REAL(r8)                 :: VAI      !ELAI+ESAI
  REAL(r8)                 :: GDIR     !average projected leaf/stem area in solar direction
  REAL(r8)                 :: EXT      !optical depth direct beam per unit leaf + stem area

! --------------------------------------------------------------------------------------------------

  NBAND = 2
  MPE = 1.E-06
  BGAP = 0.
  WGAP = 0.

! initialize output because solar radiation only done if COSZ > 0

  DO IB = 1, NBAND
    ALBD(IB) = 0.
    ALBI(IB) = 0.
    ALBGRD(IB) = 0.
    ALBGRI(IB) = 0.
    FABD(IB) = 0.
    FABI(IB) = 0.
    FTDD(IB) = 0.
    FTID(IB) = 0.
    FTII(IB) = 0.
    IF (IB.EQ.1) FSUN = 0.
  END DO

  IF(COSZ <= 0) GOTO 100

! weight reflectance/transmittance by LAI and SAI

  DO IB = 1, NBAND
    VAI = ELAI + ESAI
    WL  = ELAI / MAX(VAI,MPE)
    WS  = ESAI / MAX(VAI,MPE)
    RHO(IB) = MAX(parameters%RHOL(IB)*WL+parameters%RHOS(IB)*WS, MPE)
    TAU(IB) = MAX(parameters%TAUL(IB)*WL+parameters%TAUS(IB)*WS, MPE)
  END DO

! snow age

   CALL SNOW_AGE (parameters,DT,TG,SNEQVO,SNEQV,TAUSS,FAGE)

! snow albedos: only if COSZ > 0 and FSNO > 0

  IF(OPT_ALB == 1) &
     CALL SNOWALB_BATS (parameters,NBAND, FSNO,COSZ,FAGE,ALBSND,ALBSNI)
  IF(OPT_ALB == 2) THEN
     CALL SNOWALB_CLASS (parameters,NBAND,QSNOW,DT,ALB,ALBOLD,ALBSND,ALBSNI,ILOC,JLOC)
     ALBOLD = ALB
  END IF

! ground surface albedo

  CALL GROUNDALB (parameters,NSOIL   ,NBAND   ,ICE     ,IST     , & !in
                  FSNO    ,SMC     ,ALBSND  ,ALBSNI  ,COSZ    , & !in
                  TG      ,ILOC    ,JLOC    ,                   & !in
                  ALBGRD  ,ALBGRI  )                              !out

! loop over NBAND wavebands to calculate surface albedos and solar
! fluxes for unit incoming direct (IC=0) and diffuse flux (IC=1)

  DO IB = 1, NBAND
      IC = 0      ! direct
      CALL TWOSTREAM (parameters,IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                      FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                      TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                      FABD   ,ALBD    ,FTDD    ,FTID    ,GDIR   , &!)   !out
                      FREVD  ,FREGD   ,BGAP    ,WGAP)

      IC = 1      ! diffuse
      CALL TWOSTREAM (parameters,IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                      FWET   ,TV      ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                      TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                      FABI   ,ALBI    ,FTDI    ,FTII    ,GDIR   , & !)   !out
                      FREVI  ,FREGI   ,BGAP    ,WGAP)

  END DO

! sunlit fraction of canopy. set FSUN = 0 if FSUN < 0.01.

  EXT = GDIR/COSZ * SQRT(1.-RHO(1)-TAU(1))
  FSUN = (1.-EXP(-EXT*VAI)) / MAX(EXT*VAI,MPE)
  EXT = FSUN

  IF (EXT .LT. 0.01) THEN
     WL = 0.
  ELSE
     WL = EXT 
  END IF
  FSUN = WL

100 CONTINUE

  END SUBROUTINE ALBEDO

!== begin surrad ===================================================================================

  SUBROUTINE SURRAD (parameters,MPE     ,FSUN    ,FSHA    ,ELAI    ,VAI     , & !in
                     LAISUN  ,LAISHA  ,SOLAD   ,SOLAI   ,FABD    , & !in
                     FABI    ,FTDD    ,FTID    ,FTII    ,ALBGRD  , & !in
                     ALBGRI  ,ALBD    ,ALBI    ,ILOC    ,JLOC    , & !in
                     PARSUN  ,PARSHA  ,SAV     ,SAG     ,FSA     , & !out
                     FSR     , & !)                                       !out
                     FREVI   ,FREVD   ,FREGD   ,FREGI   ,FSRV    , &
                     FSRG) !inout

! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER, INTENT(IN)              :: ILOC
  INTEGER, INTENT(IN)              :: JLOC
  REAL(r8), INTENT(IN)                 :: MPE     !prevents underflow errors if division by zero

  REAL(r8), INTENT(IN)                 :: FSUN    !sunlit fraction of canopy
  REAL(r8), INTENT(IN)                 :: FSHA    !shaded fraction of canopy
  REAL(r8), INTENT(IN)                 :: ELAI    !leaf area, one-sided
  REAL(r8), INTENT(IN)                 :: VAI     !leaf + stem area, one-sided
  REAL(r8), INTENT(IN)                 :: LAISUN  !sunlit leaf area index, one-sided
  REAL(r8), INTENT(IN)                 :: LAISHA  !shaded leaf area index, one-sided

  REAL(r8), DIMENSION(1:2), INTENT(IN) :: SOLAD   !incoming direct solar radiation (w/m2)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: SOLAI   !incoming diffuse solar radiation (w/m2)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FABD    !flux abs by veg (per unit incoming direct flux)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FABI    !flux abs by veg (per unit incoming diffuse flux)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FTDD    !down dir flux below veg (per incoming dir flux)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FTID    !down dif flux below veg (per incoming dir flux)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FTII    !down dif flux below veg (per incoming dif flux)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: ALBGRD  !ground albedo (direct)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: ALBGRI  !ground albedo (diffuse)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: ALBD    !overall surface albedo (direct)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: ALBI    !overall surface albedo (diffuse)

  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FREVD    !overall surface albedo veg (direct)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FREVI    !overall surface albedo veg (diffuse)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FREGD    !overall surface albedo grd (direct)
  REAL(r8), DIMENSION(1:2), INTENT(IN) :: FREGI    !overall surface albedo grd (diffuse)

! output

  REAL(r8), INTENT(OUT)                :: PARSUN  !average absorbed par for sunlit leaves (w/m2)
  REAL(r8), INTENT(OUT)                :: PARSHA  !average absorbed par for shaded leaves (w/m2)
  REAL(r8), INTENT(OUT)                :: SAV     !solar radiation absorbed by vegetation (w/m2)
  REAL(r8), INTENT(OUT)                :: SAG     !solar radiation absorbed by ground (w/m2)
  REAL(r8), INTENT(OUT)                :: FSA     !total absorbed solar radiation (w/m2)
  REAL(r8), INTENT(OUT)                :: FSR     !total reflected solar radiation (w/m2)
  REAL(r8), INTENT(OUT)                :: FSRV    !reflected solar radiation by vegetation
  REAL(r8), INTENT(OUT)                :: FSRG    !reflected solar radiation by ground

! ------------------------ local variables ----------------------------------------------------
  INTEGER                          :: IB      !waveband number (1=vis, 2=nir)
  INTEGER                          :: NBAND   !number of solar radiation waveband classes

  REAL(r8)                             :: ABS     !absorbed solar radiation (w/m2)
  REAL(r8)                             :: RNIR    !reflected solar radiation [nir] (w/m2)
  REAL(r8)                             :: RVIS    !reflected solar radiation [vis] (w/m2)
  REAL(r8)                             :: LAIFRA  !leaf area fraction of canopy
  REAL(r8)                             :: TRD     !transmitted solar radiation: direct (w/m2)
  REAL(r8)                             :: TRI     !transmitted solar radiation: diffuse (w/m2)
  REAL(r8), DIMENSION(1:2)             :: CAD     !direct beam absorbed by canopy (w/m2)
  REAL(r8), DIMENSION(1:2)             :: CAI     !diffuse radiation absorbed by canopy (w/m2)
! ---------------------------------------------------------------------------------------------
   NBAND = 2

! zero summed solar fluxes

    SAG = 0.
    SAV = 0.
    FSA = 0.

! loop over nband wavebands

  DO IB = 1, NBAND

! absorbed by canopy

    CAD(IB) = SOLAD(IB)*FABD(IB)    
    CAI(IB) = SOLAI(IB)*FABI(IB)
    SAV     = SAV + CAD(IB) + CAI(IB)
    FSA     = FSA + CAD(IB) + CAI(IB)
 
! transmitted solar fluxes incident on ground

    TRD = SOLAD(IB)*FTDD(IB)
    TRI = SOLAD(IB)*FTID(IB) + SOLAI(IB)*FTII(IB)

! solar radiation absorbed by ground surface

    ABS = TRD*(1.-ALBGRD(IB)) + TRI*(1.-ALBGRI(IB))
    SAG = SAG + ABS
    FSA = FSA + ABS
  END DO

! partition visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves

     LAIFRA = ELAI / MAX(VAI,MPE)
     IF (FSUN .GT. 0.) THEN
        PARSUN = (CAD(1)+FSUN*CAI(1)) * LAIFRA / MAX(LAISUN,MPE)
        PARSHA = (FSHA*CAI(1))*LAIFRA / MAX(LAISHA,MPE)
     ELSE
        PARSUN = 0.
        PARSHA = (CAD(1)+CAI(1))*LAIFRA /MAX(LAISHA,MPE)
     ENDIF

! reflected solar radiation

     RVIS = ALBD(1)*SOLAD(1) + ALBI(1)*SOLAI(1)
     RNIR = ALBD(2)*SOLAD(2) + ALBI(2)*SOLAI(2)
     FSR  = RVIS + RNIR

! reflected solar radiation of veg. and ground (combined ground)
     FSRV = FREVD(1)*SOLAD(1)+FREVI(1)*SOLAI(1)+FREVD(2)*SOLAD(2)+FREVI(2)*SOLAI(2)
     FSRG = FREGD(1)*SOLAD(1)+FREGI(1)*SOLAI(1)+FREGD(2)*SOLAD(2)+FREGI(2)*SOLAI(2)


  END SUBROUTINE SURRAD

!== begin snow_age =================================================================================

  SUBROUTINE SNOW_AGE (parameters,DT,TG,SNEQVO,SNEQV,TAUSS,FAGE)
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ------------------------ code history ------------------------------------------------------------
! from BATS
! ------------------------ input/output variables --------------------------------------------------
!input
  type (noahmp_parameters), intent(in) :: parameters
   REAL(r8), INTENT(IN) :: DT        !main time step (s)
   REAL(r8), INTENT(IN) :: TG        !ground temperature (k)
   REAL(r8), INTENT(IN) :: SNEQVO    !snow mass at last time step(mm)
   REAL(r8), INTENT(IN) :: SNEQV     !snow water per unit ground area (mm)

!output
   REAL(r8), INTENT(OUT) :: FAGE     !snow age

!input/output
   REAL(r8), INTENT(INOUT) :: TAUSS      !non-dimensional snow age
!local
   REAL(r8)            :: TAGE       !total aging effects
   REAL(r8)            :: AGE1       !effects of grain growth due to vapor diffusion
   REAL(r8)            :: AGE2       !effects of grain growth at freezing of melt water
   REAL(r8)            :: AGE3       !effects of soot
   REAL(r8)            :: DELA       !temporary variable
   REAL(r8)            :: SGE        !temporary variable
   REAL(r8)            :: DELS       !temporary variable
   REAL(r8)            :: DELA0      !temporary variable
   REAL(r8)            :: ARG        !temporary variable
! See Yang et al. (1997) J.of Climate for detail.
!---------------------------------------------------------------------------------------------------

   IF(SNEQV.LE.0.0) THEN
          TAUSS = 0.
   ELSE IF (SNEQV.GT.800.) THEN
          TAUSS = 0.
   ELSE
          DELA0 = 1.E-6*DT
          ARG   = 5.E3*(1./TFRZ-1./TG)
          AGE1  = EXP(ARG)
          AGE2  = EXP(AMIN1(0.,10.*ARG))
          AGE3  = 0.3
          TAGE  = AGE1+AGE2+AGE3
          DELA  = DELA0*TAGE
          DELS  = AMAX1(0.0,SNEQV-SNEQVO) / parameters%SWEMX
          SGE   = (TAUSS+DELA)*(1.0-DELS)
          TAUSS = AMAX1(0.,SGE)
   ENDIF

   FAGE= TAUSS/(TAUSS+1.)

  END SUBROUTINE SNOW_AGE

!== begin snowalb_bats =============================================================================

  SUBROUTINE SNOWALB_BATS (parameters,NBAND,FSNO,COSZ,FAGE,ALBSND,ALBSNI)
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,INTENT(IN) :: NBAND  !number of waveband classes

  REAL(r8),INTENT(IN) :: COSZ    !cosine solar zenith angle
  REAL(r8),INTENT(IN) :: FSNO    !snow cover fraction (-)
  REAL(r8),INTENT(IN) :: FAGE    !snow age correction

! output

  REAL(r8), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(r8), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

  REAL(r8) :: FZEN                 !zenith angle correction
  REAL(r8) :: CF1                  !temperary variable
  REAL(r8) :: SL2                  !2.*SL
  REAL(r8) :: SL1                  !1/SL
  REAL(r8) :: SL                   !adjustable parameter
  REAL(r8), PARAMETER :: C1 = 0.2  !default in BATS 
  REAL(r8), PARAMETER :: C2 = 0.5  !default in BATS
!  REAL(r8), PARAMETER :: C1 = 0.2 * 2. ! double the default to match Sleepers River's
!  REAL(r8), PARAMETER :: C2 = 0.5 * 2. ! snow surface albedo (double aging effects)
! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

        SL=2.0
        SL1=1./SL
        SL2=2.*SL
        CF1=((1.+SL1)/(1.+SL2*COSZ)-SL1)
        FZEN=AMAX1(CF1,0.)

        ALBSNI(1)=0.95*(1.-C1*FAGE)         
        ALBSNI(2)=0.65*(1.-C2*FAGE)        

        ALBSND(1)=ALBSNI(1)+0.4*FZEN*(1.-ALBSNI(1))    !  vis direct
        ALBSND(2)=ALBSNI(2)+0.4*FZEN*(1.-ALBSNI(2))    !  nir direct

  END SUBROUTINE SNOWALB_BATS

!== begin snowalb_class ============================================================================

  SUBROUTINE SNOWALB_CLASS (parameters,NBAND,QSNOW,DT,ALB,ALBOLD,ALBSND,ALBSNI,ILOC,JLOC)
! ----------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,INTENT(IN) :: ILOC !grid index
  INTEGER,INTENT(IN) :: JLOC !grid index
  INTEGER,INTENT(IN) :: NBAND  !number of waveband classes

  REAL(r8),INTENT(IN) :: QSNOW     !snowfall (mm/s)
  REAL(r8),INTENT(IN) :: DT        !time step (sec)
  REAL(r8),INTENT(IN) :: ALBOLD    !snow albedo at last time step

! in & out

  REAL(r8),                INTENT(INOUT) :: ALB        ! 
! output

  REAL(r8), DIMENSION(1:2),INTENT(OUT) :: ALBSND !snow albedo for direct(1=vis, 2=nir)
  REAL(r8), DIMENSION(1:2),INTENT(OUT) :: ALBSNI !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  INTEGER :: IB          !waveband class

! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        ALBSND(1: NBAND) = 0.
        ALBSNI(1: NBAND) = 0.

! when cosz > 0

         ALB = 0.55 + (ALBOLD-0.55) * EXP(-0.01*DT/3600.)

! 1 mm fresh snow(SWE) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
! here assume 1cm snow depth will fully cover the old snow

         IF (QSNOW > 0.) then
           ALB = ALB + MIN(QSNOW,parameters%SWEMX/DT) * (0.84-ALB)/(parameters%SWEMX/DT)
         ENDIF

         ALBSNI(1)= ALB         ! vis diffuse
         ALBSNI(2)= ALB         ! nir diffuse
         ALBSND(1)= ALB         ! vis direct
         ALBSND(2)= ALB         ! nir direct

  END SUBROUTINE SNOWALB_CLASS

!== begin groundalb ================================================================================

  SUBROUTINE GROUNDALB (parameters,NSOIL   ,NBAND   ,ICE     ,IST     , & !in
                        FSNO    ,SMC     ,ALBSND  ,ALBSNI  ,COSZ    , & !in
                        TG      ,ILOC    ,JLOC    ,                   & !in
                        ALBGRD  ,ALBGRI  )                              !out
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                  INTENT(IN)  :: ILOC   !grid index
  INTEGER,                  INTENT(IN)  :: JLOC   !grid index
  INTEGER,                  INTENT(IN)  :: NSOIL  !number of soil layers
  INTEGER,                  INTENT(IN)  :: NBAND  !number of solar radiation waveband classes
  INTEGER,                  INTENT(IN)  :: ICE    !value of ist for land ice
  INTEGER,                  INTENT(IN)  :: IST    !surface type
  REAL(r8),                     INTENT(IN)  :: FSNO   !fraction of surface covered with snow (-)
  REAL(r8),                     INTENT(IN)  :: TG     !ground temperature (k)
  REAL(r8),                     INTENT(IN)  :: COSZ   !cosine solar zenith angle (0-1)
  REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC    !volumetric soil water content (m3/m3)
  REAL(r8), DIMENSION(1:    2), INTENT(IN)  :: ALBSND !direct beam snow albedo (vis, nir)
  REAL(r8), DIMENSION(1:    2), INTENT(IN)  :: ALBSNI !diffuse snow albedo (vis, nir)

!output

  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: ALBGRD !ground albedo (direct beam: vis, nir)
  REAL(r8), DIMENSION(1:    2), INTENT(OUT) :: ALBGRI !ground albedo (diffuse: vis, nir)

!local 

  INTEGER                               :: IB     !waveband number (1=vis, 2=nir)
  REAL(r8)                                  :: INC    !soil water correction factor for soil albedo
  REAL(r8)                                  :: ALBSOD !soil albedo (direct)
  REAL(r8)                                  :: ALBSOI !soil albedo (diffuse)
! --------------------------------------------------------------------------------------------------

  DO IB = 1, NBAND
        INC = MAX(0.11-0.40*SMC(1), 0.)
        IF (IST .EQ. 1)  THEN                     !soil
           ALBSOD = MIN(parameters%ALBSAT(IB)+INC,parameters%ALBDRY(IB))
           ALBSOI = ALBSOD
        ELSE IF (TG .GT. TFRZ) THEN               !unfrozen lake, wetland
           ALBSOD = 0.06/(MAX(0.01,COSZ)**1.7 + 0.15)
           ALBSOI = 0.06
        ELSE                                      !frozen lake, wetland
           ALBSOD = parameters%ALBLAK(IB)
           ALBSOI = ALBSOD
        END IF

! increase desert and semi-desert albedos

!        IF (IST .EQ. 1 .AND. ISC .EQ. 9) THEN
!           ALBSOD = ALBSOD + 0.10
!           ALBSOI = ALBSOI + 0.10
!        end if

        ALBGRD(IB) = ALBSOD*(1.-FSNO) + ALBSND(IB)*FSNO
        ALBGRI(IB) = ALBSOI*(1.-FSNO) + ALBSNI(IB)*FSNO
  END DO

  END SUBROUTINE GROUNDALB

!== begin twostream ================================================================================

  SUBROUTINE TWOSTREAM (parameters,IB     ,IC      ,VEGTYP  ,COSZ    ,VAI    , & !in
                        FWET   ,T       ,ALBGRD  ,ALBGRI  ,RHO    , & !in
                        TAU    ,FVEG    ,IST     ,ILOC    ,JLOC   , & !in
                        FAB    ,FRE     ,FTD     ,FTI     ,GDIR   , & !)   !out
                        FREV   ,FREG    ,BGAP    ,WGAP)

! --------------------------------------------------------------------------------------------------
! use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
   INTEGER,              INTENT(IN)  :: ILOC    !grid index
   INTEGER,              INTENT(IN)  :: JLOC    !grid index
   INTEGER,              INTENT(IN)  :: IST     !surface type
   INTEGER,              INTENT(IN)  :: IB      !waveband number
   INTEGER,              INTENT(IN)  :: IC      !0=unit incoming direct; 1=unit incoming diffuse
   INTEGER,              INTENT(IN)  :: VEGTYP  !vegetation type

   REAL(r8),                 INTENT(IN)  :: COSZ    !cosine of direct zenith angle (0-1)
   REAL(r8),                 INTENT(IN)  :: VAI     !one-sided leaf+stem area index (m2/m2)
   REAL(r8),                 INTENT(IN)  :: FWET    !fraction of lai, sai that is wetted (-)
   REAL(r8),                 INTENT(IN)  :: T       !surface temperature (k)

   REAL(r8), DIMENSION(1:2), INTENT(IN)  :: ALBGRD  !direct  albedo of underlying surface (-)
   REAL(r8), DIMENSION(1:2), INTENT(IN)  :: ALBGRI  !diffuse albedo of underlying surface (-)
   REAL(r8), DIMENSION(1:2), INTENT(IN)  :: RHO     !leaf+stem reflectance
   REAL(r8), DIMENSION(1:2), INTENT(IN)  :: TAU     !leaf+stem transmittance
   REAL(r8),                 INTENT(IN)  :: FVEG    !green vegetation fraction [0.0-1.0]

! output

   REAL(r8), DIMENSION(1:2), INTENT(OUT) :: FAB     !flux abs by veg layer (per unit incoming flux)
   REAL(r8), DIMENSION(1:2), INTENT(OUT) :: FRE     !flux refl above veg layer (per unit incoming flux)
   REAL(r8), DIMENSION(1:2), INTENT(OUT) :: FTD     !down dir flux below veg layer (per unit in flux)
   REAL(r8), DIMENSION(1:2), INTENT(OUT) :: FTI     !down dif flux below veg layer (per unit in flux)
   REAL(r8),                 INTENT(OUT) :: GDIR    !projected leaf+stem area in solar direction
   REAL(r8), DIMENSION(1:2), INTENT(OUT) :: FREV    !flux reflected by veg layer   (per unit incoming flux) 
   REAL(r8), DIMENSION(1:2), INTENT(OUT) :: FREG    !flux reflected by ground (per unit incoming flux)

! local
   REAL(r8)                              :: OMEGA   !fraction of intercepted radiation that is scattered
   REAL(r8)                              :: OMEGAL  !omega for leaves
   REAL(r8)                              :: BETAI   !upscatter parameter for diffuse radiation
   REAL(r8)                              :: BETAIL  !betai for leaves
   REAL(r8)                              :: BETAD   !upscatter parameter for direct beam radiation
   REAL(r8)                              :: BETADL  !betad for leaves
   REAL(r8)                              :: EXT     !optical depth of direct beam per unit leaf area
   REAL(r8)                              :: AVMU    !average diffuse optical depth

   REAL(r8)                              :: COSZI   !0.001 <= cosz <= 1.000
   REAL(r8)                              :: ASU     !single scattering albedo
   REAL(r8)                              :: CHIL    ! -0.4 <= xl <= 0.6

   REAL(r8)                              :: TMP0,TMP1,TMP2,TMP3,TMP4,TMP5,TMP6,TMP7,TMP8,TMP9
   REAL(r8)                              :: P1,P2,P3,P4,S1,S2,U1,U2,U3
   REAL(r8)                              :: B,C,D,D1,D2,F,H,H1,H2,H3,H4,H5,H6,H7,H8,H9,H10
   REAL(r8)                              :: PHI1,PHI2,SIGMA
   REAL(r8)                              :: FTDS,FTIS,FRES
   REAL(r8)                              :: DENFVEG
   REAL(r8)                              :: VAI_SPREAD
!jref:start
   REAL(r8)                              :: FREVEG,FREBAR,FTDVEG,FTIVEG,FTDBAR,FTIBAR
   REAL(r8)                              :: THETAZ
!jref:end   

!  variables for the modified two-stream scheme
!  Niu and Yang (2004), JGR

   REAL(r8), PARAMETER :: PAI = 3.14159265 
   REAL(r8) :: HD       !crown depth (m)
   REAL(r8) :: BB       !vertical crown radius (m)
   REAL(r8) :: THETAP   !angle conversion from SZA 
   REAL(r8) :: FA       !foliage volume density (m-1)
   REAL(r8) :: NEWVAI   !effective LSAI (-)

   REAL(r8),INTENT(INOUT) :: BGAP     !between canopy gap fraction for beam (-)
   REAL(r8),INTENT(INOUT) :: WGAP     !within canopy gap fraction for beam (-)

   REAL(r8) :: KOPEN    !gap fraction for diffue light (-)
   REAL(r8) :: GAP      !total gap fraction for beam ( <=1-shafac )

! -----------------------------------------------------------------
! compute within and between gaps
     VAI_SPREAD = VAI
     if(VAI == 0.0) THEN
         GAP     = 1.0
         KOPEN   = 1.0
     ELSE
         IF(OPT_RAD == 1) THEN
	   DENFVEG = -LOG(MAX(1.0-FVEG,0.01))/(PAI*parameters%RC**2)
           HD      = parameters%HVT - parameters%HVB
           BB      = 0.5 * HD           
           THETAP  = ATAN(BB/parameters%RC * TAN(ACOS(MAX(0.01,COSZ))) )
           ! BGAP    = EXP(-parameters%DEN * PAI * parameters%RC**2/COS(THETAP) )
           BGAP    = EXP(-DENFVEG * PAI * parameters%RC**2/COS(THETAP) )
           FA      = VAI/(1.33 * PAI * parameters%RC**3.0 *(BB/parameters%RC)*DENFVEG)
           NEWVAI  = HD*FA
           WGAP    = (1.0-BGAP) * EXP(-0.5*NEWVAI/COSZ)
           GAP     = MIN(1.0-FVEG, BGAP+WGAP)

           KOPEN   = 0.05
         END IF

         IF(OPT_RAD == 2) THEN
           GAP     = 0.0
           KOPEN   = 0.0
         END IF

         IF(OPT_RAD == 3) THEN
           GAP     = 1.0-FVEG
           KOPEN   = 1.0-FVEG
         END IF
     end if

! calculate two-stream parameters OMEGA, BETAD, BETAI, AVMU, GDIR, EXT.
! OMEGA, BETAD, BETAI are adjusted for snow. values for OMEGA*BETAD
! and OMEGA*BETAI are calculated and then divided by the new OMEGA
! because the product OMEGA*BETAI, OMEGA*BETAD is used in solution.
! also, the transmittances and reflectances (TAU, RHO) are linear
! weights of leaf and stem values.

     COSZI  = MAX(0.001, COSZ)
     CHIL   = MIN( MAX(parameters%XL, -0.4), 0.6)
     IF (ABS(CHIL) .LE. 0.01) CHIL = 0.01
     PHI1   = 0.5 - 0.633*CHIL - 0.330*CHIL*CHIL
     PHI2   = 0.877 * (1.-2.*PHI1)
     GDIR   = PHI1 + PHI2*COSZI
     EXT    = GDIR/COSZI
     AVMU   = ( 1. - PHI1/PHI2 * LOG((PHI1+PHI2)/PHI1) ) / PHI2
     OMEGAL = RHO(IB) + TAU(IB)
     TMP0   = GDIR + PHI2*COSZI
     TMP1   = PHI1*COSZI
     ASU    = 0.5*OMEGAL*GDIR/TMP0 * ( 1.-TMP1/TMP0*LOG((TMP1+TMP0)/TMP1) )
     BETADL = (1.+AVMU*EXT)/(OMEGAL*AVMU*EXT)*ASU
     BETAIL = 0.5 * ( RHO(IB)+TAU(IB) + (RHO(IB)-TAU(IB))   &
            * ((1.+CHIL)/2.)**2 ) / OMEGAL

! adjust omega, betad, and betai for intercepted snow

     IF (T .GT. TFRZ) THEN                                !no snow
        TMP0 = OMEGAL
        TMP1 = BETADL
        TMP2 = BETAIL
     ELSE
        TMP0 =   (1.-FWET)*OMEGAL        + FWET*parameters%OMEGAS(IB)
        TMP1 = ( (1.-FWET)*OMEGAL*BETADL + FWET*parameters%OMEGAS(IB)*parameters%BETADS ) / TMP0
        TMP2 = ( (1.-FWET)*OMEGAL*BETAIL + FWET*parameters%OMEGAS(IB)*parameters%BETAIS ) / TMP0
     END IF

     OMEGA = TMP0
     BETAD = TMP1
     BETAI = TMP2

! absorbed, reflected, transmitted fluxes per unit incoming radiation

     B = 1. - OMEGA + OMEGA*BETAI
     C = OMEGA*BETAI
     TMP0 = AVMU*EXT
     D = TMP0 * OMEGA*BETAD
     F = TMP0 * OMEGA*(1.-BETAD)
     TMP1 = B*B - C*C
     H = SQRT(TMP1) / AVMU
     SIGMA = TMP0*TMP0 - TMP1
     if ( ABS (SIGMA) < 1.e-6 ) SIGMA = SIGN(1.e-6_r8,SIGMA)
     P1 = B + AVMU*H
     P2 = B - AVMU*H
     P3 = B + TMP0
     P4 = B - TMP0
     S1 = EXP(-H*VAI)
     S2 = EXP(-EXT*VAI)
     IF (IC .EQ. 0) THEN
        U1 = B - C/ALBGRD(IB)
        U2 = B - C*ALBGRD(IB)
        U3 = F + C*ALBGRD(IB)
     ELSE
        U1 = B - C/ALBGRI(IB)
        U2 = B - C*ALBGRI(IB)
        U3 = F + C*ALBGRI(IB)
     END IF
     TMP2 = U1 - AVMU*H
     TMP3 = U1 + AVMU*H
     D1 = P1*TMP2/S1 - P2*TMP3*S1
     TMP4 = U2 + AVMU*H
     TMP5 = U2 - AVMU*H
     D2 = TMP4/S1 - TMP5*S1
     H1 = -D*P4 - C*F
     TMP6 = D - H1*P3/SIGMA
     TMP7 = ( D - C - H1/SIGMA*(U1+TMP0) ) * S2
     H2 = ( TMP6*TMP2/S1 - P2*TMP7 ) / D1
     H3 = - ( TMP6*TMP3*S1 - P1*TMP7 ) / D1
     H4 = -F*P3 - C*D
     TMP8 = H4/SIGMA
     TMP9 = ( U3 - TMP8*(U2-TMP0) ) * S2
     H5 = - ( TMP8*TMP4/S1 + TMP9 ) / D2
     H6 = ( TMP8*TMP5*S1 + TMP9 ) / D2
     H7 = (C*TMP2) / (D1*S1)
     H8 = (-C*TMP3*S1) / D1
     H9 = TMP4 / (D2*S1)
     H10 = (-TMP5*S1) / D2

! downward direct and diffuse fluxes below vegetation
! Niu and Yang (2004), JGR.

     IF (IC .EQ. 0) THEN
        FTDS = S2                           *(1.0-GAP) + GAP
        FTIS = (H4*S2/SIGMA + H5*S1 + H6/S1)*(1.0-GAP)
     ELSE
        FTDS = 0.
        FTIS = (H9*S1 + H10/S1)*(1.0-KOPEN) + KOPEN
     END IF
     FTD(IB) = FTDS
     FTI(IB) = FTIS

! flux reflected by the surface (veg. and ground)

     IF (IC .EQ. 0) THEN
        FRES   = (H1/SIGMA + H2 + H3)*(1.0-GAP  ) + ALBGRD(IB)*GAP        
        FREVEG = (H1/SIGMA + H2 + H3)*(1.0-GAP  ) 
        FREBAR = ALBGRD(IB)*GAP                   !jref - separate veg. and ground reflection
     ELSE
        FRES   = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN        
        FREVEG = (H7 + H8) *(1.0-KOPEN) + ALBGRI(IB)*KOPEN
        FREBAR = 0                                !jref - separate veg. and ground reflection
     END IF
     FRE(IB) = FRES

     FREV(IB) = FREVEG 
     FREG(IB) = FREBAR 

! flux absorbed by vegetation

     FAB(IB) = 1. - FRE(IB) - (1.-ALBGRD(IB))*FTD(IB) &
                            - (1.-ALBGRI(IB))*FTI(IB)

!if(iloc == 1.and.jloc ==  2) then
!  write(*,'(a7,2i2,5(a6,f8.4),2(a9,f8.4))') "ib,ic: ",ib,ic," GAP: ",GAP," FTD: ",FTD(IB)," FTI: ",FTI(IB)," FRE: ", &
!         FRE(IB)," FAB: ",FAB(IB)," ALBGRD: ",ALBGRD(IB)," ALBGRI: ",ALBGRI(IB)
!end if

  END SUBROUTINE TWOSTREAM

!== begin vege_flux ================================================================================

  SUBROUTINE VEGE_FLUX(parameters,NSNOW   ,NSOIL   ,ISNOW   ,VEGTYP  ,VEG     , & !in
                       DT      ,SAV     ,SAG     ,LWDN    ,UR      , & !in
                       UU      ,VV      ,SFCTMP  ,THAIR   ,QAIR    , & !in
                       EAIR    ,RHOAIR  ,SNOWH   ,VAI     ,GAMMAV   ,GAMMAG,  & !in
                       FWET    ,LAISUN  ,LAISHA  ,CWP     ,DZSNSO  , & !in
                       ZLVL    ,ZPD     ,Z0M     ,FVEG    , & !in
                       Z0MG    ,EMV     ,EMG     ,CANLIQ  ,FSNO,          & !in
                       CANICE  ,STC     ,DF      ,RSSUN   ,RSSHA   , & !in
                       RSURF   ,LATHEAV ,LATHEAG  ,PARSUN  ,PARSHA  ,IGS     , & !in
                       FOLN    ,CO2AIR  ,O2AIR   ,BTRAN   ,SFCPRS  , & !in
                       RHSUR   ,ILOC    ,JLOC    ,Q2      ,PAHV    ,PAHG     , & !in
                       EAH     ,TAH     ,TV      ,TG      ,CM      , & !inout
                       CH      ,DX      ,DZ8W    ,                   & !
                       TAUXV   ,TAUYV   ,IRG     ,IRC     ,SHG     , & !out
                       SHC     ,EVG     ,EVC     ,TR      ,GH      , & !out
                       T2MV    ,PSNSUN  ,PSNSHA  ,                   & !out
                       QC      ,QSFC    ,PSFC    ,                   & !in
                       Q2V     ,CAH2    ,CHLEAF  ,CHUC    )            !inout 

! --------------------------------------------------------------------------------------------------
! use newton-raphson iteration to solve for vegetation (tv) and
! ground (tg) temperatures that balance the surface energy budgets

! vegetated:
! -SAV + IRC[TV] + SHC[TV] + EVC[TV] + TR[TV] = 0
! -SAG + IRG[TG] + SHG[TG] + EVG[TG] + GH[TG] = 0
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                         INTENT(IN) :: ILOC   !grid index
  INTEGER,                         INTENT(IN) :: JLOC   !grid index
  LOGICAL,                         INTENT(IN) :: VEG    !true if vegetated surface
  INTEGER,                         INTENT(IN) :: NSNOW  !maximum no. of snow layers        
  INTEGER,                         INTENT(IN) :: NSOIL  !number of soil layers
  INTEGER,                         INTENT(IN) :: ISNOW  !actual no. of snow layers
  INTEGER,                         INTENT(IN) :: VEGTYP !vegetation physiology type
  REAL(r8),                            INTENT(IN) :: FVEG   !greeness vegetation fraction (-)
  REAL(r8),                            INTENT(IN) :: SAV    !solar rad absorbed by veg (w/m2)
  REAL(r8),                            INTENT(IN) :: SAG    !solar rad absorbed by ground (w/m2)
  REAL(r8),                            INTENT(IN) :: LWDN   !atmospheric longwave radiation (w/m2)
  REAL(r8),                            INTENT(IN) :: UR     !wind speed at height zlvl (m/s)
  REAL(r8),                            INTENT(IN) :: UU     !wind speed in eastward dir (m/s)
  REAL(r8),                            INTENT(IN) :: VV     !wind speed in northward dir (m/s)
  REAL(r8),                            INTENT(IN) :: SFCTMP !air temperature at reference height (k)
  REAL(r8),                            INTENT(IN) :: THAIR  !potential temp at reference height (k)
  REAL(r8),                            INTENT(IN) :: EAIR   !vapor pressure air at zlvl (pa)
  REAL(r8),                            INTENT(IN) :: QAIR   !specific humidity at zlvl (kg/kg)
  REAL(r8),                            INTENT(IN) :: RHOAIR !density air (kg/m**3)
  REAL(r8),                            INTENT(IN) :: DT     !time step (s)
  REAL(r8),                            INTENT(IN) :: FSNO     !snow fraction

  REAL(r8),                            INTENT(IN) :: SNOWH  !actual snow depth [m]
  REAL(r8),                            INTENT(IN) :: FWET   !wetted fraction of canopy
  REAL(r8),                            INTENT(IN) :: CWP    !canopy wind parameter

  REAL(r8),                            INTENT(IN) :: VAI    !total leaf area index + stem area index
  REAL(r8),                            INTENT(IN) :: LAISUN !sunlit leaf area index, one-sided (m2/m2)
  REAL(r8),                            INTENT(IN) :: LAISHA !shaded leaf area index, one-sided (m2/m2)
  REAL(r8),                            INTENT(IN) :: ZLVL   !reference height (m)
  REAL(r8),                            INTENT(IN) :: ZPD    !zero plane displacement (m)
  REAL(r8),                            INTENT(IN) :: Z0M    !roughness length, momentum (m)
  REAL(r8),                            INTENT(IN) :: Z0MG   !roughness length, momentum, ground (m)
  REAL(r8),                            INTENT(IN) :: EMV    !vegetation emissivity
  REAL(r8),                            INTENT(IN) :: EMG    !ground emissivity

  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC    !soil/snow temperature (k)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DF     !thermal conductivity of snow/soil (w/m/k)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !thinkness of snow/soil layers (m)
  REAL(r8),                            INTENT(IN) :: CANLIQ !intercepted liquid water (mm)
  REAL(r8),                            INTENT(IN) :: CANICE !intercepted ice mass (mm)
  REAL(r8),                            INTENT(IN) :: RSURF  !ground surface resistance (s/m)
!  REAL(r8),                            INTENT(IN) :: GAMMA  !psychrometric constant (pa/K)
!  REAL(r8),                            INTENT(IN) :: LATHEA !latent heat of vaporization/subli (j/kg)
  REAL(r8),                            INTENT(IN) :: GAMMAV  !psychrometric constant (pa/K)
  REAL(r8),                            INTENT(IN) :: LATHEAV !latent heat of vaporization/subli (j/kg)
  REAL(r8),                            INTENT(IN) :: GAMMAG  !psychrometric constant (pa/K)
  REAL(r8),                            INTENT(IN) :: LATHEAG !latent heat of vaporization/subli (j/kg)
  REAL(r8),                            INTENT(IN) :: PARSUN !par absorbed per unit sunlit lai (w/m2)
  REAL(r8),                            INTENT(IN) :: PARSHA !par absorbed per unit shaded lai (w/m2)
  REAL(r8),                            INTENT(IN) :: FOLN   !foliage nitrogen (%)
  REAL(r8),                            INTENT(IN) :: CO2AIR !atmospheric co2 concentration (pa)
  REAL(r8),                            INTENT(IN) :: O2AIR  !atmospheric o2 concentration (pa)
  REAL(r8),                            INTENT(IN) :: IGS    !growing season index (0=off, 1=on)
  REAL(r8),                            INTENT(IN) :: SFCPRS !pressure (pa)
  REAL(r8),                            INTENT(IN) :: BTRAN  !soil water transpiration factor (0 to 1)
  REAL(r8),                            INTENT(IN) :: RHSUR  !raltive humidity in surface soil/snow air space (-)

  REAL(r8)                           , INTENT(IN) :: QC     !cloud water mixing ratio
  REAL(r8)                           , INTENT(IN) :: PSFC   !pressure at lowest model layer
  REAL(r8)                           , INTENT(IN) :: DX     !grid spacing
  REAL(r8)                           , INTENT(IN) :: Q2     !mixing ratio (kg/kg)
  REAL(r8)                           , INTENT(IN) :: DZ8W   !thickness of lowest layer
  REAL(r8)                           , INTENT(INOUT) :: QSFC   !mixing ratio at lowest model layer
  REAL(r8), INTENT(IN)   :: PAHV  !precipitation advected heat - canopy net IN (W/m2)
  REAL(r8), INTENT(IN)   :: PAHG  !precipitation advected heat - ground net IN (W/m2)

! input/output
  REAL(r8),                         INTENT(INOUT) :: EAH    !canopy air vapor pressure (pa)
  REAL(r8),                         INTENT(INOUT) :: TAH    !canopy air temperature (k)
  REAL(r8),                         INTENT(INOUT) :: TV     !vegetation temperature (k)
  REAL(r8),                         INTENT(INOUT) :: TG     !ground temperature (k)
  REAL(r8),                         INTENT(INOUT) :: CM     !momentum drag coefficient
  REAL(r8),                         INTENT(INOUT) :: CH     !sensible heat exchange coefficient

! output
! -FSA + FIRA + FSH + (FCEV + FCTR + FGEV) + FCST + SSOIL = 0
  REAL(r8),                           INTENT(OUT) :: TAUXV  !wind stress: e-w (n/m2)
  REAL(r8),                           INTENT(OUT) :: TAUYV  !wind stress: n-s (n/m2)
  REAL(r8),                           INTENT(OUT) :: IRC    !net longwave radiation (w/m2) [+= to atm]
  REAL(r8),                           INTENT(OUT) :: SHC    !sensible heat flux (w/m2)     [+= to atm]
  REAL(r8),                           INTENT(OUT) :: EVC    !evaporation heat flux (w/m2)  [+= to atm]
  REAL(r8),                           INTENT(OUT) :: IRG    !net longwave radiation (w/m2) [+= to atm]
  REAL(r8),                           INTENT(OUT) :: SHG    !sensible heat flux (w/m2)     [+= to atm]
  REAL(r8),                           INTENT(OUT) :: EVG    !evaporation heat flux (w/m2)  [+= to atm]
  REAL(r8),                           INTENT(OUT) :: TR     !transpiration heat flux (w/m2)[+= to atm]
  REAL(r8),                           INTENT(OUT) :: GH     !ground heat (w/m2) [+ = to soil]
  REAL(r8),                           INTENT(OUT) :: T2MV   !2 m height air temperature (k)
  REAL(r8),                           INTENT(OUT) :: PSNSUN !sunlit leaf photosynthesis (umolco2/m2/s)
  REAL(r8),                           INTENT(OUT) :: PSNSHA !shaded leaf photosynthesis (umolco2/m2/s)
  REAL(r8),                           INTENT(OUT) :: CHLEAF !leaf exchange coefficient
  REAL(r8),                           INTENT(OUT) :: CHUC   !under canopy exchange coefficient

  REAL(r8),                           INTENT(OUT) :: Q2V
  REAL(r8) :: CAH    !sensible heat conductance, canopy air to ZLVL air (m/s)
  REAL(r8) :: U10V    !10 m wind speed in eastward dir (m/s) 
  REAL(r8) :: V10V    !10 m wind speed in eastward dir (m/s) 
  REAL(r8) :: WSPD

! ------------------------ local variables ----------------------------------------------------
  REAL(r8) :: CW           !water vapor exchange coefficient
  REAL(r8) :: FV           !friction velocity (m/s)
  REAL(r8) :: WSTAR        !friction velocity n vertical direction (m/s) (only for SFCDIF2)
  REAL(r8) :: Z0H          !roughness length, sensible heat (m)
  REAL(r8) :: Z0HG         !roughness length, sensible heat (m)
  REAL(r8) :: RB           !bulk leaf boundary layer resistance (s/m)
  REAL(r8) :: RAMC         !aerodynamic resistance for momentum (s/m)
  REAL(r8) :: RAHC         !aerodynamic resistance for sensible heat (s/m)
  REAL(r8) :: RAWC         !aerodynamic resistance for water vapor (s/m)
  REAL(r8) :: RAMG         !aerodynamic resistance for momentum (s/m)
  REAL(r8) :: RAHG         !aerodynamic resistance for sensible heat (s/m)
  REAL(r8) :: RAWG         !aerodynamic resistance for water vapor (s/m)

  REAL(r8), INTENT(OUT) :: RSSUN        !sunlit leaf stomatal resistance (s/m)
  REAL(r8), INTENT(OUT) :: RSSHA        !shaded leaf stomatal resistance (s/m)

  REAL(r8) :: MOL          !Monin-Obukhov length (m)
  REAL(r8) :: DTV          !change in tv, last iteration (k)
  REAL(r8) :: DTG          !change in tg, last iteration (k)

  REAL(r8) :: AIR,CIR      !coefficients for ir as function of ts**4
  REAL(r8) :: CSH          !coefficients for sh as function of ts
  REAL(r8) :: CEV          !coefficients for ev as function of esat[ts]
  REAL(r8) :: CGH          !coefficients for st as function of ts
  REAL(r8) :: ATR,CTR      !coefficients for tr as function of esat[ts]
  REAL(r8) :: ATA,BTA      !coefficients for tah as function of ts
  REAL(r8) :: AEA,BEA      !coefficients for eah as function of esat[ts]

  REAL(r8) :: ESTV         !saturation vapor pressure at tv (pa)
  REAL(r8) :: ESTG         !saturation vapor pressure at tg (pa)
  REAL(r8) :: DESTV        !d(es)/dt at ts (pa/k)
  REAL(r8) :: DESTG        !d(es)/dt at tg (pa/k)
  REAL(r8) :: ESATW        !es for water
  REAL(r8) :: ESATI        !es for ice
  REAL(r8) :: DSATW        !d(es)/dt at tg (pa/k) for water
  REAL(r8) :: DSATI        !d(es)/dt at tg (pa/k) for ice

  REAL(r8) :: FM           !momentum stability correction, weighted by prior iters
  REAL(r8) :: FH           !sen heat stability correction, weighted by prior iters
  REAL(r8) :: FHG          !sen heat stability correction, ground
  REAL(r8) :: HCAN         !canopy height (m) [note: hcan >= z0mg]

  REAL(r8) :: A            !temporary calculation
  REAL(r8) :: B            !temporary calculation
  REAL(r8) :: CVH          !sensible heat conductance, leaf surface to canopy air (m/s)
  REAL(r8) :: CAW          !latent heat conductance, canopy air ZLVL air (m/s)
  REAL(r8) :: CTW          !transpiration conductance, leaf to canopy air (m/s)
  REAL(r8) :: CEW          !evaporation conductance, leaf to canopy air (m/s)
  REAL(r8) :: CGW          !latent heat conductance, ground to canopy air (m/s)
  REAL(r8) :: COND         !sum of conductances (s/m)
  REAL(r8) :: UC           !wind speed at top of canopy (m/s)
  REAL(r8) :: KH           !turbulent transfer coefficient, sensible heat, (m2/s)
  REAL(r8) :: H            !temporary sensible heat flux (w/m2)
  REAL(r8) :: HG           !temporary sensible heat flux (w/m2)
  REAL(r8) :: MOZ          !Monin-Obukhov stability parameter
  REAL(r8) :: MOZG         !Monin-Obukhov stability parameter
  REAL(r8) :: MOZOLD       !Monin-Obukhov stability parameter from prior iteration
  REAL(r8) :: FM2          !Monin-Obukhov momentum adjustment at 2m
  REAL(r8) :: FH2          !Monin-Obukhov heat adjustment at 2m
  REAL(r8) :: CH2          !Surface exchange at 2m
  REAL(r8) :: THSTAR          !Surface exchange at 2m

  REAL(r8) :: THVAIR
  REAL(r8) :: THAH 
  REAL(r8) :: RAHC2        !aerodynamic resistance for sensible heat (s/m)
  REAL(r8) :: RAWC2        !aerodynamic resistance for water vapor (s/m)
  REAL(r8), INTENT(OUT):: CAH2         !sensible heat conductance for diagnostics
  REAL(r8) :: CH2V         !exchange coefficient for 2m over vegetation. 
  REAL(r8) :: CQ2V         !exchange coefficient for 2m over vegetation. 
  REAL(r8) :: EAH2         !2m vapor pressure over canopy
  REAL(r8) :: QFX        !moisture flux
  REAL(r8) :: E1           


  REAL(r8) :: VAIE         !total leaf area index + stem area index,effective
  REAL(r8) :: LAISUNE      !sunlit leaf area index, one-sided (m2/m2),effective
  REAL(r8) :: LAISHAE      !shaded leaf area index, one-sided (m2/m2),effective

  INTEGER :: K         !index
  INTEGER :: ITER      !iteration index

!jref - NITERC test from 5 to 20  
  INTEGER, PARAMETER :: NITERC = 20   !number of iterations for surface temperature
!jref - NITERG test from 3-5
  INTEGER, PARAMETER :: NITERG = 5   !number of iterations for ground temperature
  INTEGER :: MOZSGN    !number of times MOZ changes sign
  REAL(r8)    :: MPE       !prevents overflow error if division by zero

  INTEGER :: LITER     !Last iteration


  REAL(r8) :: T, TDC       !Kelvin to degree Celsius with limit -50 to +50

  character(len=80) ::  message

  TDC(T)   = MIN( 50., MAX(-50.,(T-TFRZ)) )
! ---------------------------------------------------------------------------------------------

        MPE = 1E-6
        LITER = 0
        FV = 0.1

! ---------------------------------------------------------------------------------------------
! initialization variables that do not depend on stability iteration
! ---------------------------------------------------------------------------------------------
        DTV = 0.
        DTG = 0.
        MOZ    = 0.
        MOZSGN = 0
        MOZOLD = 0.
        HG     = 0.
        H      = 0.
        QFX    = 0.

! convert grid-cell LAI to the fractional vegetated area (FVEG)

        VAIE    = MIN(6.,VAI    / FVEG)
        LAISUNE = MIN(6.,LAISUN / FVEG)
        LAISHAE = MIN(6.,LAISHA / FVEG)

! saturation vapor pressure at ground temperature

        T = TDC(TG)
        CALL ESAT(T, ESATW, ESATI, DSATW, DSATI)
        IF (T .GT. 0.) THEN
           ESTG = ESATW
        ELSE
           ESTG = ESATI
        END IF

!jref - consistent surface specific humidity for sfcdif3 and sfcdif4
!#ifdef AMIPW_PHYSICS
        QSFC = 0.622*EAIR/(PSFC-0.378*EAIR)  
!#endif
! canopy height

        HCAN = parameters%HVT
        UC = UR*LOG(HCAN/Z0M)/LOG(ZLVL/Z0M)
        UC = UR*LOG((HCAN-ZPD+Z0M)/Z0M)/LOG(ZLVL/Z0M)   ! MB: add ZPD v3.7
        IF((HCAN-ZPD) <= 0.) THEN
          WRITE(message,*) "CRITICAL PROBLEM: HCAN <= ZPD"
          print*, ( message )
          WRITE(message,*) 'i,j point=',ILOC, JLOC
          print*, ( message )
          WRITE(message,*) 'HCAN  =',HCAN
          print*, ( message )
          WRITE(message,*) 'ZPD   =',ZPD
          print*, ( message )
          write (message, *) 'SNOWH =',SNOWH
          print*, ( message )
          print*, ( "CRITICAL PROBLEM IN MODULE_SF_NOAHMPLSM:VEGEFLUX" )
        END IF

! prepare for longwave rad.

        AIR = -EMV*(1.+(1.-EMV)*(1.-EMG))*LWDN - EMV*EMG*SB*TG**4  
        CIR = (2.-EMV*(1.-EMG))*EMV*SB
! ---------------------------------------------------------------------------------------------
      loop1: DO ITER = 1, NITERC    !  begin stability iteration

       IF(ITER == 1) THEN
            Z0H  = Z0M  
            Z0HG = Z0MG
       ELSE
            Z0H  = Z0M    !* EXP(-CZIL*0.4*258.2*SQRT(FV*Z0M))
            Z0HG = Z0MG   !* EXP(-CZIL*0.4*258.2*SQRT(FV*Z0MG))
       END IF

! aerodyn resistances between heights zlvl and d+z0v

       IF(OPT_SFC == 1) THEN
          CALL SFCDIF1(parameters,ITER   ,SFCTMP ,RHOAIR ,H      ,QAIR   , & !in
                       ZLVL   ,ZPD    ,Z0M    ,Z0H    ,UR     , & !in
                       MPE    ,ILOC   ,JLOC   ,                 & !in
                       MOZ    ,MOZSGN ,FM     ,FH     ,FM2,FH2, & !inout
                       CM     ,CH     ,FV     ,CH2     )          !out
       ENDIF
     
       IF(OPT_SFC == 2) THEN
          CALL SFCDIF2(parameters,ITER   ,Z0M    ,TAH    ,THAIR  ,UR     , & !in
                       ZLVL   ,ILOC   ,JLOC   ,         & !in
                       CM     ,CH     ,MOZ    ,WSTAR  ,         & !in
                       FV     )                                   !out
          ! Undo the multiplication by windspeed that SFCDIF2 
          ! applies to exchange coefficients CH and CM:
          CH = CH / UR
          CM = CM / UR
       ENDIF

       RAMC = MAX(1.,1./(CM*UR))
       RAHC = MAX(1.,1./(CH*UR))
       RAWC = RAHC

! aerodyn resistance between heights z0g and d+z0v, RAG, and leaf
! boundary layer resistance, RB
       
       CALL RAGRB(parameters,ITER   ,VAIE   ,RHOAIR ,HG     ,TAH    , & !in
                  ZPD    ,Z0MG   ,Z0HG   ,HCAN   ,UC     , & !in
                  Z0H    ,FV     ,CWP    ,VEGTYP ,MPE    , & !in
                  TV     ,MOZG   ,FHG    ,ILOC   ,JLOC   , & !inout
                  RAMG   ,RAHG   ,RAWG   ,RB     )           !out

! es and d(es)/dt evaluated at tv

       T = TDC(TV)
       CALL ESAT(T, ESATW, ESATI, DSATW, DSATI)
       IF (T .GT. 0.) THEN
          ESTV  = ESATW
          DESTV = DSATW
       ELSE
          ESTV  = ESATI
          DESTV = DSATI
       END IF

! stomatal resistance
        
     IF(ITER == 1) THEN
        IF (OPT_CRS == 1) then  ! Ball-Berry
         CALL STOMATA (parameters,VEGTYP,MPE   ,PARSUN ,FOLN  ,ILOC  , JLOC , & !in       
                       TV    ,ESTV  ,EAH    ,SFCTMP,SFCPRS, & !in
                       O2AIR ,CO2AIR,IGS    ,BTRAN ,RB    , & !in
                       RSSUN ,PSNSUN)                         !out

         CALL STOMATA (parameters,VEGTYP,MPE   ,PARSHA ,FOLN  ,ILOC  , JLOC , & !in
                       TV    ,ESTV  ,EAH    ,SFCTMP,SFCPRS, & !in
                       O2AIR ,CO2AIR,IGS    ,BTRAN ,RB    , & !in
                       RSSHA ,PSNSHA)                         !out
        END IF

        IF (OPT_CRS == 2) then  ! Jarvis
         CALL  CANRES (parameters,PARSUN,TV    ,BTRAN ,EAH    ,SFCPRS, & !in
                       RSSUN ,PSNSUN,ILOC  ,JLOC   )          !out

         CALL  CANRES (parameters,PARSHA,TV    ,BTRAN ,EAH    ,SFCPRS, & !in
                       RSSHA ,PSNSHA,ILOC  ,JLOC   )          !out
        END IF
     END IF

! prepare for sensible heat flux above veg.

        CAH  = 1./RAHC
        CVH  = 2.*VAIE/RB
        CGH  = 1./RAHG
        COND = CAH + CVH + CGH
        ATA  = (SFCTMP*CAH + TG*CGH) / COND
        BTA  = CVH/COND
        CSH  = (1.-BTA)*RHOAIR*CPAIR*CVH

! prepare for latent heat flux above veg.

        CAW  = 1./RAWC
        CEW  = FWET*VAIE/RB
        CTW  = (1.-FWET)*(LAISUNE/(RB+RSSUN) + LAISHAE/(RB+RSSHA))
        CGW  = 1./(RAWG+RSURF)
        COND = CAW + CEW + CTW + CGW
        AEA  = (EAIR*CAW + ESTG*CGW) / COND
        BEA  = (CEW+CTW)/COND
        CEV  = (1.-BEA)*CEW*RHOAIR*CPAIR/GAMMAV   ! Barlage: change to vegetation v3.6
        CTR  = (1.-BEA)*CTW*RHOAIR*CPAIR/GAMMAV

! evaluate surface fluxes with current temperature and solve for dts

        TAH = ATA + BTA*TV               ! canopy air T.
        EAH = AEA + BEA*ESTV             ! canopy air e

        IRC = FVEG*(AIR + CIR*TV**4)
        SHC = FVEG*RHOAIR*CPAIR*CVH * (  TV-TAH)
        EVC = FVEG*RHOAIR*CPAIR*CEW * (ESTV-EAH) / GAMMAV ! Barlage: change to v in v3.6
        TR  = FVEG*RHOAIR*CPAIR*CTW * (ESTV-EAH) / GAMMAV
	IF (TV > TFRZ) THEN
          EVC = MIN(CANLIQ*LATHEAV/DT,EVC)    ! Barlage: add if block for canice in v3.6
	ELSE
          EVC = MIN(CANICE*LATHEAV/DT,EVC)
	END IF

        B   = SAV-IRC-SHC-EVC-TR+PAHV                          !additional w/m2
        A   = FVEG*(4.*CIR*TV**3 + CSH + (CEV+CTR)*DESTV) !volumetric heat capacity
        DTV = B/A

        IRC = IRC + FVEG*4.*CIR*TV**3*DTV
        SHC = SHC + FVEG*CSH*DTV
        EVC = EVC + FVEG*CEV*DESTV*DTV
        TR  = TR  + FVEG*CTR*DESTV*DTV                               

! update vegetation surface temperature
        TV  = TV + DTV
!        TAH = ATA + BTA*TV               ! canopy air T; update here for consistency

! for computing M-O length in the next iteration
        H  = RHOAIR*CPAIR*(TAH - SFCTMP) /RAHC        
        HG = RHOAIR*CPAIR*(TG  - TAH)   /RAHG

! consistent specific humidity from canopy air vapor pressure
        QSFC = (0.622*EAH)/(SFCPRS-0.378*EAH)

        IF (LITER == 1) THEN
           exit loop1 
        ENDIF
        IF (ITER >= 5 .AND. ABS(DTV) <= 0.01 .AND. LITER == 0) THEN
           LITER = 1
        ENDIF

     END DO loop1 ! end stability iteration

! under-canopy fluxes and tg

        AIR = - EMG*(1.-EMV)*LWDN - EMG*EMV*SB*TV**4
        CIR = EMG*SB
        CSH = RHOAIR*CPAIR/RAHG
        CEV = RHOAIR*CPAIR / (GAMMAG*(RAWG+RSURF))  ! Barlage: change to ground v3.6
        CGH = 2.*DF(ISNOW+1)/DZSNSO(ISNOW+1)

     loop2: DO ITER = 1, NITERG

        T = TDC(TG)
        CALL ESAT(T, ESATW, ESATI, DSATW, DSATI)
        IF (T .GT. 0.) THEN
            ESTG  = ESATW
            DESTG = DSATW
        ELSE
            ESTG  = ESATI
            DESTG = DSATI
        END IF

        IRG = CIR*TG**4 + AIR
        SHG = CSH * (TG         - TAH         )
        EVG = CEV * (ESTG*RHSUR - EAH         )
        GH  = CGH * (TG         - STC(ISNOW+1))

        B = SAG-IRG-SHG-EVG-GH+PAHG
        A = 4.*CIR*TG**3+CSH+CEV*DESTG+CGH
        DTG = B/A

        IRG = IRG + 4.*CIR*TG**3*DTG
        SHG = SHG + CSH*DTG
        EVG = EVG + CEV*DESTG*DTG
        GH  = GH  + CGH*DTG
        TG  = TG  + DTG

     END DO loop2
     
!     TAH = (CAH*SFCTMP + CVH*TV + CGH*TG)/(CAH + CVH + CGH)

! if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.

     IF(OPT_STC == 1 .OR. OPT_STC == 3) THEN
     IF (SNOWH > 0.05 .AND. TG > TFRZ) THEN
        IF(OPT_STC == 1) TG  = TFRZ
        IF(OPT_STC == 3) TG  = (1.-FSNO)*TG + FSNO*TFRZ   ! MB: allow TG>0C during melt v3.7
        IRG = CIR*TG**4 - EMG*(1.-EMV)*LWDN - EMG*EMV*SB*TV**4
        SHG = CSH * (TG         - TAH)
        EVG = CEV * (ESTG*RHSUR - EAH)
        GH  = SAG+PAHG - (IRG+SHG+EVG)
     END IF
     END IF

! wind stresses

     TAUXV = -RHOAIR*CM*UR*UU
     TAUYV = -RHOAIR*CM*UR*VV

! consistent vegetation air temperature and vapor pressure since TG is not consistent with the TAH/EAH
! calculation.
!     TAH = SFCTMP + (SHG+SHC)/(RHOAIR*CPAIR*CAH) 
!     TAH = SFCTMP + (SHG*FVEG+SHC)/(RHOAIR*CPAIR*CAH) ! ground flux need fveg
!     EAH = EAIR + (EVC+FVEG*(TR+EVG))/(RHOAIR*CAW*CPAIR/GAMMAG )
!     QFX = (QSFC-QAIR)*RHOAIR*CAW !*CPAIR/GAMMAG

! 2m temperature over vegetation ( corrected for low CQ2V values )
   IF (OPT_SFC == 1 .OR. OPT_SFC == 2) THEN
!      CAH2 = FV*1./VKC*LOG((2.+Z0H)/Z0H)
      CAH2 = FV*VKC/LOG((2.+Z0H)/Z0H)
      CAH2 = FV*VKC/(LOG((2.+Z0H)/Z0H)-FH2)
      CQ2V = CAH2
      IF (CAH2 .LT. 1.E-5 ) THEN
         T2MV = TAH
!         Q2V  = (EAH*0.622/(SFCPRS - 0.378*EAH))
         Q2V  = QSFC
      ELSE
         T2MV = TAH - (SHG+SHC/FVEG)/(RHOAIR*CPAIR) * 1./CAH2
!         Q2V = (EAH*0.622/(SFCPRS - 0.378*EAH))- QFX/(RHOAIR*FV)* 1./VKC * LOG((2.+Z0H)/Z0H)
         Q2V = QSFC - ((EVC+TR)/FVEG+EVG)/(LATHEAV*RHOAIR) * 1./CQ2V
      ENDIF
   ENDIF

! update CH for output
     CH = CAH
     CHLEAF = CVH
     CHUC = 1./RAHG

  END SUBROUTINE VEGE_FLUX

!== begin bare_flux ================================================================================

  SUBROUTINE BARE_FLUX (parameters,NSNOW   ,NSOIL   ,ISNOW   ,DT      ,SAG     , & !in
                        LWDN    ,UR      ,UU      ,VV      ,SFCTMP  , & !in
                        THAIR   ,QAIR    ,EAIR    ,RHOAIR  ,SNOWH   , & !in
                        DZSNSO  ,ZLVL    ,ZPD     ,Z0M     ,FSNO    , & !in
                        EMG     ,STC     ,DF      ,RSURF   ,LATHEA  , & !in
                        GAMMA   ,RHSUR   ,ILOC    ,JLOC    ,Q2      ,PAHB  , & !in
                        TGB     ,CM      ,CH      ,          & !inout
                        TAUXB   ,TAUYB   ,IRB     ,SHB     ,EVB     , & !out
                        GHB     ,T2MB    ,DX      ,DZ8W    ,IVGTYP  , & !out
                        QC      ,QSFC    ,PSFC    ,                   & !in
                        SFCPRS  ,Q2B     ,EHB2    )                     !in 

! --------------------------------------------------------------------------------------------------
! use newton-raphson iteration to solve ground (tg) temperature
! that balances the surface energy budgets for bare soil fraction.

! bare soil:
! -SAB + IRB[TG] + SHB[TG] + EVB[TG] + GHB[TG] = 0
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer                        , INTENT(IN) :: ILOC   !grid index
  integer                        , INTENT(IN) :: JLOC   !grid index
  INTEGER,                         INTENT(IN) :: NSNOW  !maximum no. of snow layers
  INTEGER,                         INTENT(IN) :: NSOIL  !number of soil layers
  INTEGER,                         INTENT(IN) :: ISNOW  !actual no. of snow layers
  REAL(r8),                            INTENT(IN) :: DT     !time step (s)
  REAL(r8),                            INTENT(IN) :: SAG    !solar radiation absorbed by ground (w/m2)
  REAL(r8),                            INTENT(IN) :: LWDN   !atmospheric longwave radiation (w/m2)
  REAL(r8),                            INTENT(IN) :: UR     !wind speed at height zlvl (m/s)
  REAL(r8),                            INTENT(IN) :: UU     !wind speed in eastward dir (m/s)
  REAL(r8),                            INTENT(IN) :: VV     !wind speed in northward dir (m/s)
  REAL(r8),                            INTENT(IN) :: SFCTMP !air temperature at reference height (k)
  REAL(r8),                            INTENT(IN) :: THAIR  !potential temperature at height zlvl (k)
  REAL(r8),                            INTENT(IN) :: QAIR   !specific humidity at height zlvl (kg/kg)
  REAL(r8),                            INTENT(IN) :: EAIR   !vapor pressure air at height (pa)
  REAL(r8),                            INTENT(IN) :: RHOAIR !density air (kg/m3)
  REAL(r8),                            INTENT(IN) :: SNOWH  !actual snow depth [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !thickness of snow/soil layers (m)
  REAL(r8),                            INTENT(IN) :: ZLVL   !reference height (m)
  REAL(r8),                            INTENT(IN) :: ZPD    !zero plane displacement (m)
  REAL(r8),                            INTENT(IN) :: Z0M    !roughness length, momentum, ground (m)
  REAL(r8),                            INTENT(IN) :: EMG    !ground emissivity
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC    !soil/snow temperature (k)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DF     !thermal conductivity of snow/soil (w/m/k)
  REAL(r8),                            INTENT(IN) :: RSURF  !ground surface resistance (s/m)
  REAL(r8),                            INTENT(IN) :: LATHEA !latent heat of vaporization/subli (j/kg)
  REAL(r8),                            INTENT(IN) :: GAMMA  !psychrometric constant (pa/k)
  REAL(r8),                            INTENT(IN) :: RHSUR  !raltive humidity in surface soil/snow air space (-)
  REAL(r8),                            INTENT(IN) :: FSNO     !snow fraction

!jref:start; in 
  INTEGER                        , INTENT(IN) :: IVGTYP
  REAL(r8)                           , INTENT(IN) :: QC     !cloud water mixing ratio
  REAL(r8)                           , INTENT(INOUT) :: QSFC   !mixing ratio at lowest model layer
  REAL(r8)                           , INTENT(IN) :: PSFC   !pressure at lowest model layer
  REAL(r8)                           , INTENT(IN) :: SFCPRS !pressure at lowest model layer
  REAL(r8)                           , INTENT(IN) :: DX     !horisontal grid spacing
  REAL(r8)                           , INTENT(IN) :: Q2     !mixing ratio (kg/kg)
  REAL(r8)                           , INTENT(IN) :: DZ8W   !thickness of lowest layer
!jref:end
  REAL(r8), INTENT(IN)   :: PAHB  !precipitation advected heat - ground net IN (W/m2)

! input/output
  REAL(r8),                         INTENT(INOUT) :: TGB    !ground temperature (k)
  REAL(r8),                         INTENT(INOUT) :: CM     !momentum drag coefficient
  REAL(r8),                         INTENT(INOUT) :: CH     !sensible heat exchange coefficient

! output
! -SAB + IRB[TG] + SHB[TG] + EVB[TG] + GHB[TG] = 0

  REAL(r8),                           INTENT(OUT) :: TAUXB  !wind stress: e-w (n/m2)
  REAL(r8),                           INTENT(OUT) :: TAUYB  !wind stress: n-s (n/m2)
  REAL(r8),                           INTENT(OUT) :: IRB    !net longwave rad (w/m2)   [+ to atm]
  REAL(r8),                           INTENT(OUT) :: SHB    !sensible heat flux (w/m2) [+ to atm]
  REAL(r8),                           INTENT(OUT) :: EVB    !latent heat flux (w/m2)   [+ to atm]
  REAL(r8),                           INTENT(OUT) :: GHB    !ground heat flux (w/m2)  [+ to soil]
  REAL(r8),                           INTENT(OUT) :: T2MB   !2 m height air temperature (k)
!jref:start
  REAL(r8),                           INTENT(OUT) :: Q2B    !bare ground heat conductance
  REAL(r8) :: EHB    !bare ground heat conductance
  REAL(r8) :: U10B    !10 m wind speed in eastward dir (m/s)
  REAL(r8) :: V10B    !10 m wind speed in eastward dir (m/s)
  REAL(r8) :: WSPD
!jref:end

! local variables 

  REAL(r8) :: TAUX       !wind stress: e-w (n/m2)
  REAL(r8) :: TAUY       !wind stress: n-s (n/m2)
  REAL(r8) :: FIRA       !total net longwave rad (w/m2)      [+ to atm]
  REAL(r8) :: FSH        !total sensible heat flux (w/m2)    [+ to atm]
  REAL(r8) :: FGEV       !ground evaporation heat flux (w/m2)[+ to atm]
  REAL(r8) :: SSOIL      !soil heat flux (w/m2)             [+ to soil]
  REAL(r8) :: FIRE       !emitted ir (w/m2)
  REAL(r8) :: TRAD       !radiative temperature (k)
  REAL(r8) :: TAH        !"surface" temperature at height z0h+zpd (k)

  REAL(r8) :: CW         !water vapor exchange coefficient
  REAL(r8) :: FV         !friction velocity (m/s)
  REAL(r8) :: WSTAR      !friction velocity n vertical direction (m/s) (only for SFCDIF2)
  REAL(r8) :: Z0H        !roughness length, sensible heat, ground (m)
  REAL(r8) :: RB         !bulk leaf boundary layer resistance (s/m)
  REAL(r8) :: RAMB       !aerodynamic resistance for momentum (s/m)
  REAL(r8) :: RAHB       !aerodynamic resistance for sensible heat (s/m)
  REAL(r8) :: RAWB       !aerodynamic resistance for water vapor (s/m)
  REAL(r8) :: MOL        !Monin-Obukhov length (m)
  REAL(r8) :: DTG        !change in tg, last iteration (k)

  REAL(r8) :: CIR        !coefficients for ir as function of ts**4
  REAL(r8) :: CSH        !coefficients for sh as function of ts
  REAL(r8) :: CEV        !coefficients for ev as function of esat[ts]
  REAL(r8) :: CGH        !coefficients for st as function of ts

!jref:start
  REAL(r8) :: RAHB2      !aerodynamic resistance for sensible heat 2m (s/m)
  REAL(r8) :: RAWB2      !aerodynamic resistance for water vapor 2m (s/m)
  REAL(r8),INTENT(OUT) :: EHB2       !sensible heat conductance for diagnostics
  REAL(r8) :: CH2B       !exchange coefficient for 2m temp.
  REAL(r8) :: CQ2B       !exchange coefficient for 2m temp.
  REAL(r8) :: THVAIR     !virtual potential air temp
  REAL(r8) :: THGH       !potential ground temp
  REAL(r8) :: EMB        !momentum conductance
  REAL(r8) :: QFX        !moisture flux
  REAL(r8) :: ESTG2      !saturation vapor pressure at 2m (pa)
  INTEGER :: VEGTYP     !vegetation type set to isbarren
  REAL(r8) :: E1
!jref:end

  REAL(r8) :: ESTG       !saturation vapor pressure at tg (pa)
  REAL(r8) :: DESTG      !d(es)/dt at tg (pa/K)
  REAL(r8) :: ESATW      !es for water
  REAL(r8) :: ESATI      !es for ice
  REAL(r8) :: DSATW      !d(es)/dt at tg (pa/K) for water
  REAL(r8) :: DSATI      !d(es)/dt at tg (pa/K) for ice

  REAL(r8) :: A          !temporary calculation
  REAL(r8) :: B          !temporary calculation
  REAL(r8) :: H          !temporary sensible heat flux (w/m2)
  REAL(r8) :: MOZ        !Monin-Obukhov stability parameter
  REAL(r8) :: MOZOLD     !Monin-Obukhov stability parameter from prior iteration
  REAL(r8) :: FM         !momentum stability correction, weighted by prior iters
  REAL(r8) :: FH         !sen heat stability correction, weighted by prior iters
  INTEGER :: MOZSGN  !number of times MOZ changes sign
  REAL(r8) :: FM2          !Monin-Obukhov momentum adjustment at 2m
  REAL(r8) :: FH2          !Monin-Obukhov heat adjustment at 2m
  REAL(r8) :: CH2          !Surface exchange at 2m

  INTEGER :: ITER    !iteration index
  INTEGER :: NITERB  !number of iterations for surface temperature
  REAL(r8)    :: MPE     !prevents overflow error if division by zero
!jref:start
!  DATA NITERB /3/
  DATA NITERB /5/
  SAVE NITERB
  REAL(r8) :: T, TDC     !Kelvin to degree Celsius with limit -50 to +50
  TDC(T)   = MIN( 50., MAX(-50.,(T-TFRZ)) )

! -----------------------------------------------------------------
! initialization variables that do not depend on stability iteration
! -----------------------------------------------------------------
        MPE = 1E-6
        DTG = 0.
        MOZ    = 0.
        MOZSGN = 0
        MOZOLD = 0.
        H      = 0.
        QFX    = 0.
        FV     = 0.1

        CIR = EMG*SB
        CGH = 2.*DF(ISNOW+1)/DZSNSO(ISNOW+1)
!-> xinyf guess the qsfc initial value from method from CLM
!#ifdef NEW_REGLSM
        qsfc = 0.622*eair/(psfc-0.378*eair)
!#endif
!<-
! -----------------------------------------------------------------
      loop3: DO ITER = 1, NITERB  ! begin stability iteration

        IF(ITER == 1) THEN
            Z0H = Z0M 
        ELSE
            Z0H = Z0M !* EXP(-CZIL*0.4*258.2*SQRT(FV*Z0M))
        END IF

        IF(OPT_SFC == 1) THEN
          CALL SFCDIF1(parameters,ITER   ,SFCTMP ,RHOAIR ,H      ,QAIR   , & !in
                       ZLVL   ,ZPD    ,Z0M    ,Z0H    ,UR     , & !in
                       MPE    ,ILOC   ,JLOC   ,                 & !in
                       MOZ    ,MOZSGN ,FM     ,FH     ,FM2,FH2, & !inout
                       CM     ,CH     ,FV     ,CH2     )          !out
        ENDIF

        IF(OPT_SFC == 2) THEN
          CALL SFCDIF2(parameters,ITER   ,Z0M    ,TGB    ,THAIR  ,UR     , & !in
                       ZLVL   ,ILOC   ,JLOC   ,         & !in
                       CM     ,CH     ,MOZ    ,WSTAR  ,         & !in
                       FV     )                                   !out
          ! Undo the multiplication by windspeed that SFCDIF2 
          ! applies to exchange coefficients CH and CM:
          CH = CH / UR
          CM = CM / UR
          IF(SNOWH > 0.) THEN
             CM = MIN(0.01,CM)   ! CM & CH are too large, causing
             CH = MIN(0.01,CH)   ! computational instability
          END IF

        ENDIF

        RAMB = MAX(1.,1./(CM*UR))
        RAHB = MAX(1.,1./(CH*UR))
        RAWB = RAHB

!jref - variables for diagnostics         
        EMB = 1./RAMB
        EHB = 1./RAHB

! es and d(es)/dt evaluated at tg

        T = TDC(TGB)
        CALL ESAT(T, ESATW, ESATI, DSATW, DSATI)
        IF (T .GT. 0.) THEN
            ESTG  = ESATW
            DESTG = DSATW
        ELSE
            ESTG  = ESATI
            DESTG = DSATI
        END IF

        CSH = RHOAIR*CPAIR/RAHB
        CEV = RHOAIR*CPAIR/GAMMA/(RSURF+RAWB)

! surface fluxes and dtg

        IRB   = CIR * TGB**4 - EMG*LWDN
        SHB   = CSH * (TGB        - SFCTMP      )
        EVB   = CEV * (ESTG*RHSUR - EAIR        )
        GHB   = CGH * (TGB        - STC(ISNOW+1))

        B     = SAG-IRB-SHB-EVB-GHB+PAHB
        A     = 4.*CIR*TGB**3 + CSH + CEV*DESTG + CGH
        DTG   = B/A

        IRB = IRB + 4.*CIR*TGB**3*DTG
        SHB = SHB + CSH*DTG
        EVB = EVB + CEV*DESTG*DTG
        GHB = GHB + CGH*DTG

! update ground surface temperature
        TGB = TGB + DTG

! for M-O length
        H = CSH * (TGB - SFCTMP)

        T = TDC(TGB)
        CALL ESAT(T, ESATW, ESATI, DSATW, DSATI)
        IF (T .GT. 0.) THEN
            ESTG  = ESATW
        ELSE
            ESTG  = ESATI
        END IF
        QSFC = 0.622*(ESTG*RHSUR)/(PSFC-0.378*(ESTG*RHSUR))

        QFX = (QSFC-QAIR)*CEV*GAMMA/CPAIR

     END DO loop3 ! end stability iteration
! -----------------------------------------------------------------

! if snow on ground and TG > TFRZ: reset TG = TFRZ. reevaluate ground fluxes.

     IF(OPT_STC == 1 .OR. OPT_STC == 3) THEN
     IF (SNOWH > 0.05 .AND. TGB > TFRZ) THEN
          IF(OPT_STC == 1) TGB = TFRZ
          IF(OPT_STC == 3) TGB  = (1.-FSNO)*TGB + FSNO*TFRZ  ! MB: allow TG>0C during melt v3.7
          IRB = CIR * TGB**4 - EMG*LWDN
          SHB = CSH * (TGB        - SFCTMP)
          EVB = CEV * (ESTG*RHSUR - EAIR )          !ESTG reevaluate ?
          GHB = SAG+PAHB - (IRB+SHB+EVB)
     END IF
     END IF

! wind stresses
         
     TAUXB = -RHOAIR*CM*UR*UU
     TAUYB = -RHOAIR*CM*UR*VV

!jref:start; errors in original equation corrected.
! 2m air temperature
     IF(OPT_SFC == 1 .OR. OPT_SFC ==2) THEN
       EHB2  = FV*VKC/LOG((2.+Z0H)/Z0H)
       EHB2  = FV*VKC/(LOG((2.+Z0H)/Z0H)-FH2)
       CQ2B  = EHB2
       IF (EHB2.lt.1.E-5 ) THEN
         T2MB  = TGB
         Q2B   = QSFC
       ELSE
         T2MB  = TGB - SHB/(RHOAIR*CPAIR) * 1./EHB2
         Q2B   = QSFC - EVB/(LATHEA*RHOAIR)*(1./CQ2B + RSURF)
       ENDIF
       IF (parameters%urban_flag) Q2B = QSFC
     END IF

! update CH 
     CH = EHB

  END SUBROUTINE BARE_FLUX

!== begin ragrb ====================================================================================

  SUBROUTINE RAGRB(parameters,ITER   ,VAI    ,RHOAIR ,HG     ,TAH    , & !in
                   ZPD    ,Z0MG   ,Z0HG   ,HCAN   ,UC     , & !in
                   Z0H    ,FV     ,CWP    ,VEGTYP ,MPE    , & !in
                   TV     ,MOZG   ,FHG    ,ILOC   ,JLOC   , & !inout
                   RAMG   ,RAHG   ,RAWG   ,RB     )           !out
! --------------------------------------------------------------------------------------------------
! compute under-canopy aerodynamic resistance RAG and leaf boundary layer
! resistance RB
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,              INTENT(IN) :: ILOC   !grid index
  INTEGER,              INTENT(IN) :: JLOC   !grid index
  INTEGER,              INTENT(IN) :: ITER   !iteration index
  INTEGER,              INTENT(IN) :: VEGTYP !vegetation physiology type
  REAL(r8),                 INTENT(IN) :: VAI    !total LAI + stem area index, one sided
  REAL(r8),                 INTENT(IN) :: RHOAIR !density air (kg/m3)
  REAL(r8),                 INTENT(IN) :: HG     !ground sensible heat flux (w/m2)
  REAL(r8),                 INTENT(IN) :: TV     !vegetation temperature (k)
  REAL(r8),                 INTENT(IN) :: TAH    !air temperature at height z0h+zpd (k)
  REAL(r8),                 INTENT(IN) :: ZPD    !zero plane displacement (m)
  REAL(r8),                 INTENT(IN) :: Z0MG   !roughness length, momentum, ground (m)
  REAL(r8),                 INTENT(IN) :: HCAN   !canopy height (m) [note: hcan >= z0mg]
  REAL(r8),                 INTENT(IN) :: UC     !wind speed at top of canopy (m/s)
  REAL(r8),                 INTENT(IN) :: Z0H    !roughness length, sensible heat (m)
  REAL(r8),                 INTENT(IN) :: Z0HG   !roughness length, sensible heat, ground (m)
  REAL(r8),                 INTENT(IN) :: FV     !friction velocity (m/s)
  REAL(r8),                 INTENT(IN) :: CWP    !canopy wind parameter
  REAL(r8),                 INTENT(IN) :: MPE    !prevents overflow error if division by zero

! in & out

  REAL(r8),              INTENT(INOUT) :: MOZG   !Monin-Obukhov stability parameter
  REAL(r8),              INTENT(INOUT) :: FHG    !stability correction

! outputs
  REAL(r8)                             :: RAMG   !aerodynamic resistance for momentum (s/m)
  REAL(r8)                             :: RAHG   !aerodynamic resistance for sensible heat (s/m)
  REAL(r8)                             :: RAWG   !aerodynamic resistance for water vapor (s/m)
  REAL(r8)                             :: RB     !bulk leaf boundary layer resistance (s/m)


  REAL(r8) :: KH           !turbulent transfer coefficient, sensible heat, (m2/s)
  REAL(r8) :: TMP1         !temporary calculation
  REAL(r8) :: TMP2         !temporary calculation
  REAL(r8) :: TMPRAH2      !temporary calculation for aerodynamic resistances
  REAL(r8) :: TMPRB        !temporary calculation for rb
  real :: MOLG,FHGNEW,CWPC
! --------------------------------------------------------------------------------------------------
! stability correction to below canopy resistance

       MOZG = 0.
       MOLG = 0.

       IF(ITER > 1) THEN
        TMP1 = VKC * (GRAV/TAH) * HG/(RHOAIR*CPAIR)
        IF (ABS(TMP1) .LE. MPE) TMP1 = MPE
        MOLG = -1. * FV**3 / TMP1
        MOZG = MIN( (ZPD-Z0MG)/MOLG, 1.)
       END IF

       IF (MOZG < 0.) THEN
          FHGNEW  = (1. - 15.*MOZG)**(-0.25)
       ELSE
          FHGNEW  = 1.+ 4.7*MOZG
       ENDIF

       IF (ITER == 1) THEN
          FHG = FHGNEW
       ELSE
          FHG = 0.5 * (FHG+FHGNEW)
       ENDIF

       CWPC = (CWP * VAI * HCAN * FHG)**0.5
!       CWPC = (CWP*FHG)**0.5

       TMP1 = EXP( -CWPC*Z0HG/HCAN )
       TMP2 = EXP( -CWPC*(Z0H+ZPD)/HCAN )
       TMPRAH2 = HCAN*EXP(CWPC) / CWPC * (TMP1-TMP2)

! aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.

       KH  = MAX ( VKC*FV*(HCAN-ZPD), MPE )
       RAMG = 0.
       RAHG = TMPRAH2 / KH
       RAWG = RAHG

! leaf boundary layer resistance

       TMPRB  = CWPC*50. / (1. - EXP(-CWPC/2.))
       RB     = TMPRB * SQRT(parameters%DLEAF/UC)
!       RB = 200

  END SUBROUTINE RAGRB

!== begin sfcdif1 ==================================================================================

  SUBROUTINE SFCDIF1(parameters,ITER   ,SFCTMP ,RHOAIR ,H      ,QAIR   , & !in
       &             ZLVL   ,ZPD    ,Z0M    ,Z0H    ,UR     , & !in
       &             MPE    ,ILOC   ,JLOC   ,                 & !in
       &             MOZ    ,MOZSGN ,FM     ,FH     ,FM2,FH2, & !inout
       &             CM     ,CH     ,FV     ,CH2     )          !out
! -------------------------------------------------------------------------------------------------
! computing surface drag coefficient CM for momentum and CH for heat
! -------------------------------------------------------------------------------------------------
    IMPLICIT NONE
! -------------------------------------------------------------------------------------------------
! inputs
    
  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,              INTENT(IN) :: ILOC   !grid index
    INTEGER,              INTENT(IN) :: JLOC   !grid index
    INTEGER,              INTENT(IN) :: ITER   !iteration index
    REAL(r8),                 INTENT(IN) :: SFCTMP !temperature at reference height (k)
    REAL(r8),                 INTENT(IN) :: RHOAIR !density air (kg/m**3)
    REAL(r8),                 INTENT(IN) :: H      !sensible heat flux (w/m2) [+ to atm]
    REAL(r8),                 INTENT(IN) :: QAIR   !specific humidity at reference height (kg/kg)
    REAL(r8),                 INTENT(IN) :: ZLVL   !reference height  (m)
    REAL(r8),                 INTENT(IN) :: ZPD    !zero plane displacement (m)
    REAL(r8),                 INTENT(IN) :: Z0H    !roughness length, sensible heat, ground (m)
    REAL(r8),                 INTENT(IN) :: Z0M    !roughness length, momentum, ground (m)
    REAL(r8),                 INTENT(IN) :: UR     !wind speed (m/s)
    REAL(r8),                 INTENT(IN) :: MPE    !prevents overflow error if division by zero
! in & out

    INTEGER,           INTENT(INOUT) :: MOZSGN !number of times moz changes sign
    REAL(r8),              INTENT(INOUT) :: MOZ    !Monin-Obukhov stability (z/L)
    REAL(r8),              INTENT(INOUT) :: FM     !momentum stability correction, weighted by prior iters
    REAL(r8),              INTENT(INOUT) :: FH     !sen heat stability correction, weighted by prior iters
    REAL(r8),              INTENT(INOUT) :: FM2    !sen heat stability correction, weighted by prior iters
    REAL(r8),              INTENT(INOUT) :: FH2    !sen heat stability correction, weighted by prior iters

! outputs

    REAL(r8),                INTENT(OUT) :: CM     !drag coefficient for momentum
    REAL(r8),                INTENT(OUT) :: CH     !drag coefficient for heat
    REAL(r8),                INTENT(OUT) :: FV     !friction velocity (m/s)
    REAL(r8),                INTENT(OUT) :: CH2    !drag coefficient for heat

! locals
    REAL(r8)    :: MOL                      !Monin-Obukhov length (m)
    REAL(r8)    :: TMPCM                    !temporary calculation for CM
    REAL(r8)    :: TMPCH                    !temporary calculation for CH
    REAL(r8)    :: FMNEW                    !stability correction factor, momentum, for current moz
    REAL(r8)    :: FHNEW                    !stability correction factor, sen heat, for current moz
    REAL(r8)    :: MOZOLD                   !Monin-Obukhov stability parameter from prior iteration
    REAL(r8)    :: TMP1,TMP2,TMP3,TMP4,TMP5 !temporary calculation
    REAL(r8)    :: TVIR                     !temporary virtual temperature (k)
    REAL(r8)    :: MOZ2                     !2/L
    REAL(r8)    :: TMPCM2                   !temporary calculation for CM2
    REAL(r8)    :: TMPCH2                   !temporary calculation for CH2
    REAL(r8)    :: FM2NEW                   !stability correction factor, momentum, for current moz
    REAL(r8)    :: FH2NEW                   !stability correction factor, sen heat, for current moz
    REAL(r8)    :: TMP12,TMP22,TMP32        !temporary calculation

    REAL(r8)    :: CMFM, CHFH, CM2FM2, CH2FH2
! -------------------------------------------------------------------------------------------------
! Monin-Obukhov stability parameter moz for next iteration

    MOZOLD = MOZ
  
    IF(ZLVL <= ZPD) THEN
       write(*,*) 'critical problem: ZLVL <= ZPD; model stops'
       print*,("STOP in Noah-MP")
    ENDIF

    TMPCM = LOG((ZLVL-ZPD) / Z0M)
    TMPCH = LOG((ZLVL-ZPD) / Z0H)
    TMPCM2 = LOG((2.0 + Z0M) / Z0M)
    TMPCH2 = LOG((2.0 + Z0H) / Z0H)

    IF(ITER == 1) THEN
       FV   = 0.0
       MOZ  = 0.0
       MOL  = 0.0
       MOZ2 = 0.0
    ELSE
       TVIR = (1. + 0.61*QAIR) * SFCTMP
       TMP1 = VKC * (GRAV/TVIR) * H/(RHOAIR*CPAIR)
       IF (ABS(TMP1) .LE. MPE) TMP1 = MPE
       MOL  = -1. * FV**3 / TMP1
       MOZ  = MIN( (ZLVL-ZPD)/MOL, 1.)
       MOZ2  = MIN( (2.0 + Z0H)/MOL, 1.)
    ENDIF

! accumulate number of times moz changes sign.

    IF (MOZOLD*MOZ .LT. 0.) MOZSGN = MOZSGN+1
    IF (MOZSGN .GE. 2) THEN
       MOZ = 0.
       FM = 0.
       FH = 0.
       MOZ2 = 0.
       FM2 = 0.
       FH2 = 0.
    ENDIF

! evaluate stability-dependent variables using moz from prior iteration
    IF (MOZ .LT. 0.) THEN
       TMP1 = (1. - 16.*MOZ)**0.25
       TMP2 = LOG((1.+TMP1*TMP1)/2.)
       TMP3 = LOG((1.+TMP1)/2.)
       FMNEW = 2.*TMP3 + TMP2 - 2.*ATAN(TMP1) + 1.5707963
       FHNEW = 2*TMP2

! 2-meter
       TMP12 = (1. - 16.*MOZ2)**0.25
       TMP22 = LOG((1.+TMP12*TMP12)/2.)
       TMP32 = LOG((1.+TMP12)/2.)
       FM2NEW = 2.*TMP32 + TMP22 - 2.*ATAN(TMP12) + 1.5707963
       FH2NEW = 2*TMP22
    ELSE
       FMNEW = -5.*MOZ
       FHNEW = FMNEW
       FM2NEW = -5.*MOZ2
       FH2NEW = FM2NEW
    ENDIF

! except for first iteration, weight stability factors for previous
! iteration to help avoid flip-flops from one iteration to the next

    IF (ITER == 1) THEN
       FM = FMNEW
       FH = FHNEW
       FM2 = FM2NEW
       FH2 = FH2NEW
    ELSE
       FM = 0.5 * (FM+FMNEW)
       FH = 0.5 * (FH+FHNEW)
       FM2 = 0.5 * (FM2+FM2NEW)
       FH2 = 0.5 * (FH2+FH2NEW)
    ENDIF

! exchange coefficients

    FH = MIN(FH,0.9*TMPCH)
    FM = MIN(FM,0.9*TMPCM)
    FH2 = MIN(FH2,0.9*TMPCH2)
    FM2 = MIN(FM2,0.9*TMPCM2)

    CMFM = TMPCM-FM
    CHFH = TMPCH-FH
    CM2FM2 = TMPCM2-FM2
    CH2FH2 = TMPCH2-FH2
    IF(ABS(CMFM) <= MPE) CMFM = MPE
    IF(ABS(CHFH) <= MPE) CHFH = MPE
    IF(ABS(CM2FM2) <= MPE) CM2FM2 = MPE
    IF(ABS(CH2FH2) <= MPE) CH2FH2 = MPE
    CM  = VKC*VKC/(CMFM*CMFM)
    CH  = VKC*VKC/(CMFM*CHFH)
    CH2  = VKC*VKC/(CM2FM2*CH2FH2)
        
! friction velocity

    FV = UR * SQRT(CM)
    CH2  = VKC*FV/CH2FH2

  END SUBROUTINE SFCDIF1

!== begin sfcdif2 ==================================================================================

  SUBROUTINE SFCDIF2(parameters,ITER   ,Z0     ,THZ0   ,THLM   ,SFCSPD , & !in
                     ZLM    ,ILOC   ,JLOC   ,         & !in
                     AKMS   ,AKHS   ,RLMO   ,WSTAR2 ,         & !in
                     USTAR  )                                   !out

! -------------------------------------------------------------------------------------------------
! SUBROUTINE SFCDIF (renamed SFCDIF_off to avoid clash with Eta PBL)
! -------------------------------------------------------------------------------------------------
! CALCULATE SURFACE LAYER EXCHANGE COEFFICIENTS VIA ITERATIVE PROCESS.
! SEE CHEN ET AL (1997, BLM)
! -------------------------------------------------------------------------------------------------
    IMPLICIT NONE
  type (noahmp_parameters), intent(in) :: parameters
    INTEGER, INTENT(IN) :: ILOC
    INTEGER, INTENT(IN) :: JLOC
    INTEGER, INTENT(IN) :: ITER
    REAL(r8),    INTENT(IN) :: ZLM, Z0, THZ0, THLM, SFCSPD
    REAL(r8), intent(INOUT) :: AKMS
    REAL(r8), intent(INOUT) :: AKHS
    REAL(r8), intent(INOUT) :: RLMO
    REAL(r8), intent(INOUT) :: WSTAR2
    REAL(r8),   intent(OUT) :: USTAR

    REAL(r8)     ZZ, PSLMU, PSLMS, PSLHU, PSLHS
    REAL(r8)     XX, PSPMU, YY, PSPMS, PSPHU, PSPHS
    REAL(r8)     ZILFC, ZU, ZT, RDZ, CXCH
    REAL(r8)     DTHV, DU2, BTGH, ZSLU, ZSLT, RLOGU, RLOGT
    REAL(r8)     ZETALT, ZETALU, ZETAU, ZETAT, XLU4, XLT4, XU4, XT4

    REAL(r8)     XLU, XLT, XU, XT, PSMZ, SIMM, PSHZ, SIMH, USTARK, RLMN,  &
         &         RLMA

    INTEGER  ILECH, ITR

    INTEGER, PARAMETER :: ITRMX  = 5
    REAL(r8),    PARAMETER :: WWST   = 1.2
    REAL(r8),    PARAMETER :: WWST2  = WWST * WWST
    REAL(r8),    PARAMETER :: VKRM   = 0.40
    REAL(r8),    PARAMETER :: EXCM   = 0.001
    REAL(r8),    PARAMETER :: BETA   = 1.0 / 270.0
    REAL(r8),    PARAMETER :: BTG    = BETA * GRAV
    REAL(r8),    PARAMETER :: ELFC   = VKRM * BTG
    REAL(r8),    PARAMETER :: WOLD   = 0.15
    REAL(r8),    PARAMETER :: WNEW   = 1.0 - WOLD
    REAL(r8),    PARAMETER :: PIHF   = 3.14159265 / 2.
    REAL(r8),    PARAMETER :: EPSU2  = 1.E-4
    REAL(r8),    PARAMETER :: EPSUST = 0.07
    REAL(r8),    PARAMETER :: EPSIT  = 1.E-4
    REAL(r8),    PARAMETER :: EPSA   = 1.E-8
    REAL(r8),    PARAMETER :: ZTMIN  = -5.0
    REAL(r8),    PARAMETER :: ZTMAX  = 1.0
    REAL(r8),    PARAMETER :: HPBL   = 1000.0
    REAL(r8),    PARAMETER :: SQVISC = 258.2
    REAL(r8),    PARAMETER :: RIC    = 0.183
    REAL(r8),    PARAMETER :: RRIC   = 1.0 / RIC
    REAL(r8),    PARAMETER :: FHNEU  = 0.8
    REAL(r8),    PARAMETER :: RFC    = 0.191
    REAL(r8),    PARAMETER :: RFAC   = RIC / ( FHNEU * RFC * RFC )

! ----------------------------------------------------------------------
! NOTE: THE TWO CODE BLOCKS BELOW DEFINE FUNCTIONS
! ----------------------------------------------------------------------
! LECH'S SURFACE FUNCTIONS
    PSLMU (ZZ)= -0.96* log (1.0-4.5* ZZ)
    PSLMS (ZZ)= ZZ * RRIC -2.076* (1. -1./ (ZZ +1.))
    PSLHU (ZZ)= -0.96* log (1.0-4.5* ZZ)
    PSLHS (ZZ)= ZZ * RFAC -2.076* (1. -1./ (ZZ +1.))
! PAULSON'S SURFACE FUNCTIONS
    PSPMU (XX)= -2.* log ( (XX +1.)*0.5) - log ( (XX * XX +1.)*0.5)   &
         &        +2.* ATAN (XX)                                            &
         &- PIHF
    PSPMS (YY)= 5.* YY
    PSPHU (XX)= -2.* log ( (XX * XX +1.)*0.5)
    PSPHS (YY)= 5.* YY

! THIS ROUTINE SFCDIF CAN HANDLE BOTH OVER OPEN WATER (SEA, OCEAN) AND
! OVER SOLID SURFACE (LAND, SEA-ICE).
! ----------------------------------------------------------------------
!     ZTFC: RATIO OF ZOH/ZOM  LESS OR EQUAL THAN 1
!     C......ZTFC=0.1
!     CZIL: CONSTANT C IN Zilitinkevich, S. S.1995,:NOTE ABOUT ZT
! ----------------------------------------------------------------------
    ILECH = 0

! ----------------------------------------------------------------------
    ZILFC = - parameters%CZIL * VKRM * SQVISC
    ZU = Z0
    RDZ = 1./ ZLM
    CXCH = EXCM * RDZ
    DTHV = THLM - THZ0

! BELJARS CORRECTION OF USTAR
    DU2 = MAX (SFCSPD * SFCSPD,EPSU2)
    BTGH = BTG * HPBL

    IF(ITER == 1) THEN
        IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
           WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
        ELSE
           WSTAR2 = 0.0
        END IF
        USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)
        RLMO = ELFC * AKHS * DTHV / USTAR **3
    END IF
 
! ZILITINKEVITCH APPROACH FOR ZT
    ZT = MAX(1.E-6,EXP (ZILFC * SQRT (USTAR * Z0))* Z0)
    ZSLU = ZLM + ZU
    ZSLT = ZLM + ZT
    RLOGU = log (ZSLU / ZU)
    RLOGT = log (ZSLT / ZT)

! ----------------------------------------------------------------------
! 1./MONIN-OBUKKHOV LENGTH-SCALE
! ----------------------------------------------------------------------
    ZETALT = MAX (ZSLT * RLMO,ZTMIN)
    RLMO = ZETALT / ZSLT
    ZETALU = ZSLU * RLMO
    ZETAU = ZU * RLMO
    ZETAT = ZT * RLMO

    IF (ILECH .eq. 0) THEN
       IF (RLMO .lt. 0.)THEN
          XLU4 = 1. -16.* ZETALU
          XLT4 = 1. -16.* ZETALT
          XU4  = 1. -16.* ZETAU
          XT4  = 1. -16.* ZETAT
          XLU  = SQRT (SQRT (XLU4))
          XLT  = SQRT (SQRT (XLT4))
          XU   = SQRT (SQRT (XU4))

          XT = SQRT (SQRT (XT4))
          PSMZ = PSPMU (XU)
          SIMM = PSPMU (XLU) - PSMZ + RLOGU
          PSHZ = PSPHU (XT)
          SIMH = PSPHU (XLT) - PSHZ + RLOGT
       ELSE
          ZETALU = MIN (ZETALU,ZTMAX)
          ZETALT = MIN (ZETALT,ZTMAX)
          ZETAU  = MIN (ZETAU,ZTMAX/(ZSLU/ZU))   ! Barlage: add limit on ZETAU/ZETAT
          ZETAT  = MIN (ZETAT,ZTMAX/(ZSLT/ZT))   ! Barlage: prevent SIMM/SIMH < 0
          PSMZ = PSPMS (ZETAU)
          SIMM = PSPMS (ZETALU) - PSMZ + RLOGU
          PSHZ = PSPHS (ZETAT)
          SIMH = PSPHS (ZETALT) - PSHZ + RLOGT
       END IF
! ----------------------------------------------------------------------
! LECH'S FUNCTIONS
! ----------------------------------------------------------------------
    ELSE
       IF (RLMO .lt. 0.)THEN
          PSMZ = PSLMU (ZETAU)
          SIMM = PSLMU (ZETALU) - PSMZ + RLOGU
          PSHZ = PSLHU (ZETAT)
          SIMH = PSLHU (ZETALT) - PSHZ + RLOGT
       ELSE
          ZETALU = MIN (ZETALU,ZTMAX)
          ZETALT = MIN (ZETALT,ZTMAX)
          PSMZ = PSLMS (ZETAU)
          SIMM = PSLMS (ZETALU) - PSMZ + RLOGU
          PSHZ = PSLHS (ZETAT)
          SIMH = PSLHS (ZETALT) - PSHZ + RLOGT
       END IF
! ----------------------------------------------------------------------
       END IF

! ----------------------------------------------------------------------
! BELJAARS CORRECTION FOR USTAR
! ----------------------------------------------------------------------
       USTAR = MAX (SQRT (AKMS * SQRT (DU2+ WSTAR2)),EPSUST)

! ZILITINKEVITCH FIX FOR ZT
       ZT = MAX(1.E-6,EXP (ZILFC * SQRT (USTAR * Z0))* Z0)
       ZSLT = ZLM + ZT
!-----------------------------------------------------------------------
       RLOGT = log (ZSLT / ZT)
       USTARK = USTAR * VKRM
       IF(SIMM < 1.e-6) SIMM = 1.e-6        ! Limit stability function
       AKMS = MAX (USTARK / SIMM,CXCH)
!-----------------------------------------------------------------------
! IF STATEMENTS TO AVOID TANGENT LINEAR PROBLEMS NEAR ZERO
!-----------------------------------------------------------------------
       IF(SIMH < 1.e-6) SIMH = 1.e-6        ! Limit stability function
       AKHS = MAX (USTARK / SIMH,CXCH)

       IF (BTGH * AKHS * DTHV .ne. 0.0) THEN
          WSTAR2 = WWST2* ABS (BTGH * AKHS * DTHV)** (2./3.)
       ELSE
          WSTAR2 = 0.0
       END IF
!-----------------------------------------------------------------------
       RLMN = ELFC * AKHS * DTHV / USTAR **3
!-----------------------------------------------------------------------
!     IF(ABS((RLMN-RLMO)/RLMA).LT.EPSIT)    GO TO 110
!-----------------------------------------------------------------------
       RLMA = RLMO * WOLD+ RLMN * WNEW
!-----------------------------------------------------------------------
       RLMO = RLMA

!       write(*,'(a20,10f15.6)')'SFCDIF: RLMO=',RLMO,RLMN,ELFC , AKHS , DTHV , USTAR
!    END DO
! ----------------------------------------------------------------------
  END SUBROUTINE SFCDIF2

!== begin esat =====================================================================================

  SUBROUTINE ESAT(T, ESW, ESI, DESW, DESI)
!---------------------------------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
  IMPLICIT NONE
!---------------------------------------------------------------------------------------------------
! in

  REAL(r8), intent(in)  :: T              !temperature

!out

  REAL(r8), intent(out) :: ESW            !saturation vapor pressure over water (pa)
  REAL(r8), intent(out) :: ESI            !saturation vapor pressure over ice (pa)
  REAL(r8), intent(out) :: DESW           !d(esat)/dt over water (pa/K)
  REAL(r8), intent(out) :: DESI           !d(esat)/dt over ice (pa/K)

! local

  REAL(r8) :: A0,A1,A2,A3,A4,A5,A6  !coefficients for esat over water
  REAL(r8) :: B0,B1,B2,B3,B4,B5,B6  !coefficients for esat over ice
  REAL(r8) :: C0,C1,C2,C3,C4,C5,C6  !coefficients for dsat over water
  REAL(r8) :: D0,D1,D2,D3,D4,D5,D6  !coefficients for dsat over ice

  PARAMETER (A0=6.107799961    , A1=4.436518521E-01,  &
             A2=1.428945805E-02, A3=2.650648471E-04,  &
             A4=3.031240396E-06, A5=2.034080948E-08,  &
             A6=6.136820929E-11)

  PARAMETER (B0=6.109177956    , B1=5.034698970E-01,  &
             B2=1.886013408E-02, B3=4.176223716E-04,  &
             B4=5.824720280E-06, B5=4.838803174E-08,  &
             B6=1.838826904E-10)

  PARAMETER (C0= 4.438099984E-01, C1=2.857002636E-02,  &
             C2= 7.938054040E-04, C3=1.215215065E-05,  &
             C4= 1.036561403E-07, C5=3.532421810e-10,  &
             C6=-7.090244804E-13)

  PARAMETER (D0=5.030305237E-01, D1=3.773255020E-02,  &
             D2=1.267995369E-03, D3=2.477563108E-05,  &
             D4=3.005693132E-07, D5=2.158542548E-09,  &
             D6=7.131097725E-12)

  ESW  = 100.*(A0+T*(A1+T*(A2+T*(A3+T*(A4+T*(A5+T*A6))))))
  ESI  = 100.*(B0+T*(B1+T*(B2+T*(B3+T*(B4+T*(B5+T*B6))))))
  DESW = 100.*(C0+T*(C1+T*(C2+T*(C3+T*(C4+T*(C5+T*C6))))))
  DESI = 100.*(D0+T*(D1+T*(D2+T*(D3+T*(D4+T*(D5+T*D6))))))

  END SUBROUTINE ESAT

!== begin stomata ==================================================================================

  SUBROUTINE STOMATA (parameters,VEGTYP  ,MPE     ,APAR    ,FOLN    ,ILOC    , JLOC, & !in
                      TV      ,EI      ,EA      ,SFCTMP  ,SFCPRS  , & !in
                      O2      ,CO2     ,IGS     ,BTRAN   ,RB      , & !in
                      RS      ,PSN     )                              !out
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
      INTEGER,INTENT(IN)  :: ILOC   !grid index
      INTEGER,INTENT(IN)  :: JLOC   !grid index
      INTEGER,INTENT(IN)  :: VEGTYP !vegetation physiology type

      REAL(r8), INTENT(IN)    :: IGS    !growing season index (0=off, 1=on)
      REAL(r8), INTENT(IN)    :: MPE    !prevents division by zero errors

      REAL(r8), INTENT(IN)    :: TV     !foliage temperature (k)
      REAL(r8), INTENT(IN)    :: EI     !vapor pressure inside leaf (sat vapor press at tv) (pa)
      REAL(r8), INTENT(IN)    :: EA     !vapor pressure of canopy air (pa)
      REAL(r8), INTENT(IN)    :: APAR   !par absorbed per unit lai (w/m2)
      REAL(r8), INTENT(IN)    :: O2     !atmospheric o2 concentration (pa)
      REAL(r8), INTENT(IN)    :: CO2    !atmospheric co2 concentration (pa)
      REAL(r8), INTENT(IN)    :: SFCPRS !air pressure at reference height (pa)
      REAL(r8), INTENT(IN)    :: SFCTMP !air temperature at reference height (k)
      REAL(r8), INTENT(IN)    :: BTRAN  !soil water transpiration factor (0 to 1)
      REAL(r8), INTENT(IN)    :: FOLN   !foliage nitrogen concentration (%)
      REAL(r8), INTENT(IN)    :: RB     !boundary layer resistance (s/m)

! output
      REAL(r8), INTENT(OUT)   :: RS     !leaf stomatal resistance (s/m)
      REAL(r8), INTENT(OUT)   :: PSN    !foliage photosynthesis (umol co2 /m2/ s) [always +]

! in&out
      REAL(r8)                :: RLB    !boundary layer resistance (s m2 / umol)
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
      INTEGER :: ITER     !iteration index
      INTEGER :: NITER    !number of iterations

      DATA NITER /3/
      SAVE NITER

      REAL(r8) :: AB          !used in statement functions
      REAL(r8) :: BC          !used in statement functions
      REAL(r8) :: F1          !generic temperature response (statement function)
      REAL(r8) :: F2          !generic temperature inhibition (statement function)
      REAL(r8) :: TC          !foliage temperature (degree Celsius)
      REAL(r8) :: CS          !co2 concentration at leaf surface (pa)
      REAL(r8) :: KC          !co2 Michaelis-Menten constant (pa)
      REAL(r8) :: KO          !o2 Michaelis-Menten constant (pa)
      REAL(r8) :: A,B,C,Q     !intermediate calculations for RS
      REAL(r8) :: R1,R2       !roots for RS
      REAL(r8) :: FNF         !foliage nitrogen adjustment factor (0 to 1)
      REAL(r8) :: PPF         !absorb photosynthetic photon flux (umol photons/m2/s)
      REAL(r8) :: WC          !Rubisco limited photosynthesis (umol co2/m2/s)
      REAL(r8) :: WJ          !light limited photosynthesis (umol co2/m2/s)
      REAL(r8) :: WE          !export limited photosynthesis (umol co2/m2/s)
      REAL(r8) :: CP          !co2 compensation point (pa)
      REAL(r8) :: CI          !internal co2 (pa)
      REAL(r8) :: AWC         !intermediate calculation for wc
      REAL(r8) :: VCMX        !maximum rate of carbonylation (umol co2/m2/s)
      REAL(r8) :: J           !electron transport (umol co2/m2/s)
      REAL(r8) :: CEA         !constrain ea or else model blows up
      REAL(r8) :: CF          !s m2/umol -> s/m

      F1(AB,BC) = AB**((BC-25.)/10.)
      F2(AB) = 1. + EXP((-2.2E05+710.*(AB+273.16))/(8.314*(AB+273.16)))
      REAL(r8) :: T
! ---------------------------------------------------------------------------------------------

! initialize RS=RSMAX and PSN=0 because will only do calculations
! for APAR > 0, in which case RS <= RSMAX and PSN >= 0

         CF = SFCPRS/(8.314*SFCTMP)*1.e06
         RS = 1./parameters%BP * CF
         PSN = 0.

         IF (APAR .LE. 0.) RETURN

         FNF = MIN( FOLN/MAX(MPE,parameters%FOLNMX), 1.0 )
         TC  = TV-TFRZ
         PPF = 4.6*APAR
         J   = PPF*parameters%QE25
         KC  = parameters%KC25 * F1(parameters%AKC,TC)
         KO  = parameters%KO25 * F1(parameters%AKO,TC)
         AWC = KC * (1.+O2/KO)
         CP  = 0.5*KC/KO*O2*0.21
         VCMX = parameters%VCMX25 / F2(TC) * FNF * BTRAN * F1(parameters%AVCMX,TC)

! first guess ci

         CI = 0.7*CO2*parameters%C3PSN + 0.4*CO2*(1.-parameters%C3PSN)

! rb: s/m -> s m**2 / umol

         RLB = RB/CF

! constrain ea

         CEA = MAX(0.25*EI*parameters%C3PSN+0.40*EI*(1.-parameters%C3PSN), MIN(EA,EI) )

! ci iteration
!jref: C3PSN is equal to 1 for all veg types.
       DO ITER = 1, NITER
            WJ = MAX(CI-CP,0.)*J/(CI+2.*CP)*parameters%C3PSN  + J*(1.-parameters%C3PSN)
            WC = MAX(CI-CP,0.)*VCMX/(CI+AWC)*parameters%C3PSN + VCMX*(1.-parameters%C3PSN)
            WE = 0.5*VCMX*parameters%C3PSN + 4000.*VCMX*CI/SFCPRS*(1.-parameters%C3PSN)
            PSN = MIN(WJ,WC,WE) * IGS

            CS = MAX( CO2-1.37*RLB*SFCPRS*PSN, MPE )
            A = parameters%MP*PSN*SFCPRS*CEA / (CS*EI) + parameters%BP
            B = ( parameters%MP*PSN*SFCPRS/CS + parameters%BP ) * RLB - 1.
            C = -RLB
            IF (B .GE. 0.) THEN
               Q = -0.5*( B + SQRT(B*B-4.*A*C) )
            ELSE
               Q = -0.5*( B - SQRT(B*B-4.*A*C) )
            END IF
            R1 = Q/A
            R2 = C/Q
            RS = MAX(R1,R2)
            CI = MAX( CS-PSN*SFCPRS*1.65*RS, 0. )
       END DO 

! rs, rb:  s m**2 / umol -> s/m

         RS = RS*CF

  END SUBROUTINE STOMATA

!== begin canres ===================================================================================

  SUBROUTINE CANRES (parameters,PAR   ,SFCTMP,RCSOIL ,EAH   ,SFCPRS , & !in
                     RC    ,PSN   ,ILOC   ,JLOC  )           !out

! --------------------------------------------------------------------------------------------------
! calculate canopy resistance which depends on incoming solar radiation,
! air temperature, atmospheric water vapor pressure deficit at the
! lowest model level, and soil moisture (preferably unfrozen soil
! moisture rather than total)
! --------------------------------------------------------------------------------------------------
! source:  Jarvis (1976), Noilhan and Planton (1989, MWR), Jacquemin and
! Noilhan (1990, BLM). Chen et al (1996, JGR, Vol 101(D3), 7251-7268), 
! eqns 12-14 and table 2 of sec. 3.1.2
! --------------------------------------------------------------------------------------------------
!niu    USE module_Noahlsm_utility
! --------------------------------------------------------------------------------------------------
    IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,                  INTENT(IN)  :: ILOC   !grid index
    INTEGER,                  INTENT(IN)  :: JLOC   !grid index
    REAL(r8),                     INTENT(IN)  :: PAR    !par absorbed per unit sunlit lai (w/m2)
    REAL(r8),                     INTENT(IN)  :: SFCTMP !canopy air temperature
    REAL(r8),                     INTENT(IN)  :: SFCPRS !surface pressure (pa)
    REAL(r8),                     INTENT(IN)  :: EAH    !water vapor pressure (pa)
    REAL(r8),                     INTENT(IN)  :: RCSOIL !soil moisture stress factor

!outputs

    REAL(r8),                     INTENT(OUT) :: RC     !canopy resistance per unit LAI
    REAL(r8),                     INTENT(OUT) :: PSN    !foliage photosynthesis (umolco2/m2/s)

!local

    REAL(r8)                                  :: RCQ
    REAL(r8)                                  :: RCS
    REAL(r8)                                  :: RCT
    REAL(r8)                                  :: FF
    REAL(r8)                                  :: Q2     !water vapor mixing ratio (kg/kg)
    REAL(r8)                                  :: Q2SAT  !saturation Q2
    REAL(r8)                                  :: DQSDT2 !d(Q2SAT)/d(T)

! RSMIN, RSMAX, TOPT, RGL, HS are canopy stress parameters set in REDPRM
! ----------------------------------------------------------------------
! initialize canopy resistance multiplier terms.
! ----------------------------------------------------------------------
    RC     = 0.0
    RCS    = 0.0
    RCT    = 0.0
    RCQ    = 0.0

!  compute Q2 and Q2SAT

    Q2 = 0.622 *  EAH  / (SFCPRS - 0.378 * EAH) !specific humidity [kg/kg]
    Q2 = Q2 / (1.0 + Q2)                        !mixing ratio [kg/kg]

    CALL CALHUM(parameters,SFCTMP, SFCPRS, Q2SAT, DQSDT2)

! contribution due to incoming solar radiation

    FF  = 2.0 * PAR / parameters%RGL                
    RCS = (FF + parameters%RSMIN / parameters%RSMAX) / (1.0+ FF)
    RCS = MAX (RCS,0.0001)

! contribution due to air temperature

    RCT = 1.0- 0.0016* ( (parameters%TOPT - SFCTMP)**2.0)
    RCT = MAX (RCT,0.0001)

! contribution due to vapor pressure deficit

    RCQ = 1.0/ (1.0+ parameters%HS * MAX(0.,Q2SAT-Q2))
    RCQ = MAX (RCQ,0.01)

! determine canopy resistance due to all factors

    RC  = parameters%RSMIN / (RCS * RCT * RCQ * RCSOIL)
    PSN = -999.99       ! PSN not applied for dynamic carbon

  END SUBROUTINE CANRES

!== begin calhum ===================================================================================

        SUBROUTINE CALHUM(parameters,SFCTMP, SFCPRS, Q2SAT, DQSDT2)

        IMPLICIT NONE

  type (noahmp_parameters), intent(in) :: parameters
        REAL(r8), INTENT(IN)       :: SFCTMP, SFCPRS
        REAL(r8), INTENT(OUT)      :: Q2SAT, DQSDT2
        REAL(r8), PARAMETER        :: A2=17.67,A3=273.15,A4=29.65, ELWV=2.501E6,         &
                                  A23M4=A2*(A3-A4), E0=0.611, RV=461.0,             &
                                  EPSILON=0.622
        REAL(r8)                   :: ES, SFCPRSX

! Q2SAT: saturated mixing ratio
        ES = E0 * EXP ( ELWV/RV*(1./A3 - 1./SFCTMP) )
! convert SFCPRS from Pa to KPa
        SFCPRSX = SFCPRS*1.E-3
        Q2SAT = EPSILON * ES / (SFCPRSX-ES)
! convert from  g/g to g/kg
        Q2SAT = Q2SAT * 1.E3
! Q2SAT is currently a 'mixing ratio'

! DQSDT2 is calculated assuming Q2SAT is a specific humidity
        DQSDT2=(Q2SAT/(1+Q2SAT))*A23M4/(SFCTMP-A4)**2

! DG Q2SAT needs to be in g/g when returned for SFLX
        Q2SAT = Q2SAT / 1.E3

        END SUBROUTINE CALHUM

!== begin tsnosoi ==================================================================================

  SUBROUTINE TSNOSOI (parameters,ICE     ,NSOIL   ,NSNOW   ,ISNOW   ,IST     , & !in
                      TBOT    ,ZSNSO   ,SSOIL   ,DF      ,HCPCT   , & !in
                      SAG     ,DT      ,SNOWH   ,DZSNSO  , & !in
                      TG      ,ILOC    ,JLOC    ,                   & !in
                      STC     )                                       !inout
! --------------------------------------------------------------------------------------------------
! Compute snow (up to 3L) and soil (4L) temperature. Note that snow temperatures
! during melting season may exceed melting point (TFRZ) but later in PHASECHANGE
! subroutine the snow temperatures are reset to TFRZ for melting snow.
! --------------------------------------------------------------------------------------------------
  IMPLICIT NONE
! --------------------------------------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,                         INTENT(IN)  :: ILOC
    INTEGER,                         INTENT(IN)  :: JLOC
    INTEGER,                         INTENT(IN)  :: ICE    !
    INTEGER,                         INTENT(IN)  :: NSOIL  !no of soil layers (4)
    INTEGER,                         INTENT(IN)  :: NSNOW  !maximum no of snow layers (3)
    INTEGER,                         INTENT(IN)  :: ISNOW  !actual no of snow layers
    INTEGER,                         INTENT(IN)  :: IST    !surface type

    REAL(r8),                            INTENT(IN)  :: DT     !time step (s)
    REAL(r8),                            INTENT(IN)  :: TBOT   !
    REAL(r8),                            INTENT(IN)  :: SSOIL  !ground heat flux (w/m2)
    REAL(r8),                            INTENT(IN)  :: SAG    !solar rad. absorbed by ground (w/m2)
    REAL(r8),                            INTENT(IN)  :: SNOWH  !snow depth (m)
    REAL(r8),                            INTENT(IN)  :: TG     !ground temperature (k)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: ZSNSO  !layer-bot. depth from snow surf.(m)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DZSNSO !snow/soil layer thickness (m)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DF     !thermal conductivity
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: HCPCT  !heat capacity (J/m3/k)

!input and output

    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC

!local

    INTEGER                                      :: IZ
    REAL(r8)                                         :: ZBOTSNO   !ZBOT from snow surface
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)              :: AI, BI, CI, RHSTS
    REAL(r8)                                         :: EFLXB !energy influx from soil bottom (w/m2)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)              :: PHI   !light through water (w/m2)

    REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: TBEG
    REAL(r8)                            :: ERR_EST !heat storage error  (w/m2)
    REAL(r8)                            :: SSOIL2  !ground heat flux (w/m2) (for energy check)
    REAL(r8)                            :: EFLXB2  !heat flux from the bottom (w/m2) (for energy check)
    character(len=256)              :: message
! ----------------------------------------------------------------------
! compute solar penetration through water, needs more work

    PHI(ISNOW+1:NSOIL) = 0.

! adjust ZBOT from soil surface to ZBOTSNO from snow surface

    ZBOTSNO = parameters%ZBOT - SNOWH    !from snow surface

! snow/soil heat storage for energy balance check

    DO IZ = ISNOW+1, NSOIL
       TBEG(IZ) = STC(IZ)
    ENDDO

! compute soil temperatures

      CALL HRT   (parameters,NSNOW     ,NSOIL     ,ISNOW     ,ZSNSO     , &
                  STC       ,TBOT      ,ZBOTSNO   ,DT        , &
                  DF        ,HCPCT     ,SSOIL     ,PHI       , &
                  AI        ,BI        ,CI        ,RHSTS     , &
                  EFLXB     )

      CALL HSTEP (parameters,NSNOW     ,NSOIL     ,ISNOW     ,DT        , &
                  AI        ,BI        ,CI        ,RHSTS     , &
                  STC       ) 

! update ground heat flux just for energy check, but not for final output
! otherwise, it would break the surface energy balance

    IF(OPT_TBOT == 1) THEN
       EFLXB2  = 0.
    ELSE IF(OPT_TBOT == 2) THEN
       EFLXB2  = DF(NSOIL)*(TBOT-STC(NSOIL)) / &
            (0.5*(ZSNSO(NSOIL-1)+ZSNSO(NSOIL)) - ZBOTSNO)
    END IF

    ! Skip the energy balance check for now, until we can make it work
    ! right for small time steps.
    return

! energy balance check

    ERR_EST = 0.0
    DO IZ = ISNOW+1, NSOIL
       ERR_EST = ERR_EST + (STC(IZ)-TBEG(IZ)) * DZSNSO(IZ) * HCPCT(IZ) / DT
    ENDDO

    if (OPT_STC == 1 .OR. OPT_STC == 3) THEN   ! semi-implicit
       ERR_EST = ERR_EST - (SSOIL +EFLXB)
    ELSE                     ! full-implicit
       SSOIL2 = DF(ISNOW+1)*(TG-STC(ISNOW+1))/(0.5*DZSNSO(ISNOW+1))   !M. Barlage
       ERR_EST = ERR_EST - (SSOIL2+EFLXB2)
    ENDIF

    IF (ABS(ERR_EST) > 1.) THEN    ! W/m2
       WRITE(message,*) 'TSNOSOI is losing(-)/gaining(+) false energy',ERR_EST,' W/m2'
       print*, (trim(message))
       WRITE(message,'(i6,1x,i6,1x,i3,F18.13,5F20.12)') &
            ILOC, JLOC, IST,ERR_EST,SSOIL,SNOWH,TG,STC(ISNOW+1),EFLXB
       print*, (trim(message))
       !niu      STOP
    END IF

  END SUBROUTINE TSNOSOI

!== begin hrt ======================================================================================

  SUBROUTINE HRT (parameters,NSNOW     ,NSOIL     ,ISNOW     ,ZSNSO     , &
                  STC       ,TBOT      ,ZBOT      ,DT        , &
                  DF        ,HCPCT     ,SSOIL     ,PHI       , &
                  AI        ,BI        ,CI        ,RHSTS     , &
                  BOTFLX    )
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! calculate the right hand side of the time tendency term of the soil
! thermal diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,                         INTENT(IN)  :: NSOIL  !no of soil layers (4)
    INTEGER,                         INTENT(IN)  :: NSNOW  !maximum no of snow layers (3)
    INTEGER,                         INTENT(IN)  :: ISNOW  !actual no of snow layers
    REAL(r8),                            INTENT(IN)  :: TBOT   !bottom soil temp. at ZBOT (k)
    REAL(r8),                            INTENT(IN)  :: ZBOT   !depth of lower boundary condition (m)
                                                           !from soil surface not snow surface
    REAL(r8),                            INTENT(IN)  :: DT     !time step (s)
    REAL(r8),                            INTENT(IN)  :: SSOIL  !ground heat flux (w/m2)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: ZSNSO  !depth of layer-bottom of snow/soil (m)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: STC    !snow/soil temperature (k)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: DF     !thermal conductivity [w/m/k]
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: HCPCT  !heat capacity [j/m3/k]
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)  :: PHI    !light through water (w/m2)

! output

    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: RHSTS  !right-hand side of the matrix
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: AI     !left-hand side coefficient
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: BI     !left-hand side coefficient
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: CI     !left-hand side coefficient
    REAL(r8),                            INTENT(OUT) :: BOTFLX !energy influx from soil bottom (w/m2)

! local

    INTEGER                                      :: K
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)              :: DDZ
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)              :: DZ
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)              :: DENOM
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)              :: DTSDZ
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)              :: EFLUX
    REAL(r8)                                         :: TEMP1
! ----------------------------------------------------------------------

    DO K = ISNOW+1, NSOIL
        IF (K == ISNOW+1) THEN
           DENOM(K)  = - ZSNSO(K) * HCPCT(K)
           TEMP1     = - ZSNSO(K+1)
           DDZ(K)    = 2.0 / TEMP1
           DTSDZ(K)  = 2.0 * (STC(K) - STC(K+1)) / TEMP1
           EFLUX(K)  = DF(K) * DTSDZ(K) - SSOIL - PHI(K)
        ELSE IF (K < NSOIL) THEN
           DENOM(K)  = (ZSNSO(K-1) - ZSNSO(K)) * HCPCT(K)
           TEMP1     = ZSNSO(K-1) - ZSNSO(K+1)
           DDZ(K)    = 2.0 / TEMP1
           DTSDZ(K)  = 2.0 * (STC(K) - STC(K+1)) / TEMP1
           EFLUX(K)  = (DF(K)*DTSDZ(K) - DF(K-1)*DTSDZ(K-1)) - PHI(K)
        ELSE IF (K == NSOIL) THEN
           DENOM(K)  = (ZSNSO(K-1) - ZSNSO(K)) * HCPCT(K)
           TEMP1     =  ZSNSO(K-1) - ZSNSO(K)
           IF(OPT_TBOT == 1) THEN
               BOTFLX     = 0. 
           END IF
           IF(OPT_TBOT == 2) THEN
               DTSDZ(K)  = (STC(K) - TBOT) / ( 0.5*(ZSNSO(K-1)+ZSNSO(K)) - ZBOT)
               BOTFLX    = -DF(K) * DTSDZ(K)
           END IF
           EFLUX(K)  = (-BOTFLX - DF(K-1)*DTSDZ(K-1) ) - PHI(K)
        END IF
    END DO

    DO K = ISNOW+1, NSOIL
        IF (K == ISNOW+1) THEN
           AI(K)    =   0.0
           CI(K)    = - DF(K)   * DDZ(K) / DENOM(K)
           IF (OPT_STC == 1 .OR. OPT_STC == 3 ) THEN
              BI(K) = - CI(K)
           END IF                                        
           IF (OPT_STC == 2) THEN
              BI(K) = - CI(K) + DF(K)/(0.5*ZSNSO(K)*ZSNSO(K)*HCPCT(K))
           END IF
        ELSE IF (K < NSOIL) THEN
           AI(K)    = - DF(K-1) * DDZ(K-1) / DENOM(K) 
           CI(K)    = - DF(K  ) * DDZ(K  ) / DENOM(K) 
           BI(K)    = - (AI(K) + CI (K))
        ELSE IF (K == NSOIL) THEN
           AI(K)    = - DF(K-1) * DDZ(K-1) / DENOM(K) 
           CI(K)    = 0.0
           BI(K)    = - (AI(K) + CI(K))
        END IF
           RHSTS(K)  = EFLUX(K)/ (-DENOM(K))
    END DO

  END SUBROUTINE HRT

!== begin hstep ====================================================================================

  SUBROUTINE HSTEP (parameters,NSNOW     ,NSOIL     ,ISNOW     ,DT        ,  &
                    AI        ,BI        ,CI        ,RHSTS     ,  &
                    STC       )  
! ----------------------------------------------------------------------
! CALCULATE/UPDATE THE SOIL TEMPERATURE FIELD.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,                         INTENT(IN)    :: NSOIL
    INTEGER,                         INTENT(IN)    :: NSNOW
    INTEGER,                         INTENT(IN)    :: ISNOW
    REAL(r8),                            INTENT(IN)    :: DT

! output & input
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: RHSTS
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: AI
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: BI
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: CI
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC

! local
    INTEGER                                        :: K
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)                :: RHSTSIN
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL)                :: CIIN
! ----------------------------------------------------------------------

    DO K = ISNOW+1,NSOIL
       RHSTS(K) =   RHSTS(K) * DT
       AI(K)    =      AI(K) * DT
       BI(K)    = 1. + BI(K) * DT
       CI(K)    =      CI(K) * DT
    END DO

! copy values for input variables before call to rosr12

    DO K = ISNOW+1,NSOIL
       RHSTSIN(K) = RHSTS(K)
       CIIN(K)    = CI(K)
    END DO

! solve the tri-diagonal matrix equation

    CALL ROSR12 (CI,AI,BI,CIIN,RHSTSIN,RHSTS,ISNOW+1,NSOIL,NSNOW)

! update snow & soil temperature

    DO K = ISNOW+1,NSOIL
       STC (K) = STC (K) + CI (K)
    END DO

  END SUBROUTINE HSTEP

!== begin rosr12 ===================================================================================

  SUBROUTINE ROSR12 (P,A,B,C,D,DELTA,NTOP,NSOIL,NSNOW)
! ----------------------------------------------------------------------
! SUBROUTINE ROSR12
! ----------------------------------------------------------------------
! INVERT (SOLVE) THE TRI-DIAGONAL MATRIX PROBLEM SHOWN BELOW:
! ###                                            ### ###  ###   ###  ###
! #B(1), C(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! #A(2), B(2), C(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! # 0  , A(3), B(3), C(3),  0  ,   . . .  ,    0   # #      #   # D(3) #
! # 0  ,  0  , A(4), B(4), C(4),   . . .  ,    0   # # P(4) #   # D(4) #
! # 0  ,  0  ,  0  , A(5), B(5),   . . .  ,    0   # # P(5) #   # D(5) #
! # .                                          .   # #  .   # = #   .  #
! # .                                          .   # #  .   #   #   .  #
! # .                                          .   # #  .   #   #   .  #
! # 0  , . . . , 0 , A(M-2), B(M-2), C(M-2),   0   # #P(M-2)#   #D(M-2)#
! # 0  , . . . , 0 ,   0   , A(M-1), B(M-1), C(M-1)# #P(M-1)#   #D(M-1)#
! # 0  , . . . , 0 ,   0   ,   0   ,  A(M) ,  B(M) # # P(M) #   # D(M) #
! ###                                            ### ###  ###   ###  ###
! ----------------------------------------------------------------------
    IMPLICIT NONE

    INTEGER, INTENT(IN)   :: NTOP           
    INTEGER, INTENT(IN)   :: NSOIL,NSNOW
    INTEGER               :: K, KK

    REAL(r8), DIMENSION(-NSNOW+1:NSOIL),INTENT(IN):: A, B, D
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL),INTENT(INOUT):: C,P,DELTA

! ----------------------------------------------------------------------
! INITIALIZE EQN COEF C FOR THE LOWEST SOIL LAYER
! ----------------------------------------------------------------------
    C (NSOIL) = 0.0
    P (NTOP) = - C (NTOP) / B (NTOP)
! ----------------------------------------------------------------------
! SOLVE THE COEFS FOR THE 1ST SOIL LAYER
! ----------------------------------------------------------------------
    DELTA (NTOP) = D (NTOP) / B (NTOP)
! ----------------------------------------------------------------------
! SOLVE THE COEFS FOR SOIL LAYERS 2 THRU NSOIL
! ----------------------------------------------------------------------
    DO K = NTOP+1,NSOIL
       P (K) = - C (K) * ( 1.0 / (B (K) + A (K) * P (K -1)) )
       DELTA (K) = (D (K) - A (K)* DELTA (K -1))* (1.0/ (B (K) + A (K)&
            * P (K -1)))
    END DO
! ----------------------------------------------------------------------
! SET P TO DELTA FOR LOWEST SOIL LAYER
! ----------------------------------------------------------------------
    P (NSOIL) = DELTA (NSOIL)
! ----------------------------------------------------------------------
! ADJUST P FOR SOIL LAYERS 2 THRU NSOIL
! ----------------------------------------------------------------------
    DO K = NTOP+1,NSOIL
       KK = NSOIL - K + (NTOP-1) + 1
       P (KK) = P (KK) * P (KK +1) + DELTA (KK)
    END DO
! ----------------------------------------------------------------------
  END SUBROUTINE ROSR12

!== begin phasechange ==============================================================================

  SUBROUTINE PHASECHANGE (parameters,NSNOW   ,NSOIL   ,ISNOW   ,DT      ,FACT    , & !in
                          DZSNSO  ,HCPCT   ,IST     ,ILOC    ,JLOC    , & !in
                          STC     ,SNICE   ,SNLIQ   ,SNEQV   ,SNOWH   , & !inout
                          SMC     ,SH2O    ,                            & !inout
                          QMELT   ,IMELT   ,PONDING )                     !out
! ----------------------------------------------------------------------
! melting/freezing of snow water and soil water
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER, INTENT(IN)                             :: ILOC   !grid index
  INTEGER, INTENT(IN)                             :: JLOC   !grid index
  INTEGER, INTENT(IN)                             :: NSNOW  !maximum no. of snow layers [=3]
  INTEGER, INTENT(IN)                             :: NSOIL  !No. of soil layers [=4]
  INTEGER, INTENT(IN)                             :: ISNOW  !actual no. of snow layers [<=3]
  INTEGER, INTENT(IN)                             :: IST    !surface type: 1->soil; 2->lake
  REAL(r8), INTENT(IN)                                :: DT     !land model time step (sec)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: FACT   !temporary
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: DZSNSO !snow/soil layer thickness [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)     :: HCPCT  !heat capacity (J/m3/k)

! outputs
  INTEGER, DIMENSION(-NSNOW+1:NSOIL), INTENT(OUT) :: IMELT  !phase change index
  REAL(r8),                               INTENT(OUT) :: QMELT  !snowmelt rate [mm/s]
  REAL(r8),                               INTENT(OUT) :: PONDING!snowmelt when snow has no layer [mm]

! inputs and outputs

  REAL(r8), INTENT(INOUT) :: SNEQV
  REAL(r8), INTENT(INOUT) :: SNOWH
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT)  :: STC    !snow/soil layer temperature [k]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT)  :: SH2O   !soil liquid water [m3/m3]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT)  :: SMC    !total soil water [m3/m3]
  REAL(r8), DIMENSION(-NSNOW+1:0)    , INTENT(INOUT)  :: SNICE  !snow layer ice [mm]
  REAL(r8), DIMENSION(-NSNOW+1:0)    , INTENT(INOUT)  :: SNLIQ  !snow layer liquid water [mm]

! local

  INTEGER                         :: J         !do loop index
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: HM        !energy residual [w/m2]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: XM        !melting or freezing water [kg/m2]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: WMASS0
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: WICE0 
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: WLIQ0 
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: MICE      !soil/snow ice mass [mm]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: MLIQ      !soil/snow liquid water mass [mm]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL) :: SUPERCOOL !supercooled water in soil (kg/m2)
  REAL(r8)                            :: HEATR     !energy residual or loss after melting/freezing
  REAL(r8)                            :: TEMP1     !temporary variables [kg/m2]
  REAL(r8)                            :: PROPOR
  REAL(r8)                            :: SMP       !frozen water potential (mm)
  REAL(r8)                            :: XMF       !total latent heat of phase change

! ----------------------------------------------------------------------
! Initialization

    QMELT   = 0.
    PONDING = 0.
    XMF     = 0.

    DO J = -NSNOW+1, NSOIL
         SUPERCOOL(J) = 0.0
    END DO

    DO J = ISNOW+1,0       ! all layers
         MICE(J) = SNICE(J)
         MLIQ(J) = SNLIQ(J)
    END DO

    DO J = 1, NSOIL               ! soil
         MLIQ(J) =  SH2O(J)            * DZSNSO(J) * 1000.
         MICE(J) = (SMC(J) - SH2O(J))  * DZSNSO(J) * 1000.
    END DO

    DO J = ISNOW+1,NSOIL       ! all layers
         IMELT(J)    = 0
         HM(J)       = 0.
         XM(J)       = 0.
         WICE0(J)    = MICE(J)
         WLIQ0(J)    = MLIQ(J)
         WMASS0(J)   = MICE(J) + MLIQ(J)
    ENDDO

    if(ist == 1) then
      DO J = 1,NSOIL
         IF (OPT_FRZ == 1) THEN
            IF(STC(J) < TFRZ) THEN
               SMP = HFUS*(TFRZ-STC(J))/(GRAV*STC(J))             !(m)
               SUPERCOOL(J) = parameters%SMCMAX(J)*(SMP/parameters%PSISAT(J))**(-1./parameters%BEXP(J))
               SUPERCOOL(J) = SUPERCOOL(J)*DZSNSO(J)*1000.        !(mm)
            END IF
         END IF
         IF (OPT_FRZ == 2) THEN
               CALL FRH2O (parameters,J,SUPERCOOL(J),STC(J),SMC(J),SH2O(J))
               SUPERCOOL(J) = SUPERCOOL(J)*DZSNSO(J)*1000.        !(mm)
         END IF
      ENDDO
    end if

    DO J = ISNOW+1,NSOIL
         IF (MICE(J) > 0. .AND. STC(J) >= TFRZ) THEN  !melting 
             IMELT(J) = 1
         ENDIF
         IF (MLIQ(J) > SUPERCOOL(J) .AND. STC(J) < TFRZ) THEN
             IMELT(J) = 2
         ENDIF

         ! If snow exists, but its thickness is not enough to create a layer
         IF (ISNOW == 0 .AND. SNEQV > 0. .AND. J == 1) THEN
             IF (STC(J) >= TFRZ) THEN
                IMELT(J) = 1
             ENDIF
         ENDIF
    ENDDO

! Calculate the energy surplus and loss for melting and freezing

    DO J = ISNOW+1,NSOIL
         IF (IMELT(J) > 0) THEN
             HM(J) = (STC(J)-TFRZ)/FACT(J)
             STC(J) = TFRZ
         ENDIF

         IF (IMELT(J) == 1 .AND. HM(J) < 0.) THEN
            HM(J) = 0.
            IMELT(J) = 0
         ENDIF
         IF (IMELT(J) == 2 .AND. HM(J) > 0.) THEN
            HM(J) = 0.
            IMELT(J) = 0
         ENDIF
         XM(J) = HM(J)*DT/HFUS                           
    ENDDO

! The rate of melting and freezing for snow without a layer, needs more work.

    IF (ISNOW == 0 .AND. SNEQV > 0. .AND. XM(1) > 0.) THEN  
        TEMP1  = SNEQV
        SNEQV  = MAX(0.,TEMP1-XM(1))  
        PROPOR = SNEQV/TEMP1
        SNOWH  = MAX(0.,PROPOR * SNOWH)
        HEATR  = HM(1) - HFUS*(TEMP1-SNEQV)/DT  
        IF (HEATR > 0.) THEN
              XM(1) = HEATR*DT/HFUS             
              HM(1) = HEATR                    
        ELSE
              XM(1) = 0.
              HM(1) = 0.
        ENDIF
        QMELT   = MAX(0.,(TEMP1-SNEQV))/DT
        XMF     = HFUS*QMELT
        PONDING = TEMP1-SNEQV
    ENDIF

! The rate of melting and freezing for snow and soil

    DO J = ISNOW+1,NSOIL
      IF (IMELT(J) > 0 .AND. ABS(HM(J)) > 0.) THEN

         HEATR = 0.
         IF (XM(J) > 0.) THEN                            
            MICE(J) = MAX(0., WICE0(J)-XM(J))
            HEATR = HM(J) - HFUS*(WICE0(J)-MICE(J))/DT
         ELSE IF (XM(J) < 0.) THEN                      
            IF (J <= 0) THEN                             ! snow
               MICE(J) = MIN(WMASS0(J), WICE0(J)-XM(J))  
            ELSE                                         ! soil
               IF (WMASS0(J) < SUPERCOOL(J)) THEN
                  MICE(J) = 0.
               ELSE
                  MICE(J) = MIN(WMASS0(J) - SUPERCOOL(J),WICE0(J)-XM(J))
                  MICE(J) = MAX(MICE(J),0.0)
               ENDIF
            ENDIF
            HEATR = HM(J) - HFUS*(WICE0(J)-MICE(J))/DT
         ENDIF

         MLIQ(J) = MAX(0.,WMASS0(J)-MICE(J))

         IF (ABS(HEATR) > 0.) THEN
            STC(J) = STC(J) + FACT(J)*HEATR
            IF (J <= 0) THEN                             ! snow
               IF (MLIQ(J)*MICE(J)>0.) STC(J) = TFRZ
            END IF
         ENDIF

         XMF = XMF + HFUS * (WICE0(J)-MICE(J))/DT

         IF (J < 1) THEN
            QMELT = QMELT + MAX(0.,(WICE0(J)-MICE(J)))/DT
         ENDIF
      ENDIF
    ENDDO

    DO J = ISNOW+1,0             ! snow
       SNLIQ(J) = MLIQ(J)
       SNICE(J) = MICE(J)
    END DO

    DO J = 1, NSOIL              ! soil
       SH2O(J) =  MLIQ(J)            / (1000. * DZSNSO(J))
       SMC(J)  = (MLIQ(J) + MICE(J)) / (1000. * DZSNSO(J))
    END DO
   
  END SUBROUTINE PHASECHANGE

!== begin frh2o ====================================================================================

  SUBROUTINE FRH2O (parameters,ISOIL,FREE,TKELV,SMC,SH2O)

! ----------------------------------------------------------------------
! SUBROUTINE FRH2O
! ----------------------------------------------------------------------
! CALCULATE AMOUNT OF SUPERCOOLED LIQUID SOIL WATER CONTENT IF
! TEMPERATURE IS BELOW 273.15K (TFRZ).  REQUIRES NEWTON-TYPE ITERATION
! TO SOLVE THE NONLINEAR IMPLICIT EQUATION GIVEN IN EQN 17 OF KOREN ET AL
! (1999, JGR, VOL 104(D16), 19569-19585).
! ----------------------------------------------------------------------
! NEW VERSION (JUNE 2001): MUCH FASTER AND MORE ACCURATE NEWTON
! ITERATION ACHIEVED BY FIRST TAKING LOG OF EQN CITED ABOVE -- LESS THAN
! 4 (TYPICALLY 1 OR 2) ITERATIONS ACHIEVES CONVERGENCE.  ALSO, EXPLICIT
! 1-STEP SOLUTION OPTION FOR SPECIAL CASE OF PARAMETER CK=0, WHICH
! REDUCES THE ORIGINAL IMPLICIT EQUATION TO A SIMPLER EXPLICIT FORM,
! KNOWN AS THE "FLERCHINGER EQN". IMPROVED HANDLING OF SOLUTION IN THE
! LIMIT OF FREEZING POINT TEMPERATURE TFRZ.
! ----------------------------------------------------------------------
! INPUT:

!   TKELV.........TEMPERATURE (Kelvin)
!   SMC...........TOTAL SOIL MOISTURE CONTENT (VOLUMETRIC)
!   SH2O..........LIQUID SOIL MOISTURE CONTENT (VOLUMETRIC)
!   B.............SOIL TYPE "B" PARAMETER (FROM REDPRM)
!   PSISAT........SATURATED SOIL MATRIC POTENTIAL (FROM REDPRM)

! OUTPUT:
!   FREE..........SUPERCOOLED LIQUID WATER CONTENT [m3/m3]
! ----------------------------------------------------------------------
    IMPLICIT NONE
  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,INTENT(IN)   :: ISOIL
    REAL(r8), INTENT(IN)     :: SH2O,SMC,TKELV
    REAL(r8), INTENT(OUT)    :: FREE
    REAL(r8)                 :: BX,DENOM,DF,DSWL,FK,SWL,SWLK
    INTEGER              :: NLOG,KCOUNT
!      PARAMETER(CK = 0.0)
    REAL(r8), PARAMETER      :: CK = 8.0, BLIM = 5.5, ERROR = 0.005,       &
         DICE = 920.0
    CHARACTER(LEN=80)    :: message

! ----------------------------------------------------------------------
! LIMITS ON PARAMETER B: B < 5.5  (use parameter BLIM)
! SIMULATIONS SHOWED IF B > 5.5 UNFROZEN WATER CONTENT IS
! NON-REAL(r8)ISTICALLY HIGH AT VERY LOW TEMPERATURES.
! ----------------------------------------------------------------------
    BX = parameters%BEXP(ISOIL)
! ----------------------------------------------------------------------
! INITIALIZING ITERATIONS COUNTER AND ITERATIVE SOLUTION FLAG.
! ----------------------------------------------------------------------

    IF (parameters%BEXP(ISOIL) >  BLIM) BX = BLIM
    NLOG = 0

! ----------------------------------------------------------------------
!  IF TEMPERATURE NOT SIGNIFICANTLY BELOW FREEZING (TFRZ), SH2O = SMC
! ----------------------------------------------------------------------
    KCOUNT = 0
    IF (TKELV > (TFRZ- 1.E-3)) THEN
       FREE = SMC
    ELSE

! ----------------------------------------------------------------------
! OPTION 1: ITERATED SOLUTION IN KOREN ET AL, JGR, 1999, EQN 17
! ----------------------------------------------------------------------
! INITIAL GUESS FOR SWL (frozen content)
! ----------------------------------------------------------------------
       IF (CK /= 0.0) THEN
          SWL = SMC - SH2O
! ----------------------------------------------------------------------
! KEEP WITHIN BOUNDS.
! ----------------------------------------------------------------------
          IF (SWL > (SMC -0.02)) SWL = SMC -0.02
! ----------------------------------------------------------------------
!  START OF ITERATIONS
! ----------------------------------------------------------------------
          IF (SWL < 0.) SWL = 0.
1001      Continue
          IF (.NOT.( (NLOG < 10) .AND. (KCOUNT == 0)))   goto 1002
          NLOG = NLOG +1
!          DF = ALOG ( ( parameters%PSISAT(ISOIL) * GRAV / HFUS ) * ( ( 1. + CK * SWL )**2.) * &
!               ( parameters%SMCMAX(ISOIL) / (SMC - SWL) )** BX) - ALOG ( - (               &
!               TKELV - TFRZ)/ TKELV)
          DF = LOG ( ( parameters%PSISAT(ISOIL) * GRAV / HFUS ) * ( ( 1. + CK * SWL )**2.) * &
               ( parameters%SMCMAX(ISOIL) / (SMC - SWL) )** BX) - LOG ( - (               &
               TKELV - TFRZ)/ TKELV)
          DENOM = 2. * CK / ( 1. + CK * SWL ) + BX / ( SMC - SWL )
          SWLK = SWL - DF / DENOM
! ----------------------------------------------------------------------
! BOUNDS USEFUL FOR MATHEMATICAL SOLUTION.
! ----------------------------------------------------------------------
          IF (SWLK > (SMC -0.02)) SWLK = SMC - 0.02
          IF (SWLK < 0.) SWLK = 0.

! ----------------------------------------------------------------------
! MATHEMATICAL SOLUTION BOUNDS APPLIED.
! ----------------------------------------------------------------------
          DSWL = ABS (SWLK - SWL)
! IF MORE THAN 10 ITERATIONS, USE EXPLICIT METHOD (CK=0 APPROX.)
! WHEN DSWL LESS OR EQ. ERROR, NO MORE ITERATIONS REQUIRED.
! ----------------------------------------------------------------------
          SWL = SWLK
          IF ( DSWL <= ERROR ) THEN
             KCOUNT = KCOUNT +1
          END IF
! ----------------------------------------------------------------------
!  END OF ITERATIONS
! ----------------------------------------------------------------------
! BOUNDS APPLIED WITHIN DO-BLOCK ARE VALID FOR PHYSICAL SOLUTION.
! ----------------------------------------------------------------------
          goto 1001
1002      continue
          FREE = SMC - SWL
       END IF
! ----------------------------------------------------------------------
! END OPTION 1
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! OPTION 2: EXPLICIT SOLUTION FOR FLERCHINGER EQ. i.e. CK=0
! IN KOREN ET AL., JGR, 1999, EQN 17
! APPLY PHYSICAL BOUNDS TO FLERCHINGER SOLUTION
! ----------------------------------------------------------------------
       IF (KCOUNT == 0) THEN
          write(message, '("Flerchinger used in NEW version. Iterations=", I6)') NLOG
          print*, (trim(message))
          FK = ( ( (HFUS / (GRAV * ( - parameters%PSISAT(ISOIL))))*                    &
               ( (TKELV - TFRZ)/ TKELV))** ( -1/ BX))* parameters%SMCMAX(ISOIL)
          IF (FK < 0.02) FK = 0.02
          FREE = MIN (FK, SMC)
! ----------------------------------------------------------------------
! END OPTION 2
! ----------------------------------------------------------------------
       END IF
    END IF
! ----------------------------------------------------------------------
  END SUBROUTINE FRH2O
! ----------------------------------------------------------------------
! ==================================================================================================
! **********************End of energy subroutines***********************
! ==================================================================================================

!== begin water ====================================================================================

  SUBROUTINE WATER (parameters,VEGTYP ,NSNOW  ,NSOIL  ,IMELT  ,DT     ,UU     , & !in
                    VV     ,FCEV   ,FCTR   ,QPRECC ,QPRECL ,ELAI   , & !in
                    ESAI   ,SFCTMP ,QVAP   ,QDEW   ,ZSOIL  ,BTRANI , & !in
                    FICEOLD,PONDING,TG     ,IST    ,FVEG   ,ILOC   ,JLOC ,SMCEQ , & !in
                    BDFALL ,FP     ,RAIN   ,SNOW,                    & !in  MB/AN: v3.7
		    QSNOW  ,QRAIN  ,SNOWHIN,LATHEAV,LATHEAG,frozen_canopy,frozen_ground,    & !in  MB
                    ISNOW  ,CANLIQ ,CANICE ,TV     ,SNOWH  ,SNEQV  , & !inout
                    SNICE  ,SNLIQ  ,STC    ,ZSNSO  ,SH2O   ,SMC    , & !inout
                    SICE   ,ZWT    ,WA     ,WT     ,DZSNSO ,WSLAKE , & !inout
                    SMCWTD ,DEEPRECH,RECH                          , & !inout
                    CMC    ,ECAN   ,ETRAN  ,FWET   ,RUNSRF ,RUNSUB , & !out
                    QIN    ,QDIS   ,PONDING1       ,PONDING2,        &
                    QSNBOT                                           &
#ifdef WRF_HYDRO
                        ,sfcheadrt                     &
#endif
                    )  !out
! ----------------------------------------------------------------------  
! Code history:
! Initial code: Guo-Yue Niu, Oct. 2007
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                         INTENT(IN)    :: ILOC    !grid index
  INTEGER,                         INTENT(IN)    :: JLOC    !grid index
  INTEGER,                         INTENT(IN)    :: VEGTYP  !vegetation type
  INTEGER,                         INTENT(IN)    :: NSNOW   !maximum no. of snow layers
  INTEGER                        , INTENT(IN)    :: IST     !surface type 1-soil; 2-lake
  INTEGER,                         INTENT(IN)    :: NSOIL   !no. of soil layers
  INTEGER, DIMENSION(-NSNOW+1:0) , INTENT(IN)    :: IMELT   !melting state index [1-melt; 2-freeze]
  REAL(r8),                            INTENT(IN)    :: DT      !main time step (s)
  REAL(r8),                            INTENT(IN)    :: UU      !u-direction wind speed [m/s]
  REAL(r8),                            INTENT(IN)    :: VV      !v-direction wind speed [m/s]
  REAL(r8),                            INTENT(IN)    :: FCEV    !canopy evaporation (w/m2) [+ to atm ]
  REAL(r8),                            INTENT(IN)    :: FCTR    !transpiration (w/m2) [+ to atm]
  REAL(r8),                            INTENT(IN)    :: QPRECC  !convective precipitation (mm/s)
  REAL(r8),                            INTENT(IN)    :: QPRECL  !large-scale precipitation (mm/s)
  REAL(r8),                            INTENT(IN)    :: ELAI    !leaf area index, after burying by snow
  REAL(r8),                            INTENT(IN)    :: ESAI    !stem area index, after burying by snow
  REAL(r8),                            INTENT(IN)    :: SFCTMP  !surface air temperature [k]
  REAL(r8),                            INTENT(IN)    :: QVAP    !soil surface evaporation rate[mm/s]
  REAL(r8),                            INTENT(IN)    :: QDEW    !soil surface dew rate[mm/s]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL   !depth of layer-bottom from soil surface
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)    :: BTRANI  !soil water stress factor (0 to 1)
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: FICEOLD !ice fraction at last timestep
!  REAL(r8)                           , INTENT(IN)    :: PONDING ![mm]
  REAL(r8)                           , INTENT(IN)    :: TG      !ground temperature (k)
  REAL(r8)                           , INTENT(IN)    :: FVEG    !greeness vegetation fraction (-)
  REAL(r8)                           , INTENT(IN)    :: BDFALL   !bulk density of snowfall (kg/m3) ! MB/AN: v3.7
  REAL(r8)                           , INTENT(IN)    :: FP       !fraction of the gridcell that receives precipitation ! MB/AN: v3.7
  REAL(r8)                           , INTENT(IN)    :: RAIN     !rainfall (mm/s) ! MB/AN: v3.7
  REAL(r8)                           , INTENT(IN)    :: SNOW     !snowfall (mm/s) ! MB/AN: v3.7
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)    :: SMCEQ   !equilibrium soil water content [m3/m3] (used in m-m&f groundwater dynamics)
  REAL(r8)                           , INTENT(IN)    :: QSNOW   !snow at ground srf (mm/s) [+]
  REAL(r8)                           , INTENT(IN)    :: QRAIN   !rain at ground srf (mm) [+]
  REAL(r8)                           , INTENT(IN)    :: SNOWHIN !snow depth increasing rate (m/s)

! input/output
  INTEGER,                         INTENT(INOUT) :: ISNOW   !actual no. of snow layers
  REAL(r8),                            INTENT(INOUT) :: CANLIQ  !intercepted liquid water (mm)
  REAL(r8),                            INTENT(INOUT) :: CANICE  !intercepted ice mass (mm)
  REAL(r8),                            INTENT(INOUT) :: TV      !vegetation temperature (k)
  REAL(r8),                            INTENT(INOUT) :: SNOWH   !snow height [m]
  REAL(r8),                            INTENT(INOUT) :: SNEQV   !snow water eqv. [mm]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE   !snow layer ice [mm]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ   !snow layer liquid water [mm]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC     !snow/soil layer temperature [k]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO   !depth of snow/soil layer-bottom
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO  !snow/soil layer thickness [m]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O    !soil liquid water content [m3/m3]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE    !soil ice content [m3/m3]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SMC     !total soil water content [m3/m3]
  REAL(r8),                            INTENT(INOUT) :: ZWT     !the depth to water table [m]
  REAL(r8),                            INTENT(INOUT) :: WA      !water storage in aquifer [mm]
  REAL(r8),                            INTENT(INOUT) :: WT      !water storage in aquifer 
                                                            !+ stuarated soil [mm]
  REAL(r8),                            INTENT(INOUT) :: WSLAKE  !water storage in lake (can be -) (mm)
  REAL(r8)                           , INTENT(INOUT) :: PONDING ![mm]
  REAL(r8),                            INTENT(INOUT) :: SMCWTD !soil water content between bottom of the soil and water table [m3/m3]
  REAL(r8),                            INTENT(INOUT) :: DEEPRECH !recharge to or from the water table when deep [m]
  REAL(r8),                            INTENT(INOUT) :: RECH !recharge to or from the water table when shallow [m] (diagnostic)

! output
  REAL(r8),                            INTENT(OUT)   :: CMC     !intercepted water per ground area (mm)
  REAL(r8),                            INTENT(OUT)   :: ECAN    !evap of intercepted water (mm/s) [+]
  REAL(r8),                            INTENT(OUT)   :: ETRAN   !transpiration rate (mm/s) [+]
  REAL(r8),                            INTENT(OUT)   :: FWET    !wetted/snowed fraction of canopy (-)
  REAL(r8),                            INTENT(OUT)   :: RUNSRF  !surface runoff [mm/s] 
  REAL(r8),                            INTENT(OUT)   :: RUNSUB  !baseflow (sturation excess) [mm/s]
  REAL(r8),                            INTENT(OUT)   :: QIN     !groundwater recharge [mm/s]
  REAL(r8),                            INTENT(OUT)   :: QDIS    !groundwater discharge [mm/s]
  REAL(r8),                            INTENT(OUT)   :: PONDING1
  REAL(r8),                            INTENT(OUT)   :: PONDING2
  REAL(r8),                            INTENT(OUT)   :: QSNBOT  !melting water out of snow bottom [mm/s]
  REAL(r8)                              , INTENT(IN)   :: LATHEAV !latent heat vap./sublimation (j/kg)
  REAL(r8)                              , INTENT(IN)   :: LATHEAG !latent heat vap./sublimation (j/kg)
  LOGICAL                           , INTENT(IN)   :: FROZEN_GROUND ! used to define latent heat pathway
  LOGICAL                           , INTENT(IN)   :: FROZEN_CANOPY ! used to define latent heat pathway


! local
  INTEGER                                        :: IZ
  REAL(r8)                                           :: QINSUR  !water input on soil surface [m/s]
  REAL(r8)                                           :: QSEVA   !soil surface evap rate [mm/s]
  REAL(r8)                                           :: QSDEW   !soil surface dew rate [mm/s]
  REAL(r8)                                           :: QSNFRO  !snow surface frost rate[mm/s]
  REAL(r8)                                           :: QSNSUB  !snow surface sublimation rate [mm/s]
  REAL(r8), DIMENSION(       1:NSOIL)                :: ETRANI  !transpiration rate (mm/s) [+]
  REAL(r8), DIMENSION(       1:NSOIL)                :: WCND   !hydraulic conductivity (m/s)
  REAL(r8)                                           :: QDRAIN  !soil-bottom free drainage [mm/s] 
  REAL(r8)                                           :: SNOFLOW !glacier flow [mm/s]
  REAL(r8)                                           :: FCRMAX !maximum of FCR (-)

  REAL(r8), PARAMETER ::  WSLMAX = 5000.      !maximum lake water storage (mm)

#ifdef WRF_HYDRO
  REAL(r8)                           , INTENT(INOUT)    :: sfcheadrt
#endif

! ----------------------------------------------------------------------
! initialize

   ETRANI(1:NSOIL) = 0.
   SNOFLOW         = 0.
   RUNSUB          = 0.
   QINSUR          = 0.

! canopy-intercepted snowfall/rainfall, drips, and throughfall

   CALL CANWATER (parameters,VEGTYP ,DT     , & !in
                  FCEV   ,FCTR   ,ELAI   , & !in
                  ESAI   ,TG     ,FVEG   ,ILOC   , JLOC, & !in
                  BDFALL ,FROZEN_CANOPY  , & !in     
                  CANLIQ ,CANICE ,TV     ,                 & !inout
                  CMC    ,ECAN   ,ETRAN  , & !out
                  FWET      )                           !out

! sublimation, frost, evaporation, and dew

     QSNSUB = 0.
     IF (SNEQV > 0.) THEN
       QSNSUB = MIN(QVAP, SNEQV/DT)
     ENDIF
     QSEVA = QVAP-QSNSUB

     QSNFRO = 0.
     IF (SNEQV > 0.) THEN
        QSNFRO = QDEW
     ENDIF
     QSDEW = QDEW - QSNFRO

     CALL SNOWWATER (parameters,NSNOW  ,NSOIL  ,IMELT  ,DT     ,ZSOIL  , & !in
          &          SFCTMP ,SNOWHIN,QSNOW  ,QSNFRO ,QSNSUB , & !in
          &          QRAIN  ,FICEOLD,ILOC   ,JLOC   ,         & !in
          &          ISNOW  ,SNOWH  ,SNEQV  ,SNICE  ,SNLIQ  , & !inout
          &          SH2O   ,SICE   ,STC    ,ZSNSO  ,DZSNSO , & !inout
          &          QSNBOT ,SNOFLOW,PONDING1       ,PONDING2)  !out

   IF(FROZEN_GROUND) THEN
      SICE(1) =  SICE(1) + (QSDEW-QSEVA)*DT/(DZSNSO(1)*1000.)
      QSDEW = 0.0
      QSEVA = 0.0
      IF(SICE(1) < 0.) THEN
         SH2O(1) = SH2O(1) + SICE(1)
         SICE(1) = 0.
      END IF
   END IF

! convert units (mm/s -> m/s)

    !PONDING: melting water from snow when there is no layer
    QINSUR = (PONDING+PONDING1+PONDING2)/DT * 0.001
!    QINSUR = PONDING/DT * 0.001

    IF(ISNOW == 0) THEN
       QINSUR = QINSUR+(QSNBOT + QSDEW + QRAIN) * 0.001
    ELSE
       QINSUR = QINSUR+(QSNBOT + QSDEW) * 0.001
    ENDIF

    QSEVA  = QSEVA * 0.001 

    DO IZ = 1, parameters%NROOT
       ETRANI(IZ) = ETRAN * BTRANI(IZ) * 0.001
    ENDDO

#ifdef WRF_HYDRO
       QINSUR = QINSUR+sfcheadrt/DT*0.001  !sfcheadrt units (m)
#endif

! lake/soil water balances

    IF (IST == 2) THEN                                        ! lake
       RUNSRF = 0.
       IF(WSLAKE >= WSLMAX) RUNSRF = QINSUR*1000.             !mm/s
       WSLAKE = WSLAKE + (QINSUR-QSEVA)*1000.*DT -RUNSRF*DT   !mm
    ELSE                                                      ! soil
       CALL      SOILWATER (parameters,NSOIL  ,NSNOW  ,DT     ,ZSOIL  ,DZSNSO , & !in
                            QINSUR ,QSEVA  ,ETRANI ,SICE   ,ILOC   , JLOC , & !in
                            SH2O   ,SMC    ,ZWT    ,VEGTYP , & !inout
                           SMCWTD, DEEPRECH                       , & !inout
                            RUNSRF ,QDRAIN ,RUNSUB ,WCND   ,FCRMAX )   !out
 
       IF(OPT_RUN == 1) THEN 
          CALL GROUNDWATER (parameters,NSNOW  ,NSOIL  ,DT     ,SICE   ,ZSOIL  , & !in
                            STC    ,WCND   ,FCRMAX ,ILOC   ,JLOC   , & !in
                            SH2O   ,ZWT    ,WA     ,WT     ,         & !inout
                            QIN    ,QDIS   )                           !out
          RUNSUB       = QDIS          !mm/s
       END IF

       IF(OPT_RUN == 3 .or. OPT_RUN == 4) THEN 
          RUNSUB       = RUNSUB + QDRAIN        !mm/s
       END IF

       DO IZ = 1,NSOIL
           SMC(IZ) = SH2O(IZ) + SICE(IZ)
       ENDDO
 
       IF(OPT_RUN == 5) THEN
          CALL SHALLOWWATERTABLE (parameters,NSNOW  ,NSOIL, ZSOIL, DT       , & !in
                         DZSNSO ,SMCEQ   ,ILOC , JLOC        , & !in
                         SMC    ,ZWT    ,SMCWTD ,RECH, QDRAIN  ) !inout

          SH2O(NSOIL) = SMC(NSOIL) - SICE(NSOIL)
          RUNSUB = RUNSUB + QDRAIN !it really comes from subroutine watertable, which is not called with the same frequency as the soil routines here
          WA = 0.
       ENDIF

    ENDIF

    RUNSUB       = RUNSUB + SNOFLOW         !mm/s

  END SUBROUTINE WATER

!== begin canwater =================================================================================

  SUBROUTINE CANWATER (parameters,VEGTYP ,DT     , & !in
                       FCEV   ,FCTR   ,ELAI   , & !in
                       ESAI   ,TG     ,FVEG   ,ILOC   , JLOC , & !in
                       BDFALL ,FROZEN_CANOPY  ,  & !in      
                       CANLIQ ,CANICE ,TV     ,                 & !inout
                       CMC    ,ECAN   ,ETRAN  , & !out
                       FWET      )                           !out

! ------------------------ code history ------------------------------
! canopy hydrology
! --------------------------------------------------------------------
  IMPLICIT NONE
! ------------------------ input/output variables --------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,INTENT(IN)  :: ILOC    !grid index
  INTEGER,INTENT(IN)  :: JLOC    !grid index
  INTEGER,INTENT(IN)  :: VEGTYP  !vegetation type
  REAL(r8),   INTENT(IN)  :: DT      !main time step (s)
  REAL(r8),   INTENT(IN)  :: FCEV    !canopy evaporation (w/m2) [+ = to atm]
  REAL(r8),   INTENT(IN)  :: FCTR    !transpiration (w/m2) [+ = to atm]
  REAL(r8),   INTENT(IN)  :: ELAI    !leaf area index, after burying by snow
  REAL(r8),   INTENT(IN)  :: ESAI    !stem area index, after burying by snow
  REAL(r8),   INTENT(IN)  :: TG      !ground temperature (k)
  REAL(r8),   INTENT(IN)  :: FVEG    !greeness vegetation fraction (-)
  LOGICAL                           , INTENT(IN)   :: FROZEN_CANOPY ! used to define latent heat pathway
  REAL(r8)                           , INTENT(IN)    :: BDFALL   !bulk density of snowfall (kg/m3) ! MB/AN: v3.7

! input & output
  REAL(r8), INTENT(INOUT) :: CANLIQ  !intercepted liquid water (mm)
  REAL(r8), INTENT(INOUT) :: CANICE  !intercepted ice mass (mm)
  REAL(r8), INTENT(INOUT) :: TV      !vegetation temperature (k)

! output
  REAL(r8), INTENT(OUT)   :: CMC     !intercepted water (mm)
  REAL(r8), INTENT(OUT)   :: ECAN    !evaporation of intercepted water (mm/s) [+]
  REAL(r8), INTENT(OUT)   :: ETRAN   !transpiration rate (mm/s) [+]
  REAL(r8), INTENT(OUT)   :: FWET    !wetted or snowed fraction of the canopy (-)
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  REAL(r8)                :: MAXSNO  !canopy capacity for snow interception (mm)
  REAL(r8)                :: MAXLIQ  !canopy capacity for rain interception (mm)
  REAL(r8)                :: QEVAC   !evaporation rate (mm/s)
  REAL(r8)                :: QDEWC   !dew rate (mm/s)
  REAL(r8)                :: QFROC   !frost rate (mm/s)
  REAL(r8)                :: QSUBC   !sublimation rate (mm/s)
  REAL(r8)                :: QMELTC  !melting rate of canopy snow (mm/s)
  REAL(r8)                :: QFRZC   !refreezing rate of canopy liquid water (mm/s)
  REAL(r8)                :: CANMAS  !total canopy mass (kg/m2)
! --------------------------------------------------------------------
! initialization

      ECAN    = 0.0

! --------------------------- liquid water ------------------------------
! maximum canopy water

      MAXLIQ =  parameters%CH2OP * (ELAI+ ESAI)

! evaporation, transpiration, and dew

      IF (.NOT.FROZEN_CANOPY) THEN             ! Barlage: change to frozen_canopy
        ETRAN = MAX( FCTR/HVAP, 0. )
        QEVAC = MAX( FCEV/HVAP, 0. )
        QDEWC = ABS( MIN( FCEV/HVAP, 0. ) )
        QSUBC = 0.
        QFROC = 0.
      ELSE
        ETRAN = MAX( FCTR/HSUB, 0. )
        QEVAC = 0.
        QDEWC = 0.
        QSUBC = MAX( FCEV/HSUB, 0. )
        QFROC = ABS( MIN( FCEV/HSUB, 0. ) )
      ENDIF

! canopy water balance. for convenience allow dew to bring CANLIQ above
! maxh2o or else would have to re-adjust drip

       QEVAC = MIN(CANLIQ/DT,QEVAC)
       CANLIQ=MAX(0.,CANLIQ+(QDEWC-QEVAC)*DT)
       IF(CANLIQ <= 1.E-06) CANLIQ = 0.0

! --------------------------- canopy ice ------------------------------
! for canopy ice

      MAXSNO = 6.6*(0.27+46./BDFALL) * (ELAI+ ESAI)

      QSUBC = MIN(CANICE/DT,QSUBC) 
      CANICE= MAX(0.,CANICE + (QFROC-QSUBC)*DT)
      IF(CANICE.LE.1.E-6) CANICE = 0.
     
! wetted fraction of canopy

      IF(CANICE.GT.0.) THEN
           FWET = MAX(0.,CANICE) / MAX(MAXSNO,1.E-06)
      ELSE
           FWET = MAX(0.,CANLIQ) / MAX(MAXLIQ,1.E-06)
      ENDIF
      FWET = MIN(FWET, 1.) ** 0.667

! phase change

      QMELTC = 0.
      QFRZC = 0.

      IF(CANICE.GT.1.E-6.AND.TV.GT.TFRZ) THEN
         QMELTC = MIN(CANICE/DT,(TV-TFRZ)*CICE*CANICE/DENICE/(DT*HFUS))
         CANICE = MAX(0.,CANICE - QMELTC*DT)
         CANLIQ = MAX(0.,CANLIQ + QMELTC*DT)
         TV     = FWET*TFRZ + (1.-FWET)*TV
      ENDIF

      IF(CANLIQ.GT.1.E-6.AND.TV.LT.TFRZ) THEN
         QFRZC  = MIN(CANLIQ/DT,(TFRZ-TV)*CWAT*CANLIQ/DENH2O/(DT*HFUS))
         CANLIQ = MAX(0.,CANLIQ - QFRZC*DT)
         CANICE = MAX(0.,CANICE + QFRZC*DT)
         TV     = FWET*TFRZ + (1.-FWET)*TV
      ENDIF

! total canopy water

      CMC = CANLIQ + CANICE

! total canopy evaporation

      ECAN = QEVAC + QSUBC - QDEWC - QFROC

  END SUBROUTINE CANWATER

!== begin snowwater ================================================================================

  SUBROUTINE SNOWWATER (parameters,NSNOW  ,NSOIL  ,IMELT  ,DT     ,ZSOIL  , & !in
                        SFCTMP ,SNOWHIN,QSNOW  ,QSNFRO ,QSNSUB , & !in
                        QRAIN  ,FICEOLD,ILOC   ,JLOC   ,         & !in
                        ISNOW  ,SNOWH  ,SNEQV  ,SNICE  ,SNLIQ  , & !inout
                        SH2O   ,SICE   ,STC    ,ZSNSO  ,DZSNSO , & !inout
                        QSNBOT ,SNOFLOW,PONDING1       ,PONDING2)  !out
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                         INTENT(IN)    :: ILOC   !grid index
  INTEGER,                         INTENT(IN)    :: JLOC   !grid index
  INTEGER,                         INTENT(IN)    :: NSNOW  !maximum no. of snow layers
  INTEGER,                         INTENT(IN)    :: NSOIL  !no. of soil layers
  INTEGER, DIMENSION(-NSNOW+1:0) , INTENT(IN)    :: IMELT  !melting state index [0-no melt;1-melt]
  REAL(r8),                            INTENT(IN)    :: DT     !time step (s)
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL  !depth of layer-bottom from soil surface
  REAL(r8),                            INTENT(IN)    :: SFCTMP !surface air temperature [k]
  REAL(r8),                            INTENT(IN)    :: SNOWHIN!snow depth increasing rate (m/s)
  REAL(r8),                            INTENT(IN)    :: QSNOW  !snow at ground srf (mm/s) [+]
  REAL(r8),                            INTENT(IN)    :: QSNFRO !snow surface frost rate[mm/s]
  REAL(r8),                            INTENT(IN)    :: QSNSUB !snow surface sublimation rate[mm/s]
  REAL(r8),                            INTENT(IN)    :: QRAIN  !snow surface rain rate[mm/s]
  REAL(r8), DIMENSION(-NSNOW+1:0)    , INTENT(IN)    :: FICEOLD!ice fraction at last timestep

! input & output
  INTEGER,                         INTENT(INOUT) :: ISNOW  !actual no. of snow layers
  REAL(r8),                            INTENT(INOUT) :: SNOWH  !snow height [m]
  REAL(r8),                            INTENT(INOUT) :: SNEQV  !snow water eqv. [mm]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  !snow layer ice [mm]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  !snow layer liquid water [mm]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   !soil liquid moisture (m3/m3)
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE   !soil ice moisture (m3/m3)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    !snow layer temperature [k]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO  !depth of snow/soil layer-bottom
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO !snow/soil layer thickness [m]

! output
  REAL(r8),                              INTENT(OUT) :: QSNBOT !melting water out of snow bottom [mm/s]
  REAL(r8),                              INTENT(OUT) :: SNOFLOW!glacier flow [mm]
  REAL(r8),                              INTENT(OUT) :: PONDING1
  REAL(r8),                              INTENT(OUT) :: PONDING2

! local
  INTEGER :: IZ,i
  REAL(r8)    :: BDSNOW  !bulk density of snow (kg/m3)
! ----------------------------------------------------------------------
   SNOFLOW = 0.0
   PONDING1 = 0.0
   PONDING2 = 0.0

   CALL SNOWFALL (parameters,NSOIL  ,NSNOW  ,DT     ,QSNOW  ,SNOWHIN, & !in
                  SFCTMP ,ILOC   ,JLOC   ,                 & !in
                  ISNOW  ,SNOWH  ,DZSNSO ,STC    ,SNICE  , & !inout
                  SNLIQ  ,SNEQV  )                           !inout

! MB: do each if block separately

   IF(ISNOW < 0) &        ! when multi-layer
   CALL  COMPACT (parameters,NSNOW  ,NSOIL  ,DT     ,STC    ,SNICE  , & !in
                  SNLIQ  ,ZSOIL  ,IMELT  ,FICEOLD,ILOC   , JLOC ,& !in
                  ISNOW  ,DZSNSO ,ZSNSO  )                   !inout

   IF(ISNOW < 0) &        !when multi-layer
   CALL  COMBINE (parameters,NSNOW  ,NSOIL  ,ILOC   ,JLOC   ,         & !in
                  ISNOW  ,SH2O   ,STC    ,SNICE  ,SNLIQ  , & !inout
                  DZSNSO ,SICE   ,SNOWH  ,SNEQV  ,         & !inout
                  PONDING1       ,PONDING2)                  !out

   IF(ISNOW < 0) &        !when multi-layer
   CALL   DIVIDE (parameters,NSNOW  ,NSOIL  ,                         & !in
                  ISNOW  ,STC    ,SNICE  ,SNLIQ  ,DZSNSO )   !inout

   CALL  SNOWH2O (parameters,NSNOW  ,NSOIL  ,DT     ,QSNFRO ,QSNSUB , & !in 
                  QRAIN  ,ILOC   ,JLOC   ,                 & !in
                  ISNOW  ,DZSNSO ,SNOWH  ,SNEQV  ,SNICE  , & !inout
                  SNLIQ  ,SH2O   ,SICE   ,STC    ,         & !inout
                  QSNBOT ,PONDING1       ,PONDING2)           !out

!set empty snow layers to zero

   do iz = -nsnow+1, isnow
        snice(iz) = 0.
        snliq(iz) = 0.
        stc(iz)   = 0.
        dzsnso(iz)= 0.
        zsnso(iz) = 0.
   enddo

!to obtain equilibrium state of snow in glacier region
       
   IF(SNEQV > 2000.) THEN   ! 2000 mm -> maximum water depth
      BDSNOW      = SNICE(0) / DZSNSO(0)
      SNOFLOW     = (SNEQV - 2000.)
      SNICE(0)    = SNICE(0)  - SNOFLOW 
      DZSNSO(0)   = DZSNSO(0) - SNOFLOW/BDSNOW
      SNOFLOW     = SNOFLOW / DT
   END IF

! sum up snow mass for layered snow

   IF(ISNOW < 0) THEN  ! MB: only do for multi-layer
       SNEQV = 0.
       DO IZ = ISNOW+1,0
             SNEQV = SNEQV + SNICE(IZ) + SNLIQ(IZ)
       ENDDO
   END IF

! Reset ZSNSO and layer thinkness DZSNSO

   DO IZ = ISNOW+1, 0
        DZSNSO(IZ) = -DZSNSO(IZ)
   END DO

   DZSNSO(1) = ZSOIL(1)
   DO IZ = 2,NSOIL
        DZSNSO(IZ) = (ZSOIL(IZ) - ZSOIL(IZ-1))
   END DO

   ZSNSO(ISNOW+1) = DZSNSO(ISNOW+1)
   DO IZ = ISNOW+2 ,NSOIL
       ZSNSO(IZ) = ZSNSO(IZ-1) + DZSNSO(IZ)
   ENDDO

   DO IZ = ISNOW+1 ,NSOIL
       DZSNSO(IZ) = -DZSNSO(IZ)
   END DO

  END SUBROUTINE SNOWWATER

!== begin snowfall =================================================================================

  SUBROUTINE SNOWFALL (parameters,NSOIL  ,NSNOW  ,DT     ,QSNOW  ,SNOWHIN , & !in
                       SFCTMP ,ILOC   ,JLOC   ,                  & !in
                       ISNOW  ,SNOWH  ,DZSNSO ,STC    ,SNICE   , & !inout
                       SNLIQ  ,SNEQV  )                            !inout
! ----------------------------------------------------------------------
! snow depth and density to account for the new snowfall.
! new values of snow depth & density returned.
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                            INTENT(IN) :: ILOC   !grid index
  INTEGER,                            INTENT(IN) :: JLOC   !grid index
  INTEGER,                            INTENT(IN) :: NSOIL  !no. of soil layers
  INTEGER,                            INTENT(IN) :: NSNOW  !maximum no. of snow layers
  REAL(r8),                               INTENT(IN) :: DT     !main time step (s)
  REAL(r8),                               INTENT(IN) :: QSNOW  !snow at ground srf (mm/s) [+]
  REAL(r8),                               INTENT(IN) :: SNOWHIN!snow depth increasing rate (m/s)
  REAL(r8),                               INTENT(IN) :: SFCTMP !surface air temperature [k]

! input and output

  INTEGER,                         INTENT(INOUT) :: ISNOW  !actual no. of snow layers
  REAL(r8),                            INTENT(INOUT) :: SNOWH  !snow depth [m]
  REAL(r8),                            INTENT(INOUT) :: SNEQV  !swow water equivalent [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO !thickness of snow/soil layers (m)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    !snow layer temperature [k]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE  !snow layer ice [mm]
  REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ  !snow layer liquid water [mm]

! local

  INTEGER :: NEWNODE            ! 0-no new layers, 1-creating new layers
! ----------------------------------------------------------------------
    NEWNODE  = 0

! shallow snow / no layer

    IF(ISNOW == 0 .and. QSNOW > 0.)  THEN
      SNOWH = SNOWH + SNOWHIN * DT
      SNEQV = SNEQV + QSNOW * DT
    END IF

! creating a new layer
 
    IF(ISNOW == 0  .AND. QSNOW>0. .AND. SNOWH >= 0.025) THEN !MB: change limit
!    IF(ISNOW == 0  .AND. QSNOW>0. .AND. SNOWH >= 0.05) THEN
      ISNOW    = -1
      NEWNODE  =  1
      DZSNSO(0)= SNOWH
      SNOWH    = 0.
      STC(0)   = MIN(273.16, SFCTMP)   ! temporary setup
      SNICE(0) = SNEQV
      SNLIQ(0) = 0.
    END IF

! snow with layers

    IF(ISNOW <  0 .AND. NEWNODE == 0 .AND. QSNOW > 0.) then
         SNICE(ISNOW+1)  = SNICE(ISNOW+1)   + QSNOW   * DT
         DZSNSO(ISNOW+1) = DZSNSO(ISNOW+1)  + SNOWHIN * DT
    ENDIF

! ----------------------------------------------------------------------
  END SUBROUTINE SNOWFALL

!== begin combine ==================================================================================

  SUBROUTINE COMBINE (parameters,NSNOW  ,NSOIL  ,ILOC   ,JLOC   ,         & !in
                      ISNOW  ,SH2O   ,STC    ,SNICE  ,SNLIQ  , & !inout
                      DZSNSO ,SICE   ,SNOWH  ,SNEQV  ,         & !inout
                      PONDING1       ,PONDING2)                  !out
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER, INTENT(IN)     :: ILOC
    INTEGER, INTENT(IN)     :: JLOC
    INTEGER, INTENT(IN)     :: NSNOW                        !maximum no. of snow layers
    INTEGER, INTENT(IN)     :: NSOIL                        !no. of soil layers

! input and output

    INTEGER,                         INTENT(INOUT) :: ISNOW !actual no. of snow layers
    REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O  !soil liquid moisture (m3/m3)
    REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE  !soil ice moisture (m3/m3)
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC   !snow layer temperature [k]
    REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE !snow layer ice [mm]
    REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ !snow layer liquid water [mm]
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO!snow layer depth [m]
    REAL(r8),                            INTENT(INOUT) :: sneqv !snow water equivalent [m]
    REAL(r8),                            INTENT(INOUT) :: snowh !snow depth [m]
    REAL(r8),                            INTENT(OUT) :: PONDING1
    REAL(r8),                            INTENT(OUT) :: PONDING2

! local variables:

    INTEGER :: I,J,K,L               ! node indices
    INTEGER :: ISNOW_OLD             ! number of top snow layer
    INTEGER :: MSSI                  ! node index
    INTEGER :: NEIBOR                ! adjacent node selected for combination
    REAL(r8)    :: ZWICE                 ! total ice mass in snow
    REAL(r8)    :: ZWLIQ                 ! total liquid water in snow

    REAL(r8)    :: DZMIN(3)              ! minimum of top snow layer
!    DATA DZMIN /0.045, 0.05, 0.2/
    DATA DZMIN /0.025, 0.025, 0.1/  ! MB: change limit
!-----------------------------------------------------------------------

       ISNOW_OLD = ISNOW

       DO J = ISNOW_OLD+1,0
          IF (SNICE(J) <= .1) THEN
             IF(J /= 0) THEN
                SNLIQ(J+1) = SNLIQ(J+1) + SNLIQ(J)
                SNICE(J+1) = SNICE(J+1) + SNICE(J)
             ELSE
               IF (ISNOW_OLD < -1) THEN    ! MB/KM: change to ISNOW
                SNLIQ(J-1) = SNLIQ(J-1) + SNLIQ(J)
                SNICE(J-1) = SNICE(J-1) + SNICE(J)
               ELSE
	         IF(SNICE(J) >= 0.) THEN
                  PONDING1 = SNLIQ(J)    ! ISNOW WILL GET SET TO ZERO BELOW; PONDING1 WILL GET 
                  SNEQV = SNICE(J)       ! ADDED TO PONDING FROM PHASECHANGE PONDING SHOULD BE
                  SNOWH = DZSNSO(J)      ! ZERO HERE BECAUSE IT WAS CALCULATED FOR THIN SNOW
		 ELSE   ! SNICE OVER-SUBLIMATED EARLIER
		  PONDING1 = SNLIQ(J) + SNICE(J)
		  IF(PONDING1 < 0.) THEN  ! IF SNICE AND SNLIQ SUBLIMATES REMOVE FROM SOIL
		   SICE(1) = MAX(0.0,SICE(1)+PONDING1/(DZSNSO(1)*1000.))
                   PONDING1 = 0.0
		  END IF
                  SNEQV = 0.0
                  SNOWH = 0.0
		 END IF
                 SNLIQ(J) = 0.0
                 SNICE(J) = 0.0
                 DZSNSO(J) = 0.0
               ENDIF
!                SH2O(1) = SH2O(1)+SNLIQ(J)/(DZSNSO(1)*1000.)
!                SICE(1) = SICE(1)+SNICE(J)/(DZSNSO(1)*1000.)
             ENDIF

             ! shift all elements above this down by one.
             IF (J > ISNOW+1 .AND. ISNOW < -1) THEN
                DO I = J, ISNOW+2, -1
                   STC(I)   = STC(I-1)
                   SNLIQ(I) = SNLIQ(I-1)
                   SNICE(I) = SNICE(I-1)
                   DZSNSO(I)= DZSNSO(I-1)
                END DO
             END IF
             ISNOW = ISNOW + 1
          END IF
       END DO

! to conserve water in case of too large surface sublimation

       IF(SICE(1) < 0.) THEN
          SH2O(1) = SH2O(1) + SICE(1)
          SICE(1) = 0.
       END IF

       IF(ISNOW ==0) RETURN   ! MB: get out if no longer multi-layer

       SNEQV  = 0.
       SNOWH  = 0.
       ZWICE  = 0.
       ZWLIQ  = 0.

       DO J = ISNOW+1,0
             SNEQV = SNEQV + SNICE(J) + SNLIQ(J)
             SNOWH = SNOWH + DZSNSO(J)
             ZWICE = ZWICE + SNICE(J)
             ZWLIQ = ZWLIQ + SNLIQ(J)
       END DO

! check the snow depth - all snow gone
! the liquid water assumes ponding on soil surface.

       IF (SNOWH < 0.025 .AND. ISNOW < 0 ) THEN ! MB: change limit
!       IF (SNOWH < 0.05 .AND. ISNOW < 0 ) THEN
          ISNOW  = 0
          SNEQV = ZWICE
          PONDING2 = ZWLIQ           ! LIMIT OF ISNOW < 0 MEANS INPUT PONDING
          IF(SNEQV <= 0.) SNOWH = 0. ! SHOULD BE ZERO; SEE ABOVE
       END IF

!       IF (SNOWH < 0.05 ) THEN
!          ISNOW  = 0
!          SNEQV = ZWICE
!          SH2O(1) = SH2O(1) + ZWLIQ / (DZSNSO(1) * 1000.)
!          IF(SNEQV <= 0.) SNOWH = 0.
!       END IF

! check the snow depth - snow layers combined

       IF (ISNOW < -1) THEN

          ISNOW_OLD = ISNOW
          MSSI     = 1

          DO I = ISNOW_OLD+1,0
             IF (DZSNSO(I) < DZMIN(MSSI)) THEN

                IF (I == ISNOW+1) THEN
                   NEIBOR = I + 1
                ELSE IF (I == 0) THEN
                   NEIBOR = I - 1
                ELSE
                   NEIBOR = I + 1
                   IF ((DZSNSO(I-1)+DZSNSO(I)) < (DZSNSO(I+1)+DZSNSO(I))) NEIBOR = I-1
                END IF

                ! Node l and j are combined and stored as node j.
                IF (NEIBOR > I) THEN
                   J = NEIBOR
                   L = I
                ELSE
                   J = I
                   L = NEIBOR
                END IF

                CALL COMBO (parameters,DZSNSO(J), SNLIQ(J), SNICE(J), &
                   STC(J), DZSNSO(L), SNLIQ(L), SNICE(L), STC(L) )

                ! Now shift all elements above this down one.
                IF (J-1 > ISNOW+1) THEN
                   DO K = J-1, ISNOW+2, -1
                      STC(K)   = STC(K-1)
                      SNICE(K) = SNICE(K-1)
                      SNLIQ(K) = SNLIQ(K-1)
                      DZSNSO(K) = DZSNSO(K-1)
                   END DO
                END IF

                ! Decrease the number of snow layers
                ISNOW = ISNOW + 1
                IF (ISNOW >= -1) EXIT
             ELSE

                ! The layer thickness is greater than the prescribed minimum value
                MSSI = MSSI + 1

             END IF
          END DO

       END IF

  END SUBROUTINE COMBINE

!== begin divide ===================================================================================

  SUBROUTINE DIVIDE (parameters,NSNOW  ,NSOIL  ,                         & !in
                     ISNOW  ,STC    ,SNICE  ,SNLIQ  ,DZSNSO  )  !inout
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER, INTENT(IN)                            :: NSNOW !maximum no. of snow layers [ =3]
    INTEGER, INTENT(IN)                            :: NSOIL !no. of soil layers [ =4]

! input and output

    INTEGER                        , INTENT(INOUT) :: ISNOW !actual no. of snow layers 
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC   !snow layer temperature [k]
    REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNICE !snow layer ice [mm]
    REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(INOUT) :: SNLIQ !snow layer liquid water [mm]
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO!snow layer depth [m]

! local variables:

    INTEGER                                        :: J     !indices
    INTEGER                                        :: MSNO  !number of layer (top) to MSNO (bot)
    REAL(r8)                                           :: DRR   !thickness of the combined [m]
    REAL(r8), DIMENSION(       1:NSNOW)                :: DZ    !snow layer thickness [m]
    REAL(r8), DIMENSION(       1:NSNOW)                :: SWICE !partial volume of ice [m3/m3]
    REAL(r8), DIMENSION(       1:NSNOW)                :: SWLIQ !partial volume of liquid water [m3/m3]
    REAL(r8), DIMENSION(       1:NSNOW)                :: TSNO  !node temperature [k]
    REAL(r8)                                           :: ZWICE !temporary
    REAL(r8)                                           :: ZWLIQ !temporary
    REAL(r8)                                           :: PROPOR!temporary
    REAL(r8)                                           :: DTDZ  !temporary
! ----------------------------------------------------------------------

    DO J = 1,NSNOW
          IF (J <= ABS(ISNOW)) THEN
             DZ(J)    = DZSNSO(J+ISNOW)
             SWICE(J) = SNICE(J+ISNOW)
             SWLIQ(J) = SNLIQ(J+ISNOW)
             TSNO(J)  = STC(J+ISNOW)
          END IF
    END DO

       MSNO = ABS(ISNOW)

       IF (MSNO == 1) THEN
          ! Specify a new snow layer
          IF (DZ(1) > 0.05) THEN
             MSNO = 2
             DZ(1)    = DZ(1)/2.
             SWICE(1) = SWICE(1)/2.
             SWLIQ(1) = SWLIQ(1)/2.
             DZ(2)    = DZ(1)
             SWICE(2) = SWICE(1)
             SWLIQ(2) = SWLIQ(1)
             TSNO(2)  = TSNO(1)
          END IF
       END IF

       IF (MSNO > 1) THEN
          IF (DZ(1) > 0.05) THEN
             DRR      = DZ(1) - 0.05
             PROPOR   = DRR/DZ(1)
             ZWICE    = PROPOR*SWICE(1)
             ZWLIQ    = PROPOR*SWLIQ(1)
             PROPOR   = 0.05/DZ(1)
             SWICE(1) = PROPOR*SWICE(1)
             SWLIQ(1) = PROPOR*SWLIQ(1)
             DZ(1)    = 0.05

             CALL COMBO (parameters,DZ(2), SWLIQ(2), SWICE(2), TSNO(2), DRR, &
                  ZWLIQ, ZWICE, TSNO(1))

             ! subdivide a new layer
             IF (MSNO <= 2 .AND. DZ(2) > 0.20) THEN  ! MB: change limit
!             IF (MSNO <= 2 .AND. DZ(2) > 0.10) THEN
                MSNO = 3
                DTDZ = (TSNO(1) - TSNO(2))/((DZ(1)+DZ(2))/2.)
                DZ(2)    = DZ(2)/2.
                SWICE(2) = SWICE(2)/2.
                SWLIQ(2) = SWLIQ(2)/2.
                DZ(3)    = DZ(2)
                SWICE(3) = SWICE(2)
                SWLIQ(3) = SWLIQ(2)
                TSNO(3) = TSNO(2) - DTDZ*DZ(2)/2.
                IF (TSNO(3) >= TFRZ) THEN
                   TSNO(3)  = TSNO(2)
                ELSE
                   TSNO(2) = TSNO(2) + DTDZ*DZ(2)/2.
                ENDIF

             END IF
          END IF
       END IF

       IF (MSNO > 2) THEN
          IF (DZ(2) > 0.2) THEN
             DRR = DZ(2) - 0.2
             PROPOR   = DRR/DZ(2)
             ZWICE    = PROPOR*SWICE(2)
             ZWLIQ    = PROPOR*SWLIQ(2)
             PROPOR   = 0.2/DZ(2)
             SWICE(2) = PROPOR*SWICE(2)
             SWLIQ(2) = PROPOR*SWLIQ(2)
             DZ(2)    = 0.2
             CALL COMBO (parameters,DZ(3), SWLIQ(3), SWICE(3), TSNO(3), DRR, &
                  ZWLIQ, ZWICE, TSNO(2))
          END IF
       END IF

       ISNOW = -MSNO

    DO J = ISNOW+1,0
             DZSNSO(J) = DZ(J-ISNOW)
             SNICE(J) = SWICE(J-ISNOW)
             SNLIQ(J) = SWLIQ(J-ISNOW)
             STC(J)   = TSNO(J-ISNOW)
    END DO


!    DO J = ISNOW+1,NSOIL
!    WRITE(*,'(I5,7F10.3)') J, DZSNSO(J), SNICE(J), SNLIQ(J),STC(J)
!    END DO

  END SUBROUTINE DIVIDE

!== begin combo ====================================================================================

  SUBROUTINE COMBO(parameters,DZ,  WLIQ,  WICE, T, DZ2, WLIQ2, WICE2, T2)
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------s
! input

  type (noahmp_parameters), intent(in) :: parameters
    REAL(r8), INTENT(IN)    :: DZ2   !nodal thickness of 2 elements being combined [m]
    REAL(r8), INTENT(IN)    :: WLIQ2 !liquid water of element 2 [kg/m2]
    REAL(r8), INTENT(IN)    :: WICE2 !ice of element 2 [kg/m2]
    REAL(r8), INTENT(IN)    :: T2    !nodal temperature of element 2 [k]
    REAL(r8), INTENT(INOUT) :: DZ    !nodal thickness of 1 elements being combined [m]
    REAL(r8), INTENT(INOUT) :: WLIQ  !liquid water of element 1
    REAL(r8), INTENT(INOUT) :: WICE  !ice of element 1 [kg/m2]
    REAL(r8), INTENT(INOUT) :: T     !node temperature of element 1 [k]

! local 

    REAL(r8)                :: DZC   !total thickness of nodes 1 and 2 (DZC=DZ+DZ2).
    REAL(r8)                :: WLIQC !combined liquid water [kg/m2]
    REAL(r8)                :: WICEC !combined ice [kg/m2]
    REAL(r8)                :: TC    !combined node temperature [k]
    REAL(r8)                :: H     !enthalpy of element 1 [J/m2]
    REAL(r8)                :: H2    !enthalpy of element 2 [J/m2]
    REAL(r8)                :: HC    !temporary

!-----------------------------------------------------------------------

    DZC = DZ+DZ2
    WICEC = (WICE+WICE2)
    WLIQC = (WLIQ+WLIQ2)
    H = (CICE*WICE+CWAT*WLIQ) * (T-TFRZ)+HFUS*WLIQ
    H2= (CICE*WICE2+CWAT*WLIQ2) * (T2-TFRZ)+HFUS*WLIQ2

    HC = H + H2
    IF(HC < 0.)THEN
       TC = TFRZ + HC/(CICE*WICEC + CWAT*WLIQC)
    ELSE IF (HC.LE.HFUS*WLIQC) THEN
       TC = TFRZ
    ELSE
       TC = TFRZ + (HC - HFUS*WLIQC) / (CICE*WICEC + CWAT*WLIQC)
    END IF

    DZ = DZC
    WICE = WICEC
    WLIQ = WLIQC
    T = TC

  END SUBROUTINE COMBO

!== begin compact ==================================================================================

  SUBROUTINE COMPACT (parameters,NSNOW  ,NSOIL  ,DT     ,STC    ,SNICE  , & !in
                      SNLIQ  ,ZSOIL  ,IMELT  ,FICEOLD,ILOC   , JLOC , & !in
                      ISNOW  ,DZSNSO ,ZSNSO )                    !inout
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
   INTEGER,                         INTENT(IN)    :: ILOC   !grid index
   INTEGER,                         INTENT(IN)    :: JLOC   !grid index
   INTEGER,                         INTENT(IN)    :: NSOIL  !no. of soil layers [ =4]
   INTEGER,                         INTENT(IN)    :: NSNOW  !maximum no. of snow layers [ =3]
   INTEGER, DIMENSION(-NSNOW+1:0) , INTENT(IN)    :: IMELT  !melting state index [0-no melt;1-melt]
   REAL(r8),                            INTENT(IN)    :: DT     !time step (sec)
   REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN)    :: STC    !snow layer temperature [k]
   REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: SNICE  !snow layer ice [mm]
   REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: SNLIQ  !snow layer liquid water [mm]
   REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN)    :: ZSOIL  !depth of layer-bottom from soil srf
   REAL(r8), DIMENSION(-NSNOW+1:    0), INTENT(IN)    :: FICEOLD!ice fraction at last timestep

! input and output
   INTEGER,                         INTENT(INOUT) :: ISNOW  ! actual no. of snow layers
   REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO ! snow layer thickness [m]
   REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: ZSNSO  ! depth of snow/soil layer-bottom

! local
   REAL(r8), PARAMETER     :: C2 = 21.e-3   ![m3/kg] ! default 21.e-3
   REAL(r8), PARAMETER     :: C3 = 2.5e-6   ![1/s]  
   REAL(r8), PARAMETER     :: C4 = 0.04     ![1/k]
   REAL(r8), PARAMETER     :: C5 = 2.0      !
   REAL(r8), PARAMETER     :: DM = 100.0    !upper Limit on destructive metamorphism compaction [kg/m3]
   REAL(r8), PARAMETER     :: ETA0 = 0.8e+6 !viscosity coefficient [kg-s/m2] 
                                        !according to Anderson, it is between 0.52e6~1.38e6
   REAL(r8) :: BURDEN !pressure of overlying snow [kg/m2]
   REAL(r8) :: DDZ1   !rate of settling of snow pack due to destructive metamorphism.
   REAL(r8) :: DDZ2   !rate of compaction of snow pack due to overburden.
   REAL(r8) :: DDZ3   !rate of compaction of snow pack due to melt [1/s]
   REAL(r8) :: DEXPF  !EXPF=exp(-c4*(273.15-STC)).
   REAL(r8) :: TD     !STC - TFRZ [K]
   REAL(r8) :: PDZDTC !nodal rate of change in fractional-thickness due to compaction [fraction/s]
   REAL(r8) :: VOID   !void (1 - SNICE - SNLIQ)
   REAL(r8) :: WX     !water mass (ice + liquid) [kg/m2]
   REAL(r8) :: BI     !partial density of ice [kg/m3]
   REAL(r8), DIMENSION(-NSNOW+1:0) :: FICE   !fraction of ice at current time step

   INTEGER  :: J

! ----------------------------------------------------------------------
    BURDEN = 0.0

    DO J = ISNOW+1, 0

        WX      = SNICE(J) + SNLIQ(J)
        FICE(J) = SNICE(J) / WX
        VOID    = 1. - (SNICE(J)/DENICE + SNLIQ(J)/DENH2O) / DZSNSO(J)

        ! Allow compaction only for non-saturated node and higher ice lens node.
        IF (VOID > 0.001 .AND. SNICE(J) > 0.1) THEN
           BI = SNICE(J) / DZSNSO(J)
           TD = MAX(0.,TFRZ-STC(J))
           DEXPF = EXP(-C4*TD)

           ! Settling as a result of destructive metamorphism

           DDZ1 = -C3*DEXPF

           IF (BI > DM) DDZ1 = DDZ1*EXP(-46.0E-3*(BI-DM))

           ! Liquid water term

           IF (SNLIQ(J) > 0.01*DZSNSO(J)) DDZ1=DDZ1*C5

           ! Compaction due to overburden

           DDZ2 = -(BURDEN+0.5*WX)*EXP(-0.08*TD-C2*BI)/ETA0 ! 0.5*WX -> self-burden

           ! Compaction occurring during melt

           IF (IMELT(J) == 1) THEN
              DDZ3 = MAX(0.,(FICEOLD(J) - FICE(J))/MAX(1.E-6,FICEOLD(J)))
              DDZ3 = - DDZ3/DT           ! sometimes too large
           ELSE
              DDZ3 = 0.
           END IF

           ! Time rate of fractional change in DZ (units of s-1)

           PDZDTC = (DDZ1 + DDZ2 + DDZ3)*DT
           PDZDTC = MAX(-0.5,PDZDTC)

           ! The change in DZ due to compaction

           DZSNSO(J) = DZSNSO(J)*(1.+PDZDTC)
        END IF

        ! Pressure of overlying snow

        BURDEN = BURDEN + WX

    END DO

  END SUBROUTINE COMPACT

!== begin snowh2o ==================================================================================

  SUBROUTINE SNOWH2O (parameters,NSNOW  ,NSOIL  ,DT     ,QSNFRO ,QSNSUB , & !in 
                      QRAIN  ,ILOC   ,JLOC   ,                 & !in
                      ISNOW  ,DZSNSO ,SNOWH  ,SNEQV  ,SNICE  , & !inout
                      SNLIQ  ,SH2O   ,SICE   ,STC    ,         & !inout
                      QSNBOT ,PONDING1       ,PONDING2)          !out
! ----------------------------------------------------------------------
! Renew the mass of ice lens (SNICE) and liquid (SNLIQ) of the
! surface snow layer resulting from sublimation (frost) / evaporation (dew)
! ----------------------------------------------------------------------
   IMPLICIT NONE
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
   INTEGER,                         INTENT(IN)    :: ILOC   !grid index
   INTEGER,                         INTENT(IN)    :: JLOC   !grid index
   INTEGER,                         INTENT(IN)    :: NSNOW  !maximum no. of snow layers[=3]
   INTEGER,                         INTENT(IN)    :: NSOIL  !No. of soil layers[=4]
   REAL(r8),                            INTENT(IN)    :: DT     !time step
   REAL(r8),                            INTENT(IN)    :: QSNFRO !snow surface frost rate[mm/s]
   REAL(r8),                            INTENT(IN)    :: QSNSUB !snow surface sublimation rate[mm/s]
   REAL(r8),                            INTENT(IN)    :: QRAIN  !snow surface rain rate[mm/s]

! output

   REAL(r8),                            INTENT(OUT)   :: QSNBOT !melting water out of snow bottom [mm/s]

! input and output

   INTEGER,                         INTENT(INOUT) :: ISNOW  !actual no. of snow layers
   REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: DZSNSO ! snow layer depth [m]
   REAL(r8),                            INTENT(INOUT) :: SNOWH  !snow height [m]
   REAL(r8),                            INTENT(INOUT) :: SNEQV  !snow water eqv. [mm]
   REAL(r8), DIMENSION(-NSNOW+1:0),     INTENT(INOUT) :: SNICE  !snow layer ice [mm]
   REAL(r8), DIMENSION(-NSNOW+1:0),     INTENT(INOUT) :: SNLIQ  !snow layer liquid water [mm]
   REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SH2O   !soil liquid moisture (m3/m3)
   REAL(r8), DIMENSION(       1:NSOIL), INTENT(INOUT) :: SICE   !soil ice moisture (m3/m3)
   REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(INOUT) :: STC    !snow layer temperature [k]

! local variables:

   INTEGER                     :: J         !do loop/array indices
   REAL(r8)                        :: QIN       !water flow into the element (mm/s)
   REAL(r8)                        :: QOUT      !water flow out of the element (mm/s)
   REAL(r8)                        :: WGDIF     !ice mass after minus sublimation
   REAL(r8), DIMENSION(-NSNOW+1:0) :: VOL_LIQ   !partial volume of liquid water in layer
   REAL(r8), DIMENSION(-NSNOW+1:0) :: VOL_ICE   !partial volume of ice lens in layer
   REAL(r8), DIMENSION(-NSNOW+1:0) :: EPORE     !effective porosity = porosity - VOL_ICE
   REAL(r8) :: PROPOR, TEMP
   REAL(r8) :: PONDING1, PONDING2
! ----------------------------------------------------------------------

!for the case when SNEQV becomes '0' after 'COMBINE'

   IF(SNEQV == 0.) THEN
      SICE(1) =  SICE(1) + (QSNFRO-QSNSUB)*DT/(DZSNSO(1)*1000.)  ! Barlage: SH2O->SICE v3.6
      IF(SICE(1) < 0.) THEN
         SH2O(1) = SH2O(1) + SICE(1)
         SICE(1) = 0.
      END IF
   END IF

! for shallow snow without a layer
! snow surface sublimation may be larger than existing snow mass. To conserve water,
! excessive sublimation is used to reduce soil water. Smaller time steps would tend 
! to aviod this problem.

   IF(ISNOW == 0 .and. SNEQV > 0.) THEN
      TEMP   = SNEQV
      SNEQV  = SNEQV - QSNSUB*DT + QSNFRO*DT
      PROPOR = SNEQV/TEMP
      SNOWH  = MAX(0.,PROPOR * SNOWH)

      IF(SNEQV < 0.) THEN
         SICE(1) = SICE(1) + SNEQV/(DZSNSO(1)*1000.)
         SNEQV   = 0.
         SNOWH   = 0.
      END IF
      IF(SICE(1) < 0.) THEN
         SH2O(1) = SH2O(1) + SICE(1)
         SICE(1) = 0.
      END IF
   END IF

   IF(SNOWH <= 1.E-8 .OR. SNEQV <= 1.E-6) THEN
     SNOWH = 0.0
     SNEQV = 0.0
   END IF

! for deep snow

   IF ( ISNOW < 0 ) THEN !KWM added this IF statement to prevent out-of-bounds array references

      WGDIF = SNICE(ISNOW+1) - QSNSUB*DT + QSNFRO*DT
      SNICE(ISNOW+1) = WGDIF
      IF (WGDIF < 1.e-6 .and. ISNOW <0) THEN
         CALL  COMBINE (parameters,NSNOW  ,NSOIL  ,ILOC, JLOC   , & !in
              ISNOW  ,SH2O   ,STC    ,SNICE  ,SNLIQ  , & !inout
              DZSNSO ,SICE   ,SNOWH  ,SNEQV  ,         & !inout
              PONDING1, PONDING2 )                       !out
      ENDIF
      !KWM:  Subroutine COMBINE can change ISNOW to make it 0 again?
      IF ( ISNOW < 0 ) THEN !KWM added this IF statement to prevent out-of-bounds array references
         SNLIQ(ISNOW+1) = SNLIQ(ISNOW+1) + QRAIN * DT
         SNLIQ(ISNOW+1) = MAX(0., SNLIQ(ISNOW+1))
      ENDIF
      
   ENDIF !KWM  -- Can the ENDIF be moved toward the end of the subroutine (Just set QSNBOT=0)?

! Porosity and partial volume

   !KWM Looks to me like loop index / IF test can be simplified.

   DO J = -NSNOW+1, 0
      IF (J >= ISNOW+1) THEN
         VOL_ICE(J)      = MIN(1., SNICE(J)/(DZSNSO(J)*DENICE))
         EPORE(J)        = 1. - VOL_ICE(J)
         VOL_LIQ(J)      = MIN(EPORE(J),SNLIQ(J)/(DZSNSO(J)*DENH2O))
      END IF
   END DO

   QIN = 0.
   QOUT = 0.

   !KWM Looks to me like loop index / IF test can be simplified.

   DO J = -NSNOW+1, 0
      IF (J >= ISNOW+1) THEN
         SNLIQ(J) = SNLIQ(J) + QIN
         IF (J <= -1) THEN
            IF (EPORE(J) < 0.05 .OR. EPORE(J+1) < 0.05) THEN
               QOUT = 0.
            ELSE
               QOUT = MAX(0.,(VOL_LIQ(J)-parameters%SSI*EPORE(J))*DZSNSO(J))
               QOUT = MIN(QOUT,(1.-VOL_ICE(J+1)-VOL_LIQ(J+1))*DZSNSO(J+1))
            END IF
         ELSE
            QOUT = MAX(0.,(VOL_LIQ(J) - parameters%SSI*EPORE(J))*DZSNSO(J))
         END IF
         QOUT = QOUT*1000.
         SNLIQ(J) = SNLIQ(J) - QOUT
         QIN = QOUT
      END IF
   END DO

! Liquid water from snow bottom to soil

   QSNBOT = QOUT / DT           ! mm/s

  END SUBROUTINE SNOWH2O

!== begin soilwater ================================================================================

  SUBROUTINE SOILWATER (parameters,NSOIL  ,NSNOW  ,DT     ,ZSOIL  ,DZSNSO , & !in
                        QINSUR ,QSEVA  ,ETRANI ,SICE   ,ILOC   , JLOC, & !in
                        SH2O   ,SMC    ,ZWT    ,VEGTYP ,& !inout
                        SMCWTD, DEEPRECH                       ,& !inout
                        RUNSRF ,QDRAIN ,RUNSUB ,WCND   ,FCRMAX )   !out

! ----------------------------------------------------------------------
! calculate surface runoff and soil moisture.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                     INTENT(IN) :: ILOC   !grid index
  INTEGER,                     INTENT(IN) :: JLOC   !grid index
  INTEGER,                     INTENT(IN) :: NSOIL  !no. of soil layers
  INTEGER,                     INTENT(IN) :: NSNOW  !maximum no. of snow layers
  REAL(r8),                        INTENT(IN) :: DT     !time step (sec)
  REAL(r8), INTENT(IN)                        :: QINSUR !water input on soil surface [mm/s]
  REAL(r8), INTENT(IN)                        :: QSEVA  !evap from soil surface [mm/s]
  REAL(r8), DIMENSION(1:NSOIL),    INTENT(IN) :: ZSOIL  !depth of soil layer-bottom [m]
  REAL(r8), DIMENSION(1:NSOIL),    INTENT(IN) :: ETRANI !evapotranspiration from soil layers [mm/s]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !snow/soil layer depth [m]
  REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)   :: SICE   !soil ice content [m3/m3]

  INTEGER,                     INTENT(IN) :: VEGTYP

! input & output
  REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: SH2O   !soil liquid water content [m3/m3]
  REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: SMC    !total soil water content [m3/m3]
  REAL(r8), INTENT(INOUT)                     :: ZWT    !water table depth [m]
  REAL(r8),                     INTENT(INOUT) :: SMCWTD !soil moisture between bottom of the soil and the water table [m3/m3]
  REAL(r8)                    , INTENT(INOUT) :: DEEPRECH

! output
  REAL(r8), INTENT(OUT)                       :: QDRAIN !soil-bottom free drainage [mm/s] 
  REAL(r8), INTENT(OUT)                       :: RUNSRF !surface runoff [mm/s] 
  REAL(r8), INTENT(OUT)                       :: RUNSUB !subsurface runoff [mm/s] 
  REAL(r8), INTENT(OUT)                       :: FCRMAX !maximum of FCR (-)
  REAL(r8), DIMENSION(1:NSOIL), INTENT(OUT)   :: WCND   !hydraulic conductivity (m/s)

! local
  INTEGER                                 :: K,IZ   !do-loop index
  INTEGER                                 :: ITER   !iteration index
  REAl                                    :: DTFINE !fine time step (s)
  REAL(r8), DIMENSION(1:NSOIL)                :: RHSTT  !right-hand side term of the matrix
  REAL(r8), DIMENSION(1:NSOIL)                :: AI     !left-hand side term
  REAL(r8), DIMENSION(1:NSOIL)                :: BI     !left-hand side term
  REAL(r8), DIMENSION(1:NSOIL)                :: CI     !left-hand side term

  REAL(r8)                                    :: FFF    !runoff decay factor (m-1)
  REAL(r8)                                    :: RSBMX  !baseflow coefficient [mm/s]
  REAL(r8)                                    :: PDDUM  !infiltration rate at surface (m/s)
  REAL(r8)                                    :: FICE   !ice fraction in frozen soil
  REAL(r8)                                    :: WPLUS  !saturation excess of the total soil [m]
  REAL(r8)                                    :: RSAT   !accumulation of WPLUS (saturation excess) [m]
  REAL(r8)                                    :: SICEMAX!maximum soil ice content (m3/m3)
  REAL(r8)                                    :: SH2OMIN!minimum soil liquid water content (m3/m3)
  REAL(r8)                                    :: WTSUB  !sum of WCND(K)*DZSNSO(K)
  REAL(r8)                                    :: MH2O   !water mass removal (mm)
  REAL(r8)                                    :: FSAT   !fractional saturated area (-)
  REAL(r8), DIMENSION(1:NSOIL)                :: MLIQ   !
  REAL(r8)                                    :: XS     !
  REAL(r8)                                    :: WATMIN !
  REAL(r8)                                    :: QDRAIN_SAVE !
  REAL(r8)                                    :: EPORE  !effective porosity [m3/m3]
  REAL(r8), DIMENSION(1:NSOIL)                :: FCR    !impermeable fraction due to frozen soil
  INTEGER                                 :: NITER  !iteration times soil moisture (-)
  REAL(r8)                                    :: SMCTOT !2-m averaged soil moisture (m3/m3)
  REAL(r8)                                    :: DZTOT  !2-m soil depth (m)
  REAL(r8), PARAMETER :: A = 4.0
! ----------------------------------------------------------------------
    RUNSRF = 0.0
    PDDUM  = 0.0
    RSAT   = 0.0

! for the case when snowmelt water is too large

    DO K = 1,NSOIL
       EPORE   = MAX ( 1.E-4 , ( parameters%SMCMAX(K) - SICE(K) ) )
       RSAT    = RSAT + MAX(0.,SH2O(K)-EPORE)*DZSNSO(K)  
       SH2O(K) = MIN(EPORE,SH2O(K))             
    END DO

!impermeable fraction due to frozen soil

    DO K = 1,NSOIL
       FICE    = MIN(1.0,SICE(K)/parameters%SMCMAX(K))
       FCR(K)  = MAX(0.0,EXP(-A*(1.-FICE))- EXP(-A)) /  &
                        (1.0              - EXP(-A))
    END DO

! maximum soil ice content and minimum liquid water of all layers

    SICEMAX = 0.0
    FCRMAX  = 0.0
    SH2OMIN = parameters%SMCMAX(1)
    DO K = 1,NSOIL
       IF (SICE(K) > SICEMAX) SICEMAX = SICE(K)
       IF (FCR(K)  > FCRMAX)  FCRMAX  = FCR(K)
       IF (SH2O(K) < SH2OMIN) SH2OMIN = SH2O(K)
    END DO

!subsurface runoff for runoff scheme option 2

    IF(OPT_RUN == 2) THEN 
        FFF   = 2.0
        RSBMX = 4.0
        CALL ZWTEQ (parameters,NSOIL  ,NSNOW  ,ZSOIL  ,DZSNSO ,SH2O   ,ZWT)
        RUNSUB = (1.0-FCRMAX) * RSBMX * EXP(-parameters%TIMEAN) * EXP(-FFF*ZWT)   ! mm/s
    END IF

!surface runoff and infiltration rate using different schemes

!jref impermable surface at urban
    IF ( parameters%urban_flag ) FCR(1)= 0.95

    IF(OPT_RUN == 1) THEN
       FFF = 6.0
       FSAT   = parameters%FSATMX*EXP(-0.5*FFF*(ZWT-2.0))
       IF(QINSUR > 0.) THEN
         RUNSRF = QINSUR * ( (1.0-FCR(1))*FSAT + FCR(1) )
         PDDUM  = QINSUR - RUNSRF                          ! m/s 
       END IF
    END IF

    IF(OPT_RUN == 5) THEN
       FFF = 6.0
       FSAT   = parameters%FSATMX*EXP(-0.5*FFF*MAX(-2.0-ZWT,0.))
       IF(QINSUR > 0.) THEN
         RUNSRF = QINSUR * ( (1.0-FCR(1))*FSAT + FCR(1) )
         PDDUM  = QINSUR - RUNSRF                          ! m/s
       END IF
    END IF

    IF(OPT_RUN == 2) THEN
       FFF   = 2.0
       FSAT   = parameters%FSATMX*EXP(-0.5*FFF*ZWT)
       IF(QINSUR > 0.) THEN
         RUNSRF = QINSUR * ( (1.0-FCR(1))*FSAT + FCR(1) )
         PDDUM  = QINSUR - RUNSRF                          ! m/s 
       END IF
    END IF

    IF(OPT_RUN == 3) THEN
       CALL INFIL (parameters,NSOIL  ,DT     ,ZSOIL  ,SH2O   ,SICE   , & !in
                   SICEMAX,QINSUR ,                         & !in
                   PDDUM  ,RUNSRF )                           !out
    END IF

    IF(OPT_RUN == 4) THEN
       SMCTOT = 0.
       DZTOT  = 0.
       DO K = 1,NSOIL
          DZTOT   = DZTOT  + DZSNSO(K)  
          SMCTOT  = SMCTOT + SMC(K)/parameters%SMCMAX(K)*DZSNSO(K)
          IF(DZTOT >= 2.0) EXIT
       END DO
       SMCTOT = SMCTOT/DZTOT
       FSAT   = MAX(0.01,SMCTOT) ** 4.        !BATS

       IF(QINSUR > 0.) THEN
         RUNSRF = QINSUR * ((1.0-FCR(1))*FSAT+FCR(1))  
         PDDUM  = QINSUR - RUNSRF                       ! m/s
       END IF
    END IF

! determine iteration times and finer time step

    NITER = 1

    IF(OPT_INF == 1) THEN    !OPT_INF =2 may cause water imbalance
       NITER = 3
       IF (PDDUM*DT>DZSNSO(1)*parameters%SMCMAX(1) ) THEN
          NITER = NITER*2
       END IF
    END IF                 

    DTFINE  = DT / NITER

! solve soil moisture

    QDRAIN_SAVE = 0.0
    DO ITER = 1, NITER
       CALL SRT   (parameters,NSOIL  ,ZSOIL  ,DTFINE ,PDDUM  ,ETRANI , & !in
                   QSEVA  ,SH2O   ,SMC    ,ZWT    ,FCR    , & !in
                   SICEMAX,FCRMAX ,ILOC   ,JLOC   ,SMCWTD ,         & !in
                   RHSTT  ,AI     ,BI     ,CI     ,QDRAIN , & !out
                   WCND   )                                   !out
  
       CALL SSTEP (parameters,NSOIL  ,NSNOW  ,DTFINE ,ZSOIL  ,DZSNSO , & !in
                   SICE   ,ILOC   ,JLOC   ,ZWT            ,                 & !in
                   SH2O   ,SMC    ,AI     ,BI     ,CI     , & !inout
                   RHSTT  ,SMCWTD ,QDRAIN ,DEEPRECH,                                 & !inout
                   WPLUS)                                     !out
       RSAT =  RSAT + WPLUS
       QDRAIN_SAVE = QDRAIN_SAVE + QDRAIN
    END DO

    QDRAIN = QDRAIN_SAVE/NITER

    RUNSRF = RUNSRF * 1000. + RSAT * 1000./DT  ! m/s -> mm/s
    QDRAIN = QDRAIN * 1000.

!WRF_HYDRO_DJG...
!yw    INFXSRT = RUNSRF * DT   !mm/s -> mm

! removal of soil water due to groundwater flow (option 2)

    IF(OPT_RUN == 2) THEN
         WTSUB = 0.
         DO K = 1, NSOIL
           WTSUB = WTSUB + WCND(K)*DZSNSO(K)
         END DO

         DO K = 1, NSOIL
           MH2O    = RUNSUB*DT*(WCND(K)*DZSNSO(K))/WTSUB       ! mm
           SH2O(K) = SH2O(K) - MH2O/(DZSNSO(K)*1000.)
         END DO
    END IF

! Limit MLIQ to be greater than or equal to watmin.
! Get water needed to bring MLIQ equal WATMIN from lower layer.

   IF(OPT_RUN /= 1) THEN
      DO IZ = 1, NSOIL
         MLIQ(IZ) = SH2O(IZ)*DZSNSO(IZ)*1000.
      END DO

      WATMIN = 0.01           ! mm
      DO IZ = 1, NSOIL-1
          IF (MLIQ(IZ) .LT. 0.) THEN
             XS = WATMIN-MLIQ(IZ)
          ELSE
             XS = 0.
          END IF
          MLIQ(IZ  ) = MLIQ(IZ  ) + XS
          MLIQ(IZ+1) = MLIQ(IZ+1) - XS
      END DO

        IZ = NSOIL
        IF (MLIQ(IZ) .LT. WATMIN) THEN
           XS = WATMIN-MLIQ(IZ)
        ELSE
           XS = 0.
        END IF
        MLIQ(IZ) = MLIQ(IZ) + XS
        RUNSUB   = RUNSUB - XS/DT
        IF(OPT_RUN == 5)DEEPRECH = DEEPRECH - XS*1.E-3

      DO IZ = 1, NSOIL
        SH2O(IZ)     = MLIQ(IZ) / (DZSNSO(IZ)*1000.)
      END DO
   END IF

  END SUBROUTINE SOILWATER

!== begin zwteq ====================================================================================

  SUBROUTINE ZWTEQ (parameters,NSOIL  ,NSNOW  ,ZSOIL  ,DZSNSO ,SH2O   ,ZWT)
! ----------------------------------------------------------------------
! calculate equilibrium water table depth (Niu et al., 2005)
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                         INTENT(IN) :: NSOIL  !no. of soil layers
  INTEGER,                         INTENT(IN) :: NSNOW  !maximum no. of snow layers
  REAL(r8), DIMENSION(1:NSOIL),        INTENT(IN) :: ZSOIL  !depth of soil layer-bottom [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !snow/soil layer depth [m]
  REAL(r8), DIMENSION(1:NSOIL),        INTENT(IN) :: SH2O   !soil liquid water content [m3/m3]

! output

  REAL(r8),                           INTENT(OUT) :: ZWT    !water table depth [m]

! locals

  INTEGER :: K                      !do-loop index
  INTEGER, PARAMETER :: NFINE = 100 !no. of fine soil layers of 6m soil
  REAL(r8)    :: WD1                    !water deficit from coarse (4-L) soil moisture profile
  REAL(r8)    :: WD2                    !water deficit from fine (100-L) soil moisture profile
  REAL(r8)    :: DZFINE                 !layer thickness of the 100-L soil layers to 6.0 m
  REAL(r8)    :: TEMP                   !temporary variable
  REAL(r8), DIMENSION(1:NFINE) :: ZFINE !layer-bottom depth of the 100-L soil layers to 6.0 m
! ----------------------------------------------------------------------

   WD1 = 0.
   DO K = 1,NSOIL
     WD1 = WD1 + (parameters%SMCMAX(K)-SH2O(K)) * DZSNSO(K) ! [m]
   ENDDO

   DZFINE = 3.0 * (-ZSOIL(NSOIL)) / NFINE  
   do K =1,NFINE
      ZFINE(K) = FLOAT(K) * DZFINE
   ENDDO

   ZWT = -3.*ZSOIL(NSOIL) - 0.001   ! initial value [m]

   WD2 = 0.
   DO K = 1,NFINE
     TEMP  = 1. + (ZWT-ZFINE(K))/parameters%PSISAT(K)
     WD2   = WD2 + parameters%SMCMAX(K)*(1.-TEMP**(-1./parameters%BEXP(K)))*DZFINE
     IF(ABS(WD2-WD1).LE.0.01) THEN
        ZWT = ZFINE(K)
        EXIT
     ENDIF
   ENDDO

  END SUBROUTINE ZWTEQ

!== begin infil ====================================================================================

  SUBROUTINE INFIL (parameters,NSOIL  ,DT     ,ZSOIL  ,SH2O   ,SICE   , & !in
                    SICEMAX,QINSUR ,                         & !in
                    PDDUM  ,RUNSRF )                           !out
! --------------------------------------------------------------------------------
! compute inflitration rate at soil surface and surface runoff
! --------------------------------------------------------------------------------
    IMPLICIT NONE
! --------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                  INTENT(IN) :: NSOIL  !no. of soil layers
  REAL(r8),                     INTENT(IN) :: DT     !time step (sec)
  REAL(r8), DIMENSION(1:NSOIL), INTENT(IN) :: ZSOIL  !depth of soil layer-bottom [m]
  REAL(r8), DIMENSION(1:NSOIL), INTENT(IN) :: SH2O   !soil liquid water content [m3/m3]
  REAL(r8), DIMENSION(1:NSOIL), INTENT(IN) :: SICE   !soil ice content [m3/m3]
  REAL(r8),                     INTENT(IN) :: QINSUR !water input on soil surface [mm/s]
  REAL(r8),                     INTENT(IN) :: SICEMAX!maximum soil ice content (m3/m3)

! outputs
  REAL(r8),                    INTENT(OUT) :: RUNSRF !surface runoff [mm/s] 
  REAL(r8),                    INTENT(OUT) :: PDDUM  !infiltration rate at surface

! locals
  INTEGER :: IALP1, J, JJ,  K
  REAL(r8)                     :: VAL
  REAL(r8)                     :: DDT
  REAL(r8)                     :: PX
  REAL(r8)                     :: DT1, DD, DICE
  REAL(r8)                     :: FCR
  REAL(r8)                     :: SUM
  REAL(r8)                     :: ACRT
  REAL(r8)                     :: WDF
  REAL(r8)                     :: WCND
  REAL(r8)                     :: SMCAV
  REAL(r8)                     :: INFMAX
  REAL(r8), DIMENSION(1:NSOIL) :: DMAX
  INTEGER, PARAMETER       :: CVFRZ = 3
! --------------------------------------------------------------------------------

    IF (QINSUR >  0.0) THEN
       DT1 = DT /86400.
       SMCAV = parameters%SMCMAX(1) - parameters%SMCWLT(1)

! maximum infiltration rate

       DMAX(1)= -ZSOIL(1) * SMCAV
       DICE   = -ZSOIL(1) * SICE(1)
       DMAX(1)= DMAX(1)* (1.0-(SH2O(1) + SICE(1) - parameters%SMCWLT(1))/SMCAV)

       DD = DMAX(1)

       DO K = 2,NSOIL
          DICE    = DICE + (ZSOIL(K-1) - ZSOIL(K) ) * SICE(K)
          DMAX(K) = (ZSOIL(K-1) - ZSOIL(K)) * SMCAV
          DMAX(K) = DMAX(K) * (1.0-(SH2O(K) + SICE(K) - parameters%SMCWLT(K))/SMCAV)
          DD      = DD + DMAX(K)
       END DO

       VAL = (1. - EXP ( - parameters%KDT * DT1))
       DDT = DD * VAL
       PX  = MAX(0.,QINSUR * DT)
       INFMAX = (PX * (DDT / (PX + DDT)))/ DT

! impermeable fraction due to frozen soil

       FCR = 1.
       IF (DICE >  1.E-2) THEN
          ACRT = CVFRZ * parameters%FRZX / DICE
          SUM = 1.
          IALP1 = CVFRZ - 1
          DO J = 1,IALP1
             K = 1
             DO JJ = J +1,IALP1
                K = K * JJ
             END DO
             SUM = SUM + (ACRT ** (CVFRZ - J)) / FLOAT(K)
          END DO
          FCR = 1. - EXP (-ACRT) * SUM
       END IF

! correction of infiltration limitation

       INFMAX = INFMAX * FCR

! jref for urban areas
!       IF ( parameters%urban_flag ) INFMAX == INFMAX * 0.05

       CALL WDFCND2 (parameters,WDF,WCND,SH2O(1),SICEMAX,1)
       INFMAX = MAX (INFMAX,WCND)
       INFMAX = MIN (INFMAX,PX)

       RUNSRF= MAX(0., QINSUR - INFMAX)
       PDDUM = QINSUR - RUNSRF

    END IF

  END SUBROUTINE INFIL

!== begin srt ======================================================================================

  SUBROUTINE SRT (parameters,NSOIL  ,ZSOIL  ,DT     ,PDDUM  ,ETRANI , & !in
                  QSEVA  ,SH2O   ,SMC    ,ZWT    ,FCR    , & !in
                  SICEMAX,FCRMAX ,ILOC   ,JLOC   ,SMCWTD ,         & !in
                  RHSTT  ,AI     ,BI     ,CI     ,QDRAIN , & !out
                  WCND   )                                   !out
! ----------------------------------------------------------------------
! calculate the right hand side of the time tendency term of the soil
! water diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,                  INTENT(IN)  :: ILOC   !grid index
    INTEGER,                  INTENT(IN)  :: JLOC   !grid index
    INTEGER,                  INTENT(IN)  :: NSOIL
    REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)  :: ZSOIL
    REAL(r8),                     INTENT(IN)  :: DT
    REAL(r8),                     INTENT(IN)  :: PDDUM
    REAL(r8),                     INTENT(IN)  :: QSEVA
    REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)  :: ETRANI
    REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)  :: SH2O
    REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)  :: SMC
    REAL(r8),                     INTENT(IN)  :: ZWT    ! water table depth [m]
    REAL(r8), DIMENSION(1:NSOIL), INTENT(IN)  :: FCR
    REAL(r8), INTENT(IN)                      :: FCRMAX !maximum of FCR (-)
    REAL(r8),                     INTENT(IN)  :: SICEMAX!maximum soil ice content (m3/m3)
    REAL(r8),                     INTENT(IN)  :: SMCWTD !soil moisture between bottom of the soil and the water table

! output

    REAL(r8), DIMENSION(1:NSOIL), INTENT(OUT) :: RHSTT
    REAL(r8), DIMENSION(1:NSOIL), INTENT(OUT) :: AI
    REAL(r8), DIMENSION(1:NSOIL), INTENT(OUT) :: BI
    REAL(r8), DIMENSION(1:NSOIL), INTENT(OUT) :: CI
    REAL(r8), DIMENSION(1:NSOIL), INTENT(OUT) :: WCND    !hydraulic conductivity (m/s)
    REAL(r8),                     INTENT(OUT) :: QDRAIN  !bottom drainage (m/s)

! local
    INTEGER                               :: K
    REAL(r8), DIMENSION(1:NSOIL)              :: DDZ
    REAL(r8), DIMENSION(1:NSOIL)              :: DENOM
    REAL(r8), DIMENSION(1:NSOIL)              :: DSMDZ
    REAL(r8), DIMENSION(1:NSOIL)              :: WFLUX
    REAL(r8), DIMENSION(1:NSOIL)              :: WDF
    REAL(r8), DIMENSION(1:NSOIL)              :: SMX
    REAL(r8)                                  :: TEMP1
    REAL(r8)                                  :: SMXWTD !soil moisture between bottom of the soil and water table
    REAL(r8)                                  :: SMXBOT  !soil moisture below bottom to calculate flux

! Niu and Yang (2006), J. of Hydrometeorology
! ----------------------------------------------------------------------

    IF(OPT_INF == 1) THEN
      DO K = 1, NSOIL
        CALL WDFCND1 (parameters,WDF(K),WCND(K),SMC(K),FCR(K),K)
        SMX(K) = SMC(K)
      END DO
        IF(OPT_RUN == 5)SMXWTD=SMCWTD
    END IF

    IF(OPT_INF == 2) THEN
      DO K = 1, NSOIL
        CALL WDFCND2 (parameters,WDF(K),WCND(K),SH2O(K),SICEMAX,K)
        SMX(K) = SH2O(K)
      END DO
          IF(OPT_RUN == 5)SMXWTD=SMCWTD*SH2O(NSOIL)/SMC(NSOIL)  !same liquid fraction as in the bottom layer
    END IF

    DO K = 1, NSOIL
       IF(K == 1) THEN
          DENOM(K) = - ZSOIL (K)
          TEMP1    = - ZSOIL (K+1)
          DDZ(K)   = 2.0 / TEMP1
          DSMDZ(K) = 2.0 * (SMX(K) - SMX(K+1)) / TEMP1
          WFLUX(K) = WDF(K) * DSMDZ(K) + WCND(K) - PDDUM + ETRANI(K) + QSEVA
       ELSE IF (K < NSOIL) THEN
          DENOM(k) = (ZSOIL(K-1) - ZSOIL(K))
          TEMP1    = (ZSOIL(K-1) - ZSOIL(K+1))
          DDZ(K)   = 2.0 / TEMP1
          DSMDZ(K) = 2.0 * (SMX(K) - SMX(K+1)) / TEMP1
          WFLUX(K) = WDF(K  ) * DSMDZ(K  ) + WCND(K  )         &
                   - WDF(K-1) * DSMDZ(K-1) - WCND(K-1) + ETRANI(K)
       ELSE
          DENOM(K) = (ZSOIL(K-1) - ZSOIL(K))
          IF(OPT_RUN == 1 .or. OPT_RUN == 2) THEN
             QDRAIN   = 0.
          END IF
          IF(OPT_RUN == 3) THEN
             QDRAIN   = parameters%SLOPE*WCND(K)
          END IF
          IF(OPT_RUN == 4) THEN
             QDRAIN   = (1.0-FCRMAX)*WCND(K)
          END IF
          IF(OPT_RUN == 5) THEN   !gmm new m-m&f water table dynamics formulation
             TEMP1    = 2.0 * DENOM(K)
             IF(ZWT < ZSOIL(NSOIL)-DENOM(NSOIL))THEN
!gmm interpolate from below, midway to the water table, to the middle of the auxiliary layer below the soil bottom
                SMXBOT = SMX(K) - (SMX(K)-SMXWTD) *  DENOM(K) * 2./ (DENOM(K) + ZSOIL(K) - ZWT)
             ELSE
                SMXBOT = SMXWTD
             ENDIF
             DSMDZ(K) = 2.0 * (SMX(K) - SMXBOT) / TEMP1
             QDRAIN   = WDF(K  ) * DSMDZ(K  ) + WCND(K  )
          END IF   
          WFLUX(K) = -(WDF(K-1)*DSMDZ(K-1))-WCND(K-1)+ETRANI(K) + QDRAIN
       END IF
    END DO

    DO K = 1, NSOIL
       IF(K == 1) THEN
          AI(K)    =   0.0
          BI(K)    =   WDF(K  ) * DDZ(K  ) / DENOM(K)
          CI(K)    = - BI (K)
       ELSE IF (K < NSOIL) THEN
          AI(K)    = - WDF(K-1) * DDZ(K-1) / DENOM(K)
          CI(K)    = - WDF(K  ) * DDZ(K  ) / DENOM(K)
          BI(K)    = - ( AI (K) + CI (K) )
       ELSE
          AI(K)    = - WDF(K-1) * DDZ(K-1) / DENOM(K)
          CI(K)    = 0.0
          BI(K)    = - ( AI (K) + CI (K) )
       END IF
          RHSTT(K) = WFLUX(K) / (-DENOM(K))
    END DO

! ----------------------------------------------------------------------
  END SUBROUTINE SRT

!== begin sstep ====================================================================================

  SUBROUTINE SSTEP (parameters,NSOIL  ,NSNOW  ,DT     ,ZSOIL  ,DZSNSO , & !in
                    SICE   ,ILOC   ,JLOC   ,ZWT            ,                 & !in
                    SH2O   ,SMC    ,AI     ,BI     ,CI     , & !inout
                    RHSTT  ,SMCWTD ,QDRAIN ,DEEPRECH,                                 & !inout
                    WPLUS  )                                   !out

! ----------------------------------------------------------------------
! calculate/update soil moisture content values 
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
    INTEGER,                         INTENT(IN) :: ILOC   !grid index
    INTEGER,                         INTENT(IN) :: JLOC   !grid index
    INTEGER,                         INTENT(IN) :: NSOIL  !
    INTEGER,                         INTENT(IN) :: NSNOW  !
    REAL(r8), INTENT(IN)                            :: DT
    REAL(r8), INTENT(IN)                            :: ZWT
    REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL
    REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: SICE
    REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO ! snow/soil layer thickness [m]

!input and output
    REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: SH2O
    REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: SMC
    REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: AI
    REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: BI
    REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: CI
    REAL(r8), DIMENSION(1:NSOIL), INTENT(INOUT) :: RHSTT
    REAL(r8)                    , INTENT(INOUT) :: SMCWTD
    REAL(r8)                    , INTENT(INOUT) :: QDRAIN
    REAL(r8)                    , INTENT(INOUT) :: DEEPRECH

!output
    REAL(r8), INTENT(OUT)                       :: WPLUS     !saturation excess water (m)

!local
    INTEGER                                 :: K
    REAL(r8), DIMENSION(1:NSOIL)                :: RHSTTIN
    REAL(r8), DIMENSION(1:NSOIL)                :: CIIN
    REAL(r8)                                    :: STOT
    REAL(r8)                                    :: EPORE
    REAL(r8)                                    :: WMINUS
! ----------------------------------------------------------------------
    WPLUS = 0.0

    DO K = 1,NSOIL
       RHSTT (K) =   RHSTT(K) * DT
       AI (K)    =      AI(K) * DT
       BI (K)    = 1. + BI(K) * DT
       CI (K)    =      CI(K) * DT
    END DO

! copy values for input variables before calling rosr12

    DO K = 1,NSOIL
       RHSTTIN(k) = RHSTT(K)
       CIIN(k)    = CI(K)
    END DO

! call ROSR12 to solve the tri-diagonal matrix

    CALL ROSR12 (CI,AI,BI,CIIN,RHSTTIN,RHSTT,1,NSOIL,0)

    DO K = 1,NSOIL
        SH2O(K) = SH2O(K) + CI(K)
    ENDDO

!  excessive water above saturation in a layer is moved to
!  its unsaturated layer like in a bucket

!gmmwith opt_run=5 there is soil moisture below nsoil, to the water table
  IF(OPT_RUN == 5) THEN

!update smcwtd

     IF(ZWT < ZSOIL(NSOIL)-DZSNSO(NSOIL))THEN
!accumulate qdrain to update deep water table and soil moisture later
        DEEPRECH =  DEEPRECH + DT * QDRAIN
     ELSE
        SMCWTD = SMCWTD + DT * QDRAIN  / DZSNSO(NSOIL)
        WPLUS        = MAX((SMCWTD-parameters%SMCMAX(NSOIL)), 0.0) * DZSNSO(NSOIL)
        WMINUS       = MAX((1.E-4-SMCWTD), 0.0) * DZSNSO(NSOIL)

        SMCWTD = MAX( MIN(SMCWTD,parameters%SMCMAX(NSOIL)) , 1.E-4)
        SH2O(NSOIL)    = SH2O(NSOIL) + WPLUS/DZSNSO(NSOIL)

!reduce fluxes at the bottom boundaries accordingly
        QDRAIN = QDRAIN - WPLUS/DT
        DEEPRECH = DEEPRECH - WMINUS
     ENDIF

  ENDIF

    DO K = NSOIL,2,-1
      EPORE        = MAX ( 1.E-4 , ( parameters%SMCMAX(K) - SICE(K) ) )
      WPLUS        = MAX((SH2O(K)-EPORE), 0.0) * DZSNSO(K)
      SH2O(K)      = MIN(EPORE,SH2O(K))
      SH2O(K-1)    = SH2O(K-1) + WPLUS/DZSNSO(K-1)
    END DO

    EPORE        = MAX ( 1.E-4 , ( parameters%SMCMAX(1) - SICE(1) ) )
    WPLUS        = MAX((SH2O(1)-EPORE), 0.0) * DZSNSO(1) 
    SH2O(1)      = MIN(EPORE,SH2O(1))

  END SUBROUTINE SSTEP

!== begin wdfcnd1 ==================================================================================

  SUBROUTINE WDFCND1 (parameters,WDF,WCND,SMC,FCR,ISOIL)
! ----------------------------------------------------------------------
! calculate soil water diffusivity and soil hydraulic conductivity.
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
! input 
  type (noahmp_parameters), intent(in) :: parameters
    REAL(r8),INTENT(IN)  :: SMC
    REAL(r8),INTENT(IN)  :: FCR
    INTEGER,INTENT(IN)  :: ISOIL

! output
    REAL(r8),INTENT(OUT) :: WCND
    REAL(r8),INTENT(OUT) :: WDF

! local
    REAL(r8) :: EXPON
    REAL(r8) :: FACTR
    REAL(r8) :: VKWGT
! ----------------------------------------------------------------------

! soil water diffusivity

    FACTR = MAX(0.01, SMC/parameters%SMCMAX(ISOIL))
    EXPON = parameters%BEXP(ISOIL) + 2.0
    WDF   = parameters%DWSAT(ISOIL) * FACTR ** EXPON
    WDF   = WDF * (1.0 - FCR)

! hydraulic conductivity

    EXPON = 2.0*parameters%BEXP(ISOIL) + 3.0
    WCND  = parameters%DKSAT(ISOIL) * FACTR ** EXPON
    WCND  = WCND * (1.0 - FCR)

  END SUBROUTINE WDFCND1

!== begin wdfcnd2 ==================================================================================

  SUBROUTINE WDFCND2 (parameters,WDF,WCND,SMC,SICE,ISOIL)
! ----------------------------------------------------------------------
! calculate soil water diffusivity and soil hydraulic conductivity.
! ----------------------------------------------------------------------
    IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
    REAL(r8),INTENT(IN)  :: SMC
    REAL(r8),INTENT(IN)  :: SICE
    INTEGER,INTENT(IN)  :: ISOIL

! output
    REAL(r8),INTENT(OUT) :: WCND
    REAL(r8),INTENT(OUT) :: WDF

! local
    REAL(r8) :: EXPON
    REAL(r8) :: FACTR
    REAL(r8) :: VKWGT
! ----------------------------------------------------------------------

! soil water diffusivity

    FACTR = MAX(0.01, SMC/parameters%SMCMAX(ISOIL))
    EXPON = parameters%BEXP(ISOIL) + 2.0
    WDF   = parameters%DWSAT(ISOIL) * FACTR ** EXPON

    IF (SICE > 0.0) THEN
    VKWGT = 1./ (1. + (500.* SICE)**3.)
    WDF   = VKWGT * WDF + (1.-VKWGT)*parameters%DWSAT(ISOIL)*(0.2/parameters%SMCMAX(ISOIL))**EXPON
    END IF

! hydraulic conductivity

    EXPON = 2.0*parameters%BEXP(ISOIL) + 3.0
    WCND  = parameters%DKSAT(ISOIL) * FACTR ** EXPON

  END SUBROUTINE WDFCND2

!== begin groundwater ==============================================================================

  SUBROUTINE GROUNDWATER(parameters,NSNOW  ,NSOIL  ,DT     ,SICE   ,ZSOIL  , & !in
                         STC    ,WCND   ,FCRMAX ,ILOC   ,JLOC   , & !in
                         SH2O   ,ZWT    ,WA     ,WT     ,         & !inout
                         QIN    ,QDIS   )                           !out
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                         INTENT(IN) :: ILOC  !grid index
  INTEGER,                         INTENT(IN) :: JLOC  !grid index
  INTEGER,                         INTENT(IN) :: NSNOW !maximum no. of snow layers
  INTEGER,                         INTENT(IN) :: NSOIL !no. of soil layers
  REAL(r8),                            INTENT(IN) :: DT    !timestep [sec]
  REAL(r8),                            INTENT(IN) :: FCRMAX!maximum FCR (-)
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: SICE  !soil ice content [m3/m3]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL !depth of soil layer-bottom [m]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: WCND  !hydraulic conductivity (m/s)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC   !snow/soil temperature (k)

! input and output
  REAL(r8), DIMENSION(    1:NSOIL), INTENT(INOUT) :: SH2O  !liquid soil water [m3/m3]
  REAL(r8),                         INTENT(INOUT) :: ZWT   !the depth to water table [m]
  REAL(r8),                         INTENT(INOUT) :: WA    !water storage in aquifer [mm]
  REAL(r8),                         INTENT(INOUT) :: WT    !water storage in aquifer 
                                                           !+ saturated soil [mm]
! output
  REAL(r8),                           INTENT(OUT) :: QIN   !groundwater recharge [mm/s]
  REAL(r8),                           INTENT(OUT) :: QDIS  !groundwater discharge [mm/s]

! local
  REAL(r8)                                        :: FFF   !runoff decay factor (m-1)
  REAL(r8)                                        :: RSBMX !baseflow coefficient [mm/s]
  INTEGER                                     :: IZ    !do-loop index
  INTEGER                                     :: IWT   !layer index above water table layer
  REAL(r8),  DIMENSION(    1:NSOIL)               :: DZMM  !layer thickness [mm]
  REAL(r8),  DIMENSION(    1:NSOIL)               :: ZNODE !node depth [m]
  REAL(r8),  DIMENSION(    1:NSOIL)               :: MLIQ  !liquid water mass [kg/m2 or mm]
  REAL(r8),  DIMENSION(    1:NSOIL)               :: EPORE !effective porosity [-]
  REAL(r8),  DIMENSION(    1:NSOIL)               :: HK    !hydraulic conductivity [mm/s]
  REAL(r8),  DIMENSION(    1:NSOIL)               :: SMC   !total soil water  content [m3/m3]
  REAL(KIND=8)                                    :: S_NODE!degree of saturation of IWT layer
!  REAL(r8)                                        :: S_NODE!degree of saturation of IWT layer
  REAL(r8)                                        :: DZSUM !cumulative depth above water table [m]
  REAL(r8)                                        :: SMPFZ !matric potential (frozen effects) [mm]
  REAL(r8)                                        :: KA    !aquifer hydraulic conductivity [mm/s]
  REAL(r8)                                        :: WH_ZWT!water head at water table [mm]
  REAL(r8)                                        :: WH    !water head at layer above ZWT [mm]
  REAL(r8)                                        :: WS    !water used to fill air pore [mm]
  REAL(r8)                                        :: WTSUB !sum of HK*DZMM
  REAL(r8)                                        :: WATMIN!minimum soil vol soil moisture [m3/m3]
  REAL(r8)                                        :: XS    !excessive water above saturation [mm]
  REAL(r8), PARAMETER                             :: ROUS = 0.2    !specific yield [-]
  REAL(r8), PARAMETER                             :: CMIC = 0.20   !microprore content (0.0-1.0)
                                                               !0.0-close to free drainage
! -------------------------------------------------------------
      QDIS      = 0.0
      QIN       = 0.0

! Derive layer-bottom depth in [mm]
!KWM:  Derive layer thickness in mm

      DZMM(1) = -ZSOIL(1)*1.E3
      DO IZ = 2, NSOIL
         DZMM(IZ)  = 1.E3 * (ZSOIL(IZ - 1) - ZSOIL(IZ))
      ENDDO

! Derive node (middle) depth in [m]
!KWM:  Positive number, depth below ground surface in m
      ZNODE(1) = -ZSOIL(1) / 2.
      DO IZ = 2, NSOIL
         ZNODE(IZ)  = -ZSOIL(IZ-1) + 0.5 * (ZSOIL(IZ-1) - ZSOIL(IZ))
      ENDDO

! Convert volumetric soil moisture "sh2o" to mass

      DO IZ = 1, NSOIL
         SMC(IZ)      = SH2O(IZ) + SICE(IZ)
         MLIQ(IZ)     = SH2O(IZ) * DZMM(IZ)
         EPORE(IZ)    = MAX(0.01,parameters%SMCMAX(IZ) - SICE(IZ))
         HK(IZ)       = 1.E3*WCND(IZ)
      ENDDO

! The layer index of the first unsaturated layer,
! i.e., the layer right above the water table

      IWT = NSOIL
      DO IZ = 2,NSOIL
         IF(ZWT   .LE. -ZSOIL(IZ) ) THEN
            IWT = IZ-1
            EXIT
         END IF
      ENDDO

! Groundwater discharge [mm/s]

      FFF   = 6.0
      RSBMX = 5.0

      QDIS = (1.0-FCRMAX)*RSBMX*EXP(-parameters%TIMEAN)*EXP(-FFF*(ZWT-2.0))

! Matric potential at the layer above the water table

      S_NODE = MIN(1.0,SMC(IWT)/parameters%SMCMAX(IWT) )
      S_NODE = MAX(S_NODE,REAL(0.01,KIND=8))
      SMPFZ  = -parameters%PSISAT(IWT)*1000.*S_NODE**(-parameters%BEXP(IWT))   ! m --> mm
      SMPFZ  = MAX(-120000.0,CMIC*SMPFZ)   

! Recharge rate qin to groundwater

      KA  = HK(IWT)

      WH_ZWT  = - ZWT * 1.E3                          !(mm)
      WH      = SMPFZ  - ZNODE(IWT)*1.E3              !(mm)
      QIN     = - KA * (WH_ZWT-WH)  /((ZWT-ZNODE(IWT))*1.E3)
      QIN     = MAX(-10.0/DT,MIN(10./DT,QIN))
     
! Water storage in the aquifer + saturated soil

      WT  = WT + (QIN - QDIS) * DT     !(mm)

      IF(IWT.EQ.NSOIL) THEN
         WA          = WA + (QIN - QDIS) * DT     !(mm)
         WT          = WA
         ZWT         = (-ZSOIL(NSOIL) + 25.) - WA/1000./ROUS      !(m)
         MLIQ(NSOIL) = MLIQ(NSOIL) - QIN * DT        ! [mm]

         MLIQ(NSOIL) = MLIQ(NSOIL) + MAX(0.,(WA - 5000.))
         WA          = MIN(WA, 5000.)
      ELSE
         
         IF (IWT.EQ.NSOIL-1) THEN
            ZWT = -ZSOIL(NSOIL)                   &
                 - (WT-ROUS*1000*25.) / (EPORE(NSOIL))/1000.
         ELSE
            WS = 0.   ! water used to fill soil air pores
            DO IZ = IWT+2,NSOIL
               WS = WS + EPORE(IZ) * DZMM(IZ)
            ENDDO
            ZWT = -ZSOIL(IWT+1)                  &
                  - (WT-ROUS*1000.*25.-WS) /(EPORE(IWT+1))/1000.
         ENDIF

         WTSUB = 0.
         DO IZ = 1, NSOIL
           WTSUB = WTSUB + HK(IZ)*DZMM(IZ)
         END DO

         DO IZ = 1, NSOIL           ! Removing subsurface runoff
         MLIQ(IZ) = MLIQ(IZ) - QDIS*DT*HK(IZ)*DZMM(IZ)/WTSUB
         END DO
      END IF

      ZWT = MAX(1.5,ZWT)

!
! Limit MLIQ to be greater than or equal to watmin.
! Get water needed to bring MLIQ equal WATMIN from lower layer.
!
      WATMIN = 0.01
      DO IZ = 1, NSOIL-1
          IF (MLIQ(IZ) .LT. 0.) THEN
             XS = WATMIN-MLIQ(IZ)
          ELSE
             XS = 0.
          END IF
          MLIQ(IZ  ) = MLIQ(IZ  ) + XS
          MLIQ(IZ+1) = MLIQ(IZ+1) - XS
      END DO

        IZ = NSOIL
        IF (MLIQ(IZ) .LT. WATMIN) THEN
           XS = WATMIN-MLIQ(IZ)
        ELSE
           XS = 0.
        END IF
        MLIQ(IZ) = MLIQ(IZ) + XS
        WA       = WA - XS
        WT       = WT - XS

      DO IZ = 1, NSOIL
        SH2O(IZ)     = MLIQ(IZ) / DZMM(IZ)
      END DO

  END SUBROUTINE GROUNDWATER

!== begin shallowwatertable ========================================================================

  SUBROUTINE SHALLOWWATERTABLE (parameters,NSNOW  ,NSOIL  ,ZSOIL, DT    , & !in
                         DZSNSO ,SMCEQ ,ILOC   ,JLOC         , & !in
                         SMC    ,WTD   ,SMCWTD ,RECH, QDRAIN  )  !inout
! ----------------------------------------------------------------------
!Diagnoses water table depth and computes recharge when the water table is within the resolved soil layers,
!according to the Miguez-Macho&Fan scheme
! ----------------------------------------------------------------------
  IMPLICIT NONE
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  INTEGER,                         INTENT(IN) :: NSNOW !maximum no. of snow layers
  INTEGER,                         INTENT(IN) :: NSOIL !no. of soil layers
  INTEGER,                         INTENT(IN) :: ILOC,JLOC
  REAL(r8),                            INTENT(IN) :: DT
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL !depth of soil layer-bottom [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO ! snow/soil layer thickness [m]
  REAL(r8),  DIMENSION(      1:NSOIL), INTENT(IN) :: SMCEQ  !equilibrium soil water  content [m3/m3]

! input and output
  REAL(r8),  DIMENSION(      1:NSOIL), INTENT(INOUT) :: SMC   !total soil water  content [m3/m3]
  REAL(r8),                         INTENT(INOUT) :: WTD   !the depth to water table [m]
  REAL(r8),                         INTENT(INOUT) :: SMCWTD   !soil moisture between bottom of the soil and the water table [m3/m3]
  REAL(r8),                         INTENT(OUT) :: RECH ! groundwater recharge (net vertical flux across the water table), positive up
  REAL(r8),                         INTENT(INOUT) :: QDRAIN
    
! local
  INTEGER                                     :: IZ    !do-loop index
  INTEGER                                     :: IWTD   !layer index above water table layer
  INTEGER                                     :: KWTD   !layer index where the water table layer is
  REAL(r8)                                        :: WTDOLD
  REAL(r8)                                        :: DZUP
  REAL(r8)                                        :: SMCEQDEEP
  REAL(r8),  DIMENSION(       0:NSOIL)            :: ZSOIL0
! -------------------------------------------------------------


ZSOIL0(1:NSOIL) = ZSOIL(1:NSOIL)
ZSOIL0(0) = 0.         
 
!find the layer where the water table is
     DO IZ=NSOIL,1,-1
        IF(WTD + 1.E-6 < ZSOIL0(IZ)) EXIT
     ENDDO
        IWTD=IZ

        
        KWTD=IWTD+1  !layer where the water table is
        IF(KWTD.LE.NSOIL)THEN    !wtd in the resolved layers
           WTDOLD=WTD
           IF(SMC(KWTD).GT.SMCEQ(KWTD))THEN
        
               IF(SMC(KWTD).EQ.parameters%SMCMAX(KWTD))THEN !wtd went to the layer above
                      WTD=ZSOIL0(IWTD)
                      RECH=-(WTDOLD-WTD) * (parameters%SMCMAX(KWTD)-SMCEQ(KWTD))
                      IWTD=IWTD-1
                      KWTD=KWTD-1
                   IF(KWTD.GE.1)THEN
                      IF(SMC(KWTD).GT.SMCEQ(KWTD))THEN
                      WTDOLD=WTD
                      WTD = MIN( ( SMC(KWTD)*DZSNSO(KWTD) &
                        - SMCEQ(KWTD)*ZSOIL0(IWTD) + parameters%SMCMAX(KWTD)*ZSOIL0(KWTD) ) / &
                        ( parameters%SMCMAX(KWTD)-SMCEQ(KWTD) ), ZSOIL0(IWTD))
                      RECH=RECH-(WTDOLD-WTD) * (parameters%SMCMAX(KWTD)-SMCEQ(KWTD))
                      ENDIF
                   ENDIF
               ELSE  !wtd stays in the layer
                      WTD = MIN( ( SMC(KWTD)*DZSNSO(KWTD) &
                        - SMCEQ(KWTD)*ZSOIL0(IWTD) + parameters%SMCMAX(KWTD)*ZSOIL0(KWTD) ) / &
                        ( parameters%SMCMAX(KWTD)-SMCEQ(KWTD) ), ZSOIL0(IWTD))
                      RECH=-(WTDOLD-WTD) * (parameters%SMCMAX(KWTD)-SMCEQ(KWTD))
               ENDIF
           
           ELSE    !wtd has gone down to the layer below
               WTD=ZSOIL0(KWTD)
               RECH=-(WTDOLD-WTD) * (parameters%SMCMAX(KWTD)-SMCEQ(KWTD))
               KWTD=KWTD+1
               IWTD=IWTD+1
!wtd crossed to the layer below. Now adjust it there
               IF(KWTD.LE.NSOIL)THEN
                   WTDOLD=WTD
                   IF(SMC(KWTD).GT.SMCEQ(KWTD))THEN
                   WTD = MIN( ( SMC(KWTD)*DZSNSO(KWTD) &
                   - SMCEQ(KWTD)*ZSOIL0(IWTD) + parameters%SMCMAX(KWTD)*ZSOIL0(KWTD) ) / &
                       ( parameters%SMCMAX(KWTD)-SMCEQ(KWTD) ) , ZSOIL0(IWTD) )
                   ELSE
                   WTD=ZSOIL0(KWTD)
                   ENDIF
                   RECH = RECH - (WTDOLD-WTD) * &
                                 (parameters%SMCMAX(KWTD)-SMCEQ(KWTD))

                ELSE
                   WTDOLD=WTD
!restore smoi to equilibrium value with water from the ficticious layer below
!                   SMCWTD=SMCWTD-(SMCEQ(NSOIL)-SMC(NSOIL))
!                   QDRAIN = QDRAIN - 1000 * (SMCEQ(NSOIL)-SMC(NSOIL)) * DZSNSO(NSOIL) / DT
!                   SMC(NSOIL)=SMCEQ(NSOIL)
!adjust wtd in the ficticious layer below
                   SMCEQDEEP = parameters%SMCMAX(NSOIL) * ( -parameters%PSISAT(NSOIL) / ( -parameters%PSISAT(NSOIL) - DZSNSO(NSOIL) ) ) ** (1./parameters%BEXP(NSOIL))
                   WTD = MIN( ( SMCWTD*DZSNSO(NSOIL) &
                   - SMCEQDEEP*ZSOIL0(NSOIL) + parameters%SMCMAX(NSOIL)*(ZSOIL0(NSOIL)-DZSNSO(NSOIL)) ) / &
                       ( parameters%SMCMAX(NSOIL)-SMCEQDEEP ) , ZSOIL0(NSOIL) )
                   RECH = RECH - (WTDOLD-WTD) * &
                                 (parameters%SMCMAX(NSOIL)-SMCEQDEEP)
                ENDIF
            
            ENDIF
        ELSEIF(WTD.GE.ZSOIL0(NSOIL)-DZSNSO(NSOIL))THEN
!if wtd was already below the bottom of the resolved soil crust
           WTDOLD=WTD
           SMCEQDEEP = parameters%SMCMAX(NSOIL) * ( -parameters%PSISAT(NSOIL) / ( -parameters%PSISAT(NSOIL) - DZSNSO(NSOIL) ) ) ** (1./parameters%BEXP(NSOIL))
           IF(SMCWTD.GT.SMCEQDEEP)THEN
               WTD = MIN( ( SMCWTD*DZSNSO(NSOIL) &
                 - SMCEQDEEP*ZSOIL0(NSOIL) + parameters%SMCMAX(NSOIL)*(ZSOIL0(NSOIL)-DZSNSO(NSOIL)) ) / &
                     ( parameters%SMCMAX(NSOIL)-SMCEQDEEP ) , ZSOIL0(NSOIL) )
               RECH = -(WTDOLD-WTD) * (parameters%SMCMAX(NSOIL)-SMCEQDEEP)
           ELSE
               RECH = -(WTDOLD-(ZSOIL0(NSOIL)-DZSNSO(NSOIL))) * (parameters%SMCMAX(NSOIL)-SMCEQDEEP)
               WTDOLD=ZSOIL0(NSOIL)-DZSNSO(NSOIL)
!and now even further down
               DZUP=(SMCEQDEEP-SMCWTD)*DZSNSO(NSOIL)/(parameters%SMCMAX(NSOIL)-SMCEQDEEP)
               WTD=WTDOLD-DZUP
               RECH = RECH - (parameters%SMCMAX(NSOIL)-SMCEQDEEP)*DZUP
               SMCWTD=SMCEQDEEP
           ENDIF

         
         ENDIF

IF(IWTD.LT.NSOIL .AND. IWTD.GT.0) THEN
  SMCWTD=parameters%SMCMAX(IWTD)
ELSEIF(IWTD.LT.NSOIL .AND. IWTD.LE.0) THEN
  SMCWTD=parameters%SMCMAX(1)
END IF

END  SUBROUTINE SHALLOWWATERTABLE

! ==================================================================================================
! ********************* end of water subroutines ******************************************
! ==================================================================================================

!== begin carbon ===================================================================================

  SUBROUTINE CARBON (parameters,NSNOW  ,NSOIL  ,VEGTYP ,DT     ,ZSOIL  , & !in
                     DZSNSO ,STC    ,SMC    ,TV     ,TG     ,PSN    , & !in
                     FOLN   ,BTRAN  ,APAR   ,FVEG   ,IGS    , & !in
                     TROOT  ,IST    ,LAT    ,ILOC   ,JLOC   , & !in
                     LFMASS ,RTMASS ,STMASS ,WOOD   ,STBLCP ,FASTCP , & !inout
                     GPP    ,NPP    ,NEE    ,AUTORS ,HETERS ,TOTSC  , & !out
                     TOTLB  ,XLAI   ,XSAI   )                   !out
! ------------------------------------------------------------------------------------------
      IMPLICIT NONE
! ------------------------------------------------------------------------------------------
! inputs (carbon)

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER                        , INTENT(IN) :: ILOC   !grid index
  INTEGER                        , INTENT(IN) :: JLOC   !grid index
  INTEGER                        , INTENT(IN) :: VEGTYP !vegetation type 
  INTEGER                        , INTENT(IN) :: NSNOW  !number of snow layers
  INTEGER                        , INTENT(IN) :: NSOIL  !number of soil layers
  REAL(r8)                           , INTENT(IN) :: LAT    !latitude (radians)
  REAL(r8)                           , INTENT(IN) :: DT     !time step (s)
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL  !depth of layer-bottom from soil surface
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !snow/soil layer thickness [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC    !snow/soil temperature [k]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: SMC    !soil moisture (ice + liq.) [m3/m3]
  REAL(r8)                           , INTENT(IN) :: TV     !vegetation temperature (k)
  REAL(r8)                           , INTENT(IN) :: TG     !ground temperature (k)
  REAL(r8)                           , INTENT(IN) :: FOLN   !foliage nitrogen (%)
  REAL(r8)                           , INTENT(IN) :: BTRAN  !soil water transpiration factor (0 to 1)
  REAL(r8)                           , INTENT(IN) :: PSN    !total leaf photosyn (umolco2/m2/s) [+]
  REAL(r8)                           , INTENT(IN) :: APAR   !PAR by canopy (w/m2)
  REAL(r8)                           , INTENT(IN) :: IGS    !growing season index (0=off, 1=on)
  REAL(r8)                           , INTENT(IN) :: FVEG   !vegetation greenness fraction
  REAL(r8)                           , INTENT(IN) :: TROOT  !root-zone averaged temperature (k)
  INTEGER                        , INTENT(IN) :: IST    !surface type 1->soil; 2->lake

! input & output (carbon)

  REAL(r8)                        , INTENT(INOUT) :: LFMASS !leaf mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: RTMASS !mass of fine roots [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STMASS !stem mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: WOOD   !mass of wood (incl. woody roots) [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STBLCP !stable carbon in deep soil [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: FASTCP !short-lived carbon in shallow soil [g/m2]

! outputs: (carbon)

  REAL(r8)                          , INTENT(OUT) :: GPP    !net instantaneous assimilation [g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: NPP    !net primary productivity [g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: NEE    !net ecosystem exchange [g/m2/s CO2]
  REAL(r8)                          , INTENT(OUT) :: AUTORS !net ecosystem respiration [g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: HETERS !organic respiration [g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: TOTSC  !total soil carbon [g/m2 C]
  REAL(r8)                          , INTENT(OUT) :: TOTLB  !total living carbon ([g/m2 C]
  REAL(r8)                          , INTENT(OUT) :: XLAI   !leaf area index [-]
  REAL(r8)                          , INTENT(OUT) :: XSAI   !stem area index [-]
!  REAL(r8)                          , INTENT(OUT) :: VOCFLX(5) ! voc fluxes [ug C m-2 h-1]

! local variables

  INTEGER :: J         !do-loop index
  REAL(r8)    :: WROOT     !root zone soil water [-]
  REAL(r8)    :: WSTRES    !water stress coeficient [-]  (1. for wilting )
  REAL(r8)    :: LAPM      !leaf area per unit mass [m2/g]
! ------------------------------------------------------------------------------------------

   IF ( ( VEGTYP == parameters%iswater ) .OR. ( VEGTYP == parameters%ISBARREN ) .OR. &
        ( VEGTYP == parameters%ISICE ) .or. (parameters%urban_flag) ) THEN
      XLAI   = 0.
      XSAI   = 0.
      GPP    = 0.
      NPP    = 0.
      NEE    = 0.
      AUTORS = 0.
      HETERS = 0.
      TOTSC  = 0.
      TOTLB  = 0.
      LFMASS = 0.
      RTMASS = 0.
      STMASS = 0.
      WOOD   = 0.
      STBLCP = 0.
      FASTCP = 0.

      RETURN
   END IF

      LAPM       = parameters%SLA / 1000.   ! m2/kg -> m2/g

! water stress

      WSTRES  = 1.- BTRAN

      WROOT  = 0.
      DO J=1,parameters%NROOT
        WROOT = WROOT + SMC(J)/parameters%SMCMAX(J) *  DZSNSO(J) / (-ZSOIL(parameters%NROOT))
      ENDDO

  CALL CO2FLUX (parameters,NSNOW  ,NSOIL  ,VEGTYP ,IGS    ,DT     , & !in
                DZSNSO ,STC    ,PSN    ,TROOT  ,TV     , & !in
                WROOT  ,WSTRES ,FOLN   ,LAPM   ,         & !in
                LAT    ,ILOC   ,JLOC   ,FVEG   ,         & !in
                XLAI   ,XSAI   ,LFMASS ,RTMASS ,STMASS , & !inout
                FASTCP ,STBLCP ,WOOD   ,                 & !inout
                GPP    ,NPP    ,NEE    ,AUTORS ,HETERS , & !out
                TOTSC  ,TOTLB  )                           !out

!   CALL BVOC (parameters,VOCFLX,  VEGTYP,  VEGFAC,   APAR,   TV)
!   CALL CH4

  END SUBROUTINE CARBON

!== begin co2flux ==================================================================================

  SUBROUTINE CO2FLUX (parameters,NSNOW  ,NSOIL  ,VEGTYP ,IGS    ,DT     , & !in
                      DZSNSO ,STC    ,PSN    ,TROOT  ,TV     , & !in
                      WROOT  ,WSTRES ,FOLN   ,LAPM   ,         & !in
                      LAT    ,ILOC   ,JLOC   ,FVEG   ,         & !in
                      XLAI   ,XSAI   ,LFMASS ,RTMASS ,STMASS , & !inout
                      FASTCP ,STBLCP ,WOOD   ,                 & !inout
                      GPP    ,NPP    ,NEE    ,AUTORS ,HETERS , & !out
                      TOTSC  ,TOTLB  )                           !out
! -----------------------------------------------------------------------------------------
! The original code is from RE Dickinson et al.(1998), modifed by Guo-Yue Niu, 2004
! -----------------------------------------------------------------------------------------
  IMPLICIT NONE
! -----------------------------------------------------------------------------------------

! input

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER                        , INTENT(IN) :: ILOC   !grid index
  INTEGER                        , INTENT(IN) :: JLOC   !grid index
  INTEGER                        , INTENT(IN) :: VEGTYP !vegetation physiology type
  INTEGER                        , INTENT(IN) :: NSNOW  !number of snow layers
  INTEGER                        , INTENT(IN) :: NSOIL  !number of soil layers
  REAL(r8)                           , INTENT(IN) :: DT     !time step (s)
  REAL(r8)                           , INTENT(IN) :: LAT    !latitude (radians)
  REAL(r8)                           , INTENT(IN) :: IGS    !growing season index (0=off, 1=on)
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !snow/soil layer thickness [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC    !snow/soil temperature [k]
  REAL(r8)                           , INTENT(IN) :: PSN    !total leaf photosynthesis (umolco2/m2/s)
  REAL(r8)                           , INTENT(IN) :: TROOT  !root-zone averaged temperature (k)
  REAL(r8)                           , INTENT(IN) :: TV     !leaf temperature (k)
  REAL(r8)                           , INTENT(IN) :: WROOT  !root zone soil water
  REAL(r8)                           , INTENT(IN) :: WSTRES !soil water stress
  REAL(r8)                           , INTENT(IN) :: FOLN   !foliage nitrogen (%)
  REAL(r8)                           , INTENT(IN) :: LAPM   !leaf area per unit mass [m2/g]
  REAL(r8)                           , INTENT(IN) :: FVEG   !vegetation greenness fraction

! input and output

  REAL(r8)                        , INTENT(INOUT) :: XLAI   !leaf  area index from leaf carbon [-]
  REAL(r8)                        , INTENT(INOUT) :: XSAI   !stem area index from leaf carbon [-]
  REAL(r8)                        , INTENT(INOUT) :: LFMASS !leaf mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: RTMASS !mass of fine roots [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STMASS !stem mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: FASTCP !short lived carbon [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STBLCP !stable carbon pool [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: WOOD   !mass of wood (incl. woody roots) [g/m2]

! output

  REAL(r8)                          , INTENT(OUT) :: GPP    !net instantaneous assimilation [g/m2/s]
  REAL(r8)                          , INTENT(OUT) :: NPP    !net primary productivity [g/m2]
  REAL(r8)                          , INTENT(OUT) :: NEE    !net ecosystem exchange (autors+heters-gpp)
  REAL(r8)                          , INTENT(OUT) :: AUTORS !net ecosystem resp. (maintance and growth)
  REAL(r8)                          , INTENT(OUT) :: HETERS !organic respiration
  REAL(r8)                          , INTENT(OUT) :: TOTSC  !total soil carbon (g/m2)
  REAL(r8)                          , INTENT(OUT) :: TOTLB  !total living carbon (g/m2)

! local

  REAL(r8)                   :: CFLUX    !carbon flux to atmosphere [g/m2/s]
  REAL(r8)                   :: LFMSMN   !minimum leaf mass [g/m2]
  REAL(r8)                   :: RSWOOD   !wood respiration [g/m2]
  REAL(r8)                   :: RSLEAF   !leaf maintenance respiration per timestep [g/m2]
  REAL(r8)                   :: RSROOT   !fine root respiration per time step [g/m2]
  REAL(r8)                   :: NPPL     !leaf net primary productivity [g/m2/s]
  REAL(r8)                   :: NPPR     !root net primary productivity [g/m2/s]
  REAL(r8)                   :: NPPW     !wood net primary productivity [g/m2/s]
  REAL(r8)                   :: NPPS     !wood net primary productivity [g/m2/s]
  REAL(r8)                   :: DIELF    !death of leaf mass per time step [g/m2]

  REAL(r8)                   :: ADDNPPLF !leaf assimil after resp. losses removed [g/m2]
  REAL(r8)                   :: ADDNPPST !stem assimil after resp. losses removed [g/m2]
  REAL(r8)                   :: CARBFX   !carbon assimilated per model step [g/m2]
  REAL(r8)                   :: GRLEAF   !growth respiration rate for leaf [g/m2/s]
  REAL(r8)                   :: GRROOT   !growth respiration rate for root [g/m2/s]
  REAL(r8)                   :: GRWOOD   !growth respiration rate for wood [g/m2/s]
  REAL(r8)                   :: GRSTEM   !growth respiration rate for stem [g/m2/s]
  REAL(r8)                   :: LEAFPT   !fraction of carbon allocated to leaves [-]
  REAL(r8)                   :: LFDEL    !maximum  leaf mass  available to change [g/m2/s]
  REAL(r8)                   :: LFTOVR   !stem turnover per time step [g/m2]
  REAL(r8)                   :: STTOVR   !stem turnover per time step [g/m2]
  REAL(r8)                   :: WDTOVR   !wood turnover per time step [g/m2]
  REAL(r8)                   :: RSSOIL   !soil respiration per time step [g/m2]
  REAL(r8)                   :: RTTOVR   !root carbon loss per time step by turnover [g/m2]
  REAL(r8)                   :: STABLC   !decay rate of fast carbon to slow carbon [g/m2/s]
  REAL(r8)                   :: WOODF    !calculated wood to root ratio [-]
  REAL(r8)                   :: NONLEF   !fraction of carbon to root and wood [-]
  REAL(r8)                   :: ROOTPT   !fraction of carbon flux to roots [-]
  REAL(r8)                   :: WOODPT   !fraction of carbon flux to wood [-]
  REAL(r8)                   :: STEMPT   !fraction of carbon flux to stem [-]
  REAL(r8)                   :: RESP     !leaf respiration [umol/m2/s]
  REAL(r8)                   :: RSSTEM   !stem respiration [g/m2/s]

  REAL(r8)                   :: FSW      !soil water factor for microbial respiration
  REAL(r8)                   :: FST      !soil temperature factor for microbial respiration
  REAL(r8)                   :: FNF      !foliage nitrogen adjustemt to respiration (<= 1)
  REAL(r8)                   :: TF       !temperature factor
  REAL(r8)                   :: RF       !respiration reduction factor (<= 1)
  REAL(r8)                   :: STDEL
  REAL(r8)                   :: STMSMN
  REAL(r8)                   :: SAPM     !stem area per unit mass (m2/g)
  REAL(r8)                   :: DIEST
! -------------------------- constants -------------------------------
  REAL(r8)                   :: BF       !parameter for present wood allocation [-]
  REAL(r8)                   :: RSWOODC  !wood respiration coeficient [1/s]
  REAL(r8)                   :: STOVRC   !stem turnover coefficient [1/s]
  REAL(r8)                   :: RSDRYC   !degree of drying that reduces soil respiration [-]
  REAL(r8)                   :: RTOVRC   !root turnover coefficient [1/s]
  REAL(r8)                   :: WSTRC    !water stress coeficient [-]
  REAL(r8)                   :: LAIMIN   !minimum leaf area index [m2/m2]
  REAL(r8)                   :: XSAMIN   !minimum leaf area index [m2/m2]
  REAL(r8)                   :: SC
  REAL(r8)                   :: SD
  REAL(r8)                   :: VEGFRAC

! Respiration as a function of temperature

  real :: r,x
          r(x) = exp(0.08*(x-298.16))
! ---------------------------------------------------------------------------------

! constants
    RTOVRC  = 2.0E-8        !original was 2.0e-8
    RSDRYC  = 40.0          !original was 40.0
    RSWOODC = 3.0E-10       !
    BF      = 0.90          !original was 0.90   ! carbon to roots
    WSTRC   = 100.0
    LAIMIN  = 0.05   
    XSAMIN  = 0.05     ! MB: change to prevent vegetation from not growing back in spring

    SAPM    = 3.*0.001      ! m2/kg -->m2/g
    LFMSMN  = laimin/lapm
    STMSMN  = xsamin/sapm
! ---------------------------------------------------------------------------------

! respiration

     IF(IGS .EQ. 0.) THEN
       RF = 0.5
     ELSE
       RF = 1.0
     ENDIF
            
     FNF     = MIN( FOLN/MAX(1.E-06,parameters%FOLNMX), 1.0 )
     TF      = parameters%ARM**( (TV-298.16)/10. )
     RESP    = parameters%RMF25 * TF * FNF * XLAI * RF * (1.-WSTRES) ! umol/m2/s
     RSLEAF  = MIN((LFMASS-LFMSMN)/DT,RESP*12.e-6)                         ! g/m2/s
     
     RSROOT  = parameters%RMR25*(RTMASS*1E-3)*TF *RF* 12.e-6         ! g/m2/s
     RSSTEM  = parameters%RMS25*((STMASS-STMSMN)*1E-3)*TF *RF* 12.e-6         ! g/m2/s
     RSWOOD  = RSWOODC * R(TV) * WOOD*parameters%WDPOOL

! carbon assimilation
! 1 mole -> 12 g carbon or 44 g CO2; 1 umol -> 12.e-6 g carbon;

     CARBFX  = PSN * 12.e-6              ! umol co2 /m2/ s -> g/m2/s carbon

! fraction of carbon into leaf versus nonleaf

     LEAFPT = EXP(0.01*(1.-EXP(0.75*XLAI))*XLAI)
     IF(VEGTYP == parameters%EBLFOREST) LEAFPT = EXP(0.01*(1.-EXP(0.50*XLAI))*XLAI)

     NONLEF = 1.0 - LEAFPT
     STEMPT = XLAI/10.0*LEAFPT
     LEAFPT = LEAFPT - STEMPT

!  fraction of carbon into wood versus root

     IF(WOOD > 1.e-6) THEN
        WOODF = (1.-EXP(-BF*(parameters%WRRAT*RTMASS/WOOD))/BF)*parameters%WDPOOL
     ELSE
        WOODF = parameters%WDPOOL
     ENDIF

     ROOTPT = NONLEF*(1.-WOODF)
     WOODPT = NONLEF*WOODF

! leaf and root turnover per time step

     LFTOVR = parameters%LTOVRC*5.E-7*LFMASS
     STTOVR = parameters%LTOVRC*5.E-7*STMASS
     RTTOVR = RTOVRC*RTMASS
     WDTOVR = 9.5E-10*WOOD

! seasonal leaf die rate dependent on temp and water stress
! water stress is set to 1 at permanent wilting point

     SC  = EXP(-0.3*MAX(0.,TV-parameters%TDLEF)) * (LFMASS/120.) 
     SD  = EXP((WSTRES-1.)*WSTRC)
     DIELF = LFMASS*1.E-6*(parameters%DILEFW * SD + parameters%DILEFC*SC)
     DIEST = STMASS*1.E-6*(parameters%DILEFW * SD + parameters%DILEFC*SC)

! calculate growth respiration for leaf, rtmass and wood

     GRLEAF = MAX(0.0,parameters%FRAGR*(LEAFPT*CARBFX - RSLEAF))
     GRSTEM = MAX(0.0,parameters%FRAGR*(STEMPT*CARBFX - RSSTEM))
     GRROOT = MAX(0.0,parameters%FRAGR*(ROOTPT*CARBFX - RSROOT))
     GRWOOD = MAX(0.0,parameters%FRAGR*(WOODPT*CARBFX - RSWOOD))

! Impose lower T limit for photosynthesis

     ADDNPPLF = MAX(0.,LEAFPT*CARBFX - GRLEAF-RSLEAF)
     ADDNPPST = MAX(0.,STEMPT*CARBFX - GRSTEM-RSSTEM)
!     ADDNPPLF = LEAFPT*CARBFX - GRLEAF-RSLEAF  ! MB: test Kjetil 
!     ADDNPPST = STEMPT*CARBFX - GRSTEM-RSSTEM  ! MB: test Kjetil 
     IF(TV.LT.parameters%TMIN) ADDNPPLF =0.
     IF(TV.LT.parameters%TMIN) ADDNPPST =0.

! update leaf, root, and wood carbon
! avoid reducing leaf mass below its minimum value but conserve mass

     LFDEL = (LFMASS - LFMSMN)/DT
     STDEL = (STMASS - STMSMN)/DT
     DIELF = MIN(DIELF,LFDEL+ADDNPPLF-LFTOVR)
     DIEST = MIN(DIEST,STDEL+ADDNPPST-STTOVR)

! net primary productivities

     NPPL   = MAX(ADDNPPLF,-LFDEL)
     NPPS   = MAX(ADDNPPST,-STDEL)
     NPPR   = ROOTPT*CARBFX - RSROOT - GRROOT
     NPPW   = WOODPT*CARBFX - RSWOOD - GRWOOD

! masses of plant components

     LFMASS = LFMASS + (NPPL-LFTOVR-DIELF)*DT
     STMASS = STMASS + (NPPS-STTOVR-DIEST)*DT   ! g/m2
     RTMASS = RTMASS + (NPPR-RTTOVR)      *DT

     IF(RTMASS.LT.0.0) THEN
           RTTOVR = NPPR
           RTMASS = 0.0
     ENDIF
     WOOD = (WOOD+(NPPW-WDTOVR)*DT)*parameters%WDPOOL

! soil carbon budgets

     FASTCP = FASTCP + (RTTOVR+LFTOVR+STTOVR+WDTOVR+DIELF+DIEST)*DT  ! MB: add DIEST v3.7

     FST = 2.0**( (STC(1)-283.16)/10. )
     FSW = WROOT / (0.20+WROOT) * 0.23 / (0.23+WROOT)
     RSSOIL = FSW * FST * parameters%MRP* MAX(0.,FASTCP*1.E-3)*12.E-6

     STABLC = 0.1*RSSOIL
     FASTCP = FASTCP - (RSSOIL + STABLC)*DT
     STBLCP = STBLCP + STABLC*DT

!  total carbon flux

     CFLUX  = - CARBFX + RSLEAF + RSROOT + RSWOOD + RSSTEM &  ! MB: add RSSTEM,GRSTEM,0.9*RSSOIL v3.7
          + 0.9*RSSOIL + GRLEAF + GRROOT + GRWOOD + GRSTEM    ! g/m2/s

! for outputs

     GPP    = CARBFX                                             !g/m2/s C
     NPP    = NPPL + NPPW + NPPR +NPPS                           !g/m2/s C
     AUTORS = RSROOT + RSWOOD  + RSLEAF + RSSTEM + &             !g/m2/s C  MB: add RSSTEM, GRSTEM v3.7
              GRLEAF + GRROOT + GRWOOD + GRSTEM                  !g/m2/s C  MB: add 0.9* v3.7
     HETERS = 0.9*RSSOIL                                         !g/m2/s C
     NEE    = (AUTORS + HETERS - GPP)*44./12.                    !g/m2/s CO2
     TOTSC  = FASTCP + STBLCP                                    !g/m2   C
     TOTLB  = LFMASS + RTMASS +STMASS + WOOD                     !g/m2   C  MB: add STMASS v3.7

! leaf area index and stem area index

     XLAI    = MAX(LFMASS*LAPM,LAIMIN)
     XSAI    = MAX(STMASS*SAPM,XSAMIN)
    
  END SUBROUTINE CO2FLUX

!== begin carbon_crop ==============================================================================

 SUBROUTINE CARBON_CROP (parameters,NSNOW  ,NSOIL  ,VEGTYP ,DT     ,ZSOIL  ,JULIAN , & !in
                            DZSNSO ,STC    ,SMC    ,TV     ,PSN    ,FOLN   ,BTRAN  , & !in
                            SOLDN  ,T2M    ,                                         & !in
                            LFMASS ,RTMASS ,STMASS ,WOOD   ,STBLCP ,FASTCP ,GRAIN  , & !inout
			    XLAI   ,XSAI   ,GDD    ,                                 & !inout
                            GPP    ,NPP    ,NEE    ,AUTORS ,HETERS ,TOTSC  ,TOTLB    ) !out
! ------------------------------------------------------------------------------------------
! Initial crop version created by Xing Liu
! Initial crop version added by Barlage v3.8

! ------------------------------------------------------------------------------------------
      IMPLICIT NONE
! ------------------------------------------------------------------------------------------
! inputs (carbon)

  type (noahmp_parameters), intent(in) :: parameters
  INTEGER                        , INTENT(IN) :: NSNOW  !number of snow layers
  INTEGER                        , INTENT(IN) :: NSOIL  !number of soil layers
  INTEGER                        , INTENT(IN) :: VEGTYP !vegetation type 
  REAL(r8)                           , INTENT(IN) :: DT     !time step (s)
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: ZSOIL  !depth of layer-bottomfrom soil surface
  REAL(r8)                           , INTENT(IN) :: JULIAN !Julian day of year(fractional) ( 0 <= JULIAN < YEARLEN )
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: DZSNSO !snow/soil layerthickness [m]
  REAL(r8), DIMENSION(-NSNOW+1:NSOIL), INTENT(IN) :: STC    !snow/soil temperature[k]
  REAL(r8), DIMENSION(       1:NSOIL), INTENT(IN) :: SMC    !soil moisture (ice +liq.) [m3/m3]
  REAL(r8)                           , INTENT(IN) :: TV     !vegetation temperature(k)
  REAL(r8)                           , INTENT(IN) :: PSN    !total leaf photosyn(umolco2/m2/s) [+]
  REAL(r8)                           , INTENT(IN) :: FOLN   !foliage nitrogen (%)
  REAL(r8)                           , INTENT(IN) :: BTRAN  !soil watertranspiration factor (0 to 1)
  REAL(r8)                           , INTENT(IN) :: SOLDN  !Downward solar radiation
  REAL(r8)                           , INTENT(IN) :: T2M    !air temperature

! input & output (carbon)

  REAL(r8)                        , INTENT(INOUT) :: LFMASS !leaf mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: RTMASS !mass of fine roots[g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STMASS !stem mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: WOOD   !mass of wood (incl.woody roots) [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STBLCP !stable carbon in deepsoil [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: FASTCP !short-lived carbon inshallow soil [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: GRAIN  !mass of GRAIN [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: XLAI   !leaf area index [-]
  REAL(r8)                        , INTENT(INOUT) :: XSAI   !stem area index [-]
  REAL(r8)                        , INTENT(INOUT) :: GDD    !growing degree days

! outout
  REAL(r8)                          , INTENT(OUT) :: GPP    !net instantaneous assimilation [g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: NPP    !net primary productivity [g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: NEE    !net ecosystem exchange[g/m2/s CO2]
  REAL(r8)                          , INTENT(OUT) :: AUTORS !net ecosystem respiration [g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: HETERS !organic respiration[g/m2/s C]
  REAL(r8)                          , INTENT(OUT) :: TOTSC  !total soil carbon [g/m2C]
  REAL(r8)                          , INTENT(OUT) :: TOTLB  !total living carbon ([g/m2 C]

! local variables

  INTEGER :: J         !do-loop index
  REAL(r8)    :: WROOT     !root zone soil water [-]
  REAL(r8)    :: WSTRES    !water stress coeficient [-]  (1. for wilting )
  INTEGER :: IPA       !Planting index
  INTEGER :: IHA       !Havestindex(0=on,1=off)
  INTEGER :: PGS       !Plant growth stage

  REAL(r8)    :: PSNCROP 

! ------------------------------------------------------------------------------------------
   IF ( ( VEGTYP == parameters%iswater ) .OR. ( VEGTYP == parameters%ISBARREN ) .OR. &
        ( VEGTYP == parameters%ISICE ) .or. (parameters%urban_flag) ) THEN
      XLAI   = 0.
      XSAI   = 0.
      GPP    = 0.
      NPP    = 0.
      NEE    = 0.
      AUTORS = 0.
      HETERS = 0.
      TOTSC  = 0.
      TOTLB  = 0.
      LFMASS = 0.
      RTMASS = 0.
      STMASS = 0.
      WOOD   = 0.
      STBLCP = 0.
      FASTCP = 0.
      GRAIN  = 0.
      RETURN
   END IF

! water stress


   WSTRES  = 1.- BTRAN

   WROOT  = 0.
   DO J=1,parameters%NROOT
     WROOT = WROOT + SMC(J)/parameters%SMCMAX(J) *  DZSNSO(J) / (-ZSOIL(parameters%NROOT))
   ENDDO

   CALL PSN_CROP     ( parameters,                           & !in
                       SOLDN,   XLAI,    T2M,                & !in 
                       PSNCROP                             )   !out

   CALL GROWING_GDD  (parameters,                           & !in
                      T2M ,   DT,  JULIAN,                  & !in
                      GDD ,                                 & !inout 
                      IPA ,  IHA,     PGS)                    !out                        

   CALL CO2FLUX_CROP (parameters,                              & !in
                      DT     ,STC(1) ,PSN    ,TV     ,WROOT  ,WSTRES ,FOLN   , & !in
                      IPA    ,IHA    ,PGS    ,                                 & !in XING
                      XLAI   ,XSAI   ,LFMASS ,RTMASS ,STMASS ,                 & !inout
                      FASTCP ,STBLCP ,WOOD   ,GRAIN  ,GDD    ,                 & !inout
                      GPP    ,NPP    ,NEE    ,AUTORS ,HETERS ,                 & !out
                      TOTSC  ,TOTLB  )                                           !out

  END SUBROUTINE CARBON_CROP

!== begin co2flux_crop =============================================================================

  SUBROUTINE CO2FLUX_CROP (parameters,                                              & !in
                           DT     ,STC    ,PSN    ,TV     ,WROOT  ,WSTRES ,FOLN   , & !in
                           IPA    ,IHA    ,PGS    ,                                 & !in XING
                           XLAI   ,XSAI   ,LFMASS ,RTMASS ,STMASS ,                 & !inout
                           FASTCP ,STBLCP ,WOOD   ,GRAIN  ,GDD,                     & !inout
                           GPP    ,NPP    ,NEE    ,AUTORS ,HETERS ,                 & !out
                           TOTSC  ,TOTLB  )                                           !out
! -----------------------------------------------------------------------------------------
! The original code from RE Dickinson et al.(1998) and Guo-Yue Niu(2004),
! modified by Xing Liu, 2014.
! 
! -----------------------------------------------------------------------------------------
  IMPLICIT NONE
! -----------------------------------------------------------------------------------------

! input

  type (noahmp_parameters), intent(in) :: parameters
  REAL(r8)                           , INTENT(IN) :: DT     !time step (s)
  REAL(r8)                           , INTENT(IN) :: STC    !soil temperature[k]
  REAL(r8)                           , INTENT(IN) :: PSN    !total leaf photosynthesis (umolco2/m2/s)
  REAL(r8)                           , INTENT(IN) :: TV     !leaf temperature (k)
  REAL(r8)                           , INTENT(IN) :: WROOT  !root zone soil water
  REAL(r8)                           , INTENT(IN) :: WSTRES !soil water stress
  REAL(r8)                           , INTENT(IN) :: FOLN   !foliage nitrogen (%)
  INTEGER                        , INTENT(IN) :: IPA
  INTEGER                        , INTENT(IN) :: IHA
  INTEGER                        , INTENT(IN) :: PGS

! input and output

  REAL(r8)                        , INTENT(INOUT) :: XLAI   !leaf  area index from leaf carbon [-]
  REAL(r8)                        , INTENT(INOUT) :: XSAI   !stem area index from leaf carbon [-]
  REAL(r8)                        , INTENT(INOUT) :: LFMASS !leaf mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: RTMASS !mass of fine roots [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STMASS !stem mass [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: FASTCP !short lived carbon [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: STBLCP !stable carbon pool [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: WOOD   !mass of wood (incl. woody roots) [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: GRAIN  !mass of grain (XING) [g/m2]
  REAL(r8)                        , INTENT(INOUT) :: GDD    !growing degree days (XING)

! output

  REAL(r8)                          , INTENT(OUT) :: GPP    !net instantaneous assimilation [g/m2/s]
  REAL(r8)                          , INTENT(OUT) :: NPP    !net primary productivity [g/m2]
  REAL(r8)                          , INTENT(OUT) :: NEE    !net ecosystem exchange (autors+heters-gpp)
  REAL(r8)                          , INTENT(OUT) :: AUTORS !net ecosystem resp. (maintance and growth)
  REAL(r8)                          , INTENT(OUT) :: HETERS !organic respiration
  REAL(r8)                          , INTENT(OUT) :: TOTSC  !total soil carbon (g/m2)
  REAL(r8)                          , INTENT(OUT) :: TOTLB  !total living carbon (g/m2)

! local

  REAL(r8)                   :: CFLUX    !carbon flux to atmosphere [g/m2/s]
  REAL(r8)                   :: LFMSMN   !minimum leaf mass [g/m2]
  REAL(r8)                   :: RSWOOD   !wood respiration [g/m2]
  REAL(r8)                   :: RSLEAF   !leaf maintenance respiration per timestep[g/m2]
  REAL(r8)                   :: RSROOT   !fine root respiration per time step [g/m2]
  REAL(r8)                   :: RSGRAIN  !grain respiration [g/m2]  
  REAL(r8)                   :: NPPL     !leaf net primary productivity [g/m2/s]
  REAL(r8)                   :: NPPR     !root net primary productivity [g/m2/s]
  REAL(r8)                   :: NPPW     !wood net primary productivity [g/m2/s]
  REAL(r8)                   :: NPPS     !wood net primary productivity [g/m2/s]
  REAL(r8)                   :: NPPG     !grain net primary productivity [g/m2/s] 
  REAL(r8)                   :: DIELF    !death of leaf mass per time step [g/m2]

  REAL(r8)                   :: ADDNPPLF !leaf assimil after resp. losses removed[g/m2]
  REAL(r8)                   :: ADDNPPST !stem assimil after resp. losses removed[g/m2]
  REAL(r8)                   :: CARBFX   !carbon assimilated per model step [g/m2]
  REAL(r8)                   :: CBHYDRAFX!carbonhydrate assimilated per model step [g/m2]
  REAL(r8)                   :: GRLEAF   !growth respiration rate for leaf [g/m2/s]
  REAL(r8)                   :: GRROOT   !growth respiration rate for root [g/m2/s]
  REAL(r8)                   :: GRWOOD   !growth respiration rate for wood [g/m2/s]
  REAL(r8)                   :: GRSTEM   !growth respiration rate for stem [g/m2/s]
  REAL(r8)                   :: GRGRAIN   !growth respiration rate for stem [g/m2/s]
  REAL(r8)                   :: LEAFPT   !fraction of carbon allocated to leaves [-]
  REAL(r8)                   :: LFDEL    !maximum  leaf mass  available to change[g/m2/s]
  REAL(r8)                   :: LFTOVR   !stem turnover per time step [g/m2]
  REAL(r8)                   :: STTOVR   !stem turnover per time step [g/m2]
  REAL(r8)                   :: WDTOVR   !wood turnover per time step [g/m2]
  REAL(r8)                   :: GRTOVR   !grainturnover per time step [g/m2]
  REAL(r8)                   :: RSSOIL   !soil respiration per time step [g/m2]
  REAL(r8)                   :: RTTOVR   !root carbon loss per time step by turnover[g/m2]
  REAL(r8)                   :: STABLC   !decay rate of fast carbon to slow carbon[g/m2/s]
  REAL(r8)                   :: WOODF    !calculated wood to root ratio [-]
  REAL(r8)                   :: NONLEF   !fraction of carbon to root and wood [-]
  REAL(r8)                   :: RESP     !leaf respiration [umol/m2/s]
  REAL(r8)                   :: RSSTEM   !stem respiration [g/m2/s]

  REAL(r8)                   :: FSW      !soil water factor for microbial respiration
  REAL(r8)                   :: FST      !soil temperature factor for microbialrespiration
  REAL(r8)                   :: FNF      !foliage nitrogen adjustemt to respiration(<= 1)
  REAL(r8)                   :: TF       !temperature factor
  REAL(r8)                   :: STDEL
  REAL(r8)                   :: STMSMN
  REAL(r8)                   :: SAPM     !stem area per unit mass (m2/g)
  REAL(r8)                   :: DIEST
! -------------------------- constants -------------------------------
  REAL(r8)                   :: BF       !parameter for present wood allocation [-]
  REAL(r8)                   :: RSWOODC  !wood respiration coeficient [1/s]
  REAL(r8)                   :: STOVRC   !stem turnover coefficient [1/s]
  REAL(r8)                   :: RSDRYC   !degree of drying that reduces soilrespiration [-]
  REAL(r8)                   :: RTOVRC   !root turnover coefficient [1/s]
  REAL(r8)                   :: WSTRC    !water stress coeficient [-]
  REAL(r8)                   :: LAIMIN   !minimum leaf area index [m2/m2]
  REAL(r8)                   :: XSAMIN   !minimum leaf area index [m2/m2]
  REAL(r8)                   :: SC
  REAL(r8)                   :: SD
  REAL(r8)                   :: VEGFRAC
  REAL(r8)                   :: TEMP

! Respiration as a function of temperature

  real :: r,x
          r(x) = exp(0.08*(x-298.16))
! ---------------------------------------------------------------------------------

! constants
    RSDRYC  = 40.0          !original was 40.0
    RSWOODC = 3.0E-10       !
    BF      = 0.90          !original was 0.90   ! carbon to roots
    WSTRC   = 100.0
    LAIMIN  = 0.05
    XSAMIN  = 0.01

    SAPM    = 3.*0.001      ! m2/kg -->m2/g
    LFMSMN  = laimin/0.15
    STMSMN  = xsamin/sapm
! ---------------------------------------------------------------------------------

! carbon assimilation
! 1 mole -> 12 g carbon or 44 g CO2 or 30 g CH20

     CARBFX     = PSN*12.e-6*IPA   !umol co2 /m2/ s -> g/m2/s C
     CBHYDRAFX  = PSN*30.e-6*IPA

! mainteinance respiration
     FNF     = MIN( FOLN/MAX(1.E-06,parameters%FOLN_MX), 1.0 )
     TF      = parameters%Q10MR**( (TV-298.16)/10. )
     RESP    = parameters%LFMR25 * TF * FNF * XLAI  * (1.-WSTRES)  ! umol/m2/s
     RSLEAF  = MIN(LFMASS/DT,RESP*30.e-6)                       ! g/m2/s
     RSROOT  = parameters%RTMR25*(RTMASS*1E-3)*TF * 30.e-6         ! g/m2/s
     RSSTEM  = parameters%STMR25*(STMASS*1E-3)*TF * 30.e-6         ! g/m2/s
     RSGRAIN = parameters%GRAINMR25*(GRAIN*1E-3)*TF * 30.e-6       ! g/m2/s

! calculate growth respiration for leaf, rtmass and grain

     GRLEAF  = MAX(0.0,parameters%FRA_GR*(parameters%LFPT(PGS)*CBHYDRAFX  - RSLEAF))
     GRSTEM  = MAX(0.0,parameters%FRA_GR*(parameters%STPT(PGS)*CBHYDRAFX  - RSSTEM))
     GRROOT  = MAX(0.0,parameters%FRA_GR*(parameters%RTPT(PGS)*CBHYDRAFX  - RSROOT))
     GRGRAIN = MAX(0.0,parameters%FRA_GR*(parameters%GRAINPT(PGS)*CBHYDRAFX  - RSGRAIN))

! leaf turnover, stem turnover, root turnover and leaf death caused by soil
! water and soil temperature stress

     LFTOVR  = parameters%LF_OVRC(PGS)*1.E-6*LFMASS
     RTTOVR  = parameters%RT_OVRC(PGS)*1.E-6*RTMASS
     STTOVR  = parameters%ST_OVRC(PGS)*1.E-6*STMASS
     SC  = EXP(-0.3*MAX(0.,TV-parameters%LEFREEZ)) * (LFMASS/120.)
     SD  = EXP((WSTRES-1.)*WSTRC)
     DIELF = LFMASS*1.E-6*(parameters%DILE_FW(PGS) * SD + parameters%DILE_FC(PGS)*SC)

! Allocation of CBHYDRAFX to leaf, stem, root and grain at each growth stage


     ADDNPPLF    = MAX(0.,parameters%LFPT(PGS)*CBHYDRAFX - GRLEAF-RSLEAF)
     ADDNPPST    = MAX(0.,parameters%STPT(PGS)*CBHYDRAFX - GRSTEM-RSSTEM)
    

! avoid reducing leaf mass below its minimum value but conserve mass

     LFDEL = (LFMASS - LFMSMN)/DT
     STDEL = (STMASS - STMSMN)/DT
     DIELF = MIN(DIELF,LFDEL+ADDNPPLF-LFTOVR)

! net primary productivities

     NPPL   = MAX(ADDNPPLF,-LFDEL)
     NPPS   = MAX(ADDNPPST,-STDEL)
     NPPR   = parameters%RTPT(PGS)*CBHYDRAFX - RSROOT - GRROOT
     NPPG  =  parameters%GRAINPT(PGS)*CBHYDRAFX - RSGRAIN - GRGRAIN

! masses of plant components
  
     LFMASS = LFMASS + (NPPL-LFTOVR-DIELF)*DT
     STMASS = STMASS + (NPPS-STTOVR)*DT   ! g/m2
     RTMASS = RTMASS + (NPPR-RTTOVR)*DT
     GRAIN =  GRAIN + NPPG*DT 

     GPP = CBHYDRAFX* 0.4 !!g/m2/s C  0.4=12/30, CH20 to C

     IF(PGS==6) THEN
       STMASS = STMASS - STMASS*(0.00005)
       RTMASS = RTMASS - RTMASS*(0.0005)
       GRAIN  = GRAIN + STMASS*(0.00005) + RTMASS*(0.0005) 
     END IF
    
     IF(RTMASS.LT.0.0) THEN
       RTTOVR = NPPR
       RTMASS = 0.0
     ENDIF

     IF(GRAIN.LT.0.0) THEN
       GRAIN = 0.0
     ENDIF

 ! soil carbon budgets

     IF(PGS == 1 .OR. PGS == 2 .OR. PGS == 8) THEN
       FASTCP=1000
     ELSE
       FASTCP = FASTCP + (RTTOVR+LFTOVR+STTOVR+DIELF)*DT 
     END IF
     FST = 2.0**( (STC-283.16)/10. )
     FSW = WROOT / (0.20+WROOT) * 0.23 / (0.23+WROOT)
     RSSOIL = FSW * FST * parameters%MRP* MAX(0.,FASTCP*1.E-3)*12.E-6

     STABLC = 0.1*RSSOIL
     FASTCP = FASTCP - (RSSOIL + STABLC)*DT
     STBLCP = STBLCP + STABLC*DT

!  total carbon flux

     CFLUX  = - CARBFX + RSLEAF + RSROOT  + RSSTEM &
              + RSSOIL + GRLEAF + GRROOT                  ! g/m2/s 0.4=12/30, CH20 to C

! for outputs
                                                                 !g/m2/s C

     NPP   = (NPPL + NPPS+ NPPR +NPPG)*0.4      !!g/m2/s C  0.4=12/30, CH20 to C
 
  
     AUTORS = RSROOT + RSGRAIN  + RSLEAF +  &                     !g/m2/s C
              GRLEAF + GRROOT + GRGRAIN                           !g/m2/s C

     HETERS = RSSOIL                                             !g/m2/s C
     NEE    = (AUTORS + HETERS - GPP)*44./30.                    !g/m2/s CO2
     TOTSC  = FASTCP + STBLCP                                    !g/m2   C

     TOTLB  = LFMASS + RTMASS + GRAIN         

! leaf area index and stem area index
  
     XLAI    = MAX(LFMASS*parameters%BIO2LAI,LAIMIN)
     XSAI    = MAX(STMASS*SAPM,XSAMIN)

   
!After harversting
     IF(PGS == 8 ) THEN
       LFMASS = 0.62
       STMASS = 0
       TOTLB  = 0
       GPP    = 0
       NPP    = 0
       GRAIN  = 0
       AUTORS = 0
       NEE    = 0
     END IF

    IF(PGS == 1 .OR. PGS == 2 .OR. PGS == 8) THEN
     XLAI   = 0.05
     XSAI   = 0.05
     LFMASS = LFMSMN
     STMASS = STMSMN
     RTMASS = 0
    END IF 
    
END SUBROUTINE CO2FLUX_CROP

!== begin growing_gdd ==============================================================================

  SUBROUTINE GROWING_GDD (parameters,                         & !in
                          T2M ,   DT, JULIAN,                 & !in
                          GDD ,                               & !inout 
                          IPA,   IHA,     PGS)                  !out  
!===================================================================================================

! input

  type (noahmp_parameters), intent(in) :: parameters
   REAL(r8)                     , INTENT(IN)        :: T2M     !Air temperature
   REAL(r8)                     , INTENT(IN)        :: DT      !time step (s)
   REAL(r8)                     , INTENT(IN)        :: JULIAN  !Julian day of year (fractional) ( 0 <= JULIAN < YEARLEN )

! input and output

   REAL(r8)                     , INTENT(INOUT)     :: GDD     !growing degress days

! output

   INTEGER                  , INTENT(OUT)       :: IPA     !Planting index index(0=off, 1=on)
   INTEGER                  , INTENT(OUT)       :: IHA     !Havestindex(0=on,1=off) 
   INTEGER                  , INTENT(OUT)       :: PGS     !Plant growth stage(1=S1,2=S2,3=S3)

!local 

   REAL(r8)                                         :: GDDDAY    !gap bewtween GDD and GDD8
   REAL(r8)                                         :: DAYOFS2   !DAYS in stage2
   REAL(r8)                                         :: TDIFF     !temperature difference for growing degree days calculation
   REAL(r8)                                         :: TC

   TC = T2M - 273.15

!Havestindex(0=on,1=off) 

   IPA = 1
   IHA = 1

!turn on/off the planting 
 
   IF(JULIAN < parameters%PLTDAY)  IPA = 0

!turn on/off the harvesting
    IF(JULIAN >= parameters%HSDAY) IHA = 0
   
!Calculate the growing degree days
   
    IF(TC <  parameters%GDDTBASE) THEN
      TDIFF = 0.0
    ELSEIF(TC >= parameters%GDDTCUT) THEN
      TDIFF = parameters%GDDTCUT - parameters%GDDTBASE
    ELSE
      TDIFF = TC - parameters%GDDTBASE
    END IF

    GDD     = (GDD + TDIFF) * IPA * IHA

    GDDDAY  = GDD / (86400.0 / DT)

   ! Decide corn growth stage, based on Hybrid-Maize 
   !   PGS = 1 : Before planting
   !   PGS = 2 : from tassel initiation to silking
   !   PGS = 3 : from silking to effective grain filling
   !   PGS = 4 : from effective grain filling to pysiological maturity 
   !   PGS = 5 : GDDM=1389
   !   PGS = 6 :
   !   PGS = 7 :
   !   PGS = 8 :
   !  GDDM = 1389
   !  GDDM = 1555
   ! GDDSK = 0.41*GDDM +145.4+150 !from hybrid-maize 
   ! GDDS1 = ((GDDSK-96)/38.9-4)*21
   ! GDDS1 = 0.77*GDDSK
   ! GDDS3 = GDDSK+170
   ! GDDS3 = 170

   IF(GDDDAY > 0.0) PGS = 2

   IF(GDDDAY >= parameters%GDDS1)  PGS = 3

   IF(GDDDAY >= parameters%GDDS2)  PGS = 4 

   IF(GDDDAY >= parameters%GDDS3)  PGS = 5

   IF(GDDDAY >= parameters%GDDS4)  PGS = 6

   IF(GDDDAY >= parameters%GDDS5)  PGS = 7

   IF(JULIAN >= parameters%HSDAY)  PGS = 8
 
   IF(JULIAN <  parameters%PLTDAY) PGS = 1   

END SUBROUTINE GROWING_GDD

!== begin psn_crop =================================================================================

SUBROUTINE PSN_CROP ( parameters,       & !in
                      SOLDN, XLAI,T2M,  & !in
                      PSNCROP        )    !out
!===================================================================================================

! input

  type (noahmp_parameters), intent(in) :: parameters
  REAL(r8)     , INTENT(IN)    :: SOLDN    ! downward solar radiation
  REAL(r8)     , INTENT(IN)    :: XLAI     ! LAI
  REAL(r8)     , INTENT(IN)    :: T2M      ! air temp
  REAL(r8)     , INTENT(OUT)   :: PSNCROP  !

!local

  REAL(r8)                     :: PAR      ! photosynthetically active radiation (w/m2) 1 W m-2 = 0.0864 MJ m-2 day-1
  REAL(r8)                     :: Amax     ! Maximum CO2 assimulation rate g/co2/s  
  REAL(r8)                     :: L1       ! Three Gaussian method
  REAL(r8)                     :: L2       ! Three Gaussian method
  REAL(r8)                     :: L3       ! Three Gaussian method
  REAL(r8)                     :: I1       ! Three Gaussian method
  REAL(r8)                     :: I2       ! Three Gaussian method
  REAL(r8)                     :: I3       ! Three Gaussian method
  REAL(r8)                     :: A1       ! Three Gaussian method
  REAL(r8)                     :: A2       ! Three Gaussian method
  REAL(r8)                     :: A3       ! Three Gaussian method
  REAL(r8)                     :: A        ! CO2 Assimulation 
  REAL(r8)                     :: TC

  TC = T2M - 273.15

  PAR = parameters%I2PAR * SOLDN * 0.0036  !w to MJ m-2

  IF(TC < parameters%TASSIM0) THEN
    Amax = 1E-10
  ELSEIF(TC >= parameters%TASSIM0 .and. TC < parameters%TASSIM1) THEN
    Amax = (TC - parameters%TASSIM0) * parameters%Aref / (parameters%TASSIM1 - parameters%TASSIM0)
  ELSEIF(TC >= parameters%TASSIM1 .and. TC < parameters%TASSIM2) THEN
    Amax = parameters%Aref
  ELSE
    Amax= parameters%Aref - 0.2 * (T2M - parameters%TASSIM2)
  ENDIF 
  
  Amax = max(amax,0.01)

  IF(XLAI <= 0.05) THEN
    L1 = 0.1127 * 0.05   !use initial LAI(0.05), avoid error
    L2 = 0.5    * 0.05
    L3 = 0.8873 * 0.05
  ELSE
    L1 = 0.1127 * XLAI
    L2 = 0.5    * XLAI
    L3 = 0.8873 * XLAI
  END IF

  I1 = parameters%k * PAR * exp(-parameters%k * L1)
  I2 = parameters%k * PAR * exp(-parameters%k * L2)
  I3 = parameters%k * PAR * exp(-parameters%k * L3)

  I1 = max(I1,1E-10)
  I2 = max(I2,1E-10)
  I3 = max(I3,1E-10)

  A1 = Amax * (1 - exp(-parameters%epsi * I1 / Amax))
  A2 = Amax * (1 - exp(-parameters%epsi * I2 / Amax)) * 1.6
  A3 = Amax * (1 - exp(-parameters%epsi * I3 / Amax))

  IF (XLAI <= 0.05) THEN
    A  = (A1+A2+A3) / 3.6 * 0.05
  ELSEIF (XLAI > 0.05 .and. XLAI <= 4.0) THEN
    A  = (A1+A2+A3) / 3.6 * XLAI
  ELSE
    A = (A1+A2+A3) / 3.6 * 4
  END IF

  A = A * parameters%PSNRF ! Attainable 

  PSNCROP = 6.313 * A   ! (1/44) * 1000000)/3600 = 6.313

END SUBROUTINE PSN_CROP

!== begin bvocflux =================================================================================

!  SUBROUTINE BVOCFLUX(parameters,VOCFLX,  VEGTYP,  VEGFRAC,  APAR,   TV )
!
! ------------------------------------------------------------------------------------------
!      implicit none
! ------------------------------------------------------------------------------------------
!
! ------------------------ code history ---------------------------
! source file:       BVOC
! purpose:           BVOC emissions
! DESCRIPTION:
! Volatile organic compound emission 
! This code simulates volatile organic compound emissions
! following the algorithm presented in Guenther, A., 1999: Modeling
! Biogenic Volatile Organic Compound Emissions to the Atmosphere. In
! Reactive Hydrocarbons in the Atmosphere, Ch. 3
! This model relies on the assumption that 90% of isoprene and monoterpene
! emissions originate from canopy foliage:
!    E = epsilon * gamma * density * delta
! The factor delta (longterm activity factor) applies to isoprene emission
! from deciduous plants only. We neglect this factor at the present time.
! This factor is discussed in Guenther (1997).
! Subroutine written to operate at the patch level.
! IN FINAL IMPLEMENTATION, REMEMBER:
! 1. may wish to call this routine only as freq. as rad. calculations
! 2. may wish to place epsilon values directly in pft-physiology file
! ------------------------ input/output variables -----------------
! input
!  integer                     ,INTENT(IN) :: vegtyp  !vegetation type 
!  real                        ,INTENT(IN) :: vegfrac !green vegetation fraction [0.0-1.0]
!  real                        ,INTENT(IN) :: apar    !photosynthesis active energy by canopy (w/m2)
!  real                        ,INTENT(IN) :: tv      !vegetation canopy temperature (k)
!
! output
!  real                        ,INTENT(OUT) :: vocflx(5) ! voc fluxes [ug C m-2 h-1]
!
! Local Variables
!
!  real, parameter :: R      = 8.314    ! univ. gas constant [J K-1 mol-1]
!  real, parameter :: alpha  = 0.0027   ! empirical coefficient
!  real, parameter :: cl1    = 1.066    ! empirical coefficient
!  real, parameter :: ct1    = 95000.0  ! empirical coefficient [J mol-1]
!  real, parameter :: ct2    = 230000.0 ! empirical coefficient [J mol-1]
!  real, parameter :: ct3    = 0.961    ! empirical coefficient
!  real, parameter :: tm     = 314.0    ! empirical coefficient [K]
!  real, parameter :: tstd   = 303.0    ! std temperature [K]
!  real, parameter :: bet    = 0.09     ! beta empirical coefficient [K-1]
!
!  integer ivoc        ! do-loop index
!  integer ityp        ! do-loop index
!  real epsilon(5)
!  real gamma(5)
!  real density
!  real elai
!  real par,cl,reciprod,ct
!
! epsilon :
!
!    do ivoc = 1, 5
!    epsilon(ivoc) = parameters%eps(VEGTYP,ivoc)
!    end do
!
! gamma : Activity factor. Units [dimensionless]
!
!      reciprod = 1. / (R * tv * tstd)
!      ct = exp(ct1 * (tv - tstd) * reciprod) / &
!           (ct3 + exp(ct2 * (tv - tm) * reciprod))
!
!      par = apar * 4.6 ! (multiply w/m2 by 4.6 to get umol/m2/s)
!      cl  = alpha * cl1 * par * (1. + alpha * alpha * par * par)**(-0.5)
!
!   gamma(1) = cl * ct ! for isoprenes
!
!   do ivoc = 2, 5
!   gamma(ivoc) = exp(bet * (tv - tstd))
!   end do
!
! Foliage density
!
! transform vegfrac to lai      
!
!   elai    = max(0.0,-6.5/2.5*alog((1.-vegfrac)))
!   density = elai / (parameters%slarea(VEGTYP) * 0.5)
!
! calculate the voc flux
!
!   do ivoc = 1, 5
!   vocflx(ivoc) = epsilon(ivoc) * gamma(ivoc) * density
!   end do
!
!   end subroutine bvocflux
! ==================================================================================================
! ********************************* end of carbon subroutines *****************************
! ==================================================================================================

!== begin noahmp_options ===========================================================================

  subroutine noahmp_options(idveg     ,iopt_crs  ,iopt_btr  ,iopt_run  ,iopt_sfc  ,iopt_frz , & 
                             iopt_inf  ,iopt_rad  ,iopt_alb  ,iopt_snf  ,iopt_tbot, iopt_stc, &
			     iopt_rsf )

  implicit none

  INTEGER,  INTENT(IN) :: idveg     !dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
  INTEGER,  INTENT(IN) :: iopt_crs  !canopy stomatal resistance (1-> Ball-Berry; 2->Jarvis)
  INTEGER,  INTENT(IN) :: iopt_btr  !soil moisture factor for stomatal resistance (1-> Noah; 2-> CLM; 3-> SSiB)
  INTEGER,  INTENT(IN) :: iopt_run  !runoff and groundwater (1->SIMGM; 2->SIMTOP; 3->Schaake96; 4->BATS)
  INTEGER,  INTENT(IN) :: iopt_sfc  !surface layer drag coeff (CH & CM) (1->M-O; 2->Chen97)
  INTEGER,  INTENT(IN) :: iopt_frz  !supercooled liquid water (1-> NY06; 2->Koren99)
  INTEGER,  INTENT(IN) :: iopt_inf  !frozen soil permeability (1-> NY06; 2->Koren99)
  INTEGER,  INTENT(IN) :: iopt_rad  !radiation transfer (1->gap=F(3D,cosz); 2->gap=0; 3->gap=1-Fveg)
  INTEGER,  INTENT(IN) :: iopt_alb  !snow surface albedo (1->BATS; 2->CLASS)
  INTEGER,  INTENT(IN) :: iopt_snf  !rainfall & snowfall (1-Jordan91; 2->BATS; 3->Noah)
  INTEGER,  INTENT(IN) :: iopt_tbot !lower boundary of soil temperature (1->zero-flux; 2->Noah)

  INTEGER,  INTENT(IN) :: iopt_stc  !snow/soil temperature time scheme (only layer 1)
                                    ! 1 -> semi-implicit; 2 -> full implicit (original Noah)
  INTEGER,  INTENT(IN) :: iopt_rsf  !surface resistance (1->Sakaguchi/Zeng; 2->Seller; 3->mod Sellers; 4->1+snow)

! -------------------------------------------------------------------------------------------------

  dveg = idveg
  
  opt_crs  = iopt_crs  
  opt_btr  = iopt_btr  
  opt_run  = iopt_run  
  opt_sfc  = iopt_sfc  
  opt_frz  = iopt_frz  
  opt_inf  = iopt_inf  
  opt_rad  = iopt_rad  
  opt_alb  = iopt_alb  
  opt_snf  = iopt_snf  
  opt_tbot = iopt_tbot 
  opt_stc  = iopt_stc
  opt_rsf  = iopt_rsf
  
  end subroutine noahmp_options
 
END MODULE MODULE_SF_NOAHMPLSM
