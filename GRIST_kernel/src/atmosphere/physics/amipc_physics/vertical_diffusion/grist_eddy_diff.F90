!===================================================================================================
!
!  Created by LiXiaohan on 19/06/04, adopted from CAM5
!
!  The University of Washington Moist Turbulence Scheme to compute eddy diffusion
!  coefficients associated with dry and moist turbulences in the whole atmospheric layers.
!  
!  Reference:
!  Bretherton, C., and S. Park, 2009: A New Moist Turbulence Parameterization in the Community
!  Atmosphere Model. J. Clim., 22, 3422–3448.
!
!  Park, S., and Bretherton, C., 2009: The University of Washington shallow convection and moist
!  turbulence schemes and their impact on climate simulations with the Community Atmosphere Model.
!  J. Clim., 22, 3449-3469.
!
!===================================================================================================

module grist_eddy_diff

    use grist_constants,                    only: r8, i4
    use grist_wv_saturation,                only: qsat
    use grist_handle_error,                 only: endrun
    use grist_diffusion_solver,             only: vdiff_selector
    use grist_mpi


    implicit none
    private
    public :: init_eddy_diff,       &
              compute_eddy_diff

    type(vdiff_selector)        :: fieldlist_wet                  ! Logical switches for moist mixing ratio diffusion
    type(vdiff_selector)        :: fieldlist_molec                ! Logical switches for molecular diffusion
! PBL Parameters used in the UW PBL
    character,        parameter :: sftype         = 'l'           ! Method for calculating saturation fraction

    character(len=4), parameter :: choice_evhc    = 'maxi'        ! 'orig', 'ramp', 'maxi' : recommended to be used with choice_radf 
    character(len=6), parameter :: choice_radf    = 'maxi'        ! 'orig', 'ramp', 'maxi' : recommended to be used with choice_evhc 
    character(len=6), parameter :: choice_SRCL    = 'nonamb'      ! 'origin', 'remove', 'nonamb'
 
    character(len=6), parameter :: choice_tunl    = 'rampcl'      ! 'origin', 'rampsl'(Sungsu), 'rampcl'(Chris)
    real(r8),         parameter :: ctunl          =  2._r8        !  Maximum asympt leng = ctunl*tunl when choice_tunl = 'rampsl(cl)'
    character(len=6), parameter :: choice_leng    = 'origin'      ! 'origin', 'takemn'
    real(r8),         parameter :: cleng          =  3._r8        !  Order of 'leng' when choice_leng = 'origin' [ no unit ]
    character(len=6), parameter :: choice_tkes    = 'ibprod'      ! 'ibprod' (include tkes in computing bprod), 'ebprod'(exclude)

    real(r8)                    :: lbulk_max      =  40.e3_r8     ! Maximum master length scale designed to address issues in the
                                                                  ! upper atmosphere where vertical model resolution is coarse [ m ].
                                                                  ! In order not to disturb turbulence characteristics in the lower
                                                                  ! troposphere, this should be set at least larger than ~ a few km.  
                                                                  ! atmosphere.

! Parameters for 'sedimentation-entrainment feedback' for liquid stratus
! If .false.,  no sedimentation entrainment feedback ( i.e., use default evhc )

    logical,          parameter :: id_sedfact     = .false.
    real(r8),         parameter :: ased           =  9._r8        !  Valid only when id_sedfact = .true.

! Parameters governing entrainment efficiency A = a1l(i)*evhc, evhc = 1 + a2l * a3l * L * ql / jt2slv 
! Here, 'ql' is cloud-top LWC and 'jt2slv' is the jump in 'slv' across                                
! the cloud-top entrainment zone ( across two grid layers to consider full mixture )
    real(r8),         parameter :: a1l            =   0.10_r8     ! Dry entrainment efficiency for TKE closure
                                                                  ! a1l = 0.2*tunl*erat^-1.5,
                                                                  ! where erat = <e>/wstar^2 for dry CBL =  0.3.

    real(r8),         parameter :: a1i            =   0.2_r8      ! Dry entrainment efficiency for wstar closure
    real(r8),         parameter :: ccrit          =   0.5_r8      ! Minimum allowable sqrt(tke)/wstar.
                                                                  ! Used in solving cubic equation for 'ebrk'
    real(r8),         parameter :: wstar3factcrit =   0.5_r8      ! 1/wstar3factcrit is the maximally allowed enhancement of
                                                                  ! 'wstar3' due to entrainment.

    real(r8),         parameter :: a2l            =   10._r8 !LiXH tuning, default: 30._r8  DP:10    ! Moist entrainment enhancement param (recommended range : 10~50 )
    real(r8),         parameter :: a3l            =   0.8_r8      ! Approximation to a complicated thermodynamic parameters

    real(r8),         parameter :: jbumin         =   .001_r8     ! Minimum buoyancy jump at an entrainment jump, [m/s2]
    real(r8),         parameter :: evhcmax        =   10._r8      ! Upper limit of evaporative enhancement factor

    real(r8),         parameter :: onet           =   1._r8/3._r8 ! 1/3 power in wind gradient expression [ no unit ]
    real(r8),         parameter :: qmin           =   1.e-5_r8    ! Minimum grid-mean LWC counted as clouds [kg/kg]
    real(r8),         parameter :: ntzero         =   1.e-12_r8   ! Not zero (small positive number used in 's2')
    real(r8),         parameter :: b1             =   5.8_r8      ! TKE dissipation D = e^3/(b1*leng), e = b1*W.
    real(r8)                    :: b123                           ! b1**(2/3)
    real(r8),         parameter :: tunl           =   0.085_r8    ! Asympt leng = tunl*(turb lay depth)
    real(r8),         parameter :: alph1          =   0.5562_r8   ! alph1~alph5 : Galperin instability function parameters
    real(r8),         parameter :: alph2          =  -4.3640_r8   !               These coefficients are used to calculate 
    real(r8),         parameter :: alph3          = -34.6764_r8   !               'sh' and 'sm' from 'gh'.
    real(r8),         parameter :: alph4          =  -6.1272_r8   !
    real(r8),         parameter :: alph5          =   0.6986_r8   !
    real(r8),         parameter :: ricrit         =   0.19_r8     ! Critical Richardson number for turbulence.
                                                                  ! Can be any value >= 0.19.
    real(r8),         parameter :: ae             =   1._r8       ! TKE transport efficiency [no unit]
    real(r8),         parameter :: rinc           =  -0.04_r8     ! Minimum W/<W> used for CL merging test 
    real(r8),         parameter :: wpertmin       =   1.e-6_r8    ! Minimum PBL eddy vertical velocity perturbation
    real(r8),         parameter :: wfac           =   1._r8       ! Ratio of 'wpert' to sqrt(tke) for CL.
    real(r8),         parameter :: tfac           =   1._r8       ! Ratio of 'tpert' to (w't')/wpert for CL.
                                                                  ! Same ratio also used for q
    real(r8),         parameter :: fak            =   8.5_r8      ! Constant in surface temperature excess for stable STL.
    real(r8),         parameter :: rcapmin        =   0.1_r8      ! Minimum allowable e/<e> in a CL
    real(r8),         parameter :: rcapmax        =   2.0_r8      ! Maximum allowable e/<e> in a CL
    real(r8),         parameter :: tkemax         =  20._r8       ! TKE is capped at tkemax [m2/s2]
    real(r8),         parameter :: lambda         =   0.5_r8      ! Under-relaxation factor ( 0 < lambda =< 1 )

    logical,          parameter :: use_dw_surf    =  .true.       ! Used in 'zisocl'. Default is 'true'
                                                                  ! If 'true', surface interfacial energy does not contribute
                                                                  ! to the CL mean stability functions after finishing merging.
                                                                  ! For this case, 'dl2n2_surf' is only used for a merging test
                                                                  ! based on 'l2n2'
                                                                  ! If 'false',surface interfacial enery explicitly contribute to
                                                                  ! CL mean stability functions after finishing merging.
                                                                  ! For this case, 'dl2n2_surf' and 'dl2s2_surf' are directly used
                                                                  ! for calculating surface interfacial layer energetics

    logical,          parameter :: set_qrlzero    =  .false.      ! .true. ( .false.) : turning-off ( on) radiative-turbulence
                                                                  ! interaction by setting qrl = 0.
! PBL Parameters not used in the UW PBL 
    real(r8),         parameter :: pblmaxp        =  4.e4_r8      ! PBL max depth in pressure units. 
    real(r8),         parameter :: zkmin          =  0.01_r8      ! Minimum kneutral*f(ri). 
    real(r8),         parameter :: betam          = 15.0_r8       ! Constant in wind gradient expression.
    real(r8),         parameter :: betas          =  5.0_r8       ! Constant in surface layer gradient expression.
    real(r8),         parameter :: betah          = 15.0_r8       ! Constant in temperature gradient expression.
    real(r8),         parameter :: fakn           =  7.2_r8       ! Constant in turbulent prandtl number.
    real(r8),         parameter :: ricr           =  0.3_r8       ! Critical richardson number.
    real(r8),         parameter :: sffrac         =  0.1_r8       ! Surface layer fraction of boundary layer
    real(r8),         parameter :: binm           =  betam*sffrac ! betam * sffrac
    real(r8),         parameter :: binh           =  betah*sffrac ! betah * sffrac

! PBL constants set using values from other parts of code 
    integer         :: ncvmax       ! Max numbers of CLs (good to set to 'pver')
    real(r8)        :: cpair        ! Specific heat of dry air
    real(r8)        :: rair         ! Gas const for dry air
    real(r8)        :: zvir         ! rh2o/rair - 1
    real(r8)        :: latvap       ! Latent heat of vaporization
    real(r8)        :: latice       ! Latent heat of fusion
    real(r8)        :: latsub       ! Latent heat of sublimation
    real(r8)        :: g            ! Gravitational acceleration
    real(r8)        :: vk           ! Von Karman's constant
    real(r8)        :: ccon         ! fak * sffrac * vk
    integer         :: ntop_turb    ! Top interface level to which turbulent vertical diffusion is applied ( = 1 )
    integer         :: nbot_turb    ! Bottom interface level to which turbulent vertical diffusion is applied ( = pver )

    real(r8), allocatable :: ml2(:) ! Mixing lengths squared. Not used in the UW PBL.
                                    ! Used for computing free air diffusivity.
    real(r8), allocatable :: leng_max(:)     ! Maximum length scale designed to address issues in the upper


contains
    subroutine init_eddy_diff( pver, gravx, cpairx, rairx, zvirx, latvapx, laticex,         &
                               ntop_eddy, nbot_eddy, vkx, eddy_lbulk_max, eddy_leng_max,    &
                               eddy_max_bot_pressure, pref_mid)
     use grist_diffusion_solver, only: new_fieldlist_vdiff, vdiff_select
! io
    integer(i4),  intent(in) :: pver                    ! Number of vertical layers
    integer(i4),  intent(in) :: ntop_eddy               ! Top interface level to which eddy vertical diffusivity is applied ( = 1 )
    integer(i4),  intent(in) :: nbot_eddy               ! Bottom interface level to which eddy vertical diffusivity is applied ( = pver)
    real(r8),     intent(in) :: gravx                   ! Acceleration of gravity
    real(r8),     intent(in) :: cpairx                  ! Specific heat of dry air
    real(r8),     intent(in) :: rairx                   ! Gas constant for dry air
    real(r8),     intent(in) :: zvirx                   ! rh2o/rair - 1
    real(r8),     intent(in) :: latvapx                 ! Latent heat of vaporization
    real(r8),     intent(in) :: laticex                 ! Latent heat of fusion
    real(r8),     intent(in) :: vkx                     ! Von Karman's constant
    real(r8),     intent(in) :: eddy_lbulk_max          ! Maximum master length scale
    real(r8),     intent(in) :: eddy_leng_max           ! Maximum dissipation length scale
    real(r8),     intent(in) :: eddy_max_bot_pressure   ! Bottom pressure level (hPa) at which namelist leng_max and lbulk_max
    real(r8),     intent(in) :: pref_mid(pver)

! local
    integer                  :: k

    ncvmax    = pver
    cpair     = cpairx
    rair      = rairx
    g         = gravx
    zvir      = zvirx
    latvap    = latvapx
    latice    = laticex
    latsub    = latvap + latice
    vk        = vkx
    ccon      = fak*sffrac*vk
    ntop_turb = ntop_eddy
    nbot_turb = nbot_eddy
    b123      = b1**(2._r8/3._r8)
    
    allocate(leng_max(pver));leng_max(:) = 40.e3_r8

    lbulk_max = eddy_lbulk_max
    do k = 1,pver
      if ( pref_mid(k) .le. eddy_max_bot_pressure*1.D2 ) leng_max(k)  = eddy_leng_max
    end do

    ! Set the square of the mixing lengths. Only for CAM3 HB PBL scheme.
    ! Not used for UW moist PBL. Used for free air eddy diffusivity.

    allocate(ml2(pver+1))
    ml2(1:ntop_turb) = 0._r8
    do k = ntop_turb + 1, nbot_turb
       ml2(k) = 30.0_r8**2
    end do
    ml2(nbot_turb+1:pver+1) = 0._r8
 
    ! Get fieldlists to pass to diffusion solver.
    fieldlist_wet   = new_fieldlist_vdiff(1)
    fieldlist_molec = new_fieldlist_vdiff(1)

    ! Select the fields which will be diffused 
    if(vdiff_select(fieldlist_wet,'s').ne.'')   call endrun( vdiff_select(fieldlist_wet,'s') )
    if(vdiff_select(fieldlist_wet,'q',1).ne.'') call endrun( vdiff_select(fieldlist_wet,'q',1) )
    if(vdiff_select(fieldlist_wet,'u').ne.'')   call endrun( vdiff_select(fieldlist_wet,'u') )
    if(vdiff_select(fieldlist_wet,'v').ne.'')   call endrun( vdiff_select(fieldlist_wet,'v') )

    deallocate(ml2)

!    The part below is not used in GRIST, LiXH
!    call addfld(UW_vars)  

    end subroutine init_eddy_diff

 
! Purpose: Interface to compute eddy diffusivities.
!          Eddy diffusivities are calculated in a fully implicit way through iteration process. 
    subroutine compute_eddy_diff(ncol,     pver,    t,        qv,       ztodt,  rpdel,   cldn,    qrl,     &
                                 wsedl,    z,       zi,       pmid,     pi,     u,       v,       taux,    &
                                 tauy,     shflx,   qflx,     wstarent, nturb,  rrho,    ustar,   pblh,    &
                                 kvm_in,   kvh_in,  kvm_out,  kvh_out,  kvq,    cgh,     cgs,     bprod,   &
                                 sprod,    tke,     wpert,    tpert,    qpert,  sfi,     kvinit,  tauresx, &
                                 ! for Lin Macro ------------------------>
                                 !tauresy,  ksrftms, turbtype, sm_aw,    ipbl,   kpblh,   wstarPBL )
                                 tauresy,  ksrftms, turbtype, sm_aw,    ipbl,   kpblh,  &
                                 wstarPBL, lengi,   shi,      smi )
                                 ! for Lin Macro <------------------------

    use grist_pbl_utils,                    only: calc_ustar
    use grist_diffusion_solver,             only: compute_vdiff
! io
    integer(i4),  intent(in)    :: ncol                  ! Number of columns
    integer(i4),  intent(in)    :: pver                  ! Number of vertical layers
    integer(i4),  intent(in)    :: nturb                 ! Number of iteration steps
    logical,      intent(in)    :: wstarent              ! use the 'wstar' entrainment closure
    logical,      intent(in)    :: kvinit                ! time step = 1 : initializing kvh, kvm
                                                         ! kvh, kvm = 0 in UW moist PBL scheme
    real(r8),     intent(in)    :: ztodt                 ! time step [s]
    real(r8),     intent(in)    :: t(pver, ncol)         ! Temperature [K]
    real(r8),     intent(in)    :: qv(1:3, pver, ncol)   ! Water vapor(1), liquid water(2), and ice(3) specific humidity [ kg/kg ]
    real(r8),     intent(in)    :: rpdel(pver, ncol)     ! 1./pdel [Pa]
    real(r8),     intent(in)    :: cldn(pver, ncol)      ! Stratiform cloud fraction
    real(r8),     intent(in)    :: qrl(pver, ncol)       ! Long wave radiation cooling rate
    real(r8),     intent(in)    :: wsedl(pver, ncol)     ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
    real(r8),     intent(in)    :: z(pver, ncol)         ! Layer mid-point height above surface [ m ]
    real(r8),     intent(in)    :: zi(pver+1, ncol)      ! Interface height above surface [ m ]
    real(r8),     intent(in)    :: pmid(pver, ncol)      ! Layer mid-point pressure [ Pa ]
    real(r8),     intent(in)    :: pi(pver+1, ncol)      ! Interface pressure [ Pa ]
    real(r8),     intent(in)    :: u(pver, ncol)         ! Zonal velocity [ m/s ]
    real(r8),     intent(in)    :: v(pver, ncol)         ! Meridional velocity [ m/s ]
    real(r8),     intent(in)    :: taux(ncol)            ! Zonal wind stress at surface [ N/m2 ]
    real(r8),     intent(in)    :: tauy(ncol)            ! Meridional wind stress at surface [ N/m2 ]
    real(r8),     intent(in)    :: shflx(ncol)           ! Sensible heat flux at surface [ W/m2 ]
    real(r8),     intent(in)    :: qflx(ncol)            ! Water vapor flux at surface [ kg/m2/s ]
    real(r8),     intent(in)    :: kvm_in(pver+1, ncol)  ! kvm saved from last timestep [ m2/s ]
    real(r8),     intent(in)    :: kvh_in(pver+1, ncol)  ! kvh saved from last timestep [ m2/s ]
    real(r8),     intent(in)    :: ksrftms(ncol)         ! Surface drag coefficient of turbulent mountain stress

    real(r8),     intent(inout) :: tauresx(ncol)         ! Residual stress to be added in vdiff to correct for turb
    real(r8),     intent(inout) :: tauresy(ncol)         ! Stress mismatch between sfc and atm accumulated in prior timesteps
    real(r8),     intent(out)   :: kvm_out(pver+1, ncol) ! Eddy diffusivity for momentum [ m2/s ]
    real(r8),     intent(out)   :: kvh_out(pver+1, ncol) ! Eddy diffusivity for heat [ m2/s ]
    real(r8),     intent(out)   :: kvq(pver+1, ncol)     ! Eddy diffusivity for tracers [ m2/s ]
    real(r8),     intent(out)   :: rrho(ncol)            ! Reciprocal of density at the lowest layer
    real(r8),     intent(out)   :: ustar(ncol)           ! Surface friction velocity [ m/s ]
    real(r8),     intent(out)   :: pblh(ncol)            ! PBL top height [ m ]
    real(r8),     intent(out)   :: cgh(pver+1, ncol)     ! Counter-gradient term for heat [ J/kg/m ]
    real(r8),     intent(out)   :: cgs(pver+1, ncol)     ! Counter-gradient star [ cg/flux ]
    real(r8),     intent(out)   :: tpert(ncol)           ! Convective temperature excess [ K ] 
    real(r8),     intent(out)   :: qpert(ncol)           ! Convective humidity excess [ kg/kg ]
    real(r8),     intent(out)   :: wpert(ncol)           ! Turbulent velocity excess [ m/s ]
    real(r8),     intent(out)   :: tke(pver+1, ncol)     ! Turbulent kinetic energy [ m2/s2 ]
    real(r8),     intent(out)   :: bprod(pver+1, ncol)   ! Buoyancy production [ m2/s3 ] 
    real(r8),     intent(out)   :: sprod(pver+1, ncol)   ! Shear production [ m2/s3 ]
    real(r8),     intent(out)   :: sfi(pver+1, ncol)     ! Interfacial layer saturation fraction
    integer(i4),  intent(out)   :: turbtype(pver+1, ncol)! Turbulence type identifier at all interfaces
    real(r8),     intent(out)   :: sm_aw(pver+1, ncol)   ! Normalized Galperin instability function for momentum
                                                         ! = 1 for neutral condition (Ri=0),
                                                         ! = 4.964 for maximum unstable case, 
                                                         ! = 0 for Ri > Ricrit=0.19.
    real(r8), intent(out)       :: ipbl(ncol)            ! If 1, PBL is CL,while if 0, PBL is STL.
    real(r8), intent(out)       :: kpblh(ncol)           ! Layer index containing PBL top within or at the base interface
    real(r8), intent(out)       :: wstarPBL(ncol)        ! Convective velocity within PBL [ m/s ]
    ! for Lin Macro ------------------------>
    real(r8), intent(out)   :: shi(pver+1, ncol)
    real(r8), intent(out)   :: smi(pver+1, ncol)
    real(r8), intent(out)   :: lengi(pver+1, ncol)
    ! for Lin Macro <------------------------

! local
    integer(i4)   :: icol, i, k, iturb
    real(r8)      :: kvh(pver+1, ncol)         ! Eddy diffusivity for heat [ m2/s ]
    real(r8)      :: kvm(pver+1, ncol)         ! Eddy diffusivity for momentum [ m2/s ]
    real(r8)      :: errorPBL(ncol)            ! Error function showing whether PBL produced convergent solution or not
    real(r8)      :: s2(pver, ncol)            ! Shear squared, defined at interfaces except surface [ s-2 ]
    real(r8)      :: n2(pver, ncol)            ! Buoyancy frequency, defined at interfaces except surface [ s-2 ]
    real(r8)      :: ri(pver, ncol)            ! Richardson number, 'n2/s2', defined at interfaces except surface [ s-2 ]
    real(r8)      :: pblhp(ncol)               ! PBL top pressure [ Pa ]
    real(r8)      :: minpblh(ncol)             ! Minimum PBL height
    real(r8)      :: qi(pver, ncol)            ! Ice specific humidity [ kg/kg ]
    real(r8)      :: qt(pver, ncol)            ! Total specific humidity [ kg/kg ]
    real(r8)      :: sfuh(pver, ncol)          ! Saturation fraction in upper half-layer
    real(r8)      :: sflh(pver, ncol)          ! Saturation fraction in lower half-layer
    real(r8)      :: sl(pver, ncol)            ! Liquid water static energy [ J/kg ]
    real(r8)      :: slv(pver, ncol)           ! Liquid water virtual static energy [ J/kg ]
    real(r8)      :: slslope(pver, ncol)       ! Slope of 'sl' in each layer
    real(r8)      :: qtslope(pver, ncol)       ! Slope of 'qt' in each layer
    real(r8)      :: qvfd(pver, ncol)          ! Specific humidity for diffusion [ kg/kg ]
    real(r8)      :: tfd(pver, ncol)           ! Temperature for diffusion [ K ]
    real(r8)      :: slfd(pver, ncol)          ! Liquid static energy [ J/kg ]
    real(r8)      :: qtfd(pver, ncol)          ! Total specific humidity [ kg/kg ]
    real(r8)      :: qlfd(pver, ncol)          ! Liquid water specific humidity for diffusion [ kg/kg ]
    real(r8)      :: ufd(pver, ncol)           ! U-wind for diffusion [ m/s ]
    real(r8)      :: vfd(pver, ncol)           ! V-wind for diffusion [ m/s ]
    real(r8)      :: chs(pver+1, ncol)         ! Heat buoyancy coef for sat states, defined at each interface
    real(r8)      :: chu(pver+1, ncol)         ! Heat buoyancy coef for dry states, defined at each interface
    real(r8)      :: cms(pver+1, ncol)         ! Moisture buoyancy coef for sat states, defined at each interface
    real(r8)      :: cmu(pver+1, ncol)         ! Moisture buoyancy coef for dry states, defined at each interface
    real(r8)      :: jnk1d(ncol)
    real(r8)      :: jnk2d(pver+1,ncol)
    real(r8)      :: zero(ncol)
    real(r8)      :: zero2d(pver+1, ncol)
    real(r8)      :: es                        ! Saturation vapor pressure
    real(r8)      :: qs                        ! Saturation specific humidity
    real(r8)      :: ep2, templ, temps
!-----------------------------------Lixh Test-----------------------------------
! if waccmx is avaiable, cpairv shoule be updated in subroutine physconst_update
    real(r8)      :: cpairv(pver, ncol)
!-----------------------------------Lixh Test-----------------------------------

! Variables for diagnostic output
    real(r8)      :: tkes(ncol)                ! TKE at surface interface [ m2/s2 ]
    real(r8)      :: kbase_o(ncvmax, ncol)     ! Original external base interface index of CL from 'exacol'
    real(r8)      :: ktop_o(ncvmax, ncol)      ! Original external top  interface index of CL from 'exacol'
    real(r8)      :: ncvfin_o(ncol)            ! Original number of CLs from 'exacol'
    real(r8)      :: kbase_mg(ncvmax, ncol)    ! 'kbase' after extending-merging from 'zisocl'
    real(r8)      :: ktop_mg(ncvmax, ncol)     ! 'ktop' after extending-merging from 'zisocl'
    real(r8)      :: ncvfin_mg(ncol)           ! 'ncvfin' after extending-merging from 'zisocl'
    real(r8)      :: kbase_f(ncvmax, ncol)     ! Final 'kbase' after extending-merging & including SRCL
    real(r8)      :: ktop_f(ncvmax, ncol)      ! Final 'ktop' after extending-merging & including SRCL
    real(r8)      :: ncvfin_f(ncol)            ! Final 'ncvfin' after extending-merging & including SRCL
    real(r8)      :: wet(ncvmax, ncol)         ! Entrainment rate at the CL top  [ m/s ] 
    real(r8)      :: web(ncvmax, ncol)         ! Entrainment rate at the CL base [ m/s ].
                                               ! Set to zero if CL is based at surface.
    real(r8)      :: jtbu(ncvmax, ncol)        ! Buoyancy jump across the CL top  [ m/s2 ]  
    real(r8)      :: jbbu(ncvmax, ncol)        ! Buoyancy jump across the CL base [ m/s2 ]  
    real(r8)      :: evhc(ncvmax, ncol)        ! Evaporative enhancement factor at the CL top
    real(r8)      :: jt2slv(ncvmax, ncol)      ! Jump of slv ( across two layers ) at CL top used only for evhc [ J/kg ]
    real(r8)      :: n2ht(ncvmax, ncol)        ! n2 defined at the CL top  interface but using
                                               ! sfuh(kt)   instead of sfi(kt) [ s-2 ] 
    real(r8)      :: n2hb(ncvmax, ncol)        ! n2 defined at the CL base interface but using
                                               ! sflh(kb-1) instead of sfi(kb) [ s-2 ]
    real(r8)      :: lwp(ncvmax, ncol)         ! LWP in the CL top layer [ kg/m2 ]
    real(r8)      :: opt_depth(ncvmax, ncol)   ! Optical depth of the CL top layer
    real(r8)      :: radinvfrac(ncvmax, ncol)  ! Fraction of radiative cooling confined in the top portion of CL top layer
    real(r8)      :: radf(ncvmax, ncol)        ! Buoyancy production at the CL top due to LW radiative cooling [ m2/s3 ]
    real(r8)      :: wstar(ncvmax, ncol)       ! Convective velocity in each CL [ m/s ]
    real(r8)      :: wstar3fact(ncvmax, ncol)  ! Enhancement of 'wstar3' due to entrainment (inverse) [ no unit ]
    real(r8)      :: ebrk(ncvmax, ncol)        ! Net mean TKE of CL including entrainment effect [ m2/s2 ]
    real(r8)      :: wbrk(ncvmax, ncol)        ! Net mean normalized TKE (W) of CL,
                                               ! 'ebrk/b1' including entrainment effect [ m2/s2 ]
    real(r8)      :: lbrk(ncvmax, ncol)        ! Energetic internal thickness of CL [m]
    real(r8)      :: ricl(ncvmax, ncol)        ! CL internal mean Richardson number
    real(r8)      :: ghcl(ncvmax, ncol)        ! Half of normalized buoyancy production of CL
    real(r8)      :: shcl(ncvmax, ncol)        ! Galperin instability function of heat-moisture of CL
    real(r8)      :: smcl(ncvmax, ncol)        ! Galperin instability function of mementum of CL
    real(r8)      :: ghi(pver+1, ncol)         ! Half of normalized buoyancy production at all interfaces
    ! for Lin Macro ------------------------>
    !real(r8)      :: shi(pver+1, ncol)         ! Galperin instability function of heat-moisture at all interfaces
    !real(r8)      :: smi(pver+1, ncol)         ! Galperin instability function of heat-moisture at all interfaces
    !real(r8)      :: lengi(pver+1, ncol)       ! Turbulence length scale at all interfaces [ m ]
    ! for Lin Macro <------------------------
    real(r8)      :: rii(pver+1, ncol)         ! Interfacial Richardson number defined at all interfaces
    real(r8)      :: wcap(pver+1, ncol)        ! Normalized TKE at all interfaces [ m2/s2 ]
    real(r8)      :: rairi(pver+1, ncol)       ! interface gas constant needed for compute_vdiff



    zero(:)     = 0._r8
    zero2d(:,:) = 0._r8

    ufd(:, :ncol)  = u(:,:ncol)
    vfd(:, :ncol)  = v(:,:ncol)
    tfd(:, :ncol)  = t(:,:ncol)
    qvfd(:,:ncol)  = qv(1,:,:ncol)
    qlfd(:,:ncol)  = qv(2,:,:ncol)
    qi(:,:ncol)    = qv(3,:,:ncol)

    do iturb = 1, nturb

        ! Total stress includes 'tms'.
        ! Here, in computing 'tms', we can use either iteratively changed 'ufd,vfd' or the
        ! initially given 'u,v' to the PBL scheme. Note that normal stress, 'taux, tauy'
        ! are not changed by iteration. In order to treat 'tms' in a fully implicit way,
        ! I am using updated wind, here.

        ! Compute ustar
        call calc_ustar( tfd(pver,:ncol), pmid(pver,:ncol),                  &
                         taux(:ncol) - ksrftms(:ncol) * ufd(pver,:ncol),     & ! Zonal wind stress
                         tauy(:ncol) - ksrftms(:ncol) * vfd(pver,:ncol),     & ! Meridional wind stress
                         rrho(:ncol), ustar(:ncol))
        minpblh(:ncol) = 100.0_r8 * ustar(:ncol)   ! By construction, 'minpblh' is larger than 1 [m] when 'ustar_min = 0.01'.

        ! Calculate (qt,sl,n2,s2,ri) from a given set of (t,qv,ql,qi,u,v)
        call trbintd( pver,   ncol,    z,       ufd,    vfd,    tfd,   pmid, &
                      s2,     n2,      ri,      zi,     pi,     cldn,  qtfd, &
                      qvfd,   qlfd,    qi,      sfi,    sfuh,   sflh,  slfd, &
                      slv,    slslope, qtslope, chs,    chu,    cms,   cmu   )

        ! Save initial (i.e., before iterative diffusion) profile of (qt,sl) at each iteration.         
        ! Only necessary for (qt,sl) not (u,v) because (qt,sl) are newly calculated variables. 
        if( iturb .eq. 1 ) then
            qt(:,:ncol) = qtfd(:,:ncol)
            sl(:,:ncol) = slfd(:,:ncol)
        endif

        ! Initialize kvh/kvm to send to caleddy, depending on model timestep and iteration number
        ! This is necessary for 'wstar-based' entrainment closure.
        if(iturb .eq. 1) then
        ! LiXH modified for restart
        !    if( kvinit ) then
        !        ! First iteration of first model timestep : Use free tropospheric value or zero.
        !        kvh(:,:ncol) = 0._r8
        !        kvm(:,:ncol) = 0._r8
        !    else
                ! First iteration on any model timestep except the first : Use value from previous timestep
                kvh(:,:ncol) = kvh_in(:,:ncol)
                kvm(:,:ncol) = kvm_in(:,:ncol)
        !    endif
        else
            ! Not the first iteration : Use from previous iteration
                kvh(:,:ncol) = kvh_out(:,:ncol)
                kvm(:,:ncol) = kvm_out(:,:ncol)
        end if



        ! Calculate eddy diffusivity (kvh_out,kvm_out) and (tke,bprod,sprod) using
        ! a given (kvh,kvm) which are used only for initializing (bprod,sprod)  at
        ! the first part of caleddy. (bprod,sprod) are fully updated at the end of
        ! caleddy after calculating (kvh_out,kvm_out) 
        call caleddy( pver      , ncol      ,                                   &
                      slfd      , qtfd      , qlfd      , slv      ,ufd       , &
                      vfd       , pi        , z         , zi       ,            &
                      qflx      , shflx     , slslope   , qtslope  ,            &
                      chu       , chs       , cmu       , cms      ,sfuh      , &
                      sflh      , n2        , s2        , ri       ,rrho      , &
                      pblh      , ustar     ,                                   &
                      kvh       , kvm       , kvh_out   , kvm_out  ,            &
                      tpert     , qpert     , qrl       , tke      ,            &
                      wstarent  , bprod     , sprod     , minpblh  , wpert    , &
                      tkes      , turbtype  , sm_aw     ,                       & 
                      kbase_o   , ktop_o    , ncvfin_o  ,                       &
                      kbase_mg  , ktop_mg   , ncvfin_mg ,                       &                  
                      kbase_f   , ktop_f    , ncvfin_f  ,                       &                  
                      wet       , web       , jtbu      , jbbu     ,            &
                      evhc      , jt2slv    , n2ht      , n2hb     , lwp      , &
                      opt_depth , radinvfrac, radf      , wstar    , wstar3fact,&
                      ebrk      , wbrk      , lbrk      , ricl     , ghcl     , & 
                      shcl      , smcl      ,                                   &
                      ghi       , shi       , smi       , rii      , lengi    , &
                      wcap      , pblhp     , cldn      ,                       &
                      ipbl      , kpblh     , wsedl )


        ! Calculate errorPBL to check whether PBL produced convergent solutions or not.

        if( iturb .eq. nturb ) then
           do i = 1, ncol
              errorPBL(i) = 0._r8 
              do k = 1, pver
                 errorPBL(i) = errorPBL(i) + ( kvh(k,i) - kvh_out(k,i) )**2 
              end do 
              errorPBL(i) = sqrt(errorPBL(i)/pver)
           end do
        end if

        ! Eddy diffusivities which will be used for the initialization of (bprod,
        ! sprod) in 'caleddy' at the next iteration step.

        if( iturb .gt. 1 .and. iturb .lt. nturb ) then
           kvm_out(:,:ncol) = lambda * kvm_out(:,:ncol) + ( 1._r8 - lambda ) * kvm(:,:ncol)
           kvh_out(:,:ncol) = lambda * kvh_out(:,:ncol) + ( 1._r8 - lambda ) * kvh(:,:ncol)
        endif

     ! Set nonlocal terms to zero for flux diagnostics, since not used by caleddy.

       cgh(:,:ncol) = 0._r8
       cgs(:,:ncol) = 0._r8      

       if( iturb .lt. nturb ) then

         ! Each time we diffuse the original state

           slfd(:,:ncol)  = sl(:,:ncol)
           qtfd(:,:ncol)  = qt(:,:ncol)
           ufd(:,:ncol)   = u(:,:ncol)
           vfd(:,:ncol)   = v(:,:ncol)

         ! Check to see if constituent dependent gas constant needed (WACCM-X)
!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!         if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
!           rairi(1,:ncol) = rairv(1,:ncol,lchnk)
!           do k = 2, pver
!             do i = 1, ncol
!               rairi(k,i) = 0.5_r8 * (rairv(k,i,lchnk)+rairv(k-1,i,lchnk))
!             end do
!           end do      
!         else
!           rairi(:pver+1,:ncol) = rair
!         endif
         
         rairi(:pver+1,:ncol) = rair
         cpairv(:pver,:ncol)  = cpair
!<----------------------------------Lixh close this part-----------------------------------

         ! Diffuse initial profile of each time step using a given (kvh_out,kvm_out)
         ! In the below 'compute_vdiff', (slfd,qtfd,ufd,vfd) are 'inout' variables.
         call compute_vdiff( pver    , 1        , ncol     , pmid         ,             &   
                             pi      , rpdel    , t        , ztodt        , taux      , &
                             tauy    , shflx    , qflx     , ntop_turb    , nbot_turb , &
                             kvh_out , kvm_out  , kvh_out  , cgs          , cgh       , &
                             zi      , ksrftms  , zero     , fieldlist_wet,             &
                             ufd     , vfd      , qtfd     , slfd         ,             &
                             jnk1d   , jnk1d    , jnk2d    , jnk1d        ,             &
                             tauresx , tauresy  , 0        , cpairv       , rairi , .false. )


         ! Retrieve (tfd,qvfd,qlfd) from (slfd,qtfd) in order to 
         ! use 'trbintd' at the next iteration.
          
          do k = 1, pver
             do i = 1, ncol
              ! ----------------------------------------------------- ! 
              ! Compute the condensate 'qlfd' in the updated profiles !
              ! ----------------------------------------------------- !  
              ! Option.1 : Assume grid-mean condensate is homogeneously diffused by the moist turbulence scheme.
              !            This should bs used if 'pseudodiff = .false.' in vertical_diffusion.F90.
              ! Modification : Need to be check whether below is correct in the presence of ice, qi.       
              !                I should understand why the variation of ice, qi is neglected during diffusion.
                templ     = ( slfd(k,i) - g*z(k,i) ) / cpair
                call qsat( templ, pmid(k,i), es, qs)
                ep2       =  .622_r8 
                temps     =   templ + ( qtfd(k,i) - qs ) / ( cpair / latvap + latvap * qs / ( rair * templ**2 ) )
                call qsat( temps, pmid(k,i), es, qs)
                qlfd(k,i) =   max( qtfd(k,i) - qi(k,i) - qs ,0._r8 )
              ! Option.2 : Assume condensate is not diffused by the moist turbulence scheme. 
              !            This should bs used if 'pseudodiff = .true.'  in vertical_diffusion.F90.       
              ! qlfd(k,i) = ql(k,i)
              ! ----------------------------- !
              ! Compute the other 'qvfd, tfd' ! 
              ! ----------------------------- !
                qvfd(k,i) = max( 0._r8, qtfd(k,i) - qi(k,i) - qlfd(k,i) )
                tfd(k,i)  = ( slfd(k,i) + latvap * qlfd(k,i) + latsub * qi(k,i) - g*z(k,i)) / cpair
             end do
          end do
       endif

     ! Debug 
     ! icol = phys_debug_col(lchnk) 
     ! if( icol > 0 .and. get_nstep() .ge. 1 ) then
     !     write(*,*) ' '
     !     write(*,*) 'eddy_diff debug at the end of iteration' 
     !     write(*,*) 't,     qv,     ql,     cld,     u,     v'
     !     do k = pver-3, pver
     !        write (*,*) k, tfd(k,icol), qvfd(k,icol), qlfd(k,icol), cldn(k,icol), ufd(k,icol), vfd(k,icol)
     !     end do
     ! endif
     ! Debug

    end do  ! End of 'iturb' iteration

    kvq(:,:ncol) = kvh_out(:,:ncol)

  ! Compute 'wstar' within the PBL for use in the future convection scheme.
    do i = 1, ncol
       if( ipbl(i) .eq. 1._r8 ) then 
           wstarPBL(i) = max( 0._r8, wstar(1,i) )
       else
           wstarPBL(i) = 0._r8
       endif
    end do


  ! Writing for detailed diagnostic analysis of UW moist PBL scheme 
  ! see CAM: atm/src/physics/cam/eddy_diff.F90, line782-859

    end subroutine compute_eddy_diff


! Purpose: Calculate buoyancy coefficients at all interfaces including surface. 
! Also, computes the profiles of ( sl,qt,n2,s2,ri ). Note that (n2,s2,ri) are 
! defined at each interfaces except surface.
    subroutine trbintd(pver, ncol,    z,       u,     v,     t,     pmid,  &
                       s2,   n2,      ri,      zi,    pi,    cld,   qt,    &
                       qv,   ql,      qi,      sfi,   sfuh,  sflh,  sl,    &
                       slv,  slslope, qtslope, chs,   chu,   cms,   cmu    )

    implicit none
! io
    integer,  intent(in)  :: pver                             ! Number of atmospheric layers   
    integer,  intent(in)  :: ncol                             ! Number of atmospheric columns
    real(r8), intent(in)  :: z(pver, ncol)                    ! Layer mid-point height above surface [ m ]
    real(r8), intent(in)  :: u(pver, ncol)                    ! Layer mid-point u [ m/s ]
    real(r8), intent(in)  :: v(pver, ncol)                    ! Layer mid-point v [ m/s ]
    real(r8), intent(in)  :: t(pver, ncol)                    ! Layer mid-point temperature [ K ]
    real(r8), intent(in)  :: pmid(pver, ncol)                 ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: zi(pver+1, ncol)                 ! Interface height [ m ]
    real(r8), intent(in)  :: pi(pver+1, ncol)                 ! Interface pressure [ Pa ]
    real(r8), intent(in)  :: cld(pver, ncol)                  ! Stratus fraction
    real(r8), intent(in)  :: qv(pver, ncol)                   ! Water vapor specific humidity [ kg/kg ]
    real(r8), intent(in)  :: ql(pver, ncol)                   ! Liquid water specific humidity [ kg/kg ]
    real(r8), intent(in)  :: qi(pver, ncol)                   ! Ice water specific humidity [ kg/kg ]

    real(r8), intent(out) :: s2(pver, ncol)                   ! Interfacial ( except surface ) shear squared [ s-2 ]
    real(r8), intent(out) :: n2(pver, ncol)                   ! Interfacial ( except surface ) buoyancy frequency [ s-2 ]
    real(r8), intent(out) :: ri(pver, ncol)                   ! Interfacial ( except surface ) Richardson number, 'n2/s2'
    real(r8), intent(out) :: qt(pver, ncol)                   ! Total specific humidity [ kg/kg ]
    real(r8), intent(out) :: sfi(pver+1, ncol)                ! Interfacial layer saturation fraction [ fraction ]
    real(r8), intent(out) :: sfuh(pver, ncol)                 ! Saturation fraction in upper half-layer [ fraction ]
    real(r8), intent(out) :: sflh(pver, ncol)                 ! Saturation fraction in lower half-layer [ fraction ]
    real(r8), intent(out) :: sl(pver, ncol)                   ! Liquid water static energy [ J/kg ] 
    real(r8), intent(out) :: slv(pver, ncol)                  ! Liquid water virtual static energy [ J/kg ]
    real(r8), intent(out) :: chu(pver+1, ncol)                ! Heat buoyancy coef for dry states at all interfaces, finally.
    real(r8), intent(out) :: chs(pver+1, ncol)                ! heat buoyancy coef for sat states at all interfaces, finally.
    real(r8), intent(out) :: cmu(pver+1, ncol)                ! Moisture buoyancy coef for dry states at all interfaces, finally.
    real(r8), intent(out) :: cms(pver+1, ncol)                ! Moisture buoyancy coef for sat states at all interfaces, finally.
    real(r8), intent(out) :: slslope(pver, ncol)              ! Slope of 'sl' in each layer
    real(r8), intent(out) :: qtslope(pver, ncol)              ! Slope of 'qt' in each layer

! local 
    integer               :: i                                ! Longitude index
    integer               :: k, km1                           ! Level index
    real(r8)              :: qs(pver, ncol)                   ! Saturation specific humidity
    real(r8)              :: es(pver, ncol)                   ! Saturation vapor pressure
    real(r8)              :: gam(pver, ncol)                  ! (l/cp)*(d(qs)/dT)
    real(r8)              :: rdz                              ! 1 / (delta z) between midpoints
    real(r8)              :: dsldz                            ! 'delta sl / delta z' at interface
    real(r8)              :: dqtdz                            ! 'delta qt / delta z' at interface
    real(r8)              :: ch                               ! 'sfi' weighted ch at the interface
    real(r8)              :: cm                               ! 'sfi' weighted cm at the interface
    real(r8)              :: bfact                            ! Buoyancy factor in n2 calculations
    real(r8)              :: product                          ! Intermediate vars used to find slopes
    real(r8)              :: dsldp_a, dqtdp_a                 ! Slopes across interface above 
    real(r8)              :: dsldp_b(ncol), dqtdp_b(ncol)     ! Slopes across interface below

    ! Calculate conservative scalars (qt,sl,slv) and buoyancy coefficients at the layer mid-points.
    ! Note that 'ntop_turb = 1', 'nbot_turb = pver'
    do k = ntop_turb, nbot_turb
        call qsat( t(k,:ncol), pmid(k,:ncol), es(k,:ncol), qs(k,:ncol), gam=gam(k,:ncol))
        do i = 1, ncol
            qt(k,i)  = qv(k,i) + ql(k,i) + qi(k,i) 
            sl(k,i)  = cpair * t(k,i) + g * z(k,i) - latvap * ql(k,i) - latsub * qi(k,i)
            slv(k,i) = sl(k,i) * ( 1._r8 + zvir * qt(k,i) )
            ! Thermodynamic coefficients for buoyancy flux - in this loop these are
            ! calculated at mid-points; later,  they will be averaged to interfaces,
            ! where they will ultimately be used.  At the surface, the coefficients
            ! are taken from the lowest mid point.
            bfact    = g / ( t(k,i) * ( 1._r8 + zvir * qv(k,i) - ql(k,i) - qi(k,i) ) )
            chu(k,i) = ( 1._r8 + zvir * qt(k,i) ) * bfact / cpair
            chs(k,i) = ( ( 1._r8 + ( 1._r8 + zvir ) * gam(k,i) * cpair * t(k,i) / latvap ) / ( 1._r8 + gam(k,i) ) ) * bfact / cpair
            cmu(k,i) = zvir * bfact * t(k,i)
            cms(k,i) = latvap * chs(k,i)  -  bfact * t(k,i)
        end do
    end do

    do i = 1, ncol
       chu(pver+1,i) = chu(pver,i)
       chs(pver+1,i) = chs(pver,i)
       cmu(pver+1,i) = cmu(pver,i)
       cms(pver+1,i) = cms(pver,i)
    end do

    ! Compute slopes of conserved variables sl, qt within each layer k. 
    ! 'a' indicates the 'above' gradient from layer k-1 to layer k and 
    ! 'b' indicates the 'below' gradient from layer k   to layer k+1.
    ! We take a smaller (in absolute value)  of these gradients as the
    ! slope within layer k. If they have opposite signs,   gradient in 
    ! layer k is taken to be zero. I should re-consider whether   this
    ! profile reconstruction is the best or not.
    ! This is similar to the profile reconstruction used in the UWShCu. 

    do i = 1, ncol
     ! Slopes at endpoints determined by extrapolation
       slslope(pver,i) = ( sl(pver,i) - sl(pver-1,i) ) / ( pmid(pver,i) - pmid(pver-1,i) )
       qtslope(pver,i) = ( qt(pver,i) - qt(pver-1,i) ) / ( pmid(pver,i) - pmid(pver-1,i) )
       slslope(1,i)    = ( sl(2,i) - sl(1,i) ) / ( pmid(2,i) - pmid(1,i) )
       qtslope(1,i)    = ( qt(2,i) - qt(1,i) ) / ( pmid(2,i) - pmid(1,i) )
       dsldp_b(i)      = slslope(1,i)
       dqtdp_b(i)      = qtslope(1,i)
    end do

    do k = 2, pver - 1
       do i = 1, ncol
          dsldp_a    = dsldp_b(i)
          dqtdp_a    = dqtdp_b(i)
          dsldp_b(i) = ( sl(k+1,i) - sl(k,i) ) / ( pmid(k+1,i) - pmid(k,i) )
          dqtdp_b(i) = ( qt(k+1,i) - qt(k,i) ) / ( pmid(k+1,i) - pmid(k,i) )
          product    = dsldp_a * dsldp_b(i)
          if( product .le. 0._r8 ) then 
              slslope(k,i) = 0._r8
          else if( product .gt. 0._r8 .and. dsldp_a .lt. 0._r8 ) then 
              slslope(k,i) = max( dsldp_a, dsldp_b(i) )
          else if( product .gt. 0._r8 .and. dsldp_a .gt. 0._r8 ) then 
              slslope(k,i) = min( dsldp_a, dsldp_b(i) )
          end if
          product = dqtdp_a*dqtdp_b(i)
          if( product .le. 0._r8 ) then 
              qtslope(k,i) = 0._r8
          else if( product .gt. 0._r8 .and. dqtdp_a .lt. 0._r8 ) then 
              qtslope(k,i) = max( dqtdp_a, dqtdp_b(i) )
          else if( product .gt. 0._r8 .and. dqtdp_a .gt. 0._r8 ) then 
              qtslope(k,i) = min( dqtdp_a, dqtdp_b(i) )
          end if
       end do ! i
    end do ! k

    !  Compute saturation fraction at the interfacial layers for use in buoyancy
    !  flux computation.

    call sfdiag( pver   , ncol    , qt      , ql      , sl      ,           & 
                 pi     , pmid    , zi      , cld     , sfi     , sfuh    , &
                 sflh   , slslope , qtslope )

    ! Calculate buoyancy coefficients at all interfaces (1:pver+1) and (n2,s2,ri) 
    ! at all interfaces except surface. Note 'nbot_turb = pver', 'ntop_turb = 1'.
    ! With the previous definition of buoyancy coefficients at the surface, the 
    ! resulting buoyancy coefficients at the top and surface interfaces becomes 
    ! identical to the buoyancy coefficients at the top and bottom layers. Note 
    ! that even though the dimension of (s2,n2,ri) is 'pver',  they are defined
    ! at interfaces ( not at the layer mid-points ) except the surface. 

    do k = nbot_turb, ntop_turb + 1, -1
       km1 = k - 1
       do i = 1, ncol
          rdz      = 1._r8 / ( z(km1,i) - z(k,i) )
          dsldz    = ( sl(km1,i) - sl(k,i) ) * rdz
          dqtdz    = ( qt(km1,i) - qt(k,i) ) * rdz 
          chu(k,i) = ( chu(km1,i) + chu(k,i) ) * 0.5_r8
          chs(k,i) = ( chs(km1,i) + chs(k,i) ) * 0.5_r8
          cmu(k,i) = ( cmu(km1,i) + cmu(k,i) ) * 0.5_r8
          cms(k,i) = ( cms(km1,i) + cms(k,i) ) * 0.5_r8
          ch       = chu(k,i) * ( 1._r8 - sfi(k,i) ) + chs(k,i) * sfi(k,i)
          cm       = cmu(k,i) * ( 1._r8 - sfi(k,i) ) + cms(k,i) * sfi(k,i)
          n2(k,i)  = ch * dsldz +  cm * dqtdz
          s2(k,i)  = ( ( u(km1,i) - u(k,i) )**2 + ( v(km1,i) - v(k,i) )**2) * rdz**2
          s2(k,i)  = max( ntzero, s2(k,i) )
          ri(k,i)  = n2(k,i) / s2(k,i)
       end do
    end do 
    do i = 1, ncol
       n2(1,i) = n2(2,i)
       s2(1,i) = s2(2,i)
       ri(1,i) = ri(2,i)
    end do
    return

    end subroutine trbintd


! Purpose: Interface for calculating saturation fractions  at upper and
!          lower-half layers, & interfaces for use by turbulence scheme.
!          The computed saturation fractions are repeatedly
!          used to compute buoyancy coefficients in'trbintd' & 'caleddy'.
    subroutine sfdiag( pver    , ncol    , qt      , ql      , sl      ,           &
                       pi      , pm      , zi      , cld     , sfi     , sfuh    , &
                       sflh    , slslope , qtslope )

    implicit none       

    integer,  intent(in)  :: pver                ! Number of atmospheric layers   
    integer,  intent(in)  :: ncol                ! Number of atmospheric columns   
    real(r8), intent(in)  :: sl(pver,ncol)      ! Liquid water static energy [ J/kg ]
    real(r8), intent(in)  :: qt(pver,ncol)      ! Total water specific humidity [ kg/kg ]
    real(r8), intent(in)  :: ql(pver,ncol)      ! Liquid water specific humidity [ kg/kg ]
    real(r8), intent(in)  :: pi(pver+1,ncol)    ! Interface pressures [ Pa ]
    real(r8), intent(in)  :: pm(pver,ncol)      ! Layer mid-point pressures [ Pa ]
    real(r8), intent(in)  :: zi(pver+1,ncol)    ! Interface heights [ m ]
    real(r8), intent(in)  :: cld(pver,ncol)     ! Stratiform cloud fraction [ fraction ]
    real(r8), intent(in)  :: slslope(pver,ncol) ! Slope of 'sl' in each layer
    real(r8), intent(in)  :: qtslope(pver,ncol) ! Slope of 'qt' in each layer
    real(r8), intent(out) :: sfi(pver+1,ncol)   ! Interfacial layer saturation fraction [ fraction ]
    real(r8), intent(out) :: sfuh(pver,ncol)    ! Saturation fraction in upper half-layer [ fraction ]
    real(r8), intent(out) :: sflh(pver,ncol)    ! Saturation fraction in lower half-layer [ fraction ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer               :: i                   ! Longitude index
    integer               :: k                   ! Vertical index
    integer               :: km1                 ! k-1
    integer               :: status              ! Status returned by function calls
    real(r8)              :: sltop, slbot        ! sl at top/bot of grid layer
    real(r8)              :: qttop, qtbot        ! qt at top/bot of grid layer
    real(r8)              :: tltop, tlbot  ! Liquid water temperature at top/bot of grid layer
    real(r8)              :: qxtop, qxbot        ! Sat excess at top/bot of grid layer
    real(r8)              :: qxm                 ! Sat excess at midpoint
    real(r8)              :: es               ! Saturation vapor pressure
    real(r8)              :: qs               ! Saturation spec. humidity
    real(r8)              :: cldeff(pver,ncol)  ! Effective Cloud Fraction [ fraction ]

    sfi(:,1:ncol)    = 0._r8
    sfuh(:,1:ncol)   = 0._r8
    sflh(:,1:ncol)   = 0._r8
    cldeff(:,1:ncol) = 0._r8

    select case (sftype)
    case ('d')
       ! ----------------------------------------------------------------------- !
       ! Simply use the given stratus fraction ('horizontal' cloud partitioning) !
       ! ----------------------------------------------------------------------- !
       do k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          do i = 1, ncol
             sfuh(k,i) = cld(k,i)
             sflh(k,i) = cld(k,i)
             sfi(k,i)  = 0.5_r8 * ( sflh(km1,i) + min( sflh(km1,i), sfuh(k,i) ) )
          end do
       end do
       do i = 1, ncol
          sfi(pver+1,i) = sflh(pver,i) 
       end do
    case ('l')
       ! ------------------------------------------ !
       ! Use modified stratus fraction partitioning !
       ! ------------------------------------------ !
       do k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          do i = 1, ncol
             cldeff(k,i) = cld(k,i)
             sfuh(k,i)   = cld(k,i)
             sflh(k,i)   = cld(k,i)
             if( ql(k,i) .lt. qmin ) then
                 sfuh(k,i) = 0._r8
                 sflh(k,i) = 0._r8
             end if
           ! Modification : The contribution of ice should be carefully considered.
             if( choice_evhc .eq. 'ramp' .or. choice_radf .eq. 'ramp' ) then 
                 cldeff(k,i) = cld(k,i) * min( ql(k,i) / qmin, 1._r8 )
                 sfuh(k,i)   = cldeff(k,i)
                 sflh(k,i)   = cldeff(k,i)
             elseif( choice_evhc .eq. 'maxi' .or. choice_radf .eq. 'maxi' ) then 
                 cldeff(k,i) = cld(k,i)
                 sfuh(k,i)   = cldeff(k,i)
                 sflh(k,i)   = cldeff(k,i)
             endif
           ! At the stratus top, take the minimum interfacial saturation fraction
             sfi(k,i) = 0.5_r8 * ( sflh(km1,i) + min( sfuh(k,i), sflh(km1,i) ) )
           ! Modification : Currently sfi at the top and surface interfaces are set to be zero.
           !                Also, sfuh and sflh in the top model layer is set to be zero.
           !                However, I may need to set 
           !                         do i = 1, ncol
           !                            sfi(pver+1,i) = sflh(pver,i) 
           !                         end do
           !                for treating surface-based fog. 
           ! OK. I added below block similar to the other cases.
          end do
       end do
       do i = 1, ncol
          sfi(pver+1,i) = sflh(pver,i)
       end do
    case ('u')
       ! ------------------------------------------------------------------------- !
       ! Use unsaturated buoyancy - since sfi, sfuh, sflh have already been zeroed !
       ! nothing more need be done for this case.                                  !
       ! ------------------------------------------------------------------------- !
    case ('z')
       ! ------------------------------------------------------------------------- !
       ! Calculate saturation fraction based on whether the air just above or just !
       ! below the interface is saturated, i.e. with vertical cloud partitioning.  !
       ! The saturation fraction of the interfacial layer between mid-points k and !
       ! k+1 is computed by averaging the saturation fraction   of the half-layers !
       ! above and below the interface,  with a special provision   for cloud tops !
       ! (more cloud in the half-layer below than in the half-layer above).In each !
       ! half-layer, vertical partitioning of  cloud based on the slopes diagnosed !
       ! above is used.     Loop down through the layers, computing the saturation !
       ! fraction in each half-layer (sfuh for upper half, sflh for lower half).   !
       ! Once sfuh(k,i) is computed, use with sflh(k-1,i) to determine  saturation !
       ! fraction sfi(k,i) for interfacial layer k-0.5.                            !
       ! This is 'not' chosen for full consistent treatment of stratus fraction in !
       ! all physics schemes.                                                      !
       ! ------------------------------------------------------------------------- !
       do k = ntop_turb + 1, nbot_turb
          km1 = k - 1
          do i = 1, ncol
           ! Compute saturation excess at the mid-point of layer k
             sltop    = sl(k,i) + slslope(k,i) * ( pi(k,i) - pm(k,i) )      
             qttop    = qt(k,i) + qtslope(k,i) * ( pi(k,i) - pm(k,i) )
             tltop = ( sltop - g * zi(k,i) ) / cpair 
             call qsat( tltop, pi(k,i), es, qs)
             qxtop    = qttop - qs
             slbot    = sl(k,i) + slslope(k,i) * ( pi(k+1,i) - pm(k,i) )      
             qtbot    = qt(k,i) + qtslope(k,i) * ( pi(k+1,i) - pm(k,i) )
             tlbot = ( slbot - g * zi(k+1,i) ) / cpair 
             call qsat( tlbot, pi(k+1,i), es, qs)
             qxbot    = qtbot - qs
             qxm      = qxtop + ( qxbot - qxtop ) * ( pm(k,i) - pi(k,i) ) / ( pi(k+1,i) - pi(k,i) )
           ! Find the saturation fraction sfuh(k,i) of the upper half of layer k.
             if( ( qxtop .lt. 0._r8 ) .and. ( qxm .lt. 0._r8 ) ) then
                   sfuh(k,i) = 0._r8 
             else if( ( qxtop .gt. 0._r8 ) .and. ( qxm .gt. 0._r8 ) ) then
                   sfuh(k,i) = 1._r8  
             else ! Either qxm < 0 and qxtop > 0 or vice versa
                   sfuh(k,i) = max( qxtop, qxm ) / abs( qxtop - qxm )
             end if
           ! Combine with sflh(i) (still for layer k-1) to get interfac layer saturation fraction
             sfi(k,i) = 0.5_r8 * ( sflh(k-1,i) + min( sflh(k-1,i), sfuh(k,i) ) )
           ! Update sflh to be for the lower half of layer k.             
             if( ( qxbot .lt. 0._r8 ) .and. ( qxm .lt. 0._r8 ) ) then
                   sflh(k,i) = 0._r8 
             else if( ( qxbot .gt. 0._r8 ) .and. ( qxm .gt. 0._r8 ) ) then
                   sflh(k,i) = 1._r8 
             else ! Either qxm < 0 and qxbot > 0 or vice versa
                   sflh(k,i) = max( qxbot, qxm ) / abs( qxbot - qxm )
             end if
          end do  ! i
       end do ! k
       do i = 1, ncol
          sfi(pver+1,i) = sflh(pver,i)  ! Saturation fraction in the lowest half-layer. 
       end do
    end select

  return
  end subroutine sfdiag
 

! The University of Washington Moist Turbulence Scheme:
! Purpose : This is a driver routine to compute eddy diffusion coefficients
!           for heat (sl), momentum (u, v), moisture (qt), and other  trace
!           constituents.   This scheme uses first order closure for stable
!           turbulent layers (STL). For convective layers (CL), entrainment
!           closure is used at the CL external interfaces, which is coupled
!           to the diagnosis of a CL regime mean TKE from the instantaneous
!           thermodynamic and velocity profiles.   The CLs are diagnosed by
!           extending original CL layers of moist static instability   into
!           adjacent weakly stably stratified interfaces,   stopping if the
!           stability is too strong.   This allows a realistic depiction of
!           dry convective boundary layers with a downgradient approach.
!
! NOTE:     This routine currently assumes ntop_turb = 1, nbot_turb = pver
!           ( turbulent diffusivities computed at all interior interfaces )
!           and will require modification to handle a different ntop_turb. 
    subroutine caleddy( pver         , ncol         ,                                           &
                        sl           , qt           , ql          , slv        , u            , &
                        v            , pi           , z           , zi         ,                &
                        qflx         , shflx        , slslope     , qtslope    ,                &
                        chu          , chs          , cmu         , cms        , sfuh         , &
                        sflh         , n2           , s2          , ri         , rrho         , &
                        pblh         , ustar        ,                                           &
                        kvh_in       , kvm_in       , kvh         , kvm        ,                &
                        tpert        , qpert        , qrlin       , tke        ,                & 
                        wstarent     , bprod        , sprod       , minpblh    , wpert        , &
                        tkes         , turbtype     , sm_aw       ,                             &
                        kbase_o      , ktop_o       , ncvfin_o    ,                             & 
                        kbase_mg     , ktop_mg      , ncvfin_mg   ,                             & 
                        kbase_f      , ktop_f       , ncvfin_f    ,                             & 
                        wet_CL       , web_CL       , jtbu_CL     , jbbu_CL    ,                &
                        evhc_CL      , jt2slv_CL    , n2ht_CL     , n2hb_CL    , lwp_CL       , &
                        opt_depth_CL , radinvfrac_CL, radf_CL     , wstar_CL   , wstar3fact_CL, &
                        ebrk         , wbrk         , lbrk        , ricl       , ghcl         , & 
                        shcl         , smcl         ,                                           &
                        gh_a         , sh_a         , sm_a        , ri_a       , leng         , & 
                        wcap         , pblhp        , cld         ,                             &
                        ipbl         , kpblh        , wsedl       )

    implicit none
! io
    integer,  intent(in) :: pver                      ! Number of atmospheric layers   
    integer,  intent(in) :: ncol                      ! Number of atmospheric columns   
    real(r8), intent(in) :: u(pver, ncol)             ! U wind [ m/s ]
    real(r8), intent(in) :: v(pver, ncol)             ! V wind [ m/s ]
    real(r8), intent(in) :: sl(pver, ncol)            ! Liquid water static energy, cp * T + g * z - Lv * ql - Ls * qi [ J/kg ]
    real(r8), intent(in) :: slv(pver, ncol)           ! Liquid water virtual static energy, sl * ( 1 + 0.608 * qt ) [ J/kg ]
    real(r8), intent(in) :: qt(pver, ncol)            ! Total speccific humidity  qv + ql + qi [ kg/kg ] 
    real(r8), intent(in) :: ql(pver, ncol)            ! Liquid water specific humidity [ kg/kg ]
    real(r8), intent(in) :: pi(pver+1, ncol)          ! Interface pressures [ Pa ]
    real(r8), intent(in) :: z(pver, ncol)             ! Layer midpoint height above surface [ m ]
    real(r8), intent(in) :: zi(pver+1, ncol)          ! Interface height above surface, i.e., zi(pver+1) = 0 all over the globe [ m ]
    real(r8), intent(in) :: chu(pver+1, ncol)         ! Buoyancy coeffi. unsaturated sl (heat) coef. at all interfaces.
    real(r8), intent(in) :: chs(pver+1, ncol)         ! Buoyancy coeffi. saturated sl (heat) coef. at all interfaces.
    real(r8), intent(in) :: cmu(pver+1, ncol)         ! Buoyancy coeffi. unsaturated qt (moisture) coef. at all interfaces
    real(r8), intent(in) :: cms(pver+1, ncol)         ! Buoyancy coeffi. saturated qt (moisture) coef. at all interfaces
    real(r8), intent(in) :: sfuh(pver, ncol)          ! Saturation fraction in upper half-layer [ fraction ]
    real(r8), intent(in) :: sflh(pver, ncol)          ! Saturation fraction in lower half-layer [ fraction ]
    real(r8), intent(in) :: n2(pver, ncol)            ! Interfacial (except surface) moist buoyancy frequency [ s-2 ]
    real(r8), intent(in) :: s2(pver, ncol)            ! Interfacial (except surface) shear frequency [ s-2 ]
    real(r8), intent(in) :: ri(pver, ncol)            ! Interfacial (except surface) Richardson number
    real(r8), intent(in) :: qflx(ncol)                ! Kinematic surface constituent ( water vapor ) flux [ kg/m2/s ]
    real(r8), intent(in) :: shflx(ncol)               ! Kinematic surface heat flux [ unit ? ] 
    real(r8), intent(in) :: slslope(pver, ncol)       ! Slope of 'sl' in each layer [ J/kg/Pa ]
    real(r8), intent(in) :: qtslope(pver, ncol)       ! Slope of 'qt' in each layer [ kg/kg/Pa ]
    real(r8), intent(in) :: qrlin(pver, ncol)         ! Input grid-mean LW heating rate : [ K/s ] * cpair * dp = [ W/kg*Pa ]
    real(r8), intent(in) :: wsedl(pver, ncol)         ! Sedimentation velocity of liquid stratus cloud droplet [ m/s ]
    real(r8), intent(in) :: ustar(ncol)               ! Surface friction velocity [ m/s ]
    real(r8), intent(in) :: rrho(ncol)                ! 1./bottom mid-point density. Specific volume [ m3/kg ]
    logical,  intent(in) :: wstarent                  ! Switch for choosing wstar3 entrainment parameterization
    real(r8), intent(in) :: minpblh(ncol)             ! Minimum PBL height based on surface stress [ m ]
    real(r8), intent(in) :: kvh_in(pver+1, ncol)      ! kvh saved from last timestep or last iterative step [ m2/s ] 
    real(r8), intent(in) :: kvm_in(pver+1, ncol)      ! kvm saved from last timestep or last iterative step [ m2/s ]
    real(r8), intent(in) :: cld(pver, ncol)           ! Stratus Cloud Fraction [ fraction ]

    real(r8), intent(out) :: kvh(pver+1, ncol)        ! Eddy diffusivity for heat, moisture, and tracers [ m2/s ]
    real(r8), intent(out) :: kvm(pver+1, ncol)        ! Eddy diffusivity for momentum [ m2/s ]
    real(r8), intent(out) :: pblh(ncol)               ! PBL top height [ m ]
    real(r8), intent(out) :: pblhp(ncol)              ! PBL top height pressure [ Pa ]
    real(r8), intent(out) :: tpert(ncol)              ! Convective temperature excess [ K ]
    real(r8), intent(out) :: qpert(ncol)              ! Convective humidity excess [ kg/kg ]
    real(r8), intent(out) :: wpert(ncol)              ! Turbulent velocity excess [ m/s ]
    real(r8), intent(out) :: tke(pver+1, ncol)        ! Turbulent kinetic energy [ m2/s2 ], 'tkes' at surface, pver+1.
    real(r8), intent(out) :: bprod(pver+1, ncol)      ! Buoyancy production [ m2/s3 ],     'bflxs' at surface, pver+1.
    real(r8), intent(out) :: sprod(pver+1, ncol)      ! Shear production [ m2/s3 ], (ustar(i)**3)/(vk*z(pver,i))
                                                      ! at surface, pver+1.
    integer(i4), intent(out) :: turbtype(pver+1, ncol)! Turbulence type at each interface:
                                                      ! 0. = Non turbulence interface
                                                      ! 1. = Stable turbulence interface
                                                      ! 2. = CL interior interface ( if bflxs > 0, surface is this )
                                                      ! 3. = Bottom external interface of CL
                                                      ! 4. = Top external interface of CL.
                                                      ! 5. = Double entraining CL external interface 
    real(r8), intent(out) :: sm_aw(pver+1, ncol)      ! Galperin instability function of momentum for use in the microphysics
    real(r8), intent(out) :: ipbl(ncol)               ! If 1, PBL is CL, while if 0, PBL is STL.
    real(r8), intent(out) :: kpblh(ncol)              ! Layer index containing PBL within or at the base interface

    ! Diagnostic output variables 
    real(r8) :: tkes(ncol)                            ! TKE at surface [ m2/s2 ] 
    real(r8) :: kbase_o(ncvmax, ncol)                 ! Original external base interface index of CL just after 'exacol'
    real(r8) :: ktop_o(ncvmax, ncol)                  ! Original external top  interface index of CL just after 'exacol'
    real(r8) :: ncvfin_o(ncol)                        ! Original number of CLs just after 'exacol'
    real(r8) :: kbase_mg(ncvmax, ncol)                ! kbase  just after extending-merging (after 'zisocl') but without SRCL
    real(r8) :: ktop_mg(ncvmax, ncol)                 ! ktop   just after extending-merging (after 'zisocl') but without SRCL
    real(r8) :: ncvfin_mg(ncol)                       ! ncvfin just after extending-merging (after 'zisocl') but without SRCL
    real(r8) :: kbase_f(ncvmax, ncol)                 ! Final kbase  after adding SRCL
    real(r8) :: ktop_f(ncvmax, ncol)                  ! Final ktop   after adding SRCL
    real(r8) :: ncvfin_f(ncol)                        ! Final ncvfin after adding SRCL
    real(r8) :: wet_CL(ncvmax, ncol)                  ! Entrainment rate at the CL top [ m/s ] 
    real(r8) :: web_CL(ncvmax, ncol)                  ! Entrainment rate at the CL base [ m/s ]
    real(r8) :: jtbu_CL(ncvmax, ncol)                 ! Buoyancy jump across the CL top [ m/s2 ]  
    real(r8) :: jbbu_CL(ncvmax, ncol)                 ! Buoyancy jump across the CL base [ m/s2 ]  
    real(r8) :: evhc_CL(ncvmax, ncol)                 ! Evaporative enhancement factor at the CL top
    real(r8) :: jt2slv_CL(ncvmax, ncol)               ! Jump of slv ( across two layers ) at CL top for use only in evhc [ J/kg ]
    real(r8) :: n2ht_CL(ncvmax, ncol)                 ! n2 defined at the CL top  interface
                                                      ! but using sfuh(kt)   instead of sfi(kt) [ s-2 ]
    real(r8) :: n2hb_CL(ncvmax, ncol)                 ! n2 defined at the CL base interface
                                                      ! but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
    real(r8) :: lwp_CL(ncvmax, ncol)                  ! LWP in the CL top layer [ kg/m2 ]
    real(r8) :: opt_depth_CL(ncvmax, ncol)            ! Optical depth of the CL top layer
    real(r8) :: radinvfrac_CL(ncvmax, ncol)           ! Fraction of LW radiative cooling confined in the top portion of CL
    real(r8) :: radf_CL(ncvmax, ncol)                 ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]
    real(r8) :: wstar_CL(ncvmax, ncol)                ! Convective velocity of CL including entrainment contribution finally [ m/s ]
    real(r8) :: wstar3fact_CL(ncvmax, ncol)           ! "wstar3fact" of CL. Entrainment enhancement of wstar3 (inverse)

    real(r8) :: gh_a(pver+1, ncol)                    ! Half of normalized buoyancy production, -l2n2/2e. [ no unit ]
    real(r8) :: sh_a(pver+1, ncol)                    ! Galperin instability function of heat-moisture at all interfaces [ no unit ]
    real(r8) :: sm_a(pver+1, ncol)                    ! Galperin instability function of momentum      at all interfaces [ no unit ]
    real(r8) :: ri_a(pver+1, ncol)                    ! Interfacial Richardson number                  at all interfaces [ no unit ]

    real(r8) :: ebrk(ncvmax, ncol)                    ! Net CL mean TKE [ m2/s2 ]
    real(r8) :: wbrk(ncvmax, ncol)                    ! Net CL mean normalized TKE [ m2/s2 ]
    real(r8) :: lbrk(ncvmax, ncol)                    ! Net energetic integral thickness of CL [ m ]
    real(r8) :: ricl(ncvmax, ncol)                    ! Mean Richardson number of CL ( l2n2/l2s2 )
    real(r8) :: ghcl(ncvmax, ncol)                    ! Half of normalized buoyancy production of CL                 
    real(r8) :: shcl(ncvmax, ncol)                    ! Instability function of heat and moisture of CL
    real(r8) :: smcl(ncvmax, ncol)                    ! Instability function of momentum of CL

    real(r8) :: leng(pver+1, ncol)                    ! Turbulent length scale [ m ], 0 at the surface.
    real(r8) :: wcap(pver+1, ncol)                    ! Normalized TKE [m2/s2], 'tkes/b1' at the surface and 'tke/b1' at
                                                      ! the top/bottom entrainment interfaces of CL assuming no transport.

! local
    logical :: belongcv(pver+1, ncol)                 ! True for interfaces in a CL (both interior and exterior are included)
    logical :: belongst(pver+1, ncol)                 ! True for stable turbulent layer interfaces (STL)
    logical :: in_CL                                  ! True if interfaces k,k+1 both in same CL.
    logical :: extend                                 ! True when CL is extended in zisocl
    logical :: extend_up                              ! True when CL is extended upward in zisocl
    logical :: extend_dn                              ! True when CL is extended downward in zisocl

    integer :: i                                      ! Longitude index
    integer :: k                                      ! Vertical index
    integer :: ks                                     ! Vertical index
    integer :: ncvfin(ncol)                           ! Total number of CL in column
    integer :: ncvf                                   ! Total number of CL in column prior to adding SRCL
    integer :: ncv                                    ! Index of current CL
    integer :: ncvnew                                 ! Index of added SRCL appended after regular CLs from 'zisocl'
    integer :: ncvsurf                                ! If nonzero, CL index based on surface
                                                      ! (usually 1, but can be > 1 when SRCL is based at sfc)
    integer :: kbase(ncvmax, ncol)                    ! Vertical index of CL base interface
    integer :: ktop(ncvmax, ncol)                     ! Vertical index of CL top interface
    integer :: kb, kt                                 ! kbase and ktop for current CL
    integer :: ktblw                                  ! ktop of the CL located at just below the current CL

    integer  :: ktopbl(ncol)                          ! PBL top height or interface index 
    real(r8) :: bflxs(ncol)                           ! Surface buoyancy flux [ m2/s3 ]
    real(r8) :: rcap                                  ! 'tke/ebrk' at all interfaces of CL.
                                                      ! Set to 1 at the CL entrainment interfaces
    real(r8) :: jtzm                                  ! Interface layer thickness of CL top interface [ m ]
    real(r8) :: jtsl                                  ! Jump of s_l across CL top interface [ J/kg ]
    real(r8) :: jtqt                                  ! Jump of q_t across CL top interface [ kg/kg ]
    real(r8) :: jtbu                                  ! Jump of buoyancy across CL top interface [ m/s2 ]
    real(r8) :: jtu                                   ! Jump of u across CL top interface [ m/s ]
    real(r8) :: jtv                                   ! Jump of v across CL top interface [ m/s ]
    real(r8) :: jt2slv                                ! Jump of slv ( across two layers ) at CL top for use only in evhc [ J/kg ]
    real(r8) :: radf                                  ! Buoyancy production at the CL top due to radiative cooling [ m2/s3 ]
    real(r8) :: jbzm                                  ! Interface layer thickness of CL base interface [ m ]
    real(r8) :: jbsl                                  ! Jump of s_l across CL base interface [ J/kg ]
    real(r8) :: jbqt                                  ! Jump of q_t across CL top interface [ kg/kg ]
    real(r8) :: jbbu                                  ! Jump of buoyancy across CL base interface [ m/s2 ]
    real(r8) :: jbu                                   ! Jump of u across CL base interface [ m/s ]
    real(r8) :: jbv                                   ! Jump of v across CL base interface [ m/s ]
    real(r8) :: ch                                    ! Buoyancy coefficients defined at the CL top and base interfaces
                                                      ! using CL internal
    real(r8) :: cm                                    ! sfuh(kt) and sflh(kb-1) instead of sfi(kt) and sfi(kb), respectively.
                                                      ! These are used for entrainment calculation at CL external interfaces
                                                      ! and SRCL identification.
    real(r8) :: n2ht                                  ! n2 defined at the CL top  interface
                                                      ! but using sfuh(kt)   instead of sfi(kt) [ s-2 ]
    real(r8) :: n2hb                                  ! n2 defined at the CL base interface
                                                      ! but using sflh(kb-1) instead of sfi(kb) [ s-2 ]
    real(r8) :: n2htSRCL                              ! n2 defined at the upper-half layer of SRCL.
                                                      ! This is used only for identifying SRCL.
                                                      ! n2htSRCL use SRCL internal slope sl and qt
                                                      ! as well as sfuh(kt) instead of sfi(kt) [ s-2 ]
    real(r8) :: gh                                    ! Half of normalized buoyancy production ( -l2n2/2e ) [ no unit ]
    real(r8) :: sh                                    ! Galperin instability function for heat and moisture
    real(r8) :: sm                                    ! Galperin instability function for momentum
    real(r8) :: lbulk                                 ! Depth of turbulent layer, Master length scale (not energetic length)
    real(r8) :: dzht                                  ! Thickness of top    half-layer [ m ]
    real(r8) :: dzhb                                  ! Thickness of bottom half-layer [ m ]
    real(r8) :: rootp                                 ! Sqrt(net CL-mean TKE including entrainment contribution) [ m/s ]     
    real(r8) :: evhc                                  ! Evaporative enhancement factor: (1+E)
                                                      ! with E = evap. cool. efficiency [ no unit ]
    real(r8) :: kentr                                 ! Effective entrainment diffusivity 'wet*dz', 'web*dz' [ m2/s ]
    real(r8) :: lwp                                   ! Liquid water path in the layer kt [ kg/m2 ]
    real(r8) :: opt_depth                             ! Optical depth of the layer kt [ no unit ]
    real(r8) :: radinvfrac                            ! Fraction of LW cooling in the layer kt
                                                      ! concentrated at the CL top [ no unit ]
    real(r8) :: wet                                   ! CL top entrainment rate [ m/s ]
    real(r8) :: web                                   ! CL bot entrainment rate [ m/s ]. Set to zero if CL is based at surface.
    real(r8) :: vyt                                   ! n2ht/n2 at the CL top  interface
    real(r8) :: vyb                                   ! n2hb/n2 at the CL base interface
    real(r8) :: vut                                   ! Inverse Ri (=s2/n2) at the CL top  interface
    real(r8) :: vub                                   ! Inverse Ri (=s2/n2) at the CL base interface
    real(r8) :: fact                                  ! Factor relating TKE generation to entrainment [ no unit ]
    real(r8) :: trma                                  ! Intermediate variables used for solving quadratic ( for gh from ri )
    real(r8) :: trmb                                  ! and cubic equations ( for ebrk: the net CL mean TKE )
    real(r8) :: trmc                                  !
    real(r8) :: trmp                                  !
    real(r8) :: trmq                                  !
    real(r8) :: qq                                    ! 
    real(r8) :: det                                   !
    real(r8) :: gg                                    ! Intermediate variable used for calculating stability functions of
                                                      ! SRCL or SBCL based at the surface with bflxs > 0.
    real(r8) :: dzhb5                                 ! Half thickness of the bottom-most layer of current CL regime
    real(r8) :: dzht5                                 ! Half thickness of the top-most layer of adjacent CL regime
                                                      ! just below current CL
    real(r8) :: qrlw(pver, ncol)                      ! Local grid-mean LW heating rate : [K/s] * cpair * dp = [ W/kg*Pa ]

    real(r8) :: cldeff(pver, ncol)                    ! Effective stratus fraction
    real(r8) :: qleff                                 ! Used for computing evhc
    real(r8) :: tunlramp                              ! Ramping tunl
    real(r8) :: leng_imsi                             ! For Kv = max(Kv_STL, Kv_entrain)
    real(r8) :: tke_imsi                              !
    real(r8) :: kvh_imsi                              !
    real(r8) :: kvm_imsi                              !
    real(r8) :: alph4exs                              ! For extended stability function in the stable regime
    real(r8) :: ghmin                                 !   

    real(r8) :: sedfact                               ! For 'sedimentation-entrainment feedback' 

    ! Local variables specific for 'wstar' entrainment closure
    real(r8) :: cet                                   ! Proportionality coefficient between wet and wstar3
    real(r8) :: ceb                                   ! Proportionality coefficient between web and wstar3
    real(r8) :: wstar                                 ! Convective velocity for CL [ m/s ]
    real(r8) :: wstar3                                ! Cubed convective velocity for CL [ m3/s3 ]
    real(r8) :: wstar3fact                            ! 1/(relative change of wstar^3 by entrainment)
    real(r8) :: rmin                                  ! sqrt(p)
    real(r8) :: fmin                                  ! f(rmin), where f(r) = r^3 - 3*p*r - 2q
    real(r8) :: rcrit                                 ! ccrit*wstar
    real(r8) :: fcrit                                 ! f(rcrit)
    logical     noroot                                ! True if f(r) has no root r > rcrit

    ! Option: Turn-off LW radiative-turbulence interaction in PBL scheme
    !         by setting qrlw = 0.  Logical parameter 'set_qrlzero'  was
    !         defined in the first part of 'eddy_diff.F90' module. 

    if( set_qrlzero ) then
        qrlw(:,:) = 0._r8
    else
        qrlw(:pver,:ncol) = qrlin(:pver,:ncol)
    endif

    ! Define effective stratus fraction using the grid-mean ql.
    ! Modification : The contribution of ice should be carefully considered.
    !                This should be done in combination with the 'qrlw' and
    !                overlapping assumption of liquid and ice stratus. 

    do k = 1, pver
       do i = 1, ncol
          if( choice_evhc .eq. 'ramp' .or. choice_radf .eq. 'ramp' ) then 
              cldeff(k,i) = cld(k,i) * min( ql(k,i) / qmin, 1._r8 )
          else
              cldeff(k,i) = cld(k,i)
          endif
       end do
    end do

    ! For an extended stability function in the stable regime, re-define
    ! alph4exe and ghmin. This is for future work.

    if( ricrit .eq. 0.19_r8 ) then
        alph4exs = alph4
        ghmin    = -3.5334_r8
    elseif( ricrit .gt. 0.19_r8 ) then
        alph4exs = -2._r8 * b1 * alph2 / ( alph3 - 2._r8 * b1 * alph5 ) / ricrit
        ghmin    = -1.e10_r8
    else
        call endrun('CALEDDY Error: ricrit should be larger than 0.19 in UW PBL')
    endif

    ! Initialization of Diagnostic Output
    do i = 1, ncol
       wet_CL(:ncvmax,i)        = 0._r8
       web_CL(:ncvmax,i)        = 0._r8
       jtbu_CL(:ncvmax,i)       = 0._r8
       jbbu_CL(:ncvmax,i)       = 0._r8
       evhc_CL(:ncvmax,i)       = 0._r8
       jt2slv_CL(:ncvmax,i)     = 0._r8
       n2ht_CL(:ncvmax,i)       = 0._r8
       n2hb_CL(:ncvmax,i)       = 0._r8                    
       lwp_CL(:ncvmax,i)        = 0._r8
       opt_depth_CL(:ncvmax,i)  = 0._r8
       radinvfrac_CL(:ncvmax,i) = 0._r8
       radf_CL(:ncvmax,i)       = 0._r8
       wstar_CL(:ncvmax,i)      = 0._r8          
       wstar3fact_CL(:ncvmax,i) = 0._r8
       ricl(:ncvmax,i)          = 0._r8
       ghcl(:ncvmax,i)          = 0._r8
       shcl(:ncvmax,i)          = 0._r8
       smcl(:ncvmax,i)          = 0._r8
       ebrk(:ncvmax,i)          = 0._r8
       wbrk(:ncvmax,i)          = 0._r8
       lbrk(:ncvmax,i)          = 0._r8
       gh_a(:pver+1,i)          = 0._r8
       sh_a(:pver+1,i)          = 0._r8
       sm_a(:pver+1,i)          = 0._r8
       ri_a(:pver+1,i)          = 0._r8
       sm_aw(:pver+1,i)         = 0._r8
       ipbl(i)                  = 0._r8
       kpblh(i)                 = real(pver,r8)
    end do  

    ! kvh and kvm are stored over timesteps in 'vertical_diffusion.F90' and 
    ! passed in as kvh_in and kvm_in.  However,  at the first timestep they
    ! need to be computed and these are done just before calling 'caleddy'.   
    ! kvm and kvh are also stored over iterative time step in the first part
    ! of 'eddy_diff.F90'

    ! Initialize kvh and kvm to zero
    kvh(:,:) = 0._r8
    kvm(:,:) = 0._r8
    
    ! Zero diagnostic quantities for the new diffusion step.
    wcap(:,:) = 0._r8
    leng(:,:) = 0._r8
    tke(:,:)  = 0._r8
    turbtype(:,:) = 0

    ! Initialize 'bprod' [ m2/s3 ] and 'sprod' [ m2/s3 ] at all interfaces.
    ! Note this initialization is a hybrid initialization since 'n2' [s-2] and 's2' [s-2]
    ! are calculated from the given current initial profile, while 'kvh_in' [m2/s] and 
    ! 'kvm_in' [m2/s] are from the previous iteration or previous time step.
    ! This initially guessed 'bprod' and 'sprod' will be updated at the end of this 
    ! 'caleddy' subroutine for diagnostic output.
    ! This computation of 'brpod,sprod' below is necessary for wstar-based entrainment closure.
    do k = 2, pver
       do i = 1, ncol
            bprod(k,i) = -kvh_in(k,i) * n2(k,i)
            sprod(k,i) =  kvm_in(k,i) * s2(k,i)
       end do
    end do

    ! Set 'bprod' and 'sprod' at top and bottom interface.
    ! In calculating 'surface' (actually lowest half-layer) buoyancy flux,
    ! 'chu' at surface is defined to be the same as 'chu' at the mid-point
    ! of lowest model layer (pver) at the end of 'trbind'. The same is for
    ! the other buoyancy coefficients.  'sprod(pver+1,i)'  is defined in a
    ! consistent way as the definition of 'tkes' in the original code.
    ! ( Important Option ) If I want to isolate surface buoyancy flux from
    ! the other parts of CL regimes energetically even though bflxs > 0,
    ! all I should do is to re-define 'bprod(pver+1,i)=0' in the below 'do'
    ! block. Additionally for merging test of extending SBCL based on 'l2n2'
    ! in 'zisocl', I should use 'l2n2 = - wint / sh'  for similar treatment
    ! as previous code. All other parts of the code  are fully consistently
    ! treated by these change only.
    ! My future general convection scheme will use bflxs(i).
    do i = 1, ncol
       bprod(1,i) = 0._r8 ! Top interface
       sprod(1,i) = 0._r8 ! Top interface
       ch = chu(pver+1,i) * ( 1._r8 - sflh(pver,i) ) + chs(pver+1,i) * sflh(pver,i)   
       cm = cmu(pver+1,i) * ( 1._r8 - sflh(pver,i) ) + cms(pver+1,i) * sflh(pver,i)   
       bflxs(i) = ch * shflx(i) * rrho(i) + cm * qflx(i) * rrho(i)
       if( choice_tkes .eq. 'ibprod' ) then
           bprod(pver+1,i) = bflxs(i)
       else
           bprod(pver+1,i) = 0._r8
       endif
       sprod(pver+1,i) = (ustar(i)**3)/(vk*z(pver,i))
    end do

    ! Initially identify CL regimes in 'exacol'
    !    ktop  : Interface index of the CL top  external interface
    !    kbase : Interface index of the CL base external interface
    !    ncvfin: Number of total CLs
    ! Note that if surface buoyancy flux is positive ( bflxs = bprod(pver+1,i) > 0 ),
    ! surface interface is identified as an internal interface of CL. However, even
    ! though bflxs <= 0, if 'pver' interface is a CL internal interface (ri(pver)<0),
    ! surface interface is identified as an external interface of CL. If bflxs =< 0 
    ! and ri(pver) >= 0, then surface interface is identified as a stable turbulent
    ! intereface (STL) as shown at the end of 'caleddy'. Even though a 'minpblh' is
    ! passed into 'exacol', it is not used in the 'exacol'.

    call exacol( pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin )

    ! Diagnostic output of CL interface indices before performing 'extending-merging'
    ! of CL regimes in 'zisocl'
    do i = 1, ncol
    do k = 1, ncvmax
       kbase_o(k,i) = real(kbase(k,i),r8)
       ktop_o(k,i)  = real(ktop(k,i),r8) 
       ncvfin_o(i)  = real(ncvfin(i),r8)
    end do
    end do 

    ! Perform calculation for each column
    do i = 1, ncol
       ! Define Surface Interfacial Layer TKE, 'tkes'.
       ! In the current code, 'tkes' is used as representing TKE of surface interfacial
       ! layer (low half-layer of surface-based grid layer). In the code, when bflxs>0,
       ! surface interfacial layer is assumed to be energetically  coupled to the other
       ! parts of the CL regime based at the surface. In this sense, it is conceptually
       ! more reasonable to include both 'bprod' and 'sprod' in the definition of 'tkes'.
       ! Since 'tkes' cannot be negative, it is lower bounded by small positive number. 
       ! Note that inclusion of 'bprod' in the definition of 'tkes' may increase 'ebrk'
       ! and 'wstar3', and eventually, 'wet' at the CL top, especially when 'bflxs>0'.
       ! This might help to solve the problem of too shallow PBLH over the overcast Sc
       ! regime. If I want to exclude 'bprod(pver+1,i)' in calculating 'tkes' even when
       ! bflxs > 0, all I should to do is to set 'bprod(pver+1,i) = 0' in the above 
       ! initialization 'do' loop (explained above), NOT changing the formulation of
       ! tkes(i) in the below block. This is because for consistent treatment in the 
       ! other parts of the code also.
  
     ! tkes(i) = (b1*vk*z(pver,i)*sprod(pver+1,i))**(2._r8/3._r8)
       tkes(i) = max(b1*vk*z(pver,i)*(bprod(pver+1,i)+sprod(pver+1,i)), 1.e-7_r8)**(2._r8/3._r8)
       tkes(i) = min(tkes(i), tkemax)
       tke(pver+1,i)  = tkes(i)
       wcap(pver+1,i) = tkes(i)/b1

       ! Extend and merge the initially identified CLs, relabel the CLs, and calculate
       ! CL internal mean energetics and stability functions in 'zisocl'. 
       ! The CL nearest to the surface is CL(1) and the CL index, ncv, increases 
       ! with height. The following outputs are from 'zisocl'. Here, the dimension
       ! of below outputs are (ncvmax, ncol) (except the 'ncvfin(ncol)' and 
       ! 'belongcv(pver+1,ncol)) and 'ncv' goes from 1 to 'ncvfin'. 
       ! For 'ncv = ncvfin+1, ncvmax', below output are already initialized to be zero. 
       !      ncvfin       : Total number of CLs
       !      kbase(ncv)   : Base external interface index of CL
       !      ktop         : Top  external interface index of CL
       !      belongcv     : True if the interface (either internal or external) is CL  
       !      ricl         : Mean Richardson number of internal CL
       !      ghcl         : Normalized buoyancy production '-l2n2/2e' [no unit] of internal CL
       !      shcl         : Galperin instability function of heat-moisture of internal CL
       !      smcl         : Galperin instability function of momentum of internal CL
       !      lbrk, <l>int : Thickness of (energetically) internal CL (lint, [m])
       !      wbrk, <W>int : Mean normalized TKE of internal CL  ([m2/s2])
       !      ebrk, <e>int : Mean TKE of internal CL (b1*wbrk,[m2/s2])
       ! The ncvsurf is an identifier saying which CL regime is based at the surface.
       ! If 'ncvsurf=1', then the first CL regime is based at the surface. If surface
       ! interface is not a part of CL (neither internal nor external), 'ncvsurf = 0'.
       ! After identifying and including SRCLs into the normal CL regimes (where newly
       ! identified SRCLs are simply appended to the normal CL regimes using regime 
       ! indices of 'ncvfin+1','ncvfin+2' (as will be shown in the below SRCL part),..
       ! where 'ncvfin' is the final CL regime index produced after extending-merging 
       ! in 'zisocl' but before adding SRCLs), if any newly identified SRCL (e.g., 
       ! 'ncvfin+1') is based at surface, then 'ncvsurf = ncvfin+1'. Thus 'ncvsurf' can
       ! be 0, 1, or >1. 'ncvsurf' can be a useful diagnostic output.   

       ncvsurf = 0

       if( ncvfin(i) .gt. 0 ) then 
           call zisocl( ncol   , pver     , i        ,           &
                        z      , zi       , n2       , s2      , & 
                        bprod  , sprod    , bflxs    , tkes    , &
                        ncvfin , kbase    , ktop     , belongcv, &
                        ricl   , ghcl     , shcl     , smcl    , & 
                        lbrk   , wbrk     , ebrk     ,           & 
                        extend , extend_up, extend_dn )
           if( kbase(1,i) .eq. pver + 1 ) ncvsurf = 1
       else
           belongcv(:,i) = .false.
       endif

       ! Diagnostic output after finishing extending-merging process in 'zisocl'
       ! Since we are adding SRCL additionally, we need to print out these here.

       do k = 1, ncvmax
          kbase_mg(k,i) = real(kbase(k,i),r8)
          ktop_mg(k,i)  = real(ktop(k,i),r8) 
          ncvfin_mg(i)  = real(ncvfin(i),r8)
       end do 

       ! ----------------------- !
       ! Identification of SRCLs !
       ! ----------------------- !

      ! Modification : This cannot identify the 'cirrus' layer due to the condition of
      !                ql(k,i) .gt. qmin. This should be modified in future to identify
      !                a single thin cirrus layer.  
      !                Instead of ql, we may use cldn in future, including ice 
      !                contribution.

       ! ------------------------------------------------------------------------------ !
       ! Find single-layer radiatively-driven cloud-topped convective layers (SRCLs).   !
       ! SRCLs extend through a single model layer k, with entrainment at the top and   !
       ! bottom interfaces, unless bottom interface is the surface.                     !
       ! The conditions for an SRCL is identified are:                                  ! 
       !                                                                                !
       !   1. Cloud in the layer, k : ql(k,i) .gt. qmin = 1.e-5 [ kg/kg ]               !
       !   2. No cloud in the above layer (else assuming that some fraction of the LW   !
       !      flux divergence in layer k is concentrated at just below top interface    !
       !      of layer k is invalid). Then, this condition might be sensitive to the    !
       !      vertical resolution of grid.                                              !
       !   3. LW radiative cooling (SW heating is assumed uniformly distributed through !
       !      layer k, so not relevant to buoyancy production) in the layer k. However, !
       !      SW production might also contribute, which may be considered in a future. !
       !   4. Internal stratification 'n2ht' of upper-half layer should be unstable.    !
       !      The 'n2ht' is pure internal stratification of upper half layer, obtained  !
       !      using internal slopes of sl, qt in layer k (in contrast to conventional   !
       !      interfacial slope) and saturation fraction in the upper-half layer,       !
       !      sfuh(k) (in contrast to sfi(k)).                                          !
       !   5. Top and bottom interfaces not both in the same existing convective layer. !
       !      If SRCL is within the previouisly identified CL regimes, we don't define  !
       !      a new SRCL.                                                               !
       !   6. k >= ntop_turb + 1 = 2                                                    !
       !   7. Ri at the top interface > ricrit = 0.19 (otherwise turbulent mixing will  !
       !      broadly distribute the cloud top in the vertical, preventing localized    !
       !      radiative destabilization at the top interface).                          !
       !                                                                                !
       ! Note if 'k = pver', it identifies a surface-based single fog layer, possibly,  !
       ! warm advection fog. Note also the CL regime index of SRCLs itself increases    !
       ! with height similar to the regular CLs indices identified from 'zisocl'.       !
       ! ------------------------------------------------------------------------------ !
       ncv  = 1
       ncvf = ncvfin(i)

       if( choice_SRCL .eq. 'remove' ) goto 222 

       do k = nbot_turb, ntop_turb + 1, -1 ! 'k = pver, 2, -1' is a layer index.

          if( ql(k,i) .gt. qmin .and. ql(k-1,i) .lt. qmin .and. qrlw(k,i) .lt. 0._r8 &
                                .and. ri(k,i) .ge. ricrit ) then

              ! In order to avoid any confliction with the treatment of ambiguous layer,
              ! I need to impose an additional constraint that ambiguous layer cannot be
              ! SRCL. So, I added constraint that 'k+1' interface (base interface of k
              ! layer) should not be a part of previously identified CL. Since 'belongcv'
              ! is even true for external entrainment interfaces, below constraint is
              ! fully sufficient.
 
              if( choice_SRCL .eq. 'nonamb' .and. belongcv(k+1,i) ) then
                  go to 220 
              endif

              ch = ( 1._r8 - sfuh(k,i) ) * chu(k,i) + sfuh(k,i) * chs(k,i)
              cm = ( 1._r8 - sfuh(k,i) ) * cmu(k,i) + sfuh(k,i) * cms(k,i)

              n2htSRCL = ch * slslope(k,i) + cm * qtslope(k,i)

              if( n2htSRCL .le. 0._r8 ) then

                  ! Test if bottom and top interfaces are part of the pre-existing CL. 
                  ! If not, find appropriate index for the new SRCL. Note that this
                  ! calculation makes use of 'ncv set' obtained from 'zisocl'. The 
                  ! 'in_CL' is a parameter testing whether the new SRCL is already 
                  ! within the pre-existing CLs (.true.) or not (.false.). 

                  in_CL = .false.

                  do while ( ncv .le. ncvf )
                     if( ktop(ncv,i) .le. k ) then
                        if( kbase(ncv,i) .gt. k ) then 
                            in_CL = .true.
                        endif
                        exit             ! Exit from 'do while' loop if SRCL is within the CLs.
                     else
                        ncv = ncv + 1    ! Go up one CL
                     end if
                  end do ! ncv

                  if( .not. in_CL ) then ! SRCL is not within the pre-existing CLs.

                     ! Identify a new SRCL and add it to the pre-existing CL regime group.

                     ncvfin(i)       =  ncvfin(i) + 1
                     ncvnew          =  ncvfin(i)
                     ktop(ncvnew,i)  =  k
                     kbase(ncvnew,i) =  k+1
                     belongcv(k,i)   = .true.
                     belongcv(k+1,i) = .true.

                     ! Calculate internal energy of SRCL. There is no internal energy if
                     ! SRCL is elevated from the surface. Also, we simply assume neutral 
                     ! stability function. Note that this assumption of neutral stability
                     ! does not influence numerical calculation- stability functions here
                     ! are just for diagnostic output. In general SRCLs other than a SRCL 
                     ! based at surface with bflxs <= 0, there is no other way but to use
                     ! neutral stability function.  However, in case of SRCL based at the
                     ! surface,  we can explicitly calculate non-zero stability functions            
                     ! in a consistent way.   Even though stability functions of SRCL are
                     ! just diagnostic outputs not influencing numerical calculations, it
                     ! would be informative to write out correct reasonable values rather
                     ! than simply assuming neutral stability. I am doing this right now.
                     ! Similar calculations were done for the SBCL and when surface inter
                     ! facial layer was merged by overlying CL in 'ziscol'.

                     if( k .lt. pver ) then

                         wbrk(ncvnew,i) = 0._r8
                         ebrk(ncvnew,i) = 0._r8
                         lbrk(ncvnew,i) = 0._r8
                         ghcl(ncvnew,i) = 0._r8
                         shcl(ncvnew,i) = 0._r8
                         smcl(ncvnew,i) = 0._r8
                         ricl(ncvnew,i) = 0._r8

                     else ! Surface-based fog

                         if( bflxs(i) .gt. 0._r8 ) then    ! Incorporate surface TKE into CL interior energy
                                                           ! It is likely that this case cannot exist  since
                                                           ! if surface buoyancy flux is positive,  it would
                                                           ! have been identified as SBCL in 'zisocl' ahead. 

                             ebrk(ncvnew,i) = tkes(i)
                             lbrk(ncvnew,i) = z(pver,i)
                             wbrk(ncvnew,i) = tkes(i) / b1    
        
                             print*, 'Major mistake in SRCL: bflxs > 0 for surface-based SRCL'
                             print*, 'bflxs = ', bflxs(i)
                             print*, 'ncvfin_o = ', ncvfin_o(i)
                             print*, 'ncvfin_mg = ', ncvfin_mg(i)
                             do ks = 1, ncvmax
                                print*, 'ncv =', ks, ' ', kbase_o(ks,i), ktop_o(ks,i), kbase_mg(ks,i), ktop_mg(ks,i)
                             end do
                             call endrun('CALEDDY: Major mistake in SRCL: bflxs > 0 for surface-based SRCL')

                         else                              ! Don't incorporate surface interfacial TKE into CL interior energy
                             ebrk(ncvnew,i) = 0._r8
                             lbrk(ncvnew,i) = 0._r8
                             wbrk(ncvnew,i) = 0._r8
                         endif

                         ! Calculate stability functions (ghcl, shcl, smcl, ricl) explicitly
                         ! using an reverse procedure starting from tkes(i). Note that it is
                         ! possible to calculate stability functions even when bflxs < 0.
                         ! Previous code just assumed neutral stability functions. Note that
                         ! since alph5 = 0.7 > 0, alph3 = -35 < 0, the denominator of gh  is
                         ! always positive if bflxs > 0. However, if bflxs < 0,  denominator
                         ! can be zero. For this case, we provide a possible maximum negative
                         ! value (the most stable state) to gh. Note also tkes(i) is always a
                         ! positive value by a limiter. Also, sprod(pver+1,i) > 0 by limiter.
 
                         gg = 0.5_r8 * vk * z(pver,i) * bprod(pver+1,i) / ( tkes(i)**(3._r8/2._r8) )
                         if( abs(alph5-gg*alph3) .le. 1.e-7_r8 ) then
                           ! gh = -0.28_r8
                           ! gh = -3.5334_r8
                             gh = ghmin
                         else    
                             gh = gg / ( alph5 - gg * alph3 )
                         end if 
                       ! gh = min(max(gh,-0.28_r8),0.0233_r8)
                       ! gh = min(max(gh,-3.5334_r8),0.0233_r8)
                         gh = min(max(gh,ghmin),0.0233_r8)
                         ghcl(ncvnew,i) =  gh
                         shcl(ncvnew,i) =  max(0._r8,alph5/(1._r8+alph3*gh))
                         smcl(ncvnew,i) =  max(0._r8,(alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))
                         ricl(ncvnew,i) = -(smcl(ncvnew,i)/shcl(ncvnew,i))*(bprod(pver+1,i)/sprod(pver+1,i))

                       ! 'ncvsurf' is CL regime index based at the surface. If there is no
                       ! such regime, then 'ncvsurf = 0'.
    
                         ncvsurf = ncvnew

                      end if

                  end if

              end if

          end if

   220 continue    

       end do ! End of 'k' loop where 'k' is a grid layer index running from 'pver' to 2

   222 continue

       ! -------------------------------------------------------------------------- !
       ! Up to this point, we identified all kinds of CL regimes :                  !
       !   1. A SBCL. By construction, 'bflxs > 0' for SBCL.                        !
       !   2. Surface-based CL with multiple layers and 'bflxs =< 0'                !
       !   3. Surface-based CL with multiple layers and 'bflxs > 0'                 !
       !   4. Regular elevated CL with two entraining interfaces                    ! 
       !   5. SRCLs. If SRCL is based at surface, it will be bflxs < 0.             !
       ! '1-4' were identified from 'zisocl' while '5' were identified separately   !
       ! after performing 'zisocl'. CL regime index of '1-4' increases with height  !
       ! ( e.g., CL = 1 is the CL regime nearest to the surface ) while CL regime   !
       ! index of SRCL is simply appended after the final index of CL regimes from  !
       ! 'zisocl'. However, CL regime indices of SRCLs itself increases with height !
       ! when there are multiple SRCLs, similar to the regular CLs from 'zisocl'.   !
       ! -------------------------------------------------------------------------- !

       ! Diagnostic output of final CL regimes indices
       
       do k = 1, ncvmax
          kbase_f(k,i) = real(kbase(k,i),r8)
          ktop_f(k,i)  = real(ktop(k,i),r8) 
          ncvfin_f(i)  = real(ncvfin(i),r8)
       end do 

       ! ---------------------------------------- !
       ! Perform do loop for individual CL regime !
       ! ---------------------------------------- ! -------------------------------- !
       ! For individual CLs, compute                                                 !
       !   1. Entrainment rates at the CL top and (if any) base interfaces using     !
       !      appropriate entrainment closure (current code use 'wstar' closure).    !
       !   2. Net CL mean (i.e., including entrainment contribution) TKE (ebrk)      !
       !      and normalized TKE (wbrk).                                             ! 
       !   3. TKE (tke) and normalized TKE (wcap) profiles at all CL interfaces.     !
       !   4. ( kvm, kvh ) profiles at all CL interfaces.                            !
       !   5. ( bprod, sprod ) profiles at all CL interfaces.                        !
       ! Also calculate                                                              !
       !   1. PBL height as the top external interface of surface-based CL, if any.  !
       !   2. Characteristic excesses of convective 'updraft velocity (wpert)',      !
       !      'temperature (tpert)', and 'moisture (qpert)' in the surface-based CL, !
       !      if any, for use in the separate convection scheme.                     ! 
       ! If there is no surface-based CL, 'PBL height' and 'convective excesses' are !
       ! calculated later from surface-based STL (Stable Turbulent Layer) properties.!
       ! --------------------------------------------------------------------------- !

       ktblw = 0
       do ncv = 1, ncvfin(i)

          kt = ktop(ncv,i)
          kb = kbase(ncv,i)
          ! Check whether surface interface is energetically interior or not.
          if( kb .eq. (pver+1) .and. bflxs(i) .le. 0._r8 ) then
              lbulk = zi(kt,i) - z(pver,i)
          else
              lbulk = zi(kt,i) - zi(kb,i)
          end if
          lbulk = min( lbulk, lbulk_max )

          ! Calculate 'turbulent length scale (leng)' and 'normalized TKE (wcap)'
          ! at all CL interfaces except the surface.  Note that below 'wcap' at 
          ! external interfaces are not correct. However, it does not influence 
          ! numerical calculation and correct normalized TKE at the entraining 
          ! interfaces will be re-calculated at the end of this 'do ncv' loop. 

          do k = min(kb,pver), kt, -1 
             if( choice_tunl .eq. 'rampcl' ) then
               ! In order to treat the case of 'ricl(ncv,i) >> 0' of surface-based SRCL
               ! with 'bflxs(i) < 0._r8', I changed ricl(ncv,i) -> min(0._r8,ricl(ncv,i))
               ! in the below exponential. This is necessary to prevent the model crash
               ! by too large values (e.g., 700) of ricl(ncv,i)   
                 tunlramp = ctunl*tunl*(1._r8-(1._r8-1._r8/ctunl)*exp(min(0._r8,ricl(ncv,i))))
                 tunlramp = min(max(tunlramp,tunl),ctunl*tunl)
             elseif( choice_tunl .eq. 'rampsl' ) then
                 tunlramp = ctunl*tunl
               ! tunlramp = 0.765_r8
             else
                 tunlramp = tunl
             endif
             if( choice_leng .eq. 'origin' ) then
                 leng(k,i) = ( (vk*zi(k,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
               ! leng(k,i) = vk*zi(k,i) / (1._r8+vk*zi(k,i)/(tunlramp*lbulk))
             else
                 leng(k,i) = min( vk*zi(k,i), tunlramp*lbulk )              
             endif
             leng(k,i) = min(leng_max(k), leng(k,i))
             wcap(k,i) = (leng(k,i)**2) * (-shcl(ncv,i)*n2(k,i)+smcl(ncv,i)*s2(k,i))
          end do ! k

          ! Calculate basic cross-interface variables ( jump condition ) across the 
          ! base external interface of CL.

          if( kb .lt. pver+1 ) then 

              jbzm = z(kb-1,i) - z(kb,i)                                      ! Interfacial layer thickness [m]
              jbsl = sl(kb-1,i) - sl(kb,i)                                    ! Interfacial jump of 'sl' [J/kg]
              jbqt = qt(kb-1,i) - qt(kb,i)                                    ! Interfacial jump of 'qt' [kg/kg]
              jbbu = n2(kb,i) * jbzm                                          ! Interfacial buoyancy jump [m/s2]
                                                                              ! considering saturation ( > 0 )
              jbbu = max(jbbu,jbumin)                                         ! Set minimum buoyancy jump, jbumin = 1.e-3
              jbu  = u(kb-1,i) - u(kb,i)                                      ! Interfacial jump of 'u' [m/s]
              jbv  = v(kb-1,i) - v(kb,i)                                      ! Interfacial jump of 'v' [m/s]
              ch   = (1._r8 -sflh(kb-1,i))*chu(kb,i) + sflh(kb-1,i)*chs(kb,i) ! Buoyancy coefficient just above the base interface
              cm   = (1._r8 -sflh(kb-1,i))*cmu(kb,i) + sflh(kb-1,i)*cms(kb,i) ! Buoyancy coefficient just above the base interface
              n2hb = (ch*jbsl + cm*jbqt)/jbzm                                 ! Buoyancy frequency [s-2]
                                                                              ! just above the base interface
              vyb  = n2hb*jbzm/jbbu                                           ! Ratio of 'n2hb/n2' at 'kb' interface
              vub  = min(1._r8,(jbu**2+jbv**2)/(jbbu*jbzm) )                  ! Ratio of 's2/n2 = 1/Ri' at 'kb' interface

          else 

            ! Below setting is necessary for consistent treatment when 'kb' is at the surface.
              jbbu = 0._r8
              n2hb = 0._r8
              vyb  = 0._r8
              vub  = 0._r8
              web  = 0._r8
          end if

          ! Calculate basic cross-interface variables ( jump condition ) across the 
          ! top external interface of CL. The meanings of variables are similar to
          ! the ones at the base interface.

          jtzm = z(kt-1,i) - z(kt,i)
          jtsl = sl(kt-1,i) - sl(kt,i)
          jtqt = qt(kt-1,i) - qt(kt,i)
          jtbu = n2(kt,i)*jtzm                                      ! Note : 'jtbu' is guaranteed positive by definition of CL top.
          jtbu = max(jtbu,jbumin)                                   ! But threshold it anyway to be sure.
          jtu  = u(kt-1,i) - u(kt,i)
          jtv  = v(kt-1,i) - v(kt,i)
          ch   = (1._r8 -sfuh(kt,i))*chu(kt,i) + sfuh(kt,i)*chs(kt,i) 
          cm   = (1._r8 -sfuh(kt,i))*cmu(kt,i) + sfuh(kt,i)*cms(kt,i) 
          n2ht = (ch*jtsl + cm*jtqt)/jtzm                       
          vyt  = n2ht*jtzm/jtbu                                  
          vut  = min(1._r8,(jtu**2+jtv**2)/(jtbu*jtzm))             

          ! Evaporative enhancement factor of entrainment rate at the CL top interface, evhc. 
          ! We take the full inversion strength to be 'jt2slv = slv(kt-2,i)-slv(kt,i)' 
          ! where 'kt-1' is in the ambiguous layer. However, for a cloud-topped CL overlain
          ! by another CL, it is possible that 'slv(kt-2,i) < slv(kt,i)'. To avoid negative
          ! or excessive evhc, we lower-bound jt2slv and upper-bound evhc.  Note 'jtslv' is
          ! used only for calculating 'evhc' : when calculating entrainment rate,   we will
          ! use normal interfacial buoyancy jump across CL top interface.

          evhc   = 1._r8
          jt2slv = 0._r8

        ! Modification : I should check whether below 'jbumin' produces reasonable limiting value.   
        !                In addition, our current formulation does not consider ice contribution. 

          if( choice_evhc .eq. 'orig' ) then

              if( ql(kt,i) .gt. qmin .and. ql(kt-1,i) .lt. qmin ) then 
                  jt2slv = slv(max(kt-2,1),i) - slv(kt,i)
                  jt2slv = max( jt2slv, jbumin*slv(kt-1,i)/g )
                  evhc   = 1._r8 + a2l * a3l * latvap * ql(kt,i) / jt2slv
                  evhc   = min( evhc, evhcmax )
              end if

          elseif( choice_evhc .eq. 'ramp' ) then

              jt2slv = slv(max(kt-2,1),i) - slv(kt,i)
              jt2slv = max( jt2slv, jbumin*slv(kt-1,i)/g )
              evhc   = 1._r8 + max(cldeff(kt,i)-cldeff(kt-1,i),0._r8) * a2l * a3l * latvap * ql(kt,i) / jt2slv
              evhc   = min( evhc, evhcmax )

          elseif( choice_evhc .eq. 'maxi' ) then

              qleff  = max( ql(kt-1,i), ql(kt,i) ) 
              jt2slv = slv(max(kt-2,1),i) - slv(kt,i)
              jt2slv = max( jt2slv, jbumin*slv(kt-1,i)/g )
              evhc   = 1._r8 + a2l * a3l * latvap * qleff / jt2slv
              evhc   = min( evhc, evhcmax )

          endif

          ! Calculate cloud-top radiative cooling contribution to buoyancy production.
          ! Here,  'radf' [m2/s3] is additional buoyancy flux at the CL top interface 
          ! associated with cloud-top LW cooling being mainly concentrated near the CL
          ! top interface ( just below CL top interface ).  Contribution of SW heating
          ! within the cloud is not included in this radiative buoyancy production 
          ! since SW heating is more broadly distributed throughout the CL top layer. 

          lwp        = 0._r8
          opt_depth  = 0._r8
          radinvfrac = 0._r8 
          radf       = 0._r8

          if( choice_radf .eq. 'orig' ) then

              if( ql(kt,i) .gt. qmin .and. ql(kt-1,i) .lt. qmin ) then 

                  lwp       = ql(kt,i) * ( pi(kt+1,i) - pi(kt,i) ) / g
                  opt_depth = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer

                  ! Approximate LW cooling fraction concentrated at the inversion by using
                  ! polynomial approx to exact formula 1-2/opt_depth+2/(exp(opt_depth)-1))

                  radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
                  radf        = qrlw(kt,i) / ( pi(kt,i) - pi(kt+1,i) ) ! Cp*radiative cooling = [ W/kg ] 
                  radf        = max( radinvfrac * radf * ( zi(kt,i) - zi(kt+1,i) ), 0._r8 ) * chs(kt,i)
                ! We can disable cloud LW cooling contribution to turbulence by uncommenting:
                ! radf = 0._r8

              end if

          elseif( choice_radf .eq. 'ramp' ) then

                  lwp         = ql(kt,i) * ( pi(kt+1,i) - pi(kt,i) ) / g
                  opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
                  radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
                  radinvfrac  = max(cldeff(kt,i)-cldeff(kt-1,i),0._r8) * radinvfrac 
                  radf        = qrlw(kt,i) / ( pi(kt,i) - pi(kt+1,i) ) ! Cp*radiative cooling [W/kg] 
                  radf        = max( radinvfrac * radf * ( zi(kt,i) - zi(kt+1,i) ), 0._r8 ) * chs(kt,i)

          elseif( choice_radf .eq. 'maxi' ) then

                ! Radiative flux divergence both in 'kt' and 'kt-1' layers are included 
                ! 1. From 'kt' layer
                  lwp         = ql(kt,i) * ( pi(kt+1,i) - pi(kt,i) ) / g
                  opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
                  radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth ) + opt_depth**2 )
                  radf        = max( radinvfrac * qrlw(kt,i) / ( pi(kt,i) - pi(kt+1,i) ) * ( zi(kt,i) - zi(kt+1,i) ), 0._r8 )
                ! 2. From 'kt-1' layer and add the contribution from 'kt' layer
                  lwp         = ql(kt-1,i) * ( pi(kt,i) - pi(kt-1,i) ) / g
                  opt_depth   = 156._r8 * lwp  ! Estimated LW optical depth in the CL top layer
                  radinvfrac  = opt_depth * ( 4._r8 + opt_depth ) / ( 6._r8 * ( 4._r8 + opt_depth) + opt_depth**2 )
                  radf        = radf + &
                       max( radinvfrac * qrlw(kt-1,i) / ( pi(kt-1,i) - pi(kt,i) ) * ( zi(kt-1,i) - zi(kt,i) ), 0._r8 )
                  radf        = max( radf, 0._r8 ) * chs(kt,i) 

          endif

          ! ------------------------------------------------------------------- !
          ! Calculate 'wstar3' by summing buoyancy productions within CL from   !
          !   1. Interior buoyancy production ( bprod: fcn of TKE )             !
          !   2. Cloud-top radiative cooling                                    !
          !   3. Surface buoyancy flux contribution only when bflxs > 0.        !
          !      Note that master length scale, lbulk, has already been         !
          !      corrctly defined at the first part of this 'do ncv' loop       !
          !      considering the sign of bflxs.                                 !
          ! This 'wstar3' is used for calculation of entrainment rate.          !
          ! Note that this 'wstar3' formula does not include shear production   !
          ! and the effect of drizzle, which should be included later.          !
          ! Q : Strictly speaking, in calculating interior buoyancy production, ! 
          !     the use of 'bprod' is not correct, since 'bprod' is not correct !
          !     value but initially guessed value.   More reasonably, we should ! 
          !     use '-leng(k,i)*sqrt(b1*wcap(k,i))*shcl(ncv,i)*n2(k,i)' instead !
          !     of 'bprod(k,i)', although this is still an  approximation since !
          !     tke(k,i) is not exactly 'b1*wcap(k,i)'  due to a transport term.! 
          !     However since iterative calculation will be performed after all,! 
          !     below might also be OK. But I should test this alternative.     !
          ! ------------------------------------------------------------------- !      

          dzht   = zi(kt,i)  - z(kt,i)     ! Thickness of CL top half-layer
          dzhb   = z(kb-1,i) - zi(kb,i)    ! Thickness of CL bot half-layer
          wstar3 = radf * dzht
          do k = kt + 1, kb - 1 ! If 'kt = kb - 1', this loop will not be performed. 
               wstar3 =  wstar3 + bprod(k,i) * ( z(k-1,i) - z(k,i) )
             ! Below is an alternative which may speed up convergence.
             ! However, for interfaces merged into original CL, it can
             ! be 'wcap(k,i)<0' since 'n2(k,i)>0'.  Thus, I should use
             ! the above original one.
             ! wstar3 =  wstar3 - leng(k,i)*sqrt(b1*wcap(k,i))*shcl(ncv,i)*n2(k,i)* &
             !                    (z(k-1,i) - z(k,i))
          end do      
          if( kb .eq. (pver+1) .and. bflxs(i) .gt. 0._r8 ) then
             wstar3 = wstar3 + bflxs(i) * dzhb
           ! wstar3 = wstar3 + bprod(pver+1,i) * dzhb
          end if   
          wstar3 = max( 2.5_r8 * wstar3, 0._r8 )
   
          ! -------------------------------------------------------------- !
          ! Below single block is for 'sedimentation-entrainment feedback' !
          ! -------------------------------------------------------------- !          

          if( id_sedfact ) then
            ! wsed    = 7.8e5_r8*(ql(kt,i)/ncliq(kt,i))**(2._r8/3._r8)
              sedfact = exp(-ased*wsedl(kt,i)/(wstar3**(1._r8/3._r8)+1.e-6_r8))
              if( choice_evhc .eq. 'orig' ) then
                  if (ql(kt,i).gt.qmin .and. ql(kt-1,i).lt.qmin) then
                      jt2slv = slv(max(kt-2,1),i) - slv(kt,i)
                      jt2slv = max(jt2slv, jbumin*slv(kt-1,i)/g)
                      evhc = 1._r8+sedfact*a2l*a3l*latvap*ql(kt,i) / jt2slv
                      evhc = min(evhc,evhcmax)
                  end if
              elseif( choice_evhc .eq. 'ramp' ) then
                  jt2slv = slv(max(kt-2,1),i) - slv(kt,i)
                  jt2slv = max(jt2slv, jbumin*slv(kt-1,i)/g)
                  evhc = 1._r8+max(cldeff(kt,i)-cldeff(kt-1,i),0._r8)*sedfact*a2l*a3l*latvap*ql(kt,i) / jt2slv
                  evhc = min(evhc,evhcmax)
              elseif( choice_evhc .eq. 'maxi' ) then
                  qleff  = max(ql(kt-1,i),ql(kt,i))
                  jt2slv = slv(max(kt-2,1),i) - slv(kt,i)
                  jt2slv = max(jt2slv, jbumin*slv(kt-1,i)/g)
                  evhc = 1._r8+sedfact*a2l*a3l*latvap*qleff / jt2slv
                  evhc = min(evhc,evhcmax)
              endif
          endif

          ! -------------------------------------------------------------------------- !
          ! Now diagnose CL top and bottom entrainment rates (and the contribution of  !
          ! top/bottom entrainments to wstar3) using entrainment closures of the form  !
          !                                                                            !        
          !                   wet = cet*wstar3, web = ceb*wstar3                       !
          !                                                                            !
          ! where cet and ceb depend on the entrainment interface jumps, ql, etc.      !
          ! No entrainment is diagnosed unless the wstar3 > 0. Note '1/wstar3fact' is  !
          ! a factor indicating the enhancement of wstar3 due to entrainment process.  !
          ! Q : Below setting of 'wstar3fact = max(..,0.5)'might prevent the possible  !
          !     case when buoyancy consumption by entrainment is  stronger than cloud  !
          !     top radiative cooling production. Is that OK ? No.  According to bulk  !
          !     modeling study, entrainment buoyancy consumption was always a certain  !
          !     fraction of other net productions, rather than a separate sum.  Thus,  !
          !     below max limit of wstar3fact is correct.   'wstar3fact = max(.,0.5)'  !
          !     prevents unreasonable enhancement of CL entrainment rate by cloud-top  !
          !     entrainment instability, CTEI.                                         !
          ! Q : Use of the same dry entrainment coefficient, 'a1i' both at the CL  top !
          !     and base interfaces may result in too small 'wstar3' and 'ebrk' below, !
          !     as was seen in my generalized bulk modeling study. This should be re-  !
          !     considered later                                                       !
          ! -------------------------------------------------------------------------- !
          
          if( wstar3 .gt. 0._r8 ) then
              cet = a1i * evhc / ( jtbu * lbulk )
              if( kb .eq. pver + 1 ) then 
                  wstar3fact = max( 1._r8 + 2.5_r8 * cet * n2ht * jtzm * dzht, wstar3factcrit )
              else    
                  ceb = a1i / ( jbbu * lbulk )
                  wstar3fact = max( 1._r8 + 2.5_r8 * cet * n2ht * jtzm * dzht &
                                          + 2.5_r8 * ceb * n2hb * jbzm * dzhb, wstar3factcrit )
              end if
              wstar3 = wstar3 / wstar3fact       
          else ! wstar3 == 0
              wstar3fact = 0._r8 ! This is just for dianostic output
              cet        = 0._r8
              ceb        = 0._r8
          end if 

          ! ---------------------------------------------------------------------------- !
          ! Calculate net CL mean TKE including entrainment contribution by solving a    !
          ! canonical cubic equation. The solution of cubic equ. is 'rootp**2 = ebrk'    !
          ! where 'ebrk' originally (before solving cubic eq.) was interior CL mean TKE, !
          ! but after solving cubic equation,  it is replaced by net CL mean TKE in the  !
          ! same variable 'ebrk'.                                                        !
          ! ---------------------------------------------------------------------------- !
          ! Solve cubic equation (canonical form for analytic solution)                  !
          !   r^3 - 3*trmp*r - 2*trmq = 0,   r = sqrt<e>                                 ! 
          ! to estimate <e> for CL, derived from layer-mean TKE balance:                 !
          !                                                                              !
          !   <e>^(3/2)/(b_1*<l>) \approx <B + S>   (*)                                  !
          !   <B+S> = (<B+S>_int * l_int + <B+S>_et * dzt + <B+S>_eb * dzb)/lbulk        !
          !   <B+S>_int = <e>^(1/2)/(b_1*<l>)*<e>_int                                    !
          !   <B+S>_et  = (-vyt+vut)*wet*jtbu + radf                                     !
          !   <B+S>_eb  = (-vyb+vub)*web*jbbu                                            !
          !                                                                              !
          ! where:                                                                       !
          !   <> denotes a vertical avg (over the whole CL unless indicated)             !
          !   l_int (called lbrk below) is aggregate thickness of interior CL layers     !
          !   dzt = zi(kt,i)-z(kt,i)   is thickness of top entrainment layer             !
          !   dzb = z(kb-1,i)-zi(kb,i) is thickness of bot entrainment layer             !
          !   <e>_int (called ebrk below) is the CL-mean TKE if only interior            !
          !                               interfaces contributed.                        !
          !   wet, web                  are top. bottom entrainment rates                !
          !                                                                              !
          ! For a single-level radiatively-driven convective layer, there are no         ! 
          ! interior interfaces so 'ebrk' = 'lbrk' = 0. If the CL goes to the            !
          ! surface, 'vyb' and 'vub' are set to zero before and 'ebrk' and 'lbrk'        !
          ! have already incorporated the surface interfacial layer contribution,        !
          ! so the same formulas still apply.                                            !
          !                                                                              !
          ! In the original formulation based on TKE,                                    !
          !    wet*jtbu = a1l*evhc*<e>^3/2/leng(kt,i)                                    ! 
          !    web*jbbu = a1l*<e>^3/2/leng(kt,i)                                         !
          !                                                                              !
          ! In the wstar formulation                                                     !
          !    wet*jtbu = a1i*evhc*wstar3/lbulk                                          !
          !    web*jbbu = a1i*wstar3/lbulk,                                              !
          ! ---------------------------------------------------------------------------- !

          fact = ( evhc * ( -vyt + vut ) * dzht + ( -vyb + vub ) * dzhb * leng(kb,i) / leng(kt,i) ) / lbulk

          if( wstarent ) then

              ! (Option 1) 'wstar' entrainment formulation 
              ! Here trmq can have either sign, and will usually be nonzero even for non-
              ! cloud topped CLs.  If trmq > 0, there will be two positive roots r; we take 
              ! the larger one. Why ? If necessary, we limit entrainment and wstar to prevent
              ! a solution with r < ccrit*wstar ( Why ? ) where we take ccrit = 0.5. 

              trma = 1._r8          
              trmp = ebrk(ncv,i) * ( lbrk(ncv,i) / lbulk ) / 3._r8 + ntzero
              trmq = 0.5_r8 * b1 * ( leng(kt,i)  / lbulk ) * ( radf * dzht + a1i * fact * wstar3 )

              ! Check if there is an acceptable root with r > rcrit = ccrit*wstar. 
              ! To do this, first find local minimum fmin of the cubic f(r) at sqrt(p), 
              ! and value fcrit = f(rcrit).

              rmin  = sqrt(trmp)
              fmin  = rmin * ( rmin * rmin - 3._r8 * trmp ) - 2._r8 * trmq
              wstar = wstar3**onet
              rcrit = ccrit * wstar
              fcrit = rcrit * ( rcrit * rcrit - 3._r8 * trmp ) - 2._r8 * trmq

              ! No acceptable root exists (noroot = .true.) if either:
              !    1) rmin < rcrit (in which case cubic is monotone increasing for r > rcrit)
              !       and f(rcrit) > 0.
              ! or 2) rmin > rcrit (in which case min of f(r) in r > rcrit is at rmin)
              !       and f(rmin) > 0.  
              ! In this case, we reduce entrainment and wstar3 such that r/wstar = ccrit;
              ! this changes the coefficients of the cubic.   It might be informative to
              ! check when and how many 'noroot' cases occur,  since when 'noroot',   we
              ! will impose arbitrary limit on 'wstar3, wet, web, and ebrk' using ccrit.

              noroot = ( ( rmin .lt. rcrit ) .and. ( fcrit .gt. 0._r8 ) ) &
                  .or. ( ( rmin .ge. rcrit ) .and. ( fmin  .gt. 0._r8 ) )
              if( noroot ) then ! Solve cubic for r
                  trma = 1._r8 - b1 * ( leng(kt,i) / lbulk ) * a1i * fact / ccrit**3
                  trma = max( trma, 0.5_r8 )  ! Limit entrainment enhancement of ebrk
                  trmp = trmp / trma 
                  trmq = 0.5_r8 * b1 * ( leng(kt,i) / lbulk ) * radf * dzht / trma
              end if   ! noroot

              ! Solve the cubic equation

              qq = trmq**2 - trmp**3
              if( qq .ge. 0._r8 ) then 
                  rootp = ( trmq + sqrt(qq) )**(1._r8/3._r8) + ( max( trmq - sqrt(qq), 0._r8 ) )**(1._r8/3._r8)
              else
                  rootp = 2._r8 * sqrt(trmp) * cos( acos( trmq / sqrt(trmp**3) ) / 3._r8 )
              end if
 
              ! Adjust 'wstar3' only if there is 'noroot'. 
              ! And calculate entrainment rates at the top and base interfaces.

              if( noroot )  wstar3 = ( rootp / ccrit )**3     ! Adjust wstar3 
              wet = cet * wstar3                              ! Find entrainment rates
              if( kb .lt. pver + 1 ) web = ceb * wstar3       ! When 'kb.eq.pver+1', it was set to web=0. 

          else !

              ! (Option.2) wstarentr = .false. Use original entrainment formulation.
              ! trmp > 0 if there are interior interfaces in CL, trmp = 0 otherwise.
              ! trmq > 0 if there is cloudtop radiative cooling, trmq = 0 otherwise.
             
              trma = 1._r8 - b1 * a1l * fact
              trma = max( trma, 0.5_r8 )  ! Prevents runaway entrainment instability
              trmp = ebrk(ncv,i) * ( lbrk(ncv,i) / lbulk ) / ( 3._r8 * trma )
              trmq = 0.5_r8 * b1 * ( leng(kt,i)  / lbulk ) * radf * dzht / trma

              qq = trmq**2 - trmp**3
              if( qq .ge. 0._r8 ) then 
                  rootp = ( trmq + sqrt(qq) )**(1._r8/3._r8) + ( max( trmq - sqrt(qq), 0._r8 ) )**(1._r8/3._r8)
              else ! Also part of case 3
                  rootp = 2._r8 * sqrt(trmp) * cos( acos( trmq / sqrt(trmp**3) ) / 3._r8 )
              end if   ! qq

             ! Find entrainment rates and limit them by free-entrainment values a1l*sqrt(e)

              wet = a1l * rootp * min( evhc * rootp**2 / ( leng(kt,i) * jtbu ), 1._r8 )   
              if( kb .lt. pver + 1 ) web = a1l * rootp * min( evhc * rootp**2 / ( leng(kb,i) * jbbu ), 1._r8 )

          end if ! wstarentr

          ! ---------------------------------------------------- !
          ! Finally, get the net CL mean TKE and normalized TKE  ! 
          ! ---------------------------------------------------- !

          ebrk(ncv,i) = rootp**2
          ebrk(ncv,i) = min(ebrk(ncv,i),tkemax) ! Limit CL-avg TKE used for entrainment
          wbrk(ncv,i) = ebrk(ncv,i)/b1  
        
          ! The only way ebrk = 0 is for SRCL which are actually radiatively cooled 
          ! at top interface. In this case, we remove 'convective' label from the 
          ! interfaces around this layer. This case should now be impossible, so 
          ! we flag it. Q: I can't understand why this case is impossible now. Maybe,
          ! due to various limiting procedures used in solving cubic equation ? 
          ! In case of SRCL, 'ebrk' should be positive due to cloud top LW radiative
          ! cooling contribution, although 'ebrk(internal)' of SRCL before including
          ! entrainment contribution (which include LW cooling contribution also) is
          ! zero. 

          if( ebrk(ncv,i) .le. 0._r8 ) then
              print*, 'CALEDDY: Warning, CL with zero TKE, i, kt, kb ', i, kt, kb
              belongcv(kt,i) = .false.
              belongcv(kb,i) = .false. 
          end if
 
          ! ----------------------------------------------------------------------- !
          ! Calculate complete TKE profiles at all CL interfaces, capped by tkemax. !
          ! We approximate TKE = <e> at entrainment interfaces. However when CL is  !
          ! based at surface, correct 'tkes' will be inserted to tke(pver+1,i).     !
          ! Note that this approximation at CL external interfaces do not influence !
          ! numerical calculation since 'e' at external interfaces are not used  in !
          ! actual numerical calculation afterward. In addition in order to extract !
          ! correct TKE averaged over the PBL in the cumulus scheme,it is necessary !
          ! to set e = <e> at the top entrainment interface.  Since net CL mean TKE !
          ! 'ebrk' obtained by solving cubic equation already includes tkes  ( tkes !
          ! is included when bflxs > 0 but not when bflxs <= 0 into internal ebrk ),!
          ! 'tkes' should be written to tke(pver+1,i)                               !
          ! ----------------------------------------------------------------------- !

          ! 1. At internal interfaces          
          do k = kb - 1, kt + 1, -1
             rcap = ( b1 * ae + wcap(k,i) / wbrk(ncv,i) ) / ( b1 * ae + 1._r8 )
             rcap = min( max(rcap,rcapmin), rcapmax )
             tke(k,i) = ebrk(ncv,i) * rcap
             tke(k,i) = min( tke(k,i), tkemax )
             kvh(k,i) = leng(k,i) * sqrt(tke(k,i)) * shcl(ncv,i)
             kvm(k,i) = leng(k,i) * sqrt(tke(k,i)) * smcl(ncv,i)
             bprod(k,i) = -kvh(k,i) * n2(k,i)
             sprod(k,i) =  kvm(k,i) * s2(k,i)
             turbtype(k,i) = 2                     ! CL interior interfaces.
             sm_aw(k,i) = smcl(ncv,i)/alph1        ! Diagnostic output for microphysics
          end do

          ! 2. At CL top entrainment interface
          kentr = wet * jtzm
          kvh(kt,i) = kentr
          kvm(kt,i) = kentr
          bprod(kt,i) = -kentr * n2ht + radf       ! I must use 'n2ht' not 'n2'
          sprod(kt,i) =  kentr * s2(kt,i)
          turbtype(kt,i) = 4                       ! CL top entrainment interface
          trmp = -b1 * ae / ( 1._r8 + b1 * ae )
          trmq = -(bprod(kt,i)+sprod(kt,i))*b1*leng(kt,i)/(1._r8+b1*ae)/(ebrk(ncv,i)**(3._r8/2._r8))
          rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
          rcap = min( max(rcap,rcapmin), rcapmax )
          tke(kt,i)  = ebrk(ncv,i) * rcap
          tke(kt,i)  = min( tke(kt,i), tkemax )
          sm_aw(kt,i) = smcl(ncv,i) / alph1        ! Diagnostic output for microphysics

          ! 3. At CL base entrainment interface and double entraining interfaces
          ! When current CL base is also the top interface of CL regime below,
          ! simply add the two contributions for calculating eddy diffusivity
          ! and buoyancy/shear production. Below code correctly works because
          ! we (CL regime index) always go from surface upward.

          if( kb .lt. pver + 1 ) then 

              kentr = web * jbzm

              if( kb .ne. ktblw ) then

                  kvh(kb,i) = kentr
                  kvm(kb,i) = kentr
                  bprod(kb,i) = -kvh(kb,i)*n2hb     ! I must use 'n2hb' not 'n2'
                  sprod(kb,i) =  kvm(kb,i)*s2(kb,i)
                  turbtype(kb,i) = 3                ! CL base entrainment interface
                  trmp = -b1*ae/(1._r8+b1*ae)
                  trmq = -(bprod(kb,i)+sprod(kb,i))*b1*leng(kb,i)/(1._r8+b1*ae)/(ebrk(ncv,i)**(3._r8/2._r8))
                  rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
                  rcap = min( max(rcap,rcapmin), rcapmax )
                  tke(kb,i)  = ebrk(ncv,i) * rcap
                  tke(kb,i)  = min( tke(kb,i),tkemax )

              else
                  
                  kvh(kb,i) = kvh(kb,i) + kentr 
                  kvm(kb,i) = kvm(kb,i) + kentr
                ! dzhb5 : Half thickness of the lowest  layer of  current CL regime
                ! dzht5 : Half thickness of the highest layer of adjacent CL regime just below current CL. 
                  dzhb5 = z(kb-1,i) - zi(kb,i)
                  dzht5 = zi(kb,i) - z(kb,i)
                  bprod(kb,i) = ( dzht5*bprod(kb,i) - dzhb5*kentr*n2hb )     / ( dzhb5 + dzht5 )
                  sprod(kb,i) = ( dzht5*sprod(kb,i) + dzhb5*kentr*s2(kb,i) ) / ( dzhb5 + dzht5 )
                  trmp = -b1*ae/(1._r8+b1*ae)
                  trmq = -kentr*(s2(kb,i)-n2hb)*b1*leng(kb,i)/(1._r8+b1*ae)/(ebrk(ncv,i)**(3._r8/2._r8))
                  rcap = compute_cubic(0._r8,trmp,trmq)**2._r8
                  rcap = min( max(rcap,rcapmin), rcapmax )
                  tke_imsi = ebrk(ncv,i) * rcap
                  tke_imsi = min( tke_imsi, tkemax )
                  tke(kb,i)  = ( dzht5*tke(kb,i) + dzhb5*tke_imsi ) / ( dzhb5 + dzht5 )               
                  tke(kb,i)  = min(tke(kb,i),tkemax)
                  turbtype(kb,i) = 5                ! CL double entraining interface      
                 
              end if

           else

             ! If CL base interface is surface, compute similarly using wcap(kb,i)=tkes/b1    
             ! Even when bflx < 0, use the same formula in order to impose consistency of
             ! tke(kb,i) at bflx = 0._r8
 
             rcap = (b1*ae + wcap(kb,i)/wbrk(ncv,i))/(b1*ae + 1._r8)
             rcap = min( max(rcap,rcapmin), rcapmax )
             tke(kb,i) = ebrk(ncv,i) * rcap
             tke(kb,i) = min( tke(kb,i),tkemax )

          end if

          ! For double entraining interface, simply use smcl(ncv,i) of the overlying CL. 
          ! Below 'sm_aw' is a diagnostic output for use in the microphysics.
          ! When 'kb' is surface, 'sm' will be over-written later below.

          sm_aw(kb,i) = smcl(ncv,i)/alph1             

          ! Calculate wcap at all interfaces of CL. Put a  minimum threshold on TKE
          ! to prevent possible division by zero.  'wcap' at CL internal interfaces
          ! are already calculated in the first part of 'do ncv' loop correctly.
          ! When 'kb.eq.pver+1', below formula produces the identical result to the
          ! 'tkes(i)/b1' if leng(kb,i) is set to vk*z(pver,i). Note  wcap(pver+1,i)
          ! is already defined as 'tkes(i)/b1' at the first part of caleddy.
          
          wcap(kt,i) = (bprod(kt,i)+sprod(kt,i))*leng(kt,i)/sqrt(max(tke(kt,i),1.e-6_r8))
          if( kb .lt. pver + 1 ) then
              wcap(kb,i) = (bprod(kb,i)+sprod(kb,i))*leng(kb,i)/sqrt(max(tke(kb,i),1.e-6_r8))
          end if

          ! Save the index of upper external interface of current CL-regime in order to
          ! handle the case when this interface is also the lower external interface of 
          ! CL-regime located just above. 

          ktblw = kt 

          ! Diagnostic Output

          wet_CL(ncv,i)        = wet
          web_CL(ncv,i)        = web
          jtbu_CL(ncv,i)       = jtbu
          jbbu_CL(ncv,i)       = jbbu
          evhc_CL(ncv,i)       = evhc
          jt2slv_CL(ncv,i)     = jt2slv
          n2ht_CL(ncv,i)       = n2ht
          n2hb_CL(ncv,i)       = n2hb          
          lwp_CL(ncv,i)        = lwp
          opt_depth_CL(ncv,i)  = opt_depth
          radinvfrac_CL(ncv,i) = radinvfrac
          radf_CL(ncv,i)       = radf
          wstar_CL(ncv,i)      = wstar          
          wstar3fact_CL(ncv,i) = wstar3fact          

       end do        ! ncv

       ! Calculate PBL height and characteristic cumulus excess for use in the
       ! cumulus convection shceme. Also define turbulence type at the surface
       ! when the lowest CL is based at the surface. These are just diagnostic
       ! outputs, not influencing numerical calculation of current PBL scheme.
       ! If the lowest CL is based at the surface, define the PBL depth as the
       ! CL top interface. The same rule is applied for all CLs including SRCL.

       if( ncvsurf .gt. 0 ) then

           ktopbl(i) = ktop(ncvsurf, i)
           pblh(i)   = zi(ktopbl(i), i)
           pblhp(i)  = pi(ktopbl(i), i)
           wpert(i)  = max(wfac*sqrt(ebrk(ncvsurf,i)),wpertmin)
           tpert(i)  = max(abs(shflx(i)*rrho(i)/cpair)*tfac/wpert(i),0._r8)
           qpert(i)  = max(abs(qflx(i)*rrho(i))*tfac/wpert(i),0._r8)

           if( bflxs(i) .gt. 0._r8 ) then
               turbtype(pver+1,i) = 2 ! CL interior interface
           else
               turbtype(pver+1,i) = 3 ! CL external base interface
           endif

           ipbl(i)  = 1._r8
           kpblh(i) = ktopbl(i) - 1._r8

       end if ! End of the calculationf of te properties of surface-based CL.

       ! -------------------------------------------- !
       ! Treatment of Stable Turbulent Regime ( STL ) !
       ! -------------------------------------------- !

       ! Identify top and bottom most (internal) interfaces of STL except surface.
       ! Also, calculate 'turbulent length scale (leng)' at each STL interfaces.     

       belongst(1,i) = .false.   ! k = 1 (top interface) is assumed non-turbulent
       do k = 2, pver            ! k is an interface index
          belongst(k,i) = ( ri(k,i) .lt. ricrit ) .and. ( .not. belongcv(k,i) )
          if( belongst(k,i) .and. ( .not. belongst(k-1,i) ) ) then
              kt = k             ! Top interface index of STL
          elseif( .not. belongst(k,i) .and. belongst(k-1,i) ) then
              kb = k - 1         ! Base interface index of STL
              lbulk = z(kt-1,i) - z(kb,i)
              lbulk = min( lbulk, lbulk_max )
              do ks = kt, kb
                 if( choice_tunl .eq. 'rampcl' ) then
                     tunlramp = tunl
                 elseif( choice_tunl .eq. 'rampsl' ) then
                    tunlramp = max( 1.e-3_r8, ctunl * tunl * exp(-log(ctunl)*ri(ks,i)/ricrit) )
                  ! tunlramp = 0.065_r8 + 0.7_r8 * exp(-20._r8*ri(ks,i))
                 else
                    tunlramp = tunl
                 endif
                 if( choice_leng .eq. 'origin' ) then
                     leng(ks,i) = ( (vk*zi(ks,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                   ! leng(ks,i) = vk*zi(ks,i) / (1._r8+vk*zi(ks,i)/(tunlramp*lbulk))
                 else
                     leng(ks,i) = min( vk*zi(ks,i), tunlramp*lbulk )              
                 endif
                 leng(ks,i) = min(leng_max(ks), leng(ks,i))
              end do
          end if
       end do ! k

       ! Now look whether STL extends to ground.  If STL extends to surface,
       ! re-define master length scale,'lbulk' including surface interfacial
       ! layer thickness, and re-calculate turbulent length scale, 'leng' at
       ! all STL interfaces again. Note that surface interface is assumed to
       ! always be STL if it is not CL.   
       
       belongst(pver+1,i) = .not. belongcv(pver+1,i)

       if( belongst(pver+1,i) ) then     ! kb = pver+1 (surface  STL)

           turbtype(pver+1,i) = 1        ! Surface is STL interface
          
           if( belongst(pver,i) ) then   ! STL includes interior
             ! 'kt' already defined above as the top interface of STL
               lbulk = z(kt-1,i)          
           else                          ! STL with no interior turbulence
               kt = pver+1
               lbulk = z(kt-1,i)
           end if
           lbulk = min( lbulk, lbulk_max )

           ! PBL height : Layer mid-point just above the highest STL interface
           ! Note in contrast to the surface based CL regime where  PBL height
           ! was defined at the top external interface, PBL height of  surface
           ! based STL is defined as the layer mid-point.

           ktopbl(i) = kt - 1
           pblh(i)   = z(ktopbl(i),i)
           pblhp(i)  = 0.5_r8 * ( pi(ktopbl(i),i) + pi(ktopbl(i)+1,i) )          

           ! Re-calculate turbulent length scale including surface interfacial
           ! layer contribution to lbulk.

           do ks = kt, pver
              if( choice_tunl .eq. 'rampcl' ) then
                  tunlramp = tunl
              elseif( choice_tunl .eq. 'rampsl' ) then
                  tunlramp = max(1.e-3_r8,ctunl*tunl*exp(-log(ctunl)*ri(ks,i)/ricrit))
                ! tunlramp = 0.065_r8 + 0.7_r8 * exp(-20._r8*ri(ks,i))
              else
                  tunlramp = tunl
              endif
              if( choice_leng .eq. 'origin' ) then
                  leng(ks,i) = ( (vk*zi(ks,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! leng(ks,i) = vk*zi(ks,i) / (1._r8+vk*zi(ks,i)/(tunlramp*lbulk))
              else
                  leng(ks,i) = min( vk*zi(ks,i), tunlramp*lbulk )              
              endif
             leng(ks,i) = min(leng_max(ks), leng(ks,i))
           end do ! ks

           ! Characteristic cumulus excess of surface-based STL.
           ! We may be able to use ustar for wpert.

           wpert(i) = 0._r8 
           tpert(i) = max(shflx(i)*rrho(i)/cpair*fak/ustar(i),0._r8) ! CCM stable-layer forms
           qpert(i) = max(qflx(i)*rrho(i)*fak/ustar(i),0._r8)

           ipbl(i)  = 0._r8
           kpblh(i) = ktopbl(i)

       end if

       ! Calculate stability functions and energetics at the STL interfaces
       ! except the surface. Note that tke(pver+1,i) and wcap(pver+1,i) are
       ! already calculated in the first part of 'caleddy', kvm(pver+1,i) &
       ! kvh(pver+1,i) were already initialized to be zero, bprod(pver+1,i)
       ! & sprod(pver+1,i) were direcly calculated from the bflxs and ustar.
       ! Note transport term is assumed to be negligible at STL interfaces.
           
       do k = 2, pver

          if( belongst(k,i) ) then

              turbtype(k,i) = 1    ! STL interfaces
              trma = alph3*alph4exs*ri(k,i) + 2._r8*b1*(alph2-alph4exs*alph5*ri(k,i))
              trmb = (alph3+alph4exs)*ri(k,i) + 2._r8*b1*(-alph5*ri(k,i)+alph1)
              trmc = ri(k,i)
              det = max(trmb*trmb-4._r8*trma*trmc,0._r8)
              ! Sanity Check
              if( det .lt. 0._r8 ) then
                  call endrun('CALEDDY: The det < 0. for the STL in UW eddy_diff')
              end if                  
              gh = (-trmb + sqrt(det))/(2._r8*trma)
            ! gh = min(max(gh,-0.28_r8),0.0233_r8)
            ! gh = min(max(gh,-3.5334_r8),0.0233_r8)
              gh = min(max(gh,ghmin),0.0233_r8)
              sh = max(0._r8,alph5/(1._r8+alph3*gh))
              sm = max(0._r8,(alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))

              tke(k,i)   = b1*(leng(k,i)**2)*(-sh*n2(k,i)+sm*s2(k,i))
              tke(k,i)   = min(tke(k,i),tkemax)
              wcap(k,i)  = tke(k,i)/b1
              kvh(k,i)   = leng(k,i) * sqrt(tke(k,i)) * sh
              kvm(k,i)   = leng(k,i) * sqrt(tke(k,i)) * sm
              bprod(k,i) = -kvh(k,i) * n2(k,i)
              sprod(k,i) =  kvm(k,i) * s2(k,i)

              sm_aw(k,i) = sm/alph1     ! This is diagnostic output for use in the microphysics             

          end if

       end do  ! k

       ! --------------------------------------------------- !
       ! End of treatment of Stable Turbulent Regime ( STL ) !
       ! --------------------------------------------------- !

       ! --------------------------------------------------------------- !
       ! Re-computation of eddy diffusivity at the entrainment interface !
       ! assuming that it is purely STL (0<Ri<0.19). Note even Ri>0.19,  !
       ! turbulent can exist at the entrainment interface since 'Sh,Sm'  !
       ! do not necessarily go to zero even when Ri>0.19. Since Ri can   !
       ! be fairly larger than 0.19 at the entrainment interface, I      !
       ! should set minimum value of 'tke' to be 0. in order to prevent  !
       ! sqrt(tke) from being imaginary.                                 !
       ! --------------------------------------------------------------- !

       ! goto 888

         do k = 2, pver

         if( ( turbtype(k,i) .eq. 3 ) .or. ( turbtype(k,i) .eq. 4 ) .or. &
             ( turbtype(k,i) .eq. 5 ) ) then

             trma = alph3*alph4exs*ri(k,i) + 2._r8*b1*(alph2-alph4exs*alph5*ri(k,i))
             trmb = (alph3+alph4exs)*ri(k,i) + 2._r8*b1*(-alph5*ri(k,i)+alph1)
             trmc = ri(k,i)
             det  = max(trmb*trmb-4._r8*trma*trmc,0._r8)
             gh   = (-trmb + sqrt(det))/(2._r8*trma)
           ! gh   = min(max(gh,-0.28_r8),0.0233_r8)
           ! gh   = min(max(gh,-3.5334_r8),0.0233_r8)
             gh   = min(max(gh,ghmin),0.0233_r8)
             sh   = max(0._r8,alph5/(1._r8+alph3*gh))
             sm   = max(0._r8,(alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))

             lbulk = z(k-1,i) - z(k,i)
             lbulk = min( lbulk, lbulk_max )

             if( choice_tunl .eq. 'rampcl' ) then
                 tunlramp = tunl
             elseif( choice_tunl .eq. 'rampsl' ) then
                 tunlramp = max(1.e-3_r8,ctunl*tunl*exp(-log(ctunl)*ri(k,i)/ricrit))
               ! tunlramp = 0.065_r8 + 0.7_r8*exp(-20._r8*ri(k,i))
             else
                 tunlramp = tunl
             endif
             if( choice_leng .eq. 'origin' ) then
                 leng_imsi = ( (vk*zi(k,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
               ! leng_imsi = vk*zi(k,i) / (1._r8+vk*zi(k,i)/(tunlramp*lbulk))
             else
                 leng_imsi = min( vk*zi(k,i), tunlramp*lbulk )              
             endif
             leng_imsi = min(leng_max(k), leng_imsi)

             tke_imsi = b1*(leng_imsi**2)*(-sh*n2(k,i)+sm*s2(k,i))
             tke_imsi = min(max(tke_imsi,0._r8),tkemax)
             kvh_imsi = leng_imsi * sqrt(tke_imsi) * sh
             kvm_imsi = leng_imsi * sqrt(tke_imsi) * sm

             if( kvh(k,i) .lt. kvh_imsi ) then 
                 kvh(k,i)   =  kvh_imsi
                 kvm(k,i)   =  kvm_imsi
                 leng(k,i)  = leng_imsi
                 tke(k,i)   =  tke_imsi
                 wcap(k,i)  =  tke_imsi / b1
                 bprod(k,i) = -kvh_imsi * n2(k,i)
                 sprod(k,i) =  kvm_imsi * s2(k,i)
                 sm_aw(k,i) =  sm/alph1     ! This is diagnostic output for use in the microphysics             
                 turbtype(k,i) = 1          ! This was added on Dec.10.2009 for use in microphysics.
             endif

         end if

         end do

 ! 888   continue 

       ! ------------------------------------------------------------------ !
       ! End of recomputation of eddy diffusivity at entrainment interfaces !
       ! ------------------------------------------------------------------ !

       ! As an option, we can impose a certain minimum back-ground diffusivity.

       ! do k = 1, pver+1
       !    kvh(k,i) = max(0.01_r8,kvh(k,i))
       !    kvm(k,i) = max(0.01_r8,kvm(k,i))
       ! enddo
       ! --------------------------------------------------------------------- !
       ! Diagnostic Output                                                     !
       ! Just for diagnostic purpose, calculate stability functions at  each   !
       ! interface including surface. Instead of assuming neutral stability,   !
       ! explicitly calculate stability functions using an reverse procedure   !
       ! starting from tkes(i) similar to the case of SRCL and SBCL in zisocl. !
       ! Note that it is possible to calculate stability functions even when   !
       ! bflxs < 0. Note that this inverse method allows us to define Ri even  !
       ! at the surface. Note also tkes(i) and sprod(pver+1,i) are always      !
       ! positive values by limiters (e.g., ustar_min = 0.01).                 !
       ! Dec.12.2006 : Also just for diagnostic output, re-set                 !
       ! 'bprod(pver+1,i)= bflxs(i)' here. Note that this setting does not     !
       ! influence numerical calculation at all - it is just for diagnostic    !
       ! output.                                                               !
       ! --------------------------------------------------------------------- !

       bprod(pver+1,i) = bflxs(i)
              
       gg = 0.5_r8*vk*z(pver,i)*bprod(pver+1,i)/(tkes(i)**(3._r8/2._r8))
       if( abs(alph5-gg*alph3) .le. 1.e-7_r8 ) then
         ! gh = -0.28_r8
           if( bprod(pver+1,i) .gt. 0._r8 ) then
               gh = -3.5334_r8
           else
               gh = ghmin
           endif
       else    
           gh = gg/(alph5-gg*alph3)
       end if 

     ! gh = min(max(gh,-0.28_r8),0.0233_r8)
       if( bprod(pver+1,i) .gt. 0._r8 ) then
           gh = min(max(gh,-3.5334_r8),0.0233_r8)
       else
           gh = min(max(gh,ghmin),0.0233_r8)
       endif

       gh_a(pver+1,i) = gh     
       sh_a(pver+1,i) = max(0._r8,alph5/(1._r8+alph3*gh))
       if( bprod(pver+1,i) .gt. 0._r8 ) then       
           sm_a(pver+1,i) = max(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh))
       else
           sm_a(pver+1,i) = max(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))
       endif
       sm_aw(pver+1,i) = sm_a(pver+1,i)/alph1
       ri_a(pver+1,i)  = -(sm_a(pver+1,i)/sh_a(pver+1,i))*(bprod(pver+1,i)/sprod(pver+1,i))

       do k = 1, pver
          if( ri(k,i) .lt. 0._r8 ) then
              trma = alph3*alph4*ri(k,i) + 2._r8*b1*(alph2-alph4*alph5*ri(k,i))
              trmb = (alph3+alph4)*ri(k,i) + 2._r8*b1*(-alph5*ri(k,i)+alph1)
              trmc = ri(k,i)
              det  = max(trmb*trmb-4._r8*trma*trmc,0._r8)
              gh   = (-trmb + sqrt(det))/(2._r8*trma)
              gh   = min(max(gh,-3.5334_r8),0.0233_r8)
              gh_a(k,i) = gh
              sh_a(k,i) = max(0._r8,alph5/(1._r8+alph3*gh))
              sm_a(k,i) = max(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh))
              ri_a(k,i) = ri(k,i)
          else
              if( ri(k,i) .gt. ricrit ) then
                  gh_a(k,i) = ghmin
                  sh_a(k,i) = 0._r8
                  sm_a(k,i) = 0._r8
                  ri_a(k,i) = ri(k,i)
              else
                  trma = alph3*alph4exs*ri(k,i) + 2._r8*b1*(alph2-alph4exs*alph5*ri(k,i))
                  trmb = (alph3+alph4exs)*ri(k,i) + 2._r8*b1*(-alph5*ri(k,i)+alph1)
                  trmc = ri(k,i)
                  det  = max(trmb*trmb-4._r8*trma*trmc,0._r8)
                  gh   = (-trmb + sqrt(det))/(2._r8*trma)
                  gh   = min(max(gh,ghmin),0.0233_r8)
                  gh_a(k,i) = gh
                  sh_a(k,i) = max(0._r8,alph5/(1._r8+alph3*gh))
                  sm_a(k,i) = max(0._r8,(alph1+alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4exs*gh))
                  ri_a(k,i) = ri(k,i)
              endif
          endif

       end do

    end do   ! End of column index loop, i 

    return

    end subroutine caleddy


! Purpose: Find unstable CL regimes and determine the indices
!          kbase, ktop which delimit these unstable layers :
!          ri(kbase) > 0 and ri(ktop) > 0, but ri(k) < 0 for ktop < k < kbase.
    subroutine exacol( pver, ncol, ri, bflxs, minpblh, zi, ktop, kbase, ncvfin ) 

    implicit none
! io
    integer,  intent(in) :: pver                   ! Number of atmospheric vertical layers   
    integer,  intent(in) :: ncol                   ! Number of atmospheric columns   
    real(r8), intent(in) :: ri(pver, ncol)         ! Moist gradient Richardson no.
    real(r8), intent(in) :: bflxs(ncol)            ! Buoyancy flux at surface
    real(r8), intent(in) :: minpblh(ncol)          ! Minimum PBL height based on surface stress
    real(r8), intent(in) :: zi(pver+1, ncol)       ! Interface heights

    integer, intent(out) :: kbase(ncvmax, ncol)    ! External interface index of CL base
    integer, intent(out) :: ktop(ncvmax, ncol)     ! External interface index of CL top
    integer, intent(out) :: ncvfin(ncol)           ! Total number of CLs
! local
    integer              :: i
    integer              :: k
    integer              :: ncv
    real(r8)             :: rimaxentr
    real(r8)             :: riex(pver+1)           ! Column Ri profile extended to surface

    do i = 1, ncol
       ncvfin(i) = 0
       do ncv = 1, ncvmax
          ktop(ncv,i)  = 0
          kbase(ncv,i) = 0
       end do
    end do

    ! Find CL regimes starting from the surface going upward !
    rimaxentr = 0._r8   
    do i = 1, ncol
       riex(2:pver) = ri(2:pver,i)

       ! Below allows consistent treatment of surface and other interfaces.
       ! Simply, if surface buoyancy flux is positive, Ri of surface is set to be negative.
       riex(pver+1) = rimaxentr - bflxs(i) 
       ncv = 0
       k   = pver + 1 ! Work upward from surface interface
       do while ( k .gt. ntop_turb + 1 )

        ! Below means that if 'bflxs > 0' (do not contain '=' sign), surface
        ! interface is energetically interior surface. 
          if( riex(k) .lt. rimaxentr ) then 

              ! Identify a new CL
              ncv = ncv + 1

              ! First define 'kbase' as the first interface below the lower-most unstable interface
              ! Thus, Richardson number at 'kbase' is positive.
              kbase(ncv,i) = min(k+1,pver+1)

              ! Decrement k until top unstable level
              do while( riex(k) .lt. rimaxentr .and. k .gt. ntop_turb + 1 )
                 k = k - 1
              end do

              ! ktop is the first interface above upper-most unstable interface
              ! Thus, Richardson number at 'ktop' is positive. 
              ktop(ncv,i) = k
          else

              ! Search upward for a CL.
              k = k - 1

          end if
       end do ! End of CL regime finding for each atmospheric column

       ncvfin(i) = ncv    

    end do  ! End of atmospheric column do loop

    end subroutine exacol


! Purpose: This 'zisocl' vertically extends original CLs identified from
!          'exacol' using a merging test based on either 'wint' or 'l2n2'
!          and identify new CL regimes. Similar to the case of 'exacol',
!          CL regime index increases with height.  After identifying new
!          CL regimes ( kbase, ktop, ncvfin ),calculate CL internal mean
!          energetics (lbrk : energetic thickness integral, wbrk, ebrk )
!          and stability functions (ricl, ghcl, shcl, smcl) by including
!          surface interfacial layer contribution when bflxs > 0.   Note
!          that there are two options in the treatment of the energetics
!          of surface interfacial layer (use_dw_surf= 'true' or 'false')
    subroutine zisocl( ncol   , pver  , long ,                                 & 
                       z      , zi    , n2   ,  s2      ,                      & 
                       bprod  , sprod , bflxs,  tkes    ,                      & 
                       ncvfin , kbase , ktop ,  belongcv,                      & 
                       ricl   , ghcl  , shcl ,  smcl    ,                      &
                       lbrk   , wbrk  , ebrk ,  extend  , extend_up, extend_dn )

    implicit none
! io
    integer,  intent(in)   :: long                    ! Longitude of the column
    integer,  intent(in)   :: ncol                    ! Number of atmospheric columns   
    integer,  intent(in)   :: pver                    ! Number of atmospheric vertical layers   
    real(r8), intent(in)   :: z(pver, ncol)           ! Layer mid-point height [ m ]
    real(r8), intent(in)   :: zi(pver+1, ncol)        ! Interface height [ m ]
    real(r8), intent(in)   :: n2(pver, ncol)          ! Buoyancy frequency at interfaces except surface [ s-2 ]
    real(r8), intent(in)   :: s2(pver, ncol)          ! Shear frequency at interfaces except surface [ s-2 ]
    real(r8), intent(in)   :: bprod(pver+1,ncol)      ! Buoyancy production [ m2/s3 ]. bprod(pver+1,i) = bflxs 
    real(r8), intent(in)   :: sprod(pver+1,ncol)      ! Shear production [ m2/s3 ]. sprod(pver+1,i) = usta**3/(vk*z(pver,i))
    real(r8), intent(in)   :: bflxs(ncol)             ! Surface buoyancy flux [ m2/s3 ]. bprod(pver+1,i) = bflxs 
    real(r8), intent(in)   :: tkes(ncol)              ! TKE at the surface [ s2/s2 ]

    integer, intent(inout) :: kbase(ncvmax, ncol)     ! Base external interface index of CL
    integer, intent(inout) :: ktop(ncvmax, ncol)      ! Top external interface index of CL
    integer, intent(inout) :: ncvfin(ncol)            ! Total number of CLs

    ! ---------------- !
    ! Output variables !
    ! ---------------- !
    logical,  intent(out) :: belongcv(pver+1,ncol)   ! True if interface is in a CL ( either internal or external )
    real(r8), intent(out) :: ricl(ncvmax, ncol)       ! Mean Richardson number of internal CL
    real(r8), intent(out) :: ghcl(ncvmax, ncol)       ! Half of normalized buoyancy production of internal CL
    real(r8), intent(out) :: shcl(ncvmax, ncol)       ! Galperin instability function of heat-moisture of internal CL
    real(r8), intent(out) :: smcl(ncvmax, ncol)       ! Galperin instability function of momentum of internal CL
    real(r8), intent(out) :: lbrk(ncvmax, ncol)       ! Thickness of (energetically) internal CL ( lint, [m] )
    real(r8), intent(out) :: wbrk(ncvmax, ncol)       ! Mean normalized TKE of internal CL  [ m2/s2 ]
    real(r8), intent(out) :: ebrk(ncvmax, ncol)       ! Mean TKE of internal CL ( b1*wbrk, [m2/s2] )

    ! ------------------ !
    ! Internal variables !
    ! ------------------ !
    logical               :: extend                   ! True when CL is extended in zisocl
    logical               :: extend_up                ! True when CL is extended upward in zisocl
    logical               :: extend_dn                ! True when CL is extended downward in zisocl
    logical               :: bottom                   ! True when CL base is at surface ( kb = pver + 1 )

    integer               :: i                        ! Local index for the longitude
    integer               :: ncv                      ! CL Index increasing with height
    integer               :: incv
    integer               :: k
    integer               :: kb                       ! Local index for kbase
    integer               :: kt                       ! Local index for ktop
    integer               :: ncvinit                  ! Value of ncv at routine entrance 
    integer               :: cntu                     ! Number of merged CLs during upward   extension of individual CL
    integer               :: cntd                     ! Number of merged CLs during downward extension of individual CL
    integer               :: kbinc                    ! Index for incorporating underlying CL
    integer               :: ktinc                    ! Index for incorporating  overlying CL

    real(r8)              :: wint                     ! Normalized TKE of internal CL
    real(r8)              :: dwinc                    ! Normalized TKE of CL external interfaces
    real(r8)              :: dw_surf                  ! Normalized TKE of surface interfacial layer
    real(r8)              :: dzinc
    real(r8)              :: gh
    real(r8)              :: sh
    real(r8)              :: sm
    real(r8)              :: gh_surf                  ! Half of normalized buoyancy production in surface interfacial layer 
    real(r8)              :: sh_surf                  ! Galperin instability function in surface interfacial layer  
    real(r8)              :: sm_surf                  ! Galperin instability function in surface interfacial layer 
    real(r8)              :: l2n2                     ! Vertical integral of 'l^2N^2' over CL. Include thickness product
    real(r8)              :: l2s2                     ! Vertical integral of 'l^2S^2' over CL. Include thickness product
    real(r8)              :: dl2n2                    ! Vertical integration of 'l^2*N^2' of CL external interfaces
    real(r8)              :: dl2s2                    ! Vertical integration of 'l^2*S^2' of CL external interfaces
    real(r8)              :: dl2n2_surf               ! 'dl2n2' defined in the surface interfacial layer
    real(r8)              :: dl2s2_surf               ! 'dl2s2' defined in the surface interfacial layer  
    real(r8)              :: lint                     ! Thickness of (energetically) internal CL
    real(r8)              :: dlint                    ! Interfacial layer thickness of CL external interfaces
    real(r8)              :: dlint_surf               ! Surface interfacial layer thickness 
    real(r8)              :: lbulk                    ! Master Length Scale : Whole CL thickness from top to base external interface
    real(r8)              :: lz                       ! Turbulent length scale
    real(r8)              :: ricll                    ! Mean Richardson number of internal CL 
    real(r8)              :: trma
    real(r8)              :: trmb
    real(r8)              :: trmc
    real(r8)              :: det
    real(r8)              :: zbot                     ! Height of CL base
    real(r8)              :: l2rat                    ! Square of ratio of actual to initial CL (not used)
    real(r8)              :: gg                       ! Intermediate variable used for calculating stability functions of SBCL
    real(r8)              :: tunlramp                 ! Ramping tunl

    i = long

    ! Initialize main output variables
    
    do k = 1, ncvmax
       ricl(k,i) = 0._r8
       ghcl(k,i) = 0._r8
       shcl(k,i) = 0._r8
       smcl(k,i) = 0._r8
       lbrk(k,i) = 0._r8
       wbrk(k,i) = 0._r8
       ebrk(k,i) = 0._r8
    end do
    extend    = .false.
    extend_up = .false.
    extend_dn = .false.

    ! ----------------------------------------------------------- !
    ! Loop over each CL to see if any of them need to be extended !
    ! ----------------------------------------------------------- !

    ncv = 1

    do while( ncv .le. ncvfin(i) )

       ncvinit = ncv
       cntu    = 0
       cntd    = 0
       kb      = kbase(ncv,i) 
       kt      = ktop(ncv,i)
       
       ! ---------------------------------------------------------------------------- !
       ! Calculation of CL interior energetics including surface before extension     !
       ! ---------------------------------------------------------------------------- !
       ! Note that the contribution of interior interfaces (not surface) to 'wint' is !
       ! accounted by using '-sh*l2n2 + sm*l2s2' while the contribution of surface is !
       ! accounted by using 'dwsurf = tkes/b1' when bflxs > 0. This approach is fully !
       ! reasonable. Another possible alternative,  which seems to be also consistent !
       ! is to calculate 'dl2n2_surf'  and  'dl2s2_surf' of surface interfacial layer !
       ! separately, and this contribution is explicitly added by initializing 'l2n2' !
       ! 'l2s2' not by zero, but by 'dl2n2_surf' and 'ds2n2_surf' below.  At the same !
       ! time, 'dwsurf' should be excluded in 'wint' calculation below. The only diff.!
       ! between two approaches is that in case of the latter approach, contributions !
       ! of surface interfacial layer to the CL mean stability function (ri,gh,sh,sm) !
       ! are explicitly included while the first approach is not. In this sense,  the !
       ! second approach seems to be more conceptually consistent,   but currently, I !
       ! (Sungsu) will keep the first default approach. There is a switch             !
       ! 'use_dw_surf' at the first part of eddy_diff.F90 chosing one of              !
       ! these two options.                                                           !
       ! ---------------------------------------------------------------------------- !

       ! ------------------------------------------------------ !
       ! Step 0: Calculate surface interfacial layer energetics !
       ! ------------------------------------------------------ !

       lbulk      = zi(kt,i) - zi(kb,i)
       lbulk      = min( lbulk, lbulk_max )
       dlint_surf = 0._r8
       dl2n2_surf = 0._r8
       dl2s2_surf = 0._r8
       dw_surf    = 0._r8
       if( kb .eq. pver+1 ) then

           if( bflxs(i) .gt. 0._r8 ) then

               ! Calculate stability functions of surface interfacial layer
               ! from the given 'bprod(pver+1,i)' and 'sprod(pver+1,i)' using
               ! inverse approach. Since alph5>0 and alph3<0, denominator of
               ! gg is always positive if bprod(pver+1,i)>0.               

               gg    = 0.5_r8*vk*z(pver,i)*bprod(pver+1,i)/(tkes(i)**(3._r8/2._r8))
               gh    = gg/(alph5-gg*alph3)
             ! gh    = min(max(gh,-0.28_r8),0.0233_r8)
               gh    = min(max(gh,-3.5334_r8),0.0233_r8)
               sh    = alph5/(1._r8+alph3*gh)
               sm    = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
               ricll = min(-(sm/sh)*(bprod(pver+1,i)/sprod(pver+1,i)),ricrit)

               ! Calculate surface interfacial layer contribution to CL internal
               ! energetics. By construction, 'dw_surf = -dl2n2_surf + ds2n2_surf'
               ! is exactly satisfied, which corresponds to assuming turbulent
               ! length scale of surface interfacial layer = vk * z(pver,i). Note
               ! 'dl2n2_surf','dl2s2_surf','dw_surf' include thickness product.   

               dlint_surf = z(pver,i)
               dl2n2_surf = -vk*(z(pver,i)**2)*bprod(pver+1,i)/(sh*sqrt(tkes(i)))
               dl2s2_surf =  vk*(z(pver,i)**2)*sprod(pver+1,i)/(sm*sqrt(tkes(i)))
               dw_surf    = (tkes(i)/b1)*z(pver,i) 

           else

               ! Note that this case can happen when surface is an external 
               ! interface of CL.
               lbulk = zi(kt,i) - z(pver,i)
               lbulk = min( lbulk, lbulk_max )

           end if

       end if

       ! ------------------------------------------------------ !   
       ! Step 1: Include surface interfacial layer contribution !
       ! ------------------------------------------------------ !
       
       lint = dlint_surf
       l2n2 = dl2n2_surf
       l2s2 = dl2s2_surf          
       wint = dw_surf
       if( use_dw_surf ) then
           l2n2 = 0._r8
           l2s2 = 0._r8
       else
           wint = 0._r8
       end if    
       
       ! --------------------------------------------------------------------------------- !
       ! Step 2. Include the contribution of 'pure internal interfaces' other than surface !
       ! --------------------------------------------------------------------------------- ! 
       
       if( kt .lt. kb - 1 ) then ! The case of non-SBCL.
                              
           do k = kb - 1, kt + 1, -1       
              if( choice_tunl .eq. 'rampcl' ) then
                ! Modification : I simply used the average tunlramp between the two limits.
                  tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
              elseif( choice_tunl .eq. 'rampsl' ) then
                  tunlramp = ctunl*tunl
                ! tunlramp = 0.765_r8
              else
                  tunlramp = tunl
              endif
              if( choice_leng .eq. 'origin' ) then
                  lz = ( (vk*zi(k,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! lz = vk*zi(k,i) / (1._r8+vk*zi(k,i)/(tunlramp*lbulk))
              else
                  lz = min( vk*zi(k,i), tunlramp*lbulk )              
              endif
              lz = min(leng_max(k), lz)
              dzinc = z(k-1,i) - z(k,i)
              l2n2  = l2n2 + lz*lz*n2(k,i)*dzinc
              l2s2  = l2s2 + lz*lz*s2(k,i)*dzinc
              lint  = lint + dzinc
           end do

           ! Calculate initial CL stability functions (gh,sh,sm) and net
           ! internal energy of CL including surface contribution if any. 

           ! Modification : It seems that below cannot be applied when ricrit > 0.19.
           !                May need future generalization.

           ricll = min(l2n2/max(l2s2,ntzero),ricrit) ! Mean Ri of internal CL
           trma  = alph3*alph4*ricll+2._r8*b1*(alph2-alph4*alph5*ricll)
           trmb  = ricll*(alph3+alph4)+2._r8*b1*(-alph5*ricll+alph1)
           trmc  = ricll
           det   = max(trmb*trmb-4._r8*trma*trmc,0._r8)
           gh    = (-trmb + sqrt(det))/2._r8/trma
         ! gh    = min(max(gh,-0.28_r8),0.0233_r8)
           gh    = min(max(gh,-3.5334_r8),0.0233_r8)
           sh    = alph5/(1._r8+alph3*gh)
           sm    = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
           wint  = wint - sh*l2n2 + sm*l2s2 

       else ! The case of SBCL
 
           ! If there is no pure internal interface, use only surface interfacial
           ! values. However, re-set surface interfacial values such  that it can
           ! be used in the merging tests (either based on 'wint' or 'l2n2')  and
           ! in such that surface interfacial energy is not double-counted.
           ! Note that regardless of the choise of 'use_dw_surf', below should be
           ! kept as it is below, for consistent merging test of extending SBCL. 
       
           lint = dlint_surf
           l2n2 = dl2n2_surf
           l2s2 = dl2s2_surf 
           wint = dw_surf

           ! Aug.29.2006 : Only for the purpose of merging test of extending SRCL
           ! based on 'l2n2', re-define 'l2n2' of surface interfacial layer using
           ! 'wint'. This part is designed for similar treatment of merging as in
           ! the original 'eddy_diff.F90' code,  where 'l2n2' of SBCL was defined
           ! as 'l2n2 = - wint / sh'. Note that below block is used only when (1)
           ! surface buoyancy production 'bprod(pver+1,i)' is NOT included in the
           ! calculation of surface TKE in the initialization of 'bprod(pver+1,i)'
           ! in the main subroutine ( even though bflxs > 0 ), and (2) to force 
           ! current scheme be similar to the previous scheme in the treatment of  
           ! extending-merging test of SBCL based on 'l2n2'. Otherwise below line
           ! must be commented out. Note at this stage, correct non-zero value of
           ! 'sh' has been already computed.      

           if( choice_tkes .eq. 'ebprod' ) then
               l2n2 = - wint / sh 
           endif
           
       endif

       ! Set consistent upper limits on 'l2n2' and 'l2s2'. Below limits are
       ! reasonable since l2n2 of CL interior interface is always negative.

       l2n2 = -min(-l2n2, tkemax*lint/(b1*sh))
       l2s2 =  min( l2s2, tkemax*lint/(b1*sm))
       
       ! Note that at this stage, ( gh, sh, sm )  are the values of surface
       ! interfacial layer if there is no pure internal interface, while if
       ! there is pure internal interface, ( gh, sh, sm ) are the values of
       ! pure CL interfaces or the values that include both the CL internal
       ! interfaces and surface interfaces, depending on the 'use_dw_surf'.       
       
       ! ----------------------------------------------------------------------- !
       ! Perform vertical extension-merging process                              !
       ! ----------------------------------------------------------------------- !
       ! During the merging process, we assumed ( lbulk, sh, sm ) of CL external !
       ! interfaces are the same as the ones of the original merging CL. This is !
       ! an inevitable approximation since we don't know  ( sh, sm ) of external !
       ! interfaces at this stage.     Note that current default merging test is !
       ! purely based on buoyancy production without including shear production, !
       ! since we used 'l2n2' instead of 'wint' as a merging parameter. However, !
       ! merging test based on 'wint' maybe conceptually more attractable.       !
       ! Downward CL merging process is identical to the upward merging process, !
       ! but when the base of extended CL reaches to the surface, surface inter  !
       ! facial layer contribution to the energetic of extended CL must be done  !
       ! carefully depending on the sign of surface buoyancy flux. The contribu  !
       ! tion of surface interfacial layer energetic is included to the internal !
       ! energetics of merging CL only when bflxs > 0.                           !
       ! ----------------------------------------------------------------------- !
       
       ! ---------------------------- !
       ! Step 1. Extend the CL upward !
       ! ---------------------------- !
       
       extend = .false.    ! This will become .true. if CL top or base is extended

       ! Calculate contribution of potentially incorporable CL top interface

       if( choice_tunl .eq. 'rampcl' ) then
           tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
       elseif( choice_tunl .eq. 'rampsl' ) then
           tunlramp = ctunl*tunl
         ! tunlramp = 0.765_r8
       else
           tunlramp = tunl
       endif
       if( choice_leng .eq. 'origin' ) then
           lz = ( (vk*zi(kt,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
         ! lz = vk*zi(kt,i) / (1._r8+vk*zi(kt,i)/(tunlramp*lbulk))
       else
           lz = min( vk*zi(kt,i), tunlramp*lbulk )              
       endif
       lz = min(leng_max(kt), lz)

       dzinc = z(kt-1,i)-z(kt,i)
       dl2n2 = lz*lz*n2(kt,i)*dzinc
       dl2s2 = lz*lz*s2(kt,i)*dzinc
       dwinc = -sh*dl2n2 + sm*dl2s2

       ! ------------ !
       ! Merging Test !
       ! ------------ !

       ! The part of the below test that involves kt and z has different
       ! effects based on the model top.
       ! If the model top is in the stratosphere, we want the loop to
       ! continue until it either completes normally, or kt is pushed to
       ! the top of the model. The latter case should not happen, so this
       ! causes an error.
       ! If the model top is higher, as in WACCM and WACCM-X, if kt is
       ! pushed close to the model top, this may not represent an error at
       ! all, because of very different and more variable temperature/wind
       ! profiles at the model top. Therefore we simply exit the loop early
       ! and continue with no errors.

     ! do while (  dwinc .gt. ( rinc*dzinc*wint/(lint+(1._r8-rinc)*dzinc)) )  ! Merging test based on wint
     ! do while ( -dl2n2 .gt. (-rinc*dzinc*l2n2/(lint+(1._r8-rinc)*dzinc)) )  ! Merging test based on l2n2 

       do while ( -dl2n2 .gt. (-rinc*l2n2/(1._r8-rinc)) &                     ! Integral merging test
            .and. (kt > ntop_turb+2 .or. z(kt,i) < 50000._r8) )

          ! Add contribution of top external interface to interior energy.
          ! Note even when we chose 'use_dw_surf='true.', the contribution
          ! of surface interfacial layer to 'l2n2' and 'l2s2' are included
          ! here. However it is not double counting of surface interfacial
          ! energy : surface interfacial layer energy is counted in 'wint'
          ! formula and 'l2n2' is just used for performing merging test in
          ! this 'do while' loop.     

          lint = lint + dzinc
          l2n2 = l2n2 + dl2n2
          l2n2 = -min(-l2n2, tkemax*lint/(b1*sh))
          l2s2 = l2s2 + dl2s2
          wint = wint + dwinc

          ! Extend top external interface of CL upward after merging
 
          kt        = kt - 1
          extend    = .true.
          extend_up = .true.
          if( kt .eq. ntop_turb ) then
              print*,'In caleddy: in ncvfin: in zisocl, kt=ntop_turb',kt,ntop_turb,'z=',z(kt,i)
              call endrun('zisocl: Error: Tried to extend CL to the model top')
          end if

          ! If the top external interface of extending CL is the same as the 
          ! top interior interface of the overlying CL, overlying CL will be
          ! automatically merged. Then,reduce total number of CL regime by 1. 
          ! and increase 'cntu'(number of merged CLs during upward extension)
          ! by 1.
 
          ktinc = kbase(ncv+cntu+1,i) - 1  ! Lowest interior interface of overlying CL

          if( kt .eq. ktinc ) then

              do k = kbase(ncv+cntu+1,i) - 1, ktop(ncv+cntu+1,i) + 1, -1

                 if( choice_tunl .eq. 'rampcl' ) then
                     tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
                 elseif( choice_tunl .eq. 'rampsl' ) then
                     tunlramp = ctunl*tunl
                   ! tunlramp = 0.765_r8
                 else
                     tunlramp = tunl
                 endif
                 if( choice_leng .eq. 'origin' ) then
                     lz = ( (vk*zi(k,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                   ! lz = vk*zi(k,i) / (1._r8+vk*zi(k,i)/(tunlramp*lbulk))
                 else
                     lz = min( vk*zi(k,i), tunlramp*lbulk )              
                 endif
                 lz = min(leng_max(k), lz)

                 dzinc = z(k-1,i)-z(k,i)
                 dl2n2 = lz*lz*n2(k,i)*dzinc
                 dl2s2 = lz*lz*s2(k,i)*dzinc
                 dwinc = -sh*dl2n2 + sm*dl2s2

                 lint = lint + dzinc
                 l2n2 = l2n2 + dl2n2
                 l2n2 = -min(-l2n2, tkemax*lint/(b1*sh))
                 l2s2 = l2s2 + dl2s2
                 wint = wint + dwinc

              end do 

              kt        = ktop(ncv+cntu+1,i) 
              ncvfin(i) = ncvfin(i) - 1
              cntu      = cntu + 1
        
          end if

          ! Again, calculate the contribution of potentially incorporatable CL
          ! top external interface of CL regime.

          if( choice_tunl .eq. 'rampcl' ) then
              tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
          elseif( choice_tunl .eq. 'rampsl' ) then
              tunlramp = ctunl*tunl
            ! tunlramp = 0.765_r8
          else
              tunlramp = tunl
          endif
          if( choice_leng .eq. 'origin' ) then
              lz = ( (vk*zi(kt,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
            ! lz = vk*zi(kt,i) / (1._r8+vk*zi(kt,i)/(tunlramp*lbulk))
          else
              lz = min( vk*zi(kt,i), tunlramp*lbulk )              
          endif
          lz = min(leng_max(kt), lz)

          dzinc = z(kt-1,i)-z(kt,i)
          dl2n2 = lz*lz*n2(kt,i)*dzinc
          dl2s2 = lz*lz*s2(kt,i)*dzinc
          dwinc = -sh*dl2n2 + sm*dl2s2

       end do   ! End of upward merging test 'do while' loop

       ! Update CL interface indices appropriately if any CL was merged.
       ! Note that below only updated the interface index of merged CL,
       ! not the original merging CL.  Updates of 'kbase' and 'ktop' of 
       ! the original merging CL  will be done after finishing downward
       ! extension also later.

       if( cntu .gt. 0 ) then
           do incv = 1, ncvfin(i) - ncv
              kbase(ncv+incv,i) = kbase(ncv+cntu+incv,i)
              ktop(ncv+incv,i)  = ktop(ncv+cntu+incv,i)
           end do
       end if

       ! ------------------------------ !
       ! Step 2. Extend the CL downward !
       ! ------------------------------ !

        if( kb .ne. pver + 1 ) then

           ! Calculate contribution of potentially incorporable CL base interface

           if( choice_tunl .eq. 'rampcl' ) then
               tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
           elseif( choice_tunl .eq. 'rampsl' ) then
               tunlramp = ctunl*tunl
             ! tunlramp = 0.765_r8
           else
               tunlramp = tunl
           endif
           if( choice_leng .eq. 'origin' ) then
               lz = ( (vk*zi(kb,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
             ! lz = vk*zi(kb,i) / (1._r8+vk*zi(kb,i)/(tunlramp*lbulk))
           else
               lz = min( vk*zi(kb,i), tunlramp*lbulk )              
           endif
           lz = min(leng_max(kb), lz)

           dzinc = z(kb-1,i)-z(kb,i)
           dl2n2 = lz*lz*n2(kb,i)*dzinc
           dl2s2 = lz*lz*s2(kb,i)*dzinc
           dwinc = -sh*dl2n2 + sm*dl2s2

           ! ------------ ! 
           ! Merging test !
           ! ------------ ! 

           ! In the below merging tests, I must keep '.and.(kb.ne.pver+1)',   
           ! since 'kb' is continuously updated within the 'do while' loop  
           ! whenever CL base is merged.

         ! do while( (  dwinc .gt. ( rinc*dzinc*wint/(lint+(1._r8-rinc)*dzinc)) ) &  ! Merging test based on wint
         ! do while( ( -dl2n2 .gt. (-rinc*dzinc*l2n2/(lint+(1._r8-rinc)*dzinc)) ) &  ! Merging test based on l2n2
         !             .and.(kb.ne.pver+1))
           do while( ( -dl2n2 .gt. (-rinc*l2n2/(1._r8-rinc)) ) &                     ! Integral merging test
                       .and.(kb.ne.pver+1))

              ! Add contributions from interfacial layer kb to CL interior 

              lint = lint + dzinc
              l2n2 = l2n2 + dl2n2
              l2n2 = -min(-l2n2, tkemax*lint/(b1*sh))
              l2s2 = l2s2 + dl2s2
              wint = wint + dwinc

              ! Extend the base external interface of CL downward after merging

              kb        =  kb + 1
              extend    = .true.
              extend_dn = .true.

              ! If the base external interface of extending CL is the same as the 
              ! base interior interface of the underlying CL, underlying CL  will
              ! be automatically merged. Then, reduce total number of CL by 1. 
              ! For a consistent treatment with 'upward' extension,  I should use
              ! 'kbinc = kbase(ncv-1,i) - 1' instead of 'ktop(ncv-1,i) + 1' below.
              ! However, it seems that these two methods produce the same results.
              ! Note also that in contrast to upward merging, the decrease of ncv
              ! should be performed here.
              ! Note that below formula correctly works even when upperlying CL 
              ! regime incorporates below SBCL.

              kbinc = 0
              if( ncv .gt. 1 ) kbinc = ktop(ncv-1,i) + 1
              if( kb .eq. kbinc ) then

                  do k =  ktop(ncv-1,i) + 1, kbase(ncv-1,i) - 1

                     if( choice_tunl .eq. 'rampcl' ) then
                         tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
                     elseif( choice_tunl .eq. 'rampsl' ) then
                         tunlramp = ctunl*tunl
                       ! tunlramp = 0.765_r8
                     else
                         tunlramp = tunl
                     endif
                     if( choice_leng .eq. 'origin' ) then
                         lz = ( (vk*zi(k,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                       ! lz = vk*zi(k,i) / (1._r8+vk*zi(k,i)/(tunlramp*lbulk))
                     else
                         lz = min( vk*zi(k,i), tunlramp*lbulk )              
                     endif
                     lz = min(leng_max(k), lz)

                     dzinc = z(k-1,i)-z(k,i)
                     dl2n2 = lz*lz*n2(k,i)*dzinc
                     dl2s2 = lz*lz*s2(k,i)*dzinc
                     dwinc = -sh*dl2n2 + sm*dl2s2

                     lint = lint + dzinc
                     l2n2 = l2n2 + dl2n2
                     l2n2 = -min(-l2n2, tkemax*lint/(b1*sh))
                     l2s2 = l2s2 + dl2s2
                     wint = wint + dwinc

                  end do 

                  ! We are incorporating interior of CL ncv-1, so merge
                  ! this CL into the current CL.

                  kb        = kbase(ncv-1,i)
                  ncv       = ncv - 1
                  ncvfin(i) = ncvfin(i) -1
                  cntd      = cntd + 1

              end if

              ! Calculate the contribution of potentially incorporatable CL
              ! base external interface. Calculate separately when the base
              ! of extended CL is surface and non-surface.
             
              if( kb .eq. pver + 1 ) then 

                  if( bflxs(i) .gt. 0._r8 ) then 
                      ! Calculate stability functions of surface interfacial layer
                      gg = 0.5_r8*vk*z(pver,i)*bprod(pver+1,i)/(tkes(i)**(3._r8/2._r8))
                      gh_surf = gg/(alph5-gg*alph3)
                    ! gh_surf = min(max(gh_surf,-0.28_r8),0.0233_r8)
                      gh_surf = min(max(gh_surf,-3.5334_r8),0.0233_r8)
                      sh_surf = alph5/(1._r8+alph3*gh_surf)
                      sm_surf = (alph1 + alph2*gh_surf)/(1._r8+alph3*gh_surf)/(1._r8+alph4*gh_surf)
                      ! Calculate surface interfacial layer contribution. By construction,
                      ! it exactly becomes 'dw_surf = -dl2n2_surf + ds2n2_surf'  
                      dlint_surf = z(pver,i)
                      dl2n2_surf = -vk*(z(pver,i)**2._r8)*bprod(pver+1,i)/(sh_surf*sqrt(tkes(i)))
                      dl2s2_surf =  vk*(z(pver,i)**2._r8)*sprod(pver+1,i)/(sm_surf*sqrt(tkes(i)))
                      dw_surf = (tkes(i)/b1)*z(pver,i) 
                  else
                      dlint_surf = 0._r8
                      dl2n2_surf = 0._r8
                      dl2s2_surf = 0._r8
                      dw_surf = 0._r8
                  end if
                  ! If (kb.eq.pver+1), updating of CL internal energetics should be 
                  ! performed here inside of 'do while' loop, since 'do while' loop
                  ! contains the constraint of '.and.(kb.ne.pver+1)',so updating of
                  ! CL internal energetics cannot be performed within this do while
                  ! loop when kb.eq.pver+1. Even though I updated all 'l2n2','l2s2',
                  ! 'wint' below, only the updated 'wint' is used in the following
                  ! numerical calculation.                
                  lint = lint + dlint_surf
                  l2n2 = l2n2 + dl2n2_surf
                  l2n2 = -min(-l2n2, tkemax*lint/(b1*sh))
                  l2s2 = l2s2 + dl2s2_surf 
                  wint = wint + dw_surf                

               else

                  if( choice_tunl .eq. 'rampcl' ) then
                      tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
                  elseif( choice_tunl .eq. 'rampsl' ) then
                      tunlramp = ctunl*tunl
                    ! tunlramp = 0.765_r8
                  else
                      tunlramp = tunl
                  endif
                  if( choice_leng .eq. 'origin' ) then
                      lz = ( (vk*zi(kb,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                    ! lz = vk*zi(kb,i) / (1._r8+vk*zi(kb,i)/(tunlramp*lbulk))
                  else
                      lz = min( vk*zi(kb,i), tunlramp*lbulk )              
                  endif
                  lz = min(leng_max(kb), lz)

                  dzinc = z(kb-1,i)-z(kb,i)
                  dl2n2 = lz*lz*n2(kb,i)*dzinc
                  dl2s2 = lz*lz*s2(kb,i)*dzinc
                  dwinc = -sh*dl2n2 + sm*dl2s2

              end if

          end do ! End of merging test 'do while' loop

          if( (kb.eq.pver+1) .and. (ncv.ne.1) ) then 
               call endrun('Major mistake zisocl: the CL based at surface is not indexed 1')
          end if

       end if   ! Done with bottom extension of CL 

       ! Update CL interface indices appropriately if any CL was merged.
       ! Note that below only updated the interface index of merged CL,
       ! not the original merging CL.  Updates of 'kbase' and 'ktop' of 
       ! the original merging CL  will be done later below. I should 
       ! check in detail if below index updating is correct or not.   

        if( cntd .gt. 0 ) then
           do incv = 1, ncvfin(i) - ncv
              kbase(ncv+incv,i) = kbase(ncvinit+incv,i)
              ktop(ncv+incv,i)  = ktop(ncvinit+incv,i)
           end do
       end if

       ! Sanity check for positive wint.

       if( wint .lt. 0.01_r8 ) then
           wint = 0.01_r8
       end if

       ! -------------------------------------------------------------------------- !
       ! Finally update CL mean internal energetics including surface contribution  !
       ! after finishing all the CL extension-merging process.  As mentioned above, !
       ! there are two possible ways in the treatment of surface interfacial layer, !
       ! either through 'dw_surf' or 'dl2n2_surf and dl2s2_surf' by setting logical !
       ! variable 'use_dw_surf' =.true. or .false.    In any cases, we should avoid !
       ! double counting of surface interfacial layer and one single consistent way !
       ! should be used throughout the program.                                     !
       ! -------------------------------------------------------------------------- !

       if( extend ) then

           ktop(ncv,i)  = kt
           kbase(ncv,i) = kb

           ! ------------------------------------------------------ !   
           ! Step 1: Include surface interfacial layer contribution !
           ! ------------------------------------------------------ !        
          
           lbulk      = zi(kt,i) - zi(kb,i)
           lbulk      = min( lbulk, lbulk_max )
           dlint_surf = 0._r8
           dl2n2_surf = 0._r8
           dl2s2_surf = 0._r8
           dw_surf    = 0._r8
           if( kb .eq. pver + 1 ) then
               if( bflxs(i) .gt. 0._r8 ) then
                   ! Calculate stability functions of surface interfacial layer
                   gg = 0.5_r8*vk*z(pver,i)*bprod(pver+1,i)/(tkes(i)**(3._r8/2._r8))
                   gh = gg/(alph5-gg*alph3)
                 ! gh = min(max(gh,-0.28_r8),0.0233_r8)
                   gh = min(max(gh,-3.5334_r8),0.0233_r8)
                   sh = alph5/(1._r8+alph3*gh)
                   sm = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
                   ! Calculate surface interfacial layer contribution. By construction,
                   ! it exactly becomes 'dw_surf = -dl2n2_surf + ds2n2_surf'  
                   dlint_surf = z(pver,i)
                   dl2n2_surf = -vk*(z(pver,i)**2._r8)*bprod(pver+1,i)/(sh*sqrt(tkes(i)))
                   dl2s2_surf =  vk*(z(pver,i)**2._r8)*sprod(pver+1,i)/(sm*sqrt(tkes(i)))
                   dw_surf    = (tkes(i)/b1)*z(pver,i) 
               else
                   lbulk = zi(kt,i) - z(pver,i)
                   lbulk = min( lbulk, lbulk_max )
               end if
           end if
           lint = dlint_surf
           l2n2 = dl2n2_surf
           l2s2 = dl2s2_surf
           wint = dw_surf
           if( use_dw_surf ) then
               l2n2 = 0._r8
               l2s2 = 0._r8
           else
               wint = 0._r8
           end if   
       
           ! -------------------------------------------------------------- !
           ! Step 2. Include the contribution of 'pure internal interfaces' !
           ! -------------------------------------------------------------- ! 
          
           do k = kt + 1, kb - 1
              if( choice_tunl .eq. 'rampcl' ) then
                  tunlramp = 0.5_r8*(1._r8+ctunl)*tunl
              elseif( choice_tunl .eq. 'rampsl' ) then
                  tunlramp = ctunl*tunl
                ! tunlramp = 0.765_r8
              else
                  tunlramp = tunl
              endif
              if( choice_leng .eq. 'origin' ) then
                  lz = ( (vk*zi(k,i))**(-cleng) + (tunlramp*lbulk)**(-cleng) )**(-1._r8/cleng)
                ! lz = vk*zi(k,i) / (1._r8+vk*zi(k,i)/(tunlramp*lbulk))
              else
                  lz = min( vk*zi(k,i), tunlramp*lbulk )              
              endif
              lz = min(leng_max(k), lz)
              dzinc = z(k-1,i) - z(k,i)
              lint = lint + dzinc
              l2n2 = l2n2 + lz*lz*n2(k,i)*dzinc
              l2s2 = l2s2 + lz*lz*s2(k,i)*dzinc
           end do

           ricll = min(l2n2/max(l2s2,ntzero),ricrit)
           trma = alph3*alph4*ricll+2._r8*b1*(alph2-alph4*alph5*ricll)
           trmb = ricll*(alph3+alph4)+2._r8*b1*(-alph5*ricll+alph1)
           trmc = ricll
           det = max(trmb*trmb-4._r8*trma*trmc,0._r8)
           gh = (-trmb + sqrt(det))/2._r8/trma
         ! gh = min(max(gh,-0.28_r8),0.0233_r8)
           gh = min(max(gh,-3.5334_r8),0.0233_r8)
           sh = alph5 / (1._r8+alph3*gh)
           sm = (alph1 + alph2*gh)/(1._r8+alph3*gh)/(1._r8+alph4*gh)
           ! Even though the 'wint' after finishing merging was positive, it is 
           ! possible that re-calculated 'wint' here is negative.  In this case,
           ! correct 'wint' to be a small positive number
           wint = max( wint - sh*l2n2 + sm*l2s2, 0.01_r8 )

       end if

       ! ---------------------------------------------------------------------- !
       ! Calculate final output variables of each CL (either has merged or not) !
       ! ---------------------------------------------------------------------- !

       lbrk(ncv,i) = lint
       wbrk(ncv,i) = wint/lint
       ebrk(ncv,i) = b1*wbrk(ncv,i)
       ebrk(ncv,i) = min(ebrk(ncv,i),tkemax)
       ricl(ncv,i) = ricll 
       ghcl(ncv,i) = gh 
       shcl(ncv,i) = sh
       smcl(ncv,i) = sm

       ! Increment counter for next CL. I should check if the increament of 'ncv'
       ! below is reasonable or not, since whenever CL is merged during downward
       ! extension process, 'ncv' is lowered down continuously within 'do' loop.
       ! But it seems that below 'ncv = ncv + 1' is perfectly correct.

       ncv = ncv + 1

    end do                   ! End of loop over each CL regime, ncv.

    ! ---------------------------------------------------------- !
    ! Re-initialize external interface indices which are not CLs !
    ! ---------------------------------------------------------- !

    do ncv = ncvfin(i) + 1, ncvmax
       ktop(ncv,i)  = 0
       kbase(ncv,i) = 0
    end do

    ! ------------------------------------------------ !
    ! Update CL interface identifiers, 'belongcv'      !
    ! CL external interfaces are also identified as CL !
    ! ------------------------------------------------ !

    do k = 1, pver + 1
       belongcv(k,i) = .false.
    end do

    do ncv = 1, ncvfin(i)
       do k = ktop(ncv,i), kbase(ncv,i)
          belongcv(k,i) = .true.
       end do
    end do

    return

    end subroutine zisocl

    real(r8) function compute_cubic(a,b,c)
    ! ------------------------------------------------------------------------- !
    ! Solve canonical cubic : x^3 + a*x^2 + b*x + c = 0,  x = sqrt(e)/sqrt(<e>) !
    ! Set x = max(xmin,x) at the end                                            ! 
    ! ------------------------------------------------------------------------- !
    implicit none
    real(r8), intent(in)     :: a, b, c
    real(r8)  qq, rr, dd, theta, aa, bb, x1, x2, x3
    real(r8), parameter      :: xmin = 1.e-2_r8
    
    qq = (a**2-3._r8*b)/9._r8 
    rr = (2._r8*a**3 - 9._r8*a*b + 27._r8*c)/54._r8
    
    dd = rr**2 - qq**3
    if( dd .le. 0._r8 ) then
        theta = acos(rr/qq**(3._r8/2._r8))
        x1 = -2._r8*sqrt(qq)*cos(theta/3._r8) - a/3._r8
        x2 = -2._r8*sqrt(qq)*cos((theta+2._r8*3.141592_r8)/3._r8) - a/3._r8
        x3 = -2._r8*sqrt(qq)*cos((theta-2._r8*3.141592_r8)/3._r8) - a/3._r8
        compute_cubic = max(max(max(x1,x2),x3),xmin)        
        return
    else
        if( rr .ge. 0._r8 ) then
            aa = -(sqrt(rr**2-qq**3)+rr)**(1._r8/3._r8)
        else
            aa =  (sqrt(rr**2-qq**3)-rr)**(1._r8/3._r8)
        endif
        if( aa .eq. 0._r8 ) then
            bb = 0._r8
        else
            bb = qq/aa
        endif
        compute_cubic = max((aa+bb)-a/3._r8,xmin) 
        return
    endif

    return
    end function compute_cubic

end module grist_eddy_diff
