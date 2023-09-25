!===================================================================================
!
!  Created by LiXiaohan on 19/07/10, adopted from CAM5
!
!  Zhang-McFarlane deep convection scheme
! 
!===================================================================================

 module grist_zm_conv
    use grist_constants,                    only: r8, i4
    use grist_handle_error,                 only: endrun
    use grist_mpi

    implicit none
    private
    save

    public ::   read_nml_zmconv,    &
                zm_conv_init,       &
                zm_conv_tend,       &
                zm_conv_tend_2,     &
                zm_conv_end

    ! Private module data
    real(r8), allocatable, dimension(:,:) :: mu        !(pver,ncol)
    real(r8), allocatable, dimension(:,:) :: eu        !(pver,ncol)
    real(r8), allocatable, dimension(:,:) :: du        !(pver,ncol)
    real(r8), allocatable, dimension(:,:) :: md        !(pver,ncol)
    real(r8), allocatable, dimension(:,:) :: ed        !(pver,ncol)
    real(r8), allocatable, dimension(:,:) :: dp        !(pver,ncol) 
    real(r8), allocatable, dimension(:)   :: dsubcld   !(ncol)
    integer,  allocatable, dimension(:)   :: jt        !(ncol)
    integer,  allocatable, dimension(:)   :: maxg      !(ncol)
    integer,  allocatable, dimension(:)   :: ideep     !(ncol)               

    real(r8), parameter :: unset_r8 = huge(1.0_r8)
    real(r8), parameter :: capelmt  = 70._r8                ! threshold value for cape for deep convection
    real(r8), parameter :: tiedke_add = 0.5_r8
    real(r8), parameter :: c1       = 6.112_r8
    real(r8), parameter :: c2       = 17.67_r8
    real(r8), parameter :: c3       = 243.5_r8
 
    integer                            :: lengath
    integer                            :: ncol, pver, pverp, pcnst
    logical                            :: no_deep_pbl       ! =true eliminates deep convection entirely within PBL
    integer(i4)                        :: limcnv            ! top interface level limit for convection
    real(r8)                           :: ke                ! Tunable evaporation efficiency set from namelist input zmconv_ke, 0.5e-6 - 10e-6, default:1e-6
    real(r8)                           :: c0_lnd            ! set from namelist input zmconv_c0_lnd, 1.0e-3 - 0.01, default:0.0059
    real(r8)                           :: c0_ocn            ! set from namelist input zmconv_c0_ocn),1.0e-3 - 0.1, default:0.045
    real(r8)                           :: tau               ! convective time scale set from namelist input zmconv_tau, 1800-28800, default:3600
    real(r8)                           :: tfreez, eps1, rl,     &
                                          cpres, grav, cp,      &
                                          rgrav, rgas, rh2o,    &
                                          cpliq, cpwv, latice


    contains

    subroutine read_nml_zmconv(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer             :: unitn, ierr
    real(r8)            :: zmconv_c0_lnd = unset_r8    
    real(r8)            :: zmconv_c0_ocn = unset_r8    
    real(r8)            :: zmconv_ke     = unset_r8    
    real(r8)            :: zmconv_tau    = unset_r8    

    namelist /zmconv_nl/ zmconv_c0_lnd, zmconv_c0_ocn, zmconv_ke, zmconv_tau

    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, zmconv_nl, iostat=ierr)
    if (ierr /= 0) call endrun(" error reading zmconv namelist")
    close(unitn)
 
    c0_lnd = zmconv_c0_lnd
    c0_ocn = zmconv_c0_ocn
    ke     = zmconv_ke
    tau    = zmconv_tau 
    end subroutine read_nml_zmconv


    subroutine zm_conv_init(ncol_in   , pver_in   , pverp_in  , pcnst_in  ,     & 
                            pref_edge , tfreez_in , eps_in    , rl_in     ,     & 
                            latice_in , cp_in     , cpliq_in  , cpwv_in   ,     &
                            grav_in   , rair_in   , rh2o_in   , no_deep_pbl_in)
! io
    integer(i4), intent(in)     :: ncol_in
    integer(i4), intent(in)     :: pcnst_in
    integer(i4), intent(in)     :: pver_in, pverp_in
    logical,     intent(in)     :: no_deep_pbl_in
    real(r8),    intent(in)     :: pref_edge(pverp_in)      ! reference pressures at interfaces
    real(r8),    intent(in)     :: tfreez_in, eps_in,  rl_in, latice_in,    & 
                                   cp_in, cpliq_in, cpwv_in , grav_in,      &
                                   rair_in, rh2o_in
! local
    integer                     :: k
    
    ncol  = ncol_in
    pcnst = pcnst_in
    pver  = pver_in
    pverp = pverp_in

    ! Initialization of ZM constants
    tfreez  = tfreez_in
    eps1    = eps_in
    rl      = rl_in
    latice  = latice_in
    cpres   = cp_in
    rgrav   = 1._r8/grav_in
    rgas    = rair_in
    rh2o    = rh2o_in
    grav    = grav_in
    cp      = cpres
    cpliq   = cpliq_in
    cpwv    = cpwv_in

    no_deep_pbl = no_deep_pbl_in

    allocate( mu(pver, ncol))
    allocate( eu(pver, ncol))
    allocate( du(pver, ncol))
    allocate( md(pver, ncol))
    allocate( ed(pver, ncol))
    allocate( dp(pver, ncol))
    allocate( dsubcld(ncol))
    allocate( jt(ncol))
    allocate( maxg(ncol))
    allocate( ideep(ncol))

    limcnv = 0   ! null value to check against below
    if (pref_edge(1) >= 4.e3_r8) then
       limcnv = 1
    else
       do k=1,pver
          if (pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8) then
             limcnv = k
             exit
          end if
       end do
       if ( limcnv == 0 ) limcnv = pverp
    end if
 
    if(mpi_rank()==0) print*, 'zm_conv_init: Deep convection will be capped at ',limcnv, &
                              'level, which is ', pref_edge(limcnv),' Pa'
    end subroutine zm_conv_init

 
    subroutine zm_conv_end

    if(allocated( mu))      deallocate( mu)
    if(allocated( eu))      deallocate( eu)
    if(allocated( du))      deallocate( du)
    if(allocated( md))      deallocate( md)
    if(allocated( ed))      deallocate( ed)
    if(allocated( dp))      deallocate( dp)
    if(allocated( dsubcld)) deallocate( dsubcld)
    if(allocated( jt))      deallocate( jt)
    if(allocated( maxg))    deallocate( maxg)
    if(allocated( ideep))   deallocate( ideep)
       
    end subroutine zm_conv_end


    subroutine zm_conv_tend( pblh    , mcon    , cme     , tpert   , dlf    , &
                             pflx    , zdu     , landfrac, rliq    , dt     , &
                             t       , q       , u       , v       , ps     , &
                             s       , pmid    , pint    , pdel    ,          &
                             zm      , zi      , phis    ,                    &
                             jctop   , jcbot   , ql      , fracis  , cld    , &
                             prec    , snow    , rprd    , evapcdp ,          &
                             ptend_u , ptend_v , ptend_q , ptend_s ,          &
                             flxprec , flxsnow)

    use grist_physics_update,           only: geopotential_dse
    use grist_physics_data_structure,   only: phy_tracer_info

! io
    real(r8),    intent(in)    :: dt                        ! model time increment
    real(r8),    intent(in)    :: pblh(ncol)                ! Planetary boundary layer height 
    real(r8),    intent(in)    :: tpert(ncol)               ! Thermal temperature excess 
    real(r8),    intent(in)    :: landfrac(ncol)            ! RBN - Landfrac  

    real(r8),    intent(in)    :: ps(ncol)                  ! surface pressure
    real(r8),    intent(in)    :: phis(ncol)                ! surface geopotential height
    real(r8),    intent(in)    :: pdel(pver, ncol)          ! delta pressure of model level
    real(r8),    intent(in)    :: pmid(pver, ncol)          ! pressure at model level
    real(r8),    intent(in)    :: pint(pverp, ncol)         ! pressure at interface
    real(r8),    intent(in)    :: fracis(pcnst, pver, ncol) ! fraction of tracer that is insoluble
    real(r8),    intent(in)    :: cld(pver, ncol)           ! cloud fraction
    real(r8),    intent(inout) :: zm(pver, ncol)            ! model level height
    real(r8),    intent(inout) :: zi(pverp, ncol)           ! interface height
    real(r8),    intent(inout) :: t(pver, ncol)             ! temperature at model level
    real(r8),    intent(inout) :: s(pver, ncol)             ! temperature at model level
    real(r8),    intent(inout) :: u(pver, ncol)             ! u wind speed at model level
    real(r8),    intent(inout) :: v(pver, ncol)             ! v wind speed at model level
    real(r8),    intent(inout) :: q(pcnst, pver, ncol)      ! tracer at model level

    real(r8),    intent(out)   :: mcon(pverp, ncol)         ! Convective mass flux--m sub c
    real(r8),    intent(out)   :: dlf(pver, ncol)           ! scattrd version of the detraining cld h2o tend
    real(r8),    intent(out)   :: pflx(pverp, ncol)         ! scattered precip flux at each level
    real(r8),    intent(out)   :: cme(pver, ncol)           ! cmf condensation - evaporation
    real(r8),    intent(out)   :: zdu(pver, ncol)           ! detraining mass flux
    real(r8),    intent(out)   :: ql(pver,  ncol)           ! wg grid slice of cloud liquid water
    real(r8),    intent(out)   :: rprd(pver,  ncol)         ! rain production rate
    real(r8),    intent(out)   :: evapcdp(pver,  ncol)      ! rain production rate
    real(r8),    intent(out)   :: rliq(ncol)                ! reserved liquid (not yet in cldliq) for energy integrals
    real(r8),    intent(out)   :: prec(ncol)                ! Convective-scale preciptn rate
    real(r8),    intent(out)   :: snow(ncol)                ! Convective-scale snowfall rate
    real(r8),    intent(out)   :: jctop(ncol)               ! row of top-of-deep-convection indices passed out
    real(r8),    intent(out)   :: jcbot(ncol)               ! row of base-of-cloud indices passed out

    real(r8),    intent(out)   :: ptend_u(pver, ncol)
    real(r8),    intent(out)   :: ptend_v(pver, ncol)
    real(r8),    intent(out)   :: ptend_s(pver, ncol)
    real(r8),    intent(out)   :: ptend_q(pcnst, pver, ncol)

    real(r8),    intent(out)   :: flxprec(pverp, ncol)
    real(r8),    intent(out)   :: flxsnow(pverp, ncol)

! local
    integer(i4) :: i,k,m,ii
    real(r8)    :: ptend_loc_u(pver, ncol)
    real(r8)    :: ptend_loc_v(pver, ncol)
    real(r8)    :: ptend_loc_s(pver, ncol)
    real(r8)    :: ptend_loc_q(pcnst, pver, ncol)
    real(r8)    :: ftem(pver, ncol)
    real(r8)    :: ntprprd(pver, ncol)                      ! evap outfld: net precip production in layer
    real(r8)    :: ntsnprd(pver, ncol)                      ! evap outfld: net snow production in layer
    real(r8)    :: tend_s_snwprd(pver, ncol)                ! Heating rate of snow production
    real(r8)    :: tend_s_snwevmlt(pver, ncol)              ! Heating rate of evap/melting of snow
    real(r8)    :: fake_dpdry(pver, ncol)                   ! used in convtran call
    real(r8)    :: pcont(ncol), pconb(ncol), freqzm(ncol)
    real(r8)    :: cape(ncol)
    real(r8)    :: mu_out(pver, ncol), md_out(pver, ncol)
    real(r8)    :: winds(2, pver, ncol)
    real(r8)    :: wind_tends(2, pver, ncol)
    real(r8)    :: pguall(2, pver, ncol), pgdall(2, pver, ncol)
    real(r8)    :: icwu(2, pver, ncol), icwd(2, pver, ncol)
    real(r8)    :: seten(pver, ncol)
    logical     :: l_windt(2)
    logical     :: lq(pcnst)

! initialize ptend and ptend_loc
    ptend_u     = 0._r8
    ptend_v     = 0._r8
    ptend_s     = 0._r8
    ptend_q     = 0._r8
    ptend_loc_s = 0._r8
    ptend_loc_q = 0._r8

    ftem(:,:)   = 0._r8
    mu_out(:,:) = 0._r8
    md_out(:,:) = 0._r8
    wind_tends(:,:pver,:ncol) = 0._r8

    call zm_convr( t      , q(1,:,:) , prec   , jctop  , jcbot   ,              &
                   pblh   , zm       , phis   , zi     , ptend_loc_q(1,:,:) ,   &
                   ptend_loc_s       , pmid   , pint   , pdel    , 0.5_r8*dt,   &
                   mcon   , cme      , cape   , tpert  , dlf     , pflx     ,   &
                   zdu    , rprd     , ql     , rliq   , landfrac )

! fractional occurance of ZM convection
    freqzm(:) = 0._r8
    do i = 1,lengath
       freqzm(ideep(i)) = 1.0_r8
    end do
 
! output: 'FREQZM'freqzm

! Convert mass flux from reported mb/s to kg/m^2/s
    mcon(:pver,:ncol) = mcon(:pver,:ncol) * 100._r8/grav

! Store upward and downward mass fluxes in un-gathered arrays
! + convert from mb/s to kg/m^2/s
    do i=1,lengath
       do k=1,pver
          ii = ideep(i)
          mu_out(k,ii) = mu(k,i) * 100._r8/grav
          md_out(k,ii) = md(k,i) * 100._r8/grav
       end do
    end do
 
! output: 'ZMMU'mu_out 'ZMMD'md_out
! output: 'ZMDT'ftem  'ZMDQ'ptend_loc_q(1,:,:) 
    ftem(:pver,:ncol) = ptend_loc_s(:pver,:ncol)/cp

! top and bottom pressure of cumulus  
    pcont(:ncol) = ps(:ncol)
    pconb(:ncol) = ps(:ncol)
    do i = 1,lengath
        if (maxg(i).gt.jt(i)) then
           pcont(ideep(i)) = pmid(jt(i),ideep(i))  ! gathered array (or jctop ungathered)
           pconb(ideep(i)) = pmid(maxg(i),ideep(i))! gathered array
        endif
        !     write(iulog,*) ' pcont, pconb ', pcont(i), pconb(i), cnt(i), cnb(i)
    end do
 
! output: 'PCONVT'pcont 'PCONVB'pconb

! update ptend and state
    ptend_s     = ptend_s + ptend_loc_s
    s(:,1:ncol) = s(:,1:ncol) + ptend_loc_s(:,1:ncol)*dt
    ptend_q(1,:,1:ncol) = ptend_q(1,:,1:ncol) + ptend_loc_q(1,:,1:ncol)
    q(1,:,1:ncol)       = q(1,:,1:ncol) + ptend_loc_q(1,:,1:ncol)*dt

! check negetive qv
    call qneg3('deep convection qv1', ncol, pver, 1, 1, phy_tracer_info(1)%qmin, q(1,:,1:ncol))
    !where(q(1,:,1:ncol) .lt. phy_tracer_info(1)%qmin) q(1,:,1:ncol) = phy_tracer_info(1)%qmin

    call geopotential_dse(ncol, pint, pmid, pdel, phis, s, q(1,:,1:ncol),      &
                          t, zm, zi)

! initialize ptend_loc
    ptend_loc_s = 0._r8
    ptend_loc_q = 0._r8

! Determine the phase of the precipitation produced and add latent heat of fusion
! Evaporate some of the precip directly into the environment (Sundqvist)
! Allow this to use the updated state1 and the fresh ptend_loc type
! heating and specific humidity tendencies produced
    call zm_conv_evap(t       , pmid   , pdel   , q(1,:,1:ncol)    ,            &
                      ptend_loc_s      , tend_s_snwprd   , tend_s_snwevmlt  ,   &
                      ptend_loc_q(1,:,1:ncol)   , rprd   , cld     , dt     ,   &
                      prec    , snow   , ntprprd, ntsnprd, flxprec , flxsnow)

    evapcdp(:pver,:ncol) = ptend_loc_q(1,:pver,:ncol)

    ftem(:pver,:ncol) = ptend_loc_s(:pver,:ncol)/cp
! output: 'EVAPTZM'ftem

    ftem(:pver,:ncol) = tend_s_snwprd(:pver,:ncol)/cp
! output: 'FZSNTZM'ftem

    ftem(:pver,:ncol) = tend_s_snwevmlt(:pver,:ncol)/cp
! output: 'EVSNTZM'ftem
! output: 'EVAPQZM'ptend_loc_q(1,:,:)  'ZMFLXPRC'flxprec  'ZMFLXSNW'flxsnow  'ZMNTPRPD'ntprprd 
!         'ZMNTSNPD'ntsnprd  'ZMEIHEAT'ptend_loc_s  'CMFMCDZM'mcon  'PRECCDZM'prec
         
! update ptend and state
    ptend_s     = ptend_s + ptend_loc_s
    s(:,1:ncol) = s(:,1:ncol) + ptend_loc_s(:,1:ncol)*dt
    ptend_q(1,:,1:ncol) = ptend_q(1,:,1:ncol) + ptend_loc_q(1,:,1:ncol)
    q(1,:,1:ncol)       = q(1,:,1:ncol) + ptend_loc_q(1,:,1:ncol)*dt

! check negetive qv and correct to zero, LiXH
    call qneg3('deep convection qv2', ncol, pver, 1, 1, phy_tracer_info(1)%qmin, q(1,:,1:ncol))
    !where(q(1,:,1:ncol) .lt. phy_tracer_info(1)%qmin) q(1,:,1:ncol) = phy_tracer_info(1)%qmin

    call geopotential_dse(ncol, pint, pmid, pdel, phis, s, q(1,:,1:ncol),      &
                          t, zm, zi)

! initialize ptend_loc
    ptend_loc_u = 0._r8
    ptend_loc_v = 0._r8
    ptend_loc_s = 0._r8
 
    winds(1,:pver,:ncol) = u(:pver,:ncol)
    winds(2,:pver,:ncol) = v(:pver,:ncol)
  
    l_windt(1) = .true.
    l_windt(2) = .true.

    call momtran (l_windt , winds   , 2         , maxg  ,         & 
                  1       , lengath , wind_tends, pguall, pgdall, &
                  icwu    , icwd    , dt        , seten)  

    ptend_loc_u(:pver,:ncol) = wind_tends(1,:pver,:ncol)
    ptend_loc_v(:pver,:ncol) = wind_tends(2,:pver,:ncol)
    ptend_loc_s(:pver,:ncol) = seten(:pver,:ncol)  

! update ptend and state
    ptend_s     = ptend_s + ptend_loc_s
    ptend_u     = ptend_u + ptend_loc_u
    ptend_v     = ptend_v + ptend_loc_v
    s(:,1:ncol) = s(:,1:ncol) + ptend_loc_s(:,1:ncol)*dt
    u(:,1:ncol) = u(:,1:ncol) + ptend_loc_u(:,1:ncol)*dt
    v(:,1:ncol) = v(:,1:ncol) + ptend_loc_v(:,1:ncol)*dt

    call geopotential_dse(ncol, pint, pmid, pdel, phis, s, q(1,:,1:ncol),      &
                          t, zm, zi)
    ftem(:pver,:ncol) = seten(:pver,:ncol)/cp
! output 'ZMMTT'ftem  'ZMMTU'wind_tends(1,:,:)  'ZMMTV'wind_tends(2,:,:) 
! output apparent force from  pressure gradient
!        'ZMUPGU'pguall(1,:,:)  'ZMUPGD'pgdall(1,:,:) 
!        'ZMVPGU'pguall(2,:,:)  'ZMVPGD'pgdall(2,:,:)
! output in-cloud winds
!        'ZMICUU'icwu(1,:,:)    'ZMICUD'icwd(1,:,:)
!        'ZMICVU'icwu(2,:,:)    'ZMICVD'icwd(2,:,:)    

! initialize ptend_loc
    ptend_loc_q = 0._r8

! dpdry is not used in this call to convtran since the cloud liquid and ice mixing
! ratios are moist
    fake_dpdry(:,:) = 0._r8

    lq(:) = .false.
    do m = 1, pcnst
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')lq(m) = .true.
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')   lq(m) = .true.
    end do

    call convtran (lq     , q     , pcnst      , maxg  , 1 ,   &
                   lengath, fracis, ptend_loc_q, fake_dpdry)

! output: 'ZMDICE'ptend_loc_q('cloud_liquid',:,:) 'ZMDLIQ'ptend_loc_q('cloud_ice',:,:)

! update ptend and state
    ptend_q = ptend_q + ptend_loc_q

    end subroutine zm_conv_tend


    subroutine zm_conv_tend_2(q, fracis, pdeldry, ptend_q)
    use grist_physics_data_structure,   only: phy_tracer_info
! io
    real(r8),    intent(in)    :: pdeldry(pver, ncol)
    real(r8),    intent(in)    :: fracis(pcnst, pver, ncol) ! fraction of tracer that is insoluble
    real(r8),    intent(in)    :: q(pcnst, pver, ncol)      ! tracer at model level
    real(r8),    intent(out)   :: ptend_q(pcnst, pver, ncol)

! local
    integer     :: m, i
    real(r8)    :: dpdry(pver, ncol)
    logical     :: lq(pcnst)

! initialize ptend
    ptend_q = 0._r8

    lq(:) = .false.
    do m = 2, pcnst
        if(phy_tracer_info(m)%longname .ne. 'cloud_liquid' .and. phy_tracer_info(m)%longname .ne. 'cloud_ice') lq(m) = .true.
    end do

    dpdry = 0._r8
    if(any(lq(:)))then
        do i = 1, lengath
            dpdry(:,i) = pdeldry(:,ideep(i))/100._r8
        end do

        call convtran (lq     , q     , pcnst    , maxg  , 1 ,   &
                       lengath, fracis, ptend_q  , dpdry)

     end if
    end subroutine zm_conv_tend_2


!----------------------------------------------------------------------- 
! Purpose: 
! Main driver for zhang-mcfarlane convection scheme 
! 
! Method: 
! performs deep convective adjustment based on mass-flux closure
! algorithm.
! 
! Author:guang jun zhang, m.lazare, n.mcfarlane. CAM Contact: P. Rasch
!
! This is contributed code not fully standardized by the CAM core group.
! All variables have been typed, where most are identified in comments
! The current procedure will be reimplemented in a subsequent version
! of the CAM where it will include a more straightforward formulation
! and will make use of the standard CAM nomenclature
!-----------------------------------------------------------------------

!  index of variables :
!
!  wg * alpha    array of vertical differencing used (=1. for upstream).
!  w  * cape     convective available potential energy.
!  wg * capeg    gathered convective available potential energy.
!  c  * capelmt  threshold value for cape for deep convection.
!  ic  * cpres    specific heat at constant pressure in j/kg-degk.
!  i  * dpp      
!  ic  * delt     length of model time-step in seconds.
!  wg * dp       layer thickness in mbs (between upper/lower interface).
!  wg * dqdt     mixing ratio tendency at gathered points.
!  wg * dsdt     dry static energy ("temp") tendency at gathered points.
!  wg * dudt     u-wind tendency at gathered points.
!  wg * dvdt     v-wind tendency at gathered points.
!  wg * dsubcld  layer thickness in mbs between lcl and maxi.
!  ic  * grav     acceleration due to gravy in m/sec2.
!  wg * du       detrainment in updraft. specified in mid-layer
!  wg * ed       entrainment in downdraft.
!  wg * eu       entrainment in updraft.
!  wg * hmn      moist static energy.
!  wg * hsat     saturated moist static energy.
!  w  * ideep    holds position of gathered points vs longitude index.
!  ic  * pver     number of model levels.
!  wg * j0       detrainment initiation level index.
!  wg * jd       downdraft   initiation level index.
!  ic  * jlatpr   gaussian latitude index for printing grids (if needed).
!  wg * jt       top  level index of deep cumulus convection.
!  w  * lcl      base level index of deep cumulus convection.
!  wg * lclg     gathered values of lcl.
!  w  * lel      index of highest theoretical convective plume.
!  wg * lelg     gathered values of lel.
!  w  * lon      index of onset level for deep convection.
!  w  * maxi     index of level with largest moist static energy.
!  wg * maxg     gathered values of maxi.
!  wg * mb       cloud base mass flux.
!  wg * mc       net upward (scaled by mb) cloud mass flux.
!  wg * md       downward cloud mass flux (positive up).
!  wg * mu       upward   cloud mass flux (positive up). specified
!                at interface
!  ic  * msg      number of missing moisture levels at the top of model.
!  w  * p        grid slice of ambient mid-layer pressure in mbs.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  i  * pblt     row of pbl top indices.
!  w  * pcpdh    scaled surface pressure.
!  w  * pf       grid slice of ambient interface pressure in mbs.
!  wg * pg       grid slice of gathered values of p.
!  w  * q        grid slice of mixing ratio.
!  wg * qd       grid slice of mixing ratio in downdraft.
!  wg * qg       grid slice of gathered values of q.
!  i/o * qh       grid slice of specific humidity.
!  w  * qh0      grid slice of initial specific humidity.
!  wg * qhat     grid slice of upper interface mixing ratio.
!  wg * ql       grid slice of cloud liquid water.
!  wg * qs       grid slice of saturation mixing ratio.
!  w  * qstp     grid slice of parcel temp. saturation mixing ratio.
!  wg * qstpg    grid slice of gathered values of qstp.
!  wg * qu       grid slice of mixing ratio in updraft.
!  ic  * rgas     dry air gas constant.
!  wg * rl       latent heat of vaporization.
!  w  * s        grid slice of scaled dry static energy (t+gz/cp).
!  wg * sd       grid slice of dry static energy in downdraft.
!  wg * sg       grid slice of gathered values of s.
!  wg * shat     grid slice of upper interface dry static energy.
!  wg * su       grid slice of dry static energy in updraft.
!  i/o * t       
!  o  * jctop    row of top-of-deep-convection indices passed out.
!  O  * jcbot    row of base of cloud indices passed out.
!  wg * tg       grid slice of gathered values of t.
!  w  * tl       row of parcel temperature at lcl.
!  wg * tlg      grid slice of gathered values of tl.
!  w  * tp       grid slice of parcel temperatures.
!  wg * tpg      grid slice of gathered values of tp.
!  i/o * u        grid slice of u-wind (real).
!  wg * ug       grid slice of gathered values of u.
!  i/o * utg      grid slice of u-wind tendency (real).
!  i/o * v        grid slice of v-wind (real).
!  w  * va       work array re-used by called subroutines.
!  wg * vg       grid slice of gathered values of v.
!  i/o * vtg      grid slice of v-wind tendency (real).
!  i  * w        grid slice of diagnosed large-scale vertical velocity.
!  w  * z        grid slice of ambient mid-layer height in metres.
!  w  * zf       grid slice of ambient interface height in metres.
!  wg * zfg      grid slice of gathered values of zf.
!  wg * zg       grid slice of gathered values of z.
!
!-----------------------------------------------------------------------
!
! multi-level i/o fields:
!  i      => input arrays.
!  i/o    => input/output arrays.
!  w      => work arrays.
!  wg     => work arrays operating only on gathered points.
!  ic     => input data constants.
!  c      => data constants pertaining to subroutine itself.
!-----------------------------------------------------------------------

    subroutine zm_convr( t       ,qh      ,prec    ,jctop   ,jcbot   , &
                         pblh    ,zm      ,geos    ,zi      ,qtnd    , &
                         heat    ,pap     ,paph    ,dpp     ,delt    , &
                         mcon    ,cme     ,cape    ,tpert   ,dlf     , &
                         pflx    ,zdu     ,rprd    ,ql      ,rliq    , &
                         landfrac)
! io
    real(r8), intent(in) :: delt                    ! length of model time-step in seconds.
    real(r8), intent(in) :: t(pver, ncol)           ! grid slice of temperature at mid-layer.
    real(r8), intent(in) :: qh(pver, ncol)          ! grid slice of specific humidity.
    real(r8), intent(in) :: pap(pver, ncol)     
    real(r8), intent(in) :: paph(pverp, ncol)
    real(r8), intent(in) :: dpp(pver, ncol)         ! local sigma half-level thickness (i.e. dshj).
    real(r8), intent(in) :: zm(pver, ncol)
    real(r8), intent(in) :: geos(ncol)
    real(r8), intent(in) :: zi(pverp, ncol)
    real(r8), intent(in) :: pblh(ncol)
    real(r8), intent(in) :: tpert(ncol)
    real(r8), intent(in) :: landfrac(ncol)          ! RBN Landfrac

    real(r8), intent(out) :: qtnd(pver, ncol)       ! specific humidity tendency (kg/kg/s)
    real(r8), intent(out) :: heat(pver, ncol)       ! heating rate (dry static energy tendency, W/kg)
    real(r8), intent(out) :: mcon(pverp, ncol)
    real(r8), intent(out) :: dlf(pver, ncol)        ! scattrd version of the detraining cld h2o tend
    real(r8), intent(out) :: pflx(pverp, ncol)      ! scattered precip flux at each level
    real(r8), intent(out) :: cme(pver, ncol)
    real(r8), intent(out) :: cape(ncol)             ! w  convective available potential energy.
    real(r8), intent(out) :: zdu(pver, ncol)
    real(r8), intent(out) :: rprd(pver, ncol)       ! rain production rate
    real(r8), intent(out) :: ql(pver, ncol)         ! wg grid slice of cloud liquid water.
! move these vars from local storage to output so that convective
! transports can be done in outside of conv_cam.
    real(r8), intent(out) :: jctop(ncol)            ! o row of top-of-deep-convection indices passed out.
    real(r8), intent(out) :: jcbot(ncol)            ! o row of base of cloud indices passed out.
    real(r8), intent(out) :: prec(ncol)
    real(r8), intent(out) :: rliq(ncol)             ! reserved liquid (not yet in cldliq) for energy integrals

! local
    real(r8) zs(ncol)
    real(r8) dlg(pver, ncol)            ! gathrd version of the detraining cld h2o tend
    real(r8) pflxg(pverp, ncol)         ! gather precip flux at each level
    real(r8) cug(pver, ncol)            ! gathered condensation rate
    real(r8) evpg(pver, ncol)           ! gathered evap rate of rain in downdraft
    real(r8) mumax(ncol)
!   diagnostic field used by chem/wetdep codes
    real(r8) pblt(ncol)                 ! i row of pbl top indices.

!   general work fields (local variables):
    real(r8) q(pver, ncol)              ! w  grid slice of mixing ratio.
    real(r8) p(pver, ncol)              ! w  grid slice of ambient mid-layer pressure in mbs.
    real(r8) z(pver, ncol)              ! w  grid slice of ambient mid-layer height in metres.
    real(r8) s(pver, ncol)              ! w  grid slice of scaled dry static energy (t+gz/cp).
    real(r8) tp(pver, ncol)             ! w  grid slice of parcel temperatures.
    real(r8) zf(pverp, ncol)            ! w  grid slice of ambient interface height in metres.
    real(r8) pf(pverp, ncol)            ! w  grid slice of ambient interface pressure in mbs.
    real(r8) qstp(pver, ncol)           ! w  grid slice of parcel temp. saturation mixing ratio.
    real(r8) tl(ncol)                   ! w  row of parcel temperature at lcl.
    integer lcl(ncol)                   ! w  base level index of deep cumulus convection.
    integer lel(ncol)                   ! w  index of highest theoretical convective plume.
    integer lon(ncol)                   ! w  index of onset level for deep convection.
    integer maxi(ncol)                  ! w  index of level with largest moist static energy.
    integer index(ncol)
    real(r8) precip

! gathered work fields:
    real(r8) qg(pver, ncol)             ! wg grid slice of gathered values of q.
    real(r8) tg(pver, ncol)             ! w  grid slice of temperature at interface.
    real(r8) pg(pver, ncol)             ! wg grid slice of gathered values of p.
    real(r8) zg(pver, ncol)             ! wg grid slice of gathered values of z.
    real(r8) sg(pver, ncol)             ! wg grid slice of gathered values of s.
    real(r8) tpg(pver, ncol)            ! wg grid slice of gathered values of tp.
    real(r8) zfg(pverp, ncol)           ! wg grid slice of gathered values of zf.
    real(r8) qstpg(pver, ncol)          ! wg grid slice of gathered values of qstp.
    real(r8) ug(pver, ncol)             ! wg grid slice of gathered values of u.
    real(r8) vg(pver, ncol)             ! wg grid slice of gathered values of v.
    real(r8) cmeg(pver, ncol)
    real(r8) rprdg(pver, ncol)          ! wg gathered rain production rate
    real(r8) capeg(ncol)                ! wg gathered convective available potential energy.
    real(r8) tlg(ncol)                  ! wg grid slice of gathered values of tl.
    real(r8) landfracg(ncol)            ! wg grid slice of landfrac  
    integer lclg(ncol)                  ! wg gathered values of lcl.
    integer lelg(ncol)

! work fields arising from gathered calculations.
    real(r8) dqdt(pver, ncol)           ! wg mixing ratio tendency at gathered points.
    real(r8) dsdt(pver, ncol)           ! wg dry static energy ("temp") tendency at gathered points.
!      real(r8) alpha(pver, ncol)       ! array of vertical differencing used (=1. for upstream).
    real(r8) sd(pver, ncol)             ! wg grid slice of dry static energy in downdraft.
    real(r8) qd(pver, ncol)             ! wg grid slice of mixing ratio in downdraft.
    real(r8) mc(pver, ncol)             ! wg net upward (scaled by mb) cloud mass flux.
    real(r8) qhat(pver, ncol)           ! wg grid slice of upper interface mixing ratio.
    real(r8) qu(pver, ncol)             ! wg grid slice of mixing ratio in updraft.
    real(r8) su(pver, ncol)             ! wg grid slice of dry static energy in updraft.
    real(r8) qs(pver, ncol)             ! wg grid slice of saturation mixing ratio.
    real(r8) shat(pver, ncol)           ! wg grid slice of upper interface dry static energy.
    real(r8) hmn(pver, ncol)            ! wg moist static energy.
    real(r8) hsat(pver, ncol)           ! wg saturated moist static energy.
    real(r8) qlg(pver, ncol)
    real(r8) dudt(pver, ncol)           ! wg u-wind tendency at gathered points.
    real(r8) dvdt(pver, ncol)           ! wg v-wind tendency at gathered points.
!      real(r8) ud(pver, ncol)
!      real(r8) vd(pver, ncol)
    real(r8) mb(ncol)                   ! wg cloud base mass flux.
    integer jlcl(ncol)
    integer j0(ncol)                    ! wg detrainment initiation level index.
    integer jd(ncol)                    ! wg downdraft initiation level index.
    integer i
    integer ii
    integer k
    integer msg                         !  ic number of missing moisture levels at the top of model.
    real(r8) qdifr
    real(r8) sdifr

! Set internal variable "msg" (convection limit) to "limcnv-1"
    msg = limcnv - 1

! initialize necessary arrays.
! zero out variables not used in cam
    qtnd(:,:) = 0._r8
    heat(:,:) = 0._r8
    mcon(:,:) = 0._r8
    rliq(:ncol)   = 0._r8

! initialize convective tendencies
    prec(1:ncol) = 0._r8
    do k = 1,pver
       do i = 1,ncol
          dqdt(k,i)  = 0._r8
          dsdt(k,i)  = 0._r8
          dudt(k,i)  = 0._r8
          dvdt(k,i)  = 0._r8
          pflx(k,i)  = 0._r8
          pflxg(k,i) = 0._r8
          cme(k,i)   = 0._r8
          rprd(k,i)  = 0._r8
          zdu(k,i)   = 0._r8
          ql(k,i)    = 0._r8
          qlg(k,i)   = 0._r8
          dlf(k,i)   = 0._r8
          dlg(k,i)   = 0._r8
       end do
    end do
    do i = 1,ncol
       pflx(pverp,i) = 0
       pflxg(pverp,i) = 0
    end do
 
    do i = 1,ncol
       pblt(i) = pver
       dsubcld(i) = 0._r8
 
       jctop(i) = pver
       jcbot(i) = 1
    end do

! calculate local pressure (mbs) and height (m) for both interface
! and mid-layer locations.
    do i = 1,ncol
       zs(i) = geos(i)*rgrav
       pf(pver+1,i) = paph(pver+1,i)*0.01_r8
       zf(pver+1,i) = zi(pver+1,i) + zs(i)
    end do
    do k = 1,pver
       do i = 1,ncol
          p(k,i) = pap(k,i)*0.01_r8
          pf(k,i) = paph(k,i)*0.01_r8
          z(k,i) = zm(k,i) + zs(i)
          zf(k,i) = zi(k,i) + zs(i)
       end do
    end do
 
    do k = pver - 1,msg + 1,-1
       do i = 1,ncol
          if (abs(z(k,i)-zs(i)-pblh(i)) < (zf(k,i)-zf(k+1,i))*0.5_r8) pblt(i) = k
       end do
    end do

! store incoming specific humidity field for subsequent calculation
! of precipitation (through change in storage).
! define dry static energy (normalized by cp).
    do k = 1,pver
       do i = 1,ncol
          q(k,i) = qh(k,i)
          s(k,i) = t(k,i) + (grav/cpres)*z(k,i)
          tp(k,i)=0.0_r8
!-------------LiXH Test--------------
!          shat(k,i) = s(k,i)
!          qhat(k,i) = q(k,i)
!-------------LiXH Test--------------
       end do
    end do
 
    do i = 1,ncol
       capeg(i) = 0._r8
       lclg(i) = 1
       lelg(i) = pver
       maxg(i) = 1
       tlg(i) = 400._r8
       dsubcld(i) = 0._r8
    end do

    call buoyan_dilute(q       ,t       ,p       ,z       ,pf       , &
                       tp      ,qstp    ,tl      ,rl      ,cape     , &
                       pblt    ,lcl     ,lel     ,lon     ,maxi     , &
                       rgas    ,grav    ,cpres   ,msg     ,tpert    )
 
! determine whether grid points will undergo some deep convection
! (ideep=1) or not (ideep=0), based on values of cape,lcl,lel
! (require cape.gt. 0 and lel<lcl as minimum conditions).

    lengath = 0
    do i=1,ncol
       if (cape(i) > capelmt) then
          lengath = lengath + 1
          index(lengath) = i
       end if
    end do

    if (lengath.eq.0) return

    do ii=1,lengath
       i=index(ii)
       ideep(ii)=i
    end do

! obtain gathered arrays necessary for ensuing calculations.
    do k = 1,pver
       do i = 1,lengath
          dp(k,i) = 0.01_r8*dpp(k,ideep(i))
          qg(k,i) = q(k,ideep(i))
          tg(k,i) = t(k,ideep(i))
          pg(k,i) = p(k,ideep(i))
          zg(k,i) = z(k,ideep(i))
          sg(k,i) = s(k,ideep(i))
          tpg(k,i) = tp(k,ideep(i))
          zfg(k,i) = zf(k,ideep(i))
          qstpg(k,i) = qstp(k,ideep(i))
          ug(k,i) = 0._r8
          vg(k,i) = 0._r8

!-------------LiXH Test--------------
          shat(k,i) = s(k,ideep(i))
          qhat(k,i) = q(k,ideep(i))
!-------------LiXH Test--------------
 
       end do
    end do
 
    do i = 1,lengath
       zfg(pver+1,i) = zf(pver+1,ideep(i))
    end do
    do i = 1,lengath
       capeg(i) = cape(ideep(i))
       lclg(i) = lcl(ideep(i))
       lelg(i) = lel(ideep(i))
       maxg(i) = maxi(ideep(i))
       tlg(i) = tl(ideep(i))
       landfracg(i) = landfrac(ideep(i))
    end do

! calculate sub-cloud layer pressure "thickness" for use in
! closure and tendency routines.
    do k = msg + 1,pver
       do i = 1,lengath
          if (k >= maxg(i)) then
             dsubcld(i) = dsubcld(i) + dp(k,i)
          end if
       end do
    end do

! define array of factors (alpha) which defines interfacial
! values, as well as interfacial values for (q,s) used in
! subsequent routines.
    do k = msg + 2,pver
       do i = 1,lengath
!            alpha(k,i) = 0.5
          sdifr = 0._r8
          qdifr = 0._r8
          if (sg(k,i) > 0._r8 .or. sg(k-1,i) > 0._r8) &
             sdifr = abs((sg(k,i)-sg(k-1,i))/max(sg(k-1,i),sg(k,i)))
          if (qg(k,i) > 0._r8 .or. qg(k-1,i) > 0._r8) &
             qdifr = abs((qg(k,i)-qg(k-1,i))/max(qg(k-1,i),qg(k,i)))
          if (sdifr > 1.E-6_r8) then
             shat(k,i) = log(sg(k-1,i)/sg(k,i))*sg(k-1,i)*sg(k,i)/(sg(k-1,i)-sg(k,i))
          else
             shat(k,i) = 0.5_r8* (sg(k,i)+sg(k-1,i))
          end if
          if (qdifr > 1.E-6_r8) then
             qhat(k,i) = log(qg(k-1,i)/qg(k,i))*qg(k-1,i)*qg(k,i)/(qg(k-1,i)-qg(k,i))
          else
             qhat(k,i) = 0.5_r8* (qg(k,i)+qg(k-1,i))
          end if
       end do
    end do
 
! obtain cloud properties.
    call cldprp(qg      ,tg      ,ug      ,vg      ,pg      , &
                zg      ,sg      ,                            &
                sd      ,qd      ,mc      ,                   &
                qu      ,su      ,zfg     ,qs      ,hmn     , &
                hsat    ,shat    ,qlg     ,                   &
                cmeg    ,maxg    ,lelg    ,jlcl    ,          &
                maxg    ,j0      ,jd      ,rl      ,lengath , &
                rgas    ,grav    ,cpres   ,msg     ,          &
                pflxg   ,evpg    ,cug     ,rprdg   ,limcnv  ,landfracg)

! convert detrainment from units of "1/m" to "1/mb".
    do k = msg + 1,pver
       do i = 1,lengath
          du   (k,i) = du   (k,i)* (zfg(k,i)-zfg(k+1,i))/dp(k,i)
          eu   (k,i) = eu   (k,i)* (zfg(k,i)-zfg(k+1,i))/dp(k,i)
          ed   (k,i) = ed   (k,i)* (zfg(k,i)-zfg(k+1,i))/dp(k,i)
          cug  (k,i) = cug  (k,i)* (zfg(k,i)-zfg(k+1,i))/dp(k,i)
          cmeg (k,i) = cmeg (k,i)* (zfg(k,i)-zfg(k+1,i))/dp(k,i)
          rprdg(k,i) = rprdg(k,i)* (zfg(k,i)-zfg(k+1,i))/dp(k,i)
          evpg (k,i) = evpg (k,i)* (zfg(k,i)-zfg(k+1,i))/dp(k,i)
       end do
    end do
 
    call closure(qg      ,tg      ,pg      ,zg      ,sg      , &
                 tpg     ,qs      ,qu      ,su      ,mc      , &
                 qd      ,sd      ,                            &
                 qhat    ,shat    ,qstpg   ,zfg     ,          &
                 qlg     ,mb      ,capeg   ,tlg     ,          &
                 lclg    ,lelg    ,maxg    ,1       ,          &
                 lengath ,rgas    ,grav    ,cpres   ,rl      , &
                 msg     )

! limit cloud base mass flux to theoretical upper bound.
    do i=1,lengath
       mumax(i) = 0
    end do
    do k=msg + 2,pver
       do i=1,lengath
         mumax(i) = max(mumax(i), mu(k,i)/dp(k,i))
       end do
    end do
 
    do i=1,lengath
       if (mumax(i) > 0._r8) then
          mb(i) = min(mb(i),0.5_r8/(delt*mumax(i)))
       else
          mb(i) = 0._r8
       endif
    end do
    ! If no_deep_pbl = .true., don't allow convection entirely 
    ! within PBL (suggestion of Bjorn Stevens, 8-2000)
 
    if (no_deep_pbl) then
       do i=1,lengath
          if (zm(jt(i), ideep(i)) < pblh(ideep(i))) mb(i) = 0
       end do
    end if
 
    do k=msg+1,pver
       do i=1,lengath
          mu   (k,i)  = mu   (k,i)*mb(i)
          md   (k,i)  = md   (k,i)*mb(i)
          mc   (k,i)  = mc   (k,i)*mb(i)
          du   (k,i)  = du   (k,i)*mb(i)
          eu   (k,i)  = eu   (k,i)*mb(i)
          ed   (k,i)  = ed   (k,i)*mb(i)
          cmeg (k,i)  = cmeg (k,i)*mb(i)
          rprdg(k,i)  = rprdg(k,i)*mb(i)
          cug  (k,i)  = cug  (k,i)*mb(i)
          evpg (k,i)  = evpg (k,i)*mb(i)
          pflxg(k+1,i)= pflxg(k+1,i)*mb(i)*100._r8/grav
       end do
    end do
 
    ! compute temperature and moisture changes due to convection.
     call q1q2_pjr(dqdt    ,dsdt    ,qg      ,qs      ,qu      , &
                   su      ,qhat    ,shat    ,                   &
                   sd      ,qd      ,qlg     ,                   &
                   maxg    ,1       ,lengath ,                   &
                   cpres   ,rl      ,msg     ,                   &
                   dlg     ,evpg    ,cug     )
 

    ! gather back temperature and mixing ratio.
    do k = msg + 1,pver
       do i = 1,lengath
! q is updated to compute net precip.
          q(k,ideep(i)) = qh(k,ideep(i)) + 2._r8*delt*dqdt(k,i)
          qtnd(k,ideep(i)) = dqdt (k,i)
          cme (k,ideep(i)) = cmeg (k,i)
          rprd(k,ideep(i)) = rprdg(k,i)
          zdu (k,ideep(i)) = du   (k,i)
          mcon(k,ideep(i)) = mc   (k,i)
          heat(k,ideep(i)) = dsdt (k,i)*cpres
          dlf (k,ideep(i)) = dlg  (k,i)
          pflx(k,ideep(i)) = pflxg(k,i)
          ql  (k,ideep(i)) = qlg  (k,i)
       end do
    end do
 
    do i = 1,lengath
       jctop(ideep(i)) = jt(i)
!++bee
       jcbot(ideep(i)) = maxg(i)
!--bee
       pflx(pverp,ideep(i)) = pflxg(pverp,i)
    end do
 
! Compute precip by integrating change in water vapor minus detrained cloud water
    do k = pver,msg + 1,-1
       do i = 1,ncol
          prec(i) = prec(i) - dpp(k,i)* (q(k,i)-qh(k,i)) - dpp(k,i)*dlf(k,i)*2*delt
       end do
    end do
 
! obtain final precipitation rate in m/s.
    do i = 1,ncol
       prec(i) = rgrav*max(prec(i),0._r8)/ (2._r8*delt)/1000._r8
    end do
 
! Compute reserved liquid (not yet in cldliq) for energy integrals.
! Treat rliq as flux out bottom, to be added back later.
     do k = 1, pver
        do i = 1, ncol
           rliq(i) = rliq(i) + dlf(k,i)*dpp(k,i)/grav
        end do
     end do
     rliq(:ncol) = rliq(:ncol) /1000._r8
 
    return
 
    end subroutine zm_convr


! Purpose:  Compute tendencies due to evaporation of rain from ZM scheme
!           Compute the total precipitation and snow fluxes at the surface.
!           Add in the latent heat of fusion for snow formation and melt, since it not dealt with
!           in the Zhang-MacFarlane parameterization.
!           Evaporate some of the precip directly into the environment using a Sundqvist type algorithm
    subroutine zm_conv_evap(t      , pmid   , pdel   , q      , tend_s ,            &
                            tend_s_snwprd   , tend_s_snwevmlt , tend_q , prdprec  , &
                            cldfrc , deltat , prec   , snow   , ntprprd, ntsnprd  , &
                            flxprec, flxsnow)

    use grist_wv_saturation,  only: qsat
    use cloud_fraction,       only: cldfrc_fice
! io
    real(r8), intent(in)     :: deltat                           ! time step
    real(r8), intent(in),    dimension(pver, ncol) :: t          ! temperature (K)
    real(r8), intent(in),    dimension(pver, ncol) :: pmid       ! midpoint pressure (Pa) 
    real(r8), intent(in),    dimension(pver, ncol) :: pdel       ! layer thickness (Pa)
    real(r8), intent(in),    dimension(pver, ncol) :: q          ! water vapor (kg/kg)
    real(r8), intent(in),    dimension(pver, ncol) :: prdprec    ! precipitation production (kg/ks/s)
    real(r8), intent(in),    dimension(pver, ncol) :: cldfrc     ! cloud fraction
    real(r8), intent(inout), dimension(pver, ncol) :: tend_s     ! heating rate (J/kg/s)
    real(r8), intent(inout), dimension(pver, ncol) :: tend_q     ! water vapor tendency (kg/kg/s)
    real(r8), intent(inout), dimension(ncol)      :: prec        ! Convective-scale preciptn rate
    real(r8), intent(out),   dimension(ncol)      :: snow        ! Convective-scale snowfall rate
    real(r8), intent(out),   dimension(pver, ncol) :: tend_s_snwprd   ! Heating rate of snow production
    real(r8), intent(out),   dimension(pver, ncol) :: tend_s_snwevmlt ! Heating rate of evap/melting of snow
    real(r8), intent(out),   dimension(pverp, ncol):: flxprec    ! Convective-scale flux of precip at interfaces (kg/m2/s)
    real(r8), intent(out),   dimension(pverp, ncol):: flxsnow    ! Convective-scale flux of snow   at interfaces (kg/m2/s)
    real(r8), intent(out),   dimension(pver, ncol) :: ntprprd    ! net precip production in layer
    real(r8), intent(out),   dimension(pver, ncol) :: ntsnprd    ! net snow production in layer
! local    
    real(r8) :: es    (pver, ncol)     ! Saturation vapor pressure
    real(r8) :: fice   (pver, ncol)    ! ice fraction in precip production
    real(r8) :: fsnow_conv(pver, ncol) ! snow fraction in precip production
    real(r8) :: qs   (pver, ncol)      ! saturation specific humidity
    real(r8) :: work1                  ! temp variable (pjr)
    real(r8) :: work2                  ! temp variable (pjr)
    real(r8) :: evpvint(ncol)          ! vertical integral of evaporation
    real(r8) :: evpprec(ncol)          ! evaporation of precipitation (kg/kg/s)
    real(r8) :: evpsnow(ncol)          ! evaporation of snowfall (kg/kg/s)
    real(r8) :: snowmlt(ncol)          ! snow melt tendency in layer
    real(r8) :: flxsntm(ncol)          ! flux of snow into layer, after melting
    real(r8) :: evplimit               ! temp variable for evaporation limits
    real(r8) :: rlat(ncol)
    integer :: i,k                     ! longitude,level indices

! convert input precip to kg/m2/s
    prec(:ncol) = prec(:ncol)*1000._r8

! determine saturation vapor pressure
    call qsat(t(1:pver,1:ncol), pmid(1:pver,1:ncol), es(1:pver,1:ncol), qs(1:pver,1:ncol))

! determine ice fraction in rain production (use cloud water parameterization fraction at present)
    call cldfrc_fice(ncol, t, fice, fsnow_conv)

! zero the flux integrals on the top boundary
    flxprec(1,:ncol) = 0._r8
    flxsnow(1,:ncol) = 0._r8
    evpvint(:ncol)   = 0._r8

    do k = 1, pver
       do i = 1, ncol
! Melt snow falling into layer, if necessary. 
          if (t(k,i) > tfreez) then
             flxsntm(i) = 0._r8
             snowmlt(i) = flxsnow(k,i) * grav/ pdel(k,i)
          else
             flxsntm(i) = flxsnow(k,i)
             snowmlt(i) = 0._r8
          end if

! relative humidity depression must be > 0 for evaporation
          evplimit = max(1._r8 - q(k,i)/qs(k,i), 0._r8)

! total evaporation depends on flux in the top of the layer
! flux prec is the net production above layer minus evaporation into environmet
          evpprec(i) = ke * (1._r8 - cldfrc(k,i)) * evplimit * sqrt(flxprec(k,i))
!**********************************************************
!!          evpprec(i) = 0.    ! turn off evaporation for now
!**********************************************************

! Don't let evaporation supersaturate layer (approx). Layer may already be saturated.
! Currently does not include heating/cooling change to qs
          evplimit   = max(0._r8, (qs(k,i)-q(k,i)) / deltat)

! Don't evaporate more than is falling into the layer - do not evaporate rain formed
! in this layer but if precip production is negative, remove from the available precip
! Negative precip production occurs because of evaporation in downdrafts.
!!$          evplimit   = flxprec(k,i) * grav / pdel(k,i) + min(prdprec(k,i), 0.)
          evplimit   = min(evplimit, flxprec(k,i) * grav / pdel(k,i))

! Total evaporation cannot exceed input precipitation
          evplimit   = min(evplimit, (prec(i) - evpvint(i)) * grav / pdel(k,i))

          evpprec(i) = min(evplimit, evpprec(i))

! evaporation of snow depends on snow fraction of total precipitation in the top after melting
          if (flxprec(k,i) > 0._r8) then
!            evpsnow(i) = evpprec(i) * flxsntm(i) / flxprec(k,i)
!            prevent roundoff problems
             work1 = min(max(0._r8,flxsntm(i)/flxprec(k,i)),1._r8)
             evpsnow(i) = evpprec(i) * work1
          else
             evpsnow(i) = 0._r8
          end if

! vertically integrated evaporation
          evpvint(i) = evpvint(i) + evpprec(i) * pdel(k,i)/grav

! net precip production is production - evaporation
          ntprprd(k,i) = prdprec(k,i) - evpprec(i)
! net snow production is precip production * ice fraction - evaporation - melting
!pjrworks ntsnprd(k,i) = prdprec(k,i)*fice(k,i) - evpsnow(i) - snowmlt(i)
!pjrwrks2 ntsnprd(k,i) = prdprec(k,i)*fsnow_conv(k,i) - evpsnow(i) - snowmlt(i)
! the small amount added to flxprec in the work1 expression has been increased from 
! 1e-36 to 8.64e-11 (1e-5 mm/day).  This causes the temperature based partitioning
! scheme to be used for small flxprec amounts.  This is to address error growth problems.
          if (flxprec(k,i).gt.0._r8) then
             work1 = min(max(0._r8,flxsnow(k,i)/flxprec(k,i)),1._r8)
          else
             work1 = 0._r8
          endif

          work2 = max(fsnow_conv(k,i), work1)
          if (snowmlt(i).gt.0._r8) work2 = 0._r8
!         work2 = fsnow_conv(k,i)
          ntsnprd(k,i) = prdprec(k,i)*work2 - evpsnow(i) - snowmlt(i)
          tend_s_snwprd  (k,i) = prdprec(k,i)*work2*latice
          tend_s_snwevmlt(k,i) = - ( evpsnow(i) + snowmlt(i) )*latice

! precipitation fluxes
          flxprec(k+1,i) = flxprec(k,i) + ntprprd(k,i) * pdel(k,i)/grav
          flxsnow(k+1,i) = flxsnow(k,i) + ntsnprd(k,i) * pdel(k,i)/grav

! protect against rounding error
          flxprec(k+1,i) = max(flxprec(k+1,i), 0._r8)
          flxsnow(k+1,i) = max(flxsnow(k+1,i), 0._r8)
! more protection (pjr)
!         flxsnow(k+1,i) = min(flxsnow(k+1,i), flxprec(k+1,i))

! heating (cooling) and moistening due to evaporation 
! - latent heat of vaporization for precip production has already been accounted for
! - snow is contained in prec
          tend_s(k,i)   =-evpprec(i)*rl + ntsnprd(k,i)*latice
          tend_q(k,i) = evpprec(i)
       end do
    end do

! set output precipitation rates (m/s)
    prec(:ncol) = flxprec(pver+1,:ncol) / 1000._r8
    snow(:ncol) = flxsnow(pver+1,:ncol) / 1000._r8

!**********************************************************
!!$    tend_s(:ncol,:)   = 0.      ! turn heating off
!**********************************************************

    end subroutine zm_conv_evap



! Purpose: transport of trace species
!          Mixing ratios may be with respect to either dry or moist air
! Method:  <Describe the algorithm(s) used in the routine.> 
!          <Also include any applicable external references.>
! Author: P. Rasch
    subroutine convtran(doconvtran    , q      , ncnst   , mx      ,    &
                        il1g   , il2g , fracis , dqdt    , dpdry   )
   implicit none
! io
   integer,  intent(in) :: ncnst                    ! number of tracers to transport
   integer,  intent(in) :: il1g                     ! Gathered min lon indices over which to operate
   integer,  intent(in) :: il2g                     ! Gathered max lon indices over which to operate
   logical,  intent(in) :: doconvtran(ncnst)        ! flag for doing convective transport
   real(r8), intent(in) :: q(ncnst,pver,ncol)       ! Tracer array including moisture
   real(r8), intent(in) :: fracis(ncnst,pver,ncol)  ! fraction of tracer that is insoluble
   integer,  intent(in) :: mx(ncol)                 ! Index of cloud top for each column
   real(r8), intent(in) :: dpdry(pver,ncol)         ! Delta pressure between interfaces
   real(r8), intent(out):: dqdt(ncnst,pver,ncol)    ! Tracer tendency array
! local
   integer i                    ! Work index
   integer k                    ! Work index
   integer kbm                  ! Highest altitude index of cloud base
   integer kk                   ! Work index
   integer kkp1                 ! Work index
   integer km1                  ! Work index
   integer kp1                  ! Work index
   integer ktm                  ! Highest altitude index of cloud top
   integer m                    ! Work index
   real(r8) cabv                ! Mix ratio of constituent above
   real(r8) cbel                ! Mix ratio of constituent below
   real(r8) cdifr               ! Normalized diff between cabv and cbel
   real(r8) chat(pver,ncol)     ! Mix ratio in env at interfaces
   real(r8) cond(pver,ncol)     ! Mix ratio in downdraft at interfaces
   real(r8) const(pver,ncol)    ! Gathered tracer array
   real(r8) fisg(pver,ncol)     ! gathered insoluble fraction of tracer
   real(r8) conu(pver,ncol)     ! Mix ratio in updraft at interfaces
   real(r8) dcondt(pver,ncol)   ! Gathered tend array
   real(r8) small               ! A small number
   real(r8) mbsth               ! Threshold for mass fluxes
   real(r8) mupdudp             ! A work variable
   real(r8) minc                ! A work variable
   real(r8) maxc                ! A work variable
   real(r8) fluxin              ! A work variable
   real(r8) fluxout             ! A work variable
   real(r8) netflux             ! A work variable

   real(r8) dutmp(pver,ncol)    ! Mass detraining from updraft
   real(r8) eutmp(pver,ncol)    ! Mass entraining from updraft
   real(r8) edtmp(pver,ncol)    ! Mass entraining from downdraft
   real(r8) dptmp(pver,ncol)    ! Delta pressure between interfaces


   small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
   mbsth = 1.e-15_r8

! Find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

! Loop ever each constituent
   do m = 2, ncnst
      if (doconvtran(m)) then
!------------------------------------------------------
!  Lixh closed this part, all tracers are wet.
!         if (cnst_get_type_byind(m).eq.'dry') then
!            do k = 1,pver
!               do i =il1g,il2g
!                  dptmp(k,i) = dpdry(k,i)
!                  dutmp(k,i) = du(k,i)*dp(k,i)/dpdry(k,i)
!                  eutmp(k,i) = eu(k,i)*dp(k,i)/dpdry(k,i)
!                  edtmp(k,i) = ed(k,i)*dp(k,i)/dpdry(k,i)
!               end do
!            end do
!         else
            do k = 1,pver
               do i =il1g,il2g
                  dptmp(k,i) = dp(k,i)
                  dutmp(k,i) = du(k,i)
                  eutmp(k,i) = eu(k,i)
                  edtmp(k,i) = ed(k,i)
               end do
            end do
!         endif
!------------------------------------------------------

!        dptmp = dp

! Gather up the constituent and set tend to zero
         do k = 1,pver
            do i =il1g,il2g
               const(k,i) = q(m,k,ideep(i))
               fisg(k,i) = fracis(m,k,ideep(i))
            end do
         end do

! From now on work only with gathered data

! Interpolate environment tracer values to interfaces
         do k = 1,pver
            km1 = max(1,k-1)
            do i = il1g, il2g
               minc = min(const(km1,i),const(k,i))
               maxc = max(const(km1,i),const(k,i))
               if (minc < 0) then
                  cdifr = 0._r8
               else
                  cdifr = abs(const(k,i)-const(km1,i))/max(maxc,small)
               endif

! If the two layers differ significantly use a geometric averaging
! procedure
               if (cdifr > 1.E-6_r8) then
                  cabv = max(const(km1,i),maxc*1.e-12_r8)
                  cbel = max(const(k,i),maxc*1.e-12_r8)
                  chat(k,i) = log(cabv/cbel)/(cabv-cbel)*cabv*cbel

               else             ! Small diff, so just arithmetic mean
                  chat(k,i) = 0.5_r8* (const(k,i)+const(km1,i))
               end if

! Provisional up and down draft values
               conu(k,i) = chat(k,i)
               cond(k,i) = chat(k,i)

!              provisional tends
               dcondt(k,i) = 0._r8

            end do
         end do

! Do levels adjacent to top and bottom
         k = 2
         km1 = 1
         kk = pver
         do i = il1g,il2g
            mupdudp = mu(kk,i) + dutmp(kk,i)*dptmp(kk,i)
            if (mupdudp > mbsth) then
               conu(kk,i) = (+eutmp(kk,i)*fisg(kk,i)*const(kk,i)*dptmp(kk,i))/mupdudp
            endif
            if (md(k,i) < -mbsth) then
               cond(k,i) =  (-edtmp(km1,i)*fisg(km1,i)*const(km1,i)*dptmp(km1,i))/md(k,i)
            endif
         end do

! Updraft from bottom to top
         do kk = pver-1,1,-1
            kkp1 = min(pver,kk+1)
            do i = il1g,il2g
               mupdudp = mu(kk,i) + dutmp(kk,i)*dptmp(kk,i)
               if (mupdudp > mbsth) then
                  conu(kk,i) = (  mu(kkp1,i)*conu(kkp1,i)+eutmp(kk,i)*fisg(kk,i)* &
                                  const(kk,i)*dptmp(kk,i) )/mupdudp
               endif
            end do
         end do

! Downdraft from top to bottom
         do k = 3,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (md(k,i) < -mbsth) then
                  cond(k,i) =  (  md(km1,i)*cond(km1,i)-edtmp(km1,i)*fisg(km1,i)*const(km1,i) &
                                  *dptmp(km1,i) )/md(k,i)
               endif
            end do
         end do


         do k = ktm,pver
            km1 = max(1,k-1)
            kp1 = min(pver,k+1)
            do i = il1g,il2g

! version 1 hard to check for roundoff errors
!               dcondt(k,i) =
!     $                  +(+mu(kp1,i)* (conu(kp1,i)-chat(kp1,i))
!     $                    -mu(k,i)*   (conu(k,i)-chat(k,i))
!     $                    +md(kp1,i)* (cond(kp1,i)-chat(kp1,i))
!     $                    -md(k,i)*   (cond(k,i)-chat(k,i))
!     $                   )/dp(k,i)

! version 2 hard to limit fluxes
!               fluxin =  mu(kp1,i)*conu(kp1,i) + mu(k,i)*chat(k,i)
!     $                 -(md(k,i)  *cond(k,i)   + md(kp1,i)*chat(kp1,i))
!               fluxout = mu(k,i)*conu(k,i)     + mu(kp1,i)*chat(kp1,i)
!     $                 -(md(kp1,i)*cond(kp1,i) + md(k,i)*chat(k,i))

! version 3 limit fluxes outside convection to mass in appropriate layer
! these limiters are probably only safe for positive definite quantitities
! it assumes that mu and md already satify a courant number limit of 1
               fluxin =  mu(kp1,i)*conu(kp1,i)+ mu(k,i)*min(chat(k,i),const(km1,i)) &
                         -(md(k,i)  *cond(k,i) + md(kp1,i)*min(chat(kp1,i),const(kp1,i)))
               fluxout = mu(k,i)*conu(k,i) + mu(kp1,i)*min(chat(kp1,i),const(k,i)) &
                         -(md(kp1,i)*cond(kp1,i) + md(k,i)*min(chat(k,i),const(k,i)))

               netflux = fluxin - fluxout
               if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                  netflux = 0._r8
               endif
               dcondt(k,i) = netflux/dptmp(k,i)
            end do
         end do
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
         do k = kbm,pver
            km1 = max(1,k-1)
            do i = il1g,il2g
               if (k == mx(i)) then

! version 1
!                  dcondt(k,i) = (1./dsubcld(i))*
!     $              (-mu(k,i)*(conu(k,i)-chat(k,i))
!     $               -md(k,i)*(cond(k,i)-chat(k,i))
!     $              )

! version 2
!                  fluxin =  mu(k,i)*chat(k,i) - md(k,i)*cond(k,i)
!                  fluxout = mu(k,i)*conu(k,i) - md(k,i)*chat(k,i)
! version 3
                  fluxin =  mu(k,i)*min(chat(k,i),const(km1,i)) - md(k,i)*cond(k,i)
                  fluxout = mu(k,i)*conu(k,i) - md(k,i)*min(chat(k,i),const(k,i))

                  netflux = fluxin - fluxout
                  if (abs(netflux) < max(fluxin,fluxout)*1.e-12_r8) then
                     netflux = 0._r8
                  endif
!                  dcondt(k,i) = netflux/dsubcld(i)
                  dcondt(k,i) = netflux/dptmp(k,i)
               else if (k > mx(i)) then
!                  dcondt(k,i) = dcondt(i,k-1)
                  dcondt(k,i) = 0._r8
               end if
            end do
         end do

! Initialize to zero everywhere, then scatter tendency back to full array
         dqdt(m,:,:) = 0._r8
         do k = 1,pver
            kp1 = min(pver,k+1)
            do i = il1g,il2g
               dqdt(m,k,ideep(i)) = dcondt(k,i)
            end do
         end do

      end if      ! for doconvtran

    end do

    return
    end subroutine convtran



! Purpose: Convective transport of momentum
!          Mixing ratios may be with respect to either dry or moist air
! Method:  Based on the convtran subroutine by P. Rasch
!          <Also include any applicable external references.>
! Author: J. Richter and P. Rasch
    subroutine momtran(domomtran, q      , ncnst  , mx     ,          &
                       il1g     , il2g   , dqdt   , pguall , pgdall , &
                       icwu     , icwd   , dt     , seten    )
 
    implicit none
! io
    integer,  intent(in)  :: ncnst                   ! number of tracers to transport
    real(r8), intent(in)  :: dt                      ! time step in seconds : 2*delta_t
    integer,  intent(in)  :: il1g                    ! Gathered min lon indices over which to operate
    integer,  intent(in)  :: il2g                    ! Gathered max lon indices over which to operate
    logical,  intent(in)  :: domomtran(ncnst)        ! flag for doing convective transport
    real(r8), intent(in)  :: q(ncnst,pver,ncol)      ! Wind array
    integer,  intent(in)  :: mx(ncol)                ! Index of cloud top for each column
    real(r8), intent(out) :: dqdt(ncnst,pver,ncol)   ! Tracer tendency array
    real(r8), intent(out) :: pguall(ncnst,pver,ncol) ! Apparent force from  updraft PG
    real(r8), intent(out) :: pgdall(ncnst,pver,ncol) ! Apparent force from  downdraft PG
    real(r8), intent(out) :: icwu(ncnst,pver,ncol)   ! In-cloud winds in updraft
    real(r8), intent(out) :: icwd(ncnst,pver,ncol)   ! In-cloud winds in downdraft
    real(r8), intent(out) :: seten(pver,ncol)        ! Dry static energy tendency
  
! local
    integer i                 ! Work index
    integer k                 ! Work index
    integer kbm               ! Highest altitude index of cloud base
    integer kk                ! Work index
    integer kkp1              ! Work index
    integer kkm1              ! Work index
    integer km1               ! Work index
    integer kp1               ! Work index
    integer ktm               ! Highest altitude index of cloud top
    integer m                 ! Work index
    integer ii                ! Work index
 
    real(r8) cabv                 ! Mix ratio of constituent above
    real(r8) cbel                 ! Mix ratio of constituent below
    real(r8) cdifr                ! Normalized diff between cabv and cbel
    real(r8) chat(pver,ncol)      ! Mix ratio in env at interfaces
    real(r8) cond(pver,ncol)      ! Mix ratio in downdraft at interfaces
    real(r8) const(pver,ncol)     ! Gathered wind array
    real(r8) conu(pver,ncol)      ! Mix ratio in updraft at interfaces
    real(r8) dcondt(pver,ncol)    ! Gathered tend array
    real(r8) small                ! A small number
    real(r8) mbsth                ! Threshold for mass fluxes
    real(r8) mupdudp              ! A work variable
    real(r8) minc                 ! A work variable
    real(r8) maxc                 ! A work variable
    real(r8) fluxin               ! A work variable
    real(r8) fluxout              ! A work variable
    real(r8) netflux              ! A work variable
    real(r8) momcu                ! constant for updraft pressure gradient term
    real(r8) momcd                ! constant for downdraft pressure gradient term
    real(r8) sum                  ! sum
    real(r8) sum2                 ! sum2
  
    real(r8) mududp(pver,ncol)    ! working variable
    real(r8) mddudp(pver,ncol)    ! working variable
    real(r8) pgu(pver,ncol)       ! Pressure gradient term for updraft
    real(r8) pgd(pver,ncol)       ! Pressure gradient term for downdraft
 
    real(r8) gseten(pver,ncol)         ! Gathered dry static energy tendency
    real(r8) mflux(ncnst,pverp,ncol)   ! Gathered momentum flux
    real(r8) wind0(ncnst,pver,ncol)    ! gathered  wind before time step
    real(r8) windf(ncnst,pver,ncol)    ! gathered  wind after time step
    real(r8) fkeb, fket, ketend_cons, ketend, utop, ubot, vtop, vbot, gset2
  
! Initialize outgoing fields
    pguall(:,:,:)     = 0.0_r8
    pgdall(:,:,:)     = 0.0_r8
! Initialize in-cloud winds to environmental wind
    icwu(:,:,:ncol)   = q(:,:,:ncol)
    icwd(:,:,:ncol)   = q(:,:,:ncol)
 
! Initialize momentum flux and  final winds
    mflux(:,:,:)      = 0.0_r8
    wind0(:,:,:)      = 0.0_r8
    windf(:,:,:)      = 0.0_r8
 
! Initialize dry static energy
    seten(:,:)        = 0.0_r8
    gseten(:,:)       = 0.0_r8
 
! Define constants for parameterization
    momcu = 0.4_r8
    momcd = 0.4_r8
    small = 1.e-36_r8
! mbsth is the threshold below which we treat the mass fluxes as zero (in mb/s)
    mbsth = 1.e-15_r8
 
! Find the highest level top and bottom levels of convection
    ktm = pver
    kbm = pver
    do i = il1g, il2g
       ktm = min(ktm,jt(i))
       kbm = min(kbm,mx(i))
    end do
 
! Loop ever each wind component
    do m = 1, ncnst                    !start at m = 1 to transport momentum
       if (domomtran(m)) then
 
! Gather up the winds and set tend to zero
          do k = 1,pver
             do i =il1g,il2g
                const(k,i) = q(m,k,ideep(i))
                 wind0(m,k,i) = const(k,i)
             end do
          end do
 
! From now on work only with gathered data
! Interpolate winds to interfaces
          do k = 1,pver
             km1 = max(1,k-1)
             do i = il1g, il2g
                ! use arithmetic mean
                chat(k,i) = 0.5_r8* (const(k,i)+const(km1,i))
                ! Provisional up and down draft values
                conu(k,i) = chat(k,i)
                cond(k,i) = chat(k,i)
                ! provisional tends
                dcondt(k,i) = 0._r8
             end do
          end do

! Pressure Perturbation Term
 
! Top boundary:  assume mu is zero 
          k=1
          pgu(k,:il2g) = 0.0_r8
          pgd(k,:il2g) = 0.0_r8
 
          do k=2,pver-1
             km1 = max(1,k-1)
             kp1 = min(pver,k+1)
             do i = il1g,il2g
                !interior points
                mududp(k,i) =  ( mu(k,i) * (const(k,i)- const(km1,i))/dp(km1,i) &
                            +  mu(kp1,i) * (const(kp1,i) - const(k,i))/dp(k,i))
                pgu(k,i) = - momcu * 0.5_r8 * mududp(k,i)
                mddudp(k,i) =  ( md(k,i) * (const(k,i)- const(km1,i))/dp(km1,i) &
                            +  md(kp1,i) * (const(kp1,i) - const(k,i))/dp(k,i))
                pgd(k,i) = - momcd * 0.5_r8 * mddudp(k,i)
             end do
          end do
 
        ! bottom boundary 
        k = pver
        km1 = max(1,k-1)
        do i=il1g,il2g
           mududp(k,i) =   mu(k,i) * (const(k,i)- const(km1,i))/dp(km1,i)
           pgu(k,i) = - momcu *  mududp(k,i)
           mddudp(k,i) =   md(k,i) * (const(k,i)- const(km1,i))/dp(km1,i) 
           pgd(k,i) = - momcd * mddudp(k,i)
        end do
        
 
! In-cloud velocity calculations
! Do levels adjacent to top and bottom
          k = 2
          km1 = 1
          kk = pver
          kkm1 = max(1,kk-1)
          do i = il1g,il2g
             mupdudp = mu(kk,i) + du(kk,i)*dp(kk,i)
             if (mupdudp > mbsth) then
                conu(kk,i) = (+eu(kk,i)*const(kk,i)*dp(kk,i)+pgu(kk,i)*dp(kk,i))/mupdudp
             endif
             if (md(k,i) < -mbsth) then
                cond(k,i) =  (-ed(km1,i)*const(km1,i)*dp(km1,i))-pgd(km1,i)*dp(km1,i)/md(k,i)
             endif
          end do
 
! Updraft from bottom to top
          do kk = pver-1,1,-1
             kkm1 = max(1,kk-1)
             kkp1 = min(pver,kk+1)
             do i = il1g,il2g
                mupdudp = mu(kk,i) + du(kk,i)*dp(kk,i)
                if (mupdudp > mbsth) then
                   conu(kk,i) = (  mu(kkp1,i)*conu(kkp1,i)+eu(kk,i)* &
                                   const(kk,i)*dp(kk,i)+pgu(kk,i)*dp(kk,i))/mupdudp
                endif
             end do
          end do
 
 
! Downdraft from top to bottom
          do k = 3,pver
             km1 = max(1,k-1)
             do i = il1g,il2g
                if (md(k,i) < -mbsth) then
                   cond(k,i) =  (  md(km1,i)*cond(km1,i)-ed(km1,i)*const(km1,i) &
                                   *dp(km1,i)-pgd(km1,i)*dp(km1,i) )/md(k,i)
                endif
             end do
          end do
          sum = 0._r8
          sum2 = 0._r8
          do k = ktm,pver
             km1 = max(1,k-1)
             kp1 = min(pver,k+1)
             do i = il1g,il2g
                ii = ideep(i)
                ! version 1 hard to check for roundoff errors
                dcondt(k,i) =  &
                            +(mu(kp1,i)* (conu(kp1,i)-chat(kp1,i)) &
                            -mu(k,i)*   (conu(k,i)-chat(k,i))      &
                            +md(kp1,i)* (cond(kp1,i)-chat(kp1,i)) &
                            -md(k,i)*   (cond(k,i)-chat(k,i)) &
                           )/dp(k,i)
 
             end do
          end do
 
          ! dcont for bottom layer
          do k = kbm,pver
             km1 = max(1,k-1)
             do i = il1g,il2g
                if (k == mx(i)) then
                   ! version 1
                   dcondt(k,i) = (1._r8/dp(k,i))*   &  
                        (-mu(k,i)*(conu(k,i)-chat(k,i)) &
                        -md(k,i)*(cond(k,i)-chat(k,i)) &
                        )
                end if
             end do
          end do
 
! Initialize to zero everywhere, then scatter tendency back to full array
          dqdt(m,:,:) = 0._r8
          do k = 1,pver
             do i = il1g,il2g
                ii = ideep(i)
                dqdt(m,k,ii) = dcondt(k,i)
                ! Output apparent force on the mean flow from pressure gradient
                pguall(m,k,ii) = -pgu(k,i)
                pgdall(m,k,ii) = -pgd(k,i)
                icwu(m,k,ii)   =  conu(k,i)
                icwd(m,k,ii)   =  cond(k,i)
             end do
          end do
 
          ! Calculate momentum flux in units of mb*m/s2 
          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                mflux(m,k,i) = &
                     -mu(k,i)*   (conu(k,i)-chat(k,i))      &
                     -md(k,i)*   (cond(k,i)-chat(k,i))
             end do
          end do
 
          ! Calculate winds at the end of the time step 
          do k = ktm,pver
             do i = il1g,il2g
                ii = ideep(i)
                km1 = max(1,k-1)
                kp1 = k+1
                windf(m,k,i) = const(k,i)    -   (mflux(m,kp1,i) - mflux(m,k,i)) * dt /dp(k,i)
             end do
          end do
 
        end if      ! for domomtran
    end do
 
  ! Need to add an energy fix to account for the dissipation of kinetic energy
     ! Formulation follows from Boville and Bretherton (2003)
     ! formulation by PJR
 
     do k = ktm,pver
        km1 = max(1,k-1)
        kp1 = min(pver,k+1)
        do i = il1g,il2g
 
           ii = ideep(i)
 
           ! calculate the KE fluxes at top and bot of layer 
           ! based on a discrete approximation to b&b eq(35) F_KE = u*F_u + v*F_v at interface
           utop = (wind0(1,k,i)+wind0(1,km1,i))/2._r8
           vtop = (wind0(2,k,i)+wind0(2,km1,i))/2._r8
           ubot = (wind0(1,kp1,i)+wind0(1,k,i))/2._r8
           vbot = (wind0(2,kp1,i)+wind0(2,k,i))/2._r8
           fket = utop*mflux(1,k,i)   + vtop*mflux(2,k,i)    ! top of layer
           fkeb = ubot*mflux(1,k+1,i) + vbot*mflux(2,k+1,i)  ! bot of layer
 
           ! divergence of these fluxes should give a conservative redistribution of KE
           ketend_cons = (fket-fkeb)/dp(k,i)
 
           ! tendency in kinetic energy resulting from the momentum transport
           ketend = ((windf(1,k,i)**2 + windf(2,k,i)**2) - (wind0(1,k,i)**2 + wind0(2,k,i)**2))*0.5_r8/dt
 
           ! the difference should be the dissipation
           gset2 = ketend_cons - ketend
           gseten(k,i) = gset2
        end do
     end do
 
     ! Scatter dry static energy to full array
     do k = 1,pver
        do i = il1g,il2g
           ii = ideep(i)
           seten(k,ii) = gseten(k,i)
 
        end do
     end do
 
    return
    end subroutine momtran



! Method: 
! may 09/91 - guang jun zhang, m.lazare, n.mcfarlane.
!             original version cldprop.
! 
! Author: See above, modified by P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! there are debug statements left strewn about and code segments disabled
! these are to facilitate future development. We expect to release a
! cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
    subroutine cldprp(q       ,t       ,u       ,v       ,p       , &
                      z       ,s       ,                            &
                      sd      ,qd      ,mc      ,                   &
                      qu      ,su      ,zf      ,qst     ,hmn     , &
                      hsat    ,shat    ,ql      ,                   &
                      cmeg    ,jb      ,lel     ,jlcl    ,          &
                      mx      ,j0      ,jd      ,rl      ,il2g    , &
                      rd      ,grav    ,cp      ,msg     ,          &
                      pflx    ,evp     ,cu      ,rprd    ,limcnv  ,landfrac)
   implicit none
! io
   real(r8), intent(in) :: rd                    ! gas constant for dry air
   real(r8), intent(in) :: grav                  ! gravy
   real(r8), intent(in) :: cp                    ! heat capacity of dry air
   real(r8), intent(in) :: q(pver, ncol)         ! spec. humidity of env
   real(r8), intent(in) :: t(pver, ncol)         ! temp of env
   real(r8), intent(in) :: p(pver, ncol)         ! pressure of env
   real(r8), intent(in) :: z(pver, ncol)         ! height of env
   real(r8), intent(in) :: s(pver, ncol)         ! normalized dry static energy of env
   real(r8), intent(in) :: zf(pverp, ncol)       ! height of interfaces
   real(r8), intent(in) :: u(pver, ncol)         ! zonal velocity of env
   real(r8), intent(in) :: v(pver, ncol)         ! merid. velocity of env
   real(r8), intent(in) :: landfrac(ncol)        ! RBN Landfrac
   integer, intent(in) :: jb(ncol)               ! updraft base level
   integer, intent(in) :: lel(ncol)              ! updraft launch level
   integer, intent(out) :: jlcl(ncol)            ! updraft lifting cond level
   integer, intent(in) :: mx(ncol)               ! updraft base level (same is jb)
   integer, intent(out) :: j0(ncol)              ! level where updraft begins detraining
   integer, intent(out) :: jd(ncol)              ! level of downdraft
   integer, intent(in) :: limcnv                 ! convection limiting level
   integer, intent(in) :: il2g                   !CORE GROUP REMOVE
   integer, intent(in) :: msg                    ! missing moisture vals (always 0)
   real(r8), intent(in) :: rl                    ! latent heat of vap
   real(r8), intent(in) :: shat(pver, ncol)      ! interface values of dry stat energy

   real(r8), intent(out) :: rprd(pver, ncol)     ! rate of production of precip at that layer
   real(r8), intent(out) :: hmn(pver, ncol)      ! moist stat energy of env
   real(r8), intent(out) :: hsat(pver, ncol)     ! sat moist stat energy of env
   real(r8), intent(out) :: mc(pver, ncol)       ! net mass flux
   real(r8), intent(out) :: pflx(pverp, ncol)    ! precipitation flux thru layer
   real(r8), intent(out) :: qd(pver, ncol)       ! spec humidity of downdraft
   real(r8), intent(out) :: ql(pver, ncol)       ! liq water of updraft
   real(r8), intent(out) :: qst(pver, ncol)      ! saturation mixing ratio of env.
   real(r8), intent(out) :: qu(pver, ncol)       ! spec hum of updraft
   real(r8), intent(out) :: sd(pver, ncol)       ! normalized dry stat energy of downdraft
   real(r8), intent(out) :: su(pver, ncol)       ! normalized dry stat energy of updraft
   real(r8), intent(out) :: cu(pver, ncol)
   real(r8), intent(out) :: evp(pver, ncol)
   real(r8), intent(out) :: cmeg(pver, ncol)

! Local
   real(r8) gamma(pver, ncol)
   real(r8) dz(pver, ncol)
   real(r8) iprm(pver, ncol)
   real(r8) hu(pver, ncol)
   real(r8) hd(pver, ncol)
   real(r8) eps(pver, ncol)
   real(r8) f(pver, ncol)
   real(r8) k1(pver, ncol)
   real(r8) i2(pver, ncol)
   real(r8) ihat(pver, ncol)
   real(r8) i3(pver, ncol)
   real(r8) idag(pver, ncol)
   real(r8) i4(pver, ncol)
   real(r8) qsthat(pver, ncol)
   real(r8) hsthat(pver, ncol)
   real(r8) gamhat(pver, ncol)
   real(r8) qds(pver, ncol)
! RBN For c0mask
   real(r8) c0mask(ncol)

   real(r8) hmin(ncol)
   real(r8) expdif(ncol)
   real(r8) expnum(ncol)
   real(r8) ftemp(ncol)
   real(r8) eps0(ncol)
   real(r8) rmue(ncol)
   real(r8) zuef(ncol)
   real(r8) zdef(ncol)
   real(r8) epsm(ncol)
   real(r8) ratmjb(ncol)
   real(r8) est(ncol)
   real(r8) totpcp(ncol)
   real(r8) totevp(ncol)
   real(r8) alfa(ncol)
   real(r8) ql1
   real(r8) tu
   real(r8) estu
   real(r8) qstu

   real(r8) small
   real(r8) mdt

   integer khighest
   integer klowest
   integer kount
   integer i,k

   logical doit(ncol)
   logical done(ncol)

   do i = 1,il2g
      ftemp(i) = 0._r8
      expnum(i) = 0._r8
      expdif(i) = 0._r8
      c0mask(i)  = c0_ocn * (1._r8-landfrac(i)) +   c0_lnd * landfrac(i) 
   end do

!jr Change from msg+1 to 1 to prevent blowup
   do k = 1,pver
      do i = 1,il2g
         dz(k,i) = zf(k,i) - zf(k+1,i)
      end do
   end do

! initialize many output and work variables to zero
   pflx(1,:il2g) = 0

   do k = 1,pver
      do i = 1,il2g
         k1(k,i) = 0._r8
         i2(k,i) = 0._r8
         i3(k,i) = 0._r8
         i4(k,i) = 0._r8
         mu(k,i) = 0._r8
         f(k,i) = 0._r8
         eps(k,i) = 0._r8
         eu(k,i) = 0._r8
         du(k,i) = 0._r8
         ql(k,i) = 0._r8
         cu(k,i) = 0._r8
         evp(k,i) = 0._r8
         cmeg(k,i) = 0._r8
         qds(k,i) = q(k,i)
         md(k,i) = 0._r8
         ed(k,i) = 0._r8
         sd(k,i) = s(k,i)
         qd(k,i) = q(k,i)
         mc(k,i) = 0._r8
         qu(k,i) = q(k,i)
         su(k,i) = s(k,i)
         call qmmr_hPa(t(k,i), p(k,i), est(i), qst(k,i))
!++bee
         if ( p(k,i)-est(i) <= 0._r8 ) then
            qst(k,i) = 1.0_r8
         end if
!--bee
         gamma(k,i) = qst(k,i)*(1._r8 + qst(k,i)/eps1)*eps1*rl/(rd*t(k,i)**2)*rl/cp
         hmn(k,i) = cp*t(k,i) + grav*z(k,i) + rl*q(k,i)
         hsat(k,i) = cp*t(k,i) + grav*z(k,i) + rl*qst(k,i)
         hu(k,i) = hmn(k,i)
         hd(k,i) = hmn(k,i)
         rprd(k,i) = 0._r8
      end do
   end do

!jr Set to zero things which make this routine blow up
   do k=1,msg
      do i=1,il2g
         rprd(k,i) = 0._r8
      end do
   end do

! interpolate the layer values of qst, hsat and gamma to layer interfaces
   do k = 1, msg+1
      do i = 1,il2g
         hsthat(k,i) = hsat(k,i)
         qsthat(k,i) = qst(k,i)
         gamhat(k,i) = gamma(k,i)
      end do
   end do
   do i = 1,il2g
      totpcp(i) = 0._r8
      totevp(i) = 0._r8
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (abs(qst(k-1,i)-qst(k,i)) > 1.E-6_r8) then
            qsthat(k,i) = log(qst(k-1,i)/qst(k,i))*qst(k-1,i)*qst(k,i)/ (qst(k-1,i)-qst(k,i))
         else
            qsthat(k,i) = qst(k,i)
         end if
         hsthat(k,i) = cp*shat(k,i) + rl*qsthat(k,i)
         if (abs(gamma(k-1,i)-gamma(k,i)) > 1.E-6_r8) then
            gamhat(k,i) = log(gamma(k-1,i)/gamma(k,i))*gamma(k-1,i)*gamma(k,i)/ &
                                (gamma(k-1,i)-gamma(k,i))
         else
            gamhat(k,i) = gamma(k,i)
         end if
      end do
   end do

! initialize cloud top to highest plume top.
!jr changed hard-wired 4 to limcnv+1 (not to exceed pver)
   jt(:) = pver
   do i = 1,il2g
      jt(i) = max(lel(i),limcnv+1)
      jt(i) = min(jt(i),pver)
      jd(i) = pver
      jlcl(i) = lel(i)
      hmin(i) = 1.E6_r8
   end do

! find the level of minimum hsat, where detrainment starts
   do k = msg + 1,pver
      do i = 1,il2g
         if (hsat(k,i) <= hmin(i) .and. k >= jt(i) .and. k <= jb(i)) then
            hmin(i) = hsat(k,i)
            j0(i) = k
         end if
      end do
   end do
   do i = 1,il2g
      j0(i) = min(j0(i),jb(i)-2)
      j0(i) = max(j0(i),jt(i)+2)

! Fix from Guang Zhang to address out of bounds array reference
      j0(i) = min(j0(i),pver)
   end do

! Initialize certain arrays inside cloud
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= jb(i)) then
            hu(k,i) = hmn(mx(i),i) + cp*tiedke_add
            su(k,i) = s(mx(i),i) + tiedke_add
         end if
      end do
   end do

! *********************************************************
! compute taylor series for approximate eps(z) below
! *********************************************************
   do k = pver - 1,msg + 1,-1
      do i = 1,il2g
         if (k < jb(i) .and. k >= jt(i)) then
            k1(k,i) = k1(k+1,i) + (hmn(mx(i),i)-hmn(k,i))*dz(k,i)
            ihat(k,i) = 0.5_r8* (k1(k+1,i)+k1(k,i))
            i2(k,i) = i2(k+1,i) + ihat(k,i)*dz(k,i)
            idag(k,i) = 0.5_r8* (i2(k+1,i)+i2(k,i))
            i3(k,i) = i3(k+1,i) + idag(k,i)*dz(k,i)
            iprm(k,i) = 0.5_r8* (i3(k+1,i)+i3(k,i))
            i4(k,i) = i4(k+1,i) + iprm(k,i)*dz(k,i)
         end if
      end do
   end do

! re-initialize hmin array for ensuing calculation.
   do i = 1,il2g
      hmin(i) = 1.E6_r8
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i) .and. hmn(k,i) <= hmin(i)) then
            hmin(i) = hmn(k,i)
            expdif(i) = hmn(mx(i),i) - hmin(i)
         end if
      end do
   end do

! *********************************************************
! compute approximate eps(z) using above taylor series
! *********************************************************

   do k = msg + 2,pver
      do i = 1,il2g
         expnum(i) = 0._r8
         ftemp(i) = 0._r8
         if (k < jt(i) .or. k >= jb(i)) then
            k1(k,i) = 0._r8
            expnum(i) = 0._r8
         else
            expnum(i) = hmn(mx(i),i) - (hsat(k-1,i)*(zf(k,i)-z(k,i)) + &
                        hsat(k,i)* (z(k-1,i)-zf(k,i)))/(z(k-1,i)-z(k,i))
         end if
         if ((expdif(i) > 100._r8 .and. expnum(i) > 0._r8) .and. &
            k1(k,i) > expnum(i)*dz(k,i)) then
            ftemp(i) = expnum(i)/k1(k,i)
            f(k,i) = ftemp(i) + i2(k,i)/k1(k,i)*ftemp(i)**2 + &
                     (2._r8*i2(k,i)**2-k1(k,i)*i3(k,i))/k1(k,i)**2* &
                     ftemp(i)**3 + (-5._r8*k1(k,i)*i2(k,i)*i3(k,i)+ &
                     5._r8*i2(k,i)**3+k1(k,i)**2*i4(k,i))/ &
                     k1(k,i)**3*ftemp(i)**4
            f(k,i) = max(f(k,i),0._r8)
            f(k,i) = min(f(k,i),0.0002_r8)
         end if
      end do
   end do
   do i = 1,il2g
      if (j0(i) < jb(i)) then
         if (f(j0(i),i) < 1.E-6_r8 .and. f(j0(i)+1,i) > f(j0(i),i)) j0(i) = j0(i) + 1
      end if
   end do
   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jt(i) .and. k <= j0(i)) then
            f(k,i) = max(f(k,i),f(k-1,i))
         end if
      end do
   end do
   do i = 1,il2g
      eps0(i) = f(j0(i),i)
      eps(jb(i),i) = eps0(i)
   end do

! This is set to match the Rasch and Kristjansson paper
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= j0(i) .and. k <= jb(i)) then
            eps(k,i) = f(j0(i),i)
         end if
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k < j0(i) .and. k >= jt(i)) eps(k,i) = f(k,i)
      end do
   end do

! specify the updraft mass flux mu, entrainment eu, detrainment du
! and moist static energy hu.
! here and below mu, eu,du, md and ed are all normalized by mb
   do i = 1,il2g
      if (eps0(i) > 0._r8) then
         mu(jb(i),i) = 1._r8
         eu(jb(i),i) = mu(jb(i),i)/dz(jb(i),i)
      end if
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (eps0(i) > 0._r8 .and. (k >= jt(i) .and. k < jb(i))) then
            zuef(i) = zf(k,i) - zf(jb(i),i)
            rmue(i) = (1._r8/eps0(i))* (exp(eps(k+1,i)*zuef(i))-1._r8)/zuef(i)
            mu(k,i) = (1._r8/eps0(i))* (exp(eps(k,i  )*zuef(i))-1._r8)/zuef(i)
            eu(k,i) = (rmue(i)-mu(k+1,i))/dz(k,i)
            du(k,i) = (rmue(i)-mu(k,i))/dz(k,i)
         end if
      end do
   end do

   khighest = pverp
   klowest = 1
   do i=1,il2g
      khighest = min(khighest,lel(i))
      klowest = max(klowest,jb(i))
   end do
   do k = klowest-1,khighest,-1
      do i = 1,il2g
         if (k <= jb(i)-1 .and. k >= lel(i) .and. eps0(i) > 0._r8) then
            if (mu(k,i) < 0.02_r8) then
               hu(k,i) = hmn(k,i)
               mu(k,i) = 0._r8
               eu(k,i) = 0._r8
               du(k,i) = mu(k+1,i)/dz(k,i)
            else
               hu(k,i) = mu(k+1,i)/mu(k,i)*hu(k+1,i) + &
                         dz(k,i)/mu(k,i)* (eu(k,i)*hmn(k,i)- du(k,i)*hsat(k,i))
            end if
         end if
      end do
   end do

! reset cloud top index beginning from two layers above the
! cloud base (i.e. if cloud is only one layer thick, top is not reset
   do i=1,il2g
      doit(i) = .true.
   end do
   do k=klowest-2,khighest-1,-1
      do i=1,il2g
         if (doit(i) .and. k <= jb(i)-2 .and. k >= lel(i)-1) then
           if (hu(k,i) <= hsthat(k,i) .and. hu(k+1,i) > hsthat(k+1,i) &
           .and. mu(k,i) >= 0.02_r8) then
               if (hu(k,i)-hsthat(k,i) < -2000._r8) then
                  jt(i) = k + 1
                  doit(i) = .false.
               else
                  jt(i) = k
                  doit(i) = .false.
               end if
            else if (hu(k,i) > hu(jb(i),i) .or. mu(k,i) < 0.02_r8) then
               jt(i) = k + 1
               doit(i) = .false.
            end if
         end if
      end do
   end do
   do k = pver,msg + 1,-1
      do i = 1,il2g
         if (k >= lel(i) .and. k <= jt(i) .and. eps0(i) > 0._r8) then
            mu(k,i) = 0._r8
            eu(k,i) = 0._r8
            du(k,i) = 0._r8
            hu(k,i) = hmn(k,i)
         end if
         if (k == jt(i) .and. eps0(i) > 0._r8) then
            du(k,i) = mu(k+1,i)/dz(k,i)
            eu(k,i) = 0._r8
            mu(k,i) = 0._r8
         end if
      end do
   end do

! specify downdraft properties (no downdrafts if jd.ge.jb).
! scale down downward mass flux profile so that net flux
! (up-down) at cloud base in not negative.
   do i = 1,il2g
! in normal downdraft strength run alfa=0.2.  In test4 alfa=0.1
! alfa: 0.05-0.6, LiXH
      alfa(i) = 0.1_r8
      jt(i) = min(jt(i),jb(i)-1)
      jd(i) = max(j0(i),jt(i)+1)
      jd(i) = min(jd(i),jb(i))
      hd(jd(i),i) = hmn(jd(i)-1,i)
      if (jd(i) < jb(i) .and. eps0(i) > 0._r8) then
         epsm(i) = eps0(i)
         md(jd(i),i) = -alfa(i)*epsm(i)/eps0(i)
      end if
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k > jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8) then
            zdef(i) = zf(jd(i),i) - zf(k,i)
            md(k,i) = -alfa(i)/ (2._r8*eps0(i))*(exp(2._r8*epsm(i)*zdef(i))-1._r8)/zdef(i)
         end if
      end do
   end do
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            ratmjb(i) = min(abs(mu(jb(i),i)/md(jb(i),i)),1._r8)
            md(k,i) = md(k,i)*ratmjb(i)
         end if
      end do
   end do

   small = 1.e-20_r8
   do k = msg + 1,pver
      do i = 1,il2g
         if ((k >= jt(i) .and. k <= pver) .and. eps0(i) > 0._r8) then
            ed(k-1,i) = (md(k-1,i)-md(k,i))/dz(k-1,i)
            mdt = min(md(k,i),-small)
            hd(k,i) = (md(k-1,i)*hd(k-1,i) - dz(k-1,i)*ed(k-1,i)*hmn(k-1,i))/mdt
         end if
      end do
   end do

! calculate updraft and downdraft properties.
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k >= jd(i) .and. k <= jb(i)) .and. eps0(i) > 0._r8 .and. jd(i) < jb(i)) then
            qds(k,i) = qsthat(k,i) + gamhat(k,i)*(hd(k,i)-hsthat(k,i))/ &
               (rl*(1._r8 + gamhat(k,i)))
         end if
      end do
   end do

   do i = 1,il2g
      done(i) = .false.
   end do
   kount = 0
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (k == jb(i) .and. eps0(i) > 0._r8) then
            qu(k,i) = q(mx(i),i)
            su(k,i) = (hu(k,i)-rl*qu(k,i))/cp
         end if
         if (( .not. done(i) .and. k > jt(i) .and. k < jb(i)) .and. eps0(i) > 0._r8) then
            su(k,i) = mu(k+1,i)/mu(k,i)*su(k+1,i) + &
                      dz(k,i)/mu(k,i)* (eu(k,i)-du(k,i))*s(k,i)
            qu(k,i) = mu(k+1,i)/mu(k,i)*qu(k+1,i) + dz(k,i)/mu(k,i)* (eu(k,i)*q(k,i)- &
                            du(k,i)*qst(k,i))
            tu = su(k,i) - grav/cp*zf(k,i)
            call qmmr_hPa(tu, (p(k,i)+p(k-1,i))/2._r8, estu, qstu)
            if (qu(k,i) >= qstu) then
               jlcl(i) = k
               kount = kount + 1
               done(i) = .true.
            end if
         end if
      end do
      if (kount >= il2g) goto 690
   end do
690 continue
   do k = msg + 2,pver
      do i = 1,il2g
         if ((k > jt(i) .and. k <= jlcl(i)) .and. eps0(i) > 0._r8) then
            su(k,i) = shat(k,i) + (hu(k,i)-hsthat(k,i))/(cp* (1._r8+gamhat(k,i)))
            qu(k,i) = qsthat(k,i) + gamhat(k,i)*(hu(k,i)-hsthat(k,i))/ &
                     (rl* (1._r8+gamhat(k,i)))
         end if
      end do
   end do

! compute condensation in updraft
   do k = pver,msg + 2,-1
      do i = 1,il2g
         if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8) then
            cu(k,i) = ((mu(k,i)*su(k,i)-mu(k+1,i)*su(k+1,i))/ &
                      dz(k,i)- (eu(k,i)-du(k,i))*s(k,i))/(rl/cp)
            if (k == jt(i)) cu(k,i) = 0._r8
            cu(k,i) = max(0._r8,cu(k,i))
         end if
      end do
   end do

! compute condensed liquid, rain production rate
! accumulate total precipitation (condensation - detrainment of liquid)
! Note ql1 = ql(k) + rprd(k)*dz(k)/mu(k)
! The differencing is somewhat strange (e.g. du(k,i)*ql(k+1,i)) but is
! consistently applied.
!    mu, ql are interface quantities
!    cu, du, eu, rprd are midpoint quantites
   do k = pver,msg + 2,-1
      do i = 1,il2g
         rprd(k,i) = 0._r8
         if (k >= jt(i) .and. k < jb(i) .and. eps0(i) > 0._r8 .and. mu(k,i) >= 0.0_r8) then
            if (mu(k,i) > 0._r8) then
               ql1 = 1._r8/mu(k,i)* (mu(k+1,i)*ql(k+1,i)- &
                     dz(k,i)*du(k,i)*ql(k+1,i)+dz(k,i)*cu(k,i))
               ql(k,i) = ql1/ (1._r8+dz(k,i)*c0mask(i))
            else
               ql(k,i) = 0._r8
            end if
            totpcp(i) = totpcp(i) + dz(k,i)*(cu(k,i)-du(k,i)*ql(k+1,i))
            rprd(k,i) = c0mask(i)*mu(k,i)*ql(k,i)
         end if
      end do
   end do

   do i = 1,il2g
      qd(jd(i),i) = qds(jd(i),i)
      sd(jd(i),i) = (hd(jd(i),i) - rl*qd(jd(i),i))/cp
   end do

   do k = msg + 2,pver
      do i = 1,il2g
         if (k >= jd(i) .and. k < jb(i) .and. eps0(i) > 0._r8) then
            qd(k+1,i) = qds(k+1,i)
            evp(k,i) = -ed(k,i)*q(k,i) + (md(k,i)*qd(k,i)-md(k+1,i)*qd(k+1,i))/dz(k,i)
            evp(k,i) = max(evp(k,i),0._r8)
            mdt = min(md(k+1,i),-small)
            sd(k+1,i) = ((rl/cp*evp(k,i)-ed(k,i)*s(k,i))*dz(k,i) + md(k,i)*sd(k,i))/mdt
            totevp(i) = totevp(i) - dz(k,i)*ed(k,i)*q(k,i)
         end if
      end do
   end do
   do i = 1,il2g
!*guang         totevp(i) = totevp(i) + md(jd(i),i)*q(jd(i)-1,i) -
      totevp(i) = totevp(i) + md(jd(i),i)*qd(jd(i),i) - md(jb(i),i)*qd(jb(i),i)
   end do
!!$   if (.true.) then
   if (.false.) then
      do i = 1,il2g
         k = jb(i)
         if (eps0(i) > 0._r8) then
            evp(k,i) = -ed(k,i)*q(k,i) + (md(k,i)*qd(k,i))/dz(k,i)
            evp(k,i) = max(evp(k,i),0._r8)
            totevp(i) = totevp(i) - dz(k,i)*ed(k,i)*q(k,i)
         end if
      end do
   endif

   do i = 1,il2g
      totpcp(i) = max(totpcp(i),0._r8)
      totevp(i) = max(totevp(i),0._r8)
   end do

   do k = msg + 2,pver
      do i = 1,il2g
         if (totevp(i) > 0._r8 .and. totpcp(i) > 0._r8) then
            md(k,i)  = md (k,i)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            ed(k,i)  = ed (k,i)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
            evp(k,i) = evp(k,i)*min(1._r8, totpcp(i)/(totevp(i)+totpcp(i)))
         else
            md(k,i) = 0._r8
            ed(k,i) = 0._r8
            evp(k,i) = 0._r8
         end if
! cmeg is the cloud water condensed - rain water evaporated
! rprd is the cloud water converted to rain - (rain evaporated)
         cmeg(k,i) = cu(k,i) - evp(k,i)
         rprd(k,i) = rprd(k,i)-evp(k,i)
      end do
   end do

! compute the net precipitation flux across interfaces
   pflx(1,:il2g) = 0._r8
   do k = 2,pverp
      do i = 1,il2g
         pflx(k,i) = pflx(k-1,i) + rprd(k-1,i)*dz(k-1,i)
      end do
   end do

   do k = msg + 1,pver
      do i = 1,il2g
         mc(k,i) = mu(k,i) + md(k,i)
      end do
   end do

   return
   end subroutine cldprp


! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: G. Zhang and collaborators. CCM contact:P. Rasch
! This is contributed code not fully standardized by the CCM core group.
!
! this code is very much rougher than virtually anything else in the CCM
! We expect to release cleaner code in a future release
!
! the documentation has been enhanced to the degree that we are able
    subroutine closure(q       ,t       ,p       ,z       ,s       , &
                       tp      ,qs      ,qu      ,su      ,mc      , &
                       qd      ,sd      ,                            &
                       qhat    ,shat    ,qstp    ,zf      ,          &
                       ql      ,mb      ,cape    ,tl      ,          &
                       lcl     ,lel     ,mx      ,il1g    ,          &
                       il2g    ,rd      ,grav    ,cp      ,rl      , &
                       msg     )

   implicit none
! io
   real(r8), intent(in) :: cp
   real(r8), intent(in) :: grav
   real(r8), intent(in) :: rd
   real(r8), intent(in) :: rl
   integer, intent(in)  :: msg
   integer, intent(in)  :: il1g
   integer, intent(in)  :: il2g
   real(r8), intent(inout) :: q(pver, ncol)     ! spec humidity
   real(r8), intent(inout) :: t(pver, ncol)     ! temperature
   real(r8), intent(inout) :: p(pver, ncol)     ! pressure (mb)
   real(r8), intent(inout) :: mb(ncol)         ! cloud base mass flux
   real(r8), intent(in) :: z(pver, ncol)        ! height (m)
   real(r8), intent(in) :: s(pver, ncol)        ! normalized dry static energy
   real(r8), intent(in) :: tp(pver, ncol)       ! parcel temp
   real(r8), intent(in) :: qs(pver, ncol)       ! sat spec humidity
   real(r8), intent(in) :: qu(pver, ncol)       ! updraft spec. humidity
   real(r8), intent(in) :: su(pver, ncol)       ! normalized dry stat energy of updraft
   real(r8), intent(in) :: mc(pver, ncol)       ! net convective mass flux
   real(r8), intent(in) :: qd(pver, ncol)       ! spec. humidity of downdraft
   real(r8), intent(in) :: sd(pver, ncol)       ! dry static energy of downdraft
   real(r8), intent(in) :: qhat(pver, ncol)     ! environment spec humidity at interfaces
   real(r8), intent(in) :: shat(pver, ncol)     ! env. normalized dry static energy at intrfcs
   real(r8), intent(in) :: qstp(pver, ncol)     ! spec humidity of parcel
   real(r8), intent(in) :: zf(pver+1, ncol)     ! height of interface levels
   real(r8), intent(in) :: ql(pver, ncol)       ! liquid water mixing ratio

   real(r8), intent(in) :: cape(ncol)          ! available pot. energy of column
   real(r8), intent(in) :: tl(ncol)

   integer, intent(in) :: lcl(ncol)        ! index of lcl
   integer, intent(in) :: lel(ncol)        ! index of launch leve
   integer, intent(in) :: mx(ncol)         ! base of updraft

! local
   real(r8) dtpdt(pver, ncol)
   real(r8) dqsdtp(pver, ncol)
   real(r8) dtmdt(pver, ncol)
   real(r8) dqmdt(pver, ncol)
   real(r8) dboydt(pver, ncol)
   real(r8) thetavp(pver, ncol)
   real(r8) thetavm(pver, ncol)

   real(r8) dtbdt(ncol),dqbdt(ncol),dtldt(ncol)
   real(r8) beta
   real(r8) dadt(ncol)
   real(r8) debdt
   real(r8) dltaa
   real(r8) eb

   integer i
   integer k, kmin, kmax

! change of subcloud layer properties due to convection is
! related to cumulus updrafts and downdrafts.
! mc(z)=f(z)*mb, mub=betau*mb, mdb=betad*mb are used
! to define betau, betad and f(z).
! note that this implies all time derivatives are in effect
! time derivatives per unit cloud-base mass flux, i.e. they
! have units of 1/mb instead of 1/sec.

   do i = il1g,il2g
      mb(i) = 0._r8
      eb = p(mx(i),i)*q(mx(i),i)/ (eps1+q(mx(i),i))
      dtbdt(i) = (1._r8/dsubcld(i))* (mu(mx(i),i)*(shat(mx(i),i)-su(mx(i),i))+ &
                  md(mx(i),i)* (shat(mx(i),i)-sd(mx(i),i)))
      dqbdt(i) = (1._r8/dsubcld(i))* (mu(mx(i),i)*(qhat(mx(i),i)-qu(mx(i),i))+ &
                 md(mx(i),i)* (qhat(mx(i),i)-qd(mx(i),i)))
      debdt = eps1*p(mx(i),i)/ (eps1+q(mx(i),i))**2*dqbdt(i)
      dtldt(i) = -2840._r8* (3.5_r8/t(mx(i),i)*dtbdt(i)-debdt/eb)/ &
                 (3.5_r8*log(t(mx(i),i))-log(eb)-4.805_r8)**2
   end do

!   dtmdt and dqmdt are cumulus heating and drying.
   do k = msg + 1,pver
      do i = il1g,il2g
         dtmdt(k,i) = 0._r8
         dqmdt(k,i) = 0._r8
      end do
   end do

   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k == jt(i)) then
            dtmdt(k,i) = (1._r8/dp(k,i))*(mu(k+1,i)* (su(k+1,i)-shat(k+1,i)- &
                          rl/cp*ql(k+1,i))+md(k+1,i)* (sd(k+1,i)-shat(k+1,i)))
            dqmdt(k,i) = (1._r8/dp(k,i))*(mu(k+1,i)* (qu(k+1,i)- &
                         qhat(k+1,i)+ql(k+1,i))+md(k+1,i)*(qd(k+1,i)-qhat(k+1,i)))
         end if
      end do
   end do

   beta = 0._r8
   do k = msg + 1,pver - 1
      do i = il1g,il2g
         if (k > jt(i) .and. k < mx(i)) then
            dtmdt(k,i) = (mc(k,i)* (shat(k,i)-s(k,i))+mc(k+1,i)* (s(k,i)-shat(k+1,i)))/ &
                         dp(k,i) - rl/cp*du(k,i)*(beta*ql(k,i)+ (1-beta)*ql(k+1,i))
!          dqmdt(k,i)=(mc(k,i)*(qhat(k,i)-q(k,i))
!     1                +mc(k+1,i)*(q(k,i)-qhat(k+1,i)))/dp(k,i)
!     2                +du(k,i)*(qs(k,i)-q(k,i))
!     3                +du(k,i)*(beta*ql(k,i)+(1-beta)*ql(k+1,i))

            dqmdt(k,i) = (mu(k+1,i)* (qu(k+1,i)-qhat(k+1,i)+cp/rl* (su(k+1,i)-s(k,i)))- &
                          mu(k,i)* (qu(k,i)-qhat(k,i)+cp/rl*(su(k,i)-s(k,i)))+md(k+1,i)* &
                         (qd(k+1,i)-qhat(k+1,i)+cp/rl*(sd(k+1,i)-s(k,i)))-md(k,i)* &
                         (qd(k,i)-qhat(k,i)+cp/rl*(sd(k,i)-s(k,i))))/dp(k,i) + &
                          du(k,i)* (beta*ql(k,i)+(1-beta)*ql(k+1,i))
         end if
      end do
   end do

   do k = msg + 1,pver
      do i = il1g,il2g
         if (k >= lel(i) .and. k <= lcl(i)) then
            thetavp(k,i) = tp(k,i)* (1000._r8/p(k,i))** (rd/cp)*(1._r8+1.608_r8*qstp(k,i)-q(mx(i),i))
            thetavm(k,i) = t(k,i)* (1000._r8/p(k,i))** (rd/cp)*(1._r8+0.608_r8*q(k,i))
            dqsdtp(k,i) = qstp(k,i)* (1._r8+qstp(k,i)/eps1)*eps1*rl/(rd*tp(k,i)**2)

! dtpdt is the parcel temperature change due to change of
! subcloud layer properties during convection.
            dtpdt(k,i) = tp(k,i)/ (1._r8+rl/cp* (dqsdtp(k,i)-qstp(k,i)/tp(k,i)))* &
                        (dtbdt(i)/t(mx(i),i)+rl/cp* (dqbdt(i)/tl(i)-q(mx(i),i)/ &
                         tl(i)**2*dtldt(i)))

! dboydt is the integrand of cape change.
            dboydt(k,i) = ((dtpdt(k,i)/tp(k,i)+1._r8/(1._r8+1.608_r8*qstp(k,i)-q(mx(i),i))* &
                          (1.608_r8 * dqsdtp(k,i) * dtpdt(k,i) -dqbdt(i))) - (dtmdt(k,i)/t(k,i)+0.608_r8/ &
                          (1._r8+0.608_r8*q(k,i))*dqmdt(k,i)))*grav*thetavp(k,i)/thetavm(k,i)
         end if
      end do
   end do

   do k = msg + 1,pver
      do i = il1g,il2g
         if (k > lcl(i) .and. k < mx(i)) then
            thetavp(k,i) = tp(k,i)* (1000._r8/p(k,i))** (rd/cp)*(1._r8+0.608_r8*q(mx(i),i))
            thetavm(k,i) = t(k,i)* (1000._r8/p(k,i))** (rd/cp)*(1._r8+0.608_r8*q(k,i))

! dboydt is the integrand of cape change.
            dboydt(k,i) = (dtbdt(i)/t(mx(i),i)+0.608_r8/ (1._r8+0.608_r8*q(mx(i),i))*dqbdt(i)- &
                          dtmdt(k,i)/t(k,i)-0.608_r8/ (1._r8+0.608_r8*q(k,i))*dqmdt(k,i))* &
                          grav*thetavp(k,i)/thetavm(k,i)
         end if
      end do
   end do


! buoyant energy change is set to 2/3*excess cape per 3 hours
   dadt(il1g:il2g)  = 0._r8
   kmin = minval(lel(il1g:il2g))
   kmax = maxval(mx(il1g:il2g)) - 1
   do k = kmin, kmax
      do i = il1g,il2g
         if ( k >= lel(i) .and. k <= mx(i) - 1) then
            dadt(i) = dadt(i) + dboydt(k,i)* (zf(k,i)-zf(k+1,i))
         endif
      end do
   end do
   do i = il1g,il2g
      dltaa = -1._r8* (cape(i)-capelmt)
      if (dadt(i) /= 0._r8) mb(i) = max(dltaa/tau/dadt(i),0._r8)
   end do

   return
   end subroutine closure

 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: phil rasch dec 19 1995
subroutine q1q2_pjr(dqdt    ,dsdt    ,q       ,qs      ,qu      , &
                    su      ,qhat    ,shat    ,                   &
                    sd      ,qd      ,ql      ,                   &
                    mx      ,il1g    ,il2g    ,                   &
                    cp      ,rl      ,msg     ,                   &
                    dl      ,evp     ,cu      )


   implicit none
! io
   real(r8), intent(in) :: cp
   integer, intent(in)  :: il1g
   integer, intent(in)  :: il2g
   integer, intent(in)  :: msg
   integer, intent(in)  :: mx(ncol)
   real(r8), intent(in) :: rl
   real(r8), intent(in) :: q(pver, ncol)
   real(r8), intent(in) :: qs(pver, ncol)
   real(r8), intent(in) :: qu(pver, ncol)
   real(r8), intent(in) :: su(pver, ncol)
   real(r8), intent(in) :: qhat(pver, ncol)
   real(r8), intent(in) :: shat(pver, ncol)
   real(r8), intent(in) :: sd(pver, ncol)
   real(r8), intent(in) :: qd(pver, ncol)
   real(r8), intent(in) :: ql(pver, ncol)
   real(r8), intent(in) :: evp(pver, ncol)
   real(r8), intent(in) :: cu(pver, ncol)
   real(r8),intent(out) :: dqdt(pver, ncol),dsdt(pver, ncol)
   real(r8),intent(out) :: dl(pver, ncol)
! local
   integer kbm
   integer ktm
   integer i
   integer k
   real(r8) emc

   do k = msg + 1,pver
      do i = il1g,il2g
         dsdt(k,i) = 0._r8
         dqdt(k,i) = 0._r8
         dl(k,i) = 0._r8
      end do
   end do

! find the highest level top and bottom levels of convection
   ktm = pver
   kbm = pver
   do i = il1g, il2g
      ktm = min(ktm,jt(i))
      kbm = min(kbm,mx(i))
   end do

   do k = ktm,pver-1
      do i = il1g,il2g
         emc = -cu (k,i)               &         ! condensation in updraft
               +evp(k,i)                         ! evaporating rain in downdraft

         dsdt(k,i) = -rl/cp*emc &
                     + (+mu(k+1,i)* (su(k+1,i)-shat(k+1,i)) &
                        -mu(k,i)*   (su(k,i)-shat(k,i)) &
                        +md(k+1,i)* (sd(k+1,i)-shat(k+1,i)) &
                        -md(k,i)*   (sd(k,i)-shat(k,i)) &
                       )/dp(k,i)

         dqdt(k,i) = emc + &
                    (+mu(k+1,i)* (qu(k+1,i)-qhat(k+1,i)) &
                     -mu(k,i)*   (qu(k,i)-qhat(k,i)) &
                     +md(k+1,i)* (qd(k+1,i)-qhat(k+1,i)) &
                     -md(k,i)*   (qd(k,i)-qhat(k,i)) &
                    )/dp(k,i)

         dl(k,i) = du(k,i)*ql(k+1,i)
      end do
   end do

   do k = kbm,pver
      do i = il1g,il2g
         if (k == mx(i)) then
            dsdt(k,i) = (1._r8/dsubcld(i))* &
                        (-mu(k,i)* (su(k,i)-shat(k,i)) &
                         -md(k,i)* (sd(k,i)-shat(k,i)) &
                        )
            dqdt(k,i) = (1._r8/dsubcld(i))* &
                        (-mu(k,i)*(qu(k,i)-qhat(k,i)) &
                         -md(k,i)*(qd(k,i)-qhat(k,i)) &
                        )
         else if (k > mx(i)) then
            dsdt(k,i) = dsdt(k-1,i)
            dqdt(k,i) = dqdt(k-1,i)
         end if
      end do
   end do

   return
   end subroutine q1q2_pjr


! Purpose: 
! Calculates CAPE the lifting condensation level and the convective top
! where buoyancy is first -ve.
! 
! Method: Calculates the parcel temperature based on a simple constant
! entraining plume model. CAPE is integrated from buoyancy.
! 09/09/04 - Simplest approach using an assumed entrainment rate for 
!            testing (dmpdp). 
! 08/04/05 - Swap to convert dmpdz to dmpdp  
!
! SCAM Logical Switches - DILUTE:RBN - Now Disabled 
! ---------------------
! switch(1) = .T. - Uses the dilute parcel calculation to obtain tendencies.
! switch(2) = .T. - Includes entropy/q changes due to condensate loss and freezing.
! switch(3) = .T. - Adds the PBL Tpert for the parcel temperature at all levels.
! 
! References:
! Raymond and Blythe (1992) JAS 
! 
! Author:
! Richard Neale - September 2004
    subroutine buoyan_dilute(q       ,t       ,p       ,z       ,pf      , &
                             tp      ,qstp    ,tl      ,rl      ,cape    , &
                             pblt    ,lcl     ,lel     ,lon     ,mx      , &
                             rd      ,grav    ,cp      ,msg     ,tpert   )
 
    implicit none
! io
    integer, intent(in)   :: msg
    real(r8), intent(in)  :: cp
    real(r8), intent(in)  :: rl
    real(r8), intent(in)  :: rd
    real(r8), intent(in)  :: grav
    real(r8), intent(in)  :: q(pver,ncol)         ! spec. humidity
    real(r8), intent(in)  :: t(pver,ncol)         ! temperature
    real(r8), intent(in)  :: p(pver,ncol)         ! pressure
    real(r8), intent(in)  :: z(pver,ncol)         ! height
    real(r8), intent(in)  :: pf(pver+1,ncol)      ! pressure at interfaces
    real(r8), intent(in)  :: pblt(ncol)           ! index of pbl depth
    real(r8), intent(in)  :: tpert(ncol)          ! perturbation temperature by pbl processes
 
    real(r8), intent(out) :: tp(pver,ncol)        ! parcel temperature
    real(r8), intent(out) :: qstp(pver,ncol)      ! saturation mixing ratio of parcel (only above lcl, just q below).
    real(r8), intent(out) :: tl(ncol)             ! parcel temperature at lcl
    real(r8), intent(out) :: cape(ncol)           ! convective aval. pot. energy.
    integer, intent(out)  :: lcl(ncol)        !
    integer, intent(out)  :: lel(ncol)        !
    integer, intent(out)  :: lon(ncol)        ! level of onset of deep convection
    integer, intent(out)  :: mx(ncol)         ! level of max moist static energy
 ! local
    real(r8) capeten(5,ncol)     ! provisional value of cape
    real(r8) tv(pver,ncol)       !
    real(r8) tpv(pver,ncol)      !
    real(r8) buoy(pver,ncol)
    real(r8) a1(ncol)
    real(r8) a2(ncol)
    real(r8) estp(ncol)
    real(r8) pl(ncol)
    real(r8) plexp(ncol)
    real(r8) hmax(ncol)
    real(r8) hmn(ncol)
    real(r8) y(ncol)
    logical plge600(ncol)
    integer knt(ncol)
    integer lelten(5,ncol)
    real(r8) e
    integer i
    integer k
    integer n
 
    do n = 1,5
       do i = 1,ncol
          lelten(n,i) = pver
          capeten(n,i) = 0._r8
       end do
    end do
 
    do i = 1,ncol
       lon(i) = pver
       knt(i) = 0
       lel(i) = pver
       mx(i) = lon(i)
       cape(i) = 0._r8
       hmax(i) = 0._r8
    end do
 
    tp(:,:ncol) = t(:,:ncol)
    qstp(:,:ncol) = q(:,:ncol)
 
!!! RBN - Initialize tv and buoy for output.
!!! tv=tv : tpv=tpv : qstp=q : buoy=0.
    tv(:,:ncol) = t(:,:ncol) *(1._r8+1.608_r8*q(:,:ncol))/ (1._r8+q(:,:ncol))
    tpv(:,:ncol) = tv(:,:ncol)
    buoy(:,:ncol) = 0._r8
 
! set "launching" level(mx) to be at maximum moist static energy.
! search for this level stops at planetary boundary layer top.
    do k = pver,msg + 1,-1
       do i = 1,ncol
          hmn(i) = cp*t(k,i) + grav*z(k,i) + rl*q(k,i)
          if (k >= nint(pblt(i)) .and. k <= lon(i) .and. hmn(i) > hmax(i)) then
             hmax(i) = hmn(i)
             mx(i) = k
          end if
       end do
    end do
 
 ! LCL dilute calculation - initialize to mx(i)
 ! Determine lcl in parcel_dilute and get pl,tl after parcel_dilute
 ! Original code actually sets LCL as level above wher condensate forms.
 ! Therefore in parcel_dilute lcl(i) will be at first level where qsmix < qtmix.
 
    do i = 1,ncol ! Initialise LCL variables.
       lcl(i) = mx(i)
       tl(i) = t(mx(i),i)
       pl(i) = p(mx(i),i)
    end do
 
! main buoyancy calculation.
!!! DILUTE PLUME CALCULATION USING ENTRAINING PLUME !!!
!!!   RBN 9/9/04   !!!
 
    call parcel_dilute(ncol, msg, mx, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)
 
 
! If lcl is above the nominal level of non-divergence (600 mbs),
! no deep convection is permitted (ensuing calculations
! skipped and cape retains initialized value of zero).
    do i = 1,ncol
       plge600(i) = pl(i).ge.600._r8 ! Just change to always allow buoy calculation.
    end do
 
! Main buoyancy calculation.
    do k = pver,msg + 1,-1
       do i=1,ncol
          if (k <= mx(i) .and. plge600(i)) then   ! Define buoy from launch level to cloud top.
             tv(k,i) = t(k,i)* (1._r8+1.608_r8*q(k,i))/ (1._r8+q(k,i))
             buoy(k,i) = tpv(k,i) - tv(k,i) + tiedke_add  ! +0.5K or not?
          else
             qstp(k,i) = q(k,i)
             tp(k,i)   = t(k,i)            
             tpv(k,i)  = tv(k,i)
          endif
       end do
    end do
 
    do k = msg + 2,pver
       do i = 1,ncol
          if (k < lcl(i) .and. plge600(i)) then
             if (buoy(k+1,i) > 0._r8 .and. buoy(k,i) <= 0._r8) then
                knt(i) = min(5,knt(i) + 1)
                lelten(knt(i),i) = k
             end if
          end if
       end do
    end do
! calculate convective available potential energy (cape).
    do n = 1,5
       do k = msg + 1,pver
          do i = 1,ncol
             if (plge600(i) .and. k <= mx(i) .and. k > lelten(n,i)) then
                capeten(n,i) = capeten(n,i) + rd*buoy(k,i)*log(pf(k+1,i)/pf(k,i))
             end if
          end do
       end do
    end do
! find maximum cape from all possible tentative capes from
! one sounding,
! and use it as the final cape, april 26, 1995
    do n = 1,5
       do i = 1,ncol
          if (capeten(n,i) > cape(i)) then
             cape(i) = capeten(n,i)
             lel(i) = lelten(n,i)
          end if
       end do
    end do
 !
 ! put lower bound on cape for diagnostic purposes.
 !
    do i = 1,ncol
       cape(i) = max(cape(i), 0._r8)
    end do
 !
    return
    end subroutine buoyan_dilute


    subroutine parcel_dilute (ncol, msg, klaunch, p, t, q, tpert, tp, tpv, qstp, pl, tl, lcl)

! Routine  to determine 
!   1. Tp   - Parcel temperature
!   2. qstp - Saturated mixing ratio at the parcel temperature.

    implicit none
    
    integer, intent(in) :: ncol
    integer, intent(in) :: msg
    integer, intent(in), dimension(ncol) :: klaunch(ncol)
    
    real(r8), intent(in), dimension(pver,ncol) :: p
    real(r8), intent(in), dimension(pver,ncol) :: t
    real(r8), intent(in), dimension(pver,ncol) :: q
    real(r8), intent(in), dimension(ncol) :: tpert ! PBL temperature perturbation.
    
    real(r8), intent(inout), dimension(pver,ncol) :: tp    ! Parcel temp.
    real(r8), intent(inout), dimension(pver,ncol) :: qstp  ! Parcel water vapour (sat value above lcl).
    real(r8), intent(inout), dimension(ncol) :: tl         ! Actual temp of LCL.
    real(r8), intent(inout), dimension(ncol) :: pl          ! Actual pressure of LCL. 
    
    integer, intent(inout), dimension(ncol) :: lcl ! Lifting condesation level (first model level with saturation).
    
    real(r8), intent(out), dimension(pver,ncol) :: tpv   ! Define tpv within this routine.


! Have to be careful as s is also dry static energy.

! If we are to retain the fact that CAM loops over grid-points in the internal
! loop then we need to dimension sp,atp,mp,xsh2o with ncol.

! local
    real(r8) tmix(pver,ncol)        ! Tempertaure of the entraining parcel.
    real(r8) qtmix(pver,ncol)       ! Total water of the entraining parcel.
    real(r8) qsmix(pver,ncol)       ! Saturated mixing ratio at the tmix.
    real(r8) smix(pver,ncol)        ! Entropy of the entraining parcel.
    real(r8) xsh2o(pver,ncol)       ! Precipitate lost from parcel.
    real(r8) ds_xsh2o(pver,ncol)    ! Entropy change due to loss of condensate.
    real(r8) ds_freeze(pver,ncol)   ! Entropy change sue to freezing of precip.
    
    real(r8) mp(ncol)    ! Parcel mass flux.
    real(r8) qtp(ncol)   ! Parcel total water.
    real(r8) sp(ncol)    ! Parcel entropy.
    real(r8) sp0(ncol)    ! Parcel launch entropy.
    real(r8) qtp0(ncol)   ! Parcel launch total water.
    real(r8) mp0(ncol)    ! Parcel launch relative mass flux.
    
    real(r8) lwmax      ! Maximum condesate that can be held in cloud before rainout.
    real(r8) dmpdp      ! Parcel fractional mass entrainment rate (/mb).
    !real(r8) dmpdpc     ! In cloud parcel mass entrainment rate (/mb).
    real(r8) dmpdz      ! Parcel fractional mass entrainment rate (/m)
    real(r8) dpdz,dzdp  ! Hydrstatic relation and inverse of.
    real(r8) senv       ! Environmental entropy at each grid point.
    real(r8) qtenv      ! Environmental total water "   "   ".
    real(r8) penv       ! Environmental total pressure "   "   ".
    real(r8) tenv       ! Environmental total temperature "   "   ".
    real(r8) new_s      ! Hold value for entropy after condensation/freezing adjustments.
    real(r8) new_q      ! Hold value for total water after condensation/freezing adjustments.
    real(r8) dp         ! Layer thickness (center to center)
    real(r8) tfguess    ! First guess for entropy inversion - crucial for efficiency!
    real(r8) tscool     ! Super cooled temperature offset (in degC) (eg -35).
    
    real(r8) qxsk, qxskp1        ! LCL excess water (k, k+1)
    real(r8) dsdp, dqtdp, dqxsdp ! LCL s, qt, p gradients (k, k+1)
    real(r8) slcl,qtlcl,qslcl    ! LCL s, qt, qs values.
    
    integer rcall       ! Number of ientropy call for errors recording
    integer nit_lheat     ! Number of iterations for condensation/freezing loop.
    integer i,k,ii   ! Loop counters.

!======================================================================
!    SUMMARY
!
!  9/9/04 - Assumes parcel is initiated from level of maxh (klaunch)
!           and entrains at each level with a specified entrainment rate.
!
! 15/9/04 - Calculates lcl(i) based on k where qsmix is first < qtmix.          
!
!======================================================================
!
! Set some values that may be changed frequently.
!

    nit_lheat = 2 ! iterations for ds,dq changes from condensation freezing.
    dmpdz=-1.e-3_r8        ! Entrainment rate. (-ve for /m), -2e-3 - -0.2e-3, default: -1e-3
    !dmpdpc = 3.e-2_r8   ! In cloud entrainment rate (/mb).
    lwmax = 1.e-3_r8    ! Need to put formula in for this.
    tscool = 0.0_r8   ! Temp at which water loading freezes in the cloud.
    
    qtmix=0._r8
    smix=0._r8
    
    qtenv = 0._r8
    senv = 0._r8
    tenv = 0._r8
    penv = 0._r8
    
    qtp0 = 0._r8
    sp0  = 0._r8
    mp0 = 0._r8
    
    qtp = 0._r8
    sp = 0._r8
    mp = 0._r8
    
    new_q = 0._r8
    new_s = 0._r8
    
    do k = pver, msg+1, -1
       do i=1,ncol 
    
    ! Initialize parcel values at launch level.
    
          if (k == klaunch(i)) then 
             qtp0(i) = q(k,i)   ! Parcel launch total water (assuming subsaturated) - OK????.
             sp0(i)  = entropy(t(k,i),p(k,i),qtp0(i))  ! Parcel launch entropy.
             mp0(i)  = 1._r8       ! Parcel launch relative mass (i.e. 1 parcel stays 1 parcel for dmpdp=0, undilute). 
             smix(k,i)  = sp0(i)
             qtmix(k,i) = qtp0(i)
             tfguess = t(k,i)
             rcall = 1
             call ientropy (rcall,i,smix(k,i),p(k,i),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)
          end if
    
    ! Entraining levels
          
          if (k < klaunch(i)) then 
    
    ! Set environmental values for this level.                 
             
             dp = (p(k,i)-p(k+1,i)) ! In -ve mb as p decreasing with height - difference between center of layers.
             qtenv = 0.5_r8*(q(k,i)+q(k+1,i))         ! Total water of environment.
             tenv  = 0.5_r8*(t(k,i)+t(k+1,i)) 
             penv  = 0.5_r8*(p(k,i)+p(k+1,i))
    
             senv  = entropy(tenv,penv,qtenv)  ! Entropy of environment.   
    
    ! Determine fractional entrainment rate /pa given value /m.
    
             dpdz = -(penv*grav)/(rgas*tenv) ! in mb/m since  p in mb.
             dzdp = 1._r8/dpdz                  ! in m/mb
             dmpdp = dmpdz*dzdp              ! /mb Fractional entrainment
    
    ! Sum entrainment to current level
    ! entrains q,s out of intervening dp layers, in which linear variation is assumed
    ! so really it entrains the mean of the 2 stored values.
    
             sp(i)  = sp(i)  - dmpdp*dp*senv 
             qtp(i) = qtp(i) - dmpdp*dp*qtenv 
             mp(i)  = mp(i)  - dmpdp*dp
                
    ! Entrain s and qt to next level.
    
             smix(k,i)  = (sp0(i)  +  sp(i)) / (mp0(i) + mp(i))

             qtmix(k,i) = (qtp0(i) + qtp(i)) / (mp0(i) + mp(i))
    
    ! Invert entropy from s and q to determine T and saturation-capped q of mixture.
    ! t(k,i) used as a first guess so that it converges faster.
    
             tfguess = tmix(k+1,i)
             rcall = 2
             call ientropy(rcall,i,smix(k,i),p(k,i),qtmix(k,i),tmix(k,i),qsmix(k,i),tfguess)   
    !
    ! Determine if this is lcl of this column if qsmix <= qtmix.
    ! FIRST LEVEL where this happens on ascending.
    
             if (qsmix(k,i) <= qtmix(k,i) .and. qsmix(k+1,i) > qtmix(k+1,i)) then
                lcl(i) = k
                qxsk   = qtmix(k,i) - qsmix(k,i)
                qxskp1 = qtmix(k+1,i) - qsmix(k+1,i)
                dqxsdp = (qxsk - qxskp1)/dp
                pl(i)  = p(k+1,i) - qxskp1/dqxsdp    ! pressure level of actual lcl.
                dsdp   = (smix(k,i)  - smix(k+1,i))/dp
                dqtdp  = (qtmix(k,i) - qtmix(k+1,i))/dp
                slcl   = smix(k+1,i)  +  dsdp* (pl(i)-p(k+1,i))  
                qtlcl  = qtmix(k+1,i) +  dqtdp*(pl(i)-p(k+1,i))
    
                tfguess = tmix(k,i)
                rcall = 3
                call ientropy (rcall,i,slcl,pl(i),qtlcl,tl(i),qslcl,tfguess)
    
    !            write(iulog,*)' '
    !            write(iulog,*)' p',p(k+1,i),pl(i),p(lcl(i),i)
    !            write(iulog,*)' t',tmix(k+1,i),tl(i),tmix(lcl(i),i)
    !            write(iulog,*)' s',smix(k+1,i),slcl,smix(lcl(i),i)
    !            write(iulog,*)'qt',qtmix(k+1,i),qtlcl,qtmix(lcl(i),i)
    !            write(iulog,*)'qs',qsmix(k+1,i),qslcl,qsmix(lcl(i),i)
    
             endif
          end if !  k < klaunch
     
       end do ! Levels loop
    end do ! Columns loop

!!!!!!!!!!!!!!!!!!!!!!!!!!END ENTRAINMENT LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! Could stop now and test with this as it will provide some estimate of buoyancy
!! without the effects of freezing/condensation taken into account for tmix.

!! So we now have a profile of entropy and total water of the entraining parcel
!! Varying with height from the launch level klaunch parcel=environment. To the 
!! top allowed level for the existence of convection.

!! Now we have to adjust these values such that the water held in vaopor is < or 
!! = to qsmix. Therefore, we assume that the cloud holds a certain amount of
!! condensate (lwmax) and the rest is rained out (xsh2o). This, obviously 
!! provides latent heating to the mixed parcel and so this has to be added back 
!! to it. But does this also increase qsmix as well? Also freezing processes
 

    xsh2o = 0._r8
    ds_xsh2o = 0._r8
    ds_freeze = 0._r8

!!!!!!!!!!!!!!!!!!!!!!!!!PRECIPITATION/FREEZING LOOP!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Iterate solution twice for accuracy



    do k = pver, msg+1, -1
       do i=1,ncol    
          
    ! Initialize variables at k=klaunch
          
          if (k == klaunch(i)) then
    
    ! Set parcel values at launch level assume no liquid water.            
    
             tp(k,i)    = tmix(k,i)
             qstp(k,i)  = q(k,i) 
             tpv(k,i)   =  (tp(k,i) + tpert(i)) * (1._r8+1.608_r8*qstp(k,i)) / (1._r8+qstp(k,i))
             
          end if
    
          if (k < klaunch(i)) then
                
    ! Initiaite loop if switch(2) = .T. - RBN:DILUTE - TAKEN OUT BUT COULD BE RETURNED LATER.
    
    ! Iterate nit_lheat times for s,qt changes.
    
             do ii=0,nit_lheat-1            
    
    ! Rain (xsh2o) is excess condensate, bar LWMAX (Accumulated loss from qtmix).
    
                xsh2o(k,i) = max (0._r8, qtmix(k,i) - qsmix(k,i) - lwmax)
    
    ! Contribution to ds from precip loss of condensate (Accumulated change from smix).(-ve)                     
                         
                ds_xsh2o(k,i) = ds_xsh2o(k+1,i) - cpliq * log (tmix(k,i)/tfreez) * max(0._r8,(xsh2o(k,i)-xsh2o(k+1,i)))
    !
    ! Entropy of freezing: latice times amount of water involved divided by T.
    !
     
                if (tmix(k,i) <= tfreez+tscool .and. ds_freeze(k+1,i) == 0._r8) then ! One off freezing of condensate. 
                   ds_freeze(k,i) = (latice/tmix(k,i)) * max(0._r8,qtmix(k,i)-qsmix(k,i)-xsh2o(k,i)) ! Gain of LH
                end if
                
                if (tmix(k,i) <= tfreez+tscool .and. ds_freeze(k+1,i) /= 0._r8) then ! Continual freezing of additional condensate.
                   ds_freeze(k,i) = ds_freeze(k+1,i)+(latice/tmix(k,i)) * max(0._r8,(qsmix(k+1,i)-qsmix(k,i)))
                end if
                
    ! Adjust entropy and accordingly to sum of ds (be careful of signs).
    
                new_s = smix(k,i) + ds_xsh2o(k,i) + ds_freeze(k,i) 
    
    ! Adjust liquid water and accordingly to xsh2o.
    
                new_q = qtmix(k,i) - xsh2o(k,i)
    
    ! Invert entropy to get updated Tmix and qsmix of parcel.
    
                tfguess = tmix(k,i)
                rcall =4
                call ientropy (rcall,i,new_s, p(k,i), new_q, tmix(k,i), qsmix(k,i), tfguess)
                
             end do  ! Iteration loop for freezing processes.
    
    ! tp  - Parcel temp is temp of mixture.
    ! tpv - Parcel v. temp should be density temp with new_q total water. 
    
             tp(k,i)    = tmix(k,i)
    
    ! tpv = tprho in the presence of condensate (i.e. when new_q > qsmix)
    
             if (new_q > qsmix(k,i)) then  ! Super-saturated so condensate present - reduces buoyancy.
                qstp(k,i) = qsmix(k,i)
             else                          ! Just saturated/sub-saturated - no condensate virtual effects.
                qstp(k,i) = new_q
             end if
    
             tpv(k,i) = (tp(k,i)+tpert(i))* (1._r8+1.608_r8*qstp(k,i)) / (1._r8+ new_q) 
    
          end if ! k < klaunch
          
       end do ! Loop for columns
       
    end do  ! Loop for vertical levels.
    
    
    return
    end subroutine parcel_dilute


! TK(K),p(mb),qtot(kg/kg)
! from Raymond and Blyth 1992
    real(r8) function entropy(TK,p,qtot)
    real(r8), intent(in) :: p,qtot,TK
    real(r8) :: qv,qst,e,est,L,eref,pref

    pref = 1000.0_r8           ! mb
    eref = 6.106_r8            ! sat p at tfreez (mb)
    
    L = rl - (cpliq - cpwv)*(TK-tfreez)         ! T IN CENTIGRADE
    
    ! Replace call to satmixutils.
    
    call qmmr_hPa(TK, p, est, qst)
    
    qv = min(qtot,qst)                         ! Partition qtot into vapor part only.
    e = qv*p / (eps1 +qv)
    
    entropy = (cpres + qtot*cpliq)*log( TK/tfreez) - rgas*log( (p-e)/pref ) + &
            L*qv/TK - qv*rh2o*log(qv/qst)
     
    return
    end FUNCTION entropy


! p(mb), Tfg/T(K), qt/qv(kg/kg), s(J/kg). 
! Inverts entropy, pressure and total water qt 
! for T and saturated vapor mixing ratio
    subroutine ientropy (rcall,icol,s,p,qt,T,qst,Tfg)

     integer, intent(in) :: icol, rcall
     real(r8), intent(in)  :: s, p, Tfg, qt
     real(r8), intent(out) :: qst, T
     real(r8) :: qv,Ts,dTs,fs1,fs2,est
     real(r8) :: pref,eref,L,e
     real(r8) :: this_lat,this_lon
     integer :: LOOPMAX,i


    LOOPMAX = 100                   !* max number of iteration loops 
    
    ! Values for entropy
    pref = 1000.0_r8           ! mb ref pressure.
    eref = 6.106_r8           ! sat p at tfreez (mb)
    
    ! Invert the entropy equation -- use Newton's method
    
    Ts = Tfg                  ! Better first guess based on Tprofile from conv.
    
    converge: do i=0, LOOPMAX
    
       L = rl - (cpliq - cpwv)*(Ts-tfreez) 
    
       call qmmr_hPa(Ts, p, est, qst)
       qv = min(qt,qst) 
       e = qv*p / (eps1 +qv)  ! Bolton (eq. 16)
       fs1 = (cpres + qt*cpliq)*log( Ts/tfreez ) - rgas*log( (p-e)/pref ) + &
            L*qv/Ts - qv*rh2o*log(qv/qst) - s
       
       L = rl - (cpliq - cpwv)*(Ts-1._r8-tfreez)         
    
       call qmmr_hPa(Ts-1._r8, p, est, qst)
       qv = min(qt,qst) 
       e = qv*p / (eps1 +qv)
       fs2 = (cpres + qt*cpliq)*log( (Ts-1._r8)/tfreez ) - rgas*log( (p-e)/pref ) + &
            L*qv/(Ts-1._r8) - qv*rh2o*log(qv/qst) - s 
       
       dTs = fs1/(fs2 - fs1)
       Ts  = Ts+dTs
       if (abs(dTs).lt.0.001_r8) exit converge
       if (i .eq. LOOPMAX - 1) then
          !this_lat = get_rlat_p(lchnk, icol)*57.296_r8
          !this_lon = get_rlon_p(lchnk, icol)*57.296_r8
          if(mpi_rank()==0)then
          print*, '*** ZM_CONV: IENTROPY: Failed and about to exit, info follows ****'
          print*,'rank = ',mpi_rank()
          print*, 'ZM_CONV: IENTROPY. Details: rcall,icol= ',rcall,icol, &
          ! ' lat: ',this_lat,' lon: ',this_lon, &
           ' P(mb)= ', p, ' Tfg(K)= ', Tfg, 'T(K)= ', T,                 &
           ' qt(g/kg) = ', 1000._r8*qt, ' qst(g/kg) = ', 1000._r8*qst,', s(J/kg) = ',s
          call endrun("**** ZM_CONV IENTROPY: Tmix did not converge ****")
          end if
       end if
    enddo converge
    
    ! Replace call to satmixutils.
    
    call qmmr_hPa(Ts, p, est, qst)
    
    qv = min(qt,qst)                             !       /* check for saturation */
    T = Ts 
    
!     100    format (A,I1,I4,I4,7(A,F6.2))
    
    return
    end SUBROUTINE ientropy


    elemental subroutine qmmr_hPa(t, p, es, qm)
    use grist_wv_saturation, only: qmmr

    ! Inputs
    real(r8), intent(in) :: t    ! Temperature (K)
    real(r8), intent(in) :: p    ! Pressure (hPa)
    ! Outputs
    real(r8), intent(out) :: es  ! Saturation vapor pressure (hPa)
    real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                                  ! (vapor mass over dry mass, kg/kg)

    call qmmr(t, p*100._r8, es, qm)

    es = es*0.01_r8

    end subroutine qmmr_hPa


 end module grist_zm_conv
