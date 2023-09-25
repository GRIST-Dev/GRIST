
!----------------------------------------------------------------------------
! Created on 2019
! Author: Yi Zhang
! Version 1.0
! Description: This driver mainly follows the naming convention of WRF, but 
!              the dim conversion requires a transpose from GRIST based vars
! Revision history:
!            1. YSUV381, CAMUW_V381 as two options
!--------------------------------------------------------------------------

module grist_PhysW_pbl_driver

  use grist_constants,          only: i4, r8, r_d=>rdry, r_v=>rvap, cp, gravity, zero, pi
  use grist_mpi
  use grist_nml_module,           only: doAquaPlanet
  use grist_PhysW_nml_module,     only: PhysW_bl_scheme, use_gwdo, step_bl
  use grist_PhysW_data_structure, only: p_qv, p_qc, p_qr, p_qi, p_qs, p_qg, param_first_scalar
  use grist_PhysW_data_structure, only: pstate_wrf, ptend_wrf, psurf_wrf
! module 
  use module_wrfmodel_constants,only: rcp, rovg, xls, xlv, xlf, &
                                      ep_1, ep_2, svp1, svp2, svp3, svpt0, &
                                      karman, eomeg, stbolt
  use module_bl_ysu_v381,     only: ysu_v381=>ysu
  use module_bl_gwdo_v381,    only: gwdo_v381=>gwdo
#ifdef UW_CAM5
  use module_bl_camuwpbl_driver_v381, only: camuwpblinit, camuwpbl
#endif

  implicit  none

  private
  public  ::  grist_PhysW_pbl_init, &
              grist_PhysW_pbl_run,  &
              grist_PhysW_pbl_final

   integer(i4), public   :: stepbl
   logical               :: is_CAMMGMP_used = .false.
!   real(r8),allocatable :: kvm(:,:,:)
!   real(r8),allocatable :: kvh(:,:,:)
!   real(r8),allocatable :: tauresx2d(:,:)
!   real(r8),allocatable :: tauresy2d(:,:)

contains

    subroutine grist_PhysW_pbl_init(ncell, nLevel)
      integer(i4),  intent(in) :: ncell, nLevel
  
      stepbl = step_bl

      if(.not.allocated(ptend_wrf%ruublten)) allocate(ptend_wrf%ruublten(1:ncell, 1:nLevel, 1:1))
      if(.not.allocated(ptend_wrf%rvvblten)) allocate(ptend_wrf%rvvblten(1:ncell, 1:nLevel, 1:1))
      if(.not.allocated(ptend_wrf%rthblten)) allocate(ptend_wrf%rthblten(1:ncell, 1:nLevel, 1:1))
      if(.not.allocated(ptend_wrf%rqvblten)) allocate(ptend_wrf%rqvblten(1:ncell, 1:nLevel, 1:1))
      if(.not.allocated(ptend_wrf%rqcblten)) allocate(ptend_wrf%rqcblten(1:ncell, 1:nLevel, 1:1))
      if(.not.allocated(ptend_wrf%rqiblten)) allocate(ptend_wrf%rqiblten(1:ncell, 1:nLevel, 1:1))
!
! not using "init" in their module
!
      ptend_wrf%ruublten = zero
      ptend_wrf%rvvblten = zero
      ptend_wrf%rthblten = zero
      ptend_wrf%rqvblten = zero
      ptend_wrf%rqcblten = zero
      ptend_wrf%rqiblten = zero
#ifdef UW_CAM5
!For CAMUW Scheme, LiXH add.
     if(trim(PhysW_bl_scheme) .eq. 'CAMUW_V381')then
!        if(.not.allocated(kvm)) allocate(kvm(1:ncell, 1:nLevel, 1:1))
!        if(.not.allocated(kvh)) allocate(kvh(1:ncell, 1:nLevel, 1:1))
!        if(.not.allocated(tauresx2d)) allocate(tauresx2d(1:ncell, 1:1))
!        if(.not.allocated(tauresy2d)) allocate(tauresy2d(1:ncell, 1:1))
!        kvm = 0._r8
!        kvh = 0._r8
!        tauresx2d = 0._r8
!        tauresy2d = 0._r8

        call camuwpblinit(is_CAMMGMP_used=is_CAMMGMP_used,          &
                   ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,  &
                   ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
                   its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )

     end if
#endif

      return
    end subroutine grist_PhysW_pbl_init

    subroutine grist_PhysW_pbl_final()
      if(allocated(ptend_wrf%ruublten)) deallocate(ptend_wrf%ruublten)
      if(allocated(ptend_wrf%rvvblten)) deallocate(ptend_wrf%rvvblten)
      if(allocated(ptend_wrf%rthblten)) deallocate(ptend_wrf%rthblten)
      if(allocated(ptend_wrf%rqvblten)) deallocate(ptend_wrf%rqvblten)
      if(allocated(ptend_wrf%rqcblten)) deallocate(ptend_wrf%rqcblten)
      if(allocated(ptend_wrf%rqiblten)) deallocate(ptend_wrf%rqiblten)

!For CAMUW Scheme, LiXH add.
!      if(allocated(kvm)) deallocate(kvm)
!      if(allocated(kvh)) deallocate(kvh)
!      if(allocated(tauresx2d)) deallocate(tauresx2d)
!      if(allocated(tauresy2d)) deallocate(tauresy2d)

      return
    end subroutine grist_PhysW_pbl_final
 
    subroutine grist_PhysW_pbl_run(ncell,nLevel,nspecies,itimestep,dtime)

    implicit none
!======================================================================
! grid structure in physics part of wrf
!----------------------------------------------------------------------
! the horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. all predicted
! variables are carried at half levels except w, which is at full
! levels. some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! in wrf, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  in your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         .
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! definitions
!-----------
! rho_d      dry density (kg/m^3)
! theta_m    moist potential temperature (k)
! qv         water vapor mixing ratio (kg/kg)
! qc         cloud water mixing ratio (kg/kg)
! qr         rain water mixing ratio (kg/kg)
! qi         cloud ice mixing ratio (kg/kg)
! qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- rublten       u tendency due to 
!                 pbl parameterization (m/s^2)
!-- rvblten       v tendency due to 
!                 pbl parameterization (m/s^2)
!-- rthblten      theta tendency due to 
!                 pbl parameterization (k/s)
!-- rqvblten      qv tendency due to 
!                 pbl parameterization (kg/kg/s)
!-- rqcblten      qc tendency due to 
!                 pbl parameterization (kg/kg/s)
!-- rqiblten      qi tendency due to 
!                 pbl parameterization (kg/kg/s)
!-- itimestep     number of time steps
!-- glw           downward long wave flux at ground surface (w/m^2)
!-- gsw           downward short wave flux at ground surface (w/m^2)
!-- emiss         surface emissivity (between 0 and 1)
!-- tsk           surface temperature (k)
!-- tmn           soil temperature at lower boundary (k)
!-- xland         land mask (1 for land, 2 for water)
!-- znt           roughness length (m)
!-- mavail        surface moisture availability (between 0 and 1)
!-- ust           u* in similarity theory (m/s)
!-- mol           q* (similarity theory) (kg/kg)
!-- hol           pbl height over monin-obukhov length
!-- pblh          pbl height (m)
!-- capg          heat capacity for soil (j/k/m^3)
!-- thc           thermal inertia (cal/cm/k/s^0.5)
!-- snowc         flag indicating snow coverage (1 for snow cover)
!-- hfx           upward heat flux at the surface (w/m^2)
!-- qfx           upward moisture flux at the surface (kg/m^2/s)
!-- regime        flag indicating pbl regime (stable, unstable, etc.)
!-- tke_myj       turbulence kinetic energy from mellor-yamada-janjic (myj) (m^2/s^2)
!-- akhs          sfc exchange coefficient of heat/moisture from myj
!-- akms          sfc exchange coefficient of momentum from myj
!-- thz0          potential temperature at roughness length (k)
!-- uz0           u wind component at roughness length (m/s)
!-- vz0           v wind component at roughness length (m/s)
!-- qsfc          specific humidity at lower boundary (kg/kg)
!-- th2           diagnostic 2-m theta from surface layer and lsm
!-- t2            diagnostic 2-m temperature from surface layer and lsm
!-- q2            diagnostic 2-m mixing ratio from surface layer and lsm
!-- lowlyr        index of lowest model layer above ground
!-- rr            dry air density (kg/m^3)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (k)
!-- moist         moisture array (4d - last index is species) (kg/kg)
!-- p_phy         pressure (pa)
!-- pi_phy        exner function (dimensionless)
!-- p8w           pressure at full levels (pa)
!-- t_phy         temperature (k)
!-- dz8w          dz between full levels (m)
!-- z             height above sea level (m)
!-- config_flags
!-- dx            horizontal space interval (m)
!-- dt            time step (second)
!-- nspecies       number of moisture species
!-- psfc          pressure at the surface (pa)
!-- tslb          
!-- zs
!-- dzs
!-- num_soil_layers number of soil layer
!-- ifsnow      ifsnow=1 for snow-cover effects
!
!-- p_qv          species index for water vapor
!-- p_qc          species index for cloud water
!-- p_qr          species index for rain water
!-- p_qi          species index for cloud ice
!-- p_qs          species index for snow
!-- p_qg          species index for graupel
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!
!io
    integer(i4), intent(in)    :: ncell
    integer(i4), intent(in)    :: nLevel
    integer(i4), intent(in)    :: nspecies
    integer(i4), intent(in)    :: itimestep
    real(r8),    intent(in)    :: dtime
! local
!   real(r8)     :: xland  (1:ncell, 1:1)
!   real(r8)     :: psim   (1:ncell, 1:1)
!   real(r8)     :: psih   (1:ncell, 1:1)
!   real(r8)     :: gz1oz0 (1:ncell, 1:1)
!   real(r8)     :: br     (1:ncell, 1:1)
!   real(r8)     :: tsk    (1:ncell, 1:1)
!   real(r8)     :: ust    (1:ncell, 1:1)
!   real(r8)     :: mol    (1:ncell, 1:1)
!   real(r8)     :: pblh   (1:ncell, 1:1)
!   real(r8)     :: hfx    (1:ncell, 1:1)
!   real(r8)     :: qfx    (1:ncell, 1:1)
!   real(r8)     :: regime (1:ncell, 1:1)
!   real(r8)     :: znt    (1:ncell, 1:1)
!   real(r8)     :: wspd   (1:ncell, 1:1)
! local
!   real(r8)     :: u_frame
!   real(r8)     :: v_frame
   real(r8)     :: hol    (1:ncell, 1:1)
   real(r8)     :: wstar  (1:ncell, 1:1)
   real(r8)     :: delta  (1:ncell, 1:1)
   real(r8)     :: uoce   (1:ncell, 1:1), voce(1:ncell,1:1)
   real(r8)     :: exch_h (1:ncell, 1:nLevel, 1:1)
   real(r8)     :: dtaux3d(1:ncell, 1:nLevel, 1:1)
   real(r8)     :: dtauy3d(1:ncell, 1:nLevel, 1:1)
   real(r8)     :: dusfcg (1:ncell, 1:1)
   real(r8)     :: dvsfcg (1:ncell, 1:1)

! SSW add for Shin-Hong scheme
!   real(r8)     :: tke_pbl(1:ncell, 1:nLevel, 1:1)
!   real(r8)     :: el_pbl (1:ncell, 1:nLevel, 1:1)

!   real(r8)     :: v_phytmp( 1:ncell, 1:nlev, 1:1 )
!   real(r8)     :: u_phytmp( 1:ncell, 1:nlev, 1:1 )
!   real(r8)     :: tskold  ( 1:ncell,         1:1 )
!   real(r8)     :: ustold  ( 1:ncell,         1:1 )
!   real(r8)     :: zntold  ( 1:ncell,         1:1 )
!   real(r8)     :: zol     ( 1:ncell,         1:1 )
!   real(r8)     :: psfc    ( 1:ncell,         1:1 )
   real(r8)     :: dtmin
   real(r8)     :: dtbl
! For CAMUW, LiXH add
   real(r8)     :: tpert2d(1:ncell, 1:1)
   real(r8)     :: qpert2d(1:ncell, 1:1)
   real(r8)     :: wpert2d(1:ncell, 1:1)
   real(r8)     :: smaw3d(1:ncell, 1:nLevel, 1:1)
   real(r8)     :: turbtype3d(1:ncell, 1:nLevel, 1:1)

   dtmin= dtime/60._r8
! pbl schemes need pbl time step for updates
   dtbl = dtime*stepbl

IF(itimestep.eq.0 .or. itimestep .eq. 1 .or.  mod(itimestep,stepbl) .eq. 0) then

     select case(trim(PhysW_bl_scheme))

     case('YSUV381') ! from WRF-V3.8.1
! assume no-slip sea surface
     uoce = zero
     voce = zero

     call ysu_v381(u3d         = pstate_wrf%u_phy(  1:ncell,1:nLevel,1:1), &
                   v3d         = pstate_wrf%v_phy(  1:ncell,1:nLevel,1:1), &
                   th3d        = pstate_wrf%th_phy( 1:ncell,1:nLevel,1:1), &
                   t3d         = pstate_wrf%t_phy(  1:ncell,1:nLevel,1:1), &
                   qv3d        = pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qv), &
                   qc3d        = pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qc), &
                   qi3d        = pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qi), &
                   p3d         = pstate_wrf%p_phy(  1:ncell,1:nLevel,1:1), &
                   p3di        = pstate_wrf%p8w  (  1:ncell,1:nLevel,1:1), &
                   pi3d        = pstate_wrf%pi_phy( 1:ncell,1:nLevel,1:1), &
                   rublten     = ptend_wrf%ruublten(1:ncell,1:nLevel,1:1), &
                   rvblten     = ptend_wrf%rvvblten(1:ncell,1:nLevel,1:1), &
                   rthblten    = ptend_wrf%rthblten(1:ncell,1:nLevel,1:1), &
                   rqvblten    = ptend_wrf%rqvblten(1:ncell,1:nLevel,1:1), &
                   rqcblten    = ptend_wrf%rqcblten(1:ncell,1:nLevel,1:1), &
                   rqiblten    = ptend_wrf%rqiblten(1:ncell,1:nLevel,1:1), &
                   flag_qi     = .false.,     &
                   cp          = cp,          &
                   g           = gravity,     &
                   rovcp       = rcp,         &
                   rd          = r_d,         &
                   rovg        = rovg,        &
                   ep1         = ep_1,        &
                   ep2         = ep_2,        &
                   karman      = karman,      &
                   xlv         = xlv,         &
                   rv          = r_v,         &
                   dz8w        = pstate_wrf%dz8w(1:ncell,1:nLevel,1:1), &
                   psfc        = psurf_wrf%psfc(1:ncell,1:1),           &
                   !znu,    ! not used
                   !znw,    ! --
                   !mut,    ! --
                   !p_top,  ! --
                   znt         = psurf_wrf%znt(   1:ncell,1:1), & ! inout, rst
                   ust         = psurf_wrf%ust(   1:ncell,1:1), & ! inout, rst
                   hpbl        = psurf_wrf%pblh(  1:ncell,1:1), & ! inout, rst
                   psim        = psurf_wrf%psim(  1:ncell,1:1), & ! in
                   psih        = psurf_wrf%psih(  1:ncell,1:1), & ! in
                   xland       = psurf_wrf%xland( 1:ncell,1:1), & ! in
                   hfx         = psurf_wrf%hfx(   1:ncell,1:1), & ! in
                   qfx         = psurf_wrf%qfx(   1:ncell,1:1), & ! in
                   wspd        = psurf_wrf%wspd(  1:ncell,1:1), & ! inout
                   br          = psurf_wrf%br(    1:ncell,1:1), & ! in
                   dt          = dtime,                         &
                   kpbl2d      = pstate_wrf%kpbl(1:ncell,1:1),  & ! out,   actually not used in YWU
                   exch_h      = exch_h(1:ncell,1:nLevel,1:1),  & ! out,   not used
                   wstar       = wstar(1:ncell,1:1),            & ! inout, 1st calc inside
                   delta       = delta(1:ncell,1:1),            & ! inout, 1st calc inside
                   u10         = psurf_wrf%u10(   1:ncell,1:1), & ! inout
                   v10         = psurf_wrf%v10(   1:ncell,1:1), & ! inout
                   uoce        = uoce(1:ncell,1:1),             &
                   voce        = voce(1:ncell,1:1),             &
                   rthraten    = ptend_wrf%rthraten(1:ncell,1:nLevel,1:1), & ! in
                   ysu_topdown_pblmix = 0,             & ! 1 is active, not using now
                   !ctopo,   ! not used
                   !ctopo2,  ! not used 
                   ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,  &
                   ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
                   its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )
              !optional, ! not used
                   !regime )

#ifdef UW_CAM5
        case ('CAMUW_V381')
          call camuwpbl(&
                   dt          = dtime,  &
                   u_phy       = pstate_wrf%u_phy(  1:ncell,1:nLevel,1:1), &
                   v_phy       = pstate_wrf%v_phy(  1:ncell,1:nLevel,1:1), &
                   th_phy      = pstate_wrf%th_phy( 1:ncell,1:nLevel,1:1), &
                   rho         = pstate_wrf%rhom ( 1:ncell,1:nLevel, 1:1), &
                   qv_curr     = pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qv), &
                   hfx         = psurf_wrf%hfx(   1:ncell,1:1), & ! in
                   qfx         = psurf_wrf%qfx(   1:ncell,1:1), & ! in
                   ustar       = psurf_wrf%ust(   1:ncell,1:1), & ! in 
                   p8w         = pstate_wrf%p8w  (  1:ncell,1:nLevel,1:1), &
                   p_phy       = pstate_wrf%p_phy(  1:ncell,1:nLevel,1:1), &
                   z           = pstate_wrf%zzz  (  1:ncell,1:nLevel,1:1), &
                   t_phy       = pstate_wrf%t_phy(  1:ncell,1:nLevel,1:1), &
                   qc_curr     = pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qc), &
                   qi_curr     = pstate_wrf%moist(  1:ncell,1:nLevel,1:1,p_qi), &
                   z_at_w      = pstate_wrf%z8w  (  1:ncell,1:nLevel,1:1), &
                   !CLDFRA_OLD_mp,    ! not used, CAMMG microp is not available.
                   CLDFRA      = pstate_wrf%cldstr(1:ncell,1:nLevel,1:1) , &  
                   ht          = psurf_wrf%ht(1:ncell,1:1),         &
                   rthratenlw  = ptend_wrf%rthraten_lw(1:ncell,1:nLevel, 1:1),  &
                   exner       = pstate_wrf%pi_phy( 1:ncell,1:nLevel,1:1), &
                   is_CAMMGMP_used=is_CAMMGMP_used,                 &
                   itimestep   = itimestep,                         &
                   !qnc_curr, qni_curr, wsedl3d,  !not used?
                   ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,  &
                   ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
                   its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1,&
                   tauresx2d   = pstate_wrf%tauresx2d(1:ncell,1:1),        &
                   tauresy2d   = pstate_wrf%tauresy2d(1:ncell,1:1),        &
                   rublten     = ptend_wrf%ruublten(1:ncell,1:nLevel,1:1), &
                   rvblten     = ptend_wrf%rvvblten(1:ncell,1:nLevel,1:1), &
                   rthblten    = ptend_wrf%rthblten(1:ncell,1:nLevel,1:1), &
                   rqiblten    = ptend_wrf%rqiblten(1:ncell,1:nLevel,1:1), &
                   !rqniblten,                    !not used?  WSM6 is Single Moment scheme
                   rqvblten    = ptend_wrf%rqvblten(1:ncell,1:nLevel,1:1), &
                   rqcblten    = ptend_wrf%rqcblten(1:ncell,1:nLevel,1:1), &
                   kvm3D       = pstate_wrf%kvm(1:ncell,1:nLevel,1:1),     &
                   kvh3D       = pstate_wrf%kvh(1:ncell,1:nLevel,1:1),     & 
                   tpert2d     = tpert2d(1:ncell,1:1),              &
                   qpert2d     = qpert2d(1:ncell,1:1),              &
                   wpert2d     = wpert2d(1:ncell,1:1),              &
                   smaw3d      = smaw3d(1:ncell,1:nLevel,1:1),      & 
                   turbtype3d  = turbtype3d(1:ncell,1:nLevel,1:1),  & 
                   tke_pbl     = pstate_wrf%tke(1:ncell, 1:nLevel, 1:1),   & 
                   leng_pbl    = pstate_wrf%leng(1:ncell, 1:nLevel, 1:1),  &
                   pblh2D      = psurf_wrf%pblh(  1:ncell,1:1),     &
                   kpbl2D      = pstate_wrf%kpbl(1:ncell,1:1)       )

#endif
        case default
          if(mpi_rank().eq.0) print*, "you must select a PBL scheme: YSUV381 (default), CAMUW_V381 (needs -DUW_CAM5 for compiling), stop"
          call mpi_abort()
      end select
!
! This should be used only for real case, not Aqua Planet
!
      if(use_gwdo .and. .not.doAquaPlanet)then
!
! GWD will add tendency to ruublten and rvvblten
!
      call gwdo_v381(&
                u3d     = pstate_wrf%u_phy(1:ncell,1:nLevel,1:1),      &
                v3d     = pstate_wrf%v_phy(1:ncell,1:nLevel,1:1),      &
                t3d     = pstate_wrf%t_phy(1:ncell,1:nLevel,1:1),      &
                qv3d    = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv), &
                p3d     = pstate_wrf%p_phy(1:ncell,1:nLevel,1:1),      &
                p3di    = pstate_wrf%p8w  (1:ncell,1:nLevel,1:1),      &
                pi3d    = pstate_wrf%pi_phy( 1:ncell,1:nLevel,1:1),    &
                z       = pstate_wrf%zzz  (1:ncell,1:nLevel,1:1),      &
                rublten = ptend_wrf%ruublten(1:ncell,1:nLevel,1:1),    &
                rvblten = ptend_wrf%rvvblten(1:ncell,1:nLevel,1:1),    &
                dtaux3d = dtaux3d(1:ncell,1:nLevel,1:1), & ! inout, 1st calc inside
                dtauy3d = dtauy3d(1:ncell,1:nLevel,1:1), & ! inout, 1st calc inside
                dusfcg  = dusfcg (1:ncell,         1:1), & ! inout, 1st calc inside
                dvsfcg  = dvsfcg (1:ncell,         1:1), & ! inout, 1st calc inside
                var2d   = psurf_wrf%var2d(1:ncell,1:1) , & ! filled from static
                oc12d   = psurf_wrf%oc12d(1:ncell,1:1) , & ! filled from static
                oa2d1   = psurf_wrf%oa2d1(1:ncell,1:1) , & ! filled from static
                oa2d2   = psurf_wrf%oa2d2(1:ncell,1:1) , & ! filled from static
                oa2d3   = psurf_wrf%oa2d3(1:ncell,1:1) , & ! filled from static
                oa2d4   = psurf_wrf%oa2d4(1:ncell,1:1) , & ! filled from static
                ol2d1   = psurf_wrf%ol2d1(1:ncell,1:1) , & ! filled from static
                ol2d2   = psurf_wrf%ol2d2(1:ncell,1:1) , & ! filled from static
                ol2d3   = psurf_wrf%ol2d3(1:ncell,1:1) , & ! filled from static
                ol2d4   = psurf_wrf%ol2d4(1:ncell,1:1) , & ! filled from static
              ! optional, not used now
              ! znu,
              ! znw,
              ! mut,
              ! p_top,                         &
              ! optional, not used now
                 cp     = cp ,    & ! constants
                 g      = gravity,& ! constants
                 rd     = r_d,    & ! constants
                 rv     = r_v,    & ! constants
                 ep1    = ep_1,   & ! constants
                 pi     = pi,     & ! constants
                 dt     = dtime,  &
                 dx     = pstate_wrf%dxmean(1:ncell), &
                 kpbl2d = pstate_wrf%kpbl(1:ncell,1:1), & ! from pbl scheme
                 itimestep = itimestep,      &
                 ids=1,ide = ncell, jds=1,jde=1, kds=1,kde=nLevel, &
                 ims=1,ime = ncell, jms=1,jme=1, kms=1,kme=nLevel, &
                 its=1,ite = ncell, jts=1,jte=1, kts=1,kte=nLevel-1)

     end if

END IF

     return
   end subroutine  grist_PhysW_pbl_run 

end module grist_PhysW_pbl_driver
