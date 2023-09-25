
!----------------------------------------------------------------------------
! Created on 2019
! Author: Yi Zhang
! Version 1.0
! Description: This driver mainly follows the naming convention of WRF, but 
!              the dim conversion requires a transpose from GRIST based vars, 
!              currently only support Lin1983
! Revision history:
!              1. LIN_V381, MORR_TWOM_V381, WSM6_V381, YLINV381 can be used
!--------------------------------------------------------------------------

module grist_PhysW_microphysics_driver

  use grist_constants,          only: i4, r8, r_d=>rdry, r_v=>rvap, cp, cv, gravity, zero
  use grist_mpi
  use grist_nml_module,                only: sub_physpkg
  use grist_physics_idealized_package, only: grist_idealized_physics_kessler_klemp15
  use grist_PhysW_data_structure, only: pstate_wrf, ptend_wrf, psurf_wrf
  use grist_PhysW_data_structure, only: p_qv, p_qc, p_qr, p_qi, p_qs, p_qg, p_ng, p_ni, p_ns, p_nr, param_first_scalar
! PhysW nml
  use grist_PhysW_nml_module, only: PhysW_mp_scheme, use_cond
! wrf module
  use module_wrfmodel_constants,only: xls, xlv, xlf, rhowater, rhosnow, ep_2, svp1, svp2, svp3, svpt0, cpd, cpv, &
                                      psat, cliq, cice,rhoair0,rhowater,svpt0,epsilon,ep_1,ep_2
!
! single moment scheme
!
  use module_mp_lin_v381,       only: lin_v381=>lin_et_al
  use module_mp_sbu_ylin_v381,  only: sbu_ylin_v381=>sbu_ylin
#ifdef ZHM_WSM6
  use module_mp_wsm6_v381_EXP,  only: wsm6init_v381=>wsm6init, wsm6_v381=>wsm6 ! tested for climate runs (hdc&ndc)
#else
  use module_mp_wsm6_v381,      only: wsm6init_v381=>wsm6init, wsm6_v381=>wsm6 ! tested for climate runs (hdc&ndc)
#endif
!
! two-moment scheme added here
!
  use module_mp_morr_two_moment_v381, only: morr_twom_init_v381=>MORR_TWO_MOMENT_INIT, morr_twom_run_v381=>MP_MORR_TWO_MOMENT

#ifdef NJCXY
  ! Xiaoqi Xu added this, modified some features
  ! A: aerosol, C: cloud, E: entrainment
  use module_mp_morr_two_moment_v381_ACE, only: morr_twom_init_v381_ace=>MORR_TWO_MOMENT_INIT, &
                                                morr_twom_run_v381_ace=>MP_MORR_TWO_MOMENT
#endif

! macro pcond, just for test!:)
  use cam3_cloud_fraction,      only: pcond

  implicit  none

  private

  public         ::  grist_PhysW_microphysics_init, &
                     grist_PhysW_microphysics_run,  &
                     grist_PhysW_microphysics_final

    real(r8), allocatable   :: moist_old(:,:,:,:)
    real(r8), allocatable   :: th_old(:,:,:)

contains

   subroutine grist_PhysW_microphysics_init(ncell,nLevel,nspecies)
    integer(i4), intent(in) :: ncell, nLevel,  nspecies

! These will always be needed by any microphysics, six-class at most for now
    if(.not.allocated(moist_old))            allocate(moist_old(1:ncell,1:nLevel, 1:1, nspecies))
    if(.not.allocated(th_old))               allocate(   th_old(1:ncell,1:nLevel, 1:1))
    if(.not.allocated(ptend_wrf%rthmpten))   allocate(ptend_wrf%rthmpten(1:ncell,1:nLevel, 1:1))
    if(.not.allocated(ptend_wrf%rqvmpten))   allocate(ptend_wrf%rqvmpten(1:ncell,1:nLevel, 1:1))
    if(.not.allocated(ptend_wrf%rqcmpten))   allocate(ptend_wrf%rqcmpten(1:ncell,1:nLevel, 1:1))
    if(.not.allocated(ptend_wrf%rqrmpten))   allocate(ptend_wrf%rqrmpten(1:ncell,1:nLevel, 1:1))
    if(.not.allocated(ptend_wrf%rqimpten))   allocate(ptend_wrf%rqimpten(1:ncell,1:nLevel, 1:1))
    if(.not.allocated(ptend_wrf%rqsmpten))   allocate(ptend_wrf%rqsmpten(1:ncell,1:nLevel, 1:1))
    if(.not.allocated(ptend_wrf%rqgmpten))   allocate(ptend_wrf%rqgmpten(1:ncell,1:nLevel, 1:1))

    select case(trim(PhysW_mp_scheme))
    case('WSM6V381')
       call wsm6init_v381(rhoair0,rhowater,rhosnow,cliq,cpv,0,.false.)
    case('LINV381','YLINV381')
     ! do nothing
    case('MORR_TWOM_V381')
! needed by two-moment MP like morrison
       if(.not.allocated(ptend_wrf%rnrmpten))   allocate(ptend_wrf%rnrmpten(1:ncell,1:nLevel, 1:1))
       if(.not.allocated(ptend_wrf%rnimpten))   allocate(ptend_wrf%rnimpten(1:ncell,1:nLevel, 1:1))
       if(.not.allocated(ptend_wrf%rnsmpten))   allocate(ptend_wrf%rnsmpten(1:ncell,1:nLevel, 1:1))
       if(.not.allocated(ptend_wrf%rngmpten))   allocate(ptend_wrf%rngmpten(1:ncell,1:nLevel, 1:1))

       call morr_twom_init_v381(0)
#ifdef NJCXY
    case('MORR_TWOM_V381_ACE')
! needed by two-moment MP like morrison
       if(.not.allocated(ptend_wrf%rnrmpten))   allocate(ptend_wrf%rnrmpten(1:ncell,1:nLevel, 1:1))
       if(.not.allocated(ptend_wrf%rnimpten))   allocate(ptend_wrf%rnimpten(1:ncell,1:nLevel, 1:1))
       if(.not.allocated(ptend_wrf%rnsmpten))   allocate(ptend_wrf%rnsmpten(1:ncell,1:nLevel, 1:1))
       if(.not.allocated(ptend_wrf%rngmpten))   allocate(ptend_wrf%rngmpten(1:ncell,1:nLevel, 1:1))

       call morr_twom_init_v381_ace(0)
#endif
   
    case default
       if(mpi_rank().eq.0) print*, "you must select a MP scheme, LIN, WSM6, WSM6V381 stop"
       call mpi_abort()
    end select

    return
   end subroutine grist_PhysW_microphysics_init

   subroutine grist_PhysW_microphysics_final()

       if(allocated(moist_old))            deallocate(moist_old)
       if(allocated(th_old))               deallocate(   th_old)
       if(allocated(ptend_wrf%rthmpten))   deallocate(ptend_wrf%rthmpten)
       if(allocated(ptend_wrf%rqvmpten))   deallocate(ptend_wrf%rqvmpten)
       if(allocated(ptend_wrf%rqcmpten))   deallocate(ptend_wrf%rqcmpten)
       if(allocated(ptend_wrf%rqrmpten))   deallocate(ptend_wrf%rqrmpten)
       if(allocated(ptend_wrf%rqimpten))   deallocate(ptend_wrf%rqimpten)
       if(allocated(ptend_wrf%rqsmpten))   deallocate(ptend_wrf%rqsmpten)
       if(allocated(ptend_wrf%rqgmpten))   deallocate(ptend_wrf%rqgmpten)
! two-moment
       if(allocated(ptend_wrf%rnrmpten))   deallocate(ptend_wrf%rnrmpten)
       if(allocated(ptend_wrf%rnimpten))   deallocate(ptend_wrf%rnimpten)
       if(allocated(ptend_wrf%rnsmpten))   deallocate(ptend_wrf%rnsmpten)
       if(allocated(ptend_wrf%rngmpten))   deallocate(ptend_wrf%rngmpten)

     return
   end subroutine grist_PhysW_microphysics_final

   subroutine grist_PhysW_microphysics_run(ncell,nLevel,nspecies,itimestep,dtime)

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
!----------------------------------------------------------------------
!-- th_phy        potential temperature    (k)
!-- moist_new     updated moisture array   (kg/kg)
!-- moist_old     old moisture array       (kg/kg)
!-- rho           density of air           (kg/m^3)
!-- pi            exner function           (dimensionless)
!-- p             pressure                 (pa)
!-- rainnc        grid scale precipitation (mm)
!-- rainncv       one time step grid scale precipitation (mm/step)
!!!-- sr            one time step mass ratio of snow to total precip
!-- z             height above sea level   (m)
!-- dt            time step              (s)
!-- config_flags  flag for configuration      ! change ---  ?????   
!-- nspecies       number of water substances   (integer)
!-- g             acceleration due to gravity  (m/s^2)
!-- cp            heat capacity at constant pressure for dry air (j/kg/k)
!-- r_d           gas constant for dry air (j/kg/k)
!-- r_v           gas constant for water vapor (j/kg/k)
!-- xls           latent heat of sublimation   (j/kg)
!-- xlv           latent heat of vaporization  (j/kg)
!-- xlf           latent heat of melting       (j/kg)
!-- rhowater      water density                      (kg/m^3)
!-- rhosnow       snow density               (kg/m^3)
!-- f_ice_phy     fraction of ice.
!-- f_rain_phy    fraction of rain.
!-- f_rimef_phy   mass ratio of rimed ice (rime factor)
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
!-- i_start       start indices for i in tile
!-- i_end         end indices for i in tile
!-- j_start       start indices for j in tile
!-- j_end         end indices for j in tile
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!-- num_tiles     number of tiles
!======================================================================
! io
   integer(i4), intent(in)    :: ncell
   integer(i4), intent(in)    :: nLevel
   integer(i4), intent(in)    :: nspecies
   integer(i4), intent(in)    :: itimestep
   real(r8),    intent(in)    :: dtime
!
! local
!
   integer(i4)   :: icell
   real(r8)      :: surface_rain(ncell)
   real(r8)      :: precl(ncell) ! m/s as rate
   real(r8)      :: re_cloud(ncell,nLevel,1)
   real(r8)      :: re_ice  (ncell,nLevel,1)
   real(r8)      :: re_snow (ncell,nLevel,1)
   real(r8)      :: ri3d    (ncell,nLevel,1)
   real(r8)      :: pmid(ncell,nLevel-1)     ! Pressure at layer midpoints
   real(r8)      :: pdel(ncell,nLevel-1)     ! Delta p at each model level
   real(r8)      :: tttt(ncell,nLevel-1)     ! Temperature
   real(r8)      :: qqqq(ncell,nLevel-1)     ! Specific humidity
   real(r8)      :: qclv(ncell,nLevel-1)     ! Moisture tendency due to rainout
   real(r8)      :: thpt(ncell,nLevel-1)     ! Moisture tendency due to rainout

   ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1) = zero
   ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1) = zero
   ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1) = zero

!   if(use_cond)then
!
! large-scale condensation rate, hold its tends and precl,
! modify physics state for next processes, from top to bottom as in CAM
!

!      pmid(1:ncell,1:nLevel-1) = pstate_wrf%p_phy(1:ncell,nLevel-1:1:-1,1)
!     pdel(1:ncell,1:nLevel-1) = pstate_wrf%p8w  (1:ncell,nLevel-1:1:-1,1)-pstate_wrf%p8w  (1:ncell,nLevel:2:-1,1)
!     tttt(1:ncell,1:nLevel-1) = pstate_wrf%t_phy(1:ncell,nLevel-1:1:-1,1)
!     qqqq(1:ncell,1:nLevel-1) = pstate_wrf%moist(1:ncell,nLevel-1:1:-1,1,p_qv)
!     qclv(1:ncell,1:nLevel-1) = pstate_wrf%moist(1:ncell,nLevel-1:1:-1,1,p_qc)

!     call pcond(plond = ncell, plon=ncell, plev=nLevel-1, &
!                tdt   = dtime, & ! in
!                pmid  = pmid , & ! in
!                pdel  = pdel , & ! in
!                t     = tttt , & ! inout
!                q     = qqqq , & ! inout
!                qc    = qclv , & ! inout
!                precl = precl)   ! out

! obtain th based on new t assuming pii unchanged as other physics
!      thpt(1:ncell,1:nLevel-1) = tttt(1:ncell,1:nLevel-1)/pstate_wrf%pi_phy(1:ncell,nLevel-1:1:-1,1)

!    
! store this part of explicit tend into mptend
! bottom to top
!

!     ptend_wrf%rthmpten(1:ncell,1:nLevel,1) = (thpt(1:ncell,nLevel-1:1:-1)-pstate_wrf%th_phy(1:ncell,1:nLevel-1,1))/dtime
!      ptend_wrf%rqvmpten(1:ncell,1:nLevel,1) = (qqqq(1:ncell,nLevel-1:1:-1)-pstate_wrf%moist (1:ncell,1:nLevel-1,1,p_qv))/dtime
!      ptend_wrf%rqcmpten(1:ncell,1:nLevel,1) = (qclv(1:ncell,nLevel-1:1:-1)-pstate_wrf%moist (1:ncell,1:nLevel-1,1,p_qc))/dtime

!
! let pcond modify t, th, qv, qc, then use them as input for later physics
! 

!     pstate_wrf%t_phy (1:ncell,1:nLevel-1,1)      = tttt(1:ncell,nLevel-1:1:-1)
!     pstate_wrf%th_phy(1:ncell,1:nLevel-1,1)      = thpt(1:ncell,nLevel-1:1:-1)
!     pstate_wrf%moist (1:ncell,1:nLevel-1,1,p_qv) = qqqq(1:ncell,nLevel-1:1:-1)
!     pstate_wrf%moist (1:ncell,1:nLevel-1,1,p_qc) = qclv(1:ncell,nLevel-1:1:-1)

!    end if

!---------------------------------------------------------------------
!  check for microphysics type.  we need a clean way to 
!  specify these things!
!---------------------------------------------------------------------

! run this produce same solutions as run outside in the simple physics package
! can not run now because hardcode changes
      if(trim(sub_physpkg).eq.'DCMIP2016-SC')then
        do icell = 1, ncell
         call grist_idealized_physics_kessler_klemp15(pstate_wrf%th_phy(icell,1:nLevel-1,1)     ,&
                                                      pstate_wrf%moist (icell,1:nLevel-1,1,p_qv),&
                                                      pstate_wrf%moist (icell,1:nLevel-1,1,p_qc),&
                                                      pstate_wrf%moist (icell,1:nLevel-1,1,p_qr),&
                                                      pstate_wrf%rhod  (icell,1:nLevel-1,1)     ,&
                                                      pstate_wrf%pi_phy(icell,1:nLevel-1,1)     ,&
                                                      dtime                                     ,&
                                                      pstate_wrf%zzz(icell,1:nLevel-1,1)        ,&
                                                      nLevel-1                                  ,&
                                                      psurf_wrf%rainnc (icell,1)               ,&
                                                      ptend_wrf%rthmpten(icell,1:nLevel-1,1)    ,&
                                                      ptend_wrf%rqvmpten(icell,1:nLevel-1,1)    ,&
                                                      ptend_wrf%rqcmpten(icell,1:nLevel-1,1)    ,&
                                                      ptend_wrf%rqrmpten(icell,1:nLevel-1,1))
        end do

      else

      select case(trim(PhysW_mp_scheme))

      case('LINV381')

      moist_old(1:ncell,1:nLevel,1:1,1:nspecies) = pstate_wrf%moist(1:ncell,1:nLevel,1:1,1:nspecies)
      th_old(1:ncell,1:nLevel,1:1)              = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)

        call lin_v381(th      = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1),      &
                      qv      = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv), &
                      ql      = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc), &
                      qr      = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qr), &
                      qi      = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi), &
                      qs      = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qs), &
                      rho     = pstate_wrf%rhom ( 1:ncell,1:nLevel,1:1),      & 
                      pii     = pstate_wrf%pi_phy(1:ncell,1:nLevel,1:1),      & 
                      p       = pstate_wrf%p_phy( 1:ncell,1:nLevel,1:1),      &
                      dt_in   = dtime ,                                       &
                      z       = pstate_wrf%zzz(1:ncell,1:nLevel,1:1),         &
                      ht      = psurf_wrf%ht(1:ncell,1:1),                    &
                      dz8w    = pstate_wrf%dz8w( 1:ncell,1:nLevel,1:1),       &
                      grav    = gravity, cp = cpd, Rair = r_d, rvapor= r_v,   &
                      XLS     = xls, XLV=xlv, XLF=xlf, rhowater = rhowater, rhosnow = rhosnow, &
                      EP2     = ep_2, SVP1 = svp1, SVP2 = svp2, SVP3 = svp3, SVPT0  = svpt0  , &
                      RAINNC    = psurf_wrf%rainnc( 1:ncell,1:1 ),  & 
                      RAINNCV   = psurf_wrf%rainncv(1:ncell,1:1),   &
                      SNOWNC    = psurf_wrf%snownc( 1:ncell,1:1),   &
                      SNOWNCV   = psurf_wrf%snowncv(1:ncell,1:1),   &
                      GRAUPELNC = psurf_wrf%graupelnc( 1:ncell,1:1),& 
                      GRAUPELNCV= psurf_wrf%graupelncv(1:ncell,1:1),&
                      SR = surface_rain(1:ncell),                   &
                      !refl_10cm, diagflag, do_radar_ref           &
                      ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,&
                      ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,&
                      its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )
                  ! Optional 
                  !    ,qlsink, precr, preci, precs, precg          &
                  !    , F_QG,F_QNDROP                              &
                  !    , qg, qndrop                                 &
      ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)-th_old(1:ncell,1:nLevel,1:1))/dtime
      ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv)-moist_old(1:ncell,1:nLevel,1:1,p_qv))/dtime
      ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc)-moist_old(1:ncell,1:nLevel,1:1,p_qc))/dtime
      ptend_wrf%rqrmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr)-moist_old(1:ncell,1:nLevel,1:1,p_qr))/dtime
      ptend_wrf%rqimpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi)-moist_old(1:ncell,1:nLevel,1:1,p_qi))/dtime
      ptend_wrf%rqsmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs)-moist_old(1:ncell,1:nLevel,1:1,p_qs))/dtime

      case('WSM6V381')
!
! store here for manully diagnosing tendency
!
      moist_old(1:ncell,1:nLevel,1:1,1:nspecies) = pstate_wrf%moist(1:ncell,1:nLevel,1:1,1:nspecies)
      th_old(1:ncell,1:nLevel,1:1)               = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)
      call wsm6_v381(&
             th  = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1),       &
             q   = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv),  &
             qc  = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc),  &
             qr  = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qr),  &
             qi  = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi),  &
             qs  = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qs),  &
             qg  = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qg),  &
             den = pstate_wrf%rhom(  1:ncell,1:nLevel,1:1),       &
             pii = pstate_wrf%pi_phy(1:ncell,1:nLevel,1:1),       &
             p   = pstate_wrf%p_phy( 1:ncell,1:nLevel,1:1),       &
             delz= pstate_wrf%dz8w(  1:ncell,1:nLevel,1:1),       &
             delt= dtime,       &
             g   = gravity,     &
             cpd = cpd,         &
             cpv = cpv,         &
             rd  = r_d,         &
             rv  = r_v,         &
             t0c = svpt0,       &
             ep1 = ep_1,        &
             ep2 = ep_2,        &
             qmin= epsilon,     &
             XLS = xls,         &
             XLV0= xlv,         &
             XLF0= xlf,         &
             den0= rhoair0,     &
             denr= rhowater,    &
             cliq= cliq,        &
             cice= cice,        &
             psat= psat,        &
             rain    = psurf_wrf%rainnc( 1:ncell,1:1 ), & !mm
             rainncv = psurf_wrf%rainncv(1:ncell,1:1 ), & !mm/step
             snow    = psurf_wrf%snownc( 1:ncell,1:1),  & !mm
             snowncv = psurf_wrf%snowncv(1:ncell,1:1),  & !mm/step
             sr      = surface_rain(1:ncell),           &
             refl_10cm  = pstate_wrf%refl_10cm(1:ncell,1:nLevel,1:1), & 
             diagflag   = .true., & 
             do_radar_ref=  1, &
             graupel    = psurf_wrf%graupelnc( 1:ncell ,1:1), &  ! 
             graupelncv = psurf_wrf%graupelncv( 1:ncell ,1:1),&  !
             has_reqc   = 1, has_reqi = 1 , has_reqs = 1 ,    &  ! for radiation
             re_cloud   = pstate_wrf%re_cloud(1:ncell,1:nLevel,1:1),& ! for radiation
             re_ice     = pstate_wrf%re_ice  (1:ncell,1:nLevel,1:1),& ! for radiation 
             re_snow    = pstate_wrf%re_snow (1:ncell,1:nLevel,1:1),& ! for radiation
             ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,  &
             ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
             its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )
!
! explicit tendencies
!
      ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)-th_old(1:ncell,1:nLevel,1:1))/dtime
      ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv)-moist_old(1:ncell,1:nLevel,1:1,p_qv))/dtime
      ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc)-moist_old(1:ncell,1:nLevel,1:1,p_qc))/dtime
      ptend_wrf%rqrmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr)-moist_old(1:ncell,1:nLevel,1:1,p_qr))/dtime
      ptend_wrf%rqimpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi)-moist_old(1:ncell,1:nLevel,1:1,p_qi))/dtime
      ptend_wrf%rqsmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs)-moist_old(1:ncell,1:nLevel,1:1,p_qs))/dtime
      ptend_wrf%rqgmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qg)-moist_old(1:ncell,1:nLevel,1:1,p_qg))/dtime

      case('YLINV381')
! store here for manully diagnosing tendency
      moist_old(1:ncell,1:nLevel,1:1,1:nspecies) = pstate_wrf%moist(1:ncell,1:nLevel,1:1,1:nspecies)
      th_old(1:ncell,1:nLevel,1:1)               = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)

      call sbu_ylin_v381(&
                    th     = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1),      &
                    qv     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv), &
                    ql     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc), &
                    qr     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qr), &
                    qi     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi), &
                    qs     = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qs), &
                    Ri3D   = ri3d(1:ncell,1:nLevel,1:1),                   &  ! seems only a diagnosed var
                    rho    = pstate_wrf%rhom(  1:ncell,1:nLevel,1:1),      &
                    pii    = pstate_wrf%pi_phy(1:ncell,1:nLevel,1:1),      &
                    p      = pstate_wrf%p_phy( 1:ncell,1:nLevel,1:1),      &
                    dt_in  = dtime,                                        &
                    z      = pstate_wrf%zzz(   1:ncell,1:nLevel,1:1),      &  ! height above sea level   (m)
                    ht     = psurf_wrf%ht(     1:ncell,1:1),               &  ! surface height
                    dz8w   = pstate_wrf%dz8w(  1:ncell,1:nLevel,1:1),      &
                    RAINNC = psurf_wrf%rainnc( 1:ncell,1:1 ),              &
                    RAINNCV= psurf_wrf%rainncv(1:ncell,1:1 ),              &
                    ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel, &
                    ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel, &
                    its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1)
!
! explicit tendencies
!
      ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)-th_old(1:ncell,1:nLevel,1:1))/dtime
      ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv)-moist_old(1:ncell,1:nLevel,1:1,p_qv))/dtime
      ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc)-moist_old(1:ncell,1:nLevel,1:1,p_qc))/dtime
      ptend_wrf%rqrmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr)-moist_old(1:ncell,1:nLevel,1:1,p_qr))/dtime
      ptend_wrf%rqimpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi)-moist_old(1:ncell,1:nLevel,1:1,p_qi))/dtime
      ptend_wrf%rqsmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs)-moist_old(1:ncell,1:nLevel,1:1,p_qs))/dtime

      case('MORR_TWOM_V381')
          moist_old(1:ncell,1:nLevel,1:1,1:nspecies) = pstate_wrf%moist(1:ncell,1:nLevel,1:1,1:nspecies)
          th_old(1:ncell,1:nLevel,1:1)              = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)

          call morr_twom_run_v381(itimestep = itimestep                       ,      &
                                  th        = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1),      &
                                  qv        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv), &
                                  qc        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc), &
                                  qr        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qr), &
                                  qi        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi), &
                                  qs        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qs), &
                                  qg        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qg), &
                                  ni        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_ni), &
                                  ns        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_ns), &
                                  nr        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_nr), &
                                  ng        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_ng), &
                                  rho       = pstate_wrf%rhom(  1:ncell,1:nLevel,1:1),      & ! not used actually
                                  pii       = pstate_wrf%pi_phy(1:ncell,1:nLevel,1:1), & 
                                  p         = pstate_wrf%p_phy(1:ncell,1:nLevel,1:1) , &
                                  dt_in     = dtime                                  , &
                                  dz        = pstate_wrf%dz8w(  1:ncell,1:nLevel,1:1), & 
                                  ht        = psurf_wrf%ht(1:ncell,1:1)              , &  ! surface height, not used actually 
                                  w         = pstate_wrf%www(1:ncell,1:nLevel,1),& ! not used
                                  rainnc    = psurf_wrf%rainnc( 1:ncell,1:1 ),   &
                                  rainncv   = psurf_wrf%rainncv( 1:ncell,1:1 ),  &
                                  sr        = surface_rain(1:ncell),             &
                                  snownc    = psurf_wrf%snownc( 1:ncell,1:1),    &
                                  snowncv   = psurf_wrf%snowncv(1:ncell,1:1),    &
                                  graupelnc = psurf_wrf%graupelnc( 1:ncell,1:1), &
                                  graupelncv= psurf_wrf%graupelncv(1:ncell,1:1), & ! hm added 7/13/13
                                  ! yizhang: currently we donot diagnose reflectivity
                                  !refl_10cm, 
                                  !diagflag,
                                  !do_radar_ref,      & ! GT added for reflectivity calcs
                                  ! yizhang: currently we donot diagnose reflectivity
                                  qrcuten   = ptend_wrf%rqrcuten(1:ncell, 1:nLevel, 1:1), & ! zero for some schemes
                                  qscuten   = ptend_wrf%rqscuten(1:ncell, 1:nLevel, 1:1), & ! zero for some schemes
                                  qicuten   = ptend_wrf%rqicuten(1:ncell, 1:nLevel, 1:1), & ! need restart 
                        !
                        ! In WRF, mu comes from solve_em, to scale tend from partial(pi)/paritla(eta)
                        ! We do not need this procedure, as all our physics are inside physics
                        !
                                  !mu        = ,           & ! hm added
                                  !yizhang: below are for WRF-chem, not used yet
                                  !f_qndrop, 
                                  !qndrop                        & ! hm added, wrf-chem 
                                  !yizhang: below are for WRF-chem, not used yet
                                  ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel  , & ! domain dims
                                  ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel  , & ! memory dims
                                  its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )! tile   dims            )
!jdf       ,C2PREC3D,CSED3D,ISED3D,SSED3D,GSED3D,RSED3D & ! HM ADD, WRF-CHEM
                                  !yizhang: below are for WRF-chem, not used yet
                                  !rainprod, evapprod                      &
                                  !qlsink,precr,preci,precs,precg &        ! HM ADD, WRF-CHEM

!
! explicit tendencies
!
      ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)-th_old(1:ncell,1:nLevel,1:1))/dtime
      ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv)-moist_old(1:ncell,1:nLevel,1:1,p_qv))/dtime
      ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc)-moist_old(1:ncell,1:nLevel,1:1,p_qc))/dtime
      ptend_wrf%rqrmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr)-moist_old(1:ncell,1:nLevel,1:1,p_qr))/dtime
      ptend_wrf%rqimpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi)-moist_old(1:ncell,1:nLevel,1:1,p_qi))/dtime
      ptend_wrf%rqsmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs)-moist_old(1:ncell,1:nLevel,1:1,p_qs))/dtime
      ptend_wrf%rqgmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qg)-moist_old(1:ncell,1:nLevel,1:1,p_qg))/dtime

      ptend_wrf%rnrmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_nr)-moist_old(1:ncell,1:nLevel,1:1,p_nr))/dtime
      ptend_wrf%rnimpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_ni)-moist_old(1:ncell,1:nLevel,1:1,p_ni))/dtime
      ptend_wrf%rnsmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_ns)-moist_old(1:ncell,1:nLevel,1:1,p_ns))/dtime
      ptend_wrf%rngmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_ng)-moist_old(1:ncell,1:nLevel,1:1,p_ng))/dtime

#ifdef NJCXY
      case('MORR_TWOM_V381_ACE')
          moist_old(1:ncell,1:nLevel,1:1,1:nspecies) = pstate_wrf%moist(1:ncell,1:nLevel,1:1,1:nspecies)
          th_old(1:ncell,1:nLevel,1:1)              = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)

          call morr_twom_run_v381_ace(itimestep = itimestep                       ,      &
                                  th        = pstate_wrf%th_phy(1:ncell,1:nLevel,1:1),      &
                                  qv        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qv), &
                                  qc        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qc), &
                                  qr        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qr), &
                                  qi        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qi), &
                                  qs        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qs), &
                                  qg        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_qg), &
                                  ni        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_ni), &
                                  ns        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_ns), &
                                  nr        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_nr), &
                                  ng        = pstate_wrf%moist( 1:ncell,1:nLevel,1:1,p_ng), &
                                  rho       = pstate_wrf%rhom(  1:ncell,1:nLevel,1:1),      & ! not used actually
                                  pii       = pstate_wrf%pi_phy(1:ncell,1:nLevel,1:1), & 
                                  p         = pstate_wrf%p_phy(1:ncell,1:nLevel,1:1) , &
                                  dt_in     = dtime                                  , &
                                  dz        = pstate_wrf%dz8w(  1:ncell,1:nLevel,1:1), & 
                                  ht        = psurf_wrf%ht(1:ncell,1:1)              , &  ! surface height, not used actually 
                                  w         = pstate_wrf%www(1:ncell,1:nLevel,1),& ! not used
                                  rainnc    = psurf_wrf%rainnc( 1:ncell,1:1 ),   &
                                  rainncv   = psurf_wrf%rainncv( 1:ncell,1:1 ),  &
                                  sr        = surface_rain(1:ncell),             &
                                  snownc    = psurf_wrf%snownc( 1:ncell,1:1),    &
                                  snowncv   = psurf_wrf%snowncv(1:ncell,1:1),    &
                                  graupelnc = psurf_wrf%graupelnc( 1:ncell,1:1), &
                                  graupelncv= psurf_wrf%graupelncv(1:ncell,1:1), & ! hm added 7/13/13
                                  ! yizhang: currently we donot diagnose reflectivity
                                  !refl_10cm, 
                                  !diagflag,
                                  !do_radar_ref,      & ! GT added for reflectivity calcs
                                  ! yizhang: currently we donot diagnose reflectivity
                                  qrcuten   = ptend_wrf%rqrcuten(1:ncell, 1:nLevel, 1:1), & ! zero for some schemes
                                  qscuten   = ptend_wrf%rqscuten(1:ncell, 1:nLevel, 1:1), & ! zero for some schemes
                                  qicuten   = ptend_wrf%rqicuten(1:ncell, 1:nLevel, 1:1), & ! need restart 
                        !
                        ! In WRF, mu comes from solve_em, to scale tend from partial(pi)/paritla(eta)
                        ! We do not need this procedure, as all our physics are inside physics
                        !
                                  !mu        = ,           & ! hm added
                                  !yizhang: below are for WRF-chem, not used yet
                                  !f_qndrop, 
                                  !qndrop                        & ! hm added, wrf-chem 
                                  !yizhang: below are for WRF-chem, not used yet
                                  ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel  , & ! domain dims
                                  ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel  , & ! memory dims
                                  its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1 )! tile   dims            )
!jdf       ,C2PREC3D,CSED3D,ISED3D,SSED3D,GSED3D,RSED3D & ! HM ADD, WRF-CHEM
                                  !yizhang: below are for WRF-chem, not used yet
                                  !rainprod, evapprod                      &
                                  !qlsink,precr,preci,precs,precg &        ! HM ADD, WRF-CHEM

!
! explicit tendencies
!
      ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rthmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%th_phy(1:ncell,1:nLevel,1:1)-th_old(1:ncell,1:nLevel,1:1))/dtime
      ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqvmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv)-moist_old(1:ncell,1:nLevel,1:1,p_qv))/dtime
      ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1) = ptend_wrf%rqcmpten(1:ncell,1:nLevel,1:1)+(pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qc)-moist_old(1:ncell,1:nLevel,1:1,p_qc))/dtime
      ptend_wrf%rqrmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qr)-moist_old(1:ncell,1:nLevel,1:1,p_qr))/dtime
      ptend_wrf%rqimpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qi)-moist_old(1:ncell,1:nLevel,1:1,p_qi))/dtime
      ptend_wrf%rqsmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qs)-moist_old(1:ncell,1:nLevel,1:1,p_qs))/dtime
      ptend_wrf%rqgmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qg)-moist_old(1:ncell,1:nLevel,1:1,p_qg))/dtime

      ptend_wrf%rnrmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_nr)-moist_old(1:ncell,1:nLevel,1:1,p_nr))/dtime
      ptend_wrf%rnimpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_ni)-moist_old(1:ncell,1:nLevel,1:1,p_ni))/dtime
      ptend_wrf%rnsmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_ns)-moist_old(1:ncell,1:nLevel,1:1,p_ns))/dtime
      ptend_wrf%rngmpten(1:ncell,1:nLevel,1:1) = (pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_ng)-moist_old(1:ncell,1:nLevel,1:1,p_ng))/dtime
#endif

      case default
        if(mpi_rank().eq.0) print*, "you must select a MP scheme: LIN, WSM6, stop"
        call mpi_abort()
      end select
      
      end if

! if os, also update t_phy inside
    pstate_wrf%t_phy (1:ncell,1:nLevel-1,1) = pstate_wrf%th_phy(1:ncell,1:nLevel-1,1)*pstate_wrf%pi_phy(1:ncell,1:nLevel-1,1)

    return
   end subroutine grist_PhysW_microphysics_run

end module grist_PhysW_microphysics_driver
