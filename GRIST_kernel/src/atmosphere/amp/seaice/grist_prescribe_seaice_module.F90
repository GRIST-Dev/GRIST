
!---------------------------------------------------------
! Purpose: A prescribeed sea ice module
! version: 1.0
! History: First ported to GRIST by Rongxy at ~2020-July/Aug
!    Modified to GRIST-A21 version at 2021-01-25 yizhang
!    Remove pstate dependancy with iostream
!    Add albice from CAM3 code
!---------------------------------------------------------


MODULE grist_prescribe_seaice_module
  
  use grist_data_types,   only: scalar_1d_field, scalar_2d_field
  use grist_constants,    only: gravity,i4,r8,zero
  use grist_mpi
  use grist_nml_module,   only: seaice_albedo_scheme

  implicit none
  
  public :: grist_seaice_init, grist_seaice_run

  type(scalar_1d_field), public   :: sicetemp_at_pc_surface   ! sea ice temp
!  type(scalar_1d_field), public   :: sicealb_at_pc_surface    ! sea ice albedo

  real(r8), parameter   :: tmelt   = 273.15_r8    !  melting temperature of ice/snow
  real(r8), parameter   :: ctfreez = 271.38_r8   !   temperature at which sea starts freezing/melting

  !
#ifdef AMIPW_CLIMATE
  REAL(r8), PARAMETER   ::    calbmns = 0.70_r8  ! minimum (snow)
  REAL(r8), PARAMETER   ::    calbmxs = 0.85_r8  ! maximum (snow)
  REAL(r8), PARAMETER   ::    calbmni = 0.60_r8  ! minimum (bare sea ice)
  REAL(r8), PARAMETER   ::    calbmxi = 0.75_r8  ! maximum (bare sea ice)
#else
  REAL(r8), PARAMETER   ::    calbmns = 0.70_r8  ! minimum (snow)
  REAL(r8), PARAMETER   ::    calbmxs = 0.85_r8  ! maximum (snow)
  REAL(r8), PARAMETER   ::    calbmni = 0.60_r8  ! minimum (bare sea ice)
  REAL(r8), PARAMETER   ::    calbmxi = 0.75_r8  ! maximum (bare sea ice)
#endif

CONTAINS

!
!=====================================================================================
!
SUBROUTINE grist_seaice_init(ncols)
!
!  sea ice initialization
!
    implicit none

    integer, intent(in)   :: ncols

    if(.not.allocated(sicetemp_at_pc_surface%f))     allocate(sicetemp_at_pc_surface%f(ncols))
!    if(.not.allocated(sicealb_at_pc_surface%f))      allocate(sicealb_at_pc_surface%f(ncols))
 
    sicetemp_at_pc_surface%f(:) = tmelt   ! set sea ice temperature to melting point
!    sicealb_at_pc_surface%f(:)  = 0.75_r8


END SUBROUTINE grist_seaice_init

!
!===========================================================================
!

SUBROUTINE grist_seaice_run(ncols, nlev, lats, coszrs, ocnfrac_at_pc_surface, &
                            pseaice, snowh,  plwds, pswds,  &
                            pahfs,  pahfl, pevap, plwup, ptaux, ptauy, ptsi, pts, &
                            palbdvis, palbivis, palbdnir, palbinir              , &
                            zbot, ubot, vbot, thbot, tbot, qbot, rbot )
!
!   sea ice temperature, albedo and surface atm-sea fluxes
!
!   use grist_physics_data_structure, only: pstate   
   use grist_domain_types, only : global_domain
   USE grist_nml_module,   only : model_timestep

   implicit none
!
! we need sensible, latent/qflx, ts, albedo from this module for dtp2
! argument
   integer(i4),          intent(in)            :: ncols, nlev
   real(r8),             intent(in)            :: lats(ncols), coszrs(ncols), ocnfrac_at_pc_surface(ncols)
   real(r8),             intent(in)            :: pseaice(ncols), snowh(ncols), plwds(ncols), pswds(ncols) ! ice concentration, snow depth, lw and sw down
   real(r8),             intent(inout)         :: pahfs(ncols),    pahfl(ncols)        ! air sensible and latent (hfx, lh), needed
   real(r8),             intent(inout)         :: pevap(ncols),    plwup(ncols)        ! qflx (needed), lwup
   real(r8),             intent(inout)         :: ptaux(ncols),    ptauy(ncols)        ! taux, tauy
   real(r8),             intent(inout)         :: ptsi(ncols),     pts(ncols)          ! ts-ice, ts (needed)
   real(r8),             intent(inout)         :: palbdvis(ncols), palbivis(ncols)     ! 4-band albedo, needed, merge by seaice
   real(r8),             intent(inout)         :: palbdnir(ncols), palbinir(ncols)
   real(r8), intent(in)                        :: zbot(ncols), ubot(ncols), vbot(ncols), thbot(ncols),tbot(ncols), qbot(ncols), rbot(ncols)
!local
   real(r8) :: psni(ncols)
   real(r8) :: zseni(ncols),  zlati(ncols),  zlwupi(ncols), zalbi(ncols), psiced(ncols),    &
               zevapi(ncols), ztauxi(ncols), ztauyi(ncols), zlwdn(ncols), zswdn(ncols)
   real(r8) :: wgt1, wgt2
   real(r8) :: asdir(ncols), asdif(ncols), aldir(ncols), aldif(ncols)
   integer(i4)        :: jl, mask(ncols)
   character(len=128) :: seaice_albedo

   psni(1:ncols) = snowh(1:ncols)
!
! cam3 gives higher sp albedo than np, realistic, but it does not account night
! points, changed for diagnostic
!
   seaice_albedo = trim(seaice_albedo_scheme)
#ifdef REG_AMP21
   seaice_albedo = 'CAM3'
#endif

! ocean mask
   mask(1:ncols) = nint(ocnfrac_at_pc_surface(1:ncols),i4) ! lt.0.5 is land

   do jl = 1,ncols
      if(mask(jl) .ge. 0.5 .and. pseaice(jl) .ge. 0.01_r8) then ! must be sea and have seaice
        if(lats(jl) .gt. 0._r8)then
           psiced(jl) = 2.0_r8  ! northern hemisphere
        else
           psiced(jl) = 1.0_r8   !southern hemisphere
        endif
      else
        !pseaice(jl) = 0.0_r8
        psiced(jl)  = 0.0_r8
        ptsi(jl)    = tmelt
      endif
   end do

!   psni(:) = 0.0_r8 ! snow not taken into account

! atm - sea ice fluxes
   call grist_seaice_flux(ncols, nlev, psiced, ptsi, &
                          zseni, zlati, zlwupi, zevapi, ztauxi, ztauyi, &
                          zbot, ubot, vbot, thbot, tbot, qbot, rbot)
 
! sea ice albedo  
   select case(trim(seaice_albedo))
   case('CAM3')
     call cam3_albice(ncols, ncols, tbot,  psni, coszrs, pseaice, psiced, &
                 asdir , aldir  , asdif, aldif )
     zalbi = asdir
   case('ECHAM5')
     call get_seaice_albedo (ncols, pseaice, ptsi, psni, zalbi)
   case default
     if(mpi_rank().eq.0) print*,"you must select a seaice albedo scheme, model stop"
     call mpi_abort()
   end select

   !surface net downward radiation
   zlwdn(:)=plwds(:)+zlwupi(:)         !long wave
   zswdn(:)=pswds(:)*(1._r8-zalbi(:))  !short wave

   ! sea ice temperature
   call grist_seaice_temp (ncols, model_timestep,                 &
                           psiced,     psni,      ptsi,           &
                           zlwdn,      zswdn,    zseni,    zlati )

!$omp parallel  private(jl,wgt1,wgt2) 
!$omp do schedule(static,10) 
   do jl = 1,ncols
     if(mask(jl) .ge. 0.5  .and. psiced(jl) .ge. 0.1_r8) then
        wgt1 = pseaice(jl)
        wgt2 = 1._r8-wgt1
        pahfs(jl) = pahfs(jl)*wgt2 - zseni(jl)*wgt1       ! ice flux is downward positive
        pahfl(jl) = pahfl(jl)*wgt2 - zlati(jl)*wgt1
        plwup(jl) = plwup(jl)*wgt2 - zlwupi(jl)*wgt1
        pevap(jl) = pevap(jl)*wgt2 - zevapi(jl)*wgt1
        ptaux(jl) = ptaux(jl)*wgt2 - ztauxi(jl)*wgt1
        ptauy(jl) = ptauy(jl)*wgt2 - ztauyi(jl)*wgt1
        pts(jl)   = pts(jl)*wgt2+ptsi(jl)*wgt1
!        sicealb_at_pc_surface%f(jl)=zalbi(jl)

        if(coszrs(jl).gt.zero)then
          if(trim(seaice_albedo).eq.'ECHAM5')then
            palbdvis(jl) = palbdvis(jl)*wgt2 + zalbi(jl)*wgt1
            palbivis(jl) = palbivis(jl)*wgt2 + zalbi(jl)*wgt1
            palbdnir(jl) = palbdnir(jl)*wgt2 + zalbi(jl)*wgt1
            palbinir(jl) = palbinir(jl)*wgt2 + zalbi(jl)*wgt1
          else if ( trim(seaice_albedo).eq.'CAM3' ) then
            palbdvis(jl) = palbdvis(jl)*wgt2 + asdir(jl)*wgt1
            palbivis(jl) = palbivis(jl)*wgt2 + asdif(jl)*wgt1
            palbdnir(jl) = palbdnir(jl)*wgt2 + aldir(jl)*wgt1
            palbinir(jl) = palbinir(jl)*wgt2 + aldif(jl)*wgt1
          else
! do nothing
          end if
        else
          palbdvis(jl) = zero
          palbivis(jl) = zero
          palbdnir(jl) = zero
          palbinir(jl) = zero
        endif
     endif
   end do
!$omp end do nowait
!$omp end parallel 

   return

END SUBROUTINE grist_seaice_run
!
!=====================================================================================
!
SUBROUTINE grist_seaice_flux(ncols, nlev, psiced, ptsi, &
                                          sen, lat, lwup, evap, taux, tauy,&
                                          zbot, ubot, vbot, thbot, tbot, qbot, rbot)
!
! atm-sea ice fluxes
!
!   use grist_physics_data_structure, only: pstate
   use grist_amp_shr_flux_mod,       only: shr_flux_atmIce   

   implicit none

   ! arguments
   integer,  intent(in)    :: ncols, nlev
   real(r8), intent(in)    :: psiced(ncols), ptsi(ncols)
   real(r8), intent(inout) :: sen(ncols),  lat(ncols),  lwup(ncols), &
                              evap(ncols), taux(ncols), tauy(ncols)
   real(r8), intent(in)    :: zbot(ncols), ubot(ncols), vbot(ncols), &
                              thbot(ncols),tbot(ncols), qbot(ncols), rbot(ncols)
   ! local
   
   real(r8) :: tref(ncols), qref(ncols), ts(ncols)
   integer  :: mask(ncols) 
   integer  :: jl
!
   
   ! set mask
   do jl=1,ncols
      if(psiced(jl) .ge. 0.1_r8) then
         mask(jl) = 1
      else
         mask(jl) = 0
      endif
   end do

   ! model states
   !zbot (1:ncols)  = pstate%z_at_pc_full_level%f(nlev,1:ncols)
   !ubot (1:ncols)  = pstate%u_wind_at_pc_full_level%f(nlev,1:ncols)
   !vbot (1:ncols)  = pstate%v_wind_at_pc_full_level%f(nlev,1:ncols)
   !thbot(1:ncols)  = pstate%temp_at_pc_full_level%f(nlev,1:ncols)*pstate%exner_at_pc_full_level%f(nlev,1:ncols)
   !qbot (1:ncols)  = pstate%tracer_mxrt_at_pc_full_level%f(1,nlev,1:ncols)
   !rbot (1:ncols)  = pstate%delp_at_pc_full_level%f(nlev,1:ncols)/gravity/ &
   !                           (pstate%z_at_pc_face_level%f(nlev,1:ncols)-pstate%z_at_pc_face_level%f(nlev+1,1:ncols))
   !tbot (1:ncols)  = pstate%temp_at_pc_full_level%f(nlev,1:ncols)
   ts   (1:ncols)  = ptsi(1:ncols)

   ! atm-ice fluxes
   call shr_flux_atmIce(mask  ,zbot  ,ubot  ,vbot  ,thbot    &
                         ,qbot  ,rbot  ,tbot  ,ts    ,sen    &
                         ,lat   ,lwup  ,evap  ,taux  ,tauy   &
                         ,tref  ,qref                        )

END SUBROUTINE grist_seaice_flux
!
!=====================================================================================
!
SUBROUTINE get_seaice_albedo (ncols, pseaice, ptsi, psni, palbi)
 
   !- Author:
  !
  !  Marco Giorgetta, MPI, May 2000
  !  Revised by Rongxy, Aug 10, 2020 for coupling with grist
 
   implicit none

   ! arguments
   integer, intent(in)     :: ncols
   real(r8), intent(in)    :: pseaice(ncols), ptsi(ncols), psni(ncols)
   real(r8), intent(inout) :: palbi(ncols)

   ! local
   REAL(r8) :: ztalb    ! upper temp. limit for cold snow albedo
   REAL(r8) :: ztsalb   ! upper temp. limit for cold snow albedo on sea ice
   REAL(r8) :: zalbmax  ! maximum snow albedo
   REAL(r8) :: zalbmin  ! minimum snow albedo
   REAL(r8) :: zdalb    ! snow albedo change per deg C

   INTEGER :: jl    ! loop index


   ztalb=tmelt-5.0_r8
   ztsalb=tmelt-1.0_r8

   DO jl = 1,ncols
       !IF (pseaice(jl)>0.0_r8) THEN
       IF (pseaice(jl)>0.01_r8) THEN ! yizhang modifies here to avoid modifying pseaice in the driver
          ! minimum and maximum albedo
          IF (psni(jl)>0.01_r8) THEN
             ! on snow covered sea ice
             zalbmin=calbmns
             zalbmax=calbmxs
          ELSE
             ! on bare sea ice
             zalbmin=calbmni
             zalbmax=calbmxi
          END IF
          ! temperature dependent snow albedo
          IF (ptsi(jl)>=tmelt) THEN
             palbi(jl)=zalbmin
          ELSE IF (ptsi(jl)<ztsalb) THEN
             palbi(jl)=zalbmax
          ELSE
             zdalb=(zalbmax-zalbmin)/(tmelt-ztsalb)
             palbi(jl)=zalbmin+zdalb*(tmelt-ptsi(jl))
          END IF
       ELSE
#ifdef AMIPW_CLIMATE
          palbi(jl) = 0.85_r8
#else
          palbi(jl) = 0.75_r8
#endif
       END IF
   END DO
   
   RETURN

END SUBROUTINE get_seaice_albedo

!
!=====================================================================================
!
SUBROUTINE grist_seaice_temp (ncols, delta_time                         &
                             , psiced,     psni,      ptsi              &
                             , ptrfli,     psofli,    pahfsi,    pahfli )              
  ! Description:
  !
  ! Prognostic calculation of sea-ice temperature
  !
  !
  ! Authors:
  !
  ! F. Lunkeit, MI, April 1991, original source
  ! L. Kornblueh, MPI, May 1998, f90 rewrite
  ! U. Schulzweida, MPI, May 1998, f90 rewrite
  ! A, Rhodin, MPI, Jan 1999, argument list added
  ! M. Esch, MPI, June 1999, ECHAM5-modifications
  ! R. Voss, U.Schlese, December 1999, modifications for coupling
  ! I. Kirchner, MPI, December 2000, time control
  ! U. Schlese, M. Esch, MPI, September 2002, mixed layer ocean
  ! U. Schlese, MPI December 2002, ice thickness over ocean
  !
  ! RongXY Aug 10 2020, modified for coupling with grist 
  !
  ! for more details see file AUTHORS
  !
  IMPLICIT NONE
  !

! Arguments
  INTEGER,  INTENT(IN) :: ncols
  REAL(r8), INTENT(IN) :: delta_time
  REAL(r8):: psiced(ncols),      psni(ncols),      ptsi(ncols),        & 
             ptrfli(ncols),      psofli(ncols),    pahfsi(ncols),    pahfli(ncols)    

! local

  REAL(r8):: zdtime
  REAL(r8):: zalpha, zalphas, zrho_sea, zrho_sn, ziscond, zcpice,       &
             zrhoice, zdice, zcpcon, zcpdt, zsniced, zicefl, zsflx
  INTEGER :: jl
  
!  Executable statements
!
!-- 1. Set up constants
!
  zdtime = delta_time
  zalpha = 2.1656_r8
  zalphas=0.31_r8
  zrho_sea=1025._r8
  zrho_sn=330._r8
  ziscond=zalpha/zalphas*zrho_sea/zrho_sn
  zcpice = 2106._r8
  zrhoice = 910._r8
  zdice = 0.10_r8
  zcpcon = zrhoice*zcpice*zdice
  zcpdt = zcpcon/zdtime
!
!-- 2. Compute new skin-temperature
!    

  DO jl=1,ncols
     IF (psiced(jl).GE.zdice) THEN
        zsniced=psiced(jl)+ziscond*psni(jl)
        zicefl=zalpha*ctfreez/zsniced
        zsflx=ptrfli(jl)+psofli(jl)+pahfsi(jl)+pahfli(jl)
        ptsi(jl)=(zcpdt*ptsi(jl)+zsflx+zicefl)/(zcpdt+zalpha/zsniced)    
        IF (ptsi(jl).GT.tmelt) ptsi(jl)=tmelt
    
     ELSE
        ptsi(jl)=tmelt
     END IF
  END DO

!
  RETURN

END SUBROUTINE grist_seaice_temp

! seaice albedo from cam3

subroutine cam3_albice(pcols,   ncol,    Tair,  snowh   ,coszrs  , icefrac, sicthk, &
                  asdir   ,aldir   ,asdif, aldif   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute surface albedos
!
! Method: 
! Computes surface albedos for direct/diffuse incident radiation for
! two spectral intervals:
!   s = 0.2-0.7 micro-meters
!   l = 0.7-5.0 micro-meters
!
! Albedos specified as follows:
! Ocean with      Surface albs specified; combined with overlying snow
!   sea ice       
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
! 
! Author: CCM1
! 
!-----------------------------------------------------------------------
!
! $Id: albice.F90,v 1.1.4.6 2004/04/23 13:56:25 eaton Exp $
! $Author: eaton $
!
!-----------------------------------------------------------------------

!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid,       only: pcols
!  use dycore, only: dycore_is, get_resolution
!  use comsrf,       only: icefrac, sicthk
!  use grist_constants, only: r8
  use ice_constants,only: c0, c1, p2, p33, rhofresh, rhos, Tffresh, snowpatch, timelt

  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: pcols            ! chunk identifier
  integer , intent(in) :: ncol             ! number of atmospheric columns

  real(r8), intent(in) :: Tair(pcols)      ! bottom level air temp
  real(r8), intent(in) :: snowh(pcols)     ! Snow depth (liquid water equivalent, m? we just input snow depth, so no further conversion is needed)
  real(r8), intent(in) :: coszrs(pcols)    ! Cosine solar zenith angle
  real(r8), intent(in) :: icefrac(pcols)   ! 
  real(r8), intent(in) :: sicthk(pcols)    !
  real(r8), intent(out):: asdir(pcols)     ! Srf alb for direct rad   0.2-0.7 micro-ms
  real(r8), intent(out):: aldir(pcols)     ! Srf alb for direct rad   0.7-5.0 micro-ms
  real(r8), intent(out):: asdif(pcols)     ! Srf alb for diffuse rad  0.2-0.7 micro-ms
  real(r8), intent(out):: aldif(pcols)     ! Srf alb for diffuse rad  0.7-5.0 micro-ms
      ! albedos for ice in each category
  real(r8) :: alvdrn (pcols) ! visible, direct   (fraction)
  real(r8) :: alidrn (pcols) ! near-ir, direct   (fraction)
  real(r8) :: alvdfn (pcols) ! visible, diffuse  (fraction)
  real(r8) :: alidfn (pcols) ! near-ir, diffuse  (fraction)
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i                 ! Longitude index
!-----------------------------------------------------------------------
#ifdef AMIPW_CLIMATE
  real (r8), parameter :: albocn = 0.06_r8  ! ocean albedo
#else
  real (r8), parameter :: albocn = 0.06_r8  ! ocean albedo
#endif
  real (r8), parameter :: &
       ahmax    = 1.0_r8,   &! thickns above which ice alb is const,
#ifdef AMIPW_CLIMATE
       albicev  = 0.68_r8,  &! visible ice albedo for h > ahmax
       albicei  = 0.30_r8,  &! near-ir ice albedo for h > ahmax
#else
       albicev  = 0.68_r8,  &! visible ice albedo for h > ahmax
       albicei  = 0.30_r8,  &! near-ir ice albedo for h > ahmax
#endif
       dT_mlt   = 1._r8,    &! change in temp to give dalb_mlt change
       dalb_mlt = -0.075_r8,&! albedo change per dT_mlt change
       dalb_mltv= -0.100_r8,&! albedo vis change per dT_mlt change in temp for snow
       dalb_mlti= -0.150_r8  ! albedo nir change per dT_mlt change in temp for snow

  real (r8) :: &
       albsnowv,            &! cold snow albedo, visible
       albsnowi              ! cold snow albedo, near IR

! parameter for fractional snow area 
  real(r8)  fhtan ! factor used in dependence of albedo on ice thickness
  real(r8)  vicen(pcols),vsnon(pcols),aicen(pcols),tsfcn(pcols)
  real (r8) hi    ! ice thickness  (m)
  real (r8) hs    ! snow thickness (m)
  real (r8) snw   !
  real (r8) albo  ! effective ocean albedo, function of ice thickness
  real (r8) asnow ! snow-covered area fraction
  real (r8) asnwv ! snow albedo, visible 
  real (r8) asnwi ! snow albedo, near IR
  real (r8) fh    ! piecewise linear function of thickness 
  real (r8) fT    ! piecewise linear function of surface temperature
  real (r8) dTs   ! difference of Tsfc and Timelt
!-----------------------------------------------------------------------

!
! Set snow albedos (every time, since there is no initialization procedure)
!
!  if ( dycore_is ('LR') ) then
#ifdef AMIPW_CLIMATE
     albsnowv = 0.96_r8 ! cold snow albedo, visible
     albsnowi = 0.68_r8 ! cold snow albedo, near IR
#else
     albsnowv = 0.96_r8 ! cold snow albedo, visible
     albsnowi = 0.68_r8 ! cold snow albedo, near IR
#endif
!  else 
!     if ( get_resolution() == 'T31' ) then
!        albsnowv = 0.91_r8 ! cold snow albedo, visible
!        albsnowi = 0.63_r8 ! cold snow albedo, near IR
!     else
!        albsnowv = 0.96_r8 ! cold snow albedo, visible
!        albsnowi = 0.68_r8 ! cold snow albedo, near IR
!     endif
!  endif
  
!
! Initialize all sea ice surface albedos to zero
!
  asdir(:) = 0.
  aldir(:) = 0.
  asdif(:) = 0.
  aldif(:) = 0.
  alvdrn(:) = 0.
  alidrn(:) = 0.
  alvdfn(:) = 0.
  alidfn(:) = 0.

  fhtan = atan(ahmax*5._r8) 

!$omp parallel  private(i,hi,snw,hs,fh,albo,dTs,fT,asnow,asnwv,asnwi) 
!$omp do schedule(static,10) 
  do i=1,ncol
     if (icefrac(i) > 0._r8 .and. coszrs(i)>0.0 ) then !yizhang comment for
     !diagnostic
        hi  = sicthk(i)
        snw = snowh(i) !snowh(i)*rhofresh/rhos ! yizhang
        aicen(i) = icefrac(i)
        vicen(i) = hi*aicen(i)
        !---------------------------------------------------------
        ! keep snow/ice boundary above sea level by reducing snow
        !---------------------------------------------------------
        vsnon(i) = min(snw*aicen(i),p33*vicen(i))
        Tsfcn(i) = min(Tair(i)-Tffresh,-p2)   ! deg C       
        !---------------------------------------------------------
        ! make linear temp profile and compute enthalpy
        !---------------------------------------------------------

        hi = vicen(i) / aicen(i)
        hs = vsnon(i) / aicen(i)
        
        ! bare ice, thickness dependence
        fh = min(atan(hi*5.)/fhtan,c1)
        albo = albocn*(c1-fh)
        alvdfn(i) = albicev*fh + albo
        alidfn(i) = albicei*fh + albo

        ! bare ice, temperature dependence
        dTs = Timelt - Tsfcn(i)
        fT = min(dTs/dT_mlt-c1,c0)
        alvdfn(i) = alvdfn(i) - dalb_mlt*fT
        alidfn(i) = alidfn(i) - dalb_mlt*fT

        if( hs .gt. 0._r8 ) then
           ! fractional area of snow on ice (thickness dependent)
           asnow = hs / ( hs + snowpatch ) 
           asnwv = albsnowv
           asnwi = albsnowi
           ! snow on ice, temperature dependence
           asnwv = asnwv - dalb_mltv*fT
           asnwi = asnwi - dalb_mlti*fT
           
           ! combine ice and snow albedos
           alvdfn(i) = alvdfn(i)*(c1-asnow) + asnwv*asnow
           alidfn(i) = alidfn(i)*(c1-asnow) + asnwi*asnow
        endif
        alvdrn(i) = alvdfn(i)
        alidrn(i) = alidfn(i)
     endif  ! aicen > 0._r8

     if (icefrac(i) > 0._r8 .and. coszrs(i)>0.0 ) then
        asdir(i) = alvdrn(i)
        aldir(i) = alidrn(i)
        asdif(i) = alvdfn(i)
        aldif(i) = alidfn(i)
     end if
  enddo
!$omp end do nowait
!$omp end parallel 
!
  return
end subroutine cam3_albice

END MODULE grist_prescribe_seaice_module
