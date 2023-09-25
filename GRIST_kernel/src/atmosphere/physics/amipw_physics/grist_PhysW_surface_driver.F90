
!----------------------------------------------------------------------------
! Created on 2019
! Author: Yi Zhang
! Version 1.0
! Description: This driver mainly follows the naming convention of WRF, but 
!              the dim conversion requires a transpose from GRIST based vars, 
! Revision history: 
!              1. sfclay from v3.8.1 (default)
!--------------------------------------------------------------------------

module grist_PhysW_surface_driver

  use grist_constants,         only: i4, r8, r_d=>rdry, r_v=>rvap, cp, gravity,  zero, latvap
  use grist_nml_module,        only: working_mode
  use grist_mpi
#ifdef USE_SCM
! why here?
  use grist_scm_comm_module
#endif
  use grist_PhysW_nml_module,       only: PhysW_sf_scheme
  use grist_PhysW_data_structure,   only: pstate_wrf, psurf_wrf, num_soil_layers, p_qv
  use grist_PhysW_pbl_driver,       only: stepbl
  use grist_PhysW_radiation_driver, only: call_radiation
  use module_wrfmodel_constants,    only: rcp, rovg, xls, xlv, xlf, &
                                          ep_1, ep_2, svp1, svp2, svp3, svpt0, &
                                          karman, eomeg, stbolt, p1000mb
  use module_sf_sfclay_v381,        only: sfclayinit_v381=>sfclayinit, sfclay_v381=>sfclay

  implicit  none

  private
  public   ::   grist_PhysW_surface_init, &
                grist_PhysW_surface_run,  &
                grist_PhysW_surface_final

    integer(i4)  :: ifsnow
    integer,parameter,private :: isfflx   = 1   !=1 for surface heat and moisture fluxes.
    integer,parameter,private :: isftcflx = 0   !=0,(Charnock and Carlson-Boland).
    integer,parameter,private :: iz0tlnd  = 0   !=0,(Carlson-Boland).

contains
    
    subroutine grist_PhysW_surface_init

       ifsnow = 1
 
       select case(trim(PhysW_sf_scheme))
       case('SFCLAYV381')
          call sfclayinit_v381(.false.)
       case default
          if(mpi_rank().eq.0) print*, "you must select a SFC scheme: SFCLAY, SFCLAYV341, stop"
          call mpi_abort()
       end select
       return

    end subroutine grist_PhysW_surface_init

    subroutine grist_PhysW_surface_final
       return
    end subroutine grist_PhysW_surface_final

    subroutine grist_PhysW_surface_run(ncell,nLevel,nspecies,itimestep,dtime,dxmean)

!------------------------------------------------------------------
!  implicit none
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
! theta      potential temperature (k)
! qv         water vapor mixing ratio (kg/kg)
! qc         cloud water mixing ratio (kg/kg)
! qr         rain water mixing ratio (kg/kg)
! qi         cloud ice mixing ratio (kg/kg)
! qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- itimestep     number of time steps
!-- glw           downward long wave flux at ground surface (w/m^2)
!-- gsw           downward short wave flux at ground surface (w/m^2)
!-- emiss         surface emissivity (between 0 and 1)
!-- tsk           surface temperature (k)
!-- tmn           soil temperature at lower boundary (k)
!-- xland         land mask (1 for land, 2 for water)
!-- znt           time-varying roughness length (m)
!-- z0            background roughness length (m)
!-- mavail        surface moisture availability (between 0 and 1)
!-- ust           u* in similarity theory (m/s)
!-- mol           q* (similarity theory) (kg/kg)
!-- hol           pbl height over monin-obukhov length
!-- pblh          pbl height (m)
!-- capg          heat capacity for soil (j/k/m^3)
!-- thc           thermal inertia (cal/cm/k/s^0.5)
!-- snowc         flag indicating snow coverage (1 for snow cover)
!-- hfx           net upward heat flux at the surface (w/m^2)
!-- qfx           net upward moisture flux at the surface (kg/m^2/s)
!-- lh            net upward latent heat flux at surface (w/m^2)
!-- regime        flag indicating pbl regime (stable, unstable, etc.)
!-- tke_myj       turbulence kinetic energy from mellor-yamada-janjic (myj) (m^2/s^2)
!-- akhs          sfc exchange coefficient of heat/moisture from myj
!-- akms          sfc exchange coefficient of momentum from myj
!-- thz0          potential temperature at roughness length (k)
!-- uz0           u wind component at roughness length (m/s)
!-- vz0           v wind component at roughness length (m/s)
!-- qsfc          specific humidity at lower boundary (kg/kg)
!-- u10           diagnostic 10-m u component from surface layer
!-- v10           diagnostic 10-m v component from surface layer
!-- th2           diagnostic 2-m theta from surface layer and lsm
!-- t2            diagnostic 2-m temperature from surface layer and lsm
!-- q2            diagnostic 2-m mixing ratio from surface layer and lsm
!-- tshltr        diagnostic 2-m theta from myj
!-- th10          diagnostic 10-m theta from myj
!-- qshltr        diagnostic 2-m specific humidity from myj
!-- q10           diagnostic 10-m specific humidity from myj
!-- lowlyr        index of lowest model layer above ground
!-- rr            dry air density (kg/m^3)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (k)
!-- moist         moisture array (4d - last index is species) (kg/kg)
!-- p_phy         pressure (pa)
!-- pi_phy        exner function (dimensionless)
!-- pshltr        diagnostic shelter (2m) pressure from myj (pa)
!-- p8w           pressure at full levels (pa)
!-- t_phy         temperature (k)
!-- dz8w          dz between full levels (m)
!-- z             height above sea level (m)
!-- config_flags
!-- dx            horizontal space interval (m)
!-- dt            time step (second)
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
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!
!------------------------------------------------------------------ 
! io
   integer(i4), intent(in)    :: ncell
   integer(i4), intent(in)    :: nLevel
   integer(i4), intent(in)    :: nspecies
   integer(i4), intent(in)    :: itimestep
   real(r8),    intent(in)    :: dtime
   real(r8),    intent(in)    :: dxmean
! local       
   real(r8)   :: qgh (1:ncell, 1:1)
   real(r8)   :: chs (1:ncell, 1:1)
   real(r8)   :: chs2(1:ncell, 1:1)
   real(r8)   :: cqs2(1:ncell, 1:1)
   real(r8)   :: cpm (1:ncell, 1:1)
   real(r8)   :: fm(  1:ncell, 1:1)
   real(r8)   :: fh(  1:ncell, 1:1)
   real(r8)   :: rmol(1:ncell, 1:1)
   real(r8)   :: dtmin
   real(r8)   :: dtbl

!
! surface schemes need pbl time step for updates and accumulations
! assume these schemes provide no tendencies
!
     dtmin = dtime/60._r8
     dtbl  = dtime*stepbl

     select case(trim(PhysW_sf_scheme))

     case('SFCLAYV381')

        !pstate_wrf%dxmean(1:ncell) = dxmean
     call sfclay_v381(&
                 u3d    = pstate_wrf%u_phy(1:ncell,1:nLevel,1:1),     & ! in
                 v3d    = pstate_wrf%v_phy(1:ncell,1:nLevel,1:1),     & ! in
                 t3d    = pstate_wrf%t_phy(1:ncell,1:nLevel,1:1),     & ! in
                 qv3d   = pstate_wrf%moist(1:ncell,1:nLevel,1:1,p_qv),& ! in
                 p3d    = pstate_wrf%p_phy(1:ncell,1:nLevel,1:1),     & ! in
                 dz8w   = pstate_wrf%dz8w (1:ncell,1:nLevel,1:1),     & ! in
                 cp     = cp,       &
                 g      = gravity,  &
                 rovcp  = rcp,      &
                 r      = r_d,      &
                 xlv    = xlv,      &
                 psfc   = psurf_wrf%psfc(1:ncell,1:1),   & ! in
                 chs    = chs(1:ncell,1:1),              & ! inout, 1st calc inside
                 chs2   = chs2(1:ncell,1:1),             & ! inout, 1st calc inside  
                 cqs2   = cqs2(1:ncell,1:1),             & ! inout, 1st calc inside
                 cpm    = cpm(1:ncell,1:1),              & ! inout, 1st calc inside
                 znt    = psurf_wrf%znt(   1:ncell,1:1), & ! inout
                 ust    = psurf_wrf%ust(   1:ncell,1:1), & ! inout
                 pblh   = psurf_wrf%pblh(  1:ncell,1:1), & ! in
                 mavail = psurf_wrf%mavail(1:ncell,1:1), & ! in
                 zol    = psurf_wrf%zol(   1:ncell,1:1), & ! inout, 1st-calc-inside, z/l height over monin-obukhov length
                 mol    = psurf_wrf%mol(   1:ncell,1:1), & ! inout
                 regime = psurf_wrf%regime(1:ncell,1:1), & ! inout, 1st-calc-inside
                 psim   = psurf_wrf%psim(  1:ncell,1:1), & ! inout, 1st-calc-inside
                 psih   = psurf_wrf%psih(  1:ncell,1:1), & ! inout, 1st-calc-inside
                 fm     = fm(  1:ncell,1:1),             & ! inout, 1st-calc-inside
                 fh     = fh(  1:ncell,1:1),             & ! inout, 1st-calc-inside
                 xland  = psurf_wrf%xland(1:ncell,1:1),  & ! in
                 hfx    = psurf_wrf%hfx(  1:ncell,1:1),  & ! inout, 1st calc inside
                 qfx    = psurf_wrf%qfx(  1:ncell,1:1),  & ! inout, 1st calc inside
                 lh     = psurf_wrf%lh(   1:ncell,1:1),  & ! inout, 1st calc inside
                 tsk    = psurf_wrf%tsk(  1:ncell,1:1),  & ! in
                 flhc   = psurf_wrf%flhc( 1:ncell,1:1),  & ! inout, 1st calc inside
                 flqc   = psurf_wrf%flqc( 1:ncell,1:1),  & ! inout, 1st calc inside
                 qgh    = qgh(1:ncell,1:1),              & ! inout, 1st calc inside
                 qsfc   = psurf_wrf%qsfc(  1:ncell,1:1), & ! out, only water point
                 rmol   = rmol(  1:ncell,1:1)          , & ! inout, 1st calc inside
                 u10    = psurf_wrf%u10(   1:ncell,1:1), & ! out
                 v10    = psurf_wrf%v10(   1:ncell,1:1), & ! out
                 th2    = psurf_wrf%th2(   1:ncell,1:1), & ! out
                 t2     = psurf_wrf%t2(    1:ncell,1:1), & ! out
                 q2     = psurf_wrf%q2(    1:ncell,1:1), & ! out
                 gz1oz0 = psurf_wrf%gz1oz0(1:ncell,1:1), & ! inout, 1st calc inside
                 wspd   = psurf_wrf%wspd(  1:ncell,1:1), & ! inout, 1st calc inside
                 br     = psurf_wrf%br(    1:ncell,1:1), & ! inout, 1st calc inside
                 isfflx = isfflx, &
                 dx     = pstate_wrf%dxmean(1:ncell), &
                 svp1   = svp1,   &
                 svp2   = svp2,   &
                 svp3   = svp3,   &
                 svpt0  = svpt0,  &
                 ep1    = ep_1,   &
                 ep2    = ep_2,   &
                 karman = karman, &
                 eomeg  = eomeg,  &
                 stbolt = stbolt, &
                 p1000mb= p1000mb,&
                 ids=1,ide=ncell, jds=1,jde=1, kds=1,kde=nLevel,  &
                 ims=1,ime=ncell, jms=1,jme=1, kms=1,kme=nLevel,  &
                 its=1,ite=ncell, jts=1,jte=1, kts=1,kte=nLevel-1, &
                 isftcflx =isftcflx, iz0tlnd=iz0tlnd)
                 !ustm,ck,cka,cd,cda,isftcflx,iz0tlnd,scm_force_flux ! not used

     case default
          if(mpi_rank().eq.0) print*, "you must select a SFC scheme: SFCLAY, SFCLAYV341, SFCLAYV381, stop"
          call mpi_abort()
     end select

   return
   end subroutine grist_PhysW_surface_run

end module grist_PhysW_surface_driver
