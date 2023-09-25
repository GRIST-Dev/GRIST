
!==========================================================================
! LiXH use to test physical package.
! This module will be removed when physpkg is completed.
!==========================================================================

 module lixh_test_simple_physics

   use grist_constants,     only: r8, i4, rdry, cp, cv, p00, one,       &
                                  gravity, rvap, rearth, pi, omega,     &
                                  latvap
   use grist_domain_types,  only: global_domain
   !use grist_scm_comm_module, only: tground
   use grist_physics_data_structure, only: tground
   use grist_mpi,            only: mpi_rank
 
   implicit none
   private
   
   public :: large_scale_prcp,      &
             surface_flux,          &
             grist_idealized_physics_mitc

   real(r8)     :: gravit
   real(r8)     :: rair
   real(r8)     :: cpair                       ! Specific heat of dry air
   real(r8)     :: rh2o
   real(r8)     :: epsilo
   real(r8)     :: zvir
   real(r8)     :: a

   real(r8), allocatable ::  Tsurf(:)   ! Sea Surface Temperature (K)
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
   real(r8) pbltop                      ! Pressure position (Pa) of the top of boundary layer
   real(r8) pblconst                    ! Pressure constant (Pa) for the calculation of the decay of diffusivity
 
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

   integer  i,k
   real(r8) dlnpint
   real(r8) qsat
   real(r8) qsats
   real(r8) tmp
   real(r8), allocatable ::  za(:)
   real(r8), allocatable ::  wind(:)
   real(r8), allocatable ::  Cd(:)
   real(r8), allocatable ::  pdel(:,:)
   real(r8), allocatable ::  rpdel(:,:)
   real(r8), allocatable ::  dtdt(:,:)
   real(r8), allocatable ::  dqdt(:,:)

   contains


   subroutine large_scale_prcp(pcols, pver, dtime, ps, pmid, pint,      &
                               t, q, precl)
! io
   integer, intent(in)     :: pcols
   integer, intent(in)     :: pver
   real(r8), intent(in)    :: dtime
   real(r8), intent(in)    :: ps(pcols)
   real(r8), intent(in)    :: pmid(pver,pcols)
   real(r8), intent(in)    :: pint(pver+1,pcols)
   real(r8), intent(inout) :: t(pver,pcols)
   real(r8), intent(inout) :: q(pver,pcols)
   real(r8), intent(out)   :: precl(pcols)
   if(.not.allocated(Tsurf)) allocate(Tsurf(pcols));Tsurf=0.
   if(.not.allocated(za))    allocate(za(pcols));za=0.
   if(.not.allocated(wind))  allocate(wind(pcols));wind=0.
   if(.not.allocated(Cd))    allocate(Cd(pcols));Cd=0.
   if(.not.allocated(pdel))  allocate(pdel(pver,pcols));pdel=0.
   if(.not.allocated(rpdel)) allocate(rpdel(pver,pcols));rpdel=0.
   if(.not.allocated(dtdt))  allocate(dtdt(pver,pcols));dtdt=0.
   if(.not.allocated(dqdt))  allocate(dqdt(pver,pcols));dqdt=0.
  
   gravit = gravity                       ! Gravity (9.80616 m/s^2)
   rair   = rdry                          ! Gas constant for dry air: 287 J/(kg K)
   cpair  = cp                            ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   rh2o   = rvap                          ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                     ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8           ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = rearth                        ! Reference Earth's Radius (m)

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
   t_max       = 294._r8                    ! modified maximum HS equilibrium temperature (original is 315 K)
   delta_T     = 65._r8                     ! difference in eq-polar HS equilibrium temperature (K)
   delta_theta = 10._r8                     ! parameter for vertical temperature gradient (K)

! Large-Scale Condensation and Large-Scale Precipitation Rate

      do k=1,pver
         do i=1,pcols
            pdel(k,i)  = abs(pint(k+1,i)-pint(k,i))                                 ! Pressure thickness (Pa)
            rpdel(k,i) = 1._r8 / pdel(k,i)                                          ! reciprocal of the pressure thickness
!------------
!           Calculate the moisture removal
!------------
            qsat = epsilo*e0/pmid(k,i)*exp(-latvap/rh2o*((1._r8/t(k,i))-1._r8/T0))  ! saturation specific humidity
            if (q(k,i) > qsat) then                                                 ! saturated?
               tmp  = 1._r8/dtime*(q(k,i)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(k,i)**2))) ! condensation rate
               dtdt(k,i) = dtdt(k,i)+latvap/cpair*tmp                               ! temperature tendency
               dqdt(k,i) = dqdt(k,i)-tmp                                            ! specific humidity tendency
               precl(i) = precl(i) + tmp*pdel(k,i)/(gravit*rhow)                    ! precipitation rate
            end if
         end do
      end do
!
! Update moisture and temperature fields from Large-Scale Precipitation Scheme
!
if(.false.)then
      do k=1,pver
         do i=1,pcols
            ! LiXH test qv before large scale condensation
             if(q(k,i) .lt. 0._r8)then
                print*,'Negetive qv appears! rank =', mpi_rank(), 'k=',k,'q=',q(k,i)
             end if
           
            t(k,i) =  t(k,i) + dtdt(k,i)*dtime    ! update the state variable T
            q(k,i) =  q(k,i) + dqdt(k,i)*dtime    ! update the state variable q

            ! LiXH add a simple positive fixer
            if(q(k,i) .lt. 0._r8)then
                q(k,i) =  q(k,i) - dqdt(k,i)*dtime
                t(k,i) =  t(k,i) - dtdt(k,i)*dtime
            end if
         end do
      end do
end if
   end subroutine large_scale_prcp

   subroutine surface_flux(pcols, pver, pmid, pint, u, v, q, t, ps,      &
                           taux, tauy, shflx, qflx)
   !use grist_scm_comm_module,              only: have_lhflx, have_shflx
   use grist_physics_data_structure,        only: have_lhflx, have_shflx
! io
   integer, intent(in)     :: pcols
   integer, intent(in)     :: pver
   real(r8), intent(in)    :: pmid(pver,pcols)
   real(r8), intent(in)    :: pint(pver+1,pcols)
   real(r8), intent(in)    :: u(pver,pcols)
   real(r8), intent(in)    :: v(pver,pcols)
   real(r8), intent(in)    :: t(pver,pcols)
   real(r8), intent(in)    :: q(pver,pcols)
   real(r8), intent(in)    :: ps(pcols)
   real(r8), intent(out)   :: taux(pcols), tauy(pcols), shflx(pcols), qflx(pcols)
! local
   real(r8) rho

! Set Surface Temperature
     do i=1,pcols
        Tsurf(i) = tground
     end do
! Calculate hydrostatic height za of the lowest model level
     do i=1,pcols 
        dlnpint = log(ps(i)) - log(pint(pver,i))  ! ps(i) is identical to pint(pver+1,i), 
                                                  ! note: this gives the correct positive sign and corrects
                                                  ! the typo in the Reed and Jablonowski (2012) paper
        za(i) = rair/gravit*t(pver,i)*(1._r8+zvir*q(pver,i))*0.5_r8*dlnpint
     end do

     do i = 1,pcols
        wind(i) = sqrt(u(pver,i)**2+v(pver,i)**2)
        if( wind(i) .lt. v20) then
            Cd(i) = Cd0+Cd1*wind(i)
        else
            Cd(i) = Cm
        end if
        rho   = pmid(pver,i)/(rair * t(pver,i) * (1._r8+zvir*q(pver,i)))
        qsats = epsilo*e0/ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))
        taux(i) = -rho * Cd(i) * wind(i) * u(pver,i)
        tauy(i) = -rho * Cd(i) * wind(i) * v(pver,i)

        if(.not.(have_shflx .and. have_lhflx))then
            shflx(i)     = rho * cpair  * C*wind(i)*(Tsurf(i)-t(pver,i))
            qflx(i)      = rho * C*wind(i)*(qsats-q(pver,i))
        end if
     end do

   if(allocated(Tsurf)) deallocate(Tsurf)
   if(allocated(za))    deallocate(za)
   if(allocated(wind))  deallocate(wind)
   if(allocated(Cd))    deallocate(Cd)
   if(allocated(pdel))  deallocate(pdel)
   if(allocated(rpdel)) deallocate(rpdel)
   if(allocated(dtdt))  deallocate(dtdt)
   if(allocated(dqdt))  deallocate(dqdt)
 
   end subroutine surface_flux


subroutine grist_idealized_physics_mitc(pcols, pver, dtime, lat, t, q, u, v,        &
                                        pmid, pint, etamid, ps, precl, mode)
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
   real(r8), intent(in) :: pmid(pver,pcols)   ! Pressure at full-model level (Pa)
   real(r8), intent(in) :: pint(pver+1,pcols) ! Pressure at model interfaces (Pa), position pver+1 is the surface
   real(r8), intent(in) :: ps(pcols)          ! Surface Pressue (Pa)

!
! Lixh add a 'mode' variable for SCM:
! if (mode .eq. 'scm'), then use the surface temperature of the forcing data.
!
   character(*),intent(in) :: mode

! Input/Output arguments 
!
!  pcols is the maximum number of vertical columns per 'chunk' of atmosphere
!
   real(r8), intent(inout) :: t(pver,pcols)      ! Temperature at full-model level (K)
   real(r8), intent(inout) :: q(pver,pcols)      ! Specific Humidity at full-model level (kg/kg)
   real(r8), intent(inout) :: u(pver,pcols)      ! Zonal wind at full-model level (m/s)
   real(r8), intent(inout) :: v(pver,pcols)      ! Meridional wind at full-model level (m/s)

!
! Output arguments 
!
   real(r8), intent(out) :: precl(pcols)         ! Large-scale precipitation rate (m_water / s)

!
!---------------------------Local workspace-----------------------------
!

! Integers for loops

   integer  i,k                         ! Longitude and level indices

! Physical Constants, all in SI units - Many of these may be model dependent

   real(r8) gravit                      ! Gravity
   real(r8) rair                        ! Gas constant for dry air
   real(r8) cpair                       ! Specific heat of dry air
   real(r8) rh2o                        ! Gas constant for water vapor
   real(r8) epsilo                      ! Ratio of gas constant for dry air to that for vapor
   real(r8) zvir                        ! Constant for virtual temp. calc. =(rh2o/rair) - 1
   real(r8) a                           ! Reference Earth's Radius (m)

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

   real(r8) dtdt(pver,pcols)            ! Temperature tendency
   real(r8) dqdt(pver,pcols)            ! Specific humidity tendency
   real(r8) dudt(pver,pcols)            ! Zonal wind tendency
   real(r8) dvdt(pver,pcols)            ! Meridional wind tendency

! Surface fluxes for diagnostic purposes

   real(r8) taux(pcols)                 ! surface momentum flux in the longitudinal direction
   real(r8) tauy(pcols)                 ! surface momentum flux in the latitudinal direction
   real(r8) shflx(pcols)                ! sensible heat flux at the surface
   real(r8) lhflx(pcols)                ! latent heat flux at the surface

! Temporary variables and arrays

   real(r8) tmp                         ! holds coeffiecients
   real(r8) qsat                        ! Saturation vapor pressure (Pa)
   real(r8) qsats                       ! Saturation vapor pressure (Pa) at the surface with prescribed SST
   real(r8) pdel(pver,pcols)            ! Pressure thickness (Pa) of the layer, pdel is the positive
                                        ! difference between the surrounding model interface pressures
   real(r8) rpdel(pver,pcols)           ! Reciprocal of the pressure layer thickness (1/Pa)


! Variables for Boundary Layer Calculation

   real(r8) wind(pcols)                 ! Magnitude of Wind
   real(r8) Cd(pcols)                   ! Drag coefficient for momentum
   real(r8) Km(pver+1,pcols)            ! Eddy diffusivity for boundary layer calculations: u and v
   real(r8) Ke(pver+1,pcols)            ! Eddy diffusivity for boundary layer calculations: T and q
   real(r8) rho                         ! Moist density
   real(r8) za(pcols)                   ! Heights at midpoints of first model level (m)
   real(r8) dlnpint                     ! Used for calculation of heights
   real(r8) pbltop                      ! Pressure position (Pa) of the top of boundary layer
   real(r8) pblconst                    ! Pressure constant (Pa) for the calculation of the decay of diffusivity
   real(r8) CA(pver,pcols)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CC(pver,pcols)              ! Matrix Coefficents for PBL Scheme 
   real(r8) CE(pver+1,pcols)            ! Matrix Coefficents for PBL Scheme
   real(r8) CAm(pver,pcols)             ! Matrix Coefficents for PBL Scheme
   real(r8) CCm(pver,pcols)             ! Matrix Coefficents for PBL Scheme
   real(r8) CEm(pver+1,pcols)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFu(pver+1,pcols)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFv(pver+1,pcols)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFt(pver+1,pcols)           ! Matrix Coefficents for PBL Scheme
   real(r8) CFq(pver+1,pcols)           ! Matrix Coefficents for PBL Scheme


! Save the initial specific humidity contents to be able to compute the moist pressure
! adjustments at the end of the routine, whether this is needed depends on the model

   real(r8) qini(pver,pcols)            ! Initial specific humidity

!===============================================================================
!
! Physical Constants - MAY BE MODEL DEPENDENT 
!
!===============================================================================
!   gravit = 9.80616_r8                   ! Gravity (9.80616 m/s^2)
!   rair   = 287.0_r8                     ! Gas constant for dry air: 287 J/(kg K)
!   cpair  = 1.0045e3_r8                  ! Specific heat of dry air: here we use 1004.5 J/(kg K)
!   latvap = 2.5e6_r8                     ! Latent heat of vaporization (J/kg)
!   rh2o   = 461.5_r8                     ! Gas constant for water vapor: 461.5 J/(kg K)
!   epsilo = rair/rh2o                    ! Ratio of gas constant for dry air to that for vapor
!   zvir   = (rh2o/rair) - 1._r8          ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
!   a      = 6371220.0_r8                 ! Reference Earth's Radius (m)
!   omega  = 7.29212d-5                   ! Reference rotation rate of the Earth (s^-1)
!   pi     = 4._r8*atan(1._r8)            ! pi

   gravit = gravity                       ! Gravity (9.80616 m/s^2)
   rair   = rdry                          ! Gas constant for dry air: 287 J/(kg K)
   cpair  = cp                            ! Specific heat of dry air: here we use 1004.5 J/(kg K)
   rh2o   = rvap                          ! Gas constant for water vapor: 461.5 J/(kg K)
   epsilo = rair/rh2o                     ! Ratio of gas constant for dry air to that for vapor
   zvir   = (rh2o/rair) - 1._r8           ! Constant for virtual temp. calc. =(rh2o/rair) - 1 is approx. 0.608
   a      = rearth                        ! Reference Earth's Radius (m)

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
   t_max       = 294._r8                    ! modified maximum HS equilibrium temperature (original is 315 K)
   delta_T     = 65._r8                     ! difference in eq-polar HS equilibrium temperature (K)
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
        dlnpint = log(ps(i)) - log(pint(pver,i))  ! ps(i) is identical to pint(pver+1,i), 
                                                  ! note: this gives the correct positive sign and corrects
                                                  ! the typo in the Reed and Jablonowski (2012) paper
        za(i) = rair/gravit*t(pver,i)*(1._r8+zvir*q(pver,i))*0.5_r8*dlnpint
     end do
!
! Set Initial Specific Humidity - For dry mass adjustment at the end
!
     qini(:pver,:pcols) = q(:pver,:pcols)

!--------------------------------------------------------------
! Set Sea Surface Temperature
!--------------------------------------------------------------
     if(trim(mode) .eq. 'scm')then
         do i=1,pcols
            Tsurf(i) = tground
         end do
     else
        do i=1,pcols
            Tsurf(i) = del_T*exp(-(((lat(i))**2.0_r8)/(2.0_r8*(T_width**2.0_r8)))) + T_min
        end do
     end if

!===============================================================================
!
! Set initial physics time tendencies and precipitation field to zero
!
!===============================================================================
     dtdt(:pver,:pcols)  = 0._r8            ! initialize temperature tendency with zero
     dqdt(:pver,:pcols)  = 0._r8            ! initialize specific humidity tendency with zero
     dudt(:pver,:pcols)  = 0._r8            ! initialize zonal wind tendency with zero
     dvdt(:pver,:pcols)  = 0._r8            ! initialize meridional wind tendency with zero
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
            pdel(k,i)  = abs(pint(k+1,i)-pint(k,i))                                 ! Pressure thickness (Pa)
            rpdel(k,i) = 1._r8 / pdel(k,i)                                          ! reciprocal of the pressure thickness
!------------
!           Calculate the moisture removal
!------------
            qsat = epsilo*e0/pmid(k,i)*exp(-latvap/rh2o*((1._r8/t(k,i))-1._r8/T0))  ! saturation specific humidity
            if (q(k,i) > qsat) then                                                 ! saturated?
               tmp  = 1._r8/dtime*(q(k,i)-qsat)/(1._r8+(latvap/cpair)*(epsilo*latvap*qsat/(rair*t(k,i)**2))) ! condensation rate
               dtdt(k,i) = dtdt(k,i)+latvap/cpair*tmp                               ! temperature tendency
               dqdt(k,i) = dqdt(k,i)-tmp                                            ! specific humidity tendency
               precl(i) = precl(i) + tmp*pdel(k,i)/(gravit*rhow)                    ! precipitation rate
            end if
         end do
      end do
!
! Update moisture and temperature fields from Large-Scale Precipitation Scheme
!
      do k=1,pver
         do i=1,pcols
            t(k,i) =  t(k,i) + dtdt(k,i)*dtime    ! update the state variable T
            q(k,i) =  q(k,i) + dqdt(k,i)*dtime    ! update the state variable q
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

     do i=1,pcols
        wind(i) = sqrt(u(pver,i)**2+v(pver,i)**2)    ! wind magnitude at the lowest level
     end do
     do i=1,pcols
        Ke(pver+1,i) = C*wind(i)*za(i)               ! Ke: latent and sensible heat mixing coefficient
        if( wind(i) .lt. v20) then                   ! Km: momentum mixing coefficient
           Cd(i) = Cd0+Cd1*wind(i) 
           Km(pver+1,i) = Cd(i)*wind(i)*za(i)
        else
           Cd(i) = Cm
           Km(pver+1,i) = Cm*wind(i)*za(i)
        endif
     end do

      do k=1,pver                                     ! adjustments in the vertical direction
         do i=1,pcols
            if( pint(k,i) .ge. pbltop) then
               Km(k,i) = Km(pver+1,i)                 ! constant Km below 850 hPa level
               Ke(k,i) = Ke(pver+1,i)                 ! constant Ke below 850 hPa level
            else
               Km(k,i) = Km(pver+1,i)*exp(-(pbltop-pint(k,i))**2/(pblconst)**2)  ! Km tapered to 0
               Ke(k,i) = Ke(pver+1,i)*exp(-(pbltop-pint(k,i))**2/(pblconst)**2)  ! Ke tapered to 0
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

     do i=1,pcols
        rho   = pmid(pver,i)/(rair * t(pver,i) * (1._r8+zvir*q(pver,i)))      ! moist density at lowest model level rho = p/(Rd Tv)
        qsats = epsilo*e0/ps(i)*exp(-latvap/rh2o*((1._r8/Tsurf(i))-1._r8/T0))      ! saturation specific humidity at the surface
!---------
!       Surface momentum fluxes are replaced by Rayleigh friction, commented out
!
!       dudt(pver,i) = dudt(pver,i) + (u(pver,i) &                                ! replaced by Rayleigh friction below
!                           /(1._r8+Cd(i)*wind(i)*dtime/za(i))-u(pver,i))/dtime
!       dvdt(pver,i) = dvdt(pver,i) + (v(pver,i) &                                ! replaced by Rayleigh friction below
!                          /(1._r8+Cd(i)*wind(i)*dtime/za(i))-v(pver,i))/dtime
!       taux(i)     = -rho * Cd(i) * wind(i) * u(pver,i)                          ! surface momentum flux x in N/m2
!       tauy(i)     = -rho * Cd(i) * wind(i) * v(pver,i)                          ! surface momentum flux y in N/m2
!       u(pver,i)   = u(pver,i)/(1._r8+Cd(i)*wind(i)*dtime/za(i))                 ! update of the state variable u
!       v(pver,i)   = v(pver,i)/(1._r8+Cd(i)*wind(i)*dtime/za(i))                 ! update of the state variable v
!--------

!---------------------------
!       Sensible heat fluxes
!---------------------------
        dtdt(pver,i) = dtdt(pver,i) +((t(pver,i)+C*wind(i)*Tsurf(i)*dtime/za(i)) & ! tendency due to sensible heat flux
                            /(1._r8+C*wind(i)*dtime/za(i))-t(pver,i))/dtime
        shflx(i)     = rho * cpair  * C*wind(i)*(Tsurf(i)-t(pver,i))               ! sensible heat flux in W/m2
        t(pver,i)    = (t(pver,i)+C*wind(i)*Tsurf(i)*dtime/za(i)) &                ! update of the state variable t
                            /(1._r8+C*wind(i)*dtime/za(i))  
!-------------------------
!       Latent heat fluxes
!-------------------------
        dqdt(pver,i) = dqdt(pver,i) +((q(pver,i)+C*wind(i)*qsats*dtime/za(i)) &    ! tendency due to latent heat flux
                            /(1._r8+C*wind(i)*dtime/za(i))-q(pver,i))/dtime
        lhflx(i)     = rho * latvap * C*wind(i)*(qsats-q(pver,i))                  ! latent heat flux in W/m2
        q(pver,i)    = (q(pver,i)+C*wind(i)*qsats*dtime/za(i)) &                   ! update of the state variable q
                            /(1._r8+C*wind(i)*dtime/za(i))
     end do
!===============================================================================


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
            rho = (pint(k+1,i)/(rair*(t(k+1,i)*(1._r8+zvir*q(k+1,i))+t(k,i)*(1._r8+zvir*q(k,i)))/2.0_r8))! density of moist air
            CAm(k,i)   = rpdel(k,i)*dtime*gravit*gravit*Km(k+1,i)*rho*rho   &
                         /(pmid(k+1,i)-pmid(k,i))    
            CCm(k+1,i) = rpdel(k+1,i)*dtime*gravit*gravit*Km(k+1,i)*rho*rho &
                         /(pmid(k+1,i)-pmid(k,i))
            CA(k,i)    = rpdel(k,i)*dtime*gravit*gravit*Ke(k+1,i)*rho*rho   &
                         /(pmid(k+1,i)-pmid(k,i))
            CC(k+1,i)  = rpdel(k+1,i)*dtime*gravit*gravit*Ke(k+1,i)*rho*rho &
                         /(pmid(k+1,i)-pmid(k,i))
         end do
      end do
      do i=1,pcols
         CAm(pver,i)   = 0._r8
         CCm(1,i)      = 0._r8
         CEm(pver+1,i) = 0._r8
         CA(pver,i)    = 0._r8
         CC(1,i)       = 0._r8
         CE(pver+1,i)  = 0._r8
         CFu(pver+1,i) = 0._r8
         CFv(pver+1,i) = 0._r8
         CFt(pver+1,i) = 0._r8
         CFq(pver+1,i) = 0._r8 
      end do
      do i=1,pcols
         do k=pver,1,-1
            CE(k,i)  = CC(k,i)/(1._r8+CA(k,i)+CC(k,i)-CA(k,i)*CE(k+1,i)) 
            CEm(k,i) = CCm(k,i)/(1._r8+CAm(k,i)+CCm(k,i)-CAm(k,i)*CEm(k+1,i))
            CFu(k,i) = (u(k,i)+CAm(k,i)*CFu(k+1,i)) &
                       /(1._r8+CAm(k,i)+CCm(k,i)-CAm(k,i)*CEm(k+1,i))
            CFv(k,i) = (v(k,i)+CAm(k,i)*CFv(k+1,i)) &
                       /(1._r8+CAm(k,i)+CCm(k,i)-CAm(k,i)*CEm(k+1,i))
            CFt(k,i) = ((p0/pmid(k,i))**(rair/cpair)*t(k,i)+CA(k,i)*CFt(k+1,i)) &
                       /(1._r8+CA(k,i)+CC(k,i)-CA(k,i)*CE(k+1,i)) 
            CFq(k,i) = (q(k,i)+CA(k,i)*CFq(k+1,i)) &
                       /(1._r8+CA(k,i)+CC(k,i)-CA(k,i)*CE(k+1,i))
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
!            dudt(1,i)  = dudt(1,i)+(CFu(1,i)-u(1,i))/dtime
!            dvdt(1,i)  = dvdt(1,i)+(CFv(1,i)-v(1,i))/dtime
!            u(1,i)    = CFu(1,i)
!            v(1,i)    = CFv(1,i)
!
            dtdt(1,i)  = dtdt(1,i)+(CFt(1,i)*(pmid(1,i)/p0)**(rair/cpair)-t(1,i))/dtime
            t(1,i)    = CFt(1,i)*(pmid(1,i)/p0)**(rair/cpair)
            dqdt(1,i)  = dqdt(1,i)+(CFq(1,i)-q(1,i))/dtime
            q(1,i)  = CFq(1,i)
      end do
!
! Loop over the remaining level
!
      do i=1,pcols
         do k=2,pver
!
!            Mixing of momentum is done via the Rayleigh friction
!
!            dudt(k,i)  = dudt(k,i)+(CEm(k,i)*u(k-1,i)+CFu(k,i)-u(k,i))/dtime
!            dvdt(k,i)  = dvdt(k,i)+(CEm(k,i)*v(k-1,i)+CFv(k,i)-v(k,i))/dtime
!            u(k,i)    = CEm(k,i)*u(k-1,i)+CFu(k,i)
!            v(k,i)    = CEm(k,i)*v(k-1,i)+CFv(k,i)
!
            dtdt(k,i)  = dtdt(k,i)+((CE(k,i)*t(k-1,i) &
                              *(p0/pmid(k-1,i))**(rair/cpair)+CFt(k,i)) &
                              *(pmid(k,i)/p0)**(rair/cpair)-t(k,i))/dtime 
            t(k,i)    = (CE(k,i)*t(k-1,i)*(p0/pmid(k-1,i))**(rair/cpair)+CFt(k,i)) &
                              *(pmid(k,i)/p0)**(rair/cpair)
            dqdt(k,i)  = dqdt(k,i)+(CE(k,i)*q(k-1,i)+CFq(k,i)-q(k,i))/dtime
            q(k,i)  = CE(k,i)*q(k-1,i)+CFq(k,i)
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
               dudt(k,i) = dudt(k,i) + tmp*u(k,i)          ! compute u tendency
               dvdt(k,i) = dvdt(k,i) + tmp*v(k,i)          ! compute v tendency
               u(k,i)    = u(k,i) +  tmp*u(k,i)*dtime      ! update u
               v(k,i)    = v(k,i) +  tmp*v(k,i)*dtime      ! update v
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
               trefa = (trefc - delta_theta*cos(lat(i))**2*log((pmid(k,i)/ps0)))*(pmid(k,i)/ps0)**(rair/cpair)
               trefa = max(t00,trefa)
               dtdt(k,i) = dtdt(k,i) + (trefa - t(k,i))*tmp             ! compute T tendency
               t(k,i)    = t(k,i) + (trefa - t(k,i))*tmp*dtime          ! update T
            end do
         else
            tmp = ka
            do i=1,pcols
               trefc = T_max - delta_T*sin(lat(i))**2
               trefa = (trefc - delta_theta*cos(lat(i))**2*log((pmid(k,i)/ps0)))*(pmid(k,i)/ps0)**(rair/cpair)
               trefa = max(t00,trefa)
               dtdt(k,i) = dtdt(k,i) + (trefa - t(k,i))*tmp             ! compute T tendency
               t(k,i)    = t(k,i) + (trefa - t(k,i))*tmp*dtime          ! update T
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



end module lixh_test_simple_physics
