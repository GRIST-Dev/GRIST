!======================================================
!
!  Created by LiXiaohan on 19/5/13.
!  interface to RRTMG solar radiation
!======================================================

 module radsw

    use grist_constants,                    only: i4, r8
    use grist_nml_module,                   only: nlev, nlevp
    use grist_handle_error,                 only: endrun
    use parrrsw,                            only: nbndsw, ngptsw
    use rrtmg_sw_init,                      only: rrtmg_sw_ini
    use rrtmg_sw_rad,                       only: rrtmg_sw, rrtmg_4DDA_sw
    use radconstants,                       only: idx_sw_diag
!    use grist_mpi

    implicit none
    private
    save

    real(r8) :: fractional_solar_irradiance(1:nbndsw) ! fraction of solar irradiance in each band
    real(r8) :: solar_band_irrad(1:nbndsw) ! rrtmg-assumed solar irradiance in each sw band

    public :: radsw_init ,     &! initialize constants
              rad_rrtmg_sw      ! driver for solar radiation code


 contains

! Purpose: Solar radiation code
! Method:  RRTMG, two-stream, with McICA
! Divides solar spectrum into 14 intervals from 0.2-12.2 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  Includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud, 
! and surface absorption. Computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. Adds the layers 
! assuming scattering between layers to be isotropic, and distinguishes 
! direct solar beam from scattered radiation.
! 
! Longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. coszrs > 0) computations are done.
! 
! Note that an extra layer above the model top layer is added.
! 
! mks units are used.
! 
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
    subroutine rad_rrtmg_sw(ncol     ,rrtmg_levs ,r_state      ,                    &
                            E_pmid   ,E_cld      ,                                  &
                            E_aer_tau,E_aer_tau_w,E_aer_tau_w_g,E_aer_tau_w_f,      &
                            eccf     ,E_coszrs   ,solin        ,sfac         ,      &
                            E_asdir  ,E_asdif    ,E_aldir      ,E_aldif      ,      &
                            qrs      ,qrsc       ,fsnt         ,fsntc        ,      &
                            fsntoa   ,fsutoa     ,fsntoac      ,fsnirtoa     ,      &
                            fsnrtoac ,fsnrtoaq   ,fsns         ,fsnsc        ,      &
                            fsdsc    ,fsusc      ,fsdtc        ,fsutc        ,      & 
                            fsds     ,fsus       ,fsdt         ,fsut         ,      & 
                            sols     ,soll       ,solsd        ,solld        ,      &
                            fns      ,fcns       ,Nday         ,Nnite        ,      &
                            IdxDay   ,IdxNite    ,su           ,sd           ,      &
                            E_cld_tau, E_cld_tau_w, E_cld_tau_w_g, E_cld_tau_w_f,   &
                            old_convert, flg_4DDA_sw)
    use cmparray_mod,        only: CmpDayNite, ExpDayNite
    use mcica_subcol_gen_sw, only: mcica_subcol_sw
    use grist_constants,     only: cpair => cp
    use rrtmg_state,         only: rrtmg_state_t
 
    ! io
    integer, intent(in) :: ncol                      ! number of atmospheric columns
    integer, intent(in) :: rrtmg_levs                ! number of levels rad is applied

    type(rrtmg_state_t), intent(in) :: r_state

    integer, intent(in) :: Nday                      ! Number of daylight columns
    integer, intent(in) :: Nnite                     ! Number of night columns
    integer, intent(in), dimension(ncol) :: IdxDay   ! Indicies of daylight coumns
    integer, intent(in), dimension(ncol) :: IdxNite  ! Indicies of night coumns

    real(r8), intent(in) :: E_pmid(nlev, ncol)       ! Level pressure (Pascals)
    real(r8), intent(in) :: E_cld(nlev, ncol)        ! Fractional cloud cover

    real(r8), intent(in) :: E_aer_tau    (nbndsw, 0:nlev, ncol)      ! aerosol optical depth
    real(r8), intent(in) :: E_aer_tau_w  (nbndsw, 0:nlev, ncol)      ! aerosol OD * ssa
    real(r8), intent(in) :: E_aer_tau_w_g(nbndsw, 0:nlev, ncol)      ! aerosol OD * ssa * asm
    real(r8), intent(in) :: E_aer_tau_w_f(nbndsw, 0:nlev, ncol)      ! aerosol OD * ssa * fwd

    real(r8), intent(in) :: eccf               ! Eccentricity factor (1./earth-sun dist^2)
    real(r8), intent(in) :: E_coszrs(ncol)     ! Cosine solar zenith angle
    real(r8), intent(in) :: E_asdir(ncol)      ! 0.2-0.7 micro-meter srfc alb: direct rad
    real(r8), intent(in) :: E_aldir(ncol)      ! 0.7-5.0 micro-meter srfc alb: direct rad
    real(r8), intent(in) :: E_asdif(ncol)      ! 0.2-0.7 micro-meter srfc alb: diffuse rad
    real(r8), intent(in) :: E_aldif(ncol)      ! 0.7-5.0 micro-meter srfc alb: diffuse rad
    real(r8), intent(in) :: sfac(nbndsw)       ! factor to account for solar variability in each band 

    real(r8), intent(in) :: E_cld_tau    (nbndsw, nlev, ncol)      ! cloud optical depth
    real(r8), intent(in) :: E_cld_tau_w  (nbndsw, nlev, ncol)      ! cloud optical 
    real(r8), intent(in) :: E_cld_tau_w_g(nbndsw, nlev, ncol)      ! cloud optical 
    real(r8), intent(in) :: E_cld_tau_w_f(nbndsw, nlev, ncol)      ! cloud optical 
    logical,  intent(in) :: old_convert
    logical, intent(in)  :: FLG_4DDA_SW   ! linhan

    real(r8), intent(out) :: solin(ncol)       ! Incident solar flux
    real(r8), intent(out) :: qrs (nlev, ncol)  ! Solar heating rate
    real(r8), intent(out) :: qrsc(nlev, ncol)  ! Clearsky solar heating rate
    real(r8), intent(out) :: fsns(ncol)        ! Surface absorbed solar flux
    real(r8), intent(out) :: fsnt(ncol)        ! Total column absorbed solar flux
    real(r8), intent(out) :: fsntoa(ncol)      ! Net solar flux at TOA
    real(r8), intent(out) :: fsutoa(ncol)      ! Upward solar flux at TOA
    real(r8), intent(out) :: fsds(ncol)        ! Flux shortwave downwelling surface
    real(r8), intent(out) :: fsus(ncol)        ! Flux shortwave upwelling surface
    real(r8), intent(out) :: fsdt(ncol)        ! Flux shortwave downwelling at TOA
    real(r8), intent(out) :: fsut(ncol)        ! Flux shortwave upwelling at TOA

    real(r8), intent(out) :: fsdsc(ncol)       ! Clear sky surface downwelling solar flux
    real(r8), intent(out) :: fsusc(ncol)       ! Clear sky surface upwelling solar flux
    real(r8), intent(out) :: fsdtc(ncol)       ! Clear sky downwelling solar flux at TOA
    real(r8), intent(out) :: fsutc(ncol)       ! Clear sky upwelling solar flux at TOA 
 
    real(r8), intent(out) :: fsnsc(ncol)       ! Clear sky surface absorbed solar flux
    real(r8), intent(out) :: fsntc(ncol)       ! Clear sky total column absorbed solar flx
    real(r8), intent(out) :: fsntoac(ncol)     ! Clear sky net solar flx at TOA
    real(r8), intent(out) :: sols(ncol)        ! Direct solar rad on surface (< 0.7)
    real(r8), intent(out) :: soll(ncol)        ! Direct solar rad on surface (>= 0.7)
    real(r8), intent(out) :: solsd(ncol)       ! Diffuse solar rad on surface (< 0.7)
    real(r8), intent(out) :: solld(ncol)       ! Diffuse solar rad on surface (>= 0.7)
    real(r8), intent(out) :: fsnirtoa(ncol)    ! Near-IR flux absorbed at toa
    real(r8), intent(out) :: fsnrtoac(ncol)    ! Clear sky near-IR flux absorbed at toa
    real(r8), intent(out) :: fsnrtoaq(ncol)    ! Net near-IR flux at toa >= 0.7 microns

    real(r8), intent(out) :: fns(nlevp, ncol)  ! net flux at interfaces
    real(r8), intent(out) :: fcns(nlevp, ncol) ! net clear-sky flux at interfaces

    real(r8), intent(out) :: su(nbndsw,nlevp,ncol) ! shortwave spectral flux up
    real(r8), intent(out) :: sd(nbndsw,nlevp,ncol) ! shortwave spectral flux down


    ! local
    ! Minimum cloud amount (as a fraction of the grid-box area) to 
    ! distinguish from clear sky
    real(r8), parameter :: cldmin = 1.0e-80_r8

    ! Decimal precision of cloud amount (0 -> preserve full resolution;
    ! 10^-n -> preserve n digits of cloud amount)
    real(r8), parameter :: cldeps = 0.0_r8

    real(r8) :: pmid(nlev, ncol)               ! Level pressure (Pascals)

    real(r8) :: cld(rrtmg_levs-1, ncol)        ! Fractional cloud cover
    real(r8) :: cicewp(rrtmg_levs-1, ncol)     ! in-cloud cloud ice water path
    real(r8) :: cliqwp(rrtmg_levs-1, ncol)     ! in-cloud cloud liquid water path
    real(r8) :: rel(rrtmg_levs-1, ncol)        ! Liquid effective drop size (microns)
    real(r8) :: rei(rrtmg_levs-1, ncol)        ! Ice effective drop size (microns)

    real(r8) :: coszrs(ncol)                   ! Cosine solar zenith angle
    real(r8) :: asdir(ncol)                    ! 0.2-0.7 micro-meter srfc alb: direct rad
    real(r8) :: aldir(ncol)                    ! 0.7-5.0 micro-meter srfc alb: direct rad
    real(r8) :: asdif(ncol)                    ! 0.2-0.7 micro-meter srfc alb: diffuse rad
    real(r8) :: aldif(ncol)                    ! 0.7-5.0 micro-meter srfc alb: diffuse rad

    real(r8) :: h2ovmr(rrtmg_levs, ncol)       ! h2o volume mixing ratio
    real(r8) :: o3vmr(rrtmg_levs, ncol)        ! o3 volume mixing ratio
    real(r8) :: co2vmr(rrtmg_levs, ncol)       ! co2 volume mixing ratio 
    real(r8) :: ch4vmr(rrtmg_levs, ncol)       ! ch4 volume mixing ratio 
    real(r8) :: o2vmr(rrtmg_levs, ncol)        ! o2  volume mixing ratio 
    real(r8) :: n2ovmr(rrtmg_levs, ncol)       ! n2o volume mixing ratio 

    real(r8) :: tsfc(ncol)                     ! surface temperature

    integer :: inflgsw                         ! flag for cloud parameterization method
    integer :: iceflgsw                        ! flag for ice cloud parameterization method
    integer :: liqflgsw                        ! flag for liquid cloud parameterization method
    integer :: icld                            ! Flag for cloud overlap method
                                               ! 0=clear, 1=random, 2=maximum/random, 3=maximum
    integer :: dyofyr                          ! Set to day of year for Earth/Sun distance calculation in
                                               ! rrtmg_sw, or pass in adjustment directly into adjes
    real(r8) :: solvar(nbndsw)                 ! solar irradiance variability in each band

    integer, parameter :: nsubcsw = ngptsw     ! rrtmg_sw g-point (quadrature point) dimension
    integer :: permuteseed                     ! permute seed for sub-column generator

    real(r8) :: diagnostic_od(nlev, ncol)      ! cloud optical depth - diagnostic temp variable

    real(r8) :: tauc_sw(nbndsw, rrtmg_levs-1, ncol)         ! cloud optical depth
    real(r8) :: ssac_sw(nbndsw, rrtmg_levs-1, ncol)         ! cloud single scat. albedo
    real(r8) :: asmc_sw(nbndsw, rrtmg_levs-1, ncol)         ! cloud asymmetry parameter
    real(r8) :: fsfc_sw(nbndsw, rrtmg_levs-1, ncol)         ! cloud forward scattering fraction

    real(r8) :: tau_aer_sw(nbndsw, rrtmg_levs-1, ncol)      ! aer optical depth
    real(r8) :: ssa_aer_sw(nbndsw, rrtmg_levs-1, ncol)      ! aer single scat. albedo
    real(r8) :: asm_aer_sw(nbndsw, rrtmg_levs-1, ncol)      ! aer asymmetry parameter

    real(r8) :: cld_stosw(nsubcsw, rrtmg_levs-1, ncol)      ! stochastic cloud fraction
    real(r8) :: rei_stosw(rrtmg_levs-1, ncol)               ! stochastic ice particle size 
    real(r8) :: rel_stosw(rrtmg_levs-1, ncol)               ! stochastic liquid particle size
    real(r8) :: cicewp_stosw(nsubcsw, rrtmg_levs-1, ncol)   ! stochastic cloud ice water path
    real(r8) :: cliqwp_stosw(nsubcsw, rrtmg_levs-1, ncol)   ! stochastic cloud liquid wter path
    real(r8) :: tauc_stosw(nsubcsw, rrtmg_levs-1, ncol)     ! stochastic cloud optical depth (optional)
    real(r8) :: ssac_stosw(nsubcsw, rrtmg_levs-1, ncol)     ! stochastic cloud single scat. albedo (optional)
    real(r8) :: asmc_stosw(nsubcsw, rrtmg_levs-1, ncol)     ! stochastic cloud asymmetry parameter (optional)
    real(r8) :: fsfc_stosw(nsubcsw, rrtmg_levs-1, ncol)     ! stochastic cloud forward scattering fraction (optional)

    real(r8), parameter :: dps = 1._r8/86400._r8            ! Inverse of seconds per day
 
    real(r8) :: swuflx(rrtmg_levs+1, ncol)                  ! Total sky shortwave upward flux (W/m2)
    real(r8) :: swdflx(rrtmg_levs+1, ncol)                  ! Total sky shortwave downward flux (W/m2)
    real(r8) :: swhr(rrtmg_levs, ncol)                      ! Total sky shortwave radiative heating rate (K/d)
    real(r8) :: swuflxc(rrtmg_levs+1, ncol)                 ! Clear sky shortwave upward flux (W/m2)
    real(r8) :: swdflxc(rrtmg_levs+1, ncol)                 ! Clear sky shortwave downward flux (W/m2)
    real(r8) :: swhrc(rrtmg_levs, ncol)                     ! Clear sky shortwave radiative heating rate (K/d)
    real(r8) :: swuflxs(nbndsw, rrtmg_levs+1, ncol)         ! Shortwave spectral flux up
    real(r8) :: swdflxs(nbndsw, rrtmg_levs+1, ncol)         ! Shortwave spectral flux down

    real(r8) :: dirdnuv(rrtmg_levs+1, ncol)                 ! Direct downward shortwave flux, UV/vis
    real(r8) :: difdnuv(rrtmg_levs+1, ncol)                 ! Diffuse downward shortwave flux, UV/vis
    real(r8) :: dirdnir(rrtmg_levs+1, ncol)                 ! Direct downward shortwave flux, near-IR
    real(r8) :: difdnir(rrtmg_levs+1, ncol)                 ! Diffuse downward shortwave flux, near-IR

    ! Added for net near-IR diagnostic
    real(r8) :: ninflx(rrtmg_levs+1, ncol)                  ! Net shortwave flux, near-IR
    real(r8) :: ninflxc(rrtmg_levs+1, ncol)                 ! Net clear sky shortwave flux, near-IR

    ! Other
    integer :: i, k, ns       ! indices

    ! Cloud radiative property arrays
    real(r8) :: tauxcl(0:nlev, ncol) ! water cloud extinction optical depth
    real(r8) :: tauxci(0:nlev, ncol) ! ice cloud extinction optical depth
    real(r8) :: wcl(0:nlev, ncol)    ! liquid cloud single scattering albedo
    real(r8) :: gcl(0:nlev, ncol)    ! liquid cloud asymmetry parameter
    real(r8) :: fcl(0:nlev, ncol)    ! liquid cloud forward scattered fraction
    real(r8) :: wci(0:nlev, ncol)    ! ice cloud single scattering albedo
    real(r8) :: gci(0:nlev, ncol)    ! ice cloud asymmetry parameter
    real(r8) :: fci(0:nlev, ncol)    ! ice cloud forward scattered fraction

    ! Aerosol radiative property arrays
    real(r8) :: tauxar(0:nlev, ncol) ! aerosol extinction optical depth
    real(r8) :: wa(0:nlev, ncol)     ! aerosol single scattering albedo
    real(r8) :: ga(0:nlev, ncol)     ! aerosol assymetry parameter
    real(r8) :: fa(0:nlev, ncol)     ! aerosol forward scattered fraction

    ! CRM
    real(r8) :: fus(nlevp, ncol)     ! Upward flux (added for CRM)
    real(r8) :: fds(nlevp, ncol)     ! Downward flux (added for CRM)
    real(r8) :: fusc(nlevp, ncol)    ! Upward clear-sky flux (added for CRM)
    real(r8) :: fdsc(nlevp, ncol)    ! Downward clear-sky flux (added for CRM)

    integer :: kk

    real(r8) :: pmidmb(rrtmg_levs, ncol)   ! Level pressure (hPa)
    real(r8) :: pintmb(rrtmg_levs+1, ncol) ! Model interface pressure (hPa)
    real(r8) :: tlay(rrtmg_levs, ncol)     ! mid point temperature
    real(r8) :: tlev(rrtmg_levs+1, ncol)   ! interface temperature

    ! Initialize output fields:
    fsds(1:ncol)     = 0.0_r8
    fsdsc(1:ncol)    = 0.0_r8
    fsus(1:ncol)     = 0.0_r8
    fsusc(1:ncol)    = 0.0_r8
    fsdt(1:ncol)     = 0.0_r8
    fsdtc(1:ncol)    = 0.0_r8
    fsut(1:ncol)     = 0.0_r8
    fsutc(1:ncol)    = 0.0_r8

    fsnirtoa(1:ncol) = 0.0_r8
    fsnrtoac(1:ncol) = 0.0_r8
    fsnrtoaq(1:ncol) = 0.0_r8

    fsns(1:ncol)     = 0.0_r8
    fsnsc(1:ncol)    = 0.0_r8

    fsnt(1:ncol)     = 0.0_r8
    fsntc(1:ncol)    = 0.0_r8
    fsntoa(1:ncol)   = 0.0_r8
    fsutoa(1:ncol)   = 0.0_r8
    fsntoac(1:ncol)  = 0.0_r8

    solin(1:ncol)    = 0.0_r8

    sols(1:ncol)     = 0.0_r8
    soll(1:ncol)     = 0.0_r8
    solsd(1:ncol)    = 0.0_r8
    solld(1:ncol)    = 0.0_r8

    qrs (1:nlev,1:ncol)  = 0.0_r8
    qrsc(1:nlev,1:ncol)  = 0.0_r8
    fns(1:nlevp,1:ncol)  = 0.0_r8
    fcns(1:nlevp,1:ncol) = 0.0_r8
    !--------------LiXH has not completed scm_crm_mode------------------>
    !if (single_column.and.scm_crm_mode) then 
    !   fus(1:nlevp,1:ncol)  = 0.0_r8
    !   fds(1:nlevp,1:ncol)  = 0.0_r8
    !   fusc(1:nlevp,1:ncol) = 0.0_r8
    !   fdsc(1:nlevp,1:ncol) = 0.0_r8
    !endif
    !<-------------LiXH has not completed scm_crm_mode-------------------

    su(:,:,1:ncol) = 0.0_r8
    sd(:,:,1:ncol) = 0.0_r8

    ! If night everywhere, return:
    if ( Nday == 0 ) then
      return
    endif

    ! Rearrange input arrays
    call CmpDayNite(E_pmid(nlevp-rrtmg_levs+1:nlev,1:ncol), pmid(1:rrtmg_levs-1,1:ncol), &
         Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs-1, 1, ncol)
    call CmpDayNite(E_cld(nlevp-rrtmg_levs+1:nlev,1:ncol),  cld(1:rrtmg_levs-1,1:ncol),  &
         Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs-1, 1, ncol)

    call CmpDayNite(r_state%pintmb, pintmb, Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs+1, 1, ncol)
    call CmpDayNite(r_state%pmidmb, pmidmb, Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)
    call CmpDayNite(r_state%h2ovmr, h2ovmr, Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)
    call CmpDayNite(r_state%o3vmr,  o3vmr,  Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)
    call CmpDayNite(r_state%co2vmr, co2vmr, Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)

    call CmpDayNite(E_coszrs, coszrs,    Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call CmpDayNite(E_asdir,  asdir,     Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call CmpDayNite(E_aldir,  aldir,     Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call CmpDayNite(E_asdif,  asdif,     Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call CmpDayNite(E_aldif,  aldif,     Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call CmpDayNite(r_state%tlay,   tlay,   Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)
    call CmpDayNite(r_state%tlev,   tlev,   Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs+1, 1, ncol)
    call CmpDayNite(r_state%ch4vmr, ch4vmr, Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)
    call CmpDayNite(r_state%o2vmr,  o2vmr,  Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)
    call CmpDayNite(r_state%n2ovmr, n2ovmr, Nday, IdxDay, Nnite, IdxNite, 1, rrtmg_levs, 1, ncol)

    ! These fields are no longer input by CAM.
    cicewp = 0.0_r8
    cliqwp = 0.0_r8
    rel = 0.0_r8
    rei = 0.0_r8

    ! Aerosol daylight map
    ! Also convert to optical properties of rrtmg interface, even though
    ! these quantities are later multiplied back together inside rrtmg 
    ! Why does rrtmg use the factored quantities?
    ! There are several different ways this factoring could be done.
    ! Other ways might allow for better optimization
    do ns = 1, nbndsw
       do k  = 1, rrtmg_levs-1
          kk=(nlevp-rrtmg_levs) + k
          do i  = 1, Nday
             if(E_aer_tau_w(ns, kk, IdxDay(i)) > 1.e-80_r8) then
                asm_aer_sw(ns,k,i) = E_aer_tau_w_g(ns,kk,IdxDay(i))/E_aer_tau_w(ns,kk,IdxDay(i))
             else
                asm_aer_sw(ns,k,i) = 0._r8
             endif
             if(E_aer_tau(ns,kk,IdxDay(i)) > 0._r8) then
                ssa_aer_sw(ns,k,i) = E_aer_tau_w(ns,kk,IdxDay(i))/E_aer_tau(ns,kk,IdxDay(i))
                tau_aer_sw(ns,k,i) = E_aer_tau(ns,kk,IdxDay(i))
             else
                ssa_aer_sw(ns,k,i) = 1._r8
                tau_aer_sw(ns,k,i) = 0._r8
             endif
          enddo
       enddo
    enddo

    !--------------LiXH has not completed scm_crm_mode------------------>
    !if (scm_crm_mode) then
    !   ! overwrite albedos for CRM
    !   if(have_asdir) asdir = asdirobs(1)
    !   if(have_asdif) asdif = asdifobs(1)
    !   if(have_aldir) aldir = aldirobs(1)
    !   if(have_aldif) aldif = aldifobs(1)
    !endif
    !<-------------LiXH has not completed scm_crm_mode-------------------

    ! Define solar incident radiation
    do i = 1, Nday
       solin(i)  = sum(sfac(:)*solar_band_irrad(:)) * eccf * coszrs(i)
    end do

    ! Calculate cloud optical properties here if using CAM method, or if using one of the
    ! methods in RRTMG_SW, then pass in cloud physical properties and zero out cloud optical 
    ! properties here

    ! Zero optional cloud optical property input arrays tauc_sw, ssac_sw, asmc_sw, 
    ! if inputting cloud physical properties to RRTMG_SW
    !tauc_sw(:,:,:) = 0.0_r8
    !ssac_sw(:,:,:) = 1.0_r8
    !asmc_sw(:,:,:) = 0.0_r8
    !fsfc_sw(:,:,:) = 0.0_r8
    !
    ! Or, calculate and pass in CAM cloud shortwave optical properties to RRTMG_SW
    !if (present(old_convert)) print *, 'old_convert',old_convert
    !if (present(ancientmethod)) print *, 'ancientmethod',ancientmethod

    !----------LiXH modified, ifort can not use optional--------------->
    !if (present(old_convert))then
       if (old_convert)then ! convert without limits
          do i = 1, Nday
          do k = 1, rrtmg_levs-1
          kk=(nlevp-rrtmg_levs) + k
          do ns = 1, nbndsw
            if (E_cld_tau_w(ns,kk,IdxDay(i)) > 0._r8) then
               fsfc_sw(ns,k,i)=E_cld_tau_w_f(ns,kk,IdxDay(i))/E_cld_tau_w(ns,kk,IdxDay(i))
               asmc_sw(ns,k,i)=E_cld_tau_w_g(ns,kk,IdxDay(i))/E_cld_tau_w(ns,kk,IdxDay(i))
            else
               fsfc_sw(ns,k,i) = 0._r8
               asmc_sw(ns,k,i) = 0._r8
            endif
    
            tauc_sw(ns,k,i)=E_cld_tau(ns,kk,IdxDay(i))
            if (tauc_sw(ns,k,i) > 0._r8) then
               ssac_sw(ns,k,i)=E_cld_tau_w(ns,kk,IdxDay(i))/tauc_sw(ns,k,i)
            else
               tauc_sw(ns,k,i) = 0._r8
               fsfc_sw(ns,k,i) = 0._r8
               asmc_sw(ns,k,i) = 0._r8
               ssac_sw(ns,k,i) = 1._r8
            endif
          enddo
          enddo
          enddo
       else
          ! eventually, when we are done with archaic versions, This set of code will become the default.
          do i = 1, Nday
          do k = 1, rrtmg_levs-1
          kk=(nlevp-rrtmg_levs) + k
          do ns = 1, nbndsw
            if (E_cld_tau_w(ns,kk,IdxDay(i)) > 0._r8) then
               fsfc_sw(ns,k,i)=E_cld_tau_w_f(ns,kk,IdxDay(i))/max(E_cld_tau_w(ns,kk,IdxDay(i)), 1.e-80_r8)
               asmc_sw(ns,k,i)=E_cld_tau_w_g(ns,kk,IdxDay(i))/max(E_cld_tau_w(ns,kk,IdxDay(i)), 1.e-80_r8)
            else
               fsfc_sw(ns,k,i) = 0._r8
               asmc_sw(ns,k,i) = 0._r8
            endif
    
            tauc_sw(ns,k,i)=E_cld_tau(ns,kk,IdxDay(i))
            if (tauc_sw(ns,k,i) > 0._r8) then
               ssac_sw(ns,k,i)=max(E_cld_tau_w(ns,kk,IdxDay(i)),1.e-80_r8)/max(tauc_sw(ns,k,i),1.e-80_r8)
            else
               tauc_sw(ns,k,i) = 0._r8
               fsfc_sw(ns,k,i) = 0._r8
               asmc_sw(ns,k,i) = 0._r8
               ssac_sw(ns,k,i) = 1._r8
            endif
          enddo
          enddo
          enddo
       endif
    !----------LiXH modified, ifort can not use optional--------------->
    !else
    !   do i = 1, Nday
    !   do k = 1, rrtmg_levs-1
    !   kk=(nlevp-rrtmg_levs) + k
    !   do ns = 1, nbndsw
    !     if (E_cld_tau_w(ns,kk,IdxDay(i)) > 0._r8) then
    !        fsfc_sw(ns,k,i)=E_cld_tau_w_f(ns,kk,IdxDay(i))/max(E_cld_tau_w(ns,kk,IdxDay(i)), 1.e-80_r8)
    !        asmc_sw(ns,k,i)=E_cld_tau_w_g(ns,kk,IdxDay(i))/max(E_cld_tau_w(ns,kk,IdxDay(i)), 1.e-80_r8)
    !     else
    !        fsfc_sw(ns,k,i) = 0._r8
    !        asmc_sw(ns,k,i) = 0._r8
    !     endif

    !     tauc_sw(ns,k,i)=E_cld_tau(ns,kk,IdxDay(i))
    !     if (tauc_sw(ns,k,i) > 0._r8) then
    !        ssac_sw(ns,k,i)=max(E_cld_tau_w(ns,kk,IdxDay(i)),1.e-80_r8)/max(tauc_sw(ns,k,i),1.e-80_r8)
    !     else
    !        tauc_sw(ns,k,i) = 0._r8
    !        fsfc_sw(ns,k,i) = 0._r8
    !        asmc_sw(ns,k,i) = 0._r8
    !        ssac_sw(ns,k,i) = 1._r8
    !     endif
    !   enddo
    !   enddo
    !   enddo
    !endif
    !<---------LiXH modified, ifort can not use optional----------------

    ! Call mcica sub-column generator for RRTMG_SW

    ! Call sub-column generator for McICA in radiation

    ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
    icld = 2
    ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
    ! effective randomization)
    permuteseed = 1


    call mcica_subcol_sw(Nday, rrtmg_levs-1, icld, permuteseed, pmid, &
       cld, cicewp, cliqwp, rei, rel, tauc_sw, ssac_sw, asmc_sw, fsfc_sw, &
       cld_stosw, cicewp_stosw, cliqwp_stosw, rei_stosw, rel_stosw, &
       tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw)

    ! Call RRTMG_SW for all layers for daylight columns

    ! Select parameterization of cloud ice and liquid optical depths
    ! Use CAM shortwave cloud optical properties directly
    inflgsw = 0 
    iceflgsw = 0
    liqflgsw = 0
    ! Use E&C param for ice to mimic CAM3 for now
    !   inflgsw = 2 
    !   iceflgsw = 1
    !   liqflgsw = 1
    ! Use merged Fu and E&C params for ice 
    !   inflgsw = 2 
    !   iceflgsw = 3
    !   liqflgsw = 1

    ! Set day of year for Earth/Sun distance calculation in rrtmg_sw, or
    ! set to zero and pass E/S adjustment (eccf) directly into array adjes
    dyofyr = 0

    tsfc(:ncol) = tlev(rrtmg_levs+1,:ncol)

    solvar(1:nbndsw) = sfac(1:nbndsw)

    if (FLG_4DDA_SW) then
    call rrtmg_4DDA_sw(Nday, rrtmg_levs, icld,         &
                  pmidmb, pintmb, tlay, tlev, tsfc, &
                  h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, &
                  asdir, asdif, aldir, aldif, &
                  coszrs, eccf, dyofyr, solvar, &
                  inflgsw, iceflgsw, liqflgsw, &
                  cld_stosw, tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw, &
                  cicewp_stosw, cliqwp_stosw, rei, rel, &
                  tau_aer_sw, ssa_aer_sw, asm_aer_sw, &
                  swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
                  dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc, swuflxs, swdflxs)
    else
    call rrtmg_sw(Nday, rrtmg_levs, icld,         &
                  pmidmb, pintmb, tlay, tlev, tsfc, &
                  h2ovmr, o3vmr, co2vmr, ch4vmr, o2vmr, n2ovmr, &
                  asdir, asdif, aldir, aldif, &
                  coszrs, eccf, dyofyr, solvar, &
                  inflgsw, iceflgsw, liqflgsw, &
                  cld_stosw, tauc_stosw, ssac_stosw, asmc_stosw, fsfc_stosw, &
                  cicewp_stosw, cliqwp_stosw, rei, rel, &
                  tau_aer_sw, ssa_aer_sw, asm_aer_sw, &
                  swuflx, swdflx, swhr, swuflxc, swdflxc, swhrc, &
                  dirdnuv, dirdnir, difdnuv, difdnir, ninflx, ninflxc, swuflxs, swdflxs)
    endif

    ! Flux units are in W/m2 on output from rrtmg_sw and contain output for
    ! extra layer above model top with vertical indexing from bottom to top.
    !
    ! Heating units are in J/kg/s on output from rrtmg_sw and contain output 
    ! for extra layer above model top with vertical indexing from bottom to top.  
    !
    ! Reverse vertical indexing to go from top to bottom for CAM output.

    ! Set the net absorted shortwave flux at TOA (top of extra layer)
    fsntoa(1:Nday) = swdflx(rrtmg_levs+1,1:Nday) - swuflx(rrtmg_levs+1,1:Nday)
    fsutoa(1:Nday) = swuflx(rrtmg_levs+1,1:Nday)
    fsntoac(1:Nday) = swdflxc(rrtmg_levs+1,1:Nday) - swuflxc(rrtmg_levs+1,1:Nday)

    ! Set net near-IR flux at top of the model
    fsnirtoa(1:Nday) = ninflx(rrtmg_levs,1:Nday)
    fsnrtoaq(1:Nday) = ninflx(rrtmg_levs,1:Nday)
    fsnrtoac(1:Nday) = ninflxc(rrtmg_levs,1:Nday)

    ! Set the net absorbed shortwave flux at the model top level
    fsnt(1:Nday) = swdflx(rrtmg_levs,1:Nday) - swuflx(rrtmg_levs,1:Nday)
    fsntc(1:Nday) = swdflxc(rrtmg_levs,1:Nday) - swuflxc(rrtmg_levs,1:Nday)

    ! Set the downwelling flux at the model top level
    fsdt(1:Nday) = swdflx(rrtmg_levs,1:Nday)
    fsdtc(1:Nday) = swdflxc(rrtmg_levs,1:Nday)
    ! Set the upwelling flux at the model top level
    fsut(1:Nday) = swuflx(rrtmg_levs,1:Nday)
    fsutc(1:Nday) = swuflxc(rrtmg_levs,1:Nday)

    ! Set the downwelling flux at the surface 
    fsds(1:Nday) = swdflx(1,1:Nday)
    fsdsc(1:Nday) = swdflxc(1,1:Nday)
    ! Set the upwelling flux at the surface 
    fsus(1:Nday) = swuflx(1,1:Nday)
    fsusc(1:Nday) = swuflxc(1,1:Nday)

    ! Set the net shortwave flux at the surface
    fsns(1:Nday) = swdflx(1,1:Nday) - swuflx(1,1:Nday)
    fsnsc(1:Nday) = swdflxc(1,1:Nday) - swuflxc(1,1:Nday)

    ! Set the UV/vis and near-IR direct and dirruse downward shortwave flux at surface
    sols(1:Nday) = dirdnuv(1,1:Nday)
    soll(1:Nday) = dirdnir(1,1:Nday)
    solsd(1:Nday) = difdnuv(1,1:Nday)
    solld(1:Nday) = difdnir(1,1:Nday)


    ! Set the net, up and down fluxes at model interfaces
    fns (nlevp-rrtmg_levs+1:nlevp,1:Nday) =  swdflx(rrtmg_levs:1:-1,1:Nday) -  swuflx(rrtmg_levs:1:-1,1:Nday)
    fcns(nlevp-rrtmg_levs+1:nlevp,1:Nday) = swdflxc(rrtmg_levs:1:-1,1:Nday) - swuflxc(rrtmg_levs:1:-1,1:Nday)
    fus (nlevp-rrtmg_levs+1:nlevp,1:Nday) =  swuflx(rrtmg_levs:1:-1,1:Nday)
    fusc(nlevp-rrtmg_levs+1:nlevp,1:Nday) = swuflxc(rrtmg_levs:1:-1,1:Nday)
    fds (nlevp-rrtmg_levs+1:nlevp,1:Nday) =  swdflx(rrtmg_levs:1:-1,1:Nday)
    fdsc(nlevp-rrtmg_levs+1:nlevp,1:Nday) = swdflxc(rrtmg_levs:1:-1,1:Nday)

    ! Set solar heating, reverse layering
    ! Pass shortwave heating to CAM arrays and convert from K/d to J/kg/s
    qrs (nlevp-rrtmg_levs+1:nlev,1:Nday) = swhr (rrtmg_levs-1:1:-1,1:Nday)*cpair*dps
    qrsc(nlevp-rrtmg_levs+1:nlev,1:Nday) = swhrc(rrtmg_levs-1:1:-1,1:Nday)*cpair*dps
    !linhan test
   ! if(mpi_rank()==0)then
   ! print*,'linhan test sw' !linhan
   ! print*,rrtmg_levs+1,swuflx(rrtmg_levs+1,1:Nday),swdflx(rrtmg_levs+1,1:Nday) !linhan
   !     do k=rrtmg_levs,1,-1   !linhan
   !     print*,k,
   !     print*,swhr(k,1:Nday),swuflx(k,1:Nday),swdflx(k,1:Nday) !linhan
   !     enddo  !linhan
   ! print*,'linhan test swend' !linhan
   ! end if

    ! Set spectral fluxes, reverse layering
    ! order=(/3,1,2/) maps the first index of swuflxs to the third index of su.
    su(:,nlevp-rrtmg_levs+1:nlevp,1:Nday) = swuflxs(:,rrtmg_levs:1:-1,1:Nday)

    sd(:,nlevp-rrtmg_levs+1:nlevp,1:Nday) = swdflxs(:,rrtmg_levs:1:-1,1:Nday)

    ! Rearrange output arrays.
    !
    ! intent(out)

    call ExpDayNite(solin,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(qrs,	Nday, IdxDay, Nnite, IdxNite, 1, nlev, 1, ncol)
    call ExpDayNite(qrsc,	Nday, IdxDay, Nnite, IdxNite, 1, nlev, 1, ncol)
    call ExpDayNite(fns,	Nday, IdxDay, Nnite, IdxNite, 1, nlevp, 1, ncol)
    call ExpDayNite(fcns,	Nday, IdxDay, Nnite, IdxNite, 1, nlevp, 1, ncol)
    call ExpDayNite(fsns,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsnt,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsntoa,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsutoa,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsds,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsus,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsdt,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsut,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsnsc,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsdsc,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsusc,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsdtc,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsutc,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsntc,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsntoac,Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(sols,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(soll,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(solsd,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(solld,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsnirtoa,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsnrtoac,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)
    call ExpDayNite(fsnrtoaq,	Nday, IdxDay, Nnite, IdxNite, 1, ncol)

    call ExpDayNite(su,	Nday, IdxDay, Nnite, IdxNite, 1, nbndsw, 1, nlevp, 1, ncol)
    call ExpDayNite(sd,	Nday, IdxDay, Nnite, IdxNite, 1, nbndsw, 1, nlevp, 1, ncol)

    !--------------LiXH has not completed scm_crm_mode------------------>
    !  these outfld calls don't work for spmd only outfield in scm mode (nonspmd)
    !if (single_column .and. scm_crm_mode) then 
    !   ! Following outputs added for CRM
    !   call ExpDayNite(fus,Nday, IdxDay, Nnite, IdxNite, 1, nlevp, 1, ncol)
    !   call ExpDayNite(fds,Nday, IdxDay, Nnite, IdxNite, 1, nlevp, 1, ncol)
    !   call ExpDayNite(fusc,Nday, IdxDay, Nnite, IdxNite, 1, nlevp, 1, ncol)
    !   call ExpDayNite(fdsc,Nday, IdxDay, Nnite, IdxNite, 1, nlevp, 1, ncol)
    !   ! output:
    !   ! fus, fds, fusc, fdsc
    !endif
    !<-------------LiXH has not completed scm_crm_mode-------------------

    end subroutine rad_rrtmg_sw


! Purpose : Initialize various constants for radiation scheme.
    subroutine radsw_init()
    use radconstants,  only: get_solar_band_fraction_irrad, get_ref_solar_band_irrad

    ! get the reference fractional solar irradiance in each band
    call get_solar_band_fraction_irrad(fractional_solar_irradiance)
    call get_ref_solar_band_irrad( solar_band_irrad )

    ! Initialize rrtmg_sw
    call rrtmg_sw_ini
 
    end subroutine radsw_init


 end module radsw
