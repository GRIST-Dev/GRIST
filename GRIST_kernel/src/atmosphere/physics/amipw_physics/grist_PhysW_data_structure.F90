
!----------------------------------------------------------------------------
! Created on 2019
! Author: Yi Zhang
! Version 1.0
! Description: This module contains data specific to WRF-based physics, it 
! also follows their name convention currently
! Revision history:
!--------------------------------------------------------------------------

 module grist_PhysW_data_structure

   use grist_constants,          only: r8, i4, zero
   use grist_nml_module,         only: run_type
   use grist_PhysW_nml_module,   only: PhysW_mp_scheme,   &
                                       PhysW_bl_scheme
   use grist_data_types,         only: scalar_1d_field
 
   implicit none

   private

   public :: pstate_wrf,&
             psurf_wrf, &
             ptend_wrf

   public :: grist_PhysW_data_structure_construct, &
             grist_PhysW_data_structure_destruct

    logical,     public  :: restart         = .false.
    logical,     public  :: warm_rain       = .false.
    integer(i4), public  :: num_soil_layers = 1
    integer(i4), public  :: p_qv      ! species index for water vapor
    integer(i4), public  :: p_qc      ! species index for cloud water
    integer(i4), public  :: p_qr      ! species index for rain water
    integer(i4), public  :: p_qi      ! species index for cloud ice
    integer(i4), public  :: p_qs      ! species index for snow
    integer(i4), public  :: p_qg      ! species index for graupel
    integer(i4), public  :: p_nr      ! species index for num of rain water
    integer(i4), public  :: p_ni      ! species index for num of cloud ice
    integer(i4), public  :: p_ns      ! species index for num of snow
    integer(i4), public  :: p_ng      ! species index for num of graupel
    integer(i4), public  :: param_first_scalar = 1

!-------------------------------------------------------------------------------------------

    type physics_state_wrf
!
! atmosphere state
!
      real(r8),  allocatable   ::   u_phy  (:, :, :)  ! U, eq. pstate%u
      real(r8),  allocatable   ::   v_phy  (:, :, :)  ! V, eq. pstate%v
      real(r8),  allocatable   ::   th_phy (:, :, :)  ! PT
      real(r8),  allocatable   ::   t_phy  (:, :, :)  ! T, eq. pstate%temp
      real(r8),  allocatable   ::   t8w    (:, :, :)  ! T
      real(r8),  allocatable   ::   moist  (:, :, :, :) ! tracer, eq. pstate%tracer_mxrt
      real(r8),  allocatable   ::   p_phy  (:, :, :)  ! Pres, eq. pstate%pressure_at_pc_full_level
      real(r8),  allocatable   ::   pi_phy (:, :, :)  ! exner
      real(r8),  allocatable   ::   p8w    (:, :, :)  ! Pres, eq. pstate%pressure_at_pc_face_level
      real(r8),  allocatable   ::   rhom   (:, :, :)  ! moist/full density
      real(r8),  allocatable   ::   rhod   (:, :, :)  ! dry density
      real(r8),  allocatable   ::   dz8w   (:, :, :)  ! dz between full level (w-level is WRF's full)
      real(r8),  allocatable   ::   zzz    (:, :, :)  ! height
      real(r8),  allocatable   ::   www    (:, :, :)  ! vertical speed at w-level, face level of GRIST
      real(r8),  allocatable   ::   omega  (:, :, :)  ! mpressure vertical speed at GRIST's full level, used for CU at HDC
      real(r8),  allocatable   ::   w0avg  (:, :, :)  ! averaged vertical speed, at full level of GRIST
      real(r8),  allocatable   ::   cmfmc  (:, :, :)  ! net updraft mass flux (intended for cam3 cf)
      real(r8),  allocatable   ::   cldfra (:, :, :)  ! total cf
      real(r8),  allocatable   ::   cldcum (:, :, :)  ! cumulus  cf
      real(r8),  allocatable   ::   cldstr (:, :, :)  ! stratiform cf
      real(r8),  allocatable   ::   re_cloud(:, :, :) !
      real(r8),  allocatable   ::   re_ice  (:, :, :) !
      real(r8),  allocatable   ::   re_snow (:, :, :) !
      real(r8),  allocatable   ::   refl_10cm(:, :, :) !
      !real(r8),  allocatable   ::   relhum (:, :, :) !
      real(r8),  allocatable   ::   z8w    (:, :, :)  ! height at w-level, LiXH add
#ifdef CLIMATE_O3
      real(r8),  allocatable   ::   o3     (:, :, :)  ! Ozone, LiXH add   
#endif
!
! geographical
!
      real(r8),  allocatable   ::   xlat   (:, :)
      real(r8),  allocatable   ::   xlong  (:, :)
      real(r8),  allocatable   ::   dxmean (:)
      real(r8),  allocatable   ::   coriolis (:)
      integer(i4),allocatable  ::   kpbl   (:, :)
!
! diagnosed radiation flux from CAMRAD or RRTMG
!
      real(r8), allocatable   ::    swupt(:, :),swuptc(:, :),swdnt(:, :),swdntc(:, :)
      real(r8), allocatable   ::    lwupt(:, :),lwuptc(:, :),lwdnt(:, :),lwdntc(:, :)
      real(r8), allocatable   ::    swupb(:, :),swupbc(:, :),swdnb(:, :),swdnbc(:, :)
      real(r8), allocatable   ::    lwupb(:, :),lwupbc(:, :),lwdnb(:, :),lwdnbc(:, :)
      real(r8), allocatable   ::    swcf (:, :),lwcf  (:, :),olr  (:, :),coszr (:, :)
      real(r8), allocatable   ::    cemiss(:,:,:), taucldc(:,:,:), taucldi(:,:,:)
      real(r8), allocatable   ::    cldtot(:,:), cldlow(:,:), cldmed(:,:), cldhgh(:,:)
!
! diagnosed PBL TKE and Mixing length
!
      real(r8), allocatable   ::    tke (:, :, :) !diag TKE
      real(r8), allocatable   ::    leng(:, :, :) !mixing length
      !Only for UW scheme
#ifdef UW_CAM5
      real(r8), allocatable   ::    kvm(:, :, :)
      real(r8), allocatable   ::    kvh(:, :, :)
      real(r8), allocatable   :: tauresx2d(:, :)
      real(r8), allocatable   :: tauresy2d(:, :)
#endif

    end type physics_state_wrf
!
! surface
!
    type physics_surface_wrf

      real(r8),  allocatable   ::   albedo (:, :)    ! albedo (between 0 and 1)
      real(r8),  allocatable   ::   asdir  (:, :)    ! albedo (between 0 and 1)
      real(r8),  allocatable   ::   asdif  (:, :)    ! albedo (between 0 and 1)
      real(r8),  allocatable   ::   aldir  (:, :)    ! albedo (between 0 and 1)
      real(r8),  allocatable   ::   aldif  (:, :)    ! albedo (between 0 and 1)
      real(r8),  allocatable   ::   br     (:, :)    ! bulk richardson number in surface layer
      !real(r8),  allocatable   ::   capg   (:, :)    ! heat capacity for soil (j/k/m^3)
      !real(r8),  allocatable   ::   chklowq(:, :)    ! out of slab, needed by eta pbl
      real(r8),  allocatable   ::   emiss  (:, :)    ! surface emissivity (between 0 and 1)
      real(r8),  allocatable   ::   flhc   (:, :)    ! exchange coefficient for heat (m/s)
      real(r8),  allocatable   ::   flqc   (:, :)    ! exchange coefficient for moisture (m/s)
      real(r8),  allocatable   ::   glw    (:, :)    ! downward long wave flux at ground surface (w/m^2) 
      real(r8),  allocatable   ::   gsw    (:, :)    ! net downward short wave flux at ground surface (w/m^2) 
      real(r8),  allocatable   ::   swdown (:, :)    ! downward short wave flux at ground surface (w/m^2) 
      real(r8),  allocatable   ::   gz1oz0 (:, :)    ! log(z/z0) where z0 is roughness length
      real(r8),  allocatable   ::   hfx    (:, :)    ! net upward heat flux at the surface (w/m^2)
      real(r8),  allocatable   ::   ht     (:, :)    ! net upward latent heat flux at surface (w/m^2)
      real(r8),  allocatable   ::   lh     (:, :)    ! net upward heat flux at the surface (w/m^2)
      real(r8),  allocatable   ::   mavail (:, :)    ! surface moisture availability (between 0 and 1)
      real(r8),  allocatable   ::   mol    (:, :)    ! q* (similarity theory) (kg/kg)

      real(r8),  allocatable   ::   pblh   (:, :)    ! pbl height (m)
      real(r8),  allocatable   ::   psih   (:, :)    ! similarity stability function for heat
      real(r8),  allocatable   ::   psim   (:, :)    ! similarity stability function for momentum
      real(r8),  allocatable   ::   psfc   (:, :)    ! surface pressure
      real(r8),  allocatable   ::   q2     (:, :)    ! out, diagnostic 2m mixing ratio (kg/kg)
      real(r8),  allocatable   ::   qfx    (:, :)    ! net upward moisture flux at the surface (kg/m^2/s)
      real(r8),  allocatable   ::   qsfc   (:, :)    ! specific humidity at lower boundary (kg/kg)

      real(r8),  allocatable   ::   regime (:, :)    ! flag indicating pbl regime (stable, unstable, etc.)
      !real(r8),  allocatable   ::   soil_dzs(:)
      !real(r8),  allocatable   ::   soil_zs (:)
      real(r8),  allocatable   ::   snowc  (:, :)    ! flag indicating snow coverage (1 for snow cover)
      !real(r8),  allocatable   ::   tslb   (:, :, :) ! soil temperature in 5-layer model
      real(r8),  allocatable   ::   thc    (:, :)    ! thermal inertia (cal/cm/k/s^0.5)
      real(r8),  allocatable   ::   tmn    (:, :)    ! soil temperature at lower boundary (k)
      real(r8),  allocatable   ::   t2     (:, :)    ! diagnostic 2-m temperature from surface layer and lsm
      real(r8),  allocatable   ::   th2    (:, :)    ! diagnostic 2-m theta from surface layer and lsm
      !real(r8),  allocatable   ::   thz0   (:, :)    ! potential temperature at roughness length (k)
      real(r8),  allocatable   ::   tsk    (:, :)    ! surface temperature (k)

      real(r8),  allocatable   ::   xland  (:, :)    ! land mask (1 for land, 2 for water)
      real(r8),  allocatable   ::   xice   (:, :)    ! 
      real(r8),  allocatable   ::   snow   (:, :)    ! mm snow depth at surface 
      real(r8),  allocatable   ::   ust    (:, :)    ! u* in similarity theory (m/s)
      real(r8),  allocatable   ::   u10    (:, :)    ! diagnostic 10-m u component from surface layer
      real(r8),  allocatable   ::   v10    (:, :)    ! diagnostic 10-m v component from surface layer
      real(r8),  allocatable   ::   wspd   (:, :)    ! wind speed at lowest model level (m/s)
      real(r8),  allocatable   ::   znt    (:, :)    ! time-varying roughness length (m)
      real(r8),  allocatable   ::   zol    (:, :)    ! z/l height over monin-obukhov length
! rain diagnose
      real(r8),  allocatable   ::   rainc  (:, :) ! accumulated cumulus rain (mm)
      real(r8),  allocatable   ::   raincv (:, :) ! one-step cumulus rain (mm/step)
      real(r8),  allocatable   ::   rainnc (:, :) ! accumulate non-cumulus rain (mm)
      real(r8),  allocatable   ::   rainncv(:, :) ! one-step non-cumulus rain(mm/step)
      real(r8),  allocatable   ::   raint(:, :)   ! one-step total rain(mm/step)

      real(r8),  allocatable   ::   snownc (:, :)    ! accumulate snow (mm)
      real(r8),  allocatable   ::   snowncv(:, :)    ! one-step snow (mm/step)
      real(r8),  allocatable   ::   graupelnc (:, :) ! accumulate graupel (mm)
      real(r8),  allocatable   ::   graupelncv(:, :) ! one-step graupel (mm/step)
! static, gwdo
      real(r8),  allocatable   ::   var2d(:,:)
      real(r8),  allocatable   ::   oc12d(:,:)
      real(r8),  allocatable   ::   oa2d1(:,:)
      real(r8),  allocatable   ::   oa2d2(:,:)
      real(r8),  allocatable   ::   oa2d3(:,:)
      real(r8),  allocatable   ::   oa2d4(:,:)
      real(r8),  allocatable   ::   ol2d1(:,:)
      real(r8),  allocatable   ::   ol2d2(:,:)
      real(r8),  allocatable   ::   ol2d3(:,:)
      real(r8),  allocatable   ::   ol2d4(:,:)
!
! Below are diagnosed from CESM routine, only for ocean points; 
! for taux, tauy we will use CESM's information. For other points, this can be a good check to
! see what we have modified; because they are offline with WRF-physics, we may use scalar_1d_field to define them.
!
      type(scalar_1d_field)  :: taux
      type(scalar_1d_field)  :: tauy
      type(scalar_1d_field)  :: hfx_atmOcn
      type(scalar_1d_field)  :: qfx_atmOcn
      type(scalar_1d_field)  :: lh_atmOcn
      type(scalar_1d_field)  :: lwup_atmOcn
      type(scalar_1d_field)  :: taux_atmOcn
      type(scalar_1d_field)  :: tauy_atmOcn
      type(scalar_1d_field)  :: t2m_atmOcn
      type(scalar_1d_field)  :: q2m_atmOcn
      type(scalar_1d_field)  :: uu10m_atmOcn ! square
      type(scalar_1d_field)  :: ustar_atmOcn

    end type physics_surface_wrf

    type physics_tend_wrf
!
! 3d advection tend, for cumulus
!
      real(r8), allocatable   :: rqvdyten(:, :, :)
      real(r8), allocatable   :: rthdyten(:, :, :)
!
! cumulus
!
      real(r8), allocatable   :: ruucuten(:, :, :)
      real(r8), allocatable   :: rvvcuten(:, :, :)
      real(r8), allocatable   :: rthcuten(:, :, :)
      real(r8), allocatable   :: rqvcuten(:, :, :)
      real(r8), allocatable   :: rqccuten(:, :, :)
      real(r8), allocatable   :: rqrcuten(:, :, :)
      real(r8), allocatable   :: rqicuten(:, :, :)
      real(r8), allocatable   :: rqscuten(:, :, :)
!
! micro
!
      real(r8), allocatable   :: rthmpten(:, :, :)
      real(r8), allocatable   :: rqvmpten(:, :, :)
      real(r8), allocatable   :: rqcmpten(:, :, :)
      real(r8), allocatable   :: rqrmpten(:, :, :)
      real(r8), allocatable   :: rqimpten(:, :, :)
      real(r8), allocatable   :: rqsmpten(:, :, :)
      real(r8), allocatable   :: rqgmpten(:, :, :)
! for two-morr
      real(r8), allocatable   :: rnrmpten(:, :, :)
      real(r8), allocatable   :: rnimpten(:, :, :)
      real(r8), allocatable   :: rnsmpten(:, :, :)
      real(r8), allocatable   :: rngmpten(:, :, :)
!
! pbl
!
      real(r8), allocatable   :: ruublten(:, :, :)
      real(r8), allocatable   :: rvvblten(:, :, :)
      real(r8), allocatable   :: rthblten(:, :, :)
      real(r8), allocatable   :: rqvblten(:, :, :)
      real(r8), allocatable   :: rqcblten(:, :, :)
      real(r8), allocatable   :: rqiblten(:, :, :)
!
! radiation
!
      real(r8), allocatable   :: rthraten   (:, :, :)
      real(r8), allocatable   :: rthraten_lw(:, :, :)
      real(r8), allocatable   :: rthraten_sw(:, :, :)
 
    end type physics_tend_wrf
!
! define
!
    type(physics_state_wrf)    :: pstate_wrf
    type(physics_surface_wrf)  :: psurf_wrf
    type(physics_tend_wrf)     :: ptend_wrf

  contains

    subroutine grist_PhysW_data_structure_construct(ncell, nLevel, nspecies)
      integer(i4),  intent(in)  :: ncell, nLevel, nspecies
!
! first set tracer index
! default, we only have single-moment scheme
!
        if(nspecies.eq.6.or.nspecies.eq.5)then ! e.g., WSM6, linv381
           p_qv   = 1   ! species index for water vapor
           p_qc   = 2   ! species index for cloud water
           p_qr   = 3   ! species index for rain water
           p_qi   = 4   ! species index for cloud ice
           p_qs   = 5   ! species index for snow
           p_qg   = 6   ! species index for graupel
        end if
        if(nspecies.eq.3)then ! e.g., WSM6, linv381
           p_qv   = 1   ! species index for water vapor
           p_qc   = 2   ! species index for cloud water
           p_qi   = 3   ! species index for cloud ice
           p_qr   = 999 ! not use
           p_qs   = 999 ! not use
           p_qg   = 999 ! not use
        end if
!
! for two-moment scheme, we should not use nspecies to judge whether two or single moment;
! because some single moment scheme might support 7 species, if any
!
        if(trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381' .or. &
           trim(PhysW_mp_scheme).eq.'MORR_TWOM_V381_ACE' )then
           p_qv   = 1   ! species index for water vapor
           p_qc   = 2   ! species index for cloud water
           p_qr   = 3   ! species index for rain water
           p_qi   = 4   ! species index for cloud ice
           p_qs   = 5   ! species index for snow
           p_qg   = 6   ! species index for graupel
           p_ni   = 7
           p_ns   = 8
           p_nr   = 9
           p_ng   = 10
        end if

        restart = trim(run_type).eq.'restart'

        if(.not.allocated(ptend_wrf%rqvdyten))   allocate(ptend_wrf%rqvdyten(1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(ptend_wrf%rthdyten))   allocate(ptend_wrf%rthdyten(1:ncell, 1:nLevel, 1:1),source=zero)
!
! pstate_wrf
!
        if(.not.allocated(pstate_wrf%u_phy  ))   allocate(pstate_wrf%u_phy  (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%v_phy  ))   allocate(pstate_wrf%v_phy  (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%th_phy ))   allocate(pstate_wrf%th_phy (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%t_phy  ))   allocate(pstate_wrf%t_phy  (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%t8w    ))   allocate(pstate_wrf%t8w    (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%moist  ))   allocate(pstate_wrf%moist  (1:ncell, 1:nLevel, 1:1, nspecies),source=zero)
        if(.not.allocated(pstate_wrf%p_phy  ))   allocate(pstate_wrf%p_phy  (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%pi_phy ))   allocate(pstate_wrf%pi_phy (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%p8w    ))   allocate(pstate_wrf%p8w    (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%rhom   ))   allocate(pstate_wrf%rhom   (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%rhod   ))   allocate(pstate_wrf%rhod   (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%dz8w   ))   allocate(pstate_wrf%dz8w   (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%z8w    ))   allocate(pstate_wrf%z8w    (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%zzz    ))   allocate(pstate_wrf%zzz    (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%www    ))   allocate(pstate_wrf%www    (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%omega  ))   allocate(pstate_wrf%omega  (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%w0avg  ))   allocate(pstate_wrf%w0avg  (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%cmfmc  ))   allocate(pstate_wrf%cmfmc  (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%cldfra ))   allocate(pstate_wrf%cldfra (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%cldcum ))   allocate(pstate_wrf%cldcum (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%cldstr ))   allocate(pstate_wrf%cldstr (1:ncell, 1:nLevel, 1:1),source=zero)
        !if(.not.allocated(pstate_wrf%relhum ))   allocate(pstate_wrf%relhum (1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%xlat ))     allocate(pstate_wrf%xlat   (1:ncell, 1:1))
        if(.not.allocated(pstate_wrf%xlong ))    allocate(pstate_wrf%xlong  (1:ncell, 1:1))
        if(.not.allocated(pstate_wrf%dxmean))    allocate(pstate_wrf%dxmean (1:ncell))
        if(.not.allocated(pstate_wrf%coriolis))  allocate(pstate_wrf%coriolis (1:ncell))
        if(.not.allocated(pstate_wrf%kpbl))      allocate(pstate_wrf%kpbl   (1:ncell,1:1),source=1)
        if(.not.allocated(pstate_wrf%re_cloud))  allocate(pstate_wrf%re_cloud(1:ncell, 1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%re_ice  ))  allocate(pstate_wrf%re_ice (1:ncell,  1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%re_snow ))  allocate(pstate_wrf%re_snow(1:ncell,  1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%refl_10cm )) allocate(pstate_wrf%refl_10cm(1:ncell, 1:nLevel, 1:1),source=zero)
#ifdef CLIMATE_O3
        if(.not.allocated(pstate_wrf%o3    ))    allocate(pstate_wrf%o3     (1:ncell,  1:nLevel, 1:1),source=zero)
#endif

        if(PhysW_bl_scheme .eq. 'S&HV381' .or.    &
           PhysW_bl_scheme .eq. 'CAMUW_V381' )then
        if(.not.allocated(pstate_wrf%tke))      allocate(pstate_wrf%tke     (1:ncell,  1:nLevel, 1:1),source=zero)
        if(.not.allocated(pstate_wrf%leng))     allocate(pstate_wrf%leng    (1:ncell,  1:nLevel, 1:1),source=zero)
        end if
#ifdef UW_CAM5
        if(PhysW_bl_scheme .eq. 'CAMUW_V381' )then
           if(.not.allocated(pstate_wrf%kvm))      allocate(pstate_wrf%kvm     (1:ncell,  1:nLevel, 1:1),source=zero)
           if(.not.allocated(pstate_wrf%kvh))      allocate(pstate_wrf%kvh     (1:ncell,  1:nLevel, 1:1),source=zero)
           if(.not.allocated(pstate_wrf%tauresx2d))allocate(pstate_wrf%tauresx2d(1:ncell, 1:1),source=zero)
           if(.not.allocated(pstate_wrf%tauresy2d))allocate(pstate_wrf%tauresy2d(1:ncell, 1:1),source=zero)
        end if
#endif
        
        ! type data definition
        if(.not.allocated(psurf_wrf%taux%f))           allocate(psurf_wrf%taux%f       (1:ncell) ,source=zero); psurf_wrf%taux%pos         = 0
        if(.not.allocated(psurf_wrf%tauy%f))           allocate(psurf_wrf%tauy%f       (1:ncell) ,source=zero); psurf_wrf%tauy%pos         = 0
        if(.not.allocated(psurf_wrf%hfx_atmOcn%f))     allocate(psurf_wrf%hfx_atmOcn%f (1:ncell) ,source=zero); psurf_wrf%hfx_atmOcn%pos   = 0
        if(.not.allocated(psurf_wrf%qfx_atmOcn%f))     allocate(psurf_wrf%qfx_atmOcn%f (1:ncell) ,source=zero); psurf_wrf%qfx_atmOcn%pos   = 0
        if(.not.allocated(psurf_wrf%lh_atmOcn%f))      allocate(psurf_wrf%lh_atmOcn%f  (1:ncell) ,source=zero); psurf_wrf%lh_atmOcn%pos    = 0
        if(.not.allocated(psurf_wrf%lwup_atmOcn%f))    allocate(psurf_wrf%lwup_atmOcn%f(1:ncell) ,source=zero); psurf_wrf%lwup_atmOcn%pos  = 0
        if(.not.allocated(psurf_wrf%taux_atmOcn%f))    allocate(psurf_wrf%taux_atmOcn%f(1:ncell) ,source=zero); psurf_wrf%taux_atmOcn%pos  = 0
        if(.not.allocated(psurf_wrf%tauy_atmOcn%f))    allocate(psurf_wrf%tauy_atmOcn%f(1:ncell) ,source=zero); psurf_wrf%tauy_atmOcn%pos  = 0
        if(.not.allocated(psurf_wrf%t2m_atmOcn%f))     allocate(psurf_wrf%t2m_atmOcn%f (1:ncell) ,source=zero); psurf_wrf%t2m_atmOcn%pos   = 0
        if(.not.allocated(psurf_wrf%q2m_atmOcn%f))     allocate(psurf_wrf%q2m_atmOcn%f (1:ncell) ,source=zero); psurf_wrf%q2m_atmOcn%pos   = 0
        if(.not.allocated(psurf_wrf%uu10m_atmOcn%f))   allocate(psurf_wrf%uu10m_atmOcn%f(1:ncell),source=zero); psurf_wrf%uu10m_atmOcn%pos = 0
        if(.not.allocated(psurf_wrf%ustar_atmOcn%f))   allocate(psurf_wrf%ustar_atmOcn%f(1:ncell),source=zero); psurf_wrf%ustar_atmOcn%pos = 0

! energy flux
        call wrap_allocate_2d(var1=pstate_wrf%swupt ,var2=pstate_wrf%swuptc,var3=pstate_wrf%swdnt ,var4=pstate_wrf%swdntc,ncell=ncell)
        call wrap_allocate_2d(var1=pstate_wrf%lwupt ,var2=pstate_wrf%lwuptc,var3=pstate_wrf%lwdnt ,var4=pstate_wrf%lwdntc,ncell=ncell)
        call wrap_allocate_2d(var1=pstate_wrf%swupb ,var2=pstate_wrf%swupbc,var3=pstate_wrf%swdnb ,var4=pstate_wrf%swdnbc,ncell=ncell)
        call wrap_allocate_2d(var1=pstate_wrf%lwupb ,var2=pstate_wrf%lwupbc,var3=pstate_wrf%lwdnb ,var4=pstate_wrf%lwdnbc,ncell=ncell)
        call wrap_allocate_2d(var1=pstate_wrf%swcf  ,var2=pstate_wrf%lwcf  ,var3=pstate_wrf%olr   ,var4=pstate_wrf%coszr ,ncell=ncell)
        call wrap_allocate_2d(var1=pstate_wrf%cldlow,var2=pstate_wrf%cldmed,var3=pstate_wrf%cldhgh,var4=pstate_wrf%cldtot,ncell=ncell)
        call wrap_allocate_3d(var1=pstate_wrf%cemiss,var2=pstate_wrf%taucldc,var3=pstate_wrf%taucldi,ncell=ncell, nLevel=nLevel)
!
! psurf_wrf, refactor later
!
        if(.not.allocated(psurf_wrf%albedo ))    allocate(psurf_wrf%albedo (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%asdir  ))    allocate(psurf_wrf%asdir  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%asdif  ))    allocate(psurf_wrf%asdif  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%aldir  ))    allocate(psurf_wrf%aldir  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%aldif  ))    allocate(psurf_wrf%aldif  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%br     ))    allocate(psurf_wrf%br     (1:ncell, 1:1),source=zero)
        !if(.not.allocated(psurf_wrf%capg   ))    allocate(psurf_wrf%capg   (1:ncell, 1:1),source=zero)
        !if(.not.allocated(psurf_wrf%chklowq))    allocate(psurf_wrf%chklowq(1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%emiss  ))    allocate(psurf_wrf%emiss  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%flhc   ))    allocate(psurf_wrf%flhc   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%flqc   ))    allocate(psurf_wrf%flqc   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%glw    ))    allocate(psurf_wrf%glw    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%gsw    ))    allocate(psurf_wrf%gsw    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%swdown ))    allocate(psurf_wrf%swdown (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%gz1oz0 ))    allocate(psurf_wrf%gz1oz0 (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%hfx    ))    allocate(psurf_wrf%hfx    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%ht     ))    allocate(psurf_wrf%ht     (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%lh     ))    allocate(psurf_wrf%lh     (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%mavail ))    allocate(psurf_wrf%mavail (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%mol    ))    allocate(psurf_wrf%mol    (1:ncell, 1:1),source=zero)
                                                                             
        if(.not.allocated(psurf_wrf%pblh   ))    allocate(psurf_wrf%pblh   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%psih   ))    allocate(psurf_wrf%psih   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%psim   ))    allocate(psurf_wrf%psim   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%psfc   ))    allocate(psurf_wrf%psfc   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%q2     ))    allocate(psurf_wrf%q2     (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%qfx    ))    allocate(psurf_wrf%qfx    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%qsfc   ))    allocate(psurf_wrf%qsfc   (1:ncell, 1:1),source=zero)
                                    
        if(.not.allocated(psurf_wrf%regime  ))   allocate(psurf_wrf%regime (1:ncell, 1:1),source=zero)
        !if(.not.allocated(psurf_wrf%soil_dzs))   allocate(psurf_wrf%soil_dzs(1:num_soil_layers),source=zero)
        !if(.not.allocated(psurf_wrf%soil_zs ))   allocate(psurf_wrf%soil_zs (1:num_soil_layers),source=zero)
        if(.not.allocated(psurf_wrf%snowc   ))   allocate(psurf_wrf%snowc  (1:ncell, 1:1),source=zero)
        !if(.not.allocated(psurf_wrf%tslb    ))   allocate(psurf_wrf%tslb(   1:ncell, 1:num_soil_layers, 1:1 ))
        if(.not.allocated(psurf_wrf%thc     ))   allocate(psurf_wrf%thc    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%tmn     ))   allocate(psurf_wrf%tmn    (1:ncell, 1:1),source=zero) 
        if(.not.allocated(psurf_wrf%t2      ))   allocate(psurf_wrf%t2     (1:ncell, 1:1),source=zero) 
        if(.not.allocated(psurf_wrf%th2     ))   allocate(psurf_wrf%th2    (1:ncell, 1:1),source=zero) 
        !if(.not.allocated(psurf_wrf%thz0    ))   allocate(psurf_wrf%thz0   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%tsk     ))   allocate(psurf_wrf%tsk    (1:ncell, 1:1),source=zero)   
                                   
        if(.not.allocated(psurf_wrf%xland   ))   allocate(psurf_wrf%xland  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%xice    ))   allocate(psurf_wrf%xice   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%snow    ))   allocate(psurf_wrf%snow   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%ust     ))   allocate(psurf_wrf%ust    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%u10     ))   allocate(psurf_wrf%u10    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%v10     ))   allocate(psurf_wrf%v10    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%wspd    ))   allocate(psurf_wrf%wspd   (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%znt     ))   allocate(psurf_wrf%znt    (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%zol     ))   allocate(psurf_wrf%zol    (1:ncell, 1:1),source=zero)

        if(.not.allocated(psurf_wrf%rainc  ))    allocate(psurf_wrf%rainc  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%raincv ))    allocate(psurf_wrf%raincv (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%rainnc ))    allocate(psurf_wrf%rainnc (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%rainncv))    allocate(psurf_wrf%rainncv(1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%raint  ))    allocate(psurf_wrf%raint  (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%snownc ))    allocate(psurf_wrf%snownc (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%snowncv))    allocate(psurf_wrf%snowncv(1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%graupelnc )) allocate(psurf_wrf%graupelnc (1:ncell, 1:1),source=zero)
        if(.not.allocated(psurf_wrf%graupelncv)) allocate(psurf_wrf%graupelncv(1:ncell, 1:1),source=zero)

        if(.not.allocated(psurf_wrf%var2d))      allocate(psurf_wrf%var2d(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%oc12d))      allocate(psurf_wrf%oc12d(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%oa2d1))      allocate(psurf_wrf%oa2d1(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%oa2d2))      allocate(psurf_wrf%oa2d2(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%oa2d3))      allocate(psurf_wrf%oa2d3(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%oa2d4))      allocate(psurf_wrf%oa2d4(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%ol2d1))      allocate(psurf_wrf%ol2d1(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%ol2d2))      allocate(psurf_wrf%ol2d2(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%ol2d3))      allocate(psurf_wrf%ol2d3(1:ncell,1:1),source=zero)
        if(.not.allocated(psurf_wrf%ol2d4))      allocate(psurf_wrf%ol2d4(1:ncell,1:1),source=zero)

      return
    end subroutine grist_PhysW_data_structure_construct

    subroutine grist_PhysW_data_structure_destruct()

        if(allocated(ptend_wrf%rqvdyten))   deallocate(ptend_wrf%rqvdyten)
        if(allocated(ptend_wrf%rthdyten))   deallocate(ptend_wrf%rthdyten)
!
! pstate_wrf
!
        if(allocated(pstate_wrf%u_phy  ))   deallocate(pstate_wrf%u_phy  )
        if(allocated(pstate_wrf%v_phy  ))   deallocate(pstate_wrf%v_phy  )
        if(allocated(pstate_wrf%th_phy ))   deallocate(pstate_wrf%th_phy )
        if(allocated(pstate_wrf%t_phy  ))   deallocate(pstate_wrf%t_phy  )
        if(allocated(pstate_wrf%t8w    ))   deallocate(pstate_wrf%t8w    )
        if(allocated(pstate_wrf%moist  ))   deallocate(pstate_wrf%moist  )
        if(allocated(pstate_wrf%p_phy  ))   deallocate(pstate_wrf%p_phy  )
        if(allocated(pstate_wrf%pi_phy ))   deallocate(pstate_wrf%pi_phy )
        if(allocated(pstate_wrf%p8w    ))   deallocate(pstate_wrf%p8w    )
        if(allocated(pstate_wrf%rhom   ))   deallocate(pstate_wrf%rhom   )
        if(allocated(pstate_wrf%rhod   ))   deallocate(pstate_wrf%rhod   )
        if(allocated(pstate_wrf%dz8w   ))   deallocate(pstate_wrf%dz8w   )
        if(allocated(pstate_wrf%z8w    ))   deallocate(pstate_wrf%z8w    )
        if(allocated(pstate_wrf%zzz    ))   deallocate(pstate_wrf%zzz    )
        if(allocated(pstate_wrf%www    ))   deallocate(pstate_wrf%www    )
        if(allocated(pstate_wrf%omega  ))   deallocate(pstate_wrf%omega  )
        if(allocated(pstate_wrf%w0avg  ))   deallocate(pstate_wrf%w0avg  )
        if(allocated(pstate_wrf%cmfmc  ))   deallocate(pstate_wrf%cmfmc  )
        if(allocated(pstate_wrf%cldfra ))   deallocate(pstate_wrf%cldfra )
        if(allocated(pstate_wrf%cldcum ))   deallocate(pstate_wrf%cldcum )
        if(allocated(pstate_wrf%cldstr ))   deallocate(pstate_wrf%cldstr )
        !if(allocated(pstate_wrf%relhum ))   deallocate(pstate_wrf%relhum )
        if(allocated(pstate_wrf%xlat   ))   deallocate(pstate_wrf%xlat   )
        if(allocated(pstate_wrf%xlong  ))   deallocate(pstate_wrf%xlong  )
        if(allocated(pstate_wrf%dxmean ))   deallocate(pstate_wrf%dxmean )
        if(allocated(pstate_wrf%coriolis )) deallocate(pstate_wrf%coriolis)
        if(allocated(pstate_wrf%kpbl   ))   deallocate(pstate_wrf%kpbl   )
        if(allocated(pstate_wrf%re_cloud )) deallocate(pstate_wrf%re_cloud)
        if(allocated(pstate_wrf%re_ice ))   deallocate(pstate_wrf%re_ice )
        if(allocated(pstate_wrf%re_snow ))  deallocate(pstate_wrf%re_snow)
        if(allocated(pstate_wrf%refl_10cm ))deallocate(pstate_wrf%refl_10cm)
        if(allocated(pstate_wrf%tke    ))   deallocate(pstate_wrf%tke    )
        if(allocated(pstate_wrf%leng   ))   deallocate(pstate_wrf%leng   )
#ifdef CLIMATE_O3
        if(allocated(pstate_wrf%o3     ))   deallocate(pstate_wrf%o3     )
#endif
#ifdef UW_CAM5
        if(allocated(pstate_wrf%kvm))       deallocate(pstate_wrf%kvm    )
        if(allocated(pstate_wrf%kvh))       deallocate(pstate_wrf%kvh    )
        if(allocated(pstate_wrf%tauresx2d)) deallocate(pstate_wrf%tauresx2d)
        if(allocated(pstate_wrf%tauresy2d)) deallocate(pstate_wrf%tauresy2d)
#endif
        
        if(allocated(psurf_wrf%taux%f))          deallocate(psurf_wrf%taux%f        )
        if(allocated(psurf_wrf%tauy%f))          deallocate(psurf_wrf%tauy%f        )
        if(allocated(psurf_wrf%hfx_atmOcn%f))    deallocate(psurf_wrf%hfx_atmOcn%f  )
        if(allocated(psurf_wrf%qfx_atmOcn%f))    deallocate(psurf_wrf%qfx_atmOcn%f  )
        if(allocated(psurf_wrf%lh_atmOcn%f))     deallocate(psurf_wrf%lh_atmOcn%f   )
        if(allocated(psurf_wrf%lwup_atmOcn%f))   deallocate(psurf_wrf%lwup_atmOcn%f )
        if(allocated(psurf_wrf%taux_atmOcn%f))   deallocate(psurf_wrf%taux_atmOcn%f )
        if(allocated(psurf_wrf%tauy_atmOcn%f))   deallocate(psurf_wrf%tauy_atmOcn%f )
        if(allocated(psurf_wrf%t2m_atmOcn%f))    deallocate(psurf_wrf%t2m_atmOcn%f  )
        if(allocated(psurf_wrf%q2m_atmOcn%f))    deallocate(psurf_wrf%q2m_atmOcn%f  )
        if(allocated(psurf_wrf%uu10m_atmOcn%f))  deallocate(psurf_wrf%uu10m_atmOcn%f)
        if(allocated(psurf_wrf%ustar_atmOcn%f))  deallocate(psurf_wrf%ustar_atmOcn%f)
! energy flux
        call wrap_deallocate_2d(var1=pstate_wrf%swupt ,var2=pstate_wrf%swuptc ,var3=pstate_wrf%swdnt ,var4=pstate_wrf%swdntc)
        call wrap_deallocate_2d(var1=pstate_wrf%lwupt ,var2=pstate_wrf%lwuptc ,var3=pstate_wrf%lwdnt ,var4=pstate_wrf%lwdntc)
        call wrap_deallocate_2d(var1=pstate_wrf%swupb ,var2=pstate_wrf%swupbc ,var3=pstate_wrf%swdnb ,var4=pstate_wrf%swdnbc)
        call wrap_deallocate_2d(var1=pstate_wrf%lwupb ,var2=pstate_wrf%lwupbc ,var3=pstate_wrf%lwdnb ,var4=pstate_wrf%lwdnbc)
        call wrap_deallocate_2d(var1=pstate_wrf%swcf  ,var2=pstate_wrf%lwcf   ,var3=pstate_wrf%olr   ,var4=pstate_wrf%coszr )
        call wrap_deallocate_3d(var1=pstate_wrf%cemiss,var2=pstate_wrf%taucldc,var3=pstate_wrf%taucldi)
        call wrap_deallocate_2d(var1=pstate_wrf%cldtot,var2=pstate_wrf%cldlow ,var3=pstate_wrf%cldmed,var4=pstate_wrf%cldhgh)
!
! psurf_wrf
!
        if(allocated(psurf_wrf%albedo ))    deallocate(psurf_wrf%albedo )
        if(allocated(psurf_wrf%asdir  ))    deallocate(psurf_wrf%asdir  )
        if(allocated(psurf_wrf%asdif  ))    deallocate(psurf_wrf%asdif  )
        if(allocated(psurf_wrf%aldir  ))    deallocate(psurf_wrf%aldir  )
        if(allocated(psurf_wrf%aldif  ))    deallocate(psurf_wrf%aldif  )
        if(allocated(psurf_wrf%br     ))    deallocate(psurf_wrf%br     )
        !if(allocated(psurf_wrf%capg   ))    deallocate(psurf_wrf%capg   )
        !if(allocated(psurf_wrf%chklowq))    deallocate(psurf_wrf%chklowq)
        if(allocated(psurf_wrf%emiss  ))    deallocate(psurf_wrf%emiss  )
        if(allocated(psurf_wrf%flhc   ))    deallocate(psurf_wrf%flhc   )
        if(allocated(psurf_wrf%flqc   ))    deallocate(psurf_wrf%flqc   )
        if(allocated(psurf_wrf%glw    ))    deallocate(psurf_wrf%glw    )
        if(allocated(psurf_wrf%gsw    ))    deallocate(psurf_wrf%gsw    )
        if(allocated(psurf_wrf%swdown ))    deallocate(psurf_wrf%swdown )
        if(allocated(psurf_wrf%gz1oz0 ))    deallocate(psurf_wrf%gz1oz0 )
        if(allocated(psurf_wrf%hfx    ))    deallocate(psurf_wrf%hfx    )
        if(allocated(psurf_wrf%ht     ))    deallocate(psurf_wrf%ht     )
        if(allocated(psurf_wrf%lh     ))    deallocate(psurf_wrf%lh     )
        if(allocated(psurf_wrf%mavail ))    deallocate(psurf_wrf%mavail )
        if(allocated(psurf_wrf%mol    ))    deallocate(psurf_wrf%mol    )
                                                                       
        if(allocated(psurf_wrf%pblh   ))    deallocate(psurf_wrf%pblh   )
        if(allocated(psurf_wrf%psih   ))    deallocate(psurf_wrf%psih   )
        if(allocated(psurf_wrf%psim   ))    deallocate(psurf_wrf%psim   )
        if(allocated(psurf_wrf%psfc   ))    deallocate(psurf_wrf%psfc   )
        if(allocated(psurf_wrf%q2     ))    deallocate(psurf_wrf%q2     )
        if(allocated(psurf_wrf%qfx    ))    deallocate(psurf_wrf%qfx    )
        if(allocated(psurf_wrf%qsfc   ))    deallocate(psurf_wrf%qsfc   )
                              
        if(allocated(psurf_wrf%regime  ))   deallocate(psurf_wrf%regime )
        !if(allocated(psurf_wrf%soil_dzs))   deallocate(psurf_wrf%soil_dzs)
        !if(allocated(psurf_wrf%soil_zs ))   deallocate(psurf_wrf%soil_zs)
        if(allocated(psurf_wrf%snowc   ))   deallocate(psurf_wrf%snowc  )
        !if(allocated(psurf_wrf%tslb    ))   deallocate(psurf_wrf%tslb   )
        if(allocated(psurf_wrf%thc     ))   deallocate(psurf_wrf%thc    )
        if(allocated(psurf_wrf%tmn     ))   deallocate(psurf_wrf%tmn    ) 
        if(allocated(psurf_wrf%t2      ))   deallocate(psurf_wrf%t2     ) 
        if(allocated(psurf_wrf%th2     ))   deallocate(psurf_wrf%th2    ) 
        !if(allocated(psurf_wrf%thz0    ))   deallocate(psurf_wrf%thz0   )
        if(allocated(psurf_wrf%tsk     ))   deallocate(psurf_wrf%tsk    )   
                              
        if(allocated(psurf_wrf%xland   ))   deallocate(psurf_wrf%xland  )
        if(allocated(psurf_wrf%xice    ))   deallocate(psurf_wrf%xice   )
        if(allocated(psurf_wrf%snow    ))   deallocate(psurf_wrf%snow   )
        if(allocated(psurf_wrf%ust     ))   deallocate(psurf_wrf%ust    )
        if(allocated(psurf_wrf%u10     ))   deallocate(psurf_wrf%u10    )
        if(allocated(psurf_wrf%v10     ))   deallocate(psurf_wrf%v10    )
        if(allocated(psurf_wrf%wspd    ))   deallocate(psurf_wrf%wspd   )
        if(allocated(psurf_wrf%znt     ))   deallocate(psurf_wrf%znt    )
        if(allocated(psurf_wrf%zol     ))   deallocate(psurf_wrf%zol    )
! rain
        if(allocated(psurf_wrf%rainc  ))   deallocate(psurf_wrf%rainc  ) ! accumulated
        if(allocated(psurf_wrf%raincv ))   deallocate(psurf_wrf%raincv ) ! one-step
        if(allocated(psurf_wrf%rainnc ))   deallocate(psurf_wrf%rainnc )
        if(allocated(psurf_wrf%rainncv))   deallocate(psurf_wrf%rainncv)
        if(allocated(psurf_wrf%snownc ))   deallocate(psurf_wrf%snownc )
        if(allocated(psurf_wrf%snowncv))   deallocate(psurf_wrf%snowncv)
        if(allocated(psurf_wrf%graupelnc ))deallocate(psurf_wrf%graupelnc )
        if(allocated(psurf_wrf%graupelncv))deallocate(psurf_wrf%graupelncv)
        if(allocated(psurf_wrf%raint))     deallocate(psurf_wrf%raint) ! diagnosed
! static, gwdo
        if(allocated(psurf_wrf%var2d))     deallocate(psurf_wrf%var2d)
        if(allocated(psurf_wrf%oc12d))     deallocate(psurf_wrf%oc12d)
        if(allocated(psurf_wrf%oa2d1))     deallocate(psurf_wrf%oa2d1)
        if(allocated(psurf_wrf%oa2d2))     deallocate(psurf_wrf%oa2d2)
        if(allocated(psurf_wrf%oa2d3))     deallocate(psurf_wrf%oa2d3)
        if(allocated(psurf_wrf%oa2d4))     deallocate(psurf_wrf%oa2d4)
        if(allocated(psurf_wrf%ol2d1))     deallocate(psurf_wrf%ol2d1)
        if(allocated(psurf_wrf%ol2d2))     deallocate(psurf_wrf%ol2d2)
        if(allocated(psurf_wrf%ol2d3))     deallocate(psurf_wrf%ol2d3)
        if(allocated(psurf_wrf%ol2d4))     deallocate(psurf_wrf%ol2d4)

      return
    end subroutine grist_PhysW_data_structure_destruct

    subroutine wrap_allocate_2d(var1,var2,var3,var4,ncell)
      real(r8), allocatable, optional,intent(inout) :: var1(:,:)
      real(r8), allocatable, optional,intent(inout) :: var2(:,:)
      real(r8), allocatable, optional,intent(inout) :: var3(:,:)
      real(r8), allocatable, optional,intent(inout) :: var4(:,:)
      integer(i4),  intent(in) ::  ncell
      
      if(present(var1).and..not.allocated(var1)) allocate(var1(1:ncell,1:1))
      if(present(var2).and..not.allocated(var2)) allocate(var2(1:ncell,1:1))
      if(present(var3).and..not.allocated(var3)) allocate(var3(1:ncell,1:1))
      if(present(var4).and..not.allocated(var4)) allocate(var4(1:ncell,1:1))

      if(present(var1)) var1 = zero
      if(present(var2)) var2 = zero
      if(present(var3)) var3 = zero
      if(present(var4)) var4 = zero

      return
    end subroutine wrap_allocate_2d

    subroutine wrap_allocate_3d(var1,var2,var3,var4,ncell,nLevel)
      real(r8), allocatable, optional,intent(inout) :: var1(:,:,:)
      real(r8), allocatable, optional,intent(inout) :: var2(:,:,:)
      real(r8), allocatable, optional,intent(inout) :: var3(:,:,:)
      real(r8), allocatable, optional,intent(inout) :: var4(:,:,:)
      integer(i4),  intent(in) ::  ncell, nLevel
      
      if(present(var1).and..not.allocated(var1)) allocate(var1(1:ncell,1:nLevel,1:1))
      if(present(var2).and..not.allocated(var2)) allocate(var2(1:ncell,1:nLevel,1:1))
      if(present(var3).and..not.allocated(var3)) allocate(var3(1:ncell,1:nLevel,1:1))
      if(present(var4).and..not.allocated(var4)) allocate(var4(1:ncell,1:nLevel,1:1))

      if(present(var1)) var1 = zero
      if(present(var2)) var2 = zero
      if(present(var3)) var3 = zero
      if(present(var4)) var4 = zero

      return
    end subroutine wrap_allocate_3d

    subroutine wrap_deallocate_2d(var1,var2,var3,var4)
      real(r8), allocatable, optional,intent(inout) :: var1(:,:)
      real(r8), allocatable, optional,intent(inout) :: var2(:,:)
      real(r8), allocatable, optional,intent(inout) :: var3(:,:)
      real(r8), allocatable, optional,intent(inout) :: var4(:,:)
      
      if(present(var1).and.allocated(var1)) deallocate(var1)
      if(present(var2).and.allocated(var2)) deallocate(var2)
      if(present(var3).and.allocated(var3)) deallocate(var3)
      if(present(var4).and.allocated(var4)) deallocate(var4)

      return
    end subroutine wrap_deallocate_2d

    subroutine wrap_deallocate_3d(var1,var2,var3,var4)
      real(r8), allocatable, optional,intent(inout) :: var1(:,:,:)
      real(r8), allocatable, optional,intent(inout) :: var2(:,:,:)
      real(r8), allocatable, optional,intent(inout) :: var3(:,:,:)
      real(r8), allocatable, optional,intent(inout) :: var4(:,:,:)
      
      if(present(var1).and.allocated(var1)) deallocate(var1)
      if(present(var2).and.allocated(var2)) deallocate(var2)
      if(present(var3).and.allocated(var3)) deallocate(var3)
      if(present(var4).and.allocated(var4)) deallocate(var4)

      return
    end subroutine wrap_deallocate_3d

 end module grist_PhysW_data_structure
