!======================================================
!
!  Created by LiXiaohan on 19/4/10.
!  GRIST_SCM specific routines and vars
!======================================================

 module grist_scm_comm_module
    use grist_constants,                only: i4, r8
    use grist_nml_module,               only: nlev, ntracer, start_ymd, start_tod
    use grist_physics_data_structure,   only: tground, &
                                              have_lhflx, have_shflx, have_tg

    implicit none
    public

    logical, parameter :: use_phys_vars = .true.

    ! variables in grist_scm namelist:
!    integer(i4)                :: start_ymd            ! start date in yyyymmdd format
!    integer(i4)                :: start_tod            ! time of day relative to start date
    real(r8)                   :: scm_lat              ! latitude of the single model
    real(r8)                   :: scm_lon              ! longitude of the single model
    logical                    :: scm_relaxation       ! do relaxation
    character(1024)            :: scm_file_name        ! dynamical forcing data filename
    character(1024)            :: scam_file_name       ! dynamical forcing data from 3D CAM to supplement high-levels
    character(1024)            :: scm_test_name

    ! variables in SCM Module:
    integer(i4)                :: io_time_index        ! time index of file reading
    logical                    :: doioupdate
    real(r8)                   :: tsair                ! air temperature at the surface
    logical                    :: have_tsair
!    real(r8)                   :: tground              ! ground temperature
!    logical                    :: have_tg
    real(r8)                   :: psobs                ! actual surface pressure
    real(r8)                   :: lhflxobs             ! observed surface latent heat flux
!    logical                    :: have_lhflx
    real(r8)                   :: shflxobs             ! observed surface sensible heat flux
!    logical                    :: have_shflx
    logical                    :: have_alb             ! observed surface albedo 
    real(r8)                   :: albobs               ! observed surface albedo
    real(r8)                   :: tpertobs             ! observed pertubation temperature
    real(r8)                   :: qpertobs             ! observed pertubation specific humidity
    logical                    :: have_tpert
    logical                    :: have_qpert
    real(r8), allocatable      :: tobs(:)              ! actual temperature
    real(r8), allocatable      :: qobs(:)              ! actual water vapor Mixing ratio
    real(r8), allocatable      :: cldliqobs(:)
    real(r8), allocatable      :: cldiceobs(:)
    real(r8), allocatable      :: cldobs(:)            ! observed cloud
    real(r8), allocatable      :: concldobs(:)         ! observed convective cloud
    logical                    :: have_cld
    logical                    :: have_concld
    real(r8), allocatable      :: geoobs(:)            ! geopotential_height
    logical                    :: have_geo
    real(r8), allocatable      :: clwpobs(:)           ! observed clwp
    logical                    :: have_clwp
    real(r8), allocatable      :: wfld(:)              ! vertical pressure velocity at model levels
    real(r8), allocatable      :: wfldint(:)           ! vertical pressure velocity at face levels
    real(r8), allocatable      :: uobs(:)              ! actual u wind
    real(r8), allocatable      :: vobs(:)              ! actual v wind
    real(r8), allocatable      :: divq3d(:)            ! 3D q advective tendency
    real(r8), allocatable      :: divq(:)              ! horizontal q advective tendency
    real(r8), allocatable      :: vertdivq(:)          ! vertical q advective tendency
    real(r8), allocatable      :: divt3d(:)            ! 3D temperature advective
    real(r8), allocatable      :: divt(:)              ! horizontal temperature advective tendency
    real(r8), allocatable      :: vertdivt(:)          ! vertical temperature advective tendency
    logical                    :: have_divq
    logical                    :: have_vertdivq
    logical                    :: have_divt
    logical                    :: have_vertdivt
    logical                    :: use_3dfrc            ! 3D q and T advective tendency
    real(r8)                   :: ptend                ! surface pressure tendency
#ifdef AMIPW_PHYSICS
    real(r8), allocatable      :: t_tend_dycore(:)
    real(r8), allocatable      :: q_tend_dycore(:)
#endif
    contains

    subroutine grist_scm_vars_construct
        allocate(tobs(nlev));tobs=0.
        allocate(qobs(nlev));qobs=0.
        allocate(cldliqobs(nlev));cldliqobs=0.
        allocate(cldiceobs(nlev));cldiceobs=0.
        allocate(cldobs(nlev));cldobs=0.
        allocate(concldobs(nlev));concldobs=0.
        allocate(clwpobs(nlev));clwpobs=0.
        allocate(wfld(nlev));wfld=0.
        allocate(wfldint(nlev+1));wfldint=0.
        allocate(uobs(nlev));uobs=0.
        allocate(vobs(nlev));vobs=0.
        allocate(geoobs(nlev));geoobs=0.
        allocate(divq3d(nlev));divq3d=0.
        allocate(divq(nlev));divq=0.
        allocate(vertdivq(nlev));vertdivq=0.
        allocate(divt3d(nlev));divt3d=0.
        allocate(divt(nlev));divt=0.
        allocate(vertdivt(nlev));vertdivt=0.
#ifdef AMIPW_PHYSICS
        allocate(t_tend_dycore(nlev));t_tend_dycore=0.
        allocate(q_tend_dycore(nlev));q_tend_dycore=0.
#endif
    end subroutine grist_scm_vars_construct

    subroutine grist_scm_vars_destruct
        deallocate(tobs)
        deallocate(qobs)
        deallocate(cldliqobs)
        deallocate(cldiceobs)
        deallocate(cldobs)
        deallocate(concldobs)
        deallocate(clwpobs)
        deallocate(wfld)
        deallocate(wfldint)
        deallocate(uobs)
        deallocate(vobs)
        deallocate(geoobs)
        deallocate(divq3d)
        deallocate(divq)
        deallocate(vertdivq)
        deallocate(divt3d)
        deallocate(divt)
        deallocate(vertdivt)

#ifdef AMIPW_PHYSICS
        deallocate(t_tend_dycore)
        deallocate(q_tend_dycore)
#endif
    end subroutine grist_scm_vars_destruct

 end module grist_scm_comm_module
