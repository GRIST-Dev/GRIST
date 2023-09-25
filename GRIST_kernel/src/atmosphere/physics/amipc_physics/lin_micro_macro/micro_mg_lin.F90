!===================================================================================================
!
!  Created by LiXiaohan on 20/08/10
!
!  interface for MG_LIN microphysics
!
!===================================================================================================

  module micro_mg_lin

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8,      i4,            &
                                                  gravity, rdry,   rvap,  &
                                                  cp,      rhoh2o, tmelt, &
                                                  latvap,  latice, pi
 
    use grist_nml_module,                   only: nlev, nlevp, ntracer
    use grist_mpi

    implicit none
    private
    save
    
   
    public :: read_nml_micro_mg_lin,    &
              micro_mg_init_lin,        &
              micro_mg_tend_lin,        &
              micro_mg_end_lin

    integer :: micro_mg_version     = 1      ! Version number for MG.
    integer :: micro_mg_sub_version = 5      ! Second part of version number.
    
    logical :: microp_uniform = .false.
    
    logical, public :: do_cldliq ! Prognose cldliq flag
    logical, public :: do_cldice ! Prognose cldice flag
    
    integer, parameter :: ncnst = 4       ! Number of constituents
!    character(len=8), parameter :: &      ! Constituent names
!       cnst_names(ncnst) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE'/)
    
    integer :: ixcldliq,      &! cloud liquid amount index
               ixcldice,      &! cloud ice amount index
               ixnumliq,      &! cloud liquid number index
               ixnumice        ! cloud ice water index

! if the variables below are used by other modules (pbuf in CAM), change to pstate%~, Lixh
    real(r8),   allocatable     :: nevapr(:,:)      ! 'NEVAPR' Evaporation of total precipitation (rain + snow)
    real(r8),   allocatable     :: prain(:,:)       ! 'PRAIN' Total precipitation (rain + snow)
    real(r8),   allocatable     :: tnd_qsnow(:,:)   ! 'TND_QSNOW' external tendency on snow mass (kg/kg/s)
    real(r8),   allocatable     :: tnd_nsnow(:,:)   ! 'TND_NSNOW' external tendency on snow number(#/kg/s)
    real(r8),   allocatable     :: re_ice(:,:)      ! 'RE_ICE' ice effective radius (m)
    real(r8),   allocatable     :: relvar(:,:)      ! 'RELVAR' relative variance of cloud water
    real(r8),   allocatable     :: accre_enhan(:,:) ! 'ACCRE_ENHAN' optional accretion enhancement for experimentation
    real(r8),   allocatable     :: iciwpst(:,:)     ! 'ICIWPST' Stratiform in-cloud ice water path for radiation
    real(r8),   allocatable     :: iclwpst(:,:)     ! 'ICLWPST' Stratiform in-cloud liduid water path for radiation
    real(r8),   allocatable     :: acprecl(:)       ! 'ACPRECL' accumulated precip across timesteps
    real(r8),   allocatable     :: acgcme(:)        ! 'ACGCME' accumulated condensation across timesteps
    integer,    allocatable     :: acnum(:)         ! 'ACNUM' counter for # timesteps accumulated 


! other pbuf includes:
! physics buffer fields for COSP simulator:
! mgflxprc, mgflxsnw, mgmrprc, mgmrsnw, mgreffrain, mgreffsnow, cvreffliq, cvreffice

! physics buffer fields for CARMA and AEROSOL:
! qme

    contains

    subroutine read_nml_micro_mg_lin(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer :: unitn, ierr
    logical :: micro_mg_do_cldice   = .true. ! do_cldice = .true., MG microphysics is prognosing cldice
    logical :: micro_mg_do_cldliq   = .true. ! do_cldliq = .true., MG microphysics is prognosing cldliq
  
  
    namelist /micro_mg_nl/ micro_mg_version, micro_mg_sub_version, &
                           micro_mg_do_cldice, micro_mg_do_cldliq
    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, micro_mg_nl, iostat=ierr)
    if (ierr /= 0) call endrun(" error reading micro_mg namelist")
    close(unitn)
 
    ! set local variables
    do_cldice  = micro_mg_do_cldice
    do_cldliq  = micro_mg_do_cldliq

    ! Verify that version numbers are valid.
    select case (micro_mg_version)
    case (1)
       select case (micro_mg_sub_version)
       case(5)
          ! MG version 1.5 - MG2 development
       case default
          call endrun("Invalid sub-version number specified for MG_LIN microphysics")
       end select
    case default
       call endrun("Invalid version number specified for MG_LIN microphysics")
    end select

    end subroutine read_nml_micro_mg_lin


    subroutine micro_mg_init_lin(ncol)

    use micro_mg1_0,                   only: micro_mg1_0_init
    use micro_mg1_5_lin,               only: micro_mg1_5_init_lin
    use cldwat2m_macro,                only: rhmini
    use grist_physics_data_structure,  only: phy_tracer_info, pstate 
    use grist_PhysC_data_structure,    only: pstate_cam
    use grist_physics_update,          only: ptimelevels
! io
   integer, intent(in)  :: ncol
! local
   integer :: m

    ! inquire cloud index
    do m = 1, ntracer
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')        ixcldliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')           ixcldice = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid_number') ixnumliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice_number')    ixnumice = m 
    end do


    select case (micro_mg_version)
    case (1)
       select case (micro_mg_sub_version)
       case(5)
          ! MG version 1.5 - MG2 development
          call micro_mg1_5_init_lin(gravity, rdry,  rvap,  cp,               &
                                    tmelt,   latvap,latice,rhmini,           &
                                    microp_uniform, do_cldice)
       case default
       end select
    case default
    end select

    allocate(nevapr(nlev,ncol));         nevapr      = 0._r8
    allocate(prain(nlev,ncol));          prain       = 0._r8
    allocate(tnd_qsnow(nlev,ncol));      tnd_qsnow   = 0._r8
    allocate(tnd_nsnow(nlev,ncol));      tnd_nsnow   = 0._r8
    allocate(re_ice(nlev,ncol));         re_ice      = 0._r8
    allocate(relvar(nlev,ncol));         relvar      = 2._r8
    allocate(accre_enhan(nlev,ncol));    accre_enhan = 1._r8
    allocate(iciwpst(nlev,ncol));        iciwpst     = 0._r8
    allocate(iclwpst(nlev,ncol));        iclwpst     = 0._r8
    allocate(acprecl(ncol));             acprecl     = 0._r8
    allocate(acgcme(ncol));              acgcme      = 0._r8
    allocate(acnum(ncol));               acnum       = 0

    allocate(pstate_cam%microp_wsedl_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_dei_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_des_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_res_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_mu_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_lambdac_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_rei_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_rel_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_rel_fn_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_iciwp_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_iclwp_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_icswp_at_pc_full_level%f(nlev,ncol))

    allocate(pstate_cam%microp_qrout_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_nrout_at_pc_full_level%f(nlev,ncol))

    allocate(pstate_cam%microp_cc_t_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cc_qv_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cc_ql_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cc_qi_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cc_nl_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cc_ni_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cc_qlst_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cldo_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%microp_cldfsnow_at_pc_full_level%f(ptimelevels,nlev,ncol))

    allocate(pstate_cam%microp_lwc_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%microp_iwc_at_pc_full_level%f(nlev,ncol))

! initialize position for output
    pstate_cam%microp_wsedl_at_pc_full_level%pos   = 0
    pstate_cam%microp_dei_at_pc_full_level%pos     = 0
    pstate_cam%microp_des_at_pc_full_level%pos     = 0
    pstate_cam%microp_res_at_pc_full_level%pos     = 0
    pstate_cam%microp_mu_at_pc_full_level%pos      = 0
    pstate_cam%microp_lambdac_at_pc_full_level%pos = 0
    pstate_cam%microp_rei_at_pc_full_level%pos     = 0
    pstate_cam%microp_rel_at_pc_full_level%pos     = 0
    pstate_cam%microp_rel_fn_at_pc_full_level%pos  = 0 
    pstate_cam%microp_iciwp_at_pc_full_level%pos   = 0 
    pstate_cam%microp_iclwp_at_pc_full_level%pos   = 0 
    pstate_cam%microp_icswp_at_pc_full_level%pos   = 0 
    pstate_cam%microp_cc_t_at_pc_full_level%pos    = 0
    pstate_cam%microp_cc_qv_at_pc_full_level%pos   = 0
    pstate_cam%microp_cc_ql_at_pc_full_level%pos   = 0
    pstate_cam%microp_cc_qi_at_pc_full_level%pos   = 0
    pstate_cam%microp_cc_nl_at_pc_full_level%pos   = 0
    pstate_cam%microp_cc_ni_at_pc_full_level%pos   = 0
    pstate_cam%microp_cc_qlst_at_pc_full_level%pos = 0
    pstate_cam%microp_cldo_at_pc_full_level%pos    = 0
    pstate_cam%microp_cldfsnow_at_pc_full_level%pos= 0

    pstate_cam%microp_qrout_at_pc_full_level%pos   = 0
    pstate_cam%microp_nrout_at_pc_full_level%pos   = 0

    pstate_cam%microp_lwc_at_pc_full_level%pos   = 0
    pstate_cam%microp_iwc_at_pc_full_level%pos   = 0

! initialize to 0, if we have initial data, set to initial value.
    pstate_cam%microp_wsedl_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_dei_at_pc_full_level%f     = 0._r8
    pstate_cam%microp_des_at_pc_full_level%f     = 0._r8
    pstate_cam%microp_res_at_pc_full_level%f     = 0._r8
    pstate_cam%microp_mu_at_pc_full_level%f      = 0._r8
    pstate_cam%microp_lambdac_at_pc_full_level%f = 0._r8
    pstate_cam%microp_rei_at_pc_full_level%f     = 0._r8
    pstate_cam%microp_rel_at_pc_full_level%f     = 0._r8
    pstate_cam%microp_rel_fn_at_pc_full_level%f  = 0._r8
    pstate_cam%microp_iciwp_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_iclwp_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_icswp_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_cc_t_at_pc_full_level%f    = 0._r8
    pstate_cam%microp_cc_qv_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_cc_ql_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_cc_qi_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_cc_nl_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_cc_ni_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_cc_qlst_at_pc_full_level%f = 0._r8
    pstate_cam%microp_cldo_at_pc_full_level%f    = 0._r8
    pstate_cam%microp_cldfsnow_at_pc_full_level%f= 0._r8

    pstate_cam%microp_qrout_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_nrout_at_pc_full_level%f   = 0._r8

    pstate_cam%microp_lwc_at_pc_full_level%f   = 0._r8
    pstate_cam%microp_iwc_at_pc_full_level%f   = 0._r8

    end subroutine micro_mg_init_lin


    subroutine micro_mg_end_lin

    use grist_physics_data_structure,   only: pstate
    use grist_PhysC_data_structure,     only: pstate_cam
 
    deallocate(nevapr)
    deallocate(prain)
    deallocate(tnd_qsnow)
    deallocate(tnd_nsnow)
    deallocate(re_ice)
    deallocate(relvar)
    deallocate(accre_enhan)
    deallocate(iciwpst)
    deallocate(iclwpst)
    deallocate(acprecl)
    deallocate(acgcme)
    deallocate(acnum)

    if(allocated(pstate_cam%microp_wsedl_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_wsedl_at_pc_full_level%f)
    
    if(allocated(pstate_cam%microp_dei_at_pc_full_level%f))      &
      deallocate(pstate_cam%microp_dei_at_pc_full_level%f)
  
    if(allocated(pstate_cam%microp_des_at_pc_full_level%f))      &
      deallocate(pstate_cam%microp_des_at_pc_full_level%f)
  
    if(allocated(pstate_cam%microp_res_at_pc_full_level%f))      &
      deallocate(pstate_cam%microp_res_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_mu_at_pc_full_level%f))       &
      deallocate(pstate_cam%microp_mu_at_pc_full_level%f)
  
    if(allocated(pstate_cam%microp_lambdac_at_pc_full_level%f))  &
      deallocate(pstate_cam%microp_lambdac_at_pc_full_level%f)
  
    if(allocated(pstate_cam%microp_rei_at_pc_full_level%f))      &
      deallocate(pstate_cam%microp_rei_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_rel_at_pc_full_level%f))      &
      deallocate(pstate_cam%microp_rel_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_rel_fn_at_pc_full_level%f))   &
      deallocate(pstate_cam%microp_rel_fn_at_pc_full_level%f)
  
    if(allocated(pstate_cam%microp_iciwp_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_iciwp_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_iclwp_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_iclwp_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_icswp_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_icswp_at_pc_full_level%f)

    if(allocated(pstate_cam%microp_cc_t_at_pc_full_level%f))     &
      deallocate(pstate_cam%microp_cc_t_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cc_qv_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_cc_qv_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cc_ql_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_cc_ql_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cc_qi_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_cc_qi_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cc_nl_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_cc_nl_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cc_ni_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_cc_ni_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cc_qlst_at_pc_full_level%f))  &
      deallocate(pstate_cam%microp_cc_qlst_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cldo_at_pc_full_level%f))     & 
      deallocate(pstate_cam%microp_cldo_at_pc_full_level%f)
 
    if(allocated(pstate_cam%microp_cldfsnow_at_pc_full_level%f))     & 
      deallocate(pstate_cam%microp_cldfsnow_at_pc_full_level%f)

    if(allocated(pstate_cam%microp_qrout_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_qrout_at_pc_full_level%f)

    if(allocated(pstate_cam%microp_nrout_at_pc_full_level%f))    &
      deallocate(pstate_cam%microp_nrout_at_pc_full_level%f)

    if(allocated(pstate_cam%microp_lwc_at_pc_full_level%f))    &
       deallocate(pstate_cam%microp_lwc_at_pc_full_level%f)

   if(allocated(pstate_cam%microp_iwc_at_pc_full_level%f))    &
       deallocate(pstate_cam%microp_iwc_at_pc_full_level%f)

    end subroutine micro_mg_end_lin


    subroutine micro_mg_tend_lin(ncol,dtime)

    use phys_control,                  only: phys_getopts
    use cloud_fraction,                only: top_lev=>trop_cloud_top_lev
    use grist_physics_data_structure,  only: pstate, phy_tracer_info 
    use grist_PhysC_data_structure,    only: pstate_cam, &
                                             ptend_microphysics
    use grist_physics_update,          only: old_time_level,    &
                                             geopotential_dse
    use conv_water,                    only: conv_water_4rad 

! io
    integer ,   intent(in)    :: ncol
    real(r8),   intent(in)    :: dtime
! local
    integer  :: conv_water_in_rad
    logical  :: microp_uniform = .false. ! True = configure microphysics for sub-columns
                                         ! False = use in regular mode w/o sub-columns
    integer  :: i, k, m
    real(r8) :: alst_mic(nlev, ncol)
    real(r8) :: aist_mic(nlev, ncol)
    real(r8) :: rate1cld(nlev, ncol)     ! array to hold rate1ord_cw2pr_st from microphysics
    real(r8) :: tlat(nlev, ncol)
    real(r8) :: qvlat(nlev, ncol)
    real(r8) :: qcten(nlev, ncol)
    real(r8) :: qiten(nlev, ncol)
    real(r8) :: ncten(nlev, ncol)
    real(r8) :: niten(nlev, ncol)
    real(r8) :: effc(nlev, ncol)
    real(r8) :: effc_fn(nlev, ncol)  ! Liquid effective radius at fixed number (for indirect calc)
    real(r8) :: effi(nlev, ncol)
    real(r8) :: prect(ncol)
    real(r8) :: preci(ncol)

    real(r8) :: evapsnow(nlev, ncol)                    ! Local evaporation of snow
    real(r8) :: prodsnow(nlev, ncol)                    ! Local production of snow
    real(r8) :: cmeice(nlev, ncol)                      ! Rate of cond-evap of ice within the cloud
    real(r8) :: qsout(nlev, ncol)                       ! Snow mixing ratio
    real(r8) :: rflx(nlev+1, ncol)                      ! grid-box average rain flux (kg m^-2 s^-1)
    real(r8) :: sflx(nlev+1, ncol)                      ! grid-box average snow flux (kg m^-2 s^-1)
    real(r8) :: qrout(nlev, ncol)                       ! Rain mixing ratio
    real(r8) :: reff_rain(nlev, ncol)                   ! rain effective radius (um)
    real(r8) :: reff_snow(nlev, ncol)                   ! snow effective radius (um)
    real(r8) :: qcsevap(nlev, ncol)                     ! Evaporation of falling cloud water
    real(r8) :: qisevap(nlev, ncol)                     ! Sublimation of falling cloud ice
    real(r8) :: qvres(nlev, ncol)                       ! Residual condensation term to remove excess saturation
    real(r8) :: cmeiout(nlev, ncol)                     ! Deposition/sublimation rate of cloud ice
    real(r8) :: vtrmc(nlev, ncol)                       ! Mass-weighted cloud water fallspeed
    real(r8) :: vtrmi(nlev, ncol)                       ! Mass-weighted cloud ice fallspeed
    real(r8) :: qcsedten(nlev, ncol)                    ! Cloud water mixing ratio tendency from sedimentation
    real(r8) :: qisedten(nlev, ncol)                    ! Cloud ice mixing ratio tendency from sedimentation
    real(r8) :: prao(nlev, ncol)
    real(r8) :: prco(nlev, ncol)
    real(r8) :: mnuccco(nlev, ncol)
    real(r8) :: mnuccto(nlev, ncol)
    real(r8) :: msacwio(nlev, ncol)
    real(r8) :: psacwso(nlev, ncol)
    real(r8) :: bergso(nlev, ncol)
    real(r8) :: bergo(nlev, ncol)
    real(r8) :: melto(nlev, ncol)
    real(r8) :: homoo(nlev, ncol)
    real(r8) :: qcreso(nlev, ncol)
    real(r8) :: prcio(nlev, ncol)
    real(r8) :: praio(nlev, ncol)
    real(r8) :: qireso(nlev, ncol)
    real(r8) :: mnuccro(nlev, ncol)
    real(r8) :: pracso (nlev, ncol)
    real(r8) :: meltsdt(nlev, ncol)
    real(r8) :: frzrdt (nlev, ncol)
    real(r8) :: mnuccdo(nlev, ncol)
    real(r8) :: nrout(nlev, ncol)
    real(r8) :: nsout(nlev, ncol)
    real(r8) :: refl(nlev, ncol)                        ! analytic radar reflectivity
    real(r8) :: arefl(nlev, ncol)                       ! average reflectivity will zero points outside valid range
    real(r8) :: areflz(nlev, ncol)                      ! average reflectivity in z
    real(r8) :: frefl(nlev, ncol)
    real(r8) :: csrfl(nlev, ncol)                       ! cloudsat reflectivity
    real(r8) :: acsrfl(nlev, ncol)                      ! cloudsat average
    real(r8) :: fcsrfl(nlev, ncol)
    real(r8) :: rercld(nlev, ncol)                      ! effective radius calculation for rain + cloud
    real(r8) :: ncai(nlev, ncol)                        ! output number conc of ice nuclei available (1/m3)
    real(r8) :: ncal(nlev, ncol)                        ! output number conc of CCN (1/m3)
    real(r8) :: qrout2(nlev, ncol)
    real(r8) :: qsout2(nlev, ncol)
    real(r8) :: nrout2(nlev, ncol)
    real(r8) :: nsout2(nlev, ncol)
    real(r8) :: drout2(nlev, ncol)                      ! mean rain particle diameter (m)
    real(r8) :: dsout2(nlev, ncol)                      ! mean snow particle diameter (m)
    real(r8) :: freqs(nlev, ncol)
    real(r8) :: freqr(nlev, ncol)
    real(r8) :: nfice(nlev, ncol)

    real(r8) :: mnuccdohet(nlev, ncol)

! for COSP simulator
    real(r8) :: mgflxprc(nlevp, ncol)     ! 'LS_FLXPRC' MG grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
    real(r8) :: mgflxsnw(nlevp, ncol)     ! 'LS_FLXSNW' MG grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
    real(r8) :: mgmrprc(nlev, ncol)       ! 'LS_MRPRC' MG grid-box mean mixingratio_large_scale_cloud_rain+snow at interfaces (kg/kg)
    real(r8) :: mgmrsnw(nlev, ncol)       ! 'LS_MRSNW' MG grid-box mean mixingratio_large_scale_cloud_snow at interfaces (kg/kg)
    real(r8) :: mgreffrain(nlev, ncol)    ! 'LS_REFFRAIN' MG diagnostic rain effective radius (um)
    real(r8) :: mgreffsnow(nlev, ncol)    ! 'LS_REFFSNOW' MG diagnostic snow effective radius (um)
    real(r8) :: cvreffliq(nlev, ncol)     ! 'CV_REFFLIQ' convective cloud liquid effective radius (um)
    real(r8) :: cvreffice(nlev, ncol)     ! 'CV_REFFICE' convective cloud ice effective radius (um)

   ! For rrtm optics. specificed distribution.
    real(r8) :: mucon                      ! Convective size distribution shape parameter
    real(r8) :: dcon                       ! Convective size distribution effective radius (meters)  
    real(r8) :: deicon                     ! Convective ice effective diameter (meters)

    real(r8) :: qme(nlev, ncol)            ! 'QME'
  
    real(r8) :: icecldf(nlev, ncol)        ! Ice cloud fraction
    real(r8) :: liqcldf(nlev, ncol)        ! Liquid cloud fraction (combined into cloud)

    ! in-cloud water quantities adjusted for convective water
    real(r8) :: allcld_ice(nlev, ncol)     ! All-cloud cloud ice
    real(r8) :: allcld_liq(nlev, ncol)     ! All-cloud liquid

    real(r8) :: icimr(nlev, ncol)          ! In cloud ice mixing ratio
    real(r8) :: icwmr(nlev, ncol)          ! In cloud water mixing ratio
    real(r8) :: icimrst(nlev, ncol)        ! In stratus ice mixing ratio
    real(r8) :: icwmrst(nlev, ncol)        ! In stratus water mixing ratio
    real(r8) :: icinc(nlev, ncol)          ! In cloud ice number conc
    real(r8) :: icwnc(nlev, ncol)          ! In cloud water number conc
    real(r8) :: iwc(nlev, ncol)            ! Grid box average ice water content
    real(r8) :: lwc(nlev, ncol)            ! Grid box average liquid water content  
    real(r8) :: effliq(nlev, ncol)         ! In cloud liq eff rad
    real(r8) :: effice(nlev, ncol)         ! In cloud ice eff rad
    real(r8) :: effliq_fn(nlev, ncol)      ! In cloud liq eff rad at fixed number concentration	

    real(r8) :: cdnumc(ncol)               ! Vertically-integrated droplet concentration

    ! Averaging arrays for effective radius and number....
    real(r8) :: efiout(nlev, ncol)
    real(r8) :: efcout(nlev, ncol)
    real(r8) :: ncout(nlev, ncol)
    real(r8) :: niout(nlev, ncol)
    real(r8) :: freqi(nlev, ncol)
    real(r8) :: freql(nlev, ncol)

    real(r8) :: icecldf_out(nlev, ncol)                 ! Ice cloud fraction
    real(r8) :: liqcldf_out(nlev, ncol)                 ! Liquid cloud fraction (combined into cloud)
    real(r8) :: icimrst_out(nlev, ncol)                 ! In stratus ice mixing ratio
    real(r8) :: icwmrst_out(nlev, ncol)                 ! In stratus water mixing ratio

    ! Average cloud top radius & number
    real(r8) :: ctrel(ncol)
    real(r8) :: ctrei(ncol)
    real(r8) :: ctnl(ncol)
    real(r8) :: ctni(ncol)
    real(r8) :: fcti(ncol)
    real(r8) :: fctl(ncol)

    real(r8) :: ftem(nlev, ncol)

    ! Variables for precip efficiency calculation
    real(r8) :: minlwp                     ! LWP threshold
 
    ! Variables for liquid water path and column condensation
    real(r8) :: tgliqwp(ncol)              ! column liquid
    real(r8) :: tgcmeliq(ncol)             ! column condensation rate (units)

    real(r8) :: pe(ncol)                   ! precip efficiency for output
    real(r8) :: pefrac(ncol)               ! fraction of time precip efficiency is written out
    real(r8) :: tpr(ncol)                  ! average accumulated precipitation rate in pe calculation

    ! variables for autoconversion and accretion vertical averages
    real(r8) :: vprco(ncol)     ! vertical average autoconversion
    real(r8) :: vprao(ncol)     ! vertical average accretion
    real(r8) :: racau(ncol)     ! ratio of vertical averages
    integer  :: cnt(ncol)       ! counters
 
    ! local state
    real(r8) :: state_s(nlev, ncol)
    real(r8) :: state_t(nlev, ncol)
    real(r8) :: state_zm(nlev, ncol),   state_zi(nlev+1, ncol)
    real(r8) :: state_pmid(nlev, ncol), state_pint(nlev+1, ncol)
    real(r8) :: state_pdel(nlev, ncol)
    real(r8) :: state_phis(ncol)
    real(r8) :: state_q(ntracer, nlev, ncol)

    !----------Lin micro------------
    real(r8) :: riming(nlev, ncol)
    real(r8) :: swc(nlev, ncol)
    swc=0._r8
    !----------Lin micro------------


    call phys_getopts(conv_water_in_rad_out=conv_water_in_rad)

    alst_mic = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:,:ncol)
    aist_mic = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:,:ncol)

    select case (micro_mg_version)
    case (1)
       select case (micro_mg_sub_version)
       case (5)
          call micro_mg1_5_tend_intr_lin
       end select
    end select

    !-----------LiXH add for output--------------
    pstate_cam%microp_qrout_at_pc_full_level%f = qrout
    pstate_cam%microp_nrout_at_pc_full_level%f = nrout
    !-----------LiXH add for output--------------

    mnuccdohet = 0._r8
    do k=top_lev,nlev
       do i=1,ncol
          if (pstate_cam%microp_areo_naai%f(k,i) > 0._r8) then
             mnuccdohet(k,i) = mnuccdo(k,i) - (pstate_cam%microp_areo_naai_hom%f(k,i)/pstate_cam%microp_areo_naai%f(k,i))*mnuccdo(k,i)
          end if
       end do
    end do

    mgflxprc(top_lev:nlevp,:ncol) = rflx(top_lev:nlevp,:ncol) + sflx(top_lev:nlevp,:ncol)
    mgflxsnw(top_lev:nlevp,:ncol) = sflx(top_lev:nlevp,:ncol)

    mgmrprc(top_lev:nlev,:ncol) = qrout(top_lev:nlev,:ncol) + qsout(top_lev:nlev,:ncol)
    mgmrsnw(top_lev:nlev,:ncol) = qsout(top_lev:nlev,:ncol)

    mgreffrain(top_lev:nlev,:ncol) = reff_rain(top_lev:nlev,:ncol)
    mgreffsnow(top_lev:nlev,:ncol) = reff_snow(top_lev:nlev,:ncol)

    !! calculate effective radius of convective liquid and ice using dcon and deicon (not used by code, not useful for COSP)
    !! hard-coded as average of hard-coded values used for deep/shallow convective detrainment (near line 1502/1505)
    cvreffliq(top_lev:nlev,:ncol) = 9.0_r8
    cvreffice(top_lev:nlev,:ncol) = 37.0_r8

    ! output: 'LS_REFFRAIN'mgreffrain  'LS_REFFSNOW'mgreffsnow  'CV_REFFLIQ'cvreffliq  'CV_REFFICE'cvreffice

!-------------LiXH has not completed aero model---------------->
   ! Reassign rate1 if modal aerosols
!   if (rate1_cw2pr_st_idx > 0) then
!      call pbuf_get_field(pbuf, rate1_cw2pr_st_idx, rate1ord_cw2pr_st)
!      rate1ord_cw2pr_st(top_lev:nlev,:ncol) = rate1cld(top_lev:nlev,:ncol)
!   end if
!<------------LiXH has not completed aero model-----------------

    ! Sedimentation velocity for liquid stratus cloud droplet
    pstate_cam%microp_wsedl_at_pc_full_level%f(top_lev:nlev,:ncol) = vtrmc(top_lev:nlev,:ncol)

    ! Assign default size distribution parameters for no-stratiform clouds (convection only)
    ! Also put into physics buffer for possible separate use by radiation
    dcon   = 25.e-6_r8
    mucon  = 5.3_r8
    deicon = 50._r8

    do k = top_lev, nlev
       do i = 1, ncol 
          ! Convert snow effective diameter to microns
          pstate_cam%microp_des_at_pc_full_level%f(k,i) = pstate_cam%microp_des_at_pc_full_level%f(k,i) * 1.e6_r8
          pstate_cam%microp_res_at_pc_full_level%f(k,i) = pstate_cam%microp_res_at_pc_full_level%f(k,i) * 1.e6_r8   ! Linhan for RRTMG 4DDA
          if ( pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,k,i) < 1.e-4_r8 ) then
             pstate_cam%microp_mu_at_pc_full_level%f(k,i) = mucon
             pstate_cam%microp_lambdac_at_pc_full_level%f(k,i) = (mucon + 1._r8)/dcon
             pstate_cam%microp_dei_at_pc_full_level%f(k,i) = deicon
          end if
       end do
    end do

    ! Microphysical tendencies for use in the macrophysics at the next time step
    pstate_cam%microp_cc_t_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)    = tlat(top_lev:nlev,:ncol)/cp
    pstate_cam%microp_cc_qv_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)   = qvlat(top_lev:nlev,:ncol)
    pstate_cam%microp_cc_ql_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)   = qcten(top_lev:nlev,:ncol)
    pstate_cam%microp_cc_qi_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)   = qiten(top_lev:nlev,:ncol)
    pstate_cam%microp_cc_nl_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)   = ncten(top_lev:nlev,:ncol)
    pstate_cam%microp_cc_ni_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)   = niten(top_lev:nlev,:ncol)
    pstate_cam%microp_cc_qlst_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol) = qcten(top_lev:nlev,:ncol)/max(0.01_r8,alst_mic(top_lev:nlev,:ncol))

    ! Net micro_mg_cam condensation rate
    qme(top_lev:nlev,:ncol) = pstate_cam%macrop_cmeliq_at_pc_full_level%f(top_lev:nlev,:ncol) + cmeiout(top_lev:nlev,:ncol) 

    ! For precip, accumulate only total precip in prec_pwc and snow_pwc variables.
    ! Other precip output varirables are set to 0
    pstate_cam%pcw_prec_surface%f(:ncol) = prect(:ncol)
    pstate_cam%pcw_snow_surface%f(:ncol) = preci(:ncol)
    pstate_cam%sed_prec_surface%f(:ncol) = 0._r8
    pstate_cam%sed_snow_surface%f(:ncol) = 0._r8
    pstate_cam%str_prec_surface%f(:ncol) = pstate_cam%pcw_prec_surface%f(:ncol) + pstate_cam%sed_prec_surface%f(:ncol)
    pstate_cam%str_snow_surface%f(:ncol) = pstate_cam%pcw_snow_surface%f(:ncol) + pstate_cam%sed_snow_surface%f(:ncol)

    ! re-initialize ptend
    ptend_microphysics%tend_s%f = 0._r8
    ptend_microphysics%tend_q%f = 0._r8
 
    ptend_microphysics%tend_s%f(top_lev:nlev,:ncol)          = tlat(top_lev:nlev,:ncol)
    ptend_microphysics%tend_q%f(1,top_lev:nlev,:ncol)        = qvlat(top_lev:nlev,:ncol)
    ptend_microphysics%tend_q%f(ixcldliq,top_lev:nlev,:ncol) = qcten(top_lev:nlev,:ncol)
    ptend_microphysics%tend_q%f(ixcldice,top_lev:nlev,:ncol) = qiten(top_lev:nlev,:ncol)
    ptend_microphysics%tend_q%f(ixnumliq,top_lev:nlev,:ncol) = ncten(top_lev:nlev,:ncol)
    ptend_microphysics%tend_q%f(ixnumice,top_lev:nlev,:ncol) = niten(top_lev:nlev,:ncol)
 
    ! Check to make sure that the microphysics code is respecting the flags that control
    ! whether MG should be prognosing cloud ice and cloud liquid or not.
    if (.not. do_cldice) then
       if (any(qiten(top_lev:nlev,:ncol) /= 0.0_r8)) &
            call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
            " but micro_mg_tend has ice mass tendencies.")
       if (any(niten(top_lev:nlev,:ncol) /= 0.0_r8)) &
            call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud ice,"// &
            " but micro_mg_tend has ice number tendencies.")
    end if
    if (.not. do_cldliq) then
       if (any(qcten(top_lev:nlev,:ncol) /= 0.0_r8)) &
            call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
            " but micro_mg_tend has liquid mass tendencies.")
       if (any(ncten(top_lev:nlev,:ncol) /= 0.0_r8)) &
            call endrun("micro_mg_cam:ERROR - MG microphysics is configured not to prognose cloud liquid,"// &
            " but micro_mg_tend has liquid number tendencies.")
    end if

    ! Make copy of state and update it with local tendencies
    ! Need a copy of ptend to use in physics_update since that call will
    state_t(:,1:ncol)    = pstate%temp_at_pc_full_level%f(:,1:ncol)
    state_s(:,1:ncol)    = pstate%static_energy_at_pc_full_level%f(:,1:ncol)
    state_q(:,:,1:ncol)  = pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)
    state_zm(:,1:ncol)   = pstate%z_at_pc_full_level%f(:,1:ncol)
    state_zi(:,1:ncol)   = pstate%z_at_pc_face_level%f(:,1:ncol)
    state_pmid(:,1:ncol) = pstate%pressure_at_pc_full_level%f(:,1:ncol)
    state_pint(:,1:ncol) = pstate%pressure_at_pc_face_level%f(:,1:ncol)
    state_pdel(:,1:ncol) = pstate%delp_at_pc_full_level%f(:,1:ncol)
    state_phis(1:ncol)   = pstate%geop_at_pc_surface%f(1:ncol)

    state_q(:,:,1:ncol)  = state_q(:,:,1:ncol) + ptend_microphysics%tend_q%f(:,:,1:ncol)*dtime
    state_s(:,1:ncol)    = state_s(:,1:ncol)   + ptend_microphysics%tend_s%f(:,1:ncol)*dtime

! check negetive qv
    do m = 2, 5
       if(m .eq. ixnumliq .or. m .eq. ixnumice)then
            do k = 1, nlev
               do i = 1, ncol
                  state_q(m,k,i) = max(1.e-12_r8,state_q(m,k,i))
                  state_q(m,k,i) = min(1.e10_r8,state_q(m,k,i))
                end do
            end do
        else
            call qneg3('micro physics', ncol, nlev, m, m, phy_tracer_info(m)%qmin, state_q(m,:,1:ncol))
            !where(state_q(m,:,1:ncol) .lt. phy_tracer_info(m)%qmin) state_q(m,:,1:ncol) = phy_tracer_info(m)%qmin
        end if
 
    end do

    call geopotential_dse(ncol,       state_pint, state_pmid, state_pdel,   &
                          state_phis, state_s,    state_q(1,:,1:ncol),      &
                          state_t,    state_zm,   state_zi)

    icecldf(top_lev:nlev,:ncol) = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)
    liqcldf(top_lev:nlev,:ncol) = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)

    ! Effective droplet radius
    pstate_cam%microp_rel_at_pc_full_level%f(top_lev:nlev,:ncol)     = effc(top_lev:nlev,:ncol)
    pstate_cam%microp_rel_fn_at_pc_full_level%f(top_lev:nlev,:ncol)  = effc_fn(top_lev:nlev,:ncol)
    pstate_cam%microp_rei_at_pc_full_level%f(top_lev:nlev,:ncol)     = effi(top_lev:nlev,:ncol)
 
    ! ----------------------------------------------------------- ! 
    ! Adjust in-cloud water values to take account of convective  !
    ! in-cloud water. It is used to calculate the values of       !
    ! iclwp and iciwp to pass to the radiation.                   ! 
    ! ----------------------------------------------------------- !
    if( conv_water_in_rad /= 0 ) then
        allcld_ice(:,:ncol) = 0._r8 ! Grid-avg all cloud liquid
        allcld_liq(:,:ncol) = 0._r8 ! Grid-avg all cloud ice
        call conv_water_4rad( ncol, conv_water_in_rad,                                       &
                              pstate_cam%microp_rei_at_pc_full_level%f(:,:ncol),                 &
                              state_pdel, state_q(ixcldliq,:,:), state_q(ixcldice,:,:),      &
                              allcld_liq, allcld_ice )
    else
        allcld_liq(top_lev:nlev,:ncol) = state_q(ixcldliq,top_lev:nlev,:ncol)  ! Grid-ave all cloud liquid
        allcld_ice(top_lev:nlev,:ncol) = state_q(ixcldice,top_lev:nlev,:ncol)  !           "        ice 
    end if

    ! ------------------------------------------------------------ !
    ! Compute in cloud ice and liquid mixing ratios                !
    ! Note that 'iclwp, iciwp' are used for radiation computation. !
    ! ------------------------------------------------------------ !

    icimr = 0._r8
    icwmr = 0._r8
    icimrst_out = 0._r8
    icwmrst_out = 0._r8
    icinc = 0._r8
    icwnc = 0._r8
    iwc = 0._r8
    lwc = 0._r8
    pstate_cam%microp_iciwp_at_pc_full_level%f = 0._r8
    pstate_cam%microp_iclwp_at_pc_full_level%f = 0._r8
    pstate_cam%microp_icswp_at_pc_full_level%f = 0._r8
    iciwpst = 0._r8
    iclwpst = 0._r8
    effliq = 0._r8
    effliq_fn = 0._r8
    effice = 0._r8
    pstate_cam%microp_cldfsnow_at_pc_full_level%f = 0._r8

    do k = top_lev, nlev
       do i = 1, ncol
          ! Limits for in-cloud mixing ratios consistent with MG microphysics
          ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
          icimr(k,i)     = min( allcld_ice(k,i) / max(0.0001_r8,pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i)),0.005_r8 )
          icwmr(k,i)     = min( allcld_liq(k,i) / max(0.0001_r8,pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i)),0.005_r8 )
          icimrst(k,i)   = min( state_q(ixcldice,k,i) / max(0.0001_r8,icecldf(k,i)),0.005_r8 )
          icwmrst(k,i)   = min( state_q(ixcldliq,k,i) / max(0.0001_r8,liqcldf(k,i)),0.005_r8 )
          icinc(k,i)     = state_q(ixnumice,k,i) / max(0.0001_r8,icecldf(k,i)) * &
                           state_pmid(k,i) / (287.15_r8*state_t(k,i))
          icwnc(k,i)     = state_q(ixnumliq,k,i) / max(0.0001_r8,liqcldf(k,i)) * &
                           state_pmid(k,i) / (287.15_r8*state_t(k,i))
          iwc(k,i)       = allcld_ice(k,i) * state_pmid(k,i) / (287.15_r8*state_t(k,i))
          lwc(k,i)       = allcld_liq(k,i) * state_pmid(k,i) / (287.15_r8*state_t(k,i))
          effliq(k,i)    = effc(k,i)
          effliq_fn(k,i) = effc_fn(k,i)
          effice(k,i)    = effi(k,i)

          ! grid-mean lwc and iwc for IO:
          pstate_cam%microp_lwc_at_pc_full_level%f(k,i) = allcld_liq(k,i) * state_pmid(k,i) / (287.15_r8*state_t(k,i))
          pstate_cam%microp_iwc_at_pc_full_level%f(k,i) = allcld_ice(k,i) * state_pmid(k,i) / (287.15_r8*state_t(k,i))

          ! Calculate total cloud water paths in each layer
          pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) = icimr(k,i) * state_pdel(k,i) / gravity
          pstate_cam%microp_iclwp_at_pc_full_level%f(k,i) = icwmr(k,i) * state_pdel(k,i) / gravity
          ! Calculate micro_mg_cam cloud water paths in each layer
          ! Note: uses stratiform cloud fraction!
          iciwpst(k,i)   = min(state_q(ixcldice,k,i)/                                                           &
                           max(0.0001_r8,pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,k,i)),0.005_r8)    &
                           * state_pdel(k,i) / gravity
          iclwpst(k,i)   = min(state_q(ixcldliq,k,i)/                                                           &
                           max(0.0001_r8,pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,k,i)),0.005_r8)    &
                           * state_pdel(k,i) / gravity

          ! ------------------------------ !
          ! Adjust cloud fraction for snow !
          ! ------------------------------ !
          pstate_cam%microp_cldfsnow_at_pc_full_level%f(old_time_level,k,i) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i)
          ! If cloud and only ice ( no convective cloud or ice ), then set to 0.
          if( ( pstate_cam%microp_cldfsnow_at_pc_full_level%f(old_time_level,k,i) .gt. 1.e-4_r8 ) .and. & 
              ( pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,k,i)   .lt. 1.e-4_r8 ) .and. & 
              ( state_q(ixcldliq,k,i) .lt. 1.e-10_r8 ) ) then
              pstate_cam%microp_cldfsnow_at_pc_full_level%f(old_time_level,k,i) = 0._r8
          end if
          ! If no cloud and snow, then set to 0.25
          if( ( pstate_cam%microp_cldfsnow_at_pc_full_level%f(old_time_level,k,i) .lt. 1.e-4_r8 ) .and. &
              ( qsout(k,i) .gt. 1.e-6_r8 ) ) then 
              pstate_cam%microp_cldfsnow_at_pc_full_level%f(old_time_level,k,i) = 0.25_r8
          end if
          ! Calculate in-cloud snow water path
          pstate_cam%microp_icswp_at_pc_full_level%f(k,i) = qsout(k,i)                                         & 
                     / max( 0.0001_r8, pstate_cam%microp_cldfsnow_at_pc_full_level%f(old_time_level,k,i) )     &
                     * state_pdel(k,i) / gravity
          
          !----------Lin micro------------
          swc(k,i) = qsout(k,i)*state_pdel(k,i)/gravity
          !----------Lin micro------------
       end do
    end do

    ! ------------------------------------- !
    ! Precipitation efficiency Calculation  !
    ! ------------------------------------- !

    ! Liquid water path

    ! Compute liquid water paths, and column condensation
    tgliqwp(:ncol) = 0._r8
    tgcmeliq(:ncol) = 0._r8
    do k = top_lev, nlev
       do i = 1, ncol
          tgliqwp(i)  = tgliqwp(i)                                                                          &
                      + pstate_cam%microp_iclwp_at_pc_full_level%f(k,i)                                         &
                      * pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i)

          if (pstate_cam%macrop_cmeliq_at_pc_full_level%f(k,i) > 1.e-12_r8) then
             !convert cmeliq to right units:  kgh2o/kgair/s  *  kgair/m2  / kgh2o/m3  = m/s
             tgcmeliq(i) = tgcmeliq(i)                                                                      &
                         + pstate_cam%macrop_cmeliq_at_pc_full_level%f(k,i)                                     &
                         * (state_pdel(k,i) / gravity) / rhoh2o
          end if
       end do
    end do

    ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
    ! this is 1ppmv of h2o in 10hpa
    ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9

    !-----------------------------------------------------------------------
    ! precipitation efficiency calculation  (accumulate cme and precip)

    minlwp = 0.01_r8        !minimum lwp threshold (kg/m3)

    ! zero out precip efficiency and total averaged precip
    pe(:ncol)     = 0._r8
    tpr(:ncol)    = 0._r8
    pefrac(:ncol) = 0._r8

    ! accumulate precip and condensation
    do i = 1, ncol

       acgcme(i)  = acgcme(i) + tgcmeliq(i)
       acprecl(i) = acprecl(i) + pstate_cam%str_prec_surface%f(i)
       acnum(i)   = acnum(i) + 1

       ! if LWP is zero, then 'end of cloud': calculate precip efficiency
       if (tgliqwp(i) < minlwp) then
          if (acprecl(i) > 5.e-8_r8) then
             tpr(i) = max(acprecl(i)/acnum(i), 1.e-15_r8)
             if (acgcme(i) > 1.e-10_r8) then
                pe(i) = min(max(acprecl(i)/acgcme(i), 1.e-15_r8), 1.e5_r8)
                pefrac(i) = 1._r8
             end if
          end if

          ! reset counters
!         if (pe(i) /= 0._r8 .and. (pe(i) < 1.e-8_r8 .or. pe(i) > 1.e3_r8)) then
!            write (iulog,*) 'PE:ANOMALY  pe, acprecl, acgcme, tpr, acnum ',pe(i),acprecl(i), acgcme(i), tpr(i), acnum(i)
!         endif

          acprecl(i) = 0._r8
          acgcme(i)  = 0._r8
          acnum(i)   = 0
       end if               ! end LWP zero conditional

       ! if never find any rain....(after 10^3 timesteps...)
       if (acnum(i) > 1000) then
          acnum(i)   = 0
          acprecl(i) = 0._r8
          acgcme(i)  = 0._r8
       end if

    end do

    ! output: 'PE'pe  'PEFRAC'pefrac  'APRL'tpr

    !-----------------------------------------------------------------------
    ! vertical average of non-zero accretion, autoconversion and ratio.
    ! vars: vprco(i),vprao(i),racau(i),cnt

    vprao = 0._r8
    cnt = 0
    do k = top_lev, nlev
       vprao(:ncol) = vprao(:ncol) + prao(k,:ncol)
       where (prao(k,:ncol) /= 0._r8) cnt(:ncol) = cnt(:ncol) + 1
    end do

    where (cnt > 0) vprao = vprao/cnt

    vprco = 0._r8
    cnt = 0
    do k = top_lev, nlev
       vprco(:ncol) = vprco(:ncol) + prco(k,:ncol)
       where (prco(k,:ncol) /= 0._r8) cnt(:ncol) = cnt(:ncol) + 1
    end do

    where (cnt > 0)
       vprco = vprco/cnt
       racau = vprao/vprco
    elsewhere
       racau = 0._r8
    end where

    ! output: 'VPRAO'vprao  'VPRCO'vprco  'RACAU'racau

    ! --------------------- !
    ! History Output Fields !
    ! --------------------- !

    ! Column droplet concentration
    !------------------LiXH Modified------------------>
    !cdnumc(:ncol) = sum(state_q(ixnumliq,top_lev:nlev,:ncol) * state_pdel(top_lev:nlev,:ncol)/gravity, dim=2)
    do i = 1, ncol
        cdnumc(i) = sum(state_q(ixnumliq,top_lev:nlev,i)*state_pdel(top_lev:nlev,i)/gravity)
    end do
    !<-----------------LiXH Modified-------------------

    ! Averaging for new output fields
    efcout      = 0._r8
    efiout      = 0._r8
    ncout       = 0._r8
    niout       = 0._r8
    freql       = 0._r8
    freqi       = 0._r8
    liqcldf_out = 0._r8
    icecldf_out = 0._r8
    icwmrst_out = 0._r8
    icimrst_out = 0._r8

    do k = top_lev, nlev
       do i = 1, ncol
          if ( liqcldf(k,i) > 0.01_r8 .and. icwmrst(k,i) > 5.e-5_r8 ) then
             efcout(k,i) = effc(k,i) * liqcldf(k,i)
             ncout(k,i)  = icwnc(k,i) * liqcldf(k,i)
             freql(k,i)  = liqcldf(k,i)
             liqcldf_out(k,i) = liqcldf(k,i)
             icwmrst_out(k,i) = icwmrst(k,i)
          end if
          if ( icecldf(k,i) > 0.01_r8 .and. icimrst(k,i) > 1.e-6_r8 ) then
             efiout(k,i) = effi(k,i) * icecldf(k,i)
             niout(k,i)  = icinc(k,i) * icecldf(k,i)
             freqi(k,i)  = icecldf(k,i)
             icecldf_out(k,i) = icecldf(k,i)
             icimrst_out(k,i) = icimrst(k,i)
          end if
       end do
    end do

    ! output: 'AREL'efcout  'AREI'efiout  'AWNC'ncout  'AWNI'niout  'FREQL'freql  'FREQI'freqi 

    ! Cloud top effective radius and number.
    fcti  = 0._r8
    fctl  = 0._r8
    ctrel = 0._r8
    ctrei = 0._r8
    ctnl  = 0._r8
    ctni  = 0._r8
    do i = 1, ncol
       do k = top_lev, nlev
          if ( liqcldf(k,i) > 0.01_r8 .and. icwmrst(k,i) > 1.e-7_r8 ) then
             ctrel(i) = effc(k,i) * liqcldf(k,i)
             ctnl(i)  = icwnc(k,i) * liqcldf(k,i)
             fctl(i)  = liqcldf(k,i)
             exit
          end if
          if ( icecldf(k,i) > 0.01_r8 .and. icimrst(k,i) > 1.e-7_r8 ) then
             ctrei(i) = effi(k,i) * icecldf(k,i)
             ctni(i)  = icinc(k,i) * icecldf(k,i)
             fcti(i)  = icecldf(k,i)
             exit
          end if
       end do
    end do

! output: 'ACTREL'ctrel  'ACTREI'ctrei  'ACTNL'ctnl  'ACTNI'ctni  'FCTL'fctl  'FCTI'fcti
!         'MPDT'tlat  'MPDQ'qvlat  'MPDLIQ'qcten  'MPDICE'qiten  'ICINC'icinc  'ICWNC'icwnc
!         'EFFLIQ'effliq  'EFFLIQ_IND'effliq_fn  'EFFICE'effice  'CDNUMC'cdnumc
!         'LS_FLXPRC'mgflxprc  'LS_FLXSNW'mgflxsnw  'REL'rel  'REI'rei
!         'IWC'iwc  'LWC'lwc  'ICIMR'icimr  'ICWMR'icwmr  'ICIMRST'icimrst_out  'ICWMRST'icwmrst_out
!         'CME'qme  'PRODPREC'prain  'EVAPPREC'nevapr  'EVAPSNOW'evapsnow  'QCSEVAP'qcsevap  'QISEVAP'qisevap
!         'QVRES'qvres  'CMEIOUT'cmeiout  'VTRMC'vtrmc  'VTRMI'vtrmi  'QCSEDTEN'qcsedten  'QISEDTEN'qisedten
!         'PRAO'prao  'PRCO'prco  'MNUCCCO'mnuccco  'MNUCCTO'mnuccto  'MNUCCDO'mnuccdo  'MNUCCDOhet'mnuccdohet 
!         'MSACWIO'msacwio  'PSACWSO'psacwso  'BERGSO'bergso  'BERGO'bergo  'MELTO'melto  'HOMOO'homoo
!         'QCRESO'qcreso  'PRCIO'prcio  'PRAIO'praio  'QIRESO'qireso  'MNUCCRO'mnuccro  'PRACSO'pracso
!         'MELTSDT'meltsdt  'FRZRDT'frzrdt

    ftem = 0._r8

    ftem(top_lev:nlev,:ncol) =  qcreso(top_lev:nlev,:ncol)
! output: 'MPDW2V'ftem

    ftem(top_lev:nlev,:ncol) =  melto(top_lev:nlev,:ncol) - mnuccco(top_lev:nlev,:ncol) - mnuccto(top_lev:nlev,:ncol) - &
                         bergo(top_lev:nlev,:ncol) - homoo  (top_lev:nlev,:ncol) - msacwio(top_lev:nlev,:ncol)
! output: 'MPDW2I'ftem

    ftem(top_lev:nlev,:ncol) = -prao(top_lev:nlev,:ncol) - prco(top_lev:nlev,:ncol) - psacwso(top_lev:nlev,:ncol) - &
                         bergso(top_lev:nlev,:ncol)
! output: 'MPDW2P'ftem

    ftem(top_lev:nlev,:ncol) =  cmeiout(top_lev:nlev,:ncol) + qireso (top_lev:nlev,:ncol)
! output: 'MPDI2V'ftem
    ftem(top_lev:nlev,:ncol) = -melto(top_lev:nlev,:ncol) + mnuccco(top_lev:nlev,:ncol) + mnuccto(top_lev:nlev,:ncol) + &
                         bergo(top_lev:nlev,:ncol) + homoo  (top_lev:nlev,:ncol) + msacwio(top_lev:nlev,:ncol)
! output: 'MPDI2W'ftem

    ftem(top_lev:nlev,:ncol) = -prcio(top_lev:nlev,:ncol) - praio  (top_lev:nlev,:ncol)
! output: 'MPDI2P'ftem


contains

    subroutine micro_mg1_5_tend_intr_lin
    use micro_mg1_5_lin,        only: micro_mg1_5_tend_lin, micro_mg1_5_get_cols_lin

    integer :: nlev_loc
    integer :: mgncol
    integer, allocatable :: mgcols(:)

    ! Number of total levels over which to do the microphysics
    nlev_loc  = nlev - top_lev + 1

    ! save old cloud fraction here:
    ! This seems unnecessary here, since no input variables are overwritten by MG1.5
    pstate_cam%microp_cldo_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)=pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,top_lev:nlev,:ncol)

    call micro_mg1_5_get_cols_lin(ncol, nlev_loc, top_lev, pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,:,:ncol), &
                              pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,:,:ncol), mgncol, mgcols) 


    call micro_mg1_5_tend_lin( mgncol,   mgcols,   nlev_loc,   top_lev,   dtime,                      &
                           pstate%temp_at_pc_full_level%f(:,:ncol),                               &
                           pstate%tracer_mxrt_at_pc_full_level%f(1,:,:ncol),                      &
                           pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,:,:ncol),               &
                           pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,:,:ncol),               &
                           pstate%tracer_mxrt_at_pc_full_level%f(ixnumliq,:,:ncol),               &
                           pstate%tracer_mxrt_at_pc_full_level%f(ixnumice,:,:ncol),               &  
                           relvar,   accre_enhan,                                                 &
                           pstate%pressure_at_pc_full_level%f(:,1:ncol),                          &
                           pstate%delp_at_pc_full_level%f(:,1:ncol),                              &  
                           pstate%pressure_at_pc_face_level%f(:,1:ncol),                          &
                           pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:,:ncol),      &
                           alst_mic, aist_mic, rate1cld,                                          &
                           pstate_cam%microp_areo_naai%f(:,:ncol),                                &  
                           pstate_cam%microp_areo_npccn%f(:,:ncol),                               &  
                           pstate_cam%microp_areo_rndst%f(1:4,:,:ncol),                           &
                           pstate_cam%microp_areo_nacon%f(1:4,:,:ncol),                           &  
                           tlat,     qvlat,    qcten,      qiten,     ncten,      niten,          &
                           effc,     effc_fn,  effi,                  prect,      preci,          &  
                           nevapr,   evapsnow, prain,      prodsnow,  cmeice,                     &
                           pstate_cam%microp_dei_at_pc_full_level%f(:,:ncol),                     & 
                           pstate_cam%microp_mu_at_pc_full_level%f(:,:ncol),                      &
                           pstate_cam%microp_lambdac_at_pc_full_level%f(:,:ncol), qsout,          &   
                           pstate_cam%microp_des_at_pc_full_level%f(:,:ncol),     rflx,    sflx,  &
                           qrout,              reff_rain,             reff_snow,                  & 
                           qcsevap,  qisevap,  qvres,      cmeiout,   vtrmc,      vtrmi,          &
                           qcsedten, qisedten, prao,       prco,      mnuccco,    mnuccto,        &
                           msacwio,  psacwso,  bergso,     bergo,     melto,      homoo,          &
                           qcreso,             prcio,      praio,     qireso,                     &
                           mnuccro,  pracso,   meltsdt,    frzrdt,    mnuccdo,                    &
                           nrout,    nsout,    refl,       arefl,     areflz,     frefl,          &
                           csrfl,    acsrfl,   fcsrfl,                rercld,                     &
                           ncai,     ncal,     qrout2,     qsout2,    nrout2,     nsout2,         &
                           drout2,   dsout2,   freqs,      freqr,     nfice,                      &
                           tnd_qsnow,          tnd_nsnow,             re_ice,                     &
                           !----------Lin micro------------
                           riming)
                           !----------Lin micro------------
 

    end subroutine micro_mg1_5_tend_intr_lin
    
    end subroutine micro_mg_tend_lin

  end module micro_mg_lin
