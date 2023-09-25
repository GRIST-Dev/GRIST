!===================================================================================================
!
!  Created by LiXiaohan on 19/07/10, adopted from CAM5
!
!     interface of the prognostic cloud macrophysics (Park scheme)
!
!===================================================================================================

 module grist_macrop_park

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8, i4, latice, gravity,  &
                                                  cp, tmelt
    use grist_handle_error,                 only: endrun
    use grist_mpi

    implicit none
    public          :: macrop_driver_init_park,        &
                       read_nml_macrop_park,           &
                       macrop_driver_tend_park,        &
                       end_of_macrop_park

    private
    save

! Private:
    logical :: do_cldice             ! .true., park macrophysics is prognosing cldice
    logical :: do_cldliq             ! .true., park macrophysics is prognosing cldliq
    logical :: do_detrain            ! .true., park macrophysics is detraining ice into stratiform
    logical :: use_shfrc             ! Flag, whether use sh_cld_frc from shallow_conv (LiXH)

    logical, parameter :: cu_det_st  = .false.  

    integer :: pcnst, pver, pverp
    integer :: ixcldliq, ixcldice, ixnumliq, ixnumice

! Tuning parameters, LiXH:
    !Rdet: effective radius of detrained convective condensate
    real(r8) , parameter :: Rdet_dp_liq = 8.e-6_r8 
    real(r8) , parameter :: Rdet_dp_ice = 25.e-6_r8 
    real(r8) , parameter :: Rdet_sh_liq = 8.e-6_r8  !default: 10.e-6_r8 , DP: 8.e-6_r8
    real(r8) , parameter :: Rdet_sh_ice = 25.e-6_r8 !default: 50.e-6_r8 , DP: 25.e-6_r8

 contains

    subroutine read_nml_macrop_park(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer :: unitn, ierr
    logical :: macro_do_cldice  = .true.  ! park macrophysics is prognosing cldice
    logical :: macro_do_cldliq  = .true.  ! park macrophysics is prognosing cldliq
    logical :: macro_do_detrain = .true.  ! park macrophysics is detraining ice into stratiform

    namelist /macro_nl/ macro_do_cldice, macro_do_cldliq, macro_do_detrain

    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, macro_nl, iostat=ierr)
    if (ierr /= 0) call endrun(' error reading macro namelist ')
    close(unitn)
 
    do_cldice  = macro_do_cldice
    do_cldliq  = macro_do_cldliq
    do_detrain = macro_do_detrain

    end subroutine read_nml_macrop_park


    subroutine macrop_driver_init_park(ncol)
    use grist_nml_module,               only: nlev, nlevp, ntracer
    use grist_physics_data_structure,   only: pstate, phy_tracer_info
    use grist_PhysC_data_structure,      only: pstate_cam
    use grist_physics_update,           only: ptimelevels
    use phys_control,                   only: phys_getopts
! io
    integer, intent(in) :: ncol
! local
    integer :: m
    character(len=16) :: shallow_scheme

! Initialize pstate
    allocate(pstate_cam%macrop_cmeliq_at_pc_full_level%f(nlev,ncol))
    ! multi time levels:
    allocate(pstate_cam%macrop_cld_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%macrop_concld_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%macrop_ast_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%macrop_alst_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%macrop_aist_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%macrop_qlst_at_pc_full_level%f(ptimelevels,nlev,ncol))
    allocate(pstate_cam%macrop_qist_at_pc_full_level%f(ptimelevels,nlev,ncol))

    pstate_cam%macrop_ast_at_pc_full_level%pos   = 0
    pstate_cam%macrop_alst_at_pc_full_level%pos  = 0
    pstate_cam%macrop_aist_at_pc_full_level%pos  = 0
    pstate_cam%macrop_qlst_at_pc_full_level%pos  = 0
    pstate_cam%macrop_qist_at_pc_full_level%pos  = 0
    pstate_cam%macrop_cld_at_pc_full_level%pos   = 0
    pstate_cam%macrop_concld_at_pc_full_level%pos= 0
    pstate_cam%macrop_cmeliq_at_pc_full_level%pos= 0
 
    pstate_cam%macrop_ast_at_pc_full_level%f     = 0._r8
    pstate_cam%macrop_alst_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_aist_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_qlst_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_qist_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_cld_at_pc_full_level%f     = 0._r8
    pstate_cam%macrop_concld_at_pc_full_level%f  = 0._r8
    pstate_cam%macrop_cmeliq_at_pc_full_level%f  = 0._r8

    pcnst = ntracer
    pver  = nlev
    pverp = nlevp

    allocate(pstate_cam%macrop_qcwat_at_pc_full_level%f(ptimelevels, pver, ncol))
    allocate(pstate_cam%macrop_lcwat_at_pc_full_level%f(ptimelevels, pver, ncol))
    allocate(pstate_cam%macrop_iccwat_at_pc_full_level%f(ptimelevels, pver, ncol))
    allocate(pstate_cam%macrop_nlwat_at_pc_full_level%f(ptimelevels, pver, ncol))
    allocate(pstate_cam%macrop_niwat_at_pc_full_level%f(ptimelevels, pver, ncol))
    allocate(pstate_cam%macrop_tcwat_at_pc_full_level%f(ptimelevels, pver, ncol))

    pstate_cam%macrop_qcwat_at_pc_full_level%pos = 0
    pstate_cam%macrop_lcwat_at_pc_full_level%pos = 0
    pstate_cam%macrop_iccwat_at_pc_full_level%pos= 0
    pstate_cam%macrop_nlwat_at_pc_full_level%pos = 0
    pstate_cam%macrop_niwat_at_pc_full_level%pos = 0
    pstate_cam%macrop_tcwat_at_pc_full_level%pos = 0

    pstate_cam%macrop_qcwat_at_pc_full_level%f = 0
    pstate_cam%macrop_lcwat_at_pc_full_level%f = 0
    pstate_cam%macrop_iccwat_at_pc_full_level%f= 0
    pstate_cam%macrop_nlwat_at_pc_full_level%f = 0
    pstate_cam%macrop_niwat_at_pc_full_level%f = 0
    pstate_cam%macrop_tcwat_at_pc_full_level%f = 0
 
    ! inquire cloud index
    do m = 1, pcnst
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')        ixcldliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')           ixcldice = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid_number') ixnumliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice_number')    ixnumice = m 
    end do

    call phys_getopts(shallow_scheme_out = shallow_scheme)
    if ( shallow_scheme .eq. 'UW' .or. shallow_scheme .eq. 'double_plume') then
        use_shfrc = .true.
    else
        use_shfrc = .false.
    end if

    end subroutine macrop_driver_init_park


    subroutine end_of_macrop_park

    use grist_PhysC_data_structure,   only: pstate_cam

    if(allocated(pstate_cam%macrop_qcwat_at_pc_full_level%f))      &
    deallocate(pstate_cam%macrop_qcwat_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_lcwat_at_pc_full_level%f))      &
    deallocate(pstate_cam%macrop_lcwat_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_iccwat_at_pc_full_level%f))      &
    deallocate(pstate_cam%macrop_iccwat_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_nlwat_at_pc_full_level%f))      &
    deallocate(pstate_cam%macrop_nlwat_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_niwat_at_pc_full_level%f))      &
    deallocate(pstate_cam%macrop_niwat_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_tcwat_at_pc_full_level%f))      &
    deallocate(pstate_cam%macrop_tcwat_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_ast_at_pc_full_level%f))      &
      deallocate(pstate_cam%macrop_ast_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_alst_at_pc_full_level%f))     &
      deallocate(pstate_cam%macrop_alst_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_aist_at_pc_full_level%f))     &
      deallocate(pstate_cam%macrop_aist_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_qlst_at_pc_full_level%f))     &
      deallocate(pstate_cam%macrop_qlst_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_qist_at_pc_full_level%f))     &
      deallocate(pstate_cam%macrop_qist_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_cld_at_pc_full_level%f))      &
      deallocate(pstate_cam%macrop_cld_at_pc_full_level%f)

    if(allocated(pstate_cam%macrop_concld_at_pc_full_level%f))   &
      deallocate(pstate_cam%macrop_concld_at_pc_full_level%f)
#ifndef CMAPI
    if(allocated(pstate_cam%macrop_cmeliq_at_pc_full_level%f))   &
      deallocate(pstate_cam%macrop_cmeliq_at_pc_full_level%f)
#endif

    end subroutine end_of_macrop_park


    subroutine macrop_driver_tend_park(ncol  , dtime , istep   ,     &
                                  dlf   , dlf2  ,  cmfmc2 ,     &
                                  det_s , det_ice)
    use cloud_fraction,                 only: top_lev => trop_cloud_top_lev, &
                                              cldfrc, cldfrc_fice
    use cldwat2m_macro,                 only: mmacro_pcond
    use grist_physics_data_structure,   only: pstate, phy_tracer_info
    use grist_PhysC_data_structure,      only: pstate_cam,                    &
                                              ptend_macrophysics
    use grist_physics_update,           only: geopotential_dse,              &
                                              old_time_level,                &
                                              ptimelevels
    use grist_mpi

! io
    integer,  intent(in)  :: ncol
    integer,  intent(in)  :: istep
    real(r8), intent(in)  :: dtime                   ! Timestep
    real(r8), intent(in)  :: dlf(pver, ncol)         ! Detrained water from convection schemes
    real(r8), intent(in)  :: dlf2(pver, ncol)        ! Detrained water from shallow convection scheme
    real(r8), intent(in)  :: cmfmc2(pverp, ncol)     ! Shallow convective mass flux [ kg/s/m^2 ]

    ! These two variables are needed for energy check    
    real(r8), intent(out) :: det_s(ncol)             ! Integral of detrained static energy from ice
    real(r8), intent(out) :: det_ice(ncol)           ! Integral of detrained ice for energy check

! local
    integer :: i,k,m,itim

    real(r8), allocatable :: shfrc(:,:)              ! Cloud fraction from shallow convection scheme

    ! For cldfrc
    real(r8)  cldst(pver,ncol)                       ! Stratus cloud fraction
    real(r8)  rhcloud(pver,ncol)                     ! Relative humidity cloud (last timestep)
    real(r8)  clc(ncol)                              ! Column convective cloud amount
    real(r8)  rhu00(pver,ncol)                       ! RH threshold for cloud
    real(r8)  icecldf(pver,ncol)                     ! Ice cloud fraction
    real(r8)  liqcldf(pver,ncol)                     ! Liquid cloud fraction (combined into cloud)
    real(r8)  relhum(pver,ncol)                      ! RH, output to determine drh/da

    ! For macrophysics
    real(r8)  rdtime
    real(r8)  qtend(pver, ncol)                      ! Moisture tendencies
    real(r8)  ttend(pver, ncol)                      ! Temperature tendencies
    real(r8)  ltend(pver, ncol)                      ! Cloud liquid water tendencies
    real(r8)  fice(pver,ncol)
    real(r8)  fsnow(pver,ncol)
    real(r8)  dpdlfliq(pver,ncol)
    real(r8)  dpdlfice(pver,ncol)
    real(r8)  shdlfliq(pver,ncol)
    real(r8)  shdlfice(pver,ncol)
    real(r8)  dpdlft  (pver,ncol)
    real(r8)  shdlft  (pver,ncol)

    real(r8)  dum1

    ! Output from mmacro_pcond
    real(r8)  tlat(pver, ncol)
    real(r8)  qvlat(pver, ncol)
    real(r8)  qcten(pver, ncol)
    real(r8)  qiten(pver, ncol)
    real(r8)  ncten(pver, ncol)
    real(r8)  niten(pver, ncol)
  
    ! Output from mmacro_pcond
    real(r8)  qvadj(pver, ncol)                       ! Macro-physics adjustment tendency from "positive_moisture" call (vapor)
    real(r8)  qladj(pver, ncol)                       ! Macro-physics adjustment tendency from "positive_moisture" call (liquid)
    real(r8)  qiadj(pver, ncol)                       ! Macro-physics adjustment tendency from "positive_moisture" call (ice)
    real(r8)  qllim(pver, ncol)                       ! Macro-physics tendency from "instratus_condensate" call (liquid)
    real(r8)  qilim(pver, ncol)                       ! Macro-physics tendency from "instratus_condensate" call (ice)

    ! For revised macophysics, mmacro_pcond
    real(r8)  itend(pver, ncol)
    real(r8)  lmitend(pver, ncol)
    real(r8)  t_inout(pver, ncol)
    real(r8)  qv_inout(pver, ncol)
    real(r8)  ql_inout(pver, ncol)
    real(r8)  qi_inout(pver, ncol)
    real(r8)  concld_old(pver, ncol)

    real(r8)  nl_inout(pver, ncol)
    real(r8)  ni_inout(pver, ncol)
    real(r8)  nltend(pver, ncol)
    real(r8)  nitend(pver, ncol)

    ! For detraining cumulus condensate into the 'stratus' without evaporation
    ! This is for use in mmacro_pcond
    real(r8)  dlf_T(pver,ncol)
    real(r8)  dlf_qv(pver,ncol)
    real(r8)  dlf_ql(pver,ncol)
    real(r8)  dlf_qi(pver,ncol)
    real(r8)  dlf_nl(pver,ncol)
    real(r8)  dlf_ni(pver,ncol)

    ! Local variables for CFMIP calculations
    real(r8) :: mr_lsliq(pver, ncol)            ! mixing_ratio_large_scale_cloud_liquid (kg/kg)
    real(r8) :: mr_lsice(pver, ncol)            ! mixing_ratio_large_scale_cloud_ice (kg/kg)
    real(r8) :: mr_ccliq(pver, ncol)            ! mixing_ratio_convective_cloud_liquid (kg/kg)
    real(r8) :: mr_ccice(pver, ncol)            ! mixing_ratio_convective_cloud_ice (kg/kg)

    ! CloudSat equivalent ice mass mixing ratio (kg/kg)
    real(r8) :: cldsice(pver, ncol)

    ! ptend_loc, LiXH
    real(r8)  ptend_q(pcnst, pver, ncol)
    real(r8)  ptend_s(pver, ncol)
    real(r8)  ptend_loc_q(pcnst, pver, ncol)
    real(r8)  ptend_loc_s(pver, ncol)

    ! local pstate, LiXH
    real(r8)  state_t(pver, ncol)
    real(r8)  state_s(pver, ncol)
    real(r8)  state_q(pcnst, pver, ncol)
    real(r8)  state_zm(pver, ncol)
    real(r8)  state_zi(pverp, ncol)
    real(r8)  state_pmid(pver, ncol)
    real(r8)  state_pint(pverp, ncol)
    real(r8)  state_pdel(pver, ncol)
    real(r8)  state_phis(ncol)
    real(r8)  qc(pver,ncol)
    real(r8)  qi(pver,ncol)
    real(r8)  nc(pver,ncol)
    real(r8)  ni(pver,ncol)

    real(r8)  landfrac(ncol)          ! Land fraction (fraction)
    real(r8)  snowh(ncol)             ! Snow depth over land, water equivalent (m)

    ! Initialize convective detrainment tendency
    dlf_T(:,:)  = 0._r8
    dlf_qv(:,:) = 0._r8
    dlf_ql(:,:) = 0._r8
    dlf_qi(:,:) = 0._r8
    dlf_nl(:,:) = 0._r8
    dlf_ni(:,:) = 0._r8

    ! Initialize local pstate and ptend
    landfrac(1:ncol) = pstate%landfrac_at_pc_surface%f(1:ncol)
    snowh(1:ncol)    = pstate%snowhland_at_pc_surface%f(1:ncol)


    state_s(:,1:ncol)    = pstate%static_energy_at_pc_full_level%f(:,1:ncol)
    state_t(:,1:ncol)    = pstate%temp_at_pc_full_level%f(:,1:ncol)
    state_q(:,:,1:ncol)  = pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol)
    state_zm(:,1:ncol)   = pstate%z_at_pc_full_level%f(:,1:ncol)
    state_zi(:,1:ncol)   = pstate%z_at_pc_face_level%f(:,1:ncol)
    state_pmid(:,1:ncol) = pstate%pressure_at_pc_full_level%f(:,1:ncol)
    state_pint(:,1:ncol) = pstate%pressure_at_pc_face_level%f(:,1:ncol)
    state_pdel(:,1:ncol) = pstate%delp_at_pc_full_level%f(:,1:ncol)
    state_phis(1:ncol)   = pstate%geop_at_pc_surface%f(1:ncol)

! initialize ptend and ptend_loc
    ptend_q = 0._r8
    ptend_s = 0._r8
    ptend_loc_q = 0._r8
    ptend_loc_s = 0._r8

    ! Procedures :
    ! (1) Partition detrained convective cloud water into liquid and ice based on T.
    !     This also involves heating.
    !     If convection scheme can handle this internally, this step is not necssary.
    ! (2) Assuming a certain effective droplet radius, computes number concentration
    !     of detrained convective cloud liquid and ice.
    ! (3) If 'cu_det_st = .true' ('false'), detrain convective cloud 'liquid' into 
    !     the pre-existing 'liquid' stratus ( mean environment ).  The former does
    !     not involve any macrophysical evaporation while the latter does. This is
    !     a kind of 'targetted' deposition. Then, force in-stratus LWC to be bounded 
    !     by qcst_min and qcst_max in mmacro_pcond.
    ! (4) In contrast to liquid, convective ice is detrained into the environment 
    !     and involved in the sublimation. Similar bounds as liquid stratus are imposed.
    ! This is the key procesure generating upper-level cirrus clouds.
    ! The unit of dlf : [ kg/kg/s ]

    det_s(:)   = 0._r8
    det_ice(:) = 0._r8
    dpdlfliq   = 0._r8
    dpdlfice   = 0._r8
    shdlfliq   = 0._r8
    shdlfice   = 0._r8
    dpdlft     = 0._r8
    shdlft     = 0._r8

    do k = top_lev, pver
    do i = 1, ncol
       if( state_t(k,i) > 268.15_r8 ) then
           dum1 = 0.0_r8
       elseif( state_t(k,i) < 238.15_r8 ) then
           dum1 = 1.0_r8
       else
           dum1 = ( 268.15_r8 - state_t(k,i) ) / 30._r8
       endif

       ! If detrainment was done elsewhere, still update the variables used for output
       ! assuming that the temperature split between liquid and ice is the same as assumed here.
       if (do_detrain) then
           ptend_loc_q(ixcldliq,k,i) = dlf(k,i) * ( 1._r8 - dum1 )
           ptend_loc_q(ixcldice,k,i) = dlf(k,i) * dum1
           !dum2                  = dlf(k,i) * ( 1._r8 - dum1 )
           ptend_loc_q(ixnumliq,k,i) = 3._r8 * ( max(0._r8, ( dlf(k,i) - dlf2(k,i) )) * ( 1._r8 - dum1 ) ) / &
                (4._r8*3.14_r8* Rdet_dp_liq**3*997._r8) + & ! Deep    Convection
                3._r8 * (                         dlf2(k,i)    * ( 1._r8 - dum1 ) ) / &
                (4._r8*3.14_r8*Rdet_sh_liq**3*997._r8)     ! Shallow Convection 
           !dum2                  = dlf(k,i) * dum1
           ptend_loc_q(ixnumice,k,i) = 3._r8 * ( max(0._r8, ( dlf(k,i) - dlf2(k,i) )) *  dum1 ) / &
                (4._r8*3.14_r8*Rdet_dp_ice**3*500._r8) + & ! Deep    Convection
                3._r8 * (                         dlf2(k,i)    *  dum1 ) / &
                (4._r8*3.14_r8*Rdet_sh_ice**3*500._r8)     ! Shallow Convection
           ptend_loc_s(k,i)          = dlf(k,i) * dum1 * latice
       else 
           ptend_loc_q(ixcldliq,k,i) = 0._r8
           ptend_loc_q(ixcldice,k,i) = 0._r8
           ptend_loc_q(ixnumliq,k,i) = 0._r8
           ptend_loc_q(ixnumice,k,i) = 0._r8
           ptend_loc_s(k,i)          = 0._r8
       end if
 
       ! Only rliq is saved from deep convection, which is the reserved liquid.  We need to keep
       ! track of the integrals of ice and static energy that is effected from conversion to ice
       ! so that the energy checker doesn't complain.	
       det_s(i)                  = det_s(i) + ptend_loc_s(k,i)*state_pdel(k,i)/gravity
       det_ice(i)                = det_ice(i) - ptend_loc_q(ixcldice,k,i)*state_pdel(k,i)/gravity      

       ! Targetted detrainment of convective liquid water either directly into the
       ! existing liquid stratus or into the environment. 
       if( cu_det_st ) then
           dlf_T(k,i)  = ptend_loc_s(k,i)/cp
           dlf_qv(k,i) = 0._r8
           dlf_ql(k,i) = ptend_loc_q(ixcldliq,k,i)
           dlf_qi(k,i) = ptend_loc_q(ixcldice,k,i)
           dlf_nl(k,i) = ptend_loc_q(ixnumliq,k,i)
           dlf_ni(k,i) = ptend_loc_q(ixnumice,k,i)         
           ptend_loc_q(ixcldliq,k,i) = 0._r8
           ptend_loc_q(ixcldice,k,i) = 0._r8
           ptend_loc_q(ixnumliq,k,i) = 0._r8
           ptend_loc_q(ixnumice,k,i) = 0._r8
           ptend_loc_s(k,i)          = 0._r8
           dpdlfliq(k,i)             = 0._r8
           dpdlfice(k,i)             = 0._r8
           shdlfliq(k,i)             = 0._r8
           shdlfice(k,i)             = 0._r8
           dpdlft  (k,i)             = 0._r8
           shdlft  (k,i)             = 0._r8
        else
           dpdlfliq(k,i) = ( dlf(k,i) - dlf2(k,i) ) * ( 1._r8 - dum1 )
           dpdlfice(k,i) = ( dlf(k,i) - dlf2(k,i) ) * ( dum1 )
           shdlfliq(k,i) = dlf2(k,i) * ( 1._r8 - dum1 )
           shdlfice(k,i) = dlf2(k,i) * ( dum1 )
           dpdlft  (k,i) = ( dlf(k,i) - dlf2(k,i) ) * dum1 * latice/cp
           shdlft  (k,i) = dlf2(k,i) * dum1 * latice/cp
       endif
    end do
    end do

    ! output: 'DPDLFLIQ'dpdlfliq  'DPDLFICE'dpdlfice  'SHDLFLIQ'shdlfliq  'SHDLFICE'shdlfice
    !         'DPDLFT'dpdlft  'SHDLFT'shdlft  'ZMDLF'dlf  

    det_ice(:ncol) = det_ice(:ncol)/1000._r8  ! divide by density of water
    
    ! update local state and ptend
    ptend_q(2:5,:,1:ncol) = ptend_q(2:5,:,1:ncol)+ptend_loc_q(2:5,:,1:ncol)
    ptend_s = ptend_s+ptend_loc_s

    state_q(2:5,:,1:ncol) = state_q(2:5,:,1:ncol) + ptend_loc_q(2:5,:,1:ncol)*dtime
    state_s(:,1:ncol)     = state_s(:,1:ncol) + ptend_loc_s(:,1:ncol)*dtime

    ! check negetive qv and correct to qmin
    do m = 2, 5
       if(m .eq. ixnumliq .or. m .eq. ixnumice)then
            do k = 1, pver
               do i = 1, ncol
                  state_q(m,k,i) = max(1.e-12_r8,state_q(m,k,i))
                  state_q(m,k,i) = min(1.e10_r8,state_q(m,k,i))
                end do
            end do
        else
            call qneg3('macro_physics_1', ncol, pver, m, m, phy_tracer_info(m)%qmin, state_q(m,:,1:ncol))
            !where(state_q(m,:,1:ncol) .lt. phy_tracer_info(m)%qmin) state_q(m,:,1:ncol) = phy_tracer_info(m)%qmin
        end if
    end do

    call geopotential_dse(ncol,       state_pint, state_pmid, state_pdel,   &
                          state_phis, state_s,    state_q(1,:,1:ncol),      &
                          state_t,    state_zm,   state_zi)


! Computation of Various Cloud Fractions 
!     . Cumulus AMT = Deep    Cumulus AMT ( empirical fcn of mass flux ) +      
!                     Shallow Cumulus AMT ( internally fcn of mass flux and w ) 
!     . Stratus AMT = fcn of environmental-mean RH ( no Stability Stratus )     
!     . Cumulus and Stratus are non-overlapped with higher priority on Cumulus  
!     . Cumulus ( both Deep and Shallow ) has its own LWC and IWC.              

    concld_old(top_lev:pver,:ncol) = pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)

    allocate(shfrc(pver, ncol))
    if( use_shfrc ) then
        shfrc(:,1:ncol) = pstate_cam%sh_shfrc_at_pc_full_level%f(:,1:ncol)
    else 
        shfrc(:,:)      = 0._r8
    endif

    ! CAM5 only uses 'concld' output from the below subroutine. 
    ! Stratus ('ast' = max(alst,aist)) and total cloud fraction ('cld = ast + concld')
    ! will be computed using this updated 'concld' in the stratiform macrophysics 
    ! scheme (mmacro_pcond) later below. 
    call cldfrc( ncol      ,                                                                &
                 shfrc     , use_shfrc,                                                     &
                 pstate_cam%updraft_mass_flux%f(:,1:ncol), cmfmc2    ,    &
                 pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,1:ncol)   )

   rdtime = 1._r8/dtime


! initialize ptend_loc
   ptend_loc_s = 0._r8
   ptend_loc_q = 0._r8

! Liquid Macrop_Driver Macrophysics
   qc = 0._r8; qi = 0._r8; nc = 0._r8; ni = 0._r8 
   qc(top_lev:pver,:ncol) = state_q(ixcldliq,top_lev:pver,:ncol)
   qi(top_lev:pver,:ncol) = state_q(ixcldice,top_lev:pver,:ncol)
   nc(top_lev:pver,:ncol) = state_q(ixnumliq,top_lev:pver,:ncol)
   ni(top_lev:pver,:ncol) = state_q(ixnumice,top_lev:pver,:ncol)

 ! In CAM5, 'microphysical forcing' ( CC_... ) and 'the other advective forcings' ( ttend, ... ) 
 ! are separately provided into the prognostic microp_driver macrophysics scheme. This is an
 ! attempt to resolve in-cloud and out-cloud forcings. 

   if( istep .le. 1 ) then
       pstate_cam%macrop_tcwat_at_pc_full_level%f(old_time_level,:,:ncol)   = state_t(:,:ncol)
       pstate_cam%macrop_qcwat_at_pc_full_level%f(old_time_level,:,:ncol)   = state_q(1,:,:ncol)
       pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,:,:ncol)   = qc(:,:ncol) + qi(:,:ncol)
       pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,:,:ncol)  = qi(:,:ncol)
       pstate_cam%macrop_nlwat_at_pc_full_level%f(old_time_level,:,:ncol)   = nc(:,:ncol)
       pstate_cam%macrop_niwat_at_pc_full_level%f(old_time_level,:,:ncol)   = ni(:,:ncol)
       ttend(:,:ncol)   = 0._r8
       qtend(:,:ncol)   = 0._r8
       ltend(:,:ncol)   = 0._r8
       itend(:,:ncol)   = 0._r8
       nltend(:,:ncol)  = 0._r8
       nitend(:,:ncol)  = 0._r8

       pstate_cam%microp_cc_t_at_pc_full_level%f(old_time_level,:,:ncol)    = 0._r8
       pstate_cam%microp_cc_qv_at_pc_full_level%f(old_time_level,:,:ncol)   = 0._r8
       pstate_cam%microp_cc_ql_at_pc_full_level%f(old_time_level,:,:ncol)   = 0._r8
       pstate_cam%microp_cc_qi_at_pc_full_level%f(old_time_level,:,:ncol)   = 0._r8
       pstate_cam%microp_cc_nl_at_pc_full_level%f(old_time_level,:,:ncol)   = 0._r8
       pstate_cam%microp_cc_ni_at_pc_full_level%f(old_time_level,:,:ncol)   = 0._r8
       pstate_cam%microp_cc_qlst_at_pc_full_level%f(old_time_level,:,:ncol) = 0._r8
   else
       ttend(top_lev:pver,:ncol)   = ( state_t(top_lev:pver,:ncol)                                      &
            - pstate_cam%macrop_tcwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)) * rdtime   &
            - pstate_cam%microp_cc_t_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol) 
       qtend(top_lev:pver,:ncol)   = ( state_q(1,top_lev:pver,:ncol)                                    &
            - pstate_cam%macrop_qcwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)) * rdtime   &
            - pstate_cam%microp_cc_qv_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
       ltend(top_lev:pver,:ncol)   = ( qc(top_lev:pver,:ncol) + qi(top_lev:pver,:ncol)                  &
            - pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol) ) * rdtime  &
            - (pstate_cam%microp_cc_ql_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)            &
            + pstate_cam%microp_cc_qi_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol))
       itend(top_lev:pver,:ncol)   = ( qi(top_lev:pver,:ncol)                                           &
            - pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)) * rdtime  &
            - pstate_cam%microp_cc_qi_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
       nltend(top_lev:pver,:ncol)  = ( nc(top_lev:pver,:ncol)                                           &
            - pstate_cam%macrop_nlwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)) * rdtime   &
            - pstate_cam%microp_cc_nl_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
       nitend(top_lev:pver,:ncol)  = ( ni(top_lev:pver,:ncol)                                           &      
            - pstate_cam%macrop_niwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)) * rdtime   &
            - pstate_cam%microp_cc_ni_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
   endif

   lmitend(top_lev:pver,:ncol) = ltend(top_lev:pver,:ncol) - itend(top_lev:pver,:ncol)

   t_inout(top_lev:pver,:ncol)  =  pstate_cam%macrop_tcwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol) 
   qv_inout(top_lev:pver,:ncol) =  pstate_cam%macrop_qcwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
   ql_inout(top_lev:pver,:ncol) =  pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)    &
                                -  pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
   qi_inout(top_lev:pver,:ncol) =  pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
   nl_inout(top_lev:pver,:ncol) =  pstate_cam%macrop_nlwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)
   ni_inout(top_lev:pver,:ncol) =  pstate_cam%macrop_niwat_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol)

 ! Liquid Microp_Driver Macrophysics.
 ! The main roles of this subroutines are
 ! (1) compute net condensation rate of stratiform liquid ( cmeliq )
 ! (2) compute liquid stratus and ice stratus fractions. 
 ! Note 'ttend...' are advective tendencies except microphysical process while
 !      'CC...'    are microphysical tendencies. 

   call mmacro_pcond( ncol, dtime, state_pmid, state_pdel,                             &
                      t_inout, qv_inout, ql_inout, qi_inout, nl_inout, ni_inout,       &                  
                      ttend, qtend, lmitend, itend, nltend, nitend,                    &
                      pstate_cam%microp_cc_t_at_pc_full_level%f(old_time_level,:,:ncol),   &
                      pstate_cam%microp_cc_qv_at_pc_full_level%f(old_time_level,:,:ncol),  &
                      pstate_cam%microp_cc_ql_at_pc_full_level%f(old_time_level,:,:ncol),  &
                      pstate_cam%microp_cc_qi_at_pc_full_level%f(old_time_level,:,:ncol),  &
                      pstate_cam%microp_cc_nl_at_pc_full_level%f(old_time_level,:,:ncol),  &
                      pstate_cam%microp_cc_ni_at_pc_full_level%f(old_time_level,:,:ncol),  &
                      pstate_cam%microp_cc_qlst_at_pc_full_level%f(old_time_level,:,:ncol),&
                      dlf_T, dlf_qv, dlf_ql, dlf_qi, dlf_nl, dlf_ni,                   &
                      concld_old,                                                      &
                      pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,:ncol), &
                      landfrac, snowh,                                                 &
                      tlat, qvlat, qcten, qiten, ncten, niten,                         &
                      pstate_cam%macrop_cmeliq_at_pc_full_level%f(:,:ncol),            &
                      qvadj, qladj, qiadj, qllim, qilim,                               &
                      pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,:ncol),    &
                      pstate_cam%macrop_alst_at_pc_full_level%f(old_time_level,:,:ncol),   &
                      pstate_cam%macrop_aist_at_pc_full_level%f(old_time_level,:,:ncol),   &
                      pstate_cam%macrop_qlst_at_pc_full_level%f(old_time_level,:,:ncol),   &
                      pstate_cam%macrop_qist_at_pc_full_level%f(old_time_level,:,:ncol),   &
                      do_cldice ) 

 ! Compute net stratus fraction using maximum over-lapping assumption
   pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:top_lev-1,:ncol) = 0._r8
   pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol) = &
        max( pstate_cam%macrop_alst_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol), &
             pstate_cam%macrop_aist_at_pc_full_level%f(old_time_level,top_lev:pver,:ncol) )

   do k = top_lev, pver
      do i = 1, ncol
         ptend_loc_s(k,i)          =  tlat(k,i)
         ptend_loc_q(1,k,i)        = qvlat(k,i)
         ptend_loc_q(ixcldliq,k,i) = qcten(k,i)
         ptend_loc_q(ixcldice,k,i) = qiten(k,i)
         ptend_loc_q(ixnumliq,k,i) = ncten(k,i)
         ptend_loc_q(ixnumice,k,i) = niten(k,i)

         ! Check to make sure that the macrophysics code is respecting the flags that control
         ! whether cldwat should be prognosing cloud ice and cloud liquid or not.
         if ((.not. do_cldice) .and. (qiten(k,i) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR - "// &
                 "Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice mass tendencies.")
         end if
         if ((.not. do_cldice) .and. (niten(k,i) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR -"// &
                 " Cldwat is configured not to prognose cloud ice, but mmacro_pcond has ice number tendencies.")
         end if

         if ((.not. do_cldliq) .and. (qcten(k,i) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR - "// &
                 "Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid mass tendencies.")
         end if
         if ((.not. do_cldliq) .and. (ncten(k,i) /= 0.0_r8)) then 
            call endrun("macrop_driver:ERROR - "// &
                 "Cldwat is configured not to prognose cloud liquid, but mmacro_pcond has liquid number tendencies.")
         end if
      end do
   end do

   ! update local state and ptend
   ptend_q = ptend_q+ptend_loc_q
   ptend_s = ptend_s+ptend_loc_s

   state_q(1:5,:,1:ncol) = state_q(1:5,:,1:ncol) + ptend_loc_q(1:5,:,1:ncol)*dtime
   state_s(:,1:ncol)     = state_s(:,1:ncol) + ptend_loc_s(:,1:ncol)*dtime

   ! check negetive qv and correct to qmin
   do m = 1, 5
       if(m .eq. ixnumliq .or. m .eq. ixnumice)then
           do k = 1, pver
              do i = 1, ncol
                 state_q(m,k,i) = max(1.e-12_r8,state_q(m,k,i))
                 state_q(m,k,i) = min(1.e10_r8,state_q(m,k,i))
               end do
           end do
       else
           call qneg3('macro_physics_2', ncol, pver, m, m, phy_tracer_info(m)%qmin, state_q(m,:,1:ncol))
           !where(state_q(m,:,1:ncol) .lt. phy_tracer_info(m)%qmin) state_q(m,:,1:ncol) = phy_tracer_info(m)%qmin
       end if
   end do

   call geopotential_dse(ncol,       state_pint, state_pmid, state_pdel,   &
                         state_phis, state_s,    state_q(1,:,1:ncol),      &
                         state_t,    state_zm,   state_zi)

   ! output: 'MACPDT'tlat  'MACPDQ'qvlat  'MACPDLIQ'qcten  'MACPDICE'qiten  'CLDVAPADJ'qvadj  'CLDLIQADJ'qladj
   !         'CLDICEADJ'qiadj  'CLDLIQDET'dlf_ql  'CLDICEDET'dlf_qi  'CLDLIQLIM'qllim  'CLDICELIM'qilim
   !         'AST'ast  'CONCLD'concld  'CMELIQ'cmeliq

   ! calculations and outfld calls for CLDLIQSTR, CLDICESTR, CLDLIQCON, CLDICECON for CFMIP

   ! initialize local variables
   mr_ccliq = 0._r8   !! not seen by radiation, so setting to 0 
   mr_ccice = 0._r8   !! not seen by radiation, so setting to 0
   mr_lsliq = 0._r8
   mr_lsice = 0._r8

   do k=top_lev,pver
      do i=1,ncol
         if (pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i) .gt. 0._r8) then
            mr_lsliq(k,i) = state_q(ixcldliq,k,i)
            mr_lsice(k,i) = state_q(ixcldice,k,i)
         else
            mr_lsliq(k,i) = 0._r8
            mr_lsice(k,i) = 0._r8
         end if
      end do
   end do

   ! output: 'CLDLIQSTR'mr_lsliq  'CLDICESTR'mr_lsice  'CLDLIQCON'mr_ccliq  'CLDICECON'mr_ccice

   ! Save equilibrium state variables for macrophysics at the next time step  
   cldsice = 0._r8
   do k = top_lev, pver
      pstate_cam%macrop_tcwat_at_pc_full_level%f(old_time_level,k,:ncol)  = state_t(k,:ncol)
      pstate_cam%macrop_qcwat_at_pc_full_level%f(old_time_level,k,:ncol)  = state_q(1,k,:ncol)
      pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,k,:ncol)  = state_q(ixcldliq,k,:ncol) + state_q(ixcldice,k,:ncol)
      pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,k,:ncol) = state_q(ixcldice,k,:ncol)
      pstate_cam%macrop_nlwat_at_pc_full_level%f(old_time_level,k,:ncol)  = state_q(ixnumliq,k,:ncol)
      pstate_cam%macrop_niwat_at_pc_full_level%f(old_time_level,k,:ncol)  = state_q(ixnumice,k,:ncol)
      cldsice(k,:ncol)               = pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,k,:ncol) &
            * min(1.0_r8, max(0.0_r8, (tmelt - pstate_cam%macrop_tcwat_at_pc_full_level%f(old_time_level,k,:ncol)) / 20._r8))
   end do

   !---------------LiXH add vars checking------------------- 
   if(any(isnan(pstate_cam%macrop_qcwat_at_pc_full_level%f(old_time_level,:,:))))then
       print*,'NAN appears in qcwat in macro-physics'
       print*,'rank=',mpi_rank()
       call endrun('check_non')
   end if 
   if(any(isnan(pstate_cam%macrop_lcwat_at_pc_full_level%f(old_time_level,:,:))))then
       print*,'NAN appears in lcwat in macro-physics'
       print*,'rank=',mpi_rank()
       call endrun('check_non')
   end if 
   if(any(isnan(pstate_cam%macrop_iccwat_at_pc_full_level%f(old_time_level,:,:))))then
       print*,'NAN appears in iccwat in macro-physics'
       print*,'rank=',mpi_rank()
       call endrun('check_non')
   end if 
   if(any(isnan(pstate_cam%macrop_nlwat_at_pc_full_level%f(old_time_level,:,:))))then
       print*,'NAN appears in nlwat in macro-physics'
       print*,'rank=',mpi_rank()
       call endrun('check_non')
   end if 
   if(any(isnan(pstate_cam%macrop_niwat_at_pc_full_level%f(old_time_level,:,:))))then
       print*,'NAN appears in niwat in macro-physics'
       print*,'rank=',mpi_rank()
       call endrun('check_non')
   end if 
   !---------------LiXH add vars checking------------------- 

   ! output: 'CLDSICE'cldsice

   ! re-initialize ptend
   ptend_macrophysics%tend_s%f = 0._r8
   ptend_macrophysics%tend_q%f = 0._r8
 
   ! Add the detrainment tendency to the output tendency
   ptend_macrophysics%tend_s%f(:,1:ncol)          = ptend_s(:,1:ncol)
   ptend_macrophysics%tend_q%f(1,:,1:ncol)        = ptend_q(1,:,1:ncol)
   ptend_macrophysics%tend_q%f(ixcldliq,:,1:ncol) = ptend_q(ixcldliq,:,1:ncol)
   ptend_macrophysics%tend_q%f(ixcldice,:,1:ncol) = ptend_q(ixcldice,:,1:ncol)
   ptend_macrophysics%tend_q%f(ixnumliq,:,1:ncol) = ptend_q(ixnumliq,:,1:ncol)
   ptend_macrophysics%tend_q%f(ixnumice,:,1:ncol) = ptend_q(ixnumice,:,1:ncol)
 
   deallocate(shfrc)

   end subroutine macrop_driver_tend_park

 end module grist_macrop_park
