!===================================================================================================
!
!  Created by LiXiaohan on 19/07/10, adopted from CAM5
!
!     interface of the prognostic cloud macrophysics
!
!===================================================================================================

 module grist_macrop_lin

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8, i4, latice, gravity,  &
                                                  cp, tmelt

    implicit none
    public          :: macrop_driver_init_lin,        &
                       read_nml_macrop_lin,           &
                       macrop_driver_tend_lin,        &
                       end_of_macrop_lin

    private
    save

! Private:
    logical :: do_cldice             ! .true., macrophysics is prognosing cldice
    logical :: do_cldliq             ! .true., macrophysics is prognosing cldliq
    logical :: do_detrain            ! .true., macrophysics is detraining ice into stratiform
    logical :: use_shfrc             ! Flag, whether use sh_cld_frc from shallow_conv (LiXH)

    logical, parameter :: cu_det_st  = .false.  

    integer :: pcnst, pver, pverp
    integer :: ixcldliq, ixcldice, ixnumliq, ixnumice

! Tuning parameters, LiXH:
     !Rdet: effective radius of detrained convective condensate
     real(r8) , parameter :: Rdet_dp_liq = 8.e-6_r8
     real(r8) , parameter :: Rdet_dp_ice = 25.e-6_r8
     real(r8) , parameter :: Rdet_sh_liq = 8.e-6_r8  !default: 10.e-6_r8 
     real(r8) , parameter :: Rdet_sh_ice = 25.e-6_r8 !default: 50.e-6_r8

 contains

    subroutine read_nml_macrop_lin(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer :: unitn, ierr
    logical :: macro_do_cldice  = .true.  ! macrophysics is prognosing cldice
    logical :: macro_do_cldliq  = .true.  ! macrophysics is prognosing cldliq
    logical :: macro_do_detrain = .true.  ! macrophysics is detraining ice into stratiform

    namelist /macro_nl/ macro_do_cldice, macro_do_cldliq, macro_do_detrain

    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, macro_nl, iostat=ierr)
    if (ierr /= 0) call endrun(' error reading macro namelist ')
    close(unitn)
 
    do_cldice  = macro_do_cldice
    do_cldliq  = macro_do_cldliq
    do_detrain = macro_do_detrain

    end subroutine read_nml_macrop_lin


    subroutine macrop_driver_init_lin(ncol)
    use grist_nml_module,               only: nlev, nlevp, ntracer
    use grist_physics_data_structure,   only: pstate, phy_tracer_info
    use grist_PhysC_data_structure,     only: pstate_cam
    use grist_physics_update,           only: ptimelevels

! io
    integer, intent(in) :: ncol
! local
    integer :: m

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


    allocate(pstate_cam%macrop_sgm_at_pc_full_level%f(nlev,ncol))

    pstate_cam%macrop_ast_at_pc_full_level%pos   = 0
    pstate_cam%macrop_alst_at_pc_full_level%pos  = 0
    pstate_cam%macrop_aist_at_pc_full_level%pos  = 0
    pstate_cam%macrop_qlst_at_pc_full_level%pos  = 0
    pstate_cam%macrop_qist_at_pc_full_level%pos  = 0
    pstate_cam%macrop_cld_at_pc_full_level%pos   = 0
    pstate_cam%macrop_concld_at_pc_full_level%pos= 0
    pstate_cam%macrop_cmeliq_at_pc_full_level%pos= 0
    pstate_cam%macrop_sgm_at_pc_full_level%pos   = 0
 
    pstate_cam%macrop_ast_at_pc_full_level%f     = 0._r8
    pstate_cam%macrop_alst_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_aist_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_qlst_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_qist_at_pc_full_level%f    = 0._r8
    pstate_cam%macrop_cld_at_pc_full_level%f     = 0._r8
    pstate_cam%macrop_concld_at_pc_full_level%f  = 0._r8
    pstate_cam%macrop_cmeliq_at_pc_full_level%f  = 0._r8
    pstate_cam%macrop_sgm_at_pc_full_level%f     = 0._r8

    pcnst = ntracer
    pver  = nlev
    pverp = nlevp

    ! inquire cloud index
    do m = 1, pcnst
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')        ixcldliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')           ixcldice = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid_number') ixnumliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice_number')    ixnumice = m 
    end do

!-------------------Lin Macrop------------------
        use_shfrc = .true.              !LiXH Test
!-------------------Lin Macrop------------------

    end subroutine macrop_driver_init_lin


    subroutine end_of_macrop_lin

    use grist_PhysC_data_structure,   only: pstate_cam

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

    if(allocated(pstate_cam%macrop_sgm_at_pc_full_level%f))      &
      deallocate(pstate_cam%macrop_sgm_at_pc_full_level%f)  
#endif

    end subroutine end_of_macrop_lin


    subroutine macrop_driver_tend_lin(ncol  , dtime , istep   ,     &
                                  dlf   , dlf2  ,  cmfmc2 ,     &
                                  det_s , det_ice)
    use cloud_fraction,                 only: top_lev => trop_cloud_top_lev, &
                                              cldfrc, cldfrc_fice
    use cldwat2m_macro_lin,             only: mmacro_pcond_lin
    use grist_physics_data_structure,   only: pstate, phy_tracer_info
    use grist_PhysC_data_structure,     only: pstate_cam,                    &
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

    ! For macrophysics
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
  
!-------------------Lin Macrop------------------>
   real(r8) alst_gp(pver,ncol)
   real(r8) cond_gp(pver,ncol)
   real(r8) N1(pver,ncol)
 
!   real(r8) sgm(pver,ncol)
 
   real(r8) delta_q(pver,ncol)
   real(r8) deltaq_sat(pver,ncol)
   real(r8) deltaq_uns(pver,ncol)
   real(r8) Q1_sat(pver,ncol)
   real(r8) Q1_uns(pver,ncol)
 
   real(r8) adjust_factor(pver,ncol)
!<-------------------Lin Macrop------------------

    ! For revised macophysics, mmacro_pcond
    real(r8)  t_inout(pver, ncol)
    real(r8)  qv_inout(pver, ncol)
    real(r8)  ql_inout(pver, ncol)
    real(r8)  qi_inout(pver, ncol)
    real(r8)  nl_inout(pver, ncol)
    real(r8)  ni_inout(pver, ncol)

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
    !-------------------Lin Macrop---------------------
    real(r8)  state_rpdel(pver, ncol)
    !-------------------Lin Macrop---------------------
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
    !-------------------Lin Macrop---------------------
    state_rpdel(:,1:ncol)= 1._r8/state_pdel(:,1:ncol) 
    !-------------------Lin Macrop---------------------

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
    allocate(shfrc(pver, ncol))
    if( use_shfrc ) then
        shfrc(:,1:ncol) = pstate_cam%sh_shfrc_at_pc_full_level%f(:,1:ncol)
    else 
        shfrc(:,:)      = 0._r8
    endif

    ! only uses 'concld' output from the below subroutine. 
    call cldfrc( ncol      ,                                                                &
                 shfrc     , use_shfrc,                                                     &
                 pstate_cam%updraft_mass_flux%f(:,1:ncol), cmfmc2    ,  &
                 pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,1:ncol) )

! initialize ptend_loc
   ptend_loc_s = 0._r8
   ptend_loc_q = 0._r8

! Liquid Macrop_Driver Macrophysics
   qc = 0._r8; qi = 0._r8; nc = 0._r8; ni = 0._r8 
   qc(top_lev:pver,:ncol) = state_q(ixcldliq,top_lev:pver,:ncol)
   qi(top_lev:pver,:ncol) = state_q(ixcldice,top_lev:pver,:ncol)
   nc(top_lev:pver,:ncol) = state_q(ixnumliq,top_lev:pver,:ncol)
   ni(top_lev:pver,:ncol) = state_q(ixnumice,top_lev:pver,:ncol)

   t_inout(top_lev:pver,:ncol)  =  state_t(top_lev:pver,:ncol) 
   qv_inout(top_lev:pver,:ncol) =  state_q(1,top_lev:pver,:ncol)
   ql_inout(top_lev:pver,:ncol) =  qc(top_lev:pver,:ncol)
   qi_inout(top_lev:pver,:ncol) =  qi(top_lev:pver,:ncol) 
   nl_inout(top_lev:pver,:ncol) =  nc(top_lev:pver,:ncol)
   ni_inout(top_lev:pver,:ncol) =  ni(top_lev:pver,:ncol)


   call mmacro_pcond_lin( ncol, dtime, state_pmid, state_pdel,                         &
                      t_inout, qv_inout, ql_inout, qi_inout, nl_inout, ni_inout,       &                  
                      pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,:ncol), &
                      landfrac, snowh,                                                 &
                      tlat, qvlat, qcten, qiten, ncten, niten,                         &
                      pstate_cam%macrop_cmeliq_at_pc_full_level%f(:,:ncol),            &
                      pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,:ncol),&
                      pstate_cam%macrop_alst_at_pc_full_level%f(old_time_level,:,:ncol), &
                      pstate_cam%macrop_aist_at_pc_full_level%f(old_time_level,:,:ncol), &
                      pstate_cam%macrop_qlst_at_pc_full_level%f(old_time_level,:,:ncol), &
                      pstate_cam%macrop_qist_at_pc_full_level%f(old_time_level,:,:ncol), &
                      do_cldice,                                                       &
                      state_rpdel, state_zm, state_zi,                                 &
                      pstate_cam%pbl_lengi_at_pc_face_level%f(:,:ncol),                &
                      pstate_cam%pbl_shi_at_pc_face_level%f(:,:ncol),                  &
                      pstate_cam%pbl_wstarPBL%f(:ncol),                                &
                      pstate_cam%sh_qtu_shal_at_pc_face_level%f(:,:ncol),              &
                      pstate_cam%sh_umf_shal_at_pc_face_level%f(:,:ncol),              &
                      pstate_cam%cumulus_cldtop2%f(:ncol),                             &
                      pstate_cam%cumulus_cldbot2%f(:ncol),                             &
                      pstate_cam%sh_thlu_shal_at_pc_face_level%f(:,:ncol),             &
                      alst_gp, cond_gp, N1,                                  &
                      pstate_cam%macrop_sgm_at_pc_full_level%f(:,:ncol),               &
                      delta_q, deltaq_sat, deltaq_uns, Q1_sat, Q1_uns,                 &
                      adjust_factor ) 

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

   ! output: 'MACPDT'tlat  'MACPDQ'qvlat  'MACPDLIQ'qcten  'MACPDICE'qiten
   !         'CLDLIQDET'dlf_ql  'CLDICEDET'dlf_qi
   !          'AST'ast  'CONCLD'concld  'CMELIQ'cmeliq

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

   end subroutine macrop_driver_tend_lin

 end module grist_macrop_lin
