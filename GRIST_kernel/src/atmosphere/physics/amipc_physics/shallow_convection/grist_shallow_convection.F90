!===================================================================================================
!
!  Created by LiXiaohan on 19/07/01, adopted from CAM5
!
!     interface of the shallow convection scheme 
!     only the UW scheme is avaiable now
!===================================================================================================

 module grist_shallow_convection

    use grist_handle_error,                 only: endrun
    use grist_nml_module,                   only: ntracer, nlev, nlevp
    use grist_constants,                    only: r8, i4, cp, zvir, p00, &
                                                  latvap, latice, rdry,  &
                                                  gravity, mwh2o, mwdry
    use grist_physics_update,               only: ptimelevels
    use grist_mpi

    implicit none
    private
    public          :: convect_shallow_init,        &
                       end_of_convect_shallow,      &
                       convect_shallow_tend,        &
                       convect_shallow_use_shfrc,   &
                       trigger_function_xie

! Private:
    character(len=16) :: shallow_scheme

! if the variables below are used by other modules (pbuf in CAM), change to pstate%~, Lixh
    real(r8),    allocatable     :: flxprec(:,:)    ! 'SH_FLXPRC'   uw grid-box mean rain+snow flux [kg m^-2 s^-1]
    real(r8),    allocatable     :: flxsnow(:,:)    ! 'SH_FLXSNW'   uw grid-box mean snow flux [kg m^-2 s^-1]
    real(r8),    allocatable     :: evapcsh(:,:)    ! 'NEVAPR_SHCU' Evaporation of precipitation [ kg/kg/s ]
    real(r8),    allocatable     :: snow_sh(:)      ! 'SNOW_SH'     Shallow convective snow flux at surface [m/s]
! other pbuf includes:
!   'SH_CLDLIQ', 'SH_CLDLIQ' 

 contains

    subroutine convect_shallow_init(ncol, dxmean)
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_uwshcu,                       only: init_uwshcu
    use grist_double_plume,                 only: init_double_plume
    use phys_control,                       only: phys_getopts
!------------- for Double Plume scheme, LiXH ------------>
    use grist_hpe_constants,                only: eta_face
!<------------ for Double Plume scheme, LiXH -------------
! io
    integer(i4), intent(in)     :: ncol
    real(r8),    intent(in)     :: dxmean(ncol)
!------------- for Double Plume scheme, LiXH ------------>
    ! local
    real(r8)                    :: pref_edge(nlevp)        ! Reference pressures at interfaces
    integer                     :: limcnv                  ! Top interface level limit for convection
    integer                     :: k
!<------------ for Double Plume scheme, LiXH -------------

    call phys_getopts(shallow_scheme_out = shallow_scheme)

    select case (shallow_scheme)

    case('off')
        if(mpi_rank()==0)print*,'onvect_shallow_init: shallow convection OFF'

    case('UW')
        allocate(flxprec(nlevp,ncol))   ; flxprec = 0._r8
        allocate(flxsnow(nlevp,ncol))   ; flxsnow = 0._r8 
        call init_uwshcu(latvap, cp, latice, zvir, rdry, gravity, mwh2o/mwdry)
!------------- for Double Plume scheme, LiXH ------------>
    case('double_plume')
        allocate(flxprec(nlevp,ncol))   ; flxprec = 0._r8
        allocate(flxsnow(nlevp,ncol))   ; flxsnow = 0._r8 
 
        pref_edge(:) = eta_face(:)*p00
        if( pref_edge(1) >= 4.e3_r8 ) then
           limcnv = 1
        else
             do k = 1, nlev
                if( pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8 ) then
                  limcnv = k
                  goto 11
                 endif
             enddo
             limcnv  = nlevp
        end if
        11 continue

        call init_double_plume( latvap, cp, latice, zvir, rdry, gravity, mwh2o/mwdry, limcnv, ncol, dxmean )
!<------------ for Double Plume scheme, LiXH -------------

    end select

    allocate(evapcsh(nlev,ncol))    ; evapcsh = 0._r8
    allocate(snow_sh(ncol))         ; snow_sh = 0._r8

    ! Initialize pstate:
    allocate(pstate_cam%scalar_precc_sh_surface%f(ncol))
    allocate(pstate_cam%scalar_snowc_surface%f(ncol))
    allocate(pstate_cam%sh_icwmr_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%cumulus_cldtop%f(ncol))
    allocate(pstate_cam%cumulus_cldbot%f(ncol))
    allocate(pstate_cam%cumulus_cldtop2%f(ncol))
    allocate(pstate_cam%cumulus_cldbot2%f(ncol))
    allocate(pstate_cam%cumulus_cldtop_pressure%f(ncol))
    allocate(pstate_cam%cumulus_cldbot_pressure%f(ncol))
    allocate(pstate_cam%rprdsh_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%rprdtot_at_pc_full_level%f(nlev,ncol))

!------------- for Double Plume scheme, LiXH ------------>
    allocate(pstate_cam%pcape_before_dycore%f(ncol))
    allocate(pstate_cam%pcape_previous%f(ncol))
!<------------ for Double Plume scheme, LiXH -------------

    pstate_cam%scalar_precc_sh_surface%pos       = 0
    pstate_cam%scalar_snowc_surface%pos          = 0
    pstate_cam%sh_icwmr_at_pc_full_level%pos     = 0
    pstate_cam%cumulus_cldtop%pos                = 0
    pstate_cam%cumulus_cldbot%pos                = 0
    pstate_cam%cumulus_cldtop2%pos               = 0
    pstate_cam%cumulus_cldbot2%pos               = 0
    pstate_cam%cumulus_cldtop_pressure%pos       = 0
    pstate_cam%cumulus_cldbot_pressure%pos       = 0
    pstate_cam%rprdsh_at_pc_full_level%pos       = 0
    pstate_cam%rprdtot_at_pc_full_level%pos      = 0
    pstate_cam%pcape_before_dycore%pos           = 0
    pstate_cam%pcape_previous%pos                = 0

    pstate_cam%scalar_precc_sh_surface%f         = 0._r8
    pstate_cam%scalar_snowc_surface%f            = 0._r8
    pstate_cam%sh_icwmr_at_pc_full_level%f       = 0._r8
    pstate_cam%cumulus_cldtop%f                  = nlev*1._r8
    pstate_cam%cumulus_cldbot%f                  = 1._r8
    pstate_cam%cumulus_cldtop2%f                 = nlev*1._r8
    pstate_cam%cumulus_cldbot2%f                 = 1._r8
    pstate_cam%cumulus_cldtop_pressure%f         = 0._r8
    pstate_cam%cumulus_cldbot_pressure%f         = 0._r8
    pstate_cam%rprdsh_at_pc_full_level%f         = 0._r8
    pstate_cam%rprdtot_at_pc_full_level%f        = 0._r8
    pstate_cam%pcape_before_dycore%f             = 0._r8
    pstate_cam%pcape_previous%f                  = 0._r8

    !for Lin Macro ------------------------>
    allocate(pstate_cam%sh_qtu_shal_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%sh_thlu_shal_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%sh_umf_shal_at_pc_face_level%f(nlevp,ncol))
    pstate_cam%sh_qtu_shal_at_pc_face_level%pos  = 0
    pstate_cam%sh_thlu_shal_at_pc_face_level%pos = 0
    pstate_cam%sh_umf_shal_at_pc_face_level%pos  = 0
    pstate_cam%sh_qtu_shal_at_pc_face_level%f    = 0._r8
    pstate_cam%sh_thlu_shal_at_pc_face_level%f   = 0._r8
    pstate_cam%sh_umf_shal_at_pc_face_level%f    = 0._r8
    !for Lin Macro <------------------------


    if( shallow_scheme .eq. 'UW' .or. shallow_scheme .eq. 'double_plume') then
        allocate(pstate_cam%sh_shfrc_at_pc_full_level%f(nlev,ncol))
        allocate(pstate_cam%cumulus_cush%f(ptimelevels,ncol))
        pstate_cam%sh_shfrc_at_pc_full_level%pos     = 0
        pstate_cam%cumulus_cush%pos                  = 0
        pstate_cam%sh_shfrc_at_pc_full_level%f       = 0._r8
        pstate_cam%cumulus_cush%f                    = 1000._r8
    end if

    end subroutine convect_shallow_init


    subroutine end_of_convect_shallow

    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam

     if(allocated(flxprec)) deallocate(flxprec)
     if(allocated(flxsnow)) deallocate(flxsnow)
     if(allocated(evapcsh)) deallocate(evapcsh)
     if(allocated(snow_sh)) deallocate(snow_sh)

!------------- for Double Plume scheme, LiXH ------------>
     if(allocated(pstate_cam%pcape_before_dycore%f))              &   
       deallocate(pstate_cam%pcape_before_dycore%f)

     if(allocated(pstate_cam%pcape_previous%f))                   &   
       deallocate(pstate_cam%pcape_previous%f)
!<------------ for Double Plume scheme, LiXH -------------

     if(allocated(pstate_cam%cumulus_cush%f))                     &
       deallocate(pstate_cam%cumulus_cush%f)

     if(allocated(pstate_cam%sh_shfrc_at_pc_full_level%f))        &
       deallocate(pstate_cam%sh_shfrc_at_pc_full_level%f)

     if(allocated(pstate_cam%cumulus_cldtop%f))                   &
       deallocate(pstate_cam%cumulus_cldtop%f)

     if(allocated(pstate_cam%cumulus_cldbot%f))                   &
       deallocate(pstate_cam%cumulus_cldbot%f)

     if(allocated(pstate_cam%cumulus_cldtop2%f))                  &
       deallocate(pstate_cam%cumulus_cldtop2%f)

     if(allocated(pstate_cam%cumulus_cldbot2%f))                  &
       deallocate(pstate_cam%cumulus_cldbot2%f)


     if(allocated(pstate_cam%cumulus_cldtop_pressure%f))          &
       deallocate(pstate_cam%cumulus_cldtop_pressure%f)

     if(allocated(pstate_cam%cumulus_cldbot_pressure%f))          &
       deallocate(pstate_cam%cumulus_cldbot_pressure%f)

     if(allocated(pstate_cam%sh_icwmr_at_pc_full_level%f))        &
       deallocate(pstate_cam%sh_icwmr_at_pc_full_level%f)

     if(allocated(pstate_cam%scalar_precc_sh_surface%f))          &
       deallocate(pstate_cam%scalar_precc_sh_surface%f)

     if(allocated(pstate_cam%scalar_snowc_surface%f))             &
       deallocate(pstate_cam%scalar_snowc_surface%f)

     if(allocated(pstate_cam%rprdsh_at_pc_full_level%f))          &
       deallocate(pstate_cam%rprdsh_at_pc_full_level%f)

     if(allocated(pstate_cam%rprdtot_at_pc_full_level%f))         &
       deallocate(pstate_cam%rprdtot_at_pc_full_level%f)

     !for Lin Macro ------------------------>
     if(allocated(pstate_cam%sh_qtu_shal_at_pc_face_level%f))     &
       deallocate(pstate_cam%sh_qtu_shal_at_pc_face_level%f)
     if(allocated(pstate_cam%sh_thlu_shal_at_pc_face_level%f))    &
       deallocate(pstate_cam%sh_thlu_shal_at_pc_face_level%f)  
     if(allocated(pstate_cam%sh_umf_shal_at_pc_face_level%f))     &
       deallocate(pstate_cam%sh_umf_shal_at_pc_face_level%f)  
     !for Lin Macro <------------------------
 
    end subroutine end_of_convect_shallow


! Purpose: Return true if cloud fraction should use shallow convection
!          convect_shallow_use_shfrc() = .true.   for     shallow_scheme = 'UW'
!          convect_shallow_use_shfrc() = .false.  for     all other schemes
    function convect_shallow_use_shfrc()
    
    implicit none
    logical :: convect_shallow_use_shfrc     ! Return value

    if ( shallow_scheme .eq. 'UW' .or. shallow_scheme .eq. 'double_plume') then
         convect_shallow_use_shfrc = .true.
    else
         convect_shallow_use_shfrc = .false.
    endif

    return

    end function convect_shallow_use_shfrc

    !------------------------------>
    ! LiXH add for dCAPE trigger function of Xie et al. (2018)
    subroutine trigger_function_xie(ncol)

    use grist_double_plume,                 only: cape_before_dycore
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
! io
    integer(i4), intent(in)     :: ncol

    if(shallow_scheme .eq. 'double_plume')then
    call cape_before_dycore( nlev, ncol,    &
                             pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol), &
                             pstate%temp_at_pc_full_level%f(:,1:ncol),          &
                             pstate%pressure_at_pc_full_level%f(:,1:ncol),      &
                             pstate%z_at_pc_full_level%f(:,1:ncol),             &
                             pstate%pressure_at_pc_face_level%f(:,1:ncol),      &
                             pstate%z_at_pc_face_level%f(:,1:ncol),             &
                             pstate_cam%pbl_tpert_at_pc_surface%f(1:ncol),      &
                             pstate_cam%pbl_pblh_at_pc_surface%f(1:ncol),       &
                             pstate_cam%pcape_before_dycore%f(1:ncol) )
    end if

    end subroutine trigger_function_xie
    !<------------------------------

    subroutine convect_shallow_tend(ncol, nstep, dt, cmfmc2, qc, qc2, rliq, rliq2)

    use grist_physics_data_structure,       only: pstate, phy_tracer_info 
    use grist_PhysC_data_structure,         only: pstate_cam,  &
                                                  ptend_shallow_convection
    use grist_uwshcu,                       only: compute_uwshcu_inv
!------------- for Double Plume scheme, LiXH ------------>
    use grist_double_plume,                 only: compute_double_plume_inv
!<------------ for Double Plume scheme, LiXH -------------
    use grist_physics_update,               only: old_time_level
    use grist_mpi

! io
    integer(i4), intent(in)     :: ncol
    integer(i4), intent(in)     :: nstep
    real(r8),    intent(in)     :: dt
    real(r8),    intent(inout)  :: qc(nlev, ncol)           ! dq/dt due to export of cloud water into environment by shallow and deep convection [kg/kg/s ]
    real(r8),    intent(inout)  :: rliq(ncol)               ! Vertical integral of qc [ m/s ]

    real(r8),    intent(out)    :: cmfmc2(nlevp, ncol)      ! Updraft mass flux by shallow convection [ kg/s/m2 ]
    real(r8),    intent(out)    :: qc2(nlev, ncol)          ! Same as qc but only from shallow convection scheme
    real(r8),    intent(out)    :: rliq2(ncol)              ! Vertically-integrated reserved cloud condensate [ m/s ]
! local
    integer(i4)  :: i, k, m
    integer(i4)  :: ixcldice, ixcldliq                      ! Constituent indices for cloud liquid and ice water.
    real(r8)     :: tke(nlevp, ncol)    
    real(r8)     :: pblh(ncol)
    real(r8)     :: cnt2(ncol)                              ! Top level of shallow convective activity
    real(r8)     :: cnb2(ncol)                              ! Bottom level of convective activity
    real(r8)     :: slflx(nlevp, ncol)                      ! Shallow convective liquid water static energy flux
    real(r8)     :: qtflx(nlevp, ncol)                      ! Shallow convective total water flux
    real(r8)     :: cmfdqs(nlev, ncol)                      ! Shallow convective snow production
    real(r8)     :: cbmf(ncol)                              ! Shallow cloud base mass flux [ kg/s/m2 ]
    real(r8)     :: freqsh(ncol)                            ! Frequency of shallow convection occurence
    real(r8)     :: cmfsl(nlevp, ncol)                      ! Convective flux of liquid water static energy)
    real(r8)     :: cmflq(nlevp, ncol)                      ! Convective flux of total water in energy unit)
    real(r8)     :: iccmr_UW(nlev, ncol)                    ! In-cloud Cumulus LWC+IWC [ kg/m2 ]
    real(r8)     :: icwmr_UW(nlev, ncol)                    ! In-cloud Cumulus LWC     [ kg/m2 ]
    real(r8)     :: icimr_UW(nlev, ncol)                    ! In-cloud Cumulus IWC     [ kg/m2 ] 

    real(r8)     :: ptend_qv(nlev, ncol)
    real(r8)     :: ptend_liq(nlev, ncol)
    real(r8)     :: ptend_ice(nlev, ncol)
    real(r8)     :: ptend_s(nlev, ncol)
    real(r8)     :: ptend_u(nlev, ncol)
    real(r8)     :: ptend_v(nlev, ncol)
    real(r8)     :: ptend_other_tracer(ntracer, nlev, ncol)

    do m = 1, ntracer
        if(trim(phy_tracer_info(m)%longname) .eq. 'cloud_liquid') ixcldliq = m
        if(trim(phy_tracer_info(m)%longname) .eq. 'cloud_ice')    ixcldice = m
    end do

    select case (shallow_scheme)

    case('off')

    cmfmc2    = 0._r8
    ptend_u   = 0._r8
    ptend_v   = 0._r8
    ptend_qv  = 0._r8
    ptend_liq = 0._r8
    ptend_ice = 0._r8
    ptend_other_tracer = 0._r8
    ptend_s   = 0._r8
    pstate_cam%rprdsh_at_pc_full_level%f(:,1:ncol)   = 0._r8
    cmfdqs    = 0._r8
    pstate_cam%scalar_precc_sh_surface%f(1:ncol)     = 0._r8 
    slflx     = 0._r8
    qtflx     = 0._r8
    pstate_cam%sh_icwmr_at_pc_full_level%f(:,1:ncol) = 0._r8
    rliq2     = 0._r8
    qc2       = 0._r8
    cmfsl     = 0._r8
    cmflq     = 0._r8
    cnt2      = nlev
    cnb2      = 1._r8
    evapcsh   = 0._r8
    snow_sh   = 0._r8

    case('UW')

    tke(:,1:ncol)  = pstate_cam%pbl_tke_at_pc_face_level%f(:,1:ncol)
    pblh(1:ncol)   = pstate_cam%pbl_pblh_at_pc_surface%f(1:ncol)

    call compute_uwshcu_inv(nlev      , ncol       , ntracer      , dt       ,                                                      &
                            pstate%pressure_at_pc_face_level%f(:,1:ncol)     , pstate%z_at_pc_face_level%f(:,1:ncol) ,              &
                            pstate%pressure_at_pc_full_level%f(:,1:ncol)     , pstate%z_at_pc_full_level%f(:,1:ncol) ,              &
                            pstate%delp_at_pc_full_level%f(:,1:ncol)         ,                                                      &
                            pstate%u_wind_at_pc_full_level%f(:,1:ncol)       , pstate%v_wind_at_pc_full_level%f(:,1:ncol) ,         &
                            pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol), pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,:,1:ncol),   &
                            pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,:,1:ncol), pstate%temp_at_pc_full_level%f(:,1:ncol)         ,   &
                            pstate%static_energy_at_pc_full_level%f(:,1:ncol), pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol), tke  ,   &
                            pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)        ,                                  &
                            pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,1:ncol)     ,                                  &
                            pblh      , pstate_cam%cumulus_cush%f(old_time_level,1:ncol)             ,                              &
                            cmfmc2    , slflx      , qtflx        , flxprec  , flxsnow ,                                            &
                            ptend_qv  , ptend_liq  , ptend_ice    , ptend_s  , ptend_u , ptend_v , ptend_other_tracer,              &
                            pstate_cam%rprdsh_at_pc_full_level%f(:,1:ncol)       , cmfdqs  , pstate_cam%scalar_precc_sh_surface%f(1:ncol),  &
                            snow_sh   , evapcsh , pstate_cam%sh_shfrc_at_pc_full_level%f(:,1:ncol)   , iccmr_UW, icwmr_UW, icimr_UW  ,  &
                            cbmf      , qc2     , rliq2    , cnt2 , cnb2     ,                                                          &
                            !for Lin Macro ------------------------>
                            pstate_cam%sh_qtu_shal_at_pc_face_level%f        ,                                                          &
                            pstate_cam%sh_thlu_shal_at_pc_face_level%f )

                            pstate_cam%sh_umf_shal_at_pc_face_level%f = cmfmc2
                            !for Lin Macro <------------------------

    pstate_cam%sh_icwmr_at_pc_full_level%f(:,1:ncol) = iccmr_UW(:,1:ncol)
!    ! Add shallow reserved cloud condensate to deep reserved cloud condensate
!    ! qc [ kg/kg/s] , rliq [ m/s ] 
!    qc(:nlev,:ncol)    = qc(:nlev,:ncol) + qc2(:nlev,:ncol)
!    rliq(:ncol)        = rliq(:ncol) + rliq2(:ncol)

    case('double_plume')
!------------- for Double Plume scheme, LiXH ------------>

    tke(:,1:ncol)  = pstate_cam%pbl_tke_at_pc_face_level%f(:,1:ncol)
    pblh(1:ncol)   = pstate_cam%pbl_pblh_at_pc_surface%f(1:ncol)


    !Double Plume contains deep and shallow convection
    call compute_double_plume_inv(nlev      , ncol       , ntracer      , dt , nstep  ,                                                     &
                            pstate%pressure_at_pc_face_level%f(:,1:ncol)     , pstate%z_at_pc_face_level%f(:,1:ncol) ,                      &
                            pstate%pressure_at_pc_full_level%f(:,1:ncol)     , pstate%z_at_pc_full_level%f(:,1:ncol) ,                      &
                            pstate%delp_at_pc_full_level%f(:,1:ncol)         ,                                                              &
                            pstate%u_wind_at_pc_full_level%f(:,1:ncol)       , pstate%v_wind_at_pc_full_level%f(:,1:ncol) ,                 &
                            pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol), pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,:,1:ncol),    &
                            pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,:,1:ncol), pstate%temp_at_pc_full_level%f(:,1:ncol)         ,    &
                            pstate%static_energy_at_pc_full_level%f(:,1:ncol), pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol),           &
                            pstate_cam%pbl_tpert_at_pc_surface%f(1:ncol)     , pstate%landfrac_at_pc_surface%f(1:ncol) , tke  ,             &
                            pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)        ,                                      &
                            pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,1:ncol)     ,                                      &
                            pblh      , pstate_cam%cumulus_cush%f(old_time_level,1:ncol)             , pstate_cam%pcape_previous%f(1:ncol), &
                            !---------LiXH add for (Xie) dCAPE trigger--------->
                            pstate_cam%pcape_before_dycore%f(1:ncol),                                                                       &
                            !<--------LiXH add for (Xie) dCAPE trigger----------
                            cmfmc2  , slflx   , qtflx   , flxprec , flxsnow  ,                                                              &
                            ptend_qv  , ptend_liq  , ptend_ice    , ptend_s  , ptend_u , ptend_v , ptend_other_tracer,                      &
                            pstate_cam%rprdsh_at_pc_full_level%f(:,1:ncol)   , cmfdqs  , pstate_cam%scalar_precc_sh_surface%f(1:ncol),      &
                            snow_sh   , evapcsh , pstate_cam%sh_shfrc_at_pc_full_level%f(:,1:ncol)   , iccmr_UW , icwmr_UW , icimr_UW,      &
                            cbmf      , qc2     , rliq2    , cnt2 , cnb2    ,                                                               &
                            !for Lin Macro ------------------------>
                            pstate_cam%sh_qtu_shal_at_pc_face_level%f        ,                                                              &
                            pstate_cam%sh_thlu_shal_at_pc_face_level%f )

                            pstate_cam%sh_umf_shal_at_pc_face_level%f = cmfmc2
                            !for Lin Macro <------------------------

!<------------ for Double Plume scheme, LiXH -------------

    end select
 
    ! Add shallow reserved cloud condensate to deep reserved cloud condensate
    ! qc [ kg/kg/s] , rliq [ m/s ] 
    qc(:nlev,:ncol)    = qc(:nlev,:ncol) + qc2(:nlev,:ncol)
    rliq(:ncol)        = rliq(:ncol) + rliq2(:ncol)

    pstate_cam%scalar_snowc_surface%f(1:ncol)        = pstate_cam%scalar_snowc_surface%f(1:ncol) + snow_sh(1:ncol)
!    pstate_cam%sh_icwmr_at_pc_full_level%f(:,1:ncol) = iccmr_UW(:,1:ncol)
    pstate_cam%rprdsh_at_pc_full_level%f(:,1:ncol)   = pstate_cam%rprdsh_at_pc_full_level%f(:,1:ncol) + cmfdqs(:,1:ncol)

    ! Convective fluxes of 'sl' and 'qt' in energy unit
    cmfsl(:nlevp,:ncol) = slflx(:nlevp,:ncol)
    cmflq(:nlevp,:ncol) = qtflx(:nlevp,:ncol)*latvap

    
    ! Modification : I should check whether below computation of freqsh is correct.
    freqsh(:) = 0._r8
    do i = 1, ncol
        if( maxval(cmfmc2(:nlev,i)) <= 0._r8 ) then
            freqsh(i) = 1._r8
        end if
    end do

    ! Merge shallow convection output with prior results from deep convection scheme 

    ! Combine cumulus updraft mass flux : 'cmfmc2'(shallow) + 'cmfmc'(deep)
    ! cmfmc(:nlev,:ncol) = cmfmc(:nlev,:ncol) + cmfmc2(:nlev,:ncol)
    pstate_cam%updraft_mass_flux%f(:nlev,:ncol) = pstate_cam%updraft_mass_flux%f(:nlev,:ncol) + cmfmc2(:nlev,:ncol)


    ! 'cnt2' & 'cnb2' are from shallow, 'cumulus_cldtop' & 'cumulus_cldbot' are from deep 
    ! 'cnt2' & 'cnb2' are the interface indices of cloud top & base: 
    !        cnt2 = float(kpen)
    !        cnb2 = float(krel - 1)
    ! Note that indices decreases with height.
    do i = 1, ncol
       if( cnt2(i) < pstate_cam%cumulus_cldtop%f(i)) pstate_cam%cumulus_cldtop%f(i) = cnt2(i)
       if( cnb2(i) > pstate_cam%cumulus_cldbot%f(i)) pstate_cam%cumulus_cldbot%f(i) = cnb2(i)
       pstate_cam%cumulus_cldtop_pressure%f(i) = pstate%pressure_at_pc_full_level%f(int(pstate_cam%cumulus_cldtop%f(i)),i)
       pstate_cam%cumulus_cldbot_pressure%f(i) = pstate%pressure_at_pc_full_level%f(int(pstate_cam%cumulus_cldbot%f(i)),i)     
    end do
 
    !--------------------for Lin Macrop, LiXH------------------->
    pstate_cam%cumulus_cldtop2%f(1:ncol) = cnt2(1:ncol)
    pstate_cam%cumulus_cldbot2%f(1:ncol) = cnb2(1:ncol)
    !<-------------------for Lin Macrop, LiXH--------------------

    ! This quantity was previously known as CMFDQR.
    ! Now CMFDQR is the shallow rain production only.
    pstate_cam%rprdtot_at_pc_full_level%f(:,1:ncol) =  pstate_cam%rprddp_at_pc_full_level%f(:,1:ncol) +  &
                                                       pstate_cam%rprdsh_at_pc_full_level%f(:,1:ncol)

    ! re-initialize ptend
    ptend_shallow_convection%tend_u%f  = 0._r8
    ptend_shallow_convection%tend_v%f  = 0._r8
    ptend_shallow_convection%tend_s%f  = 0._r8
    ptend_shallow_convection%tend_q%f  = 0._r8
 
    ptend_shallow_convection%tend_u%f(:,1:ncol)            = ptend_u(:,1:ncol)
    ptend_shallow_convection%tend_v%f(:,1:ncol)            = ptend_v(:,1:ncol)
    ptend_shallow_convection%tend_s%f(:,1:ncol)            = ptend_s(:,1:ncol)
    ptend_shallow_convection%tend_q%f(1,:,1:ncol)          = ptend_qv(:,1:ncol)
    ptend_shallow_convection%tend_q%f(ixcldliq,:,1:ncol)   = ptend_liq(:,1:ncol)
    ptend_shallow_convection%tend_q%f(ixcldice,:,1:ncol)   = ptend_ice(:,1:ncol)
    do m = 4, ntracer
        ptend_shallow_convection%tend_q%f(m,:,1:ncol)      = ptend_other_tracer(m,:,1:ncol)
    end do

    
    ! outfld:
    ! cmfsl, cmflq, freqsh

    end subroutine convect_shallow_tend
 end module grist_shallow_convection
