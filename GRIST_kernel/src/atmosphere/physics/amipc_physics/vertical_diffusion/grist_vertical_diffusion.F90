!===================================================================================================
!
!  Created by LiXiaohan on 19/06/01, adopted from CAM5
!
!  to compute vertical diffusion of momentum, moisture, and trace constituents
!    1. turbulent mountain stress
!    2. eddy diffusivities
!    3. molecular diffusivities
!    4. a implicit diffusion solver
!
!  Calling sequence:
!     vertical_diffusion_init:
!       init_molec_diff        Initializes molecular diffusivity module
!       init_eddy_diff         Initializes eddy diffusivity module (includes PBL)
!       init_tms               Initializes turbulent mountain stress module
!       init_vdiff             Initializes diffusion solver module
!     vertical_diffusion_ts_init:
!       Time step initialization (only used for upper boundary condition)
!     vertical_diffusion_tend:
!       compute_tms            Computes turbulent mountain stresses
!       compute_eddy_diff      Computes eddy diffusivities and countergradient terms
!       compute_vdiff          Solves vertical diffusion equations, including molecular diffusivities
!
!  Reference:
!     Bretherton, C., and S. Park, 2009: A New Moist Turbulence Parameterization in the Community 
!     Atmosphere Model. J. Clim., 22, 3422–3448.
!      
!===================================================================================================

 module grist_vertical_diffusion

!    use grist_lib
    use grist_mpi
    use grist_handle_error,                 only: endrun
    use grist_nml_module,                   only: ntracer, nlev, nlevp, sub_physpkg
    use grist_constants,                    only: r8, i4, rdry, mwdry, avogad,          &
                                                  gravity, cp, boltz, zvir,             &
                                                  latvap, latice, karman
    use grist_diffusion_solver,             only: vdiff_selector

!    use grist_fileio_0d_module_gcm,         only: wrap_bcast_0d


    implicit none
    private
    public          :: read_nml_vertical_diffusion,     &
                       vertical_diffusion_init,         &
                       vertical_diffusion_tend,         &
                       end_of_vertical_diffusion
               
! Private:
    integer, parameter   :: nturb    = 5        ! Number of iterations for solution
    logical, parameter   :: wstarent = .true.   ! Use wstar (.true.) or TKE (.false.) entrainment closure
    logical         :: kvinit                   ! Tell compute_eddy_diff/ caleddy
    logical         :: do_molec_diff            ! Switch for molecular diffusion
    type(vdiff_selector) :: fieldlist_wet       ! Logical switches for moist mixing ratio diffusion
    integer(i4)     :: ntop                     ! Top interface level to which vertical diffusion is applied ( = 1 ).
    integer(i4)     :: nbot                     ! Bottom interface level to which vertical diffusion is applied ( = nlev ).
    real(r8)        :: kv_top_pressure          ! Pressure defining the bottom of the upper atmosphere for kvh scaling (Pa)
    real(r8)        :: kv_top_scale             ! Eddy diffusivity scale factor for upper atmosphere
    real(r8)        :: kv_freetrop_scale        ! Eddy diffusivity scale factor for the free troposphere
    real(r8)        :: eddy_lbulk_max           ! Maximum master length for diag_TKE
    real(r8)        :: eddy_leng_max            ! Maximum dissipation length for diag_TKE
    real(r8)        :: eddy_max_bot_pressure    ! Bottom pressure level (hPa) for eddy_leng_max
    logical         :: diff_cnsrv_mass_check    ! do mass conservation check
    logical         :: do_iss                   ! switch for implicit turbulent surface stress
    logical         :: do_tms                   ! switch for turbulent mountain stress
    real(r8)        :: tms_orocnst
    real(r8)        :: tms_z0fac

! if the variables below are used by other modules (pbuf in CAM), change to pstate%~, Lixh
    integer(i4), allocatable     :: turbtype(:,:)       ! 'turbtype'
    real(r8),    allocatable     :: smaw(:,:)           ! 'smaw'
    
! other pbuf includes:
!

contains
    subroutine read_nml_vertical_diffusion(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer :: unitn, ierr

    namelist /vert_diff_nl/ kv_top_pressure,        &
                            kv_top_scale,           &
                            kv_freetrop_scale,      &
                            eddy_lbulk_max,         &
                            eddy_leng_max,          &
                            eddy_max_bot_pressure,  &
                            diff_cnsrv_mass_check,  &
                            do_iss,                 &
                            do_tms,                 &
                            tms_orocnst,            &
                            tms_z0fac

!    if(mpi_rank() == 0) then
        unitn = 111
        open(unitn, file=trim(nlfile), status='old')
        read(unitn, nml=vert_diff_nl, iostat=ierr)
        if(ierr /= 0) call endrun("error reading vertical_diffusion namelist")
        close(unitn)
!    end if

!    call wrap_bcast_0d(grist_comm(), 0, kv_top_pressure)
!    call wrap_bcast_0d(grist_comm(), 0, kv_top_scale)
!    call wrap_bcast_0d(grist_comm(), 0, kv_freetrop_scale)
!    call wrap_bcast_0d(grist_comm(), 0, eddy_lbulk_max)
!    call wrap_bcast_0d(grist_comm(), 0, eddy_leng_max)
!    call wrap_bcast_0d(grist_comm(), 0, eddy_max_bot_pressure)
!    call wrap_bcast_0d(grist_comm(), 0, diff_cnsrv_mass_check)
!    call wrap_bcast_0d(grist_comm(), 0, do_iss)

    end  subroutine read_nml_vertical_diffusion


    subroutine vertical_diffusion_init(ncol)
    use grist_constants,                    only: p00
    use grist_hpe_constants,                only: eta_full
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_eddy_diff,                    only: init_eddy_diff
    use grist_molec_diff,                   only: init_molec_diff
    use grist_turb_mon_stress,              only: init_tms
    use grist_diffusion_solver,             only: init_vdiff, new_fieldlist_vdiff, vdiff_select
! io
    integer(i4),    intent(in) :: ncol
! local
    integer(i4) :: i
    integer(i4) :: ntop_eddy    ! top interface level to which eddy vertical diffusion is applied (= 1)
    integer(i4) :: nbot_eddy    ! bottom interface level to which eddy vertical diffusion is applied (=nlev)
    integer(i4) :: ntop_molec   ! top interface level to which molecular vertical diffusion is applied (=1)
    integer(i4) :: nbot_molec   ! bottom interface level to which molecular vertical diffusion is applied
    real(r8)             :: lev_ref(nlev)
    real(r8), parameter  :: do_molec_pres = 0.1_r8      ! If top of model is above this pressure,
                                                        ! turn on molecular diffusion. (Pa)

    allocate(turbtype(nlevp, ncol));     turbtype = 0
    allocate(smaw(nlevp, ncol));         smaw     = 0._r8

    ! Initialize pstate_cam:
    allocate(pstate_cam%pbl_tke_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%pbl_kvm_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%pbl_kvh_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%pbl_kvt_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%pbl_pblh_at_pc_surface%f(ncol))
    allocate(pstate_cam%pbl_tauresx_at_pc_surface%f(ncol))
    allocate(pstate_cam%pbl_tauresy_at_pc_surface%f(ncol))
    if(.not.allocated(pstate_cam%pbl_tpert_at_pc_surface%f)) allocate(pstate_cam%pbl_tpert_at_pc_surface%f(ncol))
    if(.not.allocated(pstate_cam%pbl_qpert_at_pc_surface%f)) allocate(pstate_cam%pbl_qpert_at_pc_surface%f(ncol))

    pstate_cam%pbl_tke_at_pc_face_level%pos      = 0
    pstate_cam%pbl_kvm_at_pc_face_level%pos      = 0
    pstate_cam%pbl_kvh_at_pc_face_level%pos      = 0
    pstate_cam%pbl_kvt_at_pc_face_level%pos      = 0
    pstate_cam%pbl_pblh_at_pc_surface%pos        = 0
    pstate_cam%pbl_tauresx_at_pc_surface%pos     = 0
    pstate_cam%pbl_tauresy_at_pc_surface%pos     = 0
    pstate_cam%pbl_tpert_at_pc_surface%pos       = 0
    pstate_cam%pbl_qpert_at_pc_surface%pos       = 0

    pstate_cam%pbl_tke_at_pc_face_level%f        = 0.01_r8
    pstate_cam%pbl_kvm_at_pc_face_level%f        = 0._r8
    pstate_cam%pbl_kvh_at_pc_face_level%f        = 0._r8
    pstate_cam%pbl_kvt_at_pc_face_level%f        = 0._r8
    pstate_cam%pbl_pblh_at_pc_surface%f          = 0._r8
    pstate_cam%pbl_tauresx_at_pc_surface%f       = 0._r8
    pstate_cam%pbl_tauresy_at_pc_surface%f       = 0._r8
    pstate_cam%pbl_tpert_at_pc_surface%f         = 0._r8
    pstate_cam%pbl_qpert_at_pc_surface%f         = 0._r8
    ! for Lin Macro ------------------------>
    allocate(pstate_cam%pbl_lengi_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%pbl_shi_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%pbl_smi_at_pc_face_level%f(nlevp,ncol))
    allocate(pstate_cam%pbl_wstarPBL%f(ncol))
    pstate_cam%pbl_lengi_at_pc_face_level%pos    = 0
    pstate_cam%pbl_shi_at_pc_face_level%pos      = 0
    pstate_cam%pbl_smi_at_pc_face_level%pos      = 0
    pstate_cam%pbl_wstarPBL%pos                  = 0
    pstate_cam%pbl_lengi_at_pc_face_level%f      = 0._r8
    pstate_cam%pbl_shi_at_pc_face_level%f        = 0._r8
    pstate_cam%pbl_smi_at_pc_face_level%f        = 0._r8
    pstate_cam%pbl_wstarPBL%f                    = 0._r8
    ! for Lin Macro <------------------------



    ! Initialize molecular diffusion and get top and bottom molecular diffusion limits:
    ! Note that computing molecular diffusivities is a trivial expense, but constituent
    ! diffusivities depend on their molecular weights. Decomposing the diffusion matric
    ! for each constituent is a needless expense unless the diffusivity is significant.

    kvinit = .true.

    lev_ref(:) = eta_full(:)*p00
    do_molec_diff = (lev_ref(1) < do_molec_pres)
    if(do_molec_diff) then
        call init_molec_diff( ntracer, rdry, mwdry, avogad, gravity,            &
                              cp, boltz, lev_ref, ntop_molec, nbot_molec )
    else
        ntop_molec = 1
        nbot_molec = 0
    end if

    ! Initialize eddy diffusivity module
    ntop_eddy  = 1
    nbot_eddy  = nlev

!-----------------------------------Lixh close this part---------------------------------->
!    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
!       eddy_top_loop: do k = 2, pver
!          if (pref_mid(k) .gt. ntop_eddy_pres) then
!             ntop_eddy  = k-1
!             exit eddy_top_loop
!          endif
!       end do eddy_top_loop
!    end if
!<----------------------------------Lixh close this part-----------------------------------


    call init_eddy_diff( nlev, gravity, cp, rdry, zvir, latvap, latice,         &
                         ntop_eddy, nbot_eddy, karman, eddy_lbulk_max,          &
                         eddy_leng_max, eddy_max_bot_pressure, lev_ref)

    ! The vertical diffusion solver must operate over the full range of molecular and eddy diffusion
    ntop = min(ntop_molec,ntop_eddy)
    nbot = max(nbot_molec,nbot_eddy)
    
    ! Initialize turbulent mountain stress module
    if( do_tms ) then
       call init_tms( tms_orocnst, tms_z0fac, karman, gravity, rdry )
    endif

    ! Initialize diffusion solver module
    call init_vdiff( rdry, gravity, do_iss )
 
    fieldlist_wet = new_fieldlist_vdiff( ntracer )
    if(vdiff_select( fieldlist_wet, 'u' ) .ne. '' ) call endrun(vdiff_select( fieldlist_wet, 'u' ))
    if(vdiff_select( fieldlist_wet, 'v' ) .ne. '' ) call endrun(vdiff_select( fieldlist_wet, 'v' ))
    if(vdiff_select( fieldlist_wet, 's' ) .ne. '' ) call endrun(vdiff_select( fieldlist_wet, 's' ))

    IF(trim(sub_physpkg).eq.'DCMIP2016-TC')then
    do i = 1, 1
        if(vdiff_select( fieldlist_wet, 'q', i ) .ne. '') call endrun(vdiff_select( fieldlist_wet, 'q', i ))
    end do
    ELSE
    do i = 1, 3
        if(vdiff_select( fieldlist_wet, 'q', i ) .ne. '') call endrun(vdiff_select( fieldlist_wet, 'q', i ))
    end do
    END IF

    end subroutine vertical_diffusion_init

    subroutine end_of_vertical_diffusion

    use grist_PhysC_data_structure,       only: pstate_cam

     if(allocated(turbtype)) deallocate(turbtype)
     if(allocated(smaw))     deallocate(smaw)

     if(allocated(pstate_cam%pbl_tke_at_pc_face_level%f))         &
       deallocate(pstate_cam%pbl_tke_at_pc_face_level%f)

     if(allocated(pstate_cam%pbl_kvm_at_pc_face_level%f))         &
       deallocate(pstate_cam%pbl_kvm_at_pc_face_level%f)

     if(allocated(pstate_cam%pbl_kvh_at_pc_face_level%f))         &
       deallocate(pstate_cam%pbl_kvh_at_pc_face_level%f)

     if(allocated(pstate_cam%pbl_kvt_at_pc_face_level%f))         &
       deallocate(pstate_cam%pbl_kvt_at_pc_face_level%f)

     if(allocated(pstate_cam%pbl_pblh_at_pc_surface%f))           &
       deallocate(pstate_cam%pbl_pblh_at_pc_surface%f)

     if(allocated(pstate_cam%pbl_tauresx_at_pc_surface%f))        &
       deallocate(pstate_cam%pbl_tauresx_at_pc_surface%f)

     if(allocated(pstate_cam%pbl_tauresy_at_pc_surface%f))        &
       deallocate(pstate_cam%pbl_tauresy_at_pc_surface%f)

     if(allocated(pstate_cam%pbl_tpert_at_pc_surface%f))          &
       deallocate(pstate_cam%pbl_tpert_at_pc_surface%f)

     if(allocated(pstate_cam%pbl_qpert_at_pc_surface%f))          &
       deallocate(pstate_cam%pbl_qpert_at_pc_surface%f)

     ! for Lin Macro ------------------------>
     if(allocated(pstate_cam%pbl_lengi_at_pc_face_level%f))       &
       deallocate(pstate_cam%pbl_lengi_at_pc_face_level%f)  
     if(allocated(pstate_cam%pbl_shi_at_pc_face_level%f))         &
       deallocate(pstate_cam%pbl_shi_at_pc_face_level%f)
     if(allocated(pstate_cam%pbl_smi_at_pc_face_level%f))         &
       deallocate(pstate_cam%pbl_smi_at_pc_face_level%f)
     if(allocated(pstate_cam%pbl_wstarPBL%f))                     &
       deallocate(pstate_cam%pbl_wstarPBL%f)
     ! for Lin Macro <------------------------


    end subroutine end_of_vertical_diffusion

    subroutine vertical_diffusion_tend(ncol, dt)
    use grist_eddy_diff,                    only: compute_eddy_diff
    use grist_physics_data_structure,       only: pstate, phy_tracer_info
    use grist_PhysC_data_structure,         only: pstate_cam ,    &
                                                  ptend_vertical_diffusion
    use grist_pbl_utils,                    only: virtem, calc_obklen
    use grist_diffusion_solver,             only: compute_vdiff, any
    use grist_molec_diff,                   only: compute_molec_diff, vd_lu_qdecomp
    use grist_physics_update,               only: old_time_level
    use grist_turb_mon_stress,              only: compute_tms
 
! io
    integer(i4), intent(in)     :: ncol
    real(r8),    intent(in)     :: dt
! local
    integer(i4) ::  i, k
    real(r8) :: rdt
    real(r8) :: dtk(nlev, ncol)                    ! T tendency from KE dissipation
    real(r8) :: ksrftms(ncol)                      ! Turbulent mountain stress surface drag coefficient [ kg/s/m2 ]
    real(r8) :: tautmsx(ncol)                      ! U component of turbulent mountain stress [ N/m2 ]
    real(r8) :: tautmsy(ncol)                      ! V component of turbulent mountain stress [ N/m2 ]
    real(r8) :: tautotx(ncol)                      ! U component of total surface stress [ N/m2 ]
    real(r8) :: tautoty(ncol)                      ! V component of total surface stress [ N/m2 ]
    real(r8) :: rpdel(nlev, ncol)                  ! 1./pdel where 'pdel' is thickness of the layer [ Pa ]
    real(r8) :: kvq(nlevp, ncol)                   ! Eddy diffusivity for constituents [ m2/s ]
    real(r8) :: kvh(nlevp, ncol)                   ! Eddy diffusivity for heat [ m2/s ]
    real(r8) :: kvm(nlevp, ncol)                   ! Eddy diffusivity for momentum [ m2/s ]
    real(r8) :: kvm_in(nlevp, ncol)                ! kvm from previous timestep [ m2/s ]
    real(r8) :: kvh_in(nlevp, ncol)                ! kvh from previous timestep [ m2/s ]
    real(r8) :: cgh(nlevp, ncol)                   ! Counter-gradient term for heat
    real(r8) :: cgs(nlevp, ncol)                   ! Counter-gradient star  [ cg/flux ]
    real(r8) :: rrho(ncol)                         ! Reciprocal of density at surface
    real(r8) :: ustar(ncol)                        ! Surface friction velocity [ m/s ]
    real(r8) :: obklen(ncol)                       ! Obukhov length
    real(r8) :: bprod(nlevp, ncol)                 ! Buoyancy production of tke [ m2/s3]
    real(r8) :: sprod(nlevp, ncol)                 ! Shear production of tke [ m2/s3 ]
    real(r8) :: sfi(nlevp, ncol)                   ! Saturation fraction at interfaces
    real(r8) :: th(nlev, ncol)                     ! Potential temperature
    real(r8) :: topflx(ncol)                       ! Molecular heat flux at top interface
    real(r8) :: wpert(ncol)                        ! Turbulent velocity excess [ m/s ]
    ! for obklen calculation outside HB
    real(r8) :: thvs(ncol)
    real(r8) :: khfs(ncol)
    real(r8) :: kqfs(ncol)
    real(r8) :: kbfs(ncol)
 
    real(r8) :: ipbl(ncol)
    real(r8) :: kpblh(ncol)
    ! for Lin Macro ------------------------>
    !real(r8) :: wstarPBL(ncol)
    ! for Lin Macro <------------------------

    ! Copy state so we can pass to intent(inout) routines that return
    ! new state instead of a tendency.
    real(r8) :: s_tmp(nlev, ncol)
    real(r8) :: u_tmp(nlev, ncol)
    real(r8) :: v_tmp(nlev, ncol)
    real(r8) :: q_tmp(ntracer, nlev, ncol)

    ! Lixh add, rairi = rdry, cpairv = cp, this part will be modified later
    real(r8)      :: rairi(nlevp, ncol)       ! interface gas constant needed for compute_vdiff
    real(r8)      :: cpairv(nlev, ncol)
    real(r8)      :: qmincg(ntracer)

    rdt = 1._r8/dt
    rpdel(:,1:ncol) = 1./pstate%delp_at_pc_full_level%f(:,1:ncol)
    ! Computation of turbulent mountain stress
    !tautmsx = 0._r8
    !tautmsy = 0._r8
    !ksrftms(:ncol) = 0._r8
    if( do_tms ) then
        call compute_tms( nlev     , ncol      ,                         &
                          pstate%u_wind_at_pc_full_level%f(:,1:ncol)   , &
                          pstate%v_wind_at_pc_full_level%f(:,1:ncol)   , &
                          pstate%temp_at_pc_full_level%f(:,1:ncol)     , & 
                          pstate%pressure_at_pc_full_level%f(:,1:ncol) , & 
                          pstate%exner_at_pc_full_level%f(:,1:ncol)    , &
                          pstate%z_at_pc_full_level%f(:,1:ncol)        , &
                          pstate_cam%sgh30_at_pc_surface%f(1:ncol)     , &
                          ksrftms  ,  tautmsx  , tautmsy               , &
                          pstate%landfrac_at_pc_surface%f(1:ncol) )


      ! Here, both 'taux, tautmsx' are explicit surface stresses.        
      ! Note that this 'tautotx, tautoty' are different from the total stress
      ! that has been actually added into the atmosphere. This is because both
      ! taux and tautmsx are fully implicitly treated within compute_vdiff.
      ! However, 'tautotx, tautoty' are not used in the actual numerical
      ! computation in this module.   

        tautotx(:ncol) = pstate%atm_in_taux_at_pc_surface%f(1:ncol) + tautmsx(1:ncol)
        tautoty(:ncol) = pstate%atm_in_tauy_at_pc_surface%f(1:ncol) + tautmsy(1:ncol)
    else
        ksrftms(:ncol) = 0._r8
        tautotx(:ncol) = pstate%atm_in_taux_at_pc_surface%f(1:ncol)
        tautoty(:ncol) = pstate%atm_in_tauy_at_pc_surface%f(1:ncol)
    endif

    kvm_in(:,1:ncol)   = pstate_cam%pbl_kvm_at_pc_face_level%f(:,1:ncol)
    kvh_in(:,1:ncol)   = pstate_cam%pbl_kvh_at_pc_face_level%f(:,1:ncol)
    
    call compute_eddy_diff(ncol, nlev, pstate%temp_at_pc_full_level%f(:,1:ncol),            &
                           pstate%tracer_mxrt_at_pc_full_level%f(1:3,:,1:ncol), dt, rpdel,  &
                           pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:,1:ncol),   &
                           pstate_cam%lw_qrl_at_pc_full_level%f(:,1:ncol),                  &
                           pstate_cam%microp_wsedl_at_pc_full_level%f(:,1:ncol),            &
                           pstate%z_at_pc_full_level%f(:,1:ncol),                           & 
                           pstate%z_at_pc_face_level%f(:,1:ncol),                           &
                           pstate%pressure_at_pc_full_level%f(:,1:ncol),                    &
                           pstate%pressure_at_pc_face_level%f(:,1:ncol),                    &
                           pstate%u_wind_at_pc_full_level%f(:,1:ncol),                      &
                           pstate%v_wind_at_pc_full_level%f(:,1:ncol),                      &
                           pstate%atm_in_taux_at_pc_surface%f(1:ncol),                      &
                           pstate%atm_in_tauy_at_pc_surface%f(1:ncol),                      &
                           pstate%atm_in_shflx_at_pc_surface%f(1:ncol),                     &
                           pstate%atm_in_qflx_at_pc_surface%f(1,1:ncol),                    &
                           wstarent, nturb, rrho, ustar,                                    &
                           pstate_cam%pbl_pblh_at_pc_surface%f(1:ncol),                     &
                           kvm_in, kvh_in, kvm, kvh,                                        &
                           kvq, cgh, cgs, bprod, sprod,                                     &
                           pstate_cam%pbl_tke_at_pc_face_level%f(:,1:ncol), wpert,          &
                           pstate_cam%pbl_tpert_at_pc_surface%f(1:ncol),                    &
                           pstate_cam%pbl_qpert_at_pc_surface%f(1:ncol),                    &
                           sfi, kvinit,                                                     &
                           pstate_cam%pbl_tauresx_at_pc_surface%f(1:ncol),                  &
                           pstate_cam%pbl_tauresy_at_pc_surface%f(1:ncol),                  &  
                           ksrftms, turbtype, smaw,                                         &
                           ! for Lin Macro ------------------------>
                           !ipbl(1:ncol), kpblh(1:ncol), wstarPBL(1:ncol) )
                           ipbl(1:ncol), kpblh(1:ncol),                                     &
                           pstate_cam%pbl_wstarPBL%f(1:ncol),                               &
                           pstate_cam%pbl_lengi_at_pc_face_level%f(:,1:ncol),               &
                           pstate_cam%pbl_shi_at_pc_face_level%f(:,1:ncol),                 &
                           pstate_cam%pbl_smi_at_pc_face_level%f(:,1:ncol) )
                           ! for Lin Macro <------------------------

    ! The diag_TKE scheme does not calculate the Monin-Obukhov length, which is used in dry deposition calculations.
    ! Use the routines from pbl_utils to accomplish this. Assumes ustar and rrho have been set.
    th(nlev,1:ncol) = pstate%temp_at_pc_full_level%f(nlev,1:ncol)*                     &
                      pstate%exner_at_pc_full_level%f(nlev,1:ncol)
    thvs(1:ncol) = virtem(th(nlev,1:ncol),pstate%tracer_mxrt_at_pc_full_level%f(1,nlev,1:ncol))
    call calc_obklen(th(nlev,1:ncol), thvs(1:ncol), pstate%atm_in_qflx_at_pc_surface%f(1,1:ncol),        &
                     pstate%atm_in_shflx_at_pc_surface%f(1:ncol), rrho(1:ncol), ustar(1:ncol),           &
                     khfs(1:ncol), kqfs(1:ncol), kbfs(1:ncol), obklen(1:ncol))


       ! The diffusivities from diag_TKE can be much larger than from HB in the free
       ! troposphere and upper atmosphere. These seem to be larger than observations,
       ! and in WACCM the gw_drag code is already applying an eddy diffusivity in the
       ! upper atmosphere. Optionally, adjust the diffusivities in the free troposphere
       ! or the upper atmosphere.
       !
       ! NOTE: Further investigation should be done as to why the diffusivities are
       ! larger in diag_TKE.
       if ((kv_freetrop_scale /= 1._r8) .or. ((kv_top_scale /= 1._r8) .and. (kv_top_pressure > 0._r8))) then
         do i = 1, ncol
           do k = 1, nlevp
           
             ! Outside of the boundary layer?
             if (pstate%z_at_pc_face_level%f(k,i) > pstate_cam%pbl_pblh_at_pc_surface%f(i)) then

               ! In the upper atmosphere?
               if (pstate%pressure_at_pc_face_level%f(k,i) <= kv_top_pressure) then
                 kvh(k,i) = kvh(k,i) * kv_top_scale
                 kvm(k,i) = kvm(k,i) * kv_top_scale
                 kvq(k,i) = kvq(k,i) * kv_top_scale
               else
                 kvh(k,i) = kvh(k,i) * kv_freetrop_scale
                 kvm(k,i) = kvm(k,i) * kv_freetrop_scale
                 kvq(k,i) = kvq(k,i) * kv_freetrop_scale
               end if
             else
               exit
             end if
           end do
         end do
       end if

    ! kvh (in pbuf) is used by other physics parameterizations, and as an initial guess in compute_eddy_diff
    ! on the next timestep.  It is not updated by the compute_vdiff call below.
    pstate_cam%pbl_kvh_at_pc_face_level%f(:,1:ncol) = kvh(:,1:ncol)

    ! kvm (in pbuf) is only used as an initial guess in compute_eddy_diff on the next timestep.
    ! The contributions for molecular diffusion made to kvm by the call to compute_vdiff below 
    ! are not included in the pbuf as these are not needed in the initial guess by compute_eddy_diff.
    pstate_cam%pbl_kvm_at_pc_face_level%f(:,1:ncol) = kvm(:,1:ncol)

    ! Set arrays from input state.
    q_tmp(1:ntracer,:,1:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(1:ntracer,:,1:ncol)
    s_tmp(:,:ncol)   = pstate%static_energy_at_pc_full_level%f(:,1:ncol)
    u_tmp(:,:ncol)   = pstate%u_wind_at_pc_full_level%f(:,1:ncol)
    v_tmp(:,:ncol)   = pstate%v_wind_at_pc_full_level%f(:,1:ncol)

!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
!      rairi(:ncol,1) = rairv(:ncol,1,lchnk)
!      do k = 2, pver
!        do i = 1, ncol
!          rairi(i,k) = 0.5_r8 * (rairv(i,k,lchnk)+rairv(i,k-1,lchnk))
!        end do
!      end do
!    else
!      rairi(:ncol,:pver+1) = rair 
!    endif
    rairi(:nlevp,1:ncol) = rdry
    cpairv(:nlev,1:ncol) = cp
!<----------------------------------Lixh close this part-----------------------------------

    ! Modification : We may need to output 'tautotx_im,tautoty_im' from below 'compute_vdiff' and
    !                separately print out as diagnostic output, because these are different from
    !                the explicit 'tautotx, tautoty' computed above. 
    ! Note that the output 'tauresx,tauresy' from below subroutines are fully implicit ones.

    !qmincg(1:ntracer) = phy_tracer_info(1:ntracer)%qmin
    qmincg(1:ntracer) = 0._r8

    if( any(fieldlist_wet) ) then
    call compute_vdiff( nlev          , ntracer              , ncol         ,                                   &
                        pstate%pressure_at_pc_full_level%f(:,1:ncol)        ,                                   &
                        pstate%pressure_at_pc_face_level%f(:,1:ncol)        ,                                   &
                        rpdel         , pstate%temp_at_pc_full_level%f(:,1:ncol)             , dt            ,  &
                        pstate%atm_in_taux_at_pc_surface%f(1:ncol) , pstate%atm_in_tauy_at_pc_surface%f(1:ncol)          , &      
                        pstate%atm_in_shflx_at_pc_surface%f(1:ncol), pstate%atm_in_qflx_at_pc_surface%f(1:ntracer,1:ncol), & 
                        ntop          , nbot                 , kvh          , kvm            , kvq           ,  &
                        cgs           , cgh                  , pstate%z_at_pc_face_level%f(:,1:ncol)         ,  &
                        ksrftms       , qmincg               , fieldlist_wet,                                   &
                        u_tmp         , v_tmp                , q_tmp        , s_tmp          ,                  &
                        tautmsx       , tautmsy              , dtk          , topflx         ,                  &
                        pstate_cam%pbl_tauresx_at_pc_surface%f(1:ncol),pstate_cam%pbl_tauresy_at_pc_surface%f(1:ncol),     &
                        1             , cpairv(:,:ncol)      , rairi(:,:ncol),                                  &
                        do_molec_diff , compute_molec_diff   , vd_lu_qdecomp,                                   &
                        pstate_cam%pbl_kvt_at_pc_face_level%f(:,1:ncol))

    end if

    !------------------Lixh has not completed fieldlist_dry----------------->
    !if( any( fieldlist_dry ) ) then
    !
    !end if
    !<-----------------Lixh has not completed fieldlist_dry------------------

!    if (prog_modal_aero) then      !.false. LiXH

       ! Modal aerosol species not diffused, so just add the explicit surface fluxes to the
       ! lowest layer

!    tmp1(:ncol) = ztodt * gravit * state%rpdel(:ncol,pver)
!       do m = 1, pmam_ncnst
!          l = pmam_cnst_idx(m)
!          q_tmp(:ncol,pver,l) = q_tmp(:ncol,pver,l) + tmp1(:ncol) * cflx(:ncol,l)
!       enddo
!    end if
    
    ! re-initialize ptend
    ptend_vertical_diffusion%tend_u%f = 0._r8
    ptend_vertical_diffusion%tend_v%f = 0._r8
    ptend_vertical_diffusion%tend_s%f = 0._r8
    ptend_vertical_diffusion%tend_q%f = 0._r8

    ptend_vertical_diffusion%tend_u%f(:,1:ncol) = (u_tmp(:,1:ncol)-pstate%u_wind_at_pc_full_level%f(:,1:ncol))*rdt
    ptend_vertical_diffusion%tend_v%f(:,1:ncol) = (v_tmp(:,1:ncol)-pstate%v_wind_at_pc_full_level%f(:,1:ncol))*rdt
    ptend_vertical_diffusion%tend_s%f(:,1:ncol) = (s_tmp(:,1:ncol)-pstate%static_energy_at_pc_full_level%f(:,1:ncol))*rdt
    ptend_vertical_diffusion%tend_q%f(:,:,1:ncol) = (q_tmp(:,:,1:ncol)-pstate%tracer_mxrt_at_pc_full_level%f(:,:,1:ncol))*rdt

    if (diff_cnsrv_mass_check) then
        ! mass conservation check is not completed, LiXH
        call endrun('LiXH has not completed diff_cnsrv_mass_check !!')
    end if

    if(kvinit)kvinit = .false.
 
    end subroutine vertical_diffusion_tend

 end module grist_vertical_diffusion
