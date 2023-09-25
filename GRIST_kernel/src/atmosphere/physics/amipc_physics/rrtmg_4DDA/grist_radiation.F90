!======================================================
!
!  Created by LiXiaohan on 19/8/13.
!  interface to RRTMG
!======================================================

 module grist_radiation

    use grist_constants,                    only: r8, stebol, cp
    use grist_nml_module,                   only: nlev, nlevp, working_mode, start_tod, start_ymd
    use grist_handle_error,                 only: endrun
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam ,       &
                                                  ptend_radiation
    use grist_physics_update,               only: old_time_level
    use radconstants,                       only: rrtmg_sw_cloudsim_band,    &
                                                  rrtmg_lw_cloudsim_band,    &
                                                  nswbands, nlwbands
    use grist_rad_constituents,             only: N_DIAG, rad_cnst_get_info, &
                                                  rad_cnst_get_call_list
    use grist_mpi

    implicit none
    private

    public :: radiation_defaultopts, &! set default values of namelist variables in runtime_opts
              radiation_setopts,     &! set namelist values from runtime_opts
              radiation_do,          &! query which radiation calcs are done this timestep
              radiation_init,        &! calls radini
              read_nml_radiation,    &! read nml of radiation
              radiation_tend,        &! moved from radctl.F90
              end_of_radiation        ! destruct radiation variables  

    integer,public, allocatable :: cosp_cnt(:)       ! counter for cosp
    integer,public              :: cosp_cnt_init = 0 !initial value for cosp counter
    logical,public :: call_radiation
    
    ! Default values for namelist variables
    
    integer :: iradsw = -1     ! freq. of shortwave radiation calc in time steps (positive)
                               ! or hours (negative).
    integer :: iradlw = -1     ! frequency of longwave rad. calc. in time steps (positive)
                               ! or hours (negative).
    
    integer :: irad_always = 0 ! Specifies length of time in timesteps (positive)
                               ! or hours (negative) SW/LW radiation will be
                               ! run continuously from the start of an
                               ! initial or restart run
    integer :: start_step = 0  ! index of the start time step 
    logical :: spectralflux  = .false. ! calculate fluxes (up and down) per band.
    
    character(len=4) :: diag(0:N_DIAG) =(/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ','_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)
    
    logical :: dohirs = .false. ! diagnostic  brightness temperatures at the top of the
                                ! atmosphere for 7 TOVS/HIRS channels (2,4,6,8,10,11,12) and 4 TOVS/MSU 
                                ! channels (1,2,3,4).
    integer :: ihirsfq = 1      ! frequency (timesteps) of brightness temperature calcs

    logical :: FLG_4DDA_SW = .true.  ! flag for shortwave 4DDA of Feng Zhang, Linhan add
    logical :: FLG_4DDA_LW = .true.  ! flag for longwave 4DDA of Feng Zhang,  Linhan add 

    ! if the variables below are used by other modules (pbuf in CAM), change to pstate%~, Lixh

    real(r8), allocatable :: su(:,:,:)        ! shortwave upward flux (per band)
    real(r8), allocatable :: sd(:,:,:)        ! shortwave downward flux (per band)
    real(r8), allocatable :: lu(:,:,:)        ! longwave upward flux (per band)
    real(r8), allocatable :: ld(:,:,:)        ! longwave downward flux (per band)
 
    ! if the variables below are used by other modules (comsrf in CAM), change to pstate%~, Lixh
    real(r8), allocatable :: fsns(:)          ! surface absorbed solar flux
    real(r8), allocatable :: fsnt(:)          ! Net column abs solar flux at model top
    real(r8), allocatable :: flns(:)          ! Srf longwave cooling (up-down) flux
    real(r8), allocatable :: flnt(:)          ! Net outgoing lw flux at model top

    !ghg index, LiXH
    integer :: o3_idx, co2_idx


contains
 
! Purpose: Return default runtime options
    subroutine radiation_defaultopts(iradsw_out, iradlw_out, iradae_out, irad_always_out, spectralflux_out)

    integer, intent(out), optional :: iradsw_out
    integer, intent(out), optional :: iradlw_out
    integer, intent(out), optional :: iradae_out
    integer, intent(out), optional :: irad_always_out
    logical, intent(out), optional :: spectralflux_out

    if ( present(iradsw_out) )      iradsw_out = iradsw
    if ( present(iradlw_out) )      iradlw_out = iradlw
    if ( present(iradae_out) )      iradae_out = -999
    if ( present(irad_always_out) ) irad_always_out = irad_always
    if ( present(spectralflux_out) ) spectralflux_out = spectralflux

    end subroutine radiation_defaultopts


! Purpose: Set runtime options
! *** NOTE *** This routine needs information about dtime (init by dycore) 
!              and nhtfrq (init by history) to do its checking.  Being called
!              from runtime_opts provides these values possibly before they
!              have been set in the modules responsible for them.
    subroutine radiation_setopts(dtime, start_in, iradsw_in, iradlw_in, iradae_in, &
                                 irad_always_in, spectralflux_in)
    ! io
    real(r8),intent(in)           :: dtime           ! timestep size (s)
    integer, intent(in), optional :: start_in
    integer, intent(in), optional :: iradsw_in
    integer, intent(in), optional :: iradlw_in
    integer, intent(in), optional :: iradae_in
    integer, intent(in), optional :: irad_always_in
    logical, intent(in), optional :: spectralflux_in

    ! local
    integer :: ntspdy   ! no. timesteps per day
    integer :: nhtfrq1  ! local copy of input arg nhtfrq
    integer :: iradae   ! not used by RRTMG

    if ( present(start_in) )       start_step = start_in
    if ( present(iradsw_in) )      iradsw = iradsw_in
    if ( present(iradlw_in) )      iradlw = iradlw_in
    if ( present(iradae_in) )      iradae = iradae_in
    if ( present(irad_always_in) ) irad_always = irad_always_in

    ! Convert iradsw, iradlw and irad_always from hours to timesteps if necessary
    if (iradsw      < 0) iradsw      = nint((-iradsw     *3600._r8)/dtime)
    if (iradlw      < 0) iradlw      = nint((-iradlw     *3600._r8)/dtime)
    if (irad_always < 0) irad_always = nint((-irad_always*3600._r8)/dtime)

    ! Has user specified iradae?
    if (iradae /= -999) then
       call endrun('radiation_setopts: iradae not used by RRTMG.')
    end if

    end subroutine radiation_setopts


! Purpose: Returns true if the specified operation is done this timestep.
    function radiation_do(op, nstep)
    ! io
    character(len=*), intent(in)  :: op             ! name of operation
    integer, intent(in)           :: nstep
    logical                       :: radiation_do   ! return value

    select case (op)

    case ('sw') ! do a shortwave heating calc this timestep?
                    !----------LiXH add------------
                    ! gcm start step =1, scm start step =0
                    ! this part should be unified. LiXH
                    !----------LiXH add------------
       radiation_do = nstep == start_step  .or.  iradsw == 1            &
#ifdef RAD1
                     .or. nstep == 1                                    &
#endif
                     .or. (mod(nstep-1,iradsw) == 0  .and.  nstep /= 1) &
                     .or. nstep <= irad_always

    case ('lw') ! do a longwave heating calc this timestep?
       radiation_do = nstep == start_step  .or.  iradlw == 1            &
#ifdef RAD1
                     .or. nstep == 1                                    &
#endif
                     .or. (mod(nstep-1,iradlw) == 0  .and.  nstep /= 1) &
                     .or. nstep <= irad_always

    case ('aeres') ! write absorptivity/emissivity to restart file this timestep?
       ! for RRTMG there is no abs/ems restart file
       radiation_do = .false.
          
    case default
       call endrun('radiation_do: unknown operation:'//op)

    end select
    end function radiation_do


! Purpose: read radiation namelist
    subroutine read_nml_radiation(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer :: unitn, ierr

    namelist /radiation_nl/ FLG_4DDA_SW,  FLG_4DDA_LW

    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, radiation_nl, iostat=ierr)
    if (ierr /= 0) call endrun(' error reading radiation namelist ')
    close(unitn)

    end subroutine read_nml_radiation


! Purpose: Initialize the radiation parameterization
    subroutine radiation_init(ncol,dtime)

    use radsw,          only: radsw_init
    use radlw,          only: radlw_init
    use rrtmg_state,    only: rrtmg_state_init
    use modal_aer_opt,  only: modal_aer_opt_init
    use hirsbt,         only: hirsbt_init

    ! io 
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: dtime
    ! local
    integer :: icall, nmodes, i
    logical :: active_calls(0:N_DIAG)

    allocate(fsns(ncol));fsns=0._r8
    allocate(fsnt(ncol));fsnt=0._r8
    allocate(flns(ncol));flns=0._r8
    allocate(flnt(ncol));flnt=0._r8
    allocate(su(nswbands, nlevp, ncol));su=0._r8
    allocate(sd(nswbands, nlevp, ncol));sd=0._r8
    allocate(lu(nlwbands, nlevp, ncol));lu=0._r8
    allocate(ld(nlwbands, nlevp, ncol));ld=0._r8

    allocate(pstate_cam%lw_qrl_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%lw_qrlc_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%flwdsc_at_pc_surface%f(ncol))
    allocate(pstate_cam%flwusc_at_pc_surface%f(ncol))
    allocate(pstate_cam%flwdtc_at_pc_top%f(ncol))
    allocate(pstate_cam%flwutc_at_pc_top%f(ncol))
    allocate(pstate_cam%flwds_at_pc_surface%f(ncol))
    allocate(pstate_cam%flwus_at_pc_surface%f(ncol))
    allocate(pstate_cam%flwdt_at_pc_top%f(ncol))
    allocate(pstate_cam%flwut_at_pc_top%f(ncol))
    allocate(pstate_cam%lwcf_at_pc_top%f(ncol))

    pstate_cam%lw_qrl_at_pc_full_level%pos       = 0
    pstate_cam%lw_qrlc_at_pc_full_level%pos      = 0
    pstate_cam%flwdsc_at_pc_surface%pos          = 0   
    pstate_cam%flwusc_at_pc_surface%pos          = 0
    pstate_cam%flwdtc_at_pc_top%pos              = 0
    pstate_cam%flwutc_at_pc_top%pos              = 0
    pstate_cam%flwds_at_pc_surface%pos           = 0
    pstate_cam%flwus_at_pc_surface%pos           = 0
    pstate_cam%flwdt_at_pc_top%pos               = 0
    pstate_cam%flwut_at_pc_top%pos               = 0
    pstate_cam%lwcf_at_pc_top%pos                = 0
 
    pstate_cam%lw_qrl_at_pc_full_level%f         = 0._r8
    pstate_cam%lw_qrlc_at_pc_full_level%f        = 0._r8
    pstate_cam%flwdsc_at_pc_surface%f            = 0._r8   
    pstate_cam%flwusc_at_pc_surface%f            = 0._r8
    pstate_cam%flwdtc_at_pc_top%f                = 0._r8
    pstate_cam%flwutc_at_pc_top%f                = 0._r8
    pstate_cam%flwds_at_pc_surface%f             = 0._r8
    pstate_cam%flwus_at_pc_surface%f             = 0._r8
    pstate_cam%flwdt_at_pc_top%f                 = 0._r8
    pstate_cam%flwut_at_pc_top%f                 = 0._r8
    pstate_cam%lwcf_at_pc_top%f                  = 0._r8
 
    allocate(pstate_cam%sw_qrs_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%sw_qrsc_at_pc_full_level%f(nlev,ncol))
    allocate(pstate_cam%fswdsc_at_pc_surface%f(ncol))
    allocate(pstate_cam%fswusc_at_pc_surface%f(ncol))
    allocate(pstate_cam%fswdtc_at_pc_top%f(ncol))
    allocate(pstate_cam%fswutc_at_pc_top%f(ncol))
    allocate(pstate_cam%fswds_at_pc_surface%f(ncol))
    allocate(pstate_cam%fswus_at_pc_surface%f(ncol))
    allocate(pstate_cam%fswdt_at_pc_top%f(ncol))
    allocate(pstate_cam%fswut_at_pc_top%f(ncol))
    allocate(pstate_cam%fsntoa_at_pc_top%f(ncol))
    allocate(pstate_cam%fsntoac_at_pc_top%f(ncol))
    allocate(pstate_cam%swcf_at_pc_top%f(ncol))

    pstate_cam%sw_qrs_at_pc_full_level%pos       = 0
    pstate_cam%sw_qrsc_at_pc_full_level%pos      = 0
    pstate_cam%fswdsc_at_pc_surface%pos          = 0   
    pstate_cam%fswusc_at_pc_surface%pos          = 0
    pstate_cam%fswdtc_at_pc_top%pos              = 0
    pstate_cam%fswutc_at_pc_top%pos              = 0
    pstate_cam%fswds_at_pc_surface%pos           = 0
    pstate_cam%fswus_at_pc_surface%pos           = 0
    pstate_cam%fswdt_at_pc_top%pos               = 0
    pstate_cam%fswut_at_pc_top%pos               = 0
    pstate_cam%fsntoa_at_pc_top%pos              = 0
    pstate_cam%fsntoac_at_pc_top%pos             = 0
    pstate_cam%swcf_at_pc_top%pos                = 0
 
    pstate_cam%sw_qrs_at_pc_full_level%f         = 0._r8
    pstate_cam%sw_qrsc_at_pc_full_level%f        = 0._r8
    pstate_cam%fswdsc_at_pc_surface%f            = 0._r8   
    pstate_cam%fswusc_at_pc_surface%f            = 0._r8
    pstate_cam%fswdtc_at_pc_top%f                = 0._r8
    pstate_cam%fswutc_at_pc_top%f                = 0._r8
    pstate_cam%fswds_at_pc_surface%f             = 0._r8
    pstate_cam%fswus_at_pc_surface%f             = 0._r8
    pstate_cam%fswdt_at_pc_top%f                 = 0._r8
    pstate_cam%fswut_at_pc_top%f                 = 0._r8
    pstate_cam%fsntoa_at_pc_top%f                = 0._r8
    pstate_cam%fsntoac_at_pc_top%f               = 0._r8
    pstate_cam%swcf_at_pc_top%f                  = 0._r8
 
    call rrtmg_state_init()

    call radsw_init()
    call radlw_init()

    ! Determine whether modal aerosols are affecting the climate, and if so
    ! then initialize the modal aerosol optics module
    call rad_cnst_get_info(0, nmodes=nmodes)

    if (nmodes > 0) call modal_aer_opt_init()
    call hirsbt_init(dtime)

    ! "irad_always" is number of time steps to execute radiation continuously from start of
    ! initial OR restart run
!-----------LiXH modified, nstep= 0 in initial part-------------->   
!    nstep = get_nstep()
!    if ( irad_always > 0) then
!       nstep       = get_nstep()
!       irad_always = irad_always + nstep
!    end if
!<----------LiXH modified, nstep= 0 in initial part---------------   


!------------LiXH has not completed COSP------------>
!    if (docosp) call cospsimulator_intr_init
!    allocate(cosp_cnt(begchunk:endchunk))
!    if (is_first_restart_step()) then
!      cosp_cnt(begchunk:endchunk)=cosp_cnt_init
!    else
!      cosp_cnt(begchunk:endchunk)=0     
!    end if
!<------------LiXH has not completed COSP------------

    ! get list of active radiation calls
    call rad_cnst_get_call_list(active_calls)

    ! get ghg index, LiXH
    co2_idx = -1
    o3_idx  = -1
    do i = 1, pstate_cam%total_ghg_num
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'CO2')co2_idx = i
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'O3')o3_idx = i
    end do
    if(co2_idx .lt. 0 .or. o3_idx .lt. 0)then
        if(mpi_rank()==0)print*,'can not find co2 or o3 index in ghg data!'
        call endrun('radiation_init')
    end if
    end subroutine radiation_init


    subroutine end_of_radiation

    deallocate(fsns)
    deallocate(fsnt)
    deallocate(flns)
    deallocate(flnt)
    deallocate(su)
    deallocate(sd)
    deallocate(lu)
    deallocate(ld)

    if(allocated(pstate_cam%lw_qrl_at_pc_full_level%f))          &
      deallocate(pstate_cam%lw_qrl_at_pc_full_level%f)
    if(allocated(pstate_cam%lw_qrlc_at_pc_full_level%f))         &
      deallocate(pstate_cam%lw_qrlc_at_pc_full_level%f)
    if(allocated(pstate_cam%flwdsc_at_pc_surface%f))             &
      deallocate(pstate_cam%flwdsc_at_pc_surface%f)
    if(allocated(pstate_cam%flwusc_at_pc_surface%f))             & 
      deallocate(pstate_cam%flwusc_at_pc_surface%f)
    if(allocated(pstate_cam%flwdtc_at_pc_top%f))                 & 
      deallocate(pstate_cam%flwdtc_at_pc_top%f)
    if(allocated(pstate_cam%flwutc_at_pc_top%f))                 & 
      deallocate(pstate_cam%flwutc_at_pc_top%f)
    if(allocated(pstate_cam%flwds_at_pc_surface%f))              & 
      deallocate(pstate_cam%flwds_at_pc_surface%f)
    if(allocated(pstate_cam%flwus_at_pc_surface%f))              & 
      deallocate(pstate_cam%flwus_at_pc_surface%f)
    if(allocated(pstate_cam%flwdt_at_pc_top%f))                  & 
      deallocate(pstate_cam%flwdt_at_pc_top%f)
    if(allocated(pstate_cam%flwut_at_pc_top%f))                  & 
      deallocate(pstate_cam%flwut_at_pc_top%f)
    if(allocated(pstate_cam%lwcf_at_pc_top%f))                   & 
      deallocate(pstate_cam%lwcf_at_pc_top%f)

    if(allocated(pstate_cam%sw_qrs_at_pc_full_level%f))          &
      deallocate(pstate_cam%sw_qrs_at_pc_full_level%f)
    if(allocated(pstate_cam%sw_qrsc_at_pc_full_level%f))         &
      deallocate(pstate_cam%sw_qrsc_at_pc_full_level%f)
    if(allocated(pstate_cam%fswdsc_at_pc_surface%f))             &
      deallocate(pstate_cam%fswdsc_at_pc_surface%f)
    if(allocated(pstate_cam%fswusc_at_pc_surface%f))             & 
      deallocate(pstate_cam%fswusc_at_pc_surface%f)
    if(allocated(pstate_cam%fswdtc_at_pc_top%f))                 & 
      deallocate(pstate_cam%fswdtc_at_pc_top%f)
    if(allocated(pstate_cam%fswutc_at_pc_top%f))                 & 
      deallocate(pstate_cam%fswutc_at_pc_top%f)
    if(allocated(pstate_cam%fswds_at_pc_surface%f))              & 
      deallocate(pstate_cam%fswds_at_pc_surface%f)
    if(allocated(pstate_cam%fswus_at_pc_surface%f))              & 
      deallocate(pstate_cam%fswus_at_pc_surface%f)
    if(allocated(pstate_cam%fswdt_at_pc_top%f))                  & 
      deallocate(pstate_cam%fswdt_at_pc_top%f)
    if(allocated(pstate_cam%fswut_at_pc_top%f))                  & 
      deallocate(pstate_cam%fswut_at_pc_top%f)
    if(allocated(pstate_cam%fsntoa_at_pc_top%f))                 &
      deallocate(pstate_cam%fsntoa_at_pc_top%f)
    if(allocated(pstate_cam%fsntoac_at_pc_top%f))                &
      deallocate(pstate_cam%fsntoac_at_pc_top%f)
    if(allocated(pstate_cam%swcf_at_pc_top%f))                   & 
      deallocate(pstate_cam%swcf_at_pc_top%f)

    end subroutine end_of_radiation


    subroutine radiation_tend(ncol, nstep, dtime, lon, lat, net_flx)
    use grist_physics_update,           only: old_time_level
    use grist_time_manager,             only: get_curr_calday
    use grist_zenith,                   only: zenith
    use parrrtm,                        only: nbndlw
    use parrrsw,                        only: nbndsw
    use hirsbt,                         only: hirsrtm
    use hirsbtpar,                      only: pnb_hirs, pnf_msu,             &
                                              hirsname, msuname
    use rrtmg_state,                    only: rrtmg_state_t,                 &
                                              rrtmg_state_create,            &
                                              rrtmg_state_update,            &
                                              rrtmg_state_destroy,           &
                                              num_rrtmg_levs
    use grist_rad_constituents,         only: oldcldoptics,                  &
                                              liqcldoptics,                  &
                                              icecldoptics
    use radsw,                          only: rad_rrtmg_sw
    use radlw,                          only: rad_rrtmg_lw
    use aer_rad_props,                  only: aer_rad_props_sw,              &
                                              aer_rad_props_lw  
    use cloud_rad_props,                only: get_ice_optics_sw,             &
                                              get_liquid_optics_sw,          &   
                                              liquid_cloud_get_rad_props_lw, &
                                              ice_cloud_get_rad_props_lw,    &
                                              cloud_rad_props_get_lw,        &
                                              snow_cloud_get_rad_props_lw,   &
                                              get_snow_optics_sw
    use slingo,                         only: slingo_liq_get_rad_props_lw,   &
                                              slingo_liq_optics_sw
    use ebert_curry,                    only: ec_ice_optics_sw,              &
                                              ec_ice_get_rad_props_lw
    use rad_solar_var,                  only: get_variability
    !use grist_scm_comm_module,         only: have_tg, tground
    use grist_physics_data_structure,   only: have_tg, tground
    ! io
    integer, intent(in)     :: ncol
    integer, intent(in)     :: nstep
    real(r8), intent(in)    :: dtime
    real(r8), intent(in)    :: lon(ncol), lat(ncol)
    real(r8), intent(inout) :: net_flx(ncol)
    ! local
    integer  :: i, k
    integer  :: istat
    logical  :: dosw, dolw
    logical  :: cldfsnow_idx
    real(r8) :: calday                         ! current calendar day
    real(r8) :: britemp(pnf_msu,ncol)          ! Microwave brightness temperature
    real(r8) :: tb_ir(pnb_hirs,ncol)           ! Infrared brightness temperature
    real(r8) :: ts(ncol)                       ! surface temperature
    real(r8) :: pintmb(nlevp,ncol)             ! Model interface pressures (hPa)
    real(r8) :: oro(ncol)                      ! Land surface flag, sea=0, land=1
    real(r8) :: emis(nlev,ncol)                ! Cloud longwave emissivity
    real(r8) :: coszrs(ncol)                   ! Cosine solar zenith angle
    logical  :: conserve_energy = .true.       ! flag to carry (QRS,QRL)*dp across time steps
    integer  :: Nday                           ! Number of daylight columns
    integer  :: Nnite                          ! Number of night columns
    integer, dimension(ncol) :: IdxDay         ! Indicies of daylight coumns
    integer, dimension(ncol) :: IdxNite        ! Indicies of night coumns

    ! combined cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: c_cld_tau    (nbndsw,nlev,ncol)! cloud extinction optical depth
    real(r8) :: c_cld_tau_w  (nbndsw,nlev,ncol)! cloud single scattering albedo * tau
    real(r8) :: c_cld_tau_w_g(nbndsw,nlev,ncol)! cloud assymetry parameter * w * tau
    real(r8) :: c_cld_tau_w_f(nbndsw,nlev,ncol)! cloud forward scattered fraction * w * tau
    real(r8) :: c_cld_lw_abs (nbndlw,nlev,ncol)! cloud absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: cld_tau    (nbndsw,nlev,ncol)  ! cloud extinction optical depth
    real(r8) :: cld_tau_w  (nbndsw,nlev,ncol)  ! cloud single scattering albedo * tau
    real(r8) :: cld_tau_w_g(nbndsw,nlev,ncol)  ! cloud assymetry parameter * w * tau
    real(r8) :: cld_tau_w_f(nbndsw,nlev,ncol)  ! cloud forward scattered fraction * w * tau
    real(r8) :: cld_lw_abs (nbndlw,nlev,ncol)  ! cloud absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: ice_tau    (nbndsw,nlev,ncol)  ! ice extinction optical depth
    real(r8) :: ice_tau_w  (nbndsw,nlev,ncol)  ! ice single scattering albedo * tau
    real(r8) :: ice_tau_w_g(nbndsw,nlev,ncol)  ! ice assymetry parameter * tau * w
    real(r8) :: ice_tau_w_f(nbndsw,nlev,ncol)  ! ice forward scattered fraction * tau * w
    real(r8) :: ice_lw_abs (nbndlw,nlev,ncol)  ! ice absorption optics depth (LW)

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: snow_tau    (nbndsw,nlev,ncol) ! snow extinction optical depth
    real(r8) :: snow_tau_w  (nbndsw,nlev,ncol) ! snow single scattering albedo * tau
    real(r8) :: snow_tau_w_g(nbndsw,nlev,ncol) ! snow assymetry parameter * tau * w
    real(r8) :: snow_tau_w_f(nbndsw,nlev,ncol) ! snow forward scattered fraction * tau * w
    real(r8) :: snow_lw_abs (nbndlw,nlev,ncol) ! snow absorption optics depth (LW)
    real(r8) :: gb_snow_tau        (nlev,ncol) ! grid-box mean snow_tau for COSP only
    real(r8) :: gb_snow_lw         (nlev,ncol) ! grid-box mean LW snow optical depth for COSP only

    ! cloud radiative parameters are "in cloud" not "in cell"
    real(r8) :: liq_tau    (nbndsw,nlev,ncol)  ! liquid extinction optical depth
    real(r8) :: liq_tau_w  (nbndsw,nlev,ncol)  ! liquid single scattering albedo * tau
    real(r8) :: liq_tau_w_g(nbndsw,nlev,ncol)  ! liquid assymetry parameter * tau * w
    real(r8) :: liq_tau_w_f(nbndsw,nlev,ncol)  ! liquid forward scattered fraction * tau * w
    real(r8) :: liq_lw_abs (nbndlw,nlev,ncol)  ! liquid absorption optics depth (LW)

    real(r8) :: cld(nlev,ncol)                 ! cloud fraction 
    real(r8) :: cldfsnow(nlev,ncol)            ! cloud fraction of just "snow clouds- whatever they are" 
    real(r8) :: cldfprime(nlev,ncol)           ! combined cloud fraction (snow plus regular)

    real(r8) solin(ncol)         ! Solar incident flux
    real(r8) fsutoa(ncol)        ! Upwelling solar flux at TOA
    real(r8) fsnirt(ncol)        ! Near-IR flux absorbed at toa
    real(r8) fsnrtc(ncol)        ! Clear sky near-IR flux absorbed at toa
    real(r8) fsnirtsq(ncol)      ! Near-IR flux absorbed at toa >= 0.7 microns
    real(r8) fsntc(ncol)         ! Clear sky total column abs solar flux
    real(r8) fsnsc(ncol)         ! Clear sky surface abs solar flux
    real(r8) flntc(ncol)         ! Clear sky lw flux at model top
    real(r8) flnsc(ncol)         ! Clear sky lw flux at srf (up-down)
    real(r8) fln200(ncol)        ! net longwave flux interpolated to 200 mb
    real(r8) fln200c(ncol)       ! net clearsky longwave flux interpolated to 200 mb
    real(r8) fns(nlevp,ncol)     ! net shortwave flux
    real(r8) fcns(nlevp,ncol)    ! net clear-sky shortwave flux
    real(r8) fsn200(ncol)        ! fns interpolated to 200 mb
    real(r8) fsn200c(ncol)       ! fcns interpolated to 200 mb
    real(r8) fnl(nlevp,ncol)     ! net longwave flux
    real(r8) fcnl(nlevp,ncol)    ! net clear-sky longwave flux


    real(r8) pbr(nlev,ncol)                    ! Model mid-level pressures (dynes/cm2)
    real(r8) pnm(nlevp,ncol)                   ! Model interface pressures (dynes/cm2)
    real(r8) eccf                              ! Earth/sun distance factor
    real(r8) lwupcgs(ncol)                     ! Upward longwave flux in cgs units

    real(r8) :: dy                             ! Temporary layer pressure thickness
    real(r8) :: tint(nlevp,ncol)               ! Model interface temperature
    real(r8) :: sfac(1:nswbands)               ! time varying scaling factors due to Solar Spectral Irrad at 1 A.U. per band

    real(r8) :: o3(nlev,ncol)                  ! Ozone mass mixing ratio
    real(r8) :: co2(nlev,ncol)                 ! co2   mass mixing ratio
    real(r8) :: co2_col_mean(ncol)             ! co2 column mean mmr
    real(r8) :: sp_hum(nlev,ncol)              ! specific humidity

    ! Aerosol radiative properties
    real(r8) :: aer_tau    (nbndsw,0:nlev,ncol)! aerosol extinction optical depth
    real(r8) :: aer_tau_w  (nbndsw,0:nlev,ncol)! aerosol single scattering albedo * tau
    real(r8) :: aer_tau_w_g(nbndsw,0:nlev,ncol)! aerosol assymetry parameter * w * tau
    real(r8) :: aer_tau_w_f(nbndsw,0:nlev,ncol)! aerosol forward scattered fraction * w * tau
    real(r8) :: aer_lw_abs (nbndlw,nlev,ncol)  ! aerosol absorption optics depth (LW)

    integer :: icall                           ! index through climate/diagnostic radiation calls
    logical :: active_calls(0:N_DIAG)

    type(rrtmg_state_t), pointer :: r_state    ! contains the atm concentratiosn in layers needed for RRTMG

    call get_curr_calday(start_ymd, start_tod, nstep, dtime, calday)

    call zenith(calday, coszrs, ncol, lon, lat )

    cldfsnow_idx = .false.
    if(allocated(pstate_cam%microp_cldfsnow_at_pc_full_level%f))then
        cldfsnow_idx = .true.
        cldfsnow(:,1:ncol) = pstate_cam%microp_cldfsnow_at_pc_full_level%f(old_time_level,:,1:ncol)
    end if
    cld(:,1:ncol) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)
 
    ! Gather night/day column indices.
    Nday = 0
    Nnite = 0
    do i = 1, ncol
       if ( coszrs(i) > 0.0_r8 ) then
          Nday = Nday + 1
          IdxDay(Nday) = i
       else
          Nnite = Nnite + 1
          IdxNite(Nnite) = i
       end if
    end do

    dosw     = radiation_do('sw', nstep)      ! do shortwave heating calc this timestep?
    dolw     = radiation_do('lw', nstep)      ! do longwave heating calc this timestep?
    call_radiation = dosw

    if (dosw .or. dolw) then
       ! construct an RRTMG state object
       r_state => rrtmg_state_create( ncol )

       if (dosw) then
! LiXH Test: oldcldoptics = .false.
           if(oldcldoptics) then
             call ec_ice_optics_sw(ncol, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, .false.)
             call slingo_liq_optics_sw(ncol, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, .false.)
           else
             select case (icecldoptics)
! LiXH Test: icecldoptics = 'mitchell'
             case ('ebertcurry')
                call  ec_ice_optics_sw(ncol, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, .true.)
             case ('mitchell')
                call get_ice_optics_sw(ncol, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)
             case default
                call endrun('iccldoptics must be one either ebertcurry or mitchell')
             end select

             select case (liqcldoptics)
! LiXH Test: liqcldoptics = 'gammadist'
             case ('slingo')
                call slingo_liq_optics_sw(ncol, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, .true.)
             case ('gammadist')
                call get_liquid_optics_sw(ncol, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)
             case default
                call endrun('liqcldoptics must be either slingo or gammadist')
             end select
           end if 

           cld_tau    (:,:,1:ncol) =  liq_tau    (:,:,1:ncol) + ice_tau    (:,:,1:ncol)
           cld_tau_w  (:,:,1:ncol) =  liq_tau_w  (:,:,1:ncol) + ice_tau_w  (:,:,1:ncol)
           cld_tau_w_g(:,:,1:ncol) =  liq_tau_w_g(:,:,1:ncol) + ice_tau_w_g(:,:,1:ncol)
           cld_tau_w_f(:,:,1:ncol) =  liq_tau_w_f(:,:,1:ncol) + ice_tau_w_f(:,:,1:ncol)
 
           if (cldfsnow_idx) then
             ! add in snow
             call get_snow_optics_sw(ncol, snow_tau, snow_tau_w, snow_tau_w_g, snow_tau_w_f)

             do i=1, ncol
                do k=1, nlev
                   cldfprime(k,i)=max(cld(k,i),cldfsnow(k,i))
                   if(cldfprime(k,i) > 0.)then
                      c_cld_tau    (1:nbndsw,k,i) = &
                           (cldfsnow(k,i)*snow_tau    (1:nbndsw,k,i) + cld(k,i)*cld_tau    (1:nbndsw,k,i))/cldfprime(k,i)
                      c_cld_tau_w  (1:nbndsw,k,i) = &
                           (cldfsnow(k,i)*snow_tau_w  (1:nbndsw,k,i) + cld(k,i)*cld_tau_w  (1:nbndsw,k,i))/cldfprime(k,i)
                      c_cld_tau_w_g(1:nbndsw,k,i) = &
                           (cldfsnow(k,i)*snow_tau_w_g(1:nbndsw,k,i) + cld(k,i)*cld_tau_w_g(1:nbndsw,k,i))/cldfprime(k,i)
                      c_cld_tau_w_f(1:nbndsw,k,i) = &
                           (cldfsnow(k,i)*snow_tau_w_f(1:nbndsw,k,i) + cld(k,i)*cld_tau_w_f(1:nbndsw,k,i))/cldfprime(k,i)
                   else
                      c_cld_tau    (1:nbndsw,k,i)= 0._r8
                      c_cld_tau_w  (1:nbndsw,k,i)= 0._r8
                      c_cld_tau_w_g(1:nbndsw,k,i)= 0._r8
                      c_cld_tau_w_f(1:nbndsw,k,i)= 0._r8
                   endif
                enddo
             enddo
           else
             c_cld_tau    (1:nbndsw,:,1:ncol)= cld_tau    (:,:,1:ncol)
             c_cld_tau_w  (1:nbndsw,:,1:ncol)= cld_tau_w  (:,:,1:ncol)
             c_cld_tau_w_g(1:nbndsw,:,1:ncol)= cld_tau_w_g(:,:,1:ncol)
             c_cld_tau_w_f(1:nbndsw,:,1:ncol)= cld_tau_w_f(:,:,1:ncol)
           end if 
       end if
 
       if (dolw) then
          !if (FLG_4DDA_LW) then !linhan
          !   c_cld_lw_abs(1:nbndlw,:,:)= 0._r8 !linhan
          !else !linhan keep the result of (res > 91.0 um)
             if(oldcldoptics) then
                call cloud_rad_props_get_lw(ncol, cld_lw_abs, oldcloud=.true.)
             else
                select case (icecldoptics)
                case ('ebertcurry')
                   call ec_ice_get_rad_props_lw(ncol, ice_lw_abs, oldicewp=.true.)
                case ('mitchell')
                   call ice_cloud_get_rad_props_lw(ncol, ice_lw_abs)
                case default
                   call endrun('iccldoptics must be one either ebertcurry or mitchell')
                end select
                select case (liqcldoptics)
                case ('slingo')
                   call slingo_liq_get_rad_props_lw(ncol, liq_lw_abs, oldliqwp=.true.)
                case ('gammadist')
                   call liquid_cloud_get_rad_props_lw(ncol, liq_lw_abs)
                case default
                   call endrun('liqcldoptics must be either slingo or gammadist')
                end select
                cld_lw_abs(:,:,1:ncol) = liq_lw_abs(:,:,1:ncol) + ice_lw_abs(:,:,1:ncol)
             endif
             if (cldfsnow_idx) then
                ! add in snow
                call snow_cloud_get_rad_props_lw(ncol, snow_lw_abs)
                do i=1, ncol
                   do k=1, nlev
                      cldfprime(k,i)=max(cld(k,i),cldfsnow(k,i))
                      if(cldfprime(k,i) > 0.)then
                         c_cld_lw_abs(1:nbndlw,k,i)= &
                              (cldfsnow(k,i)*snow_lw_abs(1:nbndlw,k,i) + cld(k,i)*cld_lw_abs(1:nbndlw,k,i))/cldfprime(k,i)
                         if (FLG_4DDA_LW) then !linhan only use the abs of snow 
                            c_cld_lw_abs(1:nbndlw,k,i) = snow_lw_abs(1:nbndlw,k,i)*cldfsnow(k,i)/cldfprime(k,i) !linhan
                         endif !linhan
                         if (FLG_4DDA_LW) then !linhan adjust the path 
                            pstate_cam%microp_icswp_at_pc_full_level%f(k,i)=pstate_cam%microp_icswp_at_pc_full_level%f(k,i)*cldfsnow(k,i)/cldfprime(k,i) !linhan
                            pstate_cam%microp_iciwp_at_pc_full_level%f(k,i)=pstate_cam%microp_iciwp_at_pc_full_level%f(k,i)*cld(k,i)/cldfprime(k,i) !linhan
                            pstate_cam%microp_iclwp_at_pc_full_level%f(k,i)=pstate_cam%microp_iclwp_at_pc_full_level%f(k,i)*cld(k,i)/cldfprime(k,i) !linhan
                         endif !linhan
                      else
                         c_cld_lw_abs(1:nbndlw,k,i)= 0._r8
                      endif
                   enddo
                enddo
             else
                c_cld_lw_abs(1:nbndlw,:,1:ncol)=cld_lw_abs(:,:,1:ncol)
             endif
          !endif
       endif
 
       if (.not.cldfsnow_idx) then
          cldfprime(:,1:ncol)=cld(:,1:ncol)
       endif

       ! construct cgs unit reps of pmid and pint and get "eccf" - earthsundistancefactor
       call radinp(ncol, calday, pstate%pressure_at_pc_full_level%f(:,1:ncol),     &
                   pstate%pressure_at_pc_face_level%f(:,1:ncol), pbr, pnm, eccf)

       ! Calculate interface temperatures (following method
       ! used in radtpl for the longwave), using surface upward flux and
       ! stebol constant in mks units
       do i = 1, ncol
          tint(1,i)     = pstate%temp_at_pc_full_level%f(1,i)
          tint(nlevp,i) = sqrt(sqrt(pstate%atm_in_lwup_at_pc_surface%f(i)/stebol))

          !--------------LiXH Test---------------
          !if(pstate%atm_in_lwup_at_pc_surface%f(i)==0._r8 .and. working_mode=='scm' )then
          !    tint(nlevp,i) = tground
          !end if
          !--------------LiXH Test---------------

          do k = 2, nlev
             dy = (log(pstate%pressure_at_pc_face_level%f(k,i)) - log(pstate%pressure_at_pc_full_level%f(k,i))) &
                 /(log(pstate%pressure_at_pc_full_level%f(k-1,i))-log(pstate%pressure_at_pc_full_level%f(k,i)))
             tint(k,i) = pstate%temp_at_pc_full_level%f(k,i)                                                    & 
                        -dy*(pstate%temp_at_pc_full_level%f(k,i) - pstate%temp_at_pc_full_level%f(k-1,i))
          end do
       end do
       ! Solar radiation computation
       if (dosw) then

          call get_variability(sfac)

          ! Get the active climate/diagnostic shortwave calculations
          call rad_cnst_get_call_list(active_calls)

          ! The climate (icall==0) calculation must occur last.
          do icall = N_DIAG, 0, -1

              if (active_calls(icall)) then

                  ! update the concentrations in the RRTMG state object
                  call rrtmg_state_update( ncol, r_state )

                  call aer_rad_props_sw( ncol, icall, nnite, idxnite, &
                                         aer_tau, aer_tau_w, aer_tau_w_g, aer_tau_w_f)

                  call rad_rrtmg_sw( &
                       ncol,         num_rrtmg_levs, r_state,                                  &
                       pstate%pressure_at_pc_full_level%f(:,1:ncol),           cldfprime,      &
                       aer_tau,      aer_tau_w,    aer_tau_w_g,  aer_tau_w_f,                  &
                       eccf,         coszrs,       solin,        sfac,                         &
                       pstate%atm_in_asdir_at_pc_surface%f(1:ncol),                            &
                       pstate%atm_in_asdif_at_pc_surface%f(1:ncol),                            &
                       pstate%atm_in_aldir_at_pc_surface%f(1:ncol),                            &
                       pstate%atm_in_aldif_at_pc_surface%f(1:ncol),                            &
                       pstate_cam%sw_qrs_at_pc_full_level%f(:,1:ncol) ,                        &
                       pstate_cam%sw_qrsc_at_pc_full_level%f(:,1:ncol),                        & 
                       fsnt,         fsntc,        pstate_cam%fsntoa_at_pc_top%f(1:ncol),      &
                       fsutoa,       pstate_cam%fsntoac_at_pc_top%f(1:ncol),                   &
                       fsnirt,       fsnrtc,       fsnirtsq,     fsns,         fsnsc,          &
                       pstate_cam%fswdsc_at_pc_surface%f(1:ncol),                              & 
                       pstate_cam%fswusc_at_pc_surface%f(1:ncol),                              & 
                       pstate_cam%fswdtc_at_pc_top%f(1:ncol),                                  & 
                       pstate_cam%fswutc_at_pc_top%f(1:ncol),                                  & 
                       pstate_cam%fswds_at_pc_surface%f(1:ncol),                               & 
                       pstate_cam%fswus_at_pc_surface%f(1:ncol),                               & 
                       pstate_cam%fswdt_at_pc_top%f(1:ncol),                                   & 
                       pstate_cam%fswut_at_pc_top%f(1:ncol),                                   & 
                       pstate%atm_out_sols_at_pc_surface%f(1:ncol),                            &
                       pstate%atm_out_soll_at_pc_surface%f(1:ncol),                            &
                       pstate%atm_out_solsd_at_pc_surface%f(1:ncol),                           &
                       pstate%atm_out_solld_at_pc_surface%f(1:ncol), fns,      fcns,           &
                       Nday,         Nnite,        IdxDay,       IdxNite,                      &
                       su,           sd,                                                       &
                       c_cld_tau,    c_cld_tau_w,  c_cld_tau_w_g,c_cld_tau_w_f,                &
                       .false., FLG_4DDA_SW)
                       ! LiXH changes:  why can not use optional? 
                       !E_cld_tau=c_cld_tau, E_cld_tau_w=c_cld_tau_w, E_cld_tau_w_g=c_cld_tau_w_g, E_cld_tau_w_f=c_cld_tau_w_f, &
                       !old_convert = .false.)

                  do i=1,ncol
                     pstate_cam%swcf_at_pc_top%f(i)=pstate_cam%fsntoa_at_pc_top%f(i) - pstate_cam%fsntoac_at_pc_top%f(i)
                  end do
!------------LiXH Test--------------
!                  !  Output net fluxes at 200 mb
!                  call vertinterp(ncol, pcols, nlevp, state%pint, 20000._r8, fcns, fsn200c)
!                  call vertinterp(ncol, pcols, nlevp, state%pint, 20000._r8, fns, fsn200)
!
!                  ! Dump shortwave radiation information to history tape buffer (diagnostics)
!                  ftem(:ncol,:nlev) = qrs(:ncol,:nlev)/cpair
!                  call outfld('QRS'//diag(icall),ftem  ,pcols,lchnk)
!                  ftem(:ncol,:nlev) = qrsc(:ncol,:nlev)/cpair
!                  call outfld('QRSC'//diag(icall),ftem  ,pcols,lchnk)
!                  call outfld('SOLIN'//diag(icall),solin ,pcols,lchnk)
!                  call outfld('FSDS'//diag(icall),fsds  ,pcols,lchnk)
!                  call outfld('FSNIRTOA'//diag(icall),fsnirt,pcols,lchnk)
!                  call outfld('FSNRTOAC'//diag(icall),fsnrtc,pcols,lchnk)
!                  call outfld('FSNRTOAS'//diag(icall),fsnirtsq,pcols,lchnk)
!                  call outfld('FSNT'//diag(icall),fsnt  ,pcols,lchnk)
!                  call outfld('FSNS'//diag(icall),fsns  ,pcols,lchnk)
!                  call outfld('FSNTC'//diag(icall),fsntc ,pcols,lchnk)
!                  call outfld('FSNSC'//diag(icall),fsnsc ,pcols,lchnk)
!                  call outfld('FSDSC'//diag(icall),fsdsc ,pcols,lchnk)
!                  call outfld('FSNTOA'//diag(icall),fsntoa,pcols,lchnk)
!                  call outfld('FSUTOA'//diag(icall),fsutoa,pcols,lchnk)
!                  call outfld('FSNTOAC'//diag(icall),fsntoac,pcols,lchnk)
!                  call outfld('SOLS'//diag(icall),cam_out%sols  ,pcols,lchnk)
!                  call outfld('SOLL'//diag(icall),cam_out%soll  ,pcols,lchnk)
!                  call outfld('SOLSD'//diag(icall),cam_out%solsd ,pcols,lchnk)
!                  call outfld('SOLLD'//diag(icall),cam_out%solld ,pcols,lchnk)
!                  call outfld('FSN200'//diag(icall),fsn200,pcols,lchnk)
!                  call outfld('FSN200C'//diag(icall),fsn200c,pcols,lchnk)
!                  call outfld('SWCF'//diag(icall),swcf  ,pcols,lchnk)
!------------LiXH Test--------------

              end if ! (active_calls(icall))
          end do ! icall

!------------LiXH Test--------------
!          ! Output cloud optical depth fields for the visible band
!          tot_icld_vistau(:ncol,:)  = c_cld_tau(idx_sw_diag,:ncol,:)
!          liq_icld_vistau(:ncol,:)  = liq_tau(idx_sw_diag,:ncol,:)
!          ice_icld_vistau(:ncol,:)  = ice_tau(idx_sw_diag,:ncol,:)
!          if (cldfsnow_idx) then
!             snow_icld_vistau(:ncol,:) = snow_tau(idx_sw_diag,:ncol,:)
!          endif
!          ! multiply by total cloud fraction to get gridbox value
!          tot_cld_vistau(:ncol,:) = c_cld_tau(idx_sw_diag,:ncol,:)*cldfprime(:ncol,:)
!
!          ! add fillvalue for night columns
!          do i = 1, Nnite
!              tot_cld_vistau(IdxNite(i),:)   = fillvalue
!              tot_icld_vistau(IdxNite(i),:)  = fillvalue
!              liq_icld_vistau(IdxNite(i),:)  = fillvalue
!              ice_icld_vistau(IdxNite(i),:)  = fillvalue
!              if (cldfsnow_idx) then
!                 snow_icld_vistau(IdxNite(i),:) = fillvalue
!              endif
!          end do
!
!          call outfld('TOT_CLD_VISTAU', tot_cld_vistau, pcols, lchnk)       
!          call outfld('TOT_ICLD_VISTAU', tot_icld_vistau, pcols, lchnk)
!          call outfld('LIQ_ICLD_VISTAU', liq_icld_vistau, pcols, lchnk)
!          call outfld('ICE_ICLD_VISTAU', ice_icld_vistau, pcols, lchnk)
!          if (cldfsnow_idx) then
!             call outfld('SNOW_ICLD_VISTAU', snow_icld_vistau, pcols, lchnk)
!          endif
!------------LiXH Test--------------
 
       end if !dosw
 
!------------LiXH Test--------------
       ! Output aerosol mmr
!       call rad_cnst_out(0, state, pbuf)
!------------LiXH Test--------------

       ! Longwave radiation computation

       if (dolw) then
          !
          ! Convert upward longwave flux units to CGS
          !
          do i=1,ncol
             lwupcgs(i) = pstate%atm_in_lwup_at_pc_surface%f(i)*1000._r8
             !--------------LiXH has not completed scm_crm_mode------------->
             !if(trim(working_mode) .eq. 'scm' .and. have_tgi .and. scm_crm_mode) lwupcgs(i) = 1000*stebol*tground**4
             !<-------------LiXH has not completed scm_crm_mode--------------
          end do

          call rad_cnst_get_call_list(active_calls)

          ! The climate (icall==0) calculation must occur last.
          do icall = N_DIAG, 0, -1

              if (active_calls(icall)) then

                  ! update the conctrations in the RRTMG state object
                  call  rrtmg_state_update( ncol, r_state)

                  call aer_rad_props_lw(ncol, icall, aer_lw_abs)
 
                  call rad_rrtmg_lw( &
                       ncol,         num_rrtmg_levs,  r_state,                                   &
                       pstate%pressure_at_pc_full_level%f(:,1:ncol),                             &
                       aer_lw_abs,   cldfprime,       c_cld_lw_abs,                              &
                       pstate_cam%microp_iciwp_at_pc_full_level%f(:,1:ncol),                     &
                       pstate_cam%microp_icswp_at_pc_full_level%f(:,1:ncol),                     &
                       pstate_cam%microp_iclwp_at_pc_full_level%f(:,1:ncol),                     &
                       pstate_cam%microp_rei_at_pc_full_level%f(:,1:ncol),                       &
                       pstate_cam%microp_res_at_pc_full_level%f(:,1:ncol),                       &
                       pstate_cam%microp_rel_at_pc_full_level%f(:,1:ncol),                       &
                       pstate_cam%lw_qrl_at_pc_full_level%f(:,1:ncol),                           &
                       pstate_cam%lw_qrlc_at_pc_full_level%f(:,1:ncol),                          &
                       flns,         flnt,            flnsc,           flntc,                    &
                       pstate_cam%flwut_at_pc_top%f(1:ncol),                                     &
                       pstate_cam%flwdt_at_pc_top%f(1:ncol),                                     &
                       pstate_cam%flwus_at_pc_surface%f(1:ncol),                                 &
                       pstate_cam%flwds_at_pc_surface%f(1:ncol),                                 &  
                       pstate_cam%flwutc_at_pc_top%f(1:ncol),                                    &
                       pstate_cam%flwdtc_at_pc_top%f(1:ncol),                                    &
                       pstate_cam%flwusc_at_pc_surface%f(1:ncol),                                &
                       pstate_cam%flwdsc_at_pc_surface%f(1:ncol),                                &  
                       fnl,          fcnl,            lu,              ld,                       &
                       FLG_4DDA_LW,  coszrs)


                  do i=1,ncol
                     pstate_cam%lwcf_at_pc_top%f(i)=pstate_cam%flwutc_at_pc_top%f(i) - pstate_cam%flwut_at_pc_top%f(i)
                  end do

!------------LiXH Test--------------
                !output
                !  !  Output fluxes at 200 mb
                !  call vertinterp(ncol, pcols, nlevp, state%pint, 20000._r8, fnl, fln200)
                !  call vertinterp(ncol, pcols, nlevp, state%pint, 20000._r8, fcnl, fln200c)
                !  ! Dump longwave radiation information to history tape buffer (diagnostics)
                !  call outfld('QRL'//diag(icall),qrl (:ncol,:)/cpair,ncol,lchnk)
                !  call outfld('QRLC'//diag(icall),qrlc(:ncol,:)/cpair,ncol,lchnk)
                !  call outfld('FLNT'//diag(icall),flnt  ,pcols,lchnk)
                !  call outfld('FLUT'//diag(icall),flut  ,pcols,lchnk)
                !  call outfld('FLUTC'//diag(icall),flutc ,pcols,lchnk)
                !  call outfld('FLNTC'//diag(icall),flntc ,pcols,lchnk)
                !  call outfld('FLNS'//diag(icall),flns  ,pcols,lchnk)
                !  
                !  call outfld('FLDSC'//diag(icall),fldsc ,pcols,lchnk)
                !  call outfld('FLNSC'//diag(icall),flnsc ,pcols,lchnk)
                !  call outfld('LWCF'//diag(icall),lwcf  ,pcols,lchnk)
                !  call outfld('FLN200'//diag(icall),fln200,pcols,lchnk)
                !  call outfld('FLN200C'//diag(icall),fln200c,pcols,lchnk)
!------------LiXH Test--------------

              end if
          end do
       end if   ! dolw

       call rrtmg_state_destroy(r_state)

       ! mji/hirsrtm - Add call to HIRSRTM package
       ! HIRS brightness temperature calculation in 7 infra-red channels and 4 microwave
       ! channels as a diagnostic to compare to TOV/MSU satellite data.
       ! Done if dohirs set to .true. at time step frequency ihirsfq

       if ( dohirs .and. (mod(nstep-1,ihirsfq) .eq. 0) ) then

          do i= 1, ncol
             ts(i) = sqrt(sqrt(pstate%atm_in_lwup_at_pc_surface%f(i)/stebol))
             ! Set oro (land/sea flag) for compatibility with landfrac/icefrac/ocnfrac
             ! oro=0 (sea or ice); oro=1 (land)
             if (pstate%landfrac_at_pc_surface%f(i).ge.0.001) then
                oro(i)=1.
             else
                oro(i)=0.
             endif
             ! Convert pressure from Pa to hPa
             do k = 1, nlev
                pintmb(k,i) = pstate%pressure_at_pc_face_level%f(k,i)*1.e-2_r8        
             end do
             pintmb(nlevp,i) = pstate%pressure_at_pc_face_level%f(nlevp,i)*1.e-2_r8 
          end do
          
          ! Get specific humidity
          sp_hum(:,1:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol)
          ! Get ozone mass mixing ratio.
          o3(:,1:ncol)     = pstate_cam%ghg_at_pc_full_level(o3_idx)%f(:,1:ncol)
          ! Get CO2 mass mixing ratio
          co2(:,1:ncol)    = pstate_cam%ghg_at_pc_full_level(co2_idx)%f(:,1:ncol)

          call calc_col_mean(ncol, co2, co2_col_mean)

          call hirsrtm( ncol   ,pintmb      ,                     &
                        pstate%temp_at_pc_full_level%f(:,1:ncol), &
                        sp_hum ,co2_col_mean,                     &
                        o3     ,ts       ,oro    ,tb_ir  ,britemp )

          !---------------LiXH has not completed----------------->
          !do i = 1, pnb_hirs
          !   call outfld(hirsname(i),tb_ir(1,i),pcols,lchnk)
          !end do
          !do i = 1, pnf_msu
          !   call outfld(msuname(i),britemp(1,i),pcols,lchnk)
          !end do
          !<--------------LiXH has not completed------------------

       end if

       !! initialize and calculate emis
       emis(:,:) = 0._r8
       emis(:,:ncol) = 1._r8 - exp(-cld_lw_abs(rrtmg_lw_cloudsim_band,:,:ncol))
       !---------------LiXH has not completed----------------->
       !call outfld('EMIS', emis, pcols, lchnk)
       !<--------------LiXH has not completed------------------

       !! compute grid-box mean SW and LW snow optical depth for use by COSP
       gb_snow_tau(:,:) = 0._r8
       gb_snow_lw(:,:) = 0._r8
       !if (cldfsnow_idx > 0) then
       if (cldfsnow_idx ) then
          do i=1,ncol
             do k=1,nlev
                if(cldfsnow(k,i) > 0.)then
                   gb_snow_tau(k,i) = snow_tau(rrtmg_sw_cloudsim_band,k,i)*cldfsnow(k,i)
                   gb_snow_lw(k,i) = snow_lw_abs(rrtmg_lw_cloudsim_band,k,i)*cldfsnow(k,i)
                end if
             enddo
          enddo
       end if

!------------LiXH has not completed COSP------------>
!       if (docosp) then
!          !! cosp_cnt referenced for each chunk... cosp_cnt(lchnk)
!          !! advance counter for this timestep
!          cosp_cnt(lchnk) = cosp_cnt(lchnk) + 1
!
!          !! if counter is the same as cosp_nradsteps, run cosp and reset counter
!           if (cosp_nradsteps .eq. cosp_cnt(lchnk)) then
!              !call should be compatible with camrt radiation.F90 interface too, should be with (in),optional
!              ! N.B.: For snow optical properties, the GRID-BOX MEAN shortwave and longwave optical depths are passed.
!          call cospsimulator_intr_run(state,  pbuf, cam_in, emis, coszrs, &
!                   cld_swtau_in=cld_tau(rrtmg_sw_cloudsim_band,:,:),&
!                   snow_tau_in=gb_snow_tau,snow_emis_in=gb_snow_lw)
!              cosp_cnt(lchnk) = 0  !! reset counter
!           end if
!       end if
!<------------LiXH has not completed COSP------------

    else   !  if (dosw .or. dolw) then

       ! convert radiative heating rates from Q*dp to Q for energy conservation
       if (conserve_energy) then
          do k =1 , nlev
             do i = 1, ncol
                pstate_cam%sw_qrs_at_pc_full_level%f(k,i) = pstate_cam%sw_qrs_at_pc_full_level%f(k,i)   &
                                                       /pstate%delp_at_pc_full_level%f(k,i)
                pstate_cam%lw_qrl_at_pc_full_level%f(k,i) = pstate_cam%lw_qrl_at_pc_full_level%f(k,i)   &
                                                       /pstate%delp_at_pc_full_level%f(k,i)
             end do
          end do
       end if

    end if !dosw .or. dolw

    ! re-initialize ptend
    ptend_radiation%tend_s%f = 0._r8

    ptend_radiation%tend_s%f(:,1:ncol) = pstate_cam%sw_qrs_at_pc_full_level%f(:,:ncol)   &
                                       + pstate_cam%lw_qrl_at_pc_full_level%f(:,:ncol)

    net_flx(1:ncol) = fsnt(1:ncol) - fsns(1:ncol) - flnt(1:ncol) + flns(1:ncol)

    if (conserve_energy) then
       do k =1 , nlev
          do i = 1, ncol
             pstate_cam%sw_qrs_at_pc_full_level%f(k,i) = pstate_cam%sw_qrs_at_pc_full_level%f(k,i)  &
                                                    *pstate%delp_at_pc_full_level%f(k,i)
             pstate_cam%lw_qrl_at_pc_full_level%f(k,i) = pstate_cam%lw_qrl_at_pc_full_level%f(k,i)  &
                                                    *pstate%delp_at_pc_full_level%f(k,i)
          end do
       end do
    end if
 
    pstate%atm_out_netsw_at_pc_surface%f(:ncol) = fsns(:ncol)
    pstate%atm_out_fswds_at_pc_surface%f(:ncol) = pstate_cam%fswds_at_pc_surface%f(1:ncol) 
    pstate%atm_out_flwds_at_pc_surface%f(:ncol) = pstate_cam%flwds_at_pc_surface%f(1:ncol) 

    end subroutine radiation_tend


! Purpose: Set latitude and time dependent arrays for input to solar
!          and longwave radiation. Convert model pressures to cgs.
    subroutine radinp(ncol, calday, pmid, pint, pmidrd, pintrd, eccf)
    use grist_zenith,        only: orb_decl
    ! io
    integer, intent(in)   :: ncol               ! number of atmospheric columns
    real(r8), intent(in)  :: calday
    real(r8), intent(in)  :: pmid(nlev,ncol)    ! Pressure at model mid-levels (pascals)
    real(r8), intent(in)  :: pint(nlevp,ncol)   ! Pressure at model interfaces (pascals)

    real(r8), intent(out) :: pmidrd(nlev,ncol)  ! Pressure at mid-levels (dynes/cm*2)
    real(r8), intent(out) :: pintrd(nlevp,ncol) ! Pressure at interfaces (dynes/cm*2)
    real(r8), intent(out) :: eccf               ! Earth-sun distance factor
    ! local
    integer i                ! Longitude loop index
    integer k                ! Vertical loop index
    real(r8) :: delta        ! Solar declination angle

    call orb_decl(calday  ,delta   ,eccf)

    ! Convert pressure from pascals to dynes/cm2
    do k=1, nlev
       do i=1, ncol
          pmidrd(k,i) = pmid(k,i)*10.0_r8
          pintrd(k,i) = pint(k,i)*10.0_r8
       end do
    end do
    do i=1,ncol
       pintrd(nlevp,i) = pint(nlevp,i)*10.0_r8
    end do

    end subroutine radinp


! Compute the column mean mass mixing ratio.  
    subroutine calc_col_mean(ncol, mmr_pointer, mean_value)

    integer,  intent(in)                    :: ncol
    real(r8), dimension(nlev, ncol)         :: mmr_pointer  ! mass mixing ratio (lev)
    real(r8), dimension(ncol), intent(out)  :: mean_value   ! column mean mmr

    integer  :: i, k
    real(r8) :: ptot(ncol)
    !-----------------------------------------------------------------------

    mean_value   = 0.0_r8
    ptot         = 0.0_r8

    do k=1, nlev
       do i=1, ncol
          mean_value(i) = mean_value(i) + mmr_pointer(k,i)*pstate%delp_at_pc_full_level%f(k,i)
          ptot(i)       = ptot(i) + pstate%delp_at_pc_full_level%f(k,i)
       end do
    end do
    do i=1,ncol
       mean_value(i) = mean_value(i) / ptot(i)
    end do

    end subroutine calc_col_mean


 end module grist_radiation
