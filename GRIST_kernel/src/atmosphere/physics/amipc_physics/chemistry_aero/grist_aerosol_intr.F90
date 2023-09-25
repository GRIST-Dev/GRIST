module grist_aerosol_intr

    use grist_constants,                    only: r8
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_nml_module,                   only: nlev, ntracer
    use grist_rad_constituents,             only: rad_cnst_get_info
    use phys_control,                       only: phys_getopts
    use grist_modal_aero_wateruptake,       only: modal_aero_wateruptake_dr
    use grist_modal_aero_calcsize,          only: modal_aero_calcsize_diag
    use grist_handle_error,                 only: endrun

implicit none
private
save


public ::   read_nml_aerosol,     &
            aerosol_init,         &    
            aerosol_wet_intr,     &
            aerosol_drydep_intr,  &
            end_of_grist_aerosol

 ! prognostic bulk sea salt
logical :: progsslt = .false.

! prognostic bulk dust
logical :: dust = .false.

! modal aerosols
logical :: prog_modal_aero = .false. ! if true then prognostic modal aerosols are present

! modal aerosols
logical :: clim_modal_aero ! if true then modal aerosols are used in the climate calculation.

! is_any_aerosol is set to .true. to indicate that prognostic dust, sea salt, or
! modal aerosols are present.  It is only used to initialize dry deposition module
logical :: is_any_aerosol = .false.

! Namelist variables
real(r8)      :: dust_emis_fact = -1.e36_r8   ! tuning parameter for dust emissions
character(256) :: soil_erod = 'soil_erod'   ! full pathname for soil erodibility dataset


contains


    subroutine read_nml_aerosol(nlfile)

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    namelist /aerosol_nl/ dust_emis_fact, soil_erod

    ! Read namelist
    unitn = 111
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, aerosol_nl, iostat=ierr)
    if (ierr /= 0) call endrun('ERROR reading aerosol_nl namelist')
    close(unitn)

    end subroutine read_nml_aerosol


    subroutine aerosol_init(ncol)
    !io
    integer, intent(in) :: ncol
    !local
    integer :: nmodes


    !---------------LiXH has not completed chemistry------------->
    !do m=1,ndst
    !   dst_idx(m) = get_spc_ndx( dust_names(m) )
    !enddo
    !do m=1,nsst
    !   slt_idx(m) = get_spc_ndx( progseasalts_names(m) )
    !enddo
 
    !dust      = all( dst_idx > 0 )
    !progsslt  = all( slt_idx > 0 )

    dust     = .false.
    progsslt = .false.
    !<--------------LiXH has not completed chemistry--------------

    ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    ! ***N.B.*** The wet radius and water uptake calculations are currently only being done
    !            for the aerosol modes in the climate calc.  This could be extended to compute
    !            these quantities for the modes used in diagnostic radiative calcs as well.
    !
    ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
    ! The modal aerosols can be either prognostic or prescribed.
    call rad_cnst_get_info(0, nmodes=nmodes)
    clim_modal_aero = (nmodes > 0)

    is_any_aerosol = (dust .or. progsslt .or. prog_modal_aero)


    !---------------LiXH has not completed chemistry------------->
    !if (is_any_aerosol .or. carma_do_drydep) then

    !   ! tell camsrfexch to allocate fv & ram1 -- needed by prodsslts and dust
    !   call hub2atm_setopts(aero_dust_in=.true.)

    !end if
    !<--------------LiXH has not completed chemistry--------------


    ! The following fields are diagnosed for either prognostic or prescribed
    ! modal aerosols.  Currently only modes that affect the climate are
    ! accounted for.
    if (clim_modal_aero) then
        allocate(pstate_cam%aerosol_dgnum%f(nmodes,nlev,ncol))
        allocate(pstate_cam%aerosol_dgnumwet%f(nmodes,nlev,ncol))
        allocate(pstate_cam%aerosol_wetdens_ap%f(nmodes,nlev,ncol))
        allocate(pstate_cam%aerosol_qaerwat%f(nmodes,nlev,ncol))

        pstate_cam%aerosol_dgnum%f               = 0._r8
        pstate_cam%aerosol_dgnumwet%f            = 0._r8
        pstate_cam%aerosol_wetdens_ap%f          = 0._r8
        pstate_cam%aerosol_qaerwat%f             = 0._r8

        pstate_cam%aerosol_dgnum%pos             = 0
        pstate_cam%aerosol_dgnumwet%pos          = 0
        pstate_cam%aerosol_wetdens_ap%pos        = 0
        pstate_cam%aerosol_qaerwat%pos           = 0
    end if
 
    allocate(pstate_cam%aerosol_fracis%f(ntracer, nlev, ncol))
    pstate_cam%aerosol_fracis%pos                = 0
    pstate_cam%aerosol_fracis%f                  = 1._r8


    !---------------LiXH has not completed chemistry------------->
    !if (dust .or. prog_modal_aero) then
    !   call dust_initialize(dust_emis_fact, soil_erod)
    !end if
    !
    !if (progsslt .or. prog_modal_aero) then
    !   call progseasalts_initialize
    !end if
    !<--------------LiXH has not completed chemistry--------------

    if (clim_modal_aero) then

       if (.not. prog_modal_aero) then
          ! If climate calculations are affected by prescribed modal aerosols,
          ! the initialization routine for the dry mode radius calculation is called
          ! here.  For prognostic MAM the initialization is called from
          ! modal_aero_initialize
    !---------------LiXH has not completed chemistry------------->
    !      call modal_aero_calcsize_init()
    !<--------------LiXH has not completed chemistry--------------
       endif

       ! not needed in GRIST, LiXH
       !call modal_aero_wateruptake_init()      

    end if

    !---------------LiXH has not completed chemistry------------->
    ! Dry deposition needs to be initialized if any of the aerosols
    ! are running.
    !if (is_any_aerosol) then
    !   call inidrydep(rair, gravit)
    !endif
    !<--------------LiXH has not completed chemistry--------------


    end subroutine aerosol_init


    subroutine end_of_grist_aerosol ()
      if(allocated(pstate_cam%aerosol_fracis%f))                   &
        deallocate(pstate_cam%aerosol_fracis%f)

      if(allocated(pstate_cam%aerosol_dgnum%f))                   &
        deallocate(pstate_cam%aerosol_dgnum%f)

      if(allocated(pstate_cam%aerosol_dgnumwet%f))                   &
        deallocate(pstate_cam%aerosol_dgnumwet%f)

      if(allocated(pstate_cam%aerosol_wetdens_ap%f))                   &
        deallocate(pstate_cam%aerosol_wetdens_ap%f)

      if(allocated(pstate_cam%aerosol_qaerwat%f))                   &
        deallocate(pstate_cam%aerosol_qaerwat%f)

    end subroutine end_of_grist_aerosol


    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Interface to wet processing of all aerosols
    ! 
    ! Method: 
    !  use a modified version of the scavenging parameterization described in
    !     Barth et al, 2000, JGR (sulfur cycle paper)
    !     Rasch et al, 2001, JGR (INDOEX paper)
    ! 
    ! Author: Phil Rasch
    ! 
    !-----------------------------------------------------------------------
    subroutine aerosol_wet_intr (ncol)

    ! Arguments:

    integer,     intent(in)    :: ncol


    !LiXH has not completed mz_aero!!!
    if (clim_modal_aero) then

        call modal_aero_calcsize_diag(ncol)

        call modal_aero_wateruptake_dr(ncol)
    end if

    end subroutine aerosol_wet_intr


    subroutine aerosol_drydep_intr (ncol)

    ! Arguments:

    integer,     intent(in)    :: ncol


    !LiXH has not completed mz_aero!!!

    return
 
    end subroutine aerosol_drydep_intr

end module grist_aerosol_intr
