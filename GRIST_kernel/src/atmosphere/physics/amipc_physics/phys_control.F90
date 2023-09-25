!======================================================
!
!  Created by LiXiaohan on 19/7/13.
!  Provides a control interface to CAM physics packages
!======================================================

  module phys_control
    use grist_handle_error,                 only: endrun

    implicit none
    private
    save
    
    public :: read_nml_phys_ctl,   &! read namelist from file
              phys_getopts,        &! generic query method
              phys_deepconv_pbl,   &! return true if deep convection is allowed in the PBL
              lsm_scheme            ! option for LSM  !--cheyz 20200419
    
    ! private
    character(len=16), parameter :: unset_str = 'UNSET'
    integer,           parameter :: unset_int = huge(1)
    
    ! Namelist variables:
    character(len=16) :: deep_scheme          = unset_str  ! deep convection package
    character(len=16) :: shallow_scheme       = unset_str  ! shallow convection package
    character(len=16) :: eddy_scheme          = unset_str  ! vertical diffusion package
    character(len=16) :: microp_scheme        = unset_str  ! microphysics package
    character(len=16) :: macrop_scheme        = unset_str  ! macrophysics package
    character(len=16) :: radiation_scheme     = unset_str  ! radiation package
    integer           :: conv_water_in_rad    = unset_int
    logical           :: state_debug_checks   = .false.    ! check physics state  
    logical           :: prog_modal_aero      = .false.    ! check physics state  

    character(len=16) :: lsm_scheme           = unset_str  ! land surface model---cheyz 
 
    contains

    subroutine read_nml_phys_ctl(nlfile)

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    ! Local variables
    integer :: unitn, ierr

    namelist /phys_ctl_nl/deep_scheme, shallow_scheme, eddy_scheme, microp_scheme,  &
                          macrop_scheme, radiation_scheme, conv_water_in_rad,       &
                          state_debug_checks,prog_modal_aero,                       &
                          lsm_scheme  !--cheyz

    unitn = 111
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, phys_ctl_nl, iostat=ierr)
    if (ierr /= 0) call endrun("error reading cldfrc namelist")
    close(unitn)

    end subroutine read_nml_phys_ctl


! Purpose: Return runtime settings of physics package
    subroutine phys_getopts(deep_scheme_out, shallow_scheme_out, eddy_scheme_out, microp_scheme_out, &
                            radiation_scheme_out, macrop_scheme_out, conv_water_in_rad_out,          &
                            prog_modal_aero_out,state_debug_checks_out)

   character(len=16), intent(out), optional :: deep_scheme_out
   character(len=16), intent(out), optional :: shallow_scheme_out
   character(len=16), intent(out), optional :: eddy_scheme_out
   character(len=16), intent(out), optional :: microp_scheme_out
   character(len=16), intent(out), optional :: radiation_scheme_out
   character(len=16), intent(out), optional :: macrop_scheme_out
   integer,           intent(out), optional :: conv_water_in_rad_out
   logical,           intent(out), optional :: prog_modal_aero_out
   logical,           intent(out), optional :: state_debug_checks_out

   if ( present(deep_scheme_out         ) ) deep_scheme_out          = deep_scheme
   if ( present(shallow_scheme_out      ) ) shallow_scheme_out       = shallow_scheme
   if ( present(eddy_scheme_out         ) ) eddy_scheme_out          = eddy_scheme
   if ( present(microp_scheme_out       ) ) microp_scheme_out        = microp_scheme
   if ( present(radiation_scheme_out    ) ) radiation_scheme_out     = radiation_scheme
   if ( present(macrop_scheme_out       ) ) macrop_scheme_out        = macrop_scheme
   if ( present(conv_water_in_rad_out   ) ) conv_water_in_rad_out    = conv_water_in_rad
   if ( present(prog_modal_aero_out     ) ) prog_modal_aero_out      = prog_modal_aero
   if ( present(state_debug_checks_out  ) ) state_debug_checks_out   = state_debug_checks

   end subroutine phys_getopts


! Purpose: Return true if deep convection is allowed in the PBL
   function phys_deepconv_pbl()

   logical phys_deepconv_pbl

   ! Don't allow deep convection in PBL if running UW PBL scheme
   if ( (eddy_scheme .eq. 'diag_TKE' ) .or. (shallow_scheme .eq. 'UW' ) ) then
      phys_deepconv_pbl = .true.
   else
      phys_deepconv_pbl = .false.
   endif

   return

   end function phys_deepconv_pbl


 end module phys_control
