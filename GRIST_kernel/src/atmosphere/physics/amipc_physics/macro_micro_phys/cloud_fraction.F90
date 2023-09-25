!======================================================
!
!  Created by LiXiaohan on 19/7/13.
!  Cloud fraction parameterization.
!======================================================


  module cloud_fraction
      use grist_constants,                    only: i4, r8, tmelt, cappa,   &
                                                    gravity,rdry
      use grist_handle_error,                 only: endrun
      use phys_control,                       only: phys_getopts
      use grist_mpi
 
  implicit none
  private
  save

  ! Public interfaces
  public ::  read_nml_cldfrc,   &! read cldfrc_nl namelist
             cldfrc_init,       &! Inititialization of cloud_fraction run-time parameters
             end_of_cldfrc,     &! destract pstate 
             cldfrc_getparams,  &! public access of tuning parameters
             cldfrc,            &! Computation of cloud fraction
             cldfrc_fice,       &! Calculate fraction of condensate in ice phase (radiation partitioning)
             trop_cloud_top_lev,&
             clim_modal_aero_top_lev


  ! Private data
  real(r8), parameter :: unset_r8 = huge(1.0_r8)

  ! Top level
  integer :: top_lev = 1

  integer :: pver, pverp

  ! Namelist variables
  logical  :: cldfrc_freeze_dry                   ! switch for Vavrus correction
  logical  :: cldfrc_ice                          ! switch to compute ice cloud fraction
  real(r8) :: cldfrc_rhminl = unset_r8            ! minimum rh for low stable clouds
  real(r8) :: cldfrc_rhminl_adj_land = unset_r8   ! rhminl adjustment for snowfree land
  real(r8) :: cldfrc_rhminh = unset_r8            ! minimum rh for high stable clouds
  real(r8) :: cldfrc_sh1    = unset_r8            ! parameter for shallow convection cloud fraction
  real(r8) :: cldfrc_sh2    = unset_r8            ! parameter for shallow convection cloud fraction
  real(r8) :: cldfrc_dp1    = unset_r8            ! parameter for deep convection cloud fraction
  real(r8) :: cldfrc_dp2    = unset_r8            ! parameter for deep convection cloud fraction
  real(r8) :: cldfrc_premit = unset_r8            ! top pressure bound for mid level cloud
  real(r8) :: cldfrc_premib  = unset_r8           ! bottom pressure bound for mid level cloud
  integer  :: cldfrc_iceopt                       ! option for ice cloud closure
                                                  ! 1=wang & sassen 2=schiller (iciwc)
                                                  ! 3=wood & field, 4=Wilson (based on smith)
  real(r8) :: cldfrc_icecrit = unset_r8           ! Critical RH for ice clouds in Wilson & Ballard closure (smaller = more ice clouds)
  real(r8) :: trop_cloud_top_press = 0._r8        ! Pressure used to set troposphere cloud physics top (Pa)
  integer  :: trop_cloud_top_lev                  ! Top level for troposphere cloud physics
  real(r8) :: clim_modal_aero_top_press = 0._r8   ! Pressure used to set MAM process top (Pa)
  integer  :: clim_modal_aero_top_lev             ! Top level for MAM processes that impact climate

  real(r8) :: rhminl             ! set from namelist input cldfrc_rhminl
  real(r8) :: rhminl_adj_land    ! set from namelist input cldfrc_rhminl_adj_land
  real(r8) :: rhminh             ! set from namelist input cldfrc_rhminh
  real(r8) :: sh1, sh2           ! set from namelist input cldfrc_sh1, cldfrc_sh2
  real(r8) :: dp1,dp2            ! set from namelist input cldfrc_dp1, cldfrc_dp2
  real(r8) :: premit             ! set from namelist input cldfrc_premit
  real(r8) :: premib             ! set from namelist input cldfrc_premib
  integer  :: iceopt             ! set from namelist input cldfrc_iceopt
  real(r8) :: icecrit            ! set from namelist input cldfrc_icecrit

  ! query interfaces for scheme settings
  character(len=16) :: shallow_scheme, eddy_scheme, macrop_scheme
 
  ! constants
  real(r8), parameter :: pnot = 1.e5_r8         ! reference pressure
  real(r8), parameter :: lapse = 6.5e-3_r8      ! U.S. Standard Atmosphere lapse rate
  real(r8), parameter :: pretop = 1.0e2_r8      ! pressure bounding high cloud

  integer count

  logical :: inversion_cld_off    ! Turns off stratification-based cld frc

  integer :: k700   ! model level nearest 700 mb


  contains

    subroutine read_nml_cldfrc(nlfile)

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr

    namelist /cldfrc_nl/ cldfrc_dp1,     cldfrc_dp2,      cldfrc_freeze_dry,  &
                         cldfrc_ice,     cldfrc_icecrit,  cldfrc_iceopt,      &
                         cldfrc_premib,  cldfrc_premit,   cldfrc_rhminh,      &
                         cldfrc_rhminl,  cldfrc_rhminl_adj_land,              &
                         cldfrc_sh1,     cldfrc_sh2,                          &
                         trop_cloud_top_press, clim_modal_aero_top_press

    unitn = 111
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, cldfrc_nl, iostat=ierr)
    if (ierr /= 0) call endrun('error reading cldfrc namelist')
    close(unitn)

    rhminl = cldfrc_rhminl
    rhminl_adj_land = cldfrc_rhminl_adj_land
    rhminh = cldfrc_rhminh
    sh1    = cldfrc_sh1
    sh2    = cldfrc_sh2
    dp1    = cldfrc_dp1
    dp2    = cldfrc_dp2
    premit = cldfrc_premit
    premib  = cldfrc_premib
    iceopt  = cldfrc_iceopt
    icecrit = cldfrc_icecrit

    ! Check that top for modal aerosols is not lower than top for clouds.
    if (clim_modal_aero_top_press > trop_cloud_top_press) &
        call endrun('ERROR: clim_modal_aero_top press must be less than or equal to trop_cloud_top_press.')

    end subroutine read_nml_cldfrc


! Purpose: Return cldfrc tuning parameters
    subroutine cldfrc_getparams(rhminl_out, rhminl_adj_land_out, rhminh_out,  premit_out, &
                                premib_out, iceopt_out,          icecrit_out)

    real(r8),          intent(out), optional :: rhminl_out
    real(r8),          intent(out), optional :: rhminl_adj_land_out
    real(r8),          intent(out), optional :: rhminh_out
    real(r8),          intent(out), optional :: premit_out
    real(r8),          intent(out), optional :: premib_out
    integer,           intent(out), optional :: iceopt_out
    real(r8),          intent(out), optional :: icecrit_out

    if ( present(rhminl_out) )      rhminl_out = rhminl
    if ( present(rhminl_adj_land_out) ) rhminl_adj_land_out = rhminl_adj_land
    if ( present(rhminh_out) )      rhminh_out = rhminh
    if ( present(premit_out) )      premit_out = premit
    if ( present(premib_out) )      premib_out  = premib
    if ( present(iceopt_out) )      iceopt_out  = iceopt
    if ( present(icecrit_out) )     icecrit_out = icecrit

    end subroutine cldfrc_getparams

! Purpose: Initialize cloud fraction run-time parameters
    subroutine cldfrc_init(pver_in, pverp_in, ncol)
       use grist_constants,                    only: p00
       use grist_hpe_constants,                only: eta_full
       use grist_PhysC_data_structure,          only: pstate_cam
 
!io
    integer, intent(in)  :: pver_in, pverp_in, ncol

    integer  :: k
    real(r8) :: pref_mid(pver_in)
    !-----------------------------------------------------------------------------
 
    pver  = pver_in
    pverp = pverp_in

    call phys_getopts(shallow_scheme_out = shallow_scheme ,&
                      eddy_scheme_out    = eddy_scheme    ,&
                      macrop_scheme_out  = macrop_scheme  )
 
    ! trop_cloud_top_lev and clim_modal_aero_top_lev are set in ref_pres.F90 in CAM, 
    ! LiXH moves them here.
    pref_mid(:) = eta_full(:)*p00
    ! Find level corresponding to the top of troposphere clouds.
    trop_cloud_top_lev = minloc(abs(pref_mid - trop_cloud_top_press), 1)
    ! Find level corresponding to the top for MAM processes.
    clim_modal_aero_top_lev = minloc(abs(pref_mid - clim_modal_aero_top_press),1)

    ! Limit CAM5 cloud physics to below top cloud level.
    top_lev = trop_cloud_top_lev
 
 
    ! Turn off inversion_cld if any UW PBL scheme is being used
    if ( (eddy_scheme .eq. 'diag_TKE' ) .or. (shallow_scheme .eq.  'UW' ) .or. (shallow_scheme .eq.  'double_plume' )) then
       inversion_cld_off = .true.
    else
       inversion_cld_off = .false.
    endif
 
    if(mpi_rank()==0)then
       print*,'tuning parameters cldfrc_init: inversion_cld_off',inversion_cld_off
       print*,'tuning parameters cldfrc_init: dp1',dp1,'dp2',dp2,'sh1',sh1,'sh2',sh2
       if (shallow_scheme .ne. 'UW' .and. shallow_scheme .ne. 'double_plume') then
          print*,'tuning parameters cldfrc_init: rhminl',rhminl,'rhminl_adj_land',rhminl_adj_land, &
                        'rhminh',rhminh,'premit',premit,'premib',premib
          print*,'tuning parameters cldfrc_init: iceopt',iceopt,'icecrit',icecrit
       endif
    endif
 
    if (pref_mid(top_lev) > 7.e4_r8) &
         call endrun ('cldfrc_init: model levels bracketing 700 mb not found')
 
    ! Find vertical level nearest 700 mb.
    k700 = minloc(abs(pref_mid(top_lev:pver) - 7.e4_r8), 1)
 
    if(mpi_rank()==0)then
       print*,'cldfrc_init: model level nearest 700 mb is',k700,'which is',pref_mid(k700),'pascals'
    end if

! Initialize pstate
    allocate(pstate_cam%cld_sh_frac_at_pc_full_level%f(pver,ncol))
    allocate(pstate_cam%cld_dp_frac_at_pc_full_level%f(pver,ncol))

    pstate_cam%cld_sh_frac_at_pc_full_level%pos   = 0
    pstate_cam%cld_dp_frac_at_pc_full_level%pos   = 0
 
    pstate_cam%cld_sh_frac_at_pc_full_level%f     = 0._r8
    pstate_cam%cld_dp_frac_at_pc_full_level%f     = 0._r8

    end subroutine cldfrc_init


    subroutine end_of_cldfrc

    use grist_PhysC_data_structure,   only: pstate_cam

    if(allocated(pstate_cam%cld_sh_frac_at_pc_full_level%f))      &
      deallocate(pstate_cam%cld_sh_frac_at_pc_full_level%f)

    if(allocated(pstate_cam%cld_dp_frac_at_pc_full_level%f))      &
      deallocate(pstate_cam%cld_dp_frac_at_pc_full_level%f)

    end subroutine end_of_cldfrc


! Purpose: Compute cumulus fraction (concld)
    subroutine cldfrc(ncol    ,                                               &
                      shfrc   ,use_shfrc,                                     &
                      cmfmc   ,cmfmc2   ,concld  )

    use grist_PhysC_data_structure,          only: pstate_cam

! io
    integer, intent(in) :: ncol                   ! number of atmospheric columns
    
    real(r8), intent(in) :: cmfmc(pverp,ncol)    ! convective mass flux--m sub c
    real(r8), intent(in) :: cmfmc2(pverp,ncol)   ! shallow convective mass flux--m sub c
    real(r8), intent(in) :: shfrc(pver,ncol)     ! cloud fraction from convect_shallow
    logical,  intent(in)  :: use_shfrc

    real(r8), intent(out) :: concld(pver,ncol)   ! convective cloud cover

! local


    integer i, k           ! column, level indices

    concld   = 0._r8

    ! Estimate of local convective cloud cover based on convective mass flux (ZM)
    ! Modify local large-scale relative humidity to account for presence of 
    ! convective cloud when evaluating relative humidity based layered cloud amount
    
    do k=top_lev,pver-1
       do i=1,ncol
          if ( .not. use_shfrc ) then
             pstate_cam%cld_sh_frac_at_pc_full_level%f(k,i) = max(0.0_r8,min(sh1*log(1.0_r8+sh2*cmfmc2(k+1,i)),0.30_r8))
          else
             pstate_cam%cld_sh_frac_at_pc_full_level%f(k,i) = shfrc(k,i)
          endif
          !-----------LiXH modified------------->
          !Double Plume scheme contains deep and shallow convections, cld_dp has been calculated in its module, LiXH.
          if (shallow_scheme .ne. 'double_plume')then
              pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i) = max(0.0_r8,min(dp1*log(1.0_r8+dp2*(cmfmc(k+1,i)-cmfmc2(k+1,i))),0.60_r8))
          else
              !pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i) = max(max(0.0_r8,min(dp1*log(1.0_r8+dp2*(cmfmc(k+1,i)-cmfmc2(k+1,i))),0.60_r8)),pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i))
              pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i) = pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i)
          end if
          !<----------LiXH modified--------------

          concld(k,i) = min(pstate_cam%cld_sh_frac_at_pc_full_level%f(k,i)+pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i),0.80_r8)
       end do
    end do

    return
    end subroutine cldfrc


! Purpose: Compute the fraction of the total cloud water which is in ice phase.
!          The fraction depends on temperature only. 
!          This is the form that was used for radiation, the code came from cldefr originally
! 
! Author: B. A. Boville Sept 10, 2002
! modified: PJR 3/13/03 (added fsnow to ascribe snow production for convection )
    subroutine cldfrc_fice(ncol, t, fice, fsnow)
! io
    integer,  intent(in)  :: ncol                 ! number of active columns
    real(r8), intent(in)  :: t(pver,ncol)         ! temperature

    real(r8), intent(out) :: fice(pver,ncol)      ! Fractional ice content within cloud
    real(r8), intent(out) :: fsnow(pver,ncol)     ! Fractional snow content for convection

! Local variables
    real(r8) :: tmax_fice                         ! max temperature for cloud ice formation
    real(r8) :: tmin_fice                         ! min temperature for cloud ice formation
    real(r8) :: tmax_fsnow                        ! max temperature for transition to convective snow
    real(r8) :: tmin_fsnow                        ! min temperature for transition to convective snow
    integer :: i,k                                ! loop indexes

!-----------------------------------------------------------------------

    tmax_fice = tmelt - 10._r8        ! max temperature for cloud ice formation
    tmin_fice = tmax_fice - 30._r8    ! min temperature for cloud ice formation
    tmax_fsnow = tmelt                ! max temperature for transition to convective snow
    tmin_fsnow = tmelt - 5._r8        ! min temperature for transition to convective snow

    fice(:top_lev-1,:) = 0._r8
    fsnow(:top_lev-1,:) = 0._r8

! Define fractional amount of cloud that is ice
    do k=top_lev,pver
       do i=1,ncol

! If warmer than tmax then water phase
          if (t(k,i) > tmax_fice) then
             fice(k,i) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(k,i) < tmin_fice) then
             fice(k,i) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fice(k,i) =(tmax_fice - t(k,i)) / (tmax_fice - tmin_fice)
          end if

! snow fraction partitioning

! If warmer than tmax then water phase
          if (t(k,i) > tmax_fsnow) then
             fsnow(k,i) = 0.0_r8

! If colder than tmin then ice phase
          else if (t(k,i) < tmin_fsnow) then
             fsnow(k,i) = 1.0_r8

! Otherwise mixed phase, with ice fraction decreasing linearly from tmin to tmax
          else 
             fsnow(k,i) =(tmax_fsnow - t(k,i)) / (tmax_fsnow - tmin_fsnow)
          end if

       end do
    end do

    end subroutine cldfrc_fice

  end module cloud_fraction
