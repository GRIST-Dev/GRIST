!========================================================================================
!
!  Created by LiXiaohan on 19/8/13, adopted from CAM5.
!
!  Purpose: Provides greenhouse gas (ghg) values at the Earth's surface.
!           These values may be time dependent.
!========================================================================================

 module chem_surfvals

    use grist_constants,              only: r8
    use grist_nml_module,             only: ntracer, start_tod, start_ymd
    use grist_handle_error,           only: endrun
    use grist_physics_data_structure, only: pstate
    use grist_mpi
 
    implicit none
    private                   ! Make default access private
    save

! Public methods
    public :: read_nml_chem_surfvals, &! read namelist input 
              chem_surfvals_init,     &! initialize options that depend on namelist input
              chem_surfvals_get,      &! return surface values for: CO2VMR, CO2MMR, CH4VMR
                                       ! N2OVMR, F11VMR, and F12VMR
              chem_surfvals_co2_rad    ! return co2 for radiation
!   public ::&
!      chem_surfvals_readnl,  &! read namelist input
!      chem_surfvals_init,    &! initialize options that depend on namelist input
!      chem_surfvals_set,     &! set ghg surface values when ramp option is on
!      chem_surfvals_get,     &! return surface values for: CO2VMR, CO2MMR, CH4VMR
!                              ! N2OVMR, F11VMR, and F12VMR
!      chem_surfvals_co2_rad   ! return co2 for radiation

    public :: flbc_list

    ! Default values for namelist variables -- now set by build-namelist
    real(r8) :: o2mmr = .23143_r8               ! o2 mass mixing ratio
    real(r8) :: co2vmr_rad = -1.0_r8            ! co2 vmr override for radiation
    real(r8) :: co2vmr = -1.0_r8                ! co2   volume mixing ratio 
    real(r8) :: n2ovmr = -1.0_r8                ! n2o   volume mixing ratio 
    real(r8) :: ch4vmr = -1.0_r8                ! ch4   volume mixing ratio 
    real(r8) :: f11vmr = -1.0_r8                ! cfc11 volume mixing ratio 
    real(r8) :: f12vmr = -1.0_r8                ! cfc12 volume mixing ratio 
    character(len=16) :: scenario_ghg = 'FIXED' ! 'FIXED','RAMPED' or 'RAMP_CO2_ONLY'
    integer  :: rampYear_ghg = 0                ! ramped gases fixed at this year (if > 0)
    character(len=256) :: bndtvghg = ' '        ! filename for ramped data
    integer  :: ramp_co2_start_ymd = 0          ! start date for co2 ramping (yyyymmdd)
    real(r8) :: ramp_co2_annual_rate = 1.0_r8   ! % amount of co2 ramping per yr; default is 1% 
    real(r8) :: ramp_co2_cap = -9999.0_r8       ! co2 ramp cap if rate>0, floor otherwise 
                                                ! as multiple or fraction of inital value
                                                ! ex. 4.0 => cap at 4x initial co2 setting 
    integer  :: ghg_yearStart_model = 0         ! model start year
    integer  :: ghg_yearStart_data  = 0         ! data  start year   

    logical  :: ghg_use_calendar                ! true => data year = model year
    logical  :: doRamp_ghg                      ! true => turn on ramping for ghg
    logical  :: ramp_just_co2                   ! true => ramping to be done just for co2 and not other ghg's
    integer  :: fixYear_ghg                     ! year at which Ramped gases are fixed
    integer  :: co2_start                       ! date at which co2 begins ramping
    real(r8) :: co2_daily_factor                ! daily multiplier to achieve annual rate of co2 ramp
    real(r8) :: co2_limit                       ! value of co2vmr where ramping ends
    real(r8) :: co2_base                        ! initial co2 volume mixing ratio, before any ramping
    integer :: ntim = -1                        ! number of yearly data values
    integer,  allocatable :: yrdata(:)          ! yearly data values
    real(r8), allocatable :: co2(:)             ! co2 mixing ratios in ppmv 
    real(r8), allocatable :: ch4(:)             ! ppbv
    real(r8), allocatable :: n2o(:)             ! ppbv
    real(r8), allocatable :: f11(:)             ! pptv
    real(r8), allocatable :: f12(:)             ! pptv
    real(r8), allocatable :: adj(:)             ! unitless adjustment factor for f11 & f12
    
    ! fixed lower boundary 
    character(len=256) :: flbc_file = ' '
    character(len=16), allocatable  :: flbc_list(:)
!-------------LiXH has not completed waccm/cam-chem-------------->   
!    type(time_ramp)    :: flbc_timing     != time_ramp( "CYCLICAL",  19970101, 0 )
!<------------LiXH has not completed waccm/cam-chem---------------   

 contains

!Purpose :  Read chem_surfvals_nl namelist group.
    subroutine read_nml_chem_surfvals(nlfile)

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    character(len=*), parameter :: subname = 'read_nml_chem_surfvals'
    ! Local variables
    integer :: unitn, ierr, i
    
    character(len=8)   :: flbc_type = 'CYCLICAL'     ! 'CYCLICAL' | 'SERIAL' | 'FIXED'
    integer            :: flbc_cycle_yr = 0
    integer            :: flbc_fixed_ymd = 0
    integer            :: flbc_fixed_tod = 0

    namelist /chem_surfvals_nl/ co2vmr, n2ovmr, ch4vmr, f11vmr, f12vmr, &
                                co2vmr_rad, scenario_ghg, rampyear_ghg, bndtvghg, &
                                ramp_co2_start_ymd, ramp_co2_annual_rate, ramp_co2_cap, &
                                ghg_yearStart_model, ghg_yearStart_data
    ! waccm/cam-chem naemlist
    namelist /chem_surfvals_nl/ flbc_type, flbc_cycle_yr, flbc_fixed_ymd, flbc_fixed_tod, flbc_list, flbc_file

    allocate(flbc_list(ntracer))
    flbc_list = ''

    unitn = 111
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, chem_surfvals_nl, iostat=ierr)
    if (ierr /= 0) call endrun('read_nml_chem_surfvals: ERROR reading namelist')
    close(unitn)

!-------------LiXH has not completed waccm/cam-chem-------------->   
!   flbc_timing%type      = flbc_type
!   flbc_timing%cycle_yr  = flbc_cycle_yr
!   flbc_timing%fixed_ymd = flbc_fixed_ymd
!   flbc_timing%fixed_tod = flbc_fixed_tod
!<------------LiXH has not completed waccm/cam-chem---------------   

    if ( len_trim(bndtvghg) > 0 .and. len_trim(flbc_file) > 0 ) then
       call endrun('read_nml_chem_surfvals: Cannot specify both bndtvghg and flbc_file ')
    endif

    if (co2vmr_rad > 0._r8) then
       if (mpi_rank()==0) print*, trim(subname)//': co2vmr_rad override is set to ', co2vmr_rad
    end if

    end subroutine read_nml_chem_surfvals



! Purpose: Initialize the ramp options that are controlled by namelist input.
!          Set surface values at initial time.
    subroutine chem_surfvals_init()
    use grist_time_manager,     only: get_start_date
    
    ! local
    integer :: yr, mon, day, ncsec

    if (scenario_ghg == 'FIXED') then
       doRamp_ghg    = .false.
       ramp_just_co2 = .false.
       if (mpi_rank()==0) print*,'chem_surfvals_init: ghg surface values are fixed as follows'

    else if (scenario_ghg == 'RAMPED') then
       doRamp_ghg = .true.
       ramp_just_co2 = .false.
       !----------LiXH has not completed ramp read--------->
       call endrun ('chem_surfvals_init: can not do ghg_ramp_read!')
       !call ghg_ramp_read
       !<---------LiXH has not completed ramp read----------

       fixYear_ghg = rampYear_ghg     ! set private member to namelist var
       if (mpi_rank()==0) then
          if ( fixYear_ghg > 0 ) then
             print*, '  FIXED values from year ',fixYear_ghg
          else
             print*, '  RAMPED values initialized to'
          end if
       end if
       !----------LiXH has not completed chem_surfvals_set--------->
       call endrun ('chem_surfvals_init: can not call chem_surfvals_set')
       !call chem_surfvals_set()
       !<---------LiXH has not completed chem_surfvals_set----------

    else if (scenario_ghg == 'RAMP_CO2_ONLY') then
       if(ramp_co2_start_ymd == 0) then
          ! by default start the ramp at the initial run time
          call get_start_date(start_ymd, start_tod, yr, mon, day, ncsec)
          ramp_co2_start_ymd = yr*10000 + mon*100 + day
       end if
       co2_start = ramp_co2_start_ymd

       if(ramp_co2_annual_rate <= -100.0_r8) then
          print*, 'RAMP_CO2:  invalid ramp_co2_annual_rate= ',ramp_co2_annual_rate
          call endrun ('chem_surfvals_init: RAMP_CO2_ANNUAL_RATE must be greater than -100.0')
       end if

       doRamp_ghg = .true.
       ramp_just_co2 = .true.
       co2_base = co2vmr        ! save initial setting 
       if (mpi_rank()==0) print*, 'RAMPED values initialized to'

       co2_daily_factor = (ramp_co2_annual_rate*0.01_r8+1.0_r8)**(1.0_r8/365.0_r8)

       if(ramp_co2_cap > 0.0_r8) then  
          co2_limit = ramp_co2_cap * co2_base
       else                                  ! if no cap/floor specified, provide default
          if(ramp_co2_annual_rate < 0.0_r8) then
             co2_limit = 0.0_r8
          else
             co2_limit = transfer( Z'7F800000',co2_limit)
          end if
       end if
       if((ramp_co2_annual_rate<0.0_r8 .and. co2_limit>co2_base) .or. &
          (ramp_co2_annual_rate>0.0_r8 .and. co2_limit<co2_base)) then
          print*, 'RAMP_CO2: ramp_co2_cap is unreachable'
          print*, 'RAMP_CO2: ramp_co2_annual_rate= ',ramp_co2_annual_rate,' ramp_co2_cap= ',ramp_co2_cap
          call endrun('chem_surfvals_init:  ramp_co2_annual_rate and ramp_co2_cap incompatible')
       end if

       !----------LiXH has not completed chem_surfvals_set--------->
       call endrun ('chem_surfvals_init: can not call chem_surfvals_set')
       !call chem_surfvals_set()
       !<---------LiXH has not completed chem_surfvals_set----------
    else
       call endrun ('chem_surfvals_init: input namelist SCENARIO_GHG must be set to either FIXED, RAMPED or RAMP_CO2_ONLY')
    endif

!-------------LiXH has not completed waccm/cam-chem-------------->   
    ! waccm/cam-chem fixed lower boundary conditions
!    call flbc_inti( flbc_file, flbc_list, flbc_timing, co2vmr, ch4vmr, n2ovmr, f11vmr, f12vmr )
!<------------LiXH has not completed waccm/cam-chem---------------   

    if (mpi_rank()==0) then
       print*, '  co2 volume mixing ratio = ',co2vmr
       print*, '  ch4 volume mixing ratio = ',ch4vmr
       print*, '  n2o volume mixing ratio = ',n2ovmr
       print*, '  f11 volume mixing ratio = ',f11vmr
       print*, '  f12 volume mixing ratio = ',f12vmr
    end if

    end subroutine chem_surfvals_init


    function chem_surfvals_get(name)
    use grist_constants,    only: mwdry, mwco2

    character(len=*), intent(in) :: name
    ! local
    real(r8) :: rmwco2 
    real(r8) :: chem_surfvals_get

    rmwco2 = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
    select case (name)
    case ('CO2VMR')
       chem_surfvals_get = co2vmr
    case ('CO2MMR')
       chem_surfvals_get = rmwco2 * co2vmr
    case ('N2OVMR')
       chem_surfvals_get = n2ovmr
    case ('CH4VMR')
       chem_surfvals_get = ch4vmr
    case ('F11VMR')
       chem_surfvals_get = f11vmr
    case ('F12VMR')
       chem_surfvals_get = f12vmr
    case ('O2MMR')
       chem_surfvals_get = o2mmr
    case default
       call endrun('chem_surfvals_get does not know name')
    end select

    end function chem_surfvals_get


! Purpose : Return the value of CO2 (as mmr) that is radiatively active.
!           This method is used by ghg_data to set the prescribed value of CO2 in
!           the physics buffer.  If the user has set the co2vmr_rad namelist
!           variable then that value will override either the value set by the
!           co2vmr namelist variable, or the values time interpolated from a
!           dataset.
!           This method is also used by cam_history to write the radiatively active
!           CO2 to the history file.  The optional argument allows returning the
!           value as vmr.
    function chem_surfvals_co2_rad(vmr_in)
    use grist_constants,    only: mwdry, mwco2

    ! io
    logical, intent(in), optional :: vmr_in  ! return CO2 as vmr
    ! Return value
    real(r8) :: chem_surfvals_co2_rad

    ! Local variables
    real(r8) :: convert_vmr      ! convert vmr to desired output

    ! by default convert vmr to mmr
    convert_vmr = mwco2/mwdry    ! ratio of molecular weights of co2 to dry air
    if (present(vmr_in)) then
       ! if request return vmr
       if (vmr_in) convert_vmr = 1.0_r8
    end if

    if (co2vmr_rad > 0._r8) then
       chem_surfvals_co2_rad = convert_vmr * co2vmr_rad
    else                           
       chem_surfvals_co2_rad = convert_vmr * co2vmr     
    end if

    end function chem_surfvals_co2_rad


 end module chem_surfvals
