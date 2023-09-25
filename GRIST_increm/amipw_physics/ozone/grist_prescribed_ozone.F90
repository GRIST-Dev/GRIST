!=======================================================
!
!  Created by LiXiaohan on 20/07/01.
!  manages reading and interpolation of prescribed ozone
!
!=======================================================

 module grist_prescribed_ozone

    use grist_constants,                    only: i4, i8, r8, rad2deg, boltz, mwdry
    use grist_physics_data_structure,       only: pstate
    use grist_PhysW_data_structure,         only: pstate_wrf
    use grist_nml_module,                   only: nlev
    use grist_time_manager,                 only: get_curr_date
    use grist_tracer_data,                  only: trfld, trfile
    use grist_handle_error,                 only: endrun
    use grist_mpi

    implicit none
    private

    public :: prescribed_ozone_adv,         &
              prescribed_ozone_init,        &
              read_nml_prescribed_ozone

    integer            :: ozone_nlev

    type(trfld), allocatable :: fields(:)
    type(trfile)             :: file


    logical            :: has_prescribed_ozone = .false.
    character(len=8), parameter :: ozone_name = 'ozone'
    character(len=16)  :: fld_name = 'ozone'
    character(len=256) :: filename = ' '
    character(len=256) :: filelist = ' '
    character(len=256) :: datapath = ' '
    character(len=32)  :: data_type = 'SERIAL'
    logical            :: rmv_file = .false.
    integer            :: cycle_yr  = 0
    integer            :: fixed_ymd = 0
    integer            :: fixed_tod = 0

 contains

    subroutine read_nml_prescribed_ozone(nlfile)
! io
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
! local
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'prescribed_ozone_readnl'
 
    character(len=16)  :: prescribed_ozone_name
    character(len=256) :: prescribed_ozone_file
    character(len=256) :: prescribed_ozone_filelist
    character(len=256) :: prescribed_ozone_datapath
    character(len=32)  :: prescribed_ozone_type
    logical            :: prescribed_ozone_rmfile
    integer            :: prescribed_ozone_cycle_yr
    integer            :: prescribed_ozone_fixed_ymd
    integer            :: prescribed_ozone_fixed_tod
 
    namelist /prescribed_ozone_nl/ &
       prescribed_ozone_name,      &
       prescribed_ozone_file,      &
       prescribed_ozone_filelist,  &
       prescribed_ozone_datapath,  &
       prescribed_ozone_type,      &
       prescribed_ozone_rmfile,    &
       prescribed_ozone_cycle_yr,  &
       prescribed_ozone_fixed_ymd, &
       prescribed_ozone_fixed_tod      
    !-----------------------------------------------------------------------------
 
    ! Initialize namelist variables from local module variables.
    prescribed_ozone_name     = fld_name
    prescribed_ozone_file     = filename
    prescribed_ozone_filelist = filelist
    prescribed_ozone_datapath = datapath
    prescribed_ozone_type     = data_type
    prescribed_ozone_rmfile   = rmv_file
    prescribed_ozone_cycle_yr = cycle_yr
    prescribed_ozone_fixed_ymd= fixed_ymd
    prescribed_ozone_fixed_tod= fixed_tod

    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, prescribed_ozone_nl, iostat=ierr)
    if (ierr /= 0) call endrun(' error reading prescribed_ozone namelist ')
    close(unitn)

    ! Update module variables with user settings.
    fld_name   = prescribed_ozone_name
    filename   = prescribed_ozone_file
    filelist   = prescribed_ozone_filelist
    datapath   = prescribed_ozone_datapath
    data_type  = prescribed_ozone_type
    rmv_file   = prescribed_ozone_rmfile
    cycle_yr   = prescribed_ozone_cycle_yr
    fixed_ymd  = prescribed_ozone_fixed_ymd
    fixed_tod  = prescribed_ozone_fixed_tod

    ! Turn on prescribed volcanics if user has specified an input dataset.
    if (len_trim(filename) > 0 ) has_prescribed_ozone = .true.

    end subroutine read_nml_prescribed_ozone

    subroutine prescribed_ozone_init(ncol, dtime, lat, lon)
    
    use grist_tracer_data, only : trcdata_init
! io
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: dtime
    real(r8), intent(in) :: lat(:)             ! longitude in radian
    real(r8), intent(in) :: lon(:)             ! latitude in radian 
! local
    character(len=32) :: specifier(1)
    integer :: i
 
    if( has_prescribed_ozone ) then
        if(mpi_rank()==0)print*,'ozone is prescribed in :'//trim(filename)
    else
        return
    end if

    !do i = 1, pstate_cam%total_ghg_num
    !    if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'O3')then
    !        o3_idx = pstate_cam%ghg_at_pc_full_level(i)%idx
    !    end if
    !end do

    specifier(1) = trim(ozone_name)//':'//trim(fld_name)

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.       !Only aerosol use in_pbuf, LiXH
    call trcdata_init(  dtime, ncol, lon, lat,                                 & 
                        specifier, filename, filelist, datapath, fields, file, &
                        rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)
    

    end subroutine prescribed_ozone_init

    subroutine prescribed_ozone_adv(ncol, nstep, dtime, lon, lat)

    use string_utils,       only : to_lower, GLC
    use grist_tracer_data,  only : advance_trcdata

! io
    integer(i4), intent(in) :: ncol
    integer(i4), intent(in) :: nstep
    real(r8)   , intent(in) :: dtime
    real(r8)   , intent(in) :: lat(:)             ! longitude in radian
    real(r8)   , intent(in) :: lon(:)             ! latitude in radian 

! local
    integer(i4)             :: i, k, count
    real(r8)                :: zzdh,zzd1,zzd
    real(r8), allocatable   :: ozone_ver(:,:)
    real(r8), allocatable   :: log_lev_data(:), log_lev_model(:)

   ! real(r8) :: to_mmr(nlev,ncol)
    real(r8) :: to_mmr(ncol,nlev,1)
    real(r8) :: molmass
    real(r8) :: amass
    character(len=32) :: units_str

    if( .not. has_prescribed_ozone ) return

        molmass = 47.9981995_r8
        amass   = mwdry

        call advance_trcdata(nstep, dtime, ncol, lon, lat, &
                             fields, file)

        units_str = trim(to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))))

        select case ( units_str )
        case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
           ! to_mmr(:,:ncol) = (molmass*1.e6_r8*boltz*pstate%temp_at_pc_full_level%f(:,:ncol))/(amass*pstate%pressure_at_pc_full_level%f(:,:ncol))
           to_mmr(:ncol,:nlev,1) = (molmass*1.e6_r8*boltz*pstate_wrf%t_phy(:ncol,:nlev,1))/(amass*pstate_wrf%p_phy(:ncol,:nlev,1))
        case ('kg/kg','mmr')
           ! to_mmr(:,:ncol) = 1._r8
           to_mmr(:ncol,:nlev,1) = 1._r8
        case ('mol/mol','mole/mole','vmr','fraction')
           ! to_mmr(:,:ncol) = molmass/amass
           to_mmr(:ncol,:nlev,1) = molmass/amass
        case default
           if(mpi_rank()==0)print*, 'prescribed_ozone_adv: units = ',trim(fields(1)%units) ,' are not recognized'
           call endrun('prescribed_ozone_adv: units are not recognized')
        end select

        pstate_wrf%o3(1:ncol,1:nlev,1) = transpose(fields(1)%data(nlev:1:-1,1:ncol))*to_mmr(1:ncol,1:nlev,1)

    end subroutine prescribed_ozone_adv

 end module grist_prescribed_ozone 
