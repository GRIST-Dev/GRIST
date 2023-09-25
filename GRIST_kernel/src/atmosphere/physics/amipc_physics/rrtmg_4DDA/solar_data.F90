!==========================================================================
!
!  Created by LiXiaohan on 19/8/19, adopted from CAM5 (chemistry).
!
!  Solar irradiance / photon flux data
!===========================================================================
 module solar_data
    use grist_constants,          only: r8
    use grist_handle_error,       only: endrun
    use grist_wrap_nf
    use grist_mpi
    use grist_nml_module,         only: start_ymd, start_tod
 
  implicit none
  save
  private

  public :: read_nml_solar_data,   &
            solar_data_init,       &
            end_of_solar_data,     &
            solar_data_advance

  public :: nbins, we  ! number of wavelength samples of spectrum, wavelength endpoints
  public :: sol_etf
  public :: sol_irrad
  public :: sol_tsi
  public :: do_spctrl_scaling
  public :: has_spectrum
  public :: has_ref_spectrum
  public :: ssi_ref
  public :: ref_tsi

  integer :: nbins
  integer :: ntimes
  real(r8), allocatable :: sol_etf(:)
  real(r8), allocatable :: irradi(:,:)
  real(r8)              :: itsi(2)
  real(r8), allocatable :: ssi_ref(:)  ! a reference spectrum constructed from 3 solar cycles of data

  real(r8), allocatable :: we(:)
  real(r8), allocatable :: data_times(:)

  integer :: last_index = 1
  integer :: file_id
  integer :: ssi_vid
  integer :: tsi_vid
  integer :: ref_vid
  integer :: tsi_ref_vid

  logical :: initialized = .false.
  logical :: has_spectrum = .false.
  logical :: has_ref_spectrum = .false.
  logical :: has_tsi = .false.

  real(r8), allocatable :: sol_irrad(:)
  real(r8)              :: sol_tsi = -1.0_r8
  real(r8)              :: ref_tsi

  real(r8), allocatable :: irrad_fac(:)
  real(r8), allocatable :: etf_fac(:)
  real(r8), allocatable :: dellam(:)

  logical  :: fixed_scon
  logical  :: fixed_solar
  real(r8) :: offset_time

  logical            :: do_spctrl_scaling = .false.

! namelist vars
  character(len=256) :: solar_data_file = ''
  character(len=8)   :: solar_data_type = 'SERIAL'      ! "FIXED" or "SERIAL"
  integer            :: solar_data_ymd = 0              ! YYYYMMDD for "FIXED" type
  integer            :: solar_data_tod = 0              ! seconds of day for "FIXED" type
  real(r8)           :: solar_const = -9999._r8         ! constant TSI (W/m2)
  logical            :: solar_htng_spctrl_scl = .false. ! do rad heating spectral scaling

 contains

    subroutine read_nml_solar_data( nlfile )
    ! io
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    ! local
    integer :: unitn, ierr

    namelist /solar_inparm/ &
         solar_data_file, solar_data_type, solar_data_ymd, solar_data_tod, solar_const, &
         solar_htng_spctrl_scl
 
    unitn = 111
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, solar_inparm, iostat=ierr)
    if (ierr /= 0) call endrun('read_nml_solar_data: ERROR reading namelist')
    close(unitn)

    if ( (len_trim(solar_data_file) > 0) .and. (solar_const>0._r8) ) then
       call endrun('read_nml_solar_data: ERROR cannot specify both solar_data_file and solar_const')      
    endif

    if ( len_trim(solar_data_file) > 0 ) then
       fixed_scon = .false.
    else
       fixed_scon = .true.
    endif
    
    if ( solar_const>0._r8 ) then
       sol_tsi = solar_const
    endif
    
    !---LiXH did not initialize ref_tsi here--->
    !ref_tsi = nan
    !<--LiXH did not initialize ref_tsi here----
    
    if (mpi_rank()==0) then
       print*, 'read_nml_solar_data: solar_const (W/m2) = ', solar_const
       if ( .not.fixed_scon ) then
          print*, 'read_nml_solar_data: solar_data_file = ',trim(solar_data_file)
          print*, 'read_nml_solar_data: solar_data_type = ',trim(solar_data_type)
          print*, 'read_nml_solar_data: solar_data_ymd  = ',solar_data_ymd
          print*, 'read_nml_solar_data: solar_data_tod  = ',solar_data_tod
       endif
    endif

    end subroutine read_nml_solar_data


    subroutine solar_data_init(nstep, dtime)

    use grist_constants, only : c0, planck
    use grist_physics_iofile,   only: getfile

    include 'netcdf.inc'
    ! io
    integer,  intent(in)     :: nstep
    real(r8), intent(in)     :: dtime
    ! local
    integer  :: astat, dimid, vid
    character(len=256) :: filen   
    real(r8), allocatable :: lambda(:)
    integer,  allocatable :: dates(:)
    integer,  allocatable :: datesecs(:)

    integer  :: i, wvl_vid
    real(r8), parameter :: c = c0     ! speed of light (m/s)
    real(r8), parameter :: h = planck ! Planck's constant (Joule sec)
    real(r8), parameter :: fac = 1._r8/(h*c)

    real(r8) :: model_time, time
    real(r8) :: tmp(1)
    integer  :: ierr

    has_spectrum = .false.
    if ( fixed_scon ) return

    fixed_solar = trim(solar_data_type) == 'FIXED'

    call getfile( solar_data_file, filen, 0 )
    call wrap_open( trim(filen), 0, file_id )
    if(mpi_rank()==0) print*,'solar_data_init: data file = ',trim(filen)
    ierr = nf_inq_varid (file_id, 'ssi', ssi_vid)
    has_spectrum = ierr==NF_NOERR

    ierr = nf_inq_varid (file_id, 'tsi', tsi_vid)
    has_tsi      = ierr==NF_NOERR .and. solar_const<0._r8

    ierr = nf_inq_varid (file_id, 'ssi_ref', ref_vid)
    has_ref_spectrum = ierr==NF_NOERR

    call wrap_inq_dimid(file_id, 'time', dimid)
    call wrap_inq_dimlen(file_id, dimid, ntimes)

    if ( has_spectrum ) then
        ierr = nf_inq_varid (file_id, 'wavelength', wvl_vid)
        if ( ierr==NF_NOERR ) then
            call wrap_inq_dimid(file_id,'wavelength', dimid )
        else ! for backwards compatibility
            call wrap_inq_varid(file_id, 'wvl', wvl_vid)
            call wrap_inq_dimid(file_id, 'wvl', dimid )
        end if
        call wrap_inq_dimlen(file_id, dimid, nbins)
        if ( has_ref_spectrum ) then
            call wrap_inq_varid(file_id, 'tsi_ref', tsi_ref_vid)
        end if
    end if

    do_spctrl_scaling = has_spectrum .and. solar_htng_spctrl_scl

    allocate(lambda(nbins))
    allocate(dellam(nbins))
    allocate(data_times(ntimes))
    allocate(dates(ntimes))
    allocate(datesecs(ntimes))
    allocate(irrad_fac(nbins))
    allocate(etf_fac(nbins))
    allocate(sol_irrad(nbins))
    allocate(ssi_ref(nbins))

    call wrap_inq_varid(file_id, 'date', vid)
    call wrap_get_var_int(file_id, vid, dates)

    ierr = nf_inq_varid (file_id, 'datesec', vid)
    if ( ierr==NF_NOERR ) then
        call wrap_get_var_int(file_id, vid, datesecs)
    else
       datesecs(:) = 0
    endif
       
    if (has_spectrum) then
        call wrap_get_var_realx(file_id, wvl_vid, lambda)
        call wrap_inq_varid(file_id, 'band_width', vid)
        call wrap_get_var_realx(file_id, vid, dellam)
    endif

    if(mpi_rank()==0) print*, 'solar_data_init: has_ref_spectrum',has_ref_spectrum

    if ( has_ref_spectrum ) then
        call wrap_inq_varid(file_id, 'ssi_ref', vid)
        call wrap_get_var_realx(file_id, vid, ssi_ref)
        call wrap_get_var_realx(file_id, tsi_ref_vid, tmp)
        ref_tsi = tmp(1)
    endif

    if ( has_spectrum ) then
       allocate(sol_etf(nbins))
       allocate(irradi(nbins,2))
       allocate(we(nbins+1))

       we(:nbins)  = lambda(:nbins) - 0.5_r8*dellam(:nbins)
       we(nbins+1) = lambda(nbins)  + 0.5_r8*dellam(nbins)
       do i = 1,nbins
          irrad_fac(i) = 1.e-3_r8                ! mW/m2/nm --> W/m2/nm
          etf_fac(i)   = 1.e-16_r8*lambda(i)*fac ! mW/m2/nm --> photons/cm2/sec/nm
       enddo
       if(has_ref_spectrum) then
          ssi_ref = ssi_ref * 1.e-3_r8        ! mW/m2/nm --> W/m2/nm
       endif
    endif

    offset_time = 0._r8

    if ( solar_data_ymd > 0 ) then
      call get_model_time( nstep, dtime, model_time )
      call convert_date( solar_data_ymd, solar_data_tod, time )
      offset_time = time - model_time
    endif

    call convert_dates( dates, datesecs, data_times )

    data_times = data_times - offset_time

    deallocate(lambda)
    deallocate(dates)
    deallocate(datesecs)

    ! need to force data loading when the model starts at a time =/ 00:00:00.000
    ! -- may occur in restarts also
    call solar_data_advance(nstep, dtime)
    initialized = .true.

    end subroutine solar_data_init


    subroutine end_of_solar_data()
    deallocate(dellam)
    deallocate(data_times)
    deallocate(irrad_fac)
    deallocate(etf_fac)
    deallocate(sol_irrad)
    deallocate(ssi_ref)
    if(allocated(sol_etf)) deallocate(sol_etf)
    if(allocated(irradi))  deallocate(irradi)
    if(allocated(we))      deallocate(we) 

    end subroutine end_of_solar_data


! Purpose : Reads in the ETF data for the current date.  
    subroutine solar_data_advance(nstep, dtime)

    use grist_constants,   only : cday
    ! io
    integer,  intent(in)     :: nstep
    real(r8), intent(in)     :: dtime
    ! local
    integer  :: year, month, day, sec
    integer  :: index, i, nt
    integer  :: offset(2), count(2)
    logical  :: do_adv, read_data
    real(r8) :: time, delt
    real(r8) :: data(nbins)
    integer  :: ierr

    if ( fixed_scon ) return
    if ( fixed_solar .and. initialized ) return

    index = -1
    call get_model_time(nstep, dtime, time, year=year, month=month, day=day, seconds=sec )
    read_data = time > data_times(last_index) .or. .not.initialized

    if ( read_data ) then

       find_ndx: do i = last_index, ntimes
          if ( data_times(i) > time ) then
             index = i-1
             exit find_ndx
          endif
       enddo find_ndx

       last_index = index+1

       nt = 2

       if ( index < 1 ) then
          write(*,102) year,month,day,sec
          call endrun('solar_data_advance: failed to read data from '//trim(solar_data_file))
       endif

       ! get the surrounding time slices
       offset = (/ 1, index /)
       count =  (/ nbins, nt /)

   
       if (has_spectrum) then
          call wrap_get_vara_realx(file_id, ssi_vid, offset, count, irradi) 
       endif
       if (has_tsi .and. (.not.do_spctrl_scaling)) then
          call wrap_get_vara_realx(file_id, tsi_vid, (/index/), (/nt/), itsi)
          if ( any(itsi(:nt) < 0._r8) ) then
             call endrun( 'solar_data_advance: invalid or missing tsi data  ' )
          endif
       endif
    else
       index = last_index - 1
    endif

    delt = ( time - data_times(index) ) / ( data_times(index+1) - data_times(index) )

    ! this assures that FIXED data are b4b on restarts
    if ( fixed_solar ) then
       delt = dble(int(delt*cday+.5_r8))/dble(cday)
    endif

    if (has_spectrum) then
       data(:) = irradi(:,1) + delt*( irradi(:,2) - irradi(:,1) )
       do i = 1,nbins
          sol_irrad(i) = data(i)*irrad_fac(i) ! W/m2/nm
          sol_etf(i)   = data(i)*etf_fac(i)   ! photons/cm2/sec/nm 
       enddo
    endif
    if (has_tsi .and. (.not.do_spctrl_scaling)) then
       sol_tsi = itsi(1) + delt*( itsi(2) - itsi(1) )
       if ( mpi_rank()==0 ) then
          if (day == 1 .and. sec == 0) then
             write(*,101) year, month, day, sec, sol_tsi
          end if
       endif
    endif

 101 FORMAT('solar_data_advance: date, TSI : ',i4.4,'-',i2.2,'-',i2.2,'-',i5.5,',  ',f12.6)
 102 FORMAT('solar_data_advance: not able to find data for : ',i4.4,'-',i2.2,'-',i2.2,'-',i5.5)

  end subroutine solar_data_advance
 

  subroutine convert_dates( dates, secs, times )

    use grist_time_manager, only: set_time_float_from_date

    integer,  intent(in)  :: dates(:)
    integer,  intent(in)  :: secs(:)

    real(r8), intent(out) :: times(:)

    integer :: year, month, day, sec,n ,i

    n = size( dates ) 

    do i=1,n
       year = dates(i)/10000
       month = (dates(i)-year*10000)/100
       day = dates(i)-year*10000-month*100
       sec = secs(i)
       call set_time_float_from_date(start_ymd, start_tod, times(i), year, month, day, sec )
    enddo

  end subroutine convert_dates


  subroutine convert_date( date, sec, time )

    integer,  intent(in)  :: date
    integer,  intent(in)  :: sec
    real(r8), intent(out) :: time

    integer :: dates(1), secs(1)
    real(r8) :: times(1)
    dates(1) = date
    secs(1) = sec
    call convert_dates( dates, secs, times )

    time = times(1)
  end subroutine convert_date


  subroutine get_model_time( nstep, dtime, time, year, month, day, seconds )

    use grist_time_manager, only: get_curr_date
    !io
    integer,  intent(in)  :: nstep
    real(r8), intent(in)  :: dtime
    real(r8), intent(out) :: time
    integer, optional, intent(out) :: year, month, day, seconds
    ! local
    integer  :: yr, mn, dy, sc, date

    call get_curr_date(start_ymd, start_tod, nstep, dtime, yr, mn, dy, sc)

    date = yr*10000 + mn*100 + dy
    call convert_date( date, sc, time )

    if (present(year))    year = yr
    if (present(month))   month = mn
    if (present(day))     day = dy
    if (present(seconds)) seconds = sc

  end subroutine get_model_time

 end module solar_data 
