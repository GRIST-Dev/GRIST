!===========================================================================================
!
! Created by LiXiaohan on 20/12/15, adopted from CAM5
!
! module used to read (and interpolate) offline tracer data (sources and
! mixing ratios)
!
!===========================================================================================

 module grist_tracer_data
   use grist_constants,                    only: i4, r8, deg2rad, pi
   use grist_physics_data_structure,       only: pstate
   use grist_PhysC_data_structure,         only: pstate_cam
   use grist_nml_module,                   only: nlev, nlevp, start_ymd, start_tod
   use grist_time_manager,                 only: get_curr_date, &
                                                 set_time_float_from_date
   use grist_handle_error,                 only: endrun
   use grist_wrap_nf
   use grist_mpi


   implicit none
   include 'netcdf.inc'

   private  ! all unless made public
   save
 
   public :: trfld, trfile
   public :: trcdata_init,   &
             advance_trcdata


   logical :: log_print = .true.

! not used in GRIST, LiXH
!   type input3d
!      real(r8), dimension(:,:,:), allocatable :: data
!   endtype input3d     
 
   type input2d
      real(r8), dimension(:,:), allocatable :: data
   endtype input2d

   type input1d
      real(r8), dimension(:), allocatable :: data
   endtype input1d


   type trfld
      real(r8), dimension(:,:), allocatable :: data
      type(input2d), dimension(4) :: input
      character(len=32) :: srcnam
      character(len=32) :: fldnam
      character(len=32) :: units
      integer :: var_id
      integer :: coords(4) ! LATDIM | LONDIM | LEVDIM | TIMDIM
      integer :: order(4) ! LATDIM | LONDIM | LEVDIM | TIMDIM
      logical :: srf_fld = .false.
      integer :: pbuf_ndx = -1
   endtype trfld

   type trfile
      type(input1d), dimension(4) :: ps_in
      character(len=256) :: pathname = ' '
      character(len=256) :: curr_filename = ' '
      character(len=256) :: next_filename = ' '
      integer            :: curr_fileid
      integer            :: next_fileid
 
!      type(var_desc_t) :: currfnameid  ! pio restart file var id 
!      type(var_desc_t) :: nextfnameid  ! pio restart file var id
 
      character(len=256) :: filenames_list = ''
      real(r8) :: datatimem = -1.e36_r8     ! time of prv. values read in
      real(r8) :: datatimep = -1.e36_r8     ! time of nxt. values read in
      real(r8) :: datatimes(4)
      integer :: interp_recs
      real(r8), dimension(:), allocatable :: curr_data_times
      real(r8), dimension(:), allocatable :: next_data_times
      logical :: remove_trc_file = .false.  ! delete file when finished with it
      real(r8) :: offset_time
      integer :: cyc_ndx_beg
      integer :: cyc_ndx_end
      integer :: cyc_yr = 0
      real(r8) :: one_yr = 0
      real(r8) :: curr_mod_time ! model time - calendar day
      real(r8) :: next_mod_time ! model time - calendar day - next time step
      integer :: nlon
      integer :: nlat
      integer :: nlev
      integer :: nilev
      integer :: ps_coords(3) ! LATDIM | LONDIM | TIMDIM
      integer :: ps_order(3) ! LATDIM | LONDIM | TIMDIM
      real(r8), dimension(:),   allocatable :: lons
      real(r8), dimension(:),   allocatable :: lats
      real(r8), dimension(:),   allocatable :: levs
      real(r8), dimension(:),   allocatable :: ilevs
      real(r8), dimension(:),   allocatable :: hyam
      real(r8), dimension(:),   allocatable :: hybm
      real(r8), dimension(:,:), allocatable :: ps
      real(r8), dimension(:),   allocatable :: hyai
      real(r8), dimension(:),   allocatable :: hybi
      real(r8), dimension(:,:), allocatable :: weight_x, weight_y
      integer,  dimension(:,:), allocatable :: index_x, index_y
      real(r8)                        :: p0
      integer :: ps_id
      logical,  dimension(:), allocatable :: in_pbuf
      logical :: has_ps = .false.
      logical :: zonal_ave = .false.
      logical :: alt_data = .false.
      logical :: cyclical = .false.
      logical :: cyclical_list = .false.
      logical :: weight_by_lat = .false.
      logical :: conserve_column = .false.
      logical :: fill_in_months = .false.
      logical :: fixed = .false.
      logical :: initialized = .false.
      logical :: top_bndry = .false.
      logical :: stepTime = .false.  ! Do not interpolate in time, but use stepwise times
   endtype trfile
 
   integer, public, parameter :: MAXTRCRS = 100
 
   integer, parameter :: LONDIM = 1
   integer, parameter :: LATDIM = 2
   integer, parameter :: LEVDIM = 3
   integer, parameter :: TIMDIM = 4
 
   integer, parameter :: PS_TIMDIM = 3
 
   integer, parameter :: ZA_LATDIM = 1
   integer, parameter :: ZA_LEVDIM = 2
   integer, parameter :: ZA_TIMDIM = 3

   integer, parameter :: nm=1    ! array index for previous (minus) data
   integer, parameter :: np=2    ! array index for next (plus) data
 
   contains

    subroutine trcdata_init( dtime, ncol, grid_lon, grid_lat,                     &
                             specifier, filename, filelist, datapath, flds, file, &
                             rmv_file, data_cycle_yr, data_fixed_ymd, data_fixed_tod, data_type )

    use grist_horizontal_interpolate, only : xy_interp_init

    implicit none

    real(r8),            intent(in)    :: dtime
    integer,             intent(in)    :: ncol
    real(r8),            intent(in)    :: grid_lat(:)             ! longitude in radian
    real(r8),            intent(in)    :: grid_lon(:)             ! latitude in radian  
    character(len=*),    intent(in)    :: specifier(:)
    character(len=*),    intent(in)    :: filename
    character(len=*),    intent(in)    :: filelist
    character(len=*),    intent(in)    :: datapath
    type(trfld), allocatable, intent(inout) :: flds(:)
    type(trfile),        intent(inout) :: file
    logical,             intent(in)    :: rmv_file
    integer,             intent(in)    :: data_cycle_yr
    integer,             intent(in)    :: data_fixed_ymd
    integer,             intent(in)    :: data_fixed_tod
    character(len=*),    intent(in)    :: data_type

    integer :: f, mxnflds, astat
    integer :: str_yr, str_mon, str_day
    integer :: lon_dimid, lat_dimid, lev_dimid, tim_dimid
    integer :: dimids(4), did
    integer :: idx, varid
    integer :: ierr
    integer :: errcode
    real(r8) :: start_time, time1, time2
    integer :: i1,i2,j1,j2,i
    integer :: nvardims, vardimids(4)
    character(len=80) :: data_units
    real(r8):: tmp(1)

    logical :: find_aero = .false.

    call specify_fields( specifier, flds )

    file%datatimep=-1.e36_r8
    file%datatimem=-1.e36_r8

    mxnflds = 0 
    if (allocated(flds)) mxnflds = size( flds )

    if (mxnflds < 1) return
    
    file%remove_trc_file = rmv_file
    file%pathname = trim(datapath)
    file%filenames_list = trim(filelist)

    file%fill_in_months = .false.
    file%cyclical = .false.  
    file%cyclical_list = .false.  

! does not work when compiled with pathf90
!    select case ( to_upper(data_type) )
    select case ( data_type )
    case( 'FIXED' )
       file%fixed = .true.
    case( 'INTERP_MISSING_MONTHS' )
       file%fill_in_months = .true.
    case( 'CYCLICAL' )
       file%cyclical = .true.
       file%cyc_yr = data_cycle_yr
    case( 'CYCLICAL_LIST' )
       file%cyclical_list = .true.
       file%cyc_yr = data_cycle_yr
    case( 'SERIAL' )
    case default 
       if(mpi_rank()==0)then 
       print*, 'trcdata_init: invalid data type: '//trim(data_type)//' file: '//trim(filename)
       print*, 'trcdata_init: valid data types: SERIAL | CYCLICAL | CYCLICAL_LIST | FIXED | INTERP_MISSING_MONTHS '
       end if
       call endrun('trcdata_init: invalid data type: '//trim(data_type)//' file: '//trim(filename))
    endselect

    if ( (.not.file%fixed) .and. ((data_fixed_ymd>0._r8) .or.(data_fixed_tod>0._r8))) then
       call endrun('trcdata_init: Cannot specify data_fixed_ymd or data_fixed_tod if data type is not FIXED')
    endif
    if ( (.not.file%cyclical) .and. (data_cycle_yr>0._r8) ) then
       call endrun('trcdata_init: Cannot specify data_cycle_yr if data type is not CYCLICAL')
    endif

    if (mpi_rank()==0) then
       print*, 'trcdata_init: data type: '//trim(data_type)//' file: '//trim(filename)
    endif

    ! if there is no list of files (len_trim(file%filenames_list)<1) then
    !  -> set curr_filename from namelist rather from restart data
    if ( len_trim(file%curr_filename)<1 .or. len_trim(file%filenames_list)<1 .or. file%fixed ) then ! initial run
       file%curr_filename = trim(filename)

       call get_model_time(0, dtime, file)

       if ( file%fixed ) then
          str_yr = data_fixed_ymd/10000
          str_mon = (data_fixed_ymd - str_yr*10000)/100
          str_day = data_fixed_ymd - str_yr*10000 - str_mon*100
          call set_time_float_from_date(start_ymd, start_tod, start_time, str_yr, str_mon, str_day, data_fixed_tod )
          file%offset_time = start_time - file%curr_mod_time
       else
          file%offset_time = 0
       endif
    endif

    call set_time_float_from_date(start_ymd, start_tod, time2, 2, 1, 1, 0 )
    call set_time_float_from_date(start_ymd, start_tod, time1, 1, 1, 1, 0 )

    file%one_yr = time2-time1

    if ( file%cyclical .or. file%cyclical_list) then
       file%cyc_ndx_beg = -1
       file%cyc_ndx_end = -1
       if ( file%cyc_yr /= 0 ) then
          call set_time_float_from_date(start_ymd, start_tod, time1, file%cyc_yr  , 1, 1, 0 )
          call set_time_float_from_date(start_ymd, start_tod, time2, file%cyc_yr+1, 1, 1, 0 )
          file%one_yr = time2-time1
      endif

       call open_trc_datafile( file%curr_filename, file%pathname, file%curr_fileid, file%curr_data_times, &
            cyc_ndx_beg=file%cyc_ndx_beg, cyc_ndx_end=file%cyc_ndx_end, cyc_yr=file%cyc_yr )
    else
       call open_trc_datafile( file%curr_filename, file%pathname, file%curr_fileid, file%curr_data_times )
       file%curr_data_times = file%curr_data_times - file%offset_time
    endif

    ierr = nf_inq_dimid(file%curr_fileid, 'lon', idx)

    file%zonal_ave = (ierr/=NF_NOERR)

    if ( file%zonal_ave ) then
       file%nlon = 1
    else
       call get_dimension( file%curr_fileid, 'lon', file%nlon, dimid=lon_dimid, data=file%lons )

       file%lons =  file%lons * deg2rad
    endif

    call wrap_inq_dimid(file%curr_fileid, 'time', tim_dimid)

    call get_dimension( file%curr_fileid, 'lat', file%nlat, dimid=lat_dimid, data=file%lats )

    file%lats =  file%lats * deg2rad


   ! LiXH 
   ! allocate( file%ps(file%nlon,file%nlat) )
    ierr = nf_inq_varid( file%curr_fileid, 'PS', file%ps_id )
    file%has_ps = (ierr==NF_NOERR)

    ierr = nf_inq_dimid( file%curr_fileid, 'altitude', idx )
    file%alt_data = (ierr==NF_NOERR)

    if ( file%has_ps) then
       if ( file%zonal_ave ) then
          call wrap_inq_vardimid (file%curr_fileid, file%ps_id, dimids(1:2))
          do did = 1,2
             if      ( dimids(did) == lat_dimid ) then
                file%ps_coords(LATDIM) = did
                file%ps_order(did) = LATDIM
             else if ( dimids(did) == tim_dimid ) then
                file%ps_coords(PS_TIMDIM) = did
                file%ps_order(did) = PS_TIMDIM
             endif
          enddo
       else
           call wrap_inq_vardimid (file%curr_fileid, file%ps_id, dimids(1:3))
          do did = 1,3
             if      ( dimids(did) == lon_dimid ) then
                file%ps_coords(LONDIM) = did
                file%ps_order(did) = LONDIM
             else if ( dimids(did) == lat_dimid ) then
                file%ps_coords(LATDIM) = did
                file%ps_order(did) = LATDIM
             else if ( dimids(did) == tim_dimid ) then
                file%ps_coords(PS_TIMDIM) = did
                file%ps_order(did) = PS_TIMDIM
             endif
          enddo
       end if
    endif

    if (log_print .and. mpi_rank()==0)print*, 'trcdata_init: file%has_ps = ' , file%has_ps 

    if (file%alt_data) then
       call get_dimension( file%curr_fileid, 'altitude_int', file%nilev, data=file%ilevs  )

       call get_dimension( file%curr_fileid, 'altitude',     file%nlev, dimid=lev_dimid,  data=file%levs  )
    else
       call get_dimension( file%curr_fileid, 'lev',          file%nlev, dimid=lev_dimid,  data=file%levs  )

       if (lev_dimid>0) then
          file%levs =  file%levs*100._r8 ! mbar->pascals
       endif
    endif

    if (file%has_ps) then

       allocate( file%hyam(file%nlev),  file%hybm(file%nlev) )
       allocate( file%hyai(file%nlev+1),  file%hybi(file%nlev+1) )

       ierr = nf_inq_varid( file%curr_fileid, 'P0', varid)
       if ( ierr == NF_NOERR ) then
          call wrap_get_var_realx( file%curr_fileid, varid, tmp )
          file%p0 = tmp(1)
       else
          file%p0 = 100000._r8
       endif
       call wrap_inq_varid( file%curr_fileid, 'hyam', varid )
       call wrap_get_var_realx( file%curr_fileid, varid, file%hyam )
       call wrap_inq_varid( file%curr_fileid, 'hybm', varid )
       call wrap_get_var_realx( file%curr_fileid, varid, file%hybm )
       if (file%conserve_column) then
          call wrap_inq_varid( file%curr_fileid, 'hyai', varid )
          call wrap_get_var_realx( file%curr_fileid, varid, file%hyai )
          call wrap_inq_varid( file%curr_fileid, 'hybi', varid )
          call wrap_get_var_realx( file%curr_fileid, varid, file%hybi )
       endif

       !LiXH
       !allocate( file%ps(ncol) )
       allocate( file%ps_in(1)%data(ncol) )
       allocate( file%ps_in(2)%data(ncol) )
       if( file%fill_in_months ) then
          allocate( file%ps_in(3)%data(ncol) )
          allocate( file%ps_in(4)%data(ncol) )
       end if
    endif

    flds_loop: do f = 1,mxnflds

       ! get netcdf variable id for the field
       call wrap_inq_varid( file%curr_fileid, flds(f)%srcnam, flds(f)%var_id )

       ! determine if the field has a vertical dimension
       if (lev_dimid>0) then
          call wrap_inq_varndims(  file%curr_fileid, flds(f)%var_id, nvardims )
          call wrap_inq_vardimid(  file%curr_fileid, flds(f)%var_id, vardimids(:nvardims) )
          flds(f)%srf_fld = .not.any(vardimids(:nvardims)==lev_dimid)
       else
          flds(f)%srf_fld = .true.
       endif

       ! allocate memory only if not already in pbuf2d
       ! LiXH: only aerosol is set in pbuf2d
       if ( .not. file%in_pbuf(f) ) then 
       !LiXH: unify flds(f)%data, data_out, and pstate_cam%aerosol_at_pc_full_level 
       !   if ( flds(f)%srf_fld .or. file%top_bndry ) then
       !      allocate( flds(f)%data(1,ncol) )
       !   else
             allocate( flds(f)%data(nlev,ncol) ); flds(f)%data=0._r8
       !   endif
       else
           find_aero = .false.
           do i = 1, pstate_cam%total_aerosol_num
             if(trim(flds(f)%fldnam) .eq. trim(pstate_cam%aerosol_at_pc_full_level(i)%name))then
                 flds(f)%pbuf_ndx = pstate_cam%aerosol_at_pc_full_level(i)%idx
                 find_aero = .true.
                 exit
             end if
           end do
           if(.not.find_aero)then
               if(mpi_rank()==0)print*,'Error in grist_aero_init:cannot find aerosol'
               call endrun('In trcdata_init')
           end if 
       endif
 
       if (flds(f)%srf_fld) then
          allocate( flds(f)%input(1)%data(1, ncol) )
       else
          allocate( flds(f)%input(1)%data(file%nlev, ncol) )
       endif
       if (flds(f)%srf_fld) then
          allocate( flds(f)%input(2)%data(1, ncol) )
       else
          allocate( flds(f)%input(2)%data(file%nlev, ncol) )
       endif

       if( file%fill_in_months ) then
          if (flds(f)%srf_fld) then
             allocate( flds(f)%input(3)%data(1, ncol) )
          else
             allocate( flds(f)%input(3)%data(file%nlev, ncol) )
          endif
          if (flds(f)%srf_fld) then
             allocate( flds(f)%input(4)%data(1, ncol) )
          else
             allocate( flds(f)%input(4)%data(file%nlev, ncol) )
          endif
       endif

       if ( file%zonal_ave ) then
          call wrap_inq_vardimid(  file%curr_fileid, flds(f)%var_id, dimids(1:3) )
          do did = 1,3
             if      ( dimids(did) == lat_dimid ) then
                flds(f)%coords(ZA_LATDIM) = did
                flds(f)%order(did) = ZA_LATDIM
             else if ( dimids(did) == lev_dimid ) then
                flds(f)%coords(ZA_LEVDIM) = did
                flds(f)%order(did) = ZA_LEVDIM
             else if ( dimids(did) == tim_dimid ) then
                flds(f)%coords(ZA_TIMDIM) = did
                flds(f)%order(did) = ZA_TIMDIM
             endif
          enddo
       else if ( flds(f)%srf_fld ) then
          call wrap_inq_vardimid(  file%curr_fileid, flds(f)%var_id, dimids(1:3) )
          do did = 1,3
             if      ( dimids(did) == lon_dimid ) then
                flds(f)%coords(LONDIM) = did
                flds(f)%order(did) = LONDIM
             else if ( dimids(did) == lat_dimid ) then
                flds(f)%coords(LATDIM) = did
                flds(f)%order(did) = LATDIM
             else if ( dimids(did) == tim_dimid ) then
                flds(f)%coords(PS_TIMDIM) = did
                flds(f)%order(did) = PS_TIMDIM
             endif
          enddo
       else
          call wrap_inq_vardimid(  file%curr_fileid, flds(f)%var_id, dimids )
           do did = 1,4
             if      ( dimids(did) == lon_dimid ) then
                flds(f)%coords(LONDIM) = did
                flds(f)%order(did) = LONDIM
             else if ( dimids(did) == lat_dimid ) then
                flds(f)%coords(LATDIM) = did
                flds(f)%order(did) = LATDIM
             else if ( dimids(did) == lev_dimid ) then
                flds(f)%coords(LEVDIM) = did
                flds(f)%order(did) = LEVDIM
             else if ( dimids(did) == tim_dimid ) then
                flds(f)%coords(TIMDIM) = did
                flds(f)%order(did) = TIMDIM
             endif
          enddo
       endif
       call wrap_get_att_value (file%curr_fileid, flds(f)%var_id, 'units', data_units)
       data_units = trim(data_units)
       flds(f)%units = data_units(1:32)

    enddo flds_loop

! if weighting by latitude, compute weighting for horizontal interpolation
    if( file%weight_by_lat ) then
! weight_x & weight_y are weighting function for x & y interpolation
! LiXH use linear interpolation
        allocate(file%weight_x(ncol,2))
        allocate(file%weight_y(ncol,2))
        allocate(file%index_x(ncol,2))
        allocate(file%index_y(ncol,2))
        file%weight_x(:,:) = 0.0_r8
        file%weight_y(:,:) = 0.0_r8
        file%index_x(:,:) = 0
        file%index_y(:,:) = 0

        call xy_interp_init(ncol,grid_lon,grid_lat,                     &
                            file%nlon,file%nlat,file%lons,file%lats,    &
                            file%weight_x,file%weight_y,                &
                            file%index_x,file%index_y)
                            
    end if

    end subroutine trcdata_init


    subroutine advance_trcdata(nstep, dtime, ncol, grid_lon, grid_lat, &
                               flds, file )

    implicit none
    integer,             intent(in)    :: nstep
    integer,             intent(in)    :: ncol
    real(r8),            intent(in)    :: dtime
    real(r8),            intent(in)    :: grid_lon(:), grid_lat(:)
    type(trfile),        intent(inout) :: file
    type(trfld),         intent(inout) :: flds(:)

    ! local
    real(r8) :: data_time

    if ( .not.( file%fixed .and. file%initialized ) ) then

       call get_model_time(nstep, dtime, file)
       data_time = file%datatimep

       if ( file%cyclical .or. file%cyclical_list ) then
          ! wrap around
          if ( (file%datatimep<file%datatimem) .and. (file%curr_mod_time>file%datatimem) ) then
             data_time = data_time + file%one_yr 
          endif
       endif

    ! For stepTime need to advance if the times are equal
    ! Should not impact other runs?
       if ( file%curr_mod_time >= data_time ) then
          call read_next_trcdata( nstep, dtime , ncol, grid_lon, grid_lat, flds, file)
          if(mpi_rank()==0) print*, 'READ_NEXT_TRCDATA ', flds(:)%fldnam
       end if

    endif
    
    ! need to interpolate the data, regardless
    ! each mpi task needs to interpolate
    call interpolate_trcdata( ncol, flds, file )

    file%initialized = .true.

    end subroutine advance_trcdata
  
    
    subroutine get_model_time(nstep, dtime, file)

    ! io                                        
    integer,  intent(in)     :: nstep
    real(r8), intent(in)     :: dtime
    type(trfile), intent(inout) :: file

    !local
    integer  :: yr, mn, day, nsec

    call get_curr_date(start_ymd, start_tod, nstep, dtime, yr, mn, day, nsec)

    if ( file%cyclical .or. file%cyclical_list) yr = file%cyc_yr
    call set_time_float_from_date (start_ymd, start_tod, file%curr_mod_time, yr, mn, day, nsec)
    file%next_mod_time = file%curr_mod_time + dtime/86400._r8

    end subroutine get_model_time


    subroutine check_files( file, fids, itms, times_found)

    implicit none
    type(trfile),      intent(inout) :: file
    integer,           intent(out)   :: fids(2) ! ids of files that contains these recs
    integer, optional, intent(out)   :: itms(2)
    logical, optional, intent(inout) :: times_found

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    logical            :: list_cycled

    list_cycled = .false.

   !-----------------------------------------------------------------------
   !        If next time beyond the end of the time list,
   !        then increment the filename and move on to the next file
   !-----------------------------------------------------------------------
    if ((file%next_mod_time > file%curr_data_times(size(file%curr_data_times))).or.file%cyclical_list) then
       if (file%cyclical_list) then
          if ( allocated(file%next_data_times) ) then
              if ((file%curr_mod_time > file%datatimep)) then

               call advance_file(file)
     
            endif
         endif

       endif
       if ( .not. allocated(file%next_data_times) ) then
          ! open next file if not already opened...
          if (file%cyclical_list) then
              file%next_filename = incr_filename( file%curr_filename, filenames_list=file%filenames_list, datapath=file%pathname ,&
                                                  cyclical_list=file%cyclical_list, list_cycled=list_cycled)
          else
              file%next_filename = incr_filename( file%curr_filename, filenames_list=file%filenames_list, datapath=file%pathname)
          endif

          call open_trc_datafile( file%next_filename, file%pathname, file%next_fileid, file%next_data_times )
 
          file%next_data_times = file%next_data_times - file%offset_time

       endif
    endif
    
    !-----------------------------------------------------------------------
    !        If using next_data_times and the current is greater than or equal to the next, then
    !        close the current file, and set up for next file.
    !-----------------------------------------------------------------------
    if ( allocated(file%next_data_times) ) then
       if (file%cyclical_list .and. list_cycled) then    ! special case - list cycled 

          file%datatimem = file%curr_data_times(size(file%curr_data_times)) 
          itms(1)=size(file%curr_data_times)
          fids(1)=file%curr_fileid

          file%datatimep = file%next_data_times(1)
          itms(2)=1
          fids(2) = file%next_fileid

          times_found = .true.

       else if (file%curr_mod_time >= file%next_data_times(1)) then

          call advance_file(file)

       endif
    endif

    end subroutine check_files


    function incr_filename( filename, filenames_list, datapath, cyclical_list, list_cycled, abort )

    !-----------------------------------------------------------------------
    ! 	... Increment or decrement a date string withing a filename
    !           the filename date section is assumed to be of the form
    !           yyyy-dd-mm
    !-----------------------------------------------------------------------

    use string_utils, only : incstr

    implicit none

    character(len=256),           intent(in)    :: filename ! present dynamical dataset filename
    character(len=256), optional, intent(in)    :: filenames_list 
    character(len=256), optional, intent(in)    :: datapath
    logical           , optional, intent(in)    :: cyclical_list  ! If true, allow list to cycle
    logical           , optional, intent(out)   :: list_cycled
    logical           , optional, intent(in)    :: abort

    character(len=256)                :: incr_filename         ! next filename in the sequence

    ! set new next_filename ...

    !-----------------------------------------------------------------------
    !	... local variables
    !-----------------------------------------------------------------------
    integer :: pos, pos1, istat
    character(len=256) :: fn_new, line, filepath
    character(len=6)   :: seconds
    character(len=5)   :: num
    integer :: ios,unitnumber
    logical :: abort_run
 
    if (present(abort)) then 
       abort_run = abort
    else
       abort_run = .true.
    endif

    if (present(list_cycled)) list_cycled = .false.

    if (( .not. present(filenames_list)) .or.(len_trim(filenames_list) == 0)) then
       !-----------------------------------------------------------------------
       !	... ccm type filename
       !-----------------------------------------------------------------------
       pos = len_trim( filename )
       fn_new = filename(:pos)
       if ( mpi_rank()==0 ) print*, 'incr_flnm: old filename = ',trim(fn_new)
       if( fn_new(pos-2:) == '.nc' ) then
          pos = pos - 3
       end if
       istat = incstr( fn_new(:pos), 1 )
       if( istat /= 0 .and. mpi_rank()==0) then
          print*, 'incr_flnm: incstr returned ', istat
          print*, '           while trying to decrement ',trim( fn_new )
          call endrun('In incr_filename')
       end if

    else

       !-------------------------------------------------------------------
       !  ... open filenames_list
       !-------------------------------------------------------------------
       if ( mpi_rank()==0 ) print*, 'incr_flnm: old filename = ',trim(filename)
       if ( mpi_rank()==0 ) print*, 'incr_flnm: open filenames_list : ',trim(filenames_list)
       unitnumber = 123
       if ( present(datapath) ) then
         filepath = trim(datapath) //'/'// trim(filenames_list)
       else
         filepath = trim(filenames_list)
       endif

       open( unit=unitnumber, file=filepath, iostat=ios, status="OLD")
       if (ios /= 0) then
          call endrun('not able to open filenames_list file: '//trim(filepath))
       endif

       !-------------------------------------------------------------------
       !  ...  read file names
       !-------------------------------------------------------------------
       read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
       if (ios /= 0) then
          if (abort_run) then
             call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
          else 
             fn_new = 'NOT_FOUND'
             incr_filename = trim(fn_new)
             return
          endif
       endif

       !-------------------------------------------------------------------
       !      If current filename is '', then initialize with the first filename read in
       !      and skip this section.
       !-------------------------------------------------------------------
       if (filename /= '') then 

          !-------------------------------------------------------------------
          !       otherwise read until find current filename
          !-------------------------------------------------------------------
          do while( trim(line) /= trim(filename) )
             read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
             if (ios /= 0) then
                if (abort_run) then
                   call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
                else 
                   fn_new = 'NOT_FOUND'
                   incr_filename = trim(fn_new)
                   return
                endif
             endif
          enddo
   
          !-------------------------------------------------------------------
          !      Read next filename
          !-------------------------------------------------------------------
          read( unit=unitnumber, fmt='(A)', iostat=ios ) line 

          !---------------------------------------------------------------------------------
          !       If cyclical_list, then an end of file is not an error, but rather 
          !       a signal to rewind and start over
          !---------------------------------------------------------------------------------

          if (ios /= 0) then
             if (present(cyclical_list)) then
                if (cyclical_list) then
                   list_cycled=.true.
                   rewind(unitnumber)
                   read( unit=unitnumber, fmt='(A)', iostat=ios ) line 
                     ! Error here should never happen, but check just in case
                   if (ios /= 0) then
                      call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
                   endif
                else
                   call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
                endif
             else
                if (abort_run) then
                   call endrun('not able to increment file name from filenames_list file: '//trim(filenames_list))
                else 
                   fn_new = 'NOT_FOUND'
                   incr_filename = trim(fn_new)
                   return
                endif
             endif
          endif

       endif

       !---------------------------------------------------------------------------------
       !     Assign the current filename and close the filelist
       !---------------------------------------------------------------------------------
       fn_new = trim(line)

       close(unit=unitnumber)
    endif

    !---------------------------------------------------------------------------------
    !      return the current filename 
    !---------------------------------------------------------------------------------
    incr_filename = trim(fn_new)
    if ( mpi_rank()==0 ) print*, 'incr_flnm: new filename = ',trim(incr_filename)

  end function incr_filename





    subroutine find_times( itms, fids, time, file, datatimem, datatimep, times_found )

    implicit none

    type(trfile), intent(in) :: file
    real(r8), intent(out) :: datatimem, datatimep
    integer, intent(out) :: itms(2) ! record numbers that bracket time
    integer, intent(out) :: fids(2)
    real(r8), intent(in) :: time    ! time of interest
    logical, intent(inout)  :: times_found

    integer :: np1        ! current forward time index of dataset
    integer :: n,i      ! 
    integer :: curr_tsize, next_tsize, all_tsize
    integer :: astat
    integer :: cyc_tsize

    real(r8), allocatable, dimension(:):: all_data_times

    curr_tsize = size(file%curr_data_times)
    next_tsize = 0
    if ( allocated(file%next_data_times)) next_tsize = size(file%next_data_times)

    all_tsize = curr_tsize + next_tsize

    allocate( all_data_times( all_tsize ) )

    all_data_times(:curr_tsize) = file%curr_data_times(:)
    if (next_tsize > 0) all_data_times(curr_tsize+1:all_tsize) = file%next_data_times(:)

    if ( .not. file%cyclical ) then
       if ( all( all_data_times(:) > time ) .and. mpi_rank()==0) then
          print*, 'FIND_TIMES: ALL data times are after ', time
          print*, 'FIND_TIMES: data times: ',all_data_times(:)
          print*, 'FIND_TIMES: time: ',time
          call endrun('find_times: all(all_data_times(:) > time) '// trim(file%curr_filename) )
       endif

       ! find bracketing times 
       find_times_loop : do n=1, all_tsize-1
          np1 = n + 1
          datatimem = all_data_times(n)   !+ file%offset_time
          datatimep = all_data_times(np1) !+ file%offset_time
       ! When stepTime, datatimep may not equal the time (as only datatimem is used)
       ! Should not break other runs?
          if ( (time .ge. datatimem) .and. (time .lt. datatimep) ) then
             times_found = .true.
             exit find_times_loop
          endif
       enddo find_times_loop

    else  ! file%cyclical

       cyc_tsize = file%cyc_ndx_end - file%cyc_ndx_beg + 1

       if ( cyc_tsize > 1 ) then

          call findplb(all_data_times(file%cyc_ndx_beg:file%cyc_ndx_end),cyc_tsize, time, n )

          if (n == cyc_tsize) then
             np1 = 1
          else
             np1 = n+1
          endif

          datatimem = all_data_times(n  +file%cyc_ndx_beg-1)   
          datatimep = all_data_times(np1+file%cyc_ndx_beg-1) 
          times_found = .true.

       endif
    endif

    if ( .not. times_found ) then
       if (mpi_rank()==0) then
          print*,'FIND_TIMES: Failed to find dates bracketing desired time =', time
          print*,' datatimem = ',file%datatimem
          print*,' datatimep = ',file%datatimep
          print*,' all_data_times = ',all_data_times
          !call endrun()
          return
       endif
    endif

    deallocate( all_data_times )
 
    if ( .not. file%cyclical ) then
      itms(1) = n
      itms(2) = np1
    else
      itms(1) = n   +file%cyc_ndx_beg-1
      itms(2) = np1 +file%cyc_ndx_beg-1
    endif

    fids(:) = file%curr_fileid

    do i=1,2
       if ( itms(i) > curr_tsize ) then 
          itms(i) = itms(i) - curr_tsize 
          fids(i) = file%next_fileid
       endif
    enddo

    end subroutine find_times


    subroutine read_next_trcdata(nstep, dtime , ncol, grid_lon, grid_lat,   &
                                 flds, file)
    implicit none
    integer,             intent(in)    :: nstep
    integer,             intent(in)    :: ncol
    real(r8),            intent(in)    :: dtime
    real(r8),            intent(in)    :: grid_lon(:), grid_lat(:)
    type (trfile), intent(inout)       :: file
    type (trfld),intent(inout)         :: flds(:)

    integer :: recnos(4),i,f,nflds      ! 
    integer :: cnt4(4)            ! array of counts for each dimension
    integer :: strt4(4)           ! array of starting indices
    integer :: cnt3(3)            ! array of counts for each dimension
    integer :: strt3(3)           ! array of starting indices
    integer :: fids(4)
    logical :: times_found 

    integer  :: cur_yr, cur_mon, cur_day, cur_sec, yr1, yr2, mon, date, sec
    real(r8) :: series1_time, series2_time
    integer  :: fid1, fid2
   
    !----------------LiXH confirms lon varying between 0-2pi-----------
    real(r8) :: loc_lon(ncol), loc_lat(ncol)


    loc_lon = grid_lon
    loc_lat = grid_lat
    do i = 1, ncol
        if(loc_lon(i) .lt. 0.)loc_lon(i) = loc_lon(i) + 2.*pi
    end do
    !----------------LiXH confirms lon varying between 0-2pi-----------

    nflds = size(flds)
    times_found = .false.
    
    do while( .not. times_found )
    call find_times( recnos(1:2), fids(1:2), file%curr_mod_time, file, file%datatimem, file%datatimep, times_found )
       if ( .not. times_found ) then
           call check_files( file, fids(1:2), recnos(1:2), times_found )
       endif
    enddo
          
    !--------------------------------------------------------------
    !       If stepTime, then no time interpolation is to be done
    !--------------------------------------------------------------
    if (file%stepTime) then
       file%interp_recs = 1
    else
       file%interp_recs = 2
    end if

    if ( file%fill_in_months ) then     !LiXH Test: F

       if( file%datatimep-file%datatimem > file%one_yr ) then

          call get_curr_date(start_ymd, start_tod, nstep, dtime, cur_yr, cur_mon, cur_day, cur_sec)

          call set_date_from_time_float(file%datatimem, yr1, mon, date, sec )
          call set_date_from_time_float(file%datatimep, yr2, mon, date, sec )

          call set_time_float_from_date(start_ymd, start_tod, series1_time, yr1, cur_mon, cur_day, cur_sec )
          call set_time_float_from_date(start_ymd, start_tod, series2_time, yr2, cur_mon, cur_day, cur_sec )

          fid1 = fids(1)
          fid2 = fids(2)
          file%cyclical = .true.
          call set_cycle_indices( fid1, file%cyc_ndx_beg, file%cyc_ndx_end, yr1)
          call find_times( recnos(1:2), fids(1:2), series1_time, file, file%datatimes(1), file%datatimes(2), times_found )

          if ( .not. times_found ) then
              call endrun('read_next_trcdata: time not found for series1_time')
          endif
          call set_cycle_indices( fid2, file%cyc_ndx_beg, file%cyc_ndx_end, yr2)

         !--------------LiXH is not sure what fid1%fh means, seems like file handle, must check!---------------- 
         ! if ( fid1%fh /= fid2%fh ) then
          if ( fid1 /= fid2 ) then
            file%cyc_ndx_beg = file%cyc_ndx_beg + size(file%curr_data_times)
            file%cyc_ndx_end = file%cyc_ndx_end + size(file%curr_data_times)
          endif
         !--------------LiXH is not sure what fid1%fh means, seems like file handle, must check!---------------- 

          call find_times( recnos(3:4), fids(3:4), series2_time, file, file%datatimes(3), file%datatimes(4), times_found )
          if ( .not. times_found ) then
              call endrun('read_next_trcdata: time not found for series2_time')
          endif
          file%cyclical = .false.
          file%interp_recs = 4

          call set_date_from_time_float(file%datatimes(1), yr1, mon, date, sec )
          call set_time_float_from_date(start_ymd, start_tod, file%datatimem, cur_yr,  mon, date, sec )

          if (file%datatimes(1) > file%datatimes(2) ) then ! wrap around
            if ( cur_mon == 1 ) then 
               call set_time_float_from_date(start_ymd, start_tod, file%datatimem, cur_yr-1,  mon, date, sec )
            endif 
          endif

          call set_date_from_time_float(file%datatimes(2), yr1, mon, date, sec )
          call set_time_float_from_date(start_ymd, start_tod, file%datatimep, cur_yr,  mon, date, sec )

          if (file%datatimes(1) > file%datatimes(2) ) then ! wrap around
            if ( cur_mon == 12 ) then
              call set_time_float_from_date(start_ymd, start_tod, file%datatimep, cur_yr+1,  mon, date, sec )
            endif 
          endif

       endif

    endif

    !
    ! Set up hyperslab corners
    !
    do i=1,file%interp_recs

       strt4(:) = 1
       strt3(:) = 1

       do f = 1,nflds
          if ( file%zonal_ave ) then        !LiXH Test: F
             cnt3(flds(f)%coords(ZA_LATDIM)) = file%nlat
             if (flds(f)%srf_fld) then
                cnt3(flds(f)%coords(ZA_LEVDIM)) = 1
             else
                cnt3(flds(f)%coords(ZA_LEVDIM)) = file%nlev
             endif
             cnt3(flds(f)%coords(ZA_TIMDIM)) = 1
             strt3(flds(f)%coords(ZA_TIMDIM)) = recnos(i)
             call read_za_trc( fids(i), flds(f)%var_id, ncol, loc_lon, loc_lat, flds(f)%input(i)%data, strt3, cnt3, file, &
                  (/ flds(f)%order(ZA_LATDIM),flds(f)%order(ZA_LEVDIM) /) )
          else if ( flds(f)%srf_fld ) then
             cnt3( flds(f)%coords(LONDIM)) = file%nlon
             cnt3( flds(f)%coords(LATDIM)) = file%nlat
             cnt3( flds(f)%coords(PS_TIMDIM)) = 1
             strt3(flds(f)%coords(PS_TIMDIM)) = recnos(i)
             call read_2d_trc( fids(i), flds(f)%var_id, ncol, loc_lon, loc_lat, flds(f)%input(i)%data(1,:), strt3, cnt3, file, &
                 (/ flds(f)%order(LONDIM),flds(f)%order(LATDIM) /) )
          else
             cnt4(flds(f)%coords(LONDIM)) = file%nlon
             cnt4(flds(f)%coords(LATDIM)) = file%nlat
             cnt4(flds(f)%coords(LEVDIM)) = file%nlev
             cnt4(flds(f)%coords(TIMDIM)) = 1
             strt4(flds(f)%coords(TIMDIM)) = recnos(i)

             call read_3d_trc( fids(i), flds(f)%var_id, ncol, loc_lon, loc_lat, flds(f)%input(i)%data, strt4, cnt4, file, &
                  (/ flds(f)%order(LONDIM),flds(f)%order(LATDIM),flds(f)%order(LEVDIM) /))

          endif

       enddo

       if ( file%has_ps ) then      !LiXH Test: T
          cnt3 = 1
          strt3 = 1
          if (.not. file%zonal_ave) then
             cnt3(file%ps_coords(LONDIM)) = file%nlon
          end if
          cnt3(file%ps_coords(LATDIM)) = file%nlat
          cnt3(file%ps_coords(PS_TIMDIM)) = 1
          strt3(file%ps_coords(PS_TIMDIM)) = recnos(i)
          if (file%zonal_ave) then
             call read_2d_trc( fids(i), file%ps_id, ncol, loc_lon, loc_lat, file%ps_in(i)%data, strt3(1:2), cnt3(1:2), file,&
                  (/ 1, 2 /) )
          else
             call read_2d_trc( fids(i), file%ps_id, ncol, loc_lon, loc_lat, file%ps_in(i)%data, strt3, cnt3, file, &
                  (/ file%ps_order(LONDIM),file%ps_order(LATDIM) /) )
          end if
       endif

    enddo

    end subroutine read_next_trcdata


    subroutine read_2d_trc( fid, vid, loc_ncol, loc_lon, loc_lat, loc_arr, strt, cnt, file, order )

    use grist_horizontal_interpolate, only : xy_interp,     &
                                             lininterp_init, lininterp, interp_type, lininterp_finish


    implicit none
    integer, intent(in) :: fid
    integer, intent(in) :: vid
    integer, intent(in) :: loc_ncol
    integer, intent(in) :: strt(:), cnt(:), order(2)
    real(r8),intent(in) :: loc_lon(:), loc_lat(:)
    real(r8),intent(out)  :: loc_arr(:)
    type (trfile), intent(in) :: file

    !local
    real(r8), allocatable, target :: wrk2d(:,:)
    real(r8), pointer :: wrk2d_in(:,:)

    integer :: ncols
    real(r8), parameter :: zero=0_r8, twopi=2_r8*pi
    type(interp_type) :: lon_wgts, lat_wgts    
    real(r8) :: file_lats(file%nlat)

    nullify(wrk2d_in)
    allocate( wrk2d(cnt(1),cnt(2)) )
 
    if(order(1)/=1 .or. order(2)/=2 .or. cnt(1)/=file%nlon .or. cnt(2)/=file%nlat) then
       allocate( wrk2d_in(file%nlon, file%nlat) )
    end if
 
 
    call  wrap_get_vara_realx( fid, vid, strt, cnt, wrk2d )
    if(associated(wrk2d_in)) then
       wrk2d_in = reshape( wrk2d(:,:),(/file%nlon,file%nlat/), order=order )
       deallocate(wrk2d)
    else
       wrk2d_in => wrk2d
    end if

    ! PGI 13.9 bug workaround.
    file_lats = file%lats

    ! For zonal average, only interpolate along latitude.
    if (file%zonal_ave) then
       ncols = loc_ncol
       call lininterp_init(file_lats, file%nlat, loc_lat, ncols, 1, lat_wgts)
       call lininterp(wrk2d_in(1,:), file%nlat, loc_arr(1:ncols), ncols, lat_wgts)
       call lininterp_finish(lat_wgts)
    else
       ! if weighting by latitude, the perform horizontal interpolation by using weight_x, weight_y
       if(file%weight_by_lat) then
          ncols = loc_ncol
          call xy_interp(file%nlon,file%nlat,1,ncols,file%weight_x,file%weight_y,wrk2d_in,loc_arr(1:ncols),  &
                         file%index_x,file%index_y) 
       else
          ncols = loc_ncol
          call lininterp_init(file%lons, file%nlon, loc_lon, ncols, 2, lon_wgts, zero, twopi)
          call lininterp_init(file%lats, file%nlat, loc_lat, ncols, 1, lat_wgts)
          call lininterp(wrk2d_in, file%nlon, file%nlat, loc_arr(1:ncols), ncols, lon_wgts, lat_wgts)    
       
          call lininterp_finish(lon_wgts)
          call lininterp_finish(lat_wgts)
       endif

    end if

    if(allocated(wrk2d)) then
       deallocate(wrk2d)
    else
       deallocate(wrk2d_in)
    end if
  
    end subroutine read_2d_trc


    subroutine read_za_trc( fid, vid, loc_ncol, loc_lon, loc_lat, loc_arr, strt, cnt, file, order )
    use grist_horizontal_interpolate, only : lininterp_init, lininterp, interp_type, lininterp_finish

    implicit none
    integer, intent(in) :: fid
    integer, intent(in) :: vid
    integer, intent(in) :: loc_ncol
    real(r8),intent(in) :: loc_lon(:), loc_lat(:)
    integer, intent(in) :: strt(:), cnt(:)
    integer, intent(in) :: order(2)
    real(r8),intent(out):: loc_arr(:,:)
    type (trfile), intent(in) :: file

    type(interp_type) :: lat_wgts
    real(r8) :: wrk(loc_ncol)
    real(r8), allocatable, target :: wrk2d(:,:)
    real(r8), pointer :: wrk2d_in(:,:)
    integer :: ncols, k

     nullify(wrk2d_in)
     allocate( wrk2d(cnt(1),cnt(2)) )

     if(order(1)/=1 .or. order(2)/=2 .or. cnt(1)/=file%nlat .or. cnt(2)/=file%nlev) then
        allocate( wrk2d_in(file%nlat, file%nlev) )
     end if


    call wrap_get_vara_realx( fid, vid, strt, cnt, wrk2d )
    if(associated(wrk2d_in)) then
       wrk2d_in = reshape( wrk2d(:,:),(/file%nlat,file%nlev/), order=order )
       deallocate(wrk2d)
    else
       wrk2d_in => wrk2d
    end if

    ncols = loc_ncol

    call lininterp_init(file%lats, file%nlat, loc_lat, ncols, 1, lat_wgts)
    do k=1,file%nlev
       call lininterp(wrk2d_in(:,k), file%nlat, wrk(1:ncols), ncols, lat_wgts)    
       loc_arr(k,1:ncols) = wrk(1:ncols)
    end do
    call lininterp_finish(lat_wgts)

    if(allocated(wrk2d)) then
       deallocate(wrk2d)
    else
       deallocate(wrk2d_in)
    end if
    
    end subroutine read_za_trc

    subroutine read_3d_trc( fid, vid, loc_ncol, loc_lon, loc_lat, loc_arr, strt, cnt, file, order)
    use grist_horizontal_interpolate, only : xy_interp, &
                                             lininterp_init, lininterp, interp_type, lininterp_finish

    implicit none
    integer, intent(in) :: fid
    integer, intent(in) :: vid
    integer, intent(in) :: loc_ncol
    real(r8),intent(in) :: loc_lon(:), loc_lat(:)
    integer, intent(in) :: strt(:), cnt(:), order(3)
    real(r8),intent(out)  :: loc_arr(:,:)
    type (trfile),     intent(in) :: file

    integer :: i,j,k, astat, c, ncols

    integer                     :: jlim(2), jl, ju, ierr
    integer                     :: gndx

    real(r8), allocatable, target :: wrk3d(:,:,:)
    real(r8), pointer :: wrk3d_in(:,:,:)
    real(r8), parameter :: zero=0._r8, twopi=2._r8*pi
    type(interp_type) :: lon_wgts, lat_wgts    

    loc_arr(:,:) = 0._r8
    nullify(wrk3d_in)
    allocate(wrk3d(cnt(1),cnt(2),cnt(3)))

    call wrap_get_vara_realx( fid, vid, strt, cnt, wrk3d )

    if(order(1)/=1 .or. order(2)/=2 .or. order(3)/=3 .or. &
         cnt(1)/=file%nlon.or.cnt(2)/=file%nlat.or.cnt(3)/=file%nlev) then
       allocate(wrk3d_in(file%nlon,file%nlat,file%nlev))
       wrk3d_in = reshape( wrk3d(:,:,:),(/file%nlon,file%nlat,file%nlev/), order=order )
       deallocate(wrk3d)
    else
       wrk3d_in => wrk3d
    end if

    j=1

! If weighting by latitude, then perform horizontal interpolation by using weight_x, weight_y

   if(file%weight_by_lat) then

       ncols = loc_ncol
       call xy_interp(file%nlon,file%nlat,file%nlev,ncols,file%weight_x,file%weight_y,wrk3d_in, &
            loc_arr(:,:), file%index_x,file%index_y) 

   else
       ncols = loc_ncol
       call lininterp_init(file%lons, file%nlon, loc_lon, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(file%lats, file%nlat, loc_lat, ncols, 1, lat_wgts)

       call lininterp(wrk3d_in, file%nlon, file%nlat, file%nlev, loc_arr(:,:), ncols, lon_wgts, lat_wgts)

       call lininterp_finish(lon_wgts)
       call lininterp_finish(lat_wgts)
   endif


    if(allocated(wrk3d)) then
       deallocate( wrk3d )
    else
       deallocate( wrk3d_in )
    end if
    end subroutine read_3d_trc


    subroutine interpolate_trcdata( ncol, flds, file )
    use mo_util,      only : rebin

    implicit none
    !io
    integer,             intent(in)    :: ncol
    type (trfld),        intent(inout) :: flds(:)
    type (trfile),       intent(inout) :: file

    !local
    real(r8) :: fact1, fact2
    real(r8) :: deltat
    integer :: f,nflds, i,k
    real(r8) :: ps(ncol)
    real(r8) :: datain(file%nlev, ncol)
    real(r8) :: pin(file%nlev, ncol)
    real(r8)            :: model_z(nlevp)
    real(r8), parameter :: m2km  = 1.e-3_r8
    real(r8), parameter :: cday  = 86400._r8
    real(r8), dimension(:,:), allocatable :: data_out

    nflds = size(flds)

    if ( file%interp_recs == 4 ) then
       deltat = file%datatimes(3) - file%datatimes(1)
       fact1 = (file%datatimes(3) - file%datatimem)/deltat
       fact2 = 1._r8-fact1
       if ( file%has_ps ) then
          file%ps_in(1)%data(:ncol) = fact1*file%ps_in(1)%data(:ncol) + fact2*file%ps_in(3)%data(:ncol) 
       endif
       do f = 1,nflds
       flds(f)%input(1)%data(:,:ncol) = fact1*flds(f)%input(1)%data(:,:ncol) + fact2*flds(f)%input(3)%data(:,:ncol) 
       enddo

       deltat = file%datatimes(4) - file%datatimes(2)
       fact1 = (file%datatimes(4) - file%datatimep)/deltat
       fact2 = 1._r8-fact1

       if ( file%has_ps ) then
          file%ps_in(2)%data(:ncol) = fact1*file%ps_in(2)%data(:ncol) + fact2*file%ps_in(4)%data(:ncol) 
       endif
       do f = 1,nflds
       flds(f)%input(2)%data(:,:ncol) = fact1*flds(f)%input(2)%data(:,:ncol) + fact2*flds(f)%input(4)%data(:,:ncol) 
       enddo

    endif
    !-------------------------------------------------------------------------
    !       If file%interp_recs=1 then no time interpolation -- set
    !       fact1=1 and fact2=0 and will just use first value unmodified
    !-------------------------------------------------------------------------

    if (file%interp_recs == 1) then
       fact1=1._r8
       fact2=0._r8
    else
       file%interp_recs = 2

       deltat = file%datatimep - file%datatimem

       if ( file%cyclical .and. (deltat < 0._r8) ) then
          deltat = deltat+file%one_yr
          if ( file%datatimep >= file%curr_mod_time ) then
             fact1 = (file%datatimep - file%curr_mod_time)/deltat
          else
             fact1 = (file%datatimep+file%one_yr - file%curr_mod_time)/deltat
          endif
       else
             fact1 = (file%datatimep - file%curr_mod_time)/deltat
       endif

       ! this assures that FIXED data are b4b on restarts
       if ( file%fixed ) then
          fact1 = dble(int(fact1*cday+.5_r8))/dble(cday)
       endif
       fact2 = 1._r8-fact1
    endif


    fld_loop: do f = 1,nflds

       allocate(data_out(nlev,ncol)); data_out=0._r8

       if (file%alt_data) then

          if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
              datain(:,:ncol) = fact1*flds(f)%input(nm)%data(:,:ncol)
          else
              datain(:,:ncol) = fact1*flds(f)%input(nm)%data(:,:ncol) + fact2*flds(f)%input(np)%data(:,:ncol) 
          end if
          do i = 1,ncol
             model_z(1:nlevp) = m2km * pstate%z_at_pc_face_level%f(nlevp:1:-1,i)
             call rebin( file%nlev, nlev, file%ilevs, model_z, datain(:,i), data_out(:,i) )
          enddo

       else

          if ( file%nlev>1 ) then
             if ( file%has_ps ) then
                if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
                   ps(:ncol) = fact1*file%ps_in(nm)%data(:ncol)
                else
                   ps(:ncol) = fact1*file%ps_in(nm)%data(:ncol) + fact2*file%ps_in(np)%data(:ncol) 
                end if
                do i = 1,ncol
                   do k = 1,file%nlev
                      pin(k,i) = file%p0*file%hyam(k) + ps(i)*file%hybm(k)
                   enddo
                enddo
             else
                do k = 1,file%nlev
                   pin(k,:) = file%levs(k)
                enddo
             endif
          endif

          if (flds(f)%srf_fld) then
             do i = 1,ncol
                if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
                   data_out(1,i) = &
                        fact1*flds(f)%input(nm)%data(1,i)
                else
                   data_out(1,i) = &
                        fact1*flds(f)%input(nm)%data(1,i) + fact2*flds(f)%input(np)%data(1,i) 
                endif
             enddo

          else
             if (fact2 == 0) then  ! This needed as %data is not set if fact2=0 (and lahey compiler core dumps)
                 datain(:,:ncol) = fact1*flds(f)%input(nm)%data(:,:ncol)
             else
                 datain(:,:ncol) = fact1*flds(f)%input(nm)%data(:,:ncol) + fact2*flds(f)%input(np)%data(:,:ncol)
             end if
             if ( file%top_bndry ) then
                 call vert_interp_ub(ncol, file%nlev, file%levs,  datain(:,:ncol), data_out(1,:ncol) )
             else if(file%conserve_column) then
                call vert_interp_mixrat(ncol,file%nlev,nlev,pstate%pressure_at_pc_face_level%f(:,1:ncol), &
                     datain, data_out(:,:), &
                     file%p0,ps,file%hyai,file%hybi)
             else
                call vert_interp(ncol, file%nlev, pin, pstate%pressure_at_pc_full_level%f(:,1:ncol), datain, data_out(:,:) )
             endif
          endif

       endif

       if (flds(f)%pbuf_ndx<=0) then
          flds(f)%data = data_out
       else
          pstate_cam%aerosol_at_pc_full_level(flds(f)%pbuf_ndx)%f(:,1:ncol) = data_out 
       endif

       deallocate(data_out)
    enddo fld_loop

    end subroutine interpolate_trcdata


    subroutine get_dimension( fileid, dname, dsize, dimid, data )
    ! io
    character(*), intent(in) :: dname
    integer, intent(in)  :: fileid
    integer, intent(out) :: dsize
    integer, optional, intent(out) :: dimid
    real(r8), optional, allocatable :: data(:) 
    ! local
    integer :: vid, ierr, id

    ierr = nf_inq_dimid(fileid, dname, id)

    if ( ierr==NF_NOERR ) then

        call wrap_inq_dimlen(fileid, id, dsize)

       if ( present(dimid) ) then
          dimid = id
       endif

       if ( present(data) ) then
           if ( allocated(data) )deallocate(data)
           allocate( data(dsize) )

           call wrap_inq_varid (fileid, dname, vid)
           call wrap_get_var_realx (fileid, vid, data)
       end if 

    else
       dsize = 1
       if ( present(dimid) ) then
          dimid = -1
       endif
    endif

    end subroutine get_dimension


    subroutine set_cycle_indices( fileid, cyc_ndx_beg, cyc_ndx_end, cyc_yr )

    implicit none

    integer, intent(in)  :: fileid
    integer, intent(out) :: cyc_ndx_beg
    integer, intent(out) :: cyc_ndx_end
    integer, intent(in)  :: cyc_yr

    integer, allocatable , dimension(:) :: dates, datesecs
    integer :: timesize, i, astat, year, ierr
    integer :: dateid

    call get_dimension( fileid, 'time', timesize )
    cyc_ndx_beg=-1

    allocate( dates(timesize) )

    call wrap_inq_varid(   fileid, 'date',  dateid  )
    call wrap_get_var_int( fileid, dateid, dates )

    do i=1,timesize
       year = dates(i) / 10000
       if ( year == cyc_yr ) then
          if  (cyc_ndx_beg < 0)  then
             cyc_ndx_beg = i
          endif
          cyc_ndx_end = i
       endif
    enddo
    deallocate( dates )

    if (cyc_ndx_beg < 0) then
       print*, 'set_cycle_indices: cycle year not found : ' , cyc_yr
       call endrun('set_cycle_indices: cycle year not found')
    endif

    end subroutine set_cycle_indices


    subroutine open_trc_datafile( fname, path, fileid, times, cyc_ndx_beg, cyc_ndx_end, cyc_yr )
    ! io
    character(*), intent(in)  :: fname
    character(*), intent(in)  :: path
    integer,      intent(out) :: fileid
    real(r8), allocatable,  intent(out) :: times(:)

    integer, optional, intent(out) :: cyc_ndx_beg
    integer, optional, intent(out) :: cyc_ndx_end
    integer, optional, intent(in) :: cyc_yr

    character(len=256)  :: filepath
    integer :: year, month, day, dsize, i, timesize
    integer :: dateid,secid
    integer, allocatable , dimension(:) :: dates, datesecs
    integer :: astat, ierr

    integer(i4), parameter            :: omode = 0
    integer(i4)                       :: dimid

    if (len_trim(path) == 0) then
       filepath = trim(fname)
    else
       filepath = trim(path) // '/' // trim(fname)
    end if

    ! open file and get fileid
    if(log_print .and. mpi_rank()==0)print*,'open_aero_datafile: ',trim(filepath)
    call wrap_open(trim(filepath), omode, fileid)

    call get_dimension(fileid, 'time', timesize)

    if ( allocated(times) )deallocate(times)
    allocate( times(timesize) )

    allocate( dates(timesize) )
    allocate( datesecs(timesize) )

    ierr = nf_inq_varid (fileid, 'date', dateid)
    call wrap_get_var_int (fileid, dateid, dates)
    
    ierr = nf_inq_varid (fileid, 'datesec', secid)

    if (ierr==NF_NOERR)then
       call wrap_get_var_int (fileid, secid, datesecs)
    else
       datesecs=0
    end if

    do i=1,timesize
       year = dates(i) / 10000
       month = mod(dates(i),10000)/100
       day = mod(dates(i),100)
       call set_time_float_from_date(start_ymd, start_tod, times(i), year, month, day, datesecs(i) )
       if ( present(cyc_yr) ) then
          if ( year == cyc_yr ) then
             if ( present(cyc_ndx_beg) .and. (cyc_ndx_beg < 0) ) then
                cyc_ndx_beg = i
             endif
             if ( present(cyc_ndx_end) ) then
                cyc_ndx_end = i
             endif
          endif
       endif
    enddo

    deallocate( dates )
    deallocate( datesecs )       
       
    if ( present(cyc_yr) .and. present(cyc_ndx_beg) ) then
       if (cyc_ndx_beg < 0) then
           if(mpi_rank()==0)print*, 'open_trc_datafile: cycle year not found : ' , cyc_yr
          call endrun('open_trc_datafile: cycle year not found')
       endif
    endif

    end subroutine open_trc_datafile


    subroutine specify_fields( specifier, fields )

    implicit none

    character(len=*), intent(in) :: specifier(:)
    type(trfld), allocatable, intent(inout) :: fields(:)

    integer :: fld_cnt
    integer :: i,j
    character(len=256) :: str1, str2
    character(len=32), allocatable, dimension(:) :: fld_name,  src_name
    integer :: nflds

    nflds = size(specifier)

    allocate(fld_name(nflds),  src_name(nflds))

    fld_cnt = 0

    count_cnst: do i = 1, nflds

       if ( len_trim( specifier(i) ) == 0 ) then
          exit count_cnst
       endif

       j = scan( specifier(i),':')

       if (j > 0) then
          str1 = trim(adjustl( specifier(i)(:j-1) ))
          str2 = trim(adjustl( specifier(i)(j+1:) ))
          fld_name(i) = trim(adjustl( str1 ))
          src_name(i) = trim(adjustl( str2 ))
       else
          fld_name(i) = trim(adjustl( specifier(i) ))
          src_name(i) = trim(adjustl( specifier(i) ))
       endif

       fld_cnt = fld_cnt + 1

    enddo count_cnst

    if( fld_cnt < 1 ) then
       if(allocated(fields))deallocate(fields) 
       return
    end if

    !-----------------------------------------------------------------------
    ! 	... allocate field type array
    !-----------------------------------------------------------------------
    if(allocated(fields))deallocate(fields)
    allocate( fields(fld_cnt))

    do i = 1,fld_cnt
       fields(i)%fldnam = fld_name(i)
       fields(i)%srcnam = src_name(i)
    enddo

    deallocate(fld_name, src_name)

    end subroutine specify_fields

    
    subroutine vert_interp_mixrat( ncol, nsrc, ntrg, trg_x, src, trg, p0, ps, hyai, hybi)
  
    implicit none

    integer, intent(in)   :: ncol 
    integer, intent(in)   :: nsrc                  ! dimension source array
    integer, intent(in)   :: ntrg                  ! dimension target array
    real(r8)              :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)      :: trg_x(ntrg+1,ncol)         ! target coordinates
    real(r8), intent(in)      :: src(nsrc,ncol)             ! source array
    real(r8), intent(out)     :: trg(ntrg,ncol)             ! target array

    real(r8) :: ps(ncol), p0, hyai(nsrc+1), hybi(nsrc+1)
    !---------------------------------------------------------------
    !   ... local variables
    !---------------------------------------------------------------
    integer  :: i, j, n
    integer  :: sil
    real(r8)     :: tl, y
    real(r8)     :: bot, top
   
    do n = 1,ncol
    
    do i=1,nsrc+1
     src_x(i) = p0*hyai(i)+ps(n)*hybi(i)
    enddo

    do i = 1, ntrg
       tl = trg_x(i+1,n)
       if( (tl.gt.src_x(1)).and.(trg_x(i,n).lt.src_x(nsrc+1)) ) then
          do sil = 1,nsrc
             if( (tl-src_x(sil))*(tl-src_x(sil+1)).le.0.0_r8 ) then
                exit
             end if
          end do

          if( tl.gt.src_x(nsrc+1)) sil = nsrc

          y = 0.0_r8
          bot = min(tl,src_x(nsrc+1))   
          top = trg_x(i,n)
          do j = sil,1,-1
           if( top.lt.src_x(j) ) then
             y = y+(bot-src_x(j))*src(j,n)
            bot = src_x(j)
           else
            y = y+(bot-top)*src(j,n)
            exit
           endif
          enddo
          trg(i,n) = y
       else
        trg(i,n) = 0.0_r8
       end if
    end do

    if( trg_x(ntrg+1,n).lt.src_x(nsrc+1) ) then
     top = trg_x(ntrg+1,n)
     bot = src_x(nsrc+1)
     y = 0.0_r8
     do j=nsrc,1,-1
      if( top.lt.src_x(j) ) then
       y = y+(bot-src_x(j))*src(j,n)
       bot = src_x(j)
      else
       y = y+(bot-top)*src(j,n)
       exit
      endif
     enddo
     trg(ntrg,n) = trg(ntrg,n)+y
    endif

! turn mass into mixing ratio 
    do i=1,ntrg
     trg(i,n) = trg(i,n)/(trg_x(i+1,n)-trg_x(i,n))
    enddo
    
    enddo

   end subroutine vert_interp_mixrat


  subroutine vert_interp( ncol, levsiz, pin, pmid, datain, dataout )
    !-------------------------------------------------------------------------- 
    ! 
    ! Interpolate data from current time-interpolated values to model levels
    !--------------------------------------------------------------------------
    implicit none
    ! Arguments
    !
    integer,  intent(in)  :: ncol                ! number of atmospheric columns
    integer,  intent(in)  :: levsiz
    real(r8), intent(in)  :: pin(levsiz,ncol)
    real(r8), intent(in)  :: pmid(nlev,ncol)          ! level pressures 
    real(r8), intent(in)  :: datain(levsiz,ncol)
    real(r8), intent(out) :: dataout(nlev,ncol)     

    !
    ! local storage
    !

    integer ::  i                  ! longitude index
    integer ::  k, kk, kkstart     ! level indices
    integer ::  kupper(ncol)       ! Level indices for interpolation
    real(r8) :: dpu                ! upper level pressure difference
    real(r8) :: dpl                ! lower level pressure difference



    !--------------------------------------------------------------------------
    !
    ! Initialize index array
    !
    do i=1,ncol
       kupper(i) = 1
    end do

    do k=1,nlev
       !
       ! Top level we need to start looking is the top level for the previous k
       ! for all column points
       !
       kkstart = levsiz
       do i=1,ncol
          kkstart = min0(kkstart,kupper(i))
       end do
       !
       ! Store level indices for interpolation
       !
       do kk=kkstart,levsiz-1
          do i=1,ncol
             if (pin(kk,i).lt.pmid(k,i) .and. pmid(k,i).le.pin(kk+1,i)) then
                kupper(i) = kk
             end if
          end do
       end do
       ! interpolate or extrapolate...
       do i=1,ncol
          if (pmid(k,i) .lt. pin(1,i)) then
             dataout(k,i) = datain(1,i)*pmid(k,i)/pin(1,i)
          else if (pmid(k,i) .gt. pin(levsiz,i)) then
             dataout(k,i) = datain(levsiz,i)
          else
             dpu = pmid(k,i) - pin(kupper(i),i)
             dpl = pin(kupper(i)+1,i) - pmid(k,i)
             dataout(k,i) = (datain(kupper(i),i)*dpl + &
                  datain(kupper(i)+1,i)*dpu)/(dpl + dpu)
          end if
       end do
    end do


    end subroutine vert_interp


    subroutine vert_interp_ub( ncol, nlevs, plevs,  datain, dataout )
    !----------------------------------------------------------------------- 
    ! 
    ! Interpolate data from current time-interpolated values to top interface pressure
    !  -- from mo_tgcm_ubc.F90
    !--------------------------------------------------------------------------
    implicit none
    ! Arguments
    !
    integer,  intent(in)  :: ncol
    integer,  intent(in)  :: nlevs
    real(r8), intent(in)  :: plevs(nlevs)
    real(r8), intent(in)  :: datain(nlevs,ncol)
    real(r8), intent(out) :: dataout(ncol)   

    !
    ! local variables
    !
    integer  :: i,ku,kl,kk
    real(r8) :: pinterp, delp
    
    do i = 1,ncol

    pinterp = pstate%pressure_at_pc_face_level%f(1,i)

    if( pinterp <= plevs(1) ) then
       kl = 1
       ku = 1
       delp = 0._r8
    else if( pinterp >= plevs(nlevs) ) then
       kl = nlevs
       ku = nlevs
       delp = 0._r8
    else

       do kk = 2,nlevs
          if( pinterp <= plevs(kk) ) then
             ku = kk
             kl = kk - 1
             delp = log( pinterp/plevs(kk) ) / log( plevs(kk-1)/plevs(kk) )
             exit
          end if
       end do
    end if

    dataout(i) = datain(kl,i) + delp * (datain(ku,i) - datain(kl,i))

    end do

    end subroutine vert_interp_ub


    subroutine set_date_from_time_float(curr_mod_time, yr, mn, day, nsec)
    !io
    real(r8), intent(in)      :: curr_mod_time
    integer,  intent(out)     :: yr, mn, day, nsec
    ! local
    integer  :: month(12)=(/31,28,31,30,31,30,31,31,30,31,30,31/)
    integer  :: calday, i

    nsec = int((curr_mod_time - floor(curr_mod_time))*86400.)
    if(nsec .lt. 0)call endrun('Error in set_date_from_time_float')

    yr = floor(curr_mod_time)/365+1
    calday = mod(floor(curr_mod_time), 365)

    mn = 1
    do i = 1, 12
        if(calday .gt. month(i))then
            calday = calday - month(i)
            mn = mn+1
        else
            exit
        end if
    end do

    day = calday


    end subroutine set_date_from_time_float


    subroutine advance_file(file)

    !------------------------------------------------------------------------------
    !   This routine advances to the next file
    !------------------------------------------------------------------------------

    use ioFileMod, only: getfil

    implicit none

    type(trfile), intent(inout) :: file

    !-----------------------------------------------------------------------
    !   local variables
    !-----------------------------------------------------------------------
    character(len=256) :: ctmp
    character(len=256) :: loc_fname   
    integer            :: istat, astat

    !-----------------------------------------------------------------------
    !   close current file ...
    !-----------------------------------------------------------------------
    call wrap_close( file%curr_fileid )

    !-----------------------------------------------------------------------
    !   remove if requested
    !-----------------------------------------------------------------------
    if( file%remove_trc_file ) then
       call getfil( file%curr_filename, loc_fname, 0 )
       if(mpi_rank()==0)then
          print*, 'advance_file: removing file = ',trim(loc_fname) 
          ctmp = 'rm -f ' // trim(loc_fname) 
          print*, 'advance_file: fsystem issuing command - '
          print*, trim(ctmp)
       end if
     !-------LiXH has not completed shr_sys_system-------
     !  call shr_sys_system( ctmp, istat )
     !-------LiXH has not completed shr_sys_system-------
    end if
   
    !-----------------------------------------------------------------------
    !   Advance the filename and file id
    !-----------------------------------------------------------------------
    file%curr_filename = file%next_filename
    file%curr_fileid = file%next_fileid
   
    !-----------------------------------------------------------------------
    !   Advance the curr_data_times
    !-----------------------------------------------------------------------
    deallocate( file%curr_data_times )
    allocate( file%curr_data_times( size( file%next_data_times ) ) )
    file%curr_data_times(:) = file%next_data_times(:)
    
    !-----------------------------------------------------------------------
    !   delete information about next file (as was just assigned to current)
    !-----------------------------------------------------------------------
    file%next_filename = ''
    
    deallocate( file%next_data_times )

    end subroutine advance_file


    subroutine findplb( x, nx, xval, index )

    !----------------------------------------------------------------------- 
    ! Purpose: 
    ! "find periodic lower bound"
    ! Search the input array for the lower bound of the interval that
    ! contains the input value.  The returned index satifies:
    ! x(index) .le. xval .lt. x(index+1)
    ! Assume the array represents values in one cycle of a periodic coordinate.
    ! So, if xval .lt. x(1), or xval .ge. x(nx), then the index returned is nx.
    !
    ! Author: B. Eaton
    !----------------------------------------------------------------------- 

    implicit none

    integer, intent(in) ::   nx         ! size of x
    real(r8), intent(in) ::  x(nx)      ! strictly increasing array
    real(r8), intent(in) ::  xval       ! value to be searched for in x
    
    integer, intent(out) ::  index

    ! Local variables:
    integer i
    !-----------------------------------------------------------------------

    if ( xval .lt. x(1) .or. xval .ge. x(nx) ) then
       index = nx
       return
    end if

    do i = 2, nx
       if ( xval .lt. x(i) ) then
          index = i-1
          return
       end if
    end do

    end subroutine findplb


 end module grist_tracer_data
