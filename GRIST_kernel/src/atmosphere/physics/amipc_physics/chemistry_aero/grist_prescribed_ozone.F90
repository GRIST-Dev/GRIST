!=======================================================
!
!  Created by LiXiaohan on 20/07/01.
!  manages reading and interpolation of prescribed ozone
!
!=======================================================

 module grist_prescribed_ozone

    use grist_constants,                    only: i4, i8, r8, rad2deg, boltz, mwdry
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_nml_module,                   only: doAquaPlanet, nlev
    use grist_time_manager,                 only: get_curr_date
    use grist_domain_types,                 only: global_domain
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

    real(r8), allocatable   :: lev_data_ape(:)
    real(r8), allocatable   :: ozone_data_ape(:,:)


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

    integer :: o3_idx
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

    do i = 1, pstate_cam%total_ghg_num
        if(trim(pstate_cam%ghg_at_pc_full_level(i)%name) .eq. 'O3')then
            o3_idx = pstate_cam%ghg_at_pc_full_level(i)%idx
        end if
    end do

    if(doAquaPlanet)then
        call ozone_ape(ncol, lat)

    else
        specifier(1) = trim(ozone_name)//':'//trim(fld_name)

        allocate(file%in_pbuf(size(specifier)))
        file%in_pbuf(:) = .false.       !Only aerosol use in_pbuf, LiXH
        call trcdata_init(  dtime, ncol, lon, lat,                                 & 
                            specifier, filename, filelist, datapath, fields, file, &
                            rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)
    
    end if

    end subroutine prescribed_ozone_init

    subroutine ozone_ape(ncol, lat)
! io
    integer(i4), intent(in) :: ncol
    real(r8)   , intent(in) :: lat(ncol)            ! latitude in radian
! local
    integer  :: i, j, k, count, jg, kg, lg
    real(r8) :: w1,w2,dy
    real(r8) :: o3vubc,prsubc,altubc
    character:: label*120
    real(r8) :: lat_deg(ncol)
    real(r8), dimension(:),   allocatable :: log_lev_data, log_lev_model
    real(r8), dimension(:),   allocatable :: glat, gprs, galt, gpri, gali
    real(r8), dimension(:,:), allocatable :: gdat

    lat_deg = lat*rad2deg

    open (01,file=trim(filename),status='old')
    read (01,'(i3)') jg
    read (01,'(i3)') kg
    lg = 1 + kg
 
    ozone_nlev = kg
    if(.not.allocated(lev_data_ape)) allocate(lev_data_ape(ozone_nlev))
    if(.not.allocated(ozone_data_ape)) allocate(ozone_data_ape(ozone_nlev, ncol))

    allocate(glat(1:jg));glat=0._r8
    allocate(gprs(1:kg));gprs=0._r8
    allocate(galt(1:kg));galt=0._r8
    allocate(gpri(1:lg));gpri=0._r8
    allocate(gali(1:lg));gali=0._r8
    allocate(gdat(1:jg,1:kg));gdat=0._r8

    read (01,'(3(1pe12.5))') o3vubc,prsubc,altubc
    read (01,'(a)') label   ! grid latitudes (deg) --------->
    read (01,'(10(1pe12.5))') (glat(j),j=1,jg)
    read (01,'(a)') label   ! layer pressure (mb) ---------->
    read (01,'(10(1pe12.5))') (gprs(k),k=1,kg)
    read (01,'(a)') label   ! layer altitude (km) ---------->
    read (01,'(10(1pe12.5))') (galt(k),k=1,kg)
    read (01,'(a)') label   ! interface pressure (mb) ------>
    read (01,'(10(1pe12.5))') (gpri(k),k=1,lg)
    read (01,'(a)') label   ! interface altitude (km) ------>
    read (01,'(10(1pe12.5))') (gali(k),k=1,lg)
    read (01,'(a)') label   ! o3 vmr (ppmv) lat-alt
    read (01,'(10(1pe12.5))')((gdat(j,k),j=1,jg),k=1,kg)
    close(01)

    lev_data_ape = gprs


    ! horizontal interp
    do i = 1, ncol
       if(lat_deg(i) .le. glat(1))then
           ozone_data_ape(1:kg,i)=gdat(1,1:kg)
       elseif(lat_deg(i) .ge. glat(jg))then
           ozone_data_ape(1:kg,i)=gdat(jg,1:kg)
       else
loop1:     do count = 1, jg-1
              if(lat_deg(i) .ge. glat(count) .and. lat_deg(i) .lt. glat(count+1))then
                  dy = glat(count+1)-glat(count)
                  w1 = (lat_deg(i)-glat(count))/dy
                  w2 = 1._r8-w1
                  ozone_data_ape(1:kg,i) = w1*gdat(count+1,1:kg)+w2*gdat(count,1:kg)
                  exit loop1
               endif
            enddo loop1
       endif
    enddo

    ozone_data_ape=ozone_data_ape*1.E-6

    deallocate(glat)
    deallocate(gprs)
    deallocate(galt)
    deallocate(gpri)
    deallocate(gali)
    deallocate(gdat)
 
 
    end subroutine ozone_ape

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

    real(r8) :: to_mmr(nlev,ncol)
    real(r8) :: molmass
    real(r8) :: amass
    character(len=32) :: units_str


    if( .not. has_prescribed_ozone ) return

    if(doAquaPlanet)then

        ! vertical interp
        allocate(log_lev_data(1:ozone_nlev))
        allocate(log_lev_model(1:nlev))
        allocate(ozone_ver(1:nlev,1:ncol))
 
        do i = 1, ncol
            log_lev_data  = log(lev_data_ape*100.)  ! hPa to Pa
            log_lev_model = log(pstate%pressure_at_pc_full_level%f(:,i)) 
            do k = 1, nlev
                if(log_lev_model(k) .ge. log_lev_data(ozone_nlev))then
                    ozone_ver(k,i) = ozone_data_ape(ozone_nlev,i)*(log_lev_model(k)/log_lev_data(ozone_nlev))**3.3
                else if(log_lev_model(k) .le. log_lev_data(1))then
                    ozone_ver(k,i) = ozone_data_ape(1,i)*(log_lev_model(k)/log_lev_data(1))**3.3
                else
loop2:              do count = 1, ozone_nlev-1
                        if(log_lev_data(count) .le. log_lev_model(k) .and. log_lev_data(count+1) .gt. log_lev_model(k))then
                            zzdh = log_lev_data(count+1)-log_lev_data(count)
                            zzd1 = log_lev_data(count+1)-log_lev_model(k)
                            zzd=zzd1/zzdh
                            ozone_ver(k,i) = ozone_data_ape(count+1,i)*(1-zzd)+ozone_data_ape(count,i)*zzd
                            exit loop2
                        end if
                    end do loop2
                end if
            end do
        end do
        
        pstate_cam%ghg_at_pc_full_level(o3_idx)%f(:,1:ncol) = ozone_ver(:,1:ncol)
        
        deallocate(log_lev_data)
        deallocate(log_lev_model)
        deallocate(ozone_ver)

    else
        
        molmass = 47.9981995_r8
        amass   = mwdry

        call advance_trcdata(nstep, dtime, ncol, lon, lat, &
                             fields, file)

        units_str = trim(to_lower(trim(fields(1)%units(:GLC(fields(1)%units)))))

        select case ( units_str )
        case ("molec/cm3","/cm3","molecules/cm3","cm^-3","cm**-3")
            to_mmr(:,:ncol) = (molmass*1.e6_r8*boltz*pstate%temp_at_pc_full_level%f(:,:ncol))/(amass*pstate%pressure_at_pc_full_level%f(:,:ncol))
        case ('kg/kg','mmr')
            to_mmr(:,:ncol) = 1._r8
        case ('mol/mol','mole/mole','vmr','fraction')
            to_mmr(:,:ncol) = molmass/amass
        case default
           if(mpi_rank()==0)print*, 'prescribed_ozone_adv: units = ',trim(fields(1)%units) ,' are not recognized'
           call endrun('prescribed_ozone_adv: units are not recognized')
        end select

        pstate_cam%ghg_at_pc_full_level(o3_idx)%f(:,1:ncol) = fields(1)%data(:,1:ncol)*to_mmr(:,:ncol)

    end if

    end subroutine prescribed_ozone_adv

 end module grist_prescribed_ozone
