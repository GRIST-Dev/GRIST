!======================================================
!
!  Created by LiXiaohan on 20/5/08.
!  Provides the interface to surface layer
!======================================================

 module grist_scm_surface_module
    use grist_constants,                    only: i4, r8, latvap, rad2deg, zero
    use grist_nml_module,                   only: nlev, nlevp, ntracer, physpkg
    use grist_scm_comm_module
!    use lixh_test_simple_physics,           only: large_scale_prcp,     &
!                                                  surface_flux
#ifdef AMIPC_PHYSICS
    use grist_amp_surface_flux_module ,     only: grist_amp_surface_flux_run
#endif
    use grist_time_manager,                 only: get_curr_date
    use grist_mpi

    implicit none
    private
    public  :: grist_scm_surface

 !  Private
    real(r8) :: albedo(24)

contains
#ifdef AMIPC_PHYSICS
    subroutine grist_scm_surface(istep, dtime, lat, lon)
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam 
! io
    real(r8)   , intent(in)            :: dtime
    integer(i4), intent(in)            :: istep
    real(r8)   , intent(in)            :: lat(1) 
    real(r8)   , intent(in)            :: lon(1)  

! local
    integer             :: yr, mn, dy, sc, hr
    character(len=256)  :: filename, filepath
    real(r8)            :: lat_loc, lon_loc

    if(pstate%landfrac_at_pc_surface%f(1) .eq. 1._r8)then
    ! Land Model
    ! This part should be modified to call the real land model, LiXH 
    !    call large_scale_prcp(1, nlev, dtime,      &
    !                          pstate%pressure_at_pc_surface%f(1),    &
    !                          pstate%pressure_at_pc_full_level%f(:,1),   &
    !                          pstate%pressure_at_pc_face_level%f(:,1),   &
    !                          pstate%temp_at_pc_full_level%f(:,1),         &
    !                          pstate%tracer_mxrt_at_pc_full_level%f(1,:,1),&
    !                          pstate%scalar_precl_surface%f(1))

    !    call surface_flux(1, nlev, pstate%pressure_at_pc_full_level%f(:,1),  &
    !                      pstate%pressure_at_pc_face_level%f(:,1),   &
    !                      pstate%u_wind_at_pc_full_level%f(:,1),   &
    !                      pstate%v_wind_at_pc_full_level%f(:,1),   &
    !                      pstate%tracer_mxrt_at_pc_full_level%f(1,:,1),&
    !                      pstate%temp_at_pc_full_level%f(:,1),     &
    !                      pstate%pressure_at_pc_surface%f(1),    &
    !                      pstate%atm_in_taux_at_pc_surface%f(1),      &
    !                      pstate%atm_in_tauy_at_pc_surface%f(1),      &
    !                      pstate%atm_in_shflx_at_pc_surface%f(1),     &
    !                      pstate%atm_in_qflx_at_pc_surface%f(1,1))

        if(trim(scm_test_name) .eq. 'arm95' .or. trim(scm_test_name) .eq. 'ARM95' .or.        &
           trim(scm_test_name) .eq. 'arm97' .or. trim(scm_test_name) .eq. 'ARM97')then

                if(istep .eq. 1)then
                    if(trim(scm_test_name) .eq. 'arm95' .or. trim(scm_test_name) .eq. 'ARM95')then
                        filename = 'adaptor.mars.for.ARM95.nc'
                    else
                        filename = 'adaptor.mars.for.ARM97.nc'
                    end if
                    filepath = '/g13/pengxd/lixh/grist/grist_v2/grist_merge/init_data/scm_data/'//trim(filename)

                    lat_loc = lat(1)
                    lon_loc = lon(1)
                    call get_albedo_from_data(filepath, lat_loc, lon_loc)
                end if

                call get_curr_date(start_ymd, start_tod, istep, dtime, yr, mn, dy, sc)

                hr = sc/3600+1
                if(hr .gt. 24)print*,'Error in albedo !!!'

                pstate%atm_in_asdir_at_pc_surface%f       = albedo(hr)
                pstate%atm_in_asdif_at_pc_surface%f       = albedo(hr)
                pstate%atm_in_aldir_at_pc_surface%f       = albedo(hr)
                pstate%atm_in_aldif_at_pc_surface%f       = albedo(hr)

        else
            if(mpi_rank()==0)print*,'No albedo data is fit for this scm case, set to be 0 or check code!'
            pstate%atm_in_asdir_at_pc_surface%f       = 0._r8
            pstate%atm_in_asdif_at_pc_surface%f       = 0._r8
            pstate%atm_in_aldir_at_pc_surface%f       = 0._r8
            pstate%atm_in_aldif_at_pc_surface%f       = 0._r8

        end if
 
    else
    ! Ocean model
        pstate%sst_at_pc_surface%f(1)            = pstate%ts_at_pc_surface%f(1)
        call grist_amp_surface_flux_run(1,nlev,istep,dtime,lat,lon)
    end if

    if(use_phys_vars .and. trim(physpkg) .eq. 'AMIPC_PHYSICS')then
        if(have_shflx) pstate%atm_in_shflx_at_pc_surface%f(1)    = shflxobs
        if(have_lhflx) pstate%atm_in_qflx_at_pc_surface%f(1,1)   = lhflxobs/latvap
        if(have_tpert .and. have_qpert)then
            pstate_cam%pbl_tpert_at_pc_surface%f(1)           = tpertobs
            pstate_cam%pbl_qpert_at_pc_surface%f(1)           = qpertobs
        end if
        if(have_alb)then
            pstate%atm_in_asdir_at_pc_surface%f       = albobs
            pstate%atm_in_asdif_at_pc_surface%f       = albobs
            pstate%atm_in_aldir_at_pc_surface%f       = albobs
            pstate%atm_in_aldif_at_pc_surface%f       = albobs
        end if

    end if

!    if(trim(scm_test_name) .eq. 'bomex' .or. trim(scm_test_name) .eq.'BOMEX')then
!        !use the LHF,SHF,SST,ps in the UNICON paper (Shin and Park, 2020), LiXH
!        pstate%atm_in_shflx_at_pc_surface%f(1)    = 9.5_r8
!        pstate%atm_in_qflx_at_pc_surface%f(1,1)   = 153._r8/latvap
!    end if


    end subroutine grist_scm_surface
#endif

#ifdef AMIPW_PHYSICS
    subroutine grist_scm_surface(istep, dtime, lat, lon, coszrs)
    use grist_PhysW_data_structure,         only: pstate_wrf, psurf_wrf
! io
    real(r8)   , intent(in)            :: dtime
    integer(i4), intent(in)            :: istep
    real(r8)   , intent(in)            :: lat(1) 
    real(r8)   , intent(in)            :: lon(1)  
    real(r8)   , intent(in)            :: coszrs(1)  

! local
    integer             :: yr, mn, dy, sc, hr, icell
    character(len=256)  :: filename, filepath
    real(r8)            :: lat_loc, lon_loc

    if(psurf_wrf%xland(1,1).ge.1.5_r8)then      ! ocean, include ice
        do icell = 1, 1
        if(coszrs(icell) .gt. zero)then ! daytime
! CAM3 code, default
        psurf_wrf%aldir(icell,1)  = (.026_r8/(coszrs(icell)**1.7_r8 + .065_r8)) + &
                                  (.15_r8*(coszrs(icell)- 0.10_r8)*(coszrs(icell) - 0.50_r8)*(coszrs(icell) - 1._r8))
        psurf_wrf%asdir(icell,1)  = psurf_wrf%aldir(icell,1)
        psurf_wrf%asdif(icell,1)  = 0.06_r8
        psurf_wrf%aldif(icell,1)  = 0.06_r8
        else
        psurf_wrf%asdir(icell,1)  = zero
        psurf_wrf%asdif(icell,1)  = zero
        psurf_wrf%aldir(icell,1)  = zero
        psurf_wrf%aldif(icell,1)  = zero
        end if    
        end do
    else
        if(trim(scm_test_name) .eq. 'arm95' .or. trim(scm_test_name) .eq. 'ARM95' .or.        &
           trim(scm_test_name) .eq. 'arm97' .or. trim(scm_test_name) .eq. 'ARM97')then

                if(istep .eq. 1)then
                    if(trim(scm_test_name) .eq. 'arm95' .or. trim(scm_test_name) .eq. 'ARM95')then
                        filename = 'adaptor.mars.for.ARM95.nc'
                    else
                        filename = 'adaptor.mars.for.ARM97.nc'
                    end if
                    filepath = '/g13/pengxd/lixh/grist/grist_v2/grist_merge/init_data/scm_data/'//trim(filename)

                    lat_loc = lat(1)
                    lon_loc = lon(1)
                    call get_albedo_from_data(filepath, lat_loc, lon_loc)
                end if

                call get_curr_date(start_ymd, start_tod, istep, dtime, yr, mn, dy, sc)

                hr = sc/3600+1
                if(hr .gt. 24)print*,'Error in albedo !!!'

                psurf_wrf%asdir       = albedo(hr)
                psurf_wrf%asdif       = albedo(hr)
                psurf_wrf%aldir       = albedo(hr)
                psurf_wrf%aldif       = albedo(hr)

        else
            if(mpi_rank()==0)print*,'No albedo data is fit for this scm case, set to be 0 or check code!'
            psurf_wrf%asdir       = 0._r8
            psurf_wrf%asdif       = 0._r8
            psurf_wrf%aldir       = 0._r8
            psurf_wrf%aldif       = 0._r8
        end if
    end if

    psurf_wrf%albedo = psurf_wrf%asdir

    if(use_phys_vars .and. trim(physpkg) .eq. 'AMIPW_PHYSICS')then
        if(have_shflx) psurf_wrf%hfx(1,1)    = shflxobs
        if(have_lhflx) psurf_wrf%qfx(1,1)   = lhflxobs/latvap
        if(have_alb)then
            psurf_wrf%asdir       = albobs
            psurf_wrf%asdif       = albobs
            psurf_wrf%aldir       = albobs
            psurf_wrf%aldif       = albobs
        end if
    end if

!    if(trim(scm_test_name) .eq. 'bomex' .or. trim(scm_test_name) .eq.'BOMEX')then
!        !use the LHF,SHF,SST,ps in the UNICON paper (Shin and Park, 2020), LiXH
!        pstate_wrf%hfx(1,1)    = 9.5_r8
!        pstate_wrf%qfx(1,1)   = 153._r8/latvap
!    end if

    end subroutine grist_scm_surface
#endif

    subroutine get_albedo_from_data(filepath, lat, lon)

    use grist_wrap_nf
    
    ! io
    real(r8),            intent(in)    :: lon, lat  !in radian
    character(len=256),  intent(in)    :: filepath

    ! local
    integer             :: i, j
    integer             :: fileid, dimid, timeid, latid, lonid, varid
    integer             :: ntime, nlat, nlon
    integer             :: lon_i(2), lat_j(2)
    real(r8)            :: grid_lon, grid_lat
    real(r8)            :: scale_facter, add_offset
    real(r8)            :: weight_lon(2), weight_lat(2)
    real(r8),   allocatable, dimension(:)     :: data_lon ,data_lat
    integer(2), allocatable, dimension(:,:,:) :: data_albedo_shrt
    real(r8),   allocatable, dimension(:,:,:) :: data_albedo

    grid_lon = lon*rad2deg
    grid_lat = lat*rad2deg
    
    if(grid_lon .ge. 180._r8 ) grid_lon = grid_lon - 360._r8

    call wrap_open(trim(filepath), 0, fileid)
    call wrap_inq_dimid(fileid,  'time',       dimid )
    call wrap_inq_dimlen(fileid, dimid,        ntime  )
    call wrap_inq_dimid(fileid,  'latitude',   dimid  )
    call wrap_inq_dimlen(fileid, dimid,        nlat   )
    call wrap_inq_dimid(fileid,  'longitude',  dimid  )
    call wrap_inq_dimlen(fileid, dimid,        nlon   )

    if(ntime .ne. 24) print*,'ERROR in albedo data! check the code!! LiXH'

    allocate(data_lon(nlon))
    allocate(data_lat(nlat))
    allocate(data_albedo_shrt(ntime, nlat, nlon))
    allocate(data_albedo(ntime, nlat, nlon))

    call wrap_inq_varid(    fileid,  'latitude',  latid)
    call wrap_get_var_realx(fileid,  latid,       data_lat )
    call wrap_inq_varid(    fileid,  'longitude', lonid)
    call wrap_get_var_realx(fileid,  lonid,       data_lon )

    do i = 1, nlon-1
        if(grid_lon .ge. data_lon(i) .and. grid_lon .lt. data_lon(i+1))then
            lon_i(1)      = i
            lon_i(2)      = i+1
            weight_lon(1) = (grid_lon - data_lon(i))/(data_lon(i+1)-data_lon(i))
            weight_lon(2) = 1._r8 - weight_lon(1)
            exit
        end if
    end do

    if(grid_lon .gt. data_lon(nlon))then
        lon_i(1)      = nlon
        lon_i(2)      = 1
        weight_lon(1) = (grid_lon - data_lon(nlon))/(data_lon(nlon)-data_lon(nlon-1))
        weight_lon(2) = 1._r8 - weight_lon(1)
    end if

    do j = 1, nlat-1
        if(grid_lat .ge. data_lat(j) .and.  grid_lat .lt. data_lat(j+1) )then
            lat_j(1)      = j
            lat_j(2)      = j+1
            weight_lat(1) = (grid_lat - data_lat(j))/(data_lat(j+1)-data_lat(j))
            weight_lat(2) = 1._r8 - weight_lat(1)
            exit
        end if
    end do

    if(grid_lat .ge. data_lat(nlat))then
        lat_j(1)      = nlat
        lat_j(2)      = nlat
        weight_lat(1) = 1._r8
        weight_lat(2) = 0._r8
    end if
    if(grid_lat .lt. data_lat(1))then
        lat_j(1)      = 1
        lat_j(2)      = 1
        weight_lat(1) = 1._r8
        weight_lat(2) = 0._r8
    end if

    call wrap_inq_varid(fileid, 'fal', varid)
    call wrap_get_var_short(fileid, varid, data_albedo_shrt)

    call wrap_get_att_value(fileid, varid, 'scale_factor', scale_facter)
    call wrap_get_att_value(fileid, varid, 'add_offset',   add_offset)

    data_albedo = real(data_albedo_shrt,r8)*scale_facter+add_offset

    do i = 1, 24
    albedo(i)   = (data_albedo(i,lat_j(1),lon_i(1))*weight_lon(2)+data_albedo(i,lat_j(1),lon_i(2))*weight_lon(1))*weight_lat(2) &
                 +(data_albedo(i,lat_j(2),lon_i(1))*weight_lon(2)+data_albedo(i,lat_j(2),lon_i(2))*weight_lon(1))*weight_lat(1)
    end do

    if(mpi_rank()==0)then
        print*,'LiXH Test albedo:',albedo
    end if

    end subroutine get_albedo_from_data


 end module grist_scm_surface_module
