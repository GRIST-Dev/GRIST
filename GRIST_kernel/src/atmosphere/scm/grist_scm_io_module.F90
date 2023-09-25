!======================================================
!
!  Created by LiXiaohan on 19/4/10.
!  Read forcing data for GRIST_SCM
!======================================================

 module grist_scm_io_module
    use grist_scm_comm_module
    use grist_mpi
    use grist_wrap_nf
    use grist_constants,                    only: i4, r8, latvap, cp
    use grist_nml_module,                   only: outdir,                       &
                                                  fname_output,                 &
                                                  comm_group_size,              &
                                                  nlev, nlevp, model_timestep
    use grist_handle_error,                 only: endrun
    use grist_hpe_constants,                only: eta_full_a, eta_full_b, p0
    use grist_scm_pressure_update,          only: time_integration_renew_mass_state
! major refactor for here
    use grist_scm_dyn_vars,                 only: n3, n3m1, scm_u, scm_v, scm_ps3, scm_t3, scm_q3
    use grist_domain_types,                 only: global_domain, group_comm, block_structure
    use grist_lib

    implicit none
    include 'netcdf.inc'

    private
    public  :: scm_output_final,            &
               grist_scm_io_read,           &
               grist_scm_io_output,         &
               scam_data_io_read


contains

    subroutine scm_output_final(mesh)
    ! io
    type(global_domain), intent(inout), target :: mesh
    ! local
    type(group_comm),      pointer             :: comm
    type(block_structure), pointer             :: local_block

    local_block => mesh%local_block
    comm => mesh%gcomm_write

    deallocate(comm%idx_localv)
    deallocate(comm%idx_localt)
    deallocate(comm%idx_locale)
    deallocate(comm%counts_sendv)
    deallocate(comm%counts_sendt)
    deallocate(comm%counts_sende)
    deallocate(comm%displs_sendv)
    deallocate(comm%displs_sendt)
    deallocate(comm%displs_sende)
    deallocate(comm%counts_recvv)
    deallocate(comm%counts_recvt)
    deallocate(comm%counts_recve)
    deallocate(comm%displs_recvv)
    deallocate(comm%displs_recvt)
    deallocate(comm%displs_recve)
    deallocate(comm%idx_recvv)
    deallocate(comm%idx_recvt)
    deallocate(comm%idx_recve)

    return
    end subroutine scm_output_final


    subroutine grist_scm_io_read(itimestep)

! io
    integer(i4), intent(in)           :: itimestep
! local
    integer(i4), parameter            :: omode = 0
    integer(i4)                       :: i, k
    integer(i4)                       :: ncid, dimid, varid, status
    integer(i4)                       :: plev, total_lev, bdate, ntime
    real(r8)                          :: var_srf
    integer(i4)                       :: current_ymd, current_tod, interval_s
    integer(i4)                       :: next_io_ymd, next_io_tod
    real(r8), allocatable             :: levs(:)
    integer(i4), allocatable          :: tsec(:)
    real(r8)                          :: tmp(1),tmpps(1,1,1)
    integer(i4)                       :: tmpint(1)
    real(r8), allocatable             :: pint(:), pdel(:), pmid(:), pdel_int(:)
    real(r8)                          :: weight
    logical                           :: fill_ends

    fill_ends= .false.

!================================================
!       get time information
!================================================
!    if(mpi_rank() .eq. 0)then
        call wrap_open(trim(scm_file_name), omode, ncid)
        status = nf_inq_varid (ncid, 'bdate', varid)
        if(status .ne. NF_NOERR)then
            status = nf_inq_varid (ncid, 'basedate', varid)
            if(status .ne. NF_NOERR)then
                print*,'can not get var ID for bdate in file',trim(scm_file_name)
                call endrun('Error in reading scm io file')
            end if
        end if
        call wrap_get_var_int(ncid, varid, tmpint)
        bdate = tmpint(1)

        status = nf_inq_dimid (ncid, 'time', dimid)
        if(status .ne. NF_NOERR)then
            status = nf_inq_dimid (ncid, 'tsec', dimid)
            if(status .ne. NF_NOERR)then
                print*,'can not get dim ID for time in file',trim(scm_file_name)
                call endrun('Error in reading scm io file')
            end if
        end if
        call wrap_inq_dimlen(ncid, dimid, ntime)

        allocate(tsec(ntime))
        call wrap_inq_varid(ncid, 'tsec', varid)
        call wrap_get_var_int(ncid, varid, tsec)

        if(itimestep .eq. 0)then
            io_time_index = 0
            do i = 1, ntime
                call get_current_time(bdate, 0, next_io_ymd, next_io_tod, tsec(i))
                if(start_ymd .gt. next_io_ymd .or.          &
                   (start_ymd .eq. next_io_ymd .and. start_tod .ge. next_io_tod))then
                   io_time_index = i
                endif
            end do

            if(io_time_index .eq. 0 .or. io_time_index .ge. ntime)then
                call get_current_time(bdate, 0, next_io_ymd, next_io_tod, tsec(ntime))
                print*, 'Error! Start model time does not fall within io_file period'
                print*, 'Start model date is ',start_ymd,' and ',start_tod,' seconds'
                print*, 'io_file start is ',bdate, ' and ',tsec(1),' seconds'
                print*, 'io_file end is   ',next_io_ymd, ' and ',next_io_tod,' seconds'
                call endrun("Error in reading scm_data_file")
            end if
            doioupdate = .true.
        else
            interval_s = itimestep*model_timestep
            call get_current_time(start_ymd, start_tod, current_ymd, current_tod, interval_s)
            call get_current_time(bdate, 0, next_io_ymd, next_io_tod, tsec(io_time_index+1))
            if(current_ymd .gt. next_io_ymd .or.            &
              (current_ymd .eq. next_io_ymd .and. current_tod .ge. next_io_tod))then
                io_time_index = io_time_index+1

                if(io_time_index .gt. ntime)then
                    print*, 'Reached the end of the time varient dataset'
                    doioupdate = .false.
                else
                    doioupdate = .true.
                end if
            else
                doioupdate = .false.
            end if
        end if

        deallocate(tsec)

        if(.not.doioupdate)then
            call wrap_close(ncid)
            return
        end if
!    end if

!================================================
!       read netDCF file containing
!       initial or forcing data
!================================================
        if(mpi_rank() .eq. 0) print*, "Read input data file:", trim(scm_file_name)
        call wrap_inq_dimid(ncid, 'lev', dimid)
        call wrap_inq_dimlen(ncid, dimid, plev)
! pressure levels:     
        allocate(levs(plev+1)); levs=0._r8
        call wrap_inq_varid(ncid, 'lev', varid)
        call wrap_get_var_realx(ncid, varid, levs(1:plev))
        call wrap_inq_varid(ncid, 'Ps', varid)
        call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
        psobs = tmp(1)
        levs(plev+1) = psobs
        total_lev = plev+1
        do k = plev, 1, -1
            if(levs(k) .gt. psobs)then
                total_lev = k
                levs(k) = psobs
            end if
        end do
        plev = total_lev
        if(plev .eq. 1) call endrun("Ps in scm_file is too low")
        if(psobs .lt. 1200.)then
            !turn the unit of ps and lev from hPa to Pa
            psobs = psobs*100.
            levs(:plev) = levs(:plev)*100.
        end if

! air temperature at the surface:
        status = nf_inq_varid (ncid, 'Tsair', varid)
        if(status .ne. NF_NOERR)then
            have_tsair = .false.
        else
            have_tsair = .true.
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            tsair = tmp(1)
        end if

! actual temperature:   LiXH test
!        if(itimestep .gt. 1)then
!            tobs(:) = scm_t3(:,n3)
!        end if

        if(have_tsair)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'T', levs(:plev),      &
                                            plev, psobs, tobs, status, fill_ends, tsair)
        else
            call get_vertical_interp_ncdata(ncid, io_time_index, 'T', levs(:plev),      &
                                            plev, psobs, tobs, status, fill_ends )
        end if

! ground temperature:
        status = nf_inq_varid (ncid, 'Tg', varid)
        if(status .ne. NF_NOERR)then
            have_tg = .false.
            if(have_tsair)then
!                print*, "Do not have Tg, using Tsair"
                tground = tsair
            else
!                print*, "Do not have Tg, using T at lowest level"
                tground = tobs(plev)
            end if
        else
            have_tg = .true.
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            tground = tmp(1)
        end if

! geopotential height:
! this value is used in the non-hydrostatic version of GRIST:
! if no geoobs exists, it will be diagnosed under hydrostatic balance.

!        status = nf_inq_varid (ncid, 'geopotential', varid)
!        if(status .ne. NF_NOERR)then
!            have_geo = .false.
!        else
!            have_geo = .true.
!            call get_vertical_interp_ncdata(ncid, io_time_index, 'geopotential', levs(:plev), &
!                                            plev, psobs, geoobs, status, fill_ends)
!        end if
 
! actual water vapor Mixing ratio: LiXH Test
!        if(itimestep .gt. 1)then
!            qobs(:) = scm_q3(1,:,n3)
!        end if

        status = nf_inq_varid (ncid, 'qsrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'q', levs(:plev),      &
                                            plev, psobs, qobs, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid,(/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'q', levs(:plev),      &
                                            plev, psobs, qobs, status, fill_ends, var_srf)
        end if

! cldliqobs has only one level in ARM95.nc data
!        call get_vertical_interp_ncdata(ncid, io_time_index, 'cldliq', levs(:plev),     &
!                                        plev, psobs, cldliqobs, status, .false.)
!        call get_vertical_interp_ncdata(ncid, io_time_index, 'cldice', levs(:plev),     &
!                                        plev, psobs, cldiceobs, status, .false.)

! observed cloud:
        status = nf_inq_varid (ncid, 'CLOUD', varid)
        if(status .ne. NF_NOERR)then
            status = nf_inq_varid (ncid, 'cld', varid)
            if(status .ne. NF_NOERR)then
!                print*, "Could not find variable CLOUD or cld in scm_data. Setting to zero"
                cldobs   = 0.
                have_cld = .false.
            else
                call get_vertical_interp_ncdata(ncid, io_time_index, 'cld', levs(:plev),&
                                        plev, psobs, cldobs, status, fill_ends)
                have_cld = .true.
            end if
        else
            call get_vertical_interp_ncdata(ncid, io_time_index, 'CLOUD', levs(:plev),      &
                                            plev, psobs, cldobs, status, fill_ends)
            have_cld = .true.
        end if

        status = nf_inq_varid (ncid, 'CONCLD', varid)
        if(status .ne. NF_NOERR)then
            concldobs   = 0.
            have_concld = .false.
        else
            call get_vertical_interp_ncdata(ncid, io_time_index, 'CONCLD', levs(:plev),      &
                                            plev, psobs, concldobs, status, fill_ends)
            have_concld = .true.
        end if

! observed clwp:        
        call get_vertical_interp_ncdata(ncid, io_time_index, 'clwp', levs(:plev),       &
                                        plev, psobs, clwpobs, status, fill_ends)
        if(status .ne. NF_NOERR)then
            have_clwp = .false.
        else
            have_clwp = .true.
        end if


! horizontal q advection:
        status = nf_inq_varid (ncid, 'divqsrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'divq', levs(:plev),   &
                                            plev, psobs, divq, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'divq', levs(:plev),   &
                                            plev, psobs, divq, status, fill_ends, var_srf)
        end if
        if(status .ne. NF_NOERR)then
            have_divq = .false.
        else
            have_divq = .true.
        end if

! read vertical advection if available:
        status = nf_inq_varid (ncid, 'vertdivqsrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'vertdivq', levs(:plev), &
                                            plev, psobs, vertdivq, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid,(/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'vertdivq', levs(:plev), &
                                            plev, psobs, vertdivq, status, fill_ends, var_srf)
        end if
        if(status .ne. NF_NOERR)then
            have_vertdivq = .false.
        else
            have_vertdivq = .true.
        end if
! read 3D advection if available:
! Attention! The variable name of divq3d may be changed=================>
        call get_vertical_interp_ncdata(ncid, io_time_index, 'q_dten', levs(:plev),      &
                                        plev, psobs, divq3d, status, fill_ends)
! Attention! The variable name of divq3d may be changed<=================
        if(status .ne. NF_NOERR)then
            if(have_divq .and. have_vertdivq)then
!                print*,"Do not have divq3d, using divq + vertdivq"
                divq3d = divq+vertdivq
            else
                divq3d = 0._r8
            end if
        end if


! horizontal temperature advection:
        status = nf_inq_varid (ncid, 'divTsrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'divT', levs(:plev),   &
                                            plev, psobs, divt, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid,(/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'divT', levs(:plev),   &
                                            plev, psobs, divt, status, fill_ends, var_srf)
        end if
        if(status .ne. NF_NOERR)then
            have_divt = .false.
        else
            have_divt = .true.
        end if
! read vertical advection if available:
        status = nf_inq_varid (ncid, 'vertdivTsrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'vertdivT', levs(:plev), &
                                            plev, psobs, vertdivt, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'vertdivT', levs(:plev), &
                                            plev, psobs, vertdivt, status, fill_ends, var_srf)
        end if
        if(status .ne. NF_NOERR)then
            have_vertdivt = .false.
        else
            have_vertdivt = .true.
        end if
! read 3D advection if available:
        status = nf_inq_varid (ncid, 'divT3dsrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'divT3d', levs(:plev), &
                                            plev, psobs, divt3d, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'divT3d', levs(:plev), &
                                            plev, psobs, divt3d, status, fill_ends, var_srf)
        end if
        if(status .ne. NF_NOERR)then
            if(have_divt .and. have_vertdivt)then
!                print*,"Do not have divT3d, using divT + vertdivT"
                divt3d = divt + vertdivt
            else
                divt3d = 0._r8
            end if
        end if

! make sure that use_3dfrc flag is set to true if we only have 3d forcing available:
        if(.not. have_divt .or. .not. have_divq)then
            use_3dfrc = .true.
        end if

! surface pressure tendency:
        status = nf_inq_varid (ncid, 'ptend', varid)
        if(status .ne. NF_NOERR)then
!            print*, "Could not find variable Ptend in scm_data. Setting to zero"
            ptend = 0.0_r8
        else
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            ptend = tmp(1)
        end if 

! omega:
        call get_vertical_interp_ncdata(ncid, io_time_index, 'omega', levs(:plev),      &
                                        plev, psobs, wfld, status, fill_ends, ptend)

        allocate(pint(nlevp));    pint=0._r8
        allocate(pdel(nlev));     pdel=0._r8
        allocate(pmid(nlev));     pmid=0._r8
        allocate(pdel_int(nlevp));pdel_int=0._r8
        tmp(1) = psobs

        call time_integration_renew_mass_state(1, tmp, pint, pdel, pmid, pdel_int)
        wfldint(1)     = 0._r8
        wfldint(nlevp) = 0._r8
        do k = 2, nlev
            weight = (pint(k) - pmid(k-1))/(pmid(k) - pmid(k-1))
            wfldint(k)= (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
        end do

        deallocate(pint)
        deallocate(pdel)
        deallocate(pmid)
        deallocate(pdel_int)

! u wind:
        status = nf_inq_varid (ncid, 'usrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'u', levs(:plev),      &
                                            plev, psobs, uobs, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'u', levs(:plev),      &
                                            plev, psobs, uobs, status, fill_ends, var_srf)
        end if
! v wind:
        status = nf_inq_varid (ncid, 'vsrf', varid)
        if(status .ne. NF_NOERR)then
            call get_vertical_interp_ncdata(ncid, io_time_index, 'v', levs(:plev),      &
                                            plev, psobs, vobs, status, fill_ends)
        else
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            var_srf = tmp(1)
            call get_vertical_interp_ncdata(ncid, io_time_index, 'v', levs(:plev),      &
                                            plev, psobs, vobs, status, fill_ends, var_srf)
        end if

! lhflx, shflx:
        status = nf_inq_varid (ncid, 'lhflx', varid)
        if(status .ne. NF_NOERR)then
            status = nf_inq_varid (ncid, 'lh', varid)
            if(status .ne. NF_NOERR)then
                have_lhflx = .false.
            else
                have_lhflx = .true.
                call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
                lhflxobs   = tmp(1)
            end if
        else
            have_lhflx = .true.
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            lhflxobs   = tmp(1)
        end if

        status = nf_inq_varid (ncid, 'shflx', varid)
        if(status .ne. NF_NOERR)then
            status = nf_inq_varid (ncid, 'sh', varid)
            if(status .ne. NF_NOERR)then
                have_shflx = .false.
            else
                have_shflx = .true.
                call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
                shflxobs   = tmp(1)
            end if
        else
            have_shflx = .true.
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/),tmp)
            shflxobs   = tmp(1)
        end if

! surface albedo:
        status = nf_inq_varid (ncid, 'alb', varid)
        if(status .ne. NF_NOERR)then
            have_alb = .false.
        else
            have_alb = .true.
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            albobs   = tmp(1)
        end if

! pertubation temperature and specific humidity:
        status = nf_inq_varid (ncid, 'tpert', varid)
        if(status .ne. NF_NOERR)then
            have_tpert = .false.
        else
            have_tpert = .true.
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            tpertobs   = tmp(1)
        end if

        status = nf_inq_varid (ncid, 'qpert', varid)
        if(status .ne. NF_NOERR)then
            have_qpert = .false.
        else
            have_qpert = .true.
            call wrap_get_vara_realx(ncid, varid, (/1,1,io_time_index/), (/1,1,1/), tmp)
            qpertobs   = tmp(1)
        end if

        deallocate(levs)
        call wrap_close(ncid)

    end subroutine grist_scm_io_read

    subroutine grist_scm_io_output(mesh, itimestep)
        use grist_data_types,                only: scalar_1d_field
        use grist_physics_data_structure,    only: pstate
#ifdef AMIPC_PHYSICS
        use grist_PhysC_data_structure,      only: pstate_cam
        use grist_physics_update,            only: old_time_level
#endif
#ifdef AMIPW_PHYSICS
        use grist_PhysW_data_structure,        only: pstate_wrf,p_qv,p_qc,p_qr
#endif
 
        use grist_time_manager,              only: get_current_date
        use grist_fileio_list_1d_module_par, only: wrap_output_init_1d , &
                                                   wrap_add_field_1d   , &
                                                   wrap_output_1d=>wrap_output_1d_group      , &
                                                   wrap_output_clean_1d
        use grist_scm_diagnose_module,       only: diag_physics_vars
! io
    integer(i4), intent(in)           :: itimestep
    type(global_domain), intent(in)   :: mesh
! local
    integer :: ierr
    real(r8)                          :: dt
    type(scalar_1d_field)             :: hyam, hybm
    type(scalar_1d_field)             :: uwind, vwind, temp, qv, qc, qi, &
                                         rh, pblh, tke, cmfmc
    type(scalar_1d_field)             :: precc, precl, prect, tmq
    type(scalar_1d_field)             :: cld,concld
    type(scalar_1d_field)             :: qrain,nrain
    type(scalar_1d_field)             :: fsns,flnt,lwcf,swcf
    type(scalar_1d_field)             :: temp_obs, qv_obs
    type(scalar_1d_field)             :: lwc, iwc, iclwc_sh, iclwc_dp
    type(scalar_1d_field)             :: dq_deep,dq_shallow,dq_micro,dq_macro,dq_pbl
    character(len=5)                  :: day
    character(len=5)                  :: sec

    dt = model_timestep
    call get_current_date(itimestep,dt,day,sec)

    allocate(hyam%f(nlev))
    allocate(hybm%f(nlev))
    allocate(uwind%f(nlev))
    allocate(vwind%f(nlev))
    allocate(temp%f(nlev))
    allocate(qv%f(nlev))
    allocate(qc%f(nlev))
    allocate(qi%f(nlev))
    allocate(temp_obs%f(nlev))
    allocate(qv_obs%f(nlev))
    allocate(pblh%f(nlev))
    allocate(tke%f(nlevp))
    allocate(cmfmc%f(nlevp))
    allocate(precc%f(nlev))
    allocate(precl%f(nlev))
    allocate(prect%f(nlev))
    allocate(tmq%f(nlev))
    allocate(cld%f(nlev))
    allocate(concld%f(nlev))
    allocate(rh%f(nlev))
    allocate(fsns%f(nlev))
    allocate(flnt%f(nlev))
    allocate(lwcf%f(nlev))
    allocate(swcf%f(nlev))

    allocate(qrain%f(nlev))
    allocate(nrain%f(nlev))
    allocate(iclwc_sh%f(nlev))
    allocate(iclwc_dp%f(nlev))
    allocate(lwc%f(nlev))
    allocate(iwc%f(nlev))

    allocate(dq_deep%f(nlev))
    allocate(dq_shallow%f(nlev))
    allocate(dq_micro%f(nlev))
    allocate(dq_macro%f(nlev))
    allocate(dq_pbl%f(nlev))

    hyam%pos        = 8
    hybm%pos        = 8
    uwind%pos       = 8
    vwind%pos       = 8
    temp%pos        = 8
    qv%pos          = 8
    qc%pos          = 8
    qi%pos          = 8
    temp_obs%pos    = 8
    qv_obs%pos      = 8
    precc%pos       = 8
    precl%pos       = 8
    prect%pos       = 8
    tmq%pos       = 8
    pblh%pos        = 8
    tke%pos         = 9
    cmfmc%pos       = 9
    cld%pos         = 8 
    concld%pos      = 8 
    rh%pos          = 8 
    fsns%pos        = 8 
    flnt%pos        = 8 
    lwcf%pos        = 8
    swcf%pos        = 8

    qrain%pos       = 8
    nrain%pos       = 8
    iclwc_sh%pos    = 8
    iclwc_dp%pos    = 8
    lwc%pos         = 8
    iwc%pos         = 8

    dq_deep%pos     = 8
    dq_shallow%pos     = 8
    dq_micro%pos     = 8
    dq_macro%pos     = 8
    dq_pbl%pos     = 8

    hyam%f(:)  = eta_full_a
    hybm%f(:)  = eta_full_b
    uwind%f(:) = scm_u            !scalar_U_wind_at_pc_full_level_n%f(:,1)
    vwind%f(:) = scm_v            !scalar_V_wind_at_pc_full_level_n%f(:,1)
    temp%f(:)  = scm_t3(:,n3m1)   !scalar_temp_at_pc_full_level_n%f(:,1)
    temp_obs%f(:) = tobs(:)
    qv_obs%f(:)   = qobs(:)
    precc%f(:) = spread(diag_physics_vars%precc%f(1), dim=1, ncopies=nlev)
    precl%f(:) = spread(diag_physics_vars%precl%f(1), dim=1, ncopies=nlev)
    prect%f(:) = spread(diag_physics_vars%prect%f(1), dim=1, ncopies=nlev)
    tmq%f(:)   = spread(diag_physics_vars%tmq%f(1),   dim=1, ncopies=nlev)
    lwcf%f(:)  = spread(diag_physics_vars%lwcf%f(1),  dim=1, ncopies=nlev)
    swcf%f(:)  = spread(diag_physics_vars%swcf%f(1),  dim=1, ncopies=nlev)
    cld%f(:)   = diag_physics_vars%cloud%f(:,1)
    qc%f(:)    = diag_physics_vars%cldliq%f(:,1) !scalar_tracer_mxrt_at_pc_full_level_n%f(2,:,1)
    qi%f(:)    = diag_physics_vars%cldice%f(:,1) !scalar_tracer_mxrt_at_pc_full_level_n%f(2,:,1)

    dq_deep%f(:)   = diag_physics_vars%dq_deep%f(:,1)
    dq_shallow%f(:)= diag_physics_vars%dq_shallow%f(:,1)
    dq_micro%f(:)  = diag_physics_vars%dq_micro%f(:,1)
    dq_macro%f(:)  = diag_physics_vars%dq_macro%f(:,1)
    dq_pbl%f(:)    = diag_physics_vars%dq_pbl%f(:,1)
 
#ifdef AMIPC_PHYSICS
    pblh%f(:)  = spread(pstate_cam%pbl_pblh_at_pc_surface%f(1), dim=1, ncopies=nlev)
    fsns%f(:)  = spread(pstate_cam%fswds_at_pc_surface%f(1), dim=1, ncopies=nlev)   &
                -spread(pstate_cam%fswus_at_pc_surface%f(1), dim=1, ncopies=nlev)
    flnt%f(:)  = spread(pstate_cam%flwut_at_pc_top%f(1), dim=1, ncopies=nlev)   &
                -spread(pstate_cam%flwdt_at_pc_top%f(1), dim=1, ncopies=nlev)
    tke%f(:)   = pstate_cam%pbl_tke_at_pc_face_level%f(:,1)
    cmfmc%f(:) = pstate_cam%updraft_mass_flux%f(:,1)
    concld%f(:)= pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,1)
    rh%f(:)    = pstate_cam%diag_relhum%f(:,1)

    qv%f(:)    = pstate%tracer_mxrt_at_pc_full_level%f(1,:,1) !scalar_tracer_mxrt_at_pc_full_level_n%f(1,:,1)
    qrain%f(:) = pstate_cam%microp_qrout_at_pc_full_level%f(:,1)
    nrain%f(:) = pstate_cam%microp_nrout_at_pc_full_level%f(:,1)
    lwc%f(:)   = pstate_cam%microp_lwc_at_pc_full_level%f(:,1)
    iwc%f(:)   = pstate_cam%microp_iwc_at_pc_full_level%f(:,1)
    iclwc_sh%f(:) = pstate_cam%sh_icwmr_at_pc_full_level%f(:,1)
    iclwc_dp%f(:) = pstate_cam%dp_icwmr_at_pc_full_level%f(:,1)
#endif

    fname_output = "GRIST.SCM."//trim(day)//"-"//trim(sec)//".1d.nc"

    call wrap_output_init_1d(mesh)
    call wrap_add_field_1d(hyam           , "hyam" )
    call wrap_add_field_1d(hybm           , "hybm" )
    call wrap_add_field_1d(temp_obs       , "Tobs" )
    call wrap_add_field_1d(qv_obs         , "qvobs")
    call wrap_add_field_1d(uwind          , "u"    )
    call wrap_add_field_1d(vwind          , "v"    )
    call wrap_add_field_1d(temp           , "T"    )
    call wrap_add_field_1d(precc          , "precc")
    call wrap_add_field_1d(precl          , "precl")
    call wrap_add_field_1d(prect          , "prect")
    call wrap_add_field_1d(tmq            , "tmq")
    call wrap_add_field_1d(cld            , "cloud") 
    call wrap_add_field_1d(qc             , "qc"   )
    call wrap_add_field_1d(qi             , "qi"   )
    call wrap_add_field_1d(lwcf           , "lwcf" )
    call wrap_add_field_1d(swcf           , "swcf" )
#ifdef AMIPC_PHYSICS
    call wrap_add_field_1d(pblh           , "pblh" )
    call wrap_add_field_1d(tke            , "tke"  )
    call wrap_add_field_1d(cmfmc          , "cmfmc")
    call wrap_add_field_1d(concld         , "concld") 
    call wrap_add_field_1d(rh             , "rh"   ) 
    call wrap_add_field_1d(fsns           , "swns" ) 
    call wrap_add_field_1d(flnt           , "lwnt" )
    call wrap_add_field_1d(qv             , "qv"   )
    call wrap_add_field_1d(qrain          , "qrain") 
    call wrap_add_field_1d(nrain          , "nrain") 
    call wrap_add_field_1d(lwc            , "lwc") 
    call wrap_add_field_1d(iwc            , "iwc") 
    call wrap_add_field_1d(iclwc_dp       , "iclwc_dp") 
    call wrap_add_field_1d(iclwc_sh       , "iclwc_sh") 
#endif

    call wrap_add_field_1d(dq_deep        , "dq_deep")
    call wrap_add_field_1d(dq_shallow     , "dq_shallow")
    call wrap_add_field_1d(dq_micro       , "dq_micro")
    call wrap_add_field_1d(dq_macro       , "dq_macro")
    call wrap_add_field_1d(dq_pbl         , "dq_pbl")


    call wrap_output_1d(mesh, outdir, fname_output)

    call wrap_output_clean_1d()

    deallocate(hyam%f, hybm%f, uwind%f, vwind%f, temp%f, qv%f, qc%f, qi%f)
    deallocate(tke%f, pblh%f, cmfmc%f, rh%f)
    deallocate(precc%f, precl%f, prect%f)
    deallocate(temp_obs%f, qv_obs%f)
    deallocate(cld%f, concld%f)
    deallocate(fsns%f, flnt%f, lwcf%f, swcf%f)
    deallocate(qrain%f, nrain%f, lwc%f, iwc%f, iclwc_sh%f, iclwc_dp%f)
    deallocate(dq_deep%f, dq_shallow%f, dq_micro%f, dq_macro%f, dq_pbl%f)

    end subroutine grist_scm_io_output

    subroutine get_vertical_interp_ncdata(ncid, timeidx, varname, data_levs,            &
                                          data_plev, ps, var, status, fillends,         &
                                          surfdat )
    implicit none
! io
    integer(i4), intent(in)           :: ncid
    integer(i4), intent(in)           :: timeidx
    character(*), intent(in)          :: varname
    integer(i4), intent(in)           :: data_plev
    real(r8), intent(in)              :: data_levs(data_plev)
    real(r8), intent(in)              :: ps
    real(r8), optional, intent(in)    :: surfdat
    logical, intent(in)               :: fillends
    integer(i4), intent(out)          :: status
    real(r8), intent(out)             :: var(:)
! local
    integer(i4)                       :: i
    integer(i4)                       :: status2
    integer(i4)                       :: plev, ndims, varid, attid, dimlen
    real(r8)                          :: missing_val
    character(256)                    :: dimname
    integer(i4), allocatable          :: var_dimids(:)
    integer(i4), allocatable          :: start(:), count(:)
    real(r8), allocatable             :: tmp(:), tmp1(:) 
    real(r8)                          :: dx, dy, m
    logical                           :: usable_dim


    status = nf_inq_varid (ncid, trim(varname), varid)
    if(status/=NF_NOERR)return

    call wrap_inq_varndims(ncid, varid, ndims)

    if(ndims .gt. 4)then
        call endrun('Error! extract data: The Var has more than 4 dimensions')
    end if
    plev = 0
    allocate(var_dimids(ndims))
    allocate(start(ndims), count(ndims))
    call wrap_inq_vardimid(ncid, varid, var_dimids)
    do i = ndims, 1, -1
        usable_dim = .false.
        call wrap_inq_dim(ncid, var_dimids(i), dimname, dimlen)
        if(trim(dimname) .eq. 'lon' .or. trim(dimname) .eq. 'longitude')then
            start(i) = 1
            count(i) = 1
            usable_dim = .true.
        end if
        if(trim(dimname) .eq. 'lat' .or. trim(dimname) .eq. 'latitude')then
            start(i) = 1
            count(i) = 1
            usable_dim = .true.
        end if
        if(trim(dimname) .eq. 'lev' .or. trim(dimname) .eq. 'ilev')then
            start(i) = 1
            count(i) = dimlen
            plev     = dimlen
            usable_dim = .true.
        end if
        if(trim(dimname) .eq. 'time' .or. trim(dimname) .eq. 'tsec')then
            start(i) = timeidx
            count(i) = 1
            usable_dim = .true.
        end if

        if(.not.(usable_dim))then
            print*,'The input var ', trim(varname),     &
                   ' has an unusable dimension ', trim(dimname)
            call endrun('Error! Check the var dimension')
        end if
    end do

    if(plev .eq. 0)then
        allocate(tmp1(1))
        call wrap_get_vara_realx(ncid, varid, start, count, tmp1)
        var(1) = tmp1(1)
        deallocate(tmp1)
        return
    else
        allocate(tmp1(plev+1))
        call wrap_get_vara_realx(ncid, varid, start, count, tmp1(:plev))
    end if


    deallocate(var_dimids)
    deallocate(start)
    deallocate(count)

    if(plev+1 .ge. data_plev)then
        allocate(tmp(data_plev))
        tmp(1:data_plev-1) = tmp1(1:data_plev-1)
    else
        call endrun('Error! var '//trim(varname)//' has levels less than plev')
    end if

    if(present(surfdat))then
        tmp(data_plev) = surfdat
    else
        dy = data_levs(data_plev-1) - data_levs(data_plev-2)
        dx = tmp(data_plev-1) - tmp(data_plev-2)
        if(dx .ne. 0._r8)then
            m =dy/dx
            tmp(data_plev) = ((data_levs(data_plev) - data_levs(data_plev-1))/m) + tmp(data_plev-1)
        else
            tmp(data_plev) = tmp(data_plev-1)
        end if
    end if


! check data for missing value
if(.false.)then
    status2 = nf_inq_attid(ncid, varid, 'missing_value', attid)
    if(status2 .eq. NF_NOERR)then
        call wrap_get_att_value(ncid, varid, 'missing_value', missing_val)
    else
        missing_val = -9999999.0_r8
    end if
    do i = 1, data_plev
        if(tmp(i) .eq. missing_val)then
            print*,'missing value',i,tmp(i),missing_val
            call endrun('Error! extract data: missing value found in Var')
        end if
    end do
end if

! vertical interpolation
    call vertical_interp(tmp, data_levs, data_plev, ps, var(:nlev), fillends)

! check negetive qv and correct to zero, LiXH
    if(trim(varname) .eq. 'q' .or. trim(varname) .eq. 'cldliq' .or. trim(varname) .eq. 'cldice')then
        where(var .lt. 0._r8)var = 0._r8
    end if

    deallocate(tmp)
    deallocate(tmp1)

    end subroutine get_vertical_interp_ncdata


    subroutine vertical_interp(inputdata, data_levs, dplev, ps, outdata, fillends)


    implicit none
! io
    integer(i4),  intent(in)          :: dplev
    real(r8), intent(in)              :: inputdata(dplev)
    real(r8), intent(in)              :: data_levs(dplev)
    real(r8), intent(in)              :: ps
    logical, intent(in)               :: fillends
    real(r8), intent(out)             :: outdata(nlev)
! local
    integer(i4)                       :: i, j
    integer(i4)                       :: mstart_lev, mend_lev
    integer(i4)                       :: im(nlev), ip(nlev)
    real(r8)                          :: wgtm(nlev), wgtp(nlev)
    real(r8)                          :: model_levs(nlev)

    do i = 1, nlev
        model_levs(i) = p0*eta_full_a(i) + ps*eta_full_b(i)
    end do

    ! check data pressure leves
    !if(data_levs(1) .gt. model_levs(1)) call endrun("Error! interp data: data levs less than model levs")
    !if(data_levs(dplev) .lt. model_levs(dplev)) call endrun("Error! interp data: data levs less than model levs")

    mstart_lev = 1
    do i = nlev, 1, -1
        if(model_levs(i) .gt. data_levs(1)) mstart_lev = i
    end do

    mend_lev = nlev
    do i = 1, nlev
        if(model_levs(i) .le. data_levs(dplev)) mend_lev = i
    end do

    do i = mstart_lev, mend_lev 
        do j = 1, dplev-1
            if(model_levs(i) .gt. data_levs(j) .and. model_levs(i) .le. data_levs(j+1))then
                im(i)   = j
                ip(i)   = j+1
                wgtm(i) = (data_levs(j+1)-model_levs(i))/(data_levs(j+1)-data_levs(j))
                wgtp(i) = (model_levs(i)-data_levs(j))/(data_levs(j+1)-data_levs(j))
                exit
            end if
        end do
        outdata(i) = inputdata(im(i))*wgtm(i) + inputdata(ip(i))*wgtp(i)
    end do

    if(fillends)then
        do i = 1, mstart_lev
            outdata(i) = inputdata(1)
        end do
        do i = mend_lev, nlev
            outdata(i) = inputdata(dplev)
        end do
    end if

    end subroutine vertical_interp

    subroutine get_current_time(start_date, start_sec, current_date, current_sec, inc_s)

! io
    integer(i4), intent(in)           :: start_date, start_sec
    integer(i4), intent(in)           :: inc_s
    integer(i4), intent(out)          :: current_date, current_sec
! local
    integer(i4)                       :: i
    integer(i4)                       :: year, mon, day, julday
    integer(i4)                       :: month(12),             &
                                         leap_year_month(12)
    logical                           :: leap_year

    month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
    leap_year_month = (/31,29,31,30,31,30,31,31,30,31,30,31/)

    year  = start_date/10000
    mon   = (start_date-year*10000)/100
    day   = start_date-year*10000-mon*100

    current_sec = start_sec+inc_s
    do while(current_sec .ge. 86400)
        day = day+1
        current_sec = current_sec-86400
    end do

    julday = day
    if((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or.     &
       (mod(year,400) .eq. 0))then
        do i = 1, mon-1
            julday = julday+leap_year_month(i)
        end do
    else
        do i = 1, mon-1
            julday = julday+month(i)
        end do
    end if

    do while(julday .gt. 365)
        if((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or.     &
           (mod(year,400) .eq. 0))then
            leap_year = .true.
        else
            leap_year = .false.
        end if

        if(leap_year .and. julday .gt. 366)then
            year = year+1
            julday = julday-366
        else if(leap_year .and. julday .eq. 366)then
            exit
        else
            year = year+1
            julday = julday-365
        end if
    end do
   
    mon = 1
    if((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or.     &
       (mod(year,400) .eq. 0))then
        do i = 1, 12
            if(julday .gt. leap_year_month(i))then
                julday = julday-leap_year_month(i)
                mon    = mon+1
            else
                day    = julday
                exit
            end if
        end do
    else
        do i = 1, 12
            if(julday .gt. month(i))then
                julday = julday-month(i)
                mon    = mon+1
            else
                day    = julday
                exit
            end if
        end do
    end if

    current_date = year*10000+mon*100+day

    end subroutine get_current_time

    subroutine scam_data_io_read(scam_file_name)
    !io
    character(*)            :: scam_file_name
    ! local
    integer(i4), parameter            :: omode = 0
    integer(i4)                       :: ncid, dimid, varid, status
    integer(i4)                       :: plev
    real(r8), allocatable             :: tmp(:)
 
    status = nf_open (trim(scam_file_name), omode, ncid)
    if (status/=NF_NOERR) then
      return
    end if 

!    call wrap_open(trim(scam_file_name), omode, ncid)
    call wrap_inq_dimid(ncid, 'lev', dimid)
    call wrap_inq_dimlen(ncid, dimid, plev)

    if(plev .ne. nlev)call endrun('Can not use scam data! nlev of data is not equal to grist')

    allocate(tmp(plev))
    call wrap_inq_varid(ncid, 'T', varid)
    call wrap_get_vara_realx(ncid, varid, (/1,1,1,1/), (/1,1,plev,1/), tmp)
    tobs(1:nlev) = tmp(1:plev)

    call wrap_inq_varid(ncid, 'Q', varid)
    call wrap_get_vara_realx(ncid, varid, (/1,1,1,1/), (/1,1,plev,1/), tmp)
    qobs(1:nlev) = tmp(1:plev)

    deallocate(tmp)
 
    end subroutine scam_data_io_read

 end module grist_scm_io_module
