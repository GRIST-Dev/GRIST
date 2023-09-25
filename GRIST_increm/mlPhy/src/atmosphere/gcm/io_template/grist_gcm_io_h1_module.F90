
 module grist_gcm_io_h1_module

  use grist_lib
  use grist_constants,        only: i4, i8
  use grist_hpe_constants,    only: eta_full_a, eta_full_b, eta_face_a, eta_face_b
  use grist_domain_types,     only: global_domain, group_comm, block_structure
  use grist_data_types,       only: scalar_1d_field, scalar_3d_field

  use grist_nml_module,       only: fname_output,          &
                                    working_mode,          &
                                    advection_scheme,      &
                                    conserve_scheme,       &
                                    model_timestep,        &
                                    testcase,              &
                                    outdir, nlev, ntracer, &
                                    comm_group_size,       &
                                    start_ymd, start_tod, grid_info

  use grist_time_manager,     only: get_current_date, get_curr_cdate

  use grist_fileio_list_3d_module_par, only: wrap_output_init_3d     ,&
                                             wrap_add_field_3d       ,&
                                             wrap_output_3d_group    ,&
                                             wrap_output_3d_group_sp ,&
                                             wrap_output_clean_3d

  use grist_fileio_list_2d_module_par, only: wrap_output_init_2d      ,&
                                             wrap_add_field_2d        ,&
                                             wrap_output_2d_group     ,&
                                             wrap_output_2d_group_sp  ,&
                                             wrap_output_clean_2d

  use grist_fileio_list_1d_module_par, only: wrap_output_init_1d    ,&
                                             wrap_add_field_1d      ,&
                                             wrap_output_1d_group   ,&
                                             wrap_output_1d_group_sp,&
                                             wrap_output_clean_1d   ,&
                                             wrap_read_1d_group_rst
! data module
  use grist_dycore_vars_module, only: dycoreVarGeography, dycoreVarCellFull, dycoreVarCellFace, dycoreVarVertFull,dycoreVarSurface
  use grist_tracer_transport_vars_module, only: tracerVarCellFull, tracerVarCellFace, tracerVarEdgeFull
  use grist_dtp_vars_module
  use grist_physics_data_structure, only: pstate
  use grist_datam_static_data_module,  only: staticData_phis_at_pc_surface, &
                                             staticData_sst_at_pc_surface
  use grist_datam_initial_data_module, only: initialData_uuu_at_pc_full_level, &
                                             initialData_vvv_at_pc_full_level, &
                                             initialData_ttt_at_pc_full_level, &
                                             initialData_ps_at_pc_surface
#ifdef AMIPW_PHYSICS
  use grist_wrf_data_structure,    only: pstate_wrf, psurf_wrf,  ptend_wrf
#endif
#ifdef AMIPC_PHYSICS
  use grist_cam5_data_structure,   only: pstate_cam
  use grist_physics_update,        only: old_time_level
#endif
!
! gcm-diagnose
!
  use grist_gcm_diagnose_h1_module,only: gcm_h1_1d_inst_physics_variables, &
                                         gcm_h1_2d_inst_physics_variables, &
                                         gcm_h1_1d_accu_physics_variables, &
                                         gcm_h1_1d_dump_physics_variables, &
                                         gcm_h1_1d_rest_physics_variables, &
                                         gcm_h1_2d_accu_physics_variables, &
                                         gcm_h1_2d_dump_physics_variables, &
                                         gcm_h1_2d_rest_physics_variables, &
                                         diag_phys_accu_vars_h1_1d,        &
                                         diag_phys_accu_vars_h1_2d,        &
                                         diag_phys_inst_vars_h1_1d,        &
                                         diag_phys_inst_vars_h1_2d

  use grist_dycore_diagnose_module_2d,  only: dycore_diagnose_variables

  implicit none

  private

  public  :: gcm_output_init         , &
             gcm_output_final        , &
             gcm_output_history_h1   , &
             gcm_output_history_h1_1d, &
             gcm_output_history_h1_2d, &
             gcm_output_history_h1_3d, &
             gcm_output_atm_h1_file

  contains
!================================================
! [1] Construct filename according to itimestep
! [2] Calculate Diag Vars using Prog vars
! [3] Do an Output Flow
!================================================

  subroutine gcm_output_init(mesh)
    ! io
    type(global_domain),   intent(inout), target :: mesh
    ! local
    integer(i4)                                  :: i, ig, np_output
    integer(i4)                                  :: nv_avg, nv_res, vpos, &
                                                    ne_avg, ne_res, epos, &
                                                    nt_avg, nt_res, tpos
    type(set)                                    :: sv_list, se_list, st_list
    integer,               allocatable           :: buff_svi(:), buff_sti(:), buff_sei(:)

    type(block_structure), pointer               :: local_block
    type(group_comm),      pointer               :: comm
    integer                                      :: nprocs, ierr
    integer, allocatable                         :: recvvi(:), recvti(:), recvei(:)
    integer(i4)                                  :: rvtot, rttot, retot

    local_block => mesh%local_block
    comm => mesh%gcomm_write
    comm%all_comm = local_block%comm
    comm%gsize = comm_group_size
    comm%color = mod(mpi_rank(comm%all_comm),comm%gsize)
    call mpi_comm_split(comm%all_comm, comm%color, mpi_rank(comm%all_comm), comm%gcomm,ierr)
    call mpi_comm_rank(comm%gcomm,comm%grank,ierr)

    nprocs = mpi_size(comm%all_comm)
    np_output = ceiling(nprocs/real(comm_group_size, 8))
    ! for vtx
    ! sort v_index & get localv
    do i = 1, size(local_block%output%v)
        call sv_list%insert(local_block%output%v(i))
    enddo
    allocate(buff_svi(sv_list%size()))
    call sv_list%sort(buff_svi)
    call sv_list%final
    allocate(comm%idx_localv(size(buff_svi)))
    do i = 1, size(buff_svi)
       comm%idx_localv(i) = mesh%map_g2l_v%find(buff_svi(i))
    end do
    ! sendv counts & displs
    allocate(comm%counts_sendv(1:nprocs), source=0)
    allocate(comm%displs_sendv(1:nprocs), source=0)
    nv_avg = mesh%nv_all/np_output
    nv_res = mod(mesh%nv_all, np_output)
    vpos = (nv_avg+1)*nv_res
    do i = 1, size(buff_svi)
        if(buff_svi(i) .gt. vpos) then
            ig = (buff_svi(i)-vpos-1)/nv_avg + nv_res
        else
            ig = (buff_svi(i)-1)/(nv_avg+1)
        end if
        comm%counts_sendv(ig*comm_group_size+1) = comm%counts_sendv(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sendv(i) = comm%displs_sendv(i-1) + comm%counts_sendv(i-1)
    end do
    ! recvv counts & displs
    allocate(comm%counts_recvv(nprocs), source=0)
    allocate(comm%displs_recvv(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sendv, 1, MPI_INTEGER, comm%counts_recvv, 1, MPI_INTEGER, comm%all_comm, ierr)
    if(comm%color .eq. 0) then
        rvtot = comm%counts_recvv(1)
        do i = 2, nprocs
            comm%displs_recvv(i) = comm%displs_recvv(i-1) + comm%counts_recvv(i-1)
            rvtot = rvtot + comm%counts_recvv(i)
        end do
    else
        rvtot = 1
    end if
    allocate(recvvi(rvtot))
    call MPI_Alltoallv(buff_svi, comm%counts_sendv, comm%displs_sendv, MPI_INTEGER, recvvi, &
                       comm%counts_recvv, comm%displs_recvv, MPI_INTEGER, comm%all_comm, ierr)
    !open(101, file="recvv"//i2s(mpi_rank()))
    !do i = 1, size(recvvi)
    !  write(101, *) recvvi(i)
    !end do
    !close(101)
    ! recvv local idx
    allocate(comm%idx_recvv(rvtot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. nv_res) then
            comm%mystart_v = comm%grank*(nv_avg+1) + 1
        else
            comm%mystart_v = vpos + (comm%grank-nv_res)*nv_avg + 1
        end if
        if(comm%grank .lt. nv_res) then
            comm%mycount_v = nv_avg + 1
        else
            comm%mycount_v = nv_avg
        end if
        do i = 1, rvtot
            comm%idx_recvv(i) = recvvi(i)-comm%mystart_v+1
        end do
    end if
    deallocate(buff_svi)
    deallocate(recvvi)

    ! for edg
    ! sort e_index & get locale
    do i = 1, size(local_block%output%e)
        call se_list%insert(local_block%output%e(i))
    enddo
    allocate(buff_sei(se_list%size()))
    call se_list%sort(buff_sei)
    call se_list%final
    allocate(comm%idx_locale(size(buff_sei)))
    do i = 1, size(buff_sei)
       comm%idx_locale(i) = mesh%map_g2l_e%find(buff_sei(i))
    end do
    ! sende counts & displs
    allocate(comm%counts_sende(1:nprocs), source=0)
    allocate(comm%displs_sende(1:nprocs), source=0)
    ne_avg = mesh%ne_all/np_output
    ne_res = mod(mesh%ne_all, np_output)
    epos = (ne_avg+1)*ne_res
    do i = 1, size(buff_sei)
        if(buff_sei(i) .gt. epos) then
            ig = (buff_sei(i)-epos-1)/ne_avg + ne_res
        else
            ig = (buff_sei(i)-1)/(ne_avg+1)
        end if
        comm%counts_sende(ig*comm_group_size+1) = comm%counts_sende(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sende(i) = comm%displs_sende(i-1) + comm%counts_sende(i-1)
    end do
    ! recve counts & displs
    allocate(comm%counts_recve(nprocs), source=0)
    allocate(comm%displs_recve(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sende, 1, MPI_INTEGER, comm%counts_recve, 1, MPI_INTEGER, comm%all_comm, ierr)
    if(comm%color .eq. 0) then
        retot = comm%counts_recve(1)
        do i = 2, nprocs
            comm%displs_recve(i) = comm%displs_recve(i-1) + comm%counts_recve(i-1)
            retot = retot + comm%counts_recve(i)
        end do
    else
        retot = 1
    end if
    allocate(recvei(retot))
    call MPI_Alltoallv(buff_sei, comm%counts_sende, comm%displs_sende, MPI_INTEGER, recvei, &
                       comm%counts_recve, comm%displs_recve, MPI_INTEGER, comm%all_comm, ierr)
    ! recve local idx
    allocate(comm%idx_recve(retot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. ne_res) then
            comm%mystart_e = comm%grank*(ne_avg+1) + 1
        else
            comm%mystart_e = epos + (comm%grank-ne_res)*ne_avg + 1
        end if
        if(comm%grank .lt. ne_res) then
            comm%mycount_e = ne_avg + 1
        else
            comm%mycount_e = ne_avg
        end if
        do i = 1, retot
            comm%idx_recve(i) = recvei(i)-comm%mystart_e+1
        end do
    end if
    deallocate(buff_sei)
    deallocate(recvei)

    ! for tri
    ! sort t_index & get localt
    do i = 1, size(local_block%output%t)
        call st_list%insert(local_block%output%t(i))
    enddo
    allocate(buff_sti(st_list%size()))
    call st_list%sort(buff_sti)
    call st_list%final
    allocate(comm%idx_localt(size(buff_sti)))
    do i = 1, size(buff_sti)
       comm%idx_localt(i) = mesh%map_g2l_t%find(buff_sti(i))
    end do
    ! sendt counts & displs
    allocate(comm%counts_sendt(1:nprocs), source=0)
    allocate(comm%displs_sendt(1:nprocs), source=0)
    nt_avg = mesh%nt_all/np_output
    nt_res = mod(mesh%nt_all, np_output)
    tpos = (nt_avg+1)*nt_res
    do i = 1, size(buff_sti)
        if(buff_sti(i) .gt. tpos) then
            ig = (buff_sti(i)-tpos-1)/nt_avg + nt_res
        else
            ig = (buff_sti(i)-1)/(nt_avg+1)
        end if
        comm%counts_sendt(ig*comm_group_size+1) = comm%counts_sendt(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sendt(i) = comm%displs_sendt(i-1) + comm%counts_sendt(i-1)
    end do
    ! recvt counts & displs
    allocate(comm%counts_recvt(nprocs), source=0)
    allocate(comm%displs_recvt(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sendt, 1, MPI_INTEGER, comm%counts_recvt, 1, MPI_INTEGER, comm%all_comm, ierr)
    if(comm%color .eq. 0) then
        rttot = comm%counts_recvt(1)
        do i = 2, nprocs
            comm%displs_recvt(i) = comm%displs_recvt(i-1) + comm%counts_recvt(i-1)
            rttot = rttot + comm%counts_recvt(i)
        end do
    else
        rttot = 1
    end if
    allocate(recvti(rttot))
    call MPI_Alltoallv(buff_sti, comm%counts_sendt, comm%displs_sendt, MPI_INTEGER, recvti, &
                       comm%counts_recvt, comm%displs_recvt, MPI_INTEGER, comm%all_comm, ierr)
    ! recvt local idx
    allocate(comm%idx_recvt(rttot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. nt_res) then
            comm%mystart_t = comm%grank*(nt_avg+1) + 1
        else
            comm%mystart_t = tpos + (comm%grank-nt_res)*nt_avg + 1
        end if
        if(comm%grank .lt. nt_res) then
            comm%mycount_t = nt_avg + 1
        else
            comm%mycount_t = nt_avg
        end if
        do i = 1, rttot
            comm%idx_recvt(i) = recvti(i)-comm%mystart_t+1
        end do
    end if
    deallocate(buff_sti)
    deallocate(recvti)

    return
  end subroutine gcm_output_init

  subroutine gcm_output_final(mesh)
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
  end subroutine gcm_output_final

  subroutine gcm_output_history_h1(mesh,itimestep,nsteps,history_freq)
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps 
   integer(i4)        , intent(in)     :: history_freq 

    call gcm_h1_1d_accu_physics_variables
    call gcm_h1_2d_accu_physics_variables

    if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,history_freq)==0))then
        ! only for those vars not in Diag-h1, otherwise we should put this line before h1_accu as in h0-diag;
        ! currently we do not need to accumulate variables from dycore_diagnoise to h1's Diag vars
        call dycore_diagnose_variables(mesh)
        call gcm_h1_1d_dump_physics_variables
        call gcm_h1_2d_dump_physics_variables
        call gcm_h1_1d_inst_physics_variables ! time-instant vars
        call gcm_h1_2d_inst_physics_variables ! time-instant vars
        call gcm_output_atm_h1_file(mesh, itimestep)
        call gcm_h1_1d_rest_physics_variables
        call gcm_h1_2d_rest_physics_variables
    end if
    return
  end subroutine gcm_output_history_h1

  subroutine gcm_output_history_h1_1d(mesh,itimestep,nsteps,history_freq)
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps 
   integer(i4)        , intent(in)     :: history_freq 

    call gcm_h1_1d_accu_physics_variables
    if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,history_freq)==0))then
        ! only for those vars not in Diag-h1, otherwise we should put this line before h1_accu as in h0-diag;
        ! currently we do not need to accumulate variables from dycore_diagnoise to h1's Diag vars
        call dycore_diagnose_variables(mesh)
        call gcm_h1_1d_dump_physics_variables ! obtain time-averaged vars
        call gcm_h1_1d_inst_physics_variables ! diag   time-instant vars
        call gcm_output_atm_h1_1d_file(mesh, itimestep)
        call gcm_h1_1d_rest_physics_variables
    end if
    return
  end subroutine gcm_output_history_h1_1d

  subroutine gcm_output_history_h1_2d(mesh,itimestep,nsteps,history_freq)
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps 
   integer(i4)        , intent(in)     :: history_freq 

    call gcm_h1_2d_accu_physics_variables
    if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,history_freq)==0))then
        ! only for those vars not in Diag-h1, otherwise we should put this line before h1_accu as in h0-diag;
        ! currently we do not need to accumulate variables from dycore_diagnoise to h1's Diag vars
        call dycore_diagnose_variables(mesh)
        call gcm_h1_2d_dump_physics_variables
        call gcm_h1_2d_inst_physics_variables ! diag   time-instant vars
        call gcm_output_atm_h1_2d_file(mesh, itimestep)
        call gcm_h1_2d_rest_physics_variables
    end if
    return
  end subroutine gcm_output_history_h1_2d

  subroutine gcm_output_history_h1_3d(mesh,itimestep,nsteps,history_freq)
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
   real(r8)           , intent(in)     :: nsteps 
   integer(i4)        , intent(in)     :: history_freq 

    if((itimestep == 1 .or. itimestep==nsteps .or. mod(itimestep,history_freq)==0))then
        ! only for those vars not in Diag-h1, otherwise we should put this line before h1_accu as in h0-diag;
        ! currently we do not need to accumulate variables from dycore_diagnoise to h1's Diag vars
        call gcm_output_atm_h1_3d_file(mesh, itimestep)
    end if
    return
  end subroutine gcm_output_history_h1_3d

  subroutine gcm_output_atm_h1_1d_file(mesh, itimestep)
    use grist_util_module, only: write_string
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
! local
   character(128)                      :: c_glevel
   character(128)                      :: c_testcase
   character(len=4)                    :: cyear
   character(len=2)                    :: cmon
   character(len=2)                    :: cday
   character(len=5)                    :: csec
   character(len=5)                    :: day
   character(len=5)                    :: sec

!================================================
! Info for constructing filename
!================================================

   call write_string(mesh%glevel     , c_glevel)
   call write_string(testcase        , c_testcase)
! continous date from 00000-00000
   call get_current_date(itimestep,model_timestep,day,sec)
! calendar date
   call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, &
                       cyear, cmon, cday, csec)
!
! ATM 1d file
!
    call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"1d.h1",fname_output)
    call write_atm_file1d(mesh,fname_output)
!
! LND 1d file
!
    call obtain_lnd_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"1d.h1",fname_output)
    call write_lnd_file1d(mesh,fname_output)

    return
  end subroutine gcm_output_atm_h1_1d_file

 subroutine gcm_output_atm_h1_2d_file(mesh, itimestep)
    use grist_util_module, only: write_string
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
! local
   character(128)                      :: c_glevel
   character(128)                      :: c_testcase
   character(len=4)                    :: cyear
   character(len=2)                    :: cmon
   character(len=2)                    :: cday
   character(len=5)                    :: csec
   character(len=5)                    :: day
   character(len=5)                    :: sec

!================================================
! Info for constructing filename
!================================================

   call write_string(mesh%glevel     , c_glevel)
   call write_string(testcase        , c_testcase)
! continous date from 00000-00000
   call get_current_date(itimestep,model_timestep,day,sec)
! calendar date
   call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, &
                       cyear, cmon, cday, csec)

!
! ATM 2d file
!
    call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"2d.h1",fname_output)
    call write_atm_file2d(mesh,fname_output)

    return
  end subroutine gcm_output_atm_h1_2d_file

  subroutine gcm_output_atm_h1_3d_file(mesh, itimestep)
    use grist_util_module, only: write_string
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
! local
   character(128)                      :: c_glevel
   character(128)                      :: c_testcase
   character(len=4)                    :: cyear
   character(len=2)                    :: cmon
   character(len=2)                    :: cday
   character(len=5)                    :: csec
   character(len=5)                    :: day
   character(len=5)                    :: sec

!================================================
! Info for constructing filename
!================================================

   call write_string(mesh%glevel     , c_glevel)
   call write_string(testcase        , c_testcase)
! continous date from 00000-00000
   call get_current_date(itimestep,model_timestep,day,sec)
! calendar date
   call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, &
                       cyear, cmon, cday, csec)
!
! ATM 3d file, not output when using dycore mode
!
    IF(.not.trim(working_mode).eq.'dycore')then
       call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme, &
                                    cyear,cmon,cday,csec,day,sec,"3d.h1",fname_output)
       call write_atm_file3d(mesh,fname_output)
    Endif

    return
  end subroutine gcm_output_atm_h1_3d_file

  subroutine gcm_output_atm_h1_file(mesh, itimestep)
    use grist_util_module, only: write_string
! io
   type(global_domain), intent(in)     :: mesh
   integer(i4)        , intent(in)     :: itimestep
! local
   character(128)                      :: c_glevel
   character(128)                      :: c_testcase
   character(len=4)                    :: cyear
   character(len=2)                    :: cmon
   character(len=2)                    :: cday
   character(len=5)                    :: csec
   character(len=5)                    :: day
   character(len=5)                    :: sec

!================================================
! Info for constructing filename
!================================================

   call write_string(mesh%glevel     , c_glevel)
   call write_string(testcase        , c_testcase)
! continous date from 00000-00000
   call get_current_date(itimestep,model_timestep,day,sec)
! calendar date
   call get_curr_cdate(start_ymd, start_tod, itimestep, model_timestep, &
                       cyear, cmon, cday, csec)
!
! ATM 3d file, not output when using dycore mode
!
    IF(.not.trim(working_mode).eq.'dycore')then
       call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme, &
                                    cyear,cmon,cday,csec,day,sec,"3d.h1",fname_output)
       call write_atm_file3d(mesh,fname_output)
    Endif
!
! ATM 2d file
!
    call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"2d.h1",fname_output)
    call write_atm_file2d(mesh,fname_output)
!
! ATM 1d file
!
    call obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"1d.h1",fname_output)
    call write_atm_file1d(mesh,fname_output)
!
! LND 1d file
!
    call obtain_lnd_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,"1d.h1",fname_output)
    call write_lnd_file1d(mesh,fname_output)

    return
  end subroutine gcm_output_atm_h1_file

  subroutine write_atm_file3d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output
! local
   type(scalar_3d_field)               :: tracer_mxrt
   integer(i4)                         :: ilev, itracer

#ifdef AMIPW_PHYSICS
!    if(.not.allocated(tracer_mxrt%f))        allocate(tracer_mxrt%f(ntracer,nlev,mesh%nv_full))
!    tracer_mxrt%pos  = scalar_tracer_mass_at_pc_full_level_n%pos
!    do ilev = 1, nlev
!       do itracer = 1, ntracer
!          tracer_mxrt%f(itracer,ilev,1:mesh%nv_halo(1)) = (pstate_wrf%cldfra (1:mesh%nv_halo(1), nlev+1-ilev, 1))
!       end do
!    end do
#endif

    call wrap_output_init_3d(mesh)
    call wrap_add_field_3d(tracerVarCellFull%scalar_tracer_mxrt_n,"tracerMxrt","mixing ratio of tracer species","kg/kg")
!    call wrap_add_field_3d(scalar_tracer_mass_at_pc_full_level_n,"tracerMass")

#ifdef AMIPW_PHYSICS
!    call wrap_add_field_3d(tracer_mxrt,"cldfra") ! this is internal state of physics
#endif

#ifdef SPIO
    call wrap_output_3d_group_sp(mesh, outdir,fname_output)
#else
    call wrap_output_3d_group(   mesh, outdir,fname_output)
#endif
    call wrap_output_clean_3d()

    return
  end subroutine write_atm_file3d

  subroutine write_atm_file2d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output
! local
   type(scalar_2d_field)               :: data2d
   integer(i4)                         :: ilev

! start output
    call wrap_output_init_2d(mesh)

    IF(.not.trim(working_mode).eq.'tracer')then
!
! dycore vars
!
    call wrap_add_field_2d(dycoreVarCellFull%scalar_U_wind_n           ,"uPC"   ,"zonal wind speed","m/s"   )
    call wrap_add_field_2d(dycoreVarCellFull%scalar_V_wind_n           ,"vPC"   ,"meridional wind speed","m/s"   )
    call wrap_add_field_2d(dycoreVarCellFace%scalar_www_n              ,"wwwFace","height-based vertical speed (NDC-only)", "m/s"  )
    call wrap_add_field_2d(dycoreVarCellFull%scalar_temp_n             ,"temperature","Air temperature","K")
!    call wrap_add_field_2d(dycoreVarVertFull%scalar_rel_vor_n          ,"rVor", "relative vorticity of flow", "s^-1" )
!    call wrap_add_field_2d(dycoreVarCellFull%scalar_divergence_n       ,"div" , "divergence of flow" , "s^-1"  )
    call wrap_add_field_2d(dycoreVarCellFace%scalar_geopotential_n     ,"geopotentialFace", "geopotential at face level", "m^2/s^2")
    call wrap_add_field_2d(dycoreVarCellFace%scalar_pressure_n         ,"pressureFace"  ,"full-air pressure at face level","Pa" )
    call wrap_add_field_2d(dycoreVarCellFace%scalar_mpressure_n        ,"mpressureFace" ,"moist-mass pressure at face level", "Pa")
    call wrap_add_field_2d(dycoreVarCellFull%scalar_omega_n            ,"omega" ,"(m)pressure-based vertical speed","Pa/s") ! dmp/dt
#ifdef AMIPW_PHYSICS
    call wrap_add_field_2d(diag_phys_inst_vars_h1_2d%cloud            ,"cloud3D","3d cloud fraction","")
    call wrap_add_field_2d(diag_phys_inst_vars_h1_2d%refl_10cm        ,"refl_10cm","diagnosed reflectivity from microphysics routine","")
    !call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rad_thten        ,"rad_thtenDiag","radiative heating for potential temperature:A","K/s")
#ifdef TEND_IO
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rthraten            ,"rthraten","rth tendency due to RA","") ! dmp/dt
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rthcuten           ,"rthcuten","rth tendency due to CU",""   )
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rthmpten       ,"rthmpten" ,"qv tendency due to MP" , ""  )
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rthblten ,"rthblten", "rth tendency due to PBL", "")
!    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rqiblten    ,"rqiblten","qi tendency due to PBL","" )
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rqvcuten    ,"rqvcuten" ,"qv tendency due to CU","")
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rqvmpten,"rqvmpten","rqv tendency due to MP","") ! dmp/dt
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rqvblten,"rqvblten","qv tendency due to BL","")
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%u_phy,"qdyn","","")
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%v_phy,"tdyn","","")
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%omega_wrf,"omega_wrf","","")
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%t_phy,"t_phy","","")
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%moist,"moist","","")
    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%p_phy,"p_phy","","")
!    call wrap_add_field_2d(ptend_wrf%rqvdyten,"qdyn","","")
!    call wrap_add_field_2d(ptend_wrf%rthdyten,"tdyn","","")


!    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%rhom,"rhom",    "","" )
!    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%pi_phy       ,"pi_phy","" , ""  )
!    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%dz8w,"dz8w", "", "")
!    call wrap_add_field_2d(diag_phys_accu_vars_h1_2d%zzz    ,"zzz","","")
#endif
#endif

! diag means averaged quant
    !call wrap_add_field_2d(scalar_delp_at_pc_full_level_np1           ,"delp"  )
    !call wrap_add_field_2d(scalar_eta_mass_flux_at_pc_full_level_n    ,"momega") ! generalized
    !call wrap_add_field_2d(diag_physics_vars_h1%uwind                 ,"uwindDiag")
    !call wrap_add_field_2d(diag_physics_vars_h1%vwind                 ,"vwindDiag")
    !call wrap_add_field_2d(diag_physics_vars_h1%temp                  ,"tempDiag")
    !call wrap_add_field_2d(diag_physics_vars_h1%qv                    ,"qvDiag")
    !call wrap_add_field_2d(diag_physics_vars_h1%cldcu                 ,"cldcuDiag")
    !call wrap_add_field_2d(diag_physics_vars_h1%cldst                 ,"cldstDiag")
    !----------------------------------------------------------------------------------------------
    !call wrap_add_field_2d(scalar_delhp_at_pc_full_level_n            ,"delhp" )    ! diagnosed based on hps
    !call wrap_add_field_2d(scalar_www_at_pc_full_level_n              ,"wwwFull"  ) ! diagnosed based on wwwFace 
    !call wrap_add_field_2d(scalar_abs_vor_at_dc_full_level_n          ,"aVor" )     ! diagnosed based on rvor
!    call wrap_add_field_2d(scalar_potential_temp_at_pc_full_level   ,"pt")        ! diagnosed based on temp
    !call wrap_add_field_2d(scalar_potential_temp_at_pc_full_level_iniupt,"pt0")   
    !call wrap_add_field_2d(scalar_geopotential_at_pc_full_level_n     ,"geopotentialFull") ! based on face
    !call wrap_add_field_2d(scalar_pressure_at_pc_full_level_n         ,"pressureFull"  )   ! based on face
    !call wrap_add_field_2d(scalar_mpressure_at_pc_full_level_n        ,"mpressureFull" )   ! based on face
    !call wrap_add_field_2d(scalar_U_wind_at_edge_full_level_n         ,"uEdge"      )      ! no need
    !call wrap_add_field_2d(scalar_V_wind_at_edge_full_level_n         ,"vEdge"      )      ! no need
    !call wrap_add_field_2d(scalar_hpressure_at_pc_full_level_n        ,"hPressureFull" )   ! based on hps
    !call wrap_add_field_2d(scalar_hpressure_at_pc_face_level_n        ,"hPressureFace" )   ! based on hps
    !call wrap_add_field_2d(scalar_delphi_at_pc_full_level_n           ,"delphi")           ! based on geopFace
    !call wrap_add_field_2d(density                                    ,"density")          ! based on delphi/delhp
    !----------------------------------------------------------------------------------------------

    if(trim(working_mode).ne.'dycore'.and.trim(working_mode).ne.'tracer')then
       !call wrap_add_field_2d(scalar_thetav_at_pc_full_level_n      ,"ptv")
    end if
! vars specific to CAM5 physpkg
#ifdef AMIPC_PHYSICS
    !data2d = pstate_cam%diag_relhum
    !call wrap_add_field_2d(pstate_cam%diag_relhum                     ,"relhum" )
    !data2d%f = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,1:nlev,1:mesh%nv_halo(1))
    !call wrap_add_field_2d(data2d                                     ,"cloud"  )
#endif

    END IF

#ifdef SPIO
    call wrap_output_2d_group_sp(mesh, outdir, fname_output)
#else
    call wrap_output_2d_group(mesh, outdir, fname_output)
#endif
!
! clean
!
    call wrap_output_clean_2d()

    return
  end subroutine write_atm_file2d

  subroutine write_atm_file1d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output
! local
   type(scalar_1d_field)               :: scalar_io_1d_template
   integer(i4)                         :: ilev

    if(.not.allocated(scalar_io_1d_template%f)) allocate(scalar_io_1d_template%f(mesh%nv_halo(1)))

! start output
    call wrap_output_init_1d(mesh)
! add variables: field value, varname, longname, units
! geometric
    call wrap_add_field_1d(dycoreVarGeography%scalar_lon_at_pc      , "lon_nv" ,'','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_lat_at_pc      , "lat_nv" ,'','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_lon_at_dc      , "lon_nt" ,'','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_lat_at_dc      , "lat_nt" ,'','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_lon_at_edge    , "lon_ne" ,'','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_lat_at_edge    , "lat_ne" ,'','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_area_at_pc     , "areaPC" ,'spherical area of Primal Cell (Polygon)','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_area_at_dc     , "areaDC" ,'spherical area of Dual Cell (Triangle)','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_leng_at_edp    , "lengEdp",'spherical length of primal-cell edge','')
    call wrap_add_field_1d(dycoreVarGeography%scalar_leng_at_edt    , "lengEdt",'spherical length of dual-cell edge','')
! vertical
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%vert_level     , "nlev" ,'number of full levels','')
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%vertp_level    , "nlevp",'number of face levels','')
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%hyam           , "hyam" ,'full-level hybrid-coordinate coefficient, a','')
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%hybm           , "hybm" ,'full-level hybrid-coordinate coefficient, b','')
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%hyai           , "hyai" ,'face-level hybrid-coordinate coefficient, a','')
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%hybi           , "hybi" ,'face-level hybrid-coordinate coefficient, b','')
!
! dycore vars
!
    call wrap_add_field_1d(dycoreVarSurface%scalar_hpressure_n         , "hps",'surface dry-hydrostatic/mass pressure','Pa')
    call wrap_add_field_1d(dycoreVarSurface%scalar_pressure_n          , "ps" ,'surface full-air pressure','Pa')
#ifndef DYAMOND
    call wrap_add_field_1d(dycoreVarSurface%scalar_geopotential_n      , "phis",'surface geopotential','m^2/s^2')
#endif

    if(trim(working_mode).ne.'dycore'.and.trim(working_mode).ne.'tracer')then
#ifndef DYAMOND
       call wrap_add_field_1d(pstate%scalar_precl_surface             , "precl",'large-scale precipitation rate','m/s') ! dtp
#endif
       if(trim(working_mode).eq.'amipc'.or.trim(working_mode).eq.'amipw')then

#ifndef DYAMOND
         call wrap_add_field_1d(pstate%scalar_precc_surface           , "precc",'convective precipitation rate','m/s')
         call wrap_add_field_1d(pstate%scalar_prect_surface           , "prect",'total precipitation rate','m/s')
         call wrap_add_field_1d(pstate%scalar_grapl_surface           , "precGrapl",'graupel-style precipitation rate, large-scale','m/s') ! seperate diagnosed grap for one step
         call wrap_add_field_1d(pstate%scalar_snowl_surface           , "precSnowl",'snow-style precipitation rate, large-scale','m/s') ! xxxxx snow
         call wrap_add_field_1d(pstate%ts_at_pc_surface               , "ts"       ,'surface temperature','K')
#endif
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%u10m        , "u10m",'10-meter U-wind','m/s')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%v10m        , "v10m",'10-meter V-wind','m/s')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%th2m        , "th2m",'2-meter potential temperature','K')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%t2m         , "t2m" ,'2-meter temperature' ,'K')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%q2m         , "q2m" ,'2-meter mixing ratio','kg/kg')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%hfx         , "hfx" ,'sensible heat flux at surface'  ,'W/m^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%qfx         , "qfx" ,'water-flux at surface','kg/m^2/s')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%lh          , "LH"  ,'latent heat flux at surface','W/m^2')

#ifdef TEND_IO_1d

         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%rainncv  ,"rainncv",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%snowncv  ,"snowncv",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%graupelncv   ,"graupelncv")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%xlat ,"xlat")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%albedo ,"albedo",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%psfc ,"psfc",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%coszr ,"coszr",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%glw  ,"glw")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%olr  ,"olr")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%gsw  ,"gsw")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%znt   ,"znt")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%ust ,"ust")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%pblh ,"pblh",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%psim ,"psim",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%psih ,"psih",'','')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%xland  ,"xland")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%wspd  ,"wspd")
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%br  ,"br")
#endif

#ifdef AMIPW_PHYSICS
! accumulated vars
         call wrap_add_field_1d(pstate%scalar_rainc_surface           , "rainc" ,'Accumulated conv rain'       ,'mm')
         call wrap_add_field_1d(pstate%scalar_rainnc_surface          , "rainnc",'Accumulated non-conv rain'   ,'mm')
         call wrap_add_field_1d(pstate%scalar_snownc_surface          , "snownc",'Accumulated non-conv snow'   ,'mm')
         call wrap_add_field_1d(pstate%scalar_grapnc_surface          , "grapnc",'Accumulated non-conv grapuel','mm')
#endif

!================================================
! write time-averaged state over write frequency
!================================================

#ifndef DYAMOND
! no need to seperation
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%precc       , "preccDiag",'Convective Precipitation rate:A','m/s')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%precl       , "preclDiag",'Grid-scale Precipitation rate(liquid+solid):A','m/s')
#endif
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%snowl       , "precSnowlDiag",'Grid-scale Precipitating Snow rate:A','m/s')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%grapl       , "precGraplDiag",'Grid-scale Precipitating Graupel rate:A','m/s')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%prect       , "prectDiag",'Total Precipitation rate:A','m/s')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%ts          , "tsDiag"  ,'surface temperature:A','K')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%qsfc        , "qsfcDiag",'surface mixing ratio:A','kg/kg')
#ifndef DYAMOND
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%ps          , "psDiag",'surface full-air pressure:A','Pa')
#endif
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%shflx       , "shflxDiag",'sensible heat flux at surface:A','W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%qflx        , "qflxDiag" ,'water-flux at surface:A','kg/m^2/s')
  
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwut       , "flwutDiag" ,'All-sky TOA upward LW flux:A'        ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwdt       , "flwdtDiag" ,'All-sky TOA downward LW flux:A'      ,'W/m^2') 
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwus       , "flwusDiag" ,'All-sky surface upward LW flux:A'    ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwds       , "flwdsDiag" ,'All-sky surface downward LW flux:A'  ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwutc      , "flwutcDiag",'Clear-sky TOA upward LW flux:A'      ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwdtc      , "flwdtcDiag",'Clear-sky TOA downward LW flux:A'    ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwusc      , "flwuscDiag",'Clear-sky surface upward LW flux:A'  ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%flwdsc      , "flwdscDiag",'Clear-sky surface downward LW flux:A','W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%lwcf        , "lwcfDiag"  ,'Long Wave Cloud Forcing:A','W/m^2')

         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswut       , "fswutDiag" ,'All-sky TOA upward SW flux:A'        ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswdt       , "fswdtDiag" ,'All-sky TOA downward SW flux:A'      ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswus       , "fswusDiag" ,'All-sky surface upward SW flux:A'    ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswds       , "fswdsDiag" ,'All-sky surface downward SW flux:A'  ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswutc      , "fswutcDiag",'Clear-sky TOA upward SW flux:A'      ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswdtc      , "fswdtcDiag",'Clear-sky TOA downward SW flux:A'    ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswusc      , "fswuscDiag",'Clear-sky surface upward SW flux:A'  ,'W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%fswdsc      , "fswdscDiag",'Clear-sky surface downward SW flux:A','W/m^2')
         call wrap_add_field_1d(diag_phys_accu_vars_h1_1d%swcf        , "swcfDiag"  ,'Short Wave Cloud Forcing:A','W/m^2')
!================================================
! write time-instant state over write frequency
!================================================
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%viqv        , "viqv"  ,'Vertically-integrated qv(water vapor) mixing ratio','kg/m^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%viqc        , "viqc"  ,'Vertically-integrated qc(cloud water) mixing ratio','kg/m^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%viqr        , "viqr"  ,'Vertically-integrated qr(rain water) mixing ratio','kg/m^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%viqi        , "viqi"  ,'Vertically-integrated qi(cloud ice) mixing ratio','kg/m^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%viqs        , "viqs"  ,'Vertically-integrated qs(snow) mixing ratio','kg/m^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%viqg        , "viqg"  ,'Vertically-integrated qg(graupel) mixing ratio','kg/m^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%cldtot      , "cldtot",'Total cloud fraction','fraction')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%cldlow      , "cldlow",'Low-cloud fraction'  ,'fraction')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%cldmed      , "cldmed",'Mid-cloud fraction'  ,'fraction')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%cldhgh      , "cldhgh",'High-cloud fraction' ,'fraction')
#ifdef AMIPW_PHYSICS
         call wrap_add_field_1d(psurf_wrf%taux ,"taux",'surface zonal stress'     ,'N/m^2')
         call wrap_add_field_1d(psurf_wrf%tauy ,"tauy",'surface meridional stress','N/m^2')
#endif
! vertical-interpolated
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%zzz200      , "zzz200",'200-hPa geopotential','m^2/s^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%zzz500      , "zzz500",'500-hPa geopotential','m^2/s^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%zzz700      , "zzz700",'700-hPa geopotential','m^2/s^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%zzz850      , "zzz850",'850-hPa geopotential','m^2/s^2')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%relhum200   , "relhum200",'200-hPa relative humidity','percent')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%relhum500   , "relhum500",'500-hPa relative humidity','percent')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%relhum700   , "relhum700",'700-hPa relative humidity','percent')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%relhum850   , "relhum850",'850-hPa relative humidity','percent')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%omega200    , "omega200" ,'200-hPa pressure vertical speed','Pa/s')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%omega500    , "omega500" ,'500-hPa pressure vertical speed','Pa/s')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%omega700    , "omega700" ,'700-hPa pressure vertical speed','Pa/s')
         call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%omega850    , "omega850" ,'850-hPa pressure vertical speed','Pa/s')
       end if

       if(trim(working_mode).eq.'amipc'.or.trim(working_mode).eq.'amipw')then
         !call wrap_add_field_1d(pstate%sst_at_pc_surface              , "sstPstate")
         !call wrap_add_field_1d(pstate%ocnfrac_at_pc_surface          , "ocnfrac")
         call wrap_add_field_1d(staticData_sst_at_pc_surface          , "sstRaw"  ,'Raw SST from data, for debug','K') ! for check
         call wrap_add_field_1d(pstate%landfrac_at_pc_surface         , "landfrac",'Land fraction','fraction')

#ifdef USE_NOAHMP

#ifndef DYAMOND
         call wrap_add_field_1d(pstate%atm_in_asdir_at_pc_surface     , "asdir",'albedo SW direct' ,'')  ! filled in dtp2
         call wrap_add_field_1d(pstate%atm_in_asdif_at_pc_surface     , "asdif",'albedo SW diffuse','')  ! filled in dtp2
         call wrap_add_field_1d(pstate%atm_in_aldir_at_pc_surface     , "aldir",'albedo LW direct' ,'')  ! filled in dtp2
         call wrap_add_field_1d(pstate%atm_in_aldif_at_pc_surface     , "aldif",'albedo LW diffuse','')  ! filled in dtp2
#endif
        ! call wrap_add_field_1d(pstate%atm_in_lwup_at_pc_surface      , "lwup")   ! not filled by wrfPhys
        ! call wrap_add_field_1d(pstate%atm_in_shflx_at_pc_surface     , "hfx")    ! filled in dtp2, same as hfx in LND
        ! call wrap_add_field_1d(pstate%qfx_at_pc_surface              , "qfx")    ! filled in dtp2, same as qfx in LND
        ! call wrap_add_field_1d(pstate%atm_out_fswds_at_pc_surface    , "swdown") ! filled in dtp2, as fswds
        ! call wrap_add_field_1d(pstate%atm_out_flwds_at_pc_surface    , "lwdown") ! filled in dtp2, as flwds
        ! call wrap_add_field_1d(pstate%atm_out_netsw_at_pc_surface    , "netsw")  ! filled in dtp2, as fswds-fswus
#endif
       end if
    end if
! outputing
#ifdef SPIO
    call wrap_output_1d_group_sp(mesh, outdir, fname_output)
#else
    call wrap_output_1d_group(mesh, outdir, fname_output)
#endif
! clean
    call wrap_output_clean_1d()
    deallocate(scalar_io_1d_template%f)

   return
  end subroutine write_atm_file1d

  subroutine write_lnd_file1d(mesh,fname_output)
! io
   type(global_domain), intent(in)     :: mesh
   character(len=*)   , intent(in)     :: fname_output
!  
   type(scalar_1d_field) :: io_temp1,io_temp2 ,io_temp3 ,io_temp4
   type(scalar_1d_field) :: io_temp5,io_temp6 ,io_temp7 ,io_temp8
   type(scalar_1d_field) :: io_temp9,io_temp10,io_temp11,io_temp12

    call wrap_output_init_1d(mesh)

#if (defined AMIPW_PHYSICS)

#ifndef DYAMOND
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%tskin  , "tskin")
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%emiss  , "emiss")
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%qsfc   , "qsfc")
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%lh     , "lh")
#endif
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%snowh  , "snowh","land snow height","m")
    call wrap_add_field_1d(diag_phys_inst_vars_h1_1d%xice   , "xice","sea-ice fraction","fraction")
#endif

#ifdef SPIO
    call wrap_output_1d_group_sp(mesh, outdir, fname_output)
#else
    call wrap_output_1d_group(mesh, outdir, fname_output)
#endif
! clean
    call wrap_output_clean_1d()

    return
  end subroutine write_lnd_file1d

  subroutine obtain_atm_fname_output(c_glevel,working_mode,conserve_scheme,&
                                 cyear,cmon,cday,csec,day,sec,dim,fname_output)
! io
  character(len=*),  intent(in)    :: c_glevel, working_mode, conserve_scheme
  character(len=*),  intent(in)    :: cyear, cmon, cday, csec, day, sec, dim
  character(len=*),  intent(inout) :: fname_output

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_output="GRIST.ATM."//trim(grid_info)//"."//trim(working_mode)//"."//&
#ifdef CDATE
                               trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"-"//trim(csec)//&
#else
                               trim(day)//"-"//trim(sec)//&
#endif
                               "."//trim(dim)//".nc"
    return
  end subroutine obtain_atm_fname_output

  subroutine obtain_lnd_fname_output(c_glevel,working_mode,conserve_scheme,&
                                     cyear,cmon,cday,csec,day,sec,dim,fname_output)
! io
  character(len=*),  intent(in)    :: c_glevel, working_mode, conserve_scheme
  character(len=*),  intent(in)    :: cyear, cmon, cday, csec, day, sec, dim
  character(len=*),  intent(inout) :: fname_output

    if(trim(grid_info).eq.'null') grid_info = "G"//trim(c_glevel)

    fname_output="GRIST.LND."//trim(grid_info)//"."//trim(working_mode)//"."//&
#ifdef CDATE
                               trim(cyear)//"-"//trim(cmon)//"-"//trim(cday)//"-"//trim(csec)//&
#else
                               trim(day)//"-"//trim(sec)//&
#endif
                                "."//trim(dim)//".nc"
    return
  end subroutine obtain_lnd_fname_output

  end module grist_gcm_io_h1_module
