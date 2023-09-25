!-----------------------------------------------------------
! Created on 2018
! Authors/Contact: Liu Zhuang
! Version 1.0
! Description: File Par I/O module tailored for GRIST. It 
! deals with the general conditions, must be called by an 
! explicit suit of statements.
!            How to use:
!             1) wrap_output_init_1d
!             2) wrap_add_field_1d
!             3) wrap_output_1d
!             4) wrap_output_clean_1d
!-----------------------------------------------------------

module grist_fileio_list_1d_module_par
   use grist_lib
   use grist_constants,      only: i4, i8, r4, r8
   use grist_domain_types  , only: global_domain, block_structure, group_comm
   use grist_data_types    , only: scalar_1d_field, data1d_temp_dp, data1d_temp_sp
   use grist_handle_error,   only: endrun
   use grist_nml_module,     only: nlev, nlev_inidata

   use grist_wrap_pf,       only: wrap_open,                 &
                                  wrap_create,               &
                                  wrap_enddef,               &
                                  wrap_close,                &
                                  wrap_inq_dim,              &
                                  wrap_inq_dimid,            &
                                  wrap_inq_varid,            &
                                  wrap_inq_vartype,          &
                                  wrap_put_vara_realx,       &
                                  wrap_put_vara_realx_req,   &
                                  wrap_put_vara_real,        &
                                  wrap_put_vara_real_req,    &
                                  wrap_put_att_text,         &
                                  wrap_def_var,              &
                                  wrap_def_dim,              &
                                  wrap_get_vara_realx,       &
                                  wrap_get_vara_realx_req,   &
                                  wrap_get_vara_real,        &
                                  wrap_get_vara_real_req
   use grist_list_array

   implicit none
   include 'pnetcdf.inc'
    
   private

   public :: wrap_output_init_1d    ,&
             wrap_add_field_1d      ,&
             wrap_output_1d_group   ,&
             wrap_output_1d_group_sp,&
             wrap_output_clean_1d   ,&
             wrap_read_1d_group_rst


   type(array_1d_index_list), pointer   :: scalar_1d_list   ! any 1d shape
   integer(i4)                          :: number_of_scalar
   integer(i8)                          :: mesh_nt
   integer(i8)                          :: mesh_ne
   integer(i8)                          :: mesh_nv
   integer(i4), PARAMETER               :: MAXLEN = 268000000 !required by pnetcdf (< 2GB/8B)
   real(8)                              :: default_8 = 1D38
   real(4)                              :: default_4 = 1e38

 contains

   subroutine wrap_output_init_1d(mesh)

      type(global_domain), intent(in) :: mesh

      number_of_scalar            = 0

      mesh_nt                     = mesh%nt_all
      mesh_ne                     = mesh%ne_all
      mesh_nv                     = mesh%nv_all

     return
   end subroutine wrap_output_init_1d

   subroutine wrap_add_field_1d(scalar,scalar_name,longname,units)

      type(scalar_1d_field), intent(in)            :: scalar        ! 1d array
      character*(*)        , intent(in)            :: scalar_name
      character*(*)        , intent(in), optional  :: longname
      character*(*)        , intent(in), optional  :: units

      if(present(longname).and.present(units))then
         call insert_var_1d(scalar_1d_list, scalar%f, trim(scalar_name),trim(longname),trim(units), number_of_scalar, scalar%pos)
      else
         call insert_var_1d(scalar_1d_list, scalar%f, trim(scalar_name),'','', number_of_scalar, scalar%pos)
      end if

      return
   end subroutine wrap_add_field_1d

   subroutine wrap_output_1d_group(mesh, outdir, filename)

   implicit none
!io
     type(global_domain), intent(in), target :: mesh
     character*(*),       intent(in)         :: outdir
     character*(*),       intent(in)         :: filename
! local
     type(group_comm),          pointer      :: comm
     type(block_structure),     pointer      :: local_block
     type(array_1d_index_list), pointer      :: ptr
     integer(i4)                             :: cmode, nfid, ncid, ret
     integer                                 :: err, ierr, varid, omode
     integer(i4)                             :: location_id(5)           ! dim id for the 1st dim (location)
     integer(i4)                             :: i0, i1, i2 ,i, j, lidx
     integer(i8)                             :: start_list(2)
     integer(i8)                             :: len_list(2)
     integer(i8)                             :: nlevp, nlevi8
     real(r8),                  allocatable  :: buff_var(:)
     real(8),                   allocatable  :: buff_put_dp(:)
     real(4),                   allocatable  :: buff_put_sp(:)
     real(r8),                  allocatable  :: output_var(:)
     integer(i4),               pointer      :: displs_send(:), counts_send(:)
     integer(i4),               pointer      :: displs_recv(:), counts_recv(:)
     integer(i4),               pointer      :: idx_recv(:), idx_local(:)
     integer(i4),               pointer      :: mystart, mycount
     integer(i8)                             :: start(2)
     integer(i8)                             :: count(2)
     integer(i4)                             :: req(1), st(1), act_put_num
     integer(i4)                             :: max_diml, jres, jstart

     nlevi8 = nlev
     nlevp= nlevi8+1
     comm => mesh%gcomm_write
     local_block => mesh%local_block
     omode=1

!================================================
!                   set data
!================================================
     start_list(1)= 1
     start_list(2)= 1
     
     len_list(1)= nlevi8
     len_list(2)= nlevp
!================================================
!                 create file
!================================================
     cmode = IOR(NF_CLOBBER,NF_64BIT_DATA)
     if(comm%color==0) call wrap_create(comm%gcomm,trim(outdir)//trim(filename),cmode, nfid) ! large-file support 
     if(mpi_rank() .eq. 0) print*,trim(outdir)//trim(filename)
!================================================
! 1st: define dim
!================================================
     if(comm%color==0)then
       call wrap_def_dim(nfid, 'nverts', mesh_nt , location_id(1))
       call wrap_def_dim(nfid, 'nedges', mesh_ne , location_id(2))
       call wrap_def_dim(nfid, 'ncells', mesh_nv , location_id(3))
       call wrap_def_dim(nfid, 'nlev'       , nlevi8  , location_id(4))
       call wrap_def_dim(nfid, 'nlevp'      , nlevp   , location_id(5))
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_1d_list   
     
     do while (associated(ptr))     
       if( ptr%var_pos.eq.1) i1 = 1
       if( ptr%var_pos.eq.6) i1 = 2
       if( ptr%var_pos.eq.0) i1 = 3
       if( ptr%var_pos.eq.8) i1 = 4
       if( ptr%var_pos.eq.9) i1 = 5
       if(comm%color==0) call def_var_1d_pnetcdf(nfid, ptr%varname, location_id(i1), ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================

     ptr => scalar_1d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_1d_pnetcdf(nfid, ptr%var_idlist, ptr%varname,ptr%longname,ptr%units, ptr%var_pos)
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_enddef(nfid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_1d_list
     do while (associated(ptr))
       if(ptr%var_pos<7)then
         if (ptr%var_pos .eq. 0) then
           displs_send => comm%displs_sendv
           displs_recv => comm%displs_recvv
           counts_send => comm%counts_sendv
           counts_recv => comm%counts_recvv
           idx_recv    => comm%idx_recvv
           idx_local   => comm%idx_localv
           mystart     => comm%mystart_v
           mycount     => comm%mycount_v
         else if (ptr%var_pos .eq. 1) then 
           displs_send => comm%displs_sendt
           displs_recv => comm%displs_recvt
           counts_send => comm%counts_sendt
           counts_recv => comm%counts_recvt
           idx_recv    => comm%idx_recvt
           idx_local   => comm%idx_localt
           mystart     => comm%mystart_t
           mycount     => comm%mycount_t
         else if (ptr%var_pos .eq. 6) then 
           displs_send => comm%displs_sende
           displs_recv => comm%displs_recve
           counts_send => comm%counts_sende
           counts_recv => comm%counts_recve
           idx_recv    => comm%idx_recve
           idx_local   => comm%idx_locale
           mystart     => comm%mystart_e
           mycount     => comm%mycount_e
         end if

         allocate(output_var(size(idx_local)))
         call find_output_var(ptr%var,idx_local,output_var)

         if(comm%color==0)then
           allocate(buff_var(1:mycount))
         else
           allocate(buff_var(1))
         endif
         if (r8 .eq. 8) then
           call MPI_Alltoallv(output_var, counts_send, displs_send, MPI_REAL8, buff_var, &
                              counts_recv, displs_recv, MPI_REAL8, comm%all_comm, err)
         elseif (r8 .eq. 4) then
           call MPI_Alltoallv(output_var, counts_send, displs_send, MPI_REAL, buff_var, &
                              counts_recv, displs_recv, MPI_REAL, comm%all_comm, err)
         end if
         deallocate(output_var)

         if(comm%color==0)then
           ! get the actual number of outputs (due to the limitation of pnetcdf)
           if (r8 .eq. 8) then
             max_diml = MAXLEN
           elseif (r8 .eq. 4) then
             max_diml = 2*MAXLEN
           end if
           act_put_num = ceiling(mycount/real(max_diml, r8))
           if (r8 .eq. 8) then
             allocate(buff_put_dp(1:mycount), source=default_8)
             do i = 1, size(idx_recv)
                lidx = idx_recv(i)
                buff_put_dp(lidx) = buff_var(i)
             enddo
           elseif (r8 .eq. 4) then
             allocate(buff_put_sp(1:mycount), source=default_4)
             do i = 1, size(idx_recv)
                lidx = idx_recv(i)
                buff_put_sp(lidx) = buff_var(i)
             enddo
           end if

           ! may need to put multiple times (each time < 2GB)
           jres = mod(mycount, max_diml)
           do j = 1, act_put_num
             if (j .lt. act_put_num) then
                 count(1) = max_diml
             else
                 count(1) = jres
             end if
             start(1) = mystart + (j-1)*max_diml
             jstart   = 1 + (j-1)*max_diml
             if (r8 .eq. 8) then
               call wrap_put_vara_realx_req(nfid, ptr%var_idlist, start(1), count(1), &
                                            buff_put_dp(jstart:jstart+count(1)-1),req)
             elseif (r8 .eq. 4) then
               call wrap_put_vara_real_req(nfid, ptr%var_idlist, start(1), count(1), &
                                           buff_put_sp(jstart:jstart+count(1)-1),req)
             end if
             err = nfmpi_wait_all(nfid, 1, req, st)
           end do
           if (r8 .eq. 8) then
             deallocate(buff_put_dp)
           elseif (r8 .eq. 4) then
             deallocate(buff_put_sp)
           end if
         endif
         deallocate(buff_var)
       else if(ptr%var_pos>7)then
         if(mpi_rank()==0)then
            call wrap_open(MPI_COMM_SELF,trim(outdir)//trim(filename),omode,ncid)
            call wrap_inq_varid(ncid, ptr%varname, varid)
            if(ptr%var_pos.eq.8)then
               i2 = 1
               call put_var_1d_pnetcdf(ncid, start_list(i2), len_list(i2), varid, ptr%var)
            end if
            if(ptr%var_pos.eq.9)then
               i2 = 2
               call put_var_1d_pnetcdf(ncid, start_list(i2), len_list(i2), varid, ptr%var) 
            end if
            call wrap_close(ncid)
         end if
       end if
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_close(nfid)
     return
   end subroutine wrap_output_1d_group

   subroutine wrap_output_1d_group_sp(mesh, outdir, filename)

   implicit none
!io
     type(global_domain), intent(in), target :: mesh
     character*(*),       intent(in)         :: outdir
     character*(*),       intent(in)         :: filename
! local
     type(group_comm),          pointer      :: comm
     type(block_structure),     pointer      :: local_block
     type(array_1d_index_list), pointer      :: ptr
     integer(i4)                             :: cmode, nfid, ncid, ret
     integer                                 :: err, ierr, varid, omode
     integer(i4)                             :: location_id(5)           ! dim id for the 1st dim (location)
     integer(i4)                             :: i0, i1, i2 ,i, j, lidx
     integer(i8)                             :: start_list(2)
     integer(i8)                             :: len_list(2)
     integer(i8)                             :: nlevp, nlevi8
     real(r8),                  allocatable  :: buff_var(:)
     real(4),                   allocatable  :: buff_put_sp(:)
     real(r8),                  allocatable  :: output_var(:)
     integer(i4),               pointer      :: displs_send(:), counts_send(:)
     integer(i4),               pointer      :: displs_recv(:), counts_recv(:)
     integer(i4),               pointer      :: idx_recv(:), idx_local(:)
     integer(i4),               pointer      :: mystart, mycount
     integer(i8)                             :: start(2)
     integer(i8)                             :: count(2)
     integer(i4)                             :: req(1), st(1), act_put_num
     integer(i4)                             :: max_diml, jres, jstart

     nlevi8 = nlev
     nlevp= nlevi8+1
     comm => mesh%gcomm_write
     local_block => mesh%local_block
     omode=1

!================================================
!                   set data
!================================================
     start_list(1)= 1
     start_list(2)= 1
     
     len_list(1)= nlevi8
     len_list(2)= nlevp
!================================================
!                 create file
!================================================
     cmode = IOR(NF_CLOBBER,NF_64BIT_DATA)
     if(comm%color==0) call wrap_create(comm%gcomm,trim(outdir)//trim(filename),cmode, nfid) ! large-file support 
     if(mpi_rank() .eq. 0) print*,trim(outdir)//trim(filename)
!================================================
! 1st: define dim
!================================================
     if(comm%color==0)then
       call wrap_def_dim(nfid, 'nverts', mesh_nt , location_id(1))
       call wrap_def_dim(nfid, 'nedges', mesh_ne , location_id(2))
       call wrap_def_dim(nfid, 'ncells', mesh_nv , location_id(3))
       call wrap_def_dim(nfid, 'nlev'       , nlevi8  , location_id(4))
       call wrap_def_dim(nfid, 'nlevp'      , nlevp   , location_id(5))
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_1d_list   
     
     do while (associated(ptr))     
       if( ptr%var_pos.eq.1) i1 = 1
       if( ptr%var_pos.eq.6) i1 = 2
       if( ptr%var_pos.eq.0) i1 = 3
       if( ptr%var_pos.eq.8) i1 = 4
       if( ptr%var_pos.eq.9) i1 = 5
       if(comm%color==0) call def_var_1d_pnetcdf_sp(nfid, ptr%varname, location_id(i1), ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================

     ptr => scalar_1d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_1d_pnetcdf(nfid, ptr%var_idlist, ptr%varname,ptr%longname,ptr%units, ptr%var_pos)
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_enddef(nfid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_1d_list
     do while (associated(ptr))
       if(ptr%var_pos<7)then
         if (ptr%var_pos .eq. 0) then
           displs_send => comm%displs_sendv
           displs_recv => comm%displs_recvv
           counts_send => comm%counts_sendv
           counts_recv => comm%counts_recvv
           idx_recv    => comm%idx_recvv
           idx_local   => comm%idx_localv
           mystart     => comm%mystart_v
           mycount     => comm%mycount_v
         else if (ptr%var_pos .eq. 1) then 
           displs_send => comm%displs_sendt
           displs_recv => comm%displs_recvt
           counts_send => comm%counts_sendt
           counts_recv => comm%counts_recvt
           idx_recv    => comm%idx_recvt
           idx_local   => comm%idx_localt
           mystart     => comm%mystart_t
           mycount     => comm%mycount_t
         else if (ptr%var_pos .eq. 6) then 
           displs_send => comm%displs_sende
           displs_recv => comm%displs_recve
           counts_send => comm%counts_sende
           counts_recv => comm%counts_recve
           idx_recv    => comm%idx_recve
           idx_local   => comm%idx_locale
           mystart     => comm%mystart_e
           mycount     => comm%mycount_e
         end if

         allocate(output_var(size(idx_local)))
         call find_output_var(ptr%var,idx_local,output_var)

         if(comm%color==0)then
           allocate(buff_var(1:mycount))
         else
           allocate(buff_var(1))
         endif
         if (r8 .eq. 8) then
           call MPI_Alltoallv(output_var, counts_send, displs_send, MPI_REAL8, buff_var, &
                              counts_recv, displs_recv, MPI_REAL8, comm%all_comm, err)
         elseif (r8 .eq. 4) then
           call MPI_Alltoallv(output_var, counts_send, displs_send, MPI_REAL, buff_var, &
                              counts_recv, displs_recv, MPI_REAL, comm%all_comm, err)
         end if
         deallocate(output_var)

         if(comm%color==0)then
           ! get the actual number of outputs (due to the limitation of pnetcdf)
           max_diml = 2*MAXLEN
           act_put_num = ceiling(mycount/real(max_diml, r8))
           allocate(buff_put_sp(1:mycount), source=default_4)
           do i = 1, size(idx_recv)
              lidx = idx_recv(i)
              buff_put_sp(lidx) = buff_var(i)
           enddo

           ! may need to put multiple times (each time < 2GB)
           jres = mod(mycount, max_diml)
           do j = 1, act_put_num
             if (j .lt. act_put_num) then
                 count(1) = max_diml
             else
                 count(1) = jres
             end if
             start(1) = mystart + (j-1)*max_diml
             jstart   = 1 + (j-1)*max_diml
             call wrap_put_vara_real_req(nfid, ptr%var_idlist, start(1), count(1), &
                                         buff_put_sp(jstart:jstart+count(1)-1),req)
             err = nfmpi_wait_all(nfid, 1, req, st)
           end do
           deallocate(buff_put_sp)
         endif
         deallocate(buff_var)
       else if(ptr%var_pos>7)then
         if(mpi_rank()==0)then
            call wrap_open(MPI_COMM_SELF,trim(outdir)//trim(filename),omode,ncid)
            call wrap_inq_varid(ncid, ptr%varname, varid)
            if(ptr%var_pos.eq.8)then
               i2 = 1
               call put_var_1d_pnetcdf_sp(ncid, start_list(i2), len_list(i2), varid, ptr%var)
            end if
            if(ptr%var_pos.eq.9)then
               i2 = 2
               call put_var_1d_pnetcdf_sp(ncid, start_list(i2), len_list(i2), varid, ptr%var) 
            end if
            call wrap_close(ncid)
         end if
       end if
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_close(nfid)
     return
   end subroutine wrap_output_1d_group_sp

   subroutine wrap_output_clean_1d()

     if(associated(scalar_1d_list)) deallocate(scalar_1d_list)

     return
   end subroutine wrap_output_clean_1d
!
! read routines
!

   subroutine wrap_read_1d_group_rst(comm, outdir, filename, varname, option, varout, dim_1d)
! io
     type(group_comm), intent(in), target :: comm
     character*(*),    intent(in)         :: outdir
     character*(*),    intent(in)         :: filename
     character*(*),    intent(in)         :: varname
     integer                              :: option
     real(r8),         intent(inout)      :: varout(:)
     integer,          optional           :: dim_1d
! local
     integer(i4),      parameter          :: omode = 0
     integer(i4)                          :: act_get_num
     integer(i4)                          :: varid, ncid
     integer(i4),      pointer            :: displs_send(:), counts_send(:)
     integer(i4),      pointer            :: displs_recv(:), counts_recv(:)
     integer(i4),      pointer            :: idx_send(:), idx_local(:)
     integer(i4),      pointer            :: mystart, mycount
     integer(i8)                          :: start(2)
     integer(i8)                          :: count(2)
     integer(i4)                          :: req(1), st(1)
     real(r8),         allocatable        :: buff_var(:), send_var(:), temp_var(:)
     integer                              :: i, err, j, xtype
     integer                              :: max_diml, jres, jstart
     real(4),          allocatable        :: local_data_r4(:)
     real(8),          allocatable        :: local_data_r8(:)

     if (option .le. 6) then
       if (option .eq. 0) then
         displs_send => comm%displs_sendv
         displs_recv => comm%displs_recvv
         counts_send => comm%counts_sendv
         counts_recv => comm%counts_recvv
         idx_send    => comm%idx_sendv
         idx_local   => comm%idx_localv
         mystart     => comm%mystart_v
         mycount     => comm%mycount_v
       else if (option .eq. 1) then
         displs_send => comm%displs_sendt
         displs_recv => comm%displs_recvt
         counts_send => comm%counts_sendt
         counts_recv => comm%counts_recvt
         idx_send    => comm%idx_sendt
         idx_local   => comm%idx_localt
         mystart     => comm%mystart_t
         mycount     => comm%mycount_t
       else if (option .eq. 6) then
         displs_send => comm%displs_sende
         displs_recv => comm%displs_recve
         counts_send => comm%counts_sende
         counts_recv => comm%counts_recve
         idx_send    => comm%idx_sende
         idx_local   => comm%idx_locale
         mystart     => comm%mystart_e
         mycount     => comm%mycount_e
       end if

       if (comm%color == 0 ) then
         call wrap_open(MPI_COMM_SELF,trim(outdir)//trim(filename),omode,ncid)

         call wrap_inq_varid(ncid, trim(varname), varid)
         call wrap_inq_vartype(ncid, varid, xtype)
         ! get the actual number of inputs (due to the limitation of pnetcdf)
         if (xtype .eq. NF_DOUBLE) then
           max_diml = MAXLEN
         elseif (xtype .eq. NF_FLOAT) then
           max_diml = 2*MAXLEN
         end if
         act_get_num = ceiling(mycount/real(max_diml, 8))

         allocate(buff_var(mycount))
         ! may need to get multiple times (each time < 2GB)
         jstart = 1
         jres   = mod(mycount, max_diml)
         do j = 1, act_get_num
           if (j .lt. act_get_num) then
               count(1) = max_diml
           else
               count(1) = jres
           end if
           start(1) = mystart + (j-1)*max_diml
           jstart   = jstart + (j-1)*max_diml
           if (xtype .eq. NF_DOUBLE) then
             allocate(local_data_r8(count(1)))
             call wrap_get_vara_realx_req(ncid, varid, start, count, local_data_r8, req)
             err = nfmpi_wait_all(ncid, 1, req, st)
             buff_var(jstart:jstart+count(1)-1) = local_data_r8
             deallocate(local_data_r8)
           elseif (xtype .eq. NF_FLOAT) then
             allocate(local_data_r4(count(1)))
             call wrap_get_vara_real_req(ncid, varid, start, count, local_data_r4, req)
             err = nfmpi_wait_all(ncid, 1, req, st)
             buff_var(jstart:jstart+count(1)-1) = local_data_r4
             deallocate(local_data_r4)
           end if
         end do
         call wrap_close(ncid)

         allocate(send_var(size(idx_send)))
         do i = 1, size(idx_send)
           send_var(i) = buff_var((idx_send(i)))
         end do
         deallocate(buff_var)
       else
         allocate(send_var(1))
       end if

       allocate(temp_var(size(idx_local)))
       if (r8 .eq. 8) then
         call MPI_Alltoallv(send_var, counts_send, displs_send, MPI_REAL8, temp_var, &
                            counts_recv, displs_recv, MPI_REAL8, comm%all_comm, err)
       elseif (r8 .eq. 4) then
         call MPI_Alltoallv(send_var, counts_send, displs_send, MPI_REAL4, temp_var, &
                            counts_recv, displs_recv, MPI_REAL4, comm%all_comm, err)
       end if
       do i = 1, size(varout)
          varout(i) = temp_var(idx_local(i))
       end do
       deallocate(temp_var)
     else
       call wrap_open(MPI_COMM_SELF,trim(outdir)//trim(filename),omode,ncid)
       call wrap_inq_varid(ncid, trim(varname), varid)
       call wrap_inq_vartype(ncid, varid, xtype)
       start = 1
       !-------------------LiXH add for other data levels-----------------
       if(present(dim_1d))then
          count = dim_1d
       else 
          count(1) = nlev_inidata
          count(2) = nlev_inidata+1
       end if
       if (xtype .eq. NF_DOUBLE) then
       !-------------------LiXH add for other data levels------------------
         if (option .eq. 8) then
            allocate(local_data_r8(count(1)))
            call wrap_get_vara_realx(ncid, varid, start(1), count(1), local_data_r8)
         else if(option .eq. 9)then
            allocate(local_data_r8(count(2)))
            call wrap_get_vara_realx(ncid, varid, start(2), count(2), local_data_r8)
         end if
         varout = local_data_r8
       elseif (xtype .eq. NF_FLOAT) then
       !-------------------LiXH add for other data levels------------------
         if (option .eq. 8) then
            allocate(local_data_r4(count(1)))
            call wrap_get_vara_real(ncid, varid, start(1), count(1), local_data_r4)
         else if(option .eq. 9)then
            allocate(local_data_r4(count(2)))
            call wrap_get_vara_real(ncid, varid, start(2), count(2), local_data_r4)
         end if
         varout = local_data_r4
       end if
       call wrap_close(ncid)
     end if
     if(mpi_rank().eq.0) write(6,*) "Sucessfully read ",trim(varname)," from ", trim(outdir)//trim(filename)
     return
   end subroutine wrap_read_1d_group_rst

!---------------------------------------------------------------------------
!
!                         Private routines below
!
!---------------------------------------------------------------------------

   subroutine def_var_1d_pnetcdf(nfid, varname, dim_id,  var_idlist)
! io
     character(len=*),         intent(in)    :: varname
     integer(i4),              intent(in)    :: dim_id
     integer(i4),              intent(in)    :: nfid
     integer(i4),              intent(inout) :: var_idlist
! local
     integer(i4)                             :: vdims(1)

     vdims(1) = dim_id
     if (r8 .eq. 8) then
       call wrap_def_var(nfid,trim(varname),NF_DOUBLE, 1, vdims(1), var_idlist)
     elseif (r8 .eq. 4) then
       call wrap_def_var(nfid,trim(varname),NF_FLOAT, 1, vdims(1), var_idlist)
     end if

     return
   end subroutine def_var_1d_pnetcdf

   subroutine def_var_1d_pnetcdf_sp(nfid, varname, dim_id,  var_idlist)
! io
     character(len=*),         intent(in)    :: varname
     integer(i4),              intent(in)    :: dim_id
     integer(i4),              intent(in)    :: nfid
     integer(i4),              intent(inout) :: var_idlist
! local
     integer(i4)                             :: vdims(1)

     vdims(1) = dim_id
     call wrap_def_var(nfid,trim(varname),NF_FLOAT, 1, vdims(1), var_idlist)

     return
   end subroutine def_var_1d_pnetcdf_sp

   subroutine put_att_1d_pnetcdf(nfid, var_idlist, varname, longname, units, flag)
! io
     integer(i4),              intent(in)   :: nfid
     integer(i4),              intent(in)   :: var_idlist
     character(len=*),         intent(in)   :: varname
     character(len=*),         intent(in)   :: longname
     character(len=*),         intent(in)   :: units
     integer(i4),              intent(in)   :: flag

     call wrap_put_att_text(nfid,var_idlist,'long_name', trim(longname))

     if(trim(varname).eq."lon_nt".or.&
        trim(varname).eq."lon_ne".or.&
        trim(varname).eq."lon")then
        call wrap_put_att_text(nfid,var_idlist,'standard_name',"longitude")
        call wrap_put_att_text(nfid,var_idlist,'units', "degrees_east")
        call wrap_put_att_text(nfid,var_idlist,'bounds', "lon_bnds") ! for cdo remapycon
     else if (trim(varname).eq."lat_nt".or.&
              trim(varname).eq."lat_ne".or.&
              trim(varname).eq."lat")then
              call wrap_put_att_text(nfid,var_idlist,'standard_name', "latitude")
              call wrap_put_att_text(nfid,var_idlist,'units', "degrees_north")
              call wrap_put_att_text(nfid,var_idlist,'bounds',"lat_bnds") ! for cdo remapycon
     else if (flag.eq.0)then ! cell
        call wrap_put_att_text(nfid,var_idlist,'coordinates', "lon lat")
        call wrap_put_att_text(nfid,var_idlist,'units', trim(units))
     else if (flag.eq.1)then ! vert
        call wrap_put_att_text(nfid,var_idlist,'coordinates', "lon_nt lat_nt")
        call wrap_put_att_text(nfid,var_idlist,'units', trim(units))
     else if (flag.eq.6)then ! edge
        call wrap_put_att_text(nfid,var_idlist,'coordinates', "lon_ne lat_ne")
        call wrap_put_att_text(nfid,var_idlist,'units', trim(units))
     end if

     return
   end subroutine put_att_1d_pnetcdf

   subroutine put_var_1d_pnetcdf(nfid, dim_start, dim_len, var_idlist, scalar_1d)
! io
     integer(i4),              intent(in)  :: nfid
     integer(i8),              intent(in)  :: dim_start
     integer(i8),              intent(in)  :: dim_len
     integer(i4),              intent(in)  :: var_idlist
     real(r8),                 intent(in)  :: scalar_1d(:)
! local
     integer(i8)                           :: start(1)
     integer(i8)                           :: count(1)
     real(8),                  allocatable :: temp_dp(:)
     real(4),                  allocatable :: temp_sp(:)

     start(1) = dim_start
     count(1) = dim_len

     if (r8 .eq. 8) then
       allocate(temp_dp(size(scalar_1d)))
       temp_dp = scalar_1d
       call wrap_put_vara_realx(nfid, var_idlist, start(1), count(1), temp_dp)
       deallocate(temp_dp)
     elseif (r8 .eq. 4) then
       allocate(temp_sp(size(scalar_1d)))
       temp_sp = scalar_1d
       call wrap_put_vara_real(nfid, var_idlist, start(1), count(1), temp_sp)
       deallocate(temp_sp)
     end if

     return
   end subroutine put_var_1d_pnetcdf
   
   subroutine put_var_1d_pnetcdf_sp(nfid, dim_start, dim_len, var_idlist, scalar_1d)
! io
     integer(i4),              intent(in)  :: nfid
     integer(i8),              intent(in)  :: dim_start
     integer(i8),              intent(in)  :: dim_len
     integer(i4),              intent(in)  :: var_idlist
     real(r8),                 intent(in)  :: scalar_1d(:)
! local
     integer(i8)                           :: start(1)
     integer(i8)                           :: count(1)

     real(4),                  allocatable :: tmp(:)

     allocate(tmp(size(scalar_1d)))
     tmp(:) = real(scalar_1d(:), 4)
     start(1) = dim_start
     count(1) = dim_len
     call wrap_put_vara_real(nfid, var_idlist, start(1), count(1), tmp)
     deallocate(tmp)

     return
   end subroutine put_var_1d_pnetcdf_sp
   
   subroutine find_output_var(local_var,output_local_index,output_local_var)
! io
     real(r8), pointer, intent(in)  :: local_var(:)
     integer(i4),       intent(in)  :: output_local_index(:)
     real(r8),          intent(out) :: output_local_var(:)
  
! local
     integer(i8) :: idex
     
     do idex = 1,size(output_local_index)
        output_local_var(idex) = local_var(output_local_index(idex))
     end do
  
   end subroutine find_output_var

 end module grist_fileio_list_1d_module_par
