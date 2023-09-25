
!-----------------------------------------------------------
! Created on 2018
! Authors/Contact: Liu Zhuang
! Version 1.0
! Description: File Par I/O module tailored for GRIST. 
! It deals with the general conditions, must be called by an 
! explicit suit of statements.
!            How to use:
!             1) wrap_output_init_2d
!             2) wrap_add_field_2d
!             3) wrap_output_2d
!             4) wrap_output_clean_2d
!-----------------------------------------------------------

 module grist_fileio_list_2d_module_par
   use grist_lib
   use grist_constants,      only: i4, i8, r4, r8
   use grist_domain_types  , only: global_domain, block_structure, group_comm
   use grist_data_types    , only: scalar_2d_field, data2d_temp_sp, data2d_temp_dp
   use grist_handle_error,   only: endrun
   use grist_nml_module,     only: nlev
   use grist_wrap_pf,        only: wrap_open,                 &
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
                                   wrap_get_vara_realx_req,   &
                                   wrap_get_vara_realx,       &
                                   wrap_get_vara_real_req,    &
                                   wrap_get_vara_real

   use grist_list_array

   implicit none

   include 'pnetcdf.inc'

   private

   public :: wrap_output_init_2d      ,&
             wrap_add_field_2d        ,&
             wrap_output_2d_group     ,&
             wrap_output_2d_group_sp  ,&
             wrap_output_clean_2d     ,&
             wrap_read_2d_group       ,& ! inputdata: ndim*ncell, outputdata, ncell*ndim
             wrap_read_2d_group_rst   ,& ! inputdata: ncell*ndim, outputdata, ndim*ncell
             wrap_read_2d             ,&
             wrap_bcast_2d


   type(array_2d_index_list), pointer   :: scalar_2d_list   ! any 2d shape
   integer(i4)                          :: number_of_scalar
   integer(i4)                          :: mesh_maxvnb
   integer(i8)                          :: mesh_nt
   integer(i8)                          :: mesh_ne
   integer(i8)                          :: mesh_nv
   integer(i4), PARAMETER               :: MAXLEN = 268000000 !required by pnetcdf (< 2GB/8B)
   real(8)                              :: default_8 = 1D38
   real(4)                              :: default_4 = 1e38

 contains

   subroutine wrap_output_init_2d(mesh)

      type(global_domain), intent(in) :: mesh

      mesh_maxvnb                 = mesh%maxvnb
      mesh_nt                     = mesh%nt_all
      mesh_ne                     = mesh%ne_all
      mesh_nv                     = mesh%nv_all
      number_of_scalar            = 0

      return
   end subroutine wrap_output_init_2d

   subroutine wrap_add_field_2d(scalar,scalar_name,longname,units)

      type(scalar_2d_field), intent(in)            :: scalar        ! 2d array
      character*(*)        , intent(in)            :: scalar_name
      character*(*)        , intent(in), optional  :: longname
      character*(*)        , intent(in), optional  :: units

      if(present(longname).and.present(units))then
         call insert_var_2d(scalar_2d_list, scalar%f, trim(scalar_name), trim(longname),trim(units),number_of_scalar, scalar%pos)
      else
         call insert_var_2d(scalar_2d_list, scalar%f, trim(scalar_name), '','',number_of_scalar, scalar%pos)
      end if

      return
   end subroutine wrap_add_field_2d

   subroutine wrap_output_2d_group(mesh, outdir, filename)

   implicit none
!io
     type(global_domain), intent(in), target :: mesh
     character*(*),       intent(in)         :: outdir
     character*(*),       intent(in)         :: filename
! local
     type(group_comm),          pointer      :: comm
     type(block_structure),     pointer      :: local_block
     type(array_2d_index_list), pointer      :: ptr
     integer(i4)                             :: cmode, ncid, err
     integer(i4)                             :: fst_dim_id(6)           ! dim id for the 1st dim (location)
     integer(i4)                             :: sec_dim_id(3)           ! dim id for the 2nd dim
     integer(i4)                             :: len_dim1, mpi_newtype
     integer(i4)                             :: i0, i1, i2 ,i, j, lidx
     integer(i8)                             :: nlevp, nlevi8
     real(r8),                  allocatable  :: buff_var(:), temp_var(:,:)
     real(8),                   allocatable  :: buff_put_dp(:,:) 
     real(4),                   allocatable  :: buff_put_sp(:,:) 
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
!================================================
!                 create file
!================================================
     cmode = IOR(NF_CLOBBER,NF_64BIT_DATA)
     if(comm%color==0) call wrap_create(comm%gcomm,trim(outdir)//trim(filename),cmode, ncid) ! large-file support 
     if(mpi_rank() .eq. 0) print*,trim(outdir)//trim(filename)
!================================================
! 1st: define dim
!================================================
! possible dims for 1st
     if(comm%color==0) then
       call wrap_def_dim(ncid, 'nlev'       , nlevi8      , fst_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nlevp'      , nlevi8+1    , fst_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'dim12'      , int8(12)    , fst_dim_id(3)) ! gcm    
       call wrap_def_dim(ncid, 'dim3'       , int8(3)     , fst_dim_id(4)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim4'       , int8(4)     , fst_dim_id(5)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim7'       , int8(7)     , fst_dim_id(6)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'vertices'   , int8(mesh%maxvnb)   , fst_dim_id(7)) ! 2024, yiz added for lon_bnds,lat_bnds 
! possible dims for 2nd
       call wrap_def_dim(ncid, 'nverts', mesh_nt     , sec_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nedges', mesh_ne     , sec_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'ncells', mesh_nv     , sec_dim_id(3)) ! gcm
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(ptr%var_level.eq.nlevi8)then
          i1 = 1
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.nlevp)then
          i1 = 2
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.12)then
          i1 = 3
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
! lsm restart, added by yiz 2020dec
       if(ptr%var_level.eq.3)then
          i1 = 4
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.4)then
          i1 = 5
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.7)then
          i1 = 6
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       ! added for lat/lon_bnds
       if(ptr%var_level.eq.mesh%maxvnb)then
          i1 = 7
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
! lsm restart, added by yiz 2020dec
       if(comm%color==0) call def_var_2d_pnetcdf(ptr%varname, fst_dim_id(i1),sec_dim_id(i2), ncid, ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================    
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_2d_pnetcdf( ncid, ptr%var_idlist ,ptr%varname,ptr%longname,ptr%units,ptr%var_pos)
       ptr => ptr%next
     end do 
     if(comm%color==0) call wrap_enddef(ncid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_2d_list
     do while (associated(ptr))      
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

       len_dim1 = ptr%var_level
       allocate(temp_var(len_dim1, size(idx_local)))
       call find_output_var_2d(ptr%var, ptr%var_level, idx_local, temp_var)

       if (comm%color == 0) then
         allocate(buff_var(1:mycount*len_dim1))
       else
         allocate(buff_var(1))
       endif

       if (r8 .eq. 8) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL8, mpi_newtype, err)
       elseif (r8 .eq. 4) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL, mpi_newtype, err)
       end if
       call MPI_Type_commit(mpi_newtype, err)
       call MPI_Alltoallv(temp_var, counts_send, displs_send, mpi_newtype, buff_var, &
                          counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)
       deallocate(temp_var)

       if (comm%color == 0) then
         ! get the actual number of outputs (due to the limitation of pnetcdf)
         if (r8 .eq. 8) then
           max_diml = MAXLEN/len_dim1
         elseif (r8 .eq. 4) then
           max_diml = 2*MAXLEN/len_dim1
         end if
         act_put_num = ceiling(mycount/real(max_diml, r8))

         if (r8 .eq. 8) then
           allocate(buff_put_dp(len_dim1, mycount), source=default_8)
           do i = 1, size(idx_recv)
              lidx = idx_recv(i)
              buff_put_dp(1:len_dim1,lidx) = buff_var(len_dim1*(i-1)+1:len_dim1*i)
           end do
         elseif (r8 .eq. 4) then
           allocate(buff_put_sp(len_dim1, mycount), source=default_4)
           do i = 1, size(idx_recv)
              lidx = idx_recv(i)
              buff_put_sp(1:len_dim1,lidx) = buff_var(len_dim1*(i-1)+1:len_dim1*i)
           end do
         end if

         start(1) = 1
         count(1) = len_dim1
         ! may need to put multiple times (each time < 2GB)
         jres = mod(mycount, max_diml)
         do j = 1, act_put_num
           if (j .lt. act_put_num) then
               count(2) = max_diml
           else
               count(2) = jres
           end if
           start(2) = mystart + (j-1)*max_diml
           jstart   = 1 + (j-1)*max_diml
           if (r8 .eq. 8) then
             call wrap_put_vara_realx_req(ncid, ptr%var_idlist, start(1:2), count(1:2), &
                                          buff_put_dp(:, jstart:jstart+count(2)-1),req)
           elseif (r8 .eq. 4) then
             call wrap_put_vara_real_req(ncid, ptr%var_idlist, start(1:2), count(1:2), &
                                         buff_put_sp(:, jstart:jstart+count(2)-1),req)
           end if
           err = nfmpi_wait_all(ncid, 1, req, st)
         end do
         if (r8 .eq. 8) then
           deallocate(buff_put_dp)
         elseif (r8 .eq. 4) then
           deallocate(buff_put_sp)
         end if
       endif
       deallocate(buff_var)
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_close(ncid)
     return
   end subroutine wrap_output_2d_group

   subroutine wrap_output_2d_group_sp(mesh, outdir, filename)

   implicit none
!io
     type(global_domain), intent(in), target :: mesh
     character*(*),       intent(in)         :: outdir
     character*(*),       intent(in)         :: filename
! local
     type(group_comm),          pointer      :: comm
     type(block_structure),     pointer      :: local_block
     type(array_2d_index_list), pointer      :: ptr
     integer(i4)                             :: cmode, ncid, err
     integer(i4)                             :: fst_dim_id(6)           ! dim id for the 1st dim (location)
     integer(i4)                             :: sec_dim_id(3)           ! dim id for the 2nd dim
     integer(i4)                             :: len_dim1, mpi_newtype
     integer(i4)                             :: i0, i1, i2 ,i, j, lidx
     integer(i8)                             :: nlevp, nlevi8
     real(r8),                  allocatable  :: buff_var(:), temp_var(:,:)
     real(4),                   allocatable  :: buff_put_sp(:,:) 
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
!================================================
!                 create file
!================================================
     cmode = IOR(NF_CLOBBER,NF_64BIT_DATA)
     if(comm%color==0) call wrap_create(comm%gcomm,trim(outdir)//trim(filename),cmode, ncid) ! large-file support 
     if(mpi_rank() .eq. 0) print*,trim(outdir)//trim(filename)
!================================================
! 1st: define dim
!================================================
! possible dims for 1st
     if(comm%color==0) then
       call wrap_def_dim(ncid, 'nlev'       , nlevi8      , fst_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nlevp'      , nlevi8+1    , fst_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'dim12'      , int8(12)    , fst_dim_id(3)) ! gcm    
       call wrap_def_dim(ncid, 'dim3'       , int8(3)     , fst_dim_id(4)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim4'       , int8(4)     , fst_dim_id(5)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'dim7'       , int8(7)     , fst_dim_id(6)) ! lsm restart, added by yiz 2020dec
       call wrap_def_dim(ncid, 'vertices'   , int8(mesh%maxvnb), fst_dim_id(7)) ! yiz added for lat_bnds,lon_bnds 
! possible dims for 2nd
       call wrap_def_dim(ncid, 'nverts', mesh_nt     , sec_dim_id(1)) ! gcm
       call wrap_def_dim(ncid, 'nedges', mesh_ne     , sec_dim_id(2)) ! gcm
       call wrap_def_dim(ncid, 'ncells', mesh_nv     , sec_dim_id(3)) ! gcm
     endif
!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(ptr%var_level.eq.nlevi8)then
          i1 = 1
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.nlevp)then
          i1 = 2
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       if(ptr%var_level.eq.12)then
          i1 = 3
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
! lsm restart, added by yiz 2020dec
       if(ptr%var_level.eq.3)then
          i1 = 4
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.4)then
          i1 = 5
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if      
       if(ptr%var_level.eq.7)then
          i1 = 6
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
       ! added for lat/lon_bnds
       if(ptr%var_level.eq.mesh%maxvnb)then
          i1 = 7
          if(ptr%var_pos.eq.1) i2 = 1
          if(ptr%var_pos.eq.6) i2 = 2
          if(ptr%var_pos.eq.0) i2 = 3
       end if
! lsm restart, added by yiz 2020dec
       if(comm%color==0) call def_var_2d_pnetcdf_sp(ptr%varname, fst_dim_id(i1),sec_dim_id(i2), ncid, ptr%var_idlist)
       ptr => ptr%next
     end do

!================================================
! 3rd: put att
!================================================    
     ptr => scalar_2d_list
     do while (associated(ptr))
       if(comm%color==0) call put_att_2d_pnetcdf( ncid, ptr%var_idlist ,ptr%varname,ptr%longname,ptr%units,ptr%var_pos)
       ptr => ptr%next
     end do 
     if(comm%color==0) call wrap_enddef(ncid)

!================================================
! 4th, put var
!================================================

     ptr => scalar_2d_list
     do while (associated(ptr))      
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

       len_dim1 = ptr%var_level
       allocate(temp_var(len_dim1, size(idx_local)))
       call find_output_var_2d(ptr%var, ptr%var_level, idx_local, temp_var)

       if (comm%color == 0) then
         allocate(buff_var(1:mycount*len_dim1))
       else
         allocate(buff_var(1))
       endif

       if (r8 .eq. 8) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL8, mpi_newtype, err)
       elseif (r8 .eq. 4) then
         call MPI_Type_contiguous(len_dim1, MPI_REAL, mpi_newtype, err)
       end if
       call MPI_Type_commit(mpi_newtype, err)
       call MPI_Alltoallv(temp_var, counts_send, displs_send, mpi_newtype, buff_var, &
                          counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)
       deallocate(temp_var)

       if (comm%color == 0) then
         ! get the actual number of outputs (due to the limitation of pnetcdf)
         max_diml = 2*MAXLEN/len_dim1
         act_put_num = ceiling(mycount/real(max_diml, r8))

         allocate(buff_put_sp(len_dim1, mycount), source=default_4)
         do i = 1, size(idx_recv)
            lidx = idx_recv(i)
            buff_put_sp(1:len_dim1,lidx) = buff_var(len_dim1*(i-1)+1:len_dim1*i)
         end do

         start(1) = 1
         count(1) = len_dim1
         ! may need to put multiple times (each time < 2GB)
         jres = mod(mycount, max_diml)
         do j = 1, act_put_num
           if (j .lt. act_put_num) then
               count(2) = max_diml
           else
               count(2) = jres
           end if
           start(2) = mystart + (j-1)*max_diml
           jstart   = 1 + (j-1)*max_diml
           call wrap_put_vara_real_req(ncid, ptr%var_idlist, start(1:2), count(1:2), &
                                       buff_put_sp(:, jstart:jstart+count(2)-1),req)
           err = nfmpi_wait_all(ncid, 1, req, st)
         end do
         deallocate(buff_put_sp)
       endif
       deallocate(buff_var)
       ptr => ptr%next
     end do

     if(comm%color==0) call wrap_close(ncid)
     return
   end subroutine wrap_output_2d_group_sp

   subroutine wrap_output_clean_2d()

     type(array_2d_index_list), pointer :: field_list, tmpf
    
     field_list => scalar_2d_list
     do while(associated(field_list))
       tmpf => field_list%next
       deallocate(field_list)
       field_list => tmpf
     end do
     scalar_2d_list => null()

     return
   end subroutine wrap_output_clean_2d
!
! read routines
!
   subroutine wrap_read_2d_group(comm, outdir, filename, varname, dim2_len, option, varout)
! io
     type(group_comm), intent(in), target :: comm
     character*(*),    intent(in)         :: outdir
     character*(*),    intent(in)         :: filename
     character*(*),    intent(in)         :: varname
     integer(i8),      intent(in)         :: dim2_len
     integer                              :: option
     real(r8),         intent(inout)      :: varout(:,:)
! local
     integer(i4),      parameter          :: omode = 0
     integer(i4)                          :: varid, ncid
     integer(i4)                          :: dim2, act_get_num, mpi_newtype
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
     real(4),          allocatable        :: local_data_r4(:,:)
     real(8),          allocatable        :: local_data_r8(:,:)

     dim2 = dim2_len
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
         max_diml = MAXLEN/dim2
       elseif (xtype .eq. NF_FLOAT) then
         max_diml = 2*MAXLEN/dim2
       end if
       act_get_num = ceiling(mycount/real(max_diml, 8))

       allocate(buff_var(mycount*dim2))
       start(2)       = 1
       count(2)       = dim2_len
       ! may need to get multiple times (each time < 2GB)
       jres   = mod(mycount, max_diml)
       do j = 1, act_get_num
         if (j .lt. act_get_num) then
             count(1) = max_diml
         else
             count(1) = jres
         end if
         start(1) = mystart + (j-1)*max_diml
         jstart   = 1 + (j-1)*max_diml
         if (xtype .eq. NF_DOUBLE) then
           allocate(local_data_r8(count(1), dim2))
           call wrap_get_vara_realx_req(ncid, varid, start, count, local_data_r8, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(1)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r8(i, 1:dim2)
           end do
           deallocate(local_data_r8)
         elseif (xtype .eq. NF_FLOAT) then
           allocate(local_data_r4(count(1), dim2))
           call wrap_get_vara_real_req(ncid, varid, start, count, local_data_r4, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(1)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r4(i, 1:dim2)
           end do
           deallocate(local_data_r4)
         end if
       end do
       call wrap_close(ncid)
       allocate(send_var(size(idx_send)*dim2))
       do i = 1, size(idx_send)
         send_var((i-1)*dim2+1:i*dim2) = buff_var((idx_send(i)-1)*dim2+1:idx_send(i)*dim2)
       end do
       deallocate(buff_var)
     else
       allocate(send_var(1))
     end if

     allocate(temp_var(size(idx_local)*dim2))

     if (r8 .eq. 8) then
       call MPI_Type_contiguous(dim2, MPI_REAL8, mpi_newtype, err)
     elseif (r8 .eq. 4) then
       call MPI_Type_contiguous(dim2, MPI_REAL, mpi_newtype, err)
     end if
     call MPI_Type_commit(mpi_newtype, err)
     call MPI_Alltoallv(send_var, counts_send, displs_send, mpi_newtype, temp_var, &
                        counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)

     do i = 1, size(varout,1)
        varout(i, 1:dim2) = temp_var((idx_local(i)-1)*dim2+1:idx_local(i)*dim2)
     end do

     deallocate(send_var)
     deallocate(temp_var)
     if(mpi_rank().eq.0) write(6,*) "Sucessfully read ",trim(varname)," from ", trim(outdir)//trim(filename)
     return
   end subroutine wrap_read_2d_group

   subroutine wrap_read_2d_group_rst(comm, outdir, filename, varname, dim2_len, option, varout)
! io
     type(group_comm), intent(in), target :: comm
     character*(*),    intent(in)         :: outdir
     character*(*),    intent(in)         :: filename
     character*(*),    intent(in)         :: varname
     integer(i8),      intent(in)         :: dim2_len
     integer                              :: option
     real(r8),         intent(inout)      :: varout(:,:)
! local
     integer(i4),      parameter          :: omode = 0
     integer(i4)                          :: varid, ncid
     integer(i4)                          :: dim2, act_get_num, mpi_newtype
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
     real(4),          allocatable        :: local_data_r4(:,:)
     real(8),          allocatable        :: local_data_r8(:,:)

     dim2 = dim2_len
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
         max_diml = MAXLEN/dim2
       elseif (xtype .eq. NF_FLOAT) then
         max_diml = 2*MAXLEN/dim2
       end if
       act_get_num = ceiling(mycount/real(max_diml, 8))

       allocate(buff_var(mycount*dim2))
       start(1) = 1
       count(1) = dim2_len
       ! may need to get multiple times (each time < 2GB)
       jres   = mod(mycount, max_diml)
       do j = 1, act_get_num
         if (j .lt. act_get_num) then
             count(2) = max_diml
         else
             count(2) = jres
         end if
         start(2) = mystart + (j-1)*max_diml
         jstart   = 1 + (j-1)*max_diml
         if (xtype .eq. NF_DOUBLE) then
           allocate(local_data_r8(dim2, count(2)))
           call wrap_get_vara_realx_req(ncid, varid, start, count, local_data_r8, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(2)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r8(1:dim2, i)
           end do
           deallocate(local_data_r8)
         elseif (xtype .eq. NF_FLOAT) then
           allocate(local_data_r4(dim2, count(2)))
           call wrap_get_vara_real_req(ncid, varid, start, count, local_data_r4, req)
           err = nfmpi_wait_all(ncid, 1, req, st)
           do i = 1, count(2)
             buff_var((jstart+i-2)*dim2+1:(jstart+i-1)*dim2) = local_data_r4(1:dim2, i)
           end do
           deallocate(local_data_r4)
         end if
       end do
       call wrap_close(ncid)
       allocate(send_var(size(idx_send)*dim2))
       do i = 1, size(idx_send)
         send_var((i-1)*dim2+1:i*dim2) = buff_var((idx_send(i)-1)*dim2+1:idx_send(i)*dim2)
       end do
       deallocate(buff_var)
     else
       allocate(send_var(1))
     end if

     allocate(temp_var(size(idx_local)*dim2))

     if (r8 .eq. 8) then
       call MPI_Type_contiguous(dim2, MPI_REAL8, mpi_newtype, err)
     elseif (r8 .eq. 4) then
       call MPI_Type_contiguous(dim2, MPI_REAL, mpi_newtype, err)
     end if
     call MPI_Type_commit(mpi_newtype, err)
     call MPI_Alltoallv(send_var, counts_send, displs_send, mpi_newtype, temp_var, &
                        counts_recv, displs_recv, mpi_newtype, comm%all_comm, err)

     do i = 1, size(varout,2)
        varout(1:dim2, i) = temp_var((idx_local(i)-1)*dim2+1:idx_local(i)*dim2)
     end do

     deallocate(send_var)
     deallocate(temp_var)
     if(mpi_rank().eq.0) write(6,*) "Sucessfully read ",trim(varname)," from ", trim(outdir)//trim(filename)
     return
   end subroutine wrap_read_2d_group_rst

   subroutine wrap_read_2d(comm,outdir,filename,varname,dim1_len,dim2_len,varout)
! io
     integer,               intent(in)    :: comm
     character*(*),         intent(in)    :: outdir
     character*(*),         intent(in)    :: filename
     character*(*),         intent(in)    :: varname
     integer(i8),           intent(in)    :: dim1_len
     integer(i8),           intent(in)    :: dim2_len
     type(scalar_2d_field), intent(inout) :: varout
! local
     integer(i4),           parameter     :: omode = 0
     integer(i4)                          :: varid
     integer(i4)                          :: ncid
     integer(i4)                          :: xtype
     integer(i8)                          :: start(2)
     integer(i8)                          :: count(2)
     integer                              :: max_diml, act_get_num, j, jres
! temp arrays
     real(8),               allocatable   :: temp_dp(:,:)
     real(4),               allocatable   :: temp_sp(:,:)

     call wrap_open(comm,trim(outdir)//trim(filename),omode,ncid)
     call wrap_inq_varid(ncid, trim(varname), varid)
     call wrap_inq_vartype(ncid, varid, xtype)
!================================================
!              READ vars at edge
!================================================
     start(2)       = 1
     count(2)       = dim2_len
     if (xtype .eq. NF_DOUBLE) then
       max_diml = MAXLEN/dim2_len
     elseif (xtype .eq. NF_FLOAT) then
       max_diml = 2*MAXLEN/dim2_len
     end if
     act_get_num = ceiling(dim1_len/real(max_diml, 8))
     jres = mod(dim1_len, max_diml)

     do j = 1, act_get_num
       if (j .lt. act_get_num) then
           count(1) = max_diml
       else
           count(1) = jres
       end if
       start(1) = 1 + (j-1)*max_diml
       if (xtype .eq. NF_DOUBLE) then
         allocate(temp_dp(count(1), dim2_len))
         call wrap_get_vara_realx(ncid, varid, start, count, temp_dp)
         varout%f(start(1):(start(1)+count(1)-1), 1:dim2_len) = temp_dp
         deallocate(temp_dp)
       elseif (xtype .eq. NF_FLOAT) then
         allocate(temp_sp(count(1), dim2_len))
         call wrap_get_vara_real(ncid, varid, start, count, temp_sp)
         varout%f(start(1):(start(1)+count(1)-1), 1:dim2_len) = temp_sp
         deallocate(temp_sp)
       end if
     end do
     call wrap_close(ncid)

     return

   end subroutine wrap_read_2d
   
   subroutine wrap_bcast_2d(comm, irank, dim1_len, dim2_len, varout)
! io
     integer,               intent(in)    :: comm
     integer,               intent(in)    :: irank
     integer(i8),           intent(in)    :: dim1_len
     integer(i8),           intent(in)    :: dim2_len
     type(scalar_2d_field), intent(inout) :: varout

! local
     integer                              :: err
     integer                              :: var_size
     integer                              :: root_rank
     real(r8), allocatable                :: buf(:,:)
     
     allocate(buf(dim1_len, dim2_len))
     
     buf = varout%f
     var_size = dim1_len*dim2_len
     root_rank = irank
     if (r8 .eq. 8) then
       call MPI_Bcast(buf, var_size, MPI_REAL8, root_rank, comm, err)
     elseif (r8 .eq. 4) then
       call MPI_Bcast(buf, var_size, MPI_REAL, root_rank, comm, err)
     end if
     varout%f = buf
     deallocate(buf)
   end subroutine wrap_bcast_2d

!---------------------------------------------------------------------------
!
!                         Private routines below
!
!---------------------------------------------------------------------------

   subroutine def_var_2d_pnetcdf(varname, fst_dim_id, sec_dim_id, ncid, var_idlist)
! io
     character(len=*), intent(in)    :: varname
     integer(i4),      intent(in)    :: fst_dim_id
     integer(i4),      intent(in)    :: sec_dim_id
     integer(i4),      intent(in)    :: ncid
     integer(i4),      intent(inout) :: var_idlist

! local
     integer(i4)                     :: vdims(2)
     integer(i4)                     :: len_two = 2
    
     vdims(1) = fst_dim_id
     vdims(2) = sec_dim_id
     if (r8 .eq. 8) then
       call wrap_def_var(ncid,trim(varname),NF_DOUBLE, len_two, vdims(1:2), var_idlist)
     elseif (r8 .eq. 4) then
       call wrap_def_var(ncid,trim(varname),NF_FLOAT, len_two, vdims(1:2), var_idlist)
     end if

     return
   end subroutine def_var_2d_pnetcdf

   subroutine def_var_2d_pnetcdf_sp(varname, fst_dim_id, sec_dim_id, ncid, var_idlist)
! io
     character(len=*), intent(in)    :: varname
     integer(i4),      intent(in)    :: fst_dim_id
     integer(i4),      intent(in)    :: sec_dim_id
     integer(i4),      intent(in)    :: ncid
     integer(i4),      intent(inout) :: var_idlist

! local
     integer(i4)                     :: vdims(2)
     integer(i4)                     :: len_two = 2
    
     vdims(1) = fst_dim_id
     vdims(2) = sec_dim_id
     call wrap_def_var(ncid,trim(varname),NF_FLOAT, len_two, vdims(1:2), var_idlist)

     return
   end subroutine def_var_2d_pnetcdf_sp

   subroutine put_att_2d_pnetcdf(ncid, var_idlist, varname, longname,units, var_pos)
! io
     integer(i4),              intent(in)   :: ncid
     integer(i4),              intent(in)   :: var_idlist
     character(len=*),         intent(in)   :: varname
     character(len=*),         intent(in)   :: longname
     character(len=*),         intent(in)   :: units
     integer(i4),              intent(in)   :: var_pos

     call wrap_put_att_text(ncid,var_idlist,'long_name', trim(longname))

     if(trim(varname).eq."lon_nt".or.&
        trim(varname).eq."lon_ne".or.&
        trim(varname).eq."lon")then
        call wrap_put_att_text(ncid,var_idlist,'units', "degrees_east")
     else if (trim(varname).eq."lat_nt".or.&
              trim(varname).eq."lat_ne".or.&
              trim(varname).eq."lat")then
              call wrap_put_att_text(ncid,var_idlist,'units', "degrees_north")
     else if(var_pos.eq.1)then
        call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon_nt lat_nt")
        call wrap_put_att_text(ncid,var_idlist,'units', trim(units))
     else if(var_pos.eq.6)then
        call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon_ne lat_ne")
        call wrap_put_att_text(ncid,var_idlist,'units', trim(units))
     else if(var_pos.eq.0)then
        call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon lat")
        call wrap_put_att_text(ncid,var_idlist,'units', trim(units))
     end if 

     return

   end subroutine put_att_2d_pnetcdf

   subroutine find_output_var_2d(local_var,level_var,output_local_index,output_local_var)
! io
     real(r8),    pointer, intent(in)    :: local_var(:,:)
     integer(i8), intent(in)             :: level_var
     integer(i4), intent(in)             :: output_local_index(:)
     real(r8),    intent(out)            :: output_local_var(:,:)
  
! local
     integer(i8) :: idex,idey
     
     do idex = 1,size(output_local_index)
        do idey = 1,level_var 
           output_local_var(idey,idex) = local_var(idey,output_local_index(idex))
        end do
     end do

   end subroutine find_output_var_2d

 end module grist_fileio_list_2d_module_par
