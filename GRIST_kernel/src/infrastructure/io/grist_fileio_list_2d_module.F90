
!=============================================================
!  Created by zhangyi on 16/8/5.
!
!    File I/O module tailored for GRIST. 
!    It deals with the general conditions, must be called 
!    by an explicit suit of statements.
!
!    How to use:
!               1) wrap_output_init_2d
!               2) wrap_add_field_2d
!               3) wrap_output_2d
!               4) wrap_output_clean_2d
!=============================================================

 module grist_fileio_list_2d_module

   use grist_constants,      only: i4, r8
   use grist_domain_types  , only: global_domain
   use grist_data_types    , only: scalar_2d_field
   use grist_handle_error,   only: endrun
   use grist_nml_module,     only: nlev
   use grist_wrap_nf,        only: wrap_open,           &
                                   wrap_create,         &
                                   wrap_close,          &
                                   wrap_inq_dim,        &
                                   wrap_inq_dimid,      &
                                   wrap_inq_varid,      &
                                   wrap_put_vara_realx, &
                                   wrap_put_att_text,   &
                                   wrap_def_var,        &
                                   wrap_def_dim,        &
                                   wrap_get_vara_realx

   use grist_list_array

   implicit none

   include 'netcdf.inc'

   private

   public :: wrap_output_init_2d  ,&
             wrap_add_field_2d    ,&
             wrap_output_2d       ,&
             wrap_output_clean_2d ,&
             wrap_get_dim_len   ,&
             wrap_read_2d


   type(array_2d_index_list), pointer   :: scalar_2d_list   ! any 2d shape
   integer(i4)                          :: number_of_scalar
   integer(i4)                          :: mesh_maxvnb
   integer(i4)                          :: mesh_nt
   integer(i4)                          :: mesh_ne
   integer(i4)                          :: mesh_nv

 contains

   subroutine wrap_output_init_2d(mesh)

      type(global_domain), intent(in) :: mesh

      mesh_maxvnb                 = mesh%maxvnb
      mesh_nt                     = mesh%nt
      mesh_ne                     = mesh%ne
      mesh_nv                     = mesh%nv
      number_of_scalar            = 0

      return
   end subroutine wrap_output_init_2d

   subroutine wrap_add_field_2d(scalar,scalar_name)

      type(scalar_2d_field), intent(in)  :: scalar        ! 2d array
      character*(*)       ,  intent(in)  :: scalar_name

      call insert_data_2d(scalar_2d_list, scalar%f, scalar_name, number_of_scalar)

      return
   end subroutine wrap_add_field_2d

   subroutine wrap_output_2d(mesh, itimestep, dtime, outdir, filename)

   implicit none
!io
     type(global_domain),  intent(in)  :: mesh
     integer(i4),          intent(in)  :: itimestep
     real(r8),             intent(in)  :: dtime
     character*(*),        intent(in)  :: outdir
     character*(*),        intent(in)  :: filename
! local
     type(array_2d_index_list),pointer :: ptr
     integer(i4), parameter            :: omode = 0
     integer(i4)                       :: ncid
     integer                           :: ret
     integer(i4)                       :: fst_dim_id(4)           ! dim id for the 1st dim (location)
     integer(i4)                       :: sec_dim_id(7)           ! dim id for the 2nd dim
     integer(i4)                       :: len_dim1, len_dim2
     integer(i4)                       :: i1, i2 

!================================================
!                   set data
!================================================

     print*,"----------------------------------------------------------"
     print*,"    number_of_scalar= "   , number_of_scalar
     print*,"----------------------------------------------------------"

!================================================
!                 create file
!================================================

        call wrap_create (trim(outdir)//trim(filename),NF_64BIT_OFFSET, ncid) ! large-file support 

!================================================
! 1st: define dim
!================================================

        call wrap_def_dim(ncid, 'nverts', mesh_nt     , fst_dim_id(1))
        call wrap_def_dim(ncid, 'nedges', mesh_ne     , fst_dim_id(2))
        call wrap_def_dim(ncid, 'ncells', mesh_nv     , fst_dim_id(3))

        call wrap_def_dim(ncid, 'nmaxvnb'    , mesh_maxvnb , sec_dim_id(1))
        call wrap_def_dim(ncid, 'nlev'       , nlev        , sec_dim_id(2))
        call wrap_def_dim(ncid, 'none'       , 1           , sec_dim_id(3))
        call wrap_def_dim(ncid, 'ntwo'       , 2           , sec_dim_id(4))
        call wrap_def_dim(ncid, 'nthree'     , 3           , sec_dim_id(5))
        call wrap_def_dim(ncid, 'nfour'      , 4           , sec_dim_id(6))
        call wrap_def_dim(ncid, 'neight'     , 8           , sec_dim_id(7))

!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================

        ptr => scalar_2d_list

        do while (associated(ptr))

           len_dim1 = ubound(ptr%var,1)-lbound(ptr%var,1)+1
           len_dim2 = ubound(ptr%var,2)-lbound(ptr%var,2)+1
! NT
           if(len_dim1.eq.mesh_nt.and.len_dim2.eq.mesh_maxvnb)then
              i1 = 1
              i2 = 1
           end if
           if(len_dim1.eq.mesh_nt.and.len_dim2.eq.nlev)then
              i1 = 1
              i2 = 2
           end if
           if(len_dim1.eq.mesh_nt.and.len_dim2.eq.1)then
              i1 = 1
              i2 = 3
           end if
           if(len_dim1.eq.mesh_nt.and.len_dim2.eq.2)then
              i1 = 1
              i2 = 4
           end if
           if(len_dim1.eq.mesh_nt.and.len_dim2.eq.3)then
              i1 = 1
              i2 = 5
           end if
           if(len_dim1.eq.mesh_nt.and.len_dim2.eq.4)then
              i1 = 1
              i2 = 6
           end if
           if(len_dim1.eq.mesh_nt.and.len_dim2.eq.8)then
              i1 = 1
              i2 = 7
           end if
! NE
           if(len_dim1.eq.mesh_ne.and.len_dim2.eq.mesh_maxvnb)then
              i1 = 2
              i2 = 1
           end if
           if(len_dim1.eq.mesh_ne.and.len_dim2.eq.nlev)then
              i1 = 2
              i2 = 2
           end if
           if(len_dim1.eq.mesh_ne.and.len_dim2.eq.1)then
              i1 = 2
              i2 = 3
           end if
           if(len_dim1.eq.mesh_ne.and.len_dim2.eq.2)then
              i1 = 2
              i2 = 4
           end if
           if(len_dim1.eq.mesh_ne.and.len_dim2.eq.3)then
              i1 = 2
              i2 = 5
           end if
           if(len_dim1.eq.mesh_ne.and.len_dim2.eq.4)then
              i1 = 2
              i2 = 6
           end if
           if(len_dim1.eq.mesh_ne.and.len_dim2.eq.8)then
              i1 = 2
              i2 = 7
           end if
! NV
           if(len_dim1.eq.mesh_nv.and.len_dim2.eq.mesh_maxvnb)then
              i1 = 3
              i2 = 1
           end if
           if(len_dim1.eq.mesh_nv.and.len_dim2.eq.nlev)then
              i1 = 3
              i2 = 2
           end if
           if(len_dim1.eq.mesh_nv.and.len_dim2.eq.1)then
              i1 = 3
              i2 = 3
           end if
           if(len_dim1.eq.mesh_nv.and.len_dim2.eq.2)then
              i1 = 3
              i2 = 4
           end if
           if(len_dim1.eq.mesh_nv.and.len_dim2.eq.3)then
              i1 = 3
              i2 = 5
           end if
           if(len_dim1.eq.mesh_nv.and.len_dim2.eq.4)then
              i1 = 3
              i2 = 6
           end if
           if(len_dim1.eq.mesh_nv.and.len_dim2.eq.8)then
              i1 = 3
              i2 = 7
           end if

           call def_var_2d_netcdf(ptr%varname, fst_dim_id(i1),sec_dim_id(i2), ncid, ptr%var_idlist)

           ptr => ptr%next

        end do

!================================================
! 3rd: put att
!================================================

        ptr => scalar_2d_list

        do while (associated(ptr))
           call put_att_2d_netcdf( ncid, ptr%var_idlist , ptr%varname)
           ptr => ptr%next
        end do
 
        ret = nf_enddef(ncid)
        if (ret .ne. nf_noerr) call endrun("enddef wrong in wrap_output_2d@fileio")

!================================================
! 4th, put var
!================================================

        ptr => scalar_2d_list

        do while (associated(ptr))
           len_dim1 = ubound(ptr%var,1)-lbound(ptr%var,1)+1
           len_dim2 = ubound(ptr%var,2)-lbound(ptr%var,2)+1
           call put_var_2d_netcdf(ncid, len_dim1, len_dim2, ptr%var_idlist, ptr%var)
           ptr => ptr%next
        end do

        call wrap_close(ncid)

     return
   end subroutine wrap_output_2d

   subroutine wrap_output_clean_2d()

     deallocate(scalar_2d_list)

     return
   end subroutine wrap_output_clean_2d
!
! read routines
!

   subroutine wrap_get_dim_len(outdir,filename,dimname,dim_len)
! io
     character*(*),         intent(in)   :: outdir
     character*(*),         intent(in)   :: filename
     character*(*),         intent(in)   :: dimname
     integer(i4)  ,         intent(out)  :: dim_len
! local
     integer(i4), parameter              :: omode = 0
     integer(i4)                         :: ncid
     integer(i4)                         :: dimid
     character(len=128)                  :: dim_tmp

       call wrap_open(trim(outdir)//trim(filename), omode, ncid)
       call wrap_inq_dimid (ncid, trim(dimname), dimid)
       call wrap_inq_dim(ncid, dimid, dim_tmp, dim_len)

       return
   end subroutine wrap_get_dim_len

   subroutine wrap_read_2d_with_start_count(outdir,filename,varname,dim2_len,&
                   startin,countin,varout)
! io
     character*(*),         intent(in)    :: outdir
     character*(*),         intent(in)    :: filename
     character*(*),         intent(in)    :: varname
     integer      ,         intent(in)    :: dim2_len
     integer      ,         intent(in)    :: startin
     integer      ,         intent(in)    :: countin
     type(scalar_2d_field), intent(inout) :: varout
! local
     integer(i4), parameter               :: omode = 0
     integer(i4)                          :: varid
     integer(i4)                          :: ncid
     integer(i4)                          :: start(2)
     integer(i4)                          :: count(2)

        print*,"startin",startin
     call wrap_open(trim(outdir)//trim(filename),omode,ncid)
!================================================
!              READ vars at edge
!================================================
     start(1)       = startin
     start(2)       = 1
     count(1)       = startin+countin-1
     count(2)       = dim2_len

     call wrap_inq_varid(ncid, trim(varname), varid)
     call wrap_get_vara_realx(ncid, varid, start, count, varout%f)
     call wrap_close(ncid)

     return

   end subroutine wrap_read_2d_with_start_count

   subroutine wrap_read_2d(outdir,filename,varname,dim1_len,dim2_len,varout)
! io
     character*(*),         intent(in)    :: outdir
     character*(*),         intent(in)    :: filename
     character*(*),         intent(in)    :: varname
     integer      ,         intent(in)    :: dim1_len
     integer      ,         intent(in)    :: dim2_len
     type(scalar_2d_field), intent(inout) :: varout
! local
     integer(i4), parameter               :: omode = 0
     integer(i4)                          :: varid
     integer(i4)                          :: ncid
     integer(i4)                          :: start(2)
     integer(i4)                          :: count(2)

     call wrap_open(trim(outdir)//trim(filename),omode,ncid)
!================================================
!              READ vars at edge
!================================================
     start(1)       = 1
     start(2)       = 1
     count(1)       = dim1_len
     count(2)       = dim2_len

     call wrap_inq_varid(ncid, trim(varname), varid)
     call wrap_get_vara_realx(ncid, varid, start, count, varout%f)
     call wrap_close(ncid)

     return

   end subroutine wrap_read_2d

!---------------------------------------------------------------------------
!
!                         Private routines below
!
!---------------------------------------------------------------------------

   subroutine def_var_2d_netcdf(varname, fst_dim_id, sec_dim_id, ncid, var_idlist)
! io
    character(len=*),         intent(in)    :: varname
    integer(i4),              intent(in)    :: fst_dim_id
    integer(i4),              intent(in)    :: sec_dim_id
    integer(i4),              intent(in)    :: ncid
    integer(i4),              intent(inout) :: var_idlist

! local
    integer(i4)                             :: vdims(2)
    integer(i4)                             :: kk

      
       vdims(1) = fst_dim_id
       vdims(2) = sec_dim_id

       call wrap_def_var(ncid,trim(varname),NF_DOUBLE, 2, vdims(1:2), var_idlist)

    return
  end subroutine def_var_2d_netcdf

  subroutine put_att_2d_netcdf(ncid, var_idlist, varname)
! io
   integer(i4),              intent(in)   :: ncid
   integer(i4),              intent(in)   :: var_idlist
   character(len=*),         intent(in)   :: varname

     call wrap_put_att_text(ncid,var_idlist,'long_name', trim(varname))

     return
  end subroutine put_att_2d_netcdf

  subroutine put_var_2d_netcdf(ncid, dim1_len, dim2_len, var_idlist, scalar_2d)
! io
    integer(i4),              intent(in)  :: ncid
    integer(i4),              intent(in)  :: dim1_len
    integer(i4),              intent(in)  :: dim2_len
    integer(i4),              intent(in)  :: var_idlist
    real(r8),  pointer,       intent(in)  :: scalar_2d(:,:)
! local
    integer(i4)                           :: start(2)
    integer(i4)                           :: count(2)
    integer(i4)                           :: kk

       start(1) = 1
       start(2) = 1
       count(1) = dim1_len
       count(2) = dim2_len

       call wrap_put_vara_realx(ncid, var_idlist, start(1:2), count(1:2), scalar_2d)

       return

   end subroutine put_var_2d_netcdf

 end module grist_fileio_list_2d_module
