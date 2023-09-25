
 !=================================================================
 !  Created by zhangyi on 16/8/5.
 !
 !    File I/O module tailored for GRIST. 
 !    It deals with the general conditions, must be called 
 !    by an explicit suit of statements.
 !
 !    How to use:
 !             1) wrap_output_init_1d
 !             2) wrap_add_field_1d
 !             3) wrap_output_1d
 !             4) wrap_output_clean_1d
 !=================================================================
 
 module grist_fileio_list_1d_module

   use grist_constants,      only: i4, r8
   use grist_domain_types  , only: global_domain
   use grist_data_types    , only: scalar_1d_field
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

   public :: wrap_output_init_1d  ,&
             wrap_add_field_1d    ,&
             wrap_output_1d       ,&
             wrap_output_clean_1d ,&
             wrap_read_1d


   type(array_1d_index_list), pointer   :: scalar_1d_list   ! any 1d shape
   integer(i4)                          :: number_of_scalar
   integer(i4)                          :: mesh_nt
   integer(i4)                          :: mesh_ne
   integer(i4)                          :: mesh_nv

 contains

   subroutine wrap_output_init_1d(mesh)

      type(global_domain), intent(in) :: mesh

      number_of_scalar            = 0
      mesh_nt                     = mesh%nt
      mesh_ne                     = mesh%ne
      mesh_nv                     = mesh%nv

     return
   end subroutine wrap_output_init_1d

   subroutine wrap_add_field_1d(scalar,scalar_name)

      type(scalar_1d_field)  ,  intent(in)  :: scalar        ! 1d array
      character*(*)       ,  intent(in)  :: scalar_name

      call insert_data_1d(scalar_1d_list, scalar%f, scalar_name, number_of_scalar)

      return
   end subroutine wrap_add_field_1d

   subroutine wrap_output_1d(mesh, itimestep, dtime, outdir, filename)

   implicit none
!io
     type(global_domain),  intent(in)  :: mesh
     integer(i4),          intent(in)  :: itimestep
     real(r8),             intent(in)  :: dtime
     character*(*),        intent(in)  :: outdir
     character*(*),        intent(in)  :: filename
! local
     type(array_1d_index_list),pointer :: ptr
     integer(i4), parameter            :: omode = 0
     integer(i4)                       :: ncid
     integer                           :: ret
     integer(i4)                       :: location_id(3)           ! dim id for the 1st dim (location)
     integer(i4)                       :: len_dim
     integer(i4)                       :: i1
     character(len=10)                 :: flag

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
        print*,trim(outdir)//trim(filename)

!================================================
! 1st: define dim
!================================================

        call wrap_def_dim(ncid, 'nverts', mesh_nt , location_id(1))
        call wrap_def_dim(ncid, 'nedges', mesh_ne , location_id(2))
        call wrap_def_dim(ncid, 'ncells', mesh_nv , location_id(3))

!========================================================
! 2nd: define var, need lon_id, lat_id, var_idlist
!========================================================

        ptr => scalar_1d_list

        do while (associated(ptr))
           len_dim = ubound(ptr%var,1)-lbound(ptr%var,1)+1
           if(len_dim.eq.mesh_nt) i1 = 1
           if(len_dim.eq.mesh_ne) i1 = 2
           if(len_dim.eq.mesh_nv) i1 = 3
           call def_var_1d_netcdf(ptr%varname, location_id(i1), ncid, ptr%var_idlist)
           ptr => ptr%next
        end do

!================================================
! 3rd: put att
!================================================

        ptr => scalar_1d_list
        do while (associated(ptr))
           len_dim = ubound(ptr%var,1)-lbound(ptr%var,1)+1
           if(len_dim.eq.mesh_nt) flag = "nt"
           if(len_dim.eq.mesh_ne) flag = "ne"
           if(len_dim.eq.mesh_nv) flag = "nv"
           call put_att_1d_netcdf(ncid, ptr%var_idlist, ptr%varname,flag)
           ptr => ptr%next
        end do
 
        ret = nf_enddef(ncid)
        if (ret .ne. nf_noerr) call endrun("enddef wrong in wrap_output_1d@fileio")

!================================================
! 4th, put var
!================================================

        ptr => scalar_1d_list
        do while (associated(ptr))
           len_dim = ubound(ptr%var,1)-lbound(ptr%var,1)+1
           call put_var_1d_netcdf(ncid, len_dim, ptr%var_idlist, ptr%var)
           ptr => ptr%next
        end do

        call wrap_close(ncid)

     return
   end subroutine wrap_output_1d

   subroutine wrap_output_clean_1d()

     deallocate(scalar_1d_list)

     return
   end subroutine wrap_output_clean_1d
!
! read routines
!

   subroutine wrap_read_1d(outdir,filename,varname,dimlen,varout)
! io
     character*(*),      intent(in)    :: outdir
     character*(*),      intent(in)    :: filename
     character*(*),      intent(in)    :: varname
     integer      ,      intent(in)    :: dimlen
     type(scalar_1d_field), intent(inout) :: varout
! local
     integer, parameter                :: omode =0
     integer                           :: varid
     integer                           :: ncid
     integer                           :: start(1)
     integer                           :: count(1)

     call wrap_open(trim(outdir)//trim(filename),omode,ncid)
!================================================
!              READ vars at edge
!================================================
     start(1)       = 1
     count(1)       = dimlen
     call wrap_inq_varid(ncid, trim(varname), varid)
     call wrap_get_vara_realx(ncid, varid, start, count, varout%f)
     call wrap_close(ncid)

     return
   end subroutine wrap_read_1d


!---------------------------------------------------------------------------
!
!                         Private routines below
!
!---------------------------------------------------------------------------

   subroutine def_var_1d_netcdf(varname, dim_id, ncid, var_idlist)
! io
    character(len=*),         intent(in)    :: varname
    integer(i4),              intent(in)    :: dim_id
    integer(i4),              intent(in)    :: ncid
    integer(i4),              intent(inout) :: var_idlist
! local
    integer(i4)                             :: vdims(1)

       vdims(1) = dim_id
       call wrap_def_var(ncid,trim(varname),NF_DOUBLE, 1, vdims(1), var_idlist)

    return
  end subroutine def_var_1d_netcdf

  subroutine put_att_1d_netcdf(ncid, var_idlist, varname,flag)
! io
   integer(i4),              intent(in)   :: ncid
   integer(i4),              intent(in)   :: var_idlist
   character(len=*),         intent(in)   :: varname
   character(len=*),         intent(in)   :: flag

     call wrap_put_att_text(ncid,var_idlist,'long_name', trim(varname))
!
! conform for cdo remapdis
!
     if(trim(varname).eq."lon_nt".or.&
        trim(varname).eq."lon_ne".or.&
        trim(varname).eq."lon")then
        call wrap_put_att_text(ncid,var_idlist,'units', "degrees_east")
     else if (trim(varname).eq."lat_nt".or.&
              trim(varname).eq."lat_ne".or.&
              trim(varname).eq."lat")then
              call wrap_put_att_text(ncid,var_idlist,'units', "degrees_north")
     else if (trim(flag).eq."nt")then
              call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon_nt lat_nt")
     else if (trim(flag).eq."ne")then
              call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon_ne lat_ne")
     else if (trim(flag).eq."nv")then
              call wrap_put_att_text(ncid,var_idlist,'coordinates', "lon lat")
     end if

   return
  end subroutine put_att_1d_netcdf

  subroutine put_var_1d_netcdf(ncid, dim_len, var_idlist, scalar_1d)
! io
    integer(i4),              intent(in)  :: ncid
    integer(i4),              intent(in)  :: dim_len
    integer(i4),              intent(in)  :: var_idlist
    real(r8),  pointer      , intent(in)  :: scalar_1d(:)
! local
    integer(i4)                           :: start(1)
    integer(i4)                           :: count(1)

       start(1) = 1
       count(1) = dim_len
       call wrap_put_vara_realx(ncid, var_idlist, start(1), count(1), scalar_1d)

      return
   end subroutine put_var_1d_netcdf
 end module grist_fileio_list_1d_module
