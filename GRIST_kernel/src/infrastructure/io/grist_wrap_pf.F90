
!-----------------------------------------------------------
! Created on 2018
! Version 1.0
! Description: Wrapper routines for the pnetCDF library for 
! output data.
!-----------------------------------------------------------

module grist_wrap_pf
  use grist_lib

  public

CONTAINS

  subroutine wrap_create (comm, path, cmode, nfid)
    implicit none
    !------------------------------------------------------
    !
    ! Purpose:
    !
    ! Create a netCDF file for reading and/or writing
    !
    !------------------------------------------------------

    include 'pnetcdf.inc'
    integer, intent(in):: comm
    character*(*), intent(in):: path
    integer(4), intent(in):: cmode
    integer(4), intent(out):: nfid

    integer ret      ! pNetCDF return code
    integer ::mpi_info = MPI_INFO_NULL

    ret = nfmpi_create(comm, path, cmode, mpi_info, nfid)

    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_create

  !===============================================================================

  subroutine wrap_inq_dim (nfid, dimid, dimname, dimlen)
    implicit none
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets dimension name for a given dimension id
    !
    !------------------------------------------------------------
    include 'pnetcdf.inc'
    integer, intent(in)::  nfid
    integer, intent(in):: dimid
    integer(8), intent(out):: dimlen
    character*(*), intent(out):: dimname

    integer ret      ! NetCDF return code

    ret = nfmpi_inq_dim (nfid, dimid, dimname, dimlen)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_inq_dim

  !===============================================================================

  subroutine wrap_inq_dimid (nfid, dimname, dimid)
    implicit none
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets the dimension id
    !
    !------------------------------------------------------------
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(out):: dimid
    character*(*), intent(in):: dimname

    integer ret      ! NetCDF return code

    ret = nfmpi_inq_dimid (nfid, dimname, dimid)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_inq_dimid

  !===============================================================================

  subroutine wrap_inq_dimlen (nfid, dimid, dimlen)
    implicit none
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets the dimension length for a given dimension
    !
    !------------------------------------------------------------
    include 'pnetcdf.inc'

    integer, intent(in)::  nfid
    integer, intent(in)::  dimid 
    integer, intent(out):: dimlen

    integer ret      ! NetCDF return code

    ret = nfmpi_inq_dimlen (nfid, dimid, dimlen)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_inq_dimlen

  !===============================================================================

  subroutine wrap_inq_vardimid (nfid, varid, dimids)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Returns the dimension Id's from a variable
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(out)::  dimids(*)

    integer ret      ! NetCDF return code

    ret = nfmpi_inq_vardimid (nfid, varid, dimids)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_inq_vardimid

  !===============================================================================

  subroutine wrap_inq_varid (nfid, varname, varid)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Returns the variable ID
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(out):: varid
    character*(*), intent(in):: varname

    integer ret      ! NetCDF return code

    ret = nfmpi_inq_varid (nfid, varname, varid)
    if (ret/=NF_NOERR) then
       write(6,*)'wrap_inq_varid: id for ',trim(varname),' not found'
       call handle_error (ret)
    end if
  end subroutine wrap_inq_varid

  !===============================================================================

  subroutine wrap_inq_var (nfid, varid, varname, xtype, ndims, &
       dimids, natts)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Returns the variable name, type, number of dimensions, dimension ID's, and number of attributes
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(out):: xtype
    integer, intent(out):: ndims
    integer, intent(out):: dimids(*)
    integer, intent(out):: natts
    character*(*), intent(out):: varname

    integer ret      ! NetCDF return code

    ret = nfmpi_inq_var (nfid, varid, varname, xtype, ndims, dimids, &
         natts)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_inq_var

  !===============================================================================

  subroutine wrap_inq_varname (nfid, varid, varname)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Returns the variable name from the dimension ID
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer ret      ! NetCDF return code

    integer, intent(in):: nfid
    integer, intent(in):: varid
    character*(*), intent(out):: varname

    ret = nfmpi_inq_varname (nfid, varid, varname)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_inq_varname

  !===============================================================================

  subroutine wrap_get_att_text (nfid, varid, attname, atttext)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Returns the attribute text from the given variable ID and attribute name
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    character*(*), intent(in):: attname
    character*(*), intent(out):: atttext

    integer ret      ! NetCDF return code

    ret = nfmpi_get_att_text (nfid, varid, attname, atttext)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_get_att_text

  !===============================================================================

  subroutine wrap_put_att_text (nfid, varid, attname, atttext)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Puts the given attribute text to variable ID.
    !
    ! This routine violates the convetion that the wrapper codes take an identical
    ! set of arguments as the netcdf library code.  The length of the character
    ! argument is computed inside the wrapper.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    character*(*), intent(in):: attname
    character*(*), intent(in):: atttext

    integer :: ret      ! NetCDF return code
    integer(8) :: siz

    siz = len_trim(atttext)

    ret = nfmpi_put_att_text (nfid, varid, attname, siz, atttext)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_att_text

  !===============================================================================

  subroutine wrap_put_att_realx (nfid, varid, attname, xtype, len, &
       attval)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Puts the given real attribute to the variable id
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer , intent(in):: nfid
    integer , intent(in):: varid
    integer , intent(in):: xtype
    integer , intent(in):: len
    character*(*) , intent(in):: attname
    real(8) , intent(in):: attval

    integer :: ret      ! NetCDF return code

    !#ifdef SPIO
    !   ret = nf_put_att_real (nfid, varid, attname, xtype, len, attval)
    !#else
    ret = nfmpi_put_att_double (nfid, varid, attname, xtype, len, attval)
    !#endif
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_att_realx
  !===============================================================================

  subroutine wrap_def_dim (nfid, dimname, len, dimid)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Defines the input dimension
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'
    integer(4), intent(in):: nfid
    character*(*), intent(in):: dimname
    integer(8), intent(in):: len
    integer(4), intent(out):: dimid

    integer :: ret      ! NetCDF return code


    ret = nfmpi_def_dim(nfid, dimname, len, dimid)

    if (ret/=NF_NOERR) call handle_error (ret)

  end subroutine wrap_def_dim

  !===============================================================================

  subroutine wrap_def_var (nfid, name, xtype, nvdims, vdims, varid)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Defines the given variable
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer(4), intent(in)::xtype
    integer(4), intent(in)::nvdims
    integer(4), intent(out)::varid
    integer(4), intent(in):: vdims(nvdims)
    character*(*), intent(in):: name

    integer :: ret      ! NetCDF return code

    ret = nfmpi_def_var (nfid, name, xtype, nvdims, vdims, varid)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_def_var

  !===============================================================================

   subroutine wrap_get_vara_real_req (nfid, varid, start, count, arr, req)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given real variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none
!include 'netcdf.inc'
include 'pnetcdf.inc'

   integer(4), intent(in):: nfid
   integer(4), intent(in)::varid
   integer(8), intent(in)::start(*)
   integer(8), intent(in)::count(*)
   real(4),    intent(out):: arr(*)
   integer(4), intent(out):: req(*)

   integer ret      ! NetCDF return code

   ret = nfmpi_iget_vara_real (nfid, varid, start, count, arr, req)

   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VARA_REAL_REQ: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_vara_real_req

  !===============================================================================

   subroutine wrap_get_vara_realx_req (nfid, varid, start, count, arr, req)
!-------------------------------------------------------------------------------
!
! Purpose:
!
! Gets a range of the given real variable from a input file
!
!-------------------------------------------------------------------------------
   implicit none
!include 'netcdf.inc'
include 'pnetcdf.inc'

   integer(4), intent(in):: nfid
   integer(4), intent(in)::varid
   integer(8), intent(in)::start(*)
   integer(8), intent(in)::count(*)
   real(8),    intent(out):: arr(*)
   integer(4), intent(out):: req(*)

   integer ret      ! NetCDF return code

   ret = nfmpi_iget_vara_double (nfid, varid, start, count, arr, req)

   if (ret/=NF_NOERR) then
     write(6,*)'WRAP_GET_VARA_REALX: error reading varid =', varid
     call handle_error (ret)
   end if
   end subroutine wrap_get_vara_realx_req

  !===============================================================================

  subroutine wrap_get_var_real (nfid, varid, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets the given real variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    real(4), intent(out):: arr(*)

    integer ret      ! NetCDF return code

    !#ifdef SPIO
    !   ret = nf_get_var_real (nfid, varid, arr)
    !#else
        !stop
    ret = nfmpi_begin_indep_data(nfid)
    ret = nfmpi_get_var_real (nfid, varid, arr)
    ret = nfmpi_end_indep_data(nfid)
    !#endif
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VAR_REAL: error reading varid =', varid
       call handle_error (ret)
    end if
  end subroutine wrap_get_var_real

  !===============================================================================

  subroutine wrap_get_var_realx (nfid, varid, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets the given real variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    real(8), intent(out):: arr(*)

    integer ret      ! NetCDF return code

    !#ifdef SPIO
    !   ret = nf_get_var_real (nfid, varid, arr)
    !#else
        !stop
    ret = nfmpi_begin_indep_data(nfid)
    ret = nfmpi_get_var_double (nfid, varid, arr)
    ret = nfmpi_end_indep_data(nfid)
    !#endif
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VAR_REALX: error reading varid =', varid
       call handle_error (ret)
    end if
  end subroutine wrap_get_var_realx

  !===============================================================================

  subroutine wrap_get_var_int (nfid, varid, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets the given integer variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer(8), intent(out):: arr(*)

    integer ret      ! NetCDF return code

    ret = nfmpi_begin_indep_data(nfid)
    ret = nfmpi_get_var_int8 (nfid, varid, arr)
    ret = nfmpi_end_indep_data(nfid)

    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VAR_INT: error reading varid =', varid
       call handle_error (ret)
    end if
  end subroutine wrap_get_var_int

  !===============================================================================

  subroutine wrap_get_vara_real (nfid, varid,start, count, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets a range of the given real variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in)::varid
    integer(8), intent(in)::start(*)
    integer(8), intent(in)::count(*)
    real(4), intent(out):: arr(*)

    integer ret      ! NetCDF return code
    !#ifdef SPIO
    !   ret = nf_get_vara_real (nfid, varid, start, count, arr)
    !#else


    ret = nfmpi_get_vara_real_all (nfid, varid, start, count, arr)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VARA_REAL: error reading varid =', varid
       call handle_error (ret)
    end if

    !#endif
    
  end subroutine wrap_get_vara_real

  !===============================================================================

  subroutine wrap_get_vara_realx (nfid, varid,start, count, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets a range of the given real variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in)::varid
    integer(8), intent(in)::start(*)
    integer(8), intent(in)::count(*)
    real(8), intent(out):: arr(*)

    integer ret      ! NetCDF return code
    !#ifdef SPIO
    !   ret = nf_get_vara_real (nfid, varid, start, count, arr)
    !#else


    ret = nfmpi_get_vara_double_all (nfid, varid, start, count, arr)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VARA_REALX: error reading varid =', varid
       call handle_error (ret)
    end if

    !#endif
    
  end subroutine wrap_get_vara_realx

  !===============================================================================

  subroutine wrap_get_vara_int (nfid, varid, start, count, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Gets a range of the given integer variable from a input file.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(in):: start(*)
    integer, intent(in):: count(*)
    integer, intent(out):: arr(*)

    integer ret      ! NetCDF return code
    ret = nfmpi_begin_indep_data(nfid)
    ret = nfmpi_get_vara_int (nfid, varid, start, count, arr)
    ret = nfmpi_end_indep_data(nfid)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VARA_INT: error reading varid =', varid
       call handle_error (ret)
    end if
  end subroutine wrap_get_vara_int

  !===============================================================================

  subroutine wrap_open (comm, path, omode, ncid)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Open a netCDF file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in) :: comm
    character*(*), intent(in):: path
    integer, intent(in):: omode
    integer, intent(out):: ncid

    integer ret      ! NetCDF return code
    integer :: mpi_info = MPI_INFO_NULL

    ret = nfmpi_open (comm, path, omode, mpi_info, ncid)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_OPEN: nf_open failed for file ',path
       call handle_error (ret)
    end if
  end subroutine wrap_open

  !===============================================================================

  subroutine wrap_close (nfid)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Close netCDF file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid

    integer ret      ! NetCDF return code

    ret = nfmpi_close (nfid)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_CLOSE: nf_close failed for id ',nfid
       call handle_error (ret)
    end if
  end subroutine wrap_close

  !===============================================================================

   subroutine wrap_enddef (nfid)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Close netCDF file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid

    integer ret      ! NetCDF return code

    ret = nfmpi_enddef (nfid)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_ENDDEF: nfmpi_enddef failed for id ',nfid
       call handle_error (ret)
    end if
  end subroutine wrap_enddef

  !===============================================================================
  
  subroutine wrap_put_var_int (nfid, varid, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Put a integer variable on output file.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(in):: arr(*)

    integer ret      ! NetCDF return code

    ret = nfmpi_put_var_int (nfid, varid, arr)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_var_int

  !===============================================================================

  subroutine wrap_put_var1_int (nfid, varid, index, ival)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Put a variable on output file at a given index.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(in):: index(*)
    integer, intent(in):: ival

    integer ret      ! NetCDF return code

    ret = nfmpi_put_var1_int (nfid, varid, index, ival)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_var1_int

  !===============================================================================

  subroutine wrap_put_vara_int (nfid, varid, start, count, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Put a range of a integer variable on a output file.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(in):: start(*)
    integer, intent(in):: count(*)
    integer, intent(in):: arr(*)

    integer ret      ! NetCDF return code

    ret = nfmpi_put_vara_int (nfid, varid, start, count, arr)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_vara_int

  !===============================================================================

  subroutine wrap_put_vara_text (nfid, varid, start, count, text)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Put a range of the given text variable to output file.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(in):: start(*)
    integer, intent(in):: count(*)
    character*(*), intent(in):: text(*)

    integer ret      ! NetCDF return code

    ret = nfmpi_put_vara_text (nfid, varid, start, count, text)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_vara_text

  !===============================================================================

  subroutine wrap_put_var1_realx (nfid, varid, index, val)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Put the given real variable to output file at given index.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(in):: index(*)
    real(8), intent(in):: val

    integer ret      ! NetCDF return code

    !#ifdef SPIO
    !   ret = nf_put_var1_real (nfid, varid, index, val)
    !#else
    ret = nfmpi_put_var1_double (nfid, varid, index, val)
    !#endif
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_var1_realx

  !===============================================================================

  subroutine wrap_put_vara_realx (nfid, varid, start, count, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Output the given portion of the real array.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer(8), intent(in):: start(*)
    integer(8), intent(in):: count(*)
    real(8), intent(in):: arr(*)

    integer :: ret      ! NetCDF return code

    ret = nfmpi_put_vara_double_all(nfid, varid, start, count, arr)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_vara_realx

  !===============================================================================

  subroutine wrap_put_vara_realx_req (nfid, varid, start, count, arr, req) 
    !------------------------------------------------------------
    !    
    ! Purpose:
    !    
    ! Output the given portion of the real array.
    !    
    !------------------------------------------------------------
    implicit none 
    include 'pnetcdf.inc'

    integer, intent(in):: nfid 
    integer, intent(in):: varid
    integer(8), intent(in):: start(*)
    integer(8), intent(in):: count(*)
    integer(4), intent(out):: req(*)
    real(8), intent(in):: arr(*)

    integer :: ret      ! NetCDF return code

    ret = nfmpi_iput_vara_double(nfid, varid, start, count, arr,req)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_vara_realx_req

  !===============================================================================

  subroutine wrap_put_vara_real (nfid, varid, start, count, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Output the given portion of the real array.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer,    intent(in):: nfid
    integer,    intent(in):: varid
    integer(8), intent(in):: start(*)
    integer(8), intent(in):: count(*)
    real(4),    intent(in):: arr(*)

    integer ret      ! NetCDF return code
    ret = nfmpi_put_vara_real_all (nfid, varid, start, count, arr)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_vara_real

  !===============================================================================

  subroutine wrap_put_vara_real_req (nfid, varid, start, count, arr, req) 
    !------------------------------------------------------------
    !    
    ! Purpose:
    !    
    ! Output the given portion of the real array.
    !    
    !------------------------------------------------------------
    implicit none 
    include 'pnetcdf.inc'

    integer, intent(in):: nfid 
    integer, intent(in):: varid
    integer(8), intent(in):: start(*)
    integer(8), intent(in):: count(*)
    integer(4), intent(out):: req(*)
    real(4), intent(in):: arr(*)

    integer :: ret      ! NetCDF return code

    ret = nfmpi_iput_vara_real(nfid, varid, start, count, arr,req)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_vara_real_req

  !===============================================================================

  subroutine wrap_put_var_real (nfid, varid, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Put the given real variable to output file.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    real(4), intent(in):: arr(*)

    integer ret      ! NetCDF return code

    !#ifdef SPIO
    !   ret = nf_put_var_real (nfid, varid, arr)
    !#else
    ret = nfmpi_put_var_real (nfid, varid, arr)
    !#endif
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_var_real

  !===============================================================================

  subroutine wrap_put_var_realx (nfid, varid, arr)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Put the given real variable to output file.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    real(8), intent(in):: arr(*)

    integer ret      ! NetCDF return code

    !#ifdef SPIO
    !   ret = nf_put_var_real (nfid, varid, arr)
    !#else
    ret = nfmpi_put_var_double (nfid, varid, arr)
    !#endif
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_var_realx

  !===============================================================================

  subroutine handle_error(ret)
    !------------------------------------------------------------
    !
    ! Purpose:
    !
    ! Handle netCDF errors.
    !
    !------------------------------------------------------------
    !   use abortutils, only : endrun

    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: ret

    write(6,*)nfmpi_strerror(ret)
    print*,'HANDLE_ERROR'
    stop
  end subroutine handle_error

  !===============================================================================

  subroutine wrap_copy_att (nfidin, varidin, attname, nfidout, varidout)
    !------------------------------------------------------------
    ! Added by EZ
    ! 2011-03-03
    ! Purpose:
    !
    ! Copy the given attribute  to variable Out.
    !
    ! This routine violates the convetion that the wrapper codes take an identical
    ! set of arguments as the netcdf library code.  The length of the character
    ! argument is computed inside the wrapper.
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in)      :: nfidin
    integer, intent(in)      :: varidin
    character*(*), intent(in):: attname
    integer, intent(in)      :: nfidout
    integer, intent(in)      :: varidout

    integer ret      ! NetCDF return code

    ret     = nfmpi_copy_att (nfidin, varidin, attname, nfidout, varidout)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_copy_att

  !===============================================================================

  subroutine wrap_redef (nfid)
    !------------------------------------------------------------
    ! Added by EZ
    ! 2011-03-03
    ! Purpose:
    ! 
    ! Put NedCdf data into define mode
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in)      :: nfid
    integer ret      ! NetCDF return code

    ret     = nfmpi_redef (nfid)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_redef

  !===============================================================================

  subroutine wrap_inq_vartype (nfid, varid, xtype)
    !------------------------------------------------------------
    ! Added by EZ
    ! 2011-03-03
    ! Purpose:
    ! 
    ! Inquire Var's type
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in)      :: nfid, varid
    integer, intent(out)     :: xtype
    integer ret      ! NetCDF return code

    ret     = nfmpi_inq_vartype (nfid, varid, xtype)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_inq_vartype

  !===============================================================================

  subroutine wrap_get_var_text (nfid, varid, arr)
    !------------------------------------------------------------
    ! Added by EZ 
    ! 2011-03-06 
    ! Purpose:
    !
    ! Gets the given character variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    character*(*), intent(out):: arr(*)

    integer ret      ! NetCDF return code
    !ret = nfmpi_begin_indep_data(nfid) !switch into and out of independent data mode. 
    ret = nfmpi_get_var_text_all(nfid, varid, arr)
    !ret = nfmpi_end_indep_data(nfid)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VAR_TEXT: error reading varid =', varid
       call handle_error (ret)
    end if
  end subroutine wrap_get_var_text

  subroutine wrap_put_var_text (nfid, varid, text)
    !------------------------------------------------------------
    ! Added by EZ 
    ! 2011-03-06
    ! Purpose:
    !
    ! Gets the given character variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    character*(*), intent(in):: text

    integer ret      ! NetCDF return code

    ret = nfmpi_put_var_text (nfid, varid, text)
    if (ret/=NF_NOERR) call handle_error (ret)
  end subroutine wrap_put_var_text
  !===============================================================================

  subroutine wrap_get_vara_text (nfid, varid, start, count, arr)
    !------------------------------------------------------------
    ! Added by EZ 
    ! 2011-03-06 
    ! Purpose:
    !
    ! Gets the given character variable from a input file
    !
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in):: nfid
    integer, intent(in):: varid
    integer, intent(in):: start(*)
    integer, intent(in):: count(*)
    character*(*), intent(out):: arr(*)

    integer ret      ! NetCDF return code

    ret = nfmpi_get_vara_text (nfid, varid, start, count, arr)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VARA_TEXT: error reading varid =', varid
       call handle_error (ret)
    end if
  end subroutine wrap_get_vara_text

  !===============================================================================

  subroutine wrap_inq_nvars (ncid, nvars)
    !------------------------------------------------------------
    ! Added by EZ 
    ! 2011-03-06 
    ! Purpose:
    !
    !    Inquire how many variables in the ncdata
    !------------------------------------------------------------
    implicit none
    include 'pnetcdf.inc'

    integer, intent(in) :: ncid
    integer, intent(out):: nvars

    integer ret      ! NetCDF return code

    ret = nfmpi_inq_nvars (ncid, nvars)
    if (ret/=NF_NOERR) then
       write(6,*)'WRAP_GET_VARA_TEXT: error reading varid =', ncid
       call handle_error (ret)
    end if
  end subroutine wrap_inq_nvars

  !================================================================================

end module grist_wrap_pf
