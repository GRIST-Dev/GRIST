
!=====================================================================
!  Created by zhangyi on 16/8/5.
!
!    File I/O module tailored for GRIST. It deals with the general 
! conditions, and must be called by an explicit suit of statements.
!
! How to use:
!            1) wrap_output_init_0d
!            2) wrap_add_field_0d
!            3) wrap_output_0d
!            4) wrap_output_clean_0d
!=====================================================================

module grist_fileio_0d_module_gcm
  use grist_lib
  use grist_constants,      only: i4, i8, r8
  use grist_domain_types,   only: global_domain
  use grist_handle_error,   only: endrun
  use grist_wrap_pf,        only: wrap_open,           &
                                  wrap_create,         &
                                  wrap_close,          &
                                  wrap_enddef,         &
                                  wrap_inq_dimid,      &
                                  wrap_inq_varid,      &
                                  wrap_inq_vartype,    &
                                  wrap_def_var,        &
                                  wrap_def_dim,        &
                                  wrap_put_att_text,   &
                                  wrap_put_var_realx,  &
                                  wrap_put_var_real,   &
                                  wrap_put_var_int,    &
                                  wrap_put_var_text,   &
                                  wrap_get_var_realx,  &
                                  wrap_get_var_real,   &
                                  wrap_get_var_int,    &
                                  wrap_get_var_text

  implicit none
  include 'pnetcdf.inc'

  private

  public :: wrap_output_init_0d  ,&
       wrap_add_field_0d    ,&
       wrap_output_0d       ,&
       wrap_output_clean_0d ,&
       wrap_read_0d         ,&
       wrap_bcast_0d

  interface wrap_add_field_0d
     module procedure wrap_add_field_0d_real
     module procedure wrap_add_field_0d_int
     module procedure wrap_add_field_0d_char
  end interface wrap_add_field_0d

  interface wrap_read_0d
     module procedure wrap_read_0d_real
     module procedure wrap_read_0d_int
     module procedure wrap_read_0d_char
  end interface wrap_read_0d
  
  interface wrap_bcast_0d
     module procedure wrap_bcast_0d_real
     module procedure wrap_bcast_0d_int
     module procedure wrap_bcast_0d_char
  end interface wrap_bcast_0d

  type scalar_real_template
     real(r8), allocatable             :: f(:) ! the last dim is maxnum
     character(len=100), allocatable   :: varname(:)
  end type scalar_real_template

  type scalar_int_template
     integer(i4), allocatable          :: f(:) ! the last dim is maxnum
     character(len=100), allocatable   :: varname(:)
  end type scalar_int_template

  type scalar_char_template
     character(len=128), allocatable   :: f(:) ! the last dim is maxnum
     character(len=100), allocatable   :: varname(:)
  end type scalar_char_template

  type(scalar_real_template)             :: scalar_real ! maxnum
  type(scalar_int_template)              :: scalar_int  ! maxnum
  type(scalar_char_template)             :: scalar_char ! maxnum

  integer(i4)                            :: number_of_scalar_real
  integer(i4)                            :: number_of_scalar_int
  integer(i4)                            :: number_of_scalar_char

  integer(i4)                            :: maxnum

contains

  subroutine wrap_output_init_0d()

    ! maximum number of vars of each type that we can store
    maxnum                   = 20  
    number_of_scalar_real    = 0
    number_of_scalar_int     = 0
    number_of_scalar_char    = 0

    if(.not.allocated(scalar_real%f))       allocate(scalar_real%f(maxnum))
    if(.not.allocated(scalar_int%f))        allocate(scalar_int%f(maxnum))
    if(.not.allocated(scalar_char%f))       allocate(scalar_char%f(maxnum))

    if(.not.allocated(scalar_real%varname)) allocate(scalar_real%varname(maxnum))
    if(.not.allocated(scalar_int%varname))  allocate(scalar_int%varname(maxnum))
    if(.not.allocated(scalar_char%varname)) allocate(scalar_char%varname(maxnum))

    return
  end subroutine wrap_output_init_0d

  subroutine wrap_add_field_0d_real(scalar,scalar_name)
    real(r8),        intent(in)   :: scalar
    character*(*),   intent(in)   :: scalar_name
    call insert_scalar_in(scalar_name, scalar_real_in=scalar)
    return
  end subroutine wrap_add_field_0d_real

  subroutine wrap_add_field_0d_int(scalar,scalar_name)
    integer(i4),     intent(in)   :: scalar
    character*(*),   intent(in)   :: scalar_name
    call insert_scalar_in(scalar_name, scalar_int_in=scalar)
    return
  end subroutine wrap_add_field_0d_int

  subroutine wrap_add_field_0d_char(scalar,scalar_name)
    character*(*),   intent(in)   :: scalar
    character*(*),   intent(in)   :: scalar_name
    call insert_scalar_in(scalar_name, scalar_char_in=scalar)
    return
  end subroutine wrap_add_field_0d_char

  subroutine wrap_output_0d(comm, itimestep, dtime, outdir, filename)
    implicit none
    !io
    integer,              intent(in)  :: comm
    integer(i4),          intent(in)  :: itimestep
    real(r8),             intent(in)  :: dtime
    character*(*),        intent(in)  :: outdir
    character*(*),        intent(in)  :: filename
    ! local
    integer(i4), parameter            :: omode =0
    integer(i4)                       :: ncid
    integer(i4)                       :: dim_id
    integer(i4)                       :: ret

    integer(i4),allocatable           :: var_real_idlist(:)  ! nt*6 var id list
    integer(i4),allocatable           :: var_int_idlist(:)  ! ne*6 var id list
    integer(i4),allocatable           :: var_char_idlist(:)  ! nv*6 var id list

    integer(i8)                       :: dim_one   = 1
    integer(i8)                       :: dim_two   = 2
    integer(i8)                       :: dim_three = 3

    !================================================
    !                   set data
    !================================================

    print*,"----------------------------------------------------------"
    print*,"    number_of_scalar_real= ", number_of_scalar_real
    print*,"    number_of_scalar_int=  ", number_of_scalar_int
    print*,"    number_of_scalar_char= ", number_of_scalar_char
    print*,"----------------------------------------------------------"

    !================================================
    !                 create file
    !================================================

    call wrap_create (comm, trim(outdir)//trim(filename),omode, ncid)
    print*,trim(outdir)//trim(filename)

    !========================================================
    ! 1st: define dim, need dim_id
    !========================================================

    call wrap_def_dim(ncid,"scalar_dim",dim_one, dim_id)

    !========================================================
    ! 2nd: define var, var_xxx_idlist
    !========================================================

    call def_var_pnetcdf(number_of_scalar_real, ncid, dim_id, var_real_idlist,'real')
    call def_var_pnetcdf(number_of_scalar_int , ncid, dim_id, var_int_idlist ,'int')
    call def_var_pnetcdf(number_of_scalar_char, ncid, dim_id, var_char_idlist,'char')

    !================================================
    ! 3rd: put att
    !================================================

    call put_att_pnetcdf(number_of_scalar_real, ncid, var_real_idlist, 'real')
    call put_att_pnetcdf(number_of_scalar_int , ncid, var_int_idlist,  'int')
    call put_att_pnetcdf(number_of_scalar_char, ncid, var_char_idlist, 'char')

    call wrap_enddef(ncid)
    

    !================================================
    ! 4th, put var
    !================================================

    call put_var_pnetcdf(number_of_scalar_real, ncid, var_real_idlist, 'real')
    call put_var_pnetcdf(number_of_scalar_int , ncid, var_int_idlist , 'int')
    call put_var_pnetcdf(number_of_scalar_char, ncid, var_char_idlist, 'char')

    call wrap_close(ncid)

    deallocate(var_real_idlist)
    deallocate(var_int_idlist)
    deallocate(var_char_idlist)
    !================================================
    return
  end subroutine wrap_output_0d

  subroutine wrap_output_clean_0d()

    if(allocated(scalar_real%f))      deallocate(scalar_real%f)
    if(allocated(scalar_int%f))       deallocate(scalar_int%f)
    if(allocated(scalar_char%f))      deallocate(scalar_char%f)

    number_of_scalar_real = 0
    number_of_scalar_int  = 0
    number_of_scalar_char = 0

    return
  end subroutine wrap_output_clean_0d

  subroutine wrap_read_0d_real(comm, outdir,filename,varname, varout)
    
    ! io
    integer,               intent(in)    :: comm
    character*(*),         intent(in)    :: outdir
    character*(*),         intent(in)    :: filename
    character*(*),         intent(in)    :: varname
    real(r8),              intent(inout) :: varout
    ! local
    integer(i4), parameter               :: omode  =  0
    integer(i4)                          :: nfid
    integer(i4)                          :: varid
    integer(i4)                          :: xtype
    real(8)                              :: temp_dp(1)
    real(4)                              :: temp_sp(1)

    call wrap_open(comm, trim(outdir)//trim(filename), omode, nfid)
    call wrap_inq_varid(nfid, trim(varname), varid)
    call wrap_inq_vartype(nfid, varid, xtype)
  
    if (xtype .eq. NF_DOUBLE) then
       call wrap_get_var_realx(nfid, varid, temp_dp)
       varout = temp_dp(1)
    elseif (xtype .eq. NF_FLOAT) then
       call wrap_get_var_real(nfid, varid, temp_sp)
       varout = temp_sp(1)
    end if
    call wrap_close(nfid)

    return
  end subroutine wrap_read_0d_real

  subroutine wrap_read_0d_int(comm, outdir,filename,varname,varout)
    !
    ! io
    !
    integer,               intent(in)    :: comm
    character*(*),         intent(in)    :: outdir
    character*(*),         intent(in)    :: filename
    character*(*),         intent(in)    :: varname
    integer(i8),           intent(inout) :: varout
    !
    ! local
    !
    integer(i4), parameter               :: omode  =  0
    integer(i4)                          :: ncid
    integer(i4)                          :: varid
    integer(i8),allocatable              :: vartmp(:)

    allocate(vartmp(1))

    call wrap_open(comm, trim(outdir)//trim(filename), omode, ncid)
    call wrap_inq_varid(ncid, trim(varname), varid)
    call wrap_get_var_int(ncid, varid, vartmp)
    call wrap_close(ncid)
    
    varout = vartmp(1)
    deallocate(vartmp)

    return
  end subroutine wrap_read_0d_int

  subroutine wrap_read_0d_char(comm, outdir, filename, varname, varout)
    ! io
    integer,               intent(in)    :: comm
    character*(*),         intent(in)    :: outdir
    character*(*),         intent(in)    :: filename
    character*(*),         intent(in)    :: varname
    character*(*),         intent(inout) :: varout
    ! local
    integer(i4), parameter               :: omode = 0
    integer(i4)                          :: ncid
    integer(i4)                          :: varid
    character(len=32),allocatable        :: vartmp(:)

    allocate(vartmp(1))
    vartmp(1) = ""

    call wrap_open(comm, trim(outdir)//trim(filename),omode,ncid)
    call wrap_inq_varid(ncid, trim(varname), varid)
    call wrap_get_var_text(ncid, varid, vartmp)
    call wrap_close(ncid)
    varout = trim(vartmp(1))

    deallocate(vartmp)
    return
  end subroutine wrap_read_0d_char
  
  subroutine wrap_bcast_0d_real(comm, irank, varout)
    !
    ! io
    !
    integer,               intent(in)    :: comm
    integer,               intent(in)    :: irank
    real(r8),              intent(inout) :: varout
    !
    ! local
    !
    integer                              ::err
    real(r8), allocatable                :: buf(:)

    allocate(buf(1))
    
    buf(1) = varout
    if(r8 .eq. 8) then
      call MPI_Bcast(buf, 1, MPI_REAL8, 0, comm, err)
    elseif(r8 .eq. 4) then
      call MPI_Bcast(buf, 1, MPI_REAL, 0, comm, err)
    end if
    varout = buf(1)
    
    deallocate(buf)
    return
  end subroutine wrap_bcast_0d_real

  subroutine wrap_bcast_0d_int(comm, irank, varout)
    !
    ! io
    !
    integer,               intent(in)    :: comm
    integer,               intent(in)    :: irank
    integer(i8),           intent(inout) :: varout
    !
    ! local
    !
    integer                              ::err
    integer(i8),           allocatable   :: buf(:)
    allocate(buf(1))
    
    buf(1) = varout
    call MPI_Bcast(buf, 1, MPI_INTEGER8, irank, comm, err)
    varout = buf(1) 
    
    deallocate(buf)
    return
  end subroutine wrap_bcast_0d_int

    subroutine wrap_bcast_0d_char(comm, irank, varout)
    !
    ! io
    !
    integer,               intent(in)    :: comm
    integer,               intent(in)    :: irank
    character*(*),         intent(inout) :: varout
    !
    ! local
    !
    integer                              ::err
    character(len=32),     allocatable   :: buf(:)
    allocate(buf(1))
    
    buf(1) = varout
    call MPI_Bcast(buf, 32, MPI_CHARACTER, 0, comm, err)
    varout = buf(1) 
    
    deallocate(buf)
    return
  end subroutine wrap_bcast_0d_char


  !====================================================================
  !
  !                         Private routines below
  !
  !====================================================================

  subroutine check_leng(num_of_scalar,maxnum,address)
    ! io
    integer(i4),   intent(in)   :: num_of_scalar
    integer(i4),   intent(in)   :: maxnum
    character*(*), intent(in)   :: address

    if(num_of_scalar.gt.maxnum)then
       call endrun("number_of_scalar in "//trim(address)//" larger than maxnum_2d, stop and increase")
    end if

    return
  end subroutine check_leng

  subroutine insert_scalar_in(scalar_name_in, scalar_real_in, scalar_int_in, scalar_char_in)
    ! io
    character*(*),               intent(in)     :: scalar_name_in
    real(r8),        optional,   intent(in)     :: scalar_real_in
    integer(i4),     optional,   intent(in)     :: scalar_int_in
    character*(*),   optional,   intent(in)     :: scalar_char_in

    if(present(scalar_real_in))then

       number_of_scalar_real  = number_of_scalar_real + 1
       call check_leng(number_of_scalar_real, maxnum, 'real')
       scalar_real%f(number_of_scalar_real) = scalar_real_in
       scalar_real%varname(number_of_scalar_real) = scalar_name_in  

    end if

    if(present(scalar_int_in))then

       number_of_scalar_int  = number_of_scalar_int + 1
       call check_leng(number_of_scalar_int, maxnum, 'int')
       scalar_int%f(number_of_scalar_int)        = scalar_int_in
       scalar_int%varname(number_of_scalar_int)  = scalar_name_in  

    end if

    if(present(scalar_char_in))then

       number_of_scalar_char  = number_of_scalar_char + 1
       call check_leng(number_of_scalar_char, maxnum, 'char')
       scalar_char%f(number_of_scalar_char)       = scalar_char_in
       scalar_char%varname(number_of_scalar_char) = scalar_name_in

    end if

    return
  end subroutine insert_scalar_in

  subroutine def_var_pnetcdf(num_of_scalar,ncid,dim_id,var_idlist,flag)
    ! io
    integer(i4),                 intent(in)    :: num_of_scalar
    integer(i4),                 intent(in)    :: ncid
    integer(i4),                 intent(in)    :: dim_id
    integer(i4), allocatable,    intent(inout) :: var_idlist(:)
    character*(*),               intent(in)    :: flag
    ! local
    integer(i4)                                :: vdims(1)
    integer(i4)                                :: kk

    select case (trim(flag))
    case('real')
       if(num_of_scalar.ne.0)then
          allocate(var_idlist(num_of_scalar))
          vdims(1) = dim_id
          do kk = 1, num_of_scalar
             if(r8 .eq. 8) then
               call wrap_def_var(ncid,trim(scalar_real%varname(kk)),NF_DOUBLE, 1, vdims(1), var_idlist(kk))
             elseif(r8 .eq. 4) then
               call wrap_def_var(ncid,trim(scalar_real%varname(kk)),NF_REAL, 1, vdims(1), var_idlist(kk))
             end if
          end do
       end if
    case('int')
       if(num_of_scalar.ne.0)then
          allocate(var_idlist(num_of_scalar))
          vdims(1) = dim_id
          do kk = 1, num_of_scalar
             call wrap_def_var(ncid,trim(scalar_int%varname(kk)),NF_INT,     1, vdims(1), var_idlist(kk))
          end do
       end if
    case('char')
       if(num_of_scalar.ne.0)then
          allocate(var_idlist(num_of_scalar))
          vdims(1) = dim_id
          do kk = 1, num_of_scalar
             call wrap_def_var(ncid,trim(scalar_char%varname(kk)),NF_CHAR,   1, vdims(1), var_idlist(kk))
          end do
       end if
    case default
       call endrun('wrong in def_var_netcdf')
    end select

    return
  end subroutine def_var_pnetcdf

  subroutine put_att_pnetcdf(num_of_scalar, ncid, var_idlist, flag)
    ! io
    integer(i4),        intent(in)  :: num_of_scalar
    integer(i4),        intent(in)  :: ncid
    integer(i4),        intent(in)  :: var_idlist(*)
    character*(*),      intent(in)  :: flag
    ! local
    integer(i4)                     :: kk

    select case (trim(flag))

    case('real')
       if(num_of_scalar.ne.0)then
          do kk = 1, num_of_scalar
             call wrap_put_att_text(ncid, var_idlist(kk),'long_name', trim(scalar_real%varname(kk)))
          end do
       end if
    case('int')
       if(num_of_scalar.ne.0)then
          do kk = 1, num_of_scalar
             call wrap_put_att_text(ncid, var_idlist(kk),'long_name', trim(scalar_int%varname(kk)))
          end do
       end if
    case('char')
       if(num_of_scalar.ne.0)then
          do kk = 1, num_of_scalar
             call wrap_put_att_text(ncid, var_idlist(kk),'long_name', trim(scalar_char%varname(kk)))
          end do
       end if
    case default
       call endrun('wrong in put_att_netcdf')
    end select

    return
  end subroutine put_att_pnetcdf

  subroutine put_var_pnetcdf(num_of_scalar, ncid, var_idlist, flag)
    ! io
    integer(i4),        intent(in)  :: num_of_scalar
    integer(i4),        intent(in)  :: ncid
    integer(i4),        intent(in)  :: var_idlist(*)
    character*(*),      intent(in)  :: flag
    ! local
    integer(i4)                     :: kk
    real(4)                         :: temp_sp(1)
    real(8)                         :: temp_dp(1)

    select case (trim(flag))

    case('real')
       if(num_of_scalar.ne.0)then
          do kk = 1, num_of_scalar
             if (r8 .eq. 8) then
               call wrap_put_var_realx(ncid, var_idlist(kk), temp_dp)
               scalar_real%f(kk) = temp_dp(1)
             elseif (r8 .eq. 4) then
               call wrap_put_var_real(ncid, var_idlist(kk), temp_sp)
               scalar_real%f(kk) = temp_sp(1)
             end if
          end do
       end if
    case('int')
       if(num_of_scalar.ne.0)then
          do kk = 1, num_of_scalar
             call wrap_put_var_int(ncid, var_idlist(kk), scalar_int%f(kk))
          end do
       end if
    case('char')
       if(num_of_scalar.ne.0)then
          do kk = 1, num_of_scalar
             call wrap_put_var_text(ncid, var_idlist(kk), scalar_char%f(kk))
          end do
       end if
    case default
       call endrun('wrong in put_var_netcdf')
    end select

    return
  end subroutine put_var_pnetcdf

end module grist_fileio_0d_module_gcm
