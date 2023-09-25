!======================================================
!
!  Created by LiXiaohan on 19/8/14.
!  input file manipulations for physics package.
!
!======================================================

 module grist_physics_iofile


    use grist_constants,                    only: i4, r8
    use grist_handle_error,                 only: endrun
    use grist_mpi

    implicit none
    private

    public :: getfile        ! Get file from archive

 contains

! Purpose : Determine whether file is on local disk.
! . first check current working directory
! . next check full pathname[fulpath] on disk
! . by default, abort if file not found.  Setting optional iflag arg
!   to 1 overrides this behavior, and in that case the optional lexist
!   arg is used to return status of whether the file was found or not.
    subroutine getfile(fulpath, locfn, iflag, lexist)
    ! io
    character(len=*), intent(in)   :: fulpath ! full pathname on local disk
    character(len=*), intent(out)  :: locfn   ! local file name if found in working directory,
                                              ! set to fulpath if not found in working dir.
    integer, optional, intent(in)  :: iflag   ! set iflag=1 to return control to caller if
                                              ! file not found.  default is to abort.
    logical, optional, intent(out) :: lexist  ! When iflag=1 then getfil will return whether the
                                              ! file is found or not.  This flag is set .true.
                                              ! if the file is found, otherwise .false.
    ! local
    integer :: i                              ! loop index
    integer :: klen                           ! length of fulpath character string
    integer :: maxlen                         ! length of locfn input variable
    integer :: ierr                           ! error status
    logical :: lexist_in                      ! true if local file exists
    logical :: abort_on_failure
    ! --------------------------------------------------------------------
 
    abort_on_failure = .true.
    if (present(iflag)) then
       if (iflag==1) abort_on_failure = .false.
    end if
    maxlen = len(locfn)

    ! first check if file is in current working directory.
    ! get local file name from full name: start at end. look for first "/"
 
    klen = len_trim(fulpath)
    i = index(fulpath, '/', back=.true.)

    if ((klen-i) > maxlen) then
       if (abort_on_failure) then
          call endrun('physics_iofile: local filename variable is too short for path length')
       else
          if (mpi_rank()==0) print*, 'physics_iofile: local filename variable is too short for path length',klen-i,maxlen
          if (present(lexist)) lexist = .false.
          return
       end if
    end if

    locfn = fulpath(i+1:klen)
    if (len_trim(locfn) == 0) then
       call endrun ('physics_iofile: local filename has zero length')
    else if (mpi_rank()==0) then
       print*,'physics_iofile: attempting to find local file ', trim(locfn)
    end if
 
    inquire(file=locfn, exist=lexist_in)
    if (present(lexist)) lexist = lexist_in
    if (lexist_in) then
       return
    end if
 
    ! second check for full pathname on disk
 
    if (klen > maxlen) then
       if (abort_on_failure) then
          call endrun('physics_iofile: local filename variable is too short for path length')
       else
          if (mpi_rank()==0) print*, 'physics_iofile: local filename variable is too short for path length',klen,maxlen
          if (present(lexist)) lexist = .false.
          return
       end if
    end if

    locfn = trim(fulpath)
    inquire(file=locfn, exist=lexist_in)
    if (present(lexist)) lexist = lexist_in
    if (lexist_in) then
       return
    else
       if (mpi_rank()==0) print*, 'physics_iofile: all tries to get file have been unsuccessful: ',trim(fulpath)
       if (abort_on_failure) then
          call endrun ('GETFIL: FAILED to get '//trim(fulpath))
       else
          return
       endif
    endif

    end subroutine getfile
 
 end module grist_physics_iofile
