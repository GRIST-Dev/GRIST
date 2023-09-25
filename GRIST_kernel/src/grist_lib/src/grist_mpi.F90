
module grist_mpi
  include 'mpif.h'

!Rongxy Nov 7,2020 for coupling ->
  public  ::  grist_comm_set
  public  ::  grist_comm

  integer, private :: grist_comm_all
!Rongxy <-

  INTERFACE reduce
    module procedure reduce_int_1
    module procedure reduce_int_n
    module procedure reduce_real4_1
    module procedure reduce_real4_n
    module procedure reduce_real8_1
    module procedure reduce_real8_n
  END INTERFACE

contains

!Rongxy Nov 7,2020 for coupling ->
  subroutine grist_comm_set (comm)
    implicit none
    integer,intent(in) :: comm

    grist_comm_all = comm

  end subroutine grist_comm_set

  function grist_comm () result(comm)

    implicit none
    integer :: comm

    comm=grist_comm_all

  end function
!Rongxy <-

  function mpi_rank(comm) result(res)
    implicit none
    integer, optional :: comm
    integer           :: ierr
    integer           :: res

    if(present(comm)) then
      call MPI_COMM_RANK(comm, res, ierr)
    else
      call MPI_COMM_RANK(grist_comm_all, res, ierr)
    end if

  end function

  function mpi_size(comm) result(res)
    implicit none
    integer, optional :: comm
    integer           :: ierr
    integer           :: res

    if(present(comm)) then
      call MPI_COMM_SIZE(Comm, res, ierr)
    else
      call MPI_COMM_SIZE(grist_comm_all, res, ierr)
    end if

  end function

  subroutine barrier(comm)
    implicit none
    integer, optional :: comm
    integer           :: ierr
    integer           :: res

    if(present(comm)) then
      call MPI_BARRIER(Comm, ierr)
    else
      call MPI_BARRIER(grist_comm_all, ierr)
    end if

  end subroutine barrier

  subroutine reduce_int_1(sbuf, rbuf, op, comm)
    implicit none
    integer,          intent(in)  :: sbuf
    integer,          intent(out) :: rbuf
    character(len=*), intent(in)  :: op
    integer,          optional    :: comm

    ! Local
    integer                       :: comm_, ierr

    comm_ = grist_comm_all
    if(present(comm)) comm_ = comm

    if(op .eq. 'sum' .or. op .eq. 'SUM') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_INTEGER, MPI_SUM, 0, comm_, ierr)
    else if(op .eq. 'min' .or. op .eq. 'MIN') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_INTEGER, MPI_MIN, 0, comm_, ierr)
    else if(op .eq. 'max' .or. op .eq. 'MAX') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_INTEGER, MPI_MAX, 0, comm_, ierr)
    else
      print*,'The reduce option is wrong, please check it!'
    end if
 
  end subroutine reduce_int_1

  subroutine reduce_real4_1(sbuf, rbuf, op, comm)
    implicit none
    real(4),          intent(in)  :: sbuf
    real(4),          intent(out) :: rbuf
    character(len=*), intent(in)  :: op
    integer,          optional    :: comm

    ! Local
    integer                       :: comm_, ierr

    comm_ = grist_comm_all
    if(present(comm)) comm_ = comm

    if(op .eq. 'sum' .or. op .eq. 'SUM') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_REAL, MPI_SUM, 0, comm_, ierr)
    else if(op .eq. 'min' .or. op .eq. 'MIN') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_REAL, MPI_MIN, 0, comm_, ierr)
    else if(op .eq. 'max' .or. op .eq. 'MAX') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_REAL, MPI_MAX, 0, comm_, ierr)
    else
      print*,'The reduce option is wrong, please check it!'
    end if

  end subroutine reduce_real4_1

  subroutine reduce_real8_1(sbuf, rbuf, op, comm)
    implicit none
    real(8),          intent(in)  :: sbuf
    real(8),          intent(out) :: rbuf
    character(len=*), intent(in)  :: op
    integer,          optional    :: comm

    ! Local
    integer                       :: comm_, ierr

    comm_ = grist_comm_all
    if(present(comm)) comm_ = comm

    if(op .eq. 'sum' .or. op .eq. 'SUM') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_REAL8, MPI_SUM, 0, comm_, ierr)
    else if(op .eq. 'min' .or. op .eq. 'MIN') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_REAL8, MPI_MIN, 0, comm_, ierr)
    else if(op .eq. 'max' .or. op .eq. 'MAX') then
      call MPI_REDUCE(sbuf, rbuf, 1, MPI_REAL8, MPI_MAX, 0, comm_, ierr)
    else
      print*,'The reduce option is wrong, please check it!'
    end if

  end subroutine reduce_real8_1

  subroutine reduce_int_n(sbuf, rbuf, n, op, comm)
    implicit none
    integer,          intent(in)  :: sbuf(n)
    integer,          intent(out) :: rbuf(n)
    integer,          intent(in)  :: n
    character(len=*), intent(in)  :: op
    integer,          optional    :: comm

    ! Local
    integer                       :: comm_, ierr

    comm_ = grist_comm_all
    if(present(comm)) comm_ = comm

    if(op .eq. 'sum' .or. op .eq. 'SUM') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_INTEGER, MPI_SUM, 0, comm_, ierr)
    else if(op .eq. 'min' .or. op .eq. 'MIN') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_INTEGER, MPI_MIN, 0, comm_, ierr)
    else if(op .eq. 'max' .or. op .eq. 'MAX') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_INTEGER, MPI_MAX, 0, comm_, ierr)
    else
      print*,'The reduce option is wrong, please check it!'
    end if
 
  end subroutine reduce_int_n

  subroutine reduce_real4_n(sbuf, rbuf, n, op, comm)
    implicit none
    real(4),          intent(in)  :: sbuf(n)
    real(4),          intent(out) :: rbuf(n)
    integer,          intent(in)  :: n
    character(len=*), intent(in)  :: op
    integer,          optional    :: comm

    ! Local
    integer                       :: comm_, ierr

    comm_ = grist_comm_all
    if(present(comm)) comm_ = comm

    if(op .eq. 'sum' .or. op .eq. 'SUM') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_REAL, MPI_SUM, 0, comm_, ierr)
    else if(op .eq. 'min' .or. op .eq. 'MIN') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_REAL, MPI_MIN, 0, comm_, ierr)
    else if(op .eq. 'max' .or. op .eq. 'MAX') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_REAL, MPI_MAX, 0, comm_, ierr)
    else
      print*,'The reduce option is wrong, please check it!'
    end if

  end subroutine reduce_real4_n

  subroutine reduce_real8_n(sbuf, rbuf, n, op, comm)
    implicit none
    real(8),          intent(in)  :: sbuf(n)
    real(8),          intent(out) :: rbuf(n)
    integer,          intent(in)  :: n
    character(len=*), intent(in)  :: op
    integer,          optional    :: comm

    ! Local
    integer                       :: comm_, ierr

    comm_ = grist_comm_all
    if(present(comm)) comm_ = comm

    if(op .eq. 'sum' .or. op .eq. 'SUM') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_REAL8, MPI_SUM, 0, comm_, ierr)
    else if(op .eq. 'min' .or. op .eq. 'MIN') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_REAL8, MPI_MIN, 0, comm_, ierr)
    else if(op .eq. 'max' .or. op .eq. 'MAX') then
      call MPI_REDUCE(sbuf, rbuf, n, MPI_REAL8, MPI_MAX, 0, comm_, ierr)
    else
      print*,'The reduce option is wrong, please check it!'
    end if

  end subroutine reduce_real8_n

 end module grist_mpi
