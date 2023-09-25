
module grist_utils
  interface
     !> sleep for a given time period.
     !> time unit is microsecond.
     subroutine usleep(useconds) bind(C)
       use iso_c_binding 
       implicit none
       integer(c_int32_t), value :: useconds
     end subroutine
  end interface

  real(8) :: time_start, time_finish
contains

  !> convert a fortran string to C string
  function string_f2c(f_string) result(c_string)
    use iso_c_binding
    character(len=*):: f_string
    CHARACTER(LEN=LEN_TRIM(f_string)+1,KIND=C_CHAR) :: c_string

    c_string = trim(f_string) // C_NULL_CHAR
  end function


  function i2s(i) result(res)
    character(:),allocatable :: res
    integer,intent(in) :: i
    character(range(i)+2) :: tmp
    write(tmp,'(i0)') i
    res = trim(tmp)
  end function
  
  subroutine tic()
    implicit none
    
    call cpu_time(time_start)
  end subroutine

  subroutine toc()
    implicit none
    
    call cpu_time(time_finish)
    print*, "time elapsed is ", &
         time_finish-time_start, " secs."
  end subroutine
  
  subroutine assert(condition, msg)
    implicit none
    logical :: condition
    integer :: linenum
    character(len=*) :: msg

    if(.not. condition) then
       write(*,"(A)", advance="no") &
            "Error: assertation failed"
       write(*,*) "Reason: ", msg
       call abort()
    endif
  end subroutine assert

end module grist_utils

