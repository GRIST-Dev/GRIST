 module grist_handle_error

   implicit none

   private

   public :: endrun, &
             showinfo, &
             delete_file_sub
 
   character(len=200), public :: grist_message

   contains

   subroutine endrun(info, number)
!
! io
!
   character*(*), optional, intent(in)   :: info
   integer, optional, intent(in)         :: number

     if(present(info).and.present(number))then
        print*,"--------------------------------------------"
        print*,"                                            "
        print*,"          ENDRUN: ",trim(info), number
        print*,"                                            "
        print*,"--------------------------------------------"
     end if

     stop
   end subroutine endrun

   subroutine showinfo(info, number)
!
! io
!
   character*(*), optional, intent(in)   :: info
   integer, optional, intent(in)         :: number

     if(present(info).and.present(number))then
        print*,"--------------------------------------------"
        print*,"                                            "
        print*,"          SHOWINFO: ",trim(info), number
        print*,"                                            "
        print*,"--------------------------------------------"
     end if

   end subroutine showinfo

  subroutine delete_file_sub(filename)
    character(len=*), intent(in) :: filename
    logical :: deleted
    logical :: file_exists
 
    deleted = .false.
    do while (.not. deleted)
        inquire(file=filename, exist=file_exists)
        if (file_exists) then
            open(unit=10, file=filename, status='old')
            close(unit=10, status='delete')
            deleted = .true.
        else
            print *, 'File does not exist.'
            return
        end if
    end do
  end subroutine delete_file_sub

 end module grist_handle_error