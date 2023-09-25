 module grist_list

    implicit none

    type index_list
        integer  :: gidx        ! global index
        type(index_list), pointer :: next=>null()
    end type index_list

    public :: index_list   , &
              check_and_add
 contains

  subroutine check_and_add(head,new_value,n)
!
! add an index, check whether it has been
! included in the current list,
! if yes, pass, if no, add it
!
    type(index_list), pointer, intent(inout) :: head
    integer,                   intent(in)    :: new_value
    integer,                   intent(inout) :: n
!
! local
!
    type(index_list), pointer, save          :: tail ! be memorized
    type(index_list), pointer                :: ptr

    if(.not.associated(head))then
       allocate(head)
       tail=>head
       tail%gidx = new_value
       nullify(tail%next)
       n = n + 1
    else
       ptr=>head
       do while (associated(ptr))
          if(new_value.eq.ptr%gidx) return
          ptr=>ptr%next
       end do
       allocate(tail%next)      ! before, tail%next is null, after this, be associated
       tail=>tail%next          ! pointer tail points to tail%next
       tail%gidx = new_value    ! value assigned
       nullify(tail%next)       ! remove this tail's p, optional
       n = n + 1
    end if

  end subroutine check_and_add

end module grist_list