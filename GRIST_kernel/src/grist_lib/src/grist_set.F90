
#include "config.h"

module grist_set
  use iso_c_binding
  use grist_utils
  
  type, public :: set
     type(c_ptr) :: ptr = C_NULL_PTR
   contains
     procedure :: insert => set_insert
     procedure :: final  => set_final
     procedure :: init   => set_init
     procedure :: dump   => set_dump
     procedure :: size   => set_size
     procedure :: find   => set_find
     procedure :: sort   => set_sort
  end type set

contains
  
  subroutine set_init(o)
    use iso_c_binding
    implicit none
    interface
       subroutine c_set_init(p) &
            bind(C, name = "c_set_init")
         use iso_c_binding
         type(c_ptr), intent(inout) :: p
       end subroutine
    end interface
    
    class(set), intent(inout) :: o
    
    call c_set_init(o%ptr)
  end subroutine

  subroutine set_insert(o, v)
    implicit none
    interface
       subroutine c_set_insert(p, v) &
            bind(C, name = "c_set_insert")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t :: v
       end subroutine
    end interface
    
    class(set), intent(inout) :: o
    idx_t, intent(in) :: v

    if(.not. c_associated(o%ptr)) call o%init()
    
    call c_set_insert(o%ptr, v)
  end subroutine

  function set_size(o) result(v)
    implicit none
    interface
       subroutine c_set_size(p, v) &
            bind(C, name = "c_set_size")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t, intent(inout) :: v
       end subroutine
    end interface
    
    class(set), intent(in) :: o
    idx_t :: v

    call assert(c_associated(o%ptr), &
         "set must be intialized.")
    
    call c_set_size(o%ptr, v)
  end function


  function set_find(o, v) result(r)
    implicit none
    interface
       subroutine c_set_find(p, v, cr) &
            bind(C, name = "c_set_find")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t, intent(in) :: v
         integer, intent(out) :: cr
       end subroutine
    end interface
    
    class(set), intent(in) :: o
    idx_t :: v
    integer :: cr
    logical :: r

    call assert(c_associated(o%ptr), &
         "set must be intialized.")
    
    call c_set_find(o%ptr, v, cr)

    if(cr == 0) then
       r = .false.
    else
       r = .true.
    end if
  end function
  
  subroutine set_dump(o, v)
    implicit none
    interface
       subroutine c_set_dump(p, v) &
            bind(C, name = "c_set_dump")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t, intent(inout) :: v(*)
       end subroutine
    end interface
    
    class(set), intent(in) :: o
    idx_t, allocatable, intent(inout) :: v(:)

    call assert(c_associated(o%ptr), &
         "set must be intialized.")
    
    if(allocated(v)) then
       call assert(o%size() == size(v), &
            "wrong vector size.")
    else
       allocate(v(o%size()))
    end if
    
    call c_set_dump(o%ptr, v)
  end subroutine
  
  subroutine set_final(o)
    implicit none
    class(set), intent(inout) :: o
    interface
       subroutine c_set_final(p) &
            bind(C, name = "c_set_final")
         use iso_c_binding
         type(c_ptr), intent(inout) :: p
       end subroutine
    end interface

    if(c_associated(o%ptr)) then
       call c_set_final(o%ptr)
    end if
  end subroutine

  subroutine set_sort(o, v)
    implicit none
    interface
       subroutine c_set_sort(p, v) &
            bind(C, name = "c_set_sort")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t, intent(inout) :: v(*)
       end subroutine
    end interface
 
    class(set), intent(in) :: o
    idx_t, allocatable, intent(inout) :: v(:)
 
    call assert(c_associated(o%ptr), &
         "set must be intialized.")
 
    if(allocated(v)) then
       call assert(o%size() == size(v), &
            "wrong vector size.")
    else
       allocate(v(o%size()))
    end if
 
    call c_set_sort(o%ptr, v)
  end subroutine

end module
