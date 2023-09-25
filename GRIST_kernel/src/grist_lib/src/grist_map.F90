
#include "config.h"

module grist_map
  use iso_c_binding
  use grist_utils
  
  type, public :: map
     type(c_ptr) :: ptr = C_NULL_PTR
   contains
     procedure :: insert => map_insert
     procedure :: final  => map_final
     procedure :: init   => map_init
     procedure :: size   => map_size
     procedure :: find   => map_find
  end type map

contains

  !> intialize map object
  subroutine map_init(o)
    use iso_c_binding
    implicit none
    interface
       subroutine c_map_init(p) &
            bind(C, name = "c_map_init")
         use iso_c_binding
         type(c_ptr), intent(inout) :: p
       end subroutine
    end interface
    
    class(map), intent(inout) :: o
    
    call c_map_init(o%ptr)
  end subroutine

  !> insert a key-value pair
  subroutine map_insert(o, key, val)
    implicit none
    interface
       subroutine c_map_insert(p, key, val) &
            bind(C, name = "c_map_insert")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t :: key, val
       end subroutine
    end interface
    
    class(map), intent(inout) :: o
    idx_t, intent(in) :: key, val
    
    if(.not. c_associated(o%ptr)) call o%init()
    
    call c_map_insert(o%ptr, key, val)
  end subroutine

  !> get size of a map
  function map_size(o) result(v)
    implicit none
    interface
       subroutine c_map_size(p, v) &
            bind(C, name = "c_map_size")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t, intent(inout) :: v
       end subroutine
    end interface
    
    class(map), intent(in) :: o
    idx_t :: v

    call assert(c_associated(o%ptr), &
         "set must be intialized.")
    
    call c_map_size(o%ptr, v)
  end function

  !> find the value from a key
  function map_find(o, key) result(val)
    implicit none
    interface
       subroutine c_map_find(p, key, val) &
            bind(C, name = "c_map_find")
         use iso_c_binding
         type(c_ptr), intent(in) :: p
         idx_t,   intent(in) :: key
         idx_t, intent(out) :: val
       end subroutine
    end interface
    
    class(map), intent(in) :: o
    idx_t :: key, val

    call assert(c_associated(o%ptr), &
         "set must be intialized.")
    
    call c_map_find(o%ptr, key, val)

  end function

  !> destroy the map object
  subroutine map_final(o)
    implicit none
    class(map), intent(inout) :: o
    interface
       subroutine c_map_final(p) &
            bind(C, name = "c_map_final")
         use iso_c_binding
         type(c_ptr), intent(inout) :: p
       end subroutine
    end interface

    if(c_associated(o%ptr)) then
       call c_map_final(o%ptr)
    end if
  end subroutine
end module
