
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Definitions of data types
! Revision history:
!----------------------------------------------------------------------------

 module grist_data_types

  use grist_constants,     only: i4, r8, zero
  use grist_constants,     only: r4 => ns
  use grist_element_types_icosh, only: vector

  implicit none

  public :: scalar_1d_field       ,&
            scalar_2d_field       ,&
            scalar_3d_field       ,&
            exchange_field_list_1d,&
            exchange_field_list_2d,&
            exchange_field_list_3d,&
            data1d_temp           ,&
            data2d_temp           ,&
            data3d_temp           ,&
            vector_type           ,&
            wrap_allocate_data1d  ,&
            wrap_allocate_data2d  ,&
            wrap_allocate_data3d  ,&
            wrap_deallocate_data1d,&
            wrap_deallocate_data2d,&
            wrap_deallocate_data3d

  type data1d_temp
      real(r8), allocatable    :: f(:)
  end type data1d_temp

  type data1d_temp_sp
      real(4), allocatable     :: f(:)
  end type data1d_temp_sp

  type data1d_temp_dp
      real(8), allocatable     :: f(:)
  end type data1d_temp_dp

  type data2d_temp
      real(r8), allocatable    :: f(:,:)
  end type data2d_temp

  type data2d_temp_sp
      real(4), allocatable     :: f(:,:)
  end type data2d_temp_sp

  type data2d_temp_dp
      real(8), allocatable     :: f(:,:)
  end type data2d_temp_dp

  type data3d_temp
      real(r8), allocatable    :: f(:,:,:)
  end type data3d_temp

  type data3d_temp_sp
      real(4), allocatable    :: f(:,:,:)
  end type data3d_temp_sp

  type data3d_temp_dp
      real(8), allocatable    :: f(:,:,:)
  end type data3d_temp_dp

  type scalar_1d_field
     real(r8), allocatable      :: f(:)     ! field values array, ordered in the same sequence as index
     real(r4), allocatable      :: f_r4(:)
     integer(i4)                :: pos = -1 ! 0: v location; 1: t location; 6: e location
  end type scalar_1d_field

  type scalar_2d_field
     real(r8), allocatable      :: f(:,:)   ! field values array, ordered in the same sequence as index
     real(r4), allocatable      :: f_r4(:,:)
     integer(i4)                :: pos = -1 ! Position of the values relative to a mesh
  end type scalar_2d_field

  type scalar_3d_field
     real(r8), allocatable      :: f(:,:,:) ! field values array, ordered in the same sequence as index
     real(r4), allocatable      :: f_r4(:,:,:)
     integer(i4)                :: pos = -1 ! Position of the values relative to a mesh
  end type scalar_3d_field

  type vector_type
     type(vector), allocatable  :: p(:) ! vector in cartesian coord on 'pos'
     integer(i4)                :: pos  ! Position same as scalar field 
     integer(i4)                :: n    ! Number of values on the mesh
     character(len=32)          :: name
  end type vector_type
  
  type exchange_field_list_1d
      type(scalar_1d_field),pointer         :: field_data
      integer(i4)                           :: dim1
      real(r8), allocatable                 :: halo(:)
      real(r8), allocatable                 :: bdry(:)
      type(exchange_field_list_1d), pointer :: next => null()
  end type exchange_field_list_1d
  
  type exchange_field_list_2d
      type(scalar_2d_field),pointer         :: field_data
      integer(i4)                           :: dim1
      integer(i4)                           :: dim2
      type(exchange_field_list_2d), pointer :: next => null()
  end type exchange_field_list_2d

  type exchange_field_list_3d
      type(scalar_3d_field),pointer         :: field_data
      integer(i4)                           :: dim1
      integer(i4)                           :: dim2
      integer(i4)                           :: dim3
      type(exchange_field_list_3d), pointer :: next => null()
  end type exchange_field_list_3d

CONTAINS

  subroutine wrap_allocate_data1d(ncell,var)
      integer(i4),           intent(in)    :: ncell
      type(scalar_1d_field), intent(inout) :: var
      if(.not.allocated(var%f)) allocate(var%f(ncell))
      var%f    = zero
      var%pos  = 0
      return
  end subroutine wrap_allocate_data1d

  subroutine wrap_allocate_data2d(ncell,nLevel,var)
      integer(i4),           intent(in)    :: ncell
      integer(i4),           intent(in)    :: nLevel
      type(scalar_2d_field), intent(inout) :: var
      if(.not.allocated(var%f)) allocate(var%f(nLevel,ncell))
      var%f    = zero
      var%pos  = 0
      return
  end subroutine wrap_allocate_data2d

  subroutine wrap_allocate_data3d(ncell,nLevel,ntracer,var)
      integer(i4),           intent(in)    :: ncell
      integer(i4),           intent(in)    :: nLevel
      integer(i4),           intent(in)    :: ntracer
      type(scalar_3d_field), intent(inout) :: var
      if(.not.allocated(var%f)) allocate(var%f(ntracer,nLevel,ncell))
      var%f    = zero
      var%pos  = 0
      return
  end subroutine wrap_allocate_data3d

  subroutine wrap_deallocate_data1d(var)
      type(scalar_1d_field), intent(inout) :: var
      if(allocated(var%f)) deallocate(var%f)
      return
  end subroutine wrap_deallocate_data1d

  subroutine wrap_deallocate_data2d(var)
      type(scalar_2d_field), intent(inout) :: var
      if(allocated(var%f)) deallocate(var%f)
      return
  end subroutine wrap_deallocate_data2d

  subroutine wrap_deallocate_data3d(var)
      type(scalar_3d_field), intent(inout) :: var
      if(allocated(var%f)) deallocate(var%f)
      return
  end subroutine wrap_deallocate_data3d

 end module grist_data_types
