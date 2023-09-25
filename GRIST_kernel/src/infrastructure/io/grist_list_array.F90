
  module grist_list_array

    use grist_constants, only: r8, i4 ,i8

! deal with non-redundant data array
    implicit none

    type array_0d_index_list
        real(r8),  pointer                 :: var
        character(len=100)                 :: varname
        integer(i4)                        :: var_idlist  !  varid in netCDF
        type(array_0d_index_list), pointer :: next!=>null()
    end type array_0d_index_list

    type array_1d_index_list
        real(r8),  pointer                 :: var(:)
        character(len=100)                 :: varname
        character(len=999)                 :: longname
        character(len=999)                 :: units
        integer(i4)                        :: var_idlist  !  varid in netCDF
        integer(i4)                        :: var_pos
        type(array_1d_index_list), pointer :: next !=>null()
    end type array_1d_index_list

    type array_2d_index_list
        real(r8),  pointer                 :: var(:,:)
        character(len=100)                 :: varname
        character(len=999)                 :: longname
        character(len=999)                 :: units
        integer(i4)                        :: var_idlist  !  varid in netCDF
        integer(i4)                        :: var_pos
        integer(i8)                        :: var_level
        type(array_2d_index_list), pointer :: next!=>null()
    end type array_2d_index_list

    type array_3d_index_list
        real(r8),  pointer                 :: var(:,:,:)
        character(len=100)                 :: varname
        character(len=999)                 :: longname
        character(len=999)                 :: units
        integer(i4)                        :: var_idlist  !  varid in netCDF
        integer(i4)                        :: var_pos
        integer(i8)                        :: var_level
        type(array_3d_index_list), pointer :: next!=>null()
    end type array_3d_index_list

    private
    public   :: array_0d_index_list   , &
                array_1d_index_list   , &
                array_2d_index_list   , &
                array_3d_index_list   , &
                insert_data_0d        , &
                insert_data_1d        , &
                insert_data_2d        , &
                insert_var_1d         , &
                insert_var_2d         , &
                insert_var_3d

 contains

  subroutine insert_data_0d(head,new_value,new_name,count)
! io
    type(array_0d_index_list), pointer, intent(inout) :: head
    real(r8),                  target,  intent(in)    :: new_value
    character*(*)                 ,  intent(in)       :: new_name
    integer(i4),                     intent(inout)    :: count 
! local
    type(array_0d_index_list), pointer                :: tail,tmp ! be memorized

      if(.not.associated(head))then
          allocate(head)
          tail         => head
          tail%var     => new_value
          tail%varname =  new_name
          count        =  count+1 
          nullify(tail%next)
      else
          !allocate(tail%next)         ! before, tail%next is null, after this, be associated
          !tail         => tail%next   ! pointer tail points to tail%next
          !tail%var     => new_value   ! value assigned
          !tail%varname =  new_name
          !count        =  count +1 
          !nullify(tail%next)          ! remove this tail's p, optional

        tail=>head
        allocate(tmp)
        tmp%var=>new_value
        tmp%varname =  new_name
        nullify(tmp%next)
        do while(associated(tail%next))
                tail=>tail%next
        enddo
        tmp%next=>tail%next
        tail%next=>tmp

        count=count+1
      end if

    return

  end subroutine  insert_data_0d

  subroutine insert_var_1d(head,new_value,new_name,longname,units,count,pos)
! io
    type(array_1d_index_list), pointer, intent(inout) :: head
	real(r8), allocatable , target,  intent(in)       :: new_value(:)
    character*(*)                 ,  intent(in)       :: new_name
    character*(*)                 ,  intent(in)       :: longname
    character*(*)                 ,  intent(in)       :: units
    integer(i4),                     intent(inout)    :: count
    integer(i4),                     intent(in)       :: pos
! local
    type(array_1d_index_list), pointer                :: tail,tmp ! be memorized


      if(.not.associated(head))then
          allocate(head)
          tail         => head
          tail%var     => new_value
          tail%varname =  new_name
          tail%longname=  longname 
          tail%units   =  units
          tail%var_pos =  pos
          count        =  count+1 
          nullify(tail%next)
      else
          !allocate(tail%next)         ! before, tail%next is null, after this, be associated
          !tail         => tail%next   ! pointer tail points to tail%next
          !tail%var     => new_value  ! value assigned
          !tail%varname =  new_name
          !tail%var_pos =  pos
          !count        =  count +1 
          !nullify(tail%next)          ! remove this tail's p, optional
  
          !add new_value to last of list
 
        tail=>head
        allocate(tmp)
        tmp%var=>new_value
        tmp%varname = new_name
        tmp%longname= longname 
        tmp%units   = units
        tmp%var_pos = pos
        nullify(tmp%next)
        do while(associated(tail%next))
                tail=>tail%next
        enddo
        tmp%next=>tail%next
        tail%next=>tmp

        count=count+1
      end if

    return

  end subroutine  insert_var_1d
  
  subroutine insert_data_1d(head,new_value,new_name,count)
! io
    type(array_1d_index_list), pointer, intent(inout) :: head
    real(r8), allocatable , target,  intent(in)       :: new_value(:)
    character*(*)                 ,  intent(in)       :: new_name
    integer(i4),                     intent(inout)    :: count 
! local
    type(array_1d_index_list), pointer                :: tail,tmp ! be memorized

      if(.not.associated(head))then
          allocate(head)
          tail         => head
          tail%var     => new_value
          tail%varname =  new_name
          count        =  count+1 
          nullify(tail%next)
      else
          !allocate(tail%next)         ! before, tail%next is null, after this, be associated
          !tail         => tail%next   ! pointer tail points to tail%next
          !tail%var     => new_value   ! value assigned
          !tail%varname =  new_name
          !count        =  count +1 
          !nullify(tail%next)          ! remove this tail's p, optional
  
        tail=>head
        allocate(tmp)
        tmp%var=>new_value
        tmp%varname =  new_name
        nullify(tmp%next)
        do while(associated(tail%next))
                tail=>tail%next
        enddo
        tmp%next=>tail%next
        tail%next=>tmp

        count=count+1
      end if

    return

  end subroutine  insert_data_1d

  subroutine insert_data_2d(head,new_value,new_name,count)
! io
    type(array_2d_index_list), pointer, intent(inout) :: head
    real(r8), allocatable , target,  intent(in)       :: new_value(:,:)
    character*(*)                 ,  intent(in)       :: new_name
    integer(i4),                     intent(inout)    :: count 
! local
    type(array_2d_index_list), pointer                :: tail,tmp ! be memorized

      if(.not.associated(head))then
          allocate(head)
          tail         => head
          tail%var     => new_value
          tail%varname =  new_name
          count        =  count+1 
          nullify(tail%next)
      else
          !allocate(tail%next)         ! before, tail%next is null, after this, be associated
          !tail         => tail%next   ! pointer tail points to tail%next
          !tail%var     => new_value   ! value assigned
          !tail%varname =  new_name
          !count        =  count +1 
          !nullify(tail%next)          ! remove this tail's p, optional
  
        tail=>head
        allocate(tmp)
        tmp%var=>new_value
        tmp%varname =  new_name
        nullify(tmp%next)
        do while(associated(tail%next))
                tail=>tail%next
        enddo
        tmp%next=>tail%next
        tail%next=>tmp

        count=count+1
      end if

    return

  end subroutine  insert_data_2d
  
  subroutine insert_var_2d(head,new_value,new_name,longname,units,count,pos)
! io
    type(array_2d_index_list), pointer, intent(inout) :: head
    real(r8), allocatable , target,  intent(in)       :: new_value(:,:)
    character*(*)                 ,  intent(in)       :: new_name
    character*(*)                 ,  intent(in)       :: longname
    character*(*)                 ,  intent(in)       :: units
    integer(i4),                     intent(inout)    :: count 
    integer(i4),                     intent(in)       :: pos 
! local
    type(array_2d_index_list), pointer                :: tail,tmp ! be memorized

      if(.not.associated(head))then
          allocate(head)
          tail         => head
          tail%var     => new_value
          tail%varname =  new_name
          tail%longname=  longname
          tail%units   =  units
          tail%var_pos =  pos
          tail%var_level =  size(new_value(:,1))
          count        =  count+1 
          nullify(tail%next)
      else
          !allocate(tail%next)         ! before, tail%next is null, after this, be associated
          !tail         => tail%next   ! pointer tail points to tail%next
          !tail%var     => new_value   ! value assigned
          !tail%varname =  new_name
          !count        =  count +1 
          !nullify(tail%next)          ! remove this tail's p, optional
 
        tail=>head
        allocate(tmp)
        tmp%var=>new_value
        tmp%varname =  new_name
        tmp%longname=  longname
        tmp%units   =  units
        tmp%var_pos =  pos
        tmp%var_level =  size(new_value(:,1))
        nullify(tmp%next)
        do while(associated(tail%next))
                tail=>tail%next
        enddo
        tmp%next=>tail%next
        tail%next=>tmp

        count=count+1
      end if

    return

  end subroutine  insert_var_2d

  subroutine insert_var_3d(head,new_value,new_name,longname,units,count,pos)
! io
    type(array_3d_index_list), pointer, intent(inout) :: head
    real(r8), allocatable , target,  intent(in)       :: new_value(:,:,:)
    character*(*)                 ,  intent(in)       :: new_name
    character*(*)                 ,  intent(in)       :: longname
    character*(*)                 ,  intent(in)       :: units
    integer(i4),                     intent(inout)    :: count 
    integer(i4),                     intent(in)       :: pos 
! local
    type(array_3d_index_list), pointer                :: tail,tmp ! be memorized

      if(.not.associated(head))then
          allocate(head)
          tail         => head
          tail%var     => new_value
          tail%varname =  new_name
          tail%longname=  longname
          tail%units   =  units
          tail%var_pos =  pos
          tail%var_level =  size(new_value,2)
          count        =  count+1 
          nullify(tail%next)
      else
          !allocate(tail%next)         ! before, tail%next is null, after this, be associated
          !tail         => tail%next   ! pointer tail points to tail%next
          !tail%var     => new_value   ! value assigned
          !tail%varname =  new_name
          !count        =  count +1 
          !nullify(tail%next)          ! remove this tail's p, optional
 
        tail=>head
        allocate(tmp)
        tmp%var=>new_value
        tmp%varname =  new_name
        tmp%longname=  longname
        tmp%units   =  units
        tmp%var_pos =  pos
        tmp%var_level =  size(new_value,2)
        nullify(tmp%next)
        do while(associated(tail%next))
                tail=>tail%next
        enddo
        tmp%next=>tail%next
        tail%next=>tmp

        count=count+1
      end if

    return

  end subroutine  insert_var_3d
 end module grist_list_array
