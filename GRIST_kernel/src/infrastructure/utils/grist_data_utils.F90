
!-----------------------------------------------------------
! Created on Aug 08 2023
! Author: Siyuan Chen
! Version 1.0
! Description: 
!-----------------------------------------------------------

module grist_data_utils

    use grist_data_types,       only: scalar_1d_field,  &
                                      scalar_2d_field,  &
                                      scalar_3d_field

    implicit none

    interface scalar_r8_to_r4
        module procedure scalar_r8_to_r4_1d
        module procedure scalar_r8_to_r4_2d
        module procedure scalar_r8_to_r4_3d
    end interface scalar_r8_to_r4

    interface scalar_r4_to_r8
        module procedure scalar_r4_to_r8_1d
        module procedure scalar_r4_to_r8_2d
        module procedure scalar_r4_to_r8_3d
    end interface scalar_r4_to_r8

contains

    subroutine scalar_r8_to_r4_1d(scalar)

        implicit none

        type(scalar_1d_field), intent(inout)    ::  scalar
        
        if(.not.allocated(scalar%f_r4))then
            print*, "[error!] scalar_r8_to_r4_1d"
            stop
        else
            scalar%f_r4 = scalar%f
        end if

    end subroutine scalar_r8_to_r4_1d

    subroutine scalar_r8_to_r4_2d(scalar)

        implicit none

        type(scalar_2d_field), intent(inout)    ::  scalar
        
        if(.not.allocated(scalar%f_r4))then
            print*, "[error!] scalar_r8_to_r4_2d"
            stop
        else
            scalar%f_r4 = scalar%f
        end if

    end subroutine scalar_r8_to_r4_2d

    subroutine scalar_r8_to_r4_3d(scalar)

        implicit none

        type(scalar_3d_field), intent(inout)    ::  scalar
        
        if(.not.allocated(scalar%f_r4))then
            print*, "[error!] scalar_r8_to_r4_3d"
            stop
        else
            scalar%f_r4 = scalar%f
        end if

    end subroutine scalar_r8_to_r4_3d

    subroutine scalar_r4_to_r8_1d(scalar)

        implicit none

        type(scalar_1d_field), intent(inout)    ::  scalar
        
        if(.not.allocated(scalar%f_r4))then
            print*, "[error!] scalar_r4_to_r8_1d"
            stop
        else
            scalar%f = scalar%f_r4
        end if

    end subroutine scalar_r4_to_r8_1d

    subroutine scalar_r4_to_r8_2d(scalar)

        implicit none

        type(scalar_2d_field), intent(inout)    ::  scalar
        
        if(.not.allocated(scalar%f_r4))then
            print*, "[error!] scalar_r4_to_r8_2d"
            stop
        else
            scalar%f = scalar%f_r4
        end if

    end subroutine scalar_r4_to_r8_2d

    subroutine scalar_r4_to_r8_3d(scalar)

        implicit none

        type(scalar_3d_field), intent(inout)    ::  scalar
        
        if(.not.allocated(scalar%f_r4))then
            print*, "[error!] scalar_r4_to_r8_3d"
            stop
        else
            scalar%f = scalar%f_r4
        end if

    end subroutine scalar_r4_to_r8_3d

end module grist_data_utils