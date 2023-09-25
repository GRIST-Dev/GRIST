
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Additional data besides dycore and tracer that are needed for DTP mode
! Revision history:
!----------------------------------------------------------------------------

  module grist_dtp_vars_module

    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_1d_field, scalar_2d_field
    use grist_constants,    only: rearth, i4, r8, pi
    use grist_nml_module,   only: nlev

    implicit none

    public

!================================================
! Name Convention:
!   scalar_${varname}_at_$(location)_${time}
!   tend_${varname}_at_${location}_${where_from}
!================================================

!================================================
! edge, full level
!================================================

!================================================
! primal cell, full level
!================================================

    type(scalar_2d_field)  :: scalar_thetav_at_pc_full_level_n
    type(scalar_2d_field)  :: scalar_rho_at_pc_full_level_n

!================================================
! dual cell, full level
!================================================
  

!================================================
! edge, face level
!================================================


!================================================
! primal cell, face level
!================================================


!================================================
! surface
!================================================

!================================================
! geometric info
!================================================


  CONTAINS

   subroutine grist_dtp_vars_construct(mesh)
! io
     type(global_domain), intent(in) :: mesh
! local
     integer(i4)     ::  it, ie, iv

!================================================
! edge, full level
!================================================

!================================================
! primal cell, full level
!================================================
    allocate(         scalar_thetav_at_pc_full_level_n%f(nlev,mesh%nv))
    allocate(         scalar_rho_at_pc_full_level_n%f(nlev,mesh%nv))
! pos
    scalar_thetav_at_pc_full_level_n%pos        = 0
    scalar_rho_at_pc_full_level_n%pos           = 0

!================================================
! dual cell, full level
!================================================

!================================================
! edge, face level
!================================================

!================================================
! primal cell, face level
!================================================

!================================================
! surface
!================================================

     return

    end subroutine grist_dtp_vars_construct

    subroutine grist_dtp_vars_destruct()

!================================================
! edge, full level
!================================================

!================================================
! primal cell, full
!================================================
      deallocate(scalar_thetav_at_pc_full_level_n%f)      
      deallocate(scalar_rho_at_pc_full_level_n%f)

!================================================
! dual cell, full
!================================================

!================================================
! edge, face level
!================================================

!================================================
! primal cell, face level
!================================================

!================================================
! surface
!================================================

      return
    end subroutine grist_dtp_vars_destruct

   end module grist_dtp_vars_module
