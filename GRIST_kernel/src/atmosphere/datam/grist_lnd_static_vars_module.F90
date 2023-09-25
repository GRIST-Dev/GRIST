 
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Land data Module
! Revision history:
!  data directly used by other components, provided by static data
!----------------------------------------------------------------------------
  module grist_lnd_static_vars_module

    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_1d_field, scalar_2d_field
    use grist_constants,    only: rearth, i4, r8, pi
    use grist_nml_module,   only: nlev, working_mode, test_real_case
    use grist_datam_static_data_module, only:   staticData_landfrac_at_pc_surface, &
                                                staticData_phis_at_pc_surface, & !--cheyz add
                                                staticData_luIndex_at_pc_surface, &
                                                staticData_albedo_at_pc_surface, &
                                                staticData_soilTypetop_at_pc_surface, &
                                                staticData_greenfrac_at_pc_surface, &
                                                staticData_snoalb_at_pc_surface

    implicit none

    public

!================================================
! Name Convention:
!   scalar_static_${varname}_at_$(location)
!================================================
!================================================
! primal cell, full level
!================================================

   type(scalar_2d_field) :: scalar_static_albedo_at_pc
   type(scalar_2d_field) :: scalar_static_greenfrac_at_pc
   type(scalar_1d_field) :: scalar_static_lu_index_at_pc
   type(scalar_1d_field) :: scalar_static_landmask_at_pc
   type(scalar_1d_field) :: scalar_static_soilcat_top_at_pc
   type(scalar_1d_field) :: scalar_static_snoalb_at_pc
   type(scalar_1d_field) :: scalar_static_ter_at_pc
   !type(scalar_1d_field) :: scalar_static_soiltemp_at_pc
   !type(scalar_1d_field) :: scalar_static_shdmin_at_pc
   !type(scalar_1d_field) :: scalar_static_shdmax_at_pc

  CONTAINS

   subroutine grist_lnd_static_vars_construct(mesh)
! io
     type(global_domain), intent(in) :: mesh

!================================================
! primal cell, full level
!================================================

    if(.not.allocated(scalar_static_albedo_at_pc%f))     allocate(scalar_static_albedo_at_pc%f(12,mesh%nv))    ; scalar_static_albedo_at_pc%pos     = 0
    if(.not.allocated(scalar_static_greenfrac_at_pc%f))  allocate(scalar_static_greenfrac_at_pc%f(12, mesh%nv)); scalar_static_greenfrac_at_pc%pos  = 0
    if(.not.allocated(scalar_static_lu_index_at_pc%f))   allocate(scalar_static_lu_index_at_pc%f(mesh%nv))     ; scalar_static_lu_index_at_pc%pos   = 0
    if(.not.allocated(scalar_static_landmask_at_pc%f))   allocate(scalar_static_landmask_at_pc%f(mesh%nv))     ; scalar_static_landmask_at_pc%pos   = 0
    if(.not.allocated(scalar_static_soilcat_top_at_pc%f))allocate(scalar_static_soilcat_top_at_pc%f(mesh%nv))  ; scalar_static_soilcat_top_at_pc%pos= 0
    if(.not.allocated(scalar_static_snoalb_at_pc%f))     allocate(scalar_static_snoalb_at_pc%f(mesh%nv))       ; scalar_static_snoalb_at_pc%pos     = 0
    if(.not.allocated(scalar_static_ter_at_pc%f))        allocate(scalar_static_ter_at_pc%f(mesh%nv))          ; scalar_static_ter_at_pc%pos        = 0
    !if(.not.allocated(scalar_static_soiltemp_at_pc%f))   allocate(scalar_static_soiltemp_at_pc%f(mesh%nv))     ; scalar_static_soiltemp_at_pc%pos   = 0
    !if(.not.allocated(scalar_static_shdmin_at_pc%f))     allocate(scalar_static_shdmin_at_pc%f(mesh%nv))       ; scalar_static_shdmin_at_pc%pos     = 0
    !if(.not.allocated(scalar_static_shdmax_at_pc%f))     allocate(scalar_static_shdmax_at_pc%f(mesh%nv))       ; scalar_static_shdmax_at_pc%pos     = 0
 
! change this in future
! staticData are for 1:nv_full, land_vars are as physics for 1:nv_halo(1)

    scalar_static_landmask_at_pc%f  = 0
    if(test_real_case) scalar_static_landmask_at_pc%f(1:mesh%nv) = staticData_landfrac_at_pc_surface%f(1:mesh%nv)

#ifdef USE_NOAHMP
    if(test_real_case)then
!--cheyz add
       scalar_static_albedo_at_pc%f(1:12   ,1:mesh%nv)  = staticData_albedo_at_pc_surface%f(1:12   ,1:mesh%nv)
       scalar_static_greenfrac_at_pc%f(1:12,1:mesh%nv)  = staticData_greenfrac_at_pc_surface%f(1:12,1:mesh%nv)
       scalar_static_lu_index_at_pc%f(      1:mesh%nv)  = staticData_luIndex_at_pc_surface%f(       1:mesh%nv) 
       scalar_static_soilcat_top_at_pc%f(   1:mesh%nv)  = staticData_soilTypetop_at_pc_surface%f(   1:mesh%nv) 
       scalar_static_snoalb_at_pc%f(        1:mesh%nv)  = staticData_snoalb_at_pc_surface%f(        1:mesh%nv) 
       scalar_static_ter_at_pc%f(           1:mesh%nv)  = staticData_phis_at_pc_surface%f(          1:mesh%nv) 
    end if
#endif
    return
   end subroutine grist_lnd_static_vars_construct

   subroutine grist_lnd_static_vars_destruct()

      if(allocated(scalar_static_albedo_at_pc%f))      deallocate(scalar_static_albedo_at_pc%f)
      if(allocated(scalar_static_greenfrac_at_pc%f))   deallocate(scalar_static_greenfrac_at_pc%f)
      if(allocated(scalar_static_lu_index_at_pc%f))    deallocate(scalar_static_lu_index_at_pc%f)
      if(allocated(scalar_static_landmask_at_pc%f))    deallocate(scalar_static_landmask_at_pc%f)
      if(allocated(scalar_static_soilcat_top_at_pc%f)) deallocate(scalar_static_soilcat_top_at_pc%f)
      if(allocated(scalar_static_snoalb_at_pc%f))      deallocate(scalar_static_snoalb_at_pc%f)
      if(allocated(scalar_static_ter_at_pc%f))         deallocate(scalar_static_ter_at_pc%f)
      !if(allocated(scalar_static_soiltemp_at_pc%f))    deallocate(scalar_static_soiltemp_at_pc%f)
      !if(allocated(scalar_static_shdmin_at_pc%f))      deallocate(scalar_static_shdmin_at_pc%f)
      !if(allocated(scalar_static_shdmax_at_pc%f))      deallocate(scalar_static_shdmax_at_pc%f)

      return
    end subroutine grist_lnd_static_vars_destruct
   end module grist_lnd_static_vars_module
