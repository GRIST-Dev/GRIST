
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: assign mesh data
!
! Revision history:
!----------------------------------------------------------------------------

 module grist_config_mesh
 
   use grist_constants,       only: i4   
   use grist_domain_types,    only: global_domain, global_domain_data
   use grist_grid_file_vars


   implicit none

   private

   public :: config_mesh_from_file

   contains

   subroutine config_mesh_from_file(mesh)
! io
     type(global_domain_data), intent(inout) :: mesh
! local
     integer(i4)                 ::  it
     integer(i4)                 ::  ie
     integer(i4)                 ::  iv

!--------------------------
!    Global info 
!--------------------------
         mesh%nt            = mesh_nt
         mesh%ne            = mesh_ne
         mesh%nv            = mesh_nv
                           
         mesh%kind          = mesh_kind
         mesh%node          = mesh_node
         mesh%optm          = mesh_optm
         mesh%glevel        = mesh_glevel
         mesh%maxvnb        = mesh_maxvnb

         mesh%min_edt_dist  = mesh_min_edt_dist 
         mesh%max_edt_dist  = mesh_max_edt_dist 
         mesh%mean_edt_dist = mesh_mean_edt_dist
                           
         mesh%min_edp_dist  = mesh_min_edp_dist 
         mesh%max_edp_dist  = mesh_max_edp_dist 
         mesh%mean_edp_dist = mesh_mean_edp_dist
                           
         mesh%min_tri_area  = mesh_min_tri_area
         mesh%max_tri_area  = mesh_max_tri_area
         mesh%mean_tri_area = mesh_mean_tri_area
                           
         mesh%min_plg_area  = mesh_min_plg_area
         mesh%max_plg_area  = mesh_max_plg_area 
         mesh%mean_plg_area = mesh_mean_plg_area
                           
         mesh%min_tri_angle = mesh_min_tri_angle
         mesh%max_tri_angle = mesh_max_tri_angle
         mesh%mean_tri_angle= mesh_mean_tri_angle

         mesh%gridFileName  = mesh_name

     return
   end subroutine config_mesh_from_file

 end module grist_config_mesh 
