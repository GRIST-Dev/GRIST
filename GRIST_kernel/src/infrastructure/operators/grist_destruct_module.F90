
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: Destruct non-runtime mesh info       
! Revision history:        
!----------------------------------------------------------------------------

module grist_destruct_module

  use grist_constants,       only: i4   
  use grist_domain_types,    only: global_domain, &
                                   global_domain_data

  implicit none

  private

  public  :: destruct_global_mesh_var, &
             destruct_local_mesh_element

contains

  subroutine destruct_global_mesh_var(mesh)
  
  ! io
    type(global_domain_data), intent(inout) :: mesh
  ! local
    integer(i4)                 ::  iv,ie,it
 

    if (allocated(mesh%tri_v%f))  deallocate(mesh%tri_v%f)
    if (allocated(mesh%edt_v%f))  deallocate(mesh%edt_v%f)
    if (allocated(mesh%vtx_nnb))  deallocate(mesh%vtx_nnb)
    if (allocated(mesh%vtx_ltln)) deallocate(mesh%vtx_ltln)
    if (allocated(mesh%vtx_nb%f)) deallocate(mesh%vtx_nb%f)
    if (allocated(mesh%vtx_ed%f)) deallocate(mesh%vtx_ed%f)
    if (allocated(mesh%vtx_tr%f)) deallocate(mesh%vtx_tr%f)
    if (allocated(mesh%parts))    deallocate(mesh%parts)
    if (allocated(mesh%tri_nnb%f))deallocate(mesh%tri_nnb%f)
#ifdef LAM_DOMAIN
    if (allocated(mesh%plg_mark%f)) deallocate(mesh%plg_mark%f)
    if (allocated(mesh%edt_mark%f)) deallocate(mesh%edt_mark%f)
    if (allocated(mesh%tri_mark%f)) deallocate(mesh%tri_mark%f)
#endif
    return
  end subroutine destruct_global_mesh_var
!
! destruct decomposed mesh vars that have not been used;
! or already replaced by compute pattern
!

  subroutine destruct_local_mesh_element(mesh)
! io
    type(global_domain), intent(inout) :: mesh
! local

    if(allocated(mesh%vtx)) deallocate(mesh%vtx)
    if(allocated(mesh%plg)) deallocate(mesh%plg)
    if(allocated(mesh%tri)) deallocate(mesh%tri)
    if(allocated(mesh%edt)) deallocate(mesh%edt)
    if(allocated(mesh%edp)) deallocate(mesh%edp)
!
! because wrfPhys needs this info, delete it after gcm_init, not inside mesh weight
!    if(allocated(mesh%vtxCellLeng)) deallocate(mesh%vtxCellLeng)

    return
  end subroutine destruct_local_mesh_element
end module grist_destruct_module
