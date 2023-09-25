!======================================================
!
!  Created by LiXiaohan on 19/6/30.
!  contains:
!  coupling of each physical process in the physpkg
!  time level update in physics variables
!======================================================

 module grist_PhysC_diag

    use grist_constants,                only: i4, zero
    use grist_physics_data_structure,   only: pstate
    use grist_PhysC_data_structure,     only: pstate_cam
 
    implicit none
    private
    save

    public  :: grist_PhysC_diag_run

contains

  subroutine grist_PhysC_diag_run(ncol)
      integer(i4),  intent(in) :: ncol
!
! evaluate precc, precl, snowc, snowl, prect based on internal states
! after this, zero all internal states, snowc is done inside dp&sh
!
       pstate%scalar_precc_surface%f(:ncol) = pstate_cam%scalar_precc_dp_surface%f(:ncol) + pstate_cam%scalar_precc_sh_surface%f(:ncol)
       pstate%scalar_precl_surface%f(:ncol) = pstate_cam%str_prec_surface%f(:ncol)
       pstate%scalar_snowc_surface%f(:ncol) = pstate_cam%scalar_snowc_surface%f(:ncol)
       pstate%scalar_snowl_surface%f(:ncol) = pstate_cam%str_snow_surface%f(:ncol)
       pstate%scalar_prect_surface%f(:ncol) = pstate%scalar_precl_surface%f(:ncol) + pstate%scalar_precc_surface%f(:ncol)

       call physics_diag_clean

    return
  end subroutine grist_PhysC_diag_run

  subroutine physics_diag_clean

       pstate_cam%scalar_precc_dp_surface%f = zero
       pstate_cam%scalar_precc_sh_surface%f = zero 
       pstate_cam%scalar_snowc_surface%f    = zero

    return
  end subroutine physics_diag_clean

 end module grist_PhysC_diag
