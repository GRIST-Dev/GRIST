!======================================================
!
!  Created by LiXiaohan on 20/4/20.
!  Adopted from grist_gcm_diagnose_module
!  
!  Purpose:
!  currently used for accumulating vars like precip
!======================================================


 module grist_scm_diagnose_module

  use grist_constants,              only: i4, r8, zero, gravity
  use grist_domain_types,           only: global_domain
  use grist_data_types,             only: scalar_1d_field,  &
                                          scalar_2d_field
  use grist_physics_data_structure, only: pstate
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,   only: pstate_cam,                 & 
                                          ptend_vertical_diffusion,   &
                                          ptend_shallow_convection,   &
                                          ptend_deep_convection,      &
                                          ptend_microphysics,         &
                                          ptend_macrophysics

  use grist_physics_update,         only: old_time_level
#endif
#ifdef AMIPW_PHYSICS
  use grist_PhysW_data_structure,   only: pstate_wrf,ptend_wrf,p_qv,p_qc,p_qi,p_qs
#endif
  use grist_nml_module,             only: nlev
 
  implicit none

   private

   public   :: scm_diagnose_init,  &
               scm_diagnose_final, & 
               scm_accumulate_physics_variables, &
               scm_dump_physics_variables, &
               scm_reset_physics_variables, &
               diag_physics_vars

   type scm_diag_vars
        type(scalar_1d_field)  :: precc
        type(scalar_1d_field)  :: precl
        type(scalar_1d_field)  :: prect
        type(scalar_2d_field)  :: cloud
        type(scalar_2d_field)  :: cldliq
        type(scalar_2d_field)  :: cldice
        type(scalar_1d_field)  :: lwcf
        type(scalar_1d_field)  :: swcf
        type(scalar_1d_field)  :: tmq
        type(scalar_2d_field)  :: dq_deep
        type(scalar_2d_field)  :: dq_shallow
        type(scalar_2d_field)  :: dq_micro
        type(scalar_2d_field)  :: dq_macro
        type(scalar_2d_field)  :: dq_pbl

        integer(i4)            :: ncount
   end type scm_diag_vars

   type(scm_diag_vars)  :: diag_physics_vars

  contains

  subroutine scm_diagnose_init(mesh)
    type(global_domain),  intent(in)   :: mesh
    if(.not.allocated(diag_physics_vars%precc%f)) allocate(diag_physics_vars%precc%f(mesh%nv_full))
    if(.not.allocated(diag_physics_vars%precl%f)) allocate(diag_physics_vars%precl%f(mesh%nv_full))
    if(.not.allocated(diag_physics_vars%prect%f)) allocate(diag_physics_vars%prect%f(mesh%nv_full))

    if(.not.allocated(diag_physics_vars%cloud%f)) allocate(diag_physics_vars%cloud%f(nlev,mesh%nv_full))
    if(.not.allocated(diag_physics_vars%cldliq%f)) allocate(diag_physics_vars%cldliq%f(nlev,mesh%nv_full))
    if(.not.allocated(diag_physics_vars%cldice%f)) allocate(diag_physics_vars%cldice%f(nlev,mesh%nv_full))
    if(.not.allocated(diag_physics_vars%tmq%f)) allocate(diag_physics_vars%tmq%f(mesh%nv_full))
    if(.not.allocated(diag_physics_vars%lwcf%f)) allocate(diag_physics_vars%lwcf%f(mesh%nv_full))
    if(.not.allocated(diag_physics_vars%swcf%f)) allocate(diag_physics_vars%swcf%f(mesh%nv_full))
    if(.not.allocated(diag_physics_vars%dq_deep%f)) allocate(diag_physics_vars%dq_deep%f(nlev,mesh%nv_full))
    if(.not.allocated(diag_physics_vars%dq_shallow%f)) allocate(diag_physics_vars%dq_shallow%f(nlev,mesh%nv_full))
    if(.not.allocated(diag_physics_vars%dq_micro%f)) allocate(diag_physics_vars%dq_micro%f(nlev,mesh%nv_full))
    if(.not.allocated(diag_physics_vars%dq_macro%f)) allocate(diag_physics_vars%dq_macro%f(nlev,mesh%nv_full))
    if(.not.allocated(diag_physics_vars%dq_pbl%f)) allocate(diag_physics_vars%dq_pbl%f(nlev,mesh%nv_full))

    diag_physics_vars%ncount  = 0
    diag_physics_vars%precc%f = zero
    diag_physics_vars%precl%f = zero
    diag_physics_vars%prect%f = zero
    diag_physics_vars%cloud%f = zero
    diag_physics_vars%cldliq%f = zero
    diag_physics_vars%cldice%f = zero
    diag_physics_vars%tmq%f = zero
    diag_physics_vars%lwcf%f = zero
    diag_physics_vars%swcf%f = zero
    diag_physics_vars%dq_deep%f = zero
    diag_physics_vars%dq_shallow%f = zero
    diag_physics_vars%dq_micro%f = zero
    diag_physics_vars%dq_macro%f = zero
    diag_physics_vars%dq_pbl%f = zero

    diag_physics_vars%precc%pos = 0
    diag_physics_vars%precl%pos = 0
    diag_physics_vars%prect%pos = 0
    diag_physics_vars%cloud%pos = 0
    diag_physics_vars%cldliq%pos = 0
    diag_physics_vars%cldice%pos = 0
    diag_physics_vars%tmq%pos = 0
    diag_physics_vars%lwcf%pos = 0
    diag_physics_vars%swcf%pos = 0
    diag_physics_vars%dq_deep%pos = 0
    diag_physics_vars%dq_shallow%pos = 0
    diag_physics_vars%dq_micro%pos = 0
    diag_physics_vars%dq_macro%pos = 0
    diag_physics_vars%dq_pbl%pos = 0

    return
  end subroutine scm_diagnose_init

  subroutine scm_diagnose_final
    if(allocated(diag_physics_vars%precc%f)) deallocate(diag_physics_vars%precc%f)
    if(allocated(diag_physics_vars%precl%f)) deallocate(diag_physics_vars%precl%f)
    if(allocated(diag_physics_vars%prect%f)) deallocate(diag_physics_vars%prect%f)
    if(allocated(diag_physics_vars%cloud%f)) deallocate(diag_physics_vars%cloud%f)
    if(allocated(diag_physics_vars%cldliq%f))deallocate(diag_physics_vars%cldliq%f)
    if(allocated(diag_physics_vars%cldice%f))deallocate(diag_physics_vars%cldice%f)
    if(allocated(diag_physics_vars%tmq%f))   deallocate(diag_physics_vars%tmq%f)
    if(allocated(diag_physics_vars%lwcf%f)) deallocate(diag_physics_vars%lwcf%f)
    if(allocated(diag_physics_vars%swcf%f)) deallocate(diag_physics_vars%swcf%f)
    if(allocated(diag_physics_vars%dq_deep%f))   deallocate(diag_physics_vars%dq_deep%f)
    if(allocated(diag_physics_vars%dq_shallow%f))   deallocate(diag_physics_vars%dq_shallow%f)
    if(allocated(diag_physics_vars%dq_micro%f))   deallocate(diag_physics_vars%dq_micro%f)
    if(allocated(diag_physics_vars%dq_macro%f))   deallocate(diag_physics_vars%dq_macro%f)
    if(allocated(diag_physics_vars%dq_pbl%f))   deallocate(diag_physics_vars%dq_pbl%f)
    return
  end subroutine scm_diagnose_final
!
! accumulate at each model step
!
  subroutine scm_accumulate_physics_variables(ncol)
    integer, intent(in)     :: ncol
! local
    integer  :: i,k
    real(r8) :: ftem(ncol)

    diag_physics_vars%precl%f(1:ncol) = diag_physics_vars%precl%f(1:ncol)+pstate%scalar_precl_surface%f(1:ncol)
    diag_physics_vars%precc%f(1:ncol) = diag_physics_vars%precc%f(1:ncol)+pstate%scalar_precc_surface%f(1:ncol)
    diag_physics_vars%prect%f(1:ncol) = diag_physics_vars%prect%f(1:ncol)+pstate%scalar_prect_surface%f(1:ncol)
 
#ifdef AMIPC_PHYSICS
    diag_physics_vars%cloud%f(:,1:ncol) = diag_physics_vars%cloud%f(:,1:ncol)+pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)
    diag_physics_vars%cldliq%f(:,1:ncol)= diag_physics_vars%cldliq%f(:,1:ncol)+pstate%tracer_mxrt_at_pc_full_level%f(2,:,1:ncol)
    diag_physics_vars%cldice%f(:,1:ncol)= diag_physics_vars%cldice%f(:,1:ncol)+pstate%tracer_mxrt_at_pc_full_level%f(3,:,1:ncol)
    diag_physics_vars%tmq%f(1:ncol)     = diag_physics_vars%tmq%f(1:ncol)+pstate_cam%diag_tmq%f(1:ncol) 
    diag_physics_vars%lwcf%f(1:ncol)  = diag_physics_vars%lwcf%f(1:ncol)+pstate_cam%lwcf_at_pc_top%f(1:ncol)
    diag_physics_vars%swcf%f(1:ncol)  = diag_physics_vars%swcf%f(1:ncol)+pstate_cam%swcf_at_pc_top%f(1:ncol)
    diag_physics_vars%dq_deep%f(:,1:ncol) = diag_physics_vars%dq_deep%f(:,1:ncol)+ptend_deep_convection%tend_q%f(1,:,1:ncol)
    diag_physics_vars%dq_shallow%f(:,1:ncol) = diag_physics_vars%dq_shallow%f(:,1:ncol)+ptend_shallow_convection%tend_q%f(1,:,1:ncol)
    diag_physics_vars%dq_micro%f(:,1:ncol) = diag_physics_vars%dq_micro%f(:,1:ncol)+ptend_microphysics%tend_q%f(1,:,1:ncol)
    diag_physics_vars%dq_macro%f(:,1:ncol) = diag_physics_vars%dq_macro%f(:,1:ncol)+ptend_macrophysics%tend_q%f(1,:,1:ncol)
    diag_physics_vars%dq_pbl%f(:,1:ncol) = diag_physics_vars%dq_pbl%f(:,1:ncol)+ptend_vertical_diffusion%tend_q%f(1,:,1:ncol)
#endif

#ifdef AMIPW_PHYSICS
    do i = 1, ncol
      do k = 1, nlev
        diag_physics_vars%cloud%f(k,i)    = diag_physics_vars%cloud%f(k,i)+pstate_wrf%cldfra(i,nlev+1-k,1)
        diag_physics_vars%cldliq%f(k,i)   = diag_physics_vars%cldliq%f(k,i)+pstate_wrf%moist(i,nlev+1-k,1,p_qc)
        diag_physics_vars%cldice%f(k,i)   = diag_physics_vars%cldice%f(k,i)+pstate_wrf%moist(i,nlev+1-k,1,p_qi)+pstate_wrf%moist(i,nlev+1-k,1,p_qs)
        diag_physics_vars%dq_deep%f(k,i)  = diag_physics_vars%dq_deep%f(k,i)+ptend_wrf%rqvcuten(i,nlev+1-k,1)
        diag_physics_vars%dq_micro%f(k,i) = diag_physics_vars%dq_micro%f(k,i)+ptend_wrf%rqvmpten(i,nlev+1-k,1)
        diag_physics_vars%dq_pbl%f(k,i)   = diag_physics_vars%dq_pbl%f(k,i)+ptend_wrf%rqvblten(i,nlev+1-k,1)
      end do

      ftem = 0._r8
      do k = 1, nlev
          ftem(i) = ftem(i) + pstate_wrf%moist(i,nlev+1-k,1,p_qv)*pstate%delp_at_pc_full_level%f(k,i)/gravity
      end do
    end do
    diag_physics_vars%tmq%f(1:ncol)     = diag_physics_vars%tmq%f(1:ncol)+ftem(1:ncol) 

    diag_physics_vars%lwcf%f(1:ncol)    = diag_physics_vars%lwcf%f(1:ncol)+pstate_wrf%lwcf(1:ncol,1)
    diag_physics_vars%swcf%f(1:ncol)    = diag_physics_vars%swcf%f(1:ncol)+pstate_wrf%swcf(1:ncol,1)
#endif

    diag_physics_vars%ncount  = diag_physics_vars%ncount +1
    return
  end subroutine scm_accumulate_physics_variables
!
! dump depending on write_history frequency
!
  subroutine scm_dump_physics_variables

use grist_mpi

! local
    diag_physics_vars%precl%f = diag_physics_vars%precl%f/diag_physics_vars%ncount
    diag_physics_vars%precc%f = diag_physics_vars%precc%f/diag_physics_vars%ncount
    diag_physics_vars%prect%f = diag_physics_vars%prect%f/diag_physics_vars%ncount

    diag_physics_vars%cloud%f = diag_physics_vars%cloud%f/diag_physics_vars%ncount
    diag_physics_vars%cldliq%f= diag_physics_vars%cldliq%f/diag_physics_vars%ncount
    diag_physics_vars%cldice%f= diag_physics_vars%cldice%f/diag_physics_vars%ncount
    diag_physics_vars%tmq%f   = diag_physics_vars%tmq%f/diag_physics_vars%ncount
    diag_physics_vars%lwcf%f  = diag_physics_vars%lwcf%f/diag_physics_vars%ncount
    diag_physics_vars%swcf%f  = diag_physics_vars%swcf%f/diag_physics_vars%ncount
    diag_physics_vars%dq_deep%f = diag_physics_vars%dq_deep%f/diag_physics_vars%ncount
    diag_physics_vars%dq_shallow%f = diag_physics_vars%dq_shallow%f/diag_physics_vars%ncount
    diag_physics_vars%dq_micro%f = diag_physics_vars%dq_micro%f/diag_physics_vars%ncount
    diag_physics_vars%dq_macro%f = diag_physics_vars%dq_macro%f/diag_physics_vars%ncount
    diag_physics_vars%dq_pbl%f = diag_physics_vars%dq_pbl%f/diag_physics_vars%ncount
   
    return
  end subroutine scm_dump_physics_variables

  subroutine scm_reset_physics_variables
! reset to zero
    diag_physics_vars%ncount  = 0
    diag_physics_vars%precl%f = zero
    diag_physics_vars%precc%f = zero
    diag_physics_vars%prect%f = zero

    diag_physics_vars%cloud%f = zero
    diag_physics_vars%cldliq%f= zero
    diag_physics_vars%cldice%f= zero
    diag_physics_vars%tmq%f   = zero
    diag_physics_vars%lwcf%f  = zero
    diag_physics_vars%swcf%f  = zero
    diag_physics_vars%dq_deep%f = zero 
    diag_physics_vars%dq_shallow%f =  zero
    diag_physics_vars%dq_micro%f =  zero
    diag_physics_vars%dq_macro%f =  zero
    diag_physics_vars%dq_pbl%f =  zero
 
    return
  end subroutine scm_reset_physics_variables

  end module grist_scm_diagnose_module
