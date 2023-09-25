
!----------------------------------------------------------------------------
! Created on 2017-5-17
! Author: Yi Zhang
! Version 1.0
! Description: GCM Diagnostic module
! Revision history: based on the dycore diagnose module
!                   currently used for accumulating vars like precip
!----------------------------------------------------------------------------

 module grist_gcm_diagnose_h0_module

  use grist_constants,              only: i4, r8, zero, one, gravity
  use grist_domain_types,           only: global_domain
  use grist_data_types,             only: scalar_1d_field,  &
                                          scalar_2d_field, &
                                          scalar_3d_field
  use grist_dycore_vars_module,     only: dycoreVarSurface, dycoreVarCellFull, dycoreVarCellFace
  use grist_tracer_transport_vars_module, only: tracerVarCellFull
  use grist_physics_data_structure, only: pstate
  use grist_nml_module,             only: ntracer, nlev, nlevp, &
                                          write_history_h0, mif_index, nmif
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,   only: pstate_cam
  use grist_physics_update,         only: old_time_level
#endif

#ifdef AMIPW_PHYSICS
use grist_PhysW_data_structure,   only: pstate_wrf, psurf_wrf, ptend_wrf
use wv_saturation,                only: aqsat_grist
use grist_math_module,            only: lininterp
#endif

  implicit none

   private

   public   :: gcm_h0_diagnose_init,  &
               gcm_h0_diagnose_final, & 
               gcm_h0_accu_physics_variables, &
               gcm_h0_dump_physics_variables, &
               gcm_h0_rest_physics_variables, &
               diag_physics_vars_h0

   type gcm_diag_vars
        type(scalar_1d_field)  :: precc
        type(scalar_1d_field)  :: precl
        type(scalar_1d_field)  :: prect
        type(scalar_1d_field)  :: ps
        type(scalar_1d_field)  :: snowh
        type(scalar_1d_field)  :: ts
        type(scalar_1d_field)  :: shflx
        type(scalar_1d_field)  :: qflx
        type(scalar_2d_field)  :: phi
        type(scalar_2d_field)  :: uwind
        type(scalar_2d_field)  :: vwind
        type(scalar_2d_field)  :: omega
        type(scalar_2d_field)  :: temp
        type(scalar_2d_field)  :: qv
        type(scalar_2d_field)  :: qc
        type(scalar_2d_field)  :: qi
        type(scalar_2d_field)  :: rh
        type(scalar_2d_field)  :: cld
        type(scalar_2d_field)  :: concld
        type(scalar_2d_field)  :: lwc
        type(scalar_2d_field)  :: iwc
        type(scalar_1d_field)  :: cldtot
        type(scalar_1d_field)  :: cldlow
        type(scalar_1d_field)  :: cldmed
        type(scalar_1d_field)  :: cldhgh
        type(scalar_1d_field)  :: z500
        type(scalar_1d_field)  :: u200
        type(scalar_1d_field)  :: u850
        type(scalar_1d_field)  :: v850
        type(scalar_1d_field)  :: flwut
        type(scalar_1d_field)  :: flwdt
        type(scalar_1d_field)  :: flwus
        type(scalar_1d_field)  :: flwds
        type(scalar_1d_field)  :: flwutc
        type(scalar_1d_field)  :: flwdtc
        type(scalar_1d_field)  :: flwusc
        type(scalar_1d_field)  :: flwdsc
        type(scalar_1d_field)  :: fswut
        type(scalar_1d_field)  :: fswdt
        type(scalar_1d_field)  :: fswus
        type(scalar_1d_field)  :: fswds
        type(scalar_1d_field)  :: fsntoa
        type(scalar_1d_field)  :: fswutc
        type(scalar_1d_field)  :: fswdtc
        type(scalar_1d_field)  :: fswusc
        type(scalar_1d_field)  :: fswdsc
        type(scalar_1d_field)  :: fsntoac
        type(scalar_1d_field)  :: lwcf
        type(scalar_1d_field)  :: swcf
        type(scalar_1d_field)  :: tmq
        type(scalar_1d_field)  :: taux
        type(scalar_1d_field)  :: tauy
        type(scalar_1d_field)  :: lwp
        type(scalar_1d_field)  :: iwp

!-------LiXH Test for tuning GaussPDF------->
        type(scalar_2d_field)  :: sgm
!<------LiXH Test for tuning GaussPDF--------


        integer(i4)            :: ncount
   end type gcm_diag_vars

   type(gcm_diag_vars)  :: diag_physics_vars_h0
   integer(i4) :: ncell

  contains

  subroutine gcm_h0_diagnose_init(mesh)
    type(global_domain),  intent(in)   :: mesh

    diag_physics_vars_h0%ncount  = 0

    if(write_history_h0)then
    call wrap_allocate_data1d(mesh,diag_physics_vars_h0%precc)
    call wrap_allocate_data1d(mesh,diag_physics_vars_h0%precl)
    call wrap_allocate_data1d(mesh,diag_physics_vars_h0%prect)

    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%ps)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%ts)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%snowh)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%shflx)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%qflx)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%phi)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%uwind)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%vwind)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%omega)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%temp)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%qv)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%qc)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%qi)

    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%rh)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%cld)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%concld)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%lwc)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%iwc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%cldtot)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%cldlow)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%cldmed)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%cldhgh)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%z500)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%u850)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%u200)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%v850)

    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwut)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwdt)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwus)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwds)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwutc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwdtc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwusc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%flwdsc)

    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswut)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswdt)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswus)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswds)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fsntoa)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswutc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswdtc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswusc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fswdsc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%fsntoac)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%lwcf)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%swcf)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%tmq)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%taux)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%tauy)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%lwp)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%iwp)


#ifdef AMIPC_PHYSICS
    if(allocated(pstate_cam%macrop_sgm_at_pc_full_level%f))call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_h0%sgm)
#endif

    end if

    ncell = mesh%nv_halo(1)

    return
  end subroutine gcm_h0_diagnose_init

  subroutine gcm_h0_diagnose_final
    if(allocated(diag_physics_vars_h0%precc%f))  deallocate(diag_physics_vars_h0%precc%f)
    if(allocated(diag_physics_vars_h0%precl%f))  deallocate(diag_physics_vars_h0%precl%f)
    if(allocated(diag_physics_vars_h0%prect%f))  deallocate(diag_physics_vars_h0%prect%f)

    if(allocated(diag_physics_vars_h0%ps%f))     deallocate(diag_physics_vars_h0%ps%f)
    if(allocated(diag_physics_vars_h0%snowh%f))  deallocate(diag_physics_vars_h0%snowh%f)
    if(allocated(diag_physics_vars_h0%ts%f))     deallocate(diag_physics_vars_h0%ts%f)
    if(allocated(diag_physics_vars_h0%shflx%f))  deallocate(diag_physics_vars_h0%shflx%f)
    if(allocated(diag_physics_vars_h0%qflx%f))   deallocate(diag_physics_vars_h0%qflx%f)
    if(allocated(diag_physics_vars_h0%phi%f))    deallocate(diag_physics_vars_h0%phi%f)
    if(allocated(diag_physics_vars_h0%uwind%f))  deallocate(diag_physics_vars_h0%uwind%f)
    if(allocated(diag_physics_vars_h0%vwind%f))  deallocate(diag_physics_vars_h0%vwind%f)
    if(allocated(diag_physics_vars_h0%omega%f))  deallocate(diag_physics_vars_h0%omega%f)
    if(allocated(diag_physics_vars_h0%temp%f))   deallocate(diag_physics_vars_h0%temp%f)
    if(allocated(diag_physics_vars_h0%qv%f))     deallocate(diag_physics_vars_h0%qv%f)
    if(allocated(diag_physics_vars_h0%qc%f))     deallocate(diag_physics_vars_h0%qc%f)
    if(allocated(diag_physics_vars_h0%qi%f))     deallocate(diag_physics_vars_h0%qi%f)
    if(allocated(diag_physics_vars_h0%rh%f))     deallocate(diag_physics_vars_h0%rh%f)
    if(allocated(diag_physics_vars_h0%cld%f))    deallocate(diag_physics_vars_h0%cld%f)
    if(allocated(diag_physics_vars_h0%concld%f))    deallocate(diag_physics_vars_h0%concld%f)
    if(allocated(diag_physics_vars_h0%lwc%f))    deallocate(diag_physics_vars_h0%lwc%f)
    if(allocated(diag_physics_vars_h0%iwc%f))    deallocate(diag_physics_vars_h0%iwc%f)
    if(allocated(diag_physics_vars_h0%cldtot%f)) deallocate(diag_physics_vars_h0%cldtot%f)
    if(allocated(diag_physics_vars_h0%cldlow%f)) deallocate(diag_physics_vars_h0%cldlow%f)
    if(allocated(diag_physics_vars_h0%cldmed%f)) deallocate(diag_physics_vars_h0%cldmed%f)
    if(allocated(diag_physics_vars_h0%cldhgh%f)) deallocate(diag_physics_vars_h0%cldhgh%f)
    if(allocated(diag_physics_vars_h0%z500%f))   deallocate(diag_physics_vars_h0%z500%f)
    if(allocated(diag_physics_vars_h0%u850%f))   deallocate(diag_physics_vars_h0%u850%f)
    if(allocated(diag_physics_vars_h0%u200%f))   deallocate(diag_physics_vars_h0%u200%f)
    if(allocated(diag_physics_vars_h0%v850%f))   deallocate(diag_physics_vars_h0%v850%f)

    if(allocated(diag_physics_vars_h0%flwut%f))  deallocate(diag_physics_vars_h0%flwut%f )
    if(allocated(diag_physics_vars_h0%flwdt%f))  deallocate(diag_physics_vars_h0%flwdt%f )
    if(allocated(diag_physics_vars_h0%flwus%f))  deallocate(diag_physics_vars_h0%flwus%f )
    if(allocated(diag_physics_vars_h0%flwds%f))  deallocate(diag_physics_vars_h0%flwds%f )
    if(allocated(diag_physics_vars_h0%flwutc%f)) deallocate(diag_physics_vars_h0%flwutc%f)
    if(allocated(diag_physics_vars_h0%flwdtc%f)) deallocate(diag_physics_vars_h0%flwdtc%f)
    if(allocated(diag_physics_vars_h0%flwusc%f)) deallocate(diag_physics_vars_h0%flwusc%f)
    if(allocated(diag_physics_vars_h0%flwdsc%f)) deallocate(diag_physics_vars_h0%flwdsc%f)
    if(allocated(diag_physics_vars_h0%fswut%f))  deallocate(diag_physics_vars_h0%fswut%f )
    if(allocated(diag_physics_vars_h0%fswdt%f))  deallocate(diag_physics_vars_h0%fswdt%f )
    if(allocated(diag_physics_vars_h0%fswus%f))  deallocate(diag_physics_vars_h0%fswus%f )
    if(allocated(diag_physics_vars_h0%fswds%f))  deallocate(diag_physics_vars_h0%fswds%f )
    if(allocated(diag_physics_vars_h0%fsntoa%f)) deallocate(diag_physics_vars_h0%fsntoa%f )
    if(allocated(diag_physics_vars_h0%fswutc%f)) deallocate(diag_physics_vars_h0%fswutc%f)
    if(allocated(diag_physics_vars_h0%fswdtc%f)) deallocate(diag_physics_vars_h0%fswdtc%f)
    if(allocated(diag_physics_vars_h0%fswusc%f)) deallocate(diag_physics_vars_h0%fswusc%f)
    if(allocated(diag_physics_vars_h0%fswdsc%f)) deallocate(diag_physics_vars_h0%fswdsc%f)
    if(allocated(diag_physics_vars_h0%fsntoac%f))deallocate(diag_physics_vars_h0%fsntoac%f )
    if(allocated(diag_physics_vars_h0%lwcf%f))   deallocate(diag_physics_vars_h0%lwcf%f)
    if(allocated(diag_physics_vars_h0%swcf%f))   deallocate(diag_physics_vars_h0%swcf%f)
    if(allocated(diag_physics_vars_h0%tmq%f))    deallocate(diag_physics_vars_h0%tmq%f)
    if(allocated(diag_physics_vars_h0%taux%f))   deallocate(diag_physics_vars_h0%taux%f)
    if(allocated(diag_physics_vars_h0%tauy%f))   deallocate(diag_physics_vars_h0%tauy%f)
    if(allocated(diag_physics_vars_h0%lwp%f))    deallocate(diag_physics_vars_h0%lwp%f)
    if(allocated(diag_physics_vars_h0%iwp%f))    deallocate(diag_physics_vars_h0%iwp%f)

    if(allocated(diag_physics_vars_h0%sgm%f))    deallocate(diag_physics_vars_h0%sgm%f)
    return
  end subroutine gcm_h0_diagnose_final
!
! accumulate at each model step
!
  subroutine gcm_h0_accu_physics_variables
! local
#ifdef AMIPW_PHYSICS
integer  :: icell, ilev
real(r8) :: sat_specific_humidity(nlev,ncell), esvp(nlev,ncell)
real(r8) :: dpres(nlev), ftem(nlev),ftem_qc(nlev),ftem_qi(nlev)
real(r8) :: pres_out(3), vartmp_out(3)
#endif

    if(write_history_h0)then

    diag_physics_vars_h0%precl%f(1:ncell) = diag_physics_vars_h0%precl%f(1:ncell)+pstate%scalar_precl_surface%f(1:ncell)
    diag_physics_vars_h0%precc%f(1:ncell) = diag_physics_vars_h0%precc%f(1:ncell)+pstate%scalar_precc_surface%f(1:ncell)
    diag_physics_vars_h0%prect%f(1:ncell) = diag_physics_vars_h0%prect%f(1:ncell)+pstate%scalar_prect_surface%f(1:ncell)


    diag_physics_vars_h0%ps%f(1:ncell)    = diag_physics_vars_h0%ps%f(1:ncell)+dycoreVarSurface%scalar_hpressure_n%f(1:ncell)
    diag_physics_vars_h0%ts%f(1:ncell)    = diag_physics_vars_h0%ts%f(1:ncell)+pstate%ts_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%shflx%f(1:ncell) = diag_physics_vars_h0%shflx%f(1:ncell)+pstate%atm_in_shflx_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%qflx%f(1:ncell)  = diag_physics_vars_h0%qflx%f(1:ncell)+pstate%atm_in_qflx_at_pc_surface%f(1,1:ncell)
    diag_physics_vars_h0%phi%f(:,1:ncell) = diag_physics_vars_h0%phi%f(:,1:ncell)+dycoreVarCellFull%scalar_geopotential_n%f(:,1:ncell)
    diag_physics_vars_h0%uwind%f(:,1:ncell) = diag_physics_vars_h0%uwind%f(:,1:ncell)+dycoreVarCellFull%scalar_U_wind_n%f(:,1:ncell)
    diag_physics_vars_h0%vwind%f(:,1:ncell) = diag_physics_vars_h0%vwind%f(:,1:ncell)+dycoreVarCellFull%scalar_V_wind_n%f(:,1:ncell)
    diag_physics_vars_h0%omega%f(:,1:ncell) = diag_physics_vars_h0%omega%f(:,1:ncell)+dycoreVarCellFull%scalar_omega_n%f(:,1:ncell)
    diag_physics_vars_h0%temp%f(:,1:ncell)  = diag_physics_vars_h0%temp%f(:,1:ncell)+dycoreVarCellFull%scalar_temp_n%f(:,1:ncell)
    diag_physics_vars_h0%qv%f(:,1:ncell)    = diag_physics_vars_h0%qv%f(:,1:ncell)+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,1:ncell)

#ifdef AMIPC_PHYSICS
    diag_physics_vars_h0%qc%f(:,1:ncell)    = diag_physics_vars_h0%qc%f(:,1:ncell)+tracerVarCellFull%scalar_tracer_mxrt_n%f(2,:,1:ncell)
    diag_physics_vars_h0%qi%f(:,1:ncell)    = diag_physics_vars_h0%qi%f(:,1:ncell)+tracerVarCellFull%scalar_tracer_mxrt_n%f(3,:,1:ncell)
    diag_physics_vars_h0%snowh%f(1:ncell) = diag_physics_vars_h0%snowh%f(1:ncell)+pstate%snowhland_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%taux%f(1:ncell)  = diag_physics_vars_h0%taux%f(1:ncell)+pstate%atm_in_taux_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%tauy%f(1:ncell)  = diag_physics_vars_h0%tauy%f(1:ncell)+pstate%atm_in_tauy_at_pc_surface%f(1:ncell)

    diag_physics_vars_h0%rh%f(:,1:ncell)    = diag_physics_vars_h0%rh%f(:,1:ncell)+pstate_cam%diag_relhum%f(:,1:ncell)
    diag_physics_vars_h0%cld%f(:,1:ncell)   = diag_physics_vars_h0%cld%f(:,1:ncell)+pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncell)
    diag_physics_vars_h0%concld%f(:,1:ncell)   = diag_physics_vars_h0%concld%f(:,1:ncell)+pstate_cam%macrop_concld_at_pc_full_level%f(old_time_level,:,1:ncell)
    diag_physics_vars_h0%lwc%f(:,1:ncell)   = diag_physics_vars_h0%lwc%f(:,1:ncell)+pstate_cam%microp_lwc_at_pc_full_level%f(:,1:ncell)
    diag_physics_vars_h0%iwc%f(:,1:ncell)   = diag_physics_vars_h0%iwc%f(:,1:ncell)+pstate_cam%microp_iwc_at_pc_full_level%f(:,1:ncell)
    diag_physics_vars_h0%cldtot%f(1:ncell)= diag_physics_vars_h0%cldtot%f(1:ncell)+pstate_cam%diag_cloud_tot%f(1:ncell)
    diag_physics_vars_h0%cldlow%f(1:ncell)= diag_physics_vars_h0%cldlow%f(1:ncell)+pstate_cam%diag_cloud_low%f(1:ncell)
    diag_physics_vars_h0%cldmed%f(1:ncell)= diag_physics_vars_h0%cldmed%f(1:ncell)+pstate_cam%diag_cloud_med%f(1:ncell)
    diag_physics_vars_h0%cldhgh%f(1:ncell)= diag_physics_vars_h0%cldhgh%f(1:ncell)+pstate_cam%diag_cloud_hgh%f(1:ncell)
    diag_physics_vars_h0%z500%f(1:ncell)  = diag_physics_vars_h0%z500%f(1:ncell)+pstate_cam%diag_z_at_500hpa%f(1:ncell)
    diag_physics_vars_h0%u850%f(1:ncell)  = diag_physics_vars_h0%u850%f(1:ncell)+pstate_cam%diag_u_at_850hpa%f(1:ncell)
    diag_physics_vars_h0%u200%f(1:ncell)  = diag_physics_vars_h0%u200%f(1:ncell)+pstate_cam%diag_u_at_200hpa%f(1:ncell)
    diag_physics_vars_h0%v850%f(1:ncell)  = diag_physics_vars_h0%v850%f(1:ncell)+pstate_cam%diag_v_at_850hpa%f(1:ncell)
 
    diag_physics_vars_h0%flwut%f(1:ncell) = diag_physics_vars_h0%flwut%f(1:ncell)+pstate_cam%flwut_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%flwdt%f(1:ncell) = diag_physics_vars_h0%flwdt%f(1:ncell)+pstate_cam%flwdt_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%flwus%f(1:ncell) = diag_physics_vars_h0%flwus%f(1:ncell)+pstate_cam%flwus_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%flwds%f(1:ncell) = diag_physics_vars_h0%flwds%f(1:ncell)+pstate_cam%flwds_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%flwutc%f(1:ncell)= diag_physics_vars_h0%flwutc%f(1:ncell)+pstate_cam%flwutc_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%flwdtc%f(1:ncell)= diag_physics_vars_h0%flwdtc%f(1:ncell)+pstate_cam%flwdtc_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%flwusc%f(1:ncell)= diag_physics_vars_h0%flwusc%f(1:ncell)+pstate_cam%flwusc_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%flwdsc%f(1:ncell)= diag_physics_vars_h0%flwdsc%f(1:ncell)+pstate_cam%flwdsc_at_pc_surface%f(1:ncell)

    diag_physics_vars_h0%fswut%f(1:ncell) = diag_physics_vars_h0%fswut%f(1:ncell)+pstate_cam%fswut_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%fswdt%f(1:ncell) = diag_physics_vars_h0%fswdt%f(1:ncell)+pstate_cam%fswdt_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%fswus%f(1:ncell) = diag_physics_vars_h0%fswus%f(1:ncell)+pstate_cam%fswus_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%fswds%f(1:ncell) = diag_physics_vars_h0%fswds%f(1:ncell)+pstate_cam%fswds_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%fsntoa%f(1:ncell)= diag_physics_vars_h0%fsntoa%f(1:ncell)+pstate_cam%fsntoa_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%fswutc%f(1:ncell)= diag_physics_vars_h0%fswutc%f(1:ncell)+pstate_cam%fswutc_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%fswdtc%f(1:ncell)= diag_physics_vars_h0%fswdtc%f(1:ncell)+pstate_cam%fswdtc_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%fswusc%f(1:ncell)= diag_physics_vars_h0%fswusc%f(1:ncell)+pstate_cam%fswusc_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%fswdsc%f(1:ncell)= diag_physics_vars_h0%fswdsc%f(1:ncell)+pstate_cam%fswdsc_at_pc_surface%f(1:ncell)
    diag_physics_vars_h0%fsntoac%f(1:ncell)=diag_physics_vars_h0%fsntoac%f(1:ncell)+pstate_cam%fsntoac_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%lwcf%f(1:ncell)  = diag_physics_vars_h0%lwcf%f(1:ncell)+pstate_cam%lwcf_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%swcf%f(1:ncell)  = diag_physics_vars_h0%swcf%f(1:ncell)+pstate_cam%swcf_at_pc_top%f(1:ncell)
    diag_physics_vars_h0%tmq%f(1:ncell)   = diag_physics_vars_h0%tmq%f(1:ncell)+pstate_cam%diag_tmq%f(1:ncell)
    diag_physics_vars_h0%lwp%f(1:ncell)   = diag_physics_vars_h0%lwp%f(1:ncell)+pstate_cam%diag_tgliqwp%f(1:ncell)
    diag_physics_vars_h0%iwp%f(1:ncell)   = diag_physics_vars_h0%iwp%f(1:ncell)+pstate_cam%diag_tgicewp%f(1:ncell)

    if(allocated(diag_physics_vars_h0%sgm%f))then
    diag_physics_vars_h0%sgm%f(:,1:ncell) = diag_physics_vars_h0%sgm%f(:,1:ncell)+pstate_cam%macrop_sgm_at_pc_full_level%f(:,1:ncell)
    end if
#endif

#ifdef AMIPW_PHYSICS
    diag_physics_vars_h0%qc%f(:,1:ncell)    = diag_physics_vars_h0%qc%f(:,1:ncell)+tracerVarCellFull%scalar_tracer_mxrt_n%f(2,:,1:ncell)
    diag_physics_vars_h0%qi%f(:,1:ncell)    = diag_physics_vars_h0%qi%f(:,1:ncell)+tracerVarCellFull%scalar_tracer_mxrt_n%f(4,:,1:ncell)
    diag_physics_vars_h0%snowh%f(1:ncell) = diag_physics_vars_h0%snowh%f(1:ncell)+psurf_wrf%snow(1:ncell,1)
    diag_physics_vars_h0%taux%f(1:ncell)  = diag_physics_vars_h0%taux%f(1:ncell)+psurf_wrf%taux%f(1:ncell)
    diag_physics_vars_h0%tauy%f(1:ncell)  = diag_physics_vars_h0%tauy%f(1:ncell)+psurf_wrf%tauy%f(1:ncell)
    diag_physics_vars_h0%cldtot%f(1:ncell)= diag_physics_vars_h0%cldtot%f(1:ncell)+pstate_wrf%cldtot(1:ncell,1)
    diag_physics_vars_h0%cldlow%f(1:ncell)= diag_physics_vars_h0%cldlow%f(1:ncell)+pstate_wrf%cldlow(1:ncell,1)
    diag_physics_vars_h0%cldmed%f(1:ncell)= diag_physics_vars_h0%cldmed%f(1:ncell)+pstate_wrf%cldmed(1:ncell,1)
    diag_physics_vars_h0%cldhgh%f(1:ncell)= diag_physics_vars_h0%cldhgh%f(1:ncell)+pstate_wrf%cldhgh(1:ncell,1)
    pres_out(1) = 20000._r8
    pres_out(2) = 50000._r8
    pres_out(3) = 85000._r8
    do icell = 1, ncell
        call lininterp(dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,icell),&
                       log(log(dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,icell))),1,nlev,vartmp_out(1:3),log(log(pres_out)),3)
        diag_physics_vars_h0%z500%f(icell) = diag_physics_vars_h0%z500%f(icell)+vartmp_out(2)
        call lininterp(dycoreVarCellFull%scalar_U_wind_n%f(1:nlev,icell),&
                       log(log(dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,icell))),1,nlev,vartmp_out(1:3),log(log(pres_out)),3) 
        diag_physics_vars_h0%u200%f(icell) = diag_physics_vars_h0%u200%f(icell)+vartmp_out(1)
        diag_physics_vars_h0%u850%f(icell) = diag_physics_vars_h0%u850%f(icell)+vartmp_out(3)
        call lininterp(dycoreVarCellFull%scalar_V_wind_n%f(1:nlev,icell),&
                       log(log(dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,icell))),1,nlev,vartmp_out(1:3),log(log(pres_out)),3) 
        diag_physics_vars_h0%v850%f(icell) = diag_physics_vars_h0%v850%f(icell)+vartmp_out(3)
    end do
    diag_physics_vars_h0%flwut%f(1:ncell) = diag_physics_vars_h0%flwut%f(1:ncell)+pstate_wrf%lwupt(1:ncell,1)
    diag_physics_vars_h0%flwdt%f(1:ncell) = diag_physics_vars_h0%flwdt%f(1:ncell)+pstate_wrf%lwdnt(1:ncell,1)
    diag_physics_vars_h0%flwus%f(1:ncell) = diag_physics_vars_h0%flwus%f(1:ncell)+pstate_wrf%lwupb(1:ncell,1)
    diag_physics_vars_h0%flwds%f(1:ncell) = diag_physics_vars_h0%flwds%f(1:ncell)+pstate_wrf%lwdnb(1:ncell,1)
    diag_physics_vars_h0%flwutc%f(1:ncell)= diag_physics_vars_h0%flwutc%f(1:ncell)+pstate_wrf%lwuptc(1:ncell,1)
    diag_physics_vars_h0%flwdtc%f(1:ncell)= diag_physics_vars_h0%flwdtc%f(1:ncell)+pstate_wrf%lwdntc(1:ncell,1)
    diag_physics_vars_h0%flwusc%f(1:ncell)= diag_physics_vars_h0%flwusc%f(1:ncell)+pstate_wrf%lwupbc(1:ncell,1)
    diag_physics_vars_h0%flwdsc%f(1:ncell)= diag_physics_vars_h0%flwdsc%f(1:ncell)+pstate_wrf%lwdnbc(1:ncell,1)

    diag_physics_vars_h0%fswut%f(1:ncell) = diag_physics_vars_h0%fswut%f(1:ncell)+pstate_wrf%swupt(1:ncell,1)
    diag_physics_vars_h0%fswdt%f(1:ncell) = diag_physics_vars_h0%fswdt%f(1:ncell)+pstate_wrf%swdnt(1:ncell,1)
    diag_physics_vars_h0%fswus%f(1:ncell) = diag_physics_vars_h0%fswus%f(1:ncell)+pstate_wrf%swupb(1:ncell,1)
    diag_physics_vars_h0%fswds%f(1:ncell) = diag_physics_vars_h0%fswds%f(1:ncell)+pstate_wrf%swdnb(1:ncell,1)
    diag_physics_vars_h0%fsntoa%f(1:ncell)= diag_physics_vars_h0%fsntoa%f(1:ncell)+(pstate_wrf%swdnt(1:ncell,1)-pstate_wrf%swupt(1:ncell,1))
    diag_physics_vars_h0%fswutc%f(1:ncell)= diag_physics_vars_h0%fswutc%f(1:ncell)+pstate_wrf%swuptc(1:ncell,1)
    diag_physics_vars_h0%fswdtc%f(1:ncell)= diag_physics_vars_h0%fswdtc%f(1:ncell)+pstate_wrf%swdntc(1:ncell,1)
    diag_physics_vars_h0%fswusc%f(1:ncell)= diag_physics_vars_h0%fswusc%f(1:ncell)+pstate_wrf%swupbc(1:ncell,1)
    diag_physics_vars_h0%fswdsc%f(1:ncell)= diag_physics_vars_h0%fswdsc%f(1:ncell)+pstate_wrf%swdnbc(1:ncell,1)
    diag_physics_vars_h0%fsntoac%f(1:ncell)=diag_physics_vars_h0%fsntoac%f(1:ncell)+(pstate_wrf%swdntc(1:ncell,1)-pstate_wrf%swuptc(1:ncell,1))
    diag_physics_vars_h0%lwcf%f(1:ncell)  = diag_physics_vars_h0%lwcf%f(1:ncell)+pstate_wrf%lwcf(1:ncell,1)
    diag_physics_vars_h0%swcf%f(1:ncell)  = diag_physics_vars_h0%swcf%f(1:ncell)+pstate_wrf%swcf(1:ncell,1)

    do icell = 1, ncell
    do ilev = 1, nlev
       dpres(ilev)  = dycoreVarCellFace%scalar_mpressure_n%f(ilev+1,icell)-dycoreVarCellFace%scalar_mpressure_n%f(ilev,icell)
    end do
    ftem(1:nlev)    = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,icell)*dpres(1:nlev)/gravity
    ftem_qc(1:nlev) = tracerVarCellFull%scalar_tracer_mxrt_n%f(2,1:nlev,icell)*dpres(1:nlev)/gravity
    ftem_qi(1:nlev) = tracerVarCellFull%scalar_tracer_mxrt_n%f(4,1:nlev,icell)*dpres(1:nlev)/gravity
    do ilev = 2, nlev
       ftem(1) = ftem(1)+ftem(ilev) 
       ftem_qc(1) = ftem_qc(1)+ftem_qc(ilev) 
       ftem_qi(1) = ftem_qi(1)+ftem_qi(ilev) 
    end do
    diag_physics_vars_h0%tmq%f(icell)   = diag_physics_vars_h0%tmq%f(icell)+ftem(1)
    diag_physics_vars_h0%lwp%f(icell)   = diag_physics_vars_h0%lwp%f(icell)+ftem_qc(1)
    diag_physics_vars_h0%iwp%f(icell)   = diag_physics_vars_h0%iwp%f(icell)+ftem_qi(1)
    end do

    call aqsat_grist(t = dycoreVarCellFull%scalar_temp_n%f     (1:nlev,1:ncell), &
                     p = dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,1:ncell), &
                     es= esvp(1:nlev,1:ncell), &
                     qs= sat_specific_humidity(1:nlev,1:ncell), &
                     ii= ncell,ilen=ncell,kk=nlev,kstart=1,kend=nlev)
    do icell = 1, ncell
    do ilev = 1, nlev
    diag_physics_vars_h0%lwc%f(ilev,icell)   = diag_physics_vars_h0%lwc%f(ilev,icell)+&
                                               tracerVarCellFull%scalar_tracer_mxrt_n%f(2,ilev,icell)*dycoreVarCellFull%scalar_mpressure_n%f(ilev,icell)&
                                               /(287.15_r8*dycoreVarCellFull%scalar_temp_n%f(ilev,icell))
    diag_physics_vars_h0%iwc%f(ilev,icell)   = diag_physics_vars_h0%iwc%f(ilev,icell)+&
                                                tracerVarCellFull%scalar_tracer_mxrt_n%f(4,ilev,icell)*dycoreVarCellFull%scalar_mpressure_n%f(ilev,icell)&
                                               /(287.15_r8*dycoreVarCellFull%scalar_temp_n%f(ilev,icell))
 
    diag_physics_vars_h0%concld%f(ilev,icell)= diag_physics_vars_h0%concld%f(ilev,icell)+pstate_wrf%cldcum(icell,nlev+1-ilev,1)
    diag_physics_vars_h0%cld%f(ilev,icell)   = diag_physics_vars_h0%cld%f(ilev,icell)+pstate_wrf%cldfra(icell,nlev+1-ilev,1)
    diag_physics_vars_h0%rh%f(ilev,icell)    = diag_physics_vars_h0%rh%f(ilev,icell)+&
                                               tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,icell)&
                                               /(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,icell)))/sat_specific_humidity(ilev,icell)*100._r8
    end do
    end do
#endif
    end if

    diag_physics_vars_h0%ncount  = diag_physics_vars_h0%ncount + 1
    return
  end subroutine gcm_h0_accu_physics_variables

!
! dump depending on write_history frequency
!

  subroutine gcm_h0_dump_physics_variables
! local
    if(write_history_h0)then

    diag_physics_vars_h0%precl%f = diag_physics_vars_h0%precl%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%precc%f = diag_physics_vars_h0%precc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%prect%f = diag_physics_vars_h0%prect%f/diag_physics_vars_h0%ncount

    diag_physics_vars_h0%ps%f    = diag_physics_vars_h0%ps%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%snowh%f = diag_physics_vars_h0%snowh%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%ts%f    = diag_physics_vars_h0%ts%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%shflx%f = diag_physics_vars_h0%shflx%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%qflx%f  = diag_physics_vars_h0%qflx%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%phi%f   = diag_physics_vars_h0%phi%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%uwind%f = diag_physics_vars_h0%uwind%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%vwind%f = diag_physics_vars_h0%vwind%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%omega%f = diag_physics_vars_h0%omega%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%temp%f  = diag_physics_vars_h0%temp%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%qv%f    = diag_physics_vars_h0%qv%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%qc%f    = diag_physics_vars_h0%qc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%qi%f    = diag_physics_vars_h0%qi%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%rh%f    = diag_physics_vars_h0%rh%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%cld%f   = diag_physics_vars_h0%cld%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%concld%f   = diag_physics_vars_h0%concld%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%lwc%f   = diag_physics_vars_h0%lwc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%iwc%f   = diag_physics_vars_h0%iwc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%cldtot%f= diag_physics_vars_h0%cldtot%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%cldlow%f= diag_physics_vars_h0%cldlow%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%cldmed%f= diag_physics_vars_h0%cldmed%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%cldhgh%f= diag_physics_vars_h0%cldhgh%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%z500%f  = diag_physics_vars_h0%z500%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%u850%f  = diag_physics_vars_h0%u850%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%u200%f  = diag_physics_vars_h0%u200%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%v850%f  = diag_physics_vars_h0%v850%f/diag_physics_vars_h0%ncount

    diag_physics_vars_h0%flwut%f = diag_physics_vars_h0%flwut%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%flwdt%f = diag_physics_vars_h0%flwdt%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%flwus%f = diag_physics_vars_h0%flwus%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%flwds%f = diag_physics_vars_h0%flwds%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%flwutc%f= diag_physics_vars_h0%flwutc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%flwdtc%f= diag_physics_vars_h0%flwdtc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%flwusc%f= diag_physics_vars_h0%flwusc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%flwdsc%f= diag_physics_vars_h0%flwdsc%f/diag_physics_vars_h0%ncount

    diag_physics_vars_h0%fswut%f = diag_physics_vars_h0%fswut%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fswdt%f = diag_physics_vars_h0%fswdt%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fswus%f = diag_physics_vars_h0%fswus%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fswds%f = diag_physics_vars_h0%fswds%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fsntoa%f= diag_physics_vars_h0%fsntoa%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fswutc%f= diag_physics_vars_h0%fswutc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fswdtc%f= diag_physics_vars_h0%fswdtc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fswusc%f= diag_physics_vars_h0%fswusc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fswdsc%f= diag_physics_vars_h0%fswdsc%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%fsntoac%f=diag_physics_vars_h0%fsntoac%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%lwcf%f  = diag_physics_vars_h0%lwcf%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%swcf%f  = diag_physics_vars_h0%swcf%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%tmq%f   = diag_physics_vars_h0%tmq%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%lwp%f   = diag_physics_vars_h0%lwp%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%iwp%f   = diag_physics_vars_h0%iwp%f/diag_physics_vars_h0%ncount

    diag_physics_vars_h0%taux%f   = diag_physics_vars_h0%taux%f/diag_physics_vars_h0%ncount
    diag_physics_vars_h0%tauy%f   = diag_physics_vars_h0%tauy%f/diag_physics_vars_h0%ncount

    if(allocated(diag_physics_vars_h0%sgm%f))then
    diag_physics_vars_h0%sgm%f   = diag_physics_vars_h0%sgm%f/diag_physics_vars_h0%ncount
    end if

    end if
 
    return
  end subroutine gcm_h0_dump_physics_variables

  subroutine gcm_h0_rest_physics_variables
! reset to zero
    diag_physics_vars_h0%ncount  = 0

    if(write_history_h0)then
    diag_physics_vars_h0%precl%f = zero
    diag_physics_vars_h0%precc%f = zero
    diag_physics_vars_h0%prect%f = zero

    diag_physics_vars_h0%ps%f    = zero 
    diag_physics_vars_h0%ts%f    = zero 
    diag_physics_vars_h0%snowh%f = zero 
    diag_physics_vars_h0%shflx%f = zero 
    diag_physics_vars_h0%qflx%f  = zero 
    diag_physics_vars_h0%phi%f   = zero 
    diag_physics_vars_h0%uwind%f = zero 
    diag_physics_vars_h0%vwind%f = zero 
    diag_physics_vars_h0%omega%f = zero 
    diag_physics_vars_h0%temp%f  = zero 
    diag_physics_vars_h0%qv%f    = zero 
    diag_physics_vars_h0%qc%f    = zero 
    diag_physics_vars_h0%qi%f    = zero 
    diag_physics_vars_h0%rh%f    = zero 
    diag_physics_vars_h0%cld%f   = zero 
    diag_physics_vars_h0%concld%f   = zero 
    diag_physics_vars_h0%lwc%f   = zero 
    diag_physics_vars_h0%iwc%f   = zero 
    diag_physics_vars_h0%z500%f  = zero 
    diag_physics_vars_h0%u850%f  = zero 
    diag_physics_vars_h0%u200%f  = zero 
    diag_physics_vars_h0%v850%f  = zero 

    diag_physics_vars_h0%flwut%f = zero 
    diag_physics_vars_h0%flwdt%f = zero 
    diag_physics_vars_h0%flwus%f = zero 
    diag_physics_vars_h0%flwds%f = zero 
    diag_physics_vars_h0%flwutc%f= zero 
    diag_physics_vars_h0%flwdtc%f= zero 
    diag_physics_vars_h0%flwusc%f= zero 
    diag_physics_vars_h0%flwdsc%f= zero 
    diag_physics_vars_h0%fswut%f = zero 
    diag_physics_vars_h0%fswdt%f = zero 
    diag_physics_vars_h0%fswus%f = zero 
    diag_physics_vars_h0%fswds%f = zero 
    diag_physics_vars_h0%fsntoa%f= zero 
    diag_physics_vars_h0%fswutc%f= zero 
    diag_physics_vars_h0%fswdtc%f= zero 
    diag_physics_vars_h0%fswusc%f= zero 
    diag_physics_vars_h0%fswdsc%f= zero 
    diag_physics_vars_h0%fsntoac%f=zero 
    diag_physics_vars_h0%lwcf%f  = zero 
    diag_physics_vars_h0%swcf%f  = zero 
    diag_physics_vars_h0%tmq%f   = zero 
    diag_physics_vars_h0%lwp%f   = zero 
    diag_physics_vars_h0%iwp%f   = zero 
    diag_physics_vars_h0%taux%f  = zero 
    diag_physics_vars_h0%tauy%f  = zero 

    if(allocated(diag_physics_vars_h0%sgm%f))diag_physics_vars_h0%sgm%f   = zero
    end if
 
    return
  end subroutine gcm_h0_rest_physics_variables

!----------------------------------------------------------
! private wrap routines
!----------------------------------------------------------

    subroutine wrap_allocate_data1d(mesh,var)
       type(global_domain),   intent(in)    :: mesh
       type(scalar_1d_field), intent(inout) :: var
       if(.not.allocated(var%f)) allocate(var%f(mesh%nv))
       var%f    = zero
       var%pos  = 0
       return
    end subroutine wrap_allocate_data1d

    subroutine wrap_allocate_data2d(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4),           intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var
       if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%nv))
       var%f    = zero
       var%pos  = 0
       return
    end subroutine wrap_allocate_data2d

    subroutine wrap_allocate_data3d(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4),           intent(in)    :: nLevel
       type(scalar_3d_field), intent(inout) :: var
       if(.not.allocated(var%f)) allocate(var%f(ntracer,nLevel,mesh%nv))
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

  end module grist_gcm_diagnose_h0_module
