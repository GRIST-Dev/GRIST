
!----------------------------------------------------------------------------
! Created on 2017-5-17
! Author: Yi Zhang
! Version 1.0
! Description: GCM Diagnostic module
! Revision history: based on the dycore diagnose module
!                   currently used for accumulating vars like precip
!----------------------------------------------------------------------------

 module grist_gcm_diagnose_h1_module

  use grist_constants,              only: i4, r8, zero 
  use grist_domain_types,           only: global_domain
  use grist_data_types,             only: scalar_1d_field,  &
                                          scalar_2d_field
  use grist_dycore_vars_module
  use grist_tracer_transport_vars_module
  use grist_physics_data_structure, only: pstate
  use grist_nml_module,             only: nlev, nlevp,      &
                                          write_history_h1
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,   only: pstate_cam
  use grist_physics_update,         only: old_time_level
#endif

  implicit none

   private

   public   :: gcm_h1_diagnose_init,  &
               gcm_h1_diagnose_final, & 
               gcm_h1_accu_physics_variables, &
               gcm_h1_dump_physics_variables, &
               gcm_h1_rest_physics_variables, &
               diag_physics_vars_h1

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
        type(scalar_2d_field)  :: rh
        type(scalar_2d_field)  :: cld
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
        type(scalar_1d_field)  :: fswutc
        type(scalar_1d_field)  :: fswdtc
        type(scalar_1d_field)  :: fswusc
        type(scalar_1d_field)  :: fswdsc
        type(scalar_1d_field)  :: lwcf
        type(scalar_1d_field)  :: swcf

        type(scalar_1d_field)  :: asdir
        type(scalar_1d_field)  :: asdif
        type(scalar_1d_field)  :: aldir
        type(scalar_1d_field)  :: aldif
 
        integer(i4)            :: ncount
   end type gcm_diag_vars

   type(gcm_diag_vars)  :: diag_physics_vars_h1
   integer(i4) :: ncell

  contains

  subroutine gcm_h1_diagnose_init(mesh)
    type(global_domain),  intent(in)   :: mesh

    diag_physics_vars_h1%ncount  = 0

    if(write_history_h1)then
!    call wrap_allocate_data1d(mesh,diag_physics_vars_h1%precc)
!    call wrap_allocate_data1d(mesh,diag_physics_vars_h1%precl)
    call wrap_allocate_data1d(mesh,diag_physics_vars_h1%prect)

!    call wrap_allocate_data1d(mesh,     diag_physics_vars%ps)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%ts)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%snowh)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%shflx)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%qflx)
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%phi)
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%uwind)
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%vwind)
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%omega)
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%temp)
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%qv)

#ifdef AMIPC_PHYSICS
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%rh)
!    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars%cld)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%cldtot)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%cldlow)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%cldmed)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%cldhgh)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%z500)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars_h1%u850)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars_h1%u200)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars_h1%v850)

    call wrap_allocate_data1d(mesh,     diag_physics_vars_h1%flwut)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%flwdt)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%flwus)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%flwds)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%flwutc)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%flwdtc)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%flwusc)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%flwdsc)

!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswut)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswdt)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswus)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswds)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswutc)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswdtc)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswusc)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%fswdsc)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%lwcf)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%swcf)

!    call wrap_allocate_data1d(mesh,     diag_physics_vars%asdir)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%asdif)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%aldir)
!    call wrap_allocate_data1d(mesh,     diag_physics_vars%aldif)
#endif 
    end if

    ncell = mesh%nv_halo(1)

    return
  end subroutine gcm_h1_diagnose_init

  subroutine gcm_h1_diagnose_final
!    if(allocated(diag_physics_vars_h1%precc%f))  deallocate(diag_physics_vars_h1%precc%f)
!    if(allocated(diag_physics_vars_h1%precl%f))  deallocate(diag_physics_vars_h1%precl%f)
    if(allocated(diag_physics_vars_h1%prect%f))  deallocate(diag_physics_vars_h1%prect%f)

!    if(allocated(diag_physics_vars%ps%f))     deallocate(diag_physics_vars%ps%f)
!    if(allocated(diag_physics_vars%snowh%f))  deallocate(diag_physics_vars%snowh%f)
!    if(allocated(diag_physics_vars%ts%f))     deallocate(diag_physics_vars%ts%f)
!    if(allocated(diag_physics_vars%shflx%f))  deallocate(diag_physics_vars%shflx%f)
!    if(allocated(diag_physics_vars%qflx%f))   deallocate(diag_physics_vars%qflx%f)
!    if(allocated(diag_physics_vars%phi%f))    deallocate(diag_physics_vars%phi%f)
!    if(allocated(diag_physics_vars%uwind%f))  deallocate(diag_physics_vars%uwind%f)
!    if(allocated(diag_physics_vars%vwind%f))  deallocate(diag_physics_vars%vwind%f)
!    if(allocated(diag_physics_vars%omega%f))  deallocate(diag_physics_vars%omega%f)
!    if(allocated(diag_physics_vars%temp%f))   deallocate(diag_physics_vars%temp%f)
!    if(allocated(diag_physics_vars%qv%f))     deallocate(diag_physics_vars%qv%f)
!    if(allocated(diag_physics_vars%rh%f))     deallocate(diag_physics_vars%rh%f)
!    if(allocated(diag_physics_vars%cld%f))    deallocate(diag_physics_vars%cld%f)
!    if(allocated(diag_physics_vars%cldtot%f)) deallocate(diag_physics_vars%cldtot%f)
!    if(allocated(diag_physics_vars%cldlow%f)) deallocate(diag_physics_vars%cldlow%f)
!    if(allocated(diag_physics_vars%cldmed%f)) deallocate(diag_physics_vars%cldmed%f)
!    if(allocated(diag_physics_vars%cldhgh%f)) deallocate(diag_physics_vars%cldhgh%f)
!    if(allocated(diag_physics_vars%z500%f))   deallocate(diag_physics_vars%z500%f)
!    if(allocated(diag_physics_vars_h1%u850%f))   deallocate(diag_physics_vars_h1%u850%f)
!    if(allocated(diag_physics_vars_h1%u200%f))   deallocate(diag_physics_vars_h1%u200%f)
!    if(allocated(diag_physics_vars_h1%v850%f))   deallocate(diag_physics_vars_h1%v850%f)

    if(allocated(diag_physics_vars_h1%flwut%f))  deallocate(diag_physics_vars_h1%flwut%f )
!    if(allocated(diag_physics_vars%flwdt%f))  deallocate(diag_physics_vars%flwdt%f )
!    if(allocated(diag_physics_vars%flwus%f))  deallocate(diag_physics_vars%flwus%f )
!    if(allocated(diag_physics_vars%flwds%f))  deallocate(diag_physics_vars%flwds%f )
!    if(allocated(diag_physics_vars%flwutc%f)) deallocate(diag_physics_vars%flwutc%f)
!    if(allocated(diag_physics_vars%flwdtc%f)) deallocate(diag_physics_vars%flwdtc%f)
!    if(allocated(diag_physics_vars%flwusc%f)) deallocate(diag_physics_vars%flwusc%f)
!    if(allocated(diag_physics_vars%flwdsc%f)) deallocate(diag_physics_vars%flwdsc%f)
!    if(allocated(diag_physics_vars%fswut%f))  deallocate(diag_physics_vars%fswut%f )
!    if(allocated(diag_physics_vars%fswdt%f))  deallocate(diag_physics_vars%fswdt%f )
!    if(allocated(diag_physics_vars%fswus%f))  deallocate(diag_physics_vars%fswus%f )
!    if(allocated(diag_physics_vars%fswds%f))  deallocate(diag_physics_vars%fswds%f )
!    if(allocated(diag_physics_vars%fswutc%f)) deallocate(diag_physics_vars%fswutc%f)
!    if(allocated(diag_physics_vars%fswdtc%f)) deallocate(diag_physics_vars%fswdtc%f)
!    if(allocated(diag_physics_vars%fswusc%f)) deallocate(diag_physics_vars%fswusc%f)
!    if(allocated(diag_physics_vars%fswdsc%f)) deallocate(diag_physics_vars%fswdsc%f)
!    if(allocated(diag_physics_vars%lwcf%f))   deallocate(diag_physics_vars%lwcf%f)
!    if(allocated(diag_physics_vars%swcf%f))   deallocate(diag_physics_vars%swcf%f)

!    if(allocated(diag_physics_vars%asdir%f))   deallocate(diag_physics_vars%asdir%f)
!    if(allocated(diag_physics_vars%asdif%f))   deallocate(diag_physics_vars%asdif%f)
!    if(allocated(diag_physics_vars%aldir%f))   deallocate(diag_physics_vars%aldir%f)
!    if(allocated(diag_physics_vars%aldif%f))   deallocate(diag_physics_vars%aldif%f)
 

    return
  end subroutine gcm_h1_diagnose_final
!
! accumulate at each model step
!
  subroutine gcm_h1_accu_physics_variables
! local
    if(write_history_h1)then

!    diag_physics_vars_h1%precl%f(1:ncell) = diag_physics_vars_h1%precl%f(1:ncell)+pstate%scalar_precl_surface%f(1:ncell)
!    diag_physics_vars_h1%precc%f(1:ncell) = diag_physics_vars_h1%precc%f(1:ncell)+pstate%scalar_precc_surface%f(1:ncell)
    diag_physics_vars_h1%prect%f(1:ncell) = diag_physics_vars_h1%prect%f(1:ncell)+pstate%scalar_prect_surface%f(1:ncell)

!    diag_physics_vars%ps%f    = diag_physics_vars%ps%f+scalar_pressure_at_pc_surface_n%f
!    diag_physics_vars%ts%f    = diag_physics_vars%ts%f+pstate%ts_at_pc_surface%f
!    diag_physics_vars%snowh%f = diag_physics_vars%snowh%f+pstate%snowhland_at_pc_surface%f
!    diag_physics_vars%shflx%f = diag_physics_vars%shflx%f+pstate%atm_in_shflx_at_pc_surface%f
!    diag_physics_vars%qflx%f  = diag_physics_vars%qflx%f+pstate%atm_in_qflx_at_pc_surface%f(1,:)
!    diag_physics_vars%phi%f   = diag_physics_vars%phi%f+scalar_geopotential_at_pc_full_level_n%f
!    diag_physics_vars%uwind%f = diag_physics_vars%uwind%f+scalar_U_wind_at_pc_full_level_n%f 
!    diag_physics_vars%vwind%f = diag_physics_vars%vwind%f+scalar_V_wind_at_pc_full_level_n%f
!    diag_physics_vars%omega%f = diag_physics_vars%omega%f+scalar_omega_at_pc_full_level_n%f
!    diag_physics_vars%temp%f  = diag_physics_vars%temp%f+scalar_temp_at_pc_full_level_n%f
!    diag_physics_vars%qv%f    = diag_physics_vars%qv%f+scalar_tracer_mxrt_at_pc_full_level_n%f(1,:,:)
#ifdef AMIPC_PHYSICS
!    diag_physics_vars%rh%f    = diag_physics_vars%rh%f+pstate_cam%diag_relhum%f
!    diag_physics_vars%cld%f   = diag_physics_vars%cld%f+pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,:)
!    diag_physics_vars%cldtot%f= diag_physics_vars%cldtot%f+pstate_cam%diag_cloud_tot%f(:)
!    diag_physics_vars%cldlow%f= diag_physics_vars%cldlow%f+pstate_cam%diag_cloud_low%f(:)
!    diag_physics_vars%cldmed%f= diag_physics_vars%cldmed%f+pstate_cam%diag_cloud_med%f(:)
!    diag_physics_vars%cldhgh%f= diag_physics_vars%cldhgh%f+pstate_cam%diag_cloud_hgh%f(:)
!    diag_physics_vars%z500%f  = diag_physics_vars%z500%f+pstate_cam%diag_z_at_500hpa%f
!    diag_physics_vars_h1%u850%f(1:ncell)  = diag_physics_vars_h1%u850%f(1:ncell)+pstate_cam%diag_u_at_850hpa%f(1:ncell)
!    diag_physics_vars_h1%u200%f(1:ncell)  = diag_physics_vars_h1%u200%f(1:ncell)+pstate_cam%diag_u_at_200hpa%f(1:ncell)
!    diag_physics_vars_h1%v850%f(1:ncell)  = diag_physics_vars_h1%v850%f(1:ncell)+pstate_cam%diag_v_at_850hpa%f(1:ncell)
 
    diag_physics_vars_h1%flwut%f(1:ncell) = diag_physics_vars_h1%flwut%f(1:ncell)+pstate_cam%flwut_at_pc_top%f(1:ncell)
!    diag_physics_vars%flwdt%f(1:ncell) = diag_physics_vars%flwdt%f(1:ncell)+pstate_cam%flwdt_at_pc_top%f(1:ncell)
!    diag_physics_vars%flwus%f(1:ncell) = diag_physics_vars%flwus%f(1:ncell)+pstate_cam%flwus_at_pc_surface%f(1:ncell)
!    diag_physics_vars%flwds%f(1:ncell) = diag_physics_vars%flwds%f(1:ncell)+pstate_cam%flwds_at_pc_surface%f(1:ncell)
!    diag_physics_vars%flwutc%f(1:ncell)= diag_physics_vars%flwutc%f(1:ncell)+pstate_cam%flwutc_at_pc_top%f(1:ncell)
!    diag_physics_vars%flwdtc%f(1:ncell)= diag_physics_vars%flwdtc%f(1:ncell)+pstate_cam%flwdtc_at_pc_top%f(1:ncell)
!    diag_physics_vars%flwusc%f(1:ncell)= diag_physics_vars%flwusc%f(1:ncell)+pstate_cam%flwusc_at_pc_surface%f(1:ncell)
!    diag_physics_vars%flwdsc%f(1:ncell)= diag_physics_vars%flwdsc%f(1:ncell)+pstate_cam%flwdsc_at_pc_surface%f(1:ncell)
!
!    diag_physics_vars%fswut%f(1:ncell) = diag_physics_vars%fswut%f(1:ncell)+pstate_cam%fswut_at_pc_top%f(1:ncell)
!    diag_physics_vars%fswdt%f(1:ncell) = diag_physics_vars%fswdt%f(1:ncell)+pstate_cam%fswdt_at_pc_top%f(1:ncell)
!    diag_physics_vars%fswus%f(1:ncell) = diag_physics_vars%fswus%f(1:ncell)+pstate_cam%fswus_at_pc_surface%f(1:ncell)
!    diag_physics_vars%fswds%f(1:ncell) = diag_physics_vars%fswds%f(1:ncell)+pstate_cam%fswds_at_pc_surface%f(1:ncell)
!    diag_physics_vars%fswutc%f(1:ncell)= diag_physics_vars%fswutc%f(1:ncell)+pstate_cam%fswutc_at_pc_top%f(1:ncell)
!    diag_physics_vars%fswdtc%f(1:ncell)= diag_physics_vars%fswdtc%f(1:ncell)+pstate_cam%fswdtc_at_pc_top%f(1:ncell)
!    diag_physics_vars%fswusc%f(1:ncell)= diag_physics_vars%fswusc%f(1:ncell)+pstate_cam%fswusc_at_pc_surface%f(1:ncell)
!    diag_physics_vars%fswdsc%f(1:ncell)= diag_physics_vars%fswdsc%f(1:ncell)+pstate_cam%fswdsc_at_pc_surface%f(1:ncell)
!    diag_physics_vars%lwcf%f(1:ncell)  = diag_physics_vars%lwcf%f(1:ncell)+pstate_cam%lwcf_at_pc_top%f(1:ncell)
!    diag_physics_vars%swcf%f(1:ncell)  = diag_physics_vars%swcf%f(1:ncell)+pstate_cam%swcf_at_pc_top%f(1:ncell)
!
!    diag_physics_vars%asdir%f(1:ncell)  = diag_physics_vars%asdir%f(1:ncell)+pstate%atm_in_asdir_at_pc_surface%f(1:ncell)
!    diag_physics_vars%asdif%f(1:ncell)  = diag_physics_vars%asdif%f(1:ncell)+pstate%atm_in_asdif_at_pc_surface%f(1:ncell)
!    diag_physics_vars%aldir%f(1:ncell)  = diag_physics_vars%aldir%f(1:ncell)+pstate%atm_in_aldir_at_pc_surface%f(1:ncell)
!    diag_physics_vars%aldif%f(1:ncell)  = diag_physics_vars%aldif%f(1:ncell)+pstate%atm_in_aldif_at_pc_surface%f(1:ncell)
#endif
    
    end if

    diag_physics_vars_h1%ncount  = diag_physics_vars_h1%ncount + 1
    return
  end subroutine gcm_h1_accu_physics_variables

!
! dump depending on write_history frequency
!

  subroutine gcm_h1_dump_physics_variables
! local
    if(write_history_h1)then
!    diag_physics_vars_h1%precl%f = diag_physics_vars_h1%precl%f/diag_physics_vars_h1%ncount
!    diag_physics_vars_h1%precc%f = diag_physics_vars_h1%precc%f/diag_physics_vars_h1%ncount
    diag_physics_vars_h1%prect%f = diag_physics_vars_h1%prect%f/diag_physics_vars_h1%ncount

!    diag_physics_vars%ps%f    = diag_physics_vars%ps%f/diag_physics_vars%ncount
!    diag_physics_vars%snowh%f = diag_physics_vars%snowh%f/diag_physics_vars%ncount
!    diag_physics_vars%ts%f    = diag_physics_vars%ts%f/diag_physics_vars%ncount
!    diag_physics_vars%shflx%f = diag_physics_vars%shflx%f/diag_physics_vars%ncount
!    diag_physics_vars%qflx%f  = diag_physics_vars%qflx%f/diag_physics_vars%ncount
!    diag_physics_vars%phi%f   = diag_physics_vars%phi%f/diag_physics_vars%ncount
!    diag_physics_vars%uwind%f = diag_physics_vars%uwind%f/diag_physics_vars%ncount
!    diag_physics_vars%vwind%f = diag_physics_vars%vwind%f/diag_physics_vars%ncount
!    diag_physics_vars%omega%f = diag_physics_vars%omega%f/diag_physics_vars%ncount
!    diag_physics_vars%temp%f  = diag_physics_vars%temp%f/diag_physics_vars%ncount
!    diag_physics_vars%qv%f    = diag_physics_vars%qv%f/diag_physics_vars%ncount
!    diag_physics_vars%rh%f    = diag_physics_vars%rh%f/diag_physics_vars%ncount
!    diag_physics_vars%cld%f   = diag_physics_vars%cld%f/diag_physics_vars%ncount
!    diag_physics_vars%cldtot%f= diag_physics_vars%cldtot%f/diag_physics_vars%ncount
!    diag_physics_vars%cldlow%f= diag_physics_vars%cldlow%f/diag_physics_vars%ncount
!    diag_physics_vars%cldmed%f= diag_physics_vars%cldmed%f/diag_physics_vars%ncount
!    diag_physics_vars%cldhgh%f= diag_physics_vars%cldhgh%f/diag_physics_vars%ncount
!    diag_physics_vars%z500%f  = diag_physics_vars%z500%f/diag_physics_vars%ncount
!    diag_physics_vars_h1%u850%f  = diag_physics_vars_h1%u850%f/diag_physics_vars_h1%ncount
!    diag_physics_vars_h1%u200%f  = diag_physics_vars_h1%u200%f/diag_physics_vars_h1%ncount
!    diag_physics_vars_h1%v850%f  = diag_physics_vars_h1%v850%f/diag_physics_vars_h1%ncount

    diag_physics_vars_h1%flwut%f = diag_physics_vars_h1%flwut%f/diag_physics_vars_h1%ncount
!    diag_physics_vars%flwdt%f = diag_physics_vars%flwdt%f/diag_physics_vars%ncount
!    diag_physics_vars%flwus%f = diag_physics_vars%flwus%f/diag_physics_vars%ncount
!    diag_physics_vars%flwds%f = diag_physics_vars%flwds%f/diag_physics_vars%ncount
!    diag_physics_vars%flwutc%f= diag_physics_vars%flwutc%f/diag_physics_vars%ncount
!    diag_physics_vars%flwdtc%f= diag_physics_vars%flwdtc%f/diag_physics_vars%ncount
!    diag_physics_vars%flwusc%f= diag_physics_vars%flwusc%f/diag_physics_vars%ncount
!    diag_physics_vars%flwdsc%f= diag_physics_vars%flwdsc%f/diag_physics_vars%ncount

!    diag_physics_vars%fswut%f = diag_physics_vars%fswut%f/diag_physics_vars%ncount
!    diag_physics_vars%fswdt%f = diag_physics_vars%fswdt%f/diag_physics_vars%ncount
!    diag_physics_vars%fswus%f = diag_physics_vars%fswus%f/diag_physics_vars%ncount
!    diag_physics_vars%fswds%f = diag_physics_vars%fswds%f/diag_physics_vars%ncount
!    diag_physics_vars%fswutc%f= diag_physics_vars%fswutc%f/diag_physics_vars%ncount
!    diag_physics_vars%fswdtc%f= diag_physics_vars%fswdtc%f/diag_physics_vars%ncount
!    diag_physics_vars%fswusc%f= diag_physics_vars%fswusc%f/diag_physics_vars%ncount
!    diag_physics_vars%fswdsc%f= diag_physics_vars%fswdsc%f/diag_physics_vars%ncount
!    diag_physics_vars%lwcf%f  = diag_physics_vars%lwcf%f/diag_physics_vars%ncount
!    diag_physics_vars%swcf%f  = diag_physics_vars%swcf%f/diag_physics_vars%ncount

!    diag_physics_vars%asdir%f  = diag_physics_vars%asdir%f/diag_physics_vars%ncount
!    diag_physics_vars%asdif%f  = diag_physics_vars%asdif%f/diag_physics_vars%ncount
!    diag_physics_vars%aldir%f  = diag_physics_vars%aldir%f/diag_physics_vars%ncount
!    diag_physics_vars%aldif%f  = diag_physics_vars%aldif%f/diag_physics_vars%ncount
    end if

    return
  end subroutine gcm_h1_dump_physics_variables

  subroutine gcm_h1_rest_physics_variables
! reset to zero
    diag_physics_vars_h1%ncount  = 0

    if(write_history_h1)then
!    diag_physics_vars_h1%precl%f = zero
!    diag_physics_vars_h1%precc%f = zero
    diag_physics_vars_h1%prect%f = zero

!    diag_physics_vars%ps%f    = zero 
!    diag_physics_vars%ts%f    = zero 
!    diag_physics_vars%snowh%f    = zero 
!    diag_physics_vars%shflx%f = zero 
!    diag_physics_vars%qflx%f  = zero 
!    diag_physics_vars%phi%f   = zero 
!    diag_physics_vars%uwind%f = zero 
!    diag_physics_vars%vwind%f = zero 
!    diag_physics_vars%omega%f = zero 
!    diag_physics_vars%temp%f  = zero 
!    diag_physics_vars%qv%f    = zero 
!    diag_physics_vars%rh%f    = zero 
!    diag_physics_vars%cld%f   = zero 
!    diag_physics_vars%z500%f  = zero 
!    diag_physics_vars_h1%u850%f  = zero 
!    diag_physics_vars_h1%u200%f  = zero 
!    diag_physics_vars_h1%v850%f  = zero 

    diag_physics_vars_h1%flwut%f = zero 
!    diag_physics_vars%flwdt%f = zero 
!    diag_physics_vars%flwus%f = zero 
!    diag_physics_vars%flwds%f = zero 
!    diag_physics_vars%flwutc%f= zero 
!    diag_physics_vars%flwdtc%f= zero 
!    diag_physics_vars%flwusc%f= zero 
!    diag_physics_vars%flwdsc%f= zero 
!    diag_physics_vars%fswut%f = zero 
!    diag_physics_vars%fswdt%f = zero 
!    diag_physics_vars%fswus%f = zero 
!    diag_physics_vars%fswds%f = zero 
!    diag_physics_vars%fswutc%f= zero 
!    diag_physics_vars%fswdtc%f= zero 
!    diag_physics_vars%fswusc%f= zero 
!    diag_physics_vars%fswdsc%f= zero 
!    diag_physics_vars%lwcf%f  = zero 
!    diag_physics_vars%swcf%f  = zero 
!
!    diag_physics_vars%asdir%f  = zero 
!    diag_physics_vars%asdif%f  = zero 
!    diag_physics_vars%aldir%f  = zero 
!    diag_physics_vars%aldif%f  = zero 
 
    end if

    return
  end subroutine gcm_h1_rest_physics_variables

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

  end module grist_gcm_diagnose_h1_module
