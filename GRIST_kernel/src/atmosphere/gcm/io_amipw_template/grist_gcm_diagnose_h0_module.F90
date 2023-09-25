
!----------------------------------------------------------------------------
! Created on 2017-5-17
! Author: Yi Zhang
! Version 1.0
! Description: GCM Diagnostic module
! Revision history: based on the dycore diagnose module
!                   currently used for accumulating vars like precip
!                   1. seperate diagnostic vars to h0 and h1 as well
!                   2. All h0 vars are accumulated, so no sepration for inst
!                   and accu
!----------------------------------------------------------------------------

 module grist_gcm_diagnose_h0_module

  use grist_constants,              only: i4, r8, zero, one 
  use grist_domain_types,           only: global_domain
  use grist_data_types,             only: scalar_1d_field, scalar_2d_field, scalar_3d_field
  use grist_physics_data_structure, only: pstate
  use grist_nml_module,             only: ntracer, nlev, nlevp, write_history_h0, mif_index, nmif
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,   only: pstate_cam
  use grist_physics_update,         only: old_time_level
#endif
#ifdef AMIPW_PHYSICS
  use grist_PhysW_data_structure,   only: pstate_wrf, psurf_wrf, ptend_wrf
  use wv_saturation,                only: aqsat_grist
#endif
  use grist_dycore_vars_module,     only: dycoreVarCellFace, dycoreVarCellFull
  use grist_tracer_transport_vars_module, only: tracerVarCellFull

  implicit none

   private

   public   :: gcm_h0_diagnose_init         , &
               gcm_h0_diagnose_final        , & 
               gcm_h0_accu_physics_variables, &
               gcm_h0_dump_physics_variables, &
               gcm_h0_rest_physics_variables, &
               diag_physics_vars_1d_h0      , &
               diag_physics_vars_2d_h0

   type gcm_diag_vars_2d
! 2d
        type(scalar_2d_field)  :: uwind
        type(scalar_2d_field)  :: vwind
        type(scalar_2d_field)  :: temp
        type(scalar_2d_field)  :: mpressureFace
        type(scalar_2d_field)  :: omegaFull
#ifdef DYCORE_NHD
        type(scalar_2d_field)  :: wwwFace
#endif
        type(scalar_2d_field)  :: qv
        type(scalar_2d_field)  :: qc ! liqud
        type(scalar_2d_field)  :: qr ! rain 
        type(scalar_2d_field)  :: qi ! ice, snow, gra
        type(scalar_2d_field)  :: qs ! ice, snow, gra
        type(scalar_2d_field)  :: qg ! ice, snow, gra
        type(scalar_2d_field)  :: phiFace
        type(scalar_2d_field)  :: cloud
        type(scalar_2d_field)  :: cldcu
        type(scalar_2d_field)  :: cldst
        type(scalar_2d_field)  :: relhum
        type(scalar_2d_field)  :: rad_thten
        integer(i4)            :: ncount
    end type gcm_diag_vars_2d

    type gcm_diag_vars_1d
! 1d
        type(scalar_1d_field)  :: precc
        type(scalar_1d_field)  :: precl
        type(scalar_1d_field)  :: snowl
        type(scalar_1d_field)  :: grapl
        type(scalar_1d_field)  :: prect

        type(scalar_1d_field)  :: ps
        type(scalar_1d_field)  :: ts
        type(scalar_1d_field)  :: shflx
        type(scalar_1d_field)  :: qflx

!       type(scalar_1d_field)  :: z500
!       type(scalar_1d_field)  :: u200
!       type(scalar_1d_field)  :: u850
!       type(scalar_1d_field)  :: v850

        ! verticall integrated
        type(scalar_1d_field)  :: cldtot
        type(scalar_1d_field)  :: cldlow
        type(scalar_1d_field)  :: cldmed
        type(scalar_1d_field)  :: cldhgh

        type(scalar_1d_field)  :: flwut
        type(scalar_1d_field)  :: flwdt
        type(scalar_1d_field)  :: flwus
        type(scalar_1d_field)  :: flwds
        type(scalar_1d_field)  :: flwutc
        type(scalar_1d_field)  :: flwdtc
        type(scalar_1d_field)  :: flwusc
        type(scalar_1d_field)  :: flwdsc
        type(scalar_1d_field)  :: lwcf

        type(scalar_1d_field)  :: fswut
        type(scalar_1d_field)  :: fswdt
        type(scalar_1d_field)  :: fswus
        type(scalar_1d_field)  :: fswds
        type(scalar_1d_field)  :: fswutc
        type(scalar_1d_field)  :: fswdtc
        type(scalar_1d_field)  :: fswusc
        type(scalar_1d_field)  :: fswdsc
        type(scalar_1d_field)  :: swcf
! SURFACE-LSM
        type(scalar_1d_field)  :: asdir
        type(scalar_1d_field)  :: asdif
        type(scalar_1d_field)  :: aldir
        type(scalar_1d_field)  :: aldif
        type(scalar_1d_field)  :: emiss
        type(scalar_1d_field)  :: snowhland
        type(scalar_1d_field)  :: qsfc

        integer(i4)            :: ncount
   end type gcm_diag_vars_1d

   type(gcm_diag_vars_1d)   ::  diag_physics_vars_1d_h0
   type(gcm_diag_vars_2d)   ::  diag_physics_vars_2d_h0

   integer(i4) :: ncell

  contains

  subroutine gcm_h0_diagnose_init(mesh) ! till nv_full
    type(global_domain),  intent(in)   :: mesh

    diag_physics_vars_1d_h0%ncount  = 0
    diag_physics_vars_2d_h0%ncount  = 0

    if(write_history_h0)then

    call wrap_allocate_data1d(mesh,      diag_physics_vars_1d_h0%precc)
    call wrap_allocate_data1d(mesh,      diag_physics_vars_1d_h0%precl)
    call wrap_allocate_data1d(mesh,      diag_physics_vars_1d_h0%prect)
    call wrap_allocate_data1d(mesh,      diag_physics_vars_1d_h0%ps)

    call wrap_allocate_data2d(mesh,nlevp,diag_physics_vars_2d_h0%phiFace)
    call wrap_allocate_data2d(mesh,nlev ,diag_physics_vars_2d_h0%uwind)
    call wrap_allocate_data2d(mesh,nlev ,diag_physics_vars_2d_h0%vwind)
    call wrap_allocate_data2d(mesh,nlev ,diag_physics_vars_2d_h0%temp)
    call wrap_allocate_data2d(mesh,nlevp,diag_physics_vars_2d_h0%mpressureFace)
    call wrap_allocate_data2d(mesh,nlev ,diag_physics_vars_2d_h0%omegaFull)
#ifdef DYCORE_NHD
    call wrap_allocate_data2d(mesh,nlevp,diag_physics_vars_2d_h0%wwwFace)
#endif
    call wrap_allocate_data2d(mesh,nlev ,diag_physics_vars_2d_h0%qv)

#if (defined AMIPC_PHYSICS) || (defined AMIPW_PHYSICS)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%ts)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%shflx)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%qflx)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%cldtot)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%cldlow)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%cldmed)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%cldhgh)
#endif

#ifdef AMIPC_PHYSICS
    !call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%z500)
    !call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%u850)
    !call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%u200)
    !call wrap_allocate_data1d(mesh,     diag_physics_vars_h0%v850)
#endif

#if (defined AMIPC_PHYSICS) || (defined AMIPW_PHYSICS)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwut)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwdt)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwus)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwds)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwutc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwdtc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwusc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%flwdsc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%lwcf)

    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswut)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswdt)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswus)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswds)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswutc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswdtc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswusc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%fswdsc)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%swcf)

    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%cloud)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%cldcu)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%cldst)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%relhum)
#endif

#ifdef AMIPW_PHYSICS
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%snowl)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%grapl)

    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%asdir)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%asdif)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%aldir)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%aldif)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%emiss)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%snowhland)
    call wrap_allocate_data1d(mesh,     diag_physics_vars_1d_h0%qsfc)

    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%rad_thten)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%qc)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%qr)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%qi)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%qs)
    call wrap_allocate_data2d(mesh,nlev,diag_physics_vars_2d_h0%qg)
#endif

    end if

    ncell = mesh%nv_halo(1)

    return
  end subroutine gcm_h0_diagnose_init

  subroutine gcm_h0_diagnose_final
    if(allocated(diag_physics_vars_1d_h0%precc%f))         deallocate(diag_physics_vars_1d_h0%precc%f)
    if(allocated(diag_physics_vars_1d_h0%precl%f))         deallocate(diag_physics_vars_1d_h0%precl%f)
    if(allocated(diag_physics_vars_1d_h0%prect%f))         deallocate(diag_physics_vars_1d_h0%prect%f)
    if(allocated(diag_physics_vars_1d_h0%shflx%f))         deallocate(diag_physics_vars_1d_h0%shflx%f)
    if(allocated(diag_physics_vars_1d_h0%qflx%f))          deallocate(diag_physics_vars_1d_h0%qflx%f)
    if(allocated(diag_physics_vars_1d_h0%ps%f))            deallocate(diag_physics_vars_1d_h0%ps%f)
    if(allocated(diag_physics_vars_1d_h0%ts%f))            deallocate(diag_physics_vars_1d_h0%ts%f)

    if(allocated(diag_physics_vars_2d_h0%uwind%f))         deallocate(diag_physics_vars_2d_h0%uwind%f)
    if(allocated(diag_physics_vars_2d_h0%vwind%f))         deallocate(diag_physics_vars_2d_h0%vwind%f)
    if(allocated(diag_physics_vars_2d_h0%phiFace%f))       deallocate(diag_physics_vars_2d_h0%phiFace%f)
    if(allocated(diag_physics_vars_2d_h0%temp%f))          deallocate(diag_physics_vars_2d_h0%temp%f)
    if(allocated(diag_physics_vars_2d_h0%mpressureFace%f)) deallocate(diag_physics_vars_2d_h0%mpressureFace%f)
    if(allocated(diag_physics_vars_2d_h0%omegaFull%f))     deallocate(diag_physics_vars_2d_h0%omegaFull%f)
#ifdef DYCORE_NHD
    if(allocated(diag_physics_vars_2d_h0%wwwFace%f))       deallocate(diag_physics_vars_2d_h0%wwwFace%f)
#endif
    if(allocated(diag_physics_vars_2d_h0%qv%f))            deallocate(diag_physics_vars_2d_h0%qv%f)

#ifdef AMIPC_PHYSICS
    !if(allocated(diag_physics_vars_h0%z500%f))   deallocate(diag_physics_vars_h0%z500%f)
    !if(allocated(diag_physics_vars_h0%u850%f))   deallocate(diag_physics_vars_h0%u850%f)
    !if(allocated(diag_physics_vars_h0%u200%f))   deallocate(diag_physics_vars_h0%u200%f)
    !if(allocated(diag_physics_vars_h0%v850%f))   deallocate(diag_physics_vars_h0%v850%f)
#endif

#if (defined AMIPC_PHYSICS) || (AMIPW_PHYSICS)

    if(allocated(diag_physics_vars_1d_h0%flwut%f))  deallocate(diag_physics_vars_1d_h0%flwut%f )
    if(allocated(diag_physics_vars_1d_h0%flwdt%f))  deallocate(diag_physics_vars_1d_h0%flwdt%f )
    if(allocated(diag_physics_vars_1d_h0%flwus%f))  deallocate(diag_physics_vars_1d_h0%flwus%f )
    if(allocated(diag_physics_vars_1d_h0%flwds%f))  deallocate(diag_physics_vars_1d_h0%flwds%f )
    if(allocated(diag_physics_vars_1d_h0%flwutc%f)) deallocate(diag_physics_vars_1d_h0%flwutc%f)
    if(allocated(diag_physics_vars_1d_h0%flwdtc%f)) deallocate(diag_physics_vars_1d_h0%flwdtc%f)
    if(allocated(diag_physics_vars_1d_h0%flwusc%f)) deallocate(diag_physics_vars_1d_h0%flwusc%f)
    if(allocated(diag_physics_vars_1d_h0%flwdsc%f)) deallocate(diag_physics_vars_1d_h0%flwdsc%f)
    if(allocated(diag_physics_vars_1d_h0%lwcf%f))   deallocate(diag_physics_vars_1d_h0%lwcf%f)

    if(allocated(diag_physics_vars_1d_h0%fswut%f))  deallocate(diag_physics_vars_1d_h0%fswut%f )
    if(allocated(diag_physics_vars_1d_h0%fswdt%f))  deallocate(diag_physics_vars_1d_h0%fswdt%f )
    if(allocated(diag_physics_vars_1d_h0%fswus%f))  deallocate(diag_physics_vars_1d_h0%fswus%f )
    if(allocated(diag_physics_vars_1d_h0%fswds%f))  deallocate(diag_physics_vars_1d_h0%fswds%f )
    if(allocated(diag_physics_vars_1d_h0%fswutc%f)) deallocate(diag_physics_vars_1d_h0%fswutc%f)
    if(allocated(diag_physics_vars_1d_h0%fswdtc%f)) deallocate(diag_physics_vars_1d_h0%fswdtc%f)
    if(allocated(diag_physics_vars_1d_h0%fswusc%f)) deallocate(diag_physics_vars_1d_h0%fswusc%f)
    if(allocated(diag_physics_vars_1d_h0%fswdsc%f)) deallocate(diag_physics_vars_1d_h0%fswdsc%f)
    if(allocated(diag_physics_vars_1d_h0%swcf%f))   deallocate(diag_physics_vars_1d_h0%swcf%f)

    if(allocated(diag_physics_vars_2d_h0%cloud%f))  deallocate(diag_physics_vars_2d_h0%cloud%f)
    if(allocated(diag_physics_vars_2d_h0%cldcu%f))  deallocate(diag_physics_vars_2d_h0%cldcu%f)
    if(allocated(diag_physics_vars_2d_h0%cldst%f))  deallocate(diag_physics_vars_2d_h0%cldst%f)
    if(allocated(diag_physics_vars_2d_h0%relhum%f)) deallocate(diag_physics_vars_2d_h0%relhum%f)
#endif
#ifdef AMIPW_PHYSICS
    if(allocated(diag_physics_vars_1d_h0%cldtot%f)) deallocate(diag_physics_vars_1d_h0%cldtot%f)
    if(allocated(diag_physics_vars_1d_h0%cldlow%f)) deallocate(diag_physics_vars_1d_h0%cldlow%f)
    if(allocated(diag_physics_vars_1d_h0%cldmed%f)) deallocate(diag_physics_vars_1d_h0%cldmed%f)
    if(allocated(diag_physics_vars_1d_h0%cldhgh%f)) deallocate(diag_physics_vars_1d_h0%cldhgh%f)
    if(allocated(diag_physics_vars_1d_h0%snowl%f))  deallocate(diag_physics_vars_1d_h0%snowl%f)
    if(allocated(diag_physics_vars_1d_h0%grapl%f))  deallocate(diag_physics_vars_1d_h0%grapl%f)

    if(allocated(diag_physics_vars_1d_h0%asdir%f))  deallocate(diag_physics_vars_1d_h0%asdir%f)
    if(allocated(diag_physics_vars_1d_h0%asdif%f))  deallocate(diag_physics_vars_1d_h0%asdif%f)
    if(allocated(diag_physics_vars_1d_h0%aldir%f))  deallocate(diag_physics_vars_1d_h0%aldir%f)
    if(allocated(diag_physics_vars_1d_h0%aldif%f))  deallocate(diag_physics_vars_1d_h0%aldif%f)
    if(allocated(diag_physics_vars_1d_h0%emiss%f))  deallocate(diag_physics_vars_1d_h0%emiss%f)
    if(allocated(diag_physics_vars_1d_h0%snowhland%f)) deallocate(diag_physics_vars_1d_h0%snowhland%f)
    if(allocated(diag_physics_vars_1d_h0%qsfc%f))   deallocate(diag_physics_vars_1d_h0%qsfc%f)

    if(allocated(diag_physics_vars_2d_h0%rad_thten%f))  deallocate(diag_physics_vars_2d_h0%rad_thten%f)
    if(allocated(diag_physics_vars_2d_h0%qc%f))         deallocate(diag_physics_vars_2d_h0%qc%f)
    if(allocated(diag_physics_vars_2d_h0%qr%f))         deallocate(diag_physics_vars_2d_h0%qr%f)
    if(allocated(diag_physics_vars_2d_h0%qi%f))         deallocate(diag_physics_vars_2d_h0%qi%f)
    if(allocated(diag_physics_vars_2d_h0%qs%f))         deallocate(diag_physics_vars_2d_h0%qs%f)
    if(allocated(diag_physics_vars_2d_h0%qg%f))         deallocate(diag_physics_vars_2d_h0%qg%f)
#endif

    return
  end subroutine gcm_h0_diagnose_final
!
! accumulate at each model step
!
  subroutine gcm_h0_accu_physics_variables
   integer(i4) :: icell, ilev
   real(r8) :: sat_specific_humidity(nlev,ncell), esvp(nlev,ncell)
! local

    if(write_history_h0)then

    diag_physics_vars_1d_h0%precl%f = diag_physics_vars_1d_h0%precl%f+pstate%scalar_precl_surface%f
    diag_physics_vars_1d_h0%precc%f = diag_physics_vars_1d_h0%precc%f+pstate%scalar_precc_surface%f
    diag_physics_vars_1d_h0%prect%f = diag_physics_vars_1d_h0%prect%f+pstate%scalar_prect_surface%f

    diag_physics_vars_1d_h0%ps%f             = diag_physics_vars_1d_h0%ps%f            +dycoreVarCellFace%scalar_pressure_n%f(nlevp,:)
    diag_physics_vars_2d_h0%uwind%f          = diag_physics_vars_2d_h0%uwind%f         +dycoreVarCellFull%scalar_U_wind_n%f
    diag_physics_vars_2d_h0%vwind%f          = diag_physics_vars_2d_h0%vwind%f         +dycoreVarCellFull%scalar_V_wind_n%f
    diag_physics_vars_2d_h0%temp%f           = diag_physics_vars_2d_h0%temp%f          +dycoreVarCellFull%scalar_temp_n%f
    diag_physics_vars_2d_h0%mpressureFace%f  = diag_physics_vars_2d_h0%mpressureFace%f +dycoreVarCellFace%scalar_mpressure_n%f
    diag_physics_vars_2d_h0%omegaFull%f      = diag_physics_vars_2d_h0%omegaFull%f     +dycoreVarCellFull%scalar_omega_n%f 
#ifdef DYCORE_NHD
    diag_physics_vars_2d_h0%wwwFace%f        = diag_physics_vars_2d_h0%wwwFace%f       +dycoreVarCellFace%scalar_www_n%f
#endif
    diag_physics_vars_2d_h0%qv%f             = diag_physics_vars_2d_h0%qv%f            +tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:)
    diag_physics_vars_2d_h0%phiFace%f        = diag_physics_vars_2d_h0%phiFace%f       +dycoreVarCellFace%scalar_geopotential_n%f

#if (defined AMIPC_PHYSICS) || (defined AMIPW_PHYSICS)
    diag_physics_vars_1d_h0%ts%f    = diag_physics_vars_1d_h0%ts%f   +pstate%ts_at_pc_surface%f
    diag_physics_vars_1d_h0%shflx%f = diag_physics_vars_1d_h0%shflx%f+pstate%atm_in_shflx_at_pc_surface%f
    diag_physics_vars_1d_h0%qflx%f  = diag_physics_vars_1d_h0%qflx%f +pstate%atm_in_qflx_at_pc_surface%f(1,:)
#endif

#ifdef AMIPC_PHYSICS
    diag_physics_vars_2d_h0%cloud%f(:,1:ncell) = diag_physics_vars_2d_h0%cloud%f(:,1:ncell)+pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncell)
    !diag_physics_vars_h0%z500%f  = diag_physics_vars_h0%z500%f+pstate_cam%diag_z_at_500hpa%f
    !diag_physics_vars_h0%u850%f  = diag_physics_vars_h0%u850%f+pstate_cam%diag_u_at_850hpa%f
    !diag_physics_vars_h0%u200%f  = diag_physics_vars_h0%u200%f+pstate_cam%diag_u_at_200hpa%f
    !diag_physics_vars_h0%v850%f  = diag_physics_vars_h0%v850%f+pstate_cam%diag_v_at_850hpa%f
 
    diag_physics_vars_1d_h0%flwut%f(1:ncell) = diag_physics_vars_1d_h0%flwut%f(1:ncell) +pstate_cam%flwut_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%flwdt%f(1:ncell) = diag_physics_vars_1d_h0%flwdt%f(1:ncell) +pstate_cam%flwdt_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%flwus%f(1:ncell) = diag_physics_vars_1d_h0%flwus%f(1:ncell) +pstate_cam%flwus_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%flwds%f(1:ncell) = diag_physics_vars_1d_h0%flwds%f(1:ncell) +pstate_cam%flwds_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%flwutc%f(1:ncell)= diag_physics_vars_1d_h0%flwutc%f(1:ncell)+pstate_cam%flwutc_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%flwdtc%f(1:ncell)= diag_physics_vars_1d_h0%flwdtc%f(1:ncell)+pstate_cam%flwdtc_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%flwusc%f(1:ncell)= diag_physics_vars_1d_h0%flwusc%f(1:ncell)+pstate_cam%flwusc_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%flwdsc%f(1:ncell)= diag_physics_vars_1d_h0%flwdsc%f(1:ncell)+pstate_cam%flwdsc_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%lwcf%f  (1:ncell)= diag_physics_vars_1d_h0%lwcf%f(1:ncell)  +pstate_cam%lwcf_at_pc_top%f(1:ncell)

    diag_physics_vars_1d_h0%fswut%f(1:ncell) = diag_physics_vars_1d_h0%fswut%f(1:ncell) +pstate_cam%fswut_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%fswdt%f(1:ncell) = diag_physics_vars_1d_h0%fswdt%f(1:ncell) +pstate_cam%fswdt_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%fswus%f(1:ncell) = diag_physics_vars_1d_h0%fswus%f(1:ncell) +pstate_cam%fswus_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%fswds%f(1:ncell) = diag_physics_vars_1d_h0%fswds%f(1:ncell) +pstate_cam%fswds_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%fswutc%f(1:ncell)= diag_physics_vars_1d_h0%fswutc%f(1:ncell)+pstate_cam%fswutc_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%fswdtc%f(1:ncell)= diag_physics_vars_1d_h0%fswdtc%f(1:ncell)+pstate_cam%fswdtc_at_pc_top%f(1:ncell)
    diag_physics_vars_1d_h0%fswusc%f(1:ncell)= diag_physics_vars_1d_h0%fswusc%f(1:ncell)+pstate_cam%fswusc_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%fswdsc%f(1:ncell)= diag_physics_vars_1d_h0%fswdsc%f(1:ncell)+pstate_cam%fswdsc_at_pc_surface%f(1:ncell)
    diag_physics_vars_1d_h0%swcf%f  (1:ncell)= diag_physics_vars_1d_h0%swcf%f(1:ncell)  +pstate_cam%swcf_at_pc_top%f(1:ncell)
#endif

#ifdef AMIPW_PHYSICS
    diag_physics_vars_1d_h0%cldtot%f(1:ncell) = diag_physics_vars_1d_h0%cldtot%f(1:ncell) + pstate_wrf%cldtot(1:ncell,1)
    diag_physics_vars_1d_h0%cldlow%f(1:ncell) = diag_physics_vars_1d_h0%cldlow%f(1:ncell) + pstate_wrf%cldlow(1:ncell,1)
    diag_physics_vars_1d_h0%cldmed%f(1:ncell) = diag_physics_vars_1d_h0%cldmed%f(1:ncell) + pstate_wrf%cldmed(1:ncell,1)
    diag_physics_vars_1d_h0%cldhgh%f(1:ncell) = diag_physics_vars_1d_h0%cldhgh%f(1:ncell) + pstate_wrf%cldhgh(1:ncell,1)
    diag_physics_vars_1d_h0%snowl%f(1:ncell)  = diag_physics_vars_1d_h0%snowl%f (1:ncell) + pstate%scalar_snowl_surface%f(1:ncell)
    diag_physics_vars_1d_h0%grapl%f(1:ncell)  = diag_physics_vars_1d_h0%grapl%f (1:ncell) + pstate%scalar_grapl_surface%f(1:ncell)

    diag_physics_vars_1d_h0%flwut%f( 1:ncell) = diag_physics_vars_1d_h0%flwut%f (1:ncell) + pstate_wrf%lwupt(1:ncell,1)
    diag_physics_vars_1d_h0%flwdt%f( 1:ncell) = diag_physics_vars_1d_h0%flwdt%f (1:ncell) + pstate_wrf%lwdnt(1:ncell,1)
    diag_physics_vars_1d_h0%flwus%f( 1:ncell) = diag_physics_vars_1d_h0%flwus%f (1:ncell) + pstate_wrf%lwupb(1:ncell,1)
    diag_physics_vars_1d_h0%flwds%f( 1:ncell) = diag_physics_vars_1d_h0%flwds%f (1:ncell) + pstate_wrf%lwdnb(1:ncell,1)
    diag_physics_vars_1d_h0%flwutc%f(1:ncell) = diag_physics_vars_1d_h0%flwutc%f(1:ncell) + pstate_wrf%lwuptc(1:ncell,1)
    diag_physics_vars_1d_h0%flwdtc%f(1:ncell) = diag_physics_vars_1d_h0%flwdtc%f(1:ncell) + pstate_wrf%lwdntc(1:ncell,1)
    diag_physics_vars_1d_h0%flwusc%f(1:ncell) = diag_physics_vars_1d_h0%flwusc%f(1:ncell) + pstate_wrf%lwupbc(1:ncell,1)
    diag_physics_vars_1d_h0%flwdsc%f(1:ncell) = diag_physics_vars_1d_h0%flwdsc%f(1:ncell) + pstate_wrf%lwdnbc(1:ncell,1)
    diag_physics_vars_1d_h0%lwcf%f  (1:ncell) = diag_physics_vars_1d_h0%lwcf%f(1:ncell)   + pstate_wrf%lwcf(1:ncell,1)

    diag_physics_vars_1d_h0%fswut%f( 1:ncell) = diag_physics_vars_1d_h0%fswut%f (1:ncell) + pstate_wrf%swupt(1:ncell,1)
    diag_physics_vars_1d_h0%fswdt%f( 1:ncell) = diag_physics_vars_1d_h0%fswdt%f (1:ncell) + pstate_wrf%swdnt(1:ncell,1)
    diag_physics_vars_1d_h0%fswus%f( 1:ncell) = diag_physics_vars_1d_h0%fswus%f (1:ncell) + pstate_wrf%swupb(1:ncell,1)
    diag_physics_vars_1d_h0%fswds%f( 1:ncell) = diag_physics_vars_1d_h0%fswds%f (1:ncell) + pstate_wrf%swdnb(1:ncell,1)
    diag_physics_vars_1d_h0%fswutc%f(1:ncell) = diag_physics_vars_1d_h0%fswutc%f(1:ncell) + pstate_wrf%swuptc(1:ncell,1)
    diag_physics_vars_1d_h0%fswdtc%f(1:ncell) = diag_physics_vars_1d_h0%fswdtc%f(1:ncell) + pstate_wrf%swdntc(1:ncell,1)
    diag_physics_vars_1d_h0%fswusc%f(1:ncell) = diag_physics_vars_1d_h0%fswusc%f(1:ncell) + pstate_wrf%swupbc(1:ncell,1)
    diag_physics_vars_1d_h0%fswdsc%f(1:ncell) = diag_physics_vars_1d_h0%fswdsc%f(1:ncell) + pstate_wrf%swdnbc(1:ncell,1)
    diag_physics_vars_1d_h0%swcf%f  (1:ncell) = diag_physics_vars_1d_h0%swcf%f(1:ncell)   + pstate_wrf%swcf(1:ncell,1)

    diag_physics_vars_1d_h0%asdir%f(1:ncell)  = diag_physics_vars_1d_h0%asdir%f(1:ncell)  + psurf_wrf%asdir(1:ncell,1)
    diag_physics_vars_1d_h0%asdif%f(1:ncell)  = diag_physics_vars_1d_h0%asdif%f(1:ncell)  + psurf_wrf%asdif(1:ncell,1)
    diag_physics_vars_1d_h0%aldir%f(1:ncell)  = diag_physics_vars_1d_h0%aldir%f(1:ncell)  + psurf_wrf%aldir(1:ncell,1)
    diag_physics_vars_1d_h0%aldif%f(1:ncell)  = diag_physics_vars_1d_h0%aldif%f(1:ncell)  + psurf_wrf%aldif(1:ncell,1)
    diag_physics_vars_1d_h0%emiss%f(1:ncell)  = diag_physics_vars_1d_h0%emiss%f(1:ncell)  + psurf_wrf%emiss(1:ncell,1)
    diag_physics_vars_1d_h0%snowhland%f(1:ncell) = diag_physics_vars_1d_h0%snowhland%f(1:ncell) + psurf_wrf%snow(1:ncell,1)
    diag_physics_vars_1d_h0%qsfc%f(1:ncell)   = diag_physics_vars_1d_h0%qsfc%f(1:ncell)   + psurf_wrf%qsfc (1:ncell,1)

!================================================================================================
! obtain RELHUM. We use state from dynamics, not physics, because some physics
! may update 
! t but not p, so a little bit inconsistency (one-step impact is trivial). Use
! dynamical state
! is comparible with other vars, e.g., geop, omega 
!=================================================================================================

    call aqsat_grist(t = dycoreVarCellFull%scalar_temp_n%f     (1:nlev,1:ncell), &
                     p = dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,1:ncell), &
                     es= esvp(1:nlev,1:ncell)              , &
                     qs= sat_specific_humidity(1:nlev,1:ncell), &
                     ii= ncell, ilen=ncell, kk=nlev, kstart=1, kend=nlev)

    do icell = 1, ncell
       do ilev = 1, nlev
! not saved in rst, and because cldfra is called at rad step, this diag will be incorrect for the first restart-run output variable
          diag_physics_vars_2d_h0%cloud%f(ilev,icell)     = diag_physics_vars_2d_h0%cloud%f(ilev,icell)+pstate_wrf%cldfra(icell,nlev+1-ilev,1) 
          diag_physics_vars_2d_h0%cldcu%f(ilev,icell)     = diag_physics_vars_2d_h0%cldcu%f(ilev,icell)+pstate_wrf%cldcum(icell,nlev+1-ilev,1)
          diag_physics_vars_2d_h0%cldst%f(ilev,icell)     = diag_physics_vars_2d_h0%cldst%f(ilev,icell)+pstate_wrf%cldstr(icell,nlev+1-ilev,1)
          diag_physics_vars_2d_h0%rad_thten%f(ilev,icell) = diag_physics_vars_2d_h0%rad_thten%f(ilev,icell)+ptend_wrf%rthraten(icell,nlev+1-ilev,1)
          diag_physics_vars_2d_h0%relhum%f(ilev,icell)    = diag_physics_vars_2d_h0%relhum%f(ilev,icell)+&
                                                 tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,icell)&
                                               /(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,icell)))/sat_specific_humidity(ilev,icell)*100._r8
       end do
    end do
! hard-coded here
    diag_physics_vars_2d_h0%qc%f    = diag_physics_vars_2d_h0%qc%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(2,:,:) ! 6 is default, vcrisg
    diag_physics_vars_2d_h0%qr%f    = diag_physics_vars_2d_h0%qr%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(3,:,:) ! 6 is default, vcrisg
    diag_physics_vars_2d_h0%qi%f    = diag_physics_vars_2d_h0%qi%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(4,:,:)
    diag_physics_vars_2d_h0%qs%f    = diag_physics_vars_2d_h0%qs%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(5,:,:)
    diag_physics_vars_2d_h0%qg%f    = diag_physics_vars_2d_h0%qg%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(6,:,:)
#endif

    end if

    diag_physics_vars_2d_h0%ncount  = diag_physics_vars_2d_h0%ncount + 1
    diag_physics_vars_1d_h0%ncount  = diag_physics_vars_1d_h0%ncount + 1

    return
  end subroutine gcm_h0_accu_physics_variables

!
! dump depending on write_history frequency
!

  subroutine gcm_h0_dump_physics_variables
! local

    if(write_history_h0)then

    diag_physics_vars_1d_h0%precl%f = diag_physics_vars_1d_h0%precl%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%precc%f = diag_physics_vars_1d_h0%precc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%prect%f = diag_physics_vars_1d_h0%prect%f/diag_physics_vars_1d_h0%ncount

    diag_physics_vars_1d_h0%ps%f             = diag_physics_vars_1d_h0%ps%f     /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_2d_h0%uwind%f          = diag_physics_vars_2d_h0%uwind%f  /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%vwind%f          = diag_physics_vars_2d_h0%vwind%f  /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%phiFace%f        = diag_physics_vars_2d_h0%phiFace%f/diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%temp%f           = diag_physics_vars_2d_h0%temp%f   /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%mpressureFace%f  = diag_physics_vars_2d_h0%mpressureFace%f/diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%omegaFull%f      = diag_physics_vars_2d_h0%omegaFull%f    /diag_physics_vars_2d_h0%ncount
#ifdef DYCORE_NHD
    diag_physics_vars_2d_h0%wwwFace%f        = diag_physics_vars_2d_h0%wwwFace%f      /diag_physics_vars_2d_h0%ncount
#endif
    diag_physics_vars_2d_h0%qv%f             = diag_physics_vars_2d_h0%qv%f           /diag_physics_vars_2d_h0%ncount

#if (defined AMIPC_PHYSICS) || (defined AMIPW_PHYSICS)
    diag_physics_vars_1d_h0%ts%f    = diag_physics_vars_1d_h0%ts%f   /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%shflx%f = diag_physics_vars_1d_h0%shflx%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%qflx%f  = diag_physics_vars_1d_h0%qflx%f /diag_physics_vars_1d_h0%ncount
#endif

#ifdef AMIPC_PHYSICS
    !diag_physics_vars_h0%z500%f  = diag_physics_vars_h0%z500%f/diag_physics_vars_h0%ncount
    !diag_physics_vars_h0%u850%f  = diag_physics_vars_h0%u850%f/diag_physics_vars_h0%ncount
    !diag_physics_vars_h0%u200%f  = diag_physics_vars_h0%u200%f/diag_physics_vars_h0%ncount
    !diag_physics_vars_h0%v850%f  = diag_physics_vars_h0%v850%f/diag_physics_vars_h0%ncount
#endif

#if (defined AMIPC_PHYSICS) || (defined AMIPW_PHYSICS)
    diag_physics_vars_1d_h0%cldtot%f= diag_physics_vars_1d_h0%cldtot%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%cldlow%f= diag_physics_vars_1d_h0%cldlow%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%cldmed%f= diag_physics_vars_1d_h0%cldmed%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%cldhgh%f= diag_physics_vars_1d_h0%cldhgh%f/diag_physics_vars_1d_h0%ncount

    diag_physics_vars_1d_h0%snowl%f = diag_physics_vars_1d_h0%snowl%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%grapl%f = diag_physics_vars_1d_h0%grapl%f /diag_physics_vars_1d_h0%ncount

    diag_physics_vars_1d_h0%flwut%f = diag_physics_vars_1d_h0%flwut%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%flwdt%f = diag_physics_vars_1d_h0%flwdt%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%flwus%f = diag_physics_vars_1d_h0%flwus%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%flwds%f = diag_physics_vars_1d_h0%flwds%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%flwutc%f= diag_physics_vars_1d_h0%flwutc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%flwdtc%f= diag_physics_vars_1d_h0%flwdtc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%flwusc%f= diag_physics_vars_1d_h0%flwusc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%flwdsc%f= diag_physics_vars_1d_h0%flwdsc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%lwcf%f  = diag_physics_vars_1d_h0%lwcf%f  /diag_physics_vars_1d_h0%ncount

    diag_physics_vars_1d_h0%fswut%f = diag_physics_vars_1d_h0%fswut%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%fswdt%f = diag_physics_vars_1d_h0%fswdt%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%fswus%f = diag_physics_vars_1d_h0%fswus%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%fswds%f = diag_physics_vars_1d_h0%fswds%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%fswutc%f= diag_physics_vars_1d_h0%fswutc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%fswdtc%f= diag_physics_vars_1d_h0%fswdtc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%fswusc%f= diag_physics_vars_1d_h0%fswusc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%fswdsc%f= diag_physics_vars_1d_h0%fswdsc%f/diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%swcf%f  = diag_physics_vars_1d_h0%swcf%f  /diag_physics_vars_1d_h0%ncount

    diag_physics_vars_2d_h0%cloud%f = diag_physics_vars_2d_h0%cloud%f /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%cldcu%f = diag_physics_vars_2d_h0%cldcu%f /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%cldst%f = diag_physics_vars_2d_h0%cldst%f /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%relhum%f= diag_physics_vars_2d_h0%relhum%f/diag_physics_vars_2d_h0%ncount
#endif

#ifdef AMIPW_PHYSICS
    diag_physics_vars_2d_h0%rad_thten%f = diag_physics_vars_2d_h0%rad_thten%f /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_1d_h0%asdir%f     = diag_physics_vars_1d_h0%asdir%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%asdif%f     = diag_physics_vars_1d_h0%asdif%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%aldir%f     = diag_physics_vars_1d_h0%aldir%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%aldif%f     = diag_physics_vars_1d_h0%aldif%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%emiss%f     = diag_physics_vars_1d_h0%emiss%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%snowhland%f = diag_physics_vars_1d_h0%snowhland%f /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_1d_h0%qsfc%f      = diag_physics_vars_1d_h0%qsfc%f  /diag_physics_vars_1d_h0%ncount
    diag_physics_vars_2d_h0%qc%f        = diag_physics_vars_2d_h0%qc%f    /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%qr%f        = diag_physics_vars_2d_h0%qr%f    /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%qi%f        = diag_physics_vars_2d_h0%qi%f    /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%qs%f        = diag_physics_vars_2d_h0%qs%f    /diag_physics_vars_2d_h0%ncount
    diag_physics_vars_2d_h0%qg%f        = diag_physics_vars_2d_h0%qg%f    /diag_physics_vars_2d_h0%ncount
#endif
    end if
 
    return
  end subroutine gcm_h0_dump_physics_variables

  subroutine gcm_h0_rest_physics_variables
! reset to zero
    diag_physics_vars_1d_h0%ncount  = 0
    diag_physics_vars_2d_h0%ncount  = 0

    if(write_history_h0)then

    diag_physics_vars_1d_h0%precl%f  = zero
    diag_physics_vars_1d_h0%precc%f  = zero
    diag_physics_vars_1d_h0%prect%f  = zero
    diag_physics_vars_1d_h0%ps%f     = zero

    diag_physics_vars_2d_h0%uwind%f     = zero
    diag_physics_vars_2d_h0%vwind%f     = zero
    diag_physics_vars_2d_h0%phiFace%f   = zero
    diag_physics_vars_2d_h0%temp%f      = zero
    diag_physics_vars_2d_h0%mpressureFace%f  = zero
    diag_physics_vars_2d_h0%omegaFull%f = zero
#ifdef DYCORE_NHD
    diag_physics_vars_2d_h0%wwwFace%f   = zero
#endif
    diag_physics_vars_2d_h0%qv%f        = zero

#if (defined AMIPC_PHYSICS) || (defined AMIPW_PHYSICS)
    diag_physics_vars_1d_h0%ts%f     = zero
    diag_physics_vars_1d_h0%shflx%f  = zero
    diag_physics_vars_1d_h0%qflx%f   = zero
#endif

#ifdef AMIPC_PHYSICS
    !diag_physics_vars_1d_h0%z500%f   = zero 
    !diag_physics_vars_1d_h0%u850%f   = zero 
    !diag_physics_vars_1d_h0%u200%f   = zero 
    !diag_physics_vars_1d_h0%v850%f   = zero 
#endif

#if (defined AMIPC_PHYSICS) || (defined AMIPW_PHYSICS)
    diag_physics_vars_1d_h0%cldtot%f = zero
    diag_physics_vars_1d_h0%cldlow%f = zero
    diag_physics_vars_1d_h0%cldmed%f = zero
    diag_physics_vars_1d_h0%cldhgh%f = zero

    diag_physics_vars_1d_h0%flwut%f = zero
    diag_physics_vars_1d_h0%flwdt%f = zero
    diag_physics_vars_1d_h0%flwus%f = zero
    diag_physics_vars_1d_h0%flwds%f = zero
    diag_physics_vars_1d_h0%flwutc%f= zero
    diag_physics_vars_1d_h0%flwdtc%f= zero
    diag_physics_vars_1d_h0%flwusc%f= zero
    diag_physics_vars_1d_h0%flwdsc%f= zero
    diag_physics_vars_1d_h0%lwcf%f  = zero

    diag_physics_vars_1d_h0%fswut%f = zero
    diag_physics_vars_1d_h0%fswdt%f = zero
    diag_physics_vars_1d_h0%fswus%f = zero
    diag_physics_vars_1d_h0%fswds%f = zero
    diag_physics_vars_1d_h0%fswutc%f= zero
    diag_physics_vars_1d_h0%fswdtc%f= zero
    diag_physics_vars_1d_h0%fswusc%f= zero
    diag_physics_vars_1d_h0%fswdsc%f= zero
    diag_physics_vars_1d_h0%swcf%f  = zero

    diag_physics_vars_2d_h0%cloud%f = zero
    diag_physics_vars_2d_h0%cldcu%f = zero
    diag_physics_vars_2d_h0%cldst%f = zero
    diag_physics_vars_2d_h0%relhum%f = zero 
#endif

#ifdef AMIPW_PHYSICS
    diag_physics_vars_1d_h0%snowl%f  = zero
    diag_physics_vars_1d_h0%grapl%f  = zero

    diag_physics_vars_1d_h0%asdir%f  = zero
    diag_physics_vars_1d_h0%asdif%f  = zero
    diag_physics_vars_1d_h0%aldir%f  = zero
    diag_physics_vars_1d_h0%aldif%f  = zero
    diag_physics_vars_1d_h0%emiss%f  = zero
    diag_physics_vars_1d_h0%snowhland%f = zero
    diag_physics_vars_1d_h0%qsfc%f   = zero

    diag_physics_vars_2d_h0%rad_thten%f = zero
    diag_physics_vars_2d_h0%qc%f     = zero
    diag_physics_vars_2d_h0%qr%f     = zero
    diag_physics_vars_2d_h0%qi%f     = zero
    diag_physics_vars_2d_h0%qs%f     = zero
    diag_physics_vars_2d_h0%qg%f     = zero
#endif

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
