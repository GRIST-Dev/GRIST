
!----------------------------------------------------------------------------
! Created on 2017-5-17
! Author: Yi Zhang
! Version 1.0
! Description: GCM Diagnostic module
! Revision history: based on the dycore diagnose module
!                   currently used for accumulating vars like precip
!                   1. seperate diagnostic vars to h0 and h1 as well
!                   2. seperate diag vars into accu and inst; accu needs to be
!                   accumulated->dump->rest; inst only filled when needed
!----------------------------------------------------------------------------

 module grist_gcm_diagnose_h1_module

  use grist_constants,              only: i4, r8, zero, gravity, one 
  use grist_hpe_constants,          only: eta_face_a, eta_face_b, eta_full_a, eta_full_b
  use grist_domain_types,           only: global_domain
  use grist_data_types,             only: scalar_1d_field, scalar_2d_field, scalar_3d_field
  use grist_dycore_vars_module,     only: dycoreVarCellFull, dycoreVarCellFace
  use grist_tracer_transport_vars_module, only: tracerVarCellFull
  use grist_physics_data_structure, only: pstate
  use grist_nml_module,             only: ntracer, nlev, nlevp, mif_index, nmif, write_history_h1
#ifdef AMIPC_PHYSICS
  use grist_PhysC_data_structure,   only: pstate_cam
  use grist_physics_update,         only: old_time_level
#endif
#ifdef AMIPW_PHYSICS
  use grist_PhysW_data_structure,   only: pstate_wrf, psurf_wrf, ptend_wrf
  use grist_vi_cldfrac_diagnostics, only: cloud_diagnostics_calc
  use wv_saturation,                only: aqsat_grist
  use grist_math_module,            only: lininterp
#endif

  implicit none

   private

   public   :: gcm_h1_diagnose_init,  &
               gcm_h1_diagnose_final, &
               gcm_h1_1d_inst_physics_variables, &
               gcm_h1_2d_inst_physics_variables, &
               gcm_h1_1d_accu_physics_variables, &
               gcm_h1_1d_dump_physics_variables, &
               gcm_h1_1d_rest_physics_variables, &
               gcm_h1_2d_accu_physics_variables, &
               gcm_h1_2d_dump_physics_variables, &
               gcm_h1_2d_rest_physics_variables, &
               diag_phys_accu_vars_h1_1d,        &
               diag_phys_accu_vars_h1_2d,        &
               diag_phys_inst_vars_h1_1d,        &
               diag_phys_inst_vars_h1_2d

   type gcm_diag_accu_vars_2d
#ifndef DYAMOND
        type(scalar_2d_field)  :: uwind
        type(scalar_2d_field)  :: vwind
        type(scalar_2d_field)  :: temp
        type(scalar_2d_field)  :: mpressureFace
        type(scalar_2d_field)  :: omegaFull
#ifdef DYCORE_NHD
        type(scalar_2d_field)  :: wwwFace
#endif
        type(scalar_2d_field)  :: phiFace
        type(scalar_2d_field)  :: qv ! as in wsm6 and first six of morr
        type(scalar_2d_field)  :: qc ! as in wsm6 and first six of morr
        type(scalar_2d_field)  :: qr ! as in wsm6 and first six of morr
        type(scalar_2d_field)  :: qi ! as in wsm6 and first six of morr
        type(scalar_2d_field)  :: qs ! as in wsm6 and first six of morr
        type(scalar_2d_field)  :: qg ! as in wsm6 and first six of morr
#endif
        type(scalar_2d_field)  :: rad_thten ! tendency
        integer(i4)            :: ncount
    end type gcm_diag_accu_vars_2d

    type gcm_diag_accu_vars_1d
!
! These vars are either accumulated (via _accu) or instant 
!
        type(scalar_1d_field)  :: precc ! rate
        type(scalar_1d_field)  :: precl ! rate
        type(scalar_1d_field)  :: snowl ! rate
        type(scalar_1d_field)  :: grapl ! rate
        type(scalar_1d_field)  :: prect ! rate
        type(scalar_1d_field)  :: ps
        type(scalar_1d_field)  :: ts
        type(scalar_1d_field)  :: shflx
        type(scalar_1d_field)  :: qflx

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
! surface-lsm
#ifndef DYAMOND
        type(scalar_1d_field)  :: asdir
        type(scalar_1d_field)  :: asdif
        type(scalar_1d_field)  :: aldir
        type(scalar_1d_field)  :: aldif
        type(scalar_1d_field)  :: emiss
#endif
        type(scalar_1d_field)  :: qsfc
        integer(i4)            :: ncount
   end type gcm_diag_accu_vars_1d

   type gcm_diag_inst_vars_1d
        type(scalar_1d_field)  :: vert_level, vertp_level
        type(scalar_1d_field)  :: hyam, hybm, hyai, hybi
! verticall integrated
        type(scalar_1d_field)  :: viqv
        type(scalar_1d_field)  :: viqc
        type(scalar_1d_field)  :: viqr
        type(scalar_1d_field)  :: viqi
        type(scalar_1d_field)  :: viqs
        type(scalar_1d_field)  :: viqg
        type(scalar_1d_field)  :: cldtot
        type(scalar_1d_field)  :: cldlow
        type(scalar_1d_field)  :: cldmed
        type(scalar_1d_field)  :: cldhgh
! skin
        type(scalar_1d_field)  :: u10m
        type(scalar_1d_field)  :: v10m
        type(scalar_1d_field)  :: th2m
        type(scalar_1d_field)  :: t2m
        type(scalar_1d_field)  :: q2m
        type(scalar_1d_field)  :: tskin
        type(scalar_1d_field)  :: snowh
        type(scalar_1d_field)  :: xice
        type(scalar_1d_field)  :: hfx
        type(scalar_1d_field)  :: qfx
        type(scalar_1d_field)  :: lh
        type(scalar_1d_field)  :: qsfc
        type(scalar_1d_field)  :: emiss
! verticall-interpolated vars
        type(scalar_1d_field)  :: relhum200
        type(scalar_1d_field)  :: relhum500
        type(scalar_1d_field)  :: relhum700
        type(scalar_1d_field)  :: relhum850
        type(scalar_1d_field)  :: omega200
        type(scalar_1d_field)  :: omega500
        type(scalar_1d_field)  :: omega700
        type(scalar_1d_field)  :: omega850
        type(scalar_1d_field)  :: zzz200
        type(scalar_1d_field)  :: zzz500
        type(scalar_1d_field)  :: zzz700
        type(scalar_1d_field)  :: zzz850
   end type gcm_diag_inst_vars_1d

   type gcm_diag_inst_vars_2d
        type(scalar_2d_field)  :: cloud
#ifndef DYAMOND
        type(scalar_2d_field)  :: cldcu
        type(scalar_2d_field)  :: cldst
#endif
        type(scalar_2d_field)  :: relhum
        type(scalar_2d_field)  :: refl_10cm
   end type gcm_diag_inst_vars_2d

   type(gcm_diag_accu_vars_1d)   ::  diag_phys_accu_vars_h1_1d
   type(gcm_diag_accu_vars_2d)   ::  diag_phys_accu_vars_h1_2d
   type(gcm_diag_inst_vars_1d)   ::  diag_phys_inst_vars_h1_1d
   type(gcm_diag_inst_vars_2d)   ::  diag_phys_inst_vars_h1_2d

   integer(i4) :: ncell

  contains

  subroutine gcm_h1_diagnose_init(mesh) ! till nv_full
    type(global_domain),  intent(in)   :: mesh
    integer(i4) :: ilev

    diag_phys_accu_vars_h1_1d%ncount  = 0
    diag_phys_accu_vars_h1_2d%ncount  = 0

    if(write_history_h1)then
       if(.not.allocated(diag_phys_inst_vars_h1_1d%vert_level%f))  allocate(diag_phys_inst_vars_h1_1d%vert_level%f(nlev))
       if(.not.allocated(diag_phys_inst_vars_h1_1d%vertp_level%f)) allocate(diag_phys_inst_vars_h1_1d%vertp_level%f(nlevp))
       if(.not.allocated(diag_phys_inst_vars_h1_1d%hyam%f))        allocate(diag_phys_inst_vars_h1_1d%hyam%f(nlev))
       if(.not.allocated(diag_phys_inst_vars_h1_1d%hybm%f))        allocate(diag_phys_inst_vars_h1_1d%hybm%f(nlev))
       if(.not.allocated(diag_phys_inst_vars_h1_1d%hyai%f))        allocate(diag_phys_inst_vars_h1_1d%hyai%f(nlevp))
       if(.not.allocated(diag_phys_inst_vars_h1_1d%hybi%f))        allocate(diag_phys_inst_vars_h1_1d%hybi%f(nlevp))

       diag_phys_inst_vars_h1_1d%vert_level%pos  = 8; diag_phys_inst_vars_h1_1d%vertp_level%pos = 9
       diag_phys_inst_vars_h1_1d%hyam%pos        = 8; diag_phys_inst_vars_h1_1d%hybm%pos        = 8
       diag_phys_inst_vars_h1_1d%hyai%pos        = 9; diag_phys_inst_vars_h1_1d%hybi%pos        = 9

       do ilev = 1, nlev
          diag_phys_inst_vars_h1_1d%vert_level%f(ilev) = ilev
          diag_phys_inst_vars_h1_1d%vertp_level%f(ilev)= ilev
       end do
       diag_phys_inst_vars_h1_1d%vertp_level%f(nlev+1) = nlev+1
       diag_phys_inst_vars_h1_1d%hyam%f(1:nlev)  = eta_full_a; diag_phys_inst_vars_h1_1d%hybm%f(1:nlev)  = eta_full_b
       diag_phys_inst_vars_h1_1d%hyai%f(1:nlevp) = eta_face_a; diag_phys_inst_vars_h1_1d%hybi%f(1:nlevp) = eta_face_b

       call wrap_allocate_accu_vars_h1_1d
       call wrap_allocate_accu_vars_h1_2d
       call wrap_allocate_inst_vars_h1_1d
       call wrap_allocate_inst_vars_h1_2d
    end if
    ncell = mesh%nv_halo(1)

contains
    subroutine wrap_allocate_accu_vars_h1_1d
! 1d
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%precc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%precl)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%snowl)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%grapl)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%prect)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%ps)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%ts)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%shflx)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%qflx)

       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwut)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwdt)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwus)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwds)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwutc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwdtc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwusc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%flwdsc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%lwcf)

       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswut)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswdt)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswus)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswds)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswutc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswdtc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswusc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%fswdsc)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%swcf)

#ifdef AMIPW_PHYSICS
! accumulated
#ifndef DYAMOND
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%asdir)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%asdif)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%aldir)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%aldif)
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%emiss)
#endif
       call wrap_allocate_data1d(mesh,     diag_phys_accu_vars_h1_1d%qsfc)
#endif
       return
    end subroutine wrap_allocate_accu_vars_h1_1d

    subroutine wrap_allocate_accu_vars_h1_2d
! dynamics 2d
#ifndef DYAMOND
       call wrap_allocate_data2d(mesh,nlev ,diag_phys_accu_vars_h1_2d%uwind)
       call wrap_allocate_data2d(mesh,nlev ,diag_phys_accu_vars_h1_2d%vwind)
       call wrap_allocate_data2d(mesh,nlev ,diag_phys_accu_vars_h1_2d%omegaFull)
       call wrap_allocate_data2d(mesh,nlev ,diag_phys_accu_vars_h1_2d%temp)
       call wrap_allocate_data2d(mesh,nlev ,diag_phys_accu_vars_h1_2d%qv)
       call wrap_allocate_data2d(mesh,nlevp,diag_phys_accu_vars_h1_2d%mpressureFace)
       call wrap_allocate_data2d(mesh,nlevp,diag_phys_accu_vars_h1_2d%phiFace)
#ifdef DYCORE_NHD
       call wrap_allocate_data2d(mesh,nlevp,diag_phys_accu_vars_h1_2d%wwwFace)
#endif
#endif
#ifdef AMIPW_PHYSICS
       call wrap_allocate_data2d(mesh,nlev,diag_phys_accu_vars_h1_2d%rad_thten)
#ifndef DYAMOND
       call wrap_allocate_data2d(mesh,nlev,diag_phys_accu_vars_h1_2d%qc)
       call wrap_allocate_data2d(mesh,nlev,diag_phys_accu_vars_h1_2d%qr)
       call wrap_allocate_data2d(mesh,nlev,diag_phys_accu_vars_h1_2d%qi)
       call wrap_allocate_data2d(mesh,nlev,diag_phys_accu_vars_h1_2d%qs)
       call wrap_allocate_data2d(mesh,nlev,diag_phys_accu_vars_h1_2d%qg)
#endif
#endif
       return   
    end subroutine wrap_allocate_accu_vars_h1_2d

    subroutine wrap_allocate_inst_vars_h1_1d
#ifdef AMIPW_PHYSICS
! instant
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%viqv)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%viqc)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%viqr)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%viqi)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%viqs)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%viqg)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%cldtot)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%cldlow)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%cldmed)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%cldhgh)

       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%u10m )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%v10m )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%th2m )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%t2m  )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%q2m  )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%tskin)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%snowh)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%xice )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%hfx  )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%qfx  )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%lh   )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%qsfc )
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%emiss)

       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%zzz200)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%zzz500)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%zzz700)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%zzz850)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%omega200)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%omega500)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%omega700)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%omega850)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%relhum200)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%relhum500)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%relhum700)
       call wrap_allocate_data1d(mesh,  diag_phys_inst_vars_h1_1d%relhum850)
#endif

     return
    end subroutine wrap_allocate_inst_vars_h1_1d

    subroutine wrap_allocate_inst_vars_h1_2d
#ifdef AMIPW_PHYSICS
       call wrap_allocate_data2d(mesh,nlev,diag_phys_inst_vars_h1_2d%cloud)
#ifndef DYAMOND
       call wrap_allocate_data2d(mesh,nlev,diag_phys_inst_vars_h1_2d%cldcu)
       call wrap_allocate_data2d(mesh,nlev,diag_phys_inst_vars_h1_2d%cldst)
#endif
       call wrap_allocate_data2d(mesh,nlev,diag_phys_inst_vars_h1_2d%relhum)
       call wrap_allocate_data2d(mesh,nlev,diag_phys_inst_vars_h1_2d%refl_10cm)
#endif
       return
    end subroutine wrap_allocate_inst_vars_h1_2d

  end subroutine gcm_h1_diagnose_init

  subroutine gcm_h1_diagnose_final

    if(allocated(diag_phys_inst_vars_h1_1d%vert_level%f))  deallocate(diag_phys_inst_vars_h1_1d%vert_level%f)
    if(allocated(diag_phys_inst_vars_h1_1d%vertp_level%f)) deallocate(diag_phys_inst_vars_h1_1d%vertp_level%f)
    if(allocated(diag_phys_inst_vars_h1_1d%hyam%f))        deallocate(diag_phys_inst_vars_h1_1d%hyam%f)
    if(allocated(diag_phys_inst_vars_h1_1d%hybm%f))        deallocate(diag_phys_inst_vars_h1_1d%hybm%f)
    if(allocated(diag_phys_inst_vars_h1_1d%hyai%f))        deallocate(diag_phys_inst_vars_h1_1d%hyai%f)
    if(allocated(diag_phys_inst_vars_h1_1d%hybi%f))        deallocate(diag_phys_inst_vars_h1_1d%hybi%f)

    call wrap_deallocate_accu_vars_h1_1d
    call wrap_deallocate_accu_vars_h1_2d
    call wrap_deallocate_inst_vars_h1_1d
    call wrap_deallocate_inst_vars_h1_2d

contains

    subroutine wrap_deallocate_accu_vars_h1_1d

    if(allocated(diag_phys_accu_vars_h1_1d%precc%f))  deallocate(diag_phys_accu_vars_h1_1d%precc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%precl%f))  deallocate(diag_phys_accu_vars_h1_1d%precl%f)
    if(allocated(diag_phys_accu_vars_h1_1d%snowl%f))  deallocate(diag_phys_accu_vars_h1_1d%snowl%f)
    if(allocated(diag_phys_accu_vars_h1_1d%grapl%f))  deallocate(diag_phys_accu_vars_h1_1d%grapl%f)
    if(allocated(diag_phys_accu_vars_h1_1d%prect%f))  deallocate(diag_phys_accu_vars_h1_1d%prect%f)
    if(allocated(diag_phys_accu_vars_h1_1d%ps%f))     deallocate(diag_phys_accu_vars_h1_1d%ps%f)
    if(allocated(diag_phys_accu_vars_h1_1d%ts%f))     deallocate(diag_phys_accu_vars_h1_1d%ts%f)
    if(allocated(diag_phys_accu_vars_h1_1d%shflx%f))  deallocate(diag_phys_accu_vars_h1_1d%shflx%f)
    if(allocated(diag_phys_accu_vars_h1_1d%qflx%f))   deallocate(diag_phys_accu_vars_h1_1d%qflx%f)

    if(allocated(diag_phys_accu_vars_h1_1d%flwut%f))  deallocate(diag_phys_accu_vars_h1_1d%flwut%f )
    if(allocated(diag_phys_accu_vars_h1_1d%flwdt%f))  deallocate(diag_phys_accu_vars_h1_1d%flwdt%f )
    if(allocated(diag_phys_accu_vars_h1_1d%flwus%f))  deallocate(diag_phys_accu_vars_h1_1d%flwus%f )
    if(allocated(diag_phys_accu_vars_h1_1d%flwds%f))  deallocate(diag_phys_accu_vars_h1_1d%flwds%f )
    if(allocated(diag_phys_accu_vars_h1_1d%flwutc%f)) deallocate(diag_phys_accu_vars_h1_1d%flwutc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%flwdtc%f)) deallocate(diag_phys_accu_vars_h1_1d%flwdtc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%flwusc%f)) deallocate(diag_phys_accu_vars_h1_1d%flwusc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%flwdsc%f)) deallocate(diag_phys_accu_vars_h1_1d%flwdsc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%lwcf%f))   deallocate(diag_phys_accu_vars_h1_1d%lwcf%f)

    if(allocated(diag_phys_accu_vars_h1_1d%fswut%f))  deallocate(diag_phys_accu_vars_h1_1d%fswut%f )
    if(allocated(diag_phys_accu_vars_h1_1d%fswdt%f))  deallocate(diag_phys_accu_vars_h1_1d%fswdt%f )
    if(allocated(diag_phys_accu_vars_h1_1d%fswus%f))  deallocate(diag_phys_accu_vars_h1_1d%fswus%f )
    if(allocated(diag_phys_accu_vars_h1_1d%fswds%f))  deallocate(diag_phys_accu_vars_h1_1d%fswds%f )
    if(allocated(diag_phys_accu_vars_h1_1d%fswutc%f)) deallocate(diag_phys_accu_vars_h1_1d%fswutc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%fswdtc%f)) deallocate(diag_phys_accu_vars_h1_1d%fswdtc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%fswusc%f)) deallocate(diag_phys_accu_vars_h1_1d%fswusc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%fswdsc%f)) deallocate(diag_phys_accu_vars_h1_1d%fswdsc%f)
    if(allocated(diag_phys_accu_vars_h1_1d%swcf%f))   deallocate(diag_phys_accu_vars_h1_1d%swcf%f)

#ifdef AMIPW_PHYSICS
#ifndef DYAMOND
    if(allocated(diag_phys_accu_vars_h1_1d%asdir%f))     deallocate(diag_phys_accu_vars_h1_1d%asdir%f)
    if(allocated(diag_phys_accu_vars_h1_1d%asdif%f))     deallocate(diag_phys_accu_vars_h1_1d%asdif%f)
    if(allocated(diag_phys_accu_vars_h1_1d%aldir%f))     deallocate(diag_phys_accu_vars_h1_1d%aldir%f)
    if(allocated(diag_phys_accu_vars_h1_1d%aldif%f))     deallocate(diag_phys_accu_vars_h1_1d%aldif%f)
    if(allocated(diag_phys_accu_vars_h1_1d%emiss%f))     deallocate(diag_phys_accu_vars_h1_1d%emiss%f)
#endif
    if(allocated(diag_phys_accu_vars_h1_1d%qsfc%f))      deallocate(diag_phys_accu_vars_h1_1d%qsfc%f)
#endif
       return
    end subroutine wrap_deallocate_accu_vars_h1_1d

    subroutine wrap_deallocate_accu_vars_h1_2d
#ifndef DYAMOND 
    if(allocated(diag_phys_accu_vars_h1_2d%uwind%f))         deallocate(diag_phys_accu_vars_h1_2d%uwind%f)
    if(allocated(diag_phys_accu_vars_h1_2d%vwind%f))         deallocate(diag_phys_accu_vars_h1_2d%vwind%f)
    if(allocated(diag_phys_accu_vars_h1_2d%omegaFull%f))     deallocate(diag_phys_accu_vars_h1_2d%omegaFull%f)
    if(allocated(diag_phys_accu_vars_h1_2d%temp%f))          deallocate(diag_phys_accu_vars_h1_2d%temp%f)
    if(allocated(diag_phys_accu_vars_h1_2d%qv%f))            deallocate(diag_phys_accu_vars_h1_2d%qv%f)
    if(allocated(diag_phys_accu_vars_h1_2d%phiFace%f))       deallocate(diag_phys_accu_vars_h1_2d%phiFace%f)
    if(allocated(diag_phys_accu_vars_h1_2d%mpressureFace%f)) deallocate(diag_phys_accu_vars_h1_2d%mpressureFace%f)
#ifdef DYCORE_NHD
    if(allocated(diag_phys_accu_vars_h1_2d%wwwFace%f))       deallocate(diag_phys_accu_vars_h1_2d%wwwFace%f)
#endif
#endif

#ifdef AMIPW_PHYSICS
    if(allocated(diag_phys_accu_vars_h1_2d%rad_thten%f)) deallocate(diag_phys_accu_vars_h1_2d%rad_thten%f)
#ifndef DYAMOND
    if(allocated(diag_phys_accu_vars_h1_2d%qc%f))        deallocate(diag_phys_accu_vars_h1_2d%qc%f)
    if(allocated(diag_phys_accu_vars_h1_2d%qr%f))        deallocate(diag_phys_accu_vars_h1_2d%qr%f)
    if(allocated(diag_phys_accu_vars_h1_2d%qi%f))        deallocate(diag_phys_accu_vars_h1_2d%qi%f)
    if(allocated(diag_phys_accu_vars_h1_2d%qs%f))        deallocate(diag_phys_accu_vars_h1_2d%qs%f)
    if(allocated(diag_phys_accu_vars_h1_2d%qg%f))        deallocate(diag_phys_accu_vars_h1_2d%qg%f)
#endif
#endif
       return
    end subroutine wrap_deallocate_accu_vars_h1_2d

    subroutine wrap_deallocate_inst_vars_h1_1d
       if(allocated(diag_phys_inst_vars_h1_1d%viqv%f))   deallocate(diag_phys_inst_vars_h1_1d%viqv%f)
       if(allocated(diag_phys_inst_vars_h1_1d%viqc%f))   deallocate(diag_phys_inst_vars_h1_1d%viqc%f)
       if(allocated(diag_phys_inst_vars_h1_1d%viqr%f))   deallocate(diag_phys_inst_vars_h1_1d%viqr%f)
       if(allocated(diag_phys_inst_vars_h1_1d%viqi%f))   deallocate(diag_phys_inst_vars_h1_1d%viqi%f)
       if(allocated(diag_phys_inst_vars_h1_1d%viqs%f))   deallocate(diag_phys_inst_vars_h1_1d%viqs%f)
       if(allocated(diag_phys_inst_vars_h1_1d%viqg%f))   deallocate(diag_phys_inst_vars_h1_1d%viqg%f)
       if(allocated(diag_phys_inst_vars_h1_1d%cldtot%f)) deallocate(diag_phys_inst_vars_h1_1d%cldtot%f)
       if(allocated(diag_phys_inst_vars_h1_1d%cldlow%f)) deallocate(diag_phys_inst_vars_h1_1d%cldlow%f)
       if(allocated(diag_phys_inst_vars_h1_1d%cldmed%f)) deallocate(diag_phys_inst_vars_h1_1d%cldmed%f)
       if(allocated(diag_phys_inst_vars_h1_1d%cldhgh%f)) deallocate(diag_phys_inst_vars_h1_1d%cldhgh%f)

       if(allocated(diag_phys_inst_vars_h1_1d%u10m%f))   deallocate(diag_phys_inst_vars_h1_1d%u10m%f)
       if(allocated(diag_phys_inst_vars_h1_1d%v10m%f))   deallocate(diag_phys_inst_vars_h1_1d%v10m%f)
       if(allocated(diag_phys_inst_vars_h1_1d%th2m%f))   deallocate(diag_phys_inst_vars_h1_1d%th2m%f)
       if(allocated(diag_phys_inst_vars_h1_1d%t2m%f))    deallocate(diag_phys_inst_vars_h1_1d%t2m%f)
       if(allocated(diag_phys_inst_vars_h1_1d%q2m%f))    deallocate(diag_phys_inst_vars_h1_1d%q2m%f)
       if(allocated(diag_phys_inst_vars_h1_1d%tskin%f))  deallocate(diag_phys_inst_vars_h1_1d%tskin%f)
       if(allocated(diag_phys_inst_vars_h1_1d%snowh%f))  deallocate(diag_phys_inst_vars_h1_1d%snowh%f)
       if(allocated(diag_phys_inst_vars_h1_1d%xice%f))   deallocate(diag_phys_inst_vars_h1_1d%xice%f)
       if(allocated(diag_phys_inst_vars_h1_1d%hfx%f))    deallocate(diag_phys_inst_vars_h1_1d%hfx%f)
       if(allocated(diag_phys_inst_vars_h1_1d%qfx%f))    deallocate(diag_phys_inst_vars_h1_1d%qfx%f)
       if(allocated(diag_phys_inst_vars_h1_1d%lh%f))     deallocate(diag_phys_inst_vars_h1_1d%lh%f)
       if(allocated(diag_phys_inst_vars_h1_1d%qsfc%f))   deallocate(diag_phys_inst_vars_h1_1d%qsfc%f)
       if(allocated(diag_phys_inst_vars_h1_1d%emiss%f))  deallocate(diag_phys_inst_vars_h1_1d%emiss%f) 

       if(allocated(diag_phys_inst_vars_h1_1d%zzz200%f))     deallocate(diag_phys_inst_vars_h1_1d%zzz200%f)
       if(allocated(diag_phys_inst_vars_h1_1d%zzz500%f))     deallocate(diag_phys_inst_vars_h1_1d%zzz500%f)
       if(allocated(diag_phys_inst_vars_h1_1d%zzz700%f))     deallocate(diag_phys_inst_vars_h1_1d%zzz700%f)
       if(allocated(diag_phys_inst_vars_h1_1d%zzz850%f))     deallocate(diag_phys_inst_vars_h1_1d%zzz850%f)
       if(allocated(diag_phys_inst_vars_h1_1d%omega200%f))   deallocate(diag_phys_inst_vars_h1_1d%omega200%f)
       if(allocated(diag_phys_inst_vars_h1_1d%omega500%f))   deallocate(diag_phys_inst_vars_h1_1d%omega500%f)
       if(allocated(diag_phys_inst_vars_h1_1d%omega700%f))   deallocate(diag_phys_inst_vars_h1_1d%omega700%f)
       if(allocated(diag_phys_inst_vars_h1_1d%omega850%f))   deallocate(diag_phys_inst_vars_h1_1d%omega850%f)
       if(allocated(diag_phys_inst_vars_h1_1d%relhum200%f))  deallocate(diag_phys_inst_vars_h1_1d%relhum200%f)
       if(allocated(diag_phys_inst_vars_h1_1d%relhum500%f))  deallocate(diag_phys_inst_vars_h1_1d%relhum500%f)
       if(allocated(diag_phys_inst_vars_h1_1d%relhum700%f))  deallocate(diag_phys_inst_vars_h1_1d%relhum700%f)
       if(allocated(diag_phys_inst_vars_h1_1d%relhum850%f))  deallocate(diag_phys_inst_vars_h1_1d%relhum850%f)

       return
    end subroutine wrap_deallocate_inst_vars_h1_1d

    subroutine wrap_deallocate_inst_vars_h1_2d
       if(allocated(diag_phys_inst_vars_h1_2d%cloud%f))      deallocate(diag_phys_inst_vars_h1_2d%cloud%f)
#ifndef DYAMOND
       if(allocated(diag_phys_inst_vars_h1_2d%cldcu%f))      deallocate(diag_phys_inst_vars_h1_2d%cldcu%f)
       if(allocated(diag_phys_inst_vars_h1_2d%cldst%f))      deallocate(diag_phys_inst_vars_h1_2d%cldst%f)
#endif
       if(allocated(diag_phys_inst_vars_h1_2d%relhum%f))     deallocate(diag_phys_inst_vars_h1_2d%relhum%f)
       if(allocated(diag_phys_inst_vars_h1_2d%refl_10cm%f))  deallocate(diag_phys_inst_vars_h1_2d%refl_10cm%f)
     return
    end subroutine wrap_deallocate_inst_vars_h1_2d

  end subroutine gcm_h1_diagnose_final
!
! diag instant var for the current call, physics does not mean "model physics"
!
  subroutine gcm_h1_1d_inst_physics_variables
   integer(i4) :: icell, ilev
   real(r8) :: sat_specific_humidity(nlev,ncell), esvp(nlev,ncell)
   real(r8) :: pres_out(4), vartmp_out(4)

   ! diag_phys_inst_vars_h1_1d%vert_level%pos  = 8; diag_phys_inst_vars_h1_1d%vertp_level%pos = 9
   ! diag_phys_inst_vars_h1_1d%hyam%pos        = 8; diag_phys_inst_vars_h1_1d%hybm%pos        = 8
   ! diag_phys_inst_vars_h1_1d%hyai%pos        = 9; diag_phys_inst_vars_h1_1d%hybi%pos        = 9

   ! do ilev = 1, nlev
   !    diag_phys_inst_vars_h1_1d%vert_level%f(ilev)  = ilev
   !    diag_phys_inst_vars_h1_1d%vertp_level%f(ilev) = ilev
   ! end do
   ! diag_phys_inst_vars_h1_1d%vertp_level%f(nlev+1) = nlev+1

   ! diag_phys_inst_vars_h1_1d%hyam%f(1:nlev)  = eta_full_a; diag_phys_inst_vars_h1_1d%hybm%f(1:nlev)  = eta_full_b
   ! diag_phys_inst_vars_h1_1d%hyai%f(1:nlevp) = eta_face_a; diag_phys_inst_vars_h1_1d%hybi%f(1:nlevp) = eta_face_b
#ifdef AMIPW_PHYSICS
! ATM 1D
    do icell = 1, ncell
       diag_phys_inst_vars_h1_1d%viqv%f(icell)     = zero
       diag_phys_inst_vars_h1_1d%viqc%f(icell)     = zero
       diag_phys_inst_vars_h1_1d%viqr%f(icell)     = zero
       diag_phys_inst_vars_h1_1d%viqi%f(icell)     = zero
       diag_phys_inst_vars_h1_1d%viqs%f(icell)     = zero
       diag_phys_inst_vars_h1_1d%viqg%f(icell)     = zero
       do ilev = 1, nlev
          diag_phys_inst_vars_h1_1d%viqv%f(icell)  = diag_phys_inst_vars_h1_1d%viqv%f(icell) + tracerVarCellFull%scalar_tracer_mass_n%f(1,ilev,icell)/gravity
          diag_phys_inst_vars_h1_1d%viqc%f(icell)  = diag_phys_inst_vars_h1_1d%viqc%f(icell) + tracerVarCellFull%scalar_tracer_mass_n%f(2,ilev,icell)/gravity
          diag_phys_inst_vars_h1_1d%viqr%f(icell)  = diag_phys_inst_vars_h1_1d%viqr%f(icell) + tracerVarCellFull%scalar_tracer_mass_n%f(3,ilev,icell)/gravity
          diag_phys_inst_vars_h1_1d%viqi%f(icell)  = diag_phys_inst_vars_h1_1d%viqi%f(icell) + tracerVarCellFull%scalar_tracer_mass_n%f(4,ilev,icell)/gravity
          diag_phys_inst_vars_h1_1d%viqs%f(icell)  = diag_phys_inst_vars_h1_1d%viqs%f(icell) + tracerVarCellFull%scalar_tracer_mass_n%f(5,ilev,icell)/gravity
          diag_phys_inst_vars_h1_1d%viqg%f(icell)  = diag_phys_inst_vars_h1_1d%viqg%f(icell) + tracerVarCellFull%scalar_tracer_mass_n%f(6,ilev,icell)/gravity
       end do
    end do
! cloud fraction

    diag_phys_inst_vars_h1_1d%cldtot%f(1:ncell)     = pstate_wrf%cldtot(1:ncell,1)
    diag_phys_inst_vars_h1_1d%cldlow%f(1:ncell)     = pstate_wrf%cldlow(1:ncell,1)
    diag_phys_inst_vars_h1_1d%cldmed%f(1:ncell)     = pstate_wrf%cldmed(1:ncell,1)
    diag_phys_inst_vars_h1_1d%cldhgh%f(1:ncell)     = pstate_wrf%cldhgh(1:ncell,1)

! LND 1d
    diag_phys_inst_vars_h1_1d%u10m%f(1:ncell) = psurf_wrf%u10 (1:ncell,1)
    diag_phys_inst_vars_h1_1d%v10m%f(1:ncell) = psurf_wrf%v10 (1:ncell,1)
    diag_phys_inst_vars_h1_1d%th2m%f(1:ncell) = psurf_wrf%th2 (1:ncell,1)
    diag_phys_inst_vars_h1_1d%t2m%f (1:ncell) = psurf_wrf%t2  (1:ncell,1)
    diag_phys_inst_vars_h1_1d%q2m%f (1:ncell) = psurf_wrf%q2  (1:ncell,1)
    diag_phys_inst_vars_h1_1d%tskin%f(1:ncell)= psurf_wrf%tsk (1:ncell,1)
    diag_phys_inst_vars_h1_1d%snowh%f(1:ncell)= psurf_wrf%snow(1:ncell,1) ! already an accumulated var,just output, no need to accumulated over time
    diag_phys_inst_vars_h1_1d%xice%f(1:ncell) = psurf_wrf%xice(1:ncell,1)
    diag_phys_inst_vars_h1_1d%hfx%f (1:ncell) = psurf_wrf%hfx (1:ncell,1)
    diag_phys_inst_vars_h1_1d%qfx%f (1:ncell) = psurf_wrf%qfx (1:ncell,1)
    diag_phys_inst_vars_h1_1d%lh%f  (1:ncell) = psurf_wrf%lh  (1:ncell,1)
    diag_phys_inst_vars_h1_1d%qsfc%f(1:ncell) = psurf_wrf%qsfc(1:ncell,1)
    diag_phys_inst_vars_h1_1d%emiss%f(1:ncell)= psurf_wrf%emiss(1:ncell,1)

!================================================================================================
! obtain RELHUM. We use state from dynamics, not physics, because some physics may update 
! t but not p, so a little bit inconsistency (one-step impact is trivial). Use dynamical state
! is comparible with other vars, e.g., geop, omega 
!=================================================================================================

    call aqsat_grist(t= dycoreVarCellFull%scalar_temp_n%f(1:nlev,1:ncell), &
                     p= dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,1:ncell), &
                     es=esvp(1:nlev,1:ncell)              , &
                     qs=sat_specific_humidity(1:nlev,1:ncell), &
                     ii=ncell, ilen=ncell, kk=nlev, kstart=1, kend=nlev)

    do icell = 1, ncell
       do ilev = 1, nlev
          diag_phys_inst_vars_h1_2d%relhum%f(ilev,icell) = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,icell)&
                                               /(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,icell)))/sat_specific_humidity(ilev,icell)*100._r8
       end do
    end do

! verticall interpolte omega, relhum, zzz
   pres_out(1) = 20000._r8
   pres_out(2) = 50000._r8
   pres_out(3) = 70000._r8
   pres_out(4) = 85000._r8

   do icell = 1, ncell
       call lininterp(dycoreVarCellFull%scalar_omega_n%f(1:nlev,icell),&
                      log(log(dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,icell))),1,nlev,vartmp_out(1:4),log(log(pres_out)),4)
       diag_phys_inst_vars_h1_1d%omega200%f(icell) = vartmp_out(1)
       diag_phys_inst_vars_h1_1d%omega500%f(icell) = vartmp_out(2)
       diag_phys_inst_vars_h1_1d%omega700%f(icell) = vartmp_out(3)
       diag_phys_inst_vars_h1_1d%omega850%f(icell) = vartmp_out(4)

       call lininterp(dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,icell),&
                      log(log(dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,icell))),1,nlev,vartmp_out(1:4),log(log(pres_out)),4)
       diag_phys_inst_vars_h1_1d%zzz200%f(icell) = vartmp_out(1)
       diag_phys_inst_vars_h1_1d%zzz500%f(icell) = vartmp_out(2)
       diag_phys_inst_vars_h1_1d%zzz700%f(icell) = vartmp_out(3)
       diag_phys_inst_vars_h1_1d%zzz850%f(icell) = vartmp_out(4)

       call lininterp(diag_phys_inst_vars_h1_2d%relhum%f(1:nlev,icell),&
                      log(log(dycoreVarCellFull%scalar_mpressure_n%f(1:nlev,icell))),1,nlev,vartmp_out(1:4),log(log(pres_out)),4)
       diag_phys_inst_vars_h1_1d%relhum200%f(icell) = vartmp_out(1)
       diag_phys_inst_vars_h1_1d%relhum500%f(icell) = vartmp_out(2)
       diag_phys_inst_vars_h1_1d%relhum700%f(icell) = vartmp_out(3)
       diag_phys_inst_vars_h1_1d%relhum850%f(icell) = vartmp_out(4)
   end do

#endif
    return
  end subroutine gcm_h1_1d_inst_physics_variables

  subroutine gcm_h1_2d_inst_physics_variables
   integer(i4) :: icell, ilev
#ifdef AMIPW_PHYSICS
    do icell = 1, ncell
       do ilev = 1, nlev
          diag_phys_inst_vars_h1_2d%cloud%f(ilev,icell)  = pstate_wrf%cldfra(icell,nlev+1-ilev,1) 
#ifndef DYAMOND
          diag_phys_inst_vars_h1_2d%cldcu%f(ilev,icell)  = pstate_wrf%cldcum(icell,nlev+1-ilev,1)
          diag_phys_inst_vars_h1_2d%cldst%f(ilev,icell)  = pstate_wrf%cldstr(icell,nlev+1-ilev,1)
#endif
          diag_phys_inst_vars_h1_2d%refl_10cm%f(ilev,icell)  = pstate_wrf%refl_10cm(icell,nlev+1-ilev,1)
       end do
    end do
#endif
    return
  end subroutine gcm_h1_2d_inst_physics_variables
!
! accumulate at each model step
!
  subroutine gcm_h1_1d_accu_physics_variables
   integer(i4) :: icell, ilev
! local

    if(write_history_h1)then

    diag_phys_accu_vars_h1_1d%precl%f = diag_phys_accu_vars_h1_1d%precl%f+pstate%scalar_precl_surface%f
    diag_phys_accu_vars_h1_1d%snowl%f = diag_phys_accu_vars_h1_1d%snowl%f+pstate%scalar_snowl_surface%f
    diag_phys_accu_vars_h1_1d%grapl%f = diag_phys_accu_vars_h1_1d%grapl%f+pstate%scalar_grapl_surface%f
    diag_phys_accu_vars_h1_1d%precc%f = diag_phys_accu_vars_h1_1d%precc%f+pstate%scalar_precc_surface%f
    diag_phys_accu_vars_h1_1d%prect%f = diag_phys_accu_vars_h1_1d%prect%f+pstate%scalar_prect_surface%f
    diag_phys_accu_vars_h1_1d%ps%f    = diag_phys_accu_vars_h1_1d%ps%f   +dycoreVarCellFace%scalar_pressure_n%f(nlevp,:)
    diag_phys_accu_vars_h1_1d%ts%f    = diag_phys_accu_vars_h1_1d%ts%f   +pstate%ts_at_pc_surface%f
#ifdef AMIPW_PHYSICS
    diag_phys_accu_vars_h1_1d%shflx%f = diag_phys_accu_vars_h1_1d%shflx%f+pstate%atm_in_shflx_at_pc_surface%f
    diag_phys_accu_vars_h1_1d%qflx%f  = diag_phys_accu_vars_h1_1d%qflx%f +pstate%atm_in_qflx_at_pc_surface%f(1,:)
#endif

#ifdef AMIPC_PHYSICS
    !diag_physics_vars_h1%cldtot%f= diag_physics_vars_h1%cldtot%f+pstate_cam%diag_cloud_tot%f(:)
    !diag_physics_vars_h1%cldlow%f= diag_physics_vars_h1%cldlow%f+pstate_cam%diag_cloud_low%f(:)
    !diag_physics_vars_h1%cldmed%f= diag_physics_vars_h1%cldmed%f+pstate_cam%diag_cloud_med%f(:)
    !diag_physics_vars_h1%cldhgh%f= diag_physics_vars_h1%cldhgh%f+pstate_cam%diag_cloud_hgh%f(:)
    !diag_physics_vars_h1%z500%f  = diag_physics_vars_h1%z500%f+pstate_cam%diag_z_at_500hpa%f
    !diag_physics_vars_h1%u850%f  = diag_physics_vars_h1%u850%f+pstate_cam%diag_u_at_850hpa%f
    !diag_physics_vars_h1%u200%f  = diag_physics_vars_h1%u200%f+pstate_cam%diag_u_at_200hpa%f
    !diag_physics_vars_h1%v850%f  = diag_physics_vars_h1%v850%f+pstate_cam%diag_v_at_850hpa%f
 
    diag_phys_accu_vars_h1_1d%flwut%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwut%f(1:ncell) +pstate_cam%flwut_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%flwdt%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwdt%f(1:ncell) +pstate_cam%flwdt_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%flwus%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwus%f(1:ncell) +pstate_cam%flwus_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%flwds%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwds%f(1:ncell) +pstate_cam%flwds_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%flwutc%f(1:ncell)= diag_phys_accu_vars_h1_1d%flwutc%f(1:ncell)+pstate_cam%flwutc_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%flwdtc%f(1:ncell)= diag_phys_accu_vars_h1_1d%flwdtc%f(1:ncell)+pstate_cam%flwdtc_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%flwusc%f(1:ncell)= diag_phys_accu_vars_h1_1d%flwusc%f(1:ncell)+pstate_cam%flwusc_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%flwdsc%f(1:ncell)= diag_phys_accu_vars_h1_1d%flwdsc%f(1:ncell)+pstate_cam%flwdsc_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%lwcf%f  (1:ncell)= diag_phys_accu_vars_h1_1d%lwcf%f(1:ncell)  +pstate_cam%lwcf_at_pc_top%f(1:ncell)

    diag_phys_accu_vars_h1_1d%fswut%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswut%f(1:ncell) +pstate_cam%fswut_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%fswdt%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswdt%f(1:ncell) +pstate_cam%fswdt_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%fswus%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswus%f(1:ncell) +pstate_cam%fswus_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%fswds%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswds%f(1:ncell) +pstate_cam%fswds_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%fswutc%f(1:ncell)= diag_phys_accu_vars_h1_1d%fswutc%f(1:ncell)+pstate_cam%fswutc_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%fswdtc%f(1:ncell)= diag_phys_accu_vars_h1_1d%fswdtc%f(1:ncell)+pstate_cam%fswdtc_at_pc_top%f(1:ncell)
    diag_phys_accu_vars_h1_1d%fswusc%f(1:ncell)= diag_phys_accu_vars_h1_1d%fswusc%f(1:ncell)+pstate_cam%fswusc_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%fswdsc%f(1:ncell)= diag_phys_accu_vars_h1_1d%fswdsc%f(1:ncell)+pstate_cam%fswdsc_at_pc_surface%f(1:ncell)
    diag_phys_accu_vars_h1_1d%swcf%f  (1:ncell)= diag_phys_accu_vars_h1_1d%swcf%f(1:ncell)  +pstate_cam%swcf_at_pc_top%f(1:ncell)
#endif

#ifdef AMIPW_PHYSICS
    diag_phys_accu_vars_h1_1d%flwut%f( 1:ncell) = diag_phys_accu_vars_h1_1d%flwut%f (1:ncell) + pstate_wrf%lwupt(1:ncell,1)
    diag_phys_accu_vars_h1_1d%flwdt%f( 1:ncell) = diag_phys_accu_vars_h1_1d%flwdt%f (1:ncell) + pstate_wrf%lwdnt(1:ncell,1)
    diag_phys_accu_vars_h1_1d%flwus%f( 1:ncell) = diag_phys_accu_vars_h1_1d%flwus%f (1:ncell) + pstate_wrf%lwupb(1:ncell,1)
    diag_phys_accu_vars_h1_1d%flwds%f( 1:ncell) = diag_phys_accu_vars_h1_1d%flwds%f (1:ncell) + pstate_wrf%lwdnb(1:ncell,1)
    diag_phys_accu_vars_h1_1d%flwutc%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwutc%f(1:ncell) + pstate_wrf%lwuptc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%flwdtc%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwdtc%f(1:ncell) + pstate_wrf%lwdntc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%flwusc%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwusc%f(1:ncell) + pstate_wrf%lwupbc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%flwdsc%f(1:ncell) = diag_phys_accu_vars_h1_1d%flwdsc%f(1:ncell) + pstate_wrf%lwdnbc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%lwcf%f  (1:ncell) = diag_phys_accu_vars_h1_1d%lwcf%f(1:ncell)   + pstate_wrf%lwcf(1:ncell,1)

    diag_phys_accu_vars_h1_1d%fswut%f( 1:ncell) = diag_phys_accu_vars_h1_1d%fswut%f (1:ncell) + pstate_wrf%swupt(1:ncell,1)
    diag_phys_accu_vars_h1_1d%fswdt%f( 1:ncell) = diag_phys_accu_vars_h1_1d%fswdt%f (1:ncell) + pstate_wrf%swdnt(1:ncell,1)
    diag_phys_accu_vars_h1_1d%fswus%f( 1:ncell) = diag_phys_accu_vars_h1_1d%fswus%f (1:ncell) + pstate_wrf%swupb(1:ncell,1)
    diag_phys_accu_vars_h1_1d%fswds%f( 1:ncell) = diag_phys_accu_vars_h1_1d%fswds%f (1:ncell) + pstate_wrf%swdnb(1:ncell,1)
    diag_phys_accu_vars_h1_1d%fswutc%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswutc%f(1:ncell) + pstate_wrf%swuptc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%fswdtc%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswdtc%f(1:ncell) + pstate_wrf%swdntc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%fswusc%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswusc%f(1:ncell) + pstate_wrf%swupbc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%fswdsc%f(1:ncell) = diag_phys_accu_vars_h1_1d%fswdsc%f(1:ncell) + pstate_wrf%swdnbc(1:ncell,1)
    diag_phys_accu_vars_h1_1d%swcf%f  (1:ncell) = diag_phys_accu_vars_h1_1d%swcf%f(1:ncell)   + pstate_wrf%swcf(1:ncell,1)

#ifndef DYAMOND
    diag_phys_accu_vars_h1_1d%asdir%f(1:ncell)  = diag_phys_accu_vars_h1_1d%asdir%f(1:ncell)  + psurf_wrf%asdir(1:ncell,1)
    diag_phys_accu_vars_h1_1d%asdif%f(1:ncell)  = diag_phys_accu_vars_h1_1d%asdif%f(1:ncell)  + psurf_wrf%asdif(1:ncell,1)
    diag_phys_accu_vars_h1_1d%aldir%f(1:ncell)  = diag_phys_accu_vars_h1_1d%aldir%f(1:ncell)  + psurf_wrf%aldir(1:ncell,1)
    diag_phys_accu_vars_h1_1d%aldif%f(1:ncell)  = diag_phys_accu_vars_h1_1d%aldif%f(1:ncell)  + psurf_wrf%aldif(1:ncell,1)
    diag_phys_accu_vars_h1_1d%emiss%f(1:ncell)  = diag_phys_accu_vars_h1_1d%emiss%f(1:ncell)  + psurf_wrf%emiss(1:ncell,1)
#endif
    diag_phys_accu_vars_h1_1d%qsfc%f(1:ncell)   = diag_phys_accu_vars_h1_1d%qsfc%f(1:ncell)   + psurf_wrf%qsfc (1:ncell,1)
#endif

    end if

    diag_phys_accu_vars_h1_1d%ncount  = diag_phys_accu_vars_h1_1d%ncount + 1
    return
  end subroutine gcm_h1_1d_accu_physics_variables

  subroutine gcm_h1_2d_accu_physics_variables
   integer(i4) :: icell, ilev
! local

    if(write_history_h1)then
#ifndef DYAMOND
    diag_phys_accu_vars_h1_2d%uwind%f          = diag_phys_accu_vars_h1_2d%uwind%f         +dycoreVarCellFull%scalar_U_wind_n%f
    diag_phys_accu_vars_h1_2d%vwind%f          = diag_phys_accu_vars_h1_2d%vwind%f         +dycoreVarCellFull%scalar_V_wind_n%f
    diag_phys_accu_vars_h1_2d%temp%f           = diag_phys_accu_vars_h1_2d%temp%f          +dycoreVarCellFull%scalar_temp_n%f
    diag_phys_accu_vars_h1_2d%mpressureFace%f  = diag_phys_accu_vars_h1_2d%mpressureFace%f +dycoreVarCellFace%scalar_mpressure_n%f
    diag_phys_accu_vars_h1_2d%omegaFull%f      = diag_phys_accu_vars_h1_2d%omegaFull%f     +dycoreVarCellFull%scalar_omega_n%f
#ifdef DYCORE_NHD
    diag_phys_accu_vars_h1_2d%wwwFace%f        = diag_phys_accu_vars_h1_2d%wwwFace%f       +dycoreVarCellFace%scalar_www_n%f
#endif
    diag_phys_accu_vars_h1_2d%qv%f             = diag_phys_accu_vars_h1_2d%qv%f            +tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:)
    diag_phys_accu_vars_h1_2d%phiFace%f        = diag_phys_accu_vars_h1_2d%phiFace%f       +dycoreVarCellFace%scalar_geopotential_n%f
#endif

#ifdef AMIPW_PHYSICS
    do icell = 1, ncell
       do ilev = 1, nlev
! not saved in rst, and because cldfra is called at rad step, this diag will be incorrect for the first restart-run output variable
! yizhang, 20210919, we move cloud fraction to inst diagnostics
          diag_phys_accu_vars_h1_2d%rad_thten%f(ilev,icell)= diag_phys_accu_vars_h1_2d%rad_thten%f(ilev,icell)+ptend_wrf%rthraten(icell,nlev+1-ilev,1)
       end do
    end do
! hard-coded here
#ifndef DYAMOND
    diag_phys_accu_vars_h1_2d%qc%f    = diag_phys_accu_vars_h1_2d%qc%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(2,:,:) ! 6 is default, vcrisg
    diag_phys_accu_vars_h1_2d%qr%f    = diag_phys_accu_vars_h1_2d%qr%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(3,:,:) ! 6 is default, vcrisg
    diag_phys_accu_vars_h1_2d%qi%f    = diag_phys_accu_vars_h1_2d%qi%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(4,:,:)
    diag_phys_accu_vars_h1_2d%qs%f    = diag_phys_accu_vars_h1_2d%qs%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(5,:,:)
    diag_phys_accu_vars_h1_2d%qg%f    = diag_phys_accu_vars_h1_2d%qg%f  +tracerVarCellFull%scalar_tracer_mxrt_n%f(6,:,:)
#endif
#endif
    end if
    diag_phys_accu_vars_h1_2d%ncount  = diag_phys_accu_vars_h1_2d%ncount + 1
    return
  end subroutine gcm_h1_2d_accu_physics_variables

!
! dump depending on write_history frequency
!

  subroutine gcm_h1_1d_dump_physics_variables
! local

    if(write_history_h1)then

    diag_phys_accu_vars_h1_1d%precl%f = diag_phys_accu_vars_h1_1d%precl%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%snowl%f = diag_phys_accu_vars_h1_1d%snowl%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%grapl%f = diag_phys_accu_vars_h1_1d%grapl%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%precc%f = diag_phys_accu_vars_h1_1d%precc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%prect%f = diag_phys_accu_vars_h1_1d%prect%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%ps%f    = diag_phys_accu_vars_h1_1d%ps%f   /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%ts%f    = diag_phys_accu_vars_h1_1d%ts%f   /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%shflx%f = diag_phys_accu_vars_h1_1d%shflx%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%qflx%f  = diag_phys_accu_vars_h1_1d%qflx%f /diag_phys_accu_vars_h1_1d%ncount

    diag_phys_accu_vars_h1_1d%flwut%f = diag_phys_accu_vars_h1_1d%flwut%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%flwdt%f = diag_phys_accu_vars_h1_1d%flwdt%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%flwus%f = diag_phys_accu_vars_h1_1d%flwus%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%flwds%f = diag_phys_accu_vars_h1_1d%flwds%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%flwutc%f= diag_phys_accu_vars_h1_1d%flwutc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%flwdtc%f= diag_phys_accu_vars_h1_1d%flwdtc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%flwusc%f= diag_phys_accu_vars_h1_1d%flwusc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%flwdsc%f= diag_phys_accu_vars_h1_1d%flwdsc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%lwcf%f  = diag_phys_accu_vars_h1_1d%lwcf%f  /diag_phys_accu_vars_h1_1d%ncount

    diag_phys_accu_vars_h1_1d%fswut%f = diag_phys_accu_vars_h1_1d%fswut%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%fswdt%f = diag_phys_accu_vars_h1_1d%fswdt%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%fswus%f = diag_phys_accu_vars_h1_1d%fswus%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%fswds%f = diag_phys_accu_vars_h1_1d%fswds%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%fswutc%f= diag_phys_accu_vars_h1_1d%fswutc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%fswdtc%f= diag_phys_accu_vars_h1_1d%fswdtc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%fswusc%f= diag_phys_accu_vars_h1_1d%fswusc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%fswdsc%f= diag_phys_accu_vars_h1_1d%fswdsc%f/diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%swcf%f  = diag_phys_accu_vars_h1_1d%swcf%f  /diag_phys_accu_vars_h1_1d%ncount

#ifdef AMIPC_PHYSICS
    !diag_phys_accu_vars_h1_1d%cldlow%f= diag_phys_accu_vars_h1_1d%cldlow%f/diag_phys_accu_vars_h1_1d%ncount
    !diag_phys_accu_vars_h1_1d%cldmed%f= diag_phys_accu_vars_h1_1d%cldmed%f/diag_phys_accu_vars_h1_1d%ncount
    !diag_phys_accu_vars_h1_1d%cldhgh%f= diag_phys_accu_vars_h1_1d%cldhgh%f/diag_phys_accu_vars_h1_1d%ncount
    !diag_phys_accu_vars_h1_1d%z500%f  = diag_phys_accu_vars_h1_1d%z500%f  /diag_phys_accu_vars_h1_1d%ncount
    !diag_phys_accu_vars_h1_1d%u850%f  = diag_phys_accu_vars_h1_1d%u850%f  /diag_phys_accu_vars_h1_1d%ncount
    !diag_phys_accu_vars_h1_1d%u200%f  = diag_phys_accu_vars_h1_1d%u200%f  /diag_phys_accu_vars_h1_1d%ncount
    !diag_phys_accu_vars_h1_1d%v850%f  = diag_phys_accu_vars_h1_1d%v850%f  /diag_phys_accu_vars_h1_1d%ncount
#endif

#ifdef AMIPW_PHYSICS
#ifndef DYAMOND
    diag_phys_accu_vars_h1_1d%asdir%f = diag_phys_accu_vars_h1_1d%asdir%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%asdif%f = diag_phys_accu_vars_h1_1d%asdif%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%aldir%f = diag_phys_accu_vars_h1_1d%aldir%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%aldif%f = diag_phys_accu_vars_h1_1d%aldif%f /diag_phys_accu_vars_h1_1d%ncount
    diag_phys_accu_vars_h1_1d%emiss%f = diag_phys_accu_vars_h1_1d%emiss%f /diag_phys_accu_vars_h1_1d%ncount
#endif
    diag_phys_accu_vars_h1_1d%qsfc%f  = diag_phys_accu_vars_h1_1d%qsfc%f  /diag_phys_accu_vars_h1_1d%ncount
#endif
    end if
 
    return
  end subroutine gcm_h1_1d_dump_physics_variables

  subroutine gcm_h1_2d_dump_physics_variables
! local

    if(write_history_h1)then

#ifndef DYAMOND
    diag_phys_accu_vars_h1_2d%uwind%f          = diag_phys_accu_vars_h1_2d%uwind%f        /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%vwind%f          = diag_phys_accu_vars_h1_2d%vwind%f        /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%phiFace%f        = diag_phys_accu_vars_h1_2d%phiFace%f      /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%temp%f           = diag_phys_accu_vars_h1_2d%temp%f         /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%mpressureFace%f  = diag_phys_accu_vars_h1_2d%mpressureFace%f/diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%omegaFull%f      = diag_phys_accu_vars_h1_2d%omegaFull%f    /diag_phys_accu_vars_h1_2d%ncount
#ifdef DYCORE_NHD
    diag_phys_accu_vars_h1_2d%wwwFace%f        = diag_phys_accu_vars_h1_2d%wwwFace%f      /diag_phys_accu_vars_h1_2d%ncount
#endif
    diag_phys_accu_vars_h1_2d%qv%f             = diag_phys_accu_vars_h1_2d%qv%f           /diag_phys_accu_vars_h1_2d%ncount
#endif

#ifdef AMIPW_PHYSICS
    diag_phys_accu_vars_h1_2d%rad_thten%f      = diag_phys_accu_vars_h1_2d%rad_thten%f /diag_phys_accu_vars_h1_2d%ncount
#ifndef DYAMOND
    diag_phys_accu_vars_h1_2d%qc%f             = diag_phys_accu_vars_h1_2d%qc%f        /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%qr%f             = diag_phys_accu_vars_h1_2d%qr%f        /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%qi%f             = diag_phys_accu_vars_h1_2d%qi%f        /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%qs%f             = diag_phys_accu_vars_h1_2d%qs%f        /diag_phys_accu_vars_h1_2d%ncount
    diag_phys_accu_vars_h1_2d%qg%f             = diag_phys_accu_vars_h1_2d%qg%f        /diag_phys_accu_vars_h1_2d%ncount
#endif
#endif
    end if
 
    return
  end subroutine gcm_h1_2d_dump_physics_variables

  subroutine gcm_h1_1d_rest_physics_variables
! reset to zero
    diag_phys_accu_vars_h1_1d%ncount  = 0

    if(write_history_h1)then

    diag_phys_accu_vars_h1_1d%precl%f = zero; diag_phys_accu_vars_h1_1d%precc%f = zero; diag_phys_accu_vars_h1_1d%prect%f = zero
    diag_phys_accu_vars_h1_1d%snowl%f = zero; diag_phys_accu_vars_h1_1d%grapl%f = zero
    diag_phys_accu_vars_h1_1d%ps%f    = zero
    diag_phys_accu_vars_h1_1d%ts%f    = zero
    diag_phys_accu_vars_h1_1d%shflx%f = zero
    diag_phys_accu_vars_h1_1d%qflx%f  = zero

    diag_phys_accu_vars_h1_1d%flwut%f = zero; diag_phys_accu_vars_h1_1d%flwdt%f = zero; diag_phys_accu_vars_h1_1d%flwus%f = zero; diag_phys_accu_vars_h1_1d%flwds%f = zero
    diag_phys_accu_vars_h1_1d%flwutc%f= zero; diag_phys_accu_vars_h1_1d%flwdtc%f= zero; diag_phys_accu_vars_h1_1d%flwusc%f= zero; diag_phys_accu_vars_h1_1d%flwdsc%f= zero
    diag_phys_accu_vars_h1_1d%lwcf%f  = zero

    diag_phys_accu_vars_h1_1d%fswut%f = zero; diag_phys_accu_vars_h1_1d%fswdt%f = zero; diag_phys_accu_vars_h1_1d%fswus%f = zero; diag_phys_accu_vars_h1_1d%fswds%f = zero
    diag_phys_accu_vars_h1_1d%fswutc%f= zero; diag_phys_accu_vars_h1_1d%fswdtc%f= zero; diag_phys_accu_vars_h1_1d%fswusc%f= zero; diag_phys_accu_vars_h1_1d%fswdsc%f= zero
    diag_phys_accu_vars_h1_1d%swcf%f  = zero

#ifdef AMIPW_PHYSICS
#ifndef DYAMOND
    diag_phys_accu_vars_h1_1d%asdir%f  = zero
    diag_phys_accu_vars_h1_1d%asdif%f  = zero
    diag_phys_accu_vars_h1_1d%aldir%f  = zero
    diag_phys_accu_vars_h1_1d%aldif%f  = zero
    diag_phys_accu_vars_h1_1d%emiss%f  = zero
#endif
    diag_phys_accu_vars_h1_1d%qsfc%f   = zero
#endif

    end if
 
    return
  end subroutine gcm_h1_1d_rest_physics_variables

  subroutine gcm_h1_2d_rest_physics_variables
! reset to zero
    diag_phys_accu_vars_h1_2d%ncount  = 0


    if(write_history_h1)then

#ifndef DYAMOND
    diag_phys_accu_vars_h1_2d%uwind%f     = zero
    diag_phys_accu_vars_h1_2d%vwind%f     = zero
    diag_phys_accu_vars_h1_2d%phiFace%f   = zero
    diag_phys_accu_vars_h1_2d%temp%f      = zero
    diag_phys_accu_vars_h1_2d%mpressureFace%f  = zero
    diag_phys_accu_vars_h1_2d%omegaFull%f = zero
#ifdef DYCORE_NHD
    diag_phys_accu_vars_h1_2d%wwwFace%f   = zero
#endif
    diag_phys_accu_vars_h1_2d%qv%f        = zero
#endif

#ifdef AMIPW_PHYSICS
    diag_phys_accu_vars_h1_2d%rad_thten%f = zero
#ifndef DYAMOND
    diag_phys_accu_vars_h1_2d%qc%f     = zero
    diag_phys_accu_vars_h1_2d%qr%f     = zero
    diag_phys_accu_vars_h1_2d%qi%f     = zero
    diag_phys_accu_vars_h1_2d%qs%f     = zero
    diag_phys_accu_vars_h1_2d%qg%f     = zero
#endif
#endif

    end if
 
    return
  end subroutine gcm_h1_2d_rest_physics_variables

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
