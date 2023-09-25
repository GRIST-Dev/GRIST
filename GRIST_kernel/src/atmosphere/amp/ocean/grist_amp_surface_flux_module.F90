
!----------------------------------------------------------------------------
! Created on 2020
! Author: Yi Zhang
! Version 1.0
! Description:  This module is used as a column model; scheme from cesm
! Revision history:
!----------------------------------------------------------------------------

 module grist_amp_surface_flux_module

   use grist_constants,                  only: r8, i4, pi, zero, one, gravity
   use grist_amp_shr_flux_mod,           only: shr_flux_atmOcn
   use grist_nml_module,                 only: start_ymd, start_tod
   use grist_time_manager,               only: get_curr_calday
   use grist_zenith,                     only: zenith
   use grist_mpi

   implicit none
   private
  
   public :: grist_amp_surface_flux_run

   contains

   subroutine grist_amp_surface_flux_run(ncell, nlev, itimestep, dtime, lat, lon, call_radiation)
     use grist_physics_data_structure, only: pstate
! iostream
     integer(i4),  intent(in)   :: ncell, nlev
     integer(i4),  intent(in)   :: itimestep
     real(r8),     intent(in)   :: dtime           ! time step
     logical, optional, intent(in)   :: call_radiation
     real(r8),     intent(in)   :: lat(:)
     real(r8),     intent(in)   :: lon(:)
! local
     integer(i4):: scalar_mask_at_pc_surface(ncell),  icell
     real(r8)  ::  scalar_zbot_at_pc_surface(ncell)
     real(r8)  ::  scalar_ubot_at_pc_surface(ncell)
     real(r8)  ::  scalar_vbot_at_pc_surface(ncell)
     real(r8)  ::  scalar_thbot_at_pc_surface(ncell)
     real(r8)  ::  scalar_qbot_at_pc_surface(ncell)
     real(r8)  ::  scalar_rbot_at_pc_surface(ncell)
     real(r8)  ::  scalar_tbot_at_pc_surface(ncell)
     real(r8)  ::  scalar_us_at_pc_surface(ncell)
     real(r8)  ::  scalar_vs_at_pc_surface(ncell)
     real(r8)  ::  scalar_ts_at_pc_surface(ncell)

!     real(r8)  ::  scalar_shflx_at_pc_surface(ncell)
     real(r8)  ::  scalar_lhflx_at_pc_surface(ncell)
!     real(r8)  ::  scalar_qqflx_at_pc_surface(ncell)
!     real(r8)  ::  scalar_lwup_at_pc_surface(ncell)
!     real(r8)  ::  scalar_taux_at_pc_surface(ncell)
!     real(r8)  ::  scalar_tauy_at_pc_surface(ncell)

     real(r8)  ::  scalar_tref_at_pc_surface(ncell)   ! diag:  2m ref height T (K)
     real(r8)  ::  scalar_qref_at_pc_surface(ncell)   ! diag:  2m ref humidity (kg/kg)
     real(r8)  ::  scalar_duu10n_at_pc_surface(ncell) ! diag: 10m wind speed squared (m/s)^2

     real(r8)  ::  calday                         ! current calendar day
     real(r8)  ::  coszrs(ncell)                  ! Cosine solar zenith angle

       scalar_mask_at_pc_surface(1:ncell) = int(pstate%ocnfrac_at_pc_surface%f(1:ncell),i4)
!
! air state at the lowest level
!
#ifdef AMIPC_PHYSICS
       scalar_zbot_at_pc_surface(1:ncell) = pstate%z_at_pc_full_level%f(nlev,1:ncell)
       scalar_ubot_at_pc_surface(1:ncell) = pstate%u_wind_at_pc_full_level%f(nlev,1:ncell)
       scalar_vbot_at_pc_surface(1:ncell) = pstate%v_wind_at_pc_full_level%f(nlev,1:ncell)
       scalar_thbot_at_pc_surface(1:ncell)= pstate%temp_at_pc_full_level%f(nlev,1:ncell)*pstate%exner_at_pc_full_level%f(nlev,1:ncell)
       scalar_qbot_at_pc_surface(1:ncell) = pstate%tracer_mxrt_at_pc_full_level%f(1,nlev,1:ncell)
       scalar_rbot_at_pc_surface(1:ncell) = pstate%delp_at_pc_full_level%f(nlev,1:ncell)/gravity/&
                                            (pstate%z_at_pc_face_level%f(nlev,1:ncell)-pstate%z_at_pc_face_level%f(nlev+1,1:ncell))
       scalar_tbot_at_pc_surface(1:ncell) = pstate%temp_at_pc_full_level%f(nlev,1:ncell)
#endif
!
! surface state
!
       scalar_us_at_pc_surface(1:ncell)   = zero ! assume sea surface does not move as in CAM5
       scalar_vs_at_pc_surface(1:ncell)   = zero
       where(scalar_mask_at_pc_surface/=0) pstate%ts_at_pc_surface%f(1:ncell) = pstate%sst_at_pc_surface%f(1:ncell)
       scalar_ts_at_pc_surface(1:ncell)   = pstate%ts_at_pc_surface%f(1:ncell)
!
! run
!
       call shr_flux_atmOcn(nmax=ncell,  & ! in 
                            zbot = scalar_zbot_at_pc_surface , & ! in
                            ubot = scalar_ubot_at_pc_surface , & ! in
                            vbot = scalar_vbot_at_pc_surface , & ! in
                            thbot= scalar_thbot_at_pc_surface, & ! in
                            qbot = scalar_qbot_at_pc_surface , & ! in
                            rbot = scalar_rbot_at_pc_surface , & ! in
                            tbot = scalar_tbot_at_pc_surface , & ! in
                            us   = scalar_us_at_pc_surface   , & ! in
                            vs   = scalar_vs_at_pc_surface   , & ! in
                            ts   = scalar_ts_at_pc_surface   , & ! in
                            mask = scalar_mask_at_pc_surface , & ! in
                            sen  = pstate%atm_in_shflx_at_pc_surface%f(1:ncell), & ! out
                            lat  = scalar_lhflx_at_pc_surface(1:ncell)         , & ! out
                            lwup = pstate%atm_in_lwup_at_pc_surface%f(1:ncell) , & ! out
                            evap = pstate%atm_in_qflx_at_pc_surface%f(1,1:ncell),& ! out
                            taux = pstate%atm_in_taux_at_pc_surface%f(1:ncell) , & ! out
                            tauy = pstate%atm_in_tauy_at_pc_surface%f(1:ncell) , & ! out
                            tref = scalar_tref_at_pc_surface , & ! out
                            qref = scalar_qref_at_pc_surface , & ! out
                            duu10n=scalar_duu10n_at_pc_surface, & ! out
                            ustar_sv=pstate%ustar_at_pc_surface%f(1:ncell))

       call get_curr_calday(start_ymd, start_tod, itimestep, dtime, calday)
       call zenith(calday, coszrs, ncell, lon, lat )

       do icell=1, ncell
         if (scalar_mask_at_pc_surface(icell) /=0 ) then
              !if(call_radiation)then
                 pstate%atm_in_lwup_at_pc_surface%f(icell)= -one*pstate%atm_in_lwup_at_pc_surface%f(icell)
              !else
              !   pstate%atm_in_lwup_at_pc_surface%f(1:ncell)= zero
              !end if
              pstate%atm_in_shflx_at_pc_surface%f(icell)  = -one*pstate%atm_in_shflx_at_pc_surface%f(icell)
              pstate%atm_in_qflx_at_pc_surface%f(1,icell) = -one*pstate%atm_in_qflx_at_pc_surface%f(1,icell)
              pstate%atm_in_qflx_at_pc_surface%f(2:,icell)= zero
              pstate%atm_in_taux_at_pc_surface%f(icell)   = -one*pstate%atm_in_taux_at_pc_surface%f(icell)
              pstate%atm_in_tauy_at_pc_surface%f(icell )  = -one*pstate%atm_in_tauy_at_pc_surface%f(icell)
       
              ! Only for Ocean, snow and seaice have not been considered! Lixh
#ifdef OCNABD
! use this compiling option currently for regression old Aqua-planet solutions
              if(coszrs(icell) .ge. 0._r8)then
! LiXH implementation
                pstate%atm_in_asdir_at_pc_surface%f(icell) = 0.037_r8/(1.1_r8*coszrs(icell)**1.4_r8 + 0.15_r8)
                pstate%atm_in_asdif_at_pc_surface%f(icell) = pstate%atm_in_asdir_at_pc_surface%f(icell)
                pstate%atm_in_aldir_at_pc_surface%f(icell) = 0.06_r8
                pstate%atm_in_aldif_at_pc_surface%f(icell) = 0.06_r8
! CAM3 version
#ifdef CAM3OCNABD
                pstate%atm_in_aldir_at_pc_surface%f(icell)  = (.026_r8/(coszrs(icell)**1.7_r8 + .065_r8)) + &
                                                              (.15_r8*(coszrs(icell) - 0.10_r8)*(coszrs(icell) - 0.50_r8)*(coszrs(icell) - 1._r8))
                pstate%atm_in_asdir_at_pc_surface%f(icell)  = pstate%atm_in_aldir_at_pc_surface%f(icell)  
                pstate%atm_in_aldif_at_pc_surface%f(icell)  = 0.06_r8
                pstate%atm_in_asdif_at_pc_surface%f(icell)  = 0.06_r8
#endif
              else
                pstate%atm_in_asdir_at_pc_surface%f(icell) = 0._r8
                pstate%atm_in_asdif_at_pc_surface%f(icell) = 0._r8
                pstate%atm_in_aldir_at_pc_surface%f(icell) = 0._r8
                pstate%atm_in_aldif_at_pc_surface%f(icell) = 0._r8
              end if
#endif
          end if
       end do

      return
   end subroutine grist_amp_surface_flux_run
 
 end module grist_amp_surface_flux_module
