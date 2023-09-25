
!----------------------------------------------------------------------------
! Created on 2019
! Author: Yi Zhang
! Version 1.0
! Description: This driver is the top driver of the WRF-based model physics
!              used as physpkg2.
! Revision history:
!--------------------------------------------------------------------------

 module grist_PhysW_pkg

   use grist_constants,             only: r8, i4, one, zero, rad2deg
   use grist_PhysW_physics_update,  only: PhysW_update
   use grist_PhysW_data_structure,  only: psurf_wrf
   use grist_PhysW_nml_module,      only: PhysW_sf_scheme, PhysW_mp_scheme, PhysW_bl_scheme, PhysW_cu_scheme
#ifdef USE_NOAHMP
   use grist_datam_initial_data_module,  only: initialData_xice_at_pc_surface, & ! use as in LSM initial
                                               initialData_snowh_at_pc_surface   ! m
   use grist_datam_static_data_module,   only: staticData_sic_at_pc_surface 
#endif
   implicit none
   private
   
   public ::  grist_PhysW_pkg_init,   &
              grist_PhysW_pkg_run_bc, &
              grist_PhysW_pkg_run_ac, &
              grist_PhysW_pkg_run_bc1,&
              grist_PhysW_pkg_run_ac1,&
              grist_PhysW_pkg_final

   CONTAINS
 
   subroutine grist_PhysW_pkg_run_bc(ncell,nLevel,nspecies,itimestep,dtime,dxmean)
    use grist_PhysW_cumulus_driver,      only: grist_PhysW_cumulus_run
    use grist_PhysW_microphysics_driver, only: grist_PhysW_microphysics_run
    use grist_PhysW_radiation_driver,    only: grist_PhysW_radiation_run
    use grist_PhysW_surface_driver,      only: grist_PhysW_surface_run
    use grist_PhysW_pbl_driver,          only: grist_PhysW_pbl_run
    use grist_PhysW_physics_diag,        only: grist_PhysW_diag
! io
    integer(i4), intent(in)    :: ncell
    integer(i4), intent(in)    :: nLevel
    integer(i4), intent(in)    :: nspecies
    integer(i4), intent(in)    :: itimestep
    real(r8),    intent(in)    :: dtime
    real(r8),    intent(in)    :: dxmean
! local
    logical :: first_called

      first_called = .true.

      call grist_PhysW_microphysics_run(ncell,nLevel,nspecies,itimestep,dtime)
      call PhysW_update(dtime, ncell, 'MP',trim(PhysW_mp_scheme),first_called)

      call grist_PhysW_surface_run     (ncell,nLevel,nspecies,itimestep,dtime,dxmean)
! for lsm
      psurf_wrf%raint(1:ncell,1) = psurf_wrf%raincv(1:ncell,1) + psurf_wrf%rainncv(1:ncell,1) ! mm/step
      psurf_wrf%snow (1:ncell,1) = psurf_wrf%snow (1:ncell,1)  + psurf_wrf%snowncv(1:ncell,1) ! mm
#ifndef REG_AMP21
      where(psurf_wrf%xland(1:ncell,1).ge.1.5_r8) psurf_wrf%snow (1:ncell,1) = zero
#endif

     return
   end subroutine grist_PhysW_pkg_run_bc

   subroutine grist_PhysW_pkg_run_ac(ncell,nLevel,nspecies,itimestep,dtime,dxmean,coszrs)
    use grist_PhysW_cumulus_driver,      only: grist_PhysW_cumulus_run
    use grist_PhysW_microphysics_driver, only: grist_PhysW_microphysics_run
    use grist_PhysW_radiation_driver,    only: grist_PhysW_radiation_run
    use grist_PhysW_surface_driver,      only: grist_PhysW_surface_run
    use grist_PhysW_pbl_driver,          only: grist_PhysW_pbl_run
    use grist_PhysW_physics_diag,        only: grist_PhysW_diag
! io
    integer(i4), intent(in)    :: ncell
    integer(i4), intent(in)    :: nLevel
    integer(i4), intent(in)    :: nspecies
    integer(i4), intent(in)    :: itimestep
    real(r8),    intent(in)    :: dtime
    real(r8),    intent(in)    :: dxmean
    real(r8),    intent(in)    :: coszrs(ncell)
! local
    logical :: first_called

      first_called = .true.

      call grist_PhysW_pbl_run         (ncell,nLevel,nspecies,itimestep,dtime)
      call PhysW_update(dtime, ncell, 'PBL',trim(PhysW_bl_scheme))
!
! slow physics is called here, but not update model state
! CU->RAD
!
      call grist_PhysW_cumulus_run     (ncell,nLevel,nspecies,itimestep,dtime)
      call PhysW_update(dtime, ncell, 'CU',trim(PhysW_cu_scheme))

      call grist_PhysW_radiation_run   (ncell,nLevel,nspecies,itimestep,dtime,coszrs)
      call PhysW_update(dtime, ncell, 'RAD','RRTMG')

      call grist_PhysW_diag(dtime,ncell)

    return
   end subroutine grist_PhysW_pkg_run_ac
!
! bc1 and ac1 mimic cam's time-split and sequence, for testing
!
   subroutine grist_PhysW_pkg_run_bc1(ncell,nLevel,nspecies,itimestep,dtime,dxmean,coszrs)
    use grist_PhysW_cumulus_driver,      only: grist_PhysW_cumulus_run
    use grist_PhysW_microphysics_driver, only: grist_PhysW_microphysics_run
    use grist_PhysW_radiation_driver,    only: grist_PhysW_radiation_run
    use grist_PhysW_surface_driver,      only: grist_PhysW_surface_run
    use grist_PhysW_pbl_driver,          only: grist_PhysW_pbl_run
    use grist_PhysW_physics_diag,        only: grist_PhysW_diag
! io
    integer(i4), intent(in)    :: ncell
    integer(i4), intent(in)    :: nLevel
    integer(i4), intent(in)    :: nspecies
    integer(i4), intent(in)    :: itimestep
    real(r8),    intent(in)    :: dtime
    real(r8),    intent(in)    :: dxmean
    real(r8),    intent(in)    :: coszrs(ncell)
! local
    logical :: first_called

      first_called = .true.

      call grist_PhysW_cumulus_run     (ncell,nLevel,nspecies,itimestep,dtime)
      call PhysW_update(dtime, ncell, 'CU','TDK',first_called)

      call grist_PhysW_microphysics_run(ncell,nLevel,nspecies,itimestep,dtime)
      call PhysW_update(dtime, ncell, 'MP','LIN')

      call grist_PhysW_radiation_run   (ncell,nLevel,nspecies,itimestep,dtime,coszrs)
      call PhysW_update(dtime, ncell, 'RAD','RRTMG')

      call grist_PhysW_surface_run     (ncell,nLevel,nspecies,itimestep,dtime,dxmean)
! for lsm
      psurf_wrf%raint(1:ncell,1) = psurf_wrf%raincv(1:ncell,1) + psurf_wrf%rainncv(1:ncell,1)

     return
   end subroutine grist_PhysW_pkg_run_bc1

   subroutine grist_PhysW_pkg_run_ac1(ncell,nLevel,nspecies,itimestep,dtime,dxmean)
    use grist_PhysW_cumulus_driver,      only: grist_PhysW_cumulus_run
    use grist_PhysW_microphysics_driver, only: grist_PhysW_microphysics_run
    use grist_PhysW_radiation_driver,    only: grist_PhysW_radiation_run
    use grist_PhysW_surface_driver,      only: grist_PhysW_surface_run
    use grist_PhysW_pbl_driver,          only: grist_PhysW_pbl_run
    use grist_PhysW_physics_diag,        only: grist_PhysW_diag
! io
    integer(i4), intent(in)    :: ncell
    integer(i4), intent(in)    :: nLevel
    integer(i4), intent(in)    :: nspecies
    integer(i4), intent(in)    :: itimestep
    real(r8),    intent(in)    :: dtime
    real(r8),    intent(in)    :: dxmean
! local
    logical :: first_called

      call grist_PhysW_pbl_run         (ncell,nLevel,nspecies,itimestep,dtime)
      call PhysW_update(dtime, ncell, 'PBL','YSU')

      call grist_PhysW_diag(dtime,ncell)

    return
   end subroutine grist_PhysW_pkg_run_ac1

   subroutine grist_PhysW_pkg_init(ncell,ncell_full, nLevel,nspecies,dtime,lats,lons)
    use grist_physics_data_structure,  only: get_tracer_information
    use grist_PhysW_data_structure,      only: grist_PhysW_data_structure_construct, pstate_wrf
    use grist_PhysW_cumulus_driver,      only: grist_PhysW_cumulus_init
    use grist_PhysW_microphysics_driver, only: grist_PhysW_microphysics_init
    use grist_PhysW_radiation_driver,    only: grist_PhysW_radiation_init
    use grist_PhysW_surface_driver,      only: grist_PhysW_surface_init
    use grist_PhysW_pbl_driver,          only: grist_PhysW_pbl_init
#ifdef USE_NOAHMP
    use grist_lsm_noahmp_init,         only: lsm_noahmp_init
    use grist_nml_module,              only: levsoil, test_real_case
    use grist_lsm_noahmp_resVars,      only: grist_lsm_resVars_construct
#endif
    use grist_physics_data_structure,  only: pstate
    use grist_nml_module,              only: physpkg, sub_physpkg
!
! io
!
      integer(i4),   intent(in)   ::  ncell, ncell_full
      integer(i4),   intent(in)   ::  nLevel
      integer(i4),   intent(in)   ::  nspecies
      real(r8),      intent(in)   ::  dtime, lats(ncell), lons(ncell) ! radian
! local
      real(r8)    :: zzlnd, zzwtr, thinld, alblnd, xmava
      integer(i4) :: iv

!
! do initilization of all WRF-based modules
!
      call grist_PhysW_data_structure_construct(ncell,nLevel,nspecies)


      pstate_wrf%xlat(1:ncell,1)  = lats(1:ncell)*rad2deg
      pstate_wrf%xlong(1:ncell,1) = lons(1:ncell)*rad2deg

 !     call get_tracer_information(physpkg, sub_physpkg)

      call grist_PhysW_cumulus_init(ncell, nLevel)
      call grist_PhysW_microphysics_init(ncell,nLevel,nspecies)
      call grist_PhysW_radiation_init(ncell, nLevel, dtime)
      call grist_PhysW_surface_init
      call grist_PhysW_pbl_init(ncell, nLevel)
#ifdef USE_NOAHMP
      if(test_real_case)then 
           call lsm_noahmp_init(ncell, nLevel, levsoil, lats, lons)
           call grist_lsm_resVars_construct(ncell_full) ! nv_full
      end if
#endif
!
! initial some state
!
      zzlnd  = 0.1_r8
      zzwtr  = 0.0001_r8
      thinld = 0.04_r8
      alblnd = 0.2_r8
      xmava  = 0.3_r8
    
      do iv = 1, ncell
         psurf_wrf%xland(iv,1) = 2._r8-pstate%landfrac_at_pc_surface%f(iv)
         if(psurf_wrf%xland(iv,1).eq.0)then
            print*,"bad ",pstate%landfrac_at_pc_surface%f(iv),2._r8-pstate%landfrac_at_pc_surface%f(iv)
         end if
      end do
#ifdef USE_NOAHMP
      if(test_real_case)then
#ifdef REG_AMP21
         psurf_wrf%xice(1:ncell,1)  = initialData_xice_at_pc_surface%f(1:ncell)
#else
         psurf_wrf%xice(1:ncell,1)  = staticData_sic_at_pc_surface%f(1:ncell)
#endif
         psurf_wrf%snow(1:ncell,1)  = initialData_snowh_at_pc_surface%f(1:ncell)*1e3_r8 ! m->mm, lnd init file must have unit in meter
      end if
#endif
#ifndef REG_AMP21
! only valid for land, yizhang202109, mind for regression
      where(psurf_wrf%xland(1:ncell,1).ge.1.5_r8) psurf_wrf%snow (1:ncell,1) = zero
#endif

      psurf_wrf%gsw(    1:ncell,1)  = zero
      psurf_wrf%glw(    1:ncell,1)  = zero
      psurf_wrf%swdown( 1:ncell,1)  = zero
      psurf_wrf%ust(    1:ncell,1)  = zero
      psurf_wrf%mol(    1:ncell,1)  = zero
      psurf_wrf%pblh(   1:ncell,1)  = zero
      psurf_wrf%hfx(    1:ncell,1)  = zero
      psurf_wrf%qfx(    1:ncell,1)  = zero
      psurf_wrf%lh(    1:ncell,1)   = zero
      psurf_wrf%rainc(1:ncell,1)    = zero ! accumulateed precc
      psurf_wrf%rainnc(1:ncell,1)   = zero ! accumulated  rain 
      psurf_wrf%snownc(1:ncell,1)   = zero ! accumulated  
      psurf_wrf%graupelnc(1:ncell,1)= zero ! accumulated  
      psurf_wrf%raincv(1:ncell,1)   = zero
      psurf_wrf%rainncv(1:ncell,1)  = zero
      if(trim(PhysW_sf_scheme).eq.'SFCLAYV381'.or.trim(PhysW_sf_scheme).eq.'SFCLAYV341') &
      psurf_wrf%ust(    1:ncell,1)  = 0.0001_r8

      do iv = 1, ncell
         if(psurf_wrf%xland(iv,1) .lt. 1.5)then
           !albbck(iv,1) = alblnd
           psurf_wrf%albedo(iv,1) = alblnd
           psurf_wrf%emiss(iv,1)  = 0.85_r8
           psurf_wrf%thc(iv,1)    = thinld
           psurf_wrf%znt(iv,1)    = zzlnd
           !z0(iv,1)     = zzlnd
           psurf_wrf%mavail(iv,1) = xmava
           psurf_wrf%asdir(iv,1)  = psurf_wrf%albedo(iv,1)
           psurf_wrf%asdif(iv,1)  = psurf_wrf%albedo(iv,1)
           psurf_wrf%aldir(iv,1)  = psurf_wrf%albedo(iv,1)
           psurf_wrf%aldif(iv,1)  = psurf_wrf%albedo(iv,1)
         else
           !albbck(iv,1) = 0.08_r8
           psurf_wrf%albedo(iv,1) = 0.08_r8
           psurf_wrf%emiss(iv,1)  = 0.98_r8
           psurf_wrf%thc(iv,1)    = thinld
           psurf_wrf%znt(iv,1)    = zzwtr
           !z0(iv,1)     = zzwtr
           psurf_wrf%mavail(iv,1) = 1.0_r8
           psurf_wrf%asdir(iv,1) = psurf_wrf%albedo(iv,1)
           psurf_wrf%asdif(iv,1) = psurf_wrf%albedo(iv,1)
           psurf_wrf%aldir(iv,1) = psurf_wrf%albedo(iv,1)
           psurf_wrf%aldif(iv,1) = psurf_wrf%albedo(iv,1)
         endif
      enddo

      return
   end subroutine grist_PhysW_pkg_init

   subroutine grist_PhysW_pkg_final

    use grist_PhysW_data_structure,      only: grist_PhysW_data_structure_destruct
    use grist_PhysW_cumulus_driver,      only: grist_PhysW_cumulus_final
    use grist_PhysW_microphysics_driver, only: grist_PhysW_microphysics_final
    use grist_PhysW_radiation_driver,    only: grist_PhysW_radiation_final
    use grist_PhysW_surface_driver,      only: grist_PhysW_surface_final
    use grist_PhysW_pbl_driver,          only: grist_PhysW_pbl_final
#ifdef USE_NOAHMP
    use grist_lsm_noahmp_vars,         only: deallocate_lsm_noahmp_1d, deallocate_lsm_noahmp_2d
    use grist_lsm_noahmp_resVars,      only: grist_lsm_resVars_destruct
#endif

      call grist_PhysW_cumulus_final
      call grist_PhysW_microphysics_final
      call grist_PhysW_radiation_final
      call grist_PhysW_surface_final
      call grist_PhysW_pbl_final
      call grist_PhysW_data_structure_destruct
#ifdef USE_NOAHMP
      call deallocate_lsm_noahmp_2d
      call deallocate_lsm_noahmp_1d
      call grist_lsm_resVars_destruct
#endif 
      return
   end subroutine grist_PhysW_pkg_final

 end module grist_PhysW_pkg
