!======================================================
!  Yi Zhang, evaluate wrf's physics diag
!======================================================

 module grist_PhysW_physics_diag

    use grist_constants,                only: r8,i4, zero
    use grist_nml_module,               only: nlev, nlevp
    use grist_physics_data_structure,   only: pstate
    use grist_PhysW_data_structure,       only: pstate_wrf, psurf_wrf
    use grist_vi_cldfrac_diagnostics,   only: cloud_diagnostics_calc
 
    implicit none
    private
    save

    public  :: grist_PhysW_diag

contains

  subroutine grist_PhysW_diag(dtime,ncell)
      real(r8),     intent(in) :: dtime
      integer(i4),  intent(in) :: ncell

!================================================
! precipitation
!================================================
       pstate%scalar_precc_surface%f(1:ncell) = psurf_wrf%raincv (1:ncell,1)/dtime/1000._r8    ! from mm/step to m/s (step is fastPhys step, model_step)
! Seperate account for rain and solid
       pstate%scalar_precl_surface%f(1:ncell) = psurf_wrf%rainncv(1:ncell,1)   /dtime/1000._r8 ! from mm/step to m/s
       pstate%scalar_snowl_surface%f(1:ncell) = psurf_wrf%snowncv(1:ncell,1)   /dtime/1000._r8 ! from mm/step to m/s
       pstate%scalar_grapl_surface%f(1:ncell) = psurf_wrf%graupelncv(1:ncell,1)/dtime/1000._r8 ! from mm/step to m/s
       pstate%scalar_prect_surface%f(1:ncell) = pstate%scalar_precl_surface%f(1:ncell) + pstate%scalar_precc_surface%f(1:ncell)

!================================================
! cloud
!================================================
       call cloud_diagnostics_calc(ncell, nlev, nlevp, &
                                   pstate_wrf%cldfra(1:ncell,nlev:1:-1,1),pstate_wrf%p8w (1:ncell,nlevp:1:-1,1),pstate_wrf%p_phy (1:ncell,nlev:1:-1,1),&
                                   pstate_wrf%cldtot(1:ncell,1), pstate_wrf%cldlow(1:ncell,1), pstate_wrf%cldmed(1:ncell,1),pstate_wrf%cldhgh(1:ncell,1))

!================================================
! surface
! fill some pstate vars for model output, and some needed for restart
!================================================

      pstate%scalar_rainc_surface%f(1:ncell)         = psurf_wrf%rainc (1:ncell,1)! restart
      pstate%scalar_rainnc_surface%f(1:ncell)        = psurf_wrf%rainnc(1:ncell,1)! restart
      pstate%scalar_snownc_surface%f(1:ncell)        = psurf_wrf%snownc(1:ncell,1)! restart
      pstate%scalar_grapnc_surface%f(1:ncell)        = psurf_wrf%graupelnc(1:ncell,1)! restart
      pstate%scalar_raincv_surface%f(1:ncell)        = psurf_wrf%raincv(1:ncell,1)! restart, as CU is SLOW physics

      pstate%ts_at_pc_surface%f(1:ncell)             = psurf_wrf%tsk   (1:ncell,1)! restart, and h1 history
      pstate%atm_in_shflx_at_pc_surface%f(1:ncell)   = psurf_wrf%hfx   (1:ncell,1)! restart
      pstate%qfx_at_pc_surface%f (1:ncell)           = psurf_wrf%qfx   (1:ncell,1)! restart
      pstate%atm_in_qflx_at_pc_surface%f(1,1:ncell)  = psurf_wrf%qfx   (1:ncell,1)! output
! rad flux
      pstate%atm_out_flwds_at_pc_surface%f(1:ncell)  = psurf_wrf%glw   (1:ncell,1)! restart, down lw
      pstate%atm_out_netsw_at_pc_surface%f(1:ncell)  = psurf_wrf%gsw   (1:ncell,1)! restart, net sw
      pstate%atm_out_fswds_at_pc_surface%f(1:ncell)  = psurf_wrf%swdown(1:ncell,1)! restart, down sw
! surface layer+pbl
      pstate%ustar_at_pc_surface%f(1:ncell)          = psurf_wrf%ust(1:ncell,1)   ! restart
      pstate%znt_at_pc_surface%f(1:ncell)            = psurf_wrf%znt(1:ncell,1)   ! restart
      pstate%mol_at_pc_surface%f(1:ncell)            = psurf_wrf%mol(1:ncell,1)   ! restart
      pstate%pblh_at_pc_surface%f(1:ncell)           = psurf_wrf%pblh(1:ncell,1)  ! restart
      pstate%mavail_at_pc_surface%f(1:ncell)         = psurf_wrf%mavail(1:ncell,1)! restart
! albedo
      pstate%snowhland_at_pc_surface%f(1:ncell)      = psurf_wrf%snow  (1:ncell,1)! rst
      pstate%atm_in_asdir_at_pc_surface%f(1:ncell)   = psurf_wrf%asdir (1:ncell,1)! rst
      pstate%atm_in_asdif_at_pc_surface%f(1:ncell)   = psurf_wrf%asdif (1:ncell,1)! rst
      pstate%atm_in_aldir_at_pc_surface%f(1:ncell)   = psurf_wrf%aldir (1:ncell,1)! rst
      pstate%atm_in_aldif_at_pc_surface%f(1:ncell)   = psurf_wrf%aldif (1:ncell,1)! rst
      
      call physics_diag_clean

    return
  end subroutine grist_PhysW_diag

  subroutine physics_diag_clean
    return
  end subroutine physics_diag_clean

 end module grist_PhysW_physics_diag
