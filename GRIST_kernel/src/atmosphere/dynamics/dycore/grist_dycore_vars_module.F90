 
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Data Entry for GCM. This module contains almost all dycore vars 
!              for core dynamics; I will commnet some vars that have not been 
!              used till now (20181025, Yi Zhang)
!
! Revision history: 1. Feb, 2020, use new format to make code concise
!----------------------------------------------------------------------------

  module grist_dycore_vars_module

    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_1d_field, scalar_2d_field
    use grist_constants,    only: rearth, i4, r8, pi, zero
    use grist_nml_module,   only: nlev, nlevp

    implicit none

    public

!================================================
! Name Convention:
!   scalar_${varname}_at_$(location)_${time}
!   tend_${varname}_at_${location}_${where_from}
!   change to:
!   type%scalar_${varname}_${time}
!   type%tend_${varname}_${where_from}
!================================================

!================================================
! edge, full level
!================================================
   TYPE DYCORE_VAR_EDGE_FULL

   type(scalar_2d_field)   :: scalar_normal_velocity_n
   type(scalar_2d_field)   :: scalar_normal_mass_flux_n
   type(scalar_2d_field)   :: scalar_normal_pt_mass_flux_n
   type(scalar_2d_field)   :: scalar_U_wind_n 
   type(scalar_2d_field)   :: scalar_V_wind_n
   type(scalar_2d_field)   :: scalar_grad_delp_n
   type(scalar_2d_field)   :: tend_normal_velocity_vadv   ! vertical advection
   type(scalar_2d_field)   :: tend_normal_velocity_pgf    ! PGF 
   type(scalar_2d_field)   :: tend_normal_velocity_ke     ! kinetic energy
   type(scalar_2d_field)   :: tend_normal_velocity_nct    ! nonlinear coriolis term
   type(scalar_2d_field)   :: tend_hwind_laplacian_2nd
   type(scalar_2d_field)   :: tend_hwind_laplacian_4th
#ifdef GRIST_DYCORE_TEST
   type(scalar_2d_field)   :: tend_hwind_laplacian_6th
#endif

   END TYPE DYCORE_VAR_EDGE_FULL

!================================================
! primal cell, full level
!================================================
   TYPE DYCORE_VAR_CELL_FULL

   type(scalar_2d_field)   :: scalar_U_wind_n
   type(scalar_2d_field)   :: scalar_V_wind_n
   type(scalar_2d_field)   :: scalar_potential_temp_n     ! dry or modified pt
   type(scalar_2d_field)   :: scalar_potential_temp_iniupt! for diagnose
   type(scalar_2d_field)   :: scalar_mass_pt_n
   type(scalar_2d_field)   :: scalar_temp_n
   type(scalar_2d_field)   :: scalar_delhp_n
   type(scalar_2d_field)   :: scalar_delhp_np1
   type(scalar_2d_field)   :: scalar_divergence_n
   type(scalar_2d_field)   :: scalar_geopotential_n
   type(scalar_2d_field)   :: scalar_eta_mass_flux_n 
   type(scalar_2d_field)   :: scalar_hpressure_n
   type(scalar_2d_field)   :: scalar_mpressure_n     ! moist mass=pressure in hdc
   type(scalar_2d_field)   :: scalar_pressure_n      ! nh-diag-var
   type(scalar_2d_field)   :: scalar_pressure_rk     ! nh-diag-var
   type(scalar_2d_field)   :: scalar_pressure_np1    ! nh-diag-var
   type(scalar_2d_field)   :: scalar_alpha_np1       ! nh-diag-var
#ifdef DYCORE_NHD
   type(scalar_2d_field)   :: scalar_www_n           ! nh-diag-var
   type(scalar_2d_field)   :: scalar_delp_n          ! nh-temp-var
#endif
   type(scalar_2d_field)   :: scalar_delp_np1        ! nh-temp-var
   type(scalar_2d_field)   :: scalar_omega_n         ! dpim/dt
   type(scalar_2d_field)   :: scalar_omega_timavg    ! for physics 
   type(scalar_2d_field)   :: tend_mass_hori
   type(scalar_2d_field)   :: tend_mass_pt_hori
   type(scalar_2d_field)   :: tend_mass_pt_vert
   type(scalar_2d_field)   :: tend_mass_pt_laplacian_2nd
   type(scalar_2d_field)   :: tend_pt_n              ! total-adv tend across dycore steps, for cumulus

   END TYPE DYCORE_VAR_CELL_FULL

!================================================
! dual cell, full level
!================================================
   TYPE DYCORE_VAR_VERT_FULL

   type(scalar_2d_field)   :: scalar_abs_vor_n
   type(scalar_2d_field)   :: scalar_rel_vor_n
   type(scalar_2d_field)   :: scalar_pot_vor_n

   END TYPE DYCORE_VAR_VERT_FULL

!================================================
! edge, face level
!================================================
   TYPE DYCORE_VAR_EDGE_FACE

   type(scalar_2d_field)   :: scalar_normal_mass_flux_n

   END TYPE DYCORE_VAR_EDGE_FACE
!================================================
! primal cell, face level
!================================================
   TYPE DYCORE_VAR_CELL_FACE

   type(scalar_2d_field)   :: scalar_geopotential_n
   type(scalar_2d_field)   :: scalar_hpressure_n
   type(scalar_2d_field)   :: scalar_eta_mass_flux_n
   type(scalar_2d_field)   :: tend_mass_hori
   type(scalar_2d_field)   :: scalar_delhp_n
   type(scalar_2d_field)   :: scalar_delp_n
#ifdef DYCORE_NHD
   type(scalar_2d_field)   :: scalar_www_n   ! nh-prog-var
   type(scalar_2d_field)   :: scalar_phi_n   ! nh-prog-var
   type(scalar_2d_field)   :: tend_www_laplacian_2nd
   type(scalar_2d_field)   :: tend_www_laplacian_4th
#endif
   type(scalar_2d_field)   :: scalar_www_timavg   ! nh-diag-var for ourput
   type(scalar_2d_field)   :: scalar_pressure_n   ! nh-diag-var
   type(scalar_2d_field)   :: scalar_mpressure_n  ! moist mass=pressure in hdc

   END TYPE DYCORE_VAR_CELL_FACE


!================================================
! surface
!================================================
   TYPE DYCORE_VAR_SURFACE

   type(scalar_1d_field)   :: scalar_hpressure_n
   type(scalar_1d_field)   :: scalar_pressure_n      ! for output
   type(scalar_1d_field)   :: tend_hpressure_cnty
   type(scalar_1d_field)   :: scalar_geopotential_n

   END TYPE DYCORE_VAR_SURFACE
!================================================
! geometric info
!================================================
   TYPE DYCORE_VAR_GEOGRAPHY

   type(scalar_1d_field)   :: scalar_leng_at_edp
   type(scalar_1d_field)   :: scalar_leng_at_edt
   type(scalar_1d_field)   :: scalar_lon_at_edge
   type(scalar_1d_field)   :: scalar_lat_at_edge

   type(scalar_1d_field)   :: scalar_area_at_pc
   type(scalar_1d_field)   :: scalar_lon_at_pc
   type(scalar_1d_field)   :: scalar_lat_at_pc

   type(scalar_1d_field)   :: scalar_area_at_dc
   type(scalar_1d_field)   :: scalar_lon_at_dc
   type(scalar_1d_field)   :: scalar_lat_at_dc

   type(scalar_2d_field)   :: scalar_lon_bnds ! CDO convention for hori interp, lon is implicitly on cell
   type(scalar_2d_field)   :: scalar_lat_bnds !

   END TYPE DYCORE_VAR_GEOGRAPHY

   type(DYCORE_VAR_EDGE_FULL) :: dycoreVarEdgeFull
   type(DYCORE_VAR_EDGE_FACE) :: dycoreVarEdgeFace
   type(DYCORE_VAR_CELL_FULL) :: dycoreVarCellFull
   type(DYCORE_VAR_CELL_FACE) :: dycoreVarCellFace
   type(DYCORE_VAR_VERT_FULL) :: dycoreVarVertFull
   type(DYCORE_VAR_SURFACE)   :: dycoreVarSurface
   type(DYCORE_VAR_GEOGRAPHY) :: dycoreVarGeography

  CONTAINS

   subroutine grist_dycore_vars_construct(mesh)
! io
     type(global_domain), intent(in) :: mesh
!  local
     integer(i4) :: it, ie, iv
!================================================
! edge, full level
!================================================
      call dycore_fill_edgeVar_2d_0(mesh,nlev,dycoreVarEdgeFull%scalar_normal_velocity_n)
      call dycore_fill_edgeVar_2d_0(mesh,nlev,dycoreVarEdgeFull%scalar_normal_mass_flux_n) 
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n)
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%scalar_U_wind_n)    
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%scalar_V_wind_n)
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%scalar_grad_delp_n)    
      call dycore_fill_edgeVar_2d_0(mesh,nlev,dycoreVarEdgeFull%tend_normal_velocity_vadv)   
      call dycore_fill_edgeVar_2d_0(mesh,nlev,dycoreVarEdgeFull%tend_normal_velocity_pgf)           
      call dycore_fill_edgeVar_2d_0(mesh,nlev,dycoreVarEdgeFull%tend_normal_velocity_ke)     
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%tend_normal_velocity_nct)    
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%tend_hwind_laplacian_2nd)    
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%tend_hwind_laplacian_4th)
#ifdef GRIST_DYCORE_TEST
      call dycore_fill_edgeVar_2d(mesh,nlev,dycoreVarEdgeFull%tend_hwind_laplacian_6th)
#endif
!================================================
! primal cell, full level
!================================================
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_U_wind_n)     ! reconstucted U at cell
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_V_wind_n)     ! reconstucted V at cell
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_potential_temp_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_potential_temp_iniupt)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_mass_pt_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_temp_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_delhp_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_delhp_np1)
      call dycore_fill_primVar_2d_0(mesh,nlev,dycoreVarCellFull%scalar_divergence_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_geopotential_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_eta_mass_flux_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_hpressure_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_mpressure_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_pressure_n)   ! nh-diag-var
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_pressure_rk)  ! nh-diag-var
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_pressure_np1) ! nh-diag-var
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_alpha_np1)    ! nh-diag-var
#ifdef DYCORE_NHD
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_www_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_delp_n)       ! nh-temp-var
#endif
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_delp_np1)     ! nh-temp-var
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_omega_n)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%scalar_omega_timavg)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%tend_mass_hori)
      !call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%tend_mass_pt_hori)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%tend_mass_pt_vert)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%tend_mass_pt_laplacian_2nd)
      call dycore_fill_primVar_2d(mesh,nlev,dycoreVarCellFull%tend_pt_n)
!================================================
! dual cell, full level
!================================================
      call dycore_fill_dualVar_2d(mesh,nlev,dycoreVarVertFull%scalar_abs_vor_n)
      call dycore_fill_dualVar_2d_0(mesh,nlev,dycoreVarVertFull%scalar_rel_vor_n)
      call dycore_fill_dualVar_2d(mesh,nlev,dycoreVarVertFull%scalar_pot_vor_n)
!================================================
! edge, face level
!================================================
      !call dycore_fill_edgeVar_2d(mesh,nlevp, dycoreVarEdgeFace%scalar_normal_velocity_n)
#ifdef DYCORE_NHD
      call dycore_fill_edgeVar_2d(mesh,nlevp,dycoreVarEdgeFace%scalar_normal_mass_flux_n)
#endif
!================================================
! primal cell, face level
!================================================
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_geopotential_n)
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_hpressure_n)
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_eta_mass_flux_n)
#ifdef DYCORE_NHD
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%tend_mass_hori)
#endif
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_delhp_n)
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_delp_n)
#ifdef DYCORE_NHD
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_www_n)     ! nh-prog-var
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_phi_n)     ! nh-prog-var
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%tend_www_laplacian_2nd) ! diffusion tend 
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%tend_www_laplacian_4th) ! diffusion tend 
#endif
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_www_timavg)      ! nh-diag-var
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_pressure_n)      ! nh-diag-var
      call dycore_fill_primVar_2d(mesh,nlevp,dycoreVarCellFace%scalar_mpressure_n)     ! moist presssure 
!================================================
! surface
!================================================
    call dycore_fill_primVar_1d(mesh,dycoreVarSurface%scalar_hpressure_n)
    call dycore_fill_primVar_1d(mesh,dycoreVarSurface%scalar_geopotential_n)
    call dycore_fill_primVar_1d(mesh,dycoreVarSurface%scalar_pressure_n)
    call dycore_fill_primVar_1d(mesh,dycoreVarSurface%tend_hpressure_cnty)

!================================================
! geometric
!================================================
    call dycore_fill_edgeVar_1d(mesh,dycoreVarGeography%scalar_leng_at_edp)
    call dycore_fill_edgeVar_1d(mesh,dycoreVarGeography%scalar_leng_at_edt)
    call dycore_fill_edgeVar_1d(mesh,dycoreVarGeography%scalar_lon_at_edge)
    call dycore_fill_edgeVar_1d(mesh,dycoreVarGeography%scalar_lat_at_edge)

    call dycore_fill_primVar_1d(mesh,dycoreVarGeography%scalar_area_at_pc)
    call dycore_fill_primVar_1d(mesh,dycoreVarGeography%scalar_lon_at_pc)
    call dycore_fill_primVar_1d(mesh,dycoreVarGeography%scalar_lat_at_pc)

    call dycore_fill_dualVar_1d(mesh,dycoreVarGeography%scalar_area_at_dc)
    call dycore_fill_dualVar_1d(mesh,dycoreVarGeography%scalar_lon_at_dc)
    call dycore_fill_dualVar_1d(mesh,dycoreVarGeography%scalar_lat_at_dc)

    call dycore_fill_primVar_2d(mesh,mesh%maxvnb,dycoreVarGeography%scalar_lon_bnds)
    call dycore_fill_primVar_2d(mesh,mesh%maxvnb,dycoreVarGeography%scalar_lat_bnds)

!
! asign some constant values
! 
     do it = 1, mesh%nt
        dycoreVarGeography%scalar_area_at_dc%f(it)   = (rearth**2)*mesh%tri_areag(it)
        dycoreVarGeography%scalar_lon_at_dc%f(it)    = mesh%tri_c_lon(it)*180._r8/pi
        dycoreVarGeography%scalar_lat_at_dc%f(it)    = mesh%tri_c_lat(it)*180._r8/pi
     end do

     do ie = 1, mesh%ne
        dycoreVarGeography%scalar_leng_at_edp%f(ie)  = rearth*mesh%edp_leng(ie)
        dycoreVarGeography%scalar_leng_at_edt%f(ie)  = rearth*mesh%edt_leng(ie)
        dycoreVarGeography%scalar_lon_at_edge%f(ie)  = mesh%edt_c_lon(ie)*180._r8/pi
        dycoreVarGeography%scalar_lat_at_edge%f(ie)  = mesh%edt_c_lat(ie)*180._r8/pi
     end do


     do iv = 1, mesh%nv
        dycoreVarGeography%scalar_area_at_pc%f(iv)   = (rearth**2)*mesh%plg_areag(iv)
        dycoreVarGeography%scalar_lon_at_pc%f(iv)    = mesh%vtx_lon(iv)*180._r8/pi
        dycoreVarGeography%scalar_lat_at_pc%f(iv)    = mesh%vtx_lat(iv)*180._r8/pi

        do it = 1, mesh%vtx_nnb(iv)
           dycoreVarGeography%scalar_lon_bnds%f(it,iv) = mesh%tri_c_lon(mesh%vtx_tr(it,iv))*180._r8/pi
           dycoreVarGeography%scalar_lat_bnds%f(it,iv) = mesh%tri_c_lat(mesh%vtx_tr(it,iv))*180._r8/pi
        end do

        do it = mesh%vtx_nnb(iv)+1,mesh%maxvnb
           dycoreVarGeography%scalar_lon_bnds%f(it,iv) = mesh%tri_c_lon(mesh%vtx_tr(mesh%vtx_nnb(iv),iv))*180._r8/pi
           dycoreVarGeography%scalar_lat_bnds%f(it,iv) = mesh%tri_c_lat(mesh%vtx_tr(mesh%vtx_nnb(iv),iv))*180._r8/pi
        end do

     end do

     return

    end subroutine grist_dycore_vars_construct

    subroutine grist_dycore_vars_destruct()

!================================================
! edge, full level
!================================================
      call wrap_deallocate_2d(dycoreVarEdgeFull%scalar_normal_velocity_n)
      call wrap_deallocate_2d(dycoreVarEdgeFull%scalar_normal_mass_flux_n)
      call wrap_deallocate_2d(dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n)
      call wrap_deallocate_2d(dycoreVarEdgeFull%scalar_U_wind_n)   
      call wrap_deallocate_2d(dycoreVarEdgeFull%scalar_V_wind_n)
      call wrap_deallocate_2d(dycoreVarEdgeFull%scalar_grad_delp_n)
      call wrap_deallocate_2d(dycoreVarEdgeFull%tend_normal_velocity_vadv)
      call wrap_deallocate_2d(dycoreVarEdgeFull%tend_normal_velocity_pgf)
      call wrap_deallocate_2d(dycoreVarEdgeFull%tend_normal_velocity_ke)
      call wrap_deallocate_2d(dycoreVarEdgeFull%tend_normal_velocity_nct)
      call wrap_deallocate_2d(dycoreVarEdgeFull%tend_hwind_laplacian_2nd)
      call wrap_deallocate_2d(dycoreVarEdgeFull%tend_hwind_laplacian_4th)
#ifdef GRIST_DYCORE_TEST
      call wrap_deallocate_2d(dycoreVarEdgeFull%tend_hwind_laplacian_6th)
#endif
!================================================
! primal cell, full
!================================================
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_U_wind_n)      
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_V_wind_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_potential_temp_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_potential_temp_iniupt)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_mass_pt_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_temp_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_delhp_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_delhp_np1)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_divergence_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_geopotential_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_hpressure_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_eta_mass_flux_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_mpressure_n)      ! moist mass
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_pressure_n)       ! nh-diag-var
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_pressure_rk)      ! nh-diag-var
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_pressure_np1)     ! nh-diag-var
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_alpha_np1)        ! nh-diag-var

#ifdef DYCORE_NHD
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_www_n)
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_delp_n)       ! nh-temp-var
#endif

      call wrap_deallocate_2d(dycoreVarCellFull%scalar_delp_np1)     ! nh-temp-var
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_omega_n)       ! nh-diag-var
      call wrap_deallocate_2d(dycoreVarCellFull%scalar_omega_timavg)  ! nh-diag-var
      call wrap_deallocate_2d(dycoreVarCellFull%tend_mass_hori)
      !call wrap_deallocate_2d(dycoreVarCellFull%tend_mass_pt_hori)
      call wrap_deallocate_2d(dycoreVarCellFull%tend_mass_pt_vert)
      call wrap_deallocate_2d(dycoreVarCellFull%tend_mass_pt_laplacian_2nd)
      call wrap_deallocate_2d(dycoreVarCellFull%tend_pt_n)

!================================================
! dual cell, full
!================================================

      call wrap_deallocate_2d(dycoreVarVertFull%scalar_abs_vor_n)
      call wrap_deallocate_2d(dycoreVarVertFull%scalar_rel_vor_n)
      call wrap_deallocate_2d(dycoreVarVertFull%scalar_pot_vor_n)

!================================================
! edge, face level
!================================================

#ifdef DYCORE_NHD
      call wrap_deallocate_2d(dycoreVarEdgeFace%scalar_normal_mass_flux_n)
#endif

!================================================
! primal cell, face level 
!================================================
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_geopotential_n)
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_hpressure_n)
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_eta_mass_flux_n)
#ifdef DYCORE_NHD
      call wrap_deallocate_2d(dycoreVarCellFace%tend_mass_hori)
#endif
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_delhp_n)
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_delp_n)
#ifdef DYCORE_NHD
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_www_n)   ! nh-prog-var
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_phi_n)   ! nh-prog-var
      call wrap_deallocate_2d(dycoreVarCellFace%tend_www_laplacian_2nd)
      call wrap_deallocate_2d(dycoreVarCellFace%tend_www_laplacian_4th)
#endif
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_www_timavg)! nh-prog-var
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_pressure_n)   ! nh-diag-var
      call wrap_deallocate_2d(dycoreVarCellFace%scalar_mpressure_n)  ! nh-diag-var
!================================================
! surface
!================================================
      call wrap_deallocate_1d(dycoreVarSurface%scalar_hpressure_n)
      call wrap_deallocate_1d(dycoreVarSurface%scalar_pressure_n)
      call wrap_deallocate_1d(dycoreVarSurface%tend_hpressure_cnty)
      call wrap_deallocate_1d(dycoreVarSurface%scalar_geopotential_n)
!================================================
! geometric 
!================================================
      call wrap_deallocate_1d(dycoreVarGeography%scalar_area_at_pc)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_area_at_dc)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_leng_at_edp)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_leng_at_edt)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_lon_at_pc)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_lat_at_pc)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_lon_at_dc)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_lat_at_dc)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_lon_at_edge)
      call wrap_deallocate_1d(dycoreVarGeography%scalar_lat_at_edge)

      call wrap_deallocate_2d(dycoreVarGeography%scalar_lon_bnds)
      call wrap_deallocate_2d(dycoreVarGeography%scalar_lat_bnds)

      return
    end subroutine grist_dycore_vars_destruct

!================================================
! BELOW are private
!================================================

    subroutine dycore_fill_edgeVar_2d(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4)        ,   intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%ne))
        var%f    = zero
        var%pos  = 6
        return
    end subroutine dycore_fill_edgeVar_2d

    subroutine dycore_fill_edgeVar_2d_0(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4)        ,   intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%ne))
        var%f    = zero
        var%pos  = 6
        return
    end subroutine dycore_fill_edgeVar_2d_0

    subroutine dycore_fill_primVar_2d(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4)        ,   intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%nv))
        var%f    = zero
        var%pos  = 0
        return
    end subroutine dycore_fill_primVar_2d

    subroutine dycore_fill_primVar_2d_0(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4)        ,   intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%nv))
        var%f    = zero
        var%pos  = 0
        return
    end subroutine dycore_fill_primVar_2d_0

    subroutine dycore_fill_dualVar_2d(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4)        ,   intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%nt))
        var%f    = zero
        var%pos  = 1
        return
    end subroutine dycore_fill_dualVar_2d

    subroutine dycore_fill_dualVar_2d_0(mesh,nLevel,var)
       type(global_domain),   intent(in)    :: mesh
       integer(i4)        ,   intent(in)    :: nLevel
       type(scalar_2d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(nLevel,mesh%nt))
        var%f    = zero
        var%pos  = 1
        return
    end subroutine dycore_fill_dualVar_2d_0

    subroutine dycore_fill_edgeVar_1d(mesh,var)
       type(global_domain),   intent(in)    :: mesh
       type(scalar_1d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(mesh%ne))
        var%f    = zero
        var%pos  = 6
        return
    end subroutine dycore_fill_edgeVar_1d

    subroutine dycore_fill_primVar_1d(mesh,var)
       type(global_domain),   intent(in)    :: mesh
       type(scalar_1d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(mesh%nv))
        var%f    = zero
        var%pos  = 0
        return
    end subroutine dycore_fill_primVar_1d

    subroutine dycore_fill_dualVar_1d(mesh,var)
       type(global_domain),   intent(in)    :: mesh
       type(scalar_1d_field), intent(inout) :: var

        if(.not.allocated(var%f)) allocate(var%f(mesh%nt))
        var%f    = zero
        var%pos  = 1
        return
    end subroutine dycore_fill_dualVar_1d

    subroutine wrap_deallocate_2d(var)
       type(scalar_2d_field), intent(inout)  :: var
       if(allocated(var%f)) deallocate(var%f)
       return
    end subroutine wrap_deallocate_2d

    subroutine wrap_deallocate_1d(var)
       type(scalar_1d_field), intent(inout)  :: var
       if(allocated(var%f)) deallocate(var%f)
       return
    end subroutine wrap_deallocate_1d

   end module grist_dycore_vars_module
