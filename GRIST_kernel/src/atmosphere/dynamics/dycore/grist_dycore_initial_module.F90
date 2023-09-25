
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Initialize dycore with various idealized initials
!    (1) 3D RH wave, tested
!    (2) 3D Colliding MODONs, tested
!    (3) Jablonowski and Williamson Steady State (JWSS) test, tested
!    (4) Jablonowski and Williamson Perturbed (JWBW) test, tested
!    (5) Held and Suarez long term dry forcing, tested
!    (6) Circular mountain in KSP15, tested
!    (7) SCHAR mountain in KSP15, tested
!    (8) DCMIP3-1, tested
!    (9) DCMIP2-0, tested
!    (10) MIRW, tested
! Revision history:
!----------------------------------------------------------------------------

 module grist_dycore_initial_module

  use grist_constants_dbl,    only: gravity, pi, i4, r8, rearth,rdry,cp, omega,p00
  use grist_domain_types,     only: global_domain
  use grist_data_types,       only: scalar_1d_field
  use grist_math_module,      only: arcdistll, convert_vector_sph2cart
  use grist_nml_module,       only: nlev, nh_dynamics, modon_radius
  use grist_hpe_constants,    only: p0, eta_full_a, eta_full_b, eta_face_a, eta_face_b 
  use grist_hpe_hydro_pgf,    only: calc_hpe_hpressure_face_level
#ifndef SEQ_GRIST
  use grist_util_module  ,    only: exchange_init_data
  use grist_config_partition, only: debug_data_2d
#endif
  use grist_mpi,              only: mpi_rank
! data
  use grist_dycore_vars_module, only: dycoreVarCellFull, &
                                      dycoreVarSurface,  &
                                      dycoreVarEdgeFull, &
                                      dycoreVarCellFace

!  use grist_datam_static_data_module,  only: staticData_phis_at_pc_surface
                                       

  implicit none

  private

  public  :: grist_dycore_initial

  CONTAINS

  subroutine grist_dycore_initial(mesh, testcase, nlev, nlevp)
!
! io
!
   type(global_domain), intent(inout) :: mesh
   character*(*)      , intent(in)    :: testcase
   integer(i4)        , intent(in)    :: nlev, nlevp
!
! local
!
   real(r8)                            :: vector_velocity(3)
   real(r8)                            :: scalar_u
   real(r8)                            :: scalar_v
   real(r8)                            :: eta
   integer(i4)                         :: ilev
   integer(i4)                         :: it, ie, iv
! for RH3D
   real(r8), allocatable               :: stream_function_dual_cell(:)
   integer(i4)                         :: v0,v1
   real(r8)                            :: v0v1(3),flag
   real(r8), dimension(nlev)           :: scalar_template_1d_nlev_a
   real(r8), dimension(nlev)           :: scalar_template_1d_nlev_b
   real(r8), dimension(nlev)           :: scalar_template_1d_nlev_c
   real(r8), dimension(nlevp)          :: scalar_template_1d_nlevp_a
   real(r8), dimension(nlevp)          :: scalar_template_1d_nlevp_b
   real(r8)                            :: scalar_template_a
   real(r8), dimension(nlev,mesh%nv_full) ::  scalar_qqq_at_pc_full_level_n

   select case(trim(testcase))

   case('JWSS','JWBW','JWBW2')

      dycoreVarSurface%scalar_hpressure_n%f = 1e5_r8  ! Pa

      DO ilev = 1, nlev
 
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
  
         do ie = 1, mesh%ne

            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          !
          ! for normal velocity
          !
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do

         do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_temp_n%f(ilev,iv) = ini_temperature(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), &
                                                                        eta, eta_full_a(ilev), eta_full_b(ilev),&
                                                                        1._r8, 0.2_r8, 'JWSS')
            dycoreVarSurface%scalar_geopotential_n%f(iv) = ini_phis(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'JWSS')
         end do

      END DO
#ifdef DYCORE_NHD
      if(nh_dynamics)then
         do ilev = 1, nlev+1
            eta  = eta_face_a(ilev)+eta_face_b(ilev)
            do iv = 1, mesh%nv
               dycoreVarCellFace%scalar_phi_n%f(ilev,iv) = ini_geop(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8),eta,1._r8,0.2_r8,eta_face_a(ilev),eta_face_b(ilev),'JWSS')
               dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv) = dycoreVarCellFace%scalar_phi_n%f(ilev,iv)
            end do
         end do
         dycoreVarCellFace%scalar_www_n%f   = 0._r8
      end if
#endif

   case('RH3D')
! define stream function at dual cell
      allocate(stream_function_dual_cell(mesh%nt))
      do it = 1, mesh%nt
         stream_function_dual_cell(it)  = ini_stream_function(real(mesh%tri(it)%c%lon,r8),real(mesh%tri(it)%c%lat,r8),'RH3D')
      end do

      DO ilev = 1, nlev
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
! normal velocity
         do ie = 1, mesh%ne
            v0   = mesh%edp(ie)%v(1) ! tri point
            v1   = mesh%edp(ie)%v(2) ! tri point
            v0v1 = mesh%tri(v1)%c%p-mesh%tri(v0)%c%p
            flag = sign(1._r8, dot_product(v0v1,real(mesh%edp(ie)%tg,r8)))
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = &
            flag*(stream_function_dual_cell(v1)-stream_function_dual_cell(v0))/(rearth*mesh%edp(ie)%leng) 
         end do
! scalar temperature and phis
         do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_temp_n%f(ilev,iv) = ini_temperature(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), &
                                                                        eta, eta_full_a(ilev), eta_full_b(ilev), &
                                                                        1._r8, 0.2_r8, 'RH3D')
            dycoreVarSurface%scalar_geopotential_n%f(iv) =             ini_phis(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'RH3D')
            dycoreVarSurface%scalar_hpressure_n%f(iv)    = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'RH3D')
         end do
      END DO

      deallocate(stream_function_dual_cell)

   case('MIRW')

      do ilev = 1, nlev
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
         do ie = 1, mesh%ne

            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          ! for normal velocity
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do
      end do

      dycoreVarCellFull%scalar_temp_n%f   = 288._r8 

      do iv = 1, mesh%nv
         dycoreVarSurface%scalar_geopotential_n%f(iv) = ini_phis(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'MIRW')
         dycoreVarSurface%scalar_hpressure_n%f(iv)    = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'MIRW',&
                                                                          real(dycoreVarSurface%scalar_geopotential_n%f(iv),r8))
      end do
#ifdef MIRWTOPO
   case('MIRWTOPO')

      do ilev = 1, nlev
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
         do ie = 1, mesh%ne

            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          ! for normal velocity
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do
      end do

      dycoreVarCellFull%scalar_temp_n%f   = 288._r8 

      do iv = 1, mesh%nv
         dycoreVarSurface%scalar_geopotential_n%f(iv) = staticData_phis_at_pc_surface%f(iv)
         dycoreVarSurface%scalar_hpressure_n%f(iv)    = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'MIRWTOPO',&
                                                                          real(dycoreVarSurface%scalar_geopotential_n%f(iv),r8))
      end do
#endif
   case('SCHAR') ! only for nh-dynamics

!      if(.not.allocated(scalar_template_1d_nlevp_a%f)) allocate(scalar_template_1d_nlevp_a%f(nlev+1))

      DO ilev = 1, nlev
 
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
         do ie = 1, mesh%ne
            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)
            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          !
          ! for normal velocity
          !
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do

         do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_temp_n%f(ilev,iv) = ini_temperature(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), &
                                                                        eta, eta_full_a(ilev), eta_full_b(ilev),&
                                                                        1._r8, 0.2_r8, 'SCHAR')
            dycoreVarSurface%scalar_geopotential_n%f(iv) = ini_phis(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'SCHAR')
            dycoreVarSurface%scalar_hpressure_n%f(iv)    = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'SCHAR')
         end do
      END DO
!
! construct face level pressure, note this will be done again in rkfb_init,
! but effects are the same
!
       do iv = 1, mesh%nv
         !call calc_hpe_hpressure_face_level(dycoreVarSurface%scalar_hpressure_n%f(iv), & ! in
         !                                   scalar_template_1d_nlevp_a)               ! out
         do ilev = 1, nlevp
            dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do
         !dycoreVarCellFace%scalar_hpressure_n%f(:,iv)  = scalar_template_1d_nlevp_a
       end do

       do iv = 1, mesh%nv
          do ilev = 1, nlev
             dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv) = (rdry*300._r8)*&
                                                              log(1.e5_r8/dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv))-&
                                                              200._r8*(sin(mesh%vtx_lat(iv))**2)
#ifdef DYCORE_NHD
             dycoreVarCellFace%scalar_phi_n%f(ilev,iv)       = dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)
#endif
          end do
      end do
#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_www_n%f    = 0._r8
      dycoreVarCellFace%scalar_phi_n%f(nlev+1,:) = dycoreVarSurface%scalar_geopotential_n%f
! for init p-energy
      dycoreVarCellFull%scalar_geopotential_n%f(:,:) = (dycoreVarCellFace%scalar_phi_n%f(1:nlev,:)+&
                                                        dycoreVarCellFace%scalar_phi_n%f(2:nlev+1,:))*0.5_r8
#endif

!      deallocate(scalar_template_1d_nlevp_a%f)

   case('DCMIP2-0')

!      if(.not.allocated(scalar_template_1d_nlevp_a%f)) allocate(scalar_template_1d_nlevp_a%f(nlev+1))

      DO ilev = 1, nlev
 
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
  
         do ie = 1, mesh%ne

            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          !
          ! for normal velocity
          !
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = 0._r8 !dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do

         do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_temp_n%f(ilev,iv) = ini_temperature(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), &
                                                                        eta, eta_full_a(ilev), eta_full_b(ilev),&
                                                                        1._r8, 0.2_r8, 'DCMIP2-0',999._r8, &
                                                                        eta_face_a(ilev), eta_face_b(ilev), eta_face_a(ilev+1), eta_face_b(ilev+1))
            dycoreVarSurface%scalar_geopotential_n%f(iv) = ini_phis(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'DCMIP2-0')
            dycoreVarSurface%scalar_hpressure_n%f(iv)    = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'DCMIP2-0')
         end do
      END DO
#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_www_n%f           = 0._r8
#endif

! h-pressure at face
       do iv = 1, mesh%nv
         !call calc_hpe_hpressure_face_level(dycoreVarSurface%scalar_hpressure_n%f(iv), & ! in
         !                                   scalar_template_1d_nlevp_a)               ! out
         do ilev = 1, nlevp
            dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do
         !dycoreVarCellFace%scalar_hpressure_n%f(:,iv)  = scalar_template_1d_nlevp_a
       end do

       do iv = 1, mesh%nv
          do ilev = 1, nlev
             dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv) = gravity*(300._r8/0.0065_r8)*(1._r8-&
                                                                 (dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv)/1.e5_r8)**(rdry*0.0065/gravity))
#ifdef DYCORE_NHD
             dycoreVarCellFace%scalar_phi_n%f(ilev,iv)       = dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)
#endif
          end do
      end do
#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_phi_n%f(nlev+1,:) = dycoreVarSurface%scalar_geopotential_n%f
! for init p-energy
      dycoreVarCellFull%scalar_geopotential_n%f(:,:) = (dycoreVarCellFace%scalar_phi_n%f(1:nlev,:)+&
                                                       dycoreVarCellFace%scalar_phi_n%f(2:nlev+1,:))*0.5_r8

#endif
!      deallocate(scalar_template_1d_nlevp_a%f)

   case('DCMIP2-1','DCMIP210') ! only for nh-dynamics

!      if(.not.allocated(scalar_template_1d_nlevp_a%f)) allocate(scalar_template_1d_nlevp_a%f(nlev+1))

      DO ilev = 1, nlev
 
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
  
         do ie = 1, mesh%ne

            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          !
          ! for normal velocity
          !
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do

         do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_temp_n%f(ilev,iv) = ini_temperature(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), &
                                                                        eta, eta_full_a(ilev), eta_full_b(ilev),&
                                                                        1._r8, 0.2_r8, 'DCMIP2-1')
            dycoreVarSurface%scalar_geopotential_n%f(iv) = ini_phis(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),trim(testcase))
            dycoreVarSurface%scalar_hpressure_n%f(iv)    = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'DCMIP2-1')
         end do
      END DO
!
! construct face level pressure, note this will be done again in rkfb_init, 
! but effects are the same
!
       do iv = 1, mesh%nv
         !call calc_hpe_hpressure_face_level(dycoreVarSurface%scalar_hpressure_n%f(iv), & ! in
         !                                   scalar_template_1d_nlevp_a)               ! out
         do ilev = 1, nlevp
            dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do
         !dycoreVarCellFace%scalar_hpressure_n%f(:,iv)  = scalar_template_1d_nlevp_a
       end do

       do iv = 1, mesh%nv
          do ilev = 1, nlev
             dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv) = (rdry*300._r8)*&
                                                              log(1e5_r8/dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv))-&
                                                              20._r8*20._r8*(sin(mesh%vtx_lat(iv))**2)*0.5_r8
#ifdef DYCORE_NHD
             dycoreVarCellFace%scalar_phi_n%f(ilev,iv)       = dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)
#endif
          end do
      end do
#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_www_n%f    = 0._r8
      dycoreVarCellFace%scalar_phi_n%f(nlev+1,:) = dycoreVarSurface%scalar_geopotential_n%f
! for init p-energy
      dycoreVarCellFull%scalar_geopotential_n%f(:,:) = (dycoreVarCellFace%scalar_phi_n%f(1:nlev,:)+&
                                                       dycoreVarCellFace%scalar_phi_n%f(2:nlev+1,:))*0.5_r8
#endif

!      deallocate(scalar_template_1d_nlevp_a%f)

   case('DCMIP3-1') ! only for nh-dynamics

!      if(.not.allocated(scalar_template_1d_nlevp_a%f)) allocate(scalar_template_1d_nlevp_a%f(nlev+1))

      DO ilev = 1, nlev
 
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
         do ie = 1, mesh%ne

            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          !
          ! for normal velocity
          !
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do

      END DO
!
! construct face level pressure, note this will be done again in rkfb_init,
! but effects are the same
!
      do iv = 1, mesh%nv
         dycoreVarSurface%scalar_geopotential_n%f(iv) = ini_phis(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'DCMIP3-1')
         
         dycoreVarSurface%scalar_hpressure_n%f(iv)    = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),'DCMIP3-1')
         
         !call calc_hpe_hpressure_face_level(dycoreVarSurface%scalar_hpressure_n%f(iv), & ! in
         !                                   scalar_template_1d_nlevp_a)               ! out
         do ilev = 1, nlevp
            dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do
         !dycoreVarCellFace%scalar_hpressure_n%f(:,iv)  = scalar_template_1d_nlevp_a
      
      end do

#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_www_n%f               = 0._r8
#endif
  
      do ilev = 1, nlev   ! face level
          eta  = eta_face_a(ilev)+eta_face_b(ilev)
          do iv = 1, mesh%nv
            dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv) = ini_geop(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8),eta,1._r8,0.2_r8,eta_face_a(ilev),eta_face_b(ilev),'DCMIP3-1')
#ifdef DYCORE_NHD
            dycoreVarCellFace%scalar_phi_n%f(ilev,iv)          = dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)
#endif
          end do
      end do
      dycoreVarCellFace%scalar_geopotential_n%f(nlev+1,:)     = dycoreVarSurface%scalar_geopotential_n%f

#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_phi_n%f(nlev+1,:)              = dycoreVarSurface%scalar_geopotential_n%f
      dycoreVarCellFull%scalar_geopotential_n%f(:,:) = (dycoreVarCellFace%scalar_phi_n%f(1:nlev,:)+&
                                                       dycoreVarCellFace%scalar_phi_n%f(2:nlev+1,:))*0.5_r8
#endif

! for temperature
      do ilev = 1, nlev
        eta  = eta_full_a(ilev)+eta_full_b(ilev)
        do iv = 1, mesh%nv
          dycoreVarCellFull%scalar_temp_n%f(ilev,iv) = ini_temperature(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), &
                                                                      eta, eta_full_a(ilev), eta_full_b(ilev),&
                                                                      1._r8, 0.2_r8, 'DCMIP3-1', &
                                                                      dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv)/gravity)
          dycoreVarCellFull%scalar_potential_temp_iniupt%f(ilev,iv) = ini_temperature(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), &
                                                                      eta, eta_full_a(ilev), eta_full_b(ilev),&
                                                                      1._r8, 0.2_r8, 'DCMIP31-diag', &
                                                                      dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv)/gravity)
        end do
      end do

!      deallocate(scalar_template_1d_nlevp_a%f)

   case('MODON')

      dycoreVarSurface%scalar_hpressure_n%f = 1e5_r8  ! Pa

      DO ilev = 1, nlev
 
         eta  = eta_full_a(ilev)+eta_full_b(ilev)
  
         do ie = 1, mesh%ne

            call ini_vector_velocity(real(mesh%edt(ie)%c%p,r8)   , &
                                     real(mesh%edt(ie)%c%lon,r8) , &
                                     real(mesh%edt(ie)%c%lat,r8) , &
                                     eta, trim(testcase), &
                                     vector_velocity    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie)  =  scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie)  =  scalar_v
          !
          ! for normal velocity
          !
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, mesh%edp(ie)%nr)

         end do

         do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_temp_n%f(ilev,iv) = 300._r8
            dycoreVarSurface%scalar_geopotential_n%f(iv) = 0._r8
         end do

      END DO

!      if(nh_dynamics)then
!         do ilev = 1, nlev+1
!            do iv = 1, mesh%nv
!               dycoreVarCellFace%scalar_phi_n%f(ilev,iv) = (nlev-ilev+1)*2000._r8
!               dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv) = dycoreVarCellFace%scalar_phi_n%f(ilev,iv)
!            end do
!         end do
!         dycoreVarCellFace%scalar_www_n%f   = 0._r8
!      end if
   case('DCMIP1-1','DCMIP1-2','DCMIP1-3','DCMIP1-1-H','DCMIP11H')
      dycoreVarSurface%scalar_geopotential_n%f = 0._r8
      dycoreVarCellFull%scalar_temp_n%f      = 0._r8
      if (mpi_rank() .eq. 0) print*, "this is a isolated 3D advection test, no dyn state is initialized in dycore"

   case('real-WRFDA','real-GFS','real-ERAIP','real-ERAIM',&
        'DCMIP2016-BW','DCMIP2016-SC-simple', 'DCMIP2016-SC',&
        'DCMIP2016-TC','DCMIP2016-SC-A','DCMIP2016-SC-B','DCMIP2016-SC1','DCMIP2016-SCXX')
      if (mpi_rank() .eq. 0) print*, "this is a dtp test, no dyn state is initialized in dycore"

   case default
      print*,"you must select a test case,  dycore initial"
      stop
   end select

   return
  end subroutine grist_dycore_initial

!================================================
! convert UV wind component to vector velocity 
!================================================

  subroutine ini_vector_velocity(point_cart,      &
                                 lon, lat, eta, testcase, &
                                 vector_velocity, &
                                 scalar_u,        &
                                 scalar_v)
!
! io
!
    real(r8),        intent(in)    :: point_cart(3)       ! in 3D cartisian coordinate
    real(r8),        intent(in)    :: lon
    real(r8),        intent(in)    :: lat
    real(r8),        intent(in)    :: eta
    character*(*),   intent(in)    :: testcase
    real(r8),        intent(inout) :: vector_velocity(3)
    real(r8),        intent(inout) :: scalar_u
    real(r8),        intent(inout) :: scalar_v     

      scalar_u = ini_uwind(lon,lat,eta,testcase)
      scalar_v = ini_vwind(lon,lat,eta,testcase)
      call convert_vector_sph2cart(scalar_u,scalar_v,point_cart, vector_velocity)

    return
  end subroutine ini_vector_velocity

!================================================
!                 Functions
!================================================

  real(r8) function ini_temperature(lon, lat, eta, eta_a,eta_b, eta_s, eta_t, testcase, height, &
                                    eta_face_up_a,eta_face_up_b,eta_face_do_a,eta_face_do_b)
!================================================
!  initial temperature at a point
!   Lat in [-pi/2,pi/2], Lon in [-pi,pi]
!================================================
! io
    real(r8),       intent(in)   :: lon
    real(r8),       intent(in)   :: lat
    real(r8),       intent(in)   :: eta       ! at any level
    real(r8),       intent(in)   :: eta_a     ! eta's a component  
    real(r8),       intent(in)   :: eta_b     ! eta's b component
    real(r8),       intent(in)   :: eta_s     ! at surface
    real(r8),       intent(in)   :: eta_t     ! at top
    character*(*),  intent(in)   :: testcase
    real(r8),optional, intent(in) :: height
    real(r8),optional, intent(in) :: eta_face_up_a,eta_face_up_b,eta_face_do_a,eta_face_do_b
! local
    real(r8)                     :: eta_0, eta_v
    real(r8)                     :: temp_eta          ! only function of eta
    real(r8)                     :: temp0
    real(r8)                     :: u0
    real(r8)                     :: dtemp
    real(r8)                     :: pres
    real(r8)                     :: temp_lapse_rate
    real(r8)                     :: part1 
    real(r8)                     :: part2
    real(r8)                     :: part3
    real(r8)                     :: surface_pressure, surface_temp, bag_temp, ptb_temp, GG, dd, rr
    real(r8)                     :: lonc
    real(r8)                     :: pref ! for RH3D
    real(r8)                     :: pres_up, pres_do ! for RH3D

    select case(trim(testcase))

    case('JWSS','JWBW','JWBW2')
        eta_0           = 0.252
        eta_v           = (eta-eta_0)*pi*0.5_r8
        u0              = 35._r8
        temp0           = 288._r8
        dtemp           = 4.8e5_r8
        temp_lapse_rate = 0.005_r8

        if(eta.ge.eta_t.and.eta.le.eta_s)then
           temp_eta        = temp0*eta**(rdry*temp_lapse_rate/gravity)
        else  if( eta.lt.eta_t) then
           temp_eta        = temp0*eta**(rdry*temp_lapse_rate/gravity)+&
                             dtemp*(eta_t-eta)**5
        end if

        part1 = 0.75_r8*(eta*pi*u0/rdry)*sin(eta_v)*(cos(eta_v)**0.5_r8)

        part2 = ((-2*(sin(lat)**6)*(cos(lat)**2+1._r8/3._r8))+10._r8/63._r8)*2._r8*u0*(cos(eta_v)**1.5)

        part3 = (8._r8/5._r8*(cos(lat)**3)*((sin(lat)**2)+2._r8/3._r8)-pi/4._r8)*rearth*omega

        ini_temperature = temp_eta+part1*(part2+part3)

    case ('RH3D')
         temp0            = 288._r8
         temp_lapse_rate  = 0.0065_r8
         pref             = 95500._r8
         surface_pressure = ini_surface_pressure(lon,lat,testcase)
         pres             = eta_a*1e5_r8+eta_b*surface_pressure                 ! eq68
         ini_temperature  = temp0*(pres/pref)**(temp_lapse_rate*rdry/gravity)   ! eq66
    case('SCHAR','DCMIP2-1')
         ini_temperature  = 300._r8
    case('MIRW','MIRWTOPO')
         ini_temperature  = 288._r8
    case('DCMIP3-1')
         lonc             = 2._r8*pi/3._r8
         dd               = 5000._r8
         rr               = rearth*acos(cos(lat)*cos(lon-lonc))
         GG               = (gravity**2)/(1.e-4_r8*cp)
         surface_temp     = GG+(300._r8-GG)*exp(-0.25_r8*1.e-4_r8*20._r8*20._r8/(gravity**2)*(cos(2._r8*lat)-1))
         surface_pressure = ini_surface_pressure(lon,lat,'DCMIP3-1')
         pres             = eta_a*1.e5_r8+eta_b*surface_pressure
         part1            = (pres/surface_pressure)**(rdry/cp)
         bag_temp         = surface_temp*part1/((surface_temp*(part1-1._r8)/GG)+1._r8)
         ptb_temp         = (dd*dd)/(dd*dd+rr*rr)*sin(2._r8*pi*height/20000._r8)*((pres/1.e5_r8)**(rdry/cp))
         ini_temperature  = bag_temp+ptb_temp
    case('DCMIP31-diag')
         lonc             = 2._r8*pi/3._r8
         dd               = 5000._r8
         rr               = rearth*acos(cos(lat)*cos(lon-lonc))
         GG               = (gravity**2)/(1.e-4_r8*cp)
         surface_temp     = GG+(300._r8-GG)*exp(-0.25_r8*1.e-4_r8*20._r8*20._r8/(gravity**2)*(cos(2._r8*lat)-1))
         surface_pressure = ini_surface_pressure(lon,lat,'DCMIP3-1')
         pres             = eta_a*1.e5_r8+eta_b*surface_pressure
         part1            = (1.e5_r8/surface_pressure)**(rdry/cp)
         part2            =    (pres/surface_pressure)**(rdry/cp)
         bag_temp         = surface_temp*part1/((surface_temp*(part2-1._r8)/GG)+1._r8) ! already pt
         ini_temperature  = bag_temp
    case('DCMIP2-0')
         surface_pressure = ini_surface_pressure(lon,lat,'DCMIP2-0')
         pres             = eta_a*1.e5_r8+eta_b*surface_pressure
         ini_temperature  = 300._r8*(pres/1.e5_r8)**(rdry*0.0065_r8/gravity)
! FV init
       !
       !  pres_up         = eta_face_up_a*1.e5_r8+eta_face_up_b*surface_pressure
       !  pres_do         = eta_face_do_a*1.e5_r8+eta_face_do_b*surface_pressure
       !  temp0           = rdry*0.0065_r8/gravity+1._r8
       !  ini_temperature = gravity*300._r8*p00*((pres_up/p00)**temp0-(pres_do/p00)**temp0)/((rdry*0.0065_r8+gravity)*(pres_up-pres_do))
       !

    case default
       print*, "you must select a test case, please check!, ini_temperature"
       stop
    end select

    return
  end function ini_temperature

!================================================
! surface geopotential at a point
!================================================

  real(r8) function ini_phis(lon, lat, testcase)

! io
    real(r8),       intent(in)   :: lon
    real(r8),       intent(in)   :: lat
    character*(*),  intent(in)   :: testcase
! local
    real(r8)                     :: eta_0
    real(r8)                     :: eta_v
    real(r8)                     :: u0, h0, d0, kesi0, r0, r_m, rr_m
    real(r8)                     :: lonc, latc
    real(r8)                     :: part1
    real(r8)                     :: part2
    real(r8)                     :: part3

    select case(trim(testcase))
    case('JWSS','JWBW','JWBW2')
        eta_0            = 0.252
        eta_v            = (1._r8-eta_0)*pi*0.5_r8
        u0               = 35._r8
        part1            = u0*(cos(eta_v)**1.5)
        part2            = ((-2*(sin(lat)**6)*(cos(lat)**2+1._r8/3._r8))+10._r8/63._r8)*u0*(cos(eta_v)**1.5)
        part3            = (8._r8/5._r8*(cos(lat)**3)*((sin(lat)**2)+2._r8/3._r8)-pi/4._r8)*rearth*omega
        ini_phis = part1*(part2+part3)
    case('RH3D','DCMIP3-1')
        ini_phis = 0._r8
    case('MIRW')
        h0       = 2000._r8   ! meter
        d0       = 1500000._r8
        latc     = pi/6._r8
        lonc     = pi/2._r8
        r0       = rearth*acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc))
        ini_phis = gravity*h0*exp(-(r0/d0)**2)
    case('SCHAR')
        lonc     = 0._r8
        h0       = 250._r8    ! meter
        d0       = 5000._r8
        kesi0    = 4000._r8
        r0       = rearth*(lon-lonc) ! as in K15
        ini_phis = gravity*h0*exp(-(r0**2)/(d0**2))*(cos(pi*r0/kesi0)**2)*cos(lat)
    case('DCMIP2-1')
        lonc     = 0._r8 
        h0       = 250._r8    ! meter
        d0       = 5000._r8
        kesi0    = 4000._r8
        r0       = rearth*acos(cos(lat)*cos(lon-lonc))
        ini_phis = gravity*h0*exp(-(r0**2)/(d0**2))*(cos(pi*r0/kesi0)**2)
    case('DCMIP210')
        lonc     = 0._r8 
        h0       = 3000._r8    ! meter
        d0       = 2000._r8
        kesi0    = 4000._r8
        r0       = rearth*acos(cos(lat)*cos(lon-lonc))
        ini_phis = gravity*h0*exp(-(r0**2)/(d0**2)) !*(cos(pi*r0/kesi0)**2)
    case('DCMIP2-0')
        lonc     = 3._r8*pi/2._r8
        latc     = 0._r8
        r_m      = acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc))
        rr_m     = 0.75_r8*pi
        if(r_m.lt.rr_m)then
           ini_phis = gravity*1000._r8*(1._r8+cos(pi*r_m/rr_m))*(cos(pi*r_m/(pi/16._r8))**2)
           !ini_phis = gravity*2000._r8*(1._r8+cos(pi*r_m/rr_m))*(cos(pi*r_m/(pi/16._r8))**2)
        else
           ini_phis = 0._r8
        end if
    case default
       print*, "you must select a test case, please check!, ini_phis"
       stop
    end select

    return
  end function ini_phis

!================================================
! surface geopotential at a point
!================================================

  real(r8) function ini_geop(lon, lat, eta, eta_s, eta_t, eta_a, eta_b, testcase)

! io
    real(r8),       intent(in)   :: lon
    real(r8),       intent(in)   :: lat
    real(r8),       intent(in)   :: eta       ! at any level
    real(r8),       intent(in)   :: eta_s     ! at surface
    real(r8),       intent(in)   :: eta_t     ! at top
    real(r8),       intent(in)   :: eta_a
    real(r8),       intent(in)   :: eta_b
    character*(*),  intent(in)   :: testcase
! local
    real(r8)                     :: eta_0
    real(r8)                     :: eta_v
    real(r8)                     :: u0, h0, d0, kesi0, r0
    real(r8)                     :: lonc, latc
    real(r8)                     :: part0, tmp, dtemp, temp0, temp_lapse_rate
    real(r8)                     :: part1
    real(r8)                     :: part2
    real(r8)                     :: part3
    real(r8)                     :: GG, pres,surface_temp, surface_pressure

    select case(trim(testcase))
! not necessary needed for hydro case
    case('JWSS','JWBW','JWBW2')
        eta_0            = 0.252
        eta_v            = (1._r8-eta_0)*pi*0.5_r8
        u0               = 35._r8
        dtemp            = 4.8e5_r8
        temp0            = 288._r8
        temp_lapse_rate  = 0.005_r8
        part1            = u0*(cos(eta_v)**1.5)
        part2            = ((-2*(sin(lat)**6)*(cos(lat)**2+1._r8/3._r8))+10._r8/63._r8)*u0*(cos(eta_v)**1.5)
        part3            = (8._r8/5._r8*(cos(lat)**3)*((sin(lat)**2)+2._r8/3._r8)-pi/4._r8)*rearth*omega
        if(eta.le.eta_s.and.eta.ge.eta_t)then
           part0  = (temp0*gravity/temp_lapse_rate)*(1._r8-eta**(rdry*temp_lapse_rate/gravity))
        end if
        if(eta.lt.eta_t)then
           tmp = (log(eta/eta_t)+137._r8/60._r8)*(eta_t**5)-&
                  5._r8*(eta_t**4)*eta+5._r8*(eta_t**3)*(eta**2)-&
                  (10._r8/3._r8)*(eta_t**2)*(eta**3)+&
                  (5._r8/4._r8)*(eta_t)*(eta**4)-&
                  (eta**5)/5._r8
           part0  = (temp0*gravity/temp_lapse_rate)*(1._r8-eta**(rdry*temp_lapse_rate/gravity))-rdry*dtemp*tmp
        end if
        ini_geop = part0+part1*(part2+part3)
    case('DCMIP3-1')
        GG               = (gravity**2)/(1.e-4_r8*cp)
        surface_temp     = GG+(300._r8-GG)*exp(-0.25_r8*1.e-4_r8*20._r8*20._r8/(gravity**2)*(cos(2._r8*lat)-1._r8))
        surface_pressure = ini_surface_pressure(lon,lat,testcase)
        pres             = eta_a*1.e5_r8+eta_b*surface_pressure
        ini_geop         = -1.e4_r8*gravity*log((surface_temp/GG)*((pres/surface_pressure)**(rdry/cp)-1._r8)+1._r8)
        ini_geop         = ini_geop*gravity
    case default
       print*, "you must select a test case, please check!, ini_geop"
       stop
    end select

    return
  end function ini_geop

!================================================
!  U in West-East direction at a point
!  Lat in [-pi/2,pi/2], Lon in [-pi,pi]
!================================================
  real(r8) function ini_uwind(lon, lat, eta, testcase)
! io
    real(r8),       intent(in) :: lon
    real(r8),       intent(in) :: lat
    real(r8),       intent(in) :: eta
    character*(*),  intent(in) :: testcase
!
! local
!
    real(r8)                   :: u0
    real(r8)                   :: eta_0
    real(r8)                   :: eta_v
    real(r8)                   :: lonc, latc
    real(r8)                   :: up, up1
    real(r8)                   :: radius
    real(r8)                   :: center_a_lon, center_a_lat, center_b_lon, center_b_lat
    real(r8)                   :: dist_a, dist_b
    real(r8)                   :: modon_a, modon_b


! Velocity for each testcase
    select case(trim(testcase))

    case('JWSS') 
!
! jablonowski and williamson 2006 steady state test 
!
       u0        = 35._r8
       eta_0     = 0.252_r8
       eta_v     = (eta-eta_0)*pi*0.5_r8
       ini_uwind = u0*(cos(eta_v)**1.5_r8)*(sin(2._r8*lat)**2)

    case('JWBW') 
!
! jablonowski and williamson 2006 baroclinic test 
!
       u0        = 35._r8
       eta_0     = 0.252_r8
       eta_v     = (eta-eta_0)*pi*0.5_r8
       ini_uwind = u0*(cos(eta_v)**1.5_r8)*(sin(2._r8*lat)**2)
!
! pertubation of u
!
       lonc      = pi/9._r8
       latc      = 2._r8*pi/9._r8
       radius    = acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc))
       up        = exp(-(10._r8*radius)**2)
     
       ini_uwind = ini_uwind+up

    case('JWBW2') 
!
! jablonowski and williamson 2006 baroclinic test 
!
       u0        = 35._r8
       eta_0     = 0.252_r8
       eta_v     = (eta-eta_0)*pi*0.5_r8
       ini_uwind = u0*(cos(eta_v)**1.5_r8)*(sin(2._r8*lat)**2)
!
! pertubation of u
!
       lonc      = pi/9._r8
       latc      = 2._r8*pi/9._r8
       radius    = acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc))
       up        = exp(-(10._r8*radius)**2)
     
       lonc      = pi/9._r8
       latc      =-2._r8*pi/9._r8
       radius    = acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc))
       up1       = exp(-(10._r8*radius)**2)

       ini_uwind = ini_uwind+up+up1

    case('SCHAR','DCMIP2-1','DCMIP210','DCMIP3-1','MIRW','MIRWTOPO')
       u0        = 20._r8
       ini_uwind = u0*cos(lat)
    case('MODON')
       center_a_lat = 0._r8
       center_b_lat = 0._r8
       center_a_lon = pi*0.5_r8
       center_b_lon = pi*1.5_r8
       dist_a       = rearth*arcdistll(lon,lat,center_a_lon,center_a_lat)
       dist_b       = rearth*arcdistll(lon,lat,center_b_lon,center_b_lat)
       modon_a      = 40._r8*exp(-(dist_a/modon_radius)**2)
       modon_b      = 40._r8*exp(-(dist_b/modon_radius)**2)
       ini_uwind    = modon_a-modon_b
    case('DCMIP2-0')
       ini_uwind    = 0._r8
    case default
       print*, "you must select a test case, please check!, ini_uwind"
       stop
    end select

    return
  end function ini_uwind
!================================================
!  V in South-North direction at a point
!================================================
  real(r8) function ini_vwind(lon, lat, eta, testcase)
! io
    real(r8),       intent(in) :: lon
    real(r8),       intent(in) :: lat
    real(r8),       intent(in) :: eta
    character*(*),  intent(in) :: testcase

    select case(trim(testcase))

!================================================
!      Test Cases in Williamson et al. (1992)
!================================================
    case('JWSS','JWBW','JWBW2','SCHAR','DCMIP2-0','DCMIP2-1','DCMIP210','DCMIP3-1','MODON','MIRW','MIRWTOPO') ! jablonowski and williamson 2006 baroclinic test 
       ini_vwind = 0._r8
    case default
       print*, "you must select a test case, please check!, ini_vwind"
       stop
    end select

    return
  end function ini_vwind

  real(r8) function ini_stream_function(lon, lat, testcase)
!================================================
!  initial stream function defined at a point
!  see JW08
!================================================
! io
    real(r8),        intent(in)    :: lon
    real(r8),        intent(in)    :: lat
    character*(*),   intent(in)    :: testcase
! local
    real(r8)        :: mmm, kkk, nnn, u0
    
    select case(trim(testcase))
    case('RH3D')
        nnn = 4.
        u0  = 50._r8
        mmm = u0/(nnn*rearth)
        kkk = mmm
        ini_stream_function = -(rearth**2)*mmm*sin(lat)+&
                               (rearth**2)*kkk*(cos(lat)**nnn)*sin(lat)*cos(nnn*lon)
    case default
        ini_stream_function = 0._r8
    end select

    return
  end function ini_stream_function

  real(r8) function ini_surface_pressure(lon,lat,testcase,phisin)
! io
     real(r8),      intent(in)    :: lon
     real(r8),      intent(in)    :: lat
     character*(*), intent(in)    :: testcase
     real(r8), optional, intent(in)    :: phisin
! local
     real(r8)                     :: u0, temp0, temp_lapse_rate, pref
     real(r8)                     :: mmm,nnn,kkk, AAA,BBB,CCC 
     real(r8)                     :: geopotential_pert
     real(r8)                     :: ueq, lonc, h0, d0, kesi0, r0, zs
     real(r8)                     :: latc, r_m, rr_m

     select case(trim(testcase))
     case('RH3D')
         temp0           = 288._r8
         temp_lapse_rate = 0.0065_r8
         pref            = 95500._r8
         nnn             = 4.
         u0              = 50._r8
         mmm             = u0/(nnn*rearth)
         kkk             = mmm
! be careful
         AAA             = 0.5_r8*mmm*(2._r8*omega+mmm)*(cos(lat)**2)+0.25_r8*(kkk**2)*(cos(lat)**(2*nnn))*&
                           ((nnn+1)*(cos(lat)**2)+(2*nnn*nnn-nnn-2))-0.5*(nnn**2)*(kkk**2)*(cos(lat)**(2*nnn-2))

         BBB             = ((2._r8*(omega+mmm)*kkk)/((nnn+1)*(nnn+2)))*(cos(lat)**nnn)*((nnn*nnn+2*nnn+2)-((nnn+1)**2)*(cos(lat)**2))

         CCC             = 0.25*kkk*kkk*(cos(lat)**(2*nnn))*((nnn+1)*(cos(lat)**2)-(nnn+2))

         geopotential_pert = (rearth**2)*(AAA+BBB*cos(nnn*lon)+CCC*cos(2._r8*nnn*lon)) ! eq71
         ini_surface_pressure  = pref*(1._r8+temp_lapse_rate*geopotential_pert/(gravity*temp0))**(gravity/(temp_lapse_rate*rdry)) ! eq70
         !ini_surface_pressure  = (1._r8+temp_lapse_rate*geopotential_pert/(gravity*temp0))**(gravity/(temp_lapse_rate*rdry)) ! eq70
      case('SCHAR')
        ueq      = 20._r8
        lonc     = 0._r8
        h0       = 250._r8    ! meter
        d0       = 5000._r8
        kesi0    = 4000._r8
        r0       = rearth*(lon-lonc) ! as K15
        zs       = h0*exp(-(r0**2)/(d0**2))*(cos(pi*r0/kesi0)**2)*cos(lat)
        ini_surface_pressure = 1e5_r8*exp(-ueq*ueq*(sin(lat)**2)/(2._r8*rdry*300._r8)-gravity*zs/(rdry*300._r8))

      case('DCMIP2-1')
        ueq      = 20._r8
        lonc     = 0._r8
        h0       = 250._r8    ! meter
        d0       = 5000._r8
        kesi0    = 4000._r8
        r0       = rearth*acos(cos(lat)*cos(lon-lonc))
        zs       = h0*exp(-(r0**2)/(d0**2))*(cos(pi*r0/kesi0)**2)
        ini_surface_pressure = 1e5_r8*exp(-ueq*ueq*(sin(lat)**2)/(2._r8*rdry*300._r8)-gravity*zs/(rdry*300._r8))

      case('DCMIP3-1')
        AAA              = (gravity**2)/(1.e-4_r8*cp)
        temp0            = AAA+(300._r8-AAA)*exp(-0.25_r8*1.e-4_r8*20._r8*20._r8/(gravity**2)*(cos(2._r8*lat)-1._r8))
        ini_surface_pressure = 1.e5_r8*exp(20._r8*20._r8/(4._r8*AAA*rdry)*(cos(2._r8*lat)-1._r8))*&
                               (temp0/300._r8)**(cp/rdry)
      case('DCMIP2-0')
        zs  = ini_phis(lon,lat,'DCMIP2-0')
        zs  = zs/gravity
        ini_surface_pressure = 1.e5_r8*((1._r8-0.0065_r8*zs/300._r8)**(gravity/(rdry*0.0065_r8)))
      case('MIRW','MIRWTOPO')
        !zs  = ini_phis(lon,lat,'MIRW')
        zs  = phisin
        u0  = 20._r8
        nnn = 0.0182_r8
        ini_surface_pressure = 93000._r8*exp(-rearth*nnn*nnn*u0*0.5_r8/(2._r8/7._r8*gravity**2)*&
                                             (u0/rearth+2._r8*omega)*&
                                             (sin(lat)**2-1._r8)-nnn*nnn/(2._r8/7._r8*(gravity**2))*zs)
     case default

         print*,"you must select a test case"
     end select

     return
   end function ini_surface_pressure

  end module grist_dycore_initial_module
