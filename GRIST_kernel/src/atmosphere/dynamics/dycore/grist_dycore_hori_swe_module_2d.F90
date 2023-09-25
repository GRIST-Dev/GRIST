
!-------------------------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
!
! Description: Module for computation needed for SWE, 2D version of hori swe
! module
! Revision history:
!   1. A revised workflow so that PV2 is only working as PV2
!   2. If PV2 used, must be used with a proper hyperdiffusion configuration (otherwise JWBW solution 
!      is a disaster $_$
!   3. Recommended setup (based on testing efforts and experience):
!      (1) USE PV3 for G5-G8 multiyear climate runs and donot use any hyperdiffusion
!         (PV2+HYPER will damp EMF and diffuse KE at the tail; PV3's time increase is acceptable
!          because at these resolutions, dynamics is not well dominating computation)
!          PV3 also improves HDC/NDC real-world multiyear modeling stability than PV2 (G5,G6 tested)
!          at certain configurations (tested with WRF_physics).
!      (2) USE PV2+hyperdiffusion for >=G9 and high-res VR modeling at seasonal or duration slightly beyond/below
!          (impact of hyper is small for such modeling; dynamics is dominating and we should do 
!          "everything possible" to reduce its cost for the time being; PV2+hyperdiffusion is 
!          better than using NCT-PV3 once per step (only for TC test, need further check in realworld).
!          hybrid PV2&PV3 in each RK pass is also included for testing now.
!           
!          IN SHORT, IF BEST PERFORMANCE IS WANTED, JUST USE PV3 (can be active at the final RK stage; similar results as full-call in JWBW); 
!          IF SAVING TIME IS IN AN URGENT NEED, use PV2+hyperdiffusion as a poor-man version. (Note that using PV3 will increase 
!          the speed of other components (cache issue?), so the time increase is lower than that seen by time proportion statistics)
!          See Zhang (2018, JAMES) for PV details.
!
!          Test in future: further reduce PV3's stencil from 9 to 6: seems not saving time and not worsing performance...
!      (3) Please always use USE_HALO2 as default.
!-------------------------------------------------------------------------------------------

 module grist_dycore_hori_swe_module_2d

   use grist_domain_types  ,  only: global_domain
   use grist_data_types,      only: scalar_2d_field
   use grist_constants,       only: i4, r8, rearth, omega, one, zero
   use grist_nml_module,      only: pv_order, mass_wtd_vor, conserve_scheme, ke_method, ke_alpha, stencil_width
   use grist_base_flux_module,only: flux_ooo33
#ifdef USE_HALO2
   use grist_data_types      ,only: exchange_field_list_2d, scalar_2d_field
   use grist_config_partition,only: exchange_data_2d_add ,  exchange_data_2d
#endif
   use grist_lib

   implicit none

   private

   public :: calc_tend_nct_at_edge_full_level, &
             calc_grad_kinetic_energy        , &
             calc_coriolis_term

   contains

  subroutine calc_tend_nct_at_edge_full_level(mesh,nlev,irk_step,dtime                 ,&
                                              scalar_height_at_pc_full_level           ,&
                                              scalar_normal_velocity_at_edge_full_level,&
                                              scalar_normal_flux_at_edge_full_level    ,&
                                              tend_nct_at_edge_full_level             ) 
!io
   type(global_domain)  , intent(inout) :: mesh
   integer(i4)          , intent(in)    :: nlev, irk_step
   real(r8)             , intent(in)    :: dtime
   real(r8)             , intent(in)    :: scalar_height_at_pc_full_level(:,:)
   real(r8)             , intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8)             , intent(in)    :: scalar_normal_flux_at_edge_full_level(:,:)
   real(r8)             , intent(inout) :: tend_nct_at_edge_full_level(:,:)
!local-exchange
   type(scalar_2d_field)                :: scalar_pv_at_dc_full_level
   type(scalar_2d_field)                :: scalar_pv_at_edge_full_level
!local-nonexchange
   real(r8)                             :: scalar_tangen_velocity_at_edge_full_level(nlev,1:mesh%ne_full)
   real(r8)                             :: scalar_height_at_dc_full_level(nlev,1:mesh%nt_full)
   real(r8)                             :: scalar_avor_at_dc_full_level(  nlev,1:mesh%nt_full)
   real(r8)                             :: scalar_rvor_at_dc_full_level(  nlev,1:mesh%nt_full)
   real(r8)                             :: scalar_coriolis_term_at_edge(  nlev,1:mesh%ne_full)   ! what we need
   real(r8)                             ::                          tmp(  nlev,1:mesh%ne_full)
   integer(i4)                          :: it, ie, iv
   real(r8)                             :: sum_tmp(nlev)
   real(r8)                             :: length_of_edge
   integer(i4)                          :: ilev, index_of_edge, t_edge_v, cell_index
#ifdef USE_HALO2
   type(exchange_field_list_2d),pointer :: field_head_2d
   field_head_2d=>null()
#endif

   if(.not.allocated(scalar_pv_at_edge_full_level%f))    allocate(scalar_pv_at_edge_full_level%f(nlev,1:mesh%ne_full))
   if(.not.allocated(scalar_pv_at_dc_full_level%f))      allocate(scalar_pv_at_dc_full_level%f(nlev,1:mesh%nt_full))

   scalar_avor_at_dc_full_level       = zero
   scalar_height_at_dc_full_level     = one
   scalar_coriolis_term_at_edge       = zero
   scalar_pv_at_edge_full_level%f     = zero  ! this line is necessary for LAM correct initialization

!
! calculate vorticity at dual cell using Stokes' Theorem
! only until stencil_width-1 because calc_scalar_at_dual_cell needs one-more
! layer
!
     mesh%ne = mesh%ne_halo(stencil_width-1)
     mesh%nt = mesh%nt_halo(stencil_width-1)
! inline_opt_here
!     call calc_vorticity_at_dual_cell(mesh, scalar_normal_velocity_at_edge_full_level%f ,&
!                                            scalar_avor_at_dc_full_level                ,&
!                                            scalar_rvor_at_dc_full_level                ,&
!                                            nlev)

!$omp parallel  private(it,sum_tmp,ie,ilev,index_of_edge,t_edge_v)
!$omp do schedule(dynamic,5)
      do it = 1, mesh%nt
         sum_tmp(:) = zero
         do ie = 1, mesh%tri_nnb(it)
            index_of_edge          = mesh%tri_ed(ie,it)
            t_edge_v               = mesh%edt_edpTg_edtNr(index_of_edge)*mesh%tri_nr(ie,it)
            do ilev = 1, nlev
#ifndef DYCORE_OPT
               sum_tmp(ilev)       = sum_tmp(ilev)+t_edge_v*rearth*mesh%edt_leng(index_of_edge)*scalar_normal_velocity_at_edge_full_level(ilev,index_of_edge)
#else
               sum_tmp(ilev)       = sum_tmp(ilev)+t_edge_v*mesh%edt_leng(index_of_edge)*scalar_normal_velocity_at_edge_full_level(ilev,index_of_edge)
#endif
            end do
         end do
#ifndef DYCORE_OPT
         scalar_rvor_at_dc_full_level(1:nlev,it) = sum_tmp(1:nlev)/((rearth**2)*mesh%tri_areag(it))
#else
         scalar_rvor_at_dc_full_level(1:nlev,it) = sum_tmp(1:nlev)/(rearth*mesh%tri_areag(it))
#endif
         scalar_avor_at_dc_full_level(1:nlev,it) = scalar_rvor_at_dc_full_level(1:nlev,it)+2._r8*omega*sin(mesh%tri_c_lat(it))
      end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here

     if(mass_wtd_vor)then
! inline_opt_here
      !  call calc_scalar_at_dual_cell(mesh,scalar_height_at_pc_full_level%f, &
      !                                     scalar_height_at_dc_full_level  , &
      !                                     nlev)

!$omp parallel  private(it,sum_tmp,ie,ilev,cell_index)
!$omp do schedule(dynamic,5)
      do it = 1, mesh%nt
         sum_tmp(:) = 0._r8
         do ie = 1, mesh%tri_nnb(it)
            cell_index       = mesh%tri_v(ie,it)
            do ilev = 1, nlev
               sum_tmp(ilev) = sum_tmp(ilev)+(mesh%tri_kite_area(ie,it))*scalar_height_at_pc_full_level(ilev,cell_index)
            end do
         end do
         scalar_height_at_dc_full_level(1:nlev,it) = sum_tmp(1:nlev)/(mesh%tri_areag(it))
      end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here inline
        scalar_pv_at_dc_full_level%f = scalar_avor_at_dc_full_level/scalar_height_at_dc_full_level
     else
        scalar_pv_at_dc_full_level%f = scalar_avor_at_dc_full_level
     end if
!
! reconstruct tangent wind using nct operator, from stencil_width to stencil_width-1
!
     if(pv_order(irk_step).eq.3)then
! pv-dc only having until stencil_width-1 layer, so exchange data for PV3
#ifdef USE_HALO2
        call exchange_data_2d_add(mesh,field_head_2d,scalar_pv_at_dc_full_level)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

#ifdef OLDNCT
        tmp(:,:) = 1._r8
        call calc_coriolis_term(mesh,tmp  ,scalar_normal_velocity_at_edge_full_level, & ! alto both vel and mass_flux can use (only sign is needed), 
                                                                                          ! we cannot use normal mass flux which is only updated at i
                                                                                          ! ne_halo(1), and thus tangent can only be at ne_compute, 
                                                                                          ! which is not enough
                                           scalar_tangen_velocity_at_edge_full_level, &
                                           nlev)
#else
        call calc_coriolis_term_opt(mesh,tmp,scalar_normal_velocity_at_edge_full_level,scalar_tangen_velocity_at_edge_full_level,nlev,.true.)
#endif

     end if
!
! tangent flux along hx's edge
! We can just compute PV-edge until ne_halo(1) without exch
!
     mesh%ne = mesh%ne_halo(1)
     mesh%nt = mesh%nt_halo(1)

! old code compute until ne_compute for pv3 and use_halo2, which is overly conservative
#if (defined PV3_OLDCODE || defined REG_QJLAM)
#ifdef USE_HALO2
     if(pv_order(irk_step).eq.3)then
        mesh%ne = mesh%ne_compute
        mesh%nt = mesh%nt_compute
     end if
#endif
#endif

!
! we need pv flux at ne_halo(1)
!
     call calc_normal_flux_tr_edge(mesh, scalar_pv_at_dc_full_level%f              ,&
                                         scalar_tangen_velocity_at_edge_full_level ,&
                                         scalar_normal_velocity_at_edge_full_level ,&
                                         scalar_pv_at_edge_full_level%f              ,&
                                         .false.                                     ,&
                                         pv_order(irk_step), dtime, nlev)

! old code compute until ne_compute for pv3 and use_halo2, which is overly conservative
#if (defined PV3_OLDCODE || defined REG_QJLAM)
     if(pv_order(irk_step).eq.3)then
!exchange to get ne_halo(1), if halo3, no need exchange
#ifdef USE_HALO2
        call exchange_data_2d_add(mesh,field_head_2d,scalar_pv_at_edge_full_level)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
     end if
#endif

!================================================
!       The Final nonlinear Coriolis tendency
!================================================

!
! from ne_halo(1) to ne_compute
! imp: set this to ne is necessary for next KE-grad evaluation
!
     mesh%ne = mesh%ne_compute
     mesh%nt = mesh%nt_compute

     if(mass_wtd_vor)then ! default

#ifdef OLDNCT
        call calc_coriolis_term(mesh,scalar_pv_at_edge_full_level%f,scalar_normal_flux_at_edge_full_level,scalar_coriolis_term_at_edge,nlev)
#else
        call calc_coriolis_term_opt(mesh,scalar_pv_at_edge_full_level%f,scalar_normal_flux_at_edge_full_level,scalar_coriolis_term_at_edge,nlev,.false.)
#endif

     else
#ifdef OLDNCT
        call calc_coriolis_term(mesh,scalar_pv_at_edge_full_level%f,scalar_normal_velocity_at_edge_full_level,scalar_coriolis_term_at_edge,nlev)
#else
        call calc_coriolis_term_opt(mesh,scalar_pv_at_edge_full_level%f,scalar_normal_flux_at_edge_full_level,scalar_coriolis_term_at_edge,nlev,.false.)
#endif
     end if

! send back
      tend_nct_at_edge_full_level  = -1._r8*scalar_coriolis_term_at_edge

! clean 
      deallocate(scalar_pv_at_edge_full_level%f)
      deallocate(scalar_pv_at_dc_full_level%f)

    return
  end subroutine calc_tend_nct_at_edge_full_level

!==============================================
!   Calculate grad of kinetic energy at edge
!==============================================

  subroutine calc_grad_kinetic_energy(mesh, nlev, &
                                      scalar_normal_velocity_at_edge_full_level,&
                                      tend_grad_ke_at_edge_full_level )
! io
   use omp_lib
   type(global_domain),   intent(inout) :: mesh
   integer(i4),           intent(in)    :: nlev
   real(r8),              intent(in)    :: scalar_normal_velocity_at_edge_full_level(nlev,mesh%ne_full)
   real(r8),              intent(inout) :: tend_grad_ke_at_edge_full_level(:,:)
! local
   !real(r8), allocatable                :: kinetic_energy_hx(:,:)
   !real(r8), allocatable                :: kinetic_energy_tr(:,:)
   !real(r8), allocatable                :: kinetic_energy(:,:)
   real(r8)                             :: kinetic_energy_hx(nlev,mesh%nv)
   real(r8)                             :: kinetic_energy_tr(nlev,mesh%nt)
   real(r8)                             :: kinetic_energy(   nlev,mesh%nv)
   real(r8)                             :: length_of_voronoi_edge
   real(r8)                             :: length_of_triangle_edge
   real(r8)                             :: tmp_sum(nlev)
   real(r8)                             :: wind, flag
   real(r8)                             :: v1v2(3)
   integer                              :: ilev, it, ie, iv
   integer(i4)                          :: index_edge, v1, v2

!   if(.not.allocated(kinetic_energy)) allocate(kinetic_energy(nlev,mesh%nv))

   select case(ke_method)
   case (0)   ! R10's formulation

      do iv = 1, mesh%nv
         tmp_sum(:)  = 0._r8
         do ie = 1, mesh%vtx_nnb(iv)
            index_edge              = mesh%vtx_ed(ie, iv)   ! global index of this edge pair
#ifndef DYCORE_OPT
            length_of_triangle_edge = rearth*mesh%edt_leng(index_edge)
            length_of_voronoi_edge  = rearth*mesh%edp_leng(index_edge)
#else
            length_of_triangle_edge = mesh%edt_leng(index_edge)
            length_of_voronoi_edge  = mesh%edp_leng(index_edge)
#endif
            do ilev = 1, nlev
               wind                 = scalar_normal_velocity_at_edge_full_level(ilev,index_edge)
               tmp_sum(ilev)        = tmp_sum(ilev)+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
            end do
         end do
#ifndef DYCORE_OPT
         kinetic_energy(:,iv) = tmp_sum(:)/((rearth**2)*mesh%plg_areag(iv))
#else
         kinetic_energy(:,iv) = tmp_sum(:)/(mesh%plg_areag(iv))
#endif
      end do

   case(1)   ! s12/g13's formulation

!     if(.not.allocated(kinetic_energy_hx)) allocate(kinetic_energy_hx(nlev,mesh%nv))
!     if(.not.allocated(kinetic_energy_tr)) allocate(kinetic_energy_tr(nlev,mesh%nt))
!
! r10's formulation, again
!
!$omp parallel  private(iv,tmp_sum,ie,index_edge,length_of_triangle_edge,length_of_voronoi_edge,ilev,wind) 
!$omp do schedule(dynamic,5)
        do iv = 1, mesh%nv
           tmp_sum(:)  = 0._r8
           do ie = 1, mesh%vtx_nnb(iv)
              index_edge              = mesh%vtx_ed(ie, iv)   ! global index of this edge pair 
#ifndef DYCORE_OPT
              length_of_triangle_edge = rearth*mesh%edt_leng(index_edge)
              length_of_voronoi_edge  = rearth*mesh%edp_leng(index_edge)
#else
              length_of_triangle_edge = mesh%edt_leng(index_edge)
              length_of_voronoi_edge  = mesh%edp_leng(index_edge)
#endif
              do ilev = 1, nlev
                 wind                 = scalar_normal_velocity_at_edge_full_level(ilev,index_edge)
                 tmp_sum(ilev)        = tmp_sum(ilev)+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
              end do
           end do
#ifndef DYCORE_OPT
           kinetic_energy(:,iv) = tmp_sum(:)/((rearth**2)*mesh%plg_areag(iv))
#else
           kinetic_energy(:,iv) = tmp_sum(:)/(mesh%plg_areag(iv))
#endif
        end do
!$omp end do nowait
!$omp end parallel 

!
! compute ke at cell vertex/triangular center
!
!$omp parallel  private(it,tmp_sum,ie,index_edge,length_of_triangle_edge,length_of_voronoi_edge,ilev,wind)  
!$omp do schedule(dynamic,5)
        do it = 1, mesh%nt
           tmp_sum(:)  = 0._r8
           do ie = 1, mesh%tri_nnb(it)
              index_edge              = mesh%tri_ed(ie, it)
              length_of_voronoi_edge  = mesh%edp_leng(index_edge)
              length_of_triangle_edge = mesh%edt_leng(index_edge)
              do ilev = 1, nlev
                 wind                 = scalar_normal_velocity_at_edge_full_level(ilev,index_edge)
                 tmp_sum(ilev)        = tmp_sum(ilev) + 0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*(wind**2)
              end do
           end do
              kinetic_energy_tr(:,it) = tmp_sum(:)/mesh%tri_areag(it)
        end do
!$omp end do nowait
!$omp end parallel
!
! remap ke from cell vertex to cell center
!
!$omp parallel  private(iv,tmp_sum,ie,ilev)
!$omp do schedule(dynamic,5)
        do iv = 1, mesh%nv
           tmp_sum(:) = 0._r8
           do ie = 1, mesh%vtx_nnb(iv)  ! local index of each corner/tr
              do ilev = 1, nlev
                 tmp_sum(ilev) = tmp_sum(ilev)+kinetic_energy_tr(ilev,mesh%vtx_tr(ie, iv))*mesh%plg_kite_area(ie, iv)
              end do
           end do
           kinetic_energy_hx(:,iv) = tmp_sum(:)
! merge
           kinetic_energy(:,iv)    = ke_alpha*kinetic_energy(:,iv)+(1._r8-ke_alpha)*kinetic_energy_hx(:,iv)
        end do
!$omp end do nowait
!$omp end parallel

!        deallocate(kinetic_energy_hx)
!        deallocate(kinetic_energy_tr)

   case default

       print*, "you must select a ke_method, stop"
       stop
   end select
!
! gradient of KE
!
   mesh%nv = mesh%nv_compute
   mesh%nt = mesh%nt_compute

   !do ie = 1, mesh%ne
   !   v1      = mesh%edt(ie)%v(1)
   !   v2      = mesh%edt(ie)%v(2)
   !   v1v2    = mesh%vtx(v2)%p-mesh%vtx(v1)%p
   !   flag    = sign(1._r8,dot_product(v1v2,mesh%edp(ie)%nr))
   !   do ilev = 1, nlev
   !      tend_grad_ke_at_edge_full_level(ilev,ie) = -flag*(kinetic_energy(ilev,v2)-kinetic_energy(ilev,v1))/(rearth*mesh%edt(ie)%leng)
   !   end do
   !end do
!$omp parallel  private(ie,v1,v2,ilev)  
!$omp do schedule(dynamic,5)
   do ie = 1, mesh%ne
      v1      = mesh%edt_v(1, ie)
      v2      = mesh%edt_v(2, ie)
      do ilev = 1, nlev
         tend_grad_ke_at_edge_full_level(ilev,ie) = -mesh%edt_edpNr_edtTg(ie)*(kinetic_energy(ilev,v2)-kinetic_energy(ilev,v1))/(rearth*mesh%edt_leng(ie))
      end do
   end do
!$omp end do nowait
!$omp end parallel 

!   deallocate(kinetic_energy)
   return
  end subroutine calc_grad_kinetic_energy

!================================================
!
!            PRIVATE SUBROUTINES
!
!================================================

  subroutine calc_vorticity_at_dual_cell(mesh, scalar_normal_velocity_at_edge_full_level,&
                                               scalar_avor_at_dc_full_level,&
                                               scalar_rvor_at_dc_full_level,&
                                               nlev)
! io
   use omp_lib
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8),              intent(inout) :: scalar_avor_at_dc_full_level(nlev,mesh%nt_full)
   real(r8),              intent(inout) :: scalar_rvor_at_dc_full_level(nlev,mesh%nt_full)
   integer(i4)          , intent(in)    :: nlev
! local
   real(r8)                             :: coriolis
   real(r8)                             :: sum_tmp(nlev)
   real(r8)                             :: length_of_edge
   integer(i4)                          :: it, ie, ilev
   integer(i4)                          :: index_of_edge
   integer(i4)                          :: t_edge_v


   !do it = 1, mesh%nt
   !   sum_tmp(:) = 0._r8
   !   do ie = 1, 3
   !      do ilev = 1, nlev 
   !         index_of_edge          = mesh%tri(it)%ed(ie)
   !         t_edge_v               = mesh%edt(index_of_edge)%edpTg_edtNr*mesh%tri(it)%nr(ie)
   !         sum_tmp(ilev)          = sum_tmp(ilev)+t_edge_v*rearth*mesh%edt(index_of_edge)%leng*&
   !                                  scalar_normal_velocity_at_edge_full_level(ilev,index_of_edge)
   !      end do
   !   end do
   !   scalar_rvor_at_dc_full_level(:,it) = sum_tmp(:)/((rearth**2)*mesh%tri(it)%areag)
   !   coriolis                           = 2._r8*omega*sin(mesh%tri(it)%c%lat)
   !   scalar_avor_at_dc_full_level(:,it) = scalar_rvor_at_dc_full_level(:,it)+coriolis
   !end do
!$omp parallel  private(it,sum_tmp,ie,ilev,index_of_edge,t_edge_v,coriolis) 
!$omp do schedule(dynamic,5)
   do it = 1, mesh%nt
      sum_tmp(:) = 0._r8
      do ie = 1, mesh%tri_nnb(it)
         do ilev = 1, nlev
            index_of_edge          = mesh%tri_ed(ie, it)
            t_edge_v               = mesh%edt_edpTg_edtNr(index_of_edge)*mesh%tri_nr(ie, it)
            sum_tmp(ilev)          = sum_tmp(ilev)+t_edge_v*rearth*mesh%edt_leng(index_of_edge)*&
                                     scalar_normal_velocity_at_edge_full_level(ilev,index_of_edge)
         end do
      end do
      scalar_rvor_at_dc_full_level(:,it) = sum_tmp(:)/((rearth**2)*mesh%tri_areag(it))
      coriolis                           = 2._r8*omega*sin(mesh%tri_c_lat(it))
      scalar_avor_at_dc_full_level(:,it) = scalar_rvor_at_dc_full_level(:,it)+coriolis
    end do
!$omp end do nowait
!$omp end parallel 

    return
  end subroutine calc_vorticity_at_dual_cell

  subroutine calc_scalar_at_dual_cell(mesh,scalar_height_at_pc_full_level,&
                                           scalar_height_at_dc_full_level,&
                                           nlev)
! io
   use omp_lib
   type(global_domain),   intent(in)   :: mesh
   real(r8),              intent(in)   :: scalar_height_at_pc_full_level(nlev,mesh%nv_full)
   real(r8),              intent(inout):: scalar_height_at_dc_full_level(nlev,mesh%nt_full)
   integer(i4)          , intent(in)   :: nlev
! local
   integer(i4)                         :: it, ie, ilev
   integer(i4)                         :: cell_index
   real(r8)                            :: sum_tmp(nlev)


   !do it = 1, mesh%nt
   !  sum_tmp(:) = 0._r8
   !  do ie = 1, 3
   !     do ilev = 1, nlev 
   !        cell_index       = mesh%tri(it)%v(ie)
   !        sum_tmp(ilev)    = sum_tmp(ilev)+(mesh%tri(it)%kite_area(ie))*scalar_height_at_pc_full_level(ilev,cell_index)
   !     end do
   !  end do
   !     scalar_height_at_dc_full_level(:,it) = sum_tmp(:)/(mesh%tri(it)%areag)
   !end do
!$omp parallel  private(it,sum_tmp,ie,ilev,cell_index) 
!$omp do schedule(dynamic,5)
   do it = 1, mesh%nt
     sum_tmp(:) = 0._r8
     do ie = 1, mesh%tri_nnb(it)
        do ilev = 1, nlev 
           cell_index       = mesh%tri_v(ie,it)
           sum_tmp(ilev)    = sum_tmp(ilev)+(mesh%tri_kite_area(ie,it))*scalar_height_at_pc_full_level(ilev,cell_index)
        end do
     end do
        scalar_height_at_dc_full_level(:,it) = sum_tmp(:)/(mesh%tri_areag(it))
   end do
!$omp end do nowait
!$omp end parallel 

   return
  end subroutine calc_scalar_at_dual_cell

  subroutine calc_coriolis_term(mesh,scalar_pv_at_edge_full_level           ,&
                                     scalar_normal_flux_at_edge_full_level  ,&
                                     scalar_coriolis_term_at_edge_full_level,&
                                     nlev)

! io 
   use omp_lib
   type(global_domain),     intent(in)    :: mesh
   real(r8),                intent(in)    :: scalar_pv_at_edge_full_level(:,:)
   real(r8),                intent(in)    :: scalar_normal_flux_at_edge_full_level(:,:)
   real(r8),                intent(inout) :: scalar_coriolis_term_at_edge_full_level(:,:)
   integer(i4)            , intent(in)    :: nlev
! local
   real(r8)                               :: length_of_triangle
   real(r8)                               :: pv_at_edge(16)
   real(r8)                               :: cell_sum(16)
   integer(i4)                            :: ie, ilev 

!
! te conserving scheme
!
   select case(trim(conserve_scheme)) 
   case('te')

!$omp parallel  private(ie,length_of_triangle,ilev,pv_at_edge,cell_sum)
!$omp do schedule(dynamic,5)
     do ie = 1, mesh%ne  ! for each edge e
        length_of_triangle = rearth*mesh%edt_leng(ie)
        do ilev = 1, nlev

           pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie)
           cell_sum(1:mesh%edp_nedge(ie))   = mesh%edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
   scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                                (pv_at_edge(1:mesh%edp_nedge(ie))+scalar_pv_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie)))/2._r8
           scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
        end do
     end do
!$omp end do nowait
!$omp end parallel 

!
! pe conserving scheme
!
   case('pe')

     do ie = 1, mesh%ne  ! for each edge e
        length_of_triangle = rearth*mesh%edt_leng(ie)
        do ilev = 1, nlev
           pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie)
           cell_sum(1:mesh%edp_nedge(ie))   = mesh%edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
   scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                                 pv_at_edge(1:mesh%edp_nedge(ie))
           scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
        end do
     end do
   case default
     if(mpi_rank().eq.0) print*,"conserve_scheme in swe_para not set, must set to 'te' or 'pe'" 
   end select

   return
  end subroutine calc_coriolis_term

  subroutine calc_normal_flux_tr_edge(mesh,scalar_pv_at_dc_full_level                ,& ! dual-cell value 
                                           scalar_normal_velocity_at_edge_full_level ,& ! normal wind at tr edge
                                           scalar_tangen_velocity_at_edge_full_level ,& ! tangent wind at tr edge, for ffsl
                                           scalar_pv_at_edge_full_level              ,& ! edge-based value mapped from dual cell
                                           called_by_tr                              ,& ! logical
                                           pv_flag                                   ,&
                                           dtime ,nlev )
! io
   use omp_lib
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_pv_at_dc_full_level(:,:)
   real(r8),              intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8),              intent(in)    :: scalar_tangen_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(inout) :: scalar_pv_at_edge_full_level(:,:)
   logical,               intent(in)    :: called_by_tr
   integer(i4),           intent(in)    :: pv_flag
   real(r8),              intent(in)    :: dtime
   integer(i4),           intent(in)    :: nlev
! local
   real(r8)                             :: v1v2(3)           ! direction 0->1
   real(r8)                             :: der0, der1
   real(r8)                             :: flag
   real(r8)                             :: wind_edge
   integer(i4)                          :: ie, ilev
   integer(i4)                          :: v1, v2
   integer(i4)                          :: mark, v_upwind
   integer(i4)                          :: ii

!================================================
!       potential vorticity at hx's edge
!================================================

      select case(pv_flag)

      case(2) ! PV2 
        do ie = 1, mesh%ne
             !v1 = mesh%edp(ie)%v(1)
             !v2 = mesh%edp(ie)%v(2)
             v1 = mesh%edp_v(1, ie)
             v2 = mesh%edp_v(2, ie)
             do ilev = 1, nlev
                scalar_pv_at_edge_full_level(ilev,ie) = 0.5*(scalar_pv_at_dc_full_level(ilev,v1)+scalar_pv_at_dc_full_level(ilev,v2))
             end do
        end do
      case(3) ! PV3

!$omp parallel  private(ii,ie,ilev,v1,v2,wind_edge,mark,v_upwind,der0) 
!$omp do schedule(dynamic,100) 
     do ii = 1, mesh%ne*nlev,1
        ie=ceiling(ii/real(nlev,r8))
        ilev=ii-(ie-1)*nlev
!       do ie = 1, mesh%ne

           !v1        = mesh%edp(ie)%v(1)       ! 1st end of edge
           !v2        = mesh%edp(ie)%v(2)       ! 2nd end of edge
! Make sure! wind follows C0C1 positive
           !v1v2      = mesh%tri(v2)%c%p-mesh%tri(v1)%c%p
! called by! SWE solver on HX for PV
           !flag      = sign(one,dot_product(v1v2,mesh%edp(ie)%tg))
! called by! ADV solver on TR
!          ! if(called_by_tr) flag  = sign(one,dot_product(v1v2,mesh%edt(ie)%nr))

           !do ilev = 1, nlev

           !    wind_edge = flag*scalar_normal_velocity_at_edge_full_level(ilev, ie)
           !    mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
           !    v_upwind  = mesh%edp(ie)%v(mark)    ! global index

           !    call calc_der_2nd_at_tr(mesh,scalar_pv_at_dc_full_level(ilev,mesh%tri(v_upwind)%stencil_index_2nd(:)),&
           !                                 scalar_pv_at_dc_full_level(ilev,v_upwind), v_upwind, &
           !                                 ie,mark-1,mesh%tri(v_upwind)%stencil_number_2nd,der0)

           !    call flux_ooo33(scalar_pv_at_dc_full_level(ilev,v1),&
           !                   scalar_pv_at_dc_full_level(ilev,v2),&
           !                   der0 ,&
           !                   scalar_pv_at_edge_full_level(ilev,ie), &
           !                   mesh%edp(ie)%leng)
           !end do

           v1        = mesh%edp_v(1, ie)       ! 1st end of edge
           v2        = mesh%edp_v(2, ie)       ! 2nd end of edge
! Make sure wind follows C0C1 positive
           !v1v2      = mesh%tri(v2)%c%p-mesh%tri(v1)%c%p
! called by! SWE solver on HX for PV
           !flag      = sign(one,dot_product(v1v2,mesh%edp(ie)%tg))
! called by ADV solver on TR
!           if(called_by_tr) flag  = sign(one,dot_product(v1v2,mesh%edt(ie)%nr))

!           do ilev = 1, nlev

               wind_edge = mesh%edp12_edp_tg(ie)*scalar_normal_velocity_at_edge_full_level(ilev, ie)
               mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)  ! 1 or 2
               v_upwind  = mesh%edp_v(mark, ie)    ! global index

               !call calc_der_2nd_at_tr(mesh,scalar_pv_at_dc_full_level(ilev,mesh%tri(v_upwind)%stencil_index_2nd(:)),&
! inline_opt_here
             !  call calc_der_2nd_at_tr(mesh,scalar_pv_at_dc_full_level(ilev,mesh%tri_stencil_index_2nd(1:mesh%tri_stencil_number_2nd(v_upwind),v_upwind)),&
             !                               scalar_pv_at_dc_full_level(ilev,v_upwind), v_upwind, &
             !                               ie,mark-1,mesh%tri_stencil_number_2nd(v_upwind),der0)
             !  call flux_ooo33(scalar_pv_at_dc_full_level(ilev,v1),&
             !                 scalar_pv_at_dc_full_level(ilev,v2),&
             !                 der0 ,&
             !                 scalar_pv_at_edge_full_level(ilev,ie), &
             !                 real(mesh%edp_leng(ie),r8))
! inline written here
               der0   = 2._r8*dot_product(mesh%edp_v_weight_dual_cell(1:mesh%tri_stencil_number_2nd(v_upwind), mark, ie), &
                              scalar_pv_at_dc_full_level(ilev,mesh%tri_stencil_index_2nd(1:mesh%tri_stencil_number_2nd(v_upwind),v_upwind))-&
                              scalar_pv_at_dc_full_level(ilev,v_upwind))
               scalar_pv_at_edge_full_level(ilev,ie) = (scalar_pv_at_dc_full_level(ilev,v1)+scalar_pv_at_dc_full_level(ilev,v2))*0.5_r8-&
                            ((real(mesh%edp_leng(ie),r8)**2)/6._r8)*der0
! inline_opt_here

!           end do

!          call calc_der_2nd_at_tr(mesh,scalar_pv_at_dc_full_level,v2,ie,1,der1)
!================================================
!             Compute edge-based value 
!================================================
!          call flux_ooo34(wind_edge,scalar_pv_at_dc_full_level(v1),&
!                                    scalar_pv_at_dc_full_level(v2),&
!                          der0,der1,scalar_pv_at_edge(ie),     &
!                          mesh%edp(ie)%leng,beta)

       !end do
     end do
!$omp end do nowait
!$omp end parallel 

      case default
         print*,"PLEASE SET PV_ORDER IN NAMELIST!!! GRIST STOPS..."
         stop
      end select

     return
  end subroutine calc_normal_flux_tr_edge

  subroutine calc_der_2nd_at_tr(mesh,s_array,s0,it,ie,flag,nlength,der)
! io
    type(global_domain),   intent(in) :: mesh
    real(r8),              intent(in) :: s_array(:)
    real(r8),              intent(in) :: s0
    integer(i4),           intent(in) :: it
    integer(i4),           intent(in) :: ie
    integer(i4),           intent(in) :: flag
    integer(i4),           intent(in) :: nlength
    real(r8),              intent(out):: der
! local
    real(r8)                       :: m_array(1:nlength)
    real(r8)                       :: f_array

      m_array(:) = s_array-s0
      !f_array    = dot_product(mesh%edp(ie)%v_weight_dual_cell(flag+1,1:mesh%tri(it)%stencil_number_2nd), m_array(:))
      f_array    = dot_product(mesh%edp_v_weight_dual_cell(1:mesh%tri_stencil_number_2nd(it), flag+1, ie), m_array(:))
      der        = 2._r8*f_array

     return
   end subroutine calc_der_2nd_at_tr

!=================================================================
! Default now: an optimized calc_coriolis_term (it knows itself
! being used twice for different purposes)
!=================================================================

  subroutine calc_coriolis_term_opt(mesh,scalar_pv_at_edge_full_level           ,&
                                         scalar_normal_flux_at_edge_full_level  ,&
                                         scalar_coriolis_term_at_edge_full_level,&
                                         nlev, flag)

! io 
   use omp_lib
   type(global_domain),     intent(in)    :: mesh
   real(r8),                intent(in)    :: scalar_pv_at_edge_full_level(:,:)
   real(r8),                intent(in)    :: scalar_normal_flux_at_edge_full_level(:,:)
   real(r8),                intent(inout) :: scalar_coriolis_term_at_edge_full_level(:,:)
   integer(i4)            , intent(in)    :: nlev
   logical, intent(in)                    :: flag
! local
   !real(r8)                               :: length_of_triangle
   real(r8)                               :: pv_at_edge(16)
   real(r8)                               :: cell_sum(16)
   integer(i4)                            :: ie, ilev

!
! te conserving scheme
!
   select case(trim(conserve_scheme))

   case('te')

   IF(flag)then ! 1st-time called in this module, for reconstructing tangent wind

!$omp parallel  private(ie,length_of_triangle,ilev,pv_at_edge,cell_sum)
!$omp do schedule(dynamic,5)
     do ie = 1, mesh%ne  ! for each edge e
        !length_of_triangle = rearth*mesh%edt_leng(ie)
        do ilev = 1, nlev
           cell_sum(1:mesh%edp_nedge(ie))   = mesh%edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie),ie)*&
                                              scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie),ie))
           scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/(rearth*mesh%edt_leng(ie)) ! tanget flux/vel, just for checking sign
        end do
     end do
!$omp end do nowait
!$omp end parallel

   ELSE ! 2nd time called, for real nct term

!$omp parallel  private(ie,length_of_triangle,ilev,pv_at_edge,cell_sum)
!$omp do schedule(dynamic,5)
     do ie = 1, mesh%ne  ! for each edge e
        !length_of_triangle = rearth*mesh%edt_leng(ie)
        do ilev = 1, nlev
           pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie) ! real PV
           cell_sum(1:mesh%edp_nedge(ie))   = mesh%edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie),ie)*&
                                              scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                                (pv_at_edge(1:mesh%edp_nedge(ie))+scalar_pv_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie),ie)))/2._r8
           scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/(rearth*mesh%edt_leng(ie)) ! real nct
        end do
     end do
!$omp end do nowait
!$omp end parallel

   END IF
!
! pe conserving scheme, not default
!
   case('pe')

     do ie = 1, mesh%ne  ! for each edge e
        !length_of_triangle = rearth*mesh%edt_leng(ie)
        do ilev = 1, nlev
           pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie)
           cell_sum(1:mesh%edp_nedge(ie))   = mesh%edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
                                              scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                                              pv_at_edge(1:mesh%edp_nedge(ie))
           scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/(rearth*mesh%edt_leng(ie))
        end do
     end do

   case default
     if(mpi_rank().eq.0) print*,"conserve_scheme in swe_para not set, must set to 'te' or 'pe'" 
   end select

   return
  end subroutine calc_coriolis_term_opt

  end module grist_dycore_hori_swe_module_2d
