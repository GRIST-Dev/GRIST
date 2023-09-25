module grist_dycore_hori_swe_module_2d_mixed

    use grist_constants,        only: i4, r8, rearth, omega, rearth_r4 => rearth_ns, omega_r4 => omega_ns
    use grist_constants,        only: r4 => ns
    use grist_data_types,       only: scalar_2d_field, exchange_field_list_2d
    use grist_domain_types,     only: global_domain
    use grist_nml_module,       only: pv_order, mass_wtd_vor, conserve_scheme, ke_method, ke_alpha, stencil_width
    use grist_mpi,              only: mpi_rank
#ifdef USE_HALO2
    use grist_config_partition, only: exchange_data_2d_add ,  exchange_data_2d, exchange_data_2d_r4
#endif
!data
    use grist_dycore_module_vars_control_mixed, only: edt_leng, edp_leng, edp_trsk_edpl_on_edge, tri_c_lat, tri_kite_area

    use grist_clocks,                 only: clock_id, clock_begin, clock_end

    implicit none

    integer(i4)      :: clock_nct_prep, clock_nct_vor, clock_nct_mas, clock_nct_cor, clock_nct_flux, clock_nct_mpi, clock_nct_fnl

contains

!==============================================
!   Calculate grad of kinetic energy at edge
!==============================================

  subroutine calc_grad_kinetic_energy(mesh, &
                                      scalar_normal_velocity_at_edge_full_level,&
                                      tend_grad_ke_at_edge_full_level, &
                                      nlev)
! io
   use omp_lib
   type(global_domain),   intent(inout) :: mesh
   real(r4), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r4), allocatable, intent(inout) :: tend_grad_ke_at_edge_full_level(:,:)
   integer(i4)          , intent(in)    :: nlev
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
            length_of_triangle_edge = rearth*mesh%edt_leng(index_edge)
            length_of_voronoi_edge  = rearth*mesh%edp_leng(index_edge)
            do ilev = 1, nlev
               wind                 = scalar_normal_velocity_at_edge_full_level(ilev,index_edge)
               tmp_sum(ilev)        = tmp_sum(ilev)+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
            end do
         end do
         kinetic_energy(:,iv) = tmp_sum(:)/((rearth**2)*mesh%plg_areag(iv))
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
              length_of_triangle_edge = rearth*mesh%edt_leng(index_edge)
              length_of_voronoi_edge  = rearth*mesh%edp_leng(index_edge)
              do ilev = 1, nlev
                 wind                 = scalar_normal_velocity_at_edge_full_level(ilev,index_edge)
                 tmp_sum(ilev)        = tmp_sum(ilev)+0.25_r8*length_of_voronoi_edge*length_of_triangle_edge*wind**2
              end do
           end do
           kinetic_energy(:,iv) = tmp_sum(:)/((rearth**2)*mesh%plg_areag(iv))
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
         tend_grad_ke_at_edge_full_level(ilev,ie) = -mesh%edt_edpNr_edtTg(ie)*(kinetic_energy(ilev,v2)-kinetic_energy(ilev,v1))/(rearth_r4*edt_leng(ie))
         ! print *, mesh%edt_edpNr_edtTg(ie), kinetic_energy(ilev,v2), rearth, mesh%edt_leng(ie)
      end do
   end do
   ! stop
!$omp end do nowait
!$omp end parallel 

!   deallocate(kinetic_energy)
   return
  end subroutine calc_grad_kinetic_energy

  subroutine calc_tend_nct_at_edge_full_level(mesh,&
                                              scalar_height_at_pc_full_level   ,&
                                              scalar_normal_velocity_at_edge_full_level,&
                                              scalar_normal_flux_at_edge_full_level    ,&
                                              tend_nct_at_edge_full_level        ,&
                                              dtime, nlev,irk_step)
!io
   type(global_domain),   intent(inout) :: mesh
   type(scalar_2d_field), intent(in)    :: scalar_height_at_pc_full_level
   type(scalar_2d_field), intent(in)    :: scalar_normal_velocity_at_edge_full_level
   type(scalar_2d_field), intent(in)    :: scalar_normal_flux_at_edge_full_level
   type(scalar_2d_field), intent(inout) :: tend_nct_at_edge_full_level
   real(r8)             , intent(in)    :: dtime
   integer(i4)          , intent(in)    :: nlev, irk_step
!local-exchange
   type(scalar_2d_field)                :: scalar_pv_at_dc_full_level
   type(scalar_2d_field)                :: scalar_pv_at_edge_full_level
!local-nonexchange
   real(r4)                             :: scalar_tangen_velocity_at_edge_full_level(nlev,1:mesh%ne_full)
   real(r8)                             :: scalar_height_at_dc_full_level(nlev,1:mesh%nt_full)
   real(r4)                             :: scalar_avor_at_dc_full_level(  nlev,1:mesh%nt_full)
   real(r4)                             :: scalar_rvor_at_dc_full_level(  nlev,1:mesh%nt_full)
   real(r4)                             :: scalar_coriolis_term_at_edge(  nlev,1:mesh%ne_full)   ! what we need
   real(r4)                             ::                          tmp(  nlev,1:mesh%ne_full)
   integer(i4)                          :: it, ie, iv
   real(r4)                             :: sum_tmp(nlev)
   ! real(r8)                             :: length_of_edge
   integer(i4)                          :: ilev, index_of_edge, t_edge_v, cell_index
#ifdef USE_HALO2
   type(exchange_field_list_2d),pointer :: field_head_2d

   field_head_2d=>null()
#endif

   call clock_begin(clock_nct_prep)
   if(.not.allocated(scalar_pv_at_edge_full_level%f_r4))    allocate(scalar_pv_at_edge_full_level%f_r4(nlev,1:mesh%ne_full))
   if(.not.allocated(scalar_pv_at_dc_full_level%f_r4))      allocate(scalar_pv_at_dc_full_level%f_r4(nlev,1:mesh%nt_full))

   scalar_avor_at_dc_full_level       = 0._r4
   scalar_height_at_dc_full_level     = 1._r4
   scalar_coriolis_term_at_edge       = 0._r4

!
! calculate vorticity at dual cell using Stokes' Theorem
! only until stencil_width-1 because calc_scalar_at_dual_cell needs one-more
! layer
!
     mesh%ne = mesh%ne_halo(stencil_width-1)
     mesh%nt = mesh%nt_halo(stencil_width-1)
     call clock_end(clock_nct_prep)
! inline_opt_here
!     call calc_vorticity_at_dual_cell(mesh, scalar_normal_velocity_at_edge_full_level%f ,&
!                                            scalar_avor_at_dc_full_level                ,&
!                                            scalar_rvor_at_dc_full_level                ,&
!                                            nlev)
      call clock_begin(clock_nct_vor)
!$omp parallel  private(it,sum_tmp,ie,ilev,index_of_edge,t_edge_v)
!$omp do schedule(dynamic,5)
      do it = 1, mesh%nt
         sum_tmp(:) = 0._r4
         do ie = 1, mesh%tri_nnb(it)
            index_of_edge          = mesh%tri_ed(ie,it)
            t_edge_v               = mesh%edt_edpTg_edtNr(index_of_edge)*mesh%tri_nr(ie,it)
            do ilev = 1, nlev
               sum_tmp(ilev)       = sum_tmp(ilev)+t_edge_v*rearth_r4*edt_leng(index_of_edge)*scalar_normal_velocity_at_edge_full_level%f_r4(ilev,index_of_edge)
            end do
         end do
         scalar_rvor_at_dc_full_level(1:nlev,it) = sum_tmp(1:nlev)/((rearth**2)*mesh%tri_areag(it))
         scalar_avor_at_dc_full_level(1:nlev,it) = scalar_rvor_at_dc_full_level(1:nlev,it)+2._r4*omega_r4*sin(tri_c_lat(it))
      end do
!$omp end do nowait
!$omp end parallel 
      call clock_end(clock_nct_vor)
! inline_opt_here

     call clock_begin(clock_nct_mas)
     if(mass_wtd_vor)then
! inline_opt_here
      !  call calc_scalar_at_dual_cell(mesh,scalar_height_at_pc_full_level%f, &
      !                                     scalar_height_at_dc_full_level  , &
      !                                     nlev)
!$omp parallel  private(it,sum_tmp,ie,ilev,cell_index) 
!$omp do schedule(dynamic,5)
      do it = 1, mesh%nt
         sum_tmp(:) = 0._r4
         do ie = 1, mesh%tri_nnb(it)
            cell_index       = mesh%tri_v(ie,it)
            do ilev = 1, nlev
               sum_tmp(ilev) = sum_tmp(ilev)+(tri_kite_area(ie,it))*scalar_height_at_pc_full_level%f_r4(ilev,cell_index)
            end do
         end do
         scalar_height_at_dc_full_level(1:nlev,it) = sum_tmp(1:nlev)/(mesh%tri_areag(it))
      end do
!$omp end do nowait
!$omp end parallel
! inline_opt_here inline
        scalar_pv_at_dc_full_level%f_r4 = scalar_avor_at_dc_full_level/scalar_height_at_dc_full_level
     else
        scalar_pv_at_dc_full_level%f_r4 = scalar_avor_at_dc_full_level
     end if
     call clock_end(clock_nct_mas)
!
! reconstruct tangent wind using nct operator, from stencil_width to stencil_width-1
!
     if(pv_order(irk_step).eq.3)then
#ifdef USE_HALO2
        call exchange_data_2d_add(mesh,field_head_2d,scalar_pv_at_dc_full_level,"r4")
        call exchange_data_2d_r4(mesh%local_block,field_head_2d)
#endif
        call clock_begin(clock_nct_cor)
        tmp(:,:) = 1._r4
        call calc_coriolis_term1(mesh,tmp  ,scalar_normal_velocity_at_edge_full_level%f_r4, &
                                            scalar_tangen_velocity_at_edge_full_level, &
                                            nlev)
        call clock_end(clock_nct_cor)
! pv-dc only having until stencil_width-1 layer, so exchange data for PV3
        call clock_begin(clock_nct_mpi)
        call clock_end(clock_nct_mpi)
     end if

!
! tangent flux along hx's edge
! if pv2 or pv3-halo3, we can just compute PV-edge until ne_halo(1) without exch
!
     mesh%ne = mesh%ne_halo(1)
     mesh%nt = mesh%nt_halo(1)
! old code compute until ne_compute for pv3 and use_halo2, which is overly conservative
#if (defined PV3_OLDCODE)
#ifdef USE_HALO2
     if(pv_order(irk_step).eq.3)then
        mesh%ne = mesh%ne_compute
        mesh%nt = mesh%nt_compute
     end if
#endif
#endif

! we need pv flux at ne_halo(1)
! for pv2, we can get it from nt_halo(1)
! for pv3, we can only get ne_compute
     call clock_begin(clock_nct_flux)
     call calc_normal_flux_tr_edge(mesh, scalar_pv_at_dc_full_level%f_r4         ,&
                                         scalar_tangen_velocity_at_edge_full_level ,&
                                         scalar_normal_velocity_at_edge_full_level%f_r4 ,&
                                         scalar_pv_at_edge_full_level%f_r4              ,&
                                         .false.                        ,&
                                         pv_order(irk_step), dtime, nlev)
     call clock_end(clock_nct_flux)

     call clock_begin(clock_nct_mpi)
! old code compute until ne_compute for pv3 and use_halo2, which is overly conservative
#if (defined PV3_OLDCODE)
     if(pv_order(irk_step).eq.3)then
!exchange to get ne_halo(1), if halo3, no need exchange
#ifdef USE_HALO2
        call exchange_data_2d_add(mesh,field_head_2d,scalar_pv_at_edge_full_level,"r4")
        call exchange_data_2d_r4(mesh%local_block,field_head_2d)
#endif
     end if
#endif
     call clock_end(clock_nct_mpi)

!================================================
!       The Final nonlinear Coriolis tendency
!================================================

!
! from ne_halo(1) to ne_compute
! imp: set this to ne is necessary for next KE-grad evaluation
!
     mesh%ne = mesh%ne_compute
     mesh%nt = mesh%nt_compute

     call clock_begin(clock_nct_cor)
     if(mass_wtd_vor)then
        call calc_coriolis_term(mesh,scalar_pv_at_edge_full_level%f_r4,&
                                     scalar_normal_flux_at_edge_full_level%f_r4,&
                                     scalar_coriolis_term_at_edge, &
                                     nlev)

     else
      !   call calc_coriolis_term(mesh,scalar_pv_at_edge_full_level%f_r4,&
      !                                scalar_normal_velocity_at_edge_full_level%f,&
      !                                scalar_coriolis_term_at_edge, &
      !                                nlev)
      print *, "MIXCODE ERROR! mass_wtd_vor must be set .true."
      stop
     end if
     call clock_end(clock_nct_cor)

     call clock_begin(clock_nct_fnl)
! send back
      tend_nct_at_edge_full_level%f_r4  = -1._r4*scalar_coriolis_term_at_edge

! clean 
      deallocate(scalar_pv_at_edge_full_level%f_r4)
      deallocate(scalar_pv_at_dc_full_level%f_r4)
      call clock_end(clock_nct_fnl)

    return
  end subroutine calc_tend_nct_at_edge_full_level

  subroutine calc_coriolis_term(mesh,scalar_pv_at_edge_full_level           ,&
                                     scalar_normal_flux_at_edge_full_level  ,&
                                     scalar_coriolis_term_at_edge_full_level,&
                                     nlev)

! io 
   use omp_lib
   type(global_domain),     intent(in)    :: mesh
   real(r4),                intent(in)    :: scalar_pv_at_edge_full_level(nlev,mesh%ne_full)
   real(r4), allocatable  , intent(in)    :: scalar_normal_flux_at_edge_full_level(:,:)
   real(r4),                intent(inout) :: scalar_coriolis_term_at_edge_full_level(nlev,mesh%ne_full)
   integer(i4)            , intent(in)    :: nlev
! local
   real(r4)                               :: length_of_triangle
   real(r4)                               :: pv_at_edge(16)
   real(r4)                               :: cell_sum(16)
   integer(i4)                            :: ie, ilev 

!
! te conserving scheme
!
   select case(trim(conserve_scheme)) 
   case('te')
   !  do ie = 1, mesh%ne  ! for each edge e
   !     length_of_triangle = rearth*mesh%edt(ie)%leng
   !     do ilev = 1, nlev

   !        pv_at_edge(1:mesh%edp(ie)%nedge) = scalar_pv_at_edge_full_level(ilev,ie)
   !        cell_sum(1:mesh%edp(ie)%nedge)   = mesh%edp(ie)%trsk_on_edge(1:mesh%edp(ie)%nedge)*&
   !                                           mesh%edp(ie)%edpl_on_edge(1:mesh%edp(ie)%nedge)*&
   !scalar_normal_flux_at_edge_full_level(ilev,mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge))*&
   !                             (pv_at_edge(1:mesh%edp(ie)%nedge)+scalar_pv_at_edge_full_level(ilev,mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge)))/2._r8
   !        scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp(ie)%nedge))/length_of_triangle

   !     end do
   !  end do
!$omp parallel  private(ie,length_of_triangle,ilev,pv_at_edge,cell_sum) 
!$omp do schedule(dynamic,5)
     do ie = 1, mesh%ne  ! for each edge e
        length_of_triangle = rearth_r4*edt_leng(ie)
        do ilev = 1, nlev

           pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie)
           cell_sum(1:mesh%edp_nedge(ie))   = edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
   scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                                (pv_at_edge(1:mesh%edp_nedge(ie))+scalar_pv_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie)))/2._r4
           scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
        end do
     end do
!$omp end do nowait
!$omp end parallel 

!
! pe conserving scheme
!
   ! case('pe')
   ! !  do ie = 1, mesh%ne  ! for each edge e
   ! !     length_of_triangle = rearth*mesh%edt(ie)%leng
   ! !     do ilev = 1, nlev
 
   ! !        pv_at_edge(1:mesh%edp(ie)%nedge) = scalar_pv_at_edge_full_level(ilev,ie)
   ! !        cell_sum(1:mesh%edp(ie)%nedge)   = mesh%edp(ie)%trsk_on_edge(1:mesh%edp(ie)%nedge)*&
   ! !                                           mesh%edp(ie)%edpl_on_edge(1:mesh%edp(ie)%nedge)*&
   ! !scalar_normal_flux_at_edge_full_level(ilev,mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge))*&
   ! !                              pv_at_edge(1:mesh%edp(ie)%nedge)
   ! !        scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp(ie)%nedge))/length_of_triangle
   ! !     end do
   ! !  end do

   !   do ie = 1, mesh%ne  ! for each edge e
   !      length_of_triangle = rearth*mesh%edt_leng(ie)
   !      do ilev = 1, nlev

   !         pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie)
   !         cell_sum(1:mesh%edp_nedge(ie))   = mesh%edp_trsk_on_edge(1:mesh%edp_nedge(ie), ie)*&
   !                                            mesh%edp_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
   ! scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
   !                               pv_at_edge(1:mesh%edp_nedge(ie))
   !         scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
   !      end do
   !   end do
   case default
     if(mpi_rank().eq.0) print*,"conserve_scheme in swe_para not set, must set to 'te' or 'pe'" 
   end select

   return
  end subroutine calc_coriolis_term

  subroutine calc_coriolis_term1(mesh,scalar_pv_at_edge_full_level           ,&
                                      scalar_normal_flux_at_edge_full_level  ,&
                                      scalar_coriolis_term_at_edge_full_level,&
                                      nlev)

! io 
   use omp_lib
   type(global_domain),     intent(in)    :: mesh
   real(r4),                intent(in)    :: scalar_pv_at_edge_full_level(nlev,mesh%ne_full)
   real(r4), allocatable  , intent(in)    :: scalar_normal_flux_at_edge_full_level(:,:)
   real(r4),                intent(inout) :: scalar_coriolis_term_at_edge_full_level(nlev,mesh%ne_full)
   integer(i4)            , intent(in)    :: nlev
! local
   real(r4)                               :: length_of_triangle
   real(r4)                               :: pv_at_edge(16)
   real(r4)                               :: cell_sum(16)
   integer(i4)                            :: ie, ilev 

!
! te conserving scheme
!
   select case(trim(conserve_scheme)) 
   case('te')
   !  do ie = 1, mesh%ne  ! for each edge e
   !     length_of_triangle = rearth*mesh%edt(ie)%leng
   !     do ilev = 1, nlev

   !        pv_at_edge(1:mesh%edp(ie)%nedge) = scalar_pv_at_edge_full_level(ilev,ie)
   !        cell_sum(1:mesh%edp(ie)%nedge)   = mesh%edp(ie)%trsk_on_edge(1:mesh%edp(ie)%nedge)*&
   !                                           mesh%edp(ie)%edpl_on_edge(1:mesh%edp(ie)%nedge)*&
   !scalar_normal_flux_at_edge_full_level(ilev,mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge))*&
   !                             (pv_at_edge(1:mesh%edp(ie)%nedge)+scalar_pv_at_edge_full_level(ilev,mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge)))/2._r8
   !        scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp(ie)%nedge))/length_of_triangle

   !     end do
   !  end do
!$omp parallel  private(ie,length_of_triangle,ilev,pv_at_edge,cell_sum) 
!$omp do schedule(dynamic,5)
     do ie = 1, mesh%ne  ! for each edge e
        length_of_triangle = rearth_r4*edt_leng(ie)
        do ilev = 1, nlev

           pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie)
           cell_sum(1:mesh%edp_nedge(ie))   = edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
   scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
                                (pv_at_edge(1:mesh%edp_nedge(ie))+scalar_pv_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie)))/2._r4
           scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle

        end do
     end do
!$omp end do nowait
!$omp end parallel 

!
! pe conserving scheme
!
   ! case('pe')
   ! !  do ie = 1, mesh%ne  ! for each edge e
   ! !     length_of_triangle = rearth*mesh%edt(ie)%leng
   ! !     do ilev = 1, nlev
 
   ! !        pv_at_edge(1:mesh%edp(ie)%nedge) = scalar_pv_at_edge_full_level(ilev,ie)
   ! !        cell_sum(1:mesh%edp(ie)%nedge)   = mesh%edp(ie)%trsk_on_edge(1:mesh%edp(ie)%nedge)*&
   ! !                                           mesh%edp(ie)%edpl_on_edge(1:mesh%edp(ie)%nedge)*&
   ! !scalar_normal_flux_at_edge_full_level(ilev,mesh%edp(ie)%edge_on_edge(1:mesh%edp(ie)%nedge))*&
   ! !                              pv_at_edge(1:mesh%edp(ie)%nedge)
   ! !        scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp(ie)%nedge))/length_of_triangle
   ! !     end do
   ! !  end do

   !   do ie = 1, mesh%ne  ! for each edge e
   !      length_of_triangle = rearth*mesh%edt_leng(ie)
   !      do ilev = 1, nlev

   !         pv_at_edge(1:mesh%edp_nedge(ie)) = scalar_pv_at_edge_full_level(ilev,ie)
   !         cell_sum(1:mesh%edp_nedge(ie))   = mesh%edp_trsk_on_edge(1:mesh%edp_nedge(ie), ie)*&
   !                                            mesh%edp_edpl_on_edge(1:mesh%edp_nedge(ie), ie)*&
   ! scalar_normal_flux_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie), ie))*&
   !                               pv_at_edge(1:mesh%edp_nedge(ie))
   !         scalar_coriolis_term_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
   !      end do
   !   end do
   case default
     if(mpi_rank().eq.0) print*,"conserve_scheme in swe_para not set, must set to 'te' or 'pe'" 
   end select

   return
  end subroutine calc_coriolis_term1

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
   real(r4), allocatable, intent(in)    :: scalar_pv_at_dc_full_level(:,:)
   real(r4),              intent(in)    :: scalar_normal_velocity_at_edge_full_level(nlev,mesh%ne_full)
   real(r4), allocatable, intent(in)    :: scalar_tangen_velocity_at_edge_full_level(:,:)
   real(r4), allocatable, intent(inout) :: scalar_pv_at_edge_full_level(:,:)
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
                scalar_pv_at_edge_full_level(ilev,ie) = 0.5*(scalar_pv_at_dc_full_level(ilev,v1)+&
                                                             scalar_pv_at_dc_full_level(ilev,v2))
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

           !    call flux_wrf3(scalar_pv_at_dc_full_level(ilev,v1),&
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
             !  call flux_wrf3(scalar_pv_at_dc_full_level(ilev,v1),&
             !                 scalar_pv_at_dc_full_level(ilev,v2),&
             !                 der0 ,&
             !                 scalar_pv_at_edge_full_level(ilev,ie), &
             !                 real(mesh%edp_leng(ie),r8))
! inline written here
               der0   = 2._r8*dot_product(mesh%edp_v_weight_dual_cell(1:mesh%tri_stencil_number_2nd(v_upwind), mark, ie), &
                              scalar_pv_at_dc_full_level(ilev,mesh%tri_stencil_index_2nd(1:mesh%tri_stencil_number_2nd(v_upwind),v_upwind))-&
                              scalar_pv_at_dc_full_level(ilev,v_upwind))
               scalar_pv_at_edge_full_level(ilev,ie) = (scalar_pv_at_dc_full_level(ilev,v1)+scalar_pv_at_dc_full_level(ilev,v2))*0.5_r8-&
                            ((edp_leng(ie)**2)/6._r8)*der0
! inline_opt_here

!           end do

!          call calc_der_2nd_at_tr(mesh,scalar_pv_at_dc_full_level,v2,ie,1,der1)
!================================================
!             Compute edge-based value 
!================================================
!          call flux_wrf34(wind_edge,scalar_pv_at_dc_full_level(v1),&
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

  subroutine grist_dycore_hori_swe_clock_init

     implicit none

#ifndef SEQ_GRIST
      clock_nct_prep     = clock_id('nct_prep')
      clock_nct_vor      = clock_id('nct_vor' )
      clock_nct_mas      = clock_id('nct_mas' )
      clock_nct_cor      = clock_id('nct_cor' )
      clock_nct_flux     = clock_id('nct_flux')
      clock_nct_mpi      = clock_id('nct_mpi' )
      clock_nct_fnl      = clock_id('nct_fnl' )
#endif     
     
  end subroutine grist_dycore_hori_swe_clock_init

end module grist_dycore_hori_swe_module_2d_mixed
