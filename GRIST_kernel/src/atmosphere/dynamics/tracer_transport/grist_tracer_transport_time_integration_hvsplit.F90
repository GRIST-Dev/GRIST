
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: The Time driver for Tracer Transport of All Tracer Species
!
! Revision history:
!    (1) 2021,Dec, changed data structure  within vert adv and some code clean; 
!        added IEVA as in WS2020. When no IEVA is inactive, still regress old 
!        solutions
!
! References:
!            1. Zhang et al. (2020), MWR, doi:10.1175/MWR-D-19-0305.1
!            2. Li and Zhang (2022), MAP, doi:10.1007/s00703-022-00895-5
!----------------------------------------------------------------------------

  module grist_tracer_transport_time_integration_hvsplit

    use grist_constants,       only: i4, r8, rearth, one, two, zero
    use grist_domain_types,    only: global_domain
    use grist_data_types,      only: scalar_1d_field, scalar_2d_field, scalar_3d_field, exchange_field_list_3d
    use grist_nml_module,      only: ntracer, nlev, nlevp, nsteps, do_tracer_laplacian_2nd, &
                                     tracer_hadv_flag, tracer_vadv_flag, tracer_hori_limiter,  tracer_vert_limiter, tracer_hori_timestep, tracer_vert_timestep, &
                                     nrk_hori, nrk_vert, write_verbose, use_phys, ptend_tracer_f1_on, &
                                     isolate_tracer_test, tracer_hvsplit_method,  tracer_mxrt_fixer, tracer_vert_aimp, tracer_vert_fimp
! dycore data
! vertical transport
    use grist_hpe_continuity,     only: calc_hpe_tend_continuity_2d
! tracer specific operators 
    use grist_tracer_transport_vars_module,     only: tracerVarCellFull, tracerVarCellFace, tracerVarEdgeFull ! for all
    use grist_tracer_transport_tspas_module,    only: tracer_transport_vert_tspas_flux
    use grist_tracer_transport_flux_operators,  only: tracer_transport_normal_flux_at_edge, tracer_transport_vertic_flux_at_face
    use grist_tracer_transport_limiter_module,  only: tracer_transport_hori_flux_limiter_fct, tracer_transport_vert_flux_limiter_fct
    use grist_tracer_transport_utils_module,    only: tracer_transport_check_mxrt, tracer_transport_fixer_mxrt, tracer_transport_qneg3_mxrt
    use grist_tracer_transport_prescribe_module,only: tracer_transport_prescribe_ambient
! physics
    use grist_physics_data_structure, only: ptend_f1
#ifdef AMIPW_PHYSICS
    use grist_PhysW_nml_module,     only: step_cu
#endif
! MPI
    use grist_mpi
#ifndef SEQ_GRIST
    use grist_config_partition,       only: debug_data_2d, exchange_data_3d_add, exchange_data_3d
    use grist_clocks,                 only: clock_id, clock_begin, clock_end
#endif
#ifdef LAM_DOMAIN
    use grist_datam_glam_data_module, only: update_assignValueArea_3d
#endif

    implicit none
    private
    public :: tracer_transport_time_integration_init,&
              tracer_transport_time_integration_run ,&
              tracer_transport_time_integration_final
! local varitable for template io of subroutine
!
! variables used within RK step
!
     !real(r8), allocatable              :: div_sum(:,:)   ! for main/module array, compiler says it must be allocatable
     type(scalar_3d_field)               :: scalar_tracer_mxrt_at_pc_full_level_rk
     type(scalar_3d_field)               :: scalar_normal_mxrt_at_edge_full_level_rk
     type(scalar_3d_field)               :: scalar_normal_mxrt_at_edge_full_level_lo
     type(scalar_3d_field)               :: scalar_normal_mxrt_at_edge_full_level_lm
!
     real(r8)                            :: current_time, rk_substep
     logical,  allocatable               :: do_limiter_hori(:)
     logical,  allocatable               :: do_limiter_vert(:)
     logical,  allocatable               :: do_aimp_vert(:)
     logical                             :: prescribe_tracer_transport
     integer(i4)                         :: ihori_step, ivert_step, ns_hori, ns_vert, rk_number, iblock
     integer(i4)                         :: irk_step, iv, ilev, itracer, inb, edge_index
     integer(i4)                         :: clock_tracerh, clock_tracerv
     type(exchange_field_list_3d),pointer:: field_head_3d

   CONTAINS

!-----------------------------------------------------------------------
! name convention:
! for prognostic variables
! q(n+1) = q(n)+F(rk)
! n+1, new state after each sub stage
! n, state at begining of each time step
! rk itermediate state used for explicit tendency computin
!
! for other variable, _n denotes an intermediate state
! np1 denotes next (if needed)
!-----------------------------------------------------------------------

   subroutine tracer_transport_time_integration_run(mesh, dtime, itimestep,testcase,itstep,tstep_in_mstep)
!
! io
!
     use omp_lib
     type(global_domain),  intent(inout) :: mesh
     real(r8)          ,   intent(in)    :: dtime    ! time step for tracer module 
     integer(i4)       ,   intent(in)    :: itimestep
     character(len=*)  ,   intent(in)    :: testcase
     integer(i4), optional,intent(in)    :: itstep
     integer(i4), optional,intent(in)    :: tstep_in_mstep
! local
     real(r8), dimension(nlev,mesh%nv_full) :: scalar_qv_at_pc_full_old0
     real(r8), dimension(nlev,mesh%nv_full) :: scalar_qv_at_pc_full_old
     real(r8), dimension(mesh%nv_full)      :: scalar_template_a

        field_head_3d=>null()
        iblock = mpi_rank()

#ifdef AMIPW_PHYSICS
! added for cumulus tdk
        IF(present(itstep))then
        if(itstep.eq.1.and.mod(itimestep-1,step_cu).eq.0)then
           tracerVarCellFull%tend_qv_n%f  = zero
           scalar_qv_at_pc_full_old0      = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:)
        end if
        END If
        scalar_qv_at_pc_full_old  = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:)
#endif

!-----------------------------------------------------------------------
! only limit at the full-timestep step
!-----------------------------------------------------------------------

        do_limiter_hori(:)         = .false.
        do_limiter_vert(:)         = .false.
        do_aimp_vert(:)            = .false.

        do_limiter_hori(nrk_hori)  = tracer_hori_limiter
        do_limiter_vert(nrk_vert)  = tracer_vert_limiter
        do_aimp_vert(nrk_vert)     = tracer_vert_aimp

#ifndef SEQ_GRIST
        clock_tracerh = clock_id('tracerh')
        clock_tracerv = clock_id('tracerv')
#endif
        if(tracer_hori_limiter.or.tracer_hadv_flag.eq.99)then
           mesh%nv = mesh%nv_full
           mesh%ne = mesh%ne_halo(2)
#ifdef USE_HALO2
           mesh%ne = mesh%ne_halo(1)
#endif
        end if

!-----------------------------------------------------------------------
! prescribe wind (if-any)
!-----------------------------------------------------------------------

        if(isolate_tracer_test)then
           current_time = (itimestep-1)*dtime
           call tracer_transport_prescribe_ambient(mesh, current_time, testcase)
           if(iblock .eq. 0.and.write_verbose) print*,"prescribing tracer ambient",itimestep
        end if

        if(tracer_hori_limiter.or.tracer_hadv_flag.eq.99) then
           mesh%nv = mesh%nv_compute
           mesh%ne = mesh%ne_compute
        end if

        call reset_mass_continuity_equation

        SELECT CASE(tracer_hvsplit_method)

!-----------------------------------------------------------------------
! H-V, currently default 
!-----------------------------------------------------------------------

        case(1) ! H-V
#ifndef SEQ_GRIST
        call clock_begin(clock_tracerh)
#endif

!-----------------------------------------------------------------------
! H prepared, assign intermediate value from n-level
!-----------------------------------------------------------------------

           scalar_tracer_mxrt_at_pc_full_level_rk   = tracerVarCellFull%scalar_tracer_mxrt_n

!-----------------------------------------------------------------------
! horizontal call and update
!-----------------------------------------------------------------------

           call tracer_transport_time_integration_hori(mesh,dtime)
           if(iblock .eq. 0.and.write_verbose) print*,"Tracer Transport: finish updating horizontal tracer tendency"

!-----------------------------------------------------------------------
! Time-split here (due to limiter), prepared new n-level for vert
!-----------------------------------------------------------------------

           tracerVarCellFull%scalar_tracer_mass_n   = tracerVarCellFull%scalar_tracer_mass_np1
           tracerVarCellFull%scalar_tracer_mxrt_n   = tracerVarCellFull%scalar_tracer_mxrt_np1

#ifdef GRIST_TRACER_DEBUG
! optional mxrt fix
           if(tracer_mxrt_fixer)then
              call tracer_transport_fixer_mxrt(mesh%nv_full, tracerVarCellFull%scalar_tracer_mxrt_n,tracerVarCellFull%scalar_delhp_end_adv) 
              do itracer = 1, ntracer
                 tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)*tracerVarCellFull%scalar_delhp_end_adv%f
              end do
           end if
           call tracer_transport_check_mxrt(mesh, tracerVarCellFull%scalar_tracer_mxrt_n," after_tracer_transport_inte_hori")
#endif

#ifndef SEQ_GRIST
        call clock_end(clock_tracerh)
        call clock_begin(clock_tracerv)
#endif

!-----------------------------------------------------------------------
! V prepared, assign intermediate value from new n or np1
!-----------------------------------------------------------------------

           !scalar_tracer_mxrt_at_pc_full_level_rk  = tracerVarCellFull%scalar_tracer_mxrt_np1
           scalar_tracer_mxrt_at_pc_full_level_rk  = tracerVarCellFull%scalar_tracer_mxrt_n ! same as above w/o GRIST_TRACER_DEBUG

!-----------------------------------------------------------------------
! vertical call and update
!-----------------------------------------------------------------------

           !call tracer_transport_time_integration_vert(mesh,dtime)
           call tracer_transport_time_integration_vert_aimp(mesh,dtime) ! new added as default, can regress above when fully explicit
           if(iblock .eq. 0.and.write_verbose) print*,"tracer transport: finish updating vertical tracer tendency"

!-----------------------------------------------------------------------
! final update for the next adv step
!-----------------------------------------------------------------------

           tracerVarCellFull%scalar_tracer_mass_n%f = tracerVarCellFull%scalar_tracer_mass_np1%f
           tracerVarCellFull%scalar_tracer_mxrt_n%f = tracerVarCellFull%scalar_tracer_mxrt_np1%f

#ifdef GRIST_TRACER_DEBUG
! optional mxrt fix, default not use
           if(tracer_mxrt_fixer)then
              call tracer_transport_fixer_mxrt(mesh%nv_full, tracerVarCellFull%scalar_tracer_mxrt_n,tracerVarCellFull%scalar_delhp_end_adv)
              do itracer = 1, ntracer
                 tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)*tracerVarCellFull%scalar_delhp_end_adv%f
              end do
           end if
           call tracer_transport_check_mxrt(mesh, tracerVarCellFull%scalar_tracer_mxrt_n," after_tracer_transport_inte_vert")
#endif

#ifndef SEQ_GRIST
        call clock_end(clock_tracerv)
#endif

#ifdef TRACER_VH

!-----------------------------------------------------------------------
! V-H, not tuned for operational use
! in future, we may swap (1) and (2) for different tracer steps, mimicing 
! the strang splitting but across steps
!-----------------------------------------------------------------------

        case(2)

!-----------------------------------------------------------------------
! V-prepared
!-----------------------------------------------------------------------

           scalar_tracer_mxrt_at_pc_full_level_rk   = tracerVarCellFull%scalar_tracer_mxrt_n
!-----------------------------------------------------------------------
! call vert adv and update
!-----------------------------------------------------------------------
           call tracer_transport_time_integration_vert_aimp(mesh,dtime)
           if(iblock .eq. 0.and.write_verbose) print*,"Tracer Transport: finish updating vertical tracer tendency"

!-----------------------------------------------------------------------
! time-split here, prepared new n-level for horizontal advection
!-----------------------------------------------------------------------

           tracerVarCellFull%scalar_tracer_mass_n   = tracerVarCellFull%scalar_tracer_mass_np1
           tracerVarCellFull%scalar_tracer_mxrt_n   = tracerVarCellFull%scalar_tracer_mxrt_np1

!-----------------------------------------------------------------------
! intermediate value prepared
!-----------------------------------------------------------------------

           scalar_tracer_mxrt_at_pc_full_level_rk  = tracerVarCellFull%scalar_tracer_mxrt_np1

!-----------------------------------------------------------------------
! call hori and update inside
!-----------------------------------------------------------------------

           call tracer_transport_time_integration_hori(mesh,dtime)
           if(iblock .eq. 0.and.write_verbose) print*,"tracer transport: finish updating horizontal tracer tendency"

!-----------------------------------------------------------------------
! final update for the next adv step
!-----------------------------------------------------------------------
           tracerVarCellFull%scalar_tracer_mass_n%f = tracerVarCellFull%scalar_tracer_mass_np1%f
           tracerVarCellFull%scalar_tracer_mxrt_n%f = tracerVarCellFull%scalar_tracer_mxrt_np1%f

#endif
        case default
           if(iblock .eq. 0) print*,"tracer transport: you must select a hvsplit method (1 or 2) in namelist"
           stop
        end select

#ifndef SEQ_GRIST
        if(itimestep .eq. int(nsteps) .and. .false.) then
           call debug_data_2d(1,nlev,mesh%v_index,mesh%nv,"y1y",tracerVarCellFull%scalar_tracer_mass_n%f(1,:,:))
           call debug_data_2d(1,nlev,mesh%v_index,mesh%nv,"y2y",tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:))
        end if
#endif


!-----------------------------------------------------------------------
! typically we dono tracer Laplacian unless special case (e.g.,DCMIP-SC)
!-----------------------------------------------------------------------

        if(do_tracer_laplacian_2nd)then
           tracerVarCellFull%scalar_tracer_mass_n%f = tracerVarCellFull%scalar_tracer_mass_n%f+dtime*&
                                                      tracerVarCellFull%tend_tracer_mass_laplacian_2nd%f
#ifndef SEQ_GRIST
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_n)
           call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f
           end do
        end if

#ifdef AMIPW_PHYSICS
        IF(present(itstep))then
        if(itstep.eq.tstep_in_mstep.and.mod(itimestep,step_cu).eq.0)then
!-----------------------------------------------------------------------
! only use the tend at final tracer step and final cu step as adv forcing, in
! contrast of thften in dycore, we have used some assumptions here, 
! sensitivity not big
!-----------------------------------------------------------------------
           tracerVarCellFull%tend_qv_n%f = (tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:)-scalar_qv_at_pc_full_old)/dtime
        end if
        END IF
#endif

!-----------------------------------------------------------------------
!                          Coupling of F1 physics forcing
! After each tracer step, _n is the updated state, on which model physics reside 
! add f1's moisture forcing here;
! see comments in dycore's ptend_f1
!-----------------------------------------------------------------------

#ifdef USE_PTENDF1
        if(use_phys.and.ptend_tracer_f1_on)then
           tracerVarCellFull%scalar_tracer_mass_n%f       = tracerVarCellFull%scalar_tracer_mass_n%f+&
                                                           dtime*ptend_f1%tend_tracer_mass_at_pc_full_level%f
#ifndef SEQ_GRIST
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_n)
           call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)  = tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:)/&
                                                                      tracerVarCellFull%scalar_delhp_end_adv%f
           end do

           call tracer_transport_qneg3_mxrt(mesh%nv_full, nlev, ntracer,tracerVarCellFull%scalar_tracer_mxrt_n%f,"after physics coupling in F1")
           call tracer_transport_check_mxrt(mesh,                       tracerVarCellFull%scalar_tracer_mxrt_n  ,"after physics coupling in F1")
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)*tracerVarCellFull%scalar_delhp_end_adv%f
           end do
        end if
#endif

      return
   contains

     subroutine reset_mass_continuity_equation

!-------------------------------------------------------------------------------
! Assume the normal mass flux remains unchanged during the tracer transport,
! for full-level continuity eq in the large adv step,
! evaluate hori mass flux div {GRAD.(DELP*V)} at each full level
! then compute face level vertical mass flux (m*etadot)
!-------------------------------------------------------------------------------

!$omp parallel  private(iv,inb,edge_index,ilev)
!$omp do schedule(static,50)
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
! if not, problems should be fixed
        do iv = 1, mesh%nv_full
#else
        do iv = 1, mesh%nv
#endif
           tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,iv) = zero
           do inb = 1, mesh%vtx_nnb(iv)
              edge_index  = mesh%vtx_ed(inb,iv)
              do ilev = 1, nlev
                 tracerVarCellFull%tend_mass_div_avg_adv%f(ilev,iv) = tracerVarCellFull%tend_mass_div_avg_adv%f(ilev,iv)+&
                                                   tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,edge_index)*&
                                                   mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(edge_index)
              end do
           end do
           tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,iv) = tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,iv)/((rearth**2)*mesh%plg_areag(iv))
        end do
!$omp end do nowait
!$omp end parallel 
        if(iblock .eq. 0.and.write_verbose) print*, "Tracer Transport: finish evaluating horizontal mass flux div"

!-----------------------------------------------------------------------
! obtain etadot*m based on the continuity equation
! if isolate_tracer_test, eta mass flux will be prescribed
!-----------------------------------------------------------------------

        IF(.not.isolate_tracer_test)then
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
       !    do iv = 1, mesh%nv_full
            mesh%nv=mesh%nv_full
            call calc_hpe_tend_continuity_2d(mesh%nv_full, nlev, tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,1:mesh%nv), &  ! in
                                             scalar_template_a(1:mesh%nv)      , &  ! out, ps
                                             tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,1:mesh%nv))   ! out
            mesh%nv=mesh%nv_compute
#else
       !    do iv = 1, mesh%nv
            call calc_hpe_tend_continuity_2d(mesh%nv, nlev, tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,1:mesh%nv), &  ! in
                                             scalar_template_a(1:mesh%nv)      , &  ! out, ps
                                             tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,1:mesh%nv))   ! out
#endif

           !end do
           if(iblock .eq. 0.and.write_verbose) print*, "Tracer Transport: finish evaluating vertical mass flux"
        END IF
        return
     end subroutine reset_mass_continuity_equation 

   end subroutine tracer_transport_time_integration_run

   subroutine tracer_transport_time_integration_hori(mesh,dtime)

     type(global_domain),  intent(inout) :: mesh
     real(r8)          ,   intent(in)    :: dtime

        ns_hori  = int(dtime/tracer_hori_timestep)
        DO ihori_step = 1, ns_hori
        DO irk_step = 1, nrk_hori

           rk_number  = nrk_hori+1-irk_step
           rk_substep = tracer_hori_timestep/real(rk_number,r8)
!
! compute horizontal tracer flux and tracer mass flux div
!
           if (do_limiter_hori(irk_step) .and. (tracer_hadv_flag .lt. 8 .or. tracer_hadv_flag .gt. 12) .and.&
                                                tracer_hadv_flag .ne.58) &
#ifdef USE_HALO2
mesh%ne = mesh%ne_compute
#else
mesh%ne = mesh%ne_halo(1)
#endif
           call tracer_transport_normal_flux_at_edge(mesh, tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f   ,& ! in
                                                           tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f  ,&
                                                           tracerVarCellFull%scalar_tracer_mass_n%f             ,&
                                                           scalar_tracer_mxrt_at_pc_full_level_rk%f             ,& ! in
                                                           scalar_normal_mxrt_at_edge_full_level_rk%f           ,& ! out, still mixing ratio
                                                           tracer_hadv_flag             ,&
                                                           rk_substep,nlev,ntracer)
!
! Set to tracer mass
!
           do itracer = 1, ntracer
              scalar_normal_mxrt_at_edge_full_level_rk%f(itracer,:,:) = scalar_normal_mxrt_at_edge_full_level_rk%f(itracer,:,:)*&
                                                                        tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f
           end do

           if (do_limiter_hori(irk_step) .and. (tracer_hadv_flag .lt. 8 .or. tracer_hadv_flag .gt. 12) .and.&
                                                tracer_hadv_flag.ne.58) mesh%ne = mesh%ne_compute
#ifdef USE_HALO2
           call exchange_data_3d_add(mesh,field_head_3d,scalar_normal_mxrt_at_edge_full_level_rk)
           call exchange_data_3d(mesh%local_block,field_head_3d)
#else
! note: halo2 will not use this tracer_adv_flag that need 3 halos
           if (do_limiter_hori(irk_step) .and. (tracer_hadv_flag .ge. 8 .and. tracer_hadv_flag .le. 12.or.&
                                                tracer_hadv_flag.eq.58)) then
#ifndef SEQ_GRIST
              call exchange_data_3d_add(mesh,field_head_3d,scalar_normal_mxrt_at_edge_full_level_rk)
              call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
           end if

#endif
!
! limiting mixing ratio
!
           IF(do_limiter_hori(irk_step))THEN
! 1st-order upwind flux based on q(n), need halo(2)
#ifdef USE_HALO2
              mesh%ne = mesh%ne_halo(1)
#else
              mesh%ne = mesh%ne_halo(2)
#endif
              call tracer_transport_normal_flux_at_edge(mesh, tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f  ,& ! in
                                                              tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f ,&
                                                              tracerVarCellFull%scalar_tracer_mass_n%f            ,&
                                                              tracerVarCellFull%scalar_tracer_mxrt_n%f            ,& ! in
                                                              scalar_normal_mxrt_at_edge_full_level_lo%f          ,& ! out, still mixing ratio
                                                              1                         ,&
                                                              rk_substep,nlev,ntracer)

#ifdef USE_HALO2
              call exchange_data_3d_add(mesh,field_head_3d,scalar_normal_mxrt_at_edge_full_level_lo)
              call exchange_data_3d(mesh%local_block,field_head_3d)
#endif

! unchanged delhp*V*q(n) as low-order q-mass flux
              do itracer = 1, ntracer
                 scalar_normal_mxrt_at_edge_full_level_lo%f(itracer,:,:) = scalar_normal_mxrt_at_edge_full_level_lo%f(itracer,:,:)*&
                                                                           tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f
              end do

              mesh%ne = mesh%ne_compute
              call tracer_transport_hori_flux_limiter_fct(mesh, tracerVarCellFull%scalar_tracer_mxrt_n%f   ,& !       q at time n
                                                                tracerVarCellFull%scalar_tracer_mass_n%f   ,& ! delhp*q at time n
                                                                tracerVarCellFull%scalar_delhp_end_adv%f   ,&
                                                                scalar_normal_mxrt_at_edge_full_level_rk%f,& ! high-order mass flux*q
                                                                scalar_normal_mxrt_at_edge_full_level_lo%f,& !  low-order mass flux*q
                                                                scalar_normal_mxrt_at_edge_full_level_lm%f,& !    limited mass flux*q
                                                                rk_substep, ntracer, nlev)

               scalar_normal_mxrt_at_edge_full_level_rk%f = scalar_normal_mxrt_at_edge_full_level_lm%f ! reset rk flux by limited one

           END IF
!
! DIV operator
!
!$omp parallel  private(iv,inb,edge_index,ilev,itracer) 
!$omp do schedule(dynamic,50)
! compute region
           do iv = 1, mesh%nv
              tracerVarCellFull%tend_tracer_mass_hori%f(1:ntracer,1:nlev,iv) = zero
              do inb = 1, mesh%vtx_nnb(iv)
                 edge_index  = mesh%vtx_ed(inb,iv)
                 do ilev = 1, nlev
                    do itracer = 1, ntracer
                       tracerVarCellFull%tend_tracer_mass_hori%f(itracer,ilev,iv) = tracerVarCellFull%tend_tracer_mass_hori%f(itracer,ilev,iv)+& 
                                               scalar_normal_mxrt_at_edge_full_level_rk%f(itracer,ilev,edge_index)*&
                                               mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(edge_index)
                    end do
                 end do
              end do
              tracerVarCellFull%tend_tracer_mass_hori%f(1:ntracer,1:nlev,iv) = -tracerVarCellFull%tend_tracer_mass_hori%f(1:ntracer,1:nlev,iv)/((rearth**2)*mesh%plg_areag(iv))
           end do
!$omp end do nowait
!$omp end parallel 

!
! tracer mass horizontal update
!
           tracerVarCellFull%scalar_tracer_mass_np1%f = tracerVarCellFull%scalar_tracer_mass_n%f+rk_substep*tracerVarCellFull%tend_tracer_mass_hori%f

! 
! exchange data
! Note that: this exchange is only necessary when nrk_hori.gt.1. .and.
! ns_hori.gt.1 E.g., ffsl and tspas do not need this (can do in vert_adv), when hori adv is not cycled within a full tracer transport step;
! but still code here for general purpose; one data exchange is acceptable
!
#ifndef SEQ_GRIST
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_np1)
           ! no need: call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mxrt_np1)
           call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mxrt_np1%f(itracer,1:nlev,1:mesh%nv_full) = tracerVarCellFull%scalar_tracer_mass_np1%f(itracer,1:nlev,1:mesh%nv_full)/&
                                                                                     tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,1:mesh%nv_full)
           end do

#ifdef LAM_DOMAIN
           call update_assignValueArea_3d(mesh,nlev,ntracer,"tracer",tracerVarCellFull%scalar_delhp_end_adv%f,&
                                                                     tracerVarCellFull%scalar_tracer_mass_np1,&
                                                                     tracerVarCellFull%scalar_tracer_mxrt_np1)
#endif
! prog
           scalar_tracer_mxrt_at_pc_full_level_rk    = tracerVarCellFull%scalar_tracer_mxrt_np1

        END DO
        END DO

   end subroutine tracer_transport_time_integration_hori

   subroutine tracer_transport_time_integration_vert(mesh,dtime)
     use omp_lib
! io
     type(global_domain),  intent(inout) :: mesh
     real(r8)          ,   intent(in)    :: dtime
! local
     real(r8) :: scalar_mxrt_flux_1d_nlevp_ho(nlevp)
     real(r8) :: scalar_mxrt_flux_1d_nlevp_lm(nlevp)
     real(r8) :: scalar_mxrt_flux_3d_nlevp_ho(ntracer,nlevp,mesh%nv_full)
     real(r8) :: scalar_mxrt_flux_3d_nlevp_lm(ntracer,nlevp,mesh%nv_full)
 

        ns_vert  = dtime/tracer_vert_timestep
        SELECT CASE(tracer_vadv_flag)
        CASE(1,2,3,9)
        DO ivert_step = 1, ns_vert
        DO irk_step = 1, nrk_vert

           rk_number  = nrk_vert+1-irk_step
           rk_substep = tracer_vert_timestep/rk_number
!
! compute vertical tracer mass flux tendency delta(m*etadot*q)
!
!$omp parallel  private(iv,itracer,scalar_mxrt_flux_1d_nlevp_lm,scalar_mxrt_flux_1d_nlevp_ho) 
!$omp do schedule(dynamic,50)
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
           do iv = 1, mesh%nv_full
#else
           do iv = 1, mesh%nv
#endif
              do itracer = 1, ntracer
                  call tracer_transport_vertic_flux_at_face(scalar_tracer_mxrt_at_pc_full_level_rk%f(itracer,1:nlev,iv) ,&  ! full level q
                                                            tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,iv),&  ! m*etadot
                                                            tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,iv) ,&  ! delhp full
                                                            scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlevp,iv),&  ! face level q
                                                            rk_substep                ,&
                                                            tracer_vadv_flag          ,&
                                                            nlev, nlevp)
              end do ! tracer
           end do  ! ncell
!$omp end do nowait
!$omp end parallel

!   scalar_mxrt_flux_1d_nlevp_lm = scalar_template_1d_nlevp_b
           IF(do_limiter_vert(irk_step))then  ! vertical one-dimensional FCT limiter
           scalar_mxrt_flux_3d_nlevp_ho = scalar_mxrt_flux_3d_nlevp_lm
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
           do iv = 1, mesh%nv_full
#else
           do iv = 1, mesh%nv
#endif
              do itracer = 1, ntracer
                    call tracer_transport_vert_flux_limiter_fct(rk_substep, &
                                                                tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,iv) ,&  ! m*etadot
                                                                tracerVarCellFull%scalar_tracer_mass_n%f(itracer,1:nlev ,iv) ,&  ! tracer mass n
                                                                tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,1:nlev ,iv) ,&  ! tracer mxrt n
                                                                tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,iv)          ,&  ! delhp np1
                                                                scalar_mxrt_flux_3d_nlevp_ho(itracer,1:nlevp,iv),&  ! ho flux rk,  q
                                                                scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlevp,iv),&  ! lm flux np1, q
                                                                nlev, nlevp)
              end do
            end do ! ncell
            END IF

            do itracer = 1, ntracer
               tracerVarCellFull%tend_tracer_mass_vert%f(itracer,1:nlev,:) = &
               -(tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(2:nlevp,:)*scalar_mxrt_flux_3d_nlevp_lm(itracer,2:nlevp,:)-&
                 tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlev,:) *scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlev,:))
            end do

!
! tracer mass vertical update
!
           tracerVarCellFull%scalar_tracer_mass_np1%f = tracerVarCellFull%scalar_tracer_mass_n%f+rk_substep*tracerVarCellFull%tend_tracer_mass_vert%f
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mxrt_np1%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mass_np1%f(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f
           end do
#ifndef SEQ_GRIST
#ifndef TRACER_VERTNOEXC
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_np1)
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mxrt_np1)
           call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
#endif
! prog
           scalar_tracer_mxrt_at_pc_full_level_rk    = tracerVarCellFull%scalar_tracer_mxrt_np1

        END DO
        END DO

        CASE(99)  ! TSPAS, should not be used currently

        DO ivert_step = 1, ns_vert

           rk_substep      = tracer_vert_timestep
! compute vertical tracer mass flux tendency delta(m*etadot*q)
           do iv = 1, mesh%nv
              do itracer = 1, ntracer
                 call tracer_transport_vert_tspas_flux(rk_substep                       , &
                                                       tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,iv)  , &
                                                       tracerVarCellFull%scalar_tracer_mass_n%f(itracer,1:nlev,iv)   , &
                                                       tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,1:nlev,iv)   , &
                                                       tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,iv)   , &
                                                       scalar_mxrt_flux_1d_nlevp_lm, & ! limited q*delhp
                                                       nlev, nlevp)
                 tracerVarCellFull%tend_tracer_mass_vert%f(itracer,1:nlev,iv) = &
                                      -(tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(2:nlev+1,iv)*scalar_mxrt_flux_1d_nlevp_lm(2:nlev+1)-&
                                        tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlev,iv)  *scalar_mxrt_flux_1d_nlevp_lm(1:nlev))/tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,iv)
              end do
           end do
!
! tracer mass vertical update
!
           tracerVarCellFull%scalar_tracer_mass_np1%f = tracerVarCellFull%scalar_tracer_mass_n%f+rk_substep*tracerVarCellFull%tend_tracer_mass_vert%f
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mxrt_np1%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mass_np1%f(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f
           end do
!
! exchange data
!
#ifndef SEQ_GRIST
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_np1)
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mxrt_np1)
           call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
! prog
           scalar_tracer_mxrt_at_pc_full_level_rk    = tracerVarCellFull%scalar_tracer_mxrt_np1
        END DO

        case default
           if(iblock .eq. 0) print*,"tracer transport: tracer_vadv_flag configured incorrectly, please modify"
        end select

   end subroutine tracer_transport_time_integration_vert

   subroutine tracer_transport_time_integration_init(mesh)
    use grist_dycore_vars_module,                only: dycoreVarCellFull
    use grist_tracer_transport_ffsl_module,      only: tracer_transport_ffsl_flux_sm10_init
    use grist_nml_module,                        only: tracer_transport_hori_scheme
! io
     type(global_domain),  intent(inout) :: mesh

        if(.not.allocated(do_limiter_hori)) allocate(do_limiter_hori(nrk_hori))
        if(.not.allocated(do_limiter_vert)) allocate(do_limiter_vert(nrk_vert))
        if(.not.allocated(do_aimp_vert))    allocate(do_aimp_vert(nrk_vert))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_rk%f)) allocate(scalar_normal_mxrt_at_edge_full_level_rk%f(ntracer,nlev,mesh%ne_full))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_lo%f)) allocate(scalar_normal_mxrt_at_edge_full_level_lo%f(ntracer,nlev,mesh%ne_full))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_lm%f)) allocate(scalar_normal_mxrt_at_edge_full_level_lm%f(ntracer,nlev,mesh%ne_full))
!
! when model starts, tracer mix ratio is given by initial state,
! scalar_delhp_at_pc_full_level_beg_adv is given by dyn
! adv_init only evaluates the initial state of tracer mass
! this tracer mass will be integrated since after,
! no need for init any more
!
        scalar_normal_mxrt_at_edge_full_level_rk%f = 0._r8
        scalar_normal_mxrt_at_edge_full_level_lo%f = 0._r8
        scalar_normal_mxrt_at_edge_full_level_lm%f = 0._r8
        tracerVarCellFull%scalar_delhp_end_adv%f = dycoreVarCellFull%scalar_delhp_n%f

        do itracer = 1, ntracer
           tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)*&
                                                                   dycoreVarCellFull%scalar_delhp_n%f
        end do

      if(trim(tracer_transport_hori_scheme).eq.'ffsl')then
#ifdef USE_HALO2
      mesh%ne = mesh%ne_halo(1)
#else
      mesh%ne = mesh%ne_halo(2)
#endif
      call tracer_transport_ffsl_flux_sm10_init(mesh)
      mesh%ne = mesh%ne_full

      end if

      return
   end subroutine tracer_transport_time_integration_init

   subroutine tracer_transport_time_integration_final
      use grist_tracer_transport_ffsl_module, only: tracer_transport_ffsl_flux_sm10_final
      use grist_nml_module,                   only: tracer_transport_hori_scheme

      deallocate(do_limiter_hori); deallocate(do_limiter_vert); deallocate(do_aimp_vert)
      deallocate(scalar_normal_mxrt_at_edge_full_level_rk%f)
      deallocate(scalar_normal_mxrt_at_edge_full_level_lo%f)
      deallocate(scalar_normal_mxrt_at_edge_full_level_lm%f)

      if(trim(tracer_transport_hori_scheme).eq.'ffsl') call tracer_transport_ffsl_flux_sm10_final

   end subroutine tracer_transport_time_integration_final

!--------------------------------------------------------
! BELOW added for optional IEVA advection as in WS2020
!--------------------------------------------------------

   subroutine tracer_transport_aimp_blend_mass_flux(mesh,dtime,nlev,nlevp,ncell_all, ncell_cal, nedge_all, &
                                                    scalar_eta_mass_flux_at_pc_face        , &
                                                    scalar_delhp_at_pc_full_level          , &
                                                    scalar_normal_velocity_at_pc_full_level, &
                                                    scalar_eta_mass_flux_at_pc_face_exp    , &
                                                    scalar_eta_mass_flux_at_pc_face_imp)
! io
    type(global_domain), intent(in) :: mesh
    real(r8),            intent(in) :: dtime
    integer(i4), intent(in)  :: nlev, nlevp, ncell_all, ncell_cal, nEdge_all
    real(r8),  intent(in)    :: scalar_eta_mass_flux_at_pc_face(nlevp,nCell_all)
    real(r8),  intent(in)    :: scalar_delhp_at_pc_full_level(nlev,nCell_all)
    real(r8),  intent(in)    :: scalar_normal_velocity_at_pc_full_level(nlev,nEdge_all)
    real(r8),  intent(inout) :: scalar_eta_mass_flux_at_pc_face_exp(nlevp,nCell_all) ! we want
    real(r8),  intent(inout) :: scalar_eta_mass_flux_at_pc_face_imp(nlevp,nCell_all) ! we want
!local
    real(r8), parameter :: kesi=0.9_r8, alpha_max = 1.1_r8, alpha_min = 0.8_r8
    integer(i4) :: iv, ie, inb, ilev, mark
    real(r8)    :: outflow_div_at_pc_full_level(nlev,ncell_all)
    real(r8)    :: alpha_face, alpha_h, alpha_star_max, alpha_star_min, ggg

!
! Evaluate alpha(h) at face, which means the sum of all outflow divergence,
! at the upwind level of etadot. Because alpha(h) is for g, g is with omega, 
! so only calc 2->nlev face level
!
      do iv = 1, ncell_cal
         outflow_div_at_pc_full_level(1:nlev,iv)  = zero
         do inb = 1, mesh%vtx_nnb(iv)
            ie      = mesh%vtx_ed(inb,iv)
            do ilev = 1, nlev
               outflow_div_at_pc_full_level(ilev,iv) = outflow_div_at_pc_full_level(ilev,iv) + &
               max(scalar_normal_velocity_at_pc_full_level(ilev,ie)*mesh%plg_nr(inb,iv)*mesh%edp_leng(ie),zero)
            end do
          end do
          outflow_div_at_pc_full_level(1:nlev,iv) = outflow_div_at_pc_full_level(1:nlev,iv) / (mesh%plg_areag(iv)*rearth)
      end do

      do iv = 1, ncell_cal
         do ilev = 2, nlev
            alpha_face = dtime*abs(scalar_eta_mass_flux_at_pc_face(ilev,iv)*two/(scalar_delhp_at_pc_full_level(ilev,iv)+scalar_delhp_at_pc_full_level(ilev-1,iv)))
          ! if eta_mass<0, upward motion for this face level, upwind-full is ilev
          ! if eta_mass>0, downwa motion for this face level, upwind-full is ilev-1
          ! there are not-a-few places in GRIST that selects upwind in this way, make sure single-precision does not harm it!
            mark       = int(0.5_r8*sign(one,scalar_eta_mass_flux_at_pc_face(ilev,iv))+0.5_r8)
            alpha_h    = dtime*outflow_div_at_pc_full_level(ilev-mark,iv)
            alpha_star_max = (one-kesi*alpha_h/alpha_max)
            alpha_star_min = alpha_min*alpha_star_max/alpha_max

          ! can we avoid if here?
            if(alpha_face.le.alpha_star_min)then
               ggg = one
            else if(alpha_face.gt.alpha_star_min.and.alpha_face.le.(2*alpha_star_max-alpha_star_min))then
               ggg = one/(one+((alpha_face-alpha_star_min)**2)/(4*alpha_star_max*(alpha_star_max-alpha_star_min)))
            else
               ggg = alpha_star_max/alpha_face
            end if
            if(tracer_vert_fimp) ggg = zero
            scalar_eta_mass_flux_at_pc_face_exp(ilev,iv) =      ggg *scalar_eta_mass_flux_at_pc_face(ilev,iv)
            scalar_eta_mass_flux_at_pc_face_imp(ilev,iv) = (one-ggg)*scalar_eta_mass_flux_at_pc_face(ilev,iv)
         end do
         scalar_eta_mass_flux_at_pc_face_exp(1,iv) = zero; scalar_eta_mass_flux_at_pc_face_exp(nlevp,iv) = zero
         scalar_eta_mass_flux_at_pc_face_imp(1,iv) = zero; scalar_eta_mass_flux_at_pc_face_imp(nlevp,iv) = zero
      end do

      return
   end subroutine tracer_transport_aimp_blend_mass_flux

!----------------------------------------------------------------------------
! The below part is copied from xxx_vert routine with modifications to
! allow doing aimp at selected stage of nrk_vert, can be 1, 2, or more
! dafault at the final stage; if no aimp is used at a stage, still regression
! as before; verified in idealized and real-world setup
!----------------------------------------------------------------------------

   subroutine tracer_transport_time_integration_vert_aimp(mesh,dtime)
     use omp_lib
     use grist_tridiagonal_solver, only: solve_tridiagonal
! io
     type(global_domain),  intent(inout) :: mesh
     real(r8)          ,   intent(in)    :: dtime
! local
     real(r8) :: scalar_mxrt_flux_1d_nlevp_ho(nlevp)
     real(r8) :: scalar_mxrt_flux_1d_nlevp_lm(nlevp)
     real(r8) :: scalar_mxrt_flux_3d_nlevp_ho(ntracer,nlevp,mesh%nv_full)
     real(r8) :: scalar_mxrt_flux_3d_nlevp_lm(ntracer,nlevp,mesh%nv_full)
     real(r8) :: scalar_eta_mass_flux_at_pc_face_exp (nlevp,mesh%nv_full)
     real(r8) :: scalar_eta_mass_flux_at_pc_face_imp (nlevp,mesh%nv_full)
     real(r8) :: trid_aaa(nlev),trid_bbb(nlev),trid_ccc(nlev),trid_rrr(nlev),trid_out(nlev)

        ns_vert  = dtime/tracer_vert_timestep

        SELECT CASE(tracer_vadv_flag)
        CASE(1,2,3,9)
        DO ivert_step = 1, ns_vert
        DO irk_step = 1, nrk_vert

           rk_number  = nrk_vert+1-irk_step
           rk_substep = tracer_vert_timestep/rk_number
!
! default full explicit method
!
           scalar_eta_mass_flux_at_pc_face_exp = tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f

!----------------------------------------------------
! ADD IEVA HERE; if IEVA not called, this routine
! must produce identical solutions as before
!----------------------------------------------------

           IF(do_aimp_vert(irk_step))then
#ifdef TRACER_VERTNOEXC
              call tracer_transport_aimp_blend_mass_flux(mesh,rk_substep,nlev,nlevp,mesh%nv_full, mesh%nv_full, mesh%ne_full, &
#else
              call tracer_transport_aimp_blend_mass_flux(mesh,rk_substep,nlev,nlevp,mesh%nv_full, mesh%nv, mesh%ne_full, &
#endif
                                                    tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,1:mesh%nv_full)   , &
                                                    tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,1:mesh%nv_full)            , &
                                                    tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f(1:nlev,1:mesh%ne_full), &
                                                    scalar_eta_mass_flux_at_pc_face_exp , & ! need
                                                    scalar_eta_mass_flux_at_pc_face_imp)    ! need
           end if

!----------------------------------------------------
! END IEVA HERE; if IEVA not called, this routine
! must produce identical solutions as before
!----------------------------------------------------

!
! compute vertical tracer mass flux tendency delta(m*etadot*q)
!
!$omp parallel  private(iv,itracer,scalar_mxrt_flux_1d_nlevp_lm,scalar_mxrt_flux_1d_nlevp_ho) 
!$omp do schedule(dynamic,50)
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
           do iv = 1, mesh%nv_full
#else
           do iv = 1, mesh%nv
#endif
              do itracer = 1, ntracer
                  call tracer_transport_vertic_flux_at_face(scalar_tracer_mxrt_at_pc_full_level_rk%f(itracer,1:nlev,iv) ,&  ! full level q
                                                            scalar_eta_mass_flux_at_pc_face_exp(1:nlevp,iv)    ,&  ! m*etadot
                                                            tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,iv),&  ! delhp full
                                                            scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlevp,iv)   ,&  ! face level q
                                                            rk_substep                ,&
                                                            tracer_vadv_flag          ,&
                                                            nlev, nlevp)
              end do ! tracer
           end do  ! ncell
!$omp end do nowait
!$omp end parallel

           IF(do_limiter_vert(irk_step))then  ! vertical one-dimensional FCT limiter
              scalar_mxrt_flux_3d_nlevp_ho = scalar_mxrt_flux_3d_nlevp_lm
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
           do iv = 1, mesh%nv_full
#else
           do iv = 1, mesh%nv
#endif
              do itracer = 1, ntracer
                    call tracer_transport_vert_flux_limiter_fct(rk_substep, &
                                                                scalar_eta_mass_flux_at_pc_face_exp(1:nlevp,iv) ,&  ! m*etadot
                                                                tracerVarCellFull%scalar_tracer_mass_n%f(itracer,1:nlev ,iv) ,&  ! tracer mass n
                                                                tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,1:nlev ,iv) ,&  ! tracer mxrt n
                                                                tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,iv)          ,&  ! delhp np1
                                                                scalar_mxrt_flux_3d_nlevp_ho(itracer,1:nlevp,iv),&  ! ho flux rk,  q
                                                                scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlevp,iv),&  ! lm flux np1, q
                                                                nlev, nlevp)
              end do
            end do ! ncell
            END IF

            do itracer = 1, ntracer
               tracerVarCellFull%tend_tracer_mass_vert%f(itracer,1:nlev,:) = &
               -(scalar_eta_mass_flux_at_pc_face_exp(2:nlevp,:)*scalar_mxrt_flux_3d_nlevp_lm(itracer,2:nlevp,:)-&
                 scalar_eta_mass_flux_at_pc_face_exp(1:nlev,:) *scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlev,:))
            end do
!
! tracer mass vertical update
!
           tracerVarCellFull%scalar_tracer_mass_np1%f = tracerVarCellFull%scalar_tracer_mass_n%f+rk_substep*tracerVarCellFull%tend_tracer_mass_vert%f
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mxrt_np1%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mass_np1%f(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f
           end do

!----------------------------------------------------
! ADD IEVA HERE; if IEVA not called, this routine
! must produce identical solutions as before
!----------------------------------------------------

           IF(do_aimp_vert(irk_step))then
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
           DO iv = 1, mesh%nv_full
#else
           do iv = 1, mesh%nv
#endif
              trid_aaa(1) = zero; trid_ccc(nlev) = zero;
              do ilev = 2, nlev
                 trid_aaa(ilev)   = -rk_substep*scalar_eta_mass_flux_at_pc_face_imp(ilev,iv)*(one+sign(one,scalar_eta_mass_flux_at_pc_face_imp(ilev,iv)))
                 trid_ccc(ilev-1) =  rk_substep*scalar_eta_mass_flux_at_pc_face_imp(ilev,iv)*(one-sign(one,scalar_eta_mass_flux_at_pc_face_imp(ilev,iv)))
              end do
              do ilev = 1, nlev
                 trid_bbb(ilev)  = rk_substep*scalar_eta_mass_flux_at_pc_face_imp(ilev+1,iv)*(one+sign(one,scalar_eta_mass_flux_at_pc_face_imp(ilev+1,iv)))-& 
                                   rk_substep*scalar_eta_mass_flux_at_pc_face_imp(ilev,iv)  *(one-sign(one,scalar_eta_mass_flux_at_pc_face_imp(ilev,iv)))+&
                 two*tracerVarCellFull%scalar_delhp_end_adv%f(ilev,iv)
              end do
              do itracer = 1, ntracer
                 trid_rrr(1:nlev) = two*tracerVarCellFull%scalar_tracer_mass_np1%f(itracer,1:nlev,iv)
                 call solve_tridiagonal(trid_aaa,trid_bbb,trid_ccc,trid_rrr,trid_out,nlev)
                 tracerVarCellFull%scalar_tracer_mxrt_np1%f(itracer,1:nlev,iv) = trid_out(1:nlev)
                 tracerVarCellFull%scalar_tracer_mass_np1%f(itracer,1:nlev,iv) = tracerVarCellFull%scalar_tracer_mxrt_np1%f(itracer,1:nlev,iv)*&
                                                                                tracerVarCellFull%scalar_delhp_end_adv%f(1:nlev,iv)
              end do
           END DO
           END IF

!----------------------------------------------------
! END IEVA HERE; if IEVA not called, this routine
! must produce identical solutions as before
!----------------------------------------------------

#ifndef SEQ_GRIST
#ifndef TRACER_VERTNOEXC
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_np1)
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mxrt_np1)
           call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
#endif

#ifdef LAM_DOMAIN
           call update_assignValueArea_3d(mesh,nlev,ntracer,"tracer",tracerVarCellFull%scalar_delhp_end_adv%f,&
                                                                     tracerVarCellFull%scalar_tracer_mass_np1,&
                                                                     tracerVarCellFull%scalar_tracer_mxrt_np1)
#endif
! prog
           scalar_tracer_mxrt_at_pc_full_level_rk    = tracerVarCellFull%scalar_tracer_mxrt_np1

        END DO
        END DO

        case default
           if(iblock .eq. 0) print*,"tracer transport: tracer_vadv_flag configured incorrectly, please modify"
        end select

   end subroutine tracer_transport_time_integration_vert_aimp

  end module grist_tracer_transport_time_integration_hvsplit
