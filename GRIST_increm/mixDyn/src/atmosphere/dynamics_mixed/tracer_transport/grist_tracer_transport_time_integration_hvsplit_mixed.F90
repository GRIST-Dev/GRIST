
!----------------------------------------------------------------------------
! Created on Aug 11 2023
! Author: Siyuan Chen
! Version 1.0
! Description: Building upon the grist_tracer_transport_time_integration_hvsplit module, 
!              replace selected variables with single precision while maintaining similarity 
!              to the results of double precision variables.
!----------------------------------------------------------------------------
module grist_tracer_transport_time_integration_hvsplit_mixed

    use grist_constants,        only: r8, i4, zero => zero_ns, one => one_ns, two => two_ns, rearth_r8 => rearth, rearth => rearth_ns, zero_r4 => zero_ns
    use grist_constants,        only: r4 => ns
    use grist_data_types,       only: scalar_3d_field, exchange_field_list_3d
    use grist_domain_types,     only: global_domain
    use grist_mpi,              only: mpi_rank
    use grist_nml_module,       only: ntracer, nlev, nlevp, nrk_hori, nrk_vert, isolate_tracer_test, tracer_hori_limiter,   &
                                      tracer_vert_limiter, tracer_vert_aimp, tracer_hadv_flag, write_verbose,               &
                                      tracer_hvsplit_method, tracer_hori_timestep, tracer_mxrt_fixer, tracer_vert_fimp, tracer_transport_hori_scheme,    &
                                      tracer_vert_timestep, tracer_vadv_flag, do_tracer_laplacian_2nd, doNotDiagnose
    use grist_tracer_hpe_continuity_mixed,              only: calc_hpe_tend_continuity_2d
    use grist_tracer_transport_prescribe_module_mixed,  only: tracer_transport_prescribe_ambient
    use grist_tracer_transport_flux_operators_mixed,    only: tracer_transport_normal_flux_at_edge,     &
                                                              tracer_transport_normal_flux_at_edge_r8,  &
                                                              tracer_transport_vertic_flux_at_face
    use grist_tracer_transport_limiter_module_mixed,    only: tracer_transport_hori_flux_limiter_fct, tracer_transport_vert_flux_limiter_fct
    use grist_tracer_transport_utils_module_mixed,      only: tracer_transport_check_mxrt, tracer_transport_fixer_mxrt

#ifndef SEQ_GRIST
    use grist_config_partition,       only: exchange_data_3d_add, exchange_data_3d, exchange_data_3d_r4
    use grist_clocks,                 only: clock_id, clock_begin, clock_end
#endif

! data
    use grist_tracer_module_vars_control_mixed,         only: edp_leng, plg_areag
    use grist_tracer_transport_vars_module,             only: tracerVarCellFace, tracerVarCellFull, tracerVarEdgeFull, tracerVarSurface
#ifdef AMIPW_PHYSICS
    use grist_PhysW_nml_module,     only: step_cu
#endif


    implicit none

    type(scalar_3d_field)                   ::  scalar_tracer_mxrt_at_pc_full_level_rk
    type(scalar_3d_field)                   ::  scalar_normal_mxrt_at_edge_full_level_rk
    type(scalar_3d_field)                   ::  scalar_normal_mxrt_at_edge_full_level_lo
    type(scalar_3d_field)                   ::  scalar_normal_mxrt_at_edge_full_level_lm

    real(r8)                                ::  current_time, rk_substep
    ! real(r4)                                ::  scalar_template_a
    logical,  allocatable                   ::  do_limiter_hori(:)
    logical,  allocatable                   ::  do_limiter_vert(:)
    logical,  allocatable                   ::  do_aimp_vert(:)
    logical                                 ::  prescribe_tracer_transport
    integer(i4)                             ::  ihori_step, ivert_step, ns_hori, ns_vert, rk_number, iblock
    integer(i4)                             ::  irk_step, iv, ilev, itracer, inb, edge_index
    type(exchange_field_list_3d), pointer   ::  field_head_3d

#ifndef SEQ_GRIST
    integer(i4)                             ::  clock_tracerh, clock_tracerh_mpi, clock_tracerh_flux, clock_tracerh_limiter, clock_tracerh_div, clock_tracerh_update, &
                                                clock_tracerv
#endif

    public :: tracer_transport_time_integration_init, &
              tracer_transport_time_integration_run , &
              tracer_transport_time_integration_final
contains

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
        real(r4), dimension(nlev,mesh%nv_full) :: scalar_qv_at_pc_full_old0
        real(r4), dimension(nlev,mesh%nv_full) :: scalar_qv_at_pc_full_old
#ifdef DBL_MASS_TRACER
        real(r8)                               :: div_sum(nlev)
        real(r8), dimension(mesh%nv)           :: scalar_template_a
#else
        real(r4)                               :: div_sum(nlev)
        real(r4), dimension(mesh%nv)           :: scalar_template_a
#endif

        field_head_3d=>null()
        iblock = mpi_rank()

!================================================
! added for cumulus tdk
!================================================
#ifdef AMIPW_PHYSICS
        IF(present(itstep))then
!================================================
! added for cumulus tdk
!================================================
        if(itstep.eq.1.and.mod(itimestep-1,step_cu).eq.0)then
           tracerVarCellFull%tend_qv_n%f_r4 = zero
           scalar_qv_at_pc_full_old0        = tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(1,:,:)
        end if
        END If
#endif

!
! only limit at the full-timestep step
!
        prescribe_tracer_transport = isolate_tracer_test

        do_limiter_hori(:)         = .false.
        do_limiter_vert(:)         = .false.
        do_aimp_vert(:)            = .false.
        do_limiter_hori(nrk_hori)  = tracer_hori_limiter
        do_limiter_vert(nrk_vert)  = tracer_vert_limiter
        do_aimp_vert(nrk_vert)     = tracer_vert_aimp

!
! prescribe wind
!
        if(tracer_hori_limiter.or.tracer_hadv_flag.eq.99)then
            mesh%nv = mesh%nv_full
            mesh%ne = mesh%ne_halo(2)
#ifdef USE_HALO2
            mesh%ne = mesh%ne_halo(1)
#endif
        end if
        if(prescribe_tracer_transport)then
            current_time = (itimestep-1)*dtime
            call tracer_transport_prescribe_ambient(mesh, current_time, testcase)
            if(iblock .eq. 0.and.write_verbose) print*,"prescribing tracer ambient",itimestep
        end if
        if(tracer_hori_limiter.or.tracer_hadv_flag.eq.99) then
            mesh%nv = mesh%nv_compute
            mesh%ne = mesh%ne_compute
        end if

!-------------------------------------------------------------------------------
! Assume the normal mass flux remains unchanged during the tracer transport,
! for full level continuity eq in the large adv step,
! evaluate mass flux div {GRAD.(DELP*V)} at each full level
! then compute face level vertical mass flux (m*etadot)
!-------------------------------------------------------------------------------

!$omp parallel  private(iv,inb,edge_index,ilev) 
!$omp do schedule(static,50)
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
        do iv = 1, mesh%nv_full
#else
        do iv = 1, mesh%nv
#endif
            !div_sum(1,:) = 0._r8
#ifdef DBL_MASS_TRACER
            tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,iv) = zero
#else
            tracerVarCellFull%tend_mass_div_avg_adv%f_r4(1:nlev,iv) = zero_r4
#endif
            do inb = 1, mesh%vtx_nnb(iv)
                edge_index  = mesh%vtx_ed(inb,iv)
                do ilev = 1, nlev
                ! div_sum(1,ilev) = div_sum(1,ilev)+tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,edge_index)*&
#ifdef DBL_MASS_TRACER
                tracerVarCellFull%tend_mass_div_avg_adv%f(ilev,iv) = tracerVarCellFull%tend_mass_div_avg_adv%f(ilev,iv)+&
                                                   tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,edge_index)*&
                                                   mesh%plg_nr(inb,iv)*rearth_r8*mesh%edp_leng(edge_index)
#else
                tracerVarCellFull%tend_mass_div_avg_adv%f_r4(ilev,iv) = tracerVarCellFull%tend_mass_div_avg_adv%f_r4(ilev,iv)+&
                                                   tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,edge_index)*&
                                                   mesh%plg_nr(inb,iv)*rearth*edp_leng(edge_index)
#endif
                end do
            end do
            !tracerVarCellFull%tend_mass_div_avg_adv%f(:,iv) = div_sum(1,:)/((rearth**2)*mesh%plg_areag(iv))
#ifdef DBL_MASS_TRACER
            tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,iv) = tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,iv)/((rearth_r8**2)*mesh%plg_areag(iv))
#else
            tracerVarCellFull%tend_mass_div_avg_adv%f_r4(1:nlev,iv) = tracerVarCellFull%tend_mass_div_avg_adv%f_r4(1:nlev,iv)/((rearth**2)*plg_areag(iv))
#endif
        end do
!$omp end do nowait
!$omp end parallel 
        if(iblock .eq. 0.and.write_verbose) print*, "Tracer Transport: finish evaluating horizontal mass flux div"

!         IF(.not.isolate_tracer_test)then
! ! obtain etadot*m based on the continuity equation
! #ifdef TRACER_VERTNOEXC
! ! This will produce identical solutions as using exchange inside vertical advection
!               do iv = 1, mesh%nv_full
! #else
!               do iv = 1, mesh%nv
! #endif

!               !scalar_template_1d_nlev = tracerVarCellFull%tend_mass_div_avg_adv%f(:,iv)
!               call calc_hpe_tend_continuity(tracerVarCellFull%tend_mass_div_avg_adv%f_r4(1:nlev,iv), &  ! in
!                                             scalar_template_a      , &  ! out, ps
!                                             tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f_r4(1:nlevp,iv))   ! out
!               !tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(:,iv) = scalar_template_1d_nlevp

!           end do
!           if(iblock .eq. 0.and.write_verbose) print*, "Tracer Transport: finish evaluating vertical mass flux"
!         END IF
        IF(.not.isolate_tracer_test)then
! obtain etadot*m based on the continuity equation
#ifdef TRACER_VERTNOEXC
! This will produce identical solutions as using exchange inside vertical advection
       !    do iv = 1, mesh%nv_full
#ifdef DBL_MASS_TRACER
            call calc_hpe_tend_continuity_2d(mesh%nv_full, nlev, tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,1:mesh%nv), &  ! in
#else
            call calc_hpe_tend_continuity_2d(mesh%nv_full, nlev, tracerVarCellFull%tend_mass_div_avg_adv%f_r4(1:nlev,1:mesh%nv), &  ! in
#endif

#else
       !    do iv = 1, mesh%nv
#ifdef DBL_MASS_TRACER
            call calc_hpe_tend_continuity_2d(mesh%nv, nlev, tracerVarCellFull%tend_mass_div_avg_adv%f(1:nlev,1:mesh%nv), &  ! in
#else
            call calc_hpe_tend_continuity_2d(mesh%nv, nlev, tracerVarCellFull%tend_mass_div_avg_adv%f_r4(1:nlev,1:mesh%nv), &  ! in
#endif

#endif
                                               scalar_template_a(1:mesh%nv)      , &  ! out, ps
#ifdef DBL_MASS_TRACER
                                               tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,1:mesh%nv))   ! out
#else
                                               tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f_r4(1:nlevp,1:mesh%nv))   ! out
#endif
           !end do
           If(.not.doNotDiagnose)then
            tracerVarSurface%tend_hpressure%f = scalar_template_a
           End if
           if(iblock .eq. 0.and.write_verbose) print*, "Tracer Transport: finish evaluating vertical mass flux"
        END IF

        scalar_qv_at_pc_full_old  = tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(1,:,:)

        SELECT CASE(tracer_hvsplit_method)
        case(1) ! H-V
#ifndef SEQ_GRIST
            call clock_begin(clock_tracerh)
#endif
! H prepared
            ! scalar_tracer_mxrt_at_pc_full_level_rk   = tracerVarCellFull%scalar_tracer_mxrt_n
            scalar_tracer_mxrt_at_pc_full_level_rk%f_r4= tracerVarCellFull%scalar_tracer_mxrt_n%f_r4
            scalar_tracer_mxrt_at_pc_full_level_rk%pos = tracerVarCellFull%scalar_tracer_mxrt_n%pos
! call
            call tracer_transport_time_integration_hori(mesh,dtime)
            if(iblock .eq. 0.and.write_verbose) print*,"Tracer Transport: finish updating horizontal tracer tendency"
! time-split here, prepared for vertical advection
            ! tracerVarCellFull%scalar_tracer_mass_n   = tracerVarCellFull%scalar_tracer_mass_np1
            tracerVarCellFull%scalar_tracer_mass_n%f_r4   = tracerVarCellFull%scalar_tracer_mass_np1%f_r4
            tracerVarCellFull%scalar_tracer_mass_n%pos    = tracerVarCellFull%scalar_tracer_mass_np1%pos
            ! tracerVarCellFull%scalar_tracer_mxrt_n   = tracerVarCellFull%scalar_tracer_mxrt_np1
            tracerVarCellFull%scalar_tracer_mxrt_n%f_r4   = tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4
            tracerVarCellFull%scalar_tracer_mxrt_n%pos    = tracerVarCellFull%scalar_tracer_mxrt_np1%pos
 ! optional mxrt fix
            if(tracer_mxrt_fixer)then
                call tracer_transport_fixer_mxrt(mesh%nv_full, tracerVarCellFull%scalar_tracer_mxrt_n,tracerVarCellFull%scalar_delhp_end_adv) 
                do itracer = 1, ntracer
                    tracerVarCellFull%scalar_tracer_mass_n%f_r4(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(itracer,:,:)*tracerVarCellFull%scalar_delhp_end_adv%f_r4
                end do
            end if
            call tracer_transport_check_mxrt(mesh, tracerVarCellFull%scalar_tracer_mxrt_n," after_tracer_transport_inte_hori")
#ifndef SEQ_GRIST
            call clock_end(clock_tracerh)
            call clock_begin(clock_tracerv)
#endif
! V prepared
            ! scalar_tracer_mxrt_at_pc_full_level_rk  = tracerVarCellFull%scalar_tracer_mxrt_np1
            scalar_tracer_mxrt_at_pc_full_level_rk%f_r4  = tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4
            scalar_tracer_mxrt_at_pc_full_level_rk%pos   = tracerVarCellFull%scalar_tracer_mxrt_np1%pos
! call
            !call tracer_transport_time_integration_vert(mesh,dtime)
            call tracer_transport_time_integration_vert_aimp(mesh,dtime)
            if(iblock .eq. 0.and.write_verbose) print*,"tracer transport: finish updating vertical tracer tendency"
! final update for the next adv step
            tracerVarCellFull%scalar_tracer_mass_n%f_r4 = tracerVarCellFull%scalar_tracer_mass_np1%f_r4
            tracerVarCellFull%scalar_tracer_mxrt_n%f_r4 = tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4
! optional mxrt fix
            if(tracer_mxrt_fixer)then
                call tracer_transport_fixer_mxrt(mesh%nv_full, tracerVarCellFull%scalar_tracer_mxrt_n,tracerVarCellFull%scalar_delhp_end_adv)
                do itracer = 1, ntracer
                    tracerVarCellFull%scalar_tracer_mass_n%f_r4(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(itracer,:,:)*tracerVarCellFull%scalar_delhp_end_adv%f_r4
                end do
            end if
            call tracer_transport_check_mxrt(mesh, tracerVarCellFull%scalar_tracer_mxrt_n," after_tracer_transport_inte_vert")
#ifndef SEQ_GRIST
            call clock_end(clock_tracerv)
#endif
        case default
            if(iblock .eq. 0) print*,"tracer transport: you must select a hvsplit method in namelist"
            stop
        end select

        if(do_tracer_laplacian_2nd)then
           tracerVarCellFull%scalar_tracer_mass_n%f_r4 = tracerVarCellFull%scalar_tracer_mass_n%f_r4+dtime*&
                                                         tracerVarCellFull%tend_tracer_mass_laplacian_2nd%f
#ifndef SEQ_GRIST
           call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_n,"r4")
           call exchange_data_3d_r4(mesh%local_block,field_head_3d)
#endif
           do itracer = 1, ntracer
              tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(itracer,:,:) = &
              tracerVarCellFull%scalar_tracer_mass_n%f_r4(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f_r4
           end do
        end if

#ifdef AMIPW_PHYSICS
        IF(present(itstep))then
        if(itstep.eq.tstep_in_mstep.and.mod(itimestep,step_cu).eq.0)then
! only use the tend at final tracer step and final cu step as adv forcing, in
! contrast of thften in dycore, we have used some assumptions here, sensitivity
! not big
           tracerVarCellFull%tend_qv_n%f_r4 = (tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(1,:,:)-scalar_qv_at_pc_full_old)/dtime
        end if
        END IF
#endif

      return
    end subroutine tracer_transport_time_integration_run

    subroutine tracer_transport_time_integration_hori(mesh,dtime)

        type(global_domain),  intent(inout) :: mesh
        real(r8)          ,   intent(in)    :: dtime

        ns_hori  = int(dtime/tracer_hori_timestep)
        DO ihori_step = 1, ns_hori
        DO irk_step = 1, nrk_hori

            rk_number  = nrk_hori+1-irk_step
            rk_substep = tracer_hori_timestep/real(rk_number,r4)
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
#ifndef SEQ_GRIST
            call clock_begin(clock_tracerh_flux)
#endif
            call tracer_transport_normal_flux_at_edge_r8(mesh, tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f_r4   ,& ! in
                                                               tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f  ,&
                                                               tracerVarCellFull%scalar_tracer_mass_n%f_r4               ,&
                                                               scalar_tracer_mxrt_at_pc_full_level_rk%f_r4              ,& ! in
                                                               scalar_normal_mxrt_at_edge_full_level_rk%f_r4            ,& ! out, still mixing ratio
                                                               tracer_hadv_flag             ,&
                                                               rk_substep,nlev,ntracer)
!
! Set to tracer mass
!
            do itracer = 1, ntracer
                scalar_normal_mxrt_at_edge_full_level_rk%f_r4(itracer,:,:) = scalar_normal_mxrt_at_edge_full_level_rk%f_r4(itracer,:,:)*&
                                                                            tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f
            end do
#ifndef SEQ_GRIST
            call clock_end(clock_tracerh_flux)
#endif
            if (do_limiter_hori(irk_step) .and. (tracer_hadv_flag .lt. 8 .or. tracer_hadv_flag .gt. 12) .and.&
                                                 tracer_hadv_flag.ne.58) mesh%ne = mesh%ne_compute
! note: halo2 will not use this tracer_adv_flag that need 3 halos
            if (do_limiter_hori(irk_step) .and. (tracer_hadv_flag .ge. 8 .and. tracer_hadv_flag .le. 12.or.&
                                                 tracer_hadv_flag.eq.58)) then
#ifndef SEQ_GRIST
                call exchange_data_3d_add(mesh,field_head_3d,scalar_normal_mxrt_at_edge_full_level_rk,"r4")
                call clock_begin(clock_tracerh_mpi)
                call exchange_data_3d_r4(mesh%local_block,field_head_3d)
                call clock_end(clock_tracerh_mpi)
#endif
            end if
#ifdef USE_HALO2
            call exchange_data_3d_add(mesh,field_head_3d,scalar_normal_mxrt_at_edge_full_level_rk,"r4")
            call clock_begin(clock_tracerh_mpi)
            call exchange_data_3d_r4(mesh%local_block,field_head_3d)
            call clock_end(clock_tracerh_mpi)
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
#ifndef SEQ_GRIST
                call clock_begin(clock_tracerh_flux)
#endif
                call tracer_transport_normal_flux_at_edge(mesh, tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f_r4  ,& ! in
                                                                tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f ,&
                                                                tracerVarCellFull%scalar_tracer_mass_n%f_r4              ,&
                                                                tracerVarCellFull%scalar_tracer_mxrt_n%f_r4              ,& ! in
                                                                scalar_normal_mxrt_at_edge_full_level_lo%f_r4           ,& ! out, still mixing ratio
                                                                1                         ,&
                                                                rk_substep,nlev,ntracer)
#ifndef SEQ_GRIST
                call clock_end(clock_tracerh_flux)
#endif

#ifdef USE_HALO2
                call exchange_data_3d_add(mesh,field_head_3d,scalar_normal_mxrt_at_edge_full_level_lo,"r4")
                call clock_begin(clock_tracerh_mpi)
                call exchange_data_3d_r4(mesh%local_block,field_head_3d)
                call clock_end(clock_tracerh_mpi)
#endif

! unchanged delhp*V*q(n) as low-order q-mass flux
#ifndef SEQ_GRIST
                call clock_begin(clock_tracerh_limiter)
#endif
                do itracer = 1, ntracer
                    scalar_normal_mxrt_at_edge_full_level_lo%f_r4(itracer,:,:) = scalar_normal_mxrt_at_edge_full_level_lo%f_r4(itracer,:,:)*&
                                                                              tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f
                end do
                mesh%ne = mesh%ne_compute

                call tracer_transport_hori_flux_limiter_fct(mesh, tracerVarCellFull%scalar_tracer_mxrt_n%f_r4   ,& !       q at time n
                                                                  tracerVarCellFull%scalar_tracer_mass_n%f_r4   ,& ! delhp*q at time n
                                                                  tracerVarCellFull%scalar_delhp_end_adv%f_r4   ,&
                                                                  scalar_normal_mxrt_at_edge_full_level_rk%f_r4,& ! high-order mass flux*q
                                                                  scalar_normal_mxrt_at_edge_full_level_lo%f_r4,& !  low-order mass flux*q
                                                                  scalar_normal_mxrt_at_edge_full_level_lm%f_r4,& !    limited mass flux*q
                                                                  rk_substep, ntracer, nlev)

                scalar_normal_mxrt_at_edge_full_level_rk%f_r4 = scalar_normal_mxrt_at_edge_full_level_lm%f_r4 ! reset rk flux by limited one

#ifndef SEQ_GRIST
                call clock_end(clock_tracerh_limiter)
#endif
            END IF
!
! DIV operator
!
#ifndef SEQ_GRIST
                call clock_begin(clock_tracerh_div)
#endif
!$omp parallel  private(iv,inb,edge_index,ilev,itracer) 
!$omp do schedule(dynamic,50) 
            do iv = 1, mesh%nv
                !div_sum(:,:) = 0._r8
                tracerVarCellFull%tend_tracer_mass_hori%f_r4(1:ntracer,1:nlev,iv) = zero
                do inb = 1, mesh%vtx_nnb(iv)
                    edge_index  = mesh%vtx_ed(inb,iv)
                    do ilev = 1, nlev
                        do itracer = 1, ntracer
                            !div_sum(itracer,ilev) = div_sum(itracer,ilev)+scalar_normal_mxrt_at_edge_full_level_rk%f(itracer,ilev,edge_index)*&
                            tracerVarCellFull%tend_tracer_mass_hori%f_r4(itracer,ilev,iv) = tracerVarCellFull%tend_tracer_mass_hori%f_r4(itracer,ilev,iv)+& 
                                                    scalar_normal_mxrt_at_edge_full_level_rk%f_r4(itracer,ilev,edge_index)*&
                                                    mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(edge_index)
                        end do
                    end do
                end do
                !tracerVarCellFull%tend_tracer_mass_hori%f(:,:,iv) = -div_sum(:,:)/((rearth**2)*mesh%plg_areag(iv))
                tracerVarCellFull%tend_tracer_mass_hori%f_r4(1:ntracer,1:nlev,iv) = -tracerVarCellFull%tend_tracer_mass_hori%f_r4(1:ntracer,1:nlev,iv)/((rearth**2)*mesh%plg_areag(iv))
            end do
  !$omp end do nowait
  !$omp end parallel 
#ifndef SEQ_GRIST
                call clock_end(clock_tracerh_div)
#endif
!
! tracer mass horizontal update
!
#ifndef SEQ_GRIST
                call clock_begin(clock_tracerh_update)
#endif
            tracerVarCellFull%scalar_tracer_mass_np1%f_r4 = tracerVarCellFull%scalar_tracer_mass_n%f_r4+rk_substep*tracerVarCellFull%tend_tracer_mass_hori%f_r4
#ifndef SEQ_GRIST
                call clock_end(clock_tracerh_update)
#endif
! 
! exchange data
!
#ifndef SEQ_GRIST
            call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_np1,"r4")
            !call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mxrt_np1,"r4")
            call clock_begin(clock_tracerh_mpi)
            call exchange_data_3d_r4(mesh%local_block,field_head_3d)
            call clock_end(clock_tracerh_mpi)
#endif
            do itracer = 1, ntracer
               tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4(itracer,:,1:mesh%nv_full) = tracerVarCellFull%scalar_tracer_mass_np1%f_r4(itracer,:,1:mesh%nv_full)/tracerVarCellFull%scalar_delhp_end_adv%f_r4(:,1:mesh%nv_full)
            end do
! prog
            ! scalar_tracer_mxrt_at_pc_full_level_rk    = tracerVarCellFull%scalar_tracer_mxrt_np1
            scalar_tracer_mxrt_at_pc_full_level_rk%f_r4    = tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4
            scalar_tracer_mxrt_at_pc_full_level_rk%pos     = tracerVarCellFull%scalar_tracer_mxrt_np1%pos
 
        END DO
        END DO

    end subroutine tracer_transport_time_integration_hori

    subroutine tracer_transport_time_integration_init(mesh)
        use grist_dycore_vars_module, only: dycoreVarCellFull
        use grist_tracer_transport_ffsl_module_mixed, only: tracer_transport_ffsl_flux_sm10_init
! io
        type(global_domain),  intent(inout) :: mesh

        !if(.not.allocated(div_sum))         allocate(div_sum(ntracer,nlev))
        if(.not.allocated(do_limiter_hori)) allocate(do_limiter_hori(nrk_hori))
        if(.not.allocated(do_limiter_vert)) allocate(do_limiter_vert(nrk_vert))
        if(.not.allocated(do_aimp_vert))    allocate(do_aimp_vert(nrk_vert))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_rk%f)) allocate(scalar_normal_mxrt_at_edge_full_level_rk%f(ntracer,nlev,mesh%ne_full))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_lo%f)) allocate(scalar_normal_mxrt_at_edge_full_level_lo%f(ntracer,nlev,mesh%ne_full))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_lm%f)) allocate(scalar_normal_mxrt_at_edge_full_level_lm%f(ntracer,nlev,mesh%ne_full))
! MIXCODE
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_rk%f_r4)) allocate(scalar_normal_mxrt_at_edge_full_level_rk%f_r4(ntracer,nlev,mesh%ne_full))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_lo%f_r4)) allocate(scalar_normal_mxrt_at_edge_full_level_lo%f_r4(ntracer,nlev,mesh%ne_full))
        if(.not.allocated(scalar_normal_mxrt_at_edge_full_level_lm%f_r4)) allocate(scalar_normal_mxrt_at_edge_full_level_lm%f_r4(ntracer,nlev,mesh%ne_full))
!
! when model starts, tracer mix ratio is given by initial state,
! scalar_delhp_at_pc_full_level_beg_adv is given by dyn
! adv_init only evaluates the initial state of tracer mass
! this tracer mass will be integrated since after,
! no need for init any more
!
        scalar_normal_mxrt_at_edge_full_level_rk%f    = 0._r4
        scalar_normal_mxrt_at_edge_full_level_rk%f_r4 = 0._r4
        scalar_normal_mxrt_at_edge_full_level_lo%f    = 0._r4
        scalar_normal_mxrt_at_edge_full_level_lo%f_r4 = 0._r4
        scalar_normal_mxrt_at_edge_full_level_lm%f    = 0._r4
        scalar_normal_mxrt_at_edge_full_level_lm%f_r4 = 0._r4
        tracerVarCellFull%scalar_delhp_end_adv%f_r4   = dycoreVarCellFull%scalar_delhp_n%f

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

#ifndef SEQ_GRIST
        clock_tracerh           = clock_id('tracerh')
        clock_tracerh_mpi       = clock_id('tracerh_mpi')
        clock_tracerh_flux      = clock_id('tracerh_flux')
        clock_tracerh_limiter   = clock_id('tracerh_limiter')
        clock_tracerh_div       = clock_id('tracerh_div')
        clock_tracerh_update    = clock_id('tracerh_update')
        clock_tracerv           = clock_id('tracerv')
#endif
     return
   end subroutine tracer_transport_time_integration_init

   subroutine tracer_transport_time_integration_final

        use grist_dycore_vars_module, only: dycoreVarCellFull
        use grist_tracer_transport_ffsl_module_mixed, only: tracer_transport_ffsl_flux_sm10_final

        !if(.not.allocated(div_sum))         allocate(div_sum(ntracer,nlev))
        if(allocated(do_limiter_hori)) deallocate(do_limiter_hori)
        if(allocated(do_limiter_vert)) deallocate(do_limiter_vert)
        if(allocated(do_aimp_vert))    deallocate(do_aimp_vert)
        if(allocated(scalar_normal_mxrt_at_edge_full_level_rk%f)) deallocate(scalar_normal_mxrt_at_edge_full_level_rk%f)
        if(allocated(scalar_normal_mxrt_at_edge_full_level_lo%f)) deallocate(scalar_normal_mxrt_at_edge_full_level_lo%f)
        if(allocated(scalar_normal_mxrt_at_edge_full_level_lm%f)) deallocate(scalar_normal_mxrt_at_edge_full_level_lm%f)
! MIXCODE
        if(allocated(scalar_normal_mxrt_at_edge_full_level_rk%f_r4)) deallocate(scalar_normal_mxrt_at_edge_full_level_rk%f_r4)
        if(allocated(scalar_normal_mxrt_at_edge_full_level_lo%f_r4)) deallocate(scalar_normal_mxrt_at_edge_full_level_lo%f_r4)
        if(allocated(scalar_normal_mxrt_at_edge_full_level_lm%f_r4)) deallocate(scalar_normal_mxrt_at_edge_full_level_lm%f_r4)

        if(trim(tracer_transport_hori_scheme).eq.'ffsl') call tracer_transport_ffsl_flux_sm10_final

       return
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
        real(r8), intent(in) :: dtime
        integer(i4), intent(in)  :: nlev, nlevp, ncell_all, ncell_cal, nEdge_all
#ifdef DBL_MASS_TRACER
        real(r8),  intent(in)    :: scalar_eta_mass_flux_at_pc_face(nlevp,nCell_all)
#else
        real(r4),  intent(in)    :: scalar_eta_mass_flux_at_pc_face(nlevp,nCell_all)
#endif
        real(r4),  intent(in)    :: scalar_delhp_at_pc_full_level(nlev,nCell_all)
        real(r4),  intent(in)    :: scalar_normal_velocity_at_pc_full_level(nlev,nEdge_all)
        real(r4),  intent(inout) :: scalar_eta_mass_flux_at_pc_face_exp(nlevp,nCell_all) ! we want
        real(r4),  intent(inout) :: scalar_eta_mass_flux_at_pc_face_imp(nlevp,nCell_all) ! we want
!local
        real(r4), parameter :: kesi=0.9_r4, alpha_max = 1.1_r4, alpha_min = 0.8_r4
        integer(i4) :: iv, ie, inb, ilev, mark
        real(r4)    :: outflow_div_at_pc_full_level(nlev,ncell_all)
        real(r4)    :: alpha_face, alpha_h, alpha_star_max, alpha_star_min, ggg

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
                mark       = int(0.5_r4*sign(one,scalar_eta_mass_flux_at_pc_face(ilev,iv))+0.5_r4)
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
        use grist_tridiagonal_solver, only: solve_tridiagonal1
! io
        type(global_domain),  intent(inout) :: mesh
        real(r8)          ,   intent(in)    :: dtime
! local
        real(r4) :: scalar_mxrt_flux_1d_nlevp_ho(nlevp)
        real(r4) :: scalar_mxrt_flux_1d_nlevp_lm(nlevp)
        real(r4) :: scalar_mxrt_flux_3d_nlevp_ho(ntracer,nlevp,mesh%nv_full)
        real(r4) :: scalar_mxrt_flux_3d_nlevp_lm(ntracer,nlevp,mesh%nv_full)
        real(r4) :: scalar_eta_mass_flux_at_pc_face_exp (nlevp,mesh%nv_full)
        real(r4) :: scalar_eta_mass_flux_at_pc_face_imp (nlevp,mesh%nv_full)
        real(r4) :: trid_aaa(nlev),trid_bbb(nlev),trid_ccc(nlev),trid_rrr(nlev),trid_out(nlev)

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
#ifdef DBL_MASS_TRACER
            scalar_eta_mass_flux_at_pc_face_exp = tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f
#else
            scalar_eta_mass_flux_at_pc_face_exp = tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f_r4
#endif

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
#ifdef DBL_MASS_TRACER
                                                    tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1:nlevp,1:mesh%nv_full)   , &
#else
                                                    tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f_r4(1:nlevp,1:mesh%nv_full)   , &
#endif
                                                    tracerVarCellFull%scalar_delhp_end_adv%f_r4(1:nlev,1:mesh%nv_full)            , &
                                                    tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f_r4(1:nlev,1:mesh%ne_full), &
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
                    call tracer_transport_vertic_flux_at_face(scalar_tracer_mxrt_at_pc_full_level_rk%f_r4(itracer,1:nlev,iv) ,&  ! full level q
                                                              scalar_eta_mass_flux_at_pc_face_exp(1:nlevp,iv)    ,&  ! m*etadot
                                                              tracerVarCellFull%scalar_delhp_end_adv%f_r4(1:nlev,iv) ,&  ! delhp full
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
                                                                    tracerVarCellFull%scalar_tracer_mass_n%f_r4(itracer,1:nlev ,iv) ,&  ! tracer mass n
                                                                    tracerVarCellFull%scalar_tracer_mxrt_n%f_r4(itracer,1:nlev ,iv) ,&  ! tracer mxrt n
                                                                    tracerVarCellFull%scalar_delhp_end_adv%f_r4(1:nlev,iv)          ,&  ! delhp np1
                                                                    scalar_mxrt_flux_3d_nlevp_ho(itracer,1:nlevp,iv),&  ! ho flux rk,  q
                                                                    scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlevp,iv),&  ! lm flux np1, q
                                                                    nlev, nlevp)
                    end do
                end do ! ncell
            END IF

            do itracer = 1, ntracer
                tracerVarCellFull%tend_tracer_mass_vert%f_r4(itracer,1:nlev,:) = &
                -(scalar_eta_mass_flux_at_pc_face_exp(2:nlevp,:)*scalar_mxrt_flux_3d_nlevp_lm(itracer,2:nlevp,:)-&
                  scalar_eta_mass_flux_at_pc_face_exp(1:nlev,:) *scalar_mxrt_flux_3d_nlevp_lm(itracer,1:nlev,:))
            end do
!
! tracer mass vertical update
!
            tracerVarCellFull%scalar_tracer_mass_np1%f_r4 = tracerVarCellFull%scalar_tracer_mass_n%f_r4+rk_substep*tracerVarCellFull%tend_tracer_mass_vert%f_r4
            do itracer = 1, ntracer
               tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4(itracer,:,:) = tracerVarCellFull%scalar_tracer_mass_np1%f_r4(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f_r4
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
                    two*tracerVarCellFull%scalar_delhp_end_adv%f_r4(ilev,iv)
                end do
                do itracer = 1, ntracer
                    trid_rrr(1:nlev) = two*tracerVarCellFull%scalar_tracer_mass_np1%f_r4(itracer,1:nlev,iv)
                    call solve_tridiagonal1(trid_aaa,trid_bbb,trid_ccc,trid_rrr,trid_out,nlev)
                    tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4(itracer,1:nlev,iv) = trid_out(1:nlev)
                    tracerVarCellFull%scalar_tracer_mass_np1%f_r4(itracer,1:nlev,iv) = tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4(itracer,1:nlev,iv)*&
                                                                                    tracerVarCellFull%scalar_delhp_end_adv%f_r4(1:nlev,iv)
                end do
            END DO
            END IF

!----------------------------------------------------
! END IEVA HERE; if IEVA not called, this routine
! must produce identical solutions as before
!----------------------------------------------------
#ifndef SEQ_GRIST
#ifndef TRACER_VERTNOEXC
        call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_np1,"r4")
        call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mxrt_np1,"r4")
        call exchange_data_3d_r4(mesh%local_block,field_head_3d)
#endif
#endif
! prog
           ! scalar_tracer_mxrt_at_pc_full_level_rk    = tracerVarCellFull%scalar_tracer_mxrt_np1
           scalar_tracer_mxrt_at_pc_full_level_rk%f_r4    = tracerVarCellFull%scalar_tracer_mxrt_np1%f_r4
           scalar_tracer_mxrt_at_pc_full_level_rk%pos     = tracerVarCellFull%scalar_tracer_mxrt_np1%pos

        END DO
        END DO

        case default
           if(iblock .eq. 0) print*,"tracer transport: tracer_vadv_flag configured incorrectly, please modify"
        end select

   end subroutine tracer_transport_time_integration_vert_aimp

end module grist_tracer_transport_time_integration_hvsplit_mixed
