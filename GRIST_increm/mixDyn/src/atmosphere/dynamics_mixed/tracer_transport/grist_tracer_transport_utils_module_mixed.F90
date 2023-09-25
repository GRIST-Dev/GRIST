module grist_tracer_transport_utils_module_mixed

    use grist_constants,        only: i4, r8, zero
    use grist_data_types,       only: scalar_2d_field, scalar_3d_field
    use grist_domain_types,     only: global_domain
    use grist_nml_module,       only: ntracer, nlev, write_verbose
    use grist_mpi,              only: mpi_rank, reduce

! data
    use grist_tracer_transport_vars_module, only: scalar_tracer_mxrt_min

    implicit none

contains

    subroutine tracer_transport_check_mxrt(mesh, scalar_tracer_mxrt_at_pc_full_level_n,called_env)
! io
        use omp_lib
        type(global_domain),   intent(in) :: mesh
        type(scalar_3d_field), intent(in) :: scalar_tracer_mxrt_at_pc_full_level_n
        character(len=*),      intent(in) :: called_env
! local
        integer(i4)     :: itracer, iv, ilev
        integer(i4)     :: nsum, nsum_global, ierr
        integer(i4)     :: ii

        nsum = 0
!$omp parallel  private(ii,iv,ilev,itracer) 
!$omp do schedule(static,100) reduction(+:nsum)
        do ii = 1, mesh%nv_compute*nlev*ntracer,1
            iv=ceiling(ii/real(nlev*ntracer,r8))
            ilev=ceiling((ii-(iv-1)*nlev*ntracer)/real(ntracer,r8))
            itracer=ii-(iv-1)*nlev*ntracer-(ilev-1)*ntracer
!        do iv = 1, mesh%nv_compute
!           do itracer = 1, ntracer
!              do ilev = 1, nlev
                if(scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev,iv).lt.scalar_tracer_mxrt_min(itracer))then
                    nsum = nsum + 1
                end if
!              end do
!           end do
!        end do
        end do
!$omp end do nowait
!$omp end parallel
#ifndef SEQ_GRIST
        call reduce(nsum, nsum_global, 'sum')
#else
        nsum_global = nsum
#endif
        if(mpi_rank().eq.0.and.write_verbose) print*,"nsum_global=",nsum_global,"in violation of min",trim(called_env)

        return
    end subroutine tracer_transport_check_mxrt

!----------------------------------------------------------------------------------
! This filling procedure is similar to that used by CAM-FV; Although horizontal
! transport is monotonic and refuses undershooting to machine round off, certain
! small negative values may still be generated in some ill conditions
! (diabatic/hvsplit/boundary)
!----------------------------------------------------------------------------------

    subroutine tracer_transport_fixer_mxrt(ncell,scalar_tracer_mxrt_at_pc_full_level_n,&
                                                 scalar_delhp_at_pc_full_level_n)
! io
        integer(i4),           intent(in)    :: ncell
        type(scalar_3d_field), intent(inout) :: scalar_tracer_mxrt_at_pc_full_level_n
        type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_full_level_n
! local
        integer(i4)     :: itracer, iv, ilev
        integer(i4)     :: nsum, nsum_global, ierr
        real(r8)        :: qkold, qfill_up, qfill_lo, qfill

        do iv = 1, ncell
            do itracer = 1, ntracer
                do ilev = 1, nlev
                    if(scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev,iv).lt.scalar_tracer_mxrt_min(itracer))then
             ! as long as tracer_mxrt violates minimum constraint, it will be reset to min value,
             ! regardless whether we borrow from upper or lower levels; when not
             ! borrowing, this will violate mass conservation, but such computational
             ! sources are quite small.
                        qkold = scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev,iv)
                        scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev,iv) = scalar_tracer_mxrt_min(itracer)
             ! borrow from upper or lower levels
                        qfill_up = zero
                        qfill_lo = zero
                        if(ilev.ne.1)    qfill_up = (scalar_tracer_mxrt_min(itracer)-qkold)*&
                                    scalar_delhp_at_pc_full_level_n%f_r4(ilev,iv)/scalar_delhp_at_pc_full_level_n%f_r4(ilev-1,iv)
                        if(ilev.ne.nlev) qfill_lo = (scalar_tracer_mxrt_min(itracer)-qkold)*&
                                    scalar_delhp_at_pc_full_level_n%f_r4(ilev,iv)/scalar_delhp_at_pc_full_level_n%f_r4(ilev+1,iv)
             ! check how to borrow:
             ! 1) from upper level:
#ifndef NEWFIXER
                        if(ilev.ne.1.and.(scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev-1,iv)-qfill_up).gt.scalar_tracer_mxrt_min(itracer))then
#else
                            if(ilev.ne.1)then
                                if((scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev-1,iv)-qfill_up).gt.scalar_tracer_mxrt_min(itracer))then

#endif
                                    scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev-1,iv) = scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev-1,iv)-qfill_up
                                    cycle
#ifdef NEWFIXER
                                end if
#endif
             ! 2) from lower level:
#ifndef NEWFIXER
                            else if(ilev.ne.nlev.and.(scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev+1,iv)-qfill_lo).gt.scalar_tracer_mxrt_min(itracer))then
#else
                            else if(ilev.ne.nlev)then
                                if((scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev+1,iv)-qfill_lo).gt.scalar_tracer_mxrt_min(itracer))then
#endif
                            scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev+1,iv) = scalar_tracer_mxrt_at_pc_full_level_n%f_r4(itracer,ilev+1,iv)-qfill_lo
                            cycle
#ifdef NEWFIXER
                            endif
#endif
             ! 3) not borrow:
                        else
                            cycle
                        end if
                    end if
                end do
            end do
        end do
       
        return
    end subroutine tracer_transport_fixer_mxrt

end module grist_tracer_transport_utils_module_mixed