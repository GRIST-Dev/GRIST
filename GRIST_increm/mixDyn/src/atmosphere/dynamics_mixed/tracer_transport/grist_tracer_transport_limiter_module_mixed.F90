module grist_tracer_transport_limiter_module_mixed

    use grist_constants,        only: i4, r8, rearth => rearth_ns, zero_r4 => zero_ns
    use grist_constants,        only: r4 => ns
    use grist_domain_types,     only: global_domain
! data
    use grist_tracer_module_vars_control_mixed, only: edp_leng, plg_nr, plg_areag

    implicit none

    real(r4), parameter :: one  = 1._r4
    real(r4), parameter :: zero = 0._r4
    real(r4), parameter :: half = 0.5_r4
#ifndef SPCODE
#ifdef REG_KESI
    real(r4), parameter :: kesi = 1.e-80_r8
#else
    real(r4), parameter :: kesi = TINY(1._r4)
#endif
#else
    ! real(r8), parameter :: kesi = 1.e-40_r8
    real(r8), parameter :: kesi = TINY(1._r8)
#endif
    real(r4), parameter :: bigy = 1.e16_r4

contains

!================================================
! horizontal flux limiter
!================================================

    subroutine tracer_transport_hori_flux_limiter_fct(mesh,&
                                                      scalar_tracer_mxrt_at_pc_full_level         ,&
                                                      scalar_tracer_mass_at_pc_full_level         ,&
                                                      scalar_delhp_at_pc_full_level_new           ,&
                                                      scalar_normal_tracer_hoflux_at_edge_full_level,&
                                                      scalar_normal_tracer_loflux_at_edge_full_level,&
                                                      scalar_normal_tracer_lmflux_at_edge_full_level,&
                                                      dtime, ntracer, nlev)
! io
        use omp_lib
        type(global_domain)    , intent(in)    :: mesh
        real(r4), allocatable  , intent(in)    :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)          !       q at time n
        real(r4), allocatable  , intent(in)    :: scalar_tracer_mass_at_pc_full_level(:,:,:)          ! delhp*q at time n
        real(r4), allocatable  , intent(in)    :: scalar_delhp_at_pc_full_level_new(:,:)                ! delhp at time np1
        real(r4), allocatable  , intent(in)    :: scalar_normal_tracer_hoflux_at_edge_full_level(:,:,:) ! ho flux, q*delhp*V
        real(r4), allocatable  , intent(in)    :: scalar_normal_tracer_loflux_at_edge_full_level(:,:,:) ! lo flux, q*delhp*V
        real(r4), allocatable  , intent(inout) :: scalar_normal_tracer_lmflux_at_edge_full_level(:,:,:) ! lm flux, q*delhp*V
        real(r8)               , intent(in)    :: dtime
        integer(i4)            , intent(in)    :: ntracer, nlev
! local
        real(r4), allocatable                  :: df_flux_at_edge(:,:,:)    ! ho-lo flux
        real(r4), allocatable                  :: scalar_td_at_cell(:,:,:)  ! TD solution
        real(r4), allocatable                  :: r_inflow_at_cell(:,:,:)   ! R
        real(r4), allocatable                  :: r_ouflow_at_cell(:,:,:)   ! R
        real(r4), allocatable                  :: qmax(:,:,:)
        real(r4), allocatable                  :: qmin(:,:,:)
        real(r4), allocatable                  :: p_inflow(:,:,:)
        real(r4), allocatable                  :: p_ouflow(:,:,:)
        real(r4)                               :: v0v1(3)
        real(r4)                               :: div_sum(ntracer,nlev)
        real(r4)                               :: maf_inflow, maf_ouflow
        real(r4)                               :: flag, coef
        integer(i4)                            :: cell_number, edge_number
        integer(i4)                            :: inb
        integer(i4)                            :: index_edge, index_cell, mark
        integer(i4)                            :: icell, ie, v0, v1, v_upwind,v_dowind,v_nb, ilev, itracer
        integer(i4)                            :: ii
! added for sat-limiter
        real(r4)  ::  sat_specific_humidity(nlev,mesh%nv_full), esvp(nlev,mesh%nv_full)

        cell_number = mesh%nv
        edge_number = mesh%ne

        allocate(scalar_td_at_cell(ntracer,nlev,1:mesh%nv_full))
        allocate(  df_flux_at_edge(ntracer,nlev,1:mesh%ne_full))
        allocate( r_inflow_at_cell(ntracer,nlev,1:mesh%nv_full))
        allocate( r_ouflow_at_cell(ntracer,nlev,1:mesh%nv_full))
        allocate(             qmax(ntracer,nlev,1:mesh%nv_full))
        allocate(             qmin(ntracer,nlev,1:mesh%nv_full))
        allocate(         p_inflow(ntracer,nlev,1:mesh%nv_full))
        allocate(         p_ouflow(ntracer,nlev,1:mesh%nv_full))

!$omp parallel  private(ii,ie,ilev,itracer) 
!$omp do schedule(static,100) 
        do ii = 1, mesh%ne_halo(1)*nlev*ntracer,1
            ie=ceiling(ii/real(nlev*ntracer,r4))
            ilev=ceiling((ii-(ie-1)*nlev*ntracer)/real(ntracer,r4))
            itracer=ii-(ie-1)*nlev*ntracer-(ilev-1)*ntracer
!    do ie = 1, mesh%ne_halo(1)!edge_number
!       do ilev = 1, nlev
!          do itracer = 1, ntracer 
            df_flux_at_edge(itracer,ilev,ie) = (scalar_normal_tracer_hoflux_at_edge_full_level(itracer,ilev,ie)-&
                                                scalar_normal_tracer_loflux_at_edge_full_level(itracer,ilev,ie))*rearth*edp_leng(ie)
!          end do
!       end do
!    end do
        end do
!$omp end do nowait
!$omp end parallel 
!
! [1] preupdate TD solution at cell
!
!$omp parallel  private(icell,div_sum,inb,index_edge,itracer,ilev) 
!$omp do schedule(dynamic,20) 
        do icell = 1, mesh%nv_halo(2)!cell_number     ! for each vertice in the mesh
            div_sum(:,:) = zero_r4
            do inb = 1, mesh%vtx_nnb(icell)
                index_edge = mesh%vtx_ed(inb,icell)
                do itracer = 1, ntracer
                    do ilev = 1, nlev
                        div_sum(itracer,ilev) = div_sum(itracer,ilev)+scalar_normal_tracer_loflux_at_edge_full_level(itracer,ilev,index_edge)*&
                                                                      plg_nr(inb,icell)*edp_leng(index_edge)*rearth
                    end do
                end do
            end do
            scalar_td_at_cell(:,:,icell) = scalar_tracer_mass_at_pc_full_level(:,:,icell)-dtime*div_sum(:,:)/((rearth**2)*plg_areag(icell))
        end do
!$omp end do nowait
!$omp end parallel
!
! mixing ratio
!
        do itracer = 1, ntracer
            scalar_td_at_cell(itracer,:,:) = scalar_td_at_cell(itracer,:,:)/scalar_delhp_at_pc_full_level_new(:,:)
        end do

! yizhang for sat-checker: [1.5]: calculate qsat

!     if(tracer_sat_limiter)then
!        call aqsat_grist(t=dycoreVarCellFull%scalar_temp_n%f, & ! nlev, nv_full
!                         p=scalar_pressure_at_pc_full_level_n%f, &
!                        es=esvp, &
!                        qs=sat_specific_humidity, &
!                        ii=mesh%nv_full, ilen=mesh%nv_halo(1), kk=nlev, kstart=1, kend=nlev)
       
!        ! convert to mxrt, but still use the same array; then do a limiter for the td solution
!        do icell = 1, mesh%nv_full
!           do ilev  = 1, nlev
!              sat_specific_humidity(ilev,icell) = sat_specific_humidity(ilev,icell)/tracerVarCellFull%scalar_mif_n%f(ilev,icell)
!              scalar_td_at_cell(1,ilev,icell)   = min(scalar_td_at_cell(1,ilev,icell),sat_specific_humidity(ilev,icell))
!           end do
!        end do
!     end if

!
! [2] compute R values at cell
!
        p_inflow  = zero
        p_ouflow  = zero
        qmax      =-bigy
        qmin      = bigy

!$omp parallel  private(icell,inb,ilev,itracer,index_edge,index_cell) 
!$omp do schedule(dynamic,20) 
        do icell   = 1, mesh%nv_halo(1) !cell_number
            do inb  = 1, mesh%vtx_nnb(icell)
                index_edge = mesh%vtx_ed(inb,icell)
                index_cell = mesh%vtx_nb(inb,icell)

                do ilev    = 1, nlev
                do itracer = 1, ntracer

                    qmax(itracer,ilev,icell)   = max(qmax(itracer,ilev,icell),scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,icell),scalar_td_at_cell(itracer,ilev,icell))
                    qmin(itracer,ilev,icell)   = min(qmin(itracer,ilev,icell),scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,icell),scalar_td_at_cell(itracer,ilev,icell))
                    
                    qmax(itracer,ilev,icell)   = max(qmax(itracer,ilev,icell), scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,index_cell), &
                                                                                                scalar_td_at_cell(itracer,ilev,index_cell))
                    qmin(itracer,ilev,icell)   = min(qmin(itracer,ilev,icell), scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,index_cell), &
                                                                                                    scalar_td_at_cell(itracer,ilev,index_cell))
                    p_inflow(itracer,ilev,icell)  = p_inflow(itracer,ilev,icell)+min(zero ,mesh%plg_nr(inb,icell)*df_flux_at_edge(itracer,ilev,index_edge))
                    p_ouflow(itracer,ilev,icell)  = p_ouflow(itracer,ilev,icell)+max(zero ,mesh%plg_nr(inb,icell)*df_flux_at_edge(itracer,ilev,index_edge))
                
                end do
                end do

            end do
        end do
!$omp end do nowait
!$omp end parallel 

!$omp parallel  private(ii,icell,ilev,itracer,maf_inflow,maf_ouflow) 
!$omp do schedule(dynamic,100) 
        do ii = 1, mesh%nv_halo(1)*nlev*ntracer,1
            icell=ceiling(ii/real(nlev*ntracer,r4))
            ilev=ceiling((ii-(icell-1)*nlev*ntracer)/real(ntracer,r4))
            itracer=ii-(icell-1)*nlev*ntracer-(ilev-1)*ntracer
!     do icell   = 1, mesh%nv_halo(1)!cell_number
!     do ilev    = 1, nlev
!     do itracer = 1, ntracer
            maf_inflow = scalar_delhp_at_pc_full_level_new(ilev,icell)*(qmax(itracer,ilev,icell)-scalar_td_at_cell(itracer,ilev,icell))*mesh%plg_areag(icell)*(rearth**2)/dtime
            maf_ouflow = scalar_delhp_at_pc_full_level_new(ilev,icell)*(scalar_td_at_cell(itracer,ilev,icell)-qmin(itracer,ilev,icell))*mesh%plg_areag(icell)*(rearth**2)/dtime
#ifdef REGRESSION_FCT
            if(p_inflow(itracer,ilev,icell).lt.zero)then
                r_inflow_at_cell(itracer,ilev,icell)  = min(one,maf_inflow/(abs(p_inflow(itracer,ilev,icell))))
            else
                r_inflow_at_cell(itracer,ilev,icell)  = zero
            end if

            if(p_ouflow(itracer,ilev,icell).gt.zero)then
                r_ouflow_at_cell(itracer,ilev,icell) = min(one,maf_ouflow/(abs(p_ouflow(itracer,ilev,icell))))
            else
                r_ouflow_at_cell(itracer,ilev,icell) = zero
            end if
#else
            r_inflow_at_cell(itracer,ilev,icell) = min(one,maf_inflow/(kesi+abs(p_inflow(itracer,ilev,icell))))
            r_ouflow_at_cell(itracer,ilev,icell) = min(one,maf_ouflow/(kesi+abs(p_ouflow(itracer,ilev,icell))))
#endif


!    end do
!    end do
!    end do
        end do
!$omp end do nowait
!$omp end parallel 
!
![3] compute limited flux at edge
!
!$omp parallel  private(ie,v0,v1,v0v1,flag,ilev,itracer,mark,v_upwind,v_dowind,coef) 
!$omp do schedule(dynamic,50) 
        do ie = 1, edge_number
            !v0                = mesh%edt(ie)%v(1)
            !v1                = mesh%edt(ie)%v(2)
            v0                = mesh%edt_v(1,ie)
            v1                = mesh%edt_v(2,ie)
            v0v1              = mesh%vtx_p(1:3,v1)-mesh%vtx_p(1:3,v0)
            flag              = sign(one,dot_product(v0v1,real(mesh%edp_nr(1:3,ie),r4)))
            do ilev = 1, nlev
            do itracer = 1, ntracer
                mark     = int(-0.5*sign(one,flag*df_flux_at_edge(itracer,ilev,ie))+1.5)
                !v_upwind = mesh%edt(ie)%v(mark)
                !v_dowind = mesh%edt(ie)%v(3-mark)
                v_upwind = mesh%edt_v(mark,ie)
                v_dowind = mesh%edt_v(3-mark,ie)
                coef     = min(r_inflow_at_cell(itracer,ilev,v_dowind),r_ouflow_at_cell(itracer,ilev,v_upwind))
                scalar_normal_tracer_lmflux_at_edge_full_level(itracer,ilev,ie) = scalar_normal_tracer_loflux_at_edge_full_level(itracer,ilev,ie)+&
                                                                            coef*(scalar_normal_tracer_hoflux_at_edge_full_level(itracer,ilev,ie)-&
                                                                                  scalar_normal_tracer_loflux_at_edge_full_level(itracer,ilev,ie))
            end do
            end do
        end do
!$omp end do nowait
!$omp end parallel 

        deallocate(scalar_td_at_cell)
        deallocate(df_flux_at_edge)
        deallocate(r_inflow_at_cell)
        deallocate(r_ouflow_at_cell)
        deallocate(p_inflow)
        deallocate(p_ouflow)
        deallocate(qmax)
        deallocate(qmin)

        return
    end subroutine tracer_transport_hori_flux_limiter_fct


!================================================
! vertical flux limiter
!================================================

    subroutine tracer_transport_vert_flux_limiter_fct(dtime, &
                                                      scalar_mass_eta_velocity_at_face_level, &
                                                      scalar_tracer_mass_at_full_level_n, &
                                                      scalar_tracer_mxrt_at_full_level_n,  &
                                                      scalar_delhp_at_full_level_np1,&
                                                      scalar_ho_flux_at_face_level,  &
                                                      scalar_lm_flux_at_face_level,  &
                                                      nlev, nlevp)
! io
        real(r8),              intent(in)    :: dtime
        !real(r8), allocatable, intent(in)    :: scalar_mass_eta_velocity_at_face_level(:) !downward positive
        !real(r8), allocatable, intent(in)    :: scalar_tracer_mass_at_full_level_n(:)
        !real(r8), allocatable, intent(in)    :: scalar_tracer_mxrt_at_full_level_n(:)
        !real(r8), allocatable, intent(in)    :: scalar_delhp_at_full_level_np1(:)
        !real(r8), allocatable, intent(in)    :: scalar_ho_flux_at_face_level(:)   ! q
        !real(r8), allocatable, intent(inout) :: scalar_lm_flux_at_face_level(:)   ! limited q
        real(r4), intent(in)    :: scalar_mass_eta_velocity_at_face_level(nlevp) !downward positive
        real(r4), intent(in)    :: scalar_tracer_mass_at_full_level_n(nlev)
        real(r4), intent(in)    :: scalar_tracer_mxrt_at_full_level_n(nlev)
        real(r4), intent(in)    :: scalar_delhp_at_full_level_np1(nlev)
        real(r4), intent(in)    :: scalar_ho_flux_at_face_level(nlevp)   ! q
        real(r4), intent(inout) :: scalar_lm_flux_at_face_level(nlevp)   ! limited q
        integer(i4),           intent(in)    :: nlev, nlevp
! local 
        real(r4)       :: scalar_lo_flux_at_face_level(nlevp) ! q
        real(r4)       :: scalar_ad_flux_at_face_level(nlevp) ! w.q
        real(r4)       :: scalar_tracer_mxrt_at_full_level_td(nlev)
        real(r4)       :: r_inflow(nlev), r_ouflow(nlev)
        real(r4)       :: mass_eta_velocity, coef
        real(r4)       :: qmax, qmin
        real(r4)       :: p_inflow, p_ouflow, m_inflow, m_ouflow
        integer(i4)    :: ilev, mark, v_dowind, v_upwind

!
!(1) compute 1st-order upwind flux and update scalar_td
! integration is for mass*q, limitation is for q
! all flux only contain q, not mass flux, we call it normalized flux
! mass flux remains the same as prescribed (consistency)
!
! face
        do ilev = 2, nlev
            mass_eta_velocity                  = scalar_mass_eta_velocity_at_face_level(ilev)
            scalar_lo_flux_at_face_level(ilev) = half*(                            scalar_tracer_mxrt_at_full_level_n(ilev)+scalar_tracer_mxrt_at_full_level_n(ilev-1))-&
                                                 half*sign(one,mass_eta_velocity)*(scalar_tracer_mxrt_at_full_level_n(ilev)-scalar_tracer_mxrt_at_full_level_n(ilev-1))
            scalar_ad_flux_at_face_level(ilev) = (scalar_ho_flux_at_face_level(ilev)-scalar_lo_flux_at_face_level(ilev))*mass_eta_velocity
        end do

        scalar_lo_flux_at_face_level(1)       = zero
        scalar_ad_flux_at_face_level(1)       = zero
        scalar_lo_flux_at_face_level(nlevp)   = zero
        scalar_ad_flux_at_face_level(nlevp)   = zero
! full
        do ilev = 1, nlev
            scalar_tracer_mxrt_at_full_level_td(ilev) = scalar_tracer_mass_at_full_level_n(ilev)-dtime*&
                                            (scalar_mass_eta_velocity_at_face_level(ilev+1)*scalar_lo_flux_at_face_level(ilev+1)-&
                                            scalar_mass_eta_velocity_at_face_level(ilev)  *scalar_lo_flux_at_face_level(ilev))
            scalar_tracer_mxrt_at_full_level_td(ilev) = scalar_tracer_mxrt_at_full_level_td(ilev)/scalar_delhp_at_full_level_np1(ilev)
        end do
!
! (2) compute R values at full level
!
        DO ilev = 1, nlev
            if(ilev.eq.1)then
                qmax      = max(scalar_tracer_mxrt_at_full_level_n(ilev), scalar_tracer_mxrt_at_full_level_n(ilev+1),&
                                scalar_tracer_mxrt_at_full_level_td(ilev),scalar_tracer_mxrt_at_full_level_td(ilev+1))
                qmin      = min(scalar_tracer_mxrt_at_full_level_n(ilev), scalar_tracer_mxrt_at_full_level_n(ilev+1),&
                                scalar_tracer_mxrt_at_full_level_td(ilev),scalar_tracer_mxrt_at_full_level_td(ilev+1))
            else if(ilev.eq.nlev)then
                qmax      = max(scalar_tracer_mxrt_at_full_level_n(ilev), scalar_tracer_mxrt_at_full_level_n(ilev-1),&
                                scalar_tracer_mxrt_at_full_level_td(ilev),scalar_tracer_mxrt_at_full_level_td(ilev-1))
                qmin      = min(scalar_tracer_mxrt_at_full_level_n(ilev), scalar_tracer_mxrt_at_full_level_n(ilev-1),&
                                scalar_tracer_mxrt_at_full_level_td(ilev),scalar_tracer_mxrt_at_full_level_td(ilev-1))
            else
                qmax      = max(scalar_tracer_mxrt_at_full_level_n(ilev), scalar_tracer_mxrt_at_full_level_n(ilev-1), scalar_tracer_mxrt_at_full_level_n(ilev+1),&
                                scalar_tracer_mxrt_at_full_level_td(ilev),scalar_tracer_mxrt_at_full_level_td(ilev-1),scalar_tracer_mxrt_at_full_level_td(ilev+1))
                qmin      = min(scalar_tracer_mxrt_at_full_level_n(ilev), scalar_tracer_mxrt_at_full_level_n(ilev-1), scalar_tracer_mxrt_at_full_level_n(ilev+1),&
                                scalar_tracer_mxrt_at_full_level_td(ilev),scalar_tracer_mxrt_at_full_level_td(ilev-1),scalar_tracer_mxrt_at_full_level_td(ilev+1))
            end if

            p_inflow  = max(zero,scalar_ad_flux_at_face_level(ilev))  -min(zero,scalar_ad_flux_at_face_level(ilev+1)) 
            p_ouflow  = max(zero,scalar_ad_flux_at_face_level(ilev+1))-min(zero,scalar_ad_flux_at_face_level(ilev))
            m_inflow  = (qmax-scalar_tracer_mxrt_at_full_level_td(ilev))*scalar_delhp_at_full_level_np1(ilev)/dtime
            m_ouflow  = (scalar_tracer_mxrt_at_full_level_td(ilev)-qmin)*scalar_delhp_at_full_level_np1(ilev)/dtime
#ifdef REGRESSION_FCT
            if(p_inflow.gt.zero)then
                r_inflow(ilev)  = min(one,m_inflow/(p_inflow))
            else
                r_inflow(ilev)  = zero
            end if
            if(p_ouflow.gt.zero)then
                r_ouflow(ilev)  = min(one,m_ouflow/(p_ouflow))
            else
                r_ouflow(ilev)  = zero
            end if
#else
            r_inflow(ilev)  = min(one,m_inflow/(kesi+p_inflow))
            r_ouflow(ilev)  = min(one,m_ouflow/(kesi+p_ouflow))
#endif
        END DO
!
! (3) flux limiting at face level, limiting q flux!
! 
        do ilev = 2, nlev
            mark     = int(-0.5_r4*sign(one,scalar_ad_flux_at_face_level(ilev))+0.5_r4) ! 0 or 1
            v_upwind = ilev-(1-mark)
            v_dowind = ilev-mark
            coef     = min(r_inflow(v_dowind),r_ouflow(v_upwind)) ! downwind associated with inflow, upwind associated with outflow
            scalar_lm_flux_at_face_level(ilev) = scalar_lo_flux_at_face_level(ilev)+coef*(scalar_ho_flux_at_face_level(ilev)-scalar_lo_flux_at_face_level(ilev))
        end do

        scalar_lm_flux_at_face_level(1)     = zero
        scalar_lm_flux_at_face_level(nlevp) = zero

        return
   end subroutine tracer_transport_vert_flux_limiter_fct

end module grist_tracer_transport_limiter_module_mixed
