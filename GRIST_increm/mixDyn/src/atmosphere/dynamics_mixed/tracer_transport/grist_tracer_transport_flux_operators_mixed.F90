module grist_tracer_transport_flux_operators_mixed

    use grist_constants,    only: i4, r8, ns, one => one_ns, half => half_ns, rearth => rearth_ns
    use grist_constants,    only: r4 => ns
    use grist_domain_types, only: global_domain
    use grist_tracer_transport_ffsl_module_mixed, only: tracer_transport_ffsl_flux_sm10
! data
    use grist_tracer_module_vars_control_mixed, only: edt_leng, edp_nr, vtx_p, edt_v_weight_prime_cell

    implicit none

#ifndef SPCODE
    real(r8), parameter  :: kesi = 1.e-80_r8
#else
    ! real(r8), parameter  :: kesi = 1.e-40_r8
    real(r8), parameter :: kesi = TINY(1._r8)
#endif

contains
!
! all normal flux are reconstructed as q, not multiplied by wind/mass flux
!

    subroutine tracer_transport_normal_flux_at_edge(mesh,scalar_normal_velocity_at_edge_full_level , &
                                                         scalar_normal_mass_flux_at_edge_full_level, &
                                                         scalar_tracer_mass_at_pc_full_level       , &
                                                         scalar_tracer_mxrt_at_pc_full_level       , &
                                                         scalar_normal_mxrt_at_edge_full_level     , &
                                                         adv_flag                                  , &
                                                         dtime, nlev, ntracer)
! io
        use omp_lib
        type(global_domain),   intent(in)    :: mesh
        real(r4), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:) ! nlev, ne
        real(r8), allocatable, intent(in)    :: scalar_normal_mass_flux_at_edge_full_level(:,:)! nlev, ne
        real(r4), allocatable, intent(in)    :: scalar_tracer_mass_at_pc_full_level(:,:,:)     ! ntracer, nlev, nv
        real(r4), allocatable, intent(in)    :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)     ! ntracer, nlev, nv
        real(r4), allocatable, intent(inout) :: scalar_normal_mxrt_at_edge_full_level(:,:,:)   ! ntracer, nlev, ne
        integer(i4),           intent(in)    :: adv_flag
        real(r8)   ,           intent(in)    :: dtime
        integer(i4),           intent(in)    :: nlev, ntracer
! local
        real(r4),dimension(1:3)              :: v0v1 !0->1
        real(r4)                             :: wind_edge
        real(r4)                             :: der0
        real(r4)                             :: part1, part2, part3
        integer(i4)                          :: v0, v1, v_upwind
        integer(i4)                          :: ie, itracer, ilev, inb, kk
        integer(i4)                          :: edge_index, cell_index, isum, flag, mark, kkk

        select case(adv_flag)
            
        case (1)  ! 1st-order upwind flux
!$omp parallel  private(ie,v0,v1,v0v1,flag,ilev,wind_edge,itracer)  
!$omp do schedule(dynamic,20) 
        do ie = 1, mesh%ne
            !v0        = mesh%edt(ie)%v(1)  ! 1st end of edge
            !v1        = mesh%edt(ie)%v(2)  ! 2nd end of edge
            v0        = mesh%edt_v(1,ie)   ! 1st end of edge
            v1        = mesh%edt_v(2,ie)   ! 2nd end of edge
            v0v1      = vtx_p(1:3,v1)-vtx_p(1:3,v0)
            flag      = sign(one,dot_product(v0v1,edp_nr(1:3,ie)))
            do ilev = 1, nlev
                wind_edge  = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
                do itracer = 1, ntracer
                    scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = &
                        half*(scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0)+scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1))-&
                        half*sign(one,wind_edge)*(scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1)-scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0))
                end do
            end do
        end do
!$omp end do nowait
!$omp end parallel 

        case(33)  ! purely upwind O3

!$omp parallel  private(ie,v0,v1,v0v1,flag,ilev,wind_edge,mark,v_upwind,itracer,der0,part1,part2)    
!$omp do schedule(dynamic,20) 
            do ie = 1, mesh%ne
                !v0        = mesh%edt(ie)%v(1)
                !v1        = mesh%edt(ie)%v(2)
                v0        = mesh%edt_v(1,ie)
                v1        = mesh%edt_v(2,ie)
                v0v1      = vtx_p(1:3,v1)-vtx_p(1:3,v0)
                flag      = sign(one,dot_product(v0v1,edp_nr(1:3,ie)))
                do ilev   = 1, nlev
                    wind_edge  = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
                    mark       = int(-0.5*sign(1._r4,wind_edge)+1.5)
                    !v_upwind   = mesh%edt(ie)%v(mark)    ! global index
                    v_upwind   = mesh%edt_v(mark,ie)    ! global index
                    do itracer = 1, ntracer
                        call calc_der_2nd_at_hx(mesh,scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind),v_upwind)),&
                                                    scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v_upwind),v_upwind,ie,mark-1,&
                                                    mesh%plg_stencil_number_2nd(v_upwind),der0)
                        part1 = (scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0)+&
                                scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1))*0.5_r4
                        part2 = -1._r4*((edt_leng(ie)**2)/6._r4)*der0
                        scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = part1+part2
                    end do
                end do
            end do
!$omp end do nowait
!$omp end parallel 

        case default
            print*,"TRACER TRANSPORT: you must select an advection scheme/flux divergence operator"

        end select

    end subroutine tracer_transport_normal_flux_at_edge

    subroutine tracer_transport_normal_flux_at_edge_r8(mesh,scalar_normal_velocity_at_edge_full_level , &
                                                            scalar_normal_mass_flux_at_edge_full_level, &
                                                            scalar_tracer_mass_at_pc_full_level       , &
                                                            scalar_tracer_mxrt_at_pc_full_level       , &
                                                            scalar_normal_mxrt_at_edge_full_level     , &
                                                            adv_flag                                  , &
                                                            dtime, nlev, ntracer)
    ! io
        use omp_lib
        type(global_domain),   intent(in)    :: mesh
        real(r4), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:) ! nlev, ne
        real(r8), allocatable, intent(in)    :: scalar_normal_mass_flux_at_edge_full_level(:,:)! nlev, ne
        real(r4), allocatable, intent(in)    :: scalar_tracer_mass_at_pc_full_level(:,:,:)     ! ntracer, nlev, nv
        real(r4), allocatable, intent(in)    :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)     ! ntracer, nlev, nv
        real(r4), allocatable, intent(inout) :: scalar_normal_mxrt_at_edge_full_level(:,:,:)   ! ntracer, nlev, ne
        integer(i4),           intent(in)    :: adv_flag
        real(r8)   ,           intent(in)    :: dtime
        integer(i4),           intent(in)    :: nlev, ntracer
! local
        real(r4), allocatable                :: scalar_tangen_velocity_at_edge_full_level(:,:) ! nlev, ne
        real(r4),dimension(1:3)              :: v0v1 !0->1
        real(r4)                             :: beta
        real(r4)                             :: wind_edge
        real(r4)                             :: der0, der1
        real(r4)                             :: der0_4th, der1_4th
        real(r4)                             :: part1, part2, part3
        integer(i4)                          :: v0, v1, v_upwind
        integer(i4)                          :: ie, itracer, ilev, inb, kk
        integer(i4)                          :: edge_index, cell_index, isum, flag, mark, kkk

        select case(adv_flag)
            
        case (1)  ! 1st-order upwind flux
!$omp parallel  private(ie,v0,v1,v0v1,flag,ilev,wind_edge,itracer)  
!$omp do schedule(dynamic,20) 
        do ie = 1, mesh%ne
            !v0        = mesh%edt(ie)%v(1)  ! 1st end of edge
            !v1        = mesh%edt(ie)%v(2)  ! 2nd end of edge
            v0        = mesh%edt_v(1,ie)   ! 1st end of edge
            v1        = mesh%edt_v(2,ie)   ! 2nd end of edge
            v0v1      = vtx_p(1:3,v1)-vtx_p(1:3,v0)
            flag      = sign(one,dot_product(v0v1,edp_nr(1:3,ie)))
            do ilev = 1, nlev
                wind_edge  = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
                do itracer = 1, ntracer
                    scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = &
                        half*(scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0)+scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1))-&
                        half*sign(one,wind_edge)*(scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1)-scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0))
                end do
            end do
        end do
!$omp end do nowait
!$omp end parallel 

        case(2)  ! 2nd-order center difference

            do ie = 1, mesh%ne
            !v0                               = mesh%edt(ie)%v(1)
            !v1                               = mesh%edt(ie)%v(2)
            v0                               = mesh%edt_v(1,ie)
            v1                               = mesh%edt_v(2,ie)
            do ilev = 1, nlev
                wind_edge                     = scalar_normal_velocity_at_edge_full_level(ilev,ie)
                do itracer = 1, ntracer
                    scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = &
                    half*(scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0)+&
                        scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1))
                end do
            end do
            end do

        case(3,4,5,6,7)  ! O3/4, damping depends on beta value

            if(adv_flag.eq.3) beta=1        ! 3rd order
            if(adv_flag.eq.4) beta=0        ! 4th order
            if(adv_flag.eq.5) beta=0.25     ! hybrid 1
            if(adv_flag.eq.6) beta=0.5      ! hybrid 2
            if(adv_flag.eq.7) beta=0.75     ! hybrid 3
        
!$omp parallel  private(ie,v0,v1,v0v1,flag,ilev,wind_edge,itracer,der0,der1,part1,part2,part3)    
!$omp do schedule(dynamic,50) 
            do ie = 1, mesh%ne
                !v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
                !v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
                v0        = mesh%edt_v(1,ie)       ! 1st end of edge
                v1        = mesh%edt_v(2,ie)       ! 2nd end of edge
                v0v1      = mesh%vtx_p(1:3,v1)-mesh%vtx_p(1:3,v0)
                flag      = sign(one,dot_product(v0v1,real(mesh%edp_nr(1:3,ie),r8)))
                do ilev   = 1, nlev
                    wind_edge = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
                    do itracer = 1, ntracer
                        call calc_der_2nd_at_hx(mesh,scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v0),v0)),&
                                                     scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0),v0,ie,0,mesh%plg_stencil_number_2nd(v0),der0)
                        call calc_der_2nd_at_hx(mesh,scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v1),v1)),&
                                                     scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1),v1,ie,1,mesh%plg_stencil_number_2nd(v1),der1)
                        part1 = (scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0)+&
                                 scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1))*0.5_r4
                        part2 = -1._r4*((mesh%edt_leng(ie)**2)/12._r4)*(der0+der1)
                        part3 = (sign(1._r4,wind_edge)*(mesh%edt_leng(ie)**2)*beta/12._r4)*(der1-der0)
                        scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = (part1+part2+part3)
                    end do
                end do
            end do
!$omp end do nowait
!$omp end parallel 

        ! case(8,9,10,11,12)  ! O5/6

        !     if(adv_flag.eq.8)  beta=1        ! 5rd order
        !     if(adv_flag.eq.9)  beta=0        ! 6th order
        !     if(adv_flag.eq.10) beta=0.25     ! hybrid
        !     if(adv_flag.eq.11) beta=0.5      ! hybrid
        !     if(adv_flag.eq.12) beta=0.75     ! hybrid
       
        !     do ie = 1, mesh%ne
        !           !v0        = mesh%edt(ie)%v(1)       ! 1st end of edge
        !           !v1        = mesh%edt(ie)%v(2)       ! 2nd end of edge
        !           v0        = mesh%edt_v(1,ie)       ! 1st end of edge
        !           v1        = mesh%edt_v(2,ie)       ! 2nd end of edge
        !           v0v1      = mesh%vtx_p(1:3,v1)-mesh%vtx_p(1:3,v0)
        !           flag      = sign(one,dot_product(v0v1,real(mesh%edp_nr(1:3,ie),r8)))
        !           do ilev = 1, nlev
        !              wind_edge = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
        !              do itracer = 1, ntracer
        !                 call calc_der_4th_at_hx(mesh,scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,mesh%plg_stencil_index_4th(1:mesh%plg_stencil_number_4th(v0),v0)),&
        !                                              scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0),v0,ie,0,mesh%plg_stencil_number_4th(v0),der0_4th,der0)
        !                 call calc_der_4th_at_hx(mesh,scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,mesh%plg_stencil_index_4th(1:mesh%plg_stencil_number_4th(v1),v1)),&
        !                                              scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1),v1,ie,1,mesh%plg_stencil_number_4th(v1),der1_4th,der1)
        !                 call flux_wrf56(wind_edge,scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0),&
        !                                           scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1),&
        !                                           der0,der1,der0_4th,der1_4th,scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie),&
        !                                           real(mesh%edt_leng(ie),r8),beta)
        !              end do
        !           end do
        !     end do

        case(33)  ! purely upwind O3

!$omp parallel  private(ie,v0,v1,v0v1,flag,ilev,wind_edge,mark,v_upwind,itracer,der0,part1,part2)    
!$omp do schedule(dynamic,20) 
            do ie = 1, mesh%ne
                !v0        = mesh%edt(ie)%v(1)
                !v1        = mesh%edt(ie)%v(2)
                v0        = mesh%edt_v(1,ie)
                v1        = mesh%edt_v(2,ie)
                v0v1      = vtx_p(1:3,v1)-vtx_p(1:3,v0)
                flag      = sign(one,dot_product(v0v1,edp_nr(1:3,ie)))
                do ilev   = 1, nlev
                    wind_edge  = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
                    mark       = int(-0.5*sign(1._r4,wind_edge)+1.5)
                    !v_upwind   = mesh%edt(ie)%v(mark)    ! global index
                    v_upwind   = mesh%edt_v(mark,ie)    ! global index
                    do itracer = 1, ntracer
                        call calc_der_2nd_at_hx_r8(mesh,scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind),v_upwind)),&
                                                    scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v_upwind),v_upwind,ie,mark-1,&
                                                    mesh%plg_stencil_number_2nd(v_upwind),der0)
                        part1 = (scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v0)+&
                                scalar_tracer_mxrt_at_pc_full_level(itracer,ilev,v1))*0.5_r4
                        part2 = -1._r4*((edt_leng(ie)**2)/6._r4)*der0
                        scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = part1+part2
                    end do
                end do
            end do
!$omp end do nowait
!$omp end parallel 

        case(28)  ! MS13, UQA-1
            !
            ! for FFSL flux, first reconstruct tangent wind based on normal velocity,
            ! then obtain the flux
            if(.not.allocated(scalar_tangen_velocity_at_edge_full_level)) allocate(scalar_tangen_velocity_at_edge_full_level(nlev,mesh%ne))
    
            call reconstruct_tangent_wind(mesh,scalar_normal_velocity_at_edge_full_level,&
                                               scalar_tangen_velocity_at_edge_full_level,&
                                               nlev)

            call tracer_transport_ffsl_flux_sm10(mesh,scalar_normal_velocity_at_edge_full_level , &
                                                      scalar_tangen_velocity_at_edge_full_level , &
                                                      scalar_tracer_mxrt_at_pc_full_level       , &
                                                      scalar_normal_mxrt_at_edge_full_level     , &
                                                      real(dtime,ns), nlev, ntracer)

       if(allocated(scalar_tangen_velocity_at_edge_full_level)) deallocate(scalar_tangen_velocity_at_edge_full_level)

        case default
            print*,"TRACER TRANSPORT: you must select an advection scheme/flux divergence operator"

        end select

    end subroutine tracer_transport_normal_flux_at_edge_r8

    subroutine tracer_transport_vertic_flux_at_face(scalar_at_full_level      ,&
                                                    scalar_mass_eta_velocity  ,&
                                                    scalar_delhp_at_full_level,&
                                                    scalar_at_face_level      ,&
                                                    dtime, order, nlev, nlevp )
! io
        real(r4), intent(in)    :: scalar_at_full_level(nlev)
        real(r4), intent(in)    :: scalar_mass_eta_velocity(nlevp)
        real(r4), intent(in)    :: scalar_delhp_at_full_level(nlev)
        real(r4), intent(inout) :: scalar_at_face_level(nlevp)
        real(r8), intent(in)    :: dtime
        integer(i4) ,  intent(in)            :: order
        integer(i4) ,  intent(in)            :: nlev
        integer(i4) ,  intent(in)            :: nlevp
! local 
        real(r4)                 :: scalar_delhp_at_face_level(nlev+1)
        real(r4)                 :: part1, der_k, der_kp1
        real(r4)                 :: qhat(nlev+1), cr_num
        integer(i4)              :: ilev,flag
!
! extrapolate, actually not used because etadot at boundaries are zero
! 
        scalar_at_face_level(1)        = scalar_at_full_level(1)
        scalar_at_face_level(nlev+1)   = scalar_at_full_level(nlev)
!
! compute delhp at face
!
        do ilev = 2, nlev
            scalar_delhp_at_face_level(ilev) = 0.5_r4*(scalar_delhp_at_full_level(ilev-1)+scalar_delhp_at_full_level(ilev))
        end do
        scalar_delhp_at_face_level(1)        = 0.5_r4*scalar_delhp_at_full_level(1)
        scalar_delhp_at_face_level(nlev+1)   = 0.5_r4*scalar_delhp_at_full_level(nlev)
!
! This is equivalent distance version, old version used for simplicity
! only used for bit reproducity
!
        select case(order)

        case(1)
            scalar_at_face_level(1)        = 0._r4
            scalar_at_face_level(nlev+1)   = 0._r4
            do ilev = 2, nlev
                scalar_at_face_level(ilev) = half*(scalar_at_full_level(ilev)+scalar_at_full_level(ilev-1))-&
                                             half*sign(one,scalar_mass_eta_velocity(ilev))*&
                                            (scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))
            end do
        case(2)
            scalar_at_face_level(1)       = 0._r4
            scalar_at_face_level(nlev+1)  = 0._r4
            do ilev = 2, nlev
                scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+&
                                              scalar_at_full_level(ilev-1))*half
            end do

        case(3) 
            scalar_at_face_level(1)       = 0._r4
            scalar_at_face_level(2)       = (scalar_at_full_level(2)   +scalar_at_full_level(1))*half
            scalar_at_face_level(nlev)    = (scalar_at_full_level(nlev)+scalar_at_full_level(nlev-1))*half
            scalar_at_face_level(nlev+1)  = 0._r4
            do ilev = 3, nlev-1
                scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+scalar_at_full_level(ilev-1))*(7._r4/12._r4)-&
                                             (scalar_at_full_level(ilev+1)+scalar_at_full_level(ilev-2))*(1._r4/12._r4)+&
                                   sign(1._r4,scalar_mass_eta_velocity(ilev))*&
                                            ((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev-2))-&
                                       3._r4*(scalar_at_full_level(ilev)  -scalar_at_full_level(ilev-1)))/12._r4
            end do

        case(9) !ppm
! set bdy
            qhat(1)       = 0._r4
            qhat(2)       = (scalar_at_full_level(2)+scalar_at_full_level(1))*half
            qhat(nlev)    = (scalar_at_full_level(nlev)+scalar_at_full_level(nlev-1))*half
            qhat(nlev+1)  = 0._r4
! compute qhat
            do ilev = 3, nlev-1
                qhat(ilev) = (7._r4*(scalar_at_full_level(ilev-1)+scalar_at_full_level(ilev))-&
                                    (scalar_at_full_level(ilev-2)+scalar_at_full_level(ilev+1)))/12._r4
            end do
! compute qface
            scalar_at_face_level(1)       = 0._r4
            scalar_at_face_level(nlev+1)  = 0._r4
            do ilev = 2, nlev
                cr_num = scalar_mass_eta_velocity(ilev)*dtime/scalar_delhp_at_face_level(ilev)
                flag   = int((cr_num-abs(cr_num))/(2*cr_num+kesi))
                if(flag.ne.0.and.flag.ne.1)then
                    print*,"flag in vadv 9 is wrong, stop"
                end if
                scalar_at_face_level(ilev) = qhat(ilev)-cr_num*(qhat(ilev)-scalar_at_full_level(ilev-1+flag))-&
                                                        abs(cr_num)*(1._r4-abs(cr_num))*&
                                                        (qhat(ilev-1+flag)-2._r4*scalar_at_full_level(ilev-1+flag)+qhat(ilev+flag))
            end do
        case default
            print*,"you must select a vertical order, stop"
            stop
        end select

        return
    end subroutine tracer_transport_vertic_flux_at_face

!----------------------------------------------------------
! private routines for reconstruction, adapted from
! core dynamics
!----------------------------------------------------------

    subroutine calc_der_2nd_at_hx(mesh,s_array,s0,iv,ie,flag,nlength,der)
! io
        type(global_domain) , intent(in) :: mesh
        real(r4)            , intent(in) :: s_array(:) 
        real(r4)            , intent(in) :: s0
        integer(i4)         , intent(in) :: iv
        integer(i4)         , intent(in) :: ie
        integer(i4)         , intent(in) :: flag  ! 0 or 1
        integer(i4)         , intent(in) :: nlength
        real(r4)            , intent(out):: der
! local
        real(r4)                         :: m_array(1:nlength)
        real(r4)                         :: f_array

        m_array(:) = s_array(:)-s0
        f_array = dot_product(edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(iv),flag+1,ie), m_array(:))
        der     = 2._r4*f_array
        return
    end subroutine calc_der_2nd_at_hx

    subroutine calc_der_2nd_at_hx_r8(mesh,s_array,s0,iv,ie,flag,nlength,der)
! io
        type(global_domain) , intent(in) :: mesh
        real(r4)            , intent(in) :: s_array(:) 
        real(r4)            , intent(in) :: s0
        integer(i4)         , intent(in) :: iv
        integer(i4)         , intent(in) :: ie
        integer(i4)         , intent(in) :: flag  ! 0 or 1
        integer(i4)         , intent(in) :: nlength
        real(r4)            , intent(out):: der
! local
        real(r4)                         :: m_array(1:nlength)
        real(r4)                         :: f_array

        m_array(:) = s_array(:)-s0
        f_array = dot_product(edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(iv),flag+1,ie), m_array(:))
        der     = 2._r4*f_array
        return
    end subroutine calc_der_2nd_at_hx_r8

   subroutine reconstruct_tangent_wind(mesh, scalar_normal_velocity_at_edge_full_level,&
                                             scalar_tangen_velocity_at_edge_full_level,&
                                             nlev)
! io
   type(global_domain),   intent(in)    :: mesh
   real(r4), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r4), allocatable, intent(inout) :: scalar_tangen_velocity_at_edge_full_level(:,:)
   integer(i4),           intent(in)    :: nlev
! local
   integer(i4)                          :: ie, ilev
   real(r4)                             :: length_of_triangle
   real(r4)                             :: cell_sum(16)

     do ie = 1, mesh%ne  ! for each edge e'
        length_of_triangle = rearth*mesh%edt_leng(ie)
        do ilev = 1, nlev
           cell_sum(1:mesh%edp_nedge(ie)) = mesh%edp_trsk_edpl_on_edge(1:mesh%edp_nedge(ie),ie)*&
                                            scalar_normal_velocity_at_edge_full_level(ilev,mesh%edp_edge_on_edge(1:mesh%edp_nedge(ie),ie))
           scalar_tangen_velocity_at_edge_full_level(ilev,ie) = sum(cell_sum(1:mesh%edp_nedge(ie)))/length_of_triangle
        end do
     end do

   return
  end subroutine reconstruct_tangent_wind

end module grist_tracer_transport_flux_operators_mixed
