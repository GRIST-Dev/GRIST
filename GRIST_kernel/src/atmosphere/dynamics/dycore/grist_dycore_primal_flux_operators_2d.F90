
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Version 1.0
! Description: High-level flux operators used by solvers. This is 2D version 
!              of flux operators specially developed for core dynamics
! Revision history:
!----------------------------------------------------------------------------

 module grist_dycore_primal_flux_operators_2d

   use grist_constants       , only: gravity, i4, r8, zero
   use grist_domain_types    , only: global_domain
   use grist_math_module     , only: norm
   use grist_base_flux_module, only: flux_ooo34, flux_ooo33

   implicit none
   private

   public :: calc_primal_normal_flux_at_edge, &
             calc_primal_normal_flux_at_edge_var2, &
             calc_primal_normal_flux_at_edge_1d

   real(r8), parameter  :: one = 1._r8

  contains

!================================================
! evaluation of two tendencies;
! they must have the same format
! note: 33 is our default and tuned
! bespoke for nhd www&phi
!================================================

 subroutine calc_primal_normal_flux_at_edge_var2(mesh,scalar_normal_velocity_at_edge_full_level, &
                                                      scalar_tangen_velocity_at_edge_full_level, & 
                                                      scalar1_height_at_pc_full_level          , &
                                                      scalar2_height_at_pc_full_level          , &
                                                      scalar1_normal_flux_at_edge_full_level   , &
                                                      scalar2_normal_flux_at_edge_full_level   , &
                                                      adv_flag                                 , &
                                                      dtime                                    , &
                                                      nlev, called_by_ndc)
! io
   use omp_lib
   type(global_domain),   intent(in)    :: mesh
   real(r8),              intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8),              intent(in)    :: scalar_tangen_velocity_at_edge_full_level(:,:)
   real(r8),              intent(in)    :: scalar1_height_at_pc_full_level(:,:)
   real(r8),              intent(in)    :: scalar2_height_at_pc_full_level(:,:)
   real(r8),              intent(inout) :: scalar1_normal_flux_at_edge_full_level(:,:) ! also q
   real(r8),              intent(inout) :: scalar2_normal_flux_at_edge_full_level(:,:) ! also q
   integer(i4),           intent(in)    :: adv_flag
   real(r8)   ,           intent(in)    :: dtime
   integer(i4),           intent(in)    :: nlev
   logical    ,           intent(in)    :: called_by_ndc
! local
   real(r8)                             :: beta
   real(r8)                             :: wind_edge
   real(r8)                             :: der1_0, der1_1, der2_0, der2_1
   integer(i4)                          :: v1, v2, v_upwind
   integer(i4)                          :: ilev, ie
   integer(i4)                          :: flag, mark
   integer(i4)                          :: ii

   select case(adv_flag)

   case(2)  ! 2nd-order center difference 

       do ie = 1, mesh%ne
          v1      = mesh%edt_v(1, ie)
          v2      = mesh%edt_v(2, ie)
          do ilev = 1, nlev
             scalar1_normal_flux_at_edge_full_level(ilev, ie) = 0.5*(scalar1_height_at_pc_full_level(ilev,v1)+&
                                                                     scalar1_height_at_pc_full_level(ilev,v2))
             scalar2_normal_flux_at_edge_full_level(ilev, ie) = 0.5*(scalar2_height_at_pc_full_level(ilev,v1)+&
                                                                     scalar2_height_at_pc_full_level(ilev,v2))
          end do
       end do

   case(3,4,5,6,7)  ! WRF3/4, damping depends on beta value

    if(adv_flag.eq.3) beta=1._r8    ! 3rd order
    if(adv_flag.eq.4) beta=0        ! 4th order
    if(adv_flag.eq.5) beta=0.25     ! hybrid 1
    if(adv_flag.eq.6) beta=0.5      ! hybrid 2
    if(adv_flag.eq.7) beta=0.75     ! hybrid 3

!$omp parallel  private(ii,ie,ilev, wind_edge,v1,v2,der1_0, der1_1, der2_0, der2_1) 
!$omp do schedule(dynamic,100) 
    do ii = 1, mesh%ne*nlev,1
        ie=ceiling(ii/real(nlev,r8))
        ilev=ii-(ie-1)*nlev

          v1        = mesh%edt_v(1, ie)
          v2        = mesh%edt_v(2, ie)

             wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge_full_level(ilev,ie)

! 1st cell
             call calc_der_2nd_at_hx_var2(mesh,scalar1_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v1), v1)),&
                                               scalar2_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v1), v1)),&
                                               scalar1_height_at_pc_full_level(ilev,v1), &
                                               scalar2_height_at_pc_full_level(ilev,v1), &
                                               v1,ie,1,mesh%plg_stencil_number_2nd(v1),der1_0,der2_0)
! 2nd cell
             call calc_der_2nd_at_hx_var2(mesh,scalar1_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v2), v2)),&
                                               scalar2_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v2), v2)),&
                                               scalar1_height_at_pc_full_level(ilev,v2), &
                                               scalar2_height_at_pc_full_level(ilev,v2), &
                                               v2,ie,2,mesh%plg_stencil_number_2nd(v2),der1_1,der2_1)
             call flux_ooo34(wind_edge, scalar1_height_at_pc_full_level(ilev,v1),&
                                        scalar1_height_at_pc_full_level(ilev,v2),&
                             der1_0,der1_1, scalar1_normal_flux_at_edge_full_level(ilev,ie),&
                             real(mesh%edt_leng(ie),r8),beta)

             call flux_ooo34(wind_edge, scalar2_height_at_pc_full_level(ilev,v1),&
                                        scalar2_height_at_pc_full_level(ilev,v2),&
                             der2_0,der2_1, scalar2_normal_flux_at_edge_full_level(ilev,ie),&
                             real(mesh%edt_leng(ie),r8),beta)
    end do
!$omp end do nowait
!$omp end parallel 

   case(33)  ! pure 3ord

!$omp parallel  private(ii,ie,ilev, wind_edge,v1,v2,mark,v_upwind,der1_0,der2_0) 
!$omp do schedule(dynamic,100) 
    do ii = 1, mesh%ne*nlev,1
        ie=ceiling(ii/real(nlev,r8))
        ilev=ii-(ie-1)*nlev

!    do ie = 1, mesh%ne
          v1        = mesh%edt_v(1, ie)       ! 1st end of edge
          v2        = mesh%edt_v(2, ie)       ! 2nd end of edge
!          do ilev   = 1, nlev
             wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge_full_level(ilev,ie)
             mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)
             v_upwind  = mesh%edt_v(mark, ie)
! only upwind cell 

! inline_opt_here
            ! call calc_der_2nd_at_hx_var2(mesh,scalar1_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind), v_upwind)),&
            !                                   scalar2_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind), v_upwind)),&
            !                                   scalar1_height_at_pc_full_level(ilev,v_upwind),&
            !                                   scalar2_height_at_pc_full_level(ilev,v_upwind),&
            !                                   v_upwind,ie,mark,mesh%plg_stencil_number_2nd(v_upwind),der1_0,der2_0)

            ! call flux_ooo33(scalar1_height_at_pc_full_level(ilev,v1),&
            !                scalar1_height_at_pc_full_level(ilev,v2),&
            !                der1_0, scalar1_normal_flux_at_edge_full_level(ilev,ie) ,&
            !                real(mesh%edt_leng(ie),r8))

            ! call flux_ooo33(scalar2_height_at_pc_full_level(ilev,v1),&
            !                scalar2_height_at_pc_full_level(ilev,v2),&
            !                der2_0, scalar2_normal_flux_at_edge_full_level(ilev,ie) ,&
            !                real(mesh%edt_leng(ie),r8))

               der1_0      = 2._r8*dot_product(mesh%edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(v_upwind), mark, ie), &
               scalar1_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind),v_upwind))-scalar1_height_at_pc_full_level(ilev,v_upwind))
               der2_0      = 2._r8*dot_product(mesh%edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(v_upwind), mark, ie), &
               scalar2_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind),v_upwind))-scalar2_height_at_pc_full_level(ilev,v_upwind))

               scalar1_normal_flux_at_edge_full_level(ilev,ie) = (scalar1_height_at_pc_full_level(ilev,v1)+&
                                                                  scalar1_height_at_pc_full_level(ilev,v2))*0.5_r8-&
                                                                 ((real(mesh%edt_leng(ie),r8)**2)/6._r8)*der1_0

               scalar2_normal_flux_at_edge_full_level(ilev,ie) = (scalar2_height_at_pc_full_level(ilev,v1)+&
                                                                  scalar2_height_at_pc_full_level(ilev,v2))*0.5_r8-&
                                                                 ((real(mesh%edt_leng(ie),r8)**2)/6._r8)*der2_0
! inline_opt_here

!          end do
!    end do
    end do
!$omp end do nowait
!$omp end parallel 

   case default
        print*,"you must select an advection scheme/flux divergence operator in gcm_para"
   end select

#ifndef REGRESSION
   do ie = 1, mesh%ne
      v1      = mesh%edt_v(1, ie)
      v2      = mesh%edt_v(2, ie)
      scalar2_normal_flux_at_edge_full_level(nlev,ie) = 0.5_r8*(scalar2_height_at_pc_full_level(nlev,v1)+scalar2_height_at_pc_full_level(nlev,v2))
   end do
#endif

#ifdef NHD_NOOPT 
! compute value at halo(1) for bdy-domain only
#ifdef LAM_DOMAIN
    if(mesh%hasBdyDomain)then
       do ie = mesh%ne_compute+1,mesh%ne_halo(1)
          if(mesh%edt_mark(ie).eq.1772.or.mesh%edt_mark(ie).eq.1771)then
             v1        = mesh%edt_v(1,ie)     ! 1st end of edge
             v2        = mesh%edt_v(2,ie)     ! 2nd end of edge
             do ilev = 1, nlev
                scalar1_normal_flux_at_edge_full_level(ilev, ie) = 0.5*(scalar1_height_at_pc_full_level(ilev,v1)+&
                                                                     scalar1_height_at_pc_full_level(ilev,v2))

                scalar2_normal_flux_at_edge_full_level(ilev, ie) = 0.5*(scalar2_height_at_pc_full_level(ilev,v1)+&
                                                                     scalar2_height_at_pc_full_level(ilev,v2))
             end do
          end if
       end do
    end if
#endif
#endif

    return
  end subroutine calc_primal_normal_flux_at_edge_var2

!================================================
! This edge-based routine can only be compute 
! until ne_compute, due to upwind, other wise 
! upwind-cell will be out of domain
!================================================

 subroutine calc_primal_normal_flux_at_edge(mesh,scalar_normal_velocity_at_edge_full_level, &
                                                 scalar_tangen_velocity_at_edge_full_level, &
                                                 scalar_height_at_pc_full_level           , &
                                                 scalar_normal_flux_at_edge_full_level    , &
                                                 adv_flag                                 , &
                                                 dtime                                    , &
                                                 nlev)
! io
   use omp_lib
   type(global_domain),   intent(in)    :: mesh
   real(r8),              intent(in)    :: scalar_normal_velocity_at_edge_full_level(nlev,mesh%ne_full)
   real(r8),              intent(in)    :: scalar_tangen_velocity_at_edge_full_level(nlev,mesh%ne_full)
   real(r8),              intent(in)    :: scalar_height_at_pc_full_level(nlev,mesh%nv_full)        ! q
   real(r8),              intent(inout) :: scalar_normal_flux_at_edge_full_level(nlev,mesh%ne_full) ! also q
   integer(i4),           intent(in)    :: adv_flag
   real(r8)   ,           intent(in)    :: dtime
   integer(i4),           intent(in)    :: nlev
! local
   real(r8)                             :: v1v2(3)
   real(r8)                             :: beta
   real(r8)                             :: wind_edge
   real(r8)                             :: der0, der1
   integer(i4)                          :: v1, v2, v_upwind
   integer(i4)                          :: ilev, ie
   integer(i4)                          :: flag, mark
   integer(i4)                          :: ii

   select case(adv_flag)

   case(2)  ! 2nd-order center difference
!$omp parallel  private(ie,ilev,v1,v2)
!$omp do schedule(static,100)
       do ie = 1, mesh%ne
          v1      = mesh%edt_v(1, ie)
          v2      = mesh%edt_v(2, ie)
          do ilev = 1, nlev
             scalar_normal_flux_at_edge_full_level(ilev, ie) = 0.5*(scalar_height_at_pc_full_level(ilev,v1)+&
                                                                    scalar_height_at_pc_full_level(ilev,v2))
          end do
       end do
!$omp end do nowait
!$omp end parallel 

   case(3,4,5,6,7)  ! WRF3/4, damping depends on beta value

    if(adv_flag.eq.3) beta=1._r8    ! 3rd order
    if(adv_flag.eq.4) beta=zero     ! 4th order
    if(adv_flag.eq.5) beta=0.25     ! hybrid 1
    if(adv_flag.eq.6) beta=0.5      ! hybrid 2
    if(adv_flag.eq.7) beta=0.75     ! hybrid 3

!$omp parallel  private(ii,ie,ilev,v1,v2,wind_edge,der0,der1) 
!$omp do schedule(dynamic,100) 
    do ii = 1, mesh%ne*nlev,1
        ie = ceiling(ii/real(nlev,r8))
        ilev=ii-(ie-1)*nlev

!    do ie = 1, mesh%ne
          v1        = mesh%edt_v(1, ie)
          v2        = mesh%edt_v(2, ie)
!          do ilev  = 1, nlev

             wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge_full_level(ilev,ie)
             call calc_der_2nd_at_hx(mesh,scalar_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v1), v1)),&
                                          scalar_height_at_pc_full_level(ilev,v1), &
                                          v1,ie,1,mesh%plg_stencil_number_2nd(v1),der0)

             call calc_der_2nd_at_hx(mesh,scalar_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v2), v2)),&
                                          scalar_height_at_pc_full_level(ilev,v2), &
                                          v2,ie,2,mesh%plg_stencil_number_2nd(v2),der1)

             call flux_ooo34(wind_edge, scalar_height_at_pc_full_level(ilev,v1),&
                                        scalar_height_at_pc_full_level(ilev,v2),&
                             der0,der1, scalar_normal_flux_at_edge_full_level(ilev,ie),&
                             real(mesh%edt_leng(ie),r8),beta)
!          end do
!    end do
    end do
!$omp end do nowait
!$omp end parallel

   case(33)  ! pure 3ord
!$omp parallel  private(ii,ie,ilev,v1,v2,wind_edge,der0,mark,v_upwind)
!$omp do schedule(dynamic,100)
    do ii = 1, mesh%ne*nlev,1
        ie=ceiling(ii/real(nlev,r8))
        ilev=ii-(ie-1)*nlev

!    do ie = 1, mesh%ne
          v1        = mesh%edt_v(1, ie)       ! 1st end of edge
          v2        = mesh%edt_v(2, ie)       ! 2nd end of edge
!          do ilev   = 1, nlev
             wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge_full_level(ilev,ie)
             mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)
             v_upwind  = mesh%edt_v(mark, ie)

             call calc_der_2nd_at_hx(mesh,scalar_height_at_pc_full_level(ilev,mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind), v_upwind)),&
                                          scalar_height_at_pc_full_level(ilev,v_upwind),&
                                          v_upwind,ie,mark,mesh%plg_stencil_number_2nd(v_upwind),der0)

             call flux_ooo33(scalar_height_at_pc_full_level(ilev,v1),&
                            scalar_height_at_pc_full_level(ilev,v2),&
                            der0, scalar_normal_flux_at_edge_full_level(ilev,ie),&
                            real(mesh%edt_leng(ie),r8))
!          end do
!    end do
    end do
!$omp end do nowait
!$omp end parallel

   case default
        print*,"you must select an advection scheme/flux divergence operator in gcm_para"
   end select

#ifdef OLD_MASS_PT 
! compute value at halo(1) for bdy-domain only
#ifdef LAM_DOMAIN
    if(mesh%hasBdyDomain)then
       do ie = mesh%ne_compute+1,mesh%ne_halo(1)
          if(mesh%edt_mark(ie).eq.1772.or.mesh%edt_mark(ie).eq.1771)then
             v1      = mesh%edt_v(1,ie)     ! 1st end of edge
             v2      = mesh%edt_v(2,ie)     ! 2nd end of edge
             do ilev = 1, nlev
                scalar_normal_flux_at_edge_full_level(ilev, ie) = 0.5*(scalar_height_at_pc_full_level(ilev,v1)+&
                                                                    scalar_height_at_pc_full_level(ilev,v2))
             end do
          end if
       end do
    end if
#endif
#endif

   return
  end subroutine calc_primal_normal_flux_at_edge

 subroutine calc_primal_normal_flux_at_edge_1d(mesh,scalar_normal_velocity_at_edge_full_level, &
                                                    scalar_tangen_velocity_at_edge_full_level, &
                                                    scalar_height_at_pc_full_level           , &
                                                    scalar_normal_flux_at_edge_full_level    , &
                                                    adv_flag                                )
! io
   use omp_lib
   type(global_domain),   intent(in)    :: mesh
   real(r8),              intent(in)    :: scalar_normal_velocity_at_edge_full_level(mesh%ne_full)
   real(r8),              intent(in)    :: scalar_tangen_velocity_at_edge_full_level(mesh%ne_full)
   real(r8),              intent(in)    :: scalar_height_at_pc_full_level(mesh%nv_full)        ! q
   real(r8),              intent(inout) :: scalar_normal_flux_at_edge_full_level(mesh%ne_full) ! also q
   integer(i4),           intent(in)    :: adv_flag
! local
   real(r8)                             :: v1v2(3)
   real(r8)                             :: beta
   real(r8)                             :: wind_edge
   real(r8)                             :: der0, der1
   integer(i4)                          :: v1, v2, v_upwind
   integer(i4)                          :: ilev, ie
   integer(i4)                          :: flag, mark
   integer(i4)                          :: ii

   select case(adv_flag)

   case(2)  ! 2nd-order center difference 
!$omp parallel  private(ie,ilev,v1,v2)   
!$omp do schedule(static,100) 
       do ie = 1, mesh%ne
          v1      = mesh%edt_v(1, ie)
          v2      = mesh%edt_v(2, ie)
          scalar_normal_flux_at_edge_full_level( ie) = 0.5*(scalar_height_at_pc_full_level(v1)+&
                                                            scalar_height_at_pc_full_level(v2))
       end do
!$omp end do nowait
!$omp end parallel 

   case(3,4,5,6,7)  ! WRF3/4, damping depends on beta value

    if(adv_flag.eq.3) beta=1._r8    ! 3rd order
    if(adv_flag.eq.4) beta=zero     ! 4th order
    if(adv_flag.eq.5) beta=0.25     ! hybrid 1
    if(adv_flag.eq.6) beta=0.5      ! hybrid 2
    if(adv_flag.eq.7) beta=0.75     ! hybrid 3

!$omp parallel  private(ii,ie,ilev,v1,v2,wind_edge,der0,der1) 
!$omp do schedule(dynamic,100) 
    do ie = 1, mesh%ne,1
          v1        = mesh%edt_v(1, ie)
          v2        = mesh%edt_v(2, ie)
             wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge_full_level(ie)
             call calc_der_2nd_at_hx(mesh,scalar_height_at_pc_full_level(mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v1), v1)),&
                                          scalar_height_at_pc_full_level(v1), &
                                          v1,ie,1,mesh%plg_stencil_number_2nd(v1),der0)

             call calc_der_2nd_at_hx(mesh,scalar_height_at_pc_full_level(mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v2), v2)),&
                                          scalar_height_at_pc_full_level(v2), &
                                          v2,ie,2,mesh%plg_stencil_number_2nd(v2),der1)

             call flux_ooo34(wind_edge, scalar_height_at_pc_full_level(v1),&
                                        scalar_height_at_pc_full_level(v2),&
                             der0,der1, scalar_normal_flux_at_edge_full_level(ie),&
                             real(mesh%edt_leng(ie),r8),beta)
    end do
!$omp end do nowait
!$omp end parallel 

   case(33)  ! pure 3ord 

!$omp parallel  private(ii,ie,ilev,v1,v2,wind_edge,der0,mark,v_upwind)    
!$omp do schedule(dynamic,100) 
    do ie = 1, mesh%ne,1

!    do ie = 1, mesh%ne
          v1        = mesh%edt_v(1, ie)       ! 1st end of edge
          v2        = mesh%edt_v(2, ie)       ! 2nd end of edge
!          do ilev   = 1, nlev
             wind_edge = mesh%edt_edpNr_edtTg(ie)*scalar_normal_velocity_at_edge_full_level(ie)
             mark      = int(-0.5*sign(1._r8,wind_edge)+1.5)
             v_upwind  = mesh%edt_v(mark, ie)
  
             call calc_der_2nd_at_hx(mesh,scalar_height_at_pc_full_level(mesh%plg_stencil_index_2nd(1:mesh%plg_stencil_number_2nd(v_upwind), v_upwind)),&
                                          scalar_height_at_pc_full_level(v_upwind),&
                                          v_upwind,ie,mark,mesh%plg_stencil_number_2nd(v_upwind),der0)

             call flux_ooo33(scalar_height_at_pc_full_level(v1),&
                            scalar_height_at_pc_full_level(v2),&
                            der0, scalar_normal_flux_at_edge_full_level(ie),&
                            real(mesh%edt_leng(ie),r8))
!          end do
!    end do
    end do
!$omp end do nowait
!$omp end parallel 

   case default
        print*,"you must select an advection scheme/flux divergence operator in gcm_para"
   end select

   return
  end subroutine calc_primal_normal_flux_at_edge_1d

!================================================
!
!              PRIVATE SUBROUTINES
!
!================================================

  subroutine calc_der_2nd_at_hx(mesh,s_array,s0,iv,ie,flag,nlength,der)
! io
    type(global_domain),  intent(in) :: mesh
    real(r8),             intent(in) :: s_array(:)
    real(r8),             intent(in) :: s0
    integer(i4)         , intent(in) :: iv
    integer(i4)         , intent(in) :: ie
    integer(i4)         , intent(in) :: flag
    integer(i4)         , intent(in) :: nlength
    real(r8)            , intent(out):: der
! local
    real(r8)                         :: m_array(1:nlength)
    real(r8)                         :: f_array

      m_array(:) = s_array(:)-s0
      !f_array = dot_product(mesh%edt(ie)%v_weight_prime_cell(flag,1:mesh%plg(iv)%stencil_number_2nd), m_array(:))
      f_array = dot_product(mesh%edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(iv), flag, ie), m_array(:))
      der     = 2._r8*f_array

      return
   end subroutine calc_der_2nd_at_hx

!================================================
! same as calc_der_2nd_at_hx but for 2 vars
!================================================

   subroutine calc_der_2nd_at_hx_var2(mesh,s1_array,s2_array,s1_0,s2_0,iv,ie,flag,nlength,der1,der2)
! io
    type(global_domain),  intent(in) :: mesh
    real(r8),             intent(in) :: s1_array(:)
    real(r8),             intent(in) :: s2_array(:)
    real(r8),             intent(in) :: s1_0
    real(r8),             intent(in) :: s2_0
    integer(i4)         , intent(in) :: iv
    integer(i4)         , intent(in) :: ie
    integer(i4)         , intent(in) :: flag
    integer(i4)         , intent(in) :: nlength
    real(r8)            , intent(out):: der1
    real(r8)            , intent(out):: der2
! local
    real(r8)                         :: m1_array(1:nlength)
    real(r8)                         :: m2_array(1:nlength)
    real(r8)                         :: f1_array
    real(r8)                         :: f2_array

      m1_array(:) = s1_array(:)-s1_0
      m2_array(:) = s2_array(:)-s2_0
      f1_array    = dot_product(mesh%edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(iv), flag, ie), m1_array(:))
      f2_array    = dot_product(mesh%edt_v_weight_prime_cell(1:mesh%plg_stencil_number_2nd(iv), flag, ie), m2_array(:))
      der1        = 2._r8*f1_array
      der2        = 2._r8*f2_array

      return
   end subroutine calc_der_2nd_at_hx_var2

  end module grist_dycore_primal_flux_operators_2d
