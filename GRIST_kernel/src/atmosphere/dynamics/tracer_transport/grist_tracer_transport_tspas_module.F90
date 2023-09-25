
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: An efficiently optimized implemetation of TSPAS; Only calculate 
!              flux using tspas, no div operator same as before to within 
!              machine truncation error; TSPAS should not used with FCT limiter 
!              because it directly produces shape-preserving flux
! Revision history:
!----------------------------------------------------------------------------

 module grist_tracer_transport_tspas_module

  use grist_constants,        only: i4, r8, rearth
  use grist_domain_types,     only: global_domain

  implicit none
  private
  public   :: tracer_transport_hori_tspas_flux, &
              tracer_transport_vert_tspas_flux

  real(r8), parameter  :: one  = 1._r8
#ifndef SPCODE
  real(r8), parameter  :: kesi = 1.e-80_r8
#else
  real(r8), parameter  :: kesi = 1.e-40_r8
#endif
  real(r8), parameter  :: zero = 0._r8

  contains

  subroutine tracer_transport_hori_tspas_flux(mesh,scalar_normal_velocity_at_edge_full_level ,& ! u
                                              scalar_normal_mass_flux_at_edge_full_level,& ! u*delhp
                                              scalar_tracer_mass_at_pc_full_level       ,& ! delhp*q at n
                                              scalar_tracer_mxrt_at_pc_full_level       ,& ! q at n
                                              scalar_normal_mxrt_at_edge_full_level     ,& ! q
                                              dtime, ntracer, nlev                     ) 
   implicit none
! io
   type(global_domain),   intent(in)    :: mesh
   real(r8), allocatable, intent(in)    :: scalar_normal_velocity_at_edge_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar_normal_mass_flux_at_edge_full_level(:,:)
   real(r8), allocatable, intent(in)    :: scalar_tracer_mass_at_pc_full_level(:,:,:)
   real(r8), allocatable, intent(in)    :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
   real(r8), allocatable, intent(inout) :: scalar_normal_mxrt_at_edge_full_level(:,:,:)
   real(r8),              intent(in)    :: dtime
   integer(i4)          , intent(in)    :: ntracer, nlev
! local 
   real(r8), allocatable    :: avalue(:,:,:)
   real(r8), allocatable    :: lw_flux(:,:,:)
   real(r8)                 :: scalar_value_astk(ntracer,nlev)
   real(r8)                 :: div_sum(ntracer,nlev)
   real(r8)                 :: scalar_max(ntracer,nlev), scalar_min(ntracer,nlev)
   real(r8)                 :: gama_max(nlev)
   real(r8)                 :: beta(nlev)
   real(r8)                 :: v1v2(1:3)
   real(r8)                 :: wind_edge, c_value, cc_value
   real(r8)                 :: gama
   real(r8)                 :: phi0, phi1
   integer(i4)              :: v1,v2 
   integer(i4)              :: flag
   integer(i4)              :: index_edge, index_cell
   integer                  :: inb,kk,itracer, iv, ie, ilev

   if(.not.allocated(avalue))      allocate(avalue(ntracer,nlev ,mesh%nv_halo(1)))
   if(.not.allocated(lw_flux))     allocate(lw_flux(ntracer,nlev,mesh%ne_halo(1)))

   kk      = 3   ! empirical parameter

! Compute LW flux and store at edges
   do ie = 1, mesh%ne_halo(1)
        !v1          = mesh%edt(ie)%v(1)
        !v2          = mesh%edt(ie)%v(2)
        v1          = mesh%edt_v(1,ie)
        v2          = mesh%edt_v(2,ie)
        v1v2        = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
        flag        = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))
        do ilev = 1, nlev
           wind_edge  = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
           do itracer = 1, ntracer
              phi0 = scalar_tracer_mass_at_pc_full_level(itracer,ilev,v1)
              phi1 = scalar_tracer_mass_at_pc_full_level(itracer,ilev,v2)
              lw_flux(itracer,ilev,ie) = (0.5_r8*(phi0+phi1)-0.5_r8*sign(1._r8,wind_edge)*(phi1-phi0)*&
                                         abs(wind_edge)*dtime/(rearth*(mesh%edt_leng(ie))))*&
                                         scalar_normal_velocity_at_edge_full_level(ilev,ie)
           end do
        end do
   end do

 ! First integration
 ! for each vertice
 ! get min and max around this vertice
 ! and PRECALCULATE the updated value using LW

   do iv=1, mesh%nv_halo(1)     ! for each vertice in the mesh
!
! 1st loop over nbs for gams_compo 
!
#ifndef SPCODE
      gama_max  = -1.e50_r8
#else
      gama_max  = -1.e38_r8
#endif
      do ie = 1, mesh%vtx_nnb(iv)
         index_edge = mesh%vtx_ed(ie,iv)
         do ilev = 1, nlev
            wind_edge  = scalar_normal_velocity_at_edge_full_level(ilev,index_edge)*mesh%plg_nr(ie,iv)  ! to ensure outward positive
            gama       = abs(wind_edge)*(1._r8-abs(wind_edge)*dtime/(rearth*mesh%edt_leng(index_edge)))*&
                                  rearth*mesh%edp_leng(index_edge)
            gama_max(ilev) = max(gama_max(ilev),gama)
         end do
      end do
! gama max
      do ilev = 1, nlev
         beta(ilev)  = 2._r8/(2._r8-kk*dtime*gama_max(ilev)/((rearth**2)*mesh%plg_areag(iv)))
         beta(ilev)  = max(beta(ilev),1._r8)
      end do
!
! 2nd loop over nbs for pre-update 
!
      div_sum(:,:) = 0._r8
      do ie = 1, mesh%vtx_nnb(iv)
         index_edge = mesh%vtx_ed(ie,iv)
         do ilev = 1, nlev
            do itracer = 1, ntracer
               div_sum(itracer,ilev) = div_sum(itracer,ilev)+beta(ilev)*lw_flux(itracer,ilev,index_edge)*(mesh%plg_nr(ie,iv))*&
                                                             (rearth*mesh%edp_leng(mesh%vtx_ed(ie,iv)))
            end do
         end do
      end do
      scalar_value_astk(:,:) = scalar_tracer_mass_at_pc_full_level(:,:,iv)-1._r8*dtime*div_sum(:,:)/((rearth**2)*mesh%plg_areag(iv))
!
! 3rd loop over nbs
!
      scalar_max       = scalar_tracer_mass_at_pc_full_level(:,:,iv)
      scalar_min       = scalar_tracer_mass_at_pc_full_level(:,:,iv)
      do inb = 1, mesh%vtx_nnb(iv)
         index_cell = mesh%vtx_nb(inb,iv)
         do ilev = 1, nlev
            do itracer = 1, ntracer
               scalar_max(itracer,ilev)  = max(scalar_max(itracer,ilev),scalar_tracer_mass_at_pc_full_level(itracer,ilev,index_cell))
               scalar_min(itracer,ilev)  = min(scalar_min(itracer,ilev),scalar_tracer_mass_at_pc_full_level(itracer,ilev,index_cell))
            end do
         end do
      end do

      avalue(:,:,iv)     = (scalar_value_astk-scalar_max)*(scalar_value_astk-scalar_min)

   end do
 ! End of First integration

! unchanged above, check unsigned reconstructed scalar at edge (q)

   do ie = 1, mesh%ne_compute
      !v1   = mesh%edt(ie)%v(1)
      !v2   = mesh%edt(ie)%v(2)
      v1   = mesh%edt_v(1,ie)
      v2   = mesh%edt_v(2,ie)
      v1v2 = mesh%vtx_p(1:3,v2)-mesh%vtx_p(1:3,v1)
      flag = sign(one,dot_product(v1v2,real(mesh%edp_nr(1:3,ie),r8)))  ! actually always 1, but we explicitly code here
      do ilev = 1, nlev
         wind_edge = flag*scalar_normal_velocity_at_edge_full_level(ilev,ie)
         do itracer = 1, ntracer
            c_value = 1._r8
            if(avalue(itracer,ilev,v1).lt.0._r8.and.avalue(itracer,ilev,v2).lt.0._r8) c_value = 0._r8
            cc_value  = c_value+(one-c_value)*abs(wind_edge)*dtime/(rearth*mesh%edt_leng(ie))
            phi0 = scalar_tracer_mass_at_pc_full_level(itracer,ilev,v1)
            phi1 = scalar_tracer_mass_at_pc_full_level(itracer,ilev,v2)
            scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = wind_edge*(0.5_r8*(phi0+phi1)-0.5_r8*sign(1._r8,wind_edge)*(phi1-phi0)*cc_value)
! this line is to generate q such that all operators have the same output style,
! however, division may potentially lead to nan when wind_edge(normal flux)=0
! since we do not want change the interface of other operators,
! an ad-hoc approach is to add a small parameter, although I donot like it very
! much.
            scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie) = scalar_normal_mxrt_at_edge_full_level(itracer,ilev,ie)/&
                                                                     (scalar_normal_mass_flux_at_edge_full_level(ilev,ie)+kesi)
         end do
      end do
   end do

   deallocate(avalue)
   deallocate(lw_flux)
  
   return
  end subroutine tracer_transport_hori_tspas_flux

  subroutine tracer_transport_vert_tspas_flux(dtime                                 , &
                                              scalar_mass_eta_velocity_at_face_level, &
                                              scalar_tracer_mass_at_full_level_n    , &
                                              scalar_tracer_mxrt_at_full_level_n    , &
                                              scalar_delhp_at_full_level_np1        , &
                                              scalar_tspas_flux_at_face_level       , & ! limited q*m*etadot
                                              nlev, nlevp)
! io
     real(r8), intent(in)    :: dtime
     real(r8), intent(in)    :: scalar_mass_eta_velocity_at_face_level(nlevp) ! downward positive
     real(r8), intent(in)    :: scalar_tracer_mass_at_full_level_n(nlev)
     real(r8), intent(in)    :: scalar_tracer_mxrt_at_full_level_n(nlev)
     real(r8), intent(in)    :: scalar_delhp_at_full_level_np1(nlev)
     real(r8), intent(inout) :: scalar_tspas_flux_at_face_level(nlevp)        ! limited delhp*q
     integer(i4),           intent(in)    :: nlev, nlevp
! local 
     real(r8)       :: scalar_lw_flux_at_face_level(nlevp)
     real(r8)       :: scalar_delhp_at_face_level_np1(nlevp)
     real(r8)       :: avalue(nlev)
     real(r8)       :: gama_up,gama_lo,beta,qstar,qmax, qmin
     real(r8)       :: c_value, cc_value
     integer(i4)    :: ilev
!
! (1) compute LW-based flux using m*etadot*delhp*q (so final update should be scaled by delhp)
!     for tspas, tendency generator and state updator use the same quantity
!
      scalar_lw_flux_at_face_level(1)       = zero
      scalar_lw_flux_at_face_level(nlevp)   = zero
      scalar_delhp_at_face_level_np1(1)     = 0.5_r8*scalar_delhp_at_full_level_np1(1)
      scalar_delhp_at_face_level_np1(nlevp) = 0.5_r8*scalar_delhp_at_full_level_np1(nlev)

      do ilev = 2, nlev
         scalar_lw_flux_at_face_level(ilev) = 0.5_r8*(scalar_tracer_mass_at_full_level_n(ilev)+scalar_tracer_mass_at_full_level_n(ilev-1))-&
                                              0.5_r8*sign(one,scalar_mass_eta_velocity_at_face_level(ilev))*&
                                             (scalar_tracer_mass_at_full_level_n(ilev)-scalar_tracer_mass_at_full_level_n(ilev-1))*&
                                              abs(scalar_mass_eta_velocity_at_face_level(ilev))*dtime/scalar_delhp_at_face_level_np1(ilev)
         scalar_delhp_at_face_level_np1(ilev) = 0.5_r8*scalar_delhp_at_full_level_np1(ilev-1)+0.5_r8*scalar_delhp_at_full_level_np1(ilev)
      end do

! (2) compute beta, preupdate using lw and obtain avalue
      do ilev = 1, nlev ! full
         gama_up =     abs(scalar_mass_eta_velocity_at_face_level(ilev))*dtime/scalar_delhp_at_face_level_np1(ilev)*&
                  (one-abs(scalar_mass_eta_velocity_at_face_level(ilev))*dtime/scalar_delhp_at_face_level_np1(ilev))
         gama_lo =     abs(scalar_mass_eta_velocity_at_face_level(ilev+1))*dtime/scalar_delhp_at_face_level_np1(ilev+1)*&
                  (one-abs(scalar_mass_eta_velocity_at_face_level(ilev+1))*dtime/scalar_delhp_at_face_level_np1(ilev+1))
         beta    = max(one,gama_up,gama_lo)
         qstar   = scalar_tracer_mass_at_full_level_n(ilev)+beta*dtime/scalar_delhp_at_full_level_np1(ilev)*&
                   (scalar_mass_eta_velocity_at_face_level(ilev+1)*scalar_lw_flux_at_face_level(ilev+1)-&
                    scalar_mass_eta_velocity_at_face_level(ilev)  *scalar_lw_flux_at_face_level(ilev))
         if(ilev.eq.1)then
            qmax    = max(scalar_tracer_mass_at_full_level_n(ilev),scalar_tracer_mass_at_full_level_n(ilev+1))
            qmin    = min(scalar_tracer_mass_at_full_level_n(ilev),scalar_tracer_mass_at_full_level_n(ilev+1))
         else if(ilev.eq.nlev)then
            qmax    = max(scalar_tracer_mass_at_full_level_n(ilev),scalar_tracer_mass_at_full_level_n(ilev-1))
            qmin    = min(scalar_tracer_mass_at_full_level_n(ilev),scalar_tracer_mass_at_full_level_n(ilev-1))
         else
            qmax    = max(scalar_tracer_mass_at_full_level_n(ilev),scalar_tracer_mass_at_full_level_n(ilev-1),scalar_tracer_mass_at_full_level_n(ilev+1))
            qmin    = min(scalar_tracer_mass_at_full_level_n(ilev),scalar_tracer_mass_at_full_level_n(ilev-1),scalar_tracer_mass_at_full_level_n(ilev+1))
         end if
         avalue(ilev) = (qstar-qmax)*(qstar-qmin)
      end do
! (3) compute final flux, do not update
      do ilev = 2, nlev
         c_value = 1._r8
         if(avalue(ilev).lt.0._r8.and.avalue(ilev-1).lt.0._r8) c_value = 0._r8
         cc_value  = c_value+(one-c_value)*abs(scalar_mass_eta_velocity_at_face_level(ilev))*dtime/scalar_delhp_at_face_level_np1(ilev)
         scalar_tspas_flux_at_face_level(ilev) = 0.5_r8*(scalar_tracer_mass_at_full_level_n(ilev)+scalar_tracer_mass_at_full_level_n(ilev-1))-&
                                                 0.5_r8*sign(one,scalar_mass_eta_velocity_at_face_level(ilev))*&
                                                (scalar_tracer_mass_at_full_level_n(ilev)-scalar_tracer_mass_at_full_level_n(ilev-1))*cc_value
      end do
 
      return
   end subroutine tracer_transport_vert_tspas_flux

  end module grist_tracer_transport_tspas_module
