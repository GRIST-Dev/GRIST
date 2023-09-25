
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: A reference atmosphere profile used for HPGF the part that 
!              satisfy the hydrostatic balance is removed to enhance the 
!              accuracy. It does not need to modify the prognostic equation,
!              just yet another HPGF formulation
!----------------------------------------------------------------------------

 module grist_dycore_ref_atmos

! para
  use grist_constants,      only: i4, r8, rearth, rdry, cp, one
  use grist_domain_types  , only: global_domain
  use grist_nml_module,     only: nlev,profile_type,hpgf_pfull_alpha
! vert
  use grist_hpe_constants,  only: p0, eta_full_a, eta_full_b, eta_face_a, eta_face_b
! dycore
  use grist_dycore_vars_module, only: dycoreVarSurface, dycoreVarCellFull

  implicit none
  private
  public :: grist_dycore_ref_atmos_init,   &
            grist_dycore_ref_atmos_final,  &
            grist_dycore_ref_atmos_create_ptb
! hold
     real(r8), public, allocatable  :: scalar_hpressure_at_pc_full_level_bar(:,:)
     real(r8), public, allocatable  ::    scalar_alphad_at_pc_full_level_bar(:,:)
     real(r8), public, allocatable  ::    scalar_masspt_at_pc_full_level_bar(:,:)
     real(r8), public, allocatable  ::      scalar_geop_at_pc_full_level_bar(:,:)
     real(r8), public, allocatable  ::      scalar_geop_at_pc_face_level_bar(:,:)
     real(r8), public, allocatable  :: scalar_grad_hpres_at_edge_full_level_bar(:,:)
! perturbed var used for PGF stored here
     real(r8), public, allocatable  :: scalar_pressure_at_pc_full_level_ptb(:,:)
     real(r8), public, allocatable  ::   scalar_alphad_at_pc_full_level_ptb(:,:)
     real(r8), public, allocatable  ::     scalar_geop_at_pc_full_level_ptb(:,:)
     real(r8), public, allocatable  ::     scalar_delp_at_pc_full_level_ptb(:,:)
! constant
     real(r8), parameter            :: at = 0.09923_r8
     real(r8), parameter            :: bt = 247.7874_r8
     real(r8), parameter            :: ct = -1.0385_r8

  CONTAINS

! called after rkfb init
     subroutine grist_dycore_ref_atmos_init(mesh)
! io
       type(global_domain), intent(in) :: mesh
! local
       real(r8)             :: hps0, hps1, pres
       real(r8)             :: temp, hp_face_up, hp_face_lo, vert_sum 
       integer(i4)          :: iv, ie
       integer(i4)          :: icell1, icell2 
       integer(i4)          :: ilev, ilevel
   
       hps0  = 101300._r8

       if(.not.allocated(scalar_hpressure_at_pc_full_level_bar))     allocate(scalar_hpressure_at_pc_full_level_bar(nlev,mesh%nv))
       if(.not.allocated(   scalar_alphad_at_pc_full_level_bar))        allocate(scalar_alphad_at_pc_full_level_bar(nlev,mesh%nv))
       if(.not.allocated(   scalar_masspt_at_pc_full_level_bar))        allocate(scalar_masspt_at_pc_full_level_bar(nlev,mesh%nv))
       if(.not.allocated(     scalar_geop_at_pc_full_level_bar))          allocate(scalar_geop_at_pc_full_level_bar(nlev,mesh%nv))
       if(.not.allocated(     scalar_geop_at_pc_face_level_bar))          allocate(scalar_geop_at_pc_face_level_bar(nlev+1,mesh%nv))
       if(.not.allocated(scalar_grad_hpres_at_edge_full_level_bar)) allocate(scalar_grad_hpres_at_edge_full_level_bar(nlev,mesh%ne))

       if(.not.allocated(scalar_pressure_at_pc_full_level_ptb))  allocate(scalar_pressure_at_pc_full_level_ptb(nlev,mesh%nv))
       if(.not.allocated(  scalar_alphad_at_pc_full_level_ptb))  allocate(  scalar_alphad_at_pc_full_level_ptb(nlev,mesh%nv))
       if(.not.allocated(    scalar_geop_at_pc_full_level_ptb))  allocate(    scalar_geop_at_pc_full_level_ptb(nlev,mesh%nv))
       if(.not.allocated(    scalar_delp_at_pc_full_level_ptb))  allocate(    scalar_delp_at_pc_full_level_ptb(nlev,mesh%nv))

! construct hpressure, alphad and geop full
       select case (trim(profile_type))
       case('reg') ! this profile actually recovers the unmodified HPGF, but provides as a regression
                   ! this is because in this config, both grad of scalar_hpressure_at_pc_full_level_bar and
                   ! scalar_geop_at_pc_full_level_bar are zero, leading to original unperturbed HPGF
                   ! it is deliberately wrong such that the model will blow if it works
       do iv = 1, mesh%nv
          do ilev = 1, nlev
             pres = eta_full_a(ilev)*p0+eta_full_b(ilev)*hps0
             scalar_hpressure_at_pc_full_level_bar(ilev,iv) = pres
                  scalar_geop_at_pc_full_level_bar(ilev,iv) =-rdry*(at*(pres-hps0)+bt/(1._r8+ct)*exp((1._r8+ct)*log(pres/hps0)))
                scalar_alphad_at_pc_full_level_bar(ilev,iv) = rdry*(at+bt*exp(ct*log(pres)))
          end do
       end do

       case('tp1')
! construct hpressure, alphad and geop full
       do iv = 1, mesh%nv
          hps1 = dycoreVarSurface%scalar_hpressure_n%f(iv) ! Pa
          do ilev = 1, nlev
             pres = eta_full_a(ilev)*p0+eta_full_b(ilev)*hps1 ! Pa
             scalar_hpressure_at_pc_full_level_bar(ilev,iv) = pres ! Pa
                  scalar_geop_at_pc_full_level_bar(ilev,iv) =-rdry*(at*(pres/100._r8-hps0/100.)+bt/(1._r8+ct)*((pres/100._r8)**(1._r8+ct)-(hps0/100.)**(1._r8+ct)))
                scalar_alphad_at_pc_full_level_bar(ilev,iv) = rdry*(at+bt*(pres/100)**ct)/100._r8
! for mass pt
             temp = pres*(at+bt*(pres/100)**ct)/100._r8
             hp_face_lo = eta_face_a(ilev+1)*p0+eta_face_b(ilev+1)*hps1
             hp_face_up = eta_face_a(ilev)  *p0+eta_face_b(ilev)  *hps1
             scalar_masspt_at_pc_full_level_bar(ilev,iv) = temp/(pres/p0)**(rdry/cp)*(hp_face_lo-hp_face_up)
          end do
          do ilev = 1, nlev+1
             pres = eta_face_a(ilev)*p0+eta_face_b(ilev)*hps1 ! Pa
             scalar_geop_at_pc_face_level_bar(ilev,iv) =-rdry*(at*(pres/100._r8-hps0/100.)+bt/(1._r8+ct)*((pres/100._r8)**(1._r8+ct)-(hps0/100.)**(1._r8+ct)))
          end do
       end do

       case('tp2')
! construct hpressure, alphad and geop full
       do iv = 1, mesh%nv
          do ilev = 1, nlev
             hps1 = hps0 !dycoreVarSurface%scalar_hpressure_n%f(iv) ! Pa
             pres = eta_full_a(ilev)*p0+eta_full_b(ilev)*hps1 ! Pa
             scalar_hpressure_at_pc_full_level_bar(ilev,iv) = pres ! Pa
                  scalar_geop_at_pc_full_level_bar(ilev,iv) =-rdry*(at*(pres/100._r8-hps0/100.)+bt/(1._r8+ct)*((pres/100._r8)**(1._r8+ct)-(hps0/100.)**(1._r8+ct)))
                scalar_alphad_at_pc_full_level_bar(ilev,iv) = rdry*(at+bt*(pres/100)**ct)/100._r8
          end do
       end do

       case('tp3')
! construct hpressure, alphad and geop full
       do iv = 1, mesh%nv
          hps1 = dycoreVarSurface%scalar_hpressure_n%f(iv) ! Pa
          do ilev = 1, nlev
             pres = eta_full_a(ilev)*p0+eta_full_b(ilev)*hps1 ! Pa
             scalar_hpressure_at_pc_full_level_bar(ilev,iv) = pres ! Pa
                scalar_alphad_at_pc_full_level_bar(ilev,iv) = rdry*(at+bt*(pres/100)**ct)/100._r8
          end do
! use vertical-integral to obtain phi at face
          do ilevel = 1, nlev     ! for each face ilevel except surface
             vert_sum = 0._r8
             do ilev = ilevel, nlev ! for each full level ilev
                pres       = eta_full_a(ilev)*p0+eta_full_b(ilev)*hps1 ! Pa
                temp       = pres*(at+bt*(pres/100)**ct)/100._r8
                hp_face_lo = eta_face_a(ilev+1)*p0+eta_face_b(ilev+1)*hps1
                hp_face_up = eta_face_a(ilev)  *p0+eta_face_b(ilev)  *hps1
                vert_sum   = vert_sum+rdry*temp*log(hp_face_lo/hp_face_up)
             end do
             scalar_geop_at_pc_face_level_bar(ilevel,iv) = dycoreVarSurface%scalar_geopotential_n%f(iv)+vert_sum
          end do
          scalar_geop_at_pc_face_level_bar(nlev+1,iv) = dycoreVarSurface%scalar_geopotential_n%f(iv)
       end do
       scalar_geop_at_pc_full_level_bar(1:nlev,:) = (scalar_geop_at_pc_face_level_bar(1:nlev,:)+&
                                                     scalar_geop_at_pc_face_level_bar(2:nlev+1,:))*0.5_r8
       case('default')
           print*,"you must select a reference profile type for profile_type in dycore_para"
       end select

! evaluate grad_hpres at edge
       do ie = 1, mesh%ne
          icell1 = mesh%edt_v(1,ie)
          icell2 = mesh%edt_v(2,ie)
          do ilev = 1, nlev
             scalar_grad_hpres_at_edge_full_level_bar(ilev,ie) = (scalar_hpressure_at_pc_full_level_bar(ilev,icell2)-&
                                                                  scalar_hpressure_at_pc_full_level_bar(ilev,icell1))/&
                                                                  (rearth*mesh%edt_leng(ie))
          end do
       end do

       return
     end subroutine grist_dycore_ref_atmos_init

     subroutine grist_dycore_ref_atmos_create_ptb

!        scalar_pressure_at_pc_full_level_ptb =     dycoreVarCellFull%scalar_pressure_np1%f-scalar_hpressure_at_pc_full_level_bar
! one recover the backward case
        scalar_pressure_at_pc_full_level_ptb =     (hpgf_pfull_alpha *dycoreVarCellFull%scalar_pressure_np1%f+&
                                               (one-hpgf_pfull_alpha)*dycoreVarCellFull%scalar_pressure_rk%f)-scalar_hpressure_at_pc_full_level_bar
          scalar_alphad_at_pc_full_level_ptb =        dycoreVarCellFull%scalar_alpha_np1%f-   scalar_alphad_at_pc_full_level_bar
            scalar_geop_at_pc_full_level_ptb = dycoreVarCellFull%scalar_geopotential_n%f  -     scalar_geop_at_pc_full_level_bar
           scalar_delp_at_pc_full_level_ptb  =         dycoreVarCellFull%scalar_delp_np1%f-    dycoreVarCellFull%scalar_delhp_np1%f 

     end subroutine grist_dycore_ref_atmos_create_ptb

     subroutine grist_dycore_ref_atmos_final

       if(allocated(scalar_hpressure_at_pc_full_level_bar))     deallocate(scalar_hpressure_at_pc_full_level_bar)
       if(allocated(   scalar_alphad_at_pc_full_level_bar))        deallocate(scalar_alphad_at_pc_full_level_bar)
       if(allocated(   scalar_masspt_at_pc_full_level_bar))        deallocate(scalar_masspt_at_pc_full_level_bar)
       if(allocated(     scalar_geop_at_pc_full_level_bar))          deallocate(scalar_geop_at_pc_full_level_bar)
       if(allocated(     scalar_geop_at_pc_face_level_bar))          deallocate(scalar_geop_at_pc_face_level_bar)
       if(allocated(scalar_grad_hpres_at_edge_full_level_bar))  deallocate(scalar_grad_hpres_at_edge_full_level_bar)

       if(allocated(scalar_pressure_at_pc_full_level_ptb))  deallocate(scalar_pressure_at_pc_full_level_ptb)
       if(allocated(  scalar_alphad_at_pc_full_level_ptb))  deallocate(  scalar_alphad_at_pc_full_level_ptb)
       if(allocated(    scalar_geop_at_pc_full_level_ptb))  deallocate(    scalar_geop_at_pc_full_level_ptb)
       if(allocated(    scalar_delp_at_pc_full_level_ptb))  deallocate(    scalar_delp_at_pc_full_level_ptb)

     end subroutine grist_dycore_ref_atmos_final

  end module grist_dycore_ref_atmos
