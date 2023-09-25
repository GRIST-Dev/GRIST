
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: This module contains major tracer variables used for online 
!              and offline run in the 3D model.
!              Name Convention:
!                 scalar_${varname}_at_$(location)_${time}
!                 tend_${varname}_at_${location}_${where_from}
!                 Similar to dyn_vars_module but independent data
! Revision history:
!                1. 2022, location is defined before in type, ${hori}_${vertical}
!                   hori = cell, edge, vert; vertical = face, full
!----------------------------------------------------------------------------

  module grist_tracer_transport_vars_module

    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_1d_field, scalar_2d_field, scalar_3d_field
    use grist_constants,    only: i4, r8, zero, one, half, eps
    use grist_nml_module,   only: nlev, ntracer, nmif, mif_index, working_mode, isolate_tracer_test
    use grist_hpe_constants,only: deta_full, deta_face
    use grist_mpi

    implicit none

    save
    PUBLIC

!================================================
! edge, full level
!================================================
    TYPE TRACER_VAR_EDGE_FULL

    type(scalar_2d_field)   ::  scalar_normal_mass_flux_avg_adv  ! during ndyn steps withing an adv step averaged
    type(scalar_2d_field)   ::  scalar_normal_velocity_avg_adv 
    type(scalar_2d_field)   ::  scalar_mif_n                     ! at the begining of dycore/tracer step

    END TYPE TRACER_VAR_EDGE_FULL
!================================================
! primal cell, full level
!================================================
    TYPE TRACER_VAR_CELL_FULL

    type(scalar_2d_field)   ::  scalar_delhp_avg_adv    ! dyn step averaged
    type(scalar_2d_field)   ::  scalar_delhp_end_adv    ! delhp state at the end of this adv step
    type(scalar_2d_field)   ::  tend_mass_div_avg_adv
    type(scalar_3d_field)   ::  scalar_tracer_mxrt_n    ! mixing ratio, diagnosed
    type(scalar_3d_field)   ::  scalar_tracer_mxrt_np1  ! mixing ratio, diagnosed
    type(scalar_3d_field)   ::  scalar_tracer_mass_n    ! mass, prognostic
    type(scalar_3d_field)   ::  scalar_tracer_mass_np1  ! mass, prognostic
    type(scalar_3d_field)   ::  tend_tracer_mass_hori   ! hori tend of tracers
    type(scalar_3d_field)   ::  tend_tracer_mass_vert   ! vert tend of tracers
    type(scalar_2d_field)   ::  scalar_mif_n            ! directly computed from tracer module
    type(scalar_3d_field)   ::  tend_tracer_mass_laplacian_2nd
! add for tiedtke conv, 3d adv tend for qv
    type(scalar_2d_field)   ::  tend_qv_n

    END TYPE TRACER_VAR_CELL_FULL

    character(len=50), allocatable :: scalar_tracer_name(:)
    real(r8),          allocatable :: scalar_tracer_mxrt_min(:)
!================================================
! primal cell, face level
!================================================
    TYPE TRACER_VAR_CELL_FACE

    type(scalar_2d_field)   :: scalar_eta_mass_flux_avg_adv ! diagnosed from scalar_delhp_at_pc_full_level_avg
    type(scalar_2d_field)   :: scalar_mif_n                 ! at the begining of dycore/tracer step
                                                            ! interpolated from tracerVarCellFull%scalar_mif_n
    END TYPE TRACER_VAR_CELL_FACE

    TYPE TRACER_VAR_SURFACE

    type(scalar_1d_field)   :: tend_hpressure

    END TYPE TRACER_VAR_SURFACE

    type(TRACER_VAR_CELL_FULL) :: tracerVarCellFull
    type(TRACER_VAR_CELL_FACE) :: tracerVarCellFace
    type(TRACER_VAR_EDGE_FULL) :: tracerVarEdgeFull
    type(TRACER_VAR_SURFACE)   :: tracerVarSurface

  CONTAINS

   subroutine tracer_transport_vars_construct(mesh)
! io
     type(global_domain), intent(in) :: mesh
! local
     integer(i4)     ::  it, ie, iv

      if(ntracer.le.0.or.ntracer.gt.100)then
        if(mpi_rank().eq.0) print*,"invalid ntracer value, please reset, model aborts"
        stop
      end if

!================================================
! edge, full level
!================================================
      allocate(tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(nlev,mesh%ne))
      allocate(tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f(nlev,mesh%ne))
      allocate(tracerVarEdgeFull%scalar_mif_n%f(nlev,mesh%ne))
      tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%pos = 6 ; tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f = zero  
      tracerVarEdgeFull%scalar_normal_velocity_avg_adv%pos  = 6 ; tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f  = zero
      tracerVarEdgeFull%scalar_mif_n%pos                    = 6 ; tracerVarEdgeFull%scalar_mif_n%f                    = 1._r8

!================================================
! primal cell, full level
!================================================

      allocate(tracerVarCellFull%scalar_delhp_avg_adv%f( nlev,mesh%nv))
      allocate(tracerVarCellFull%scalar_delhp_end_adv%f( nlev,mesh%nv))
      allocate(tracerVarCellFull%tend_mass_div_avg_adv%f(nlev,mesh%nv))
      allocate(tracerVarCellFull%scalar_tracer_mxrt_n%f(  ntracer,nlev,mesh%nv))
      allocate(tracerVarCellFull%scalar_tracer_mxrt_np1%f(ntracer,nlev,mesh%nv))
      allocate(tracerVarCellFull%scalar_tracer_mass_n%f(  ntracer,nlev,mesh%nv))
      allocate(tracerVarCellFull%scalar_tracer_mass_np1%f(ntracer,nlev,mesh%nv))
      allocate(tracerVarCellFull%tend_tracer_mass_hori%f( ntracer,nlev,mesh%nv))
      allocate(tracerVarCellFull%tend_tracer_mass_vert%f( ntracer,nlev,mesh%nv))
      allocate(scalar_tracer_name(ntracer))
      allocate(scalar_tracer_mxrt_min(ntracer))
      allocate(tracerVarCellFull%scalar_mif_n%f(nlev,mesh%nv))
      allocate(tracerVarCellFull%tend_tracer_mass_laplacian_2nd%f( ntracer,nlev,mesh%nv))
      allocate(tracerVarCellFull%tend_qv_n%f(nlev,mesh%nv))

      tracerVarCellFull%scalar_delhp_avg_adv%f     = zero  ; tracerVarCellFull%scalar_delhp_avg_adv%pos   = 0
      tracerVarCellFull%scalar_delhp_end_adv%f     = 1._r8 ; tracerVarCellFull%scalar_delhp_end_adv%pos   = 0
      tracerVarCellFull%scalar_tracer_mxrt_n%f     = zero  ; tracerVarCellFull%scalar_tracer_mxrt_n%pos   = 0
      tracerVarCellFull%scalar_tracer_mxrt_np1%f   = zero  ; tracerVarCellFull%scalar_tracer_mxrt_np1%pos = 0
      tracerVarCellFull%scalar_tracer_mass_n%f     = zero  ; tracerVarCellFull%scalar_tracer_mass_n%pos   = 0
      tracerVarCellFull%scalar_tracer_mass_np1%f   = zero  ; tracerVarCellFull%scalar_tracer_mass_np1%pos = 0
      tracerVarCellFull%tend_tracer_mass_hori%f    = zero  ; tracerVarCellFull%tend_tracer_mass_hori%pos  = 0
      tracerVarCellFull%tend_tracer_mass_vert%f    = zero  ; tracerVarCellFull%tend_tracer_mass_vert%pos  = 0 
      tracerVarCellFull%scalar_mif_n%f             = 1._r8 ; tracerVarCellFull%scalar_mif_n%pos           = 0
      tracerVarCellFull%tend_qv_n%f                = zero  ; tracerVarCellFull%tend_qv_n%pos              = 0
      tracerVarCellFull%tend_mass_div_avg_adv%pos          = 0
      tracerVarCellFull%tend_tracer_mass_laplacian_2nd%pos = 0

!================================================
! primal cell, face level
!================================================
      allocate(tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(nlev+1,mesh%nv))
      allocate(tracerVarCellFace%scalar_mif_n%f(nlev+1,mesh%nv))
! some initial
      tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f   = zero; tracerVarCellFace%scalar_eta_mass_flux_avg_adv%pos = 0
      tracerVarCellFace%scalar_mif_n%pos                 = 0
!
! set lower bounds for the coulped run
!
      if(.not.isolate_tracer_test)then
         scalar_tracer_mxrt_min(1)                           = 1e-12_r8 ! this is for CAMphys, WRF's lin is 1e-20
         if(ntracer.ge.2) scalar_tracer_mxrt_min(2:ntracer)  = zero
      else
         scalar_tracer_mxrt_min = zero
      end if

      allocate(tracerVarSurface%tend_hpressure%f    (mesh%nv))
      tracerVarSurface%tend_hpressure%f     = zero

      return
    end subroutine tracer_transport_vars_construct

    subroutine tracer_transport_vars_destruct
!================================================
! edge, full level
!================================================
      deallocate(tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f)
      deallocate(tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f)
      deallocate( tracerVarEdgeFull%scalar_mif_n%f)
!================================================
! primal cell, full level
!================================================
      deallocate(tracerVarCellFull%scalar_delhp_avg_adv%f)
      deallocate(tracerVarCellFull%scalar_delhp_end_adv%f)
      deallocate(tracerVarCellFull%tend_mass_div_avg_adv%f)
      deallocate(tracerVarCellFull%scalar_tracer_mxrt_n%f)
      deallocate(tracerVarCellFull%scalar_tracer_mxrt_np1%f)
      deallocate(tracerVarCellFull%scalar_tracer_mass_n%f)
      deallocate(tracerVarCellFull%scalar_tracer_mass_np1%f)
      deallocate(tracerVarCellFull%tend_tracer_mass_hori%f)
      deallocate(tracerVarCellFull%tend_tracer_mass_vert%f)
      deallocate(tracerVarCellFull%scalar_mif_n%f)
      deallocate(scalar_tracer_name)
      deallocate(scalar_tracer_mxrt_min)
      deallocate(tracerVarCellFull%tend_tracer_mass_laplacian_2nd%f)
      deallocate(tracerVarCellFull%tend_qv_n%f)
!================================================
! primal cell, face level
!================================================
      deallocate(tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f)
      deallocate(tracerVarCellFace%scalar_mif_n%f)
      deallocate(tracerVarSurface%tend_hpressure%f)

      return
     end subroutine tracer_transport_vars_destruct

     subroutine tracer_transport_evaluate_mif(mesh, nlev)
! io
     use omp_lib
     type(global_domain), intent(in) :: mesh
     integer(i4)        , intent(in) :: nlev
! local
     integer(i4)  :: iv, ie, ilev, icell1, icell2
     real(r8)     :: tmp

! simple test of dry core
       if(trim(working_mode).eq.'dycore')then
          tracerVarCellFull%scalar_tracer_mxrt_n%f = zero
          tracerVarCellFull%scalar_mif_n%f         = one
          tracerVarEdgeFull%scalar_mif_n%f         = one
          tracerVarCellFace%scalar_mif_n%f         = one
          return
       end if

! evaluate forcing at cell, which persists across dycore steps
       do iv = 1, mesh%nv_full
          do ilev = 1, nlev
             tracerVarCellFull%scalar_mif_n%f(ilev,iv) = one/(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,iv)))
          end do
       end do

! cell to edge
!$omp parallel  private(ie,icell1,icell2,ilev,tmp) 
!$omp do schedule(dynamic,100) 

       do ie = 1, mesh%ne
          icell1 = mesh%edt_v(1,ie)
          icell2 = mesh%edt_v(2,ie)
          do ilev = 1, nlev
#ifdef MIFAVG
             tracerVarEdgeFull%scalar_mif_n%f(ilev,ie) = half*(tracerVarCellFull%scalar_mif_n%f(ilev,icell1)+&
                                                               tracerVarCellFull%scalar_mif_n%f(ilev,icell2))
#else
             tmp  = half*(sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,icell1))+&
                          sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,icell2)))
             tracerVarEdgeFull%scalar_mif_n%f(ilev,ie) = one/(one+tmp)

#endif
          end do
       end do
!$omp end do nowait
!$omp end parallel 
! full to face
#ifdef MIFAVG
       do ilev = 2, nlev
          tracerVarCellFace%scalar_mif_n%f(ilev,:) = &
                                     0.5*(deta_full(ilev-1)/deta_face(ilev)*tracerVarCellFull%scalar_mif_n%f(ilev,:)+&
                                          deta_full(ilev)  /deta_face(ilev)*tracerVarCellFull%scalar_mif_n%f(ilev-1,:))
       end do
#else
       do iv = 1, mesh%nv_full
       do ilev = 2, nlev
               tmp = 0.5*(deta_full(ilev-1)/deta_face(ilev)*sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,iv))+&
                          deta_full(ilev)  /deta_face(ilev)*sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev-1,iv)))
               tracerVarCellFace%scalar_mif_n%f(ilev,iv) = one/(one+tmp)
       end do
       end do
#endif
! these bdy conditions actually donot matter, but we set here for consistency
       tracerVarCellFace%scalar_mif_n%f(1,:)      = tracerVarCellFull%scalar_mif_n%f(1,:)
       tracerVarCellFace%scalar_mif_n%f(nlev+1,:) = tracerVarCellFull%scalar_mif_n%f(nlev,:)

       return
     end subroutine tracer_transport_evaluate_mif

  end module grist_tracer_transport_vars_module
