
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description:  Diagnose module, similar to that in gcrm but only for tracer 
!              variables
! Revision history:
!----------------------------------------------------------------------------

 module grist_tracer_transport_diagnose_module
#ifndef SEQ_GRIST
  use grist_lib
#endif
  use grist_constants,      only: gravity, i4, r8, rearth, cp
  use grist_domain_types,   only: global_domain
  use grist_util_module,    only: write_string
  use grist_nml_module,     only: outdir, working_mode, &
                                  conserve_scheme,&
                                  pv_order, nlev, ntracer, doNotDiagnose, grid_info
! tracer
  use grist_tracer_transport_vars_module, only: tracerVarCellFace, tracerVarCellFull, tracerVarEdgeFull, tracerVarSurface
  use grist_dycore_vars_module, only: dycoreVarCellFull

  implicit none

   private

   public   :: tracer_transport_diagnose_global_tracer_mass, &
               tracer_transport_diagnose_global_crnum,       &
               tracer_transport_diagnose_global_mass

   character(len=128)     :: c_glevel
   character(len=128)     :: c_pv
   character(len=128)     :: c_tracer

  contains

!================================================
!           Diagnose Total Tracer MASS
!================================================

  subroutine tracer_transport_diagnose_global_tracer_mass(mesh)
! io
    type(global_domain),  intent(inout) :: mesh
! local
    integer(i4)                          :: iv, ilev, itracer, ierr
    real(r8), allocatable                :: scalar_tracer_mass_at_pc(:)
    real(r8), allocatable                :: scalar_tracer_mass_at_pc_all(:)

     if(.not.allocated(scalar_tracer_mass_at_pc))     allocate(scalar_tracer_mass_at_pc(ntracer))
     if(.not.allocated(scalar_tracer_mass_at_pc_all)) allocate(scalar_tracer_mass_at_pc_all(ntracer))
     scalar_tracer_mass_at_pc(:) = 0._r8

     do iv = 1, mesh%nv
         do ilev = 1, nlev
            do itracer = 1, ntracer
               scalar_tracer_mass_at_pc(itracer) = scalar_tracer_mass_at_pc(itracer)+&
                        tracerVarCellFull%scalar_tracer_mass_n%f(itracer,ilev,iv)*mesh%plg_areag(iv)
            end do
         end do
     end do
#ifndef SEQ_GRIST
     If(.not.doNotDiagnose)then

     call reduce(scalar_tracer_mass_at_pc, scalar_tracer_mass_at_pc_all, ntracer, 'sum')

     End if
#else
     scalar_tracer_mass_at_pc_all = scalar_tracer_mass_at_pc
     scalar_tracer_mass_at_pc_all = scalar_tracer_mass_at_pc
#endif

     call write_string(mesh%glevel     , c_glevel)
     call write_string(pv_order(3)     , c_pv)

     If(.not.doNotDiagnose)then

     do itracer = 1, ntracer
        call write_string(itracer, c_tracer)
#ifndef SEQ_GRIST
        if(mpi_rank() == 0) then
#endif
           open(1,file=trim(outdir)//&
                        "GCM.TRACER"//trim(c_tracer)//".GRIST-"//&
                        "GLEVEL"//trim(c_glevel)//&
                        "-"//trim(working_mode)//&
                        "-PV"//trim(c_pv)//&
                        "-"//trim(conserve_scheme)//&
                        "-TRACER-MASS.txt",access='append')
           write(1,*) scalar_tracer_mass_at_pc_all(itracer)
           close(1)
#ifndef SEQ_GRIST
        end if
#endif
     end do

     End if

     return
   end subroutine tracer_transport_diagnose_global_tracer_mass

   subroutine tracer_transport_diagnose_global_mass(mesh)
     implicit none
! io
     type(global_domain), intent(in)    ::   mesh
! local
     integer(i4)                          :: iv, ilev,ierr
     real(r8)                             :: scalar_mass_at_pc,scalar_mass_at_pc_all

     scalar_mass_at_pc  = 0._r8

     do iv = 1, mesh%nv
           scalar_mass_at_pc = scalar_mass_at_pc+tracerVarSurface%tend_hpressure%f(iv)*mesh%plg_areag(iv)
     end do
#ifndef SEQ_GRIST
     If(.not.doNotDiagnose)then
     call reduce(scalar_mass_at_pc, scalar_mass_at_pc_all, 'sum')
     End if
#else
     scalar_mass_at_pc_all = scalar_mass_at_pc
#endif
 
     call write_string(mesh%glevel     , c_glevel)

     If(.not.doNotDiagnose)then

     if(mpi_rank() == 0)then
     open(1,file=trim(outdir)//"GCM.GRIST."//trim(grid_info)//"-"//trim(working_mode)//"-"//trim(conserve_scheme)//"-TRACER-MASS-TEND.txt",access='append')
     write(1,*) scalar_mass_at_pc_all
     close(1)
     end if

     End if
     return
   end subroutine tracer_transport_diagnose_global_mass

  subroutine tracer_transport_diagnose_global_crnum(mesh,dtime)
! io
    type(global_domain),  intent(inout) :: mesh
    real(r8)           ,  intent(in)    :: dtime
! local
    real(r8)                :: scalar_tracer_hori_crnum
    real(r8)                :: scalar_tracer_vert_crnum
    real(r8)                :: scalar_tracer_hori_crnum_all
    real(r8)                :: scalar_tracer_vert_crnum_all
    integer(i4)             :: ie, iv, ilev, itracer, ierr

     scalar_tracer_hori_crnum     = 0._r8
     scalar_tracer_vert_crnum     = 0._r8
     scalar_tracer_hori_crnum_all = 0._r8
     scalar_tracer_vert_crnum_all = 0._r8

! hori crnum
     do ie = 1, mesh%ne
         do ilev = 1, nlev
            scalar_tracer_hori_crnum = max(scalar_tracer_hori_crnum,&
                                           abs(tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f(ilev,ie)*&
                                               dtime/(rearth*mesh%edt_leng(ie))))
         end do
     end do
! vert crnum
     do iv = 1, mesh%nv
         do ilev = 2, nlev
            scalar_tracer_vert_crnum = max(scalar_tracer_vert_crnum,&
                                           abs(tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(ilev,iv)*&
                                               dtime/(0.5_r8*(dycoreVarCellFull%scalar_delhp_n%f(ilev-1,iv)+&
                                                              dycoreVarCellFull%scalar_delhp_n%f(ilev,iv)))))
         end do
     end do
#ifndef SEQ_GRIST
     If(.not.doNotDiagnose)then

     call reduce(scalar_tracer_hori_crnum, scalar_tracer_hori_crnum_all, 'max')
     call reduce(scalar_tracer_vert_crnum, scalar_tracer_vert_crnum_all, 'max')

     End if
#else
     scalar_tracer_hori_crnum_all = scalar_tracer_hori_crnum
     scalar_tracer_vert_crnum_all = scalar_tracer_vert_crnum
#endif

     call write_string(mesh%glevel     , c_glevel)
     call write_string(pv_order(3)     , c_pv)

     If(.not.doNotDiagnose)then
#ifndef SEQ_GRIST
     if(mpi_rank() == 0) then
#endif
        open(1,file=trim(outdir)//&
                    "GCM.GRIST-"//&
                    "GLEVEL"//trim(c_glevel)//&
                    "-"//trim(working_mode)//&
                    "-PV"//trim(c_pv)//&
                    "-"//trim(conserve_scheme)//&
                    "-MAX.CRNUM.txt",access='append')
        write(1,*) scalar_tracer_hori_crnum_all, scalar_tracer_vert_crnum_all 
        close(1)
#ifndef SEQ_GRIST
     end if
#endif
   
     End if

     return
   end subroutine tracer_transport_diagnose_global_crnum

  end module grist_tracer_transport_diagnose_module
