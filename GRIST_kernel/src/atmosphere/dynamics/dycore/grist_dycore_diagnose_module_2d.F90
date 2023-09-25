
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Diagnostic module
! Revision history:
!----------------------------------------------------------------------------

 module grist_dycore_diagnose_module_2d


  use grist_constants,      only: gravity, i4, r8, rearth, one, zero, cp,cv,cvv
  use grist_domain_types,   only: global_domain
  use grist_data_types,     only: scalar_1d_field
  use grist_util_module,    only: write_string
  use grist_nml_module,     only: outdir, working_mode, nh_dynamics,&
                                  conserve_scheme, pv_order, grid_info, doNotDiagnose
  use grist_dycore_gcd_recon_module_2d, only: divergence_operator_2d        , &
                                              calc_vorticity_at_dual_cell_2d, &
                                              project_uv_2d                 , &
                                              vector_recon_perot_edge2cell_uv_2d
  use grist_nml_module,     only: nlev, nlevp
  use grist_dycore_vars_module, only: dycoreVarEdgeFull, dycoreVarVertFull, dycoreVarCellFace, dycoreVarCellFull, dycoreVarSurface
#ifndef SEQ_GRIST
  use grist_config_partition
#endif
  use grist_lib

  implicit none

   private

   public   :: dycore_diagnose_variables,         & 
               dycore_diagnose_total_energy,      &
               dycore_diagnose_total_energy_moist_e2,&
               dycore_diagnose_total_energy_moist_e3,&
               dycore_diagnose_mass,              &
               dycore_diagnose_mass_pt

   character(len=128)     :: c_glevel
   character(len=128)     :: c_adv
   character(len=128)     :: c_pv

  contains

!==========================================================================
!   Diagnose un-pronostic and un-required (by time integration) variables
!==========================================================================

  subroutine dycore_diagnose_variables(mesh)
! io
    type(global_domain),  intent(in)   :: mesh
! local
    integer(i4)                        :: ie
    integer(i4)                        :: iv
    real(r8)                           :: scalar
    real(r8)                           :: total_energy
    integer(i4)                        :: ilev


        call calc_vorticity_at_dual_cell_2d(mesh, dycoreVarEdgeFull%scalar_normal_velocity_n%f,&
                                                  dycoreVarVertFull%scalar_abs_vor_n%f,&
                                                  dycoreVarVertFull%scalar_rel_vor_n%f,&
                                                  nlev)

        call divergence_operator_2d(mesh, dycoreVarEdgeFull%scalar_normal_velocity_n%f, &
                                          dycoreVarCellFull%scalar_divergence_n%f, nlev)

        call vector_recon_perot_edge2cell_uv_2d(mesh, dycoreVarEdgeFull%scalar_normal_velocity_n%f, &
                                                      dycoreVarCellFull%scalar_U_wind_n%f, &
                                                      dycoreVarCellFull%scalar_V_wind_n%f, &
                                                      nlev)
#ifdef DYCORE_NHD
        dycoreVarCellFull%scalar_www_n%f(:,:) = 0.5_r8*(dycoreVarCellFace%scalar_www_n%f(1:nlev,:)+ dycoreVarCellFace%scalar_www_n%f(2:nlevp,:))
#endif

    return
  end subroutine dycore_diagnose_variables

!================================================
! Calculate moist "total-energy" at prime cell
! integral on a unit sphere
!================================================

  subroutine dycore_diagnose_total_energy_moist_e2(mesh)
  use grist_tracer_transport_vars_module, only: tracerVarCellFull
! io
   type(global_domain), intent(in)    :: mesh
! local
   integer                            :: iv,ie, ilev,ierr
   integer                            :: icell1, icell2
   real(r8)                           :: partx, part0, part1, part2, part3
   real(r8)                           :: scalar_energy_at_pc
   real(r8)                           :: scalar_energy_at_edge
   real(r8)                           :: scalar_potential_energy,scalar_potential_energy_all
   real(r8)                           :: scalar_internal_energy, scalar_internal_energy_all
   real(r8)                           :: scalar_kinetic_energy,scalar_kinetic_energy_all
   real(r8)                           :: scalar_total_energy,scalar_total_energy_all
   real(r8)                           :: delmp, gcv, delmp_c1, delmp_c2

!
! specified by locations
!
     scalar_energy_at_pc     = 0._r8
     scalar_energy_at_edge   = 0._r8
!
! specified by types
!
     scalar_potential_energy = 0._r8
     scalar_internal_energy  = 0._r8
     scalar_kinetic_energy   = 0._r8

     do iv = 1, mesh%nv
!
! compute phis*ps
!
         partx   = dycoreVarCellFace%scalar_geopotential_n%f(1,iv)*dycoreVarCellFace%scalar_mpressure_n%f(1,iv)
         part0   = zero
         part1   = zero
         part2   = zero

         do ilev = 1, nlev
! first evaluate delmp
            delmp = dycoreVarCellFace%scalar_mpressure_n%f(ilev+1,iv)-dycoreVarCellFace%scalar_mpressure_n%f(ilev,iv)
! evaluate generalized cv
            gcv   = (cv+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv)*cvv)/(one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
! phi part
            part1 = part1+delmp*dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv)
! temperature part
            part2 = part2+delmp*gcv*dycoreVarCellFull%scalar_temp_n%f(ilev,iv)

         end do

#ifdef DYCORE_NHD
         if(nh_dynamics)then
            do ilev = 1, nlevp
               if(ilev.eq.1.or.ilev.eq.nlevp)then
                  delmp = 0.5_r8*dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)
               else
                  delmp = dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)-dycoreVarCellFull%scalar_mpressure_n%f(ilev-1,iv)
               end if
               part0 = part0+0.5*delmp*(dycoreVarCellFace%scalar_www_n%f(ilev,iv)**2)
            end do
         end if
#endif
         scalar_energy_at_pc     = scalar_energy_at_pc    +(partx+part0+part1+part2)*mesh%plg_areag(iv) !*rearth*rearth
         scalar_kinetic_energy   = scalar_kinetic_energy  +part0*mesh%plg_areag(iv)
         scalar_potential_energy = scalar_potential_energy+partx*mesh%plg_areag(iv)         ! not a real PE
         scalar_internal_energy  = scalar_internal_energy +(part1+part2)*mesh%plg_areag(iv) ! not a real IE
 
     end do
!
! KE at edges
!
     do ie = 1, mesh%ne

        icell1   = mesh%edt_v(1,ie)
        icell2   = mesh%edt_v(2,ie)
        part3    = zero
        do ilev = 1, nlev
           delmp_c1 = dycoreVarCellFace%scalar_mpressure_n%f(ilev+1,icell1)-dycoreVarCellFace%scalar_mpressure_n%f(ilev,icell1)
           delmp_c2 = dycoreVarCellFace%scalar_mpressure_n%f(ilev+1,icell2)-dycoreVarCellFace%scalar_mpressure_n%f(ilev,icell2)
           part3 = part3+0.5_r8*(dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)**2)*&
                         0.5_r8*(delmp_c1+delmp_c2)
        end do

        if(icell1>mesh%nv_compute.or.icell2>mesh%nv_compute) then
           scalar_energy_at_edge = scalar_energy_at_edge+part3*mesh%edt_leng(ie)*mesh%edp_leng(ie)*0.5_r8  !*rearth*rearth ! rhombi area doubled to!*rearth*rearth
        else                                                                                               ! account two velocity component;  
           scalar_energy_at_edge = scalar_energy_at_edge+part3*mesh%edt_leng(ie)*mesh%edp_leng(ie)         ! see R10; resulting the same
        endif                                                                                              ! meaning as we sum over cell-based
     end do                                                                                                ! KE

     scalar_kinetic_energy = scalar_kinetic_energy+scalar_energy_at_edge
     scalar_total_energy   = scalar_energy_at_pc  +scalar_energy_at_edge

#ifndef SEQ_GRIST
     If(.not.doNotDiagnose)then

     call reduce(scalar_total_energy, scalar_total_energy_all, 'sum')
     call reduce(scalar_internal_energy, scalar_internal_energy_all, 'sum')
     call reduce(scalar_kinetic_energy, scalar_kinetic_energy_all, 'sum')
     call reduce(scalar_potential_energy, scalar_potential_energy_all, 'sum')

     End if
#else
     scalar_total_energy_all     = scalar_total_energy
     scalar_internal_energy_all  = scalar_internal_energy
     scalar_kinetic_energy_all   = scalar_kinetic_energy
     scalar_potential_energy_all = scalar_potential_energy
#endif

     call write_string(mesh%glevel     ,c_glevel)

     If(.not.doNotDiagnose)then

     if(mpi_rank() == 0)then
       open(1,file=trim(outdir)//"GCM.GRIST."//trim(grid_info)//"-"//trim(working_mode)//"-"//trim(conserve_scheme)//"-TE2.txt",access='append')
       write(1,*) "------------------------"
       write(1,*) "TE=",scalar_total_energy_all   ,"PE=",scalar_potential_energy_all
       write(1,*) "IE=",scalar_internal_energy_all,"KE=",scalar_kinetic_energy_all
       write(1,*) "------------------------"
       close(1)
       if(isnan(scalar_total_energy_all))then
         write(1,*) "nan comes out, model aborts"
         call mpi_abort()
       end if
     end if

     End if
   return
  end subroutine dycore_diagnose_total_energy_moist_e2

  subroutine dycore_diagnose_total_energy_moist_e3(mesh)
  use grist_tracer_transport_vars_module, only: tracerVarCellFull
! io
   type(global_domain), intent(in)    :: mesh
! local
   integer                            :: iv,ie, ilev,ierr
   integer                            :: icell1, icell2
   real(r8)                           :: partx, part0, part1, part2, part3
   real(r8)                           :: scalar_energy_at_pc
   real(r8)                           :: scalar_energy_at_edge
   real(r8)                           :: scalar_potential_energy,scalar_potential_energy_all
   real(r8)                           :: scalar_internal_energy, scalar_internal_energy_all
   real(r8)                           :: scalar_kinetic_energy,scalar_kinetic_energy_all
   real(r8)                           :: scalar_total_energy,scalar_total_energy_all
   real(r8)                           :: delmp, gcv, delmp_c1, delmp_c2

!
! specified by locations
!
     scalar_energy_at_pc     = 0._r8
     scalar_energy_at_edge   = 0._r8
!
! specified by types
!
     scalar_potential_energy = 0._r8
     scalar_internal_energy  = 0._r8
     scalar_kinetic_energy   = 0._r8

     do iv = 1, mesh%nv
!
! compute phis*ps
!
         partx   = dycoreVarSurface%scalar_geopotential_n%f(iv)*dycoreVarCellFace%scalar_mpressure_n%f(nlevp,iv)
         part0   = zero
         part1   = zero
         part2   = zero

         do ilev = 1, nlev
! first evaluate delmp
            delmp = dycoreVarCellFace%scalar_mpressure_n%f(ilev+1,iv)-dycoreVarCellFace%scalar_mpressure_n%f(ilev,iv)
! evaluate generalized cv
            gcv   = (cv+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv)*cvv)/(one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
! phi part
            part1 = part1+dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)*&
                    (dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)-dycoreVarCellFace%scalar_geopotential_n%f(ilev+1,iv))
! temperature part
            part2 = part2+delmp*gcv*dycoreVarCellFull%scalar_temp_n%f(ilev,iv)

         end do

#ifdef DYCORE_NHD
         if(nh_dynamics)then
            do ilev = 1, nlevp
               if(ilev.eq.1.or.ilev.eq.nlevp)then
                  delmp = 0.5_r8*dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)
               else
                  delmp = dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)-dycoreVarCellFull%scalar_mpressure_n%f(ilev-1,iv)
               end if
               part0 = part0+0.5*delmp*(dycoreVarCellFace%scalar_www_n%f(ilev,iv)**2)
            end do
         end if
#endif
         scalar_energy_at_pc     = scalar_energy_at_pc    +(partx+part0+part1+part2)*mesh%plg_areag(iv) !*rearth*rearth
         scalar_kinetic_energy   = scalar_kinetic_energy  +part0*mesh%plg_areag(iv)
         scalar_potential_energy = scalar_potential_energy+partx*mesh%plg_areag(iv)         ! not a real PE
         scalar_internal_energy  = scalar_internal_energy +(part1+part2)*mesh%plg_areag(iv) ! not a real IE
 
     end do
!
! KE at edges
!
     do ie = 1, mesh%ne

        icell1   = mesh%edt_v(1,ie)
        icell2   = mesh%edt_v(2,ie)
        part3    = zero
        do ilev = 1, nlev
           delmp_c1 = dycoreVarCellFace%scalar_mpressure_n%f(ilev+1,icell1)-dycoreVarCellFace%scalar_mpressure_n%f(ilev,icell1)
           delmp_c2 = dycoreVarCellFace%scalar_mpressure_n%f(ilev+1,icell2)-dycoreVarCellFace%scalar_mpressure_n%f(ilev,icell2)
           part3 = part3+0.5_r8*(dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)**2)*&
                         0.5_r8*(delmp_c1+delmp_c2)
        end do
        if(icell1>mesh%nv_compute.or.icell2>mesh%nv_compute) then                                         !*rearth*rearth ! rhombi area doubled to
           scalar_energy_at_edge = scalar_energy_at_edge+part3*mesh%edt_leng(ie)*mesh%edp_leng(ie)*0.5_r8 ! account two velocity component;
        else                                                                                              ! see R10; resulting the same
           scalar_energy_at_edge = scalar_energy_at_edge+part3*mesh%edt_leng(ie)*mesh%edp_leng(ie)        ! meaning as we sum over cell-based
        endif                                                                                             ! KE
     end do                                                                                               

     scalar_kinetic_energy = scalar_kinetic_energy+scalar_energy_at_edge
     scalar_total_energy   = scalar_energy_at_pc  +scalar_energy_at_edge
#ifndef SEQ_GRIST
     If(.not.doNotDiagnose)then

     call reduce(scalar_total_energy, scalar_total_energy_all, 'sum')
     call reduce(scalar_internal_energy, scalar_internal_energy_all, 'sum')
     call reduce(scalar_kinetic_energy, scalar_kinetic_energy_all, 'sum')
     call reduce(scalar_potential_energy, scalar_potential_energy_all, 'sum')
 
     End if
#else
     scalar_total_energy_all    = scalar_total_energy
     scalar_internal_energy_all = scalar_internal_energy
     scalar_kinetic_energy_all  = scalar_kinetic_energy
     scalar_potential_energy_all= scalar_potential_energy
#endif

     call write_string(mesh%glevel     ,c_glevel)

     If(.not.doNotDiagnose)then

     if(mpi_rank() == 0)then
       open(1,file=trim(outdir)//"GCM.GRIST."//trim(grid_info)//"-"//trim(working_mode)//"-"//trim(conserve_scheme)//"-TE3.txt",access='append')
            write(1,*) "------------------------"
            write(1,*) "TE=",scalar_total_energy_all   ,"PE=",scalar_potential_energy_all
            write(1,*) "IE=",scalar_internal_energy_all,"KE=",scalar_kinetic_energy_all
            write(1,*) "------------------------"
       close(1)
       if(isnan(scalar_total_energy_all))then
         write(1,*) "nan comes out, model aborts"
         call mpi_abort() 
       end if
     end if

     End if
   return
  end subroutine dycore_diagnose_total_energy_moist_e3

!================================================
! Calculate "total-energy" at prime cell
! used for dry core, moist  model nows
! has a  new one, but many moist tests are
! still based on this one for regression and check
!================================================

  subroutine dycore_diagnose_total_energy(mesh)
! io
   type(global_domain), intent(in)    :: mesh
! local
   integer                            :: iv,ie, ilev,ierr
   integer                            :: icell1, icell2
   real(r8)                           :: part0, part1, part2, part3
   real(r8)                           :: scalar_energy_at_pc
   real(r8)                           :: scalar_energy_at_edge
   real(r8)                           :: scalar_potential_energy,scalar_potential_energy_all
   real(r8)                           :: scalar_internal_energy, scalar_internal_energy_all
   real(r8)                           :: scalar_kinetic_energy,scalar_kinetic_energy_all
   real(r8)                           :: scalar_total_energy,scalar_total_energy_all

!
! specified by locations
     scalar_energy_at_pc     = 0._r8
     scalar_energy_at_edge   = 0._r8
! specified by types
     scalar_potential_energy = 0._r8
     scalar_internal_energy  = 0._r8
     scalar_kinetic_energy   = 0._r8

!
! Reminder: this routine does not account moisture and only for dry core,
!     |cpT+K+phis*hps is conserved for the current bdy condition;
!     |cpT+K+phi is not, but still used for development, so we can 
!      supervise "energy" for every term (phis=0 in most tests)
!

     do iv = 1, mesh%nv
!
! compute phis*ps
!
         part0  = 0._r8
         part1  = 0._r8
!         part1  = dycoreVarSurface%scalar_geopotential_n%f(iv)*dycoreVarSurface%scalar_hpressure_n%f(iv) ! "bit reproduce" uses this
         part2  = 0._r8
         do ilev = 1, nlev
! below for dev
            part1 = part1+dycoreVarCellFull%scalar_delhp_n%f(ilev,iv)*dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv)
            part2 = part2+dycoreVarCellFull%scalar_delhp_n%f(ilev,iv)*cp*dycoreVarCellFull%scalar_temp_n%f(ilev,iv)
         end do

#ifdef DYCORE_NHD
         if(nh_dynamics)then
            do ilev = 1, nlevp
               part0 = part0+0.5*dycoreVarCellFace%scalar_delhp_n%f(ilev,iv)*(dycoreVarCellFace%scalar_www_n%f(ilev,iv)**2)
            end do
         end if
#endif
         scalar_energy_at_pc     = scalar_energy_at_pc+(part0+part1+part2)*mesh%plg_areag(iv) !*rearth*rearth
         scalar_potential_energy = scalar_potential_energy+part1*mesh%plg_areag(iv)
         scalar_internal_energy  = scalar_internal_energy +part2*mesh%plg_areag(iv) ! actually enthalpy
         scalar_kinetic_energy   = scalar_kinetic_energy  +part0*mesh%plg_areag(iv)
 
     end do
 
!
! KE at edges
!
     do ie = 1, mesh%ne
        !icell1 = mesh%edt(ie)%v(1)
        !icell2 = mesh%edt(ie)%v(2)
        icell1 = mesh%edt_v(1,ie)
        icell2 = mesh%edt_v(2,ie)
        part3  = 0._r8
        do ilev = 1, nlev
           part3 = part3+0.5_r8*(dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)**2)*&
                         0.5_r8*(dycoreVarCellFull%scalar_delhp_n%f(ilev,icell1)+dycoreVarCellFull%scalar_delhp_n%f(ilev,icell2))
        end do
        if(icell1>mesh%nv_compute.or.icell2>mesh%nv_compute) then ! this is a bdy edge, remove double count
           scalar_energy_at_edge = scalar_energy_at_edge+part3*mesh%edt_leng(ie)*mesh%edp_leng(ie)*0.5_r8 !*rearth*rearth
        else                                                                                      
           scalar_energy_at_edge = scalar_energy_at_edge+part3*mesh%edt_leng(ie)*mesh%edp_leng(ie) !*rearth*rearth
        endif
     end do
     scalar_kinetic_energy = scalar_kinetic_energy+scalar_energy_at_edge
     scalar_total_energy   = scalar_energy_at_pc+scalar_energy_at_edge

#ifndef SEQ_GRIST
     If(.not.doNotDiagnose)then

     call reduce(scalar_total_energy, scalar_total_energy_all, 'sum')
     call reduce(scalar_internal_energy, scalar_internal_energy_all, 'sum')
     call reduce(scalar_kinetic_energy, scalar_kinetic_energy_all, 'sum')
     call reduce(scalar_potential_energy, scalar_potential_energy_all, 'sum')

     End if
#else
     scalar_total_energy_all    = scalar_total_energy
     scalar_internal_energy_all = scalar_internal_energy
     scalar_kinetic_energy_all  = scalar_kinetic_energy
     scalar_potential_energy_all= scalar_potential_energy
#endif
     call write_string(mesh%glevel     ,c_glevel)

     If(.not.doNotDiagnose)then
     if(mpi_rank() == 0)then
       open(1,file=trim(outdir)//"GCM.GRIST."//trim(grid_info)//"-"//trim(working_mode)//"-"//trim(conserve_scheme)//"-TE.txt",access='append')
       write(1,*) "------------------------"
       write(1,*) "TE=",scalar_total_energy_all   ,"PE=",scalar_potential_energy_all
       write(1,*) "IE=",scalar_internal_energy_all,"KE=",scalar_kinetic_energy_all
       write(1,*) "------------------------"
       close(1)
       if(isnan(scalar_total_energy_all))then
         write(1,*) "nan comes out, model aborts"
         call mpi_abort() 
       end if
     end if
     End if

   return
  end subroutine dycore_diagnose_total_energy

!================================================
!           Diagnose MASS
!================================================

  subroutine dycore_diagnose_mass(mesh)
! io
    type(global_domain),  intent(in) :: mesh
! local
    integer(i4)                          :: iv, ilev,ierr
    real(r8)                             :: scalar_mass_at_pc,scalar_mass_at_pc_all

     scalar_mass_at_pc  = 0._r8

     do iv = 1, mesh%nv
           scalar_mass_at_pc = scalar_mass_at_pc+dycoreVarSurface%scalar_hpressure_n%f(iv)*mesh%plg_areag(iv)
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
     open(1,file=trim(outdir)//"GCM.GRIST."//trim(grid_info)//"-"//trim(working_mode)//"-"//trim(conserve_scheme)//"-MASS.txt",access='append')
     write(1,*) scalar_mass_at_pc_all
     close(1)
     end if

     End if
     return
   end subroutine dycore_diagnose_mass

!================================================
!           Diagnose MASS*PT
!================================================

  subroutine dycore_diagnose_mass_pt(mesh)
! io
    type(global_domain),  intent(in) :: mesh
! local
    integer(i4)                          :: iv, ilev,ierr
    real(r8)                             :: scalar_mass_pt_at_pc,scalar_mass_pt_at_pc_all

     scalar_mass_pt_at_pc  = 0._r8

     do iv = 1, mesh%nv
         do ilev = 1, nlev
            scalar_mass_pt_at_pc = scalar_mass_pt_at_pc+&
                        dycoreVarCellFull%scalar_mass_pt_n%f(ilev,iv)*mesh%plg_areag(iv)
         end do
     end do

#ifndef SEQ_GRIST
     If(.not.doNotDiagnose)then

     call reduce(scalar_mass_pt_at_pc, scalar_mass_pt_at_pc_all, 'sum')

     End if
#else
     scalar_mass_pt_at_pc_all = scalar_mass_pt_at_pc
#endif
     call write_string(mesh%glevel     , c_glevel)

     If(.not.doNotDiagnose)then

     if(mpi_rank() == 0)then
     open(1,file=trim(outdir)//"GCM.GRIST."//trim(grid_info)//"-"//trim(working_mode)//"-"//trim(conserve_scheme)//"-PTMASS.txt",access='append')
     write(1,*) scalar_mass_pt_at_pc_all
     close(1)
     end if

     End if

     return
   end subroutine dycore_diagnose_mass_pt

  end module grist_dycore_diagnose_module_2d
