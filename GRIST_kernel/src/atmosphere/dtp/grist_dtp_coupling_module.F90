
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description:    A GENERAL DYCORE-TRACER-PHYSICS COUPLING STRATEGY (Z20, MWR)
!
!              This module handles DTP coupling stuff associated with DATA 
!          Structure, including data management/update and MPI. The general
!          purpose of PDC is to provide dynamics with appropriate physics
!          forcing without the need to tell dynamics what kind of physics
!          is used, but be controlled by explicit configuration in nml and
!          the selection of different physics packages.
!
!      1.  GRIST-A has three kinds of possible physics tendencies, as follows:
!          ptend_F1 is updated inside dycore and tracer after the dyn update;
!          ptend_F2 is updated here (after calling of physics) with one model/fast-physics step.
!          ptend_RK is updated inside Runge-Kutta loop of dycore
!          The influence of these three tendencies are controlled by ptend_${type}_on as a configurer.
!
!      2.  The seperation of F1, F2, and RK depends on model physics. Note that for a layer-integrated
!          scalar (e.g., pt), physics forcing is evaluated for non-mass-weighted quantity,
!          while updated for mass-weighted quantity, which is obtained by multiplying delhp with 
!          the non-mass-weighted tendency.
!
!      3.  For most cases, www/phi posses a zero tendency, but we still code it here for a general purpose 
!
! Present Package:  1. HSDRY_EXP: Held-Suarez dry forcing (only run with dycore)
!                   2. DCMIP2016-BW terminator
!                   3. DCMIP2016-SC Kessler microphysics following Klemp et al 2015 and DCMIP variants
!                   4. DCMIP2016-TC Simple Physics of RJ2012, JAMES
!                   Above are simple physics, below are full physics
!                   6. AMIPC_PHYSICS Physics, PhysPkg based on CAM5
!                   7. AMIPW_PHYSICS, from WRF
!
! Revision history:
!
!       1. Still exploring whether F1 is necessary inside dycore because RK
!          already does most jobs; F2 is used most often in a time-split manner.
!          F2 can be used with fast-physics when in a slow-fast seperate mode
!
!       2. Now the tendency evaluation and state coupling have been seperated to allow 
!          larger flexiblity
!
!       3. LAM_DOMAIN added
!
!----------------------------------------------------------------------------

 module grist_dtp_coupling_module

   use grist_constants,                  only: r8, i4, pi, zero, one
   use grist_nml_module,                 only: ptend_f2_on, TC_pbl_flag, nh_dynamics
   use grist_domain_types,               only: global_domain
   use grist_data_types,                 only: exchange_field_list_2d, exchange_field_list_3d
! data module, used vars should be explicitly stated here to avoid any potential conflication
   use grist_dycore_vars_module,         only: dycoreVarCellFull,  &
                                               dycoreVarCellFace,  &
                                               dycoreVarEdgeFull

   use grist_tracer_transport_vars_module, only: tracerVarCellFull, &
                                                 tracerVarEdgeFull

   use grist_tracer_transport_utils_module, only: tracer_transport_check_mxrt, &
                                                  tracer_transport_fixer_mxrt, &
                                                  tracer_transport_qneg3_mxrt

   use grist_physics_data_structure,     only: ptend_f1, ptend_f2, ptend_rk, &
                                               pstate
   use grist_dtp_dcmip2016_terminator,   only: tendency_terminator
#ifndef SEQ_GRIST
! mpi-comm
   use grist_config_partition,           only: exchange_data_2d_add, &
                                               exchange_data_3d_add, &
                                               exchange_data_2d,     &
                                               exchange_data_3d
#endif
   use grist_lib
! model physics interface; simple physics below
   use grist_dtp_interface_module,       only: grist_dtp_interface_dcmip2016_tc_a, &
                                               grist_dtp_interface_dcmip2016_sc  , &
                                               grist_dtp_interface_dcmip2016_sc_a, &
                                               grist_dtp_interface_dcmip2016_mitc, &
                                               grist_dtp_interface_hsdry_exp
   use grist_physics_idealized_package,  only: phys_tend_hsdry_exp
#ifdef AMIPC_PHYSICS
   use grist_dtp_interface_full_physpkg1,only: grist_full_physpkg1_PhysC_run
#endif
#ifdef AMIPW_PHYSICS
   use grist_dtp_interface_full_physpkg2,only: grist_full_physpkg2_PhysW_run
#endif

   implicit none
   private
  
   public :: grist_dtp_coupling_driver_tend, &
             grist_dtp_coupling_driver_coup, &
             grist_d2t_coupling_driver     , &
             grist_dtp_time_average

   real(r8), parameter    :: half  = 0.5_r8

   contains

   subroutine grist_dtp_coupling_driver_tend(mesh, nlev, ntracer, nmif, dtime, istep, physpkg)

     type(global_domain)  , intent(inout) :: mesh
     integer(i4),           intent(in)    :: nlev
     integer(i4),           intent(in)    :: ntracer
     integer(i4),           intent(in)    :: nmif
     real(r8)   ,           intent(in)    :: dtime
     integer(i4),           intent(in)    :: istep
     character(len=*)    ,  intent(in)    :: physpkg
! local
     type(exchange_field_list_2d),pointer :: field_head_2d
     type(exchange_field_list_3d),pointer :: field_head_3d
     integer(i4)                          :: itracer, iv, ilev, ie
     real(r8)                             :: tmp1

     field_head_2d=>null()
     field_head_3d=>null()

     select case(trim(physpkg))

     case('HSDRY_EXP')
!
! if first step initialzied
! if other, the latest state from last step/for this step
!
         mesh%nv  = mesh%nv_halo(1)
         dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/& ! overwritten each RK step
                                                      dycoreVarCellFull%scalar_delhp_n%f
!
! generate physics forcing from phys interface, which are
! held as constant during the dycore integration
! if hdc, pressure is hpressure
!
! forcing
         call grist_dtp_interface_hsdry_exp(mesh, nlev, dtime , physpkg,  &
                                            dycoreVarCellFull%scalar_pressure_n%f,           &
                                            dycoreVarCellFace%scalar_pressure_n%f,           &
                                            dycoreVarEdgeFull%scalar_normal_velocity_n%f  ,&
                                            dycoreVarCellFull%scalar_mass_pt_n%f            ,&
                                            dycoreVarCellFull%scalar_potential_temp_n%f     ,&
                                            ptend_f1%tend_normal_velocity_at_edge_full_level%f,&
                                            ptend_f1%tend_potential_temp_at_pc_full_level%f)

         ptend_f1%tend_www_at_pc_face_level%f         = zero
         ptend_f1%tend_tracer_mxrt_at_pc_full_level%f = zero

#ifdef REGRESSION_HSDRY
         call phys_tend_hsdry_exp(mesh, nlev,dtime, &
                                             dycoreVarCellFull%scalar_pressure_n%f  , &
                                             dycoreVarCellFace%scalar_pressure_n%f  , &
                                             dycoreVarCellFull%scalar_potential_temp_n%f     , &
                                             dycoreVarEdgeFull%scalar_normal_velocity_n%f  , &
                                             dycoreVarCellFace%scalar_www_n%f                , &
                                             ptend_f1%tend_normal_velocity_at_edge_full_level%f, &
                                             ptend_f1%tend_www_at_pc_face_level%f              , &
                                             ptend_f1%tend_potential_temp_at_pc_full_level%f  )
#endif
         ptend_f1%tend_mass_pt_at_pc_full_level%f     = ptend_f1%tend_potential_temp_at_pc_full_level%f*&
                                                        dycoreVarCellFull%scalar_delhp_n%f
         do itracer = 1, ntracer
            ptend_f1%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f1%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                           dycoreVarCellFull%scalar_delhp_n%f
         end do
         ptend_rk = ptend_f1
! forcing2
         ptend_f2%tend_normal_velocity_at_edge_full_level%f  = zero
         ptend_f2%tend_potential_temp_at_pc_full_level%f     = zero
         ptend_f2%tend_tracer_mxrt_at_pc_full_level%f        = zero
         ptend_f2%tend_www_at_pc_face_level%f                = zero
! forcing1
         ptend_f1 = ptend_f2

         mesh%nv = mesh%nv_compute

      case('DCMIP2016-TC')

         call grist_dtp_interface_dcmip2016_tc_a(mesh, nlev, mesh%nv_halo(1), dtime                 , &
                                                 dycoreVarCellFull%scalar_mpressure_n%f              , &
                                                 dycoreVarCellFace%scalar_mpressure_n%f              , &
                                                 !dycoreVarCellFull%scalar_delhp_n%f                 , &
                                                 dycoreVarEdgeFull%scalar_normal_velocity_n%f      , &
                                                 dycoreVarCellFull%scalar_potential_temp_n%f         , & !thetam
                                                 dycoreVarCellFull%scalar_temp_n%f                   , &
                                                 tracerVarCellFull%scalar_tracer_mxrt_n%f            , &
                                                 0, TC_pbl_flag                                     , &
                                                 ptend_f2%tend_normal_velocity_at_edge_full_level%f , &
                                                 ptend_f2%tend_potential_temp_at_pc_full_level%f    , &
                                                 ptend_f2%tend_tracer_mxrt_at_pc_full_level%f       , &
                                                 pstate%scalar_precl_surface%f)
!
! obtain tendency for compound variables
!
          ptend_f2%tend_mass_pt_at_pc_full_level%f = ptend_f2%tend_potential_temp_at_pc_full_level%f*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          do itracer = 1, ntracer
             ptend_f2%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          end do
          ptend_f2%tend_www_at_pc_face_level%f               = zero
          ptend_rk = ptend_f2
          ptend_f1 = ptend_f2

      case('DCMIP2016-BW-RJ')

         call grist_dtp_interface_dcmip2016_tc_a(mesh, nlev, mesh%nv_halo(1), dtime                 , &
                                                 dycoreVarCellFull%scalar_mpressure_n%f              , &
                                                 dycoreVarCellFace%scalar_mpressure_n%f              , &
                                                 dycoreVarEdgeFull%scalar_normal_velocity_n%f      , &
                                                 dycoreVarCellFull%scalar_potential_temp_n%f         , & !thetam
                                                 dycoreVarCellFull%scalar_temp_n%f                   , &
                                                 tracerVarCellFull%scalar_tracer_mxrt_n%f            , &
                                                 1, TC_pbl_flag                                     , &
                                                 ptend_f2%tend_normal_velocity_at_edge_full_level%f , &
                                                 ptend_f2%tend_potential_temp_at_pc_full_level%f    , &
                                                 ptend_f2%tend_tracer_mxrt_at_pc_full_level%f       , &
                                                 pstate%scalar_precl_surface%f)
!
! obtain tendency for compound variables
!
          ptend_f2%tend_mass_pt_at_pc_full_level%f = ptend_f2%tend_potential_temp_at_pc_full_level%f*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          do itracer = 1, ntracer
             ptend_f2%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          end do
          ptend_f2%tend_www_at_pc_face_level%f               = zero
          ptend_rk = ptend_f2
          ptend_f1 = ptend_f2

      case('DCMIP2016-BW-TMT')
! teminator
          do iv = 1, mesh%nv_halo(1)
             if(mesh%vtx_lon(iv).ge.0)then
                tmp1 = mesh%vtx_lon(iv)/(pi/180._r8)
             else
                tmp1 = (mesh%vtx_lon(iv)+2._r8*pi)/(pi/180._r8)
             end if
             do ilev = 1, nlev
                call tendency_Terminator(real(mesh%vtx_lat(iv),r8)/(pi/180._r8), tmp1, &
                                         tracerVarCellFull%scalar_tracer_mxrt_n%f(2,ilev,iv), &
                                         tracerVarCellFull%scalar_tracer_mxrt_n%f(3,ilev,iv), &
                                         dtime, &
                                         ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(2,ilev,iv), &
                                         ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(3,ilev,iv))
             end do
          end do
          ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(1,:,:) = zero
          do itracer = 1, ntracer
             ptend_f2%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          end do

      case('DCMIP2016-SC','DCMIP2016-BW-SC')

        call grist_dtp_interface_dcmip2016_sc_a(mesh, nlev, mesh%nv_halo(1), dtime           , &
                                             dycoreVarCellFull%scalar_mpressure_n%f           , &
                                             dycoreVarCellFace%scalar_geopotential_n%f        , &
                                             dycoreVarCellFull%scalar_geopotential_n%f        , &
                                             dycoreVarCellFull%scalar_delhp_n%f               , &
                                             dycoreVarCellFull%scalar_potential_temp_n%f      , & !thetam
                                             tracerVarCellFull%scalar_tracer_mxrt_n%f         , &
                                             ptend_f2%tend_potential_temp_at_pc_full_level%f , &
                                             ptend_f2%tend_tracer_mxrt_at_pc_full_level%f    , &
                                             pstate%scalar_precl_surface%f)

        ptend_f2%tend_mass_pt_at_pc_full_level%f =  ptend_f2%tend_potential_temp_at_pc_full_level%f*&
                                                    dycoreVarCellFull%scalar_delhp_n%f
        do itracer = 1, ntracer
             ptend_f2%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                     tracerVarCellFull%scalar_delhp_end_adv%f
        end do
        ptend_f2%tend_normal_velocity_at_edge_full_level%f   = zero
        ptend_f2%tend_www_at_pc_face_level%f                 = zero
        ptend_rk = ptend_f2
        ptend_f1 = ptend_f2

     case('DCMIP2016-MITC')

         call grist_dtp_interface_dcmip2016_mitc(mesh, nlev, mesh%nv_halo(1), dtime                 , &
                                                 dycoreVarCellFull%scalar_mpressure_n%f              , &
                                                 dycoreVarCellFace%scalar_mpressure_n%f              , &
                                                 dycoreVarEdgeFull%scalar_normal_velocity_n%f      , &
                                                 dycoreVarCellFull%scalar_potential_temp_n%f         , & !thetam
                                                 dycoreVarCellFull%scalar_temp_n%f                   , &
                                                 tracerVarCellFull%scalar_tracer_mxrt_n%f            , &
                                                 ptend_f2%tend_normal_velocity_at_edge_full_level%f , &
                                                 ptend_f2%tend_potential_temp_at_pc_full_level%f    , &
                                                 ptend_f2%tend_tracer_mxrt_at_pc_full_level%f       , &
                                                 pstate%scalar_precl_surface%f)
!
! obtain tendency for compound variables
!
          ptend_f2%tend_mass_pt_at_pc_full_level%f = ptend_f2%tend_potential_temp_at_pc_full_level%f*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          do itracer = 1, ntracer
             ptend_f2%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          end do
          ptend_f2%tend_www_at_pc_face_level%f               = zero
          ptend_rk = ptend_f2
          ptend_f1 = ptend_f2

#ifdef AMIPC_PHYSICS
      case('AMIPC_PHYSICS')
#ifndef NO_PHYSHALO
        call grist_full_physpkg1_PhysC_run(mesh, ntracer, nmif, nlev, mesh%nv_halo(1), istep, dtime , &
#else
        call grist_full_physpkg1_PhysC_run(mesh, ntracer, nmif, nlev, mesh%nv_compute, istep, dtime , &
#endif
                                                   tracerVarCellFull%scalar_mif_n%f                    , &
                                                   dycoreVarCellFull%scalar_geopotential_n%f           , &
                                                   dycoreVarCellFace%scalar_geopotential_n%f           , &
                                                   dycoreVarCellFull%scalar_mpressure_n%f              , &
                                                   dycoreVarCellFace%scalar_mpressure_n%f              , &
                                                   dycoreVarEdgeFull%scalar_normal_velocity_n%f      , &
                                                   dycoreVarCellFull%scalar_potential_temp_n%f         , & ! ptm
                                                   dycoreVarCellFull%scalar_temp_n%f                   , &
                                                   tracerVarCellFull%scalar_tracer_mxrt_n%f            , &
                                                   dycoreVarCellFull%scalar_omega_timavg%f          , &
                                                   ptend_f2%tend_normal_velocity_at_edge_full_level%f , &
                                                   ptend_f2%tend_potential_temp_at_pc_full_level%f    , &
                                                   ptend_f2%tend_tracer_mxrt_at_pc_full_level%f      )
!
! obtain tendency for compound variables
!
          ptend_f2%tend_mass_pt_at_pc_full_level%f = ptend_f2%tend_potential_temp_at_pc_full_level%f*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          do itracer = 1, ntracer
             ptend_f2%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          end do
          ptend_f2%tend_www_at_pc_face_level%f     = zero
          ptend_rk = ptend_f2
          ptend_f1 = ptend_f2
#endif

#ifdef AMIPW_PHYSICS
      case('AMIPW_PHYSICS')

        call grist_full_physpkg2_PhysW_run(mesh, ntracer, nmif, nlev, mesh%nv_halo(1), istep, dtime , &
                                                tracerVarCellFull%scalar_mif_n%f                    , &
                                                dycoreVarCellFull%tend_pt_n%f                       , &
                                                tracerVarCellFull%tend_qv_n%f                       , &
                                                dycoreVarCellFull%scalar_geopotential_n%f           , &
                                                dycoreVarCellFace%scalar_geopotential_n%f           , &
                                                dycoreVarCellFull%scalar_mpressure_n%f              , &
                                                dycoreVarCellFace%scalar_mpressure_n%f              , &
                                                dycoreVarEdgeFull%scalar_normal_velocity_n%f      , &
                                                dycoreVarCellFull%scalar_potential_temp_n%f         , & ! ptm
                                                dycoreVarCellFull%scalar_temp_n%f                   , &
                                                tracerVarCellFull%scalar_tracer_mxrt_n%f            , &
                                                dycoreVarCellFace%scalar_www_timavg%f               , &
                                                dycoreVarCellFull%scalar_omega_timavg%f          , &
                                                dycoreVarCellFull%scalar_delhp_n%f                  )
!
! set tendency at bdy1&2 cell/edge to zero (they are data zone)
!
#ifdef LAM_DOMAIN
        do iv = 1, mesh%nv_full
           if(mesh%plg_mark(iv).eq.177)then
              ptend_f2%tend_potential_temp_at_pc_full_level%f(:,iv)      = zero
              ptend_rk%tend_potential_temp_at_pc_full_level%f(:,iv)      = zero
              ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(:,:,iv)       = zero
           end if
        end do

        do ie = 1, mesh%ne_full
           if(mesh%edt_mark(ie).eq.1771.or.mesh%edt_mark(ie).eq.1772)then
              ptend_f2%tend_normal_velocity_at_edge_full_level%f(:,ie)   = zero
           end if
        end do
#endif

!
! obtain tendency for compound variables
!
          ptend_f2%tend_mass_pt_at_pc_full_level%f = ptend_f2%tend_potential_temp_at_pc_full_level%f*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          ptend_rk%tend_mass_pt_at_pc_full_level%f = ptend_rk%tend_potential_temp_at_pc_full_level%f*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          do itracer = 1, ntracer
             ptend_f2%tend_tracer_mass_at_pc_full_level%f(itracer,:,:) = ptend_f2%tend_tracer_mxrt_at_pc_full_level%f(itracer,:,:)*&
                                                     dycoreVarCellFull%scalar_delhp_n%f
          end do
          ptend_f2%tend_www_at_pc_face_level%f     = zero

#endif
      case('none')
! send to dynamics
         print*, "PHYSPKG: NONE"
         ptend_f1%tend_normal_velocity_at_edge_full_level%f  = zero
         ptend_f1%tend_mass_pt_at_pc_full_level%f            = zero
         ptend_f1%tend_www_at_pc_face_level%f                = zero
         ptend_f1%tend_tracer_mass_at_pc_full_level%f        = zero

         ptend_rk = ptend_f1
         ptend_f2 = ptend_f1
      case default
         print*, "PHYSPKG: NULL"
! send to dynamics
         ptend_f1%tend_normal_velocity_at_edge_full_level%f  = zero
         ptend_f1%tend_mass_pt_at_pc_full_level%f            = zero
         ptend_f1%tend_www_at_pc_face_level%f                = zero
         ptend_f1%tend_tracer_mass_at_pc_full_level%f        = zero

         ptend_rk = ptend_f1
         ptend_f2 = ptend_f1

      end select

      return
   end subroutine grist_dtp_coupling_driver_tend

!-----------------------------------------------------------------------------
! F2 physics coupling depending on physics package
!-----------------------------------------------------------------------------

   subroutine grist_dtp_coupling_driver_coup(mesh, nlev, ntracer, dtime, physpkg)

     type(global_domain)  , intent(inout) :: mesh
     integer(i4),           intent(in)    :: nlev
     integer(i4),           intent(in)    :: ntracer
     real(r8)   ,           intent(in)    :: dtime
     character(len=*)    ,  intent(in)    :: physpkg
! local
     integer(i4)        :: iv, ilev, itracer
     type(exchange_field_list_2d),pointer :: field_head_2d
     type(exchange_field_list_3d),pointer :: field_head_3d

     field_head_2d=>null()
     field_head_3d=>null()

     select case(trim(physpkg))

      case('HSDRY_EXP')
          return

      case('DCMIP2016-TC','DCMIP2016-BW-RJ','DCMIP2016-MITC','AMIPC_PHYSICS','AMIPW_PHYSICS')
!
! F2 physics
         dycoreVarEdgeFull%scalar_normal_velocity_n%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f+&
                                                         ptend_f2%tend_normal_velocity_at_edge_full_level%f*dtime
         dycoreVarCellFull%scalar_mass_pt_n%f    = dycoreVarCellFull%scalar_mass_pt_n%f+&
                                                  dtime*ptend_f2%tend_mass_pt_at_pc_full_level%f

!
! exchange data
!
#ifndef SEQ_GRIST
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
         call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
!
! diagnose raw-var for next-step dycore & tracer evaluation
!
         dycoreVarCellFull%scalar_potential_temp_n%f    = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
!
! tracer
!
         tracerVarCellFull%scalar_tracer_mass_n%f       = tracerVarCellFull%scalar_tracer_mass_n%f+&
                                                         dtime*ptend_f2%tend_tracer_mass_at_pc_full_level%f
#ifndef SEQ_GRIST
         call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_n)
         call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
         do itracer = 1, ntracer
            tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)  = tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:)/&
                                                                    dycoreVarCellFull%scalar_delhp_n%f
         end do
!
        call tracer_transport_qneg3_mxrt(mesh%nv_full, nlev, ntracer,tracerVarCellFull%scalar_tracer_mxrt_n%f,"after physics coupling in driver_coup")
        call tracer_transport_check_mxrt(mesh,tracerVarCellFull%scalar_tracer_mxrt_n,"after physics coupling in driver_coup")
        do itracer = 1, ntracer
           tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)*tracerVarCellFull%scalar_delhp_end_adv%f
        end do

      case('DCMIP2016-BW-TMT')

! tracer mass
        tracerVarCellFull%scalar_tracer_mass_n%f    = tracerVarCellFull%scalar_tracer_mass_n%f+dtime*ptend_f2%tend_tracer_mass_at_pc_full_level%f
#ifndef SEQ_GRIST
        call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_n)
        call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
        do itracer = 2, ntracer
           tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f
        end do

      case('DCMIP2016-SC','DCMIP2016-BW-SC')
!
! update F2 physics forcing depending on physics
!
! [1] mass_pt
        dycoreVarCellFull%scalar_mass_pt_n%f   = dycoreVarCellFull%scalar_mass_pt_n%f+dtime*ptend_f2%tend_mass_pt_at_pc_full_level%f
#ifndef SEQ_GRIST
        call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
        dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
! [2] tracer mass
        tracerVarCellFull%scalar_tracer_mass_n%f    = tracerVarCellFull%scalar_tracer_mass_n%f+dtime*ptend_f2%tend_tracer_mass_at_pc_full_level%f
#ifndef SEQ_GRIST
        call exchange_data_3d_add(mesh,field_head_3d,tracerVarCellFull%scalar_tracer_mass_n)
        call exchange_data_3d(mesh%local_block,field_head_3d)
#endif
        do itracer = 1, ntracer
           tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:)/tracerVarCellFull%scalar_delhp_end_adv%f
        end do
!
! do this for all halo regions so as to permit overlapping
!
        call tracer_transport_qneg3_mxrt(mesh%nv_full, nlev, ntracer,tracerVarCellFull%scalar_tracer_mxrt_n%f,"after physics coupling in driver_coup")
        call tracer_transport_check_mxrt(mesh,tracerVarCellFull%scalar_tracer_mxrt_n,"after physics coupling in driver_coup")
        do itracer = 1, ntracer
           tracerVarCellFull%scalar_tracer_mass_n%f(itracer,:,:) = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,:,:)*tracerVarCellFull%scalar_delhp_end_adv%f
        end do

      case('none')
          return
      case default
          print*, "no physpkg is activated, return"
          return
      end select

      return
   end subroutine grist_dtp_coupling_driver_coup

!=======================================================================
! This routine is called after each dycore mode to accumulate state
! needed for tracer transport mode
!=======================================================================

   subroutine grist_d2t_coupling_driver(mesh, nlev, idstep, dstep_in_tstep)
! io
     type(global_domain)  , intent(inout) :: mesh
     integer(i4),           intent(in)    :: nlev
     integer(i4),           intent(in)    :: idstep
     integer(i4),           intent(in)    :: dstep_in_tstep
! local
     integer(i4)   :: iv, ie, ilev, icell1, icell2
     real(r8)      :: tmp
     type(exchange_field_list_2d),pointer :: field_head_2d

      field_head_2d=>null()
! 1) set to zero if first step
      if(idstep.eq.1)then
         tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f  = zero
         tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f = zero
         tracerVarCellFull%scalar_delhp_avg_adv%f            = zero
      end if

! 2) do some accumulation
      do ie = 1, mesh%ne ! compute
         do ilev = 1, nlev
             tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie) = tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie)+&
                                                                dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,ie)
         end do
      end do

#ifdef LAM_DOMAIN
      ! at physical boundary, for tracer transport, only provide
      ! outflow, while inflow is zero, so de upwind side always falls within
      ! physical cell domain, shall?
      do ie = 1, mesh%ne ! compute
         icell1 = mesh%edt_v(1,ie)
         icell2 = mesh%edt_v(2,ie)
         if(mesh%plg_mark(icell1).eq.177.and.mesh%plg_mark(icell2).eq.1885)then
            do ilev = 1, nlev
                 if(sign(one,tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie))*mesh%edt_edpNr_edtTg(ie).gt.0)then
                    ! when V1 is out of bdy region, and flow is V1->V2 means inflow, so set to zero
                   ! tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie) = zero
                 end if
              end do
           end if
           if(mesh%plg_mark(icell1).eq.1885.and.mesh%plg_mark(icell2).eq.177)then
              do ilev = 1, nlev
                 if(sign(one,tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie))*mesh%edt_edpNr_edtTg(ie).lt.0)then
                    ! when V2 is out of bdy region, and flow is V2->V1 means inflow, so set to zero
                   ! tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie) = zero
                 end if
              end do
           end if
           if(mesh%edt_mark(ie).eq.1772)then
              do ilev = 1, nlev
                 tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie) = zero
              end do
           end if
      end do
#endif

      do iv = 1, mesh%nv_halo(1)
         do ilev = 1, nlev
            tracerVarCellFull%scalar_delhp_avg_adv%f(ilev,iv) = tracerVarCellFull%scalar_delhp_avg_adv%f(ilev,iv)+&
                                                               dycoreVarCellFull%scalar_delhp_n%f(ilev,iv)
         end do
      end do

! 3) average at the final step
      IF(idstep.eq.dstep_in_tstep)then
         tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f   = tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f/dstep_in_tstep
         tracerVarCellFull%scalar_delhp_avg_adv%f              = tracerVarCellFull%scalar_delhp_avg_adv%f/dstep_in_tstep
         tracerVarCellFull%scalar_delhp_end_adv%f              = dycoreVarCellFull%scalar_delhp_n%f

         do ie = 1, mesh%ne
            icell1 = mesh%edt_v(1,ie)
            icell2 = mesh%edt_v(2,ie) 
            do ilev = 1, nlev
               tmp    = (tracerVarCellFull%scalar_delhp_avg_adv%f(ilev,icell1)+&
                         tracerVarCellFull%scalar_delhp_avg_adv%f(ilev,icell2))*half
               tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f(ilev,ie) = tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie)/tmp
            end do
         end do
#ifndef SEQ_GRIST
! exchange data, tuned balance ok
        call exchange_data_2d_add(mesh,field_head_2d,tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv)
        call exchange_data_2d_add(mesh,field_head_2d,tracerVarEdgeFull%scalar_normal_velocity_avg_adv)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

      END IF
     return
   end subroutine grist_d2t_coupling_driver
!
! accumulate  www&omega and take an average during model/physics step
!
   subroutine grist_dtp_time_average(mesh, nlev, nlevp, idstep, dstep_in_tstep, itstep, tstep_in_mstep)
!io
     type(global_domain)  , intent(inout) :: mesh
     integer(i4),           intent(in)    :: nlev
     integer(i4),           intent(in)    :: nlevp
     integer(i4),           intent(in)    :: idstep
     integer(i4),           intent(in)    :: dstep_in_tstep
     integer(i4),           intent(in)    :: itstep
     integer(i4),           intent(in)    :: tstep_in_mstep
! local
     integer(i4)   :: iv, ie, ilev, icell1, icell2
     real(r8)      :: tmp
     type(exchange_field_list_2d),pointer :: field_head_2d

      field_head_2d=>null()
! 1) set to zero if first step
      if(idstep.eq.1.and.itstep.eq.1)then
        dycoreVarCellFace%scalar_www_timavg%f    = zero
        dycoreVarCellFull%scalar_omega_timavg%f  = zero
      end if
! 2) do some accumulation
      do iv = 1, mesh%nv ! compute
         do ilev = 1, nlev
             dycoreVarCellFull%scalar_omega_timavg%f(ilev,iv) = dycoreVarCellFull%scalar_omega_timavg%f(ilev,iv)+&
                                                               dycoreVarCellFull%scalar_omega_n%f(ilev,iv)
         end do
      end do

#ifdef DYCORE_NHD
      if(nh_dynamics)then
         do iv = 1, mesh%nv_full ! compute
            do ilev = 1, nlevp
               dycoreVarCellFace%scalar_www_timavg%f(ilev,iv) = dycoreVarCellFace%scalar_www_timavg%f(ilev,iv)+&
                                                               dycoreVarCellFace%scalar_www_n%f(ilev,iv)
            end do
         end do
      end if
#endif
! 3) average at the final step
      if(idstep.eq.dstep_in_tstep.and.itstep.eq.tstep_in_mstep)then
         dycoreVarCellFull%scalar_omega_timavg%f = dycoreVarCellFull%scalar_omega_timavg%f/(dstep_in_tstep*tstep_in_mstep)
#ifdef DYCORE_NHD
         if(nh_dynamics) dycoreVarCellFace%scalar_www_timavg%f  = dycoreVarCellFace%scalar_www_timavg%f/(dstep_in_tstep*tstep_in_mstep)
#endif
      endif
#ifndef SEQ_GRIST
! exchange data
        call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_omega_timavg)
        !if(nh_dynamics) call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFace%scalar_www_timavg)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
     return
   end subroutine grist_dtp_time_average

 end module grist_dtp_coupling_module
