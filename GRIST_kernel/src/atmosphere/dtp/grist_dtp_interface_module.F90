
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: This module handles DTP coupling stuff associated with Model 
!              Physics interface, and is not aware of the Dynamics data 
!              structure, i.e., Profile in, Profile out
! Revision history:
!            1. The full PhysPkg is now independent from this module, and 
!               only relies on some public variables here
!----------------------------------------------------------------------------

 module grist_dtp_interface_module
   use grist_mpi
   use grist_constants,                  only: r8, i4, pi, rdry, rvap, cp, p00, omega, gravity, rearth, zero, one
   use grist_domain_types,               only: global_domain
   use grist_dycore_gcd_recon_module_2d, only: vector_recon_perot_edge2cell_uv_2d
   use grist_math_module,                only: convert_vector_sph2cart
   use grist_hpe_constants,              only: deta_full, deta_face, eta_full
! model simple physics
   use grist_physics_idealized_package,  only: grist_idealized_physics_hsdry_exp, &
                                               grist_idealized_physics_hsdry_imp, &
                                               grist_idealized_physics_dcmip2016, &
                                               grist_idealized_physics_mitc     , &
                                               grist_idealized_physics_kessler  , &
                                               grist_idealized_physics_kessler_klemp15
   use grist_physics_kessler_mod,        only: kessler_cam

   implicit none
   private
  
   public :: grist_dtp_interface_hsdry_exp     , &
             grist_dtp_interface_dcmip2016_tc  , &
             grist_dtp_interface_dcmip2016_tc_a, &
             grist_dtp_interface_dcmip2016_sc  , &
             grist_dtp_interface_dcmip2016_sc_a, &
             grist_dtp_interface_dcmip2016_mitc, &
             grist_dtp_interface_init          , &
             grist_dtp_interface_destruct

    real(r8), parameter    :: half  = 0.5_r8
    real(r8), public,allocatable  :: scalar_U_wind_at_pc_full_level(:,:)
    real(r8), public,allocatable  :: scalar_V_wind_at_pc_full_level(:,:)
    real(r8), public,allocatable  :: scalar_qqq_at_pc_full_level(:,:,:)
    real(r8), allocatable  :: tend_U_wind_at_pc_full_level(:,:)
    real(r8), allocatable  :: tend_V_wind_at_pc_full_level(:,:)

  contains

  subroutine grist_dtp_interface_init(mesh,nlev,ntracer)
     type(global_domain) , intent(in) :: mesh
     integer(i4)         , intent(in) :: nlev
     integer(i4)         , intent(in) :: ntracer
     if(.not.allocated(scalar_U_wind_at_pc_full_level)) allocate(scalar_U_wind_at_pc_full_level(nlev,mesh%nv))
     if(.not.allocated(scalar_V_wind_at_pc_full_level)) allocate(scalar_V_wind_at_pc_full_level(nlev,mesh%nv))
     if(.not.allocated(scalar_qqq_at_pc_full_level))    allocate(scalar_qqq_at_pc_full_level(ntracer,nlev,mesh%nv))
     if(.not.allocated(tend_U_wind_at_pc_full_level))   allocate(tend_U_wind_at_pc_full_level(nlev,mesh%nv))
     if(.not.allocated(tend_V_wind_at_pc_full_level))   allocate(tend_V_wind_at_pc_full_level(nlev,mesh%nv))
     return
  end subroutine grist_dtp_interface_init

  subroutine grist_dtp_interface_destruct
     if(allocated(scalar_U_wind_at_pc_full_level))    deallocate(scalar_U_wind_at_pc_full_level)
     if(allocated(scalar_V_wind_at_pc_full_level))    deallocate(scalar_V_wind_at_pc_full_level)
     if(allocated(scalar_qqq_at_pc_full_level))       deallocate(scalar_qqq_at_pc_full_level)
     if(allocated(tend_U_wind_at_pc_full_level))      deallocate(tend_U_wind_at_pc_full_level)
     if(allocated(tend_V_wind_at_pc_full_level))      deallocate(tend_V_wind_at_pc_full_level)
     return
  end subroutine grist_dtp_interface_destruct

  subroutine grist_dtp_interface_hsdry_exp(mesh, nlev, dtime, physpkg, &
                                           scalar_pressure_at_pc_full_level,         &
                                           scalar_pressure_at_pc_face_level,         &
                                           scalar_normal_velocity_at_edge_full_level,&
                                           scalar_mass_pt_at_pc_full_level          ,&
                                           scalar_potential_temp_at_pc_full_level,   &
                                           tend_normal_velocity_at_edge_full_level,  &
                                           tend_potential_temp_at_pc_full_level)
! io
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: nlev
    real(r8)             , intent(in)   :: dtime
    character(len=*)     , intent(in)   :: physpkg
    real(r8),allocatable , intent(in)   :: scalar_pressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_pressure_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_mass_pt_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_potential_temp_at_pc_full_level(:,:)
! local
    integer(i4)                         :: ie, iv, ilev, icell1, icell2
    real(r8)                            :: tend_U_wind_at_edge_full_level, tend_V_wind_at_edge_full_level 
    real(r8)                            :: point_cart_coort(3), vector_velocity(3)

   
!
! 1) Normal wind at edge to UV and cell
!
        mesh%nv = mesh%nv_halo(1)
        call vector_recon_perot_edge2cell_uv_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                                      scalar_U_wind_at_pc_full_level, &
                                                      scalar_V_wind_at_pc_full_level, &
                                                      nlev)
!
! 2) Call held-suarez forcing based on u, v, potential temperature
!    generate all physics tendencies at cell
        select case(trim(physpkg))
        case('HSDRY_EXP')
        call grist_idealized_physics_hsdry_exp(mesh, nlev                              , &
                                               scalar_pressure_at_pc_full_level        , &
                                               scalar_pressure_at_pc_face_level        , &
                                               scalar_U_wind_at_pc_full_level          , &
                                               scalar_V_wind_at_pc_full_level          , &
                                               scalar_potential_temp_at_pc_full_level  , &
                                               tend_U_wind_at_pc_full_level            , &
                                               tend_V_wind_at_pc_full_level            , &
                                               tend_potential_temp_at_pc_full_level)
        case('HSDRY_IMP')
        call grist_idealized_physics_hsdry_imp(mesh, nlev, dtime                       , &
                                               scalar_pressure_at_pc_full_level        , &
                                               scalar_pressure_at_pc_face_level        , &
                                               scalar_U_wind_at_pc_full_level          , &
                                               scalar_V_wind_at_pc_full_level          , &
                                               scalar_mass_pt_at_pc_full_level         , &
                                               scalar_potential_temp_at_pc_full_level  , &
                                               tend_U_wind_at_pc_full_level            , &
                                               tend_V_wind_at_pc_full_level            , &
                                               tend_potential_temp_at_pc_full_level)
        case default
             print*, "grist_dtp_interface_hsdry_exp: you must select physpkg, stop" 
             stop
        end select

        mesh%nv = mesh%nv_compute
!
! 3) transform u,v tendency from cell to normal wind tendency at edge
!
     do ie = 1, mesh%ne
        !icell1 = mesh%edt(ie)%v(1)
        !icell2 = mesh%edt(ie)%v(2)
        icell1 = mesh%edt_v(1,ie)
        icell2 = mesh%edt_v(2,ie)
        point_cart_coort(:) = mesh%edt_c_p(1:3,ie)
        do ilev = 1, nlev
           tend_U_wind_at_edge_full_level = half*(tend_U_wind_at_pc_full_level(ilev,icell1)+tend_U_wind_at_pc_full_level(ilev,icell2))
           tend_V_wind_at_edge_full_level = half*(tend_V_wind_at_pc_full_level(ilev,icell1)+tend_V_wind_at_pc_full_level(ilev,icell2))
           call convert_vector_sph2cart(tend_U_wind_at_edge_full_level,&
                                        tend_V_wind_at_edge_full_level,&
                                        point_cart_coort, vector_velocity)
           tend_normal_velocity_at_edge_full_level(ilev,ie) = dot_product(vector_velocity, mesh%edp_nr(1:3,ie))
           
        end do
     end do

    return
  end subroutine grist_dtp_interface_hsdry_exp

!====================================================================
! Simple moist physics of RJ2012/DCMIP2016, for such an interface,
! both input and output style follow the dycore and tracer data,
! but the compound data tendency should be computed in the calling
! note that for NDC, hpressure and pressure differ, here we still
! input hpressure to simple physics because some assumptions used
! inside the simple physics, and we assume delhp remain unchanged
!====================================================================

  subroutine grist_dtp_interface_dcmip2016_tc(mesh, nlev, ncell, dtime                 , &
                                              scalar_hpressure_at_pc_full_level        , &
                                              scalar_hpressure_at_pc_face_level        , &
                                              scalar_delhp_at_pc_full_level            , &
                                              scalar_normal_velocity_at_edge_full_level, &
                                              scalar_potential_temp_at_pc_full_level   , & !thetam
                                              scalar_temp_at_pc_full_level             , &
                                              scalar_tracer_mxrt_at_pc_full_level      , &
                                              tend_normal_velocity_at_edge_full_level  , &
                                              tend_potential_temp_at_pc_full_level     , &
                                              tend_tracer_mxrt_at_pc_full_level        , &
                                              scalar_precl_surface)
! io
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: nlev
    integer(i4)          , intent(in)   :: ncell
    real(r8)             , intent(in)   :: dtime
    real(r8),allocatable , intent(in)   :: scalar_hpressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_hpressure_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_delhp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: tend_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: scalar_precl_surface(:)
! local
    real(r8)                            :: lats(ncell)
    real(r8)                            :: local_pmid(ncell,nlev)    ! hpressure plus moisture
    real(r8)                            :: local_pint(ncell,nlev+1)  ! hpressure plus moisture
    real(r8)                            :: local_pdel(ncell,nlev)    ! hpressure plus moisture
    real(r8)                            :: local_rpdel(ncell,nlev)   ! hpressure plus moisture
    real(r8)                            :: local_ps(ncell)           ! hpressure plus moisture
    real(r8)                            :: local_uuu(ncell,nlev)
    real(r8)                            :: local_vvv(ncell,nlev)
    real(r8)                            :: local_ttt(ncell,nlev)
    real(r8)                            :: local_qqq(ncell,nlev)
    real(r8)                            :: local_qqq_face(ncell,nlev+1)
    real(r8)                            :: local_precl(ncell)
    real(r8)                            :: local_dudt(ncell,nlev)
    real(r8)                            :: local_dvdt(ncell,nlev)
    real(r8)                            :: local_dtdt(ncell,nlev)
    real(r8)                            :: local_dqdt(ncell,nlev)
    real(r8)                            :: tend_U_wind_at_edge_full_level, tend_V_wind_at_edge_full_level 
    real(r8)                            :: point_cart_coort(3), vector_velocity(3)
    integer(i4)                         :: ie, iv, ilev, icell1, icell2, nlevp

        nlevp = nlev+1
!
! 1) prepare data
!
        if(ncell.ne.mesh%nv_halo(1)) then
           print*, "dim stuff wrong in grist_dtp_interface_dcmip2016_tc, model aborts"
           call mpi_abort()
        end if

        mesh%nv = mesh%nv_halo(1)
        call vector_recon_perot_edge2cell_uv_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                                      scalar_U_wind_at_pc_full_level, &
                                                      scalar_V_wind_at_pc_full_level, &
                                                      nlev)
        do iv = 1, ncell
           lats(iv) = mesh%vtx_lat(iv)
        end do

        scalar_qqq_at_pc_full_level = scalar_tracer_mxrt_at_pc_full_level
! get a local copy
        local_uuu(1:ncell,1:nlev) = transpose(scalar_U_wind_at_pc_full_level(1:nlev,1:ncell))
        local_vvv(1:ncell,1:nlev) = transpose(scalar_V_wind_at_pc_full_level(1:nlev,1:ncell))
        local_ttt(1:ncell,1:nlev) = transpose(  scalar_temp_at_pc_full_level(1:nlev,1:ncell))
        local_qqq(1:ncell,1:nlev) = transpose(   scalar_qqq_at_pc_full_level(1,1:nlev,1:ncell))
        do ilev = 2, nlev
           local_qqq_face(:,ilev) = 0.5*(deta_full(ilev-1)/deta_face(ilev)*local_qqq(:,ilev)+&
                                         deta_full(ilev)  /deta_face(ilev)*local_qqq(:,ilev-1))
        end do
        local_qqq_face(:,1)     = 2._r8*local_qqq(:,1)-local_qqq_face(:,2)
        local_qqq_face(:,nlevp) = 2._r8*local_qqq(:,nlev)-local_qqq_face(:,nlev)

! evaluate moist hydro-pressure for physics
        do iv = 1, ncell
           do ilev = 1, nlev
              local_pmid(iv,ilev) = scalar_hpressure_at_pc_full_level(ilev,iv)*(one+rvap/rdry*local_qqq(iv,ilev))
              local_pdel(iv,ilev) = scalar_delhp_at_pc_full_level(ilev,iv)*(one+rvap/rdry*local_qqq(iv,ilev))
              local_rpdel(iv,ilev)= one/local_pdel(iv,ilev)
           end do
           do ilev = 1, nlevp
              local_pint(iv,ilev)= scalar_hpressure_at_pc_face_level(ilev,iv)*(one+rvap/rdry*local_qqq_face(iv,ilev))
           end do
           local_ps(iv)  = local_pint(iv,nlevp)
        end do
!
! 2) call interface for tendency, local copies of state are ignored
! 
        ! save ini
        call grist_idealized_physics_dcmip2016(ncell, nlev, dtime, lats , &
                                               local_ttt, local_qqq, local_uuu, local_vvv, &
                                               local_pmid, local_pint, local_pdel, local_rpdel, local_ps, &
                                               local_precl, 0, .true., .false.                  , &
                                               local_dudt, local_dvdt, local_dtdt, local_dqdt)

        scalar_precl_surface(1:ncell) = local_precl(1:ncell)

        mesh%nv = mesh%nv_compute
!
! 3) transform u,v tendency from cell to normal wind tendency at edge
!
        do ie = 1, mesh%ne
           !icell1 = mesh%edt(ie)%v(1)
           !icell2 = mesh%edt(ie)%v(2)
           icell1 = mesh%edt_v(1,ie)
           icell2 = mesh%edt_v(2,ie)
           point_cart_coort(:) = mesh%edt_c_p(1:3,ie)
           do ilev = 1, nlev
              tend_U_wind_at_edge_full_level = half*(local_dudt(icell1,ilev)+local_dudt(icell2,ilev))
              tend_V_wind_at_edge_full_level = half*(local_dvdt(icell1,ilev)+local_dvdt(icell2,ilev))
              call convert_vector_sph2cart(tend_U_wind_at_edge_full_level,&
                                           tend_V_wind_at_edge_full_level,&
                                           point_cart_coort, vector_velocity)
              tend_normal_velocity_at_edge_full_level(ilev,ie) = dot_product(vector_velocity, mesh%edp_nr(1:3,ie))
           end do
        end do

        tend_tracer_mxrt_at_pc_full_level(1,1:nlev,1:ncell) = transpose(local_dqdt(1:ncell,1:nlev))

! find thetam tendency
        do iv = 1, ncell
           do ilev = 1, nlev
! old
!             tend_potential_temp_at_pc_full_level(ilev,iv) = local_dtdt(iv,ilev)*&
!                                                            ((p00/local_pmid(iv,ilev))**(rdry/cp))*&
!                                                            (one+(rvap/rdry)*scalar_qqq_at_pc_full_level(ilev,iv))
! (1) entire evaluation
!              tend_potential_temp_at_pc_full_level(ilev,iv) = (local_ttt(iv,ilev)*((p00/local_pmid(iv,ilev))**(rdry/cp))*(one+(rvap/rdry)*local_qqq(iv,ilev))-&
!                                                               scalar_potential_temp_at_pc_full_level(ilev,iv))/dtime
! (2) partial
             tend_potential_temp_at_pc_full_level(ilev,iv) = ((p00/local_pmid(iv,ilev))**(rdry/cp))*&
                                                              (local_dtdt(iv,ilev)*(one+(rvap/rdry)*scalar_qqq_at_pc_full_level(1,ilev,iv))+&
                                                               scalar_temp_at_pc_full_level(ilev,iv)*(rvap/rdry)*local_dqdt(iv,ilev))
           end do
        end do

      return
  end subroutine grist_dtp_interface_dcmip2016_tc

!====================================================================
! as above but use mpressure, theoratically better
!====================================================================

  subroutine grist_dtp_interface_dcmip2016_tc_a(mesh, nlev, ncell, dtime                 , &
                                              scalar_mpressure_at_pc_full_level        , &
                                              scalar_mpressure_at_pc_face_level        , &
                                              scalar_normal_velocity_at_edge_full_level, &
                                              scalar_potential_temp_at_pc_full_level   , & !thetam
                                              scalar_temp_at_pc_full_level             , &
                                              scalar_tracer_mxrt_at_pc_full_level      , &
                                              phys_flag                                , &
                                              TC_pbl_flag                              , &
                                              tend_normal_velocity_at_edge_full_level  , &
                                              tend_potential_temp_at_pc_full_level     , &
                                              tend_tracer_mxrt_at_pc_full_level        , &
                                              scalar_precl_surface)
! io
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: nlev
    integer(i4)          , intent(in)   :: ncell
    real(r8)             , intent(in)   :: dtime
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
    integer(i4)          , intent(in)   :: phys_flag
    logical              , intent(in)   :: TC_pbl_flag
    real(r8),allocatable , intent(inout):: tend_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: scalar_precl_surface(:)
! local
    real(r8)                            :: lats(ncell)
    real(r8)                            :: local_pmid(ncell,nlev)    ! hpressure plus moisture
    real(r8)                            :: local_pint(ncell,nlev+1)  ! hpressure plus moisture
    real(r8)                            :: local_pdel(ncell,nlev)    ! hpressure plus moisture
    real(r8)                            :: local_rpdel(ncell,nlev)   ! hpressure plus moisture
    real(r8)                            :: local_ps(ncell)           ! hpressure plus moisture
    real(r8)                            :: local_uuu(ncell,nlev)
    real(r8)                            :: local_vvv(ncell,nlev)
    real(r8)                            :: local_ttt(ncell,nlev)
    real(r8)                            :: local_qqq(ncell,nlev)
    real(r8)                            :: local_precl(ncell)
    real(r8)                            :: local_dudt(ncell,nlev)
    real(r8)                            :: local_dvdt(ncell,nlev)
    real(r8)                            :: local_dtdt(ncell,nlev)
    real(r8)                            :: local_dqdt(ncell,nlev)
    real(r8)                            :: tend_U_wind_at_edge_full_level, tend_V_wind_at_edge_full_level 
    real(r8)                            :: point_cart_coort(3), vector_velocity(3)
    integer(i4)                         :: ie, iv, ilev, icell1, icell2, nlevp

        nlevp = nlev+1
!
! 1) prepare data
!
        if(ncell.ne.mesh%nv_halo(1)) then
           print*, "dim stuff wrong in grist_dtp_interface_dcmip2016_tc, model aborts"
           call mpi_abort()
        end if
! obtain UV
        mesh%nv = mesh%nv_halo(1)
        call vector_recon_perot_edge2cell_uv_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                                      scalar_U_wind_at_pc_full_level, &
                                                      scalar_V_wind_at_pc_full_level, &
                                                      nlev)
        do iv = 1, ncell
           lats(iv) = mesh%vtx_lat(iv)
        end do

! from dry to moist
        scalar_qqq_at_pc_full_level = scalar_tracer_mxrt_at_pc_full_level(:,:,:)/(one+scalar_tracer_mxrt_at_pc_full_level(:,:,:))

! get a local copy
        local_uuu(1:ncell,1:nlev) = transpose(scalar_U_wind_at_pc_full_level(1:nlev,1:ncell))
        local_vvv(1:ncell,1:nlev) = transpose(scalar_V_wind_at_pc_full_level(1:nlev,1:ncell))
        local_ttt(1:ncell,1:nlev) = transpose(  scalar_temp_at_pc_full_level(1:nlev,1:ncell))
        local_qqq(1:ncell,1:nlev) = transpose(   scalar_qqq_at_pc_full_level(1,1:nlev,1:ncell))
! evaluate moist hydro-pressure for physics
        do iv = 1, ncell
           do ilev = 1, nlev
              local_pmid(iv,ilev) = scalar_mpressure_at_pc_full_level(ilev,iv)
              local_pdel(iv,ilev) = scalar_mpressure_at_pc_face_level(ilev+1,iv)-scalar_mpressure_at_pc_face_level(ilev,iv)
              local_rpdel(iv,ilev)= one/local_pdel(iv,ilev)
           end do
           do ilev = 1, nlevp
              local_pint(iv,ilev)= scalar_mpressure_at_pc_face_level(ilev,iv)
           end do
           local_ps(iv)  = local_pint(iv,nlevp)
        end do
!
! 2) call interface for tendency, local copies of state are ignored
! 
        ! save ini
        call grist_idealized_physics_dcmip2016(ncell, nlev, dtime, lats , &
                                               local_ttt, local_qqq, local_uuu, local_vvv, &
                                               local_pmid, local_pint, local_pdel, local_rpdel, local_ps, &
                                               local_precl, phys_flag, .true., TC_pbl_flag              , &
                                               local_dudt, local_dvdt, local_dtdt, local_dqdt)

        scalar_precl_surface(1:ncell) = local_precl(1:ncell)

        mesh%nv = mesh%nv_compute
!
! 3) transform u,v tendency from cell to normal wind tendency at edge
!
        do ie = 1, mesh%ne
           !icell1 = mesh%edt(ie)%v(1)
           !icell2 = mesh%edt(ie)%v(2)
           icell1 = mesh%edt_v(1,ie)
           icell2 = mesh%edt_v(2,ie)
           point_cart_coort(:) = mesh%edt_c_p(1:3,ie)
           do ilev = 1, nlev
              tend_U_wind_at_edge_full_level = half*(local_dudt(icell1,ilev)+local_dudt(icell2,ilev))
              tend_V_wind_at_edge_full_level = half*(local_dvdt(icell1,ilev)+local_dvdt(icell2,ilev))
              call convert_vector_sph2cart(tend_U_wind_at_edge_full_level,&
                                           tend_V_wind_at_edge_full_level,&
                                           point_cart_coort, vector_velocity)
              tend_normal_velocity_at_edge_full_level(ilev,ie) = dot_product(vector_velocity, mesh%edp_nr(1:3,ie))
           end do
        end do

!        tend_tracer_mxrt_at_pc_full_level(1,1:nlev,1:ncell) = transpose(local_dqdt(1:ncell,1:nlev))

! find mxrt and thetam tendency
        do iv = 1, ncell
           do ilev = 1, nlev
! moist to dry tend
              tend_tracer_mxrt_at_pc_full_level(1,ilev,iv) = local_dqdt(iv,ilev)/((one-scalar_qqq_at_pc_full_level(1,ilev,iv))**2)
! old
!             tend_potential_temp_at_pc_full_level(ilev,iv) = local_dtdt(iv,ilev)*&
!                                                            ((p00/local_pmid(iv,ilev))**(rdry/cp))*&
!                                                            (one+(rvap/rdry)*scalar_qqq_at_pc_full_level(ilev,iv))
! (1) entire evaluation
!              tend_potential_temp_at_pc_full_level(ilev,iv) = (local_ttt(iv,ilev)*((p00/local_pmid(iv,ilev))**(rdry/cp))*(one+(rvap/rdry)*local_qqq(iv,ilev))-&
!                                                               scalar_potential_temp_at_pc_full_level(ilev,iv))/dtime
! (2) partial
             tend_potential_temp_at_pc_full_level(ilev,iv) = ((p00/local_pmid(iv,ilev))**(rdry/cp))*&
                                                              (local_dtdt(iv,ilev)*(one+(rvap/rdry)*scalar_tracer_mxrt_at_pc_full_level(1,ilev,iv))+&
                                                               scalar_temp_at_pc_full_level(ilev,iv)*(rvap/rdry)*tend_tracer_mxrt_at_pc_full_level(1,ilev,iv))
           end do
        end do

      return
  end subroutine grist_dtp_interface_dcmip2016_tc_a

!=====================================================
! DCMIP-2016 supercell & Klemp et al. (2015), JAMES
! Kessler Micro Physics
!=====================================================

  subroutine grist_dtp_interface_dcmip2016_sc(mesh, nlev, ncell, dtime                 , &
                                             scalar_hpressure_at_pc_full_level        , &
                                             scalar_geopotential_at_pc_face_level     , &
                                             scalar_geopotential_at_pc_full_level     , &
                                             scalar_delhp_at_pc_full_level            , &
                                             scalar_potential_temp_at_pc_full_level   , & !thetam
                                             scalar_tracer_mxrt_at_pc_full_level      , &
                                             tend_potential_temp_at_pc_full_level     , &
                                             tend_tracer_mxrt_at_pc_full_level        , &
                                             scalar_precl_surface)
! iostream
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: nlev
    integer(i4)          , intent(in)   :: ncell
    real(r8)             , intent(in)   :: dtime
    real(r8),allocatable , intent(in)   :: scalar_hpressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_delhp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: tend_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: scalar_precl_surface(:)
!
! local data for driving the kessler micro physics, from bottom to top
!
    real(r8)     :: local_theta(nlev), local_theta_ini(nlev), local_thetam_next(nlev)
    real(r8)     :: local_qv(nlev), local_qv_ini(nlev)
    real(r8)     :: local_qc(nlev)
    real(r8)     :: local_qr(nlev)
    real(r8)     :: local_dryrho(nlev)
    real(r8)     :: local_exner(nlev)
    real(r8)     :: local_height(nlev)
    real(r8)     :: local_precl
    real(r8)     :: local_dtheta_dt(nlev)
    real(r8)     :: local_dqv_dt(nlev)
    real(r8)     :: local_dqc_dt(nlev)
    real(r8)     :: local_dqr_dt(nlev)
    integer(i4)  :: iv, ilev

    do iv = 1, ncell

      do ilev = 1, nlev
         local_qv(nlev+1-ilev)    = scalar_tracer_mxrt_at_pc_full_level(1,ilev,iv)
         local_qc(nlev+1-ilev)    = scalar_tracer_mxrt_at_pc_full_level(2,ilev,iv)
         local_qr(nlev+1-ilev)    = scalar_tracer_mxrt_at_pc_full_level(3,ilev,iv)
! Kessler requires raw theta
         local_theta(nlev+1-ilev) = scalar_potential_temp_at_pc_full_level(ilev,iv)/(one+rvap/rdry*local_qv(nlev+1-ilev))
         local_dryrho(nlev+1-ilev)= scalar_delhp_at_pc_full_level(ilev,iv)/&
                                   (scalar_geopotential_at_pc_face_level(ilev,iv)-&
                                    scalar_geopotential_at_pc_face_level(ilev+1,iv))

         local_exner(nlev+1-ilev) = (scalar_hpressure_at_pc_full_level(ilev,iv)*(one+rvap/rdry*local_qv(nlev+1-ilev))/p00)**(rdry/cp)
         local_height(nlev+1-ilev)= scalar_geopotential_at_pc_full_level(ilev,iv)/gravity
      end do
!
! call Kessler microphysics
!
      local_theta_ini = local_theta
      local_qv_ini    = local_qv
! dcmip2016 routine with rainsplit, why bad?
#ifdef DCMIP_KESSLER
      call grist_idealized_physics_kessler(local_theta, local_qv, local_qc, local_qr, & !
#else
! klemp 15 routine without rainsplit, good
      call grist_idealized_physics_kessler_klemp15(local_theta, local_qv, local_qc, local_qr, &
#endif
                                           local_dryrho,local_exner, dtime, local_height, &
                                           nlev, local_precl, &
                                           local_dtheta_dt, local_dqv_dt, local_dqc_dt, local_dqr_dt)
! CESM/CAM routine with rainsplit, same as dcmip's one but with modification to avoid some
! machine problems; this one should be used with BW case due to better stability
#ifdef CESM_KESSLER
      call kessler_cam(nlev, dtime, local_dryrho, local_height, local_exner, local_theta, local_qv, local_qc, local_qr, local_precl, &
                       local_dtheta_dt, local_dqv_dt, local_dqc_dt, local_dqr_dt)
#endif

!
! transform to final time tendency for grid-scale prognostics
!
      local_thetam_next(1:nlev)= local_theta(1:nlev)*(one+rvap/rdry*local_qv(1:nlev))

      do ilev = 1, nlev
![1] entire evaluation
         tend_potential_temp_at_pc_full_level(ilev,iv) = (local_thetam_next(nlev+1-ilev)-scalar_potential_temp_at_pc_full_level(ilev,iv))/dtime
![2] partial evaluation, prefer
         tend_potential_temp_at_pc_full_level(ilev,iv) = (one+rvap/rdry*local_qv_ini(nlev+1-ilev))*local_dtheta_dt(nlev+1-ilev)+&
                                                         (rvap/rdry*local_theta_ini(nlev+1-ilev)*local_dqv_dt(nlev+1-ilev))
         tend_tracer_mxrt_at_pc_full_level(1,ilev,iv)  = local_dqv_dt(nlev+1-ilev)
         tend_tracer_mxrt_at_pc_full_level(2,ilev,iv)  = local_dqc_dt(nlev+1-ilev)
         tend_tracer_mxrt_at_pc_full_level(3,ilev,iv)  = local_dqr_dt(nlev+1-ilev)
      end do
! send back for diagnose
      scalar_precl_surface(iv) = local_precl
    end do
    return
  end subroutine grist_dtp_interface_dcmip2016_sc


!=====================================================
! DCMIP-2016 supercell & Klemp et al. (2015), JAMES
! Kessler Micro Physics, use pressure
!=====================================================

  subroutine grist_dtp_interface_dcmip2016_sc_a(mesh, nlev, ncell, dtime                 , &
                                                scalar_mpressure_at_pc_full_level        , &
                                                scalar_geopotential_at_pc_face_level     , &
                                                scalar_geopotential_at_pc_full_level     , &
                                                scalar_delhp_at_pc_full_level            , &
                                                scalar_potential_temp_at_pc_full_level   , & !thetam
                                                scalar_tracer_mxrt_at_pc_full_level      , &
                                                tend_potential_temp_at_pc_full_level     , &
                                                tend_tracer_mxrt_at_pc_full_level        , &
                                                scalar_precl_surface)
! iostream
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: nlev
    integer(i4)          , intent(in)   :: ncell
    real(r8)             , intent(in)   :: dtime
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_geopotential_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_delhp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: tend_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: scalar_precl_surface(:)
!
! local data for driving the kessler micro physics, from bottom to top
!
    real(r8)     :: local_theta(nlev), local_theta_ini(nlev), local_thetam_next(nlev)
    real(r8)     :: local_qv(nlev), local_qv_ini(nlev)
    real(r8)     :: local_qc(nlev)
    real(r8)     :: local_qr(nlev)
    real(r8)     :: local_dryrho(nlev)
    real(r8)     :: local_exner(nlev)
    real(r8)     :: local_height(nlev)
    real(r8)     :: local_precl
    real(r8)     :: local_dtheta_dt(nlev)
    real(r8)     :: local_dqv_dt(nlev)
    real(r8)     :: local_dqc_dt(nlev)
    real(r8)     :: local_dqr_dt(nlev)
    integer(i4)  :: iv, ilev

    do iv = 1, ncell

      do ilev = 1, nlev
         local_qv(nlev+1-ilev)    = scalar_tracer_mxrt_at_pc_full_level(1,ilev,iv)
         local_qc(nlev+1-ilev)    = scalar_tracer_mxrt_at_pc_full_level(2,ilev,iv)
         local_qr(nlev+1-ilev)    = scalar_tracer_mxrt_at_pc_full_level(3,ilev,iv)
! Kessler requires raw theta
         local_theta(nlev+1-ilev) = scalar_potential_temp_at_pc_full_level(ilev,iv)/(one+rvap/rdry*local_qv(nlev+1-ilev))
         local_dryrho(nlev+1-ilev)= scalar_delhp_at_pc_full_level(ilev,iv)/&
                                   (scalar_geopotential_at_pc_face_level(ilev,iv)-&
                                    scalar_geopotential_at_pc_face_level(ilev+1,iv))

         local_exner(nlev+1-ilev) = (scalar_mpressure_at_pc_full_level(ilev,iv)/p00)**(rdry/cp)
         local_height(nlev+1-ilev)= scalar_geopotential_at_pc_full_level(ilev,iv)/gravity
      end do
!
! call Kessler microphysics
!
      local_theta_ini = local_theta
      local_qv_ini    = local_qv
! dcmip2016 routine with rainsplit, why bad?
!      call grist_idealized_physics_kessler(local_theta, local_qv, local_qc, local_qr, & !
! klemp 15 routine without rainsplit, good
      call grist_idealized_physics_kessler_klemp15(local_theta, local_qv, local_qc, local_qr, &
                                           local_dryrho,local_exner, dtime, local_height, &
                                           nlev, local_precl, &
                                           local_dtheta_dt, local_dqv_dt, local_dqc_dt, local_dqr_dt)
! CESM/CAM routine with rainsplit, same as dcmip's one but with modification to avoid some
! machine problems; this one should be used with BW case due to better stability
#ifdef CESM_KESSLER
      call kessler_cam(nlev, dtime, local_dryrho, local_height, local_exner, local_theta, local_qv, local_qc, local_qr, local_precl, &
                       local_dtheta_dt, local_dqv_dt, local_dqc_dt, local_dqr_dt)
#endif

!
! transform to final time tendency for grid-scale prognostics
!
!      local_thetam_next(1:nlev)= local_theta(1:nlev)*(one+rvap/rdry*local_qv(1:nlev))

      do ilev = 1, nlev
![1] entire evaluation
!         tend_potential_temp_at_pc_full_level(ilev,iv) = (local_thetam_next(nlev+1-ilev)-scalar_potential_temp_at_pc_full_level(ilev,iv))/dtime
![2] partial evaluation, prefer
         tend_potential_temp_at_pc_full_level(ilev,iv) = (one+rvap/rdry*local_qv_ini(nlev+1-ilev))*local_dtheta_dt(nlev+1-ilev)+&
                                                         (rvap/rdry*local_theta_ini(nlev+1-ilev)*local_dqv_dt(nlev+1-ilev))
         tend_tracer_mxrt_at_pc_full_level(1,ilev,iv)  = local_dqv_dt(nlev+1-ilev)
         tend_tracer_mxrt_at_pc_full_level(2,ilev,iv)  = local_dqc_dt(nlev+1-ilev)
         tend_tracer_mxrt_at_pc_full_level(3,ilev,iv)  = local_dqr_dt(nlev+1-ilev)
      end do
! send back for diagnose
      scalar_precl_surface(iv) = local_precl
    end do
    return
  end subroutine grist_dtp_interface_dcmip2016_sc_a

!====================================================================
! similar to tc_a but use mitc package
!====================================================================

  subroutine grist_dtp_interface_dcmip2016_mitc(mesh, nlev, ncell, dtime                 , &
                                                scalar_mpressure_at_pc_full_level        , &
                                                scalar_mpressure_at_pc_face_level        , &
                                                scalar_normal_velocity_at_edge_full_level, &
                                                scalar_potential_temp_at_pc_full_level   , & !thetam
                                                scalar_temp_at_pc_full_level             , &
                                                scalar_tracer_mxrt_at_pc_full_level      , &
                                                tend_normal_velocity_at_edge_full_level  , &
                                                tend_potential_temp_at_pc_full_level     , &
                                                tend_tracer_mxrt_at_pc_full_level        , &
                                                scalar_precl_surface)
! io
    type(global_domain)  , intent(inout):: mesh
    integer(i4)          , intent(in)   :: nlev
    integer(i4)          , intent(in)   :: ncell
    real(r8)             , intent(in)   :: dtime
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_mpressure_at_pc_face_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(in)   :: scalar_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: tend_normal_velocity_at_edge_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_potential_temp_at_pc_full_level(:,:)
    real(r8),allocatable , intent(inout):: tend_tracer_mxrt_at_pc_full_level(:,:,:)
    real(r8),allocatable , intent(inout):: scalar_precl_surface(:)
! local
    real(r8)                            :: lats(ncell)
    real(r8)                            :: local_pmid(ncell,nlev)    ! hpressure plus moisture
    real(r8)                            :: local_pint(ncell,nlev+1)  ! hpressure plus moisture
    real(r8)                            :: local_pdel(ncell,nlev)    ! hpressure plus moisture
    real(r8)                            :: local_rpdel(ncell,nlev)   ! hpressure plus moisture
    real(r8)                            :: local_ps(ncell)           ! hpressure plus moisture
    real(r8)                            :: local_uuu(ncell,nlev)
    real(r8)                            :: local_vvv(ncell,nlev)
    real(r8)                            :: local_ttt(ncell,nlev)
    real(r8)                            :: local_qqq(ncell,nlev)
    real(r8)                            :: local_precl(ncell)
    real(r8)                            :: local_dudt(ncell,nlev)
    real(r8)                            :: local_dvdt(ncell,nlev)
    real(r8)                            :: local_dtdt(ncell,nlev)
    real(r8)                            :: local_dqdt(ncell,nlev)
    real(r8)                            :: tend_U_wind_at_edge_full_level, tend_V_wind_at_edge_full_level 
    real(r8)                            :: point_cart_coort(3), vector_velocity(3)
    integer(i4)                         :: ie, iv, ilev, icell1, icell2, nlevp

        nlevp = nlev+1
!
! prepare data
!
        if(ncell.ne.mesh%nv_halo(1)) then
           print*, "dim stuff wrong in grist_dtp_interface_dcmip2016_mitc, model aborts"
           call mpi_abort()
        end if
! obtain UV
        mesh%nv = mesh%nv_halo(1)
        call vector_recon_perot_edge2cell_uv_2d(mesh, scalar_normal_velocity_at_edge_full_level, &
                                                      scalar_U_wind_at_pc_full_level, &
                                                      scalar_V_wind_at_pc_full_level, &
                                                      nlev)
        do iv = 1, ncell
           lats(iv) = mesh%vtx_lat(iv)
        end do

! from dry to moist
        scalar_qqq_at_pc_full_level = scalar_tracer_mxrt_at_pc_full_level(:,:,:)/(one+scalar_tracer_mxrt_at_pc_full_level(:,:,:))

! get a local copy
        local_uuu(1:ncell,1:nlev) = transpose(scalar_U_wind_at_pc_full_level(1:nlev,1:ncell))
        local_vvv(1:ncell,1:nlev) = transpose(scalar_V_wind_at_pc_full_level(1:nlev,1:ncell))
        local_ttt(1:ncell,1:nlev) = transpose(  scalar_temp_at_pc_full_level(1:nlev,1:ncell))
        local_qqq(1:ncell,1:nlev) = transpose(   scalar_qqq_at_pc_full_level(1,1:nlev,1:ncell))
! evaluate moist hydro-pressure for physics
        do iv = 1, ncell
           do ilev = 1, nlev
              local_pmid(iv,ilev) = scalar_mpressure_at_pc_full_level(ilev,iv)
              local_pdel(iv,ilev) = scalar_mpressure_at_pc_face_level(ilev+1,iv)-scalar_mpressure_at_pc_face_level(ilev,iv)
              local_rpdel(iv,ilev)= one/local_pdel(iv,ilev)
           end do
           do ilev = 1, nlevp
              local_pint(iv,ilev)= scalar_mpressure_at_pc_face_level(ilev,iv)
           end do
           local_ps(iv)  = local_pint(iv,nlevp)
        end do
!
! call interface for tendency, local copies of state are ignored
! 
        call grist_idealized_physics_mitc(ncell, nlev, dtime, lats, &
                                          local_ttt,  local_qqq,  local_uuu, local_vvv, &
                                          local_pmid, local_pint, real(eta_full,r8),  local_ps, local_precl, &
                                          local_dtdt, local_dqdt, local_dudt,local_dvdt)

        scalar_precl_surface(1:ncell) = local_precl(1:ncell)

        mesh%nv = mesh%nv_compute
!
! transform u,v tendency from cell to normal wind tendency at edge
!
        do ie = 1, mesh%ne
           !icell1 = mesh%edt(ie)%v(1)
           !icell2 = mesh%edt(ie)%v(2)
           icell1 = mesh%edt_v(1,ie)
           icell2 = mesh%edt_v(2,ie)
           point_cart_coort(:) = mesh%edt_c_p(1:3,ie)
           do ilev = 1, nlev
              tend_U_wind_at_edge_full_level = half*(local_dudt(icell1,ilev)+local_dudt(icell2,ilev))
              tend_V_wind_at_edge_full_level = half*(local_dvdt(icell1,ilev)+local_dvdt(icell2,ilev))
              call convert_vector_sph2cart(tend_U_wind_at_edge_full_level,&
                                           tend_V_wind_at_edge_full_level,&
                                           point_cart_coort, vector_velocity)
              tend_normal_velocity_at_edge_full_level(ilev,ie) = dot_product(vector_velocity, mesh%edp_nr(1:3,ie))
           end do
        end do

! find mxrt and thetam tendency
        do iv = 1, ncell
           do ilev = 1, nlev
! moist to dry tend
              tend_tracer_mxrt_at_pc_full_level(1,ilev,iv) = local_dqdt(iv,ilev)/((one-scalar_qqq_at_pc_full_level(1,ilev,iv))**2)
! (1) entire evaluation
!              tend_potential_temp_at_pc_full_level(ilev,iv) = (local_ttt(iv,ilev)*((p00/local_pmid(iv,ilev))**(rdry/cp))*(one+(rvap/rdry)*local_qqq(iv,ilev))-&
!                                                               scalar_potential_temp_at_pc_full_level(ilev,iv))/dtime
! (2) partial
             tend_potential_temp_at_pc_full_level(ilev,iv) = ((p00/local_pmid(iv,ilev))**(rdry/cp))*&
                                                              (local_dtdt(iv,ilev)*(one+(rvap/rdry)*scalar_tracer_mxrt_at_pc_full_level(1,ilev,iv))+&
                                                               scalar_temp_at_pc_full_level(ilev,iv)*(rvap/rdry)*tend_tracer_mxrt_at_pc_full_level(1,ilev,iv))
           end do
        end do

      return
  end subroutine grist_dtp_interface_dcmip2016_mitc

 end module grist_dtp_interface_module
