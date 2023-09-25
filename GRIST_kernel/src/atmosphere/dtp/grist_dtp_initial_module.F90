
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: This module initialize the dtp mode
!          (1) DCMIP2016, Baroclinic wave (Ullrich et al. 2014, QJRMS)
!          (2) DCMIP2016, Idealized tropical cyclone (Reed et al. 2010/2012, MWR/JAMES)
!          (3) DCMIP2016, Supercell thunderstorm (Klemp et al. 2015, JAMES)
! Revision history: 
!       1. 20190622: Be careful with the not-small sensitity associated with various
!                    initialization approach
!       2. 201909: for real case, the condition depends on the source of 
!                  input data, including ERAIM, ERAIP, GFSP, ...
!----------------------------------------------------------------------------

 module grist_dtp_initial_module

! infrastructure
  use grist_constants,             only: gravity, one, pi, i4, r8, rearth,rdry, rvap,cp, omega,p00, half, deg2rad, &
                                         fillvalue
  use grist_domain_types,          only: global_domain
  use grist_data_types,            only: scalar_1d_field,exchange_field_list_2d
  use grist_math_module,           only: arcdistll, convert_vector_sph2cart
  use grist_nml_module,            only: nlev, nlevp, nlev_inidata, ntracer, nmif, mif_index
! hpe
  use grist_hpe_constants,         only: eta_face_a, eta_face_b, eta_full_a, eta_full_b, p0
  use grist_hpe_hydro_pgf,         only: calc_hpe_hydro, calc_hpe_get_full_mass
! data structure
  use grist_dycore_vars_module,    only: dycoreVarCellFull, dycoreVarCellFace, dycoreVarEdgeFull,dycoreVarSurface
  use grist_tracer_transport_vars_module, only: tracerVarCellFull
  use grist_dtp_vars_module
! testcase
  use grist_dtp_dcmip2016_sc,      only: supercell_init, supercell_test,supercell_z
  use grist_dtp_dcmip2016_tc,      only: tropical_cyclone_test, dp, rp, exppr, cen_lat, cen_lon,&
                                         ztrop,Ttrop,ptrop, Rd, T0, q0, exponent, gamma
  use grist_dtp_dcmip2016_bw,      only: baroclinic_wave_test, evaluate_pressure_temperature
  use grist_dtp_dcmip2016_terminator, only: initial_value_Terminator
#ifndef SEQ_GRIST
! mpi-comm
  use grist_config_partition,      only: exchange_data_2d_add, exchange_data_2d
  use grist_config_partition,      only: debug_data_2d
#endif
  use grist_mpi
#ifndef SEQ_GRIST
! initial data
  use grist_datam_initial_data_module, only: initialData_uuu_at_pc_full_level, &
                                             initialData_vvv_at_pc_full_level, &
                                             initialData_ttt_at_pc_full_level, &
                                             initialData_qqq_at_pc_full_level, &
                                             !initialData_ppp_at_pc_full_level, &
                                             initialData_hypw_at_pc_face_level,&
                                             initialData_ps_at_pc_surface    , &
                                             initialData_hyai, initialData_hybi, & 
                                             initialData_hyam, initialData_hybm, initialData_plev
  use grist_datam_static_data_module, only: staticData_phis_at_pc_surface
#endif
  use grist_math_module,           only: lininterp
                                       
  implicit none

  private

  public  :: grist_dtp_initial, &
             evaluate_dry_mass_height

!
! Weights and arguments for Gaussian Qudrature
!
  real(r8)  :: zero = 0._r8
#ifndef GAUSS_30
  integer , parameter  :: nGauss = 20
  real(r8), parameter, dimension(nGauss), private :: gaussx = (/-0.0765265211334973, 0.0765265211334973, &
                                                                -0.2277858511416451, 0.2277858511416451, &
                                                                -0.3737060887154195, 0.3737060887154195, &
                                                                -0.5108670019508271, 0.5108670019508271, &
                                                                -0.6360536807265150, 0.6360536807265150, &
                                                                -0.7463319064601508, 0.7463319064601508, &
                                                                -0.8391169718222188, 0.8391169718222188, &
                                                                -0.9122344282513259, 0.9122344282513259, &
                                                                -0.9639719272779138, 0.9639719272779138, &
                                                                -0.9931285991850949, 0.9931285991850949/)

  real(r8), parameter, dimension(nGauss), private :: gaussw = (/0.1527533871307258 , 0.1527533871307258, &
                                                                0.1491729864726037 , 0.1491729864726037, &
                                                                0.1420961093183820 , 0.1420961093183820, &
                                                                0.1316886384491766 , 0.1316886384491766, &
                                                                0.1181945319615184 , 0.1181945319615184, &
                                                                0.1019301198172404 , 0.1019301198172404, &
                                                                0.0832767415767048 , 0.0832767415767048, &
                                                                0.0626720483341091 , 0.0626720483341091, &
                                                                0.0406014298003869 , 0.0406014298003869, &
                                                                0.0176140071391521 , 0.0176140071391521/) 
#endif
#ifdef GAUSS_30
  integer , parameter  :: nGauss = 30
  real(r8), parameter, dimension(nGauss), private :: gaussx = (/-0.0514718425553177, 0.0514718425553177, &
                                                                -0.1538699136085835, 0.1538699136085835, &
                                                                -0.2546369261678899, 0.2546369261678899, &
                                                                -0.3527047255308781, 0.3527047255308781, &
                                                                -0.4470337695380892, 0.4470337695380892, &
                                                                -0.5366241481420199, 0.5366241481420199, &
                                                                -0.6205261829892429, 0.6205261829892429, &
                                                                -0.6978504947933158, 0.6978504947933158, &
                                                                -0.7677774321048262, 0.7677774321048262, &
                                                                -0.8295657623827684, 0.8295657623827684, &
                                                                -0.8825605357920527, 0.8825605357920527, &
                                                                -0.9262000474292743, 0.9262000474292743, &
                                                                -0.9600218649683075, 0.9600218649683075, &
                                                                -0.9836681232797472, 0.9836681232797472, &
                                                                -0.9968934840746495, 0.9968934840746495/)
                                                                                                          
  real(r8), parameter, dimension(nGauss), private :: gaussw = (/0.1028526528935588, 0.1028526528935588, &
                                                                0.1017623897484055, 0.1017623897484055, &
                                                                0.0995934205867953, 0.0995934205867953, &
                                                                0.0963687371746443, 0.0963687371746443, &
                                                                0.0921225222377861, 0.0921225222377861, &
                                                                0.0868997872010830, 0.0868997872010830, &
                                                                0.0807558952294202, 0.0807558952294202, &
                                                                0.0737559747377052, 0.0737559747377052, &
                                                                0.0659742298821805, 0.0659742298821805, &
                                                                0.0574931562176191, 0.0574931562176191, &
                                                                0.0484026728305941, 0.0484026728305941, &
                                                                0.0387991925696271, 0.0387991925696271, &
                                                                0.0287847078833234, 0.0287847078833234, &
                                                                0.0184664683110910, 0.0184664683110910, &
                                                                0.0079681924961666, 0.0079681924961666/)
#endif
                                                                                                         
! temporaily value                                                                                       
   real(r8), parameter  ::  eps = 1e-12_r8                                                               

  CONTAINS

  subroutine grist_dtp_initial(mesh, testcase)
!
! io
!
   type(global_domain), intent(inout) :: mesh
   character*(*)       , intent(in)   :: testcase
!
! local
!
   real(r8)                            :: vector_velocity(3)
   real(r8)                            :: scalar_u
   real(r8)                            :: scalar_v
   real(r8)                            :: eta
   integer(i4)                         :: ilev
   integer(i4)                         :: it, ie, iv
! for RH3D
   real(r8), allocatable               :: stream_function_dual_cell(:)
   integer(i4)                         :: v0,v1, icell1, icell2
   integer(i4)                         :: pert
   integer(i4)                         :: iblock 
   real(r8)                            :: v0v1(3), flag
   real(r8)                            :: zref(nlevp), dz, gr, zface(nlevp), zfull(nlev)
   type(exchange_field_list_2d),pointer:: field_head_2d
   real(r8)                            :: utmp, vtmp, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
#ifdef TC_RANDOM_WIND
   real(r8)                            :: pertval_u, pertval_v
   type(scalar_2d_field)               :: temp_var
#endif
   real(r8)                            :: ztop, ptop, delp, delhp
   real(r8)                            :: pface(nlev_inidata+1), hpface(nlev_inidata+1)
   real(r8)                            :: pfull(nlev_inidata),   hpfull(nlev_inidata)
   real(r8)                            :: dpfull(nlev_inidata),  dhpfull(nlev_inidata)
! for real-case
   real(r8)                            :: scalar_template_1d_nlev_a(nlev)
   real(r8)                            :: scalar_template_1d_nlev_b(nlev)
   real(r8)                            :: scalar_template_1d_nlev_c(nlev)
   real(r8)                            :: scalar_template_1d_nlevp_a(nlevp)
   real(r8)                            :: scalar_template_1d_nlevp_b(nlevp)
   real(r8)                            :: scalar_template_a
   integer(i4)                         :: nlev_eff

   field_head_2d => null()
   iblock  = mpi_rank()

   select case(trim(testcase))
#ifndef SEQ_GRIST
      case('real-WRFDA')  ! real-case based on WRF-DA data, directly on DRY-AIR-MASS

! phis
     dycoreVarSurface%scalar_geopotential_n%f    = staticData_phis_at_pc_surface%f ! m2.s2

     do iv = 1, mesh%nv_full
!
! pressure at ERA model levels
!
        !pfull(1:nlev_inidata)   = initialData_hyam%f(1:nlev_inidata)  + initialData_hybm%f(1:nlev_inidata)  *initialData_ps_at_pc_surface%f(iv)
        !pface(1:nlev_inidata+1) = initialData_hyai%f(1:nlev_inidata+1)+ initialData_hybi%f(1:nlev_inidata+1)*initialData_ps_at_pc_surface%f(iv)
        !dpfull(1:nlev_inidata)  = pface(2:nlev_inidata+1)-pface(1:nlev_inidata)
        !hpface(1)               = pface(1)
!
! evaluate dry mass at ERA model level
!
        do ilev = 1, nlev_inidata+1
           hpface(ilev) = initialData_hypw_at_pc_face_level%f(ilev,iv)
        end do
        hpfull(1:nlev_inidata) = half*(hpface(1:nlev_inidata)+hpface(2:nlev_inidata+1))
!
! evaluate surface dry mass and geopotential
!
        !hpsum = hpface(1)
        !do ilev = 1, nlev_inidata
        !   tmp1  = initialData_hyai%f(ilev)  + initialData_hybi%f(ilev)  *initialData_ps_at_pc_surface%f(iv)
        !   tmp2  = initialData_hyai%f(ilev+1)+ initialData_hybi%f(ilev+1)*initialData_ps_at_pc_surface%f(iv)
        !   delp  = tmp2-tmp1
        !   hpsum = hpsum+delp*(one-initialData_qqq_at_pc_full_level%f(ilev,iv)) ! accumulate delhp
        !end do
        dycoreVarSurface%scalar_hpressure_n%f(iv)             = hpface(nlev_inidata+1)
        dycoreVarCellFace%scalar_hpressure_n%f(1:nlevp,iv)  = eta_face_a(1:nlevp)*p0+eta_face_b(1:nlevp)*hpface(nlev_inidata+1) ! Pa
        dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv)   = eta_full_a(1:nlev) *p0+eta_full_b(1:nlev) *hpface(nlev_inidata+1)
!
! interpolate from ERA model level to GRIST model level based on dry air mass
!
        call lininterp (initialData_uuu_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        dycoreVarCellFull%scalar_U_wind_n%f(1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_vvv_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        dycoreVarCellFull%scalar_V_wind_n%f(1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_ttt_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        dycoreVarCellFull%scalar_temp_n%f(1:nlev,iv)  , dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_qqq_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        ! moist to dry
        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv) = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv)/&
                                                            (one-tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv))
        do ilev = 1, nlev
           tracerVarCellFull%scalar_mif_n%f(ilev,iv) = one/(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,iv)))
        end do

     end do

! find mpressure at model levels
     dycoreVarCellFull%scalar_delhp_n%f(1:nlev,:)  = dycoreVarCellFace%scalar_hpressure_n%f(2:nlevp,:)-&
                                                    dycoreVarCellFace%scalar_hpressure_n%f(1:nlev ,:)

     call calc_hpe_get_full_mass(nlev, mesh%nv_full, 'dtp'                    , & ! in
                                       dycoreVarCellFace%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_delhp_n%f      , & ! in
                                       tracerVarCellFull%scalar_mif_n%f        , & ! in
                                       dycoreVarCellFull%scalar_mpressure_n%f  , & ! out
                                       dycoreVarCellFace%scalar_mpressure_n%f)     ! out

     dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_mpressure_n
     dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_mpressure_n
#ifdef DYCORE_NHD
     dycoreVarCellFull%scalar_delp_n%f(1:nlev,:) = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,:)-&
                                                   dycoreVarCellFace%scalar_pressure_n%f(1:nlev,:)
#endif
!
! evaluate geop based on the moist hydrostatic equation
!
     do iv = 1, mesh%nv_full

        scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_temp_n%f(1:nlev,iv)*&
                                     (one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv))))      ! Tv
        scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f(1:nlevp,iv)   ! pressure face
        !scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delp_n%f(1:nlev,iv)       ! delp full
        scalar_template_1d_nlev_b  = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,iv)-&
                                     dycoreVarCellFace%scalar_pressure_n%f(1:nlev,iv)     ! delp full

        scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)         ! geopotential at surface

        call calc_hpe_hydro(nlev, nlevp, &
                            scalar_template_1d_nlev_a  ,& ! Tv or T
                            scalar_template_1d_nlevp_a ,& ! p or hp
                            scalar_template_1d_nlev_b  ,& ! delp or delhp
                            scalar_template_a          ,& ! phis
                            scalar_template_1d_nlevp_b ,& ! face geop
                            scalar_template_1d_nlev_c )   ! full geop

        dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,iv) = scalar_template_1d_nlevp_b
        dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,iv)  = scalar_template_1d_nlev_c

     end do

#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_www_n%f           = 0._r8
     dycoreVarCellFace%scalar_phi_n%f(1:nlevp,:)= dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,:)
     dycoreVarCellFace%scalar_phi_n%f(nlevp,:)  = dycoreVarSurface%scalar_geopotential_n%f ! surface is given in ini condition
#endif

   case('real-ERAIP','real-GFS','real-NuGFS')   ! real-case based on any pressure level,
                                                ! treat them as full level, bounded by zero and ps
                                                ! can be used for any pressure-based initial conditions
! phis
     dycoreVarSurface%scalar_geopotential_n%f    = staticData_phis_at_pc_surface%f ! m2.s2

     do iv = 1, mesh%nv_full
!
! nominal pressure at ERA pressure levels
!
        pfull(1:nlev_inidata)   = initialData_plev%f  ! Pa
!
! check effective full level number
!
        nlev_eff = 99999 
        do ilev = 1, nlev_inidata
           if(pfull(ilev).ge.initialData_ps_at_pc_surface%f(iv))then
              nlev_eff = ilev
              exit
           end if
        end do
        if(nlev_eff.eq.99999) nlev_eff = nlev_inidata

        pface(2:nlev_eff)       = (pfull(1:nlev_eff-1)+pfull(2:nlev_eff))*half
        pface(1)                = zero
        pface(nlev_eff+1)       = initialData_ps_at_pc_surface%f(iv) ! if ps<pface(nlev_eff), let it be, good assumption! no problem till now...
        dpfull(1:nlev_eff)      = pface(2:nlev_eff+1)-pface(1:nlev_eff)
        hpface(1)               = pface(1)
!
! evaluate dry mass at ERA pressure level and surface
!
        do ilev = 2, nlev_eff+1
           hpface(ilev) = hpface(ilev-1)+dpfull(ilev-1)*(one-initialData_qqq_at_pc_full_level%f(ilev-1,iv))
        end do
        hpfull(1:nlev_eff) = half*(hpface(1:nlev_eff)+hpface(2:nlev_eff+1))

        dycoreVarSurface%scalar_hpressure_n%f(iv)             = hpface(nlev_eff+1) 
        dycoreVarCellFace%scalar_hpressure_n%f(1:nlevp,iv)  = eta_face_a(1:nlevp)*p0+eta_face_b(1:nlevp)*hpface(nlev_eff+1) ! Pa
        dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv)   = eta_full_a(1:nlev) *p0+eta_full_b(1:nlev) *hpface(nlev_eff+1)
!
! interpolate from ERA pressure level to GRIST model level based on dry air mass
!
        call lininterp (initialData_uuu_at_pc_full_level%f(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
                        dycoreVarCellFull%scalar_U_wind_n%f(1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_vvv_at_pc_full_level%f(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
                        dycoreVarCellFull%scalar_V_wind_n%f(1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_ttt_at_pc_full_level%f(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
                        dycoreVarCellFull%scalar_temp_n%f(1:nlev,iv)  , dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_qqq_at_pc_full_level%f(1:nlev_eff,iv), hpfull(1:nlev_eff), 1, nlev_eff, &
                        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        ! moist to dry
        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv) = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv)/&
                                                            (one-tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv))
        do ilev = 1, nlev
           tracerVarCellFull%scalar_mif_n%f(ilev,iv) = one/(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,iv)))
        end do

     end do

! find mpressure at model levels
     dycoreVarCellFull%scalar_delhp_n%f(1:nlev,:)  = dycoreVarCellFace%scalar_hpressure_n%f(2:nlevp,:)-&
                                                    dycoreVarCellFace%scalar_hpressure_n%f(1:nlev ,:)

     call calc_hpe_get_full_mass(nlev, mesh%nv_full, 'dtp'                    , & ! in
                                       dycoreVarCellFace%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_delhp_n%f      , & ! in
                                       tracerVarCellFull%scalar_mif_n%f        , & ! in
                                       dycoreVarCellFull%scalar_mpressure_n%f  , & ! out
                                       dycoreVarCellFace%scalar_mpressure_n%f)     ! out

     dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_mpressure_n
     dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_mpressure_n

#ifdef DYCORE_NHD
     dycoreVarCellFull%scalar_delp_n%f(1:nlev,:) = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,:)-&
                                                   dycoreVarCellFace%scalar_pressure_n%f(1:nlev,:)
#endif
!
! evaluate geop based on the moist hydrostatic equation
!
     do iv = 1, mesh%nv_full

        scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_temp_n%f(1:nlev,iv)*&
                                     (one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv))))      ! Tv
        scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f(1:nlevp,iv)  ! pressure face
        !scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delp_n%f(1:nlev,iv)       ! delp full
        scalar_template_1d_nlev_b  = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,iv)-&
                                     dycoreVarCellFace%scalar_pressure_n%f(1:nlev,iv)
        scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)         ! geopotential at surface

        call calc_hpe_hydro(nlev, nlevp,scalar_template_1d_nlev_a  ,& ! Tv or T
                            scalar_template_1d_nlevp_a ,& ! p or hp
                            scalar_template_1d_nlev_b  ,& ! delp or delhp
                            scalar_template_a          ,& ! phis
                            scalar_template_1d_nlevp_b ,& ! face geop
                            scalar_template_1d_nlev_c )   ! full geop

        dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,iv) = scalar_template_1d_nlevp_b
        dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,iv)  = scalar_template_1d_nlev_c

     end do

#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_www_n%f           = zero 
     dycoreVarCellFace%scalar_phi_n%f(1:nlevp,:)= dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,:)
     dycoreVarCellFace%scalar_phi_n%f(nlevp,:)  = dycoreVarSurface%scalar_geopotential_n%f ! surface is given in ini condition
#endif

   case('real-ERAIM')  ! real-case based on ERAI's model level

! phis
     dycoreVarSurface%scalar_geopotential_n%f    = staticData_phis_at_pc_surface%f ! m2.s2

     do iv = 1, mesh%nv_full
!
! pressure at ERA model levels
!
        pfull(1:nlev_inidata)   = initialData_hyam%f(1:nlev_inidata)  + initialData_hybm%f(1:nlev_inidata)  *initialData_ps_at_pc_surface%f(iv)
        pface(1:nlev_inidata+1) = initialData_hyai%f(1:nlev_inidata+1)+ initialData_hybi%f(1:nlev_inidata+1)*initialData_ps_at_pc_surface%f(iv)
        dpfull(1:nlev_inidata)  = pface(2:nlev_inidata+1)-pface(1:nlev_inidata)
        hpface(1)               = pface(1)
!
! evaluate dry mass at ERA model level
!
        do ilev = 2, nlev_inidata+1
           hpface(ilev) = hpface(ilev-1)+dpfull(ilev-1)*(one-initialData_qqq_at_pc_full_level%f(ilev-1,iv))
        end do
        hpfull(1:nlev_inidata) = half*(hpface(1:nlev_inidata)+hpface(2:nlev_inidata+1))
!
! evaluate surface dry mass and geopotential
!
        !hpsum = hpface(1)
        !do ilev = 1, nlev_inidata
        !   tmp1  = initialData_hyai%f(ilev)  + initialData_hybi%f(ilev)  *initialData_ps_at_pc_surface%f(iv)
        !   tmp2  = initialData_hyai%f(ilev+1)+ initialData_hybi%f(ilev+1)*initialData_ps_at_pc_surface%f(iv)
        !   delp  = tmp2-tmp1
        !   hpsum = hpsum+delp*(one-initialData_qqq_at_pc_full_level%f(ilev,iv)) ! accumulate delhp
        !end do
        dycoreVarSurface%scalar_hpressure_n%f(iv)             = hpface(nlev_inidata+1) 
        dycoreVarCellFace%scalar_hpressure_n%f(1:nlevp,iv)  = eta_face_a(1:nlevp)*p0+eta_face_b(1:nlevp)*hpface(nlev_inidata+1) ! Pa
        dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv)   = eta_full_a(1:nlev) *p0+eta_full_b(1:nlev) *hpface(nlev_inidata+1)
!
! interpolate from ERA model level to GRIST model level based on dry air mass
!
        call lininterp (initialData_uuu_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        dycoreVarCellFull%scalar_U_wind_n%f(1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_vvv_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        dycoreVarCellFull%scalar_V_wind_n%f(1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_ttt_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        dycoreVarCellFull%scalar_temp_n%f(1:nlev,iv)  , dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        call lininterp (initialData_qqq_at_pc_full_level%f(1:nlev_inidata,iv), hpfull, 1, nlev_inidata, &
                        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv), dycoreVarCellFull%scalar_hpressure_n%f(1:nlev,iv), nlev)
        ! moist to dry
        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv) = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv)/&
                                                            (one-tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv))
        do ilev = 1, nlev
           tracerVarCellFull%scalar_mif_n%f(ilev,iv) = one/(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,iv)))
        end do

     end do

! find mpressure at model levels
     dycoreVarCellFull%scalar_delhp_n%f(1:nlev,:)  = dycoreVarCellFace%scalar_hpressure_n%f(2:nlevp,:)-&
                                                    dycoreVarCellFace%scalar_hpressure_n%f(1:nlev ,:)

     call calc_hpe_get_full_mass(nlev, mesh%nv_full, 'dtp'                    , & ! in
                                       dycoreVarCellFace%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_delhp_n%f      , & ! in
                                       tracerVarCellFull%scalar_mif_n%f        , & ! in
                                       dycoreVarCellFull%scalar_mpressure_n%f  , & ! out
                                       dycoreVarCellFace%scalar_mpressure_n%f)     ! out

     dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_mpressure_n
     dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_mpressure_n
#ifdef DYCORE_NHD
     dycoreVarCellFull%scalar_delp_n%f(1:nlev,:) = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,:)-&
                                                  dycoreVarCellFace%scalar_pressure_n%f(1:nlev,:)
#endif
!
! evaluate geop based on the moist hydrostatic equation
!
     do iv = 1, mesh%nv_full

        scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_temp_n%f(1:nlev,iv)*&
                                     (one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv))))      ! Tv
        scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f(1:nlevp,iv)  ! pressure face
        !scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delp_n%f(1:nlev,iv)       ! delp full
        scalar_template_1d_nlev_b  = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,iv)-&
                                     dycoreVarCellFace%scalar_pressure_n%f(1:nlev,iv)
        scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)         ! geopotential at surface

        call calc_hpe_hydro(nlev,nlevp,scalar_template_1d_nlev_a  ,& ! Tv or T
                            scalar_template_1d_nlevp_a ,& ! p or hp
                            scalar_template_1d_nlev_b  ,& ! delp or delhp
                            scalar_template_a          ,& ! phis
                            scalar_template_1d_nlevp_b ,& ! face geop
                            scalar_template_1d_nlev_c )   ! full geop

        dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,iv) = scalar_template_1d_nlevp_b
        dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,iv)  = scalar_template_1d_nlev_c

     end do

#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_www_n%f             = 0._r8
     dycoreVarCellFace%scalar_phi_n%f(1:nlevp,:)  = dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,:)
     dycoreVarCellFace%scalar_phi_n%f(nlevp,:)    = dycoreVarSurface%scalar_geopotential_n%f ! surface is given in ini condition
#endif

#endif
   case('DCMIP2016-BW')
      tracerVarCellFull%scalar_tracer_mxrt_n%f  = 0._r8
!
! evaluate dry mass and geometric height at each level
!
      call evaluate_dry_mass_height(mesh, mesh%nv_full, nlev, 'DCMIP2016-BW', &
                                    dycoreVarSurface%scalar_hpressure_n%f      , &
                                    dycoreVarCellFace%scalar_hpressure_n%f   , &
                                    dycoreVarCellFull%scalar_hpressure_n%f   , &
                                    dycoreVarCellFace%scalar_geopotential_n%f, &
                                    dycoreVarCellFull%scalar_geopotential_n%f)
!
! final evaluation based on geop full & face
!
      pert = 0 ! expotional
      do iv = 1, mesh%nv_full
        do ilev = 1, nlevp
           call baroclinic_wave_test(0, 1, pert, one, real(mesh%vtx(iv)%lon,r8), real(mesh%vtx(iv)%lat,r8) ,&
                                     dycoreVarCellFace%scalar_pressure_n%f(ilev,iv)       ,&
                                     dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)   ,&
                                     1                                                   ,&
                                     tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
        end do
        do ilev = 1, nlev
           call baroclinic_wave_test(0, 1, pert, one, real(mesh%vtx(iv)%lon,r8), real(mesh%vtx(iv)%lat,r8) ,&
                                     dycoreVarCellFull%scalar_pressure_n%f(ilev,iv)       ,&
                                     dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv)   ,&
                                     1                                                   ,&
                                     dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv),&
                                     dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv),&
                                     dycoreVarCellFull%scalar_temp_n%f(ilev,iv)  ,&
                                     scalar_thetav_at_pc_full_level_n%f(ilev,iv),& ! not used
                                     dycoreVarSurface%scalar_geopotential_n%f(iv)  ,&
                                     dycoreVarSurface%scalar_pressure_n%f(iv)      ,& ! not used 
                                     scalar_rho_at_pc_full_level_n%f(ilev,iv)   ,& ! not used
                                     tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
! moist to dry
        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv) = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv)/&
                                                            (one-tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))

! re-evaluate dry mixing ratio based on density scaling
!            tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv) =((dycoreVarCellFace%scalar_pressure_n%f(ilev+1,iv)-&
!                                                                  dycoreVarCellFace%scalar_pressure_n%f(ilev,iv))/&
!                                                                 (dycoreVarCellFace%scalar_hpressure_n%f(ilev+1,iv)-&
!                                                                  dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv)))-one
! evaluate mif here for initial mpressure and energy
#ifndef MITC
           tracerVarCellFull%scalar_mif_n%f(ilev,iv) = one/(one+sum(tracerVarCellFull%scalar_tracer_mxrt_n%f(mif_index(1:nmif),ilev,iv)))
#endif
 
        end do
      end do
      dycoreVarCellFace%scalar_geopotential_n%f = gravity*dycoreVarCellFace%scalar_geopotential_n%f
      dycoreVarCellFull%scalar_geopotential_n%f = gravity*dycoreVarCellFull%scalar_geopotential_n%f
#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_phi_n%f          = dycoreVarCellFace%scalar_geopotential_n%f
      dycoreVarCellFace%scalar_www_n%f          = zero
#endif

#ifdef BW_TEST_TE
! find mpressure at model levels
     dycoreVarCellFull%scalar_delhp_n%f(1:nlev,:)  = dycoreVarCellFace%scalar_hpressure_n%f(2:nlevp,:)-&
                                                    dycoreVarCellFace%scalar_hpressure_n%f(1:nlev ,:)

     call calc_hpe_get_full_mass(nlev, mesh%nv_full, 'dtp'                    , & ! in
                                       dycoreVarCellFace%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_hpressure_n%f  , & ! in
                                       dycoreVarCellFull%scalar_delhp_n%f      , & ! in
                                       tracerVarCellFull%scalar_mif_n%f        , & ! in
                                       dycoreVarCellFull%scalar_mpressure_n%f  , & ! out
                                       dycoreVarCellFace%scalar_mpressure_n%f)     ! out

     dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_mpressure_n
     dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_mpressure_n
     dycoreVarCellFull%scalar_delp_n%f(1:nlev,:) = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,:)-&
                                                  dycoreVarCellFace%scalar_pressure_n%f(1:nlev,:)
!
! evaluate geop based on the moist hydrostatic equation
!
     do iv = 1, mesh%nv_full

        scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_temp_n%f(1:nlev,iv)*&
                                     (one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,1:nlev,iv))))      ! Tv
        scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f(1:nlevp,iv)  ! pressure face
        scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delp_n%f(1:nlev,iv)       ! delp full
        scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)         ! geopotential at surface

        call calc_hpe_hydro(nlev, nlevp,scalar_template_1d_nlev_a  ,& ! Tv or T
                            scalar_template_1d_nlevp_a ,& ! p or hp
                            scalar_template_1d_nlev_b  ,& ! delp or delhp
                            scalar_template_a          ,& ! phis
                            scalar_template_1d_nlevp_b ,& ! face geop
                            scalar_template_1d_nlev_c )   ! full geop

        dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,iv) = scalar_template_1d_nlevp_b
        dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,iv)  = scalar_template_1d_nlev_c

#ifdef DYCORE_NHD
        dycoreVarCellFace%scalar_phi_n%f(1:nlevp,iv)          = dycoreVarCellFace%scalar_geopotential_n%f(1:nlevp,iv)
#endif

     end do
#endif

! set toy terminator tracer
      if(ntracer.ge.2)then
      do iv = 1, mesh%nv_full
         if(mesh%vtx(iv)%lon.ge.0)then
            tmp1 = mesh%vtx(iv)%lon/(pi/180._r8)
         else
            tmp1 = (mesh%vtx(iv)%lon+2._r8*pi)/(pi/180._r8)
         end if
         do ilev = 1, nlev
            call  initial_value_Terminator(real(mesh%vtx(iv)%lat,r8)/(pi/180._r8), tmp1, &
                                           tracerVarCellFull%scalar_tracer_mxrt_n%f(2,ilev,iv), &
                                           tracerVarCellFull%scalar_tracer_mxrt_n%f(3,ilev,iv))
         end do
      end do
! or coupled with kessler
      tracerVarCellFull%scalar_tracer_mxrt_n%f(2:ntracer,:,:) = zero
      end if

   case('DCMIP2016-TC') ! FML

      tracerVarCellFull%scalar_tracer_mxrt_n%f  = 0._r8
!
! evaluate dry mass and geometric height at each level
!
      call evaluate_dry_mass_height(mesh, mesh%nv_full, nlev, 'DCMIP2016-TC', &
                                    dycoreVarSurface%scalar_hpressure_n%f      , &
                                    dycoreVarCellFace%scalar_hpressure_n%f   , &
                                    dycoreVarCellFull%scalar_hpressure_n%f   , &
                                    dycoreVarCellFace%scalar_geopotential_n%f, &
                                    dycoreVarCellFull%scalar_geopotential_n%f)
!
! final evaluation based on geop full & face
!
      do iv = 1, mesh%nv_full
         do ilev = 1, nlevp
            call tropical_cyclone_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8) , &
                                       dycoreVarCellFace%scalar_pressure_n%f(ilev,iv)    , &
                                       dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv), 1, &
                                       tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
        end do
        do ilev = 1, nlev
           call tropical_cyclone_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)  ,&
                                      dycoreVarCellFull%scalar_pressure_n%f(ilev,iv) ,&
                                      dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv), 1, &
                                      dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv)   ,&
                                      dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv)   ,&
                                      dycoreVarCellFull%scalar_temp_n%f(ilev,iv)     ,&
                                      scalar_thetav_at_pc_full_level_n%f(ilev,iv)   ,&  ! not used
                                      dycoreVarSurface%scalar_geopotential_n%f(iv)     ,&
                                      dycoreVarSurface%scalar_pressure_n%f(iv)         ,&  ! not used
                                      scalar_rho_at_pc_full_level_n%f(ilev,iv)      ,&  ! not used
                                      tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
! moist to dry
        tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv) = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv)/&
                                                            (one-tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
 
        end do
      end do

      dycoreVarCellFace%scalar_geopotential_n%f = gravity*dycoreVarCellFace%scalar_geopotential_n%f
      dycoreVarCellFull%scalar_geopotential_n%f = gravity*dycoreVarCellFull%scalar_geopotential_n%f
      dycoreVarSurface%scalar_geopotential_n%f  = zero
#ifdef DYCORE_NHD
      dycoreVarCellFace%scalar_phi_n%f          = dycoreVarCellFace%scalar_geopotential_n%f
      dycoreVarCellFace%scalar_www_n%f          = zero
#endif

   case('DCMIP2016-SCXX') ! hybrid 1, same to SC now

      pert = 0
      call supercell_init
      if(iblock.eq.0) print*, "supercell_init done"

      ptop = eta_face_a(1)*p0

      do iv = 1, mesh%nv_full
! ztop assume no moisture
         call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8),ptop,ztop,zero,zero,&
                             0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,pert)
! surface pressure
         call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8),dycoreVarSurface%scalar_pressure_n%f(iv),zero,zero,zero,&
                             1,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,pert)
! from moist to dry mixing ratio
!         tmp8 = tmp8/(one-tmp8)
! surface hpressure
         dycoreVarSurface%scalar_hpressure_n%f(iv) = dycoreVarSurface%scalar_pressure_n%f(iv)/(one+rvap/rdry*tmp8)

         do ilev = 1, nlevp
            dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do
         do ilev = 1, nlev
            dycoreVarCellFull%scalar_hpressure_n%f(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do
! evaluate face level height
         do ilev = 1, nlevp
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)                ,&
                                 dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv)   ,&
                                 dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv),& 
                                 ptop, ztop, 2,&
                                 tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7, tmp8, pert)
         end do
         dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,iv) =(dycoreVarCellFace%scalar_geopotential_n%f(1:nlev ,iv)+&
                                                               dycoreVarCellFace%scalar_geopotential_n%f(2:nlevp,iv))*half
     end do
     if(iblock.eq.0) print*, "height at full and face levels are done"
!
! evaluate final state based on height at full and face levels with pert=1
!
     pert = 1
     do iv = 1, mesh%nv_full
        do ilev = 1, nlevp
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8) , &
                                 dycoreVarCellFace%scalar_pressure_n%f(ilev,iv)    , &
                                 dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv), zero,zero, &
                                 1, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, pert)
        end do
        do ilev = 1, nlev
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8) , &
                                 dycoreVarCellFull%scalar_pressure_n%f(ilev,iv), &
                                 dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv), zero, zero, 1, &
                                 dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv)  , &
                                 dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv)  , &
                                 dycoreVarCellFull%scalar_temp_n%f(ilev,iv)    , &
                                 scalar_thetav_at_pc_full_level_n%f(ilev,iv)  , &   ! not used
                                 dycoreVarCellFull%scalar_potential_temp_iniupt%f(ilev,iv), & ! thetap
                                 dycoreVarSurface%scalar_pressure_n%f(iv)        , &   ! not used
                                 scalar_rho_at_pc_full_level_n%f(ilev,iv)     , &   ! not used
                                 tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv),pert)
! moist to dry mixing ratio
!              tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv) = tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv)/(one-&
!                                                                   tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
         end do
     end do

     dycoreVarCellFace%scalar_geopotential_n%f = gravity*dycoreVarCellFace%scalar_geopotential_n%f
     dycoreVarCellFull%scalar_geopotential_n%f = gravity*dycoreVarCellFull%scalar_geopotential_n%f
     dycoreVarSurface%scalar_geopotential_n%f  = zero
     tracerVarCellFull%scalar_tracer_mxrt_n%f(2:ntracer,:,:) = zero
#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_phi_n%f          = dycoreVarCellFace%scalar_geopotential_n%f
     dycoreVarCellFace%scalar_www_n%f          = zero
#endif

     if(iblock.eq.0) print*,"DCMIP2016-SC initial done"

   case('DCMIP2016-SC') ! as TC

      call supercell_init()
!
! evaluate dry mass and geometric height at each level
!
      call evaluate_dry_mass_height(mesh, mesh%nv_full, nlev, 'DCMIP2016-SC', &
                                    dycoreVarSurface%scalar_hpressure_n%f      , &
                                    dycoreVarCellFace%scalar_hpressure_n%f   , &
                                    dycoreVarCellFull%scalar_hpressure_n%f   , &
                                    dycoreVarCellFace%scalar_geopotential_n%f, &
                                    dycoreVarCellFull%scalar_geopotential_n%f)
      pert = 1
      do iv = 1, mesh%nv_full
         do ilev = 1, nlevp
              call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), &
                                  dycoreVarCellFace%scalar_pressure_n%f(ilev,iv)    , &
                                  dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv), zero, zero, 1, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, pert)
         end do
         do ilev = 1, nlev
              call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), &
                                  dycoreVarCellFull%scalar_pressure_n%f(ilev,iv), &
                                  dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv), zero, zero, 1, &
                                  dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv)  , &
                                  dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv)  , &
                                  dycoreVarCellFull%scalar_temp_n%f(ilev,iv)    , &
                                  scalar_thetav_at_pc_full_level_n%f(ilev,iv)  , &   ! not used
                                  dycoreVarCellFull%scalar_potential_temp_iniupt%f(ilev,iv), & ! thetap
                                  dycoreVarSurface%scalar_pressure_n%f(iv)        , &   ! not used
                                  scalar_rho_at_pc_full_level_n%f(ilev,iv)     , &   ! not used
                                  tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv),pert)
         end do
     end do

     dycoreVarCellFace%scalar_geopotential_n%f = gravity*dycoreVarCellFace%scalar_geopotential_n%f
     dycoreVarCellFull%scalar_geopotential_n%f = gravity*dycoreVarCellFull%scalar_geopotential_n%f
     dycoreVarSurface%scalar_geopotential_n%f  = zero
     tracerVarCellFull%scalar_tracer_mxrt_n%f(2:ntracer,:,:) = zero

#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_www_n%f          = zero
     dycoreVarCellFace%scalar_phi_n%f          = dycoreVarCellFace%scalar_geopotential_n%f
#endif

     if(iblock.eq.0) print*,"DCMIP2016-SC initial done"

   case('DCMIP2016-SC-simple') !round2 default

      pert = 0
      call supercell_init
      if(iblock.eq.0) print*, "supercell_init done"
!
! init surface pressure and hpressure
!
      do iv = 1, mesh%nv_full
         call supercell_z(real(mesh%vtx(iv)%lon,r8), real(mesh%vtx(iv)%lat,r8), 0._r8, dycoreVarSurface%scalar_pressure_n%f(iv), tmp1, tmp4, tmp2, tmp3,pert)
         dycoreVarSurface%scalar_hpressure_n%f(iv) = dycoreVarSurface%scalar_pressure_n%f(iv)/(one+rvap/rdry*tmp3)
      end do
      if(iblock.eq.0) print*, "supercell_z at surface done"
!
! evaluate full and face level hpressure and pressure,
!
     do iv = 1, mesh%nv_full
        do ilev = 1, nlevp
           dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
            dycoreVarCellFace%scalar_pressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)* dycoreVarSurface%scalar_pressure_n%f(iv)
        end do
        do ilev = 1, nlev
           dycoreVarCellFull%scalar_hpressure_n%f(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
            dycoreVarCellFull%scalar_pressure_n%f(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)* dycoreVarSurface%scalar_pressure_n%f(iv)
        end do
     end do
     if(iblock.eq.0) print*, "hpressure and pressure done"
!
! init face and full level height, unit: meter
!
     do iv = 1, mesh%nv_full
        do ilev = 1, nlevp
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)   ,&
                                 dycoreVarCellFace%scalar_pressure_n%f(ilev,iv),&
                                 dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv), zero, zero, 0,&
                                 tmp1,&
                                 tmp2,&
                                 tmp3,&
                                 tmp4,&
                                 tmp5,&
                                 tmp6,&
                                 tmp7, tmp8, pert)
        end do
        dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,iv) =(dycoreVarCellFace%scalar_geopotential_n%f(1:nlev,iv)+&
                                                              dycoreVarCellFace%scalar_geopotential_n%f(2:nlevp,iv))*half
     end do
     if(iblock.eq.0) print*, "height at face evaluation done"
!
! evaluate all initials defined at full level based on height
!
     pert = 1
     do iv = 1, mesh%nv_full
        do ilev = 1, nlev
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)          ,&
                                 tmp2, &
                                 dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv), zero, zero, 1,&
                                 dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_temp_n%f(ilev,iv)  ,&
                                 scalar_thetav_at_pc_full_level_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_potential_temp_iniupt%f(ilev,iv), & ! thetap
                                 tmp1                                       ,& ! ps repeated
                                 scalar_rho_at_pc_full_level_n%f(ilev,iv)   ,& ! this is full density
                                 tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv),pert)
             tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv) = max(0._r8,tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
        end do
     end do

     dycoreVarCellFace%scalar_geopotential_n%f = gravity*dycoreVarCellFace%scalar_geopotential_n%f
     dycoreVarCellFull%scalar_geopotential_n%f = gravity*dycoreVarCellFull%scalar_geopotential_n%f
     dycoreVarSurface%scalar_geopotential_n%f    = zero
     tracerVarCellFull%scalar_tracer_mxrt_n%f(2:ntracer,:,:) = zero
#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_www_n%f          = zero
     dycoreVarCellFace%scalar_phi_n%f          = dycoreVarCellFace%scalar_geopotential_n%f
#endif

     if(iblock.eq.0) print*,"DCMIP2016-SC-simple initial done"

   case('DCMIP2016-SC-A') !round2 default

      pert = 0
      call supercell_init
      if(iblock.eq.0) print*, "supercell_init done"
!
! init surface pressure and hpressure
!
      do iv = 1, mesh%nv_full
         call supercell_z(real(mesh%vtx(iv)%lon,r8), real(mesh%vtx(iv)%lat,r8), 0._r8, dycoreVarSurface%scalar_pressure_n%f(iv), tmp1, tmp4, tmp2, tmp3,pert)
         dycoreVarSurface%scalar_hpressure_n%f(iv) = dycoreVarSurface%scalar_pressure_n%f(iv)/(one+rvap/rdry*tmp3)
      end do
      if(iblock.eq.0) print*, "supercell_z at surface done"
!
! evaluate full and face level hpressure and pressure,
!
     do iv = 1, mesh%nv_full
        do ilev = 1, nlevp
           dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
            dycoreVarCellFace%scalar_pressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)* dycoreVarSurface%scalar_pressure_n%f(iv)
        end do
        do ilev = 1, nlev
           dycoreVarCellFull%scalar_hpressure_n%f(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
            dycoreVarCellFull%scalar_pressure_n%f(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)* dycoreVarSurface%scalar_pressure_n%f(iv)
        end do
     end do
     if(iblock.eq.0) print*, "hpressure and pressure done"
!
! init face and full level height, unit: meter
!
     do iv = 1, mesh%nv_full
        do ilev = 1, nlevp
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)             ,&
                                 dycoreVarCellFace%scalar_pressure_n%f(ilev,iv),&
                                 dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv), zero, zero, 0,&
                                 tmp1,&
                                 tmp2,&
                                 tmp3,&
                                 tmp4,&
                                 tmp5,&
                                 tmp6,&
                                 tmp7, tmp8, pert)
        end do
        dycoreVarCellFull%scalar_geopotential_n%f(1:nlev,iv) =(dycoreVarCellFace%scalar_geopotential_n%f(1:nlev,iv)+&
                                                              dycoreVarCellFace%scalar_geopotential_n%f(2:nlevp,iv))*half
     end do
     if(iblock.eq.0) print*, "height at face evaluation done"
!
! evaluate all initials defined at full level based on height
!
     pert = 1
     do iv = 1, mesh%nv_full
        do ilev = 1, nlev
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)          ,&
                                 tmp2, &
                                 dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv), zero, zero, 1,&
                                 dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_temp_n%f(ilev,iv)  ,&
                                 scalar_thetav_at_pc_full_level_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_potential_temp_iniupt%f(ilev,iv), & ! thetap
                                 tmp1                                       ,& ! ps repeated
                                 scalar_rho_at_pc_full_level_n%f(ilev,iv)   ,& ! this is full density
                                 tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv),pert)
             tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv) = max(0._r8,tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
        end do
     end do
!
! iterative step to bring perturbed pt to hydrostatic balance
! given hpressure, balance between temperature and geop under 
! dry constraint
!
     call grist_dtp_iterative_step(mesh, nlev,'DCMIP2016-SC'              ,&
                                   dycoreVarCellFace%scalar_hpressure_n    ,&
                                   dycoreVarCellFull%scalar_temp_n         ,&
                                   dycoreVarCellFull%scalar_geopotential_n ,&
                                   dycoreVarCellFace%scalar_geopotential_n)

!
! final evaluation based on geop full & face
!
     do iv = 1, mesh%nv_full
        do ilev = 1, nlevp
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)          ,&
                                 tmp1, &
                                 dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv), zero, zero, 1,&
                                 tmp2,&
                                 tmp3,&
                                 tmp4,&
                                 tmp5,&
                                 tmp6,& ! thetap
                                 tmp7,& ! ps repeated
                                 tmp8,& ! this is full density
                                 tmp9,pert)
             dycoreVarCellFace%scalar_pressure_n%f(ilev,iv) = dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv)*&
                                                             (one+rvap/rdry*tmp9)
        end do
        do ilev = 1, nlev
             call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)          ,&
                                 tmp2, &
                                 dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv), zero, zero, 1,&
                                 dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_temp_n%f(ilev,iv)  ,&
                                 scalar_thetav_at_pc_full_level_n%f(ilev,iv),&
                                 dycoreVarCellFull%scalar_potential_temp_iniupt%f(ilev,iv), & ! thetap
                                 tmp1                                       ,& ! ps repeated
                                 scalar_rho_at_pc_full_level_n%f(ilev,iv)   ,& ! this is full density
                                 tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv),pert)

             dycoreVarCellFull%scalar_pressure_n%f(ilev,iv) = dycoreVarCellFull%scalar_hpressure_n%f(ilev,iv)*&
                                                             (one+rvap/rdry*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv))
         end do
     end do

     dycoreVarCellFace%scalar_geopotential_n%f = gravity*dycoreVarCellFace%scalar_geopotential_n%f
     dycoreVarCellFull%scalar_geopotential_n%f = gravity*dycoreVarCellFull%scalar_geopotential_n%f
     dycoreVarSurface%scalar_geopotential_n%f    = zero
     tracerVarCellFull%scalar_tracer_mxrt_n%f(2:ntracer,:,:) = zero
#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_www_n%f          = zero
     dycoreVarCellFace%scalar_phi_n%f          = dycoreVarCellFace%scalar_geopotential_n%f
#endif

     if(iblock.eq.0) print*,"DCMIP2016-SC-A initial done"

   case('DCMIP2016-SC-B') ! hybrid, surface hps is evaluated as in SC-A, other as in SC

      pert = 1
      call supercell_init
      if(iblock.eq.0) print*, "supercell_init done"
!
! init surface pressure and hpressure
!
      do iv = 1, mesh%nv_full
         call supercell_z(real(mesh%vtx(iv)%lon,r8), real(mesh%vtx(iv)%lat,r8), 0._r8, dycoreVarSurface%scalar_pressure_n%f(iv), tmp1, tmp4, tmp2, tmp3,pert)
         dycoreVarSurface%scalar_hpressure_n%f(iv) = dycoreVarSurface%scalar_pressure_n%f(iv)/(one+rvap/rdry*tmp3)
      end do
      if(iblock.eq.0) print*, "supercell_z at surface done"
!
! evaluate full and face level height, hpressure
!
      ptop = eta_face_a(1)*p0

      do iv = 1, mesh%nv_full
         ztop  = get_z_from_pressure(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), ptop, ptop, 0.0_r8, 'DCMIP2016-SC', .false.)

         do ilev = 1, nlevp
            dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do
         do ilev = 1, nlev
            dycoreVarCellFull%scalar_hpressure_n%f(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
         end do

         do ilev = 1, nlev+1
            dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)  = get_z_from_pressure(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), &
                                                   dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv), ptop, ztop, 'DCMIP2016-SC', .true.)
         end do
         do ilev = 1, nlev
            dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv)  = 0.5_r8*(dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv)+dycoreVarCellFace%scalar_geopotential_n%f(ilev+1,iv))
         end do

      end do

      if(iblock.eq.0) print*, "hpressure and height done"

      do iv = 1, mesh%nv_full
         do ilev = 1, nlevp
              call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)            , &
                                  dycoreVarCellFace%scalar_pressure_n%f(ilev,iv)    , &
                                  dycoreVarCellFace%scalar_geopotential_n%f(ilev,iv), zero, zero, 1, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, pert)
         end do
         do ilev = 1, nlev
              call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)            , &
                                  dycoreVarCellFull%scalar_pressure_n%f(ilev,iv), &
                                  dycoreVarCellFull%scalar_geopotential_n%f(ilev,iv), zero, zero, 1, &
                                  dycoreVarCellFull%scalar_U_wind_n%f(ilev,iv)  , &
                                  dycoreVarCellFull%scalar_V_wind_n%f(ilev,iv)  , &
                                  dycoreVarCellFull%scalar_temp_n%f(ilev,iv)    , &
                                  scalar_thetav_at_pc_full_level_n%f(ilev,iv)  , &   ! not used
                                  dycoreVarCellFull%scalar_potential_temp_iniupt%f(ilev,iv), & ! thetap
                                  dycoreVarSurface%scalar_pressure_n%f(iv)        , &   ! not used
                                  scalar_rho_at_pc_full_level_n%f(ilev,iv)     , &   ! not used
                                  tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv),pert)
! re-evaluate dry mixing ratio based on density scaling
!            tracerVarCellFull%scalar_tracer_mxrt_n%f(1,ilev,iv) =((dycoreVarCellFace%scalar_pressure_n%f(ilev+1,iv)-&
!                                                                  dycoreVarCellFace%scalar_pressure_n%f(ilev,iv))/&
!                                                                 (dycoreVarCellFace%scalar_hpressure_n%f(ilev+1,iv)-&
!                                                                  dycoreVarCellFace%scalar_hpressure_n%f(ilev,iv)))-one
         end do
     end do

     dycoreVarCellFace%scalar_geopotential_n%f = gravity*dycoreVarCellFace%scalar_geopotential_n%f
     dycoreVarCellFull%scalar_geopotential_n%f = gravity*dycoreVarCellFull%scalar_geopotential_n%f
     dycoreVarSurface%scalar_geopotential_n%f    = zero
     tracerVarCellFull%scalar_tracer_mxrt_n%f(2:ntracer,:,:) = zero
#ifdef DYCORE_NHD
     dycoreVarCellFace%scalar_www_n%f          = zero
     dycoreVarCellFace%scalar_phi_n%f          = dycoreVarCellFace%scalar_geopotential_n%f
#endif

     if(iblock.eq.0) print*,"DCMIP2016-SC-B initial done"

   case default

      print*,"you must select a test case in grist_dtp_initial"
      call mpi_abort()

   end select
!
! common procedures, horizontal momentum on the model mesh
! evaluate normal wind at edge
!
   do ie = 1, mesh%ne_halo(1)
     !icell1 = mesh%edt(ie)%v(1)
     !icell2 = mesh%edt(ie)%v(2)
     icell1 = mesh%edt_v(1,ie)
     icell2 = mesh%edt_v(2,ie)
     do ilev = 1, nlev
        utmp = half*(dycoreVarCellFull%scalar_U_wind_n%f(ilev,icell1)+dycoreVarCellFull%scalar_U_wind_n%f(ilev,icell2))
        vtmp = half*(dycoreVarCellFull%scalar_V_wind_n%f(ilev,icell1)+dycoreVarCellFull%scalar_V_wind_n%f(ilev,icell2))
        call convert_vector_sph2cart(utmp, vtmp, real(mesh%edt(ie)%c%p,r8), vector_velocity)
        dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dot_product(vector_velocity, real(mesh%edp(ie)%nr,r8))
     end do
   end do

#ifdef DYCORE_NHD
   dycoreVarCellFace%scalar_www_n%f  =  0._r8
#endif

#ifdef TC_RANDOM_WIND

   call random_seed()

   temp_var = dycoreVarEdgeFull%scalar_normal_velocity_n
   do ilev = 1, nlev
!   if(mpi_rank().eq.0)then
!         call random_number(pertval_u)
!         pertval_u = 2.*pertval_u - 1.
!   end if
!   call mpi_bcast(pertval_u, 1, MPI_REAL8, 0, grist_comm(), ierr)
   !print*,"pertu=",pertval_u
   do ie = 1, mesh%ne_full
         call random_number(pertval_u)
         pertval_u = 2.*pertval_u - 1.
         dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) + 0.02*pertval_u*&
         dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)
   end do
   end do
   temp_var = dycoreVarEdgeFull%scalar_normal_velocity_n
   call exchange_data_2d_add(mesh,field_head_2d,temp_var)
   call exchange_data_2d(mesh%local_block,field_head_2d)
   do ilev = 1, nlev
      do ie = 1, mesh%ne_compute
         if(mesh%edt_v(1,ie).gt.mesh%nv_compute.and.mesh%edt_v(1,ie).le.mesh%nv_halo(1).and.mesh%edt_v(2,ie).le.mesh%nv_compute)then
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = 0.5_r8*(dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)+temp_var%f(ilev,ie))
         end if
         if(mesh%edt_v(2,ie).gt.mesh%nv_compute.and.mesh%edt_v(2,ie).le.mesh%nv_halo(1).and.mesh%edt_v(1,ie).le.mesh%nv_compute)then
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = 0.5_r8*(dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)+temp_var%f(ilev,ie))
         end if
      end do
   end do
#endif

#ifndef SEQ_GRIST
! exchange data
   call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
   call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

#ifdef LAM_DOMAIN
   do ie = 1, mesh%ne_full
      if(mesh%edt_mark(ie).eq.1772)then
         dycoreVarEdgeFull%scalar_normal_velocity_n%f(1:nlev,ie) = zero
      end if
   end do
#endif

   return
  end subroutine grist_dtp_initial

!------------------------------------------------------------------
!  evaluate dry mass and geometric height at each level
!------------------------------------------------------------------

  subroutine evaluate_dry_mass_height(mesh, ncell, nlev, testcase         , &
                                      scalar_hpressure_at_pc_surface      , &
                                      scalar_hpressure_at_pc_face_level   , &
                                      scalar_hpressure_at_pc_full_level   , &
                                      scalar_zzz_at_pc_face_level         , &
                                      scalar_zzz_at_pc_full_level)
!io
   type(global_domain),  intent(in)   :: mesh
   integer(i4),          intent(in)   :: ncell
   integer(i4),          intent(in)   :: nlev
   character(len=*),     intent(in)   :: testcase
   real(r8),             intent(out)  :: scalar_hpressure_at_pc_surface(:)
   real(r8),             intent(out)  :: scalar_hpressure_at_pc_face_level(:,:)
   real(r8),             intent(out)  :: scalar_hpressure_at_pc_full_level(:,:)
   real(r8),             intent(out)  :: scalar_zzz_at_pc_face_level(:,:)
   real(r8),             intent(out)  :: scalar_zzz_at_pc_full_level(:,:)
! local
   integer(i4)         :: iv, ilev
   real(r8)            :: pres, dr_mass, wv_mass, ptop, ztop
!
! obtain top level height
!
    ptop = eta_face_a(1)*p0

    DO iv = 1, ncell
! for each profile
       ztop  = get_ztop(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), ptop, testcase) ! analytic
       !ztop  = get_z_from_pressure(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat, ptop, ptop, 0.0_r8, testcase,.false.) ! iteration

       !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
       !ztop1 = ztop(iv)
       !call  supercell_test(mesh%vtx(iv)%lon, mesh%vtx(iv)%lat, ptop, ztop(iv), 0, tmp1, tmp2, tmp3, tmp4, &
       !                       tmp5, tmp6, tmp7, tmp8, 1)
       !ztop2 = ztop(iv)
       !ztop2 = ztop1-ztop2
       !print*,"ztop2=",ztop2, such difference is 1e-9 for supercell due to different configuration for interation
       !zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
!
! Gaussian quadrature for surface dry and wv mass, their sum should be close to analytic (height=0) produced pressure;
! this has been confirmed;
!
       dr_mass  = get_dryAirMass_from_z(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8),zero,ptop,ztop,testcase)
       !wv_mass  = get_waterVaporMass_from_z(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,zero,ztop,testcase)
       !pres     =   get_dcmipPressure_from_z(mesh%vtx(iv)%lon,mesh%vtx(iv)%lat,zero     ,testcase)
       !scalar_hpressure_at_pc_surface(iv) = pres-wv_mass
       scalar_hpressure_at_pc_surface(iv) = dr_mass
!
!  get surface hpressure based on analytical function
!
#ifdef DTPHPS
       scalar_hpressure_at_pc_surface(iv) = get_surface_hps(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8),testcase)
#endif
!
! below are same
!
       do ilev = 1, nlev+1
          scalar_hpressure_at_pc_face_level(ilev,iv) = eta_face_a(ilev)*p0+eta_face_b(ilev)*scalar_hpressure_at_pc_surface(iv)
       end do

       do ilev = 1, nlev
          scalar_hpressure_at_pc_full_level(ilev,iv) = eta_full_a(ilev)*p0+eta_full_b(ilev)*scalar_hpressure_at_pc_surface(iv)
       end do
!
! iterate to obtain z-face
! 
       do ilev = 1, nlev+1
          scalar_zzz_at_pc_face_level(ilev,iv)  = get_z_from_pressure(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8), &
                                                  scalar_hpressure_at_pc_face_level(ilev,iv), ptop, ztop, testcase,.true.)
       end do
       if(scalar_zzz_at_pc_face_level(nlev+1,iv).ne.0._r8)then
          !print*,"reset to zero", scalar_zzz_at_pc_face_level(nlev+1,iv)
          scalar_zzz_at_pc_face_level(nlev+1,iv) = zero
       end if

       do ilev = 1, nlev
          scalar_zzz_at_pc_full_level(ilev,iv)  = 0.5_r8*(scalar_zzz_at_pc_face_level(ilev,iv)+scalar_zzz_at_pc_face_level(ilev+1,iv))
       end do

    END DO

   return
  end subroutine evaluate_dry_mass_height

!------------------------------------------------------------------
! Evaluate model level height based on pressure or dry mass (dry 
! hydrostatic pressure) using fixed-point iteration; if dry mass 
! is used, the function is an vertical Gauss-integral function; 
! if (full) pressure is used, use the DCMIP analytic function now
!------------------------------------------------------------------

  real(r8) function get_z_from_pressure(lon, lat, pres, ptop, ztop, testcase, is_dry_mass)
! io
    real(r8), intent(in)  ::  lon
    real(r8), intent(in)  ::  lat
    real(r8), intent(in)  ::  pres     ! Pressure (Pa)
    real(r8), intent(in)  ::  ptop     ! Pressure (Pa)
    real(r8), intent(in)  ::  ztop
    character(len=*), intent(in) :: testcase
    logical,  intent(in)  :: is_dry_mass
! local
    integer    :: ix
    real(r8)   :: z0, z1, z2, za, zb, zc
    real(r8)   :: p0, p1, p2, pa, pb, pc


    select case(trim(testcase))

    case('DCMIP2016-BW','DCMIP2016-TC')
!
! some initial
!
    z0 = 0.0_r8
    z1 = 10000.0_r8

    if (is_dry_mass) then
       p0 = get_dryAirMass_from_z(lon,lat,z0,ptop,ztop,testcase)
       p1 = get_dryAirMass_from_z(lon,lat,z1,ptop,ztop,testcase)
    else
       p0 = get_dcmipPressure_from_z(lon,lat,z0,testcase)
       p1 = get_dcmipPressure_from_z(lon,lat,z1,testcase)
    endif
!
! fixed point iteration
!
    DO ix = 1, 1000
       z2 = z1 - (p1 - pres) * (z1 - z0) / (p1 - p0)
       if (is_dry_mass) then
          p2 = get_dryAirMass_from_z(lon,lat,z2,ptop,ztop,testcase)
       else
          p2 = get_dcmipPressure_from_z(lon,lat,z2,testcase)
       end if

       IF (ABS(p2 - pres)/pres < eps.or.ABS(z1-z2)<eps.or.ABS(p1-p2)<eps) THEN
          EXIT
       END IF

       z0 = z1
       p0 = p1

       z1 = z2
       p1 = p2
    END DO

    if (ix==1001) then
      write(*,*) "pres,p1,z1",pres,p1,z1
      print*, 'iteration did not converge in get_z_from_pressure'
      call mpi_abort
    end if
    get_z_from_pressure = z2

    case('DCMIP2016-SC') ! Use the style as in DCMIP-SC supercell_p

    za = 0._r8
    zb = 50000._r8

    if (is_dry_mass) then
       pa = get_dryAirMass_from_z(lon,lat,za,ptop,ztop,testcase)
       pb = get_dryAirMass_from_z(lon,lat,zb,ptop,ztop,testcase)
    else
       pa = get_dcmipPressure_from_z(lon,lat,za,testcase)
       pb = get_dcmipPressure_from_z(lon,lat,zb,testcase)
    endif

    if (pa .lt. pres) then
      write(*,*) 'Requested pressure out of range on bottom, adjust sample interval'
      write(*,*) pa, pres
      stop
    end if
    if (pb .gt. pres) then
      write(*,*) 'Requested pressure out of range on top, adjust sample interval'
      write(*,*) pb, pres
      stop
    end if

    ! Iterate using fixed point method
    !do iter = 1, 20
    zc = (za * (pb - pres) - zb * (pa - pres)) / (pb - pa)
    if (is_dry_mass) then
       pc = get_dryAirMass_from_z(lon,lat,zc,ptop,ztop,testcase)
    else
       pc = get_dcmipPressure_from_z(lon,lat,zc,testcase)
    endif

    do while(abs((pc - pres) / pres) .gt. eps)

      zc = (za * (pb - pres) - zb * (pa - pres)) / (pb - pa)

      if(is_dry_mass) then
        pc = get_dryAirMass_from_z(lon,lat,zc,ptop,ztop,testcase)
      else
        pc = get_dcmipPressure_from_z(lon,lat,zc,testcase)
      endif

      !write(*,*) pc

      !if (abs((pc - p) / p) .lt. 1.d-12) then
      !  exit
      !end if

      if (pc .gt. pres) then
        za = zc
        pa = pc
      else
        zb = zc
        pb = pc
      end if
    end do

    !if (iter .eq. 21) then
    !  write(*,*) 'Iteration failed to converge->iter=',iter
    !  stop
    !end if
    get_z_from_pressure = zc

    case default

    end select

  end function get_z_from_pressure

!------------------------------------------------------------------
!  Vertical Gaussian quadrature for integrating dry air mass,
!  based on full pressure, Rd, Tv (imply alpham using gas law),
!  specific humidity (not dry mixing ratio)
!------------------------------------------------------------------

  real(r8) function get_dryAirMass_from_z(lon, lat, z, ptop,  ztop, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    real(r8), intent(in)  :: z
    real(r8), intent(in)  :: ptop
    real(r8), intent(in)  :: ztop
    character(len=*), intent(in) :: testcase
! local
    real(r8)              :: xm, xr, integral
    real(r8)              :: qv, z1, z2, Tv, pfull, ztmp
    integer               :: jgw
    real(r8)              :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

      z1       = z
      z2       = ztop
      xm       = 0.5*(z1+z2)
      xr       = 0.5*(z2-z1)
      integral = zero

      select case (trim(testcase))
      case('DCMIP2016-BW')
         do jgw = 1, nGauss
            ztmp = xm + gaussx(jgw)*xr
            call baroclinic_wave_test(0, 1, 0, one, lon, lat, pfull, ztmp, 1, &
                                      tmp1, tmp2, tmp3, Tv, tmp5, tmp6, tmp7, qv)
            !this is actually ptv, transform to Tv
             Tv       = Tv*(pfull/p0)**(rdry/cp)
             integral = integral + gaussw(jgw)*gravity*pfull*(one-qv)/(rdry*Tv)
         end do

      case('DCMIP2016-TC')
         do jgw = 1, nGauss
            ztmp = xm + gaussx(jgw)*xr
            call tropical_cyclone_test(lon, lat, pfull , ztmp , 1, &
                                       tmp1, tmp2, tmp3, Tv,  tmp4 ,tmp5, tmp6, qv)
            !this is actually ptv, transform to Tv
            Tv       = Tv*(pfull/p0)**(rdry/cp)
            integral = integral + gaussw(jgw)*gravity*pfull*(one-qv)/(rdry*Tv)
         end do

      case('DCMIP2016-SC')
         do jgw = 1, nGauss
            ztmp = xm+gaussx(jgw)*xr
            call supercell_test(lon, lat, pfull, ztmp, zero, zero,1, &
                                tmp1, tmp2, tmp3, Tv, tmp5, tmp6, tmp7, qv, 0)
            !this is actually ptv, transform to Tv
            Tv       = Tv*(pfull/p0)**(rdry/cp)
            integral = integral + gaussw(jgw)*gravity*pfull*(one-qv)/(rdry*Tv)
         enddo

      case('DCMIP2016-SC-QV') ! treat qv as qv rather than sv, little sensitivity
         do jgw = 1, nGauss
            ztmp = xm+gaussx(jgw)*xr
            call supercell_test(lon, lat, pfull, ztmp, zero, zero,1, &
                                tmp1, tmp2, tmp3, Tv, tmp5, tmp6, tmp7, qv, 0)
            !this is actually ptv, transform to Tv
            Tv       = Tv*(pfull/p0)**(rdry/cp)
            integral = integral + gaussw(jgw)*gravity*pfull*(one-qv/(one+qv))/(rdry*Tv)
         enddo
      case default
         print*,"get_dryAirMass_from_z: this cannot happen"
      end select

      integral             = 0.5_r8*(z2-z1)*integral
      get_dryAirMass_from_z = integral+ptop
  end function get_dryAirMass_from_z

!------------------------------------------------------------------
!  same as above but for mass of water vapor
!------------------------------------------------------------------

  real(r8) function get_waterVaporMass_from_z(lon, lat, z, ztop, testcase)
! io
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in) :: z
    real(r8), intent(in) :: ztop
    character(len=*), intent(in) :: testcase
! local
    real (r8)            :: xm,xr,integral
    real(r8)             :: qv, z1, z2, Tv, pfull, ztmp
    integer              :: jgw
    real(r8)             :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

      z1       = z
      z2       = ztop
      xm       = 0.5_r8*(z1+z2)
      xr       = 0.5_r8*(z2-z1)
      integral = zero

      select case (trim(testcase))

      case('DCMIP2016-BW')
        do jgw = 1, nGauss
           ztmp = xm + gaussx(jgw)*xr
           call baroclinic_wave_test(0, 1, 0, one, lon, lat, pfull, ztmp, 1, &
                                     tmp1, tmp2, tmp3, Tv, tmp5, tmp6, tmp7, qv)
           !this is actually ptv, transform to Tv
            Tv       = Tv*(pfull/p0)**(rdry/cp)
            integral = integral + gaussw(jgw)*gravity*pfull*qv/(rdry*Tv)
        end do

      case('DCMIP2016-TC')
        do jgw = 1, nGauss
           ztmp = xm+gaussx(jgw)*xr
           call tropical_cyclone_test(lon, lat, pfull , ztmp , 1, &
                                      tmp1, tmp2, tmp3, Tv,  tmp4 ,tmp5, tmp6, qv)
           !this is actually ptv, transform to Tv
           Tv       = Tv*(pfull/p0)**(rdry/cp)
           integral = integral + gaussw(jgw)*gravity*pfull*qv/(rdry*Tv)
        enddo

      case('DCMIP2016-SC')
        do jgw = 1, nGauss
           ztmp = xm+gaussx(jgw)*xr
           call supercell_test(lon, lat, pfull, ztmp, zero, zero, 1, &
                               tmp1, tmp2, tmp3, Tv, tmp5, tmp6, tmp7, qv, 0)
           !this is actually ptv, transform to Tv
           Tv       = Tv*(pfull/p0)**(rdry/cp)
           integral = integral + gaussw(jgw)*gravity*pfull*qv/(rdry*Tv)
        enddo

      case default
        print*,"get_waterVaporMass_from_z: this cannot happen"
      end select

      integral                 = 0.5_r8*(z2-z1)*integral
      get_waterVaporMass_from_z = integral

  end function get_waterVaporMass_from_z

!------------------------------------------------------------------
!  produce pressure based on DCMIP analytic functions with height
!------------------------------------------------------------------

  real(r8) function get_dcmipPressure_from_z(lon, lat, z, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    real(r8), intent(in)  :: z
    character(len=*), intent(in) :: testcase
! local
    real(r8)  :: pfull, zlocal
    real(r8)  :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
    zlocal = z

    select case (trim(testcase))
    case('DCMIP2016-BW')
        call baroclinic_wave_test(0, 1, 0, one, lon, lat , pfull, zlocal, 1, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
    case('DCMIP2016-TC')
        call tropical_cyclone_test(lon, lat, pfull , zlocal, 1, &
                                   tmp1, tmp2, tmp3, tmp4,  tmp5 ,tmp6, tmp7, tmp8)
    case('DCMIP2016-SC')
        call supercell_test(lon, lat, pfull, zlocal, zero, zero, 1, &
                            tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, 0)
    case default
       print*,"get_dcmipPressure_from_z: this cannot happen"
    end select
    get_dcmipPressure_from_z = pfull
 
  end function get_dcmipPressure_from_z

!------------------------------------------------------------------
!  produce pressure based on DCMIP analytic functions with height
!------------------------------------------------------------------

  real(r8) function get_surface_hps(lon, lat, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    character(len=*), intent(in) :: testcase
! local
    real(r8)  :: pfull, zlocal
    real(r8)  :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8
    zlocal = zero

    select case (trim(testcase))
    case('DCMIP2016-BW')
        call baroclinic_wave_test(0, 1, 0, one, lon, lat , pfull, zlocal, 1, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
        pfull = pfull/(one+rvap/rdry*tmp8)
    case('DCMIP2016-TC')
        call tropical_cyclone_test(lon, lat, pfull , zlocal, 1, &
                                   tmp1, tmp2, tmp3, tmp4,  tmp5 ,tmp6, tmp7, tmp8)
        pfull = pfull/(one+rvap/rdry*(tmp8/(one-tmp8)))
    case('DCMIP2016-SC')
        call supercell_test(lon, lat, pfull, zlocal, zero, zero, 1, &
                            tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, 0)
        pfull = pfull/(one+rvap/rdry*tmp8)
    case default
       print*,"get_surface_hps: this cannot happen"
    end select
    get_surface_hps = pfull
 
  end function get_surface_hps

  real(r8) function get_ztop(lon, lat, ptop, testcase)
! io
    real(r8), intent(in)  :: lon
    real(r8), intent(in)  :: lat
    real(r8), intent(in)  :: ptop
    character(len=*), intent(in) :: testcase
! local
    real(r8)  :: pfull, zlocal
    real(r8)  :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8

    pfull = ptop

    select case (trim(testcase))
    case('DCMIP2016-BW')
        call baroclinic_wave_test(0, 1, 0, one, lon, lat , pfull, zlocal, 0, &
                                  tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8)
    case('DCMIP2016-TC')
        call tropical_cyclone_test(lon, lat, pfull , zlocal, 0, &
                                   tmp1, tmp2, tmp3, tmp4,  tmp5 ,tmp6, tmp7, tmp8)
    case('DCMIP2016-SC')
        call supercell_test(lon, lat, pfull, zlocal, zero, zero, 0, &
                            tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, 0)
    case default
       print*,"get_dcmipPressure_from_z: this cannot happen"
    end select
    get_ztop = zlocal
 
  end function get_ztop

!---------------------------------------------------------------------------
! A function for iterating dry and moist hydrostatic equartion toward
! the same state
!---------------------------------------------------------------------------

  subroutine grist_dtp_iterative_step(mesh, nlev, testcase                  , &
                                      scalar_hpressure_at_pc_face_level_n   , &
                                      scalar_temp_at_pc_full_level_n        , &
                                      scalar_geopotential_at_pc_full_level_n, &
                                      scalar_geopotential_at_pc_face_level_n)
! io 
   type(global_domain),    intent(in)     :: mesh
   integer(i4)        ,    intent(in)     :: nlev
   character(len=*)   ,    intent(in)     :: testcase
   type(scalar_2d_field),  intent(in)     :: scalar_hpressure_at_pc_face_level_n
   type(scalar_2d_field),  intent(inout)  :: scalar_temp_at_pc_full_level_n 
   type(scalar_2d_field),  intent(inout)  :: scalar_geopotential_at_pc_full_level_n
   type(scalar_2d_field),  intent(inout)  :: scalar_geopotential_at_pc_face_level_n 
! local
   real(r8)                               :: scalar_template_1d_nlev_a(nlev)
   real(r8)                               :: scalar_template_1d_nlev_b(nlev)
   real(r8)                               :: scalar_template_1d_nlev_c(nlev)
   real(r8)                               :: scalar_template_1d_nlevp_a(nlevp)
   real(r8)                               :: scalar_template_1d_nlevp_b(nlevp)
   real(r8)                               :: scalar_template_a
   real(r8)                               :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9
   real(r8)                               :: checker, limit
   integer(i4)                            :: iv, ilev


     do iv = 1, mesh%nv
! init
        scalar_template_1d_nlev_a  = scalar_temp_at_pc_full_level_n%f(:,iv)
        scalar_template_1d_nlevp_a = scalar_hpressure_at_pc_face_level_n%f(:,iv) ! hpressure face
        scalar_template_1d_nlev_b  = scalar_hpressure_at_pc_face_level_n%f(1:nlev,iv)-scalar_hpressure_at_pc_face_level_n%f(2:nlev+1,iv) ! delhp
        scalar_template_a          = zero

        call calc_hpe_hydro(nlev, nlevp,scalar_template_1d_nlev_a  ,&
                            scalar_template_1d_nlevp_a ,&
                            scalar_template_1d_nlev_b  ,&
                            scalar_template_a          ,&
                            scalar_template_1d_nlevp_b ,& ! geop face
                            scalar_template_1d_nlev_c )   ! geop full
        checker = zero
        do ilev = 1, nlev
           checker = checker+abs(scalar_template_1d_nlev_c(ilev)-scalar_geopotential_at_pc_full_level_n%f(ilev,iv)*gravity)
        end do

        if(trim(testcase).eq.'DCMIP2016-BW') limit = 1e-9_r8
        if(trim(testcase).eq.'DCMIP2016-TC') limit = 1e-9_r8
        if(trim(testcase).eq.'DCMIP2016-SC') limit = 1e-7_r8

        DO WHILE(checker.gt.limit)
           scalar_geopotential_at_pc_full_level_n%f(:,iv) = scalar_template_1d_nlev_c/gravity
! new temp profile based on new height full
           if(trim(testcase).eq.'DCMIP2016-BW')then
              do ilev = 1, nlev
                 call baroclinic_wave_test(0, 1, 0, one, real(mesh%vtx(iv)%lon,r8), real(mesh%vtx(iv)%lat,r8) ,&
                                           tmp1,&
                                           scalar_geopotential_at_pc_full_level_n%f(ilev,iv),&
                                           1   ,&
                                           tmp2,&
                                           tmp3,&
                                           scalar_temp_at_pc_full_level_n%f(ilev,iv),&
                                           tmp4,&
                                           tmp5,&
                                           tmp6,&
                                           tmp7,&
                                           tmp8)
              end do
           end if
           if(trim(testcase).eq.'DCMIP2016-TC')then
              do ilev = 1, nlev
                 call tropical_cyclone_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8)             ,&
                                         tmp1 ,&
                                         scalar_geopotential_at_pc_full_level_n%f(ilev,iv), 1,&
                                         tmp2 ,&
                                         tmp3 ,&
                                         scalar_temp_at_pc_full_level_n%f(ilev,iv)  ,&
                                         tmp4 ,&
                                         tmp5 ,&
                                         tmp2 ,&  ! ps repeated
                                         tmp6 ,&
                                         tmp7)
              end do
           end if
           if(trim(testcase).eq.'DCMIP2016-SC')then
              do ilev = 1, nlev
                 call supercell_test(real(mesh%vtx(iv)%lon,r8),real(mesh%vtx(iv)%lat,r8) ,&
                                  tmp1, &
                                  scalar_geopotential_at_pc_full_level_n%f(ilev,iv), zero, zero, 1,&
                                  tmp2, &
                                  tmp3, &
                                  scalar_temp_at_pc_full_level_n%f(ilev,iv),&
                                  tmp4,&
                                  tmp5,& ! thetap
                                  tmp6,& ! ps repeated
                                  tmp7,& ! this is full density
                                  tmp8, 1)
              end do
           end if
! re-evaluate hydro equation based on new temp 
           scalar_template_1d_nlev_a  = scalar_temp_at_pc_full_level_n%f(:,iv)
           call calc_hpe_hydro(nlev,nlevp,scalar_template_1d_nlev_a  ,&
                               scalar_template_1d_nlevp_a ,&
                               scalar_template_1d_nlev_b  ,&
                               scalar_template_a          ,&
                               scalar_template_1d_nlevp_b ,& ! geop face
                               scalar_template_1d_nlev_c)    ! geop full
           checker = zero
           do ilev = 1, nlev
              checker = checker+abs(scalar_template_1d_nlev_c(ilev)-scalar_geopotential_at_pc_full_level_n%f(ilev,iv)*gravity)
           end do

        END DO
        scalar_geopotential_at_pc_face_level_n%f(:,iv) = scalar_template_1d_nlevp_b/gravity
        scalar_geopotential_at_pc_full_level_n%f(:,iv) = scalar_template_1d_nlev_c/gravity

     end do
     return
  end subroutine grist_dtp_iterative_step

 end module grist_dtp_initial_module
