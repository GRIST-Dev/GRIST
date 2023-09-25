
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Prescribe a 3D passive advection environment for testing
!             the tracer transport module
! Revision history:
!----------------------------------------------------------------------------

 module grist_tracer_transport_prescribe_module

  use grist_constants,     only: gravity, pi, p00, i4, r8, rearth,rdry, day2sec, zero
  use grist_domain_types,  only: global_domain
  use grist_data_types,    only: scalar_1d_field, scalar_2d_field,  exchange_field_list_2d
  use grist_nml_module,    only: ntracer, nlev,nlevp, swe_use_streamf
  use grist_math_module,   only: convert_vector_sph2cart,arcdistll
! hpe
  use grist_hpe_constants, only: eta_face_a, eta_face_b, deta_face
  use grist_hpe_hydro_pgf, only: calc_hpe_hpressure_face_level, &
                                 calc_hpe_hpressure_full_level, &
                                 calc_hpe_delhp
! dycore
  use grist_dycore_vars_module, only: dycoreVarCellFace, dycoreVarCellFull, dycoreVarEdgeFull, dycoreVarSurface
! data pool
  use grist_tracer_transport_vars_module, only: tracerVarCellFull, tracerVarCellFace, tracerVarEdgeFull
! data exchange
#ifndef SEQ_GRIST
  use grist_config_partition,   only: exchange_data_2d_add, exchange_data_2d
#endif
! 
  use grist_mpi

  implicit none

  private

  public  :: tracer_transport_prescribe_ambient , & ! called each step
             tracer_transport_prescribe_initial     ! called only init


  CONTAINS

!
! this prescribe routine should be called before prescribe_initial, or initial will fail
!
  subroutine tracer_transport_prescribe_ambient(mesh, time, testcase)
!
! io
!
   type(global_domain), intent(inout) :: mesh
   real(r8)           , intent(in)    :: time
   character(len=*)   , intent(in)    :: testcase
!
! local
!
   real(r8), allocatable              :: stream_function_dual_cell(:) 
   real(r8)                           :: velocity_vector(3)
   real(r8)                           :: scalar_u, scalar_v
   real(r8)                           :: pres, height, ptop, ztop, delhp
   real(r8)                           :: scalar_normal_velocity_at_edge_uava
   integer(i4)                        :: ilev, ie, iv, icell1, icell2
   integer(i4)                        :: v1,v0, it, flag
   real(r8)                           :: v0v1(3)
   type(exchange_field_list_2d),pointer:: field_head_2d
   
      field_head_2d =>null()

      select case(trim(testcase))
      case('DCMIP1-1','DCMIP1-2','DCMIP1-1-H','DCMIP11H')
         dycoreVarSurface%scalar_hpressure_n%f  = p00
      case('DCMIP1-3') ! affect bit reproduce in intel
         do iv = 1, mesh%nv
           dycoreVarSurface%scalar_hpressure_n%f(iv) = ini_surface_pressure(real(mesh%vtx_lon(iv),r8),real(mesh%vtx_lat(iv),r8),trim(testcase))
         end do
      case default
         print*,"you must select a test case in tracer transport"
         stop
      end select

      call time_integration_renew_mass_state(mesh, dycoreVarSurface%scalar_hpressure_n    ,&
                                                   dycoreVarCellFace%scalar_hpressure_n ,&
                                                   dycoreVarCellFull%scalar_delhp_n     ,&
                                                   dycoreVarCellFull%scalar_hpressure_n ,&
                                                   dycoreVarCellFace%scalar_delhp_n     )

      if(trim(testcase).eq.'DCMIP11H') dycoreVarCellFull%scalar_delhp_n%f = 1._r8

      !scalar_delhp_at_pc_full_level_beg_adv%f = dycoreVarCellFull%scalar_delhp_n%f
      tracerVarCellFull%scalar_delhp_avg_adv%f = dycoreVarCellFull%scalar_delhp_n%f
      tracerVarCellFull%scalar_delhp_end_adv%f = dycoreVarCellFull%scalar_delhp_n%f

      do ie = 1, mesh%ne
          !icell1 = mesh%edt(ie)%v(1)
          !icell2 = mesh%edt(ie)%v(2)
          icell1 = mesh%edt_v(1,ie)
          icell2 = mesh%edt_v(2,ie)
         do ilev = 1, nlev

            ptop   = 0.5_r8*(dycoreVarCellFace%scalar_hpressure_n%f(1,icell1)+&
                             dycoreVarCellFace%scalar_hpressure_n%f(1,icell2))
            pres   = 0.5_r8*(dycoreVarCellFull%scalar_hpressure_n%f(ilev,icell1)+&
                             dycoreVarCellFull%scalar_hpressure_n%f(ilev,icell2))
            delhp  = 0.5_r8*(dycoreVarCellFull%scalar_delhp_n%f(ilev,icell1)+&
                             dycoreVarCellFull%scalar_delhp_n%f(ilev,icell2))
            call ini_velocity_vector(real(mesh%edt_c_p(1:3,ie),r8), &
                                     real(mesh%edt_c_lon(ie),r8) , &
                                     real(mesh%edt_c_lat(ie),r8) , &
                                     time,ptop,pres     , &
                                     trim(testcase)     , &
                                     velocity_vector    , &
                                     scalar_u           , &
                                     scalar_v)

            dycoreVarEdgeFull%scalar_U_wind_n%f(ilev,ie) = scalar_u
            dycoreVarEdgeFull%scalar_V_wind_n%f(ilev,ie) = scalar_v
          !
          ! for normal velocity
          !
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)        = dot_product(velocity_vector, mesh%edp_nr(1:3,ie))
            tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f(ilev,ie)  = dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)
            tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie) = dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)*delhp

         end do
      end do

!===================================================================
! use streamf here:
! If we use streamf to initialize the nondivergence part of winds,
! the above ini_velocity_vector will only provide the ud component and
! transform it to normal component; the ua and va contribution of 
! normal wind is given by the streamfunction; this leads to an
! exactly non-divergent ua&va component numerically;
! or it is not so due to the mesh irregularity
! If not, the above ini_velocity_vector will provide full velocity
!===================================================================

    if(swe_use_streamf)then
      if(.not.allocated(stream_function_dual_cell)) allocate(stream_function_dual_cell(mesh%nt_full))
      do it = 1, mesh%nt_full
         stream_function_dual_cell(it)  = ini_stream_function(real(mesh%tri_c_lon(it),r8),real(mesh%tri_c_lat(it),r8),time,trim(testcase))
      end do

      do ie = 1, mesh%ne
         !icell1 = mesh%edt(ie)%v(1)
         !icell2 = mesh%edt(ie)%v(2)
         icell1 = mesh%edt_v(1,ie)
         icell2 = mesh%edt_v(2,ie)
         !v0     = mesh%edp(ie)%v(1)
         !v1     = mesh%edp(ie)%v(2)
         v0     = mesh%edp_v(1,ie)
         v1     = mesh%edp_v(2,ie)
         v0v1   = mesh%tri_c_p(1:3,v1)-mesh%tri_c_p(1:3,v0)
         flag   = sign(1._r8, dot_product(v0v1,real(mesh%edp_tg(1:3,ie),r8)))
         scalar_normal_velocity_at_edge_uava = flag*(stream_function_dual_cell(v1)-stream_function_dual_cell(v0))/(rearth*mesh%edp_leng(ie))
! rewrite above code
         do ilev = 1, nlev
            dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie) = dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)+&
                                                                     scalar_normal_velocity_at_edge_uava
            delhp  = 0.5_r8*(dycoreVarCellFull%scalar_delhp_n%f(ilev,icell1)+dycoreVarCellFull%scalar_delhp_n%f(ilev,icell2))
            tracerVarEdgeFull%scalar_normal_velocity_avg_adv%f(ilev,ie)  = dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)
            tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv%f(ilev,ie) = dycoreVarEdgeFull%scalar_normal_velocity_n%f(ilev,ie)*delhp
         end do
      end do
      deallocate(stream_function_dual_cell)
   end if

! exchange data here
#ifdef USE_HALO2
    call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
    call exchange_data_2d_add(mesh,field_head_2d,tracerVarEdgeFull%scalar_normal_velocity_avg_adv)
    call exchange_data_2d_add(mesh,field_head_2d,tracerVarEdgeFull%scalar_normal_mass_flux_avg_adv)
    call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

! presribed generalized vertical velocity
      do iv = 1, mesh%nv
         do ilev = 2, nlev
            pres      = eta_face_a(ilev)*p00+eta_face_b(ilev)*dycoreVarSurface%scalar_hpressure_n%f(iv)
            ptop      = eta_face_a(1)*p00+eta_face_b(1)*dycoreVarSurface%scalar_hpressure_n%f(iv)
            height    = (rdry*300._r8/gravity)*log(p00/pres)
            ztop      = (rdry*300._r8/gravity)*log(p00/ptop)
            tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(ilev,iv) = prescribe_metadot(&
                                                                       time,            &
                                                                       real(mesh%vtx_lon(iv),r8),&
                                                                       real(mesh%vtx_lat(iv),r8),&
                                                                       real(eta_face_b(ilev),r8),&
                                                                       pres,height,ptop,ztop,&
                                                                       testcase)
         end do
      end do
      tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(1,:)      = 0._r8
      tracerVarCellFace%scalar_eta_mass_flux_avg_adv%f(nlev+1,:) = 0._r8

   return
  end subroutine tracer_transport_prescribe_ambient
!
! first prescribe ambient, then prescribe initial
!
  subroutine tracer_transport_prescribe_initial(mesh,testcase)
!
! io
!
   type(global_domain), intent(inout) :: mesh
   character(len=*)   , intent(in)    :: testcase
!
! local
!
   real(r8)                           :: ptop, pres, height
   real(r8)                           :: scalar(ntracer), scalar_min, scalar_min_global
   integer(i4)                        :: itracer, ilev, iv, ierr

   select case(trim(testcase))

   case('DCMIP1-1','DCMIP1-2','DCMIP1-3','DCMIP1-1-H','DCMIP11H')
      do iv = 1, mesh%nv
         do ilev = 1, nlev
            pres   = dycoreVarCellFull%scalar_hpressure_n%f(ilev,iv)
            height = (rdry*300._r8)/gravity*log(p00/pres)
            call ini_tracer_scalar(real(mesh%vtx_lon(iv),r8), real(mesh%vtx_lat(iv),r8), height, scalar, ntracer,trim(testcase))
            tracerVarCellFull%scalar_tracer_mxrt_n%f(:,ilev,iv) = scalar(:)
         end do
      end do
! set initial minimum
!      do itracer = 1, ntracer
!         scalar_min = tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,1,1)
!         do iv = 1, mesh%nv
!            do ilev = 1, nlev
!               scalar_min = min(scalar_min, tracerVarCellFull%scalar_tracer_mxrt_n%f(itracer,ilev,iv))
!            end do
!         end do
!         call reduce(scalar_min, scalar_min_global, 'min')
!         scalar_tracer_mxrt_min(itracer) = scalar_min_global
!      end do
   case default
      print*,"you must select a test case in tracer transport"
      stop
   end select
! set geop and temp to avoid energy initialization failure
   dycoreVarCellFull%scalar_geopotential_n%f = zero
   dycoreVarCellFace%scalar_geopotential_n%f = zero
   dycoreVarSurface%scalar_geopotential_n%f  = zero
   dycoreVarCellFull%scalar_temp_n%f         = zero

   return
  end subroutine tracer_transport_prescribe_initial

!================================================
! Private Routines 
! convert UV wind component to vector velocity
!================================================

  subroutine ini_velocity_vector(point_cart,           &
                                 lon, lat, time, ptop, pres, &
                                 testcase,             &
                                 velocity_vector,      &
                                 scalar_u,             &
                                 scalar_v)
!
! io
!
    real(r8),        intent(in)    :: point_cart(3)
    real(r8),        intent(in)    :: lon
    real(r8),        intent(in)    :: lat
    real(r8),        intent(in)    :: time
    real(r8),        intent(in)    :: ptop, pres
    character*(*),   intent(in)    :: testcase
    real(r8),        intent(inout) :: velocity_vector(3)
    real(r8),        intent(inout) :: scalar_u
    real(r8),        intent(inout) :: scalar_v

      scalar_u = ini_uwind(lon,lat,ptop,pres,time,testcase)
      scalar_v = ini_vwind(lon,lat,ptop,pres,time,testcase)
      
      call convert_vector_sph2cart(scalar_u, scalar_v, point_cart, velocity_vector)

    return
  end subroutine ini_velocity_vector

!================================================
!  U in West-East direction at a point
!================================================

  real(r8) function ini_uwind(lon, lat, ptop, pres, time, testcase)
! io
    real(r8),       intent(in) :: lon
    real(r8),       intent(in) :: lat
    real(r8),       intent(in) :: ptop
    real(r8),       intent(in) :: pres
    real(r8),       intent(in) :: time
    character*(*),  intent(in) :: testcase
!
! local
!
    real(r8)                   :: duration, kkk
    real(r8)                   :: lonp,omega0
    real(r8)                   :: ua, ud, u0, alpha

! Velocity for each testcase
    select case(trim(testcase))
    case('DCMIP1-1')
        duration = 12._r8*day2sec
        kkk      = 10._r8*rearth/duration
        lonp     = lon-(2._r8*pi*time/duration)
        ua       = kkk*(sin(lonp)**2)*(sin(2_r8*lat))*(cos(pi*time/duration))+(2._r8*pi*rearth*cos(lat))/duration
        omega0   = 23000._r8*pi/1036800._r8
        ud       = omega0*rearth/(0.2_r8*ptop)*cos(lonp)*(cos(lat)**2)*cos(2*pi*time/duration)*&
                   (-exp((pres-p00)/(0.2_r8*ptop))+exp((ptop-pres)/(0.2_r8*ptop)))
        ini_uwind= ua+ud
        if(swe_use_streamf) ini_uwind = ud
    case('DCMIP1-1-H','DCMIP11H')
        duration = 12._r8*day2sec
        kkk      = 10._r8*rearth/duration
        lonp     = lon-(2._r8*pi*time/duration)
        ua       = kkk*(sin(lonp)**2)*(sin(2_r8*lat))*(cos(pi*time/duration))+(2._r8*pi*rearth*cos(lat))/duration
        ini_uwind= ua
        if(swe_use_streamf) ini_uwind = 0._r8
    case('DCMIP1-2')
        ini_uwind = 40._r8*cos(lat)
    case('DCMIP1-3')
        duration  = 12._r8*day2sec
        u0        = 2._r8*pi*rearth/duration
        alpha     = pi/6._r8
        ini_uwind = u0*(cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
    case default
       print*, "you must select a test case in tracer transport, please check!"
       stop
    end select

    return
  end function ini_uwind

!================================================
!  V in South-North direction at a point
!================================================

  real(r8) function ini_vwind(lon, lat, ptop, pres, time, testcase)
! io
    real(r8),       intent(in) :: lon
    real(r8),       intent(in) :: lat
    real(r8),       intent(in) :: ptop
    real(r8),       intent(in) :: pres
    real(r8),       intent(in) :: time
    character*(*),  intent(in) :: testcase
! local
    real(r8)                   :: duration, kkk
    real(r8)                   :: lonp, omega0, height
    real(r8)                   :: rho, rho0, ztop, u0, alpha

    select case(trim(testcase))

    case('DCMIP1-1','DCMIP1-1-H','DCMIP11H')

       duration = 12._r8*day2sec
       kkk      = 10._r8*rearth/duration
       lonp     = lon-(2._r8*pi*time/duration)
       ini_vwind= kkk*(sin(2*lonp))*(cos(lat))*(cos(pi*time/duration))

       if(swe_use_streamf) ini_vwind = 0._r8
    case('DCMIP1-2')

       rho0   =  p00/(300._r8*rdry)
       rho    = pres/(300._r8*rdry)
       height = (rdry*300._r8)/gravity*log(p00/pres)
       ztop   = (rdry*300._r8)/gravity*log(p00/ptop)
       ini_vwind = -rearth*0.15_r8*pi*rho0*cos(lat)*sin(5._r8*lat)*&
                                           cos(pi*height/ztop)*&
                                           cos(pi*time/day2sec)/(5._r8*ztop*rho)
    case('DCMIP1-3')
        duration  = 12._r8*day2sec
        u0        = 2._r8*pi*rearth/duration
        alpha     = pi/6._r8
        ini_vwind = -u0*sin(lon)*sin(alpha)
    case default
       print*, "you must select a test case in tracer transport, please check!"
       stop
    end select

    return
  end function ini_vwind

  real(r8) function prescribe_metadot(time,lon,lat,etab,pres,height,ptop,ztop,testcase)
! io
    real(r8),       intent(in) :: time 
    real(r8),       intent(in) :: lon
    real(r8),       intent(in) :: lat
    real(r8),       intent(in) :: etab
    real(r8),       intent(in) :: pres
    real(r8),       intent(in) :: height
    real(r8),       intent(in) :: ptop
    real(r8),       intent(in) :: ztop
    character*(*),  intent(in) :: testcase
! local
    real(r8)                   :: u0, w0, alpha, lonp, sp,duration,rho0
    real(r8)                   :: z_s,lon_m,lat_m, r_m, rr_m, kesi_m
    real(r8)                   :: utmp,vtmp
    real(r8)                   :: partial_rm_lon, partial_rm_lat
    real(r8)                   :: partial_zs_lon, partial_zs_lat
    real(r8)                   :: partial_ps_lon, partial_ps_lat
    real(r8)                   :: partial_zz_lon, partial_zz_lat

    select case(trim(testcase))
    case('DCMIP1-1')
       duration = 1036800._r8
       lonp  = lon-2._r8*pi*time/duration
       sp    = 1._r8+exp(5._r8*(ptop-p00)/ptop)-exp(5._r8*(pres-p00)/ptop)-exp(5._r8*(ptop-pres)/ptop)
       prescribe_metadot = (23000._r8*pi/duration)*sin(lonp)*cos(lat)*cos(2*pi*time/duration)*sp
    case('DCMIP1-1-H','DCMIP11H')
       prescribe_metadot = 0._r8
    case('DCMIP1-2')
       duration = 86400._r8
       rho0     = p00/(rdry*300._r8)
       w0       = 0.15_r8
       prescribe_metadot = (-gravity*w0*rho0/5._r8)*(-2._r8*sin(5._r8*lat)*sin(lat)+5._r8*cos(lat)*cos(5*lat))*&
                            sin(pi*height/ztop)*cos(pi*time/duration)
    case('DCMIP1-3')

       alpha  = pi/6._r8
       lon_m  = 1.5_r8*pi
       lat_m  = 0._r8
       rr_m   = 3._r8*pi/4._r8
       kesi_m = pi/16._r8

       r_m    = acos(sin(lat_m)*sin(lat)+cos(lat_m)*cos(lat)*cos(lon-lon_m))
       partial_rm_lon  = cos(lat_m)*cos(lat)*sin(lon-lon_m)!/sqrt(1._r8-cos(r_m)**2)
       partial_rm_lat  =-sin(lat_m)*cos(lat)+cos(lat_m)*sin(lat)*cos(lon-lon_m)!/sqrt(1._r8-cos(r_m)**2)

       if(r_m.lt.rr_m.and.(1._r8-cos(r_m))**2.gt.0._r8)then
          partial_zs_lon  = partial_rm_lon/sqrt(1._r8-cos(r_m)**2)*&
                           (-1000._r8*pi/rr_m*sin(pi*r_m/rr_m)*(cos(pi*r_m/kesi_m)**2)-&
                             2000._r8*pi/kesi_m*(1._r8+cos(pi*r_m/rr_m))*cos(pi*r_m/kesi_m)*sin(pi*r_m/kesi_m)) 
          partial_zs_lat  = partial_rm_lat/sqrt(1._r8-cos(r_m)**2)*&
                           (-1000._r8*pi/rr_m*sin(pi*r_m/rr_m)*(cos(pi*r_m/kesi_m)**2)-&
                             2000._r8*pi/kesi_m*(1._r8+cos(pi*r_m/rr_m))*cos(pi*r_m/kesi_m)*sin(pi*r_m/kesi_m)) 
          z_s  = 1000._r8*(1._r8+cos(pi*r_m/rr_m))*(cos(pi*r_m/kesi_m)**2)
       else
          partial_zs_lon  = 0._r8
          partial_zs_lat  = 0._r8
          z_s             = 0._r8
       end if

       partial_ps_lon     = -(gravity*p00/(rdry*300._r8))*exp(-gravity*z_s/(rdry*300._r8))*partial_zs_lon
       partial_ps_lat     = -(gravity*p00/(rdry*300._r8))*exp(-gravity*z_s/(rdry*300._r8))*partial_zs_lat
       partial_zz_lon     = -rdry*300._r8/(gravity*pres)*etab*partial_ps_lon
       partial_zz_lat     = -rdry*300._r8/(gravity*pres)*etab*partial_ps_lat
       u0                 = 2*pi*rearth/1036800._r8
       utmp               = u0*(cos(lat)*cos(alpha)+sin(lat)*cos(lon)*sin(alpha))
       vtmp               =-u0*sin(lon)*sin(alpha)
       prescribe_metadot  = -utmp*partial_zz_lon/(rearth*cos(lat))-vtmp*partial_zz_lat/rearth
       prescribe_metadot  = -gravity*(pres/(rdry*300._r8))*prescribe_metadot

    case default
       print*, "you must select a test case in tracer transport, please check!"
       stop
    end select

    return
  end function prescribe_metadot 

  subroutine ini_tracer_scalar(lon, lat, height, scalar, ntracer,testcase)
! io
    real(r8),       intent(in) :: lon
    real(r8),       intent(in) :: lat
    real(r8),       intent(in) :: height
    real(r8),       intent(out):: scalar(ntracer)
    integer(i4),    intent(in) :: ntracer
    character(len=*)   , intent(in)    :: testcase
! local
    real(r8)                   :: d1, d2
    real(r8)                   :: r1, r2
    real(r8)                   :: z0, z1, z2
    real(r8)                   :: r_p, rr_p, lon_p, lat_p
    real(r8)                   :: r_z1, r_z2, r_z3, z_p1, z_p2, z_p3, dz_p1, dz_p2, dz_p3
    real(r8)                   :: height_t, height_c, radius_t
    real(r8)                   :: lonc1, lonc2, latc ! center position
! for dcmip11h
    real(r8)                   :: lon1,lat1,lon2,lat2
    real(r8)                   :: x0,y0,x1,y1,x2,y2
    real(r8)                   :: hmax,b_value,h1,h2
    real(r8)                   :: base_radius,c_value

    select case(trim(testcase))
    case('DCMIP1-1','DCMIP1-1-H')

        if(ntracer.ne.4)then
           print*,"DCMIP1-1 must choose ntracer = 4, stop"
           stop
        end if
! parameter
        lonc1    = 5._r8*pi/6._r8
        lonc2    =-5._r8*pi/6._r8
        latc     = 0._r8
        height_c = 5000._r8
        height_t = 1000._r8
        radius_t = rearth/2._r8
 
        r1  = rearth*acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc1))
        r2  = rearth*acos(sin(latc)*sin(lat)+cos(latc)*cos(lat)*cos(lon-lonc2))
        d1  = min(1._r8,(r1/radius_t)**2+((height-height_c)/height_t)**2)
        d2  = min(1._r8,(r2/radius_t)**2+((height-height_c)/height_t)**2)
        scalar(1) = 1._r8+0.5_r8*(cos(pi*d1)+cos(pi*d2))
        scalar(2) = 0.9_r8-0.8_r8*scalar(1)*scalar(1)
        scalar(3) = 0.1_r8
        if(d1.lt.0.5_r8)then
           scalar(3) = 1._r8
        else if(d2.lt.0.5_r8)then
           scalar(3) = 1._r8
        else
           scalar(3) = 0.1_r8
        end if
        if(height.gt.height_c.and. &
           lat.lt.(latc+1._r8/8._r8).and. &
           lat.gt.(latc-1._r8/8._r8)) scalar(3) = 0.1_r8
        scalar(4) = 1._r8-0.3_r8*(scalar(1)+scalar(2)+scalar(3))

    case('DCMIP11H') ! initial in Lauritzen et al. (2012)

       if(ntracer.ne.4)then
          print*,"DCMIP11H must choose ntracer = 4, stop"
          stop
       end if
! gaussian hill 
       lon1     = -5_r8*pi/6._r8
       lat1     = 0._r8
       lon2     = 5._r8*pi/6._r8
       lat2     = 0._r8

       hmax     = 0.95_r8
       b_value  = 5._r8

       x0   = cos(lat)*cos(lon)
       y0   = cos(lat)*sin(lon)
       z0   = sin(lat)

       x1   = cos(lat1)*cos(lon1)
       y1   = cos(lat1)*sin(lon1)
       z1   = sin(lat1)

       x2   = cos(lat2)*cos(lon2)
       y2   = cos(lat2)*sin(lon2)
       z2   = sin(lat2)

       h1   = hmax*exp(-b_value*((x0-x1)**2+(y0-y1)**2+(z0-z1)**2))
       h2   = hmax*exp(-b_value*((x0-x2)**2+(y0-y2)**2+(z0-z2)**2))
       scalar(1) = h1+h2
! cosine bell
       hmax        = 0.5_r8
       base_radius = 0.5_r8
       b_value     = 0.1_r8
       c_value     = 0.9_r8

       r1 = arcdistll(lon, lat, lon1, lat1)
       r2 = arcdistll(lon, lat, lon2, lat2)
       h1 = hmax*(1._r8+cos(pi*r1/base_radius))
       h2 = hmax*(1._r8+cos(pi*r2/base_radius))
       if(r1.lt.base_radius)then
          scalar(2) = b_value+c_value*h1
       else if(r2.lt.base_radius)then
          scalar(2) = b_value+c_value*h2
       else
          scalar(2) = b_value
       end if
! slotted cylinder
       base_radius = 0.5_r8
       c_value     = 1._r8
       b_value     = 0.1_r8
       r1          = arcdistll(lon, lat, lon1, lat1)
       r2          = arcdistll(lon, lat, lon2, lat2)
       if(r1.le.base_radius.and.abs(lon-lon1).ge.(base_radius/6_r8))then
         scalar(3) = c_value
       else if(r2.le.base_radius.and.abs(lon-lon2).ge.(base_radius/6_r8))then
         scalar(3) = c_value
       else if(r1.le.base_radius.and.abs(lon-lon1).lt.(base_radius/6_r8).and.&
               (lat-lat1).lt.(-5_r8*base_radius/12_r8))then
         scalar(3) = c_value
       else if(r2.le.base_radius.and.abs(lon-lon2).lt.(base_radius/6_r8).and.&
               (lat-lat2).gt.(5_r8*base_radius/12_r8))then
         scalar(3) = c_value
       else
         scalar(3) = b_value
       end if
! corrleted cosine bell
       scalar(4) = -0.8_r8*scalar(2)*scalar(2)+0.9_r8

    case('DCMIP1-2')

       if(ntracer.ne.1)then
          print*,"DCMIP1-2 must choose ntracer = 1, stop"
          stop
       end if

       z1 = 2000._r8
       z2 = 5000._r8
       z0 = 0.5_r8*(z1+z2)
       if(height.gt.z1.and.height.lt.z2)then
         scalar(1) = 0.5_r8*(1._r8+cos(2._r8*pi*(height-z0)/(z2-z1))) 
       else
         scalar(1) = 0._r8
       end if
       scalar(1:ntracer) = scalar(1)

    case('DCMIP1-3')

       if(ntracer.ne.4)then
          print*,"DCMIP1-3 must choose ntracer = 4, stop"
          stop
       end if

       lon_p = 0.5_r8*pi
       lat_p = 0._r8
       z_p1  = 3050._r8
       z_p2  = 5050._r8
       z_p3  = 8200._r8
       dz_p1 = 1000._r8
       dz_p2 = 1000._r8
       dz_p3 = 400._r8
       rr_p  = 0.25_r8*pi

       r_z1  = abs(height-z_p1)
       r_z2  = abs(height-z_p2)
       r_z3  = abs(height-z_p3)
       r_p   = acos(sin(lat_p)*sin(lat)+cos(lat_p)*cos(lat)*cos(lon-lon_p))
! q1
       if(r_z1.lt.0.5_r8*dz_p1.and.r_p.lt.rr_p)then
          scalar(1) = 0.25_r8*(1._r8+cos(2*pi*r_z1/dz_p1))*(1._r8+cos(pi*r_p/rr_p))
       else
          scalar(1) = 0._r8
       end if
! q2
       if(r_z2.lt.0.5_r8*dz_p2.and.r_p.lt.rr_p)then
          scalar(2) = 0.25_r8*(1._r8+cos(2*pi*r_z2/dz_p2))*(1._r8+cos(pi*r_p/rr_p))
       else
          scalar(2) = 0._r8
       end if
! q3
       if(r_z3.lt.0.5_r8*dz_p3.and.r_p.lt.rr_p)then
          scalar(3) = 1._r8
       else
          scalar(3) = 0._r8
       end if
! q4
       scalar(4) = scalar(1)+scalar(2)+scalar(3)

    case default
       print*, "you must select a test case in tracer transport, please check!"
       stop
    end select

    return
  end subroutine ini_tracer_scalar

  real(r8) function ini_surface_pressure(lon,lat,testcase)
! io
     real(r8)           , intent(in)    :: lon, lat 
     character(len=*)   , intent(in)    :: testcase
! local
     real(r8)             :: z_s,lon_m,lat_m, r_m, rr_m, kesi_m
 
 
     select case(trim(testcase))
     case('DCMIP1-3')
       lon_m  = 1.5_r8*pi
       lat_m  = 0._r8
       rr_m   = 3._r8*pi/4._r8
       kesi_m = pi/16._r8
       r_m    = acos(sin(lat_m)*sin(lat)+cos(lat_m)*cos(lat)*cos(lon-lon_m))
       if(r_m.lt.rr_m)then
          z_s  = 1000._r8*(1._r8+cos(pi*r_m/rr_m))*(cos(pi*r_m/kesi_m)**2)
       else
          z_s  = 0._r8
       end if
       ini_surface_pressure = p00*exp(-gravity*z_s/(rdry*300._r8))
     case default
       print*,"you must ini surface pressure, stop!"
       stop
     end select
      
     return
  end function ini_surface_pressure
!
! only used in prescribed testing of tracer transport 
!
  subroutine time_integration_renew_mass_state(mesh, scalar_hpressure_at_pc_surface    ,&
                                                      scalar_hpressure_at_pc_face_level ,&
                                                      scalar_delhp_at_pc_full_level     ,&
                                                      scalar_hpressure_at_pc_full_level ,&
                                                      scalar_delhp_at_pc_face_level     )
! io
     type(global_domain),   intent(in)    :: mesh
     type(scalar_1d_field), intent(in)    :: scalar_hpressure_at_pc_surface
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_face_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_face_level
! local
     integer(i4)                          :: iv, ilev
     real(r8)                             :: scalar_template_a
     real(r8), dimension(nlevp)           :: scalar_template_1d_nlevp_a
     real(r8), dimension(nlev)            :: scalar_template_1d_nlev_a
     real(r8), dimension(nlev)            :: scalar_template_1d_nlev_b

!        if(.not.allocated(scalar_template_1d_nlevp_a%f)) allocate(scalar_template_1d_nlevp_a%f(nlev+1))
!        if(.not.allocated(scalar_template_1d_nlev_a%f))  allocate(scalar_template_1d_nlev_a%f(nlev))
!        if(.not.allocated(scalar_template_1d_nlev_b%f))  allocate(scalar_template_1d_nlev_b%f(nlev))

!
! renew mass state, compute hpressure at face level, full level, 
! delhp, based on input surface hpressure
!
        do iv = 1, mesh%nv

           scalar_template_a   =  scalar_hpressure_at_pc_surface%f(iv)

           call calc_hpe_hpressure_face_level(nlev,nlevp,scalar_template_a, &       ! in
                                              scalar_template_1d_nlevp_a) ! out

           call calc_hpe_delhp(nlev,nlevp,scalar_template_1d_nlevp_a, &  ! in
                               scalar_template_1d_nlev_a)     ! out

           call calc_hpe_hpressure_full_level(nlev,nlevp,scalar_template_a         , & ! in
                                              scalar_template_1d_nlevp_a, & ! in
                                              scalar_template_1d_nlev_a , & ! in
                                              scalar_template_1d_nlev_b)    ! out

           scalar_hpressure_at_pc_face_level%f(:,iv)  = scalar_template_1d_nlevp_a
               scalar_delhp_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_a
           scalar_hpressure_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_b

           do ilev = 2, nlev
              scalar_delhp_at_pc_face_level%f(ilev,iv) = scalar_hpressure_at_pc_full_level%f(ilev,iv)-&
                                                         scalar_hpressure_at_pc_full_level%f(ilev-1,iv)
           end do

           scalar_delhp_at_pc_face_level%f(1,iv)      = 0.5_r8*scalar_delhp_at_pc_full_level%f(1,iv)
           scalar_delhp_at_pc_face_level%f(nlev+1,iv) = 0.5_r8*scalar_delhp_at_pc_full_level%f(nlev,iv)

        end do

        return
  end subroutine time_integration_renew_mass_state

  function ini_stream_function(lon, lat, time,testcase)
! io
    real(r8), intent(in) :: lon
    real(r8), intent(in) :: lat
    real(r8), intent(in), optional :: time
    character(len=*)   , intent(in)    :: testcase
    real(r8)             :: ini_stream_function
! local
    real(r8)             :: duration, lonp
    real(r8)             :: ooo, kkk, rrr

    select case(trim(testcase))
    case('DCMIP1-1','DCMIP1-1-H','DCMIP11H')! used by deformatinoal flow of two cosine bells and gaussian hills 
        duration = 12._r8*day2sec
        kkk      = 10._r8*rearth/duration
        lonp     = lon-(2._r8*pi*time/duration)
        ini_stream_function = rearth*(kkk*((sin(lonp))**2)*(cos(lat)**2)*(cos(pi*time/duration))-&
                              (2._r8*pi*rearth*sin(lat))/duration)
    case default
        ini_stream_function = 0._r8
    end select
    return
  end function ini_stream_function

  end module grist_tracer_transport_prescribe_module
