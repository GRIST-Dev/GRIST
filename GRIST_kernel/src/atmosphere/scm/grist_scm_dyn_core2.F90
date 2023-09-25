!======================================================
!
!  Created by LiXiaohan on 20/4/20.
!  
!  Single Column Model for GRIST.
!  Adopted from the GRIST dynamical core.
!
!  Purpose:
!  SCM dynamical core
!  Note: Remove dyn vars to an isolated module
!======================================================

 module grist_scm_dyn_core2

    use grist_constants,                    only: i4, r8, rdry, cp
    use grist_nml_module,                   only: nlev, nlevp, ntracer, &
                                                  nrk,                  &
                                                  ver_adv_flag,         &
                                                  tracer_vadv_flag,     &
                                                  eqs_vert_diff
    use grist_scm_pressure_update,          only: time_integration_renew_mass_state
    use grist_physics_data_structure,       only: phy_tracer_info
    use grist_scm_io_module
    use grist_scm_comm_module
    use grist_handle_error,                 only: endrun
    use grist_mpi
    use grist_scm_dyn_vars ! all

    implicit none

    private

    public  :: grist_scm_dyn_init2,  &
               grist_scm_dyn_end2,   &
               grist_scm_dyn_pkg2

 contains
 
    subroutine grist_scm_dyn_init2

! local
    integer  :: m,k

    !two time levels for GRIST_SCM
    n3   = 1
    n3m1 = n3
    n3m2 = n3m1

    allocate(scm_u(nlev));            scm_u      = 0._r8
    allocate(scm_v(nlev));            scm_v      = 0._r8
    allocate(scm_omega(nlev));        scm_omega  = 0._r8
    allocate(scm_ps3(1));             scm_ps3    = 0._r8
    allocate(scm_t3(nlev,1));         scm_t3     = 0._r8
    allocate(scm_q3(ntracer,nlev,1)); scm_q3     = 0._r8

    allocate(scm_t_tend(nlev));       scm_t_tend = 0._r8
    allocate(scm_u_tend(nlev));       scm_u_tend = 0._r8
    allocate(scm_v_tend(nlev));       scm_v_tend = 0._r8
    allocate(scm_q_tend(ntracer,nlev));  scm_q_tend = 0._r8
  
!    call scam_data_io_read('/GFPS8p/cess31/lixh/grist/init_data/scm_data/O3_set0/S6.cam.h1.2003-07-15-00000_CTL_relax.nc')
    call scam_data_io_read(scam_file_name)

    call grist_scm_io_read(0)

    scm_u(:)        = uobs
    scm_v(:)        = vobs
    scm_omega(:)    = wfld
    scm_ps3(n3)     = psobs
    scm_t3(:,n3)    = tobs
    scm_q3(1,:,n3)  = qobs

    if(mpi_rank()==0)then
        print*,'uobs:',uobs
        print*,'vobs:',vobs
        print*,'wfld:',wfld
        print*,'psobs:',psobs
        print*,'tobs:',tobs
        print*,'qobs:',qobs
    end if

    !-----------------LiXH Test------------
    ! cldliq and cldice in ARM95.nc has only one level 
    !scm_q3(2,:,n3)  = cldliqobs
    !scm_q3(3,:,n3)  = cldiceobs
    !-----------------LiXH Test------------

#ifdef AMIPC_PHYSICS
    do m = 1, ntracer
        where(scm_q3(m,:,n3) .lt. phy_tracer_info(m)%qmin)   &
              scm_q3(m,:,n3) = phy_tracer_info(m)%qmin
    end do
#endif

    end subroutine grist_scm_dyn_init2


    subroutine grist_scm_dyn_end2

        deallocate(scm_u)
        deallocate(scm_v)
        deallocate(scm_omega)
        deallocate(scm_ps3)
        deallocate(scm_t3)
        deallocate(scm_q3)
        deallocate(scm_t_tend)
        deallocate(scm_u_tend)
        deallocate(scm_v_tend)
        deallocate(scm_q_tend)

    end subroutine grist_scm_dyn_end2


    subroutine grist_scm_dyn_pkg2(nstep, dtime)

! io
    integer,  intent(in)  :: nstep
    real(r8), intent(in)  :: dtime
! local
    integer     :: k, m
    integer     :: rk_number, irk_step
    real(r8)    :: fac, dt
    real(r8)    :: tfcst(nlev),tfcst_rk(nlev)
    real(r8)    :: qfcst(ntracer,nlev),qfcst_rk(ntracer,nlev)
    real(r8)    :: scm_t3_face(nlevp)
    real(r8)    :: scm_q3_face(ntracer,nlevp)
    real(r8)    :: pmid(nlev), pint(nlevp), pdel(nlev), pdel_int(nlevp)
    real(r8)    :: rtau
    real(r8)    :: relaxt(nlev), relaxq(nlev)

    call time_integration_renew_mass_state(1, scm_ps3(n3), pint, pdel, pmid, pdel_int)

    tfcst_rk(1:nlev)           = scm_t3(1:nlev,n3)
    qfcst_rk(1:ntracer,1:nlev) = scm_q3(1:ntracer,1:nlev,n3)

    do irk_step = 1, nrk
        
        rk_number = nrk+1-irk_step
        dt        = dtime/rk_number

        ! update temperature
        call calc_scm_vert_flux_operator(tfcst_rk, wfldint, ver_adv_flag, scm_t3_face)

!--------------LiXH Test---------------
        do k = 1, nlev
            tfcst(k) = scm_t3(k,n3) + dt*scm_omega(k)*tfcst_rk(k)*rdry/(cp*pmid(k))      &
                                    + dt*(scm_t_tend(k)+divt(k))
        end do

!        do k = 1, nlev
!            tfcst(k) = scm_t3(k,n3) + dt*scm_omega(k)*tfcst_rk(k)*rdry/(cp*pmid(k))      &
!                                    + dt*divt(k)
!        end do
!--------------LiXH Test---------------

        do k = 2, nlev-1
            tfcst(k) = tfcst(k) - dt*scm_omega(k)/pdel(k)*(scm_t3_face(k+1)-scm_t3_face(k))
        end do
 
        k = 1
        fac = 0.5_r8*dt/pdel_int(k+1)
        tfcst(k) = tfcst(k) - fac*(wfldint(k+1)*(tfcst_rk(k+1) - tfcst_rk(k)))
        k = nlev
        fac = 0.5_r8*dt/pdel_int(k)
        tfcst(k) = tfcst(k) - fac*(wfldint(k)*(tfcst_rk(k) - tfcst_rk(k-1)))

        ! update tracer
        call calc_scm_tracer_at_face(qfcst_rk(1,:), wfldint, pdel, scm_q3_face(1,:), dt, tracer_vadv_flag)
        
!--------------LiXH Test---------------
        do m = 1, ntracer
            do k = 1, nlev
                qfcst(m,k) = scm_q3(m,k,n3) + dt*scm_q_tend(m,k)
            end do
        end do

!        do m = 1, ntracer
!            do k = 1, nlev
!                qfcst(m,k) = scm_q3(m,k,n3)
!            end do
!        end do
!--------------LiXH Test---------------

        do k = 1, nlev
            qfcst(1,k) = qfcst(1,k) + dt*divq(k)
        end do

        do k = 2, nlev-1
            qfcst(1,k) = qfcst(1,k) - dt*scm_omega(k)/pdel(k)*(scm_q3_face(1,k+1)-scm_q3_face(1,k))
        end do

        k = 1
        fac = 0.5_r8*dt/pdel_int(k+1)
        qfcst(1,k) = qfcst(1,k) - fac*(wfldint(k+1)*(qfcst_rk(1,k+1) - qfcst_rk(1,k)))
        k = nlev
        fac = 0.5_r8*dt/pdel_int(k)
        qfcst(1,k) = qfcst(1,k) - fac*(wfldint(k)*(qfcst_rk(1,k) - qfcst_rk(1,k-1)))

        tfcst_rk = tfcst
        qfcst_rk = qfcst

    end do

!--------------LiXH Test---------------
!        do k = 1, nlev
!            tfcst(k) = tfcst(k) + dtime*scm_t_tend(k)
!        end do
!
!        do m = 1, ntracer
!            do k = 1, nlev
!                qfcst(m,k) = qfcst(m,k) + dtime*scm_q_tend(m,k)
!            end do
!        end do
!--------------LiXH Test---------------

    if (scm_relaxation) then
    ! This is where we relax the solution if requested.
    ! The relaxation can be thought of as a part of the "adjustment" physics
        do k=1,nlev
            rtau       = 10800._r8          ! 3-hr adj. time scale
            rtau       = max(dtime,rtau)
            relaxt(k)  = -(tfcst(k)   - tobs(k))/rtau
            relaxq(k)  = -(qfcst(1,k) - qobs(k))/rtau
            
            tfcst(k)   = tfcst(k)   + relaxt(k)*dtime
            if (trim(scm_test_name) .eq. 'cgilsS6' .or. trim(scm_test_name) .eq. 'CGILSS6')then

            else
            qfcst(1,k) = qfcst(1,k) + relaxq(k)*dtime
            end if
        end do
    end if

#ifdef AMIPC_PHYSICS
    do m = 1, ntracer
        where(qfcst(m,:) .lt. phy_tracer_info(m)%qmin)   &
              qfcst(m,:) = phy_tracer_info(m)%qmin
    end do
#endif

    call grist_scm_io_read(nstep+1)

#ifdef AMIPW_PHYSICS
    !t_tend_dycore and q_tend_dycore is used for NTDKV381 convective parameterization.
    !Assume dt_cu = dt_SCM = dt_dycore. 
    !If use dt_cu > dt_dycore as in 3D model, modify this part to calculate Avg(t_tend) and Avg(q_tend)
    !LiXH
    t_tend_dycore(:) = (tfcst(:)-scm_t3(:,n3))/dtime
    q_tend_dycore(:) = (qfcst(1,:)-scm_q3(1,:,n3))/dtime
#endif

    scm_q3(1:ntracer,:,n3) = qfcst
    scm_t3(:,n3)           = tfcst
    scm_u       = uobs
    scm_v       = vobs
    scm_omega   = wfld
    scm_ps3(n3) = psobs

    end subroutine grist_scm_dyn_pkg2

    subroutine calc_scm_vert_flux_operator(scalar_at_full_level    ,&
                                           scalar_mass_eta_velocity,&
                                           order                   ,&
                                           scalar_at_face_level    ,&
                                           called_by_uwind)
    
    use grist_hpe_constants,  only: deta_face, deta_full

    real(r8), dimension(nlev),   intent(in)     :: scalar_at_full_level
    real(r8), dimension(nlevp),  intent(in)     :: scalar_mass_eta_velocity
    integer                   ,  intent(in)     :: order
    real(r8), dimension(nlevp),  intent(inout)  :: scalar_at_face_level
    logical, optional ,          intent(in)     :: called_by_uwind
! 
    integer(i4)                                 :: ilev
    real(r8)                                    :: part1, der_k, der_kp1
    logical                                     :: local_flag
!
! extrapolate, actually not used because etadot at boundaries are zero
! 
     scalar_at_face_level(1)        = scalar_at_full_level(1)
     scalar_at_face_level(nlev+1)   = scalar_at_full_level(nlev)

     if(present(called_by_uwind)) then
        local_flag = called_by_uwind
     else
        local_flag = .false.
     end if

     IF(eqs_vert_diff.or.local_flag)THEN
!
! This is equivalent distance version, old version used for simplicity
! only used for bit reproducity
!
     select case(order)
     case(2) 
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 2, nlev
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+&
                                           scalar_at_full_level(ilev-1))*0.5_r8
        end do

     case(3) 
        scalar_at_face_level(1)       = 0._r8       !Do not use, LiXH
        scalar_at_face_level(2)       = (scalar_at_full_level(2)+scalar_at_full_level(1))*0.5_r8
        scalar_at_face_level(nlev)    = (scalar_at_full_level(nlev)+scalar_at_full_level(nlev-1))*0.5_r8
        scalar_at_face_level(nlev+1)  = 0._r8       !Do not use, LiXH
        do ilev = 3, nlev-1
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+scalar_at_full_level(ilev-1))*(7._r8/12._r8)-&
                                          (scalar_at_full_level(ilev+1)+scalar_at_full_level(ilev-2))*(1._r8/12._r8)+&
                                sign(1._r8,scalar_mass_eta_velocity(ilev))*&
                                         ((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev-2))-&
                                    3._r8*(scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1)))/12._r8
        end do

     case default
        print*,"you must select a vertical order, stop"
        stop
     end select

     ELSE  ! non-equivalent distance version

     select case(order)
     case(2)
        do ilev = 2, nlev
           scalar_at_face_level(ilev)= (scalar_at_full_level(ilev)*deta_full(ilev-1)+&
                                          scalar_at_full_level(ilev-1)*deta_full(ilev))*0.5_r8/deta_face(ilev)
        end do
     case(3) 
        scalar_at_face_level(2)      = (scalar_at_full_level(2)*deta_full(1)+scalar_at_full_level(1)*deta_full(2))*0.5_r8/deta_face(2)
        scalar_at_face_level(nlev)   = (scalar_at_full_level(nlev)*deta_full(nlev-1)+scalar_at_full_level(nlev-1)*deta_full(nlev))*0.5_r8/deta_face(nlev)
        do ilev = 3, nlev-1 ! face level
! face level
           part1 = (scalar_at_full_level(ilev)*deta_full(ilev-1)+&
                    scalar_at_full_level(ilev-1)*deta_full(ilev))*0.5_r8/deta_face(ilev)
! full level for this face level
           der_k  = 2._r8*(deta_full(ilev-1)**2)*((scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))/deta_face(ilev)-&
                                                 (scalar_at_full_level(ilev-1)-scalar_at_full_level(ilev-2))/deta_face(ilev-1))/&
                                                 (deta_face(ilev)+deta_face(ilev-1))
           der_kp1= 2._r8*(deta_full(ilev)**2)*((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev))/deta_face(ilev+1)-&
                                                 (scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))/deta_face(ilev))/&
                                                 (deta_face(ilev+1)+deta_face(ilev))
           scalar_at_face_level(ilev)= part1-(der_k+der_kp1)/12._r8+sign(1._r8,scalar_mass_eta_velocity(ilev))*(der_kp1-der_k)
        end do
     case default
        print*, "you must select a vertical order, stop"
        stop
     end select

     END IF

     return
 
    end subroutine calc_scm_vert_flux_operator

    subroutine calc_scm_tracer_at_face(scalar_at_full_level      ,&
                                       scalar_mass_eta_velocity  ,&
                                       scalar_delhp_at_full_level,&
                                       scalar_at_face_level      ,&
                                       dt, order )

! io
    real(r8), dimension(nlev), intent(in)    :: scalar_at_full_level
    real(r8), dimension(nlevp),intent(in)    :: scalar_mass_eta_velocity
    real(r8), dimension(nlev), intent(in)    :: scalar_delhp_at_full_level
    real(r8), dimension(nlevp),intent(inout) :: scalar_at_face_level
    real(r8)                 , intent(in)    :: dt
    integer(i4)              , intent(in)    :: order
! local 
    real(r8)                 :: scalar_delhp_at_face_level(nlev+1)
    real(r8)                 :: part1, der_k, der_kp1
    real(r8)                 :: qhat(nlev+1), cr_num
    integer(i4)              :: ilev,flag
!
! extrapolate, actually not used because etadot at boundaries are zero
! 
     scalar_at_face_level(1)        = scalar_at_full_level(1)
     scalar_at_face_level(nlev+1)   = scalar_at_full_level(nlev)
!
! compute delhp at face
!
     do ilev = 2, nlev
        scalar_delhp_at_face_level(ilev) = 0.5_r8*(scalar_delhp_at_full_level(ilev-1)+scalar_delhp_at_full_level(ilev))
     end do
     scalar_delhp_at_face_level(1)       = 0.5_r8*scalar_delhp_at_full_level(1)
     scalar_delhp_at_face_level(nlev+1)  = 0.5_r8*scalar_delhp_at_full_level(nlev)
!
! This is equivalent distance version, old version used for simplicity
! only used for bit reproducity
!
     select case(order)

     case(1)
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 2, nlev
           scalar_at_face_level(ilev) = 0.5_r8*(scalar_at_full_level(ilev)+scalar_at_full_level(ilev-1))-&
                                        0.5_r8*sign(1._r8,scalar_mass_eta_velocity(ilev))*&
                                        (scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))
        end do
     case(2)
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 2, nlev
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+&
                                         scalar_at_full_level(ilev-1))*0.5_r8
        end do

     case(3) 
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(2)       = (scalar_at_full_level(2)   +scalar_at_full_level(1))*0.5_r8
        scalar_at_face_level(nlev)    = (scalar_at_full_level(nlev)+scalar_at_full_level(nlev-1))*0.5_r8
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 3, nlev-1
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+scalar_at_full_level(ilev-1))*(7._r8/12._r8)-&
                                        (scalar_at_full_level(ilev+1)+scalar_at_full_level(ilev-2))*(1._r8/12._r8)+&
                                sign(1._r8,scalar_mass_eta_velocity(ilev))*&
                                        ((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev-2))-&
                                   3._r8*(scalar_at_full_level(ilev)  -scalar_at_full_level(ilev-1)))/12._r8
        end do

     case(9) !ppm
! set bdy
        qhat(1)       = 0._r8
        qhat(2)       = (scalar_at_full_level(2)+scalar_at_full_level(1))*0.5_r8
        qhat(nlev)    = (scalar_at_full_level(nlev)+scalar_at_full_level(nlev-1))*0.5_r8
        qhat(nlev+1)  = 0._r8
! compute qhat
        do ilev = 3, nlev-1
           qhat(ilev) = (7._r8*(scalar_at_full_level(ilev-1)+scalar_at_full_level(ilev))-&
                               (scalar_at_full_level(ilev-2)+scalar_at_full_level(ilev+1)))/12._r8
        end do
! compute qface
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 2, nlev
           cr_num = scalar_mass_eta_velocity(ilev)*dt/scalar_delhp_at_face_level(ilev)
           flag   = idint((cr_num-abs(cr_num))/(2*cr_num+1.e-80_r8))
           if(flag.ne.0.and.flag.ne.1)then
              print*,"flag in vadv 9 is wrong, stop"
           end if
           scalar_at_face_level(ilev) = qhat(ilev)-cr_num*(qhat(ilev)-scalar_at_full_level(ilev-1+flag))-&
                                                  abs(cr_num)*(1._r8-abs(cr_num))*&
                                                 (qhat(ilev-1+flag)-2._r8*scalar_at_full_level(ilev-1+flag)+qhat(ilev+flag))
        end do
     case default
        print*,"you must select a vertical order, stop"
        stop
     end select

     return
 
    end subroutine calc_scm_tracer_at_face

 end module grist_scm_dyn_core2
