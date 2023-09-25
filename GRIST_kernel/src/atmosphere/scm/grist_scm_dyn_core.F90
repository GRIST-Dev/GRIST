!======================================================
!
!  Created by LiXiaohan on 19/5/20.
!  Adopted from CAM_EUL
!  
!  Purpose:
!  SCM dynamical core
!  Note: Remove dyn vars to an isolated module
!======================================================

 module grist_scm_dyn_core

    use grist_constants,                    only: i4, r8, rdry, cp
    use grist_nml_module,                   only: nlev, nlevp, ntracer
    use grist_scm_pressure_update,          only: time_integration_renew_mass_state
    use grist_scm_io_module
    use grist_physics_data_structure,       only: phy_tracer_info
    use grist_scm_comm_module
    use grist_handle_error,                 only: endrun
    use grist_mpi
    use grist_scm_dyn_vars ! all

    implicit none

    private

    public  :: grist_scm_dyn_init,  &
               grist_scm_dyn_end,   &
               grist_scm_dyn_pkg

 contains
 
    subroutine grist_scm_dyn_init

! local
    integer  :: m,k
    real(r8) :: pmid(nlev), pdel(nlev)
    real(r8) :: pint(nlevp), pdel_int(nlevp)

    n3m2 = 1
    n3m1 = 2
    n3   = 3
    allocate(scm_u(nlev));            scm_u      = 0._r8
    allocate(scm_v(nlev));            scm_v      = 0._r8
    allocate(scm_omega(nlev));        scm_omega  = 0._r8
    allocate(scm_ps3(3));             scm_ps3    = 0._r8
    allocate(scm_t3(nlev,3));         scm_t3     = 0._r8
    allocate(scm_q3(ntracer,nlev,3)); scm_q3     = 0._r8

    allocate(scm_t_tend(nlev));       scm_t_tend = 0._r8
    allocate(scm_u_tend(nlev));       scm_u_tend = 0._r8
    allocate(scm_v_tend(nlev));       scm_v_tend = 0._r8
    allocate(scm_qminus(ntracer,nlev)); scm_qminus = 0._r8

    call scam_data_io_read('/g13/pengxd/lixh/cesm1_2_2/SCMtest/run/O3_set0/dycoms.cam.h1.1999-07-11-00000.nc')

    call grist_scm_io_read(0)

    scm_u(:)        = uobs
    scm_v(:)        = vobs
    scm_omega(:)    = wfld
    scm_ps3(n3)     = psobs
    scm_t3(:,n3)    = tobs
    scm_q3(1,:,n3)  = qobs
    !-----------------LiXH Test------------
    ! cldliq and cldice in ARM95.nc has only one level 
    !scm_q3(2,:,n3)  = cldliqobs
    !scm_q3(3,:,n3)  = cldiceobs
    !-----------------LiXH Test------------

    !--------------shoule be modulated to neg3 in CAM, LiXH--------------
    do m = 1, ntracer
        where(scm_q3(m,:,n3) .lt. phy_tracer_info(m)%qmin)   &
              scm_q3(m,:,n3) = phy_tracer_info(m)%qmin
    end do
    !--------------shoule be modulated to neg3 in CAM, LiXH--------------


    ! copytimelevels
    scm_ps3(n3m2)    = scm_ps3(n3)
    scm_t3(:,n3m2)   = scm_t3(:,n3)
    scm_q3(1:ntracer,:,n3m2) = scm_q3(1:ntracer,:,n3)

    scm_ps3(n3m1)    = scm_ps3(n3)
    scm_t3(:,n3m1)   = scm_t3(:,n3)
    scm_q3(1:ntracer,:,n3m1) = scm_q3(1:ntracer,:,n3)

    ! for initial output
    !scalar_hpressure_at_pc_surface_n%f(1)          = scm_ps3(n3m1)
    !scalar_temp_at_pc_full_level_n%f(:,1)          = scm_t3(:,n3m1)
    !scalar_U_wind_at_pc_full_level_n%f(:,1)        = scm_u
    !scalar_V_wind_at_pc_full_level_n%f(:,1)        = scm_v
    !scalar_www_at_pc_full_level_n%f(:,1)           = scm_omega
    !scalar_tracer_mxrt_at_pc_full_level_n%f(1:ntracer,:,1) = scm_q3(1:ntracer,:,n3m1)

    allocate(etamid(nlev));           etamid  = 0._r8
    allocate(etaint(nlevp));          etaint  = 0._r8

    call time_integration_renew_mass_state(1, scm_ps3(n3), pint, pdel, pmid, pdel_int)
    etamid = pmid
    etaint = pint

    end subroutine grist_scm_dyn_init


    subroutine grist_scm_dyn_end

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
        deallocate(scm_qminus)

        deallocate(etamid)
        deallocate(etaint)
    end subroutine grist_scm_dyn_end


    subroutine grist_scm_dyn_pkg(nstep, dtime)
    use grist_scm_scanslt,       only: scanslt_run
! io
    integer,  intent(in)  :: nstep
    real(r8), intent(in)  :: dtime
! local
    integer     :: k, itmp, m
    real(r8)    :: om2eps
    real(r8)    :: etadot(nlevp)
    real(r8)    :: detam(nlev)
    real(r8)    :: qfcst(ntracer,nlev)


    call grist_scm_io_read(nstep)

    etadot(:) = wfldint(:)
    detam     = 0._r8
    qfcst     = 0._r8

    call scanslt_run(nstep, scm_qminus, dtime,  etadot,     &
                     detam, etamid,     etaint, qfcst)
 
    call forecast(scm_ps3(n3m1) , scm_ps3(n3m2)   , scm_ps3(n3)     ,   &
                  scm_u         , scm_v           ,                     &
                  scm_t3(:,n3)  , scm_t3(:,n3m1)  , scm_t3(:,n3m2)  ,   &
                  scm_q3(:,:,n3), scm_q3(:,:,n3m1), scm_q3(:,:,n3m2),   &
                  dtime         , scm_t_tend(:)   , scm_u_tend(:)   ,   &
                  scm_v_tend(:) , scm_omega(:)    ,                     &
                  scm_qminus    , qfcst(1:ntracer,:)  )

    !--------------shoule be modulated to neg3 in CAM, LiXH--------------
    do m = 1, ntracer
        where(scm_q3(m,:,n3) .lt. phy_tracer_info(m)%qmin)   &
              scm_q3(m,:,n3) = phy_tracer_info(m)%qmin
    end do
    !--------------shoule be modulated to neg3 in CAM, LiXH--------------

    !call tfilt_massfixrun in CAM:
    if(nstep.ge.2) then
        om2eps = 1._r8 - 2._r8*eps

        scm_t3(:,n3m1) = om2eps*scm_t3(:,n3m1)+eps*scm_t3(:,n3m2)+eps*scm_t3(:,n3)
        do m = 1,ntracer
            scm_q3(m,:,n3m1) = om2eps*scm_q3(m,:,n3m1)+eps*scm_q3(m,:,n3m2)+eps*scm_q3(m,:,n3)
        end do
        scm_ps3(n3m1)  = om2eps*scm_ps3(n3m1) + eps*scm_ps3(n3m2) + eps*scm_ps3(n3)
    end if

    itmp = n3m2
    n3m2 = n3m1
    n3m1 = n3
    n3   = itmp

    end subroutine grist_scm_dyn_pkg


    subroutine forecast(psm1, psm2, ps  , u3,  v3  ,        &
                        t3  , t3m1, t3m2, q3,  q3m1, q3m2,  &
                        dt  , t2  , fu  , fv,  omega,       &
                        qminus    , qfcst)
! io
    real(r8), intent(in) :: dt
    real(r8), intent(in) :: t2(nlev)                ! temperature tendency in physics
    real(r8), intent(in) :: fu(nlev)                ! u wind tendency in physics
    real(r8), intent(in) :: fv(nlev)                ! v wind tendency in physics
    real(r8), intent(in) :: ps(1)
    real(r8), intent(in) :: psm1(1)
    real(r8), intent(in) :: psm2(1)
    real(r8), intent(in) :: t3m1(nlev)
    real(r8), intent(in) :: t3m2(nlev)
    real(r8), intent(in) :: q3m1(ntracer,nlev)
    real(r8), intent(in) :: q3m2(ntracer,nlev)
    real(r8), intent(in) :: qminus(ntracer,nlev)
    real(r8), intent(inout) :: qfcst(ntracer,nlev)
    real(r8), intent(inout) :: q3(ntracer,nlev)
    real(r8), intent(out) :: u3(nlev)
    real(r8), intent(out) :: v3(nlev)
    real(r8), intent(out) :: t3(nlev)
    real(r8), intent(out) :: omega(nlev)
! local
    integer  :: k, m
    real(r8) :: weight, fac
    real(r8) :: wfld_int(nlevp)
    real(r8) :: tfcst(nlev)
    real(r8) :: pmid_m1(nlev), pint_m1(nlevp), pdel_m1(nlev), pdel_int_m1(nlevp)
    real(r8) :: rtau
    real(r8) :: relaxt(nlev), relaxq(nlev)


    call time_integration_renew_mass_state(1, psm1, pint_m1, pdel_m1, pmid_m1, pdel_int_m1)

    wfld_int(1)     = 0._r8
    wfld_int(nlevp) = 0._r8
    do k = 2, nlev
        weight = (pint_m1(k) - pmid_m1(k-1))/(pmid_m1(k) - pmid_m1(k-1))
        wfld_int(k) = (1.0_r8 - weight)*wfld(k-1) + weight*wfld(k)
    end do

    do k = 2, nlev-1
        fac = 0.5_r8*dt/pdel_m1(k)
        tfcst(k) = t3m2(k)-fac*(wfld_int(k+1)*(t3m1(k+1) - t3m1(k))+wfld_int(k)*(t3m1(k) - t3m1(k-1)))
    end do
    k = 1
    fac = 0.5_r8*dt/pdel_m1(k)
    tfcst(k) = t3m2(k) - fac*(wfld_int(k+1)*(t3m1(k+1) - t3m1(k)))
    k = nlev
    fac = 0.5_r8*dt/pdel_m1(k)
    tfcst(k) = t3m2(k) - fac*(wfld_int(k)*(t3m1(k) - t3m1(k-1)))
    
    do k = 1, nlev
        tfcst(k) = tfcst(k) + dt*wfld(k)*t3m1(k)*rdry/(cp*pmid_m1(k))+dt*(t2(k)+divt(k))

        qfcst(1,k) = qfcst(1,k)+dt*divq(k)
    end do

    if (scm_relaxation) then
    ! This is where we relax the solution if requested.
    ! The relaxation can be thought of as a part of the "adjustment" physics
        do k=1,nlev
            rtau       = 10800._r8          ! 3-hr adj. time scale
            rtau       = max(dt,rtau)
            relaxt(k)  = -(tfcst(k)   - tobs(k))/rtau
            relaxq(k)  = -(qfcst(1,k) - qobs(k))/rtau
            
            tfcst(k)   = tfcst(k)   + relaxt(k)*dt
    !        qfcst(1,k) = qfcst(1,k) + relaxq(k)*dt
        end do
    end if


    u3 = uobs
    v3 = vobs
    t3 = tfcst
    q3(1:ntracer,:) = qfcst(1:ntracer,:)

    omega = wfld
    end subroutine forecast

 end module grist_scm_dyn_core
