!===================================================================================================
!
!  Created by LiXiaohan on 19/07/17, adopted from CAM5
!
!     interface of the GFDL style dry adiabatic adjustment
!
!     Method:
!     if stratification is unstable, adjustment to the dry adiabatic lapse
!     rate is forced subject to the condition that enthalpy is conserved.

!===================================================================================================

 module grist_dry_adjustment

    use grist_handle_error,                 only: endrun
    use grist_nml_module,                   only: ntracer, nlev, nlevp
    use grist_constants,                    only: r8, i4, cappa, cp
    use grist_mpi

    implicit none
    private
    public          :: dadadj

! private
    integer, parameter ::  niter  = 15
    integer, parameter ::  nlvdry = 3

 contains

    subroutine dadadj(ncol, dt)
     use grist_physics_data_structure,       only: pstate
     use grist_PhysC_data_structure,         only: ptend_dry_adjustment
! io
    integer, intent(in)  :: ncol
    real(r8), intent(in) :: dt 
! local
    integer i,k             ! longitude, level indices
    integer jiter           ! iteration index

    real(r8) c1dad(nlev)        ! intermediate constant
    real(r8) c2dad(nlev)        ! intermediate constant
    real(r8) c3dad(nlev)        ! intermediate constant
    real(r8) c4dad(nlev)        ! intermediate constant
    real(r8) gammad             ! dry adiabatic lapse rate (deg/Pa)
    real(r8) zeps               ! convergence criterion (deg/Pa)
    real(r8) rdenom             ! reciprocal of denominator of expression
    real(r8) dtdp               ! delta-t/delta-p
    real(r8) zepsdp             ! zeps*delta-p
    real(r8) zgamma             ! intermediate constant
    real(r8) qave               ! mean q between levels

    logical ilconv          ! .TRUE. ==> convergence was attained
    logical dodad(ncol)    ! .TRUE. ==> do dry adjustment

    real(r8)                 :: t(nlev, ncol)
    real(r8)                 :: q(nlev, ncol)
    real(r8)                 :: pdel(nlev, ncol)
    real(r8)                 :: pmid(nlev, ncol)
    real(r8)                 :: pint(nlevp, ncol)
 
    t(:,1:ncol)        = pstate%temp_at_pc_full_level%f(:,1:ncol)
    q(:,1:ncol)        = pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol)
    pmid(:,1:ncol)     = pstate%pressure_at_pc_full_level%f(:,1:ncol)
    pint(:,1:ncol)     = pstate%pressure_at_pc_face_level%f(:,1:ncol)
    pdel(:,1:ncol)     = pstate%delp_at_pc_full_level%f(:,1:ncol)
 
    zeps = 2.0e-5_r8           ! set convergence criteria

!  Find gridpoints with unstable stratification
    do i=1,ncol
       gammad = cappa*0.5_r8*(t(2,i) + t(1,i))/pint(2,i)
       dtdp = (t(2,i) - t(1,i))/(pmid(2,i) - pmid(1,i))
       dodad(i) = (dtdp + zeps) .gt. gammad
    end do

    do k=2,nlvdry
       do i=1,ncol
          gammad = cappa*0.5_r8*(t(k+1,i) + t(k,i))/pint(k+1,i)
          dtdp = (t(k+1,i) - t(k,i))/(pmid(k+1,i) - pmid(k,i))
          dodad(i) = dodad(i) .or. (dtdp + zeps).gt.gammad
       end do
    end do

!  Make a dry adiabatic adjustment
!  Note: nlvdry ****MUST**** be < pver

    do 80 i=1,ncol
       if (dodad(i)) then
!--------------LiXH Test--------------
if(mpi_rank()==0)then
    print*,'do the dry adjustment'
end if
!--------------LiXH Test--------------
          zeps = 2.0e-5_r8
          do k=1,nlvdry
             c1dad(k) = cappa*0.5_r8*(pmid(k+1,i)-pmid(k,i))/pint(k+1,i)
             c2dad(k) = (1._r8 - c1dad(k))/(1._r8 + c1dad(k))
             rdenom = 1._r8/(pdel(k,i)*c2dad(k) + pdel(k+1,i))
             c3dad(k) = rdenom*pdel(k,i)
             c4dad(k) = rdenom*pdel(k+1,i)
          end do
50        do jiter=1,niter
             ilconv = .true.
             do k=1,nlvdry
                zepsdp = zeps*(pmid(k+1,i) - pmid(k,i))
                zgamma = c1dad(k)*(t(k,i) + t(k+1,i))
                if ((t(k+1,i)-t(k,i)) >= (zgamma+zepsdp)) then
                   ilconv = .false.
                   t(k+1,i) = t(k,i)*c3dad(k) + t(k+1,i)*c4dad(k)
                   t(k,i) = c2dad(k)*t(k+1,i)
                   qave = (pdel(k+1,i)*q(k+1,i) + pdel(k,i)*q(k,i))/(pdel(k+1,i)+ pdel(k,i))
                   q(k+1,i) = qave
                   q(k,i) = qave
                end if
             end do
             if (ilconv) go to 80 ! convergence => next longitude
          end do

!  Double convergence criterion if no convergence in niter iterations
          zeps = zeps + zeps
          if (zeps > 1.e-4_r8) then
             print*,'DADADJ: No convergence in dry adiabatic adjustment'
             print*,'Rank=',mpi_rank()
             call endrun("grist dry sdjustment")
          else
             print*,'DADADJ: Convergence criterion doubled to EPS=', zeps
             go to 50
          end if
       end if
80  continue

    ! re-initialize ptend
    ptend_dry_adjustment%tend_s%f = 0._r8
    ptend_dry_adjustment%tend_q%f = 0._r8
 
    ptend_dry_adjustment%tend_s%f(:,1:ncol)   = (t(:,1:ncol)-pstate%temp_at_pc_full_level%f(:,1:ncol))/dt*cp
    ptend_dry_adjustment%tend_q%f(1,:,1:ncol) = (q(:,1:ncol)-pstate%tracer_mxrt_at_pc_full_level%f(1,:,1:ncol))/dt

       return

    end subroutine dadadj
 end module grist_dry_adjustment
