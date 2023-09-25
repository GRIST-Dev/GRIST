!======================================================
!
!  Created by LiXiaohan on 19/5/20.
!  Adopted from CAM_EUL
!  
!  Purpose:
!  SCM dynamical core's variables
!======================================================

 module grist_scm_dyn_vars

    use grist_constants,                    only: r8

    implicit none

    private

    public  :: n3, n3m1, n3m2, scm_omega,               &
               scm_u, scm_v, scm_ps3, scm_t3, scm_q3,   &
               scm_t_tend, scm_u_tend, scm_v_tend,      &
               scm_q_tend, scm_qminus,                  &
               etamid, etaint, eps

    integer :: n3, n3m1, n3m2
    real(r8), allocatable :: etamid(:)
    real(r8), allocatable :: etaint(:)

    real(r8), allocatable :: scm_u(:)
    real(r8), allocatable :: scm_v(:)
    real(r8), allocatable :: scm_omega(:)
    real(r8), allocatable :: scm_ps3(:)
    real(r8), allocatable :: scm_t3(:,:)
    real(r8), allocatable :: scm_q3(:,:,:)
 
    real(r8), allocatable :: scm_t_tend(:)
    real(r8), allocatable :: scm_u_tend(:)
    real(r8), allocatable :: scm_v_tend(:)
    real(r8), allocatable :: scm_q_tend(:,:)
    real(r8), allocatable :: scm_qminus(:,:)

    real(r8), parameter   :: eps = 0.06_r8    !time filter coefficient. Defaults to 0.06.

 end module grist_scm_dyn_vars
