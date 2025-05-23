
!----------------------------------------------------------------------------
! Created on 2018
! Author: Yi Zhang
! Version 1.0
! Description: A tridiagnoal solver for grist; https://www.quantstart.com/
!              articles/Tridiagonal-Matrix-Solver-via-Thomas-Algorithm
! Revision history:
!----------------------------------------------------------------------------

 module grist_tridiagonal_solver_mixed
    use grist_constants,      only: i4, r8

   implicit none
   public :: solve_tridiagonal
  contains

  subroutine solve_tridiagonal(a_coef,b_coef,c_coef,r_coef,x_value,nsize)
! io
    real(r8),  intent(in)   :: a_coef(nsize)
    real(r8),  intent(in)   :: b_coef(nsize)
    real(r8),  intent(in)   :: c_coef(nsize)
    real(r8),  intent(in)   :: r_coef(nsize)
    real(r8),  intent(inout):: x_value(nsize)
    integer(i4), intent(in) :: nsize
! local
    real(r8)                :: gama(nsize)
    real(r8)                :: rho(nsize)
    integer(i4)             :: idx
!
! compute gama
!
    gama(1) = c_coef(1)/b_coef(1)
    do idx = 2, nsize-1 
        gama(idx) = c_coef(idx)/(b_coef(idx)-a_coef(idx)*gama(idx-1))
    end do
    gama(nsize) = 0.

    rho(1)  = r_coef(1)/b_coef(1)
    do idx  = 2, nsize
       rho(idx) = (r_coef(idx)-a_coef(idx)*rho(idx-1))/(b_coef(idx)-a_coef(idx)*gama(idx-1))
    end do

    x_value(nsize) = rho(nsize)
    do idx = nsize-1, 1, -1
       x_value(idx) = rho(idx)-gama(idx)*x_value(idx+1)
    end do

   return
  end subroutine solve_tridiagonal

 end module grist_tridiagonal_solver_mixed 
