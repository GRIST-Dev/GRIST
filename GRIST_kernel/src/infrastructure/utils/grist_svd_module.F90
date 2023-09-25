
   module grist_svd_module
   !use grist_constants, only: i4, r8
   use grist_constants_dbl, only: i4, r8
!============================================================================
!    This SVD routine is adopted online:
!    https://people.sc.fsu.edu/~jburkardt/f_src/svd_demo/svd_demo.html
!
!    Command parameter, integer ( i4 ) M, N, the number of rows and
!    columns of the matrix.  If M or N is not supplied on the command line,
!    the user is prompted to supply them.
!
!    Parameters:
!
!    Local, real ( r8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Local, real ( r8 ) S(M,N), the diagonal factor
!    in the singular value decomposition of A.
!
!    Output, real ( r8 ) U(M,M), the first orthogonal factor
!    in the singular value decomposition of A.
!
!    Output, real ( r8 ) V(N,N), the second orthogonal factor
!    in the singular value decomposition of A.
!=============================================================================

  implicit none

  private
  public  :: pseudo_inverse,&
             svd_lapack, migs

  contains

  subroutine pseudo_inverse( m, n, u, s, v, a_pseudo )
!===============================================================================
!
!! PSEUDO_INVERSE computes the pseudoinverse.
!
!  Discussion:
!
!    Given the singular value decomposition of a real MxN matrix A:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    the pseudo inverse is the NxM matrix A+ with the form
!
!      A+ = V * S+ * U'
!
!    where
!
!      S+ is the NxM matrix whose nonzero diagonal elements are
!      the inverses of the corresponding diagonal elements of S.
!
!  Parameters:
!
!    Input, integer ( i4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( r8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Input, real ( r8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
!    Output, real ( r8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
!
!===============================================================================

    implicit none
! io
    integer(i4), intent(in)   :: m
    integer(i4), intent(in)   :: n
    real(r8),    intent(in)   :: u(m,m)
    real(r8),    intent(in)   :: s(m,n)
    real(r8),    intent(in)   :: v(n,n)
    real(r8),    intent(inout):: a_pseudo(n,m)
! local
    integer(i4)  :: i
    real(r8)     :: sp(n,m)

    sp(1:n,1:m) = 0._r8
    do i = 1, min ( m, n )
      if ( s(i,i) /= 0._r8 ) then
        sp(i,i) = 1._r8 / s(i,i)
      end if
    end do

    a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n), &
    matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )

   return
   end subroutine pseudo_inverse

   subroutine svd_lapack ( m, n, array, u, s, v )
!===============================================================================
!
!! R8MAT_SVD_LAPACK gets the SVD of a matrix using a call to LAPACK.
!
!  Discussion:
!
!    The singular value decomposition of a real MxN matrix A has the form:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    Moreover, the nonzero entries of S are positive, and appear
!    in order, from largest magnitude to smallest.
!
!    This routine calls the LAPACK routine DGESVD to compute the
!    factorization.
!
!  Parameters:
!
!    Input, integer ( i4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( r8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Output, real ( r8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
!===============================================================================

   implicit none
! io
   integer(i4), intent(in)    :: m
   integer(i4), intent(in)    :: n
   real(r8),    intent(in)    :: array(m,n)
   real(r8),    intent(inout) :: u(m,m)
   real(r8),    intent(inout) :: s(m,n)
   real(r8),    intent(inout) :: v(n,n)
! local
   real(r8)  :: a_copy(m,n)
   real(r8)  :: sdiag(min(m,n))
   real(r8), allocatable, dimension ( : ) :: work
   integer(i4) :: i
   integer(i4) :: info
   integer(i4) :: lda
   integer(i4) :: ldu
   integer(i4) :: ldv
   integer(i4) :: lwork
   character   :: jobu
   character   :: jobv

   lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

   allocate ( work(1:lwork) )
!
!  Compute the eigenvalues and eigenvectors.
!
   jobu = 'A'
   jobv = 'A'
   lda = m
   ldu = m
   ldv = n
!
!  The input matrix is destroyed by the routine.  Since we need to keep
!  it around, we only pass a copy to the routine.
!
   a_copy(:,:) = array(:,:)

   if (r8 .eq. 8) then
     call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
                   lwork, info )
   elseif (r8 .eq. 4) then
     call sgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
                   lwork, info )
   end if

   if ( info /= 0 ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a)' ) 'R8MAT_SVD_LAPACK - Failure!'
     write ( *, '(a)' ) '  The SVD could not be calculated.'
     write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
     write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
     return
   end if
!
!  Make the MxN matrix S from the diagonal values in SDIAG.
!
   s(1:m,1:n) = 0.0D+00
   do i = 1, min ( m, n )
     s(i,i) = sdiag(i)
   end do
!
!  Transpose V.
!
   v = transpose ( v )

   deallocate ( work )

   return
  end subroutine svd_lapack

   !
   ! check the inverse routine used by MPAS
   ! we still use SVD by default
   ! Updated 10/24/2001.
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!   Program 4.4   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !                                                                       !
   ! Please Note:                                                          !
   !                                                                       !
   ! (1) This computer program is written by Tao Pang in conjunction with  !
   !     his book, "An Introduction to Computational Physics," published   !
   !     by Cambridge University Press in 1997.                            !
   !                                                                       !
   ! (2) No warranties, express or implied, are made for this program.     !
   !                                                                       !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   
   SUBROUTINE MIGS (A,N,X,INDX)
   !
   ! Subroutine to invert matrix A(N,N) with the inverse stored
   ! in X(N,N) in the output.  Copyright (c) Tao Pang 2001.
   !
     IMPLICIT NONE
     INTEGER, INTENT (IN) :: N
     INTEGER :: I,J,K
     INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
     REAL (r8), INTENT (INOUT), DIMENSION (N,N):: A
     REAL (r8), INTENT (OUT), DIMENSION (N,N):: X
     REAL (r8), DIMENSION (N,N) :: B
   !
     DO I = 1, N
       DO J = 1, N
         B(I,J) = 0.0
       END DO
     END DO
     DO I = 1, N
       B(I,I) = 1.0
     END DO
   !
     CALL ELGS (A,N,INDX)
   !
     DO I = 1, N-1
       DO J = I+1, N
         DO K = 1, N
           B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
         END DO
       END DO
     END DO
   !
     DO I = 1, N
       X(N,I) = B(INDX(N),I)/A(INDX(N),N)
       DO J = N-1, 1, -1
         X(J,I) = B(INDX(J),I)
         DO K = J+1, N
           X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
         END DO
         X(J,I) =  X(J,I)/A(INDX(J),J)
       END DO
     END DO
   END SUBROUTINE MIGS

   SUBROUTINE ELGS (A,N,INDX)
   !
   ! Subroutine to perform the partial-pivoting Gaussian elimination.
   ! A(N,N) is the original matrix in the input and transformed matrix
   ! plus the pivoting element ratios below the diagonal in the output.
   ! INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.
   !
     IMPLICIT NONE
     INTEGER, INTENT (IN) :: N
     INTEGER :: I,J,K,ITMP
     INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
     REAL (r8) :: C1,PI,PI1,PJ
     REAL (r8), INTENT (INOUT), DIMENSION (N,N) :: A
     REAL (r8), DIMENSION (N) :: C
   !
   ! Initialize the index
   !
     DO I = 1, N
       INDX(I) = I
     END DO
   !
   ! Find the rescaling factors, one from each row
   !
     DO I = 1, N
       C1= 0.0
       DO J = 1, N
         C1 = MAX(C1,ABS(A(I,J)))
       END DO
       C(I) = C1
     END DO
   !
   ! Search the pivoting (largest) element from each column
   !
     DO J = 1, N-1
       PI1 = 0.0
       DO I = J, N
         PI = ABS(A(INDX(I),J))/C(INDX(I))
         IF (PI.GT.PI1) THEN
           PI1 = PI
           K   = I
         ENDIF
       END DO
   !
   ! Interchange the rows via INDX(N) to record pivoting order
   !
       ITMP    = INDX(J)
       INDX(J) = INDX(K)
       INDX(K) = ITMP
       DO I = J+1, N
         PJ  = A(INDX(I),J)/A(INDX(J),J)
   !
   ! Record pivoting ratios below the diagonal
   !
         A(INDX(I),J) = PJ
   !
   ! Modify other elements accordingly
   !
         DO K = J+1, N
           A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
         END DO
       END DO
     END DO

   END SUBROUTINE ELGS
  end module grist_svd_module
