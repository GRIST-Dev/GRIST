
!----------------------------------------------------------------------------
! Version 1.0
! Description: This module contains some math functions collected from public 
!              domain/self-added
! Revision history: 
!              (1) use doulbe precision for all geometric calculation
!----------------------------------------------------------------------------
  module grist_math_module

   !use grist_constants ,     only: r8, i4, pi, eps
   use grist_constants_dbl , only: r4, r8, i4, pi, eps
   use grist_domain_types  , only: global_domain
   use grist_handle_error  , only: endrun
   use grist_mpi

   implicit none

    interface norm 
#ifdef SPCODE
        module procedure norm_r4
#else
        module procedure norm_r4, norm_r8
#endif
    end interface norm 
    interface cross_product 
#ifdef SPCODE
        module procedure cross_product_r4
#else
        module procedure cross_product_r4, cross_product_r8
#endif
    end interface cross_product 
    interface arcdistll
#ifdef SPCODE
        module procedure arcdistll_r4
#else
        module procedure arcdistll_r4, arcdistll_r8
#endif
    end interface arcdistll
    interface convert_vector_sph2cart
#ifdef SPCODE
        module procedure convert_vector_sph2cart_r4
#else
        module procedure convert_vector_sph2cart_r4, convert_vector_sph2cart_r8
#endif
    end interface convert_vector_sph2cart
    interface convert_vector_cart2sph
#ifdef SPCODE
        module procedure convert_vector_cart2sph_r4
#else
        module procedure convert_vector_cart2sph_r4, convert_vector_cart2sph_r8
#endif
    end interface convert_vector_cart2sph

   public

   contains

!================================================
! Given values at three consecutive full levels,
! extrapolate upper/lower boundary value based
! on a quadratic reconstruction
!================================================

    real(r8) function extrapolate_bdy(q1, q2, q3, eta1, eta2, eta3, eta_bdy)
!io
       real(r8),  intent(in)    :: q1, q2, q3
       real(r8),  intent(in)    :: eta1, eta2, eta3
       real(r8),  intent(in)    :: eta_bdy
! local
       real(r8)  :: a0, a1, a2

       a2 = ((q3-q2)/(eta3-eta2)-(q2-q1)/(eta2-eta1))/(eta3-eta1)
       a1 = (q2-q1)/(eta2-eta1)-a2*(eta2+eta1)
       a0 = q1-a1*eta1-a2*eta1**2
       extrapolate_bdy = a0+a1*eta_bdy+a2*eta_bdy**2
       return
    end function extrapolate_bdy

subroutine lininterp (arrin, yin, ncell, nlevin, arrout, &
                      yout, nlevout)
!----------------------------------------------------------------------- 
! Purpose: Do a linear interpolation from input mesh defined by yin to output
!          mesh defined by yout.  Where extrapolation is necessary, values will
!          be copied from the extreme edge of the input grid.  Vectorization is over
!          the vertical (ncell) dimension. 
! Method: Check validity of input, then determine weights, then do the N-S interpolation.
! Author: Jim Rosinski
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
! Arguments
!
   real(r8), intent(in)  :: arrin(nlevin,ncell)    ! input array of values to interpolate
   real(r8), intent(in)  :: yin(nlevin)            ! input mesh
   integer,  intent(in)  :: ncell                  ! number of vertical levels
   integer,  intent(in)  :: nlevin                 ! number of input latitudes
   real(r8), intent(out) :: arrout(nlevout,ncell)  ! interpolated array
   real(r8), intent(in)  :: yout(nlevout)          ! output mesh
   integer,  intent(in)  :: nlevout                ! number of output latitudes

!
! Local workspace
!
   integer j, jj              ! latitude indices
   integer js, jn, jjprev     ! latitude indices
   integer k                  ! level index
   integer icount             ! number of values

   real(r8) extrap            ! percent grid non-overlap
!
! Dynamic
!
   integer :: jjm(nlevout)
   integer :: jjp(nlevout)

   real(r8) :: wgts(nlevout)
   real(r8) :: wgtn(nlevout)
!
! Check validity of input coordinate arrays: must be monotonically increasing,
! and have a total of at least 2 elements
!
   if (nlevin.lt.2) then
      call endrun ('LININTERP: Must have at least 2 input points for interpolation')
   end if

   icount = 0
   do j=1,nlevin-1
      if (yin(j).gt.yin(j+1)) icount = icount + 1
   end do

   do j=1,nlevout-1
      if (yout(j).gt.yout(j+1)) icount = icount + 1
   end do

   if (icount.gt.0) then
      call endrun ('LININTERP: Non-monotonic coordinate array(s) found')
   end if
!
! Initialize index arrays for later checking
!
   do j=1,nlevout
      jjm(j) = 0
      jjp(j) = 0
   end do
!
! For values which extend beyond N and S boundaries, set weights
! such that values will just be copied.
!
   do js=1,nlevout
      if (yout(js).gt.yin(1)) goto 10
      jjm(js) = 1
      jjp(js) = 1
      wgts(js) = 1.
      wgtn(js) = 0.
   end do
10 do jn=nlevout,1,-1
      if (yout(jn).le.yin(nlevin)) goto 20
      jjm(jn) = nlevin
      jjp(jn) = nlevin
      wgts(jn) = 1.
      wgtn(jn) = 0.
   end do
!
! Loop though output indices finding input indices and weights
!
20 jjprev = 1
   do j=js,jn
      do jj=jjprev,nlevin-1
         if (yout(j).gt.yin(jj) .and. yout(j).le.yin(jj+1)) then
            jjm(j) = jj
            jjp(j) = jj + 1
            wgts(j) = (yin(jj+1)-yout(j))/(yin(jj+1)-yin(jj))
            wgtn(j) = (yout(j)-yin(jj))/(yin(jj+1)-yin(jj))
            goto 30
         end if
      end do
      write(6,*)'LININTERP: Failed to find interp values'
30    jjprev = jj
   end do
!
! Check grid overlap
!
   extrap = 100.*((js - 1.) + float(nlevout - jn))/nlevout
   if (extrap.gt.30.) then
      if(mpi_rank().eq.0) write(6,*)'LININTERP WARNING:',extrap,' % of output grid will have to be extrapolated'
   end if
!
! Check that interp/extrap points have been found for all outputs
!
   icount = 0
   do j=1,nlevout
      if (jjm(j).eq.0 .or. jjp(j).eq.0) icount = icount + 1
   end do
   if (icount.gt.0) then
      call endrun ('LININTERP: Point found without interp indices')
   end if
!
! Do the interpolation
!
   do j=1,nlevout
      do k=1,ncell
         arrout(j,k) = arrin(jjm(j),k)*wgts(j) + arrin(jjp(j),k)*wgtn(j)
      end do
   end do

   return
   end subroutine lininterp
    
   function arclen(p, q)
    !-----------------------------------------------------------
    !   This function computes the arc-length (angle in radians)            
    !    between a pair of points on the unit sphere. It is similar to
    !    arcdistxyz, but uses a calculation to avoid inverse cosine function
    !    p,q = arrays of length 3 containing the x, y, and z             
    !             coordinates (in that order) of points on the              
    !             unit sphere.                                                 
    !   Returns the angle in radians between the unit vectors              
    !       p and q.  0 <= arclen <= pi.
    !---------------------------------------------------------
! io
    real (r8), intent(in) :: p(3)
    real (r8), intent(in) :: q(3)
! local
    real (r8)             :: arclen
    real (r8)             :: dot

    dot = dot_product(p+q, p+q)

    if (dot.eq.0.) then 
       ! p and q are separated by 180 degrees.
       arclen = pi
    elseif (dot.ge.4._r8) then 
       ! p and q coincide.
       arclen = 0._r8
    else 
       arclen = 2._r8 * atan (sqrt ( (4._r8 - dot) / dot) )
    endif

    return 
   end function arclen

!----------------------------------------------------
! Arc length between p1 and p2 points on the sphere
!     Lat [-pi/2, pi/2] , lon [-pi, pi]
!-----------------------------------------------------

   function arcdistll_r8(lon1, lat1, lon2 , lat2)
    real (r8), intent(in) :: lat1
    real (r8), intent(in) :: lon1
    real (r8), intent(in) :: lat2
    real (r8), intent(in) :: lon2
    real (r8)             :: arcdistll_r8
    real (r8)             :: dlat
    real (r8)             :: dlon

    arcdistll_r8 = 0
    arcdistll_r8 = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2))
    arcdistll_r8 = abs(arcdistll_r8)

    return
  end function arcdistll_r8

  function arcdistll_r4(lon1, lat1, lon2 , lat2)
    real (r4), intent(in) :: lat1
    real (r4), intent(in) :: lon1
    real (r4), intent(in) :: lat2
    real (r4), intent(in) :: lon2
    real (r4)             :: arcdistll_r4
    real (r4)             :: dlat
    real (r4)             :: dlon

    arcdistll_r4 = 0
    arcdistll_r4 = acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon1-lon2))
    arcdistll_r4 = abs(arcdistll_r4)

    return
  end function arcdistll_r4

!------------------------------------------------------------
! The area of a spherical triangle on the unit sphere given 
!------------------------------------------------------------
  function sphtriarea(p1, p2, p3)

    !Cartesian coordinates of the nodes of spherical triangle
    real(r8), intent(in) :: p1(1:3)
    real(r8), intent(in) :: p2(1:3)
    real(r8), intent(in) :: p3(1:3)

    !Variables for spherical triangle area calculation
    ! a, b, c are lengths of the the 3 sides
    ! s is the semi-perimiter s=(a+b+c)/2
    ! e is the spherical excess of triangle
    !   e = a + b + c - 180, but we will use L'Huilier's Theorem
    real (r8)            :: a, b ,c, s, e
    real (r8)            :: tmp
    !Spherical triangle area
    real (r8)            :: sphtriarea

    !Check for degeneration
    if(norm(p1-p2)<eps/100 .or. norm(p2-p3)<eps/100 .or. norm(p1-p3)<eps/100)then
       sphtriarea=0._r8
       return
    end if

    !Calculate the sides length's and the semiperimiter
    s = 0.
    a = arclen(p1,p2)
    b = arclen(p2,p3)
    c = arclen(p3,p1)
    s = (a+b+c)/2

    !Calculate spherical triangle excess using L'Huilier's Theorem
    tmp = tan(s/2)*tan((s-a)/2)*tan((s-b)/2)*tan((s-c)/2)

    !Round off error might give almost zero negative numbers => assume zero
    if(tmp<0)then
       e=0.
    else
       e = 4*atan(sqrt(tmp))
    end if
    sphtriarea=e

    return
  end function sphtriarea

  function areas ( v1, v2, v3 )

!*****************************************************************************80
!
!  AREAS computes the area of a spherical triangle on the unit sphere.
!  stripack version
!
!  Parameters:
!
!    Input, real ( kind = 8 ) V1(3), V2(3), V3(3), the Cartesian coordinates
!    of unit vectors (the three triangle vertices in any order).  These 
!    vectors, if nonzero, are implicitly scaled to have length 1.
!
!    Output, real ( kind = 8 ) AREAS, the area of the spherical triangle 
!    defined by V1, V2, and V3, in the range 0 to 2*PI (the area of a
!    hemisphere).  AREAS = 0 (or 2*PI) if and only if V1, V2, and V3 lie in (or
!    close to) a plane containing the origin.
!
!  Local parameters:
!
!    A1,A2,A3 =    Interior angles of the spherical triangle.
!
!    CA1,CA2,CA3 = cos(A1), cos(A2), and cos(A3), respectively.
!
!    DV1,DV2,DV3 = copies of V1, V2, and V3.
!
!    I =           DO-loop index and index for Uij.
!
!    S12,S23,S31 = Sum of squared components of U12, U23, U31.
!
!    U12,U23,U31 = Unit normal vectors to the planes defined by
!                 pairs of triangle vertices.
!
  implicit none

  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) a3
  real ( kind = 8 ) areas
  real ( kind = 8 ) ca1
  real ( kind = 8 ) ca2
  real ( kind = 8 ) ca3
  real ( kind = 8 ) dv1(3)
  real ( kind = 8 ) dv2(3)
  real ( kind = 8 ) dv3(3)
  real ( kind = 8 ) s12
  real ( kind = 8 ) s23
  real ( kind = 8 ) s31
  real ( kind = 8 ) u12(3)
  real ( kind = 8 ) u23(3)
  real ( kind = 8 ) u31(3)
  real ( kind = 8 ) v1(3)
  real ( kind = 8 ) v2(3)
  real ( kind = 8 ) v3(3)

  dv1(1:3) = v1(1:3)
  dv2(1:3) = v2(1:3)
  dv3(1:3) = v3(1:3)
!
!  Compute cross products Uij = Vi X Vj.
!
  u12(1) = dv1(2) * dv2(3) - dv1(3) * dv2(2)
  u12(2) = dv1(3) * dv2(1) - dv1(1) * dv2(3)
  u12(3) = dv1(1) * dv2(2) - dv1(2) * dv2(1)

  u23(1) = dv2(2) * dv3(3) - dv2(3) * dv3(2)
  u23(2) = dv2(3) * dv3(1) - dv2(1) * dv3(3)
  u23(3) = dv2(1) * dv3(2) - dv2(2) * dv3(1)

  u31(1) = dv3(2) * dv1(3) - dv3(3) * dv1(2)
  u31(2) = dv3(3) * dv1(1) - dv3(1) * dv1(3)
  u31(3) = dv3(1) * dv1(2) - dv3(2) * dv1(1)
!
!  Normalize Uij to unit vectors.
!
  s12 = dot_product ( u12(1:3), u12(1:3) )
  s23 = dot_product ( u23(1:3), u23(1:3) )
  s31 = dot_product ( u31(1:3), u31(1:3) )
!
!  Test for a degenerate triangle associated with collinear vertices.
!
  if ( s12 == 0.0D+00 .or. s23 == 0.0D+00 .or. s31 == 0.0D+00 ) then
    areas = 0.0D+00
    return
  end if

  s12 = sqrt ( s12 )
  s23 = sqrt ( s23 )
  s31 = sqrt ( s31 )

  u12(1:3) = u12(1:3) / s12
  u23(1:3) = u23(1:3) / s23
  u31(1:3) = u31(1:3) / s31
!
!  Compute interior angles Ai as the dihedral angles between planes:
!  CA1 = cos(A1) = -<U12,U31>
!  CA2 = cos(A2) = -<U23,U12>
!  CA3 = cos(A3) = -<U31,U23>
!
  ca1 = - dot_product ( u12(1:3), u31(1:3) )
  ca2 = - dot_product ( u23(1:3), u12(1:3) )
  ca3 = - dot_product ( u31(1:3), u23(1:3) )

  ca1 = max ( ca1, -1.0D+00 )
  ca1 = min ( ca1, +1.0D+00 )
  ca2 = max ( ca2, -1.0D+00 )
  ca2 = min ( ca2, +1.0D+00 )
  ca3 = max ( ca3, -1.0D+00 )
  ca3 = min ( ca3, +1.0D+00 )

  a1 = acos ( ca1 )
  a2 = acos ( ca2 )
  a3 = acos ( ca3 )
!
!  Compute AREAS = A1 + A2 + A3 - PI.
!
  areas = a1 + a2 + a3 - acos ( -1.0D+00 )

  if ( areas < 0.0D+00 ) then
    areas = 0.0D+00
  end if

  return
end function areas

!--------------------------------------
! i if 1<=i<=n
! i-n if n< i <2n
! i+n if -n < i <= 0
! i for other
!--------------------------------------

  function modint(i, n)
    integer (i4), intent(in) :: i
    integer (i4), intent(in) :: n
    integer (i4):: modint

    modint=i

    if( i> n .and. i < 2*n )then
       modint = i-n
    elseif (-n < i .and. i <=0 )then
       modint = i + n
    end if

    return
  end function modint

!-------------------------------------------------------
! calculate the determinant of the matrix by p1, p2, p3
!-------------------------------------------------------
  function det(p1, p2, p3)
    real (r8), intent(in) :: p1(1:3)
    real (r8), intent(in) :: p2(1:3)
    real (r8), intent(in) :: p3(1:3)
    real (r8):: det

    det=dot_product(cross_product(p1,p2),p3)

    return
  end function det

!-------------------------------------------------------------------------------
!  Spherical coordinates (lat,lon) to cartesian coordinates.
!  lat [-pi/2,pi/2], lon [-pi,pi]
!  output:  X, Y, Z, the coordinates in the range -1 to 1. X^2 + Y^2 + Z^2 = 1 
!-------------------------------------------------------------------------------
  subroutine sph2cart (lon, lat, x, y, z )
! io
    real (r8), intent(in)  :: lon
    real (r8), intent(in)  :: lat
    real (r8), intent(out) :: x
    real (r8), intent(out) :: y
    real (r8), intent(out) :: z
! local
    real (r8)              :: coslat

    coslat = cos (lat)
    x = coslat * cos (lon)
    y = coslat * sin (lon)
    z = sin (lat)

    return
  end subroutine sph2cart

!---------------------------------------------------------------------
! Cartesian coordinates to spherical (lat,lon) coordinates.
! X, Y, Z, the coordinates in the range -1 to 1. X^2 + Y^2 + Z^2 = 1
! lat [-pi/2,pi/2], lon [-pi,pi]
!---------------------------------------------------------------------

  subroutine cart2sph ( x, y, z, lon, lat )
! io
    real(r8), intent(in)  :: x
    real(r8), intent(in)  :: y
    real(r8), intent(in)  :: z
    real(r8), intent(out) :: lat
    real(r8), intent(out) :: lon
! local
    real(r8)              :: pnrm

    pnrm = sqrt ( x **2 + y**2 + z**2 )
    if ( x /= 0.0D+00 .or. y /= 0.0D+00 ) then
       lon = atan2 ( y, x )
    else
       lon = 0.0D+00
    end if

    if ( pnrm /= 0.0D+00 ) then
       lat = asin ( z / pnrm )
    else
       print*, "cart2sph warning: point not in the unit sphere. norm= ", pnrm
       lat = 0.0D+00
    end if

    return
  end subroutine cart2sph

  function cross_product_r8(a,b)
    implicit none

    real(r8), intent(in) :: a(1:3)
    real(r8), intent(in) :: b(1:3)
! local    
    real(r8)             :: cross_product_r8(1:3)

    cross_product_r8(1) = a(2)*b(3) - a(3)*b(2)                                    
    cross_product_r8(2) = a(3)*b(1) - a(1)*b(3)
    cross_product_r8(3) = a(1)*b(2) - a(2)*b(1)

    return
  end function cross_product_r8

  function cross_product_r4(a,b)
    implicit none

    real(r4), intent(in) :: a(1:3)
    real(r4), intent(in) :: b(1:3)
! local    
    real(r4)             :: cross_product_r4(1:3)

    cross_product_r4(1) = a(2)*b(3) - a(3)*b(2)                                    
    cross_product_r4(2) = a(3)*b(1) - a(1)*b(3)
    cross_product_r4(3) = a(1)*b(2) - a(2)*b(1)

    return
  end function cross_product_r4

!---------------------------------------------------------------------
!   In: a point (x,y,z) and a vector at this point
!   Out: Returns the vector at this point in (lon,lat) reference
!   The vector must be tangent to the sphere, that is, V.(x,y,z)=0
!---------------------------------------------------------------------

  subroutine convert_vector_cart2sph_r8(point, vector , vector_u_component, vector_v_component)
!
! io
!
    real(r8), intent(in)  :: point(1:3)   ! Point cartesian coords
    real(r8), intent(in)  :: vector(1:3)  ! Cartesian vector on point
    real(r8), intent(out) :: vector_u_component   ! Spherical coord vector component
    real(r8), intent(out) :: vector_v_component   ! Spherical coord vector component
!
! local
!
    real(r8)              :: radius
    real(r8)              :: rho
    real(r8)              :: eel(1:3)
    real(r8)              :: eaz(1:3)
    real(r8)              :: rvec(1:3)
    real(r8)              :: zero
    real(r8)              :: test

    zero      = 0
    radius    = sqrt(point(1)**2+point(2)**2+point(3)**2)
    rho       = sqrt(point(1)**2+point(2)**2)

!Case where the point is in the north or south pole
    if(rho<10*eps)then
       vector_u_component=vector(2)
       vector_v_component=vector(1)
       return
    end if

    rvec   = (/point(1),point(2),point(3)/)
    rvec   = rvec/radius
!
! see mathworks "spherical basis representation of vectors", 
! assume a unit sphere
!
    eaz    = (/-point(2)/rho, point(1)/rho, zero/)
    eel    = (/-point(1)*point(3)/rho,-point(2)*point(3)/rho,rho/) ! old one is (rho**2)/rho

    test   = dot_product(vector,rvec)  ! should be zero

    if(abs(test)>10e-5)then
    !if(abs(test)>10e-3)then
       print *, "convert_vector_cart2sph warning: vector not tangent to sphere."
       print *, "vector:",vector(1:3)
       print *, "point:" ,rvec(1:3)
       print *, "dot product:", test
       stop
    end if

    vector_u_component   = dot_product(vector,eaz)
    vector_v_component   = dot_product(vector,eel)

    return    
  end subroutine convert_vector_cart2sph_r8

  subroutine convert_vector_cart2sph_r4(point, vector , vector_u_component, vector_v_component)
!
! io
!
    real(r4), intent(in)  :: point(1:3)   ! Point cartesian coords
    real(r4), intent(in)  :: vector(1:3)  ! Cartesian vector on point
    real(r4), intent(out) :: vector_u_component   ! Spherical coord vector component
    real(r4), intent(out) :: vector_v_component   ! Spherical coord vector component
!
! local
!
    real(r4)              :: radius
    real(r4)              :: rho
    real(r4)              :: eel(1:3)
    real(r4)              :: eaz(1:3)
    real(r4)              :: rvec(1:3)
    real(r4)              :: zero
    real(r4)              :: test

    zero      = 0
    radius    = sqrt(point(1)**2+point(2)**2+point(3)**2)
    rho       = sqrt(point(1)**2+point(2)**2)

!Case where the point is in the north or south pole
    if(rho<10*eps)then
       vector_u_component=vector(2)
       vector_v_component=vector(1)
       return
    end if

    rvec   = (/point(1),point(2),point(3)/)
    rvec   = rvec/radius
!
! see mathworks "spherical basis representation of vectors", assume a unit sphere
!
    eaz    = (/-point(2)/rho, point(1)/rho, zero/)
    eel    = (/-point(1)*point(3)/rho,-point(2)*point(3)/rho,rho/) ! old one is (rho**2)/rho

    test   = dot_product(vector,rvec)  ! should be zero
!
! using single precision is very sensitive to this line!
!
    !if(abs(test)>10e-5)then
    if(abs(test)>10e-3)then
       print *, "convert_vector_cart2sph warning: vector not tangent to sphere."
       print *, "vector:",vector(1:3)
       print *, "point:" ,rvec(1:3)
       print *, "dot product:", test
      !  stop
    end if

    vector_u_component   = dot_product(vector,eaz)
    vector_v_component   = dot_product(vector,eel)

    return    
  end subroutine convert_vector_cart2sph_r4

!---------------------------------------------------------------------
!   Input: a point(x,y,z) and a vector at this point in lon, lat system
!   in radians (vector_u_component, vector_v_component), ie, 
!   Output: the vector at this point in cartesian reference (v)
!---------------------------------------------------------------------

  subroutine convert_vector_sph2cart_r8(vector_u_component, vector_v_component, point, vector)
! io
    real(r8), intent(in)  :: vector_u_component
    real(r8), intent(in)  :: vector_v_component
    real(r8), intent(in)  :: point(1:3)  !Point cartesian coords
    real(r8), intent(out) :: vector(1:3)
! local
    real(r8)              :: radius
    real(r8)              :: rho

    radius = sqrt(point(1)**2+point(2)**2+point(3)**2)
    rho    = sqrt(point(1)**2+point(2)**2)
!
! Case where the point is in the north or south pole
!
    if(rho==0)then
       vector(1) = vector_v_component
       vector(2) = vector_u_component
       vector(3) = 0
       return
    else

       vector(1) = -vector_u_component*(point(2)/rho) - (vector_v_component*point(1)*point(3))/(rho)
       vector(2) =  vector_u_component*(point(1)/rho) - (vector_v_component*point(2)*point(3))/(rho)
       vector(3) =  vector_v_component*rho
       
    end if

    return    
  end subroutine convert_vector_sph2cart_r8
!
! r4 version of above
!
  subroutine convert_vector_sph2cart_r4(vector_u_component, vector_v_component, point, vector)
! io
    real(r4), intent(in)  :: vector_u_component
    real(r4), intent(in)  :: vector_v_component
    real(r4), intent(in)  :: point(1:3)  !Point cartesian coords
    real(r4), intent(out) :: vector(1:3)
! local
    real(r4)              :: radius
    real(r4)              :: rho

    radius = sqrt(point(1)**2+point(2)**2+point(3)**2)
    rho    = sqrt(point(1)**2+point(2)**2)
!
! case where the point is in the north or south pole
!
    if(rho==0)then
       vector(1) = vector_v_component
       vector(2) = vector_u_component
       vector(3) = 0
       return
    else

       vector(1) = -vector_u_component*(point(2)/rho) - (vector_v_component*point(1)*point(3))/(rho)
       vector(2) =  vector_u_component*(point(1)/rho) - (vector_v_component*point(2)*point(3))/(rho)
       vector(3) =  vector_v_component*rho
       
    end if

    return    
  end subroutine convert_vector_sph2cart_r4
!
! norm of a vector
!
  function norm_r8(pos)
    real(r8), intent(in) :: pos(:)
    real(r8):: norm_r8
    norm_r8 = dot_product( pos, pos)
    norm_r8 = sqrt(norm_r8)
    return
  end function norm_r8
  
  function norm_r4(pos)
    real(r4), intent(in) :: pos(:)
    real(r4):: norm_r4
    norm_r4 = dot_product( pos, pos)
    norm_r4 = sqrt(norm_r4)
    return
  end function norm_r4

!---------------------------------------------------------
! Area of triangle intersection of Voronoi/hexagonal cell
!---------------------------------------------------------
  function plg_kite_areas(cell, nnb, mesh)
    integer(i4), intent(in)         :: cell
    integer(i4), intent(in)         :: nnb
    type(global_domain), intent(in) :: mesh

    real(r8), dimension(1:nnb)::plg_kite_areas
    real(r8), dimension(1:3)  :: pc
    real(r8), dimension(1:3)  :: ped1
    real(r8), dimension(1:3)  :: pv
    real(r8), dimension(1:3)  :: ped2
    !Areas
    real(r8)    :: area, area1, area2
    !Indexes for points
    integer(i4) :: ed1, ed2, tr, inb

    pc = mesh%vtx(cell)%p
    do inb = 1, nnb
       ed1   = mesh%vtx(cell)%ed(inb)
       ed2   = mesh%vtx(cell)%ed(modint(inb+1,nnb))
       ped1  = mesh%edp(ed1)%c%p
       ped2  = mesh%edp(ed2)%c%p
       tr    = mesh%vtx(cell)%tr(inb)
       pv    = mesh%tri(tr)%c%p
       area1 = sphtriarea(pc, ped1, pv)
       area2 = sphtriarea(pc, pv, ped2)
       if(det(pc, ped1, pv)<0) area1=-area1
       if(det(pc, pv, ped2)<0) area2=-area2
       plg_kite_areas(inb) = area1+area2
    end do

    return
  end function plg_kite_areas

!-------------------------------------------------
! Kite area for all trinalge in mesh
! Area of triangle intersection with Voronoi
!-------------------------------------------------
  subroutine tri_kite_areas(mesh)
    type(global_domain), intent(inout) :: mesh

    real(r8), dimension(1:3) :: pc
    real(r8), dimension(1:3) :: pv
    real(r8), dimension(1:3) :: ped1
    real(r8), dimension(1:3) :: ped2
    real(r8)    :: area1, area2
    integer(i4) :: ed1, ed2, it, inb

    do it = 1, mesh%nt
       pc = mesh%tri(it)%c%p
       do inb = 1, 3
          ed1  = mesh%tri(it)%ed(modint(inb-1,3))
          ed2  = mesh%tri(it)%ed(inb)
          ped1 = mesh%edt(ed1)%c%p
          ped2 = mesh%edt(ed2)%c%p
          pv   = mesh%vtx(mesh%tri(it)%v(inb))%p
          area1= sphtriarea(pc, ped1, pv)
          area2= sphtriarea(pc, pv, ped2)
! normally will not happen for a quality mesh
          if(det(pc, ped1, pv)<0) area1=-area1
          if(det(pc, pv, ped2)<0) area2=-area2
          mesh%tri(it)%kite_area(inb) = area1+area2
       end do
    end do

    return
  end subroutine tri_kite_areas
end module grist_math_module
