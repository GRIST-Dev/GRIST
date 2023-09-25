!================================================================================================
! This is the "toy" chemistry module.
!================================================================================================

MODULE grist_dtp_dcmip2016_terminator

  use grist_constants,   only: pi, r8

implicit none
 
private
save

public :: initial_value_Terminator       ! initialize cl and cl2
public :: tendency_Terminator            ! interface to tendency computation

!integer, parameter :: 8 = selected_real_kind (12)

real(r8), parameter :: cly_constant = 4.e-6_8
!real(r8), parameter :: pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164_8
real(r8), parameter :: half_pi = pi*0.5_8
real(r8), parameter :: degrees_to_radians = pi/180.0_8
real(r8), parameter :: k1_lat_center =   20.d0*degrees_to_radians
real(r8), parameter :: k1_lon_center =  300.d0*degrees_to_radians

contains

!===============================================================================
!  Solar photolysis rate and recombination rate
!===============================================================================

subroutine k_vals( lat, lon, k1, k2 )

!-----------------------------------------------------------------------
! Arguments:
!-----------------------------------------------------------------------
real(r8), intent(in)    :: lat, lon  ! latitude and longitude, radians
real(r8), intent(out)   :: k1, k2    ! reaction rates

k1 = 1._8*max(0.d0,sin(lat)*sin(k1_lat_center) + cos(lat)*cos(k1_lat_center)*cos(lon-k1_lon_center))
k2 = 1._8

return

end subroutine k_vals

!===============================================================================
!  Tendencies of cl and cl2
!===============================================================================

subroutine tendency_Terminator( lat, lon, cl, cl2, dt, cl_f, cl2_f )

!-----------------------------------------------------------------------
! Arguments:
!-----------------------------------------------------------------------

real(r8), intent(in)    :: lat, lon  ! latitude and longitude, degrees
real(r8), intent(in)    :: cl, cl2   ! molar mixing ratio of cl and cl2
real(r8), intent(in)    :: dt        ! size of physics time step

real(r8), intent(out)   :: cl_f, cl2_f  ! time rate of change of cl and cl2

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

real(r8) :: r, det, expdt, el ! useful algebaic quantities used in the computation
real(r8) :: k1, k2            ! reaction rates
real(r8) :: cly               ! quantity that should be conseved

call k_vals( lat*degrees_to_radians, lon*degrees_to_radians, k1, k2 )

r = k1 / (4._8*k2)
cly = cl + 2._8* cl2

det = sqrt( r*r + 2._8*r*cly )
expdt = exp( -4._8*k2*det*dt )

if ( abs(det * k2 * dt) .gt. 1e-16 ) then
el = (1._8 - expdt) /det /dt
else
el = 4._8*k2
endif

cl_f  = -el * (cl - det + r)*(cl + det + r) / (1._8 + expdt + dt*el*(cl + r))
cl2_f = -cl_f / 2._8

return

end subroutine tendency_Terminator

!===============================================================================
!  Compute initial values
!===============================================================================

subroutine initial_value_Terminator( lat, lon, cl, cl2 )

!-----------------------------------------------------------------------
! Arguments:
!-----------------------------------------------------------------------

real(r8), intent(in)  :: lat, lon  ! latitude and longitude, degrees
real(r8), intent(out) :: cl, cl2   ! molar mixing ratio of cl and cl2

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

real(r8) :: r, det  ! useful algebraic forms
real(r8) :: k1, k2  ! reaction rates

call k_vals( lat*degrees_to_radians, lon*degrees_to_radians, k1, k2 )

r = k1 / (4._8*k2)
det = sqrt(r*r + 2._8*cly_constant*r)

cl  = (det-r)
cl2 = cly_constant/2._8 - (det-r)/2._8

return

end subroutine initial_value_Terminator

end module grist_dtp_dcmip2016_terminator



