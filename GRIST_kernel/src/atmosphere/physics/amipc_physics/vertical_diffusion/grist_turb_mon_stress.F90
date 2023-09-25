!===================================================================================================
!
!  Created by LiXiaohan on 19/06/04, adopted from CAM5
!
!===================================================================================================

module grist_turb_mon_stress

    use grist_constants,                    only: r8, i4


    implicit none
    private
    public :: init_tms,     &
              compute_tms

    real(r8), parameter :: horomin= 1._r8       ! Minimum value of subgrid orographic height for mountain stress [ m ]
    real(r8), parameter :: z0max  = 100._r8     ! Maximum value of z_0 for orography [ m ]
    real(r8), parameter :: dv2min = 0.01_r8     ! Minimum shear squared [ m2/s2 ]
    real(r8)            :: orocnst              ! Converts from standard deviation to height [ no unit ]
    real(r8)            :: z0fac                ! Factor determining z_0 from orographic standard deviation [ no unit ] 
    real(r8)            :: karman               ! von Karman constant
    real(r8)            :: gravit               ! Acceleration due to gravity
    real(r8)            :: rair                 ! Gas constant for dry air

contains

    subroutine init_tms( oro_in, z0fac_in, karman_in, gravit_in, rair_in )
! io
    real(r8), intent(in) :: oro_in, z0fac_in, karman_in, gravit_in, rair_in

    orocnst  = oro_in
    z0fac    = z0fac_in
    karman   = karman_in
    gravit   = gravit_in
    rair     = rair_in
    
  end subroutine init_tms

  subroutine compute_tms( pver    , ncol    ,                     &
                          u        , v       , t       , pmid    , exner   , &
                          zm       , sgh     , ksrf    , taux    , tauy    , & 
                          landfrac )

    !------------------------------------------------------------------------------ !
    ! Turbulent mountain stress parameterization                                    !  
    !                                                                               !
    ! Returns surface drag coefficient and stress associated with subgrid mountains !
    ! For points where the orographic variance is small ( including ocean ),        !
    ! the returned surface drag coefficient and stress is zero.                     !
    !                                                                               !
    ! Lastly arranged : Sungsu Park. Jan. 2010.                                     !
    !------------------------------------------------------------------------------ !

    ! ---------------------- !
    ! Input-Output Arguments ! 
    ! ---------------------- !

    integer,  intent(in)  :: pver                  ! Number of model layers
    integer,  intent(in)  :: ncol                  ! Number of columns actually used

    real(r8), intent(in)  :: u(pver,ncol)         ! Layer mid-point zonal wind [ m/s ]
    real(r8), intent(in)  :: v(pver,ncol)         ! Layer mid-point meridional wind [ m/s ]
    real(r8), intent(in)  :: t(pver,ncol)         ! Layer mid-point temperature [ K ]
    real(r8), intent(in)  :: pmid(pver,ncol)      ! Layer mid-point pressure [ Pa ]
    real(r8), intent(in)  :: exner(pver,ncol)     ! Layer mid-point exner function [ no unit ]
    real(r8), intent(in)  :: zm(pver,ncol)        ! Layer mid-point height [ m ]
    real(r8), intent(in)  :: sgh(ncol)            ! Standard deviation of orography [ m ]
    real(r8), intent(in)  :: landfrac(ncol)       ! Land fraction [ fraction ]
    
    real(r8), intent(out) :: ksrf(ncol)           ! Surface drag coefficient [ kg/s/m2 ]
    real(r8), intent(out) :: taux(ncol)           ! Surface zonal      wind stress [ N/m2 ]
    real(r8), intent(out) :: tauy(ncol)           ! Surface meridional wind stress [ N/m2 ]

    ! --------------- !
    ! Local Variables !
    ! --------------- !

    integer  :: i                                  ! Loop index
    integer  :: kb, kt                             ! Bottom and top of source region
    
    real(r8) :: horo                               ! Orographic height [ m ]
    real(r8) :: z0oro                              ! Orographic z0 for momentum [ m ]
    real(r8) :: dv2                                ! (delta v)**2 [ m2/s2 ]
    real(r8) :: ri                                 ! Richardson number [ no unit ]
    real(r8) :: stabfri                            ! Instability function of Richardson number [ no unit ]
    real(r8) :: rho                                ! Density [ kg/m3 ]
    real(r8) :: cd                                 ! Drag coefficient [ no unit ]
    real(r8) :: vmag                               ! Velocity magnitude [ m /s ]

    ! ----------------------- !
    ! Main Computation Begins !
    ! ----------------------- !
 
    do i = 1, ncol

     ! determine subgrid orgraphic height ( mean to peak )

       horo = orocnst * sgh(i)

     ! No mountain stress if horo is too small

       if( horo < horomin ) then

           ksrf(i) = 0._r8
           taux(i) = 0._r8
           tauy(i) = 0._r8

       else

         ! Determine z0m for orography

           z0oro = min( z0fac * horo, z0max )

         ! Calculate neutral drag coefficient

           cd = ( karman / log( ( zm(pver,i) + z0oro ) / z0oro) )**2

         ! Calculate the Richardson number over the lowest 2 layers

           kt  = pver - 1
           kb  = pver
           dv2 = max( ( u(kt,i) - u(kb,i) )**2 + ( v(kt,i) - v(kb,i) )**2, dv2min )

         ! Modification : Below computation of Ri is wrong. Note that 'Exner' function here is
         !                inverse exner function. Here, exner function is not multiplied in
         !                the denominator. Also, we should use moist Ri not dry Ri.
         !                Also, this approach using the two lowest model layers can be potentially
         !                sensitive to the vertical resolution.  
         ! OK. I only modified the part associated with exner function.

           ri  = 2._r8 * gravit * ( t(kt,i) * exner(kt,i) - t(kb,i) * exner(kb,i) ) * ( zm(kt,i) - zm(kb,i) ) &
                                / ( ( t(kt,i) * exner(kt,i) + t(kb,i) * exner(kb,i) ) * dv2 )

         ! ri  = 2._r8 * gravit * ( t(kt,i) * exner(kt,i) - t(kb,i) * exner(kb,i) ) * ( zm(kt,i) - zm(kb,i) ) &
         !                      / ( ( t(kt,i) + t(kb,i) ) * dv2 )

         ! Calculate the instability function and modify the neutral drag cofficient.
         ! We should probably follow more elegant approach like Louis et al (1982) or Bretherton and Park (2009) 
         ! but for now we use very crude approach : just 1 for ri < 0, 0 for ri > 1, and linear ramping.

           stabfri = max( 0._r8, min( 1._r8, 1._r8 - ri ) )
           cd      = cd * stabfri

         ! Compute density, velocity magnitude and stress using bottom level properties

           rho     = pmid(pver,i) / ( rair * t(pver,i) ) 
           vmag    = sqrt( u(pver,i)**2 + v(pver,i)**2 )
           ksrf(i) = rho * cd * vmag * landfrac(i)
           taux(i) = -ksrf(i) * u(pver,i)
           tauy(i) = -ksrf(i) * v(pver,i)

       end if

    end do
    
    return
  end subroutine compute_tms


end module grist_turb_mon_stress
