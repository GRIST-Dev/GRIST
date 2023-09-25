!===================================================================================================
!
!  Created by LiXiaohan on 19/06/04, adopted from CAM5
!
!  Module to hold PBL-related subprograms that may be used with multiple
!  different vertical diffusion schemes
!
!===================================================================================================

module grist_pbl_utils

    use grist_constants,                    only: r8, i4,           &
                                                  gravity, karman,  &
                                                  cp, rdry, zvir

    implicit none
    private

! Public Procedures:
! Excepting the initialization procedure, these are elemental
! procedures, so they can accept scalars or any dimension of array as
! arguments, as long as all arguments have the same number of
! elements.

    public :: calc_ustar,       &
              calc_obklen,      &
              virtem

    real(r8), parameter :: ustar_min = 0.01_r8

contains

! Purpose: Calculate ustar and bottom level density, necessary for Obukhov length calculation
    elemental subroutine calc_ustar(t, pmid, taux, tauy, rrho, ustar)
! io
    real(r8), intent(in) :: t         ! surface temperature
    real(r8), intent(in) :: pmid      ! midpoint pressure (bottom level)
    real(r8), intent(in) :: taux      ! surface u stress [N/m2]
    real(r8), intent(in) :: tauy      ! surface v stress [N/m2]

    real(r8), intent(out) :: rrho     ! 1./bottom level density
    real(r8), intent(out) :: ustar    ! surface friction velocity [m/s]

    rrho  = rdry * t / pmid
    ustar = max( sqrt( sqrt(taux**2 + tauy**2)*rrho ), ustar_min )
 
    end subroutine calc_ustar


! Purpose: Calculate Obukhov length and kinematic fluxes.
    elemental subroutine calc_obklen( ths,  thvs, qflx, shflx, rrho, ustar, &
                                      khfs, kqfs, kbfs, obklen)
! io
    real(r8), intent(in)  :: ths           ! potential temperature at surface [K]
    real(r8), intent(in)  :: thvs          ! virtual potential temperature at surface
    real(r8), intent(in)  :: qflx          ! water vapor flux (kg/m2/s)
    real(r8), intent(in)  :: shflx         ! surface heat flux (W/m2)
    real(r8), intent(in)  :: rrho          ! 1./bottom level density [ m3/kg ]
    real(r8), intent(in)  :: ustar         ! Surface friction velocity [ m/s ]
  
    real(r8), intent(out) :: khfs          ! sfc kinematic heat flux [mK/s]
    real(r8), intent(out) :: kqfs          ! sfc kinematic water vapor flux [m/s]
    real(r8), intent(out) :: kbfs          ! sfc kinematic buoyancy flux [m^2/s^3]
    real(r8), intent(out) :: obklen        ! Obukhov length
  
    ! Need kinematic fluxes for Obukhov:
    khfs = shflx*rrho/cp
    kqfs = qflx*rrho
    kbfs = khfs + zvir*ths*kqfs
  
    ! Compute Obukhov length:
    obklen = -thvs * ustar**3 / (gravity*karman*(kbfs + sign(1.e-10_r8,kbfs)))

    end subroutine calc_obklen


! Purpose: Calculate virtual temperature from temperature and specific humidity
    elemental real(r8) function virtem(t,q)

    real(r8), intent(in) :: t, q

    virtem = t * (1.0_r8 + zvir*q)

    end function virtem

end module grist_pbl_utils
