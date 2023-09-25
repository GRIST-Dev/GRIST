
!-----------------------------------------------------------
! Created on 2016
! Author: Yi Zhang
! Version 1.0
! Description: Module for 1D Flux operators
!              Common usage:
!              1) wind: wind speed follows 0->1 direction
!              2) phi0: phi at the upwind direction
!              3) phi1: phi at the downwind direction
!-----------------------------------------------------------

 module grist_base_flux_module

   use grist_constants,  only: r8, rearth, zero

   implicit none

    private

    public  ::    flux_laxwen, &
                  flux_upwind, &
                  flux_upwind_nowind, &
                  flux_laxwen_nowind, &
                  flux_diff,   & ! lax-upwind
                  flux_ooo34,  &
                  flux_ooo33,  &
                  flux_ooo56,  &
                  flux_ooo55,  &
                  flux_ppm

    real(r8), parameter :: one = 1._r8

  contains

    subroutine flux_laxwen(wind,phi0,phi1,dt,length,flux)
      implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1
      real(r8), intent(in)  :: dt
      real(r8), intent(in)  :: length
      real(r8), intent(out) :: flux

      flux = 0.5_r8*wind*(phi0+phi1)-0.5_r8*abs(wind)*(phi1-phi0)*abs(wind)*dt/length

      return
    end subroutine flux_laxwen

    subroutine flux_laxwen_nowind(wind,phi0,phi1,dt,length,flux)
      implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1
      real(r8), intent(in)  :: dt
      real(r8), intent(in)  :: length
      real(r8), intent(out) :: flux

      flux = 0.5_r8*(phi0+phi1)-0.5_r8*sign(1._r8,wind)*(phi1-phi0)*abs(wind)*dt/length

      return
    end subroutine flux_laxwen_nowind
    subroutine flux_upwind(wind,phi0,phi1,flux)
      implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1
      real(r8), intent(out) :: flux

      flux = wind*(0.5_r8*(phi0+phi1)-0.5_r8*sign(one,wind)*(phi1-phi0))
      return
    end subroutine flux_upwind

    subroutine flux_ooo34(wind,phi0,phi1,der0,der1,flux,length,beta)
    implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1  
      real(r8), intent(in)  :: der0   ! 2nd order derivative at upwind direction
      real(r8), intent(in)  :: der1   ! 2nd order derivative at downwind direction
      real(r8), intent(out) :: flux  
      real(r8), intent(in)  :: length
      real(r8), intent(in)  :: beta   ! beta in SG11
! local
      real(r8)              :: part1
      real(r8)              :: part2
      real(r8)              :: part3

      part1 = (phi0+phi1)*0.5_r8
      part2 = -1._r8*((length**2)/12._r8)*(der0+der1)
      part3 = (sign(1._r8,wind)*(length**2)*beta/12._r8)*(der1-der0)
      if(beta.ne.zero)then
         flux  = (part1+part2+part3)
      else
         flux  = part1+part2
      end if

      return
    end subroutine flux_ooo34

    subroutine flux_ooo33(phi0,phi1,der0,flux,length)
    implicit none
!io
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1  
      real(r8), intent(in)  :: der0   ! 2nd order derivative at upwind direction
      real(r8), intent(out) :: flux  
      real(r8), intent(in)  :: length
! local
      real(r8)              :: part1
      real(r8)              :: part2

      part1 = (phi0+phi1)*0.5_r8
      part2 = -1._r8*((length**2)/6._r8)*der0
      flux  = part1+part2

      return
    end subroutine flux_ooo33

    subroutine flux_ooo56(wind,phi0,phi1,der_2nd_0,der_2nd_1,der_4th_0,der_4th_1,flux,length,beta)

      implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1
      real(r8), intent(in)  :: der_2nd_0   ! 2nd order derivative at upwind direction
      real(r8), intent(in)  :: der_2nd_1   ! 2nd order derivative at downwind direction
      real(r8), intent(in)  :: der_4th_0   ! 4th order derivative at upwind direction
      real(r8), intent(in)  :: der_4th_1   ! 4th order derivative at downwind direction
      real(r8), intent(out) :: flux
      real(r8), intent(in)  :: length
      real(r8), intent(in)  :: beta        ! see my paper
! local
      real(r8)              :: part1
      real(r8)              :: part2
      real(r8)              :: part3

      part1 = (phi0+phi1)*0.5_r8
      part2 = ((length**4)/60._r8)*(der_4th_1+der_4th_0)-((length**2)/12._r8)*(der_2nd_1+der_2nd_0)
      part3 = -1._r8*(sign(1._r8,wind))*((length**4)*beta/60._r8)*(der_4th_1-der_4th_0)

      flux  = (part1+part2+part3)

      return
    end subroutine flux_ooo56

    subroutine flux_ooo55(wind,phi0,phi1,der_2nd_0,der_2nd_1,der_4th_0,flux,length,beta)
! only the fifth-order operator
! depending on the upwind fourth-order derivative
      implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1
      real(r8), intent(in)  :: der_2nd_0   ! 2nd order derivative at upwind direction
      real(r8), intent(in)  :: der_2nd_1   ! 2nd order derivative at downwind direction
      real(r8), intent(in)  :: der_4th_0   ! 4th order derivative at upwind direction
      real(r8), intent(out) :: flux
      real(r8), intent(in)  :: length
      real(r8), intent(in)  :: beta        ! see my paper
! local
      real(r8)              :: part1
      real(r8)              :: part2
      real(r8)              :: part3

      part1 = (phi0+phi1)*0.5_r8
      part2 = ((length**4)/30._r8)*der_4th_0-((length**2)/12._r8)*(der_2nd_1+der_2nd_0)
      flux  = (part1+part2)

      return
    end subroutine flux_ooo55

    subroutine flux_diff(wind,phi0,phi1,flux,dt,length)
      implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1   
      real(r8), intent(in)  :: dt
      real(r8), intent(in)  :: length
      real(r8), intent(out) :: flux

      flux = -0.5_r8*abs(wind)*(phi1-phi0)*(abs(wind)*dt/length-1)

      return
    end subroutine flux_diff

! no-wind version of 1st-order upwind
    subroutine flux_upwind_nowind(wind,phi0,phi1,flux)
      implicit none
!io
      real(r8), intent(in)  :: wind
      real(r8), intent(in)  :: phi0
      real(r8), intent(in)  :: phi1
      real(r8), intent(out) :: flux

      flux = 0.5_r8*(phi0+phi1)-0.5_r8*sign(one,wind)*(phi1-phi0)
      return
    end subroutine flux_upwind_nowind

    subroutine flux_ppm(wind,phi0,phi1,der0,der1,flux,length,dtime)
    implicit none
!io
      real(r8), intent(in)  :: wind   ! wind follows 0->1 positive
      real(r8), intent(in)  :: phi0   ! upwind point
      real(r8), intent(in)  :: phi1   ! downwind point
      real(r8), intent(in)  :: der0   ! 2nd order derivative at upwind direction
      real(r8), intent(in)  :: der1   ! 2nd order derivative at downwind direction
      real(r8), intent(out) :: flux
      real(r8), intent(in)  :: length
      real(r8), intent(in)  :: dtime
! local
      real(r8)              :: part1
      real(r8)              :: part2
      real(r8)              :: phi_arrow
      real(r8)              :: phi_upwind
      real(r8)              :: der_upwind
      real(r8)              :: cr_num  ! courant number
      real(r8)              :: ppm_coef

       ppm_coef  = 0.25_r8

       part1     = (phi0+phi1)*0.5_r8
       part2     = -1._r8*((length**2)/12._r8)*(der0+der1)
       phi_arrow = part1+part2
       cr_num    = wind*dtime/(rearth*length)

       call  flux_upwind_nowind(wind,phi0,phi1,phi_upwind)
       call  flux_upwind_nowind(wind,der0,der1,der_upwind)
!---------------------------------
! IF SELECTION VERSION
! same to machine precision
!       if(wind.gt.0)then
!          phi_upwind = phi0
!          der_upwind = der0
!       else
!          phi_upwind = phi1
!          der_upwind = der1
!       end if
!---------------------------------
       flux      = phi_arrow-abs(cr_num)*(phi_arrow-phi_upwind)+((cr_num**2)-abs(cr_num))*ppm_coef*der_upwind*(length**2)
       flux      = wind*flux

       return
    end subroutine flux_ppm

  end module grist_base_flux_module
