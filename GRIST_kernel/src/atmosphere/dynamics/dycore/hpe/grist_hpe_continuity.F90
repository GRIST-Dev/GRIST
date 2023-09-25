
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description:  Evaluate continuity equation
! Revision history:
!----------------------------------------------------------------------------

 module grist_hpe_continuity

   use grist_constants,      only: i4, r8, zero
   use grist_hpe_constants,  only: eta_face_b, eta_full_b

   implicit none

   private

   public    :: calc_hpe_tend_continuity_2d
#ifdef REFCODE
               &,calc_hpe_tend_continuity_face
               &,calc_hpe_tend_continuity
#endif

   contains

!-----------------------------------------------------------------------
! old not used codes for reference only
!-----------------------------------------------------------------------

#ifdef REFCODE
    subroutine calc_hpe_tend_continuity(nlev, nlevp, &
                                        scalar_flux_div_full_level , &
                                        scalar_ps_tendency         , &
                                        scalar_mass_eta_velocity)

      integer(i4),                intent(in)    :: nlev, nlevp 
      real(r8), dimension(nlev),  intent(in)    :: scalar_flux_div_full_level ! DEL.(V*delhp)
      real(r8)                 ,  intent(inout) :: scalar_ps_tendency         ! -PARTIAL(PS)/PARTIAL(t)
      real(r8), dimension(nlevp), intent(inout) :: scalar_mass_eta_velocity   ! m*etadot at face

      real(r8)                          :: vert_sum
      integer(i4)                       :: ilev
      integer(i4)                       :: ikkk

!
! ps tendency (3.2 of SS81)
!
         vert_sum = zero
         do ilev = 1, nlev
            vert_sum = vert_sum+scalar_flux_div_full_level(ilev)
         end do
        
         scalar_ps_tendency  = -1._r8*vert_sum
!
! m*eta velocity (3.3 of SS81)
!
         scalar_mass_eta_velocity(1) = zero

         do ilev = 1, nlev
            vert_sum = zero
            do ikkk = 1, ilev
              vert_sum = vert_sum+scalar_flux_div_full_level(ikkk)
            end do
            scalar_mass_eta_velocity(ilev+1) = eta_face_b(ilev+1)*&
                                              (-1._r8*scalar_ps_tendency)-vert_sum
         end do

!
! this should be, or the above routine is bugged
!
        if( abs(scalar_mass_eta_velocity(nlev+1)).gt.1e-15) print*,"sth wrong inside calc_hpe_tend_continuity m*etadot"

      return
   
    end subroutine calc_hpe_tend_continuity

    subroutine calc_hpe_tend_continuity_face(nlev, nlevp,&
                                             scalar_flux_div_face_level , &
                                             scalar_ps_tendency         , &
                                             scalar_mass_eta_velocity)
  
      integer(i4),               intent(in)    :: nlev, nlevp 
      real(r8), dimension(nlevp),intent(in)    :: scalar_flux_div_face_level ! DEL.(V*delhp)
      real(r8)                 , intent(in)    :: scalar_ps_tendency         ! -PARTIAL(PS)/PARTIAL(t)
      real(r8), dimension(nlev), intent(inout) :: scalar_mass_eta_velocity   ! m*etadot at full level

      real(r8)                          :: vert_sum
      integer(i4)                       :: ilev
      integer(i4)                       :: ikkk

         do ilev = 1, nlev
            vert_sum = zero
            do ikkk = 1, ilev
              vert_sum = vert_sum+scalar_flux_div_face_level(ikkk)
            end do
            scalar_mass_eta_velocity(ilev) = eta_full_b(ilev)*(-1._r8*scalar_ps_tendency)-vert_sum
         end do

      return
    end subroutine calc_hpe_tend_continuity_face
#endif

    subroutine calc_hpe_tend_continuity_2d(ncell, nlev                , &
                                           scalar_flux_div_full_level , &
                                           scalar_ps_tendency         , &
                                           scalar_mass_eta_velocity   )
! io
      use omp_lib
      integer(i4)             , intent(in)    :: ncell, nlev
      real(r8), dimension(:,:), intent(in)    :: scalar_flux_div_full_level ! DEL.(V*delhp)
      real(r8), dimension(:)  , intent(inout) :: scalar_ps_tendency         ! -PARTIAL(PS)/PARTIAL(t)
      real(r8), dimension(:,:), intent(inout) :: scalar_mass_eta_velocity   ! m*etadot at face
! local
      real(r8)                                :: vert_sum
      integer(i4)                             :: iv, ilev
      integer(i4)                             :: ikkk
!
! ps tendency
!
!$omp parallel private(iv,ilev,vert_sum,ikkk)
!$omp do schedule(dynamic,50)
       Do iv = 1, ncell
         vert_sum = zero
         do ilev = 1, nlev
            vert_sum = vert_sum+scalar_flux_div_full_level(ilev,iv)
         end do

         scalar_ps_tendency(iv)  = -1._r8*vert_sum
!
! m*eta velocity
!
         scalar_mass_eta_velocity(1,iv) = zero

         do ilev = 1, nlev
            vert_sum = zero
            do ikkk = 1, ilev
              vert_sum = vert_sum+scalar_flux_div_full_level(ikkk,iv)
            end do
            scalar_mass_eta_velocity(ilev+1,iv) = eta_face_b(ilev+1)*(-1._r8*scalar_ps_tendency(iv))-vert_sum
         end do
!
! this should be, or the above routine is bugged
!
        if (r8 .eq. 8) then
          if( abs(scalar_mass_eta_velocity(nlev+1,iv)).gt.1e-15) print*,"sth wrong inside calc_hpe_tend_continuity_2d m*etadot"
        elseif(r8 .eq. 4) then
          if( abs(scalar_mass_eta_velocity(nlev+1,iv)).gt.1e-7) print*,"sth wrong inside calc_hpe_tend_continuity_2d m*etadot"
        end if

      End do
!$omp end do nowait
!$omp end parallel 

      return
    end subroutine calc_hpe_tend_continuity_2d

 end module grist_hpe_continuity
