module grist_hpe_continuity_mixed

    use grist_constants,        only: i4, r8, zero, r4 => ns
! data
#ifdef DBL_MASS
    use grist_hpe_constants,    only: eta_face_b_r8 => eta_face_b
#endif
    use grist_dycore_module_vars_control_mixed, only: eta_face_b

    implicit none

contains

    subroutine calc_hpe_tend_continuity_2d(ncell , nlev               , &
                                           scalar_flux_div_full_level , &
                                           scalar_ps_tendency         , &
                                           scalar_mass_eta_velocity)
! io
      use omp_lib
      integer(i4)             , intent(in)    :: ncell, nlev
#ifdef DBL_MASS
      real(r8), dimension(:,:), intent(in)    :: scalar_flux_div_full_level ! DEL.(V*delhp)
      real(r8), dimension(:)  , intent(inout) :: scalar_ps_tendency         ! -PARTIAL(PS)/PARTIAL(t)
      real(r8), dimension(:,:), intent(inout) :: scalar_mass_eta_velocity   ! m*etadot at face
#else
      real(r4), dimension(:,:), intent(in)    :: scalar_flux_div_full_level ! DEL.(V*delhp)
      real(r4), dimension(:)  , intent(inout) :: scalar_ps_tendency         ! -PARTIAL(PS)/PARTIAL(t)
      real(r4), dimension(:,:), intent(inout) :: scalar_mass_eta_velocity   ! m*etadot at face
#endif

      real(r8)                          :: vert_sum
      integer(i4)                       :: iv, ilev
      integer(i4)                       :: ikkk

!
! ps tendency (3.2 of SS81)
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
! m*eta velocity (3.3 of SS81)
!
         scalar_mass_eta_velocity(1,iv) = zero

         do ilev = 1, nlev
            vert_sum = zero
            do ikkk = 1, ilev
              vert_sum = vert_sum+scalar_flux_div_full_level(ikkk,iv)
            end do
#ifdef DBL_MASS
            scalar_mass_eta_velocity(ilev+1,iv) = eta_face_b_r8(ilev+1)*&
                                              (-1._r8*scalar_ps_tendency(iv))-vert_sum
#else
            scalar_mass_eta_velocity(ilev+1,iv) = eta_face_b(ilev+1)*&
                                              (-1._r4*scalar_ps_tendency(iv))-real(vert_sum,r4)
#endif
         end do
!
! this should be, or the above routine is bugged
!
        if (r4 .eq. 8) then
          ! if( abs(scalar_mass_eta_velocity(nlev+1,iv)).gt.1e-15)then
          if( abs(scalar_mass_eta_velocity(nlev+1,iv)).gt.1e-13)then
            print*,"sth wrong inside calc_hpe_tend_continuity_2d m*etadot",r4,abs(scalar_mass_eta_velocity(nlev+1,iv))
          endif
        elseif(r4 .eq. 4) then
          if( abs(scalar_mass_eta_velocity(nlev+1,iv)).gt.1e-7)then
            print*,"sth wrong inside calc_hpe_tend_continuity_2d m*etadot",r4,abs(scalar_mass_eta_velocity(nlev+1,iv))
          end if
        end if

      End do
!$omp end do nowait
!$omp end parallel 

      return
    end subroutine calc_hpe_tend_continuity_2d

end module grist_hpe_continuity_mixed