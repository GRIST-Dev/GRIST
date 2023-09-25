module grist_tracer_hpe_continuity_mixed

    use grist_constants,                        only: i4, zero_r8 => zero, zero => zero_ns, r8, r4 => ns
    
    use grist_nml_module,                       only: nlev, nlevp
! data
#ifdef DBL_MASS_TRACER
    use grist_hpe_constants,                    only: eta_face_b
#else
    use grist_tracer_module_vars_control_mixed, only: eta_face_b
#endif

    implicit none

contains

    subroutine calc_hpe_tend_continuity_2d(ncell , nlev               , &
                                           scalar_flux_div_full_level , &
                                           scalar_ps_tendency         , &
                                           scalar_mass_eta_velocity   )
! io
      use omp_lib
      integer(i4)             , intent(in)    :: ncell, nlev
#ifdef DBL_MASS_TRACER
      real(r8), dimension(:,:), intent(in)    :: scalar_flux_div_full_level ! DEL.(V*delhp)
      real(r8), dimension(:)  , intent(inout) :: scalar_ps_tendency         ! -PARTIAL(PS)/PARTIAL(t)
      real(r8), dimension(:,:), intent(inout) :: scalar_mass_eta_velocity   ! m*etadot at face
      real(r8)                                :: vert_sum
#else
      real(r4), dimension(:,:), intent(in)    :: scalar_flux_div_full_level ! DEL.(V*delhp)
      real(r4), dimension(:)  , intent(inout) :: scalar_ps_tendency         ! -PARTIAL(PS)/PARTIAL(t)
      real(r4), dimension(:,:), intent(inout) :: scalar_mass_eta_velocity   ! m*etadot at face
      real(r4)                                :: vert_sum
#endif
      integer(i4)                             :: iv, ilev
      integer(i4)                             :: ikkk
!
! ps tendency
!
!$omp parallel private(iv,ilev,vert_sum,ikkk)
!$omp do schedule(dynamic,50)
       Do iv = 1, ncell
#ifdef DBL_MASS_TRACER
         vert_sum = zero_r8
#else
         vert_sum = zero
#endif
         do ilev = 1, nlev
            vert_sum = vert_sum+scalar_flux_div_full_level(ilev,iv)
         end do

#ifdef DBL_MASS_TRACER
         scalar_ps_tendency(iv)  = -1._r8*vert_sum
#else
         scalar_ps_tendency(iv)  = -1._r4*vert_sum
#endif
!
! m*eta velocity
!
#ifdef DBL_MASS_TRACER
         scalar_mass_eta_velocity(1,iv) = zero_r8
#else
         scalar_mass_eta_velocity(1,iv) = zero
#endif

         do ilev = 1, nlev
#ifdef DBL_MASS_TRACER
            vert_sum = zero_r8
#else
            vert_sum = zero
#endif
            do ikkk = 1, ilev
              vert_sum = vert_sum+scalar_flux_div_full_level(ikkk,iv)
            end do
#ifdef DBL_MASS_TRACER
            scalar_mass_eta_velocity(ilev+1,iv) = eta_face_b(ilev+1)*(-1._r8*scalar_ps_tendency(iv))-vert_sum
#else
            scalar_mass_eta_velocity(ilev+1,iv) = eta_face_b(ilev+1)*(-1._r4*scalar_ps_tendency(iv))-vert_sum
#endif
         end do
!
! this should be, or the above routine is bugged
!
        if (r4 .eq. 8) then
          if( abs(scalar_mass_eta_velocity(nlev+1,iv)).gt.1e-15) print*,"sth wrong inside calc_hpe_tend_continuity_2d m*etadot"
        elseif(r4 .eq. 4) then
          if( abs(scalar_mass_eta_velocity(nlev+1,iv)).gt.1e-7) print*,"sth wrong inside calc_hpe_tend_continuity_2d m*etadot"
        end if

      End do
!$omp end do nowait
!$omp end parallel 

      return
    end subroutine calc_hpe_tend_continuity_2d

end module grist_tracer_hpe_continuity_mixed