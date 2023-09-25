module grist_tracer_hpe_hydro_pgf_mixed

    use grist_constants,                        only: r8, i4
    use grist_hpe_constants,                    only: p0
    use grist_nml_module,                       only: nlev, nlevp
! data
    use grist_tracer_module_vars_control_mixed, only: eta_face_a, eta_face_b, eta_full_a, eta_full_b

    implicit none

contains

    subroutine calc_hpe_hpressure_face_level(ps,scalar_hpressure_face_level)

        real(r8),                   intent(in)     :: ps
        real(r8), dimension(nlevp), intent(inout)  :: scalar_hpressure_face_level
        integer(i4)                        :: ilev

        do ilev = 1, nlev+1
            scalar_hpressure_face_level(ilev) = eta_face_a(ilev)*p0+eta_face_b(ilev)*ps
        end do

        return
    end subroutine calc_hpe_hpressure_face_level

    subroutine calc_hpe_delhp(scalar_hpressure_face_level,scalar_delhp_full_level)

        real(r8), dimension(nlevp),  intent(in)    :: scalar_hpressure_face_level
        real(r8), dimension(nlev),  intent(inout)  :: scalar_delhp_full_level
        integer(i4)   :: ilev

            do ilev = 1, nlev     ! for each full level
            scalar_delhp_full_level(ilev) = scalar_hpressure_face_level(ilev+1)-&
                                            scalar_hpressure_face_level(ilev) 
            end do
        return
    end subroutine calc_hpe_delhp

    subroutine calc_hpe_hpressure_full_level(ps, &
                                             scalar_hpressure_face_level,&
                                             scalar_delhp_full_level    ,&
                                             scalar_hpressure_full_level)

        real(r8),                  intent(in)     ::  ps
        real(r8), dimension(nlevp),intent(in)     ::  scalar_hpressure_face_level
        real(r8), dimension(nlev), intent(in)     ::  scalar_delhp_full_level
        real(r8), dimension(nlev), intent(inout)  ::  scalar_hpressure_full_level
!
! local
!
        integer(i4)                        ::  ilev
!
! 3.17 of SS81
!
        scalar_hpressure_full_level(1) = 0.5_r8*scalar_delhp_full_level(1)

        do ilev = 2, nlev
            scalar_hpressure_full_level(ilev) = scalar_delhp_full_level(ilev)/&
                                            log(scalar_hpressure_face_level(ilev+1)/&
                                                scalar_hpressure_face_level(ilev))
        end do
!
! original hdc version use above, comment below for bit reproduce of hdc
!
        do ilev = 1, nlev
            scalar_hpressure_full_level(ilev) = eta_full_a(ilev)*p0+eta_full_b(ilev)*ps
        end do

        return
    end subroutine calc_hpe_hpressure_full_level

end module grist_tracer_hpe_hydro_pgf_mixed
