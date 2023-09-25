!======================================================
!
!  Created by LiXiaohan on 19/5/20.
!  
!  Purpose:
!  update pressure at full/face level
!  
!  adopted from 
!  atmosphere/dycore/grist_dycore_time_integration_2d.F90
!======================================================

 module grist_scm_pressure_update

    use grist_constants,            only: i4, r8
    use grist_nml_module,           only: nlev, nlevp
    use grist_hpe_hydro_pgf,        only: calc_hpe_hpressure_face_level,            &
                                          calc_hpe_hpressure_full_level,            &
                                          calc_hpe_delhp
 
    implicit none

    private

    public  :: time_integration_renew_mass_state

 contains
 
    subroutine time_integration_renew_mass_state(ncol, hpressure_at_pc_surface,     &
                                                 hpressure_at_pc_face_level,        &
                                                 delhp_at_pc_full_level,            &
                                                 hpressure_at_pc_full_level,        &
                                                 delhp_at_pc_face_level     )


! io
    integer,   intent(in)    :: ncol
    real(r8),  intent(in)    :: hpressure_at_pc_surface(ncol)
    real(r8),  intent(inout) :: hpressure_at_pc_face_level(nlevp, ncol)
    real(r8),  intent(inout) :: delhp_at_pc_full_level(nlev, ncol)
    real(r8),  intent(inout) :: hpressure_at_pc_full_level(nlev, ncol)
    real(r8),  intent(inout) :: delhp_at_pc_face_level(nlevp, ncol)
! local
    integer(i4)                          :: iv, ilev
    real(r8)                             :: scalar_template_a
    real(r8), allocatable                :: scalar_template_1d_nlevp_a(:), &
                                            scalar_template_1d_nlev_a(:),  &
                                            scalar_template_1d_nlev_b(:)

    allocate(scalar_template_1d_nlev_a(nlev))
    allocate(scalar_template_1d_nlev_b(nlev))
    allocate(scalar_template_1d_nlevp_a(nlevp))

!
! renew mass state, compute hpressure at face level, full level, 
! delhp, based on input surface hpressure
!
     
    do iv = 1, ncol
        scalar_template_a   =   hpressure_at_pc_surface(iv)

        call calc_hpe_hpressure_face_level(scalar_template_a, &       ! in
                                           scalar_template_1d_nlevp_a) ! out

        call calc_hpe_delhp(scalar_template_1d_nlevp_a, &  ! in
                            scalar_template_1d_nlev_a)     ! out

        call calc_hpe_hpressure_full_level(scalar_template_a         , & ! in
                                           scalar_template_1d_nlevp_a, & ! in
                                           scalar_template_1d_nlev_a , & ! in
                                           scalar_template_1d_nlev_b)    ! out

        hpressure_at_pc_face_level(:,iv)  = scalar_template_1d_nlevp_a
        delhp_at_pc_full_level(:,iv)  = scalar_template_1d_nlev_a
        hpressure_at_pc_full_level(:,iv)  = scalar_template_1d_nlev_b

        do ilev = 2, nlev
            delhp_at_pc_face_level(ilev,iv) = hpressure_at_pc_full_level(ilev,iv)-&
                                              hpressure_at_pc_full_level(ilev-1,iv)
        end do

        delhp_at_pc_face_level(1,iv)      = 0.5_r8*delhp_at_pc_full_level(1,iv)
        delhp_at_pc_face_level(nlevp,iv) = 0.5_r8*delhp_at_pc_full_level(nlev,iv)
    end do

        deallocate(scalar_template_1d_nlev_a)
        deallocate(scalar_template_1d_nlev_b)
        deallocate(scalar_template_1d_nlevp_a)

        return
   end subroutine time_integration_renew_mass_state


 end module grist_scm_pressure_update
