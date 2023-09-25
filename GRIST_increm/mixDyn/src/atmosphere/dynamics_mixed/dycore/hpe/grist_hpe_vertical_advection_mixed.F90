module grist_hpe_vertical_advection_mixed

    use grist_constants,    only: i4, r8, half, r4 => ns
    use grist_nml_module,   only: eqs_vert_diff, nlev, nlevp
! data
   !  use grist_hpe_constants,only: deta_face, deta_full
    use grist_dycore_module_vars_control_mixed, only: deta_face, deta_full

    implicit none

contains
!
! only used by uwind vert adv
!

   subroutine calc_hpe_vert_advection( scalar_mass_eta_velocity_face_level,&
                                       scalar_delhp_full_level            ,&
                                       scalar_at_full_level    ,&
                                       scalar_vert_advection_full_level)
   implicit none
!
! io
!
     real(r8), dimension(nlevp),   intent(in)    :: scalar_mass_eta_velocity_face_level ! m*etadot 
     real(r8), dimension(nlev),    intent(in)    :: scalar_delhp_full_level
     real(r8), dimension(nlev),    intent(in)    :: scalar_at_full_level
     real(r8), dimension(nlev),    intent(inout) :: scalar_vert_advection_full_level    ! with minus sign
!
! local
!
     integer(i4)                          :: ilev
     real(r8), dimension(nlevp)           :: scalar_at_face_level
!     real(r8), dimension(nlevp)           :: scalar_eta_velocity

!     scalar_eta_velocity  = scalar_mass_eta_velocity_face_level

!     if(.not.allocated(scalar_at_face_level)) allocate(scalar_at_face_level(nlev+1))
!
! a more generalized way
!
     call calc_hpe_vert_flux_operator1(scalar_at_full_level   ,&
                                  scalar_mass_eta_velocity_face_level    ,&
                                  2                      ,&
                                  scalar_at_face_level   ,&
                                  .true.)

     do ilev = 1, nlev

        scalar_vert_advection_full_level(ilev) =-(scalar_mass_eta_velocity_face_level(ilev+1)*&
                                       (scalar_at_face_level(ilev+1)-&
                                        scalar_at_full_level(ilev))+&
                                        scalar_mass_eta_velocity_face_level(ilev)*&
                                       (scalar_at_full_level(ilev)-&
                                        scalar_at_face_level(ilev)))/scalar_delhp_full_level(ilev)
     end do
!
! directly follow that way in SS81, same as above to machine precision when using eqs_vert_diff
!
     !do ilev = 1, nlev
     !   scalar_vert_advection(ilev) =-(scalar_mass_eta_velocity(ilev+1)*&
     !                                    (scalar_at_full_level(ilev+1)-&
     !                                     scalar_at_full_level(ilev))+&
     !                                    scalar_mass_eta_velocity(ilev)*&
     !                                    (scalar_at_full_level(ilev)-&
     !                                     scalar_at_full_level(ilev-1)))/(2._r8*scalar_delhp(ilev))
       ! if(ilev.eq.10) print*,scalar_vert_advection(ilev)
     !end do

!     deallocate(scalar_at_face_level)
!     deallocate(scalar_eta_velocity)

     return 
   end subroutine calc_hpe_vert_advection

   subroutine calc_hpe_tend_vert_mass_flux_2d(ncell, ncell_do, nlev, &
                                              scalar_at_full_level , &
                                              scalar_eta_velocity  , &
                                              order                , &
                                              tend_vert_mass_flux)
! io
    use omp_lib
    integer(i4),  intent(in)                         :: ncell, ncell_do, nlev
    real(r8), dimension(nlev,  ncell), intent(in)    :: scalar_at_full_level
#ifdef DBL_MASS
    real(r8), dimension(nlevp, ncell), intent(in)    :: scalar_eta_velocity  ! m*etadot
#else
    real(r4), dimension(nlevp, ncell), intent(in)    :: scalar_eta_velocity  ! m*etadot
#endif
    integer(i4)                      , intent(in)    :: order
    real(r8), dimension(nlev,ncell)  , intent(inout) :: tend_vert_mass_flux  ! with minus sign
! local
    real(r8), dimension(nlevp,ncell)                 :: scalar_at_face_level
    integer(i4)  :: icell, ilev

!$omp parallel private(icell, ilev)
!$omp do schedule(dynamic,5) 
     do icell = 1, ncell_do
        call calc_hpe_vert_flux_operator(scalar_at_full_level(1:nlev,icell)   ,&
                                         scalar_eta_velocity(1:nlevp,icell)    ,&
                                         order                  ,&
                                         scalar_at_face_level(1:nlevp,icell)   )
        do ilev =  1, nlev
           tend_vert_mass_flux(ilev,icell) = -(scalar_eta_velocity(ilev+1,icell)*scalar_at_face_level(ilev+1,icell)-&
                                               scalar_eta_velocity(ilev,icell)  *scalar_at_face_level(ilev,icell))

        end do
     end do
!$omp end do nowait
!$omp end parallel 

     return
   end subroutine calc_hpe_tend_vert_mass_flux_2d

!
! from full level to face level
!

   subroutine calc_hpe_vert_flux_operator(scalar_at_full_level    ,&
                                      scalar_mass_eta_velocity,&
                                      order                   ,&
                                      scalar_at_face_level    ,&
                                      called_by_uwind)

    real(r8), dimension(nlev),   intent(in)     :: scalar_at_full_level
#ifdef DBL_MASS
    real(r8), dimension(nlevp),  intent(in)     :: scalar_mass_eta_velocity
#else
    real(r4), dimension(nlevp),  intent(in)     :: scalar_mass_eta_velocity
#endif
    integer                   ,  intent(in)     :: order
    real(r8), dimension(nlevp),  intent(inout)  :: scalar_at_face_level
    logical, optional ,          intent(in)     :: called_by_uwind
! 
    integer(i4)                                 :: ilev
    real(r8)                                    :: part1, der_k, der_kp1
    logical                                     :: local_flag
!
! extrapolate, actually not used because etadot at boundaries are zero
! 
     scalar_at_face_level(1)        = scalar_at_full_level(1)
     scalar_at_face_level(nlev+1)   = scalar_at_full_level(nlev)

     if(present(called_by_uwind)) then
        local_flag = called_by_uwind
     else
        local_flag = .false.
     end if

     IF(EQS_VERT_DIFF.or.local_flag)THEN
!
! This is equivalent distance version, used as default
!
     select case(order)
     case(2) 
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 2, nlev
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+&
                                         scalar_at_full_level(ilev-1))*half
        end do

     case(3) 
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(2)       = (scalar_at_full_level(2)+scalar_at_full_level(1))*half
        scalar_at_face_level(nlev)    = (scalar_at_full_level(nlev)+scalar_at_full_level(nlev-1))*half
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 3, nlev-1
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+scalar_at_full_level(ilev-1))*(7._r8/12._r8)-&
                                          (scalar_at_full_level(ilev+1)+scalar_at_full_level(ilev-2))*(1._r8/12._r8)+&
                                sign(1._r8,scalar_mass_eta_velocity(ilev))*&
                                         ((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev-2))-&
                                    3._r8*(scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1)))/12._r8
        end do

     case default
        print*,"you must select a vertical order, stop"
        stop
     end select

     ELSE  ! non-equivalent distance version

     select case(order)
     case(2)
        do ilev = 2, nlev
           scalar_at_face_level(ilev)= (scalar_at_full_level(ilev)*deta_full(ilev-1)+&
                                          scalar_at_full_level(ilev-1)*deta_full(ilev))*half/deta_face(ilev)
        end do
     case(3) 
        scalar_at_face_level(2)      = (scalar_at_full_level(2)*deta_full(1)+scalar_at_full_level(1)*deta_full(2))*half/deta_face(2)
        scalar_at_face_level(nlev)   = (scalar_at_full_level(nlev)*deta_full(nlev-1)+scalar_at_full_level(nlev-1)*deta_full(nlev))*half/deta_face(nlev)
        do ilev = 3, nlev-1 ! face level
! face level
           part1 = (scalar_at_full_level(ilev)*deta_full(ilev-1)+&
                    scalar_at_full_level(ilev-1)*deta_full(ilev))*half/deta_face(ilev)
! full level for this face level
           der_k  = 2._r8*(deta_full(ilev-1)**2)*((scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))/deta_face(ilev)-&
                                                 (scalar_at_full_level(ilev-1)-scalar_at_full_level(ilev-2))/deta_face(ilev-1))/&
                                                 (deta_face(ilev)+deta_face(ilev-1))
           der_kp1= 2._r8*(deta_full(ilev)**2)*((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev))/deta_face(ilev+1)-&
                                                 (scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))/deta_face(ilev))/&
                                                 (deta_face(ilev+1)+deta_face(ilev))
           scalar_at_face_level(ilev)= part1-(der_k+der_kp1)/12._r8+sign(1._r8,scalar_mass_eta_velocity(ilev))*(der_kp1-der_k)
        end do
     case default
        print*, "you must select a vertical order, stop"
        stop
     end select

     END IF

     return
   end subroutine calc_hpe_vert_flux_operator

   subroutine calc_hpe_vert_flux_operator1(scalar_at_full_level    ,&
                                      scalar_mass_eta_velocity,&
                                      order                   ,&
                                      scalar_at_face_level    ,&
                                      called_by_uwind)

    real(r8), dimension(nlev),   intent(in)     :: scalar_at_full_level
    real(r8), dimension(nlevp),  intent(in)     :: scalar_mass_eta_velocity
    integer                   ,  intent(in)     :: order
    real(r8), dimension(nlevp),  intent(inout)  :: scalar_at_face_level
    logical, optional ,          intent(in)     :: called_by_uwind
! 
    integer(i4)                                 :: ilev
    real(r8)                                    :: part1, der_k, der_kp1
    logical                                     :: local_flag
!
! extrapolate, actually not used because etadot at boundaries are zero
! 
     scalar_at_face_level(1)        = scalar_at_full_level(1)
     scalar_at_face_level(nlev+1)   = scalar_at_full_level(nlev)

     if(present(called_by_uwind)) then
        local_flag = called_by_uwind
     else
        local_flag = .false.
     end if

     IF(EQS_VERT_DIFF.or.local_flag)THEN
!
! This is equivalent distance version, used as default
!
     select case(order)
     case(2) 
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 2, nlev
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+&
                                         scalar_at_full_level(ilev-1))*half
        end do

     case(3) 
        scalar_at_face_level(1)       = 0._r8
        scalar_at_face_level(2)       = (scalar_at_full_level(2)+scalar_at_full_level(1))*half
        scalar_at_face_level(nlev)    = (scalar_at_full_level(nlev)+scalar_at_full_level(nlev-1))*half
        scalar_at_face_level(nlev+1)  = 0._r8
        do ilev = 3, nlev-1
           scalar_at_face_level(ilev) = (scalar_at_full_level(ilev)+scalar_at_full_level(ilev-1))*(7._r8/12._r8)-&
                                          (scalar_at_full_level(ilev+1)+scalar_at_full_level(ilev-2))*(1._r8/12._r8)+&
                                sign(1._r8,scalar_mass_eta_velocity(ilev))*&
                                         ((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev-2))-&
                                    3._r8*(scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1)))/12._r8
        end do

     case default
        print*,"you must select a vertical order, stop"
        stop
     end select

     ELSE  ! non-equivalent distance version

     select case(order)
     case(2)
        do ilev = 2, nlev
           scalar_at_face_level(ilev)= (scalar_at_full_level(ilev)*deta_full(ilev-1)+&
                                          scalar_at_full_level(ilev-1)*deta_full(ilev))*half/deta_face(ilev)
        end do
     case(3) 
        scalar_at_face_level(2)      = (scalar_at_full_level(2)*deta_full(1)+scalar_at_full_level(1)*deta_full(2))*half/deta_face(2)
        scalar_at_face_level(nlev)   = (scalar_at_full_level(nlev)*deta_full(nlev-1)+scalar_at_full_level(nlev-1)*deta_full(nlev))*half/deta_face(nlev)
        do ilev = 3, nlev-1 ! face level
! face level
           part1 = (scalar_at_full_level(ilev)*deta_full(ilev-1)+&
                    scalar_at_full_level(ilev-1)*deta_full(ilev))*half/deta_face(ilev)
! full level for this face level
           der_k  = 2._r8*(deta_full(ilev-1)**2)*((scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))/deta_face(ilev)-&
                                                 (scalar_at_full_level(ilev-1)-scalar_at_full_level(ilev-2))/deta_face(ilev-1))/&
                                                 (deta_face(ilev)+deta_face(ilev-1))
           der_kp1= 2._r8*(deta_full(ilev)**2)*((scalar_at_full_level(ilev+1)-scalar_at_full_level(ilev))/deta_face(ilev+1)-&
                                                 (scalar_at_full_level(ilev)-scalar_at_full_level(ilev-1))/deta_face(ilev))/&
                                                 (deta_face(ilev+1)+deta_face(ilev))
           scalar_at_face_level(ilev)= part1-(der_k+der_kp1)/12._r8+sign(1._r8,scalar_mass_eta_velocity(ilev))*(der_kp1-der_k)
        end do
     case default
        print*, "you must select a vertical order, stop"
        stop
     end select

     END IF

     return
   end subroutine calc_hpe_vert_flux_operator1

end module grist_hpe_vertical_advection_mixed