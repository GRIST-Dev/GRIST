
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: vertical advection
! Revision history:
!----------------------------------------------------------------------------

 module grist_hpe_vertical_advection

   use grist_constants,      only: i4, r8, zero
   use grist_hpe_constants,  only: deta_face, deta_full, eta_face_a, eta_face_b
   use grist_nml_module,     only: eqs_vert_diff

!
!               vertical index
!-(1)--------------------------------- 1/2
!  ..............[1]..................
!-(2)--------------------------------- 1+1/2
!  ..............[2]..................
!-(3)--------------------------------- 2+1/2
!  ..............[3]..................
!-(4)--------------------------------- 3+1/2
!  ..............[4]..................
!-(5)--------------------------------- 4+1/2
!  ..............[5]..................
!-(6)------------------ilevel-1------- 5+1/2
!  ..............[6]...ilevel-1.......
!-(7)------------------ilevel--------- 6+1/2
!  ..............[7]...ilevel.........
!-(8)------------------ilevel+1------- 7+1/2
!  ..............[8]...ilevel+1.......
!-(9)------------------ilevel+2------- 8+1/2
! the positive direction of eta velocity is also
! downward positve
!

   implicit none

   private

   public :: calc_hpe_vert_advection      , &
             calc_hpe_tend_vert_mass_flux_2d
#ifdef REFCODE
            &,calc_hpe_tend_vert_mass_flux
            &,calc_hpe_tend_vert_mass_flux_face
#endif

   real(r8), parameter  :: half = 0.5_r8

   contains

!-----------------------------------------------------------------------
! only used by uwind vert adv
!-----------------------------------------------------------------------

   subroutine calc_hpe_vert_advection(nlev, nlevp                        ,&
                                      scalar_mass_eta_velocity_face_level,&
                                      scalar_delhp_full_level            ,&
                                      scalar_at_full_level               ,&
                                      scalar_vert_advection_full_level)
   implicit none
!
! io
!
     integer(i4),intent(in)    :: nlev, nlevp
     real(r8),   intent(in)    :: scalar_mass_eta_velocity_face_level(:) ! m*etadot 
     real(r8),   intent(in)    :: scalar_delhp_full_level(:)
     real(r8),   intent(in)    :: scalar_at_full_level(:)
     real(r8),   intent(inout) :: scalar_vert_advection_full_level(:)    ! with minus sign
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
     call calc_hpe_vert_flux_operator(nlev, nlevp, scalar_at_full_level   ,&
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

!
! partial m*etadot[]
!
#ifdef REFCODE
   subroutine calc_hpe_tend_vert_mass_flux(nlev, nlevp,scalar_at_full_level , &
                                           scalar_eta_velocity  , &
                                           order                , &
                                           tend_vert_mass_flux)
! io
    integer(i4),                intent(in)    :: nlev, nlevp
    real(r8), dimension(nlev),  intent(in)    :: scalar_at_full_level
    real(r8), dimension(nlevp), intent(in)    :: scalar_eta_velocity  ! m*etadot
    integer(i4)              ,  intent(in)    :: order
    real(r8), dimension(nlev),  intent(inout) :: tend_vert_mass_flux  ! with minus sign
! local
    real(r8), dimension(nlevp)                :: scalar_at_face_level
    integer(i4)                               :: ilev


     call calc_hpe_vert_flux_operator(nlev, nlevp,scalar_at_full_level   ,&
                                  scalar_eta_velocity    ,&
                                  order                  ,&
                                  scalar_at_face_level   )

     do ilev =  1, nlev

        tend_vert_mass_flux(ilev) = -(scalar_eta_velocity(ilev+1)*scalar_at_face_level(ilev+1)-&
                                        scalar_eta_velocity(ilev)*scalar_at_face_level(ilev))

     end do

     return
   end subroutine calc_hpe_tend_vert_mass_flux
#endif

   subroutine calc_hpe_tend_vert_mass_flux_2d(ncell, ncell_do, nlev, nlevp, &
                                              scalar_at_full_level , &
                                              scalar_eta_velocity  , &
                                              order                , &
                                              tend_vert_mass_flux)
! io
    use omp_lib
    integer(i4),  intent(in)                         :: ncell, ncell_do, nlev, nlevp
    real(r8), dimension(nlev,  ncell), intent(in)    :: scalar_at_full_level
    real(r8), dimension(nlevp, ncell), intent(in)    :: scalar_eta_velocity  ! m*etadot
    integer(i4)                      , intent(in)    :: order
    real(r8), dimension(nlev,ncell)  , intent(inout) :: tend_vert_mass_flux  ! with minus sign
! local
    real(r8), dimension(nlevp,ncell)                 :: scalar_at_face_level
    integer(i4)  :: icell, ilev

!$omp parallel private(icell, ilev)
!$omp do schedule(dynamic,5) 
     do icell = 1, ncell_do
        call calc_hpe_vert_flux_operator(nlev, nlevp,scalar_at_full_level(1:nlev,icell)   ,&
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

   subroutine calc_hpe_vert_flux_operator(nlev, nlevp, scalar_at_full_level    ,&
                                      scalar_mass_eta_velocity,&
                                      order                   ,&
                                      scalar_at_face_level    ,&
                                      called_by_uwind)

    integer(i4),                 intent(in)     :: nlev, nlevp
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
   end subroutine calc_hpe_vert_flux_operator
!
! BELOW FOR FACE LEVEL
!
!-----------------------------------------------------------------------
! partial m*etadot
! NHD now has bespoke version
!-----------------------------------------------------------------------

#ifdef REFCODE
   subroutine calc_hpe_tend_vert_mass_flux_face(nlev, nlevp, scalar_at_face_level      , &
                                                scalar_eta_mass_velocity_full  , &
                                                order                     , &
                                                tend_vert_mass_flux_face)
! io
    integer(i4),                 intent(in)   :: nlev, nlevp
    real(r8), dimension(nlevp),  intent(in)   :: scalar_at_face_level
    real(r8), dimension(nlev),   intent(in)   :: scalar_eta_mass_velocity_full  ! m*etadot
    integer(i4)               ,  intent(in)   :: order
    real(r8), dimension(nlevp),  intent(inout):: tend_vert_mass_flux_face  ! with minus sign
! local
    real(r8), dimension(nlev)                 :: scalar_at_full_level
    integer(i4)                               :: ilev

     call calc_hpe_vert_flux_operator_face(nlev, nlevp, scalar_at_face_level     ,&
                                       scalar_eta_mass_velocity_full ,&
                                       order                    ,&
                                       scalar_at_full_level )

     do ilev =  2, nlev
! no minus sign here
        tend_vert_mass_flux_face(ilev) = (  scalar_eta_mass_velocity_full(ilev)*scalar_at_full_level(ilev)-&
                                            scalar_eta_mass_velocity_full(ilev-1)*scalar_at_full_level(ilev-1))

     end do

     tend_vert_mass_flux_face(1)      = 0._r8
     tend_vert_mass_flux_face(nlev+1) = 0._r8

     return
   end subroutine calc_hpe_tend_vert_mass_flux_face

!
! private, from face level to full level
!

   subroutine calc_hpe_vert_flux_operator_face(nlev, nlevp, &
                                               scalar_at_face_level    ,&
                                               scalar_eta_mass_velocity_full,&
                                               order                   ,&
                                               scalar_at_full_level    )

    integer(i4),                intent(in)     :: nlev, nlevp
    real(r8), dimension(nlevp), intent(in)     :: scalar_at_face_level
    real(r8), dimension(nlev),  intent(in)     :: scalar_eta_mass_velocity_full
    integer                  ,  intent(in)     :: order
    real(r8), dimension(nlev),  intent(inout)  :: scalar_at_full_level
! local
    integer(i4)                         :: ilev
    real(r8)                            :: part1, der_k, der_kp1
    
     scalar_at_full_level(1)      = half*(scalar_at_face_level(1)+scalar_at_face_level(2))
     scalar_at_full_level(nlev)   = half*(scalar_at_face_level(nlev)+scalar_at_face_level(nlev+1))

     IF(EQS_VERT_DIFF)THEN ! default

     select case(order)
     case(2) 
        do ilev = 2, nlev-1
           scalar_at_full_level(ilev) = half*(scalar_at_face_level(ilev)+scalar_at_face_level(ilev+1))
        end do
     case(3) 
        do ilev = 2, nlev-1
           scalar_at_full_level(ilev) = (scalar_at_face_level(ilev)+scalar_at_face_level(ilev+1))*(7._r8/12._r8)-&
                                          (scalar_at_face_level(ilev+2)+scalar_at_face_level(ilev-1))*(1._r8/12._r8)+&
                                sign(1._r8,scalar_eta_mass_velocity_full(ilev))*&
                                         ((scalar_at_face_level(ilev+2)-scalar_at_face_level(ilev-1))-&
                                    3._r8*(scalar_at_face_level(ilev+1)-scalar_at_face_level(ilev)))/12._r8
        end do
     case default
        print*," you must select a vert order in calc_hpe_vert_flux_operator_face"
        stop
     end select

     ELSE
     select case(order)
     case(2) 
        do ilev = 2, nlev-1
           scalar_at_full_level(ilev) = half*(scalar_at_face_level(ilev)+scalar_at_face_level(ilev+1))
        end do
     case(3)
        do ilev = 2, nlev-1
           part1 = half*(scalar_at_face_level(ilev)+scalar_at_face_level(ilev+1))

           der_k = 2._r8*(deta_face(ilev)**2)*((scalar_at_face_level(ilev+1)-scalar_at_face_level(ilev))/deta_full(ilev)-&
                                               (scalar_at_face_level(ilev)-scalar_at_face_level(ilev-1))/deta_full(ilev-1))/&
                                               (deta_full(ilev)+deta_full(ilev-1))

           der_kp1 = 2._r8*(deta_face(ilev+1)**2)*((scalar_at_face_level(ilev+2)-scalar_at_face_level(ilev+1))/deta_full(ilev+1)-&
                                                   (scalar_at_face_level(ilev+1)-scalar_at_face_level(ilev))/deta_full(ilev))/&
                                                   (deta_full(ilev+1)+deta_full(ilev))

           scalar_at_full_level(ilev)= part1-(der_k+der_kp1)/12._r8+sign(1._r8,scalar_eta_mass_velocity_full(ilev))*(der_kp1-der_k)
        end do
     case default
        print*," you must select a vert order in calc_hpe_vert_flux_operator_face"
        stop
     end select

     END IF
     
     return
   end subroutine calc_hpe_vert_flux_operator_face
#endif

 end module grist_hpe_vertical_advection
