module grist_hpe_hydro_pgf_mixed

    use grist_constants,      only: i4, r8, rdry, r4 => ns
    use grist_hpe_constants,  only: p0
    use grist_nml_module,     only: nlev, nlevp
! data
    ! use grist_hpe_constants,  only: eta_face_a, eta_face_b, eta_full_a,eta_full_b
    use grist_dycore_module_vars_control_mixed, only: eta_face_a, eta_face_b, eta_full_a,eta_full_b

    implicit none

contains
!
!  vertical index
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
!
!
    subroutine calc_hpe_hydro(scalar_temperature_full_level  , &
                              scalar_hpressure_face_level    , &
                              scalar_delhp_full_level        , &
                              scalar_geopotential_surface    , &
                              scalar_geopotential_face_level , &
                              scalar_geopotential_full_level )
!
! io
!
      real(r8), dimension(nlev),    intent(in)    :: scalar_temperature_full_level
      real(r8), dimension(nlevp),   intent(in)    :: scalar_hpressure_face_level
      real(r8), dimension(nlev),    intent(in)    :: scalar_delhp_full_level
      real(r8),                     intent(in)    :: scalar_geopotential_surface
      real(r8), dimension(nlevp),   intent(inout) :: scalar_geopotential_face_level
      real(r8), dimension(nlev),    intent(inout) :: scalar_geopotential_full_level
!
! local
!
      real(r8)                             :: vert_sum 
      integer(i4)                          :: ilev
      integer(i4)                          :: ilevel


      do ilevel = 1, nlev     ! for each face ilevel except surface
         vert_sum = 0._r8        
         do ilev = ilevel, nlev ! for each full level ilev
            vert_sum = vert_sum+rdry*scalar_temperature_full_level(ilev)*&
                                 log(scalar_hpressure_face_level(ilev+1)/&
                                     scalar_hpressure_face_level(ilev))
         end do
         scalar_geopotential_face_level(ilevel) = scalar_geopotential_surface+vert_sum
      end do

         scalar_geopotential_face_level(nlev+1) = scalar_geopotential_surface
!
! original hdc use above, comment below for bit reproduce
! 
      do ilev = 1, nlev   ! for each full level
         scalar_geopotential_full_level(ilev) = 0.5_r8*(scalar_geopotential_face_level(ilev+1)+&
                                                        scalar_geopotential_face_level(ilev))
      end do

      return
    end subroutine calc_hpe_hydro

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

    subroutine calc_hpe_get_full_mass(nlev, ncell, working_mode          , & ! in
                                          scalar_hpressure_at_pc_face_level  , & ! in
                                          scalar_hpressure_at_pc_full_level  , & ! in
                                          scalar_delhp_at_pc_full_level      , & ! in
                                          scalar_mif_at_pc_full_level        , & ! in
                                          scalar_pressure_at_pc_full_level   , & ! out
                                          scalar_pressure_at_pc_face_level)      ! out
! io
     integer(i4)        , intent(in)    ::  nlev, ncell
     character(len=*),    intent(in)    ::  working_mode
     real(r4),            intent(in)    ::  scalar_hpressure_at_pc_face_level(:,:) 
     real(r8),            intent(in)    ::  scalar_hpressure_at_pc_full_level(:,:) 
     real(r8),            intent(in)    ::  scalar_delhp_at_pc_full_level(:,:)
     real(r8),            intent(in)    ::  scalar_mif_at_pc_full_level(:,:)
     real(r8),            intent(out)   ::  scalar_pressure_at_pc_full_level(:,:)
     real(r4),            intent(out)   ::  scalar_pressure_at_pc_face_level(:,:)
! local
     integer(i4)  :: iv, ilev, ilevel

     select case(trim(working_mode))

     case('dycore') ! for reference
       scalar_pressure_at_pc_face_level  =  scalar_hpressure_at_pc_face_level
       scalar_pressure_at_pc_full_level  =  scalar_hpressure_at_pc_full_level
       return
     case('tracer') ! for reference
       return
     case default ! this should produce very similar  pressure as hpressure in a dry mode;
                 ! but not bit reproducity due to machine round off;
                 ! at the first top, assume p=hp

       scalar_pressure_at_pc_face_level(1,1:ncell) = scalar_hpressure_at_pc_face_level(1,1:ncell)
       do ilev = 2, nlev+1 ! for each full level ilev
          scalar_pressure_at_pc_face_level(ilev,1:ncell) = scalar_pressure_at_pc_face_level(ilev-1,1:ncell)+&
                                                     scalar_delhp_at_pc_full_level(ilev-1,1:ncell)/scalar_mif_at_pc_full_level(ilev-1,1:ncell)
       end do
       scalar_pressure_at_pc_full_level(1:nlev,1:ncell)  = 0.5_r4*(scalar_pressure_at_pc_face_level(1:nlev,1:ncell)+scalar_pressure_at_pc_face_level(2:nlev+1,1:ncell))

     end select
    
     return
   end subroutine calc_hpe_get_full_mass

    subroutine calc_hpe_get_full_mass1(nlev, ncell, working_mode          , & ! in
                                          scalar_hpressure_at_pc_face_level  , & ! in
                                          scalar_hpressure_at_pc_full_level  , & ! in
                                          scalar_delhp_at_pc_full_level      , & ! in
                                          scalar_mif_at_pc_full_level        , & ! in
                                          scalar_pressure_at_pc_full_level   , & ! out
                                          scalar_pressure_at_pc_face_level)      ! out
! io
     integer(i4)        , intent(in)    ::  nlev, ncell
     character(len=*),    intent(in)    ::  working_mode
     real(r8),            intent(in)    ::  scalar_hpressure_at_pc_face_level(:,:) 
     real(r8),            intent(in)    ::  scalar_hpressure_at_pc_full_level(:,:) 
     real(r8),            intent(in)    ::  scalar_delhp_at_pc_full_level(:,:)
     real(r8),            intent(in)    ::  scalar_mif_at_pc_full_level(:,:)
     real(r8),            intent(out)   ::  scalar_pressure_at_pc_full_level(:,:)
     real(r4),            intent(out)   ::  scalar_pressure_at_pc_face_level(:,:)
! local
     integer(i4)  :: iv, ilev, ilevel

     select case(trim(working_mode))

     case('dycore') ! for reference
       scalar_pressure_at_pc_face_level  =  scalar_hpressure_at_pc_face_level
       scalar_pressure_at_pc_full_level  =  scalar_hpressure_at_pc_full_level
       return
     case('tracer') ! for reference
       return
     case default ! this should produce very similar  pressure as hpressure in a dry mode;
                 ! but not bit reproducity due to machine round off;
                 ! at the first top, assume p=hp

       scalar_pressure_at_pc_face_level(1,1:ncell) = scalar_hpressure_at_pc_face_level(1,1:ncell)
       do ilev = 2, nlev+1 ! for each full level ilev
          scalar_pressure_at_pc_face_level(ilev,1:ncell) = scalar_pressure_at_pc_face_level(ilev-1,1:ncell)+&
                                                     scalar_delhp_at_pc_full_level(ilev-1,1:ncell)/scalar_mif_at_pc_full_level(ilev-1,1:ncell)
       end do
       scalar_pressure_at_pc_full_level(1:nlev,1:ncell)  = 0.5_r4*(scalar_pressure_at_pc_face_level(1:nlev,1:ncell)+scalar_pressure_at_pc_face_level(2:nlev+1,1:ncell))

     end select
    
     return
   end subroutine calc_hpe_get_full_mass1

end module grist_hpe_hydro_pgf_mixed
