 !======================================================
 !         NAMELIST module for WRF-based physics
 !======================================================

 module grist_PhysW_nml_module

   use grist_constants,   only: i4, r8
   use grist_mpi
   use grist_nml_module,  only: model_timestep

   implicit none
   
   private

   public  :: set_PhysW_nml,   &
!
! WRF physics
!
              PhysW_cu_scheme,  &
              PhysW_mp_scheme,  &
              PhysW_bl_scheme,  &
              PhysW_cf_scheme,  &
              PhysW_ra_scheme,  &
              PhysW_rasw_scheme,&
              PhysW_ralw_scheme,&
              PhysW_sf_scheme,  &
              PhysW_lm_scheme,  &
              step_cu,            &
              step_ra,            &
              step_bl,            &
              wphys_has_req,      &
              use_gwdo,           &
              use_cond,           &
              unuse_cu,           &
              force_read_thompson,&
              write_thompson_tables,&
              para1_zoentr_test    ,&
              para2_zoentr_test    ,&
              para1_eta_test       ,&
              para2_eta_test
!
! WRF physics package related
!
   character(100)     :: PhysW_cu_scheme
   character(100)     :: PhysW_mp_scheme
   character(100)     :: PhysW_bl_scheme
   character(100)     :: PhysW_cf_scheme
   character(100)     :: PhysW_ra_scheme
   character(100)     :: PhysW_rasw_scheme
   character(100)     :: PhysW_ralw_scheme
   character(100)     :: PhysW_sf_scheme
   character(100)     :: PhysW_lm_scheme
!
! these steps control those slow physics (PDC using tendency method)
! Fast physics step will be controlled by model_timestep (PDC using operator splitting)
!
   integer(i4)        :: step_cu  = 1  ! how many model_timestep to call cumulus, default every step
   integer(i4)        :: step_ra   ! how many model_timestep to call radiation, default: 3h
   integer(i4)        :: step_bl  = 1  ! not used as bl is currently a fast physics
   integer(i4)        :: wphys_has_req  = 0
   logical            :: use_gwdo
   logical            :: use_cond
   logical            :: unuse_cu = .false.
   logical            :: force_read_thompson = .false.
   logical            :: write_thompson_tables = .false.
! sensitivity for entrainment
   real(r8)           :: para1_zoentr_test = 1._r8
   real(r8)           :: para2_zoentr_test = 1._r8
   real(r8)           :: para1_eta_test    = 1._r8
   real(r8)           :: para2_eta_test    = 1._r8

  contains

  subroutine set_PhysW_nml()

!================================================
! global vars have been defined in the header
!================================================

! local
  character(len=300) :: filename
  integer (i4)       :: fileunit

   namelist /PhysW_para/  PhysW_cu_scheme,  & ! CUmulus
                          PhysW_mp_scheme,  & ! MicroPhysics
                          PhysW_bl_scheme,  & ! Boundary Layer
                          PhysW_cf_scheme,  & ! RAdiation
                          PhysW_ra_scheme,  & ! RAdiation
                          PhysW_rasw_scheme,& ! RAdiation
                          PhysW_ralw_scheme,& ! RAdiation
                          PhysW_sf_scheme,  & ! Surface Flux
                          PhysW_lm_scheme,  & ! Land Model
                          step_cu,            &
                          step_ra,            &
                          step_bl,            &
                          wphys_has_req,      &
                          use_gwdo,           &
                          use_cond,           &
                          unuse_cu,           &
                          force_read_thompson,&
                          write_thompson_tables, &
                          para1_zoentr_test, &
                          para2_zoentr_test, &
                          para1_eta_test,    &
                          para2_eta_test

    filename = "grist_amipw_phys.nml"

    fileunit = 1

    open  (fileunit, status='old',file=filename)
    read  (fileunit, nml=PhysW_para)
    close (fileunit)

  !  step_ra = 10800/idint(model_timestep)

    if(mpi_rank() .eq. 0)then
       print*,"**********************************************************"
       print*,"     The WRFPhys is used following: ", trim(filename)
       print*,"     step_ra, step_cu, step_bl are:", step_ra, step_cu, step_bl
       print*,"**********************************************************"
    end if

    return
  end subroutine set_PhysW_nml

  end module grist_PhysW_nml_module
