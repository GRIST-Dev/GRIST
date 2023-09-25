!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Data Module for GRIST
!              Name Convention:
!              initial_data_${varname}_at_$(location)
! Revision history:
!              This entry contains data from initial file, which contains
! both atmosphere (U,V,T,Q,PS) and land-surface (SoilTemp,SoilMoist,SkinTemp,
! XICE,SNOW,SNOWH,SOILH) as standard initial input. If one wants to do "warm"
! start, please use the "restart" capability by modifying a restart file.
!
! initialDataSorc: (1) ERAIP/GFS, any pressure-level data from (re)analysis
!                  (2) ERAIM, any hybrid-sigma-moistPres model-level data
!                  (3) WRFDA, from WRFDA-like data on a dry-mass coordinate
! "Keep this part simple and stupid"
!----------------------------------------------------------------------------

  module grist_datam_initial_data_module

    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_1d_field, scalar_2d_field
    use grist_constants,    only: rearth, i4, i8, r8, pi, zero
    use grist_nml_module,   only: nlev_inidata, levsoil, initialAtmFilePath, initialDataSorc, initialLndFilePath, &
                                  initialAtmUFilePath,initialAtmVFilePath,initialAtmTFilePath,initialAtmQFilePath, large_atm_file_on
    use grist_fileio_list_1d_module_par, only: wrap_read_1d_group_rst
    use grist_fileio_list_2d_module_par, only: wrap_read_2d_group_rst
    use grist_mpi

   implicit none

   public

!================================================
! primal cell, full level
!================================================
! ATM
   type(scalar_2d_field) :: initialData_uuu_at_pc_full_level
   type(scalar_2d_field) :: initialData_vvv_at_pc_full_level
   type(scalar_2d_field) :: initialData_ttt_at_pc_full_level
   type(scalar_2d_field) :: initialData_qqq_at_pc_full_level  ! in most reanalysis, this is specific humidity
   !type(scalar_2d_field) :: initialData_ppp_at_pc_full_level ! deprecated
! yizhang added for WRFDA initial from sunwei
   type(scalar_2d_field) :: initialData_hypw_at_pc_face_level ! already on dry-mass level

   type(scalar_1d_field) :: initialData_ps_at_pc_surface
   type(scalar_1d_field) :: initialData_hyai
   type(scalar_1d_field) :: initialData_hybi
   type(scalar_1d_field) :: initialData_hyam
   type(scalar_1d_field) :: initialData_hybm
   type(scalar_1d_field) :: initialData_plev
! LND
!---cheyz add
   type(scalar_1d_field) :: initialData_skintemp_at_pc_surface
   type(scalar_2d_field) :: initialData_soiltemp_at_pc_soil_level
   type(scalar_2d_field) :: initialData_soilmoist_at_pc_soil_level
   type(scalar_1d_field) :: initialData_soilh_at_pc_surface
   type(scalar_1d_field) :: initialData_xice_at_pc_surface
   type(scalar_1d_field) :: initialData_snow_at_pc_surface
   type(scalar_1d_field) :: initialData_snowh_at_pc_surface

  CONTAINS

    subroutine grist_initial_data_construct(mesh)
! io
      type(global_domain), intent(in), target :: mesh
! local
      integer(i4)  :: it, ie, iv
      integer(i8)  :: dim2_len
      integer(i8)  :: err
      real(r8), allocatable :: buf1(:)
      real(r8), allocatable :: buf2(:)
!----------------
! ATM initial
!----------------
!
! mesh%nv is actually nv_full by call
!
    if(.not.allocated(initialData_uuu_at_pc_full_level%f))allocate(initialData_uuu_at_pc_full_level%f(nlev_inidata,mesh%nv));initialData_uuu_at_pc_full_level%pos = 0
    if(.not.allocated(initialData_vvv_at_pc_full_level%f))allocate(initialData_vvv_at_pc_full_level%f(nlev_inidata,mesh%nv));initialData_vvv_at_pc_full_level%pos = 0
    if(.not.allocated(initialData_ttt_at_pc_full_level%f))allocate(initialData_ttt_at_pc_full_level%f(nlev_inidata,mesh%nv));initialData_ttt_at_pc_full_level%pos = 0
    if(.not.allocated(initialData_qqq_at_pc_full_level%f))allocate(initialData_qqq_at_pc_full_level%f(nlev_inidata,mesh%nv));initialData_qqq_at_pc_full_level%pos = 0
    if(.not.allocated(initialData_ps_at_pc_surface%f))    allocate(initialData_ps_at_pc_surface%f(mesh%nv))                 ;initialData_ps_at_pc_surface%pos     = 0
    !if(.not.allocated(initialData_ppp_at_pc_full_level%f))allocate(initialData_ppp_at_pc_full_level%f(nlev_inidata,mesh%nv));initialData_ppp_at_pc_full_level%pos = 0 !cheyz 5/26!--cheyz
    if(trim(initialDataSorc).eq.'WRFDA')then
       if(.not.allocated(initialData_hypw_at_pc_face_level%f)) allocate(initialData_hypw_at_pc_face_level%f(nlev_inidata+1,mesh%nv));initialData_hypw_at_pc_face_level%pos = 0 ! wrfda
    end if
    if(trim(initialDataSorc).eq.'ERAIM')then
       if(.not.allocated(initialData_hyai%f))  allocate(initialData_hyai%f(nlev_inidata+1))
       if(.not.allocated(initialData_hybi%f))  allocate(initialData_hybi%f(nlev_inidata+1))
       if(.not.allocated(initialData_hyam%f))  allocate(initialData_hyam%f(nlev_inidata))
       if(.not.allocated(initialData_hybm%f))  allocate(initialData_hybm%f(nlev_inidata))
    end if
    if(trim(initialDataSorc).eq.'ERAIP'.or.(trim(initialDataSorc).eq.'GFS'))then
       if(.not.allocated(initialData_plev%f))  allocate(initialData_plev%f(nlev_inidata))
    end if

!----------------
! LND initial
!----------------

#ifdef USE_NOAHMP
    if (.not.allocated(initialData_skintemp_at_pc_surface%f))    allocate(initialData_skintemp_at_pc_surface%f(mesh%nv), source=zero) ;initialData_skintemp_at_pc_surface%pos=0
    if (.not.allocated(initialData_snow_at_pc_surface%f))        allocate(initialData_snow_at_pc_surface%f(mesh%nv), source=zero)     ;initialData_snow_at_pc_surface%pos    =0
    if (.not.allocated(initialData_soiltemp_at_pc_soil_level%f)) allocate(initialData_soiltemp_at_pc_soil_level%f(levsoil,mesh%nv), source=zero); initialData_soiltemp_at_pc_soil_level%pos=0
    if (.not.allocated(initialData_soilmoist_at_pc_soil_level%f))allocate(initialData_soilmoist_at_pc_soil_level%f(levsoil,mesh%nv), source=zero); initialData_soilmoist_at_pc_soil_level%pos=0
    if (.not.allocated(initialData_xice_at_pc_surface%f))        allocate(initialData_xice_at_pc_surface%f(mesh%nv),  source=zero) ;initialData_xice_at_pc_surface%pos   = 0
    if (.not.allocated(initialData_snowh_at_pc_surface%f))       allocate(initialData_snowh_at_pc_surface%f(mesh%nv), source=zero) ;initialData_snowh_at_pc_surface%pos  = 0
    if (.not.allocated(initialData_soilh_at_pc_surface%f))       allocate(initialData_soilh_at_pc_surface%f(mesh%nv), source=zero) ;initialData_soilh_at_pc_surface%pos  = 0
#endif

! All five info will be needed
     dim2_len = int(nlev_inidata,i8)
! for very large atm_file, we may read in from each file
     if(.not.large_atm_file_on)then
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'U', dim2_len, 0, initialData_uuu_at_pc_full_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'V', dim2_len, 0, initialData_vvv_at_pc_full_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'T', dim2_len, 0, initialData_ttt_at_pc_full_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'Q', dim2_len, 0, initialData_qqq_at_pc_full_level%f)
     else
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmUFilePath), '', 'U', dim2_len, 0, initialData_uuu_at_pc_full_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmVFilePath), '', 'V', dim2_len, 0, initialData_vvv_at_pc_full_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmTFilePath), '', 'T', dim2_len, 0, initialData_ttt_at_pc_full_level%f)
        call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmQFilePath), '', 'Q', dim2_len, 0, initialData_qqq_at_pc_full_level%f)
     end if

     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'PS',          0, initialData_ps_at_pc_surface%f)
    
     if(trim(initialDataSorc).eq.'WRFDA')then
       dim2_len = int(nlev_inidata+1,i8)
       call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'P_HYD_W', dim2_len, 0, initialData_hypw_at_pc_face_level%f)
     end if

     if(trim(initialDataSorc).eq.'ERAIM')then
        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'hyai',       9, initialData_hyai%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'hybi',       9, initialData_hybi%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'hyam',       8, initialData_hyam%f)
        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'hybm',       8, initialData_hybm%f)
     end if
     if(trim(initialDataSorc).eq.'ERAIP'.or.(trim(initialDataSorc).eq.'GFS'))then
        call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialAtmFilePath), '', 'plev',       8, initialData_plev%f)
     end if

#ifdef USE_NOAHMP
     dim2_len = int(levsoil,i8)
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialLndFilePath), '', 'SoilTemp' , dim2_len, 0, initialData_soiltemp_at_pc_soil_level%f)
     call wrap_read_2d_group_rst(mesh%gcomm_read, trim(initialLndFilePath), '', 'SoilMoist', dim2_len, 0, initialData_soilmoist_at_pc_soil_level%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialLndFilePath), '', 'SKINTEMP'           , 0, initialData_skintemp_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialLndFilePath), '', 'XICE'               , 0, initialData_xice_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialLndFilePath), '', 'SNOW'               , 0, initialData_snow_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialLndFilePath), '', 'SNOWH'              , 0, initialData_snowh_at_pc_surface%f)
     call wrap_read_1d_group_rst(mesh%gcomm_read, trim(initialLndFilePath), '', 'SOILH'              , 0, initialData_soilh_at_pc_surface%f)
#endif

     if(mpi_rank().eq.0) print*,"initial ATM and LND data constructed sucessfully in grist_initial_data_construct"

   return
   end subroutine grist_initial_data_construct

   subroutine grist_initial_data_destruct
!----------------
! ATM initial
!----------------
      if(allocated(initialData_uuu_at_pc_full_level%f)) deallocate(initialData_uuu_at_pc_full_level%f)
      if(allocated(initialData_vvv_at_pc_full_level%f)) deallocate(initialData_vvv_at_pc_full_level%f)
      if(allocated(initialData_ttt_at_pc_full_level%f)) deallocate(initialData_ttt_at_pc_full_level%f)
      if(allocated(initialData_qqq_at_pc_full_level%f)) deallocate(initialData_qqq_at_pc_full_level%f)
      if(allocated(initialData_ps_at_pc_surface%f))     deallocate(initialData_ps_at_pc_surface%f)
      !if(allocated(initialData_ppp_at_pc_full_level%f)) deallocate(initialData_ppp_at_pc_full_level%f)
      if(allocated(initialData_hypw_at_pc_face_level%f))deallocate(initialData_hypw_at_pc_face_level%f)

      if(allocated(initialData_hyai%f))                 deallocate(initialData_hyai%f)
      if(allocated(initialData_hybi%f))                 deallocate(initialData_hybi%f)
      if(allocated(initialData_hyam%f))                 deallocate(initialData_hyam%f)
      if(allocated(initialData_hybm%f))                 deallocate(initialData_hybm%f)
      if(allocated(initialData_plev%f))                 deallocate(initialData_plev%f)
!----------------
! LND initial
!----------------
      if(allocated(initialData_skintemp_at_pc_surface%f))     deallocate(initialData_skintemp_at_pc_surface%f)
      if(allocated(initialData_snow_at_pc_surface%f))         deallocate(initialData_snow_at_pc_surface%f)
      if(allocated(initialData_soiltemp_at_pc_soil_level%f))  deallocate(initialData_soiltemp_at_pc_soil_level%f)
      if(allocated(initialData_soilmoist_at_pc_soil_level%f)) deallocate(initialData_soilmoist_at_pc_soil_level%f)
      if(allocated(initialData_xice_at_pc_surface%f))         deallocate(initialData_xice_at_pc_surface%f)
      if(allocated(initialData_snowh_at_pc_surface%f))        deallocate(initialData_snowh_at_pc_surface%f)
      if(allocated(initialData_soilh_at_pc_surface%f))        deallocate(initialData_soilh_at_pc_surface%f)

      return
   end subroutine grist_initial_data_destruct
  end module grist_datam_initial_data_module
