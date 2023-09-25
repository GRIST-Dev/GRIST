
  module noahmp_nml_control
  use grist_constants,      only: i4, r8
!    use noahmp_options_config

    implicit none
    
    public :: read_nml_noahmp     ! read namelist from file
  
    integer(i4)             ::  idveg      
    integer(i4)             ::  iopt_crs   ! canopy stomatal resistance (1-> ball-berry; 2->jarvis)
    integer(i4)             ::  iopt_btr   ! soil moisture factor for stomatal resistance (1-> noah; 2-> clm; 3-> ssib)
    integer(i4)             ::  iopt_run   ! runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
    integer(i4)             ::  iopt_sfc   ! surface layer drag coeff (ch & cm) (1->m-o; 2->chen97)
    integer(i4)             ::  iopt_frz   ! supercooled liquid water (1-> ny06; 2->koren99)
    integer(i4)             ::  iopt_inf   ! frozen soil permeability (1-> ny06; 2->koren99)
    integer(i4)             ::  iopt_rad   ! radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
    integer(i4)             ::  iopt_alb   ! snow surface albedo (1->bats; 2->class)
    integer(i4)             ::  iopt_snf   ! rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
    integer(i4)             ::  iopt_tbot  ! lower boundary of soil temperature (1->zero-flux; 2->noah)
    integer(i4)             ::  iopt_stc   ! snow/soil temperature time scheme
    integer(i4)             ::  iopt_gla   ! glacier option (1->phase change; 2->simple)
    integer(i4)             ::  iopt_rsf   ! surface resistance option (1->zeng; 2->simple)
   !---cheyz
    integer(i4), public     ::  nsoil         ! number of soil layers
    real(r8),allocatable    ::  dzs1(:)       ! thickness of soil layers [m]
    integer                 ::  iz0tlnd 

   ! Namelist variables:             
    integer(i4)            :: dynamic_veg_option
    integer(i4)            :: canopy_stomatal_resistance_option
    integer(i4)            :: btr_option
    integer(i4)            :: runoff_option
    integer(i4)            :: surface_drag_option
    integer(i4)            :: supercooled_water_option
    integer(i4)            :: frozen_soil_option
    integer(i4)            :: radiative_transfer_option
    integer(i4)            :: snow_albedo_option
    integer(i4)            :: pcp_partition_option
    integer(i4)            :: tbot_option
    integer(i4)            :: temp_time_scheme_option
    integer(i4)            :: glacier_option
    integer(i4)            :: surface_resistance_option 
    integer(i4)	         :: num_soil_layers  
	 integer(i4), parameter :: max_soil_levels = 4   ! maximum soil levels in namelist
    real(r8)               :: soil_thick_input(max_soil_levels)       ! depth to soil interfaces from namelist [m]
 
    contains

   subroutine read_nml_noahmp(filename)
      implicit none
	integer :: ierr
	character(len=300),intent(in) :: filename ! filepath for file containing namelist input

    namelist /lsm_noahmp_nl/ dynamic_veg_option, canopy_stomatal_resistance_option, &
       btr_option, runoff_option, surface_drag_option, supercooled_water_option, &
       frozen_soil_option, radiative_transfer_option, snow_albedo_option, &
       pcp_partition_option, tbot_option, temp_time_scheme_option, &
       glacier_option, surface_resistance_option
	namelist /lsm_soil_layers_nl/ num_soil_layers, soil_thick_input
	
! Initialize namelist variables to dummy values, so we can tell
! if they have not been set properly.

  nsoil                   = 0
  soil_thick_input        = 0._r8
  
  open(30, file=trim(filename), status='old', iostat=ierr)
   if (ierr /= 0) then
     write(*,'(/," ***** ERROR: Problem reading namelist for noahmp lsm",/)')
     stop " ***** ERROR: Problem reading namelist for noahmp lsm"
   endif
    read(30, lsm_noahmp_nl)
	 read(30, lsm_soil_layers_nl)
   close(30)

   nsoil = num_soil_layers     
   idveg = dynamic_veg_option 
   iopt_crs = canopy_stomatal_resistance_option
   iopt_btr = btr_option
   iopt_run = runoff_option
   iopt_sfc = surface_drag_option
   iopt_frz = supercooled_water_option
   iopt_inf = frozen_soil_option
   iopt_rad = radiative_transfer_option
   iopt_alb = snow_albedo_option
   iopt_snf = pcp_partition_option
   iopt_tbot = tbot_option
   iopt_stc = temp_time_scheme_option
   iopt_gla = glacier_option
   iopt_rsf = surface_resistance_option
   dzs1     =  soil_thick_input(1:nsoil)

!---------------------------------------------------------------------
!  NAMELIST end
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!  NAMELIST check begin
!---------------------------------------------------------------------

  if (nsoil < 0) then
     stop " ***** ERROR: NSOIL must be set in the namelist."
  endif

  if (dynamic_veg_option == 2 .or. dynamic_veg_option == 5 .or. dynamic_veg_option == 6) then
     if ( canopy_stomatal_resistance_option /= 1) then
        write(*, *)
        write(*, '(" ***** Namelist error: ******************************************************")')
        write(*, '(" ***** ")')
        write(*, '(" *****       CANOPY_STOMATAL_RESISTANCE_OPTION must be 1 when DYNAMIC_VEG_OPTION == 2/5/6")')
        write(*, *)
        stop
     endif
  endif

!---------------------------------------------------------------------
!  NAMELIST check end
!---------------------------------------------------------------------

  end subroutine read_nml_noahmp
 end module noahmp_nml_control
