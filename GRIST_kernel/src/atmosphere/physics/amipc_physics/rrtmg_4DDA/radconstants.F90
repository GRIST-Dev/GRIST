!===================================================================================================
!
!  Created by LiXiaohan on 19/08/10, adopted from CAM5
!
!  This module contains constants that are specific to the radiative transfer
!  code used in the RRTMG model.
!  
!===================================================================================================

 module radconstants

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8, i4

    implicit none
    private
    save
    
    ! SHORTWAVE DATA
    
    ! number of shorwave spectral intervals
    integer, parameter, public :: nswbands = 14
    integer, parameter, public :: nbndsw = 14
    
    ! Wavenumbers of band boundaries
    !
    ! Note: Currently rad_solar_var extends the lowest band down to
    ! 100 cm^-1 if it is too high to cover the far-IR. Any changes meant
    ! to affect IR solar variability should take note of this.
    
    real(r8),parameter :: wavenum_low(nbndsw) = & ! in cm^-1
      (/2600._r8, 3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, &
        8050._r8,12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,  820._r8/)
    real(r8),parameter :: wavenum_high(nbndsw) = & ! in cm^-1
      (/3250._r8, 4000._r8, 4650._r8, 5150._r8, 6150._r8, 7700._r8, 8050._r8, &
       12850._r8,16000._r8,22650._r8,29000._r8,38000._r8,50000._r8, 2600._r8/)
    
    ! Solar irradiance at 1 A.U. in W/m^2 assumed by radiation code
    ! Rescaled so that sum is precisely 1368.22 and fractional amounts sum to 1.0
    real(r8), parameter :: solar_ref_band_irradiance(nbndsw) = & 
       (/ &
        12.11_r8,  20.3600000000001_r8, 23.73_r8, &
        22.43_r8,  55.63_r8, 102.93_r8, 24.29_r8, &
       345.74_r8, 218.19_r8, 347.20_r8, &
       129.49_r8,  50.15_r8,   3.08_r8, 12.89_r8 &
       /)
    
    ! None of the following comment appears to be the case any more? This
    ! should be reevalutated and/or removed.
    
    ! rrtmg (coarse) reference solar flux in rrtmg is initialized as the following
    ! reference data inside rrtmg seems to indicate 1366.44 instead
    !  This data references 1366.442114152342
    !real(r8), parameter :: solar_ref_band_irradiance(nbndsw) = & 
    !   (/ &
    !   12.10956827000000_r8, 20.36508467999999_r8, 23.72973826333333_r8, &
    !   22.42769644333333_r8, 55.62661262000000_r8, 102.9314315544444_r8, 24.29361887666667_r8, &
    !   345.7425138000000_r8, 218.1870300666667_r8, 347.1923147000001_r8, &
    !   129.4950181200000_r8, 48.37217043000000_r8, 3.079938997898001_r8, 12.88937733000000_r8 &
    !   /)
    !  Kurucz (fine) reference would seem to imply the following but the above values are from rrtmg_sw_init
    !  (/12.109559, 20.365097, 23.729752, 22.427697, 55.626622, 102.93142, 24.293593, &
    !    345.73655, 218.18416, 347.18406, 129.49407, 50.147238, 3.1197130, 12.793834 /)
    
    ! These are indices to the band for diagnostic output
    integer, parameter, public :: idx_sw_diag = 10 ! index to sw visible band
    integer, parameter, public :: idx_nir_diag = 8 ! index to sw near infrared (778-1240 nm) band
    integer, parameter, public :: idx_uv_diag = 11 ! index to sw uv (345-441 nm) band
    
    integer, parameter, public :: rrtmg_sw_cloudsim_band = 9  ! rrtmg band for .67 micron
    
    ! Number of evenly spaced intervals in rh
    ! The globality of this mesh may not be necessary
    ! Perhaps it could be specific to the aerosol
    ! But it is difficult to see how refined it must be
    ! for lookup.  This value was found to be sufficient
    ! for Sulfate and probably necessary to resolve the
    ! high variation near rh = 1.  Alternative methods
    ! were found to be too slow.
    ! Optimal approach would be for cam to specify size of aerosol
    ! based on each aerosol's characteristics.  Radiation 
    ! should know nothing about hygroscopic growth!
    integer, parameter, public :: nrh = 1000  
    
    ! LONGWAVE DATA
    
    ! These are indices to the band for diagnostic output
    integer, parameter, public :: idx_lw_diag = 7 ! index to (H20 window) LW band
    
    integer, parameter, public :: rrtmg_lw_cloudsim_band = 6  ! rrtmg band for 10.5 micron
    
    ! number of lw bands
    integer, parameter, public :: nlwbands = 16
    integer, parameter, public :: nbndlw = 16
    
    real(r8), parameter :: wavenumber1_longwave(nlwbands) = &! Longwave spectral band limits (cm-1)
        (/   10._r8,  350._r8, 500._r8,   630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, &
           1180._r8, 1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2600._r8 /)
    
    real(r8), parameter :: wavenumber2_longwave(nlwbands) = &! Longwave spectral band limits (cm-1)
        (/  350._r8,  500._r8,  630._r8,  700._r8,  820._r8,  980._r8, 1080._r8, 1180._r8, &
           1390._r8, 1480._r8, 1800._r8, 2080._r8, 2250._r8, 2390._r8, 2600._r8, 3250._r8 /)
    
    !These can go away when old camrt disappears
    ! Index of volc. abs., H2O non-window
    integer, public, parameter :: idx_LW_H2O_NONWND=1
    ! Index of volc. abs., H2O window
    integer, public, parameter :: idx_LW_H2O_WINDOW=2
    ! Index of volc. cnt. abs. 0500--0650 cm-1
    integer, public, parameter :: idx_LW_0500_0650=3
    ! Index of volc. cnt. abs. 0650--0800 cm-1
    integer, public, parameter :: idx_LW_0650_0800=4
    ! Index of volc. cnt. abs. 0800--1000 cm-1
    integer, public, parameter :: idx_LW_0800_1000=5
    ! Index of volc. cnt. abs. 1000--1200 cm-1
    integer, public, parameter :: idx_LW_1000_1200=6
    ! Index of volc. cnt. abs. 1200--2000 cm-1
    integer, public, parameter :: idx_LW_1200_2000=7
    
    ! GASES TREATED BY RADIATION (line spectrae)
    
    ! gasses required by radiation
    integer, public, parameter :: gasnamelength = 5
    integer, public, parameter :: nradgas = 8
    character(len=gasnamelength), public, parameter :: gaslist(nradgas) &
       = (/'H2O  ','O3   ', 'O2   ', 'CO2  ', 'N2O  ', 'CH4  ', 'CFC11', 'CFC12'/)
    
    ! what is the minimum mass mixing ratio that can be supported by radiation implementation?
    real(r8), public, parameter :: minmmr(nradgas) &
       = epsilon(1._r8)
    
    ! Length of "optics type" string specified in optics files.
    integer, parameter, public :: ot_length = 32

    public :: rad_gas_index,                    &
              get_number_sw_bands,              &
              get_sw_spectral_boundaries,       &
              get_lw_spectral_boundaries,       &
              get_ref_solar_band_irrad,         &
              get_ref_total_solar_irrad,        &
              get_solar_band_fraction_irrad

 contains

!Purpose : provide Solar Irradiance for each band in RRTMG
    subroutine get_solar_band_fraction_irrad(fractional_irradiance)

    ! fraction of solar irradiance in each band
    real(r8), intent(out) :: fractional_irradiance(1:nswbands)
    real(r8) :: tsi ! total solar irradiance

    tsi = sum(solar_ref_band_irradiance)
    fractional_irradiance = solar_ref_band_irradiance / tsi

    end subroutine get_solar_band_fraction_irrad


!Purpose : provide Total Solar Irradiance assumed by RRTMG
    subroutine get_ref_total_solar_irrad(tsi)

    real(r8), intent(out) :: tsi

    tsi = sum(solar_ref_band_irradiance)

    end subroutine get_ref_total_solar_irrad


!Purpose : get solar irradiance in each band (W/m^2)
    subroutine get_ref_solar_band_irrad( band_irrad )

    real(r8), intent(out) :: band_irrad(nswbands)
 
    band_irrad = solar_ref_band_irradiance

    end subroutine get_ref_solar_band_irrad


!Purpose : get number of solar (shortwave) bands in the rrtmg code
    subroutine get_number_sw_bands(number_of_bands)

    integer, intent(out) :: number_of_bands

    number_of_bands = nswbands

    end subroutine get_number_sw_bands


!Purpose : provide spectral boundaries of each longwave band
    subroutine get_lw_spectral_boundaries(low_boundaries, high_boundaries, units)

    real(r8), intent(out) :: low_boundaries(nlwbands), high_boundaries(nlwbands)
    character(*), intent(in) :: units ! requested units

    select case (units)
    case ('inv_cm','cm^-1','cm-1')
       low_boundaries  = wavenumber1_longwave
       high_boundaries = wavenumber2_longwave
    case('m','meter','meters')
       low_boundaries  = 1.e-2_r8/wavenumber2_longwave
       high_boundaries = 1.e-2_r8/wavenumber1_longwave
    case('nm','nanometer','nanometers')
       low_boundaries  = 1.e7_r8/wavenumber2_longwave
       high_boundaries = 1.e7_r8/wavenumber1_longwave
    case('um','micrometer','micrometers','micron','microns')
       low_boundaries  = 1.e4_r8/wavenumber2_longwave
       high_boundaries = 1.e4_r8/wavenumber1_longwave
    case('cm','centimeter','centimeters')
       low_boundaries  = 1._r8/wavenumber2_longwave
       high_boundaries = 1._r8/wavenumber1_longwave
    case default
       call endrun('get_lw_spectral_boundaries: spectral units not acceptable'//units)
    end select

    end subroutine get_lw_spectral_boundaries


!Purpose : provide spectral boundaries of each shortwave band
    subroutine get_sw_spectral_boundaries(low_boundaries, high_boundaries, units)

    real(r8), intent(out) :: low_boundaries(nswbands), high_boundaries(nswbands)
    character(*), intent(in) :: units ! requested units

    select case (units)
    case ('inv_cm','cm^-1','cm-1')
       low_boundaries = wavenum_low
       high_boundaries = wavenum_high
    case('m','meter','meters')
       low_boundaries = 1.e-2_r8/wavenum_high
       high_boundaries = 1.e-2_r8/wavenum_low
    case('nm','nanometer','nanometers')
       low_boundaries = 1.e7_r8/wavenum_high
       high_boundaries = 1.e7_r8/wavenum_low
    case('um','micrometer','micrometers','micron','microns')
       low_boundaries = 1.e4_r8/wavenum_high
       high_boundaries = 1.e4_r8/wavenum_low
    case('cm','centimeter','centimeters')
       low_boundaries  = 1._r8/wavenum_high
       high_boundaries = 1._r8/wavenum_low
    case default
       call endrun('rad_constants.F90: spectral units not acceptable'//units)
    end select

    end subroutine get_sw_spectral_boundaries

      
!Purpose : return the index in the gaslist array of the specified gasname
    integer function rad_gas_index(gasname)
    
    character(len=*),intent(in) :: gasname
    integer :: igas
    
    rad_gas_index = -1
    do igas = 1, nradgas
       if (trim(gaslist(igas)).eq.trim(gasname)) then
          rad_gas_index = igas
          return
       endif
    enddo
    call endrun ("rad_gas_index: can not find gas with name "//gasname)
    end function rad_gas_index

 end module radconstants
