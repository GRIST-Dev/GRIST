!==========================================================================
!
!  Created by LiXiaohan on 19/8/19, adopted from CAM5.
!
!  This module uses the Lean solar irradiance data to provide a solar cycle
!  scaling factor used in heating rate calculations 
!===========================================================================
 module rad_solar_var

    use grist_constants,          only: r8
    use grist_handle_error,       only: endrun
    use solar_data,               only: sol_irrad, we, nbins, has_spectrum, sol_tsi, &
                                        do_spctrl_scaling
    use grist_mpi


    implicit none
    save

    private
    public :: rad_solar_var_init,   &
              get_variability,      &
              end_of_rad_solar_var

    real(r8), allocatable :: ref_band_irrad(:)  ! scaling will be relative to ref_band_irrad in each band
    real(r8), allocatable :: irrad(:)           ! solar irradiance at model timestep in each band
    real(r8)              :: tsi_ref            ! total solar irradiance assumed by rrtmg                                                 

    real(r8), allocatable :: radbinmax(:)
    real(r8), allocatable :: radbinmin(:)
    integer :: nradbins
contains

    subroutine rad_solar_var_init( )
    use radconstants,  only : get_number_sw_bands,          &
                              get_sw_spectral_boundaries,   &
                              get_ref_solar_band_irrad,     &
                              get_ref_total_solar_irrad

    integer :: i
    integer :: ierr
    integer :: yr, mon, tod
    integer :: radmax_loc

    call get_number_sw_bands(nradbins)

    if ( do_spctrl_scaling ) then

       if ( .not.has_spectrum ) then
          call endrun('rad_solar_var_init: solar input file must have irradiance spectrum')
       endif

       allocate (radbinmax(nradbins))
       allocate (radbinmin(nradbins))
       allocate (ref_band_irrad(nradbins))
       allocate (irrad(nradbins))

       call get_sw_spectral_boundaries(radbinmin, radbinmax, 'nm')

       ! Make sure that the far-IR is included, even if RRTMG does not
       ! extend that far down. 10^5 nm corresponds to a wavenumber of
       ! 100 cm^-1.
       radmax_loc = maxloc(radbinmax,1)
       radbinmax(radmax_loc) = max(100000._r8,radbinmax(radmax_loc))

       ! for rrtmg, reference spectrum from rrtmg
       call get_ref_solar_band_irrad( ref_band_irrad )

    else

       call get_ref_total_solar_irrad(tsi_ref)

    endif

    endsubroutine rad_solar_var_init


    subroutine end_of_rad_solar_var
    if(allocated (radbinmax))       deallocate (radbinmax)
    if(allocated (radbinmin))       deallocate (radbinmin)
    if(allocated (ref_band_irrad))  deallocate (ref_band_irrad)
    if(allocated (irrad))           deallocate (irrad) 

    end subroutine end_of_rad_solar_var


    subroutine get_variability( sfac )

    real(r8), intent(out) :: sfac(nradbins)       ! scaling factors for CAM heating

    integer :: yr, mon, day, tod

    if ( do_spctrl_scaling ) then
      call integrate_spectrum( nbins, nradbins, we, radbinmin, radbinmax, sol_irrad, irrad)

      sfac(:nradbins) = irrad(:nradbins)/ref_band_irrad(:nradbins)

    else

       sfac(:nradbins) = sol_tsi/tsi_ref

    endif

    endsubroutine get_variability


    subroutine integrate_spectrum( nsrc, ntrg, src_x, min_trg, max_trg, src, trg )

    use mo_util, only : rebin

    implicit none

    !---------------------------------------------------------------
    !	... dummy arguments
    !---------------------------------------------------------------
    integer,  intent(in)  :: nsrc                  ! dimension source array
    integer,  intent(in)  :: ntrg                  ! dimension target array
    real(r8), intent(in)  :: src_x(nsrc+1)         ! source coordinates
    real(r8), intent(in)  :: max_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: min_trg(ntrg)         ! target coordinates
    real(r8), intent(in)  :: src(nsrc)             ! source array
    real(r8), intent(out) :: trg(ntrg)             ! target array
 
    !---------------------------------------------------------------
    !	... local variables
    !---------------------------------------------------------------
    real(r8) :: trg_x(2), targ(1)         ! target coordinates
    integer  :: i

    do i = 1, ntrg

       trg_x(1) = min_trg(i)
       trg_x(2) = max_trg(i)

       call rebin( nsrc, 1, src_x, trg_x, src(1:nsrc), targ(:) )
       ! W/m2/nm --> W/m2
       trg( i ) = targ(1)*(trg_x(2)-trg_x(1))

    enddo


  end subroutine integrate_spectrum

 end module rad_solar_var
