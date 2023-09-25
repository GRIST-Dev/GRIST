!========================================================================================
!
!  Created by LiXiaohan on 19/8/13, adopted from CAM5.
!  Provide default distributions of CH4, N2O, CFC11 and CFC12 to the radiation routines.
!  **NOTE** CO2 is assumed by the radiation to a be constant value.  This value is
!           currently supplied directly by the chem_surfvals module.
!
!========================================================================================

 module ghg_data

    use grist_constants,              only: r8,   &
                                            mwdry, mwch4, mwn2o, mwf11, mwf12, mwco2
    use grist_nml_module,             only: nlev, nlevp
    use grist_handle_error,           only: endrun
    use grist_physics_data_structure, only: pstate
    use grist_PhysC_data_structure,   only: pstate_cam
    use grist_mpi
 
 implicit none
 private
 save

 public :: ghg_data_init,    &
           end_of_ghg_data,  &
           ghg_data_timestep_init

 real(r8) :: rmwn2o ! = mwn2o/mwdry ! ratio of molecular weight n2o   to dry air
 real(r8) :: rmwch4 ! = mwch4/mwdry ! ratio of molecular weight ch4   to dry air
 real(r8) :: rmwf11 ! = mwf11/mwdry ! ratio of molecular weight cfc11 to dry air
 real(r8) :: rmwf12 ! = mwf12/mwdry ! ratio of molecular weight cfc12 to dry air
 real(r8) :: rmwco2 ! = mwco2/mwdry ! ratio of molecular weights of co2 to dry air
 
! integer, parameter :: ncnst = 6                        ! number of constituents
! character(len=8), dimension(ncnst), parameter :: &
!   cnst_names = (/'N2O  ', 'CH4  ', 'CFC11', 'CFC12', 'CO2  ', 'O2   '/) ! constituent names
!------------LiXH add ozone-----------
 integer, parameter :: ncnst = 7                        ! number of constituents
 character(len=8), dimension(ncnst), parameter :: &
   cnst_names = (/'N2O  ', 'CH4  ', 'CFC11', 'CFC12', 'CO2  ', 'O2   ','O3   '/) ! constituent names
!------------LiXH add ozone-----------


 contains

    subroutine ghg_data_init(ncol)
    ! io
    integer, intent(in) :: ncol
    ! local
    integer :: i

    pstate_cam%total_ghg_num = ncnst
    allocate(pstate_cam%ghg_at_pc_full_level(ncnst))
    do i = 1, ncnst
       pstate_cam%ghg_at_pc_full_level(i)%idx = i
       pstate_cam%ghg_at_pc_full_level(i)%name = trim(cnst_names(i))
       allocate(pstate_cam%ghg_at_pc_full_level(i)%f(nlev,ncol))
       pstate_cam%ghg_at_pc_full_level(i)%f = 0._r8
    end do

    end subroutine ghg_data_init


    subroutine end_of_ghg_data

    integer :: i

    do i = 1, ncnst
        deallocate(pstate_cam%ghg_at_pc_full_level(i)%f)
    end do

    deallocate(pstate_cam%ghg_at_pc_full_level)

    end subroutine end_of_ghg_data


    subroutine ghg_data_timestep_init(ncol, lat)
    ! io
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: lat(ncol)
    ! local
    integer :: iconst

    rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
    rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
    rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
    rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
    rmwco2 = mwco2/mwdry      ! ratio of molecular weights of co2 to dry air

    do iconst = 1,ncnst
        call trcmix(trim(cnst_names(iconst)), ncol, lat,                               &
                    pstate%pressure_at_pc_full_level%f(:,1:ncol),   &
                    pstate_cam%ghg_at_pc_full_level(iconst)%f(:,1:ncol))
    end do

    end subroutine ghg_data_timestep_init


! Purpose: Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and CFC12
! 
! Method: Distributions assume constant mixing ratio in the troposphere
!         and a decrease of mixing ratio in the stratosphere. Tropopause
!         defined by ptrop. The scale height of the particular trace gas
!         depends on latitude. This assumption produces a more realistic
!         stratospheric distribution of the various trace gases.
    subroutine trcmix(name, ncol, clat, pmid, q)
    use chem_surfvals,  only: chem_surfvals_get, chem_surfvals_co2_rad
    ! io
    character(len=*), intent(in)  :: name             ! constituent name
    integer,          intent(in)  :: ncol             ! number of columns
    real(r8),         intent(in)  :: clat(ncol)       ! latitude in radians for columns
    real(r8),         intent(in)  :: pmid(nlev,ncol)  ! model pressures
    real(r8),         intent(inout) :: q(nlev,ncol)   ! constituent mass mixing ratio
    ! local
    integer :: i             ! longitude loop index
    integer :: k             ! level index
    real(r8) coslat(ncol)    ! cosine of latitude
    real(r8) dlat            ! latitude in degrees
    real(r8) ptrop           ! pressure level of tropopause
    real(r8) pratio          ! pressure divided by ptrop
    real(r8) trop_mmr        ! tropospheric mass mixing ratio
    real(r8) scale           ! pressure scale height

    do i = 1, ncol
       coslat(i) = cos(clat(i))
    end do

    if (name == 'O2') then

      q = chem_surfvals_get('O2MMR')

    else if (name == 'CO2') then

       q = chem_surfvals_co2_rad()

    else if (name == 'CH4') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwch4 * chem_surfvals_get('CH4VMR')

       do k = 1, nlev
          do i = 1, ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_r8 * clat(i))
             if(dlat.le.45.0_r8) then
                scale = 0.2353_r8
             else
                scale = 0.2353_r8 + 0.0225489_r8 * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

             ! determine output mass mixing ratios
             if (pmid(k,i) >= ptrop) then
                q(k,i) = trop_mmr
             else
                pratio = pmid(k,i)/ptrop
                q(k,i) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    else if (name == 'N2O') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwn2o * chem_surfvals_get('N2OVMR')

       do k = 1,nlev
          do i = 1,ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_r8 * clat(i))
             if(dlat.le.45.0_r8) then
                scale = 0.3478_r8 + 0.00116_r8 * dlat
             else
                scale = 0.4000_r8 + 0.013333_r8 * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

             ! determine output mass mixing ratios
             if (pmid(k,i) >= ptrop) then
                q(k,i) = trop_mmr
             else
                pratio = pmid(k,i)/ptrop
                q(k,i) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    else if (name == 'CFC11') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwf11 * chem_surfvals_get('F11VMR')

       do k = 1,nlev
          do i = 1,ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_r8 * clat(i))
             if(dlat.le.45.0_r8) then
                scale = 0.7273_r8 + 0.00606_r8 * dlat
             else
                scale = 1.00_r8 + 0.013333_r8 * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

             ! determine output mass mixing ratios
             if (pmid(k,i) >= ptrop) then
                q(k,i) = trop_mmr
             else
                pratio = pmid(k,i)/ptrop
                q(k,i) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    else if (name == 'CFC12') then

       ! set tropospheric mass mixing ratios
       trop_mmr = rmwf12 * chem_surfvals_get('F12VMR')

       do k = 1,nlev
          do i = 1,ncol
             ! set stratospheric scale height factor for gases
             dlat = abs(57.2958_r8 * clat(i))
             if(dlat.le.45.0_r8) then
                scale = 0.4000_r8 + 0.00222_r8 * dlat
             else
                scale = 0.50_r8 + 0.024444_r8 * (dlat - 45)
             end if

             ! pressure of tropopause
             ptrop = 250.0e2_r8 - 150.0e2_r8*coslat(i)**2.0_r8

             ! determine output mass mixing ratios
             if (pmid(k,i) >= ptrop) then
                q(k,i) = trop_mmr
             else
                pratio = pmid(k,i)/ptrop
                q(k,i) = trop_mmr * (pratio)**scale
             end if
          end do
       end do

    else
        return

    end if

    end subroutine trcmix

 end module ghg_data
