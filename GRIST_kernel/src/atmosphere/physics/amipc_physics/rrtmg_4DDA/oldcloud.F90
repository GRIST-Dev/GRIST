!=================================================================
!
!  Created by LiXiaohan on 19/8/20.
!  
!=================================================================

 module oldcloud

    use grist_constants,                    only: r8, gravity
    use grist_nml_module,                   only: nlev, nlevp
    use grist_handle_error,                 only: endrun
    use radconstants,                       only: nswbands, nlwbands, ot_length,  &
                                                  idx_sw_diag, idx_lw_diag,       &
                                                  get_sw_spectral_boundaries
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_physics_update,               only: old_time_level
    use ebert_curry,                        only: scalefactor

    implicit none
    private
    save
    
    public :: oldcloud_init,                &
              oldcloud_lw,                  &
              old_liq_get_rad_props_lw,     &
              old_ice_get_rad_props_lw

    integer :: nmu, nlambda
    real(r8), allocatable :: g_mu(:)           ! mu samples on grid
    real(r8), allocatable :: g_lambda(:,:)     ! lambda scale samples on grid
    real(r8), allocatable :: ext_sw_liq(:,:,:)
    real(r8), allocatable :: ssa_sw_liq(:,:,:)
    real(r8), allocatable :: asm_sw_liq(:,:,:)
    real(r8), allocatable :: abs_lw_liq(:,:,:)
    
    integer :: n_g_d
    real(r8), allocatable :: g_d_eff(:)        ! radiative effective diameter samples on grid
    real(r8), allocatable :: ext_sw_ice(:,:)
    real(r8), allocatable :: ssa_sw_ice(:,:)
    real(r8), allocatable :: asm_sw_ice(:,:)
    real(r8), allocatable :: abs_lw_ice(:,:)

! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
 
    real(r8), parameter :: cldmin = 1.0e-80_r8

! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
 
    real(r8), parameter :: cldeps = 0._r8

! indexes into constituents for old optics
    integer :: ixcldice,           & ! cloud ice water index
               ixcldliq              ! cloud liquid water index

 contains

    subroutine oldcloud_init()
    use grist_nml_module,               only: ntracer
    use grist_physics_data_structure,   only: phy_tracer_info

! local
    integer :: m

    ! inquire cloud index
    do m = 1, ntracer
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')        ixcldliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')           ixcldice = m
    end do

    end subroutine oldcloud_init


    subroutine oldcloud_lw(ncol,cld_abs_od,oldwp)
    ! io
    integer, intent(in)  :: ncol
    logical,intent(in)   :: oldwp ! use old definition of waterpath
    real(r8),intent(out) :: cld_abs_od(nlwbands,nlev,ncol) ! [fraction] absorption optical depth, per layer
    ! local
    real(r8) :: gicewp(nlev,ncol)
    real(r8) :: gliqwp(nlev,ncol)
    real(r8) :: cicewp(nlev,ncol)
    real(r8) :: cliqwp(nlev,ncol)
    real(r8) :: ficemr(nlev,ncol)
    real(r8) :: cwp(nlev,ncol)
    real(r8) :: cldtau(nlev,ncol)

    real(r8) :: cldn(nlev,ncol)
    real(r8) :: rei(nlev,ncol) 
    real(r8) :: iclwpth(nlev,ncol), iciwpth(nlev,ncol)
    integer  :: lwband, i, k

    real(r8) :: kabs, kabsi
    real(r8), parameter :: kabsl = 0.090361_r8     ! longwave liquid absorption coeff (m**2/g)

    rei(:,1:ncol)  = pstate_cam%microp_rei_at_pc_full_level%f(:,1:ncol)
    cldn(:,1:ncol) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)

    if (oldwp) then
      do k=1, nlev
          do i = 1, ncol
             gicewp(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)    &
                          *pstate%delp_at_pc_full_level%f(k,i)/gravity*1000.0_r8  ! Grid box ice water path.
             gliqwp(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,k,i)    &
                          *pstate%delp_at_pc_full_level%f(k,i)/gravity*1000.0_r8  ! Grid box liquid water path.
             cicewp(k,i) = gicewp(k,i) / max(0.01_r8,cldn(k,i))                   ! In-cloud ice water path.
             cliqwp(k,i) = gliqwp(k,i) / max(0.01_r8,cldn(k,i))                   ! In-cloud liquid water path.
             ficemr(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)                   &
                          /max(1.e-10_r8,(pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)    &
                          +pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,k,i)))
          end do
      end do
      cwp(:nlev,:ncol) = cicewp(:nlev,:ncol) + cliqwp(:nlev,:ncol)
    else
      iclwpth(:,1:ncol) = pstate_cam%microp_iclwp_at_pc_full_level%f(:,1:ncol)
      iciwpth(:,1:ncol) = pstate_cam%microp_iciwp_at_pc_full_level%f(:,1:ncol)
      do k=1, nlev
         do i = 1, ncol
            cwp(k,i) = 1000.0_r8 *iclwpth(k,i) + 1000.0_r8 *iciwpth(k,i)
            ficemr(k,i) = 1000.0_r8 * iciwpth(k,i)/(max(1.e-18_r8,cwp(k,i)))
         end do
      end do
    endif

    do k=1, nlev
        do i=1, ncol
           !note that optical properties for ice valid only
           !in range of 13 > rei > 130 micron (Ebert and Curry 92)
           kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(k,i)),130._r8)
           kabs = kabsl*(1._r8-ficemr(k,i)) + kabsi*ficemr(k,i)
           !emis(k,i) = 1._r8 - exp(-1.66_r8*kabs*clwp(k,i))
           cldtau(k,i) = kabs*cwp(k,i)
        end do
    end do

    do lwband = 1,nlwbands
       cld_abs_od(lwband,1:nlev,1:ncol)=cldtau(1:nlev,1:ncol)
    enddo

    end subroutine oldcloud_lw


    subroutine old_liq_get_rad_props_lw(ncol, abs_od, oldliqwp)
    ! io
    integer, intent(in)   :: ncol
    logical, intent(in)   :: oldliqwp
    real(r8), intent(out) :: abs_od(nlwbands,nlev,ncol)
    ! local
    real(r8) :: gicewp(nlev,ncol)
    real(r8) :: gliqwp(nlev,ncol)
    real(r8) :: cicewp(nlev,ncol)
    real(r8) :: cliqwp(nlev,ncol)
    real(r8) :: ficemr(nlev,ncol)
    real(r8) :: cwp(nlev,ncol)
    real(r8) :: cldtau(nlev,ncol)
   
    real(r8) :: cldn(nlev,ncol)
    real(r8) :: rei(nlev,ncol)
    real(r8) :: iclwpth(nlev,ncol)
    real(r8) :: iciwpth(nlev,ncol)
    integer  :: lwband, i, k

    real(r8) :: kabs, kabsi
    real(r8), parameter :: kabsl = 0.090361_r8    ! longwave liquid absorption coeff (m**2/g)

    rei(:,1:ncol)  = pstate_cam%microp_rei_at_pc_full_level%f(:,1:ncol)
    cldn(:,1:ncol) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)

    if (oldliqwp) then
      do k=1, nlev
          do i = 1, ncol
             gicewp(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)    &
                          *pstate%delp_at_pc_full_level%f(k,i)/gravity*1000.0_r8  ! Grid box ice water path.
             gliqwp(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,k,i)    &
                          *pstate%delp_at_pc_full_level%f(k,i)/gravity*1000.0_r8  ! Grid box liquid water path.
             cicewp(k,i) = gicewp(k,i) / max(0.01_r8,cldn(k,i))                   ! In-cloud ice water path.
             cliqwp(k,i) = gliqwp(k,i) / max(0.01_r8,cldn(k,i))                   ! In-cloud liquid water path.
             ficemr(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)                   &
                          /max(1.e-10_r8,(pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)    &
                          +pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,k,i)))
          end do
      end do
      cwp(:nlev,:ncol) = cicewp(:nlev,:ncol) + cliqwp(:nlev,:ncol)
    else
      iclwpth(:,1:ncol) = pstate_cam%microp_iclwp_at_pc_full_level%f(:,1:ncol)
      iciwpth(:,1:ncol) = pstate_cam%microp_iciwp_at_pc_full_level%f(:,1:ncol)
 
      do k=1, nlev
         do i = 1, ncol
            cwp(k,i) = 1000.0_r8 *iclwpth(k,i) + 1000.0_r8 *iciwpth(k,i)
            ficemr(k,i) = 1000.0 * iciwpth(k,i)/(max(1.e-18_r8,cwp(k,i)))
         end do
      end do
    endif

    do k=1, nlev
        do i=1, ncol
           ! Note from Andrew Conley:
           !  Optics for RK no longer supported, This is constructed to get
           !  close to bit for bit.  Otherwise we could simply use liquid water path
           !note that optical properties for ice valid only
           !in range of 13 > rei > 130 micron (Ebert and Curry 92)
           kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(k,i)),130._r8)
           kabs = kabsl*(1._r8-ficemr(k,i)) ! + kabsi*ficemr(k,i)
           !emis(k,i) = 1._r8 - exp(-1.66_r8*kabs*clwp(k,i))
           cldtau(k,i) = kabs*cwp(k,i)
        end do
    end do

    do lwband = 1,nlwbands
       abs_od(lwband,1:nlev,1:ncol)=cldtau(1:nlev,1:ncol)
    enddo

    end subroutine old_liq_get_rad_props_lw


    subroutine old_ice_get_rad_props_lw(ncol, abs_od, oldicewp)
    ! io
    integer, intent(in)   :: ncol
    logical, intent(in) :: oldicewp
    real(r8), intent(out) :: abs_od(nlwbands,nlev,ncol)

    real(r8) :: gicewp(nlev,ncol)
    real(r8) :: gliqwp(nlev,ncol)
    real(r8) :: cicewp(nlev,ncol)
    real(r8) :: cliqwp(nlev,ncol)
    real(r8) :: ficemr(nlev,ncol)
    real(r8) :: cwp(nlev,ncol)
    real(r8) :: cldtau(nlev,ncol)

    real(r8) :: cldn(nlev,ncol)
    real(r8) :: rei(nlev,ncol)
    real(r8) :: iclwpth(nlev,ncol)
    real(r8) :: iciwpth(nlev,ncol)
 
    integer  :: lwband, i, k
    real(r8) :: kabs, kabsi

    real(r8), parameter :: kabsl = 0.090361_r8  ! longwave liquid absorption coeff (m**2/g)

    rei(:,1:ncol)  = pstate_cam%microp_rei_at_pc_full_level%f(:,1:ncol)
    cldn(:,1:ncol) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)

    if(oldicewp) then
      do k=1, nlev
          do i = 1, ncol
             gicewp(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)    &
                          *pstate%delp_at_pc_full_level%f(k,i)/gravity*1000.0_r8  ! Grid box ice water path.
             gliqwp(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,k,i)    &  
                          *pstate%delp_at_pc_full_level%f(k,i)/gravity*1000.0_r8  ! Grid box liquid water path.
             cicewp(k,i) = gicewp(k,i) / max(0.01_r8,cldn(k,i))                   ! In-cloud ice water path.
             cliqwp(k,i) = gliqwp(k,i) / max(0.01_r8,cldn(k,i))                   ! In-cloud liquid water path.
             ficemr(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)                 &
                          /max(1.e-10_r8,(pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)  &
                          +pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,k,i)))
          end do
      end do
      cwp(:nlev,:ncol)  = cicewp(:nlev,:ncol)  + cliqwp(:nlev,:ncol) 
    else
      iclwpth(:,1:ncol) = pstate_cam%microp_iclwp_at_pc_full_level%f(:,1:ncol)
      iciwpth(:,1:ncol) = pstate_cam%microp_iciwp_at_pc_full_level%f(:,1:ncol)
 
      do k=1, nlev
         do i = 1, ncol
            cwp(k,i) = 1000.0_r8 *iciwpth(k,i) + 1000.0_r8 *iclwpth(k,i)
            ficemr(k,i) = 1000.0 * iciwpth(k,i)/(max(1.e-18_r8,cwp(k,i)))
         end do
      end do
    endif

    do k=1, nlev
        do i=1, ncol

           ! Note from Andrew Conley:
           !  Optics for RK no longer supported, This is constructed to get
           !  close to bit for bit.  Otherwise we could simply use ice water path
           !note that optical properties for ice valid only
           !in range of 13 > rei > 130 micron (Ebert and Curry 92)
           kabsi = 0.005_r8 + 1._r8/min(max(13._r8,scalefactor*rei(k,i)),130._r8)
           kabs =  kabsi*ficemr(k,i) ! kabsl*(1._r8-ficemr(k,i)) + kabsi*ficemr(k,i)
           !emis(k,i) = 1._r8 - exp(-1.66_r8*kabs*clwp(k,i))
           cldtau(k,i) = kabs*cwp(k,i)
        end do
    end do

    do lwband = 1,nlwbands
       abs_od(lwband,1:nlev,1:ncol)=cldtau(1:nlev,1:ncol)
    enddo

    ! output: cicewp, iciwpth, cldtau
    end subroutine old_ice_get_rad_props_lw

 end module oldcloud
