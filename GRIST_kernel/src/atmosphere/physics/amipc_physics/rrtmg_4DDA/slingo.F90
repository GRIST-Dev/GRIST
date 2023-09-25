!=================================================================
!
!  Created by LiXiaohan on 19/8/20.
!  
!  Implements Slingo Optics for MG/RRTMG for liquid clouds and
!  a copy of the old cloud routine for reference 
!=================================================================

 module slingo

    use grist_constants,                    only: r8, gravity
    use grist_nml_module,                   only: nlev, nlevp
    use grist_handle_error,                 only: endrun
    use radconstants,                       only: nswbands, nlwbands, ot_length,  &
                                                  idx_sw_diag, idx_lw_diag,       &
                                                  get_sw_spectral_boundaries
    use grist_physics_update,               only: old_time_level
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam

    implicit none
    private
    save
    
    public :: slingo_rad_props_init,        &
!   cloud_rad_props_get_sw, & ! return SW optical props of total bulk aerosols
!   cloud_rad_props_get_lw,  & ! return LW optical props of total bulk aerosols
              slingo_liq_get_rad_props_lw,  &
              slingo_liq_optics_sw


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

    subroutine slingo_rad_props_init()
    use grist_nml_module,               only: ntracer
    use grist_physics_data_structure,   only: phy_tracer_info

! local
    integer :: m

    ! inquire cloud index
    do m = 1, ntracer
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')        ixcldliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')           ixcldice = m
    end do

    end subroutine slingo_rad_props_init


    subroutine slingo_liq_optics_sw(ncol, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f, oldliqwp)
    ! io
    integer, intent(in) :: ncol
    real(r8),intent(out) :: liq_tau    (nswbands,nlev,ncol) ! extinction optical depth
    real(r8),intent(out) :: liq_tau_w  (nswbands,nlev,ncol) ! single scattering albedo * tau
    real(r8),intent(out) :: liq_tau_w_g(nswbands,nlev,ncol) ! assymetry parameter * tau * w
    real(r8),intent(out) :: liq_tau_w_f(nswbands,nlev,ncol) ! forward scattered fraction * tau * w
    logical, intent(in) :: oldliqwp
    ! local
    real(r8) :: rel(nlev,ncol), cldn(nlev,ncol), tmpptr(nlev,ncol), iclwpth(nlev,ncol)
    real(r8), dimension(nlev,ncol) :: cliqwp
    real(r8), dimension(nswbands) :: wavmin
    real(r8), dimension(nswbands) :: wavmax

    ! Minimum cloud amount (as a fraction of the grid-box area) to 
    ! distinguish from clear sky
    real(r8), parameter :: cldmin = 1.0e-80_r8

    ! Decimal precision of cloud amount (0 -> preserve full resolution;
    ! 10^-n -> preserve n digits of cloud amount)
    real(r8), parameter :: cldeps = 0.0_r8

    ! A. Slingo's data for cloud particle radiative properties (from 'A GCM
    ! Parameterization for the Shortwave Properties of Water Clouds' JAS
    ! vol. 46 may 1989 pp 1419-1427)
    real(r8) :: abarl(4) = &  ! A coefficient for extinction optical depth
       (/ 2.817e-02_r8, 2.682e-02_r8,2.264e-02_r8,1.281e-02_r8/)
    real(r8) :: bbarl(4) = &  ! B coefficient for extinction optical depth
       (/ 1.305_r8    , 1.346_r8    ,1.454_r8    ,1.641_r8    /)
    real(r8) :: cbarl(4) = &  ! C coefficient for single scat albedo
       (/-5.62e-08_r8 ,-6.94e-06_r8 ,4.64e-04_r8 ,0.201_r8    /)
    real(r8) :: dbarl(4) = &  ! D coefficient for single  scat albedo
       (/ 1.63e-07_r8 , 2.35e-05_r8 ,1.24e-03_r8 ,7.56e-03_r8 /)
    real(r8) :: ebarl(4) = &  ! E coefficient for asymmetry parameter
       (/ 0.829_r8    , 0.794_r8    ,0.754_r8    ,0.826_r8    /)
    real(r8) :: fbarl(4) = &  ! F coefficient for asymmetry parameter
       (/ 2.482e-03_r8, 4.226e-03_r8,6.560e-03_r8,4.353e-03_r8/)

    real(r8) :: abarli        ! A coefficient for current spectral band
    real(r8) :: bbarli        ! B coefficient for current spectral band
    real(r8) :: cbarli        ! C coefficient for current spectral band
    real(r8) :: dbarli        ! D coefficient for current spectral band
    real(r8) :: ebarli        ! E coefficient for current spectral band
    real(r8) :: fbarli        ! F coefficient for current spectral band

    ! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
    ! greater than 20 micro-meters

    integer :: ns, i, k, indxsl, Nday
    real(r8) :: tmp1l, tmp2l, tmp3l, g
    real(r8) :: kext(nlev,ncol)

    Nday = ncol
    rel(:,1:ncol)  = pstate_cam%microp_rel_at_pc_full_level%f(:,1:ncol)
    cldn(:,1:ncol) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,1:ncol)

    if (oldliqwp) then
      do k=1, nlev
         do i = 1, Nday
            cliqwp(k,i) = pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,k,i)                   & 
                         *pstate%delp_at_pc_full_level%f(k,i)/(gravity*max(0.01_r8,cldn(k,i)))
         end do
      end do
    else
      tmpptr(:,1:ncol) = pstate_cam%microp_iclwp_at_pc_full_level%f(:,1:ncol)
      ! The following is the eventual target specification for in cloud liquid water path.
      cliqwp = tmpptr
    endif
   
    call get_sw_spectral_boundaries(wavmin,wavmax,'microns')
  
    do ns = 1, nswbands
       ! Set index for cloud particle properties based on the wavelength,
       ! according to A. Slingo (1989) equations 1-3:
       ! Use index 1 (0.25 to 0.69 micrometers) for visible
       ! Use index 2 (0.69 - 1.19 micrometers) for near-infrared
       ! Use index 3 (1.19 to 2.38 micrometers) for near-infrared
       ! Use index 4 (2.38 to 4.00 micrometers) for near-infrared
       if(wavmax(ns) <= 0.7_r8) then
          indxsl = 1
       else if(wavmax(ns) <= 1.25_r8) then
          indxsl = 2
       else if(wavmax(ns) <= 2.38_r8) then
          indxsl = 3
       else if(wavmax(ns) > 2.38_r8) then
          indxsl = 4
       end if

       ! Set cloud extinction optical depth, single scatter albedo,
       ! asymmetry parameter, and forward scattered fraction:
       abarli = abarl(indxsl)
       bbarli = bbarl(indxsl)
       cbarli = cbarl(indxsl)
       dbarli = dbarl(indxsl)
       ebarli = ebarl(indxsl)
       fbarli = fbarl(indxsl)

       do k=1, nlev
          do i=1, Nday

             ! note that optical properties for liquid valid only
             ! in range of 4.2 > rel > 16 micron (Slingo 89)
             if (cldn(k,i) >= cldmin .and. cldn(k,i) >= cldeps) then
                tmp1l = abarli + bbarli/min(max(4.2_r8,rel(k,i)),16._r8)
                liq_tau(ns,k,i) = 1000._r8*cliqwp(k,i)*tmp1l
             else
                liq_tau(ns,k,i) = 0.0_r8
             endif

             tmp2l = 1._r8 - cbarli - dbarli*min(max(4.2_r8,rel(k,i)),16._r8)
             tmp3l = fbarli*min(max(4.2_r8,rel(k,i)),16._r8)
             ! Do not let single scatter albedo be 1.  Delta-eddington solution
             ! for non-conservative case has different analytic form from solution
             ! for conservative case, and raddedmx is written for non-conservative case.
             liq_tau_w(ns,k,i) = liq_tau(ns,k,i) * min(tmp2l,.999999_r8)
             g = ebarli + tmp3l
             liq_tau_w_g(ns,k,i) = liq_tau_w(ns,k,i) * g
             liq_tau_w_f(ns,k,i) = liq_tau_w(ns,k,i) * g * g

          end do ! End do i=1, Nday
       end do    ! End do k=1, nlev
    end do ! nswbands

    ! output: liq_tau(idx_sw_diag,:,:), rel, cliqwp, kext
    end subroutine slingo_liq_optics_sw


    subroutine slingo_liq_get_rad_props_lw(ncol, abs_od, oldliqwp)
    ! io
    integer, intent(in) :: ncol
    real(r8), intent(out) :: abs_od(nlwbands,nlev,ncol)
    logical, intent(in) :: oldliqwp

    real(r8) :: gicewp(nlev,ncol)
    real(r8) :: gliqwp(nlev,ncol)
    real(r8) :: cicewp(nlev,ncol)
    real(r8) :: cliqwp(nlev,ncol)
    real(r8) :: ficemr(nlev,ncol)
    real(r8) :: cwp(nlev,ncol)
    real(r8) :: cldtau(nlev,ncol)

    real(r8) :: cldn(nlev,ncol), rei(nlev,ncol), iclwpth(nlev,ncol), iciwpth(nlev,ncol)
    integer :: lwband, i, k

    real(r8) :: kabs, kabsi
    real(r8), parameter :: kabsl = 0.090361_r8       ! longwave liquid absorption coeff (m**2/g)

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
               cwp   (k,i) = 1000.0_r8 * iclwpth(k,i) + 1000.0_r8 * iciwpth(k,i)
               ficemr(k,i) = 1000.0_r8 * iciwpth(k,i)/(max(1.e-18_r8, cwp(k,i)))
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
           kabs = kabsl*(1._r8-ficemr(k,i))
           cldtau(k,i) = kabs*cwp(k,i)
        end do
    end do

    do lwband = 1,nlwbands
       abs_od(lwband,1:nlev,1:ncol)=cldtau(1:nlev,1:ncol)
    enddo

    end subroutine slingo_liq_get_rad_props_lw



 end module slingo
