!======================================================
!
!  Created by LiXiaohan on 19/8/20.
!  
!======================================================

 module ebert_curry

    use grist_constants,                    only: r8, gravity
    use grist_nml_module,                   only: nlev, nlevp
    use grist_handle_error,                 only: endrun
    use radconstants,                       only: nswbands, nlwbands, ot_length,  &
                                                  idx_sw_diag, idx_lw_diag,       &
                                                  get_sw_spectral_boundaries
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_physics_update,               only: old_time_level
 
    implicit none
    private
    save
    
    public :: ec_rad_props_init,      &
!       cloud_rad_props_get_sw, & ! return SW optical props of total bulk aerosols
!       cloud_rad_props_get_lw, & ! return LW optical props of total bulk aerosols
              ec_ice_optics_sw,       &
              ec_ice_get_rad_props_lw
    
    real(r8), public, parameter:: scalefactor = 1._r8 !500._r8/917._r8

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

    subroutine ec_rad_props_init()
    use grist_nml_module,               only: ntracer
    use grist_physics_data_structure,   only: phy_tracer_info

! local
    integer :: m

    ! inquire cloud index
    do m = 1, ntracer
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')        ixcldliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')           ixcldice = m
    end do

    end subroutine ec_rad_props_init


    subroutine ec_ice_optics_sw( ncol, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f, oldicewp )

    ! io
    logical, intent(in)  :: oldicewp
    integer, intent(in)  :: ncol
    real(r8),intent(out) :: ice_tau    (nswbands,nlev,ncol) ! extinction optical depth
    real(r8),intent(out) :: ice_tau_w  (nswbands,nlev,ncol) ! single scattering albedo * tau
    real(r8),intent(out) :: ice_tau_w_g(nswbands,nlev,ncol) ! assymetry parameter * tau * w
    real(r8),intent(out) :: ice_tau_w_f(nswbands,nlev,ncol) ! forward scattered fraction * tau * w
    ! local
    real(r8), dimension(nlev,ncol) :: cicewp
    real(r8), dimension(nswbands)  :: wavmin
    real(r8), dimension(nswbands)  :: wavmax
 
    ! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
    real(r8) :: abari(4) = &     ! a coefficient for extinction optical depth
       (/ 3.448e-03_r8, 3.448e-03_r8,3.448e-03_r8,3.448e-03_r8/)
    real(r8) :: bbari(4) = &     ! b coefficient for extinction optical depth
       (/ 2.431_r8    , 2.431_r8    ,2.431_r8    ,2.431_r8    /)
    real(r8) :: cbari(4) = &     ! c coefficient for single scat albedo
       (/ 1.00e-05_r8 , 1.10e-04_r8 ,1.861e-02_r8,.46658_r8   /)
    real(r8) :: dbari(4) = &     ! d coefficient for single scat albedo
       (/ 0.0_r8      , 1.405e-05_r8,8.328e-04_r8,2.05e-05_r8 /)
    real(r8) :: ebari(4) = &     ! e coefficient for asymmetry parameter
       (/ 0.7661_r8   , 0.7730_r8   ,0.794_r8    ,0.9595_r8   /)
    real(r8) :: fbari(4) = &     ! f coefficient for asymmetry parameter
       (/ 5.851e-04_r8, 5.665e-04_r8,7.267e-04_r8,1.076e-04_r8/)

    real(r8) :: abarii           ! A coefficient for current spectral band
    real(r8) :: bbarii           ! B coefficient for current spectral band
    real(r8) :: cbarii           ! C coefficient for current spectral band
    real(r8) :: dbarii           ! D coefficient for current spectral band
    real(r8) :: ebarii           ! E coefficient for current spectral band
    real(r8) :: fbarii           ! F coefficient for current spectral band

    ! Minimum cloud amount (as a fraction of the grid-box area) to 
    ! distinguish from clear sky
    real(r8), parameter :: cldmin = 1.0e-80_r8

    ! Decimal precision of cloud amount (0 -> preserve full resolution;
    ! 10^-n -> preserve n digits of cloud amount)
    real(r8), parameter :: cldeps = 0.0_r8

    integer  :: ns, i, k, indxsl, lchnk, Nday
    integer  :: itim
    real(r8) :: tmp1i, tmp2i, tmp3i, g

    Nday = ncol

    if(oldicewp) then
      do k=1, nlev
         do i = 1, Nday
            cicewp(k,i) = 1000.0_r8*pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,k,i)     &
                         *pstate%delp_at_pc_full_level%f(k,i)                               &
                         /(gravity * max(0.01,pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i)))   
         end do
      end do
    else
      if(.not.(allocated(pstate_cam%microp_iciwp_at_pc_full_level%f)))then
          call endrun('ec_ice_optics_sw: oldicewp must be set to true since ICIWP was not found')
      end if  
      cicewp(1:nlev,1:ncol) =  1000.0_r8*pstate_cam%microp_iciwp_at_pc_full_level%f(1:nlev,1:ncol)
    endif
 
    call get_sw_spectral_boundaries(wavmin,wavmax,'microns')

    do ns = 1, nswbands

       if(wavmax(ns) <= 0.7_r8) then
          indxsl = 1
       else if(wavmax(ns) <= 1.25_r8) then
          indxsl = 2
       else if(wavmax(ns) <= 2.38_r8) then
          indxsl = 3
       else if(wavmax(ns) > 2.38_r8) then
          indxsl = 4
       end if

       abarii = abari(indxsl)
       bbarii = bbari(indxsl)
       cbarii = cbari(indxsl)
       dbarii = dbari(indxsl)
       ebarii = ebari(indxsl)
       fbarii = fbari(indxsl)

       do k=1, nlev
          do i=1, Nday
 
             ! note that optical properties for ice valid only
             ! in range of 13 > rei > 130 micron (Ebert and Curry 92)
             if (  pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i) >= cldmin     &  
             .and. pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,k,i) >= cldeps) then
                tmp1i = abarii +                                                          &
                        bbarii/max(13._r8,min(scalefactor*pstate_cam%microp_rei_at_pc_full_level%f(k,i),130._r8))
                ice_tau(ns,k,i) = cicewp(k,i)*tmp1i
             else
                ice_tau(ns,k,i) = 0.0_r8
             endif
 
             tmp2i = 1._r8 - cbarii -                                                     &  
                     dbarii*min(max(13._r8,scalefactor*pstate_cam%microp_rei_at_pc_full_level%f(k,i)),130._r8)
             tmp3i = fbarii*min(max(13._r8,scalefactor*pstate_cam%microp_rei_at_pc_full_level%f(k,i)),130._r8)
             ! Do not let single scatter albedo be 1.  Delta-eddington solution
             ! for non-conservative case has different analytic form from solution
             ! for conservative case, and raddedmx is written for non-conservative case.
             ice_tau_w(ns,k,i) = ice_tau(ns,k,i) * min(tmp2i,.999999_r8)
             g = ebarii + tmp3i
             ice_tau_w_g(ns,k,i) = ice_tau_w(ns,k,i) * g
             ice_tau_w_f(ns,k,i) = ice_tau_w(ns,k,i) * g * g
 
          end do ! End do i=1, Nday
       end do    ! End do k=1, nlev
    end do ! nswbands

    end subroutine ec_ice_optics_sw


    subroutine ec_ice_get_rad_props_lw(ncol, abs_od, oldicewp)
    ! io
    integer, intent(in)   :: ncol
    logical, intent(in)   :: oldicewp
    real(r8), intent(out) :: abs_od(nlwbands,nlev,ncol)

    real(r8) :: gicewp(nlev,ncol)
    real(r8) :: gliqwp(nlev,ncol)
    real(r8) :: cicewp(nlev,ncol)
    real(r8) :: cliqwp(nlev,ncol)
    real(r8) :: ficemr(nlev,ncol)
    real(r8) :: cwp(nlev,ncol)
    real(r8) :: cldtau(nlev,ncol)

    real(r8) :: cldn(nlev,ncol), rei(nlev,ncol)
    real(r8) :: iclwpth(nlev,ncol), iciwpth(nlev,ncol)
    
    integer  :: lwband, i, k
    real(r8) :: kabs, kabsi
    real(r8), parameter :: kabsl= 0.090361_r8       ! longwave liquid absorption coeff (m**2/g)

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
            cwp(k,i) = 1000.0_r8 *iciwpth(k,i) + 1000.0_r8 *iclwpth(k,i)
            ficemr(k,i) = 1000.0_r8*iciwpth(k,i)/(max(1.e-18_r8,cwp(k,i)))
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

    !output: cicewp(:,:)/1000, iciwpth, cldtau

    end subroutine ec_ice_get_rad_props_lw

 end module ebert_curry   
