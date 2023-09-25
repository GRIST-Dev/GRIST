!===================================================================================================
!
!  Created by LiXiaohan on 19/07/10, adopted from CAM5
!
!     interface of the aerosol activation
!
!     This module is currently hardcoded to recognize only the aerosols/modes that
!     affect the climate calculation.  This is implemented by using list
!     index 0 in all the calls to rad_constituent interfaces.
!
!     Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                     Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
!     Modifications: A. Gettelman Nov 2010  - changed to support separation of 
!                    microphysics and macrophysics and concentrate aerosol information here
!
!===================================================================================================

 module grist_microp_aero

    use grist_handle_error,             only: endrun
    use grist_constants,                only: r8, i4,       &
                                              rdry, tmelt
    use grist_physics_data_structure,   only: pstate
    use grist_PhysC_data_structure,      only: pstate_cam
    use grist_nml_module,               only: nlev, nlevp, ntracer
    use cloud_fraction,                 only: top_lev=>trop_cloud_top_lev
    use grist_rad_constituents,         only: rad_cnst_get_info,        &
                                              rad_cnst_get_aer_mmr,     &
                                              rad_cnst_get_aer_props,   &
                                              rad_cnst_get_mode_num,    &
                                              rad_cnst_get_mode_props  
    use grist_ndrop,                    only: ndrop_init, dropmixnuc,   &
                                              end_of_ndrop  
    use grist_ndrop_bam,                only: ndrop_bam_init,   &
                                              ndrop_bam_run,    &
                                              ndrop_bam_ccn
    use nucleate_ice,                   only: nucleati

    implicit none
    private
    save

    public          ::  microp_aero_init,               &
                        microp_aero_driver,             &
                        end_of_microp_aero,             &
                        read_nml_microp_aero,           &
                        bulk_scale

    ! Private:
    character(len=16)   :: eddy_scheme  ! eddy scheme

    ! contact freezing due to dust
    ! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
    real(r8), parameter :: rn_dst1 = 0.258e-6_r8
    real(r8), parameter :: rn_dst2 = 0.717e-6_r8
    real(r8), parameter :: rn_dst3 = 1.576e-6_r8
    real(r8), parameter :: rn_dst4 = 3.026e-6_r8

    real(r8) :: bulk_scale                          ! prescribed aerosol bulk sulfur scale factor

    ! smallest mixing ratio considered in microphysics
    real(r8), parameter :: qsmall = 1.e-18_r8
    
    ! minimum allowed cloud fraction
    real(r8), parameter :: mincld = 0.0001_r8

    ! indices in state%q
    integer :: ixcldliq, ixcldice, ixnumliq, ixnumice

    ! Bulk aerosols
    character(len=20), allocatable :: aername(:)
    real(r8), allocatable :: num_to_mass_aer(:)

    integer :: naer_all      ! number of aerosols affecting climate
    integer :: idxsul   = -1 ! index in aerosol list for sulfate
    integer :: idxdst1  = -1 ! index in aerosol list for dust1
    integer :: idxdst2  = -1 ! index in aerosol list for dust2
    integer :: idxdst3  = -1 ! index in aerosol list for dust3
    integer :: idxdst4  = -1 ! index in aerosol list for dust4
    integer :: idxbcphi = -1 ! index in aerosol list for Soot (BCPHIL)


    ! modal aerosols
    logical :: clim_modal_aero

    integer :: mode_accum_idx  = -1  ! index of accumulation mode
    integer :: mode_aitken_idx = -1  ! index of aitken mode
    integer :: mode_coarse_idx = -1  ! index of coarse mode
    integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
    integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
    integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
    integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode

    real(r8) :: sigmag_aitken
    logical  :: separate_dust = .false.



 contains

! Purpose: Initialize constants for aerosols needed by microphysics
    subroutine microp_aero_init(ncol)
    use grist_physics_data_structure,   only: phy_tracer_info
    use phys_control,                   only: phys_getopts

    use grist_mpi
! io
    integer , intent(in)  :: ncol
! local
    integer :: iaer
    integer :: m, n, nmodes, nspec
    character(len=32) :: str32
    character(len=*), parameter :: routine = 'microp_aero_init'

    call phys_getopts(eddy_scheme_out = eddy_scheme)

    ! inquire cloud index
    do m = 1, ntracer
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid')        ixcldliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice')           ixcldice = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_liquid_number') ixnumliq = m
        if(phy_tracer_info(m)%longname .eq. 'cloud_ice_number')    ixnumice = m 
    end do

    allocate(pstate_cam%microp_areo_naai%f(nlev,ncol))
    allocate(pstate_cam%microp_areo_naai_hom%f(nlev,ncol))
    allocate(pstate_cam%microp_areo_npccn%f(nlev,ncol))
    allocate(pstate_cam%microp_areo_rndst%f(4,nlev,ncol))
    allocate(pstate_cam%microp_areo_nacon%f(4,nlev,ncol))

    pstate_cam%microp_areo_naai%pos              = 0
    pstate_cam%microp_areo_naai_hom%pos          = 0
    pstate_cam%microp_areo_npccn%pos             = 0
    pstate_cam%microp_areo_rndst%pos             = 0
    pstate_cam%microp_areo_nacon%pos             = 0

    pstate_cam%microp_areo_naai%f                = 0._r8
    pstate_cam%microp_areo_naai_hom%f            = 0._r8
    pstate_cam%microp_areo_npccn%f               = 0._r8
    pstate_cam%microp_areo_rndst%f               = 0._r8
    pstate_cam%microp_areo_nacon%f               = 0._r8

    ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
    ! The modal aerosols can be either prognostic or prescribed.
    call rad_cnst_get_info(0, nmodes=nmodes)

    clim_modal_aero = (nmodes > 0)

!---------------LiXH Test---------------
! LiXH remove the effect of aerosol
! clim_modal_aero = .false.
!---------------LiXH Test---------------

    if (clim_modal_aero) then
      call ndrop_init()

      ! mode index for specified mode types
      do m = 1, nmodes
         call rad_cnst_get_info(0, m, mode_type=str32)
         select case (trim(str32))
         case ('accum')
            mode_accum_idx = m
         case ('aitken')
            mode_aitken_idx = m
         case ('coarse')
            mode_coarse_idx = m
         case ('coarse_dust')
            mode_coarse_dst_idx = m
         case ('coarse_seasalt')
            mode_coarse_slt_idx = m
         end select
      end do

      ! check if coarse dust is in separate mode
      separate_dust = mode_coarse_dst_idx > 0

      ! for 3-mode 
      if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
      if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_dst_idx == -1.or. mode_coarse_slt_idx == -1) then
          if(mpi_rank()==0)print*, routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

      ! species indices for specified types
      ! find indices for the dust and seasalt species in the coarse mode
      call rad_cnst_get_info(0, mode_coarse_dst_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust')
            coarse_dust_idx = n
         end select
      end do
      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do

      ! Check that required mode specie types were found
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1) then
         if(mpi_rank()==0)print*, routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

      ! get specific mode properties
      call rad_cnst_get_mode_props(0, mode_aitken_idx, sigmag=sigmag_aitken)

    else

        ! Props needed for BAM number concentration calcs.
        call rad_cnst_get_info(0, naero=naer_all)

!---------------LiXH Test---------------
! LiXH remove the effect of aerosol
!naer_all = 0
!---------------LiXH Test---------------


        !-----------LiXH 2020-10-21-------------
        if(naer_all .gt. 0)then
        allocate( aername(naer_all),        &
                  num_to_mass_aer(naer_all) )
        end if
        !-----------LiXH 2020-10-21-------------
        
        do iaer = 1, naer_all
           call rad_cnst_get_aer_props(0, iaer, &
              aername         = aername(iaer), &
              num_to_mass_aer = num_to_mass_aer(iaer) )

           ! Look for sulfate, dust, and soot in this list (Bulk aerosol only)
           if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
           if (trim(aername(iaer)) == 'DUST1') idxdst1 = iaer
           if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
           if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
           if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
           if (trim(aername(iaer)) == 'BCPHIL') idxbcphi = iaer
        end do

        call ndrop_bam_init()

    end if

    end subroutine microp_aero_init



    subroutine read_nml_microp_aero(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer :: unitn, ierr
    ! Namelist variables
    real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor
  
  
    namelist /microp_aero_nl/ microp_aero_bulk_scale
    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, microp_aero_nl, iostat=ierr)
    if (ierr /= 0) call endrun(" error reading micro_mg namelist")
    close(unitn)
 
    bulk_scale = microp_aero_bulk_scale

    end subroutine read_nml_microp_aero


    subroutine end_of_microp_aero

    if(allocated(pstate_cam%microp_areo_naai%f))          &
      deallocate(pstate_cam%microp_areo_naai%f)
 
    if(allocated(pstate_cam%microp_areo_naai_hom%f))      &
      deallocate(pstate_cam%microp_areo_naai_hom%f)
 
    if(allocated(pstate_cam%microp_areo_npccn%f))         &
      deallocate(pstate_cam%microp_areo_npccn%f)
 
    if(allocated(pstate_cam%microp_areo_rndst%f))         &
      deallocate(pstate_cam%microp_areo_rndst%f)
 
    if(allocated(pstate_cam%microp_areo_nacon%f))         &
      deallocate(pstate_cam%microp_areo_nacon%f)
 
    if(allocated(aername))          deallocate(aername)
    if(allocated(num_to_mass_aer))  deallocate(num_to_mass_aer)

    call end_of_ndrop
    end subroutine end_of_microp_aero

    
    subroutine microp_aero_driver(ncol, deltatin)
    use grist_physics_update,           only: old_time_level,           &
                                              ptimelevels
    use grist_wv_saturation,            only: qsat_water
    use grist_PhysC_data_structure,     only: ptend_microp_aero
    use grist_mpi
! io
    integer , intent(in)  :: ncol
    real(r8), intent(in)  :: deltatin
! local
    integer  :: i, k, m
    real(r8) :: icecldf(nlev,ncol)    ! ice cloud fraction   
    real(r8) :: liqcldf(nlev,ncol)    ! liquid cloud fraction
    real(r8) :: t(nlev,ncol)          ! input temperature (K)
    real(r8) :: qn(nlev,ncol)         ! input water vapor mixing ratio (kg/kg)
    ! note: all input cloud variables are grid-averaged
    real(r8) :: qc(nlev,ncol)         ! cloud water mixing ratio (kg/kg)
    real(r8) :: qi(nlev,ncol)         ! cloud ice mixing ratio (kg/kg)
    real(r8) :: nc(nlev,ncol)         ! cloud water number conc (1/kg)
    real(r8) :: ni(nlev,ncol)         ! cloud ice number conc (1/kg)
    real(r8) :: pmid(nlev,ncol)       ! pressure at layer midpoints (pa)
    real(r8) :: pdel(nlev,ncol)       ! pressure difference across level (pa)
    real(r8) :: pint(nlevp,ncol)      ! air pressure layer interfaces (pa)
    real(r8) :: rpdel(nlev,ncol)      ! inverse pressure difference across level (pa)
    real(r8) :: zm(nlev,ncol)         ! geopotential height of model levels (m)
    real(r8) :: omega(nlev,ncol)      ! vertical velocity (Pa/s)

    real(r8) :: num_accum(nlev,ncol)  ! number m.r. of accumulation mode
    real(r8) :: num_aitken(nlev,ncol) ! number m.r. of aitken mode
    real(r8) :: num_coarse(nlev,ncol) ! number m.r. of coarse mode
    real(r8) :: coarse_dust(nlev,ncol) ! mass m.r. of coarse dust
    real(r8) :: coarse_nacl(nlev,ncol) ! mass m.r. of coarse nacl

    real(r8) :: kvh(nlevp,ncol)       ! vertical eddy diff coef (m2 s-1)
    real(r8) :: tke(nlevp,ncol)       ! TKE from the UW PBL scheme (m2 s-2)
 
    real(r8) :: cldn(nlev,ncol)       ! cloud fraction
    real(r8) :: cldo(nlev,ncol)       ! old cloud fraction

    real(r8) :: aer_mmr(nlev,ncol) ! aerosol mass mixing ratio

    real(r8) :: rho(nlev,ncol)        ! air density (kg m-3)
    real(r8) :: relhum(nlev,ncol)     ! relative humidity
    real(r8) :: icldm(nlev,ncol)      ! ice cloud fraction
    real(r8) :: lcldm(nlev,ncol)      ! liq cloud fraction
    real(r8) :: nfice(nlev,ncol)      ! fice variable
    real(r8) :: dumfice               ! dummy var in fice calc
    real(r8) :: lcldn(nlev,ncol)      ! fractional coverage of new liquid cloud
    real(r8) :: lcldo(nlev,ncol)      ! fractional coverage of old liquid cloud
    real(r8) :: qcld                  ! total cloud water
    real(r8) :: nctend_mixnuc(nlev,ncol)
 
    real(r8) :: dum, dum2             ! temporary dummy variable  
    real(r8) :: dmc, ssmc             ! variables for modal scheme.

    real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
    real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
    real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
    real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)

    real(r8) :: qs(ncol)              ! liquid-ice weighted sat mixing rat (kg/kg)
    real(r8) :: es(ncol)              ! liquid-ice weighted sat vapor press (pa)
    real(r8) :: gammas(ncol)          ! parameter for cond/evap of cloud water

    ! bulk aerosol variables
    real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
    real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

    real(r8) :: wsub(nlev,ncol)       ! diagnosed sub-grid vertical velocity st. dev. (m/s)
    real(r8) :: wsubi(nlev,ncol)      ! diagnosed sub-grid vertical velocity ice (m/s)

    ! history output for ice nucleation
    real(r8) :: nihf(nlev,ncol)       ! output number conc of ice nuclei due to heterogenous freezing (1/m3)
    real(r8) :: niimm(nlev,ncol)      ! output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
    real(r8) :: nidep(nlev,ncol)      ! output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
    real(r8) :: nimey(nlev,ncol)      ! output number conc of ice nuclei due to meyers deposition (1/m3)

    real(r8) :: wght


    t(:,:ncol)     = pstate%temp_at_pc_full_level%f(:,:ncol)
    qn(:,:ncol)    = pstate%tracer_mxrt_at_pc_full_level%f(1,:,:ncol)
    qc(:,:ncol)    = pstate%tracer_mxrt_at_pc_full_level%f(ixcldliq,:,:ncol)
    qi(:,:ncol)    = pstate%tracer_mxrt_at_pc_full_level%f(ixcldice,:,:ncol)
    nc(:,:ncol)    = pstate%tracer_mxrt_at_pc_full_level%f(ixnumliq,:,:ncol)
    ni(:,:ncol)    = pstate%tracer_mxrt_at_pc_full_level%f(ixnumice,:,:ncol)
    pmid(:,:ncol)  = pstate%pressure_at_pc_full_level%f(:,:ncol)
    pdel(:,:ncol)  = pstate%delp_at_pc_full_level%f(:,:ncol)
    pint(:,:ncol)  = pstate%pressure_at_pc_face_level%f(:,:ncol)
    rpdel(:,:ncol) = 1._r8/pstate%delp_at_pc_full_level%f(:,:ncol)
    zm(:,:ncol)    = pstate%z_at_pc_full_level%f(:,:ncol)
    omega(:,:ncol) = pstate%omega_at_pc_full_level%f(:,:ncol)

    liqcldf(:nlev,:ncol) = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:nlev,:ncol)
    icecldf(:nlev,:ncol) = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:nlev,:ncol)


    if (clim_modal_aero) then
        cldn(:nlev,:ncol) = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,:nlev,:ncol)
        cldo(:nlev,:ncol) = pstate_cam%microp_cldo_at_pc_full_level%f(old_time_level,:nlev,:ncol)
    end if

    ! initialize output
    pstate_cam%microp_areo_naai%f(1:nlev,1:ncol)     = 0._r8  
    pstate_cam%microp_areo_naai_hom%f(1:nlev,1:ncol) = 0._r8  
    pstate_cam%microp_areo_npccn%f(1:nlev,1:ncol)    = 0._r8  
    pstate_cam%microp_areo_nacon%f(:,1:nlev,1:ncol)  = 0._r8

    ! set default or fixed dust bins for contact freezing
    pstate_cam%microp_areo_rndst%f(1,1:nlev,1:ncol)  = rn_dst1
    pstate_cam%microp_areo_rndst%f(2,1:nlev,1:ncol)  = rn_dst2
    pstate_cam%microp_areo_rndst%f(3,1:nlev,1:ncol)  = rn_dst3
    pstate_cam%microp_areo_rndst%f(4,1:nlev,1:ncol)  = rn_dst4

    ! initialize history output fields for ice nucleation
    nihf(1:nlev,1:ncol)  = 0._r8  
    niimm(1:nlev,1:ncol) = 0._r8  
    nidep(1:nlev,1:ncol) = 0._r8 
    nimey(1:nlev,1:ncol) = 0._r8 

    ! initialize time-varying parameters
    do k = top_lev, nlev
       do i = 1, ncol
          rho(k,i) = pmid(k,i)/(rdry*t(k,i))
       end do
    end do

    if (clim_modal_aero) then
      ! mode number mixing ratios
      call rad_cnst_get_mode_num(0, mode_accum_idx,  'a', num_accum)
      call rad_cnst_get_mode_num(0, mode_aitken_idx, 'a', num_aitken)
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', num_coarse)

      ! mode specie mass m.r.
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', coarse_dust)
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', coarse_nacl)

    else
       ! init number/mass arrays for bulk aerosols
       !-----------LiXH 2020-10-21-------------
       if(naer_all .gt. 0)then
       allocate(naer2(naer_all,nlev,ncol), &
                maerosol(naer_all,nlev,ncol))
       end if
       !-----------LiXH 2020-10-21-------------
 
       do m = 1, naer_all
          call rad_cnst_get_aer_mmr(0, m, aer_mmr)
          maerosol(m,:,:ncol) = aer_mmr(:,:ncol)*rho(:,:ncol)

          if (m .eq. idxsul) then
              naer2(m,:,:ncol) = maerosol(m,:,:ncol)*num_to_mass_aer(m)*bulk_scale
          else
              naer2(m,:,:ncol) = maerosol(m,:,:ncol)*num_to_mass_aer(m)
          end if
       end do
    end if

    ! More refined computation of sub-grid vertical velocity
    ! Set to be zero at the surface by initialization.

    select case (trim(eddy_scheme))
    case ('diag_TKE')
       tke(1:nlevp,1:ncol) = pstate_cam%pbl_tke_at_pc_face_level%f(1:nlevp,1:ncol)
    case ('CLUBB_SGS')
        print*,'In microp_aero_driver'
        call endrun('ERROR: CLUBB is not avaiable!!') 
!--------------LiXH has not completed CLUBB_SGS----------------->
!       itim = pbuf_old_tim_idx()
!       call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim/),kount=(/pcols,nlevp,1/))
!       allocate(tke(pcols,nlevp))
!       tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)
!<-------------LiXH has not completed CLUBB_SGS------------------
    case default
       kvh(1:nlevp,1:ncol) = pstate_cam%pbl_kvh_at_pc_face_level%f(1:nlevp,1:ncol) 
    end select

    ! Set minimum values above top_lev.
    wsub(:top_lev-1,:ncol)  = 0.20_r8
    wsubi(:top_lev-1,:ncol) = 0.001_r8

    do k = top_lev, nlev
       do i = 1, ncol

          select case (trim(eddy_scheme))
          case ('diag_TKE', 'CLUBB_SGS')
                wsub(k,i) = sqrt(0.5_r8*(tke(k,i) + tke(k+1,i))*(2._r8/3._r8))
                wsub(k,i) = min(wsub(k,i),10._r8)
          case default 
             ! get sub-grid vertical velocity from diff coef.
             ! following morrison et al. 2005, JAS
             ! assume mixing length of 30 m
                dum = (kvh(k,i) + kvh(k+1,i))/2._r8/30._r8
             ! use maximum sub-grid vertical vel of 10 m/s
                dum = min(dum, 10._r8)
             ! set wsub to value at current vertical level
                wsub(k,i)  = dum
          end select

          wsubi(k,i) = max(0.001_r8, wsub(k,i))
          wsubi(k,i) = min(wsubi(k,i), 0.2_r8)  !wsubimax: 0-1, default: 0.2 

!-------------LiXH has not completed CLUBB_SGS------------->
!          if (wsubi(k,i) .le. 0.04_r8) then
!             nucboast=100._r8
!             wsubi(k,i)=nucboast*wsubi(k,i)  ! boost ice SGS vertical velocity in CAM-CLUBB
!             ! to force nucleation in upper-level stratiform 
!             ! clouds.  Temporary fix until cloud-top radiative
!             ! cooling parameterization is added to CLUBB similar
!             ! to the one of appendix C of Bretherton and Park (2009).  
!          endif
!<------------LiXH has not completed CLUBB_SGS--------------
          wsub(k,i)  = max(0.20_r8, wsub(k,i))  !wsubmin: 0-1, default:0.2
       end do
    end do

    !Get humidity and saturation vapor pressures
 
    ! find wet bulk temperature and saturation value for provisional t and q without
    ! condensation
 
    do k = top_lev, nlev
    call qsat_water(t(k,:ncol), pmid(k,:ncol), es(:ncol), qs(:ncol), gam=gammas(:ncol))
 
       do i = 1, ncol
          relhum(k,i) = qn(k,i)/qs(i)
 
          ! get cloud fraction, check for minimum
          icldm(k,i) = max(icecldf(k,i), mincld)
          lcldm(k,i) = max(liqcldf(k,i), mincld)
 
          ! calculate nfice based on liquid and ice mmr (no rain and snow mmr available yet)
          nfice(k,i) = 0._r8
          dumfice    = qc(k,i) + qi(k,i)
          if (dumfice > qsmall .and. qi(k,i) > qsmall) then
             nfice(k,i) = qi(k,i)/dumfice
          end if
       end do
    end do

    !ICE Nucleation

    do k = top_lev, nlev
       do i = 1, ncol

          if (t(k,i).lt.tmelt - 5._r8) then

             ! compute aerosol number for so4, soot, and dust with units #/cm^3
             so4_num  = 0._r8
             soot_num = 0._r8
             dst1_num = 0._r8
             dst2_num = 0._r8
             dst3_num = 0._r8
             dst4_num = 0._r8
             dst_num  = 0._r8
             if (clim_modal_aero) then
               !For modal aerosols, assume for the upper troposphere:
               ! soot = accumulation mode
               ! sulfate = aiken mode
               ! dust = coarse mode
               ! since modal has internal mixtures.
               soot_num = num_accum(k,i)*rho(k,i)*1.0e-6_r8
               dmc  = coarse_dust(k,i)*rho(k,i)
               ssmc = coarse_nacl(k,i)*rho(k,i)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting 
                  wght = 1._r8
               else
                  ! 3-mode -- needs weighting for dust since dust and seasalt are combined in the "coarse" mode type
                  wght = dmc/(ssmc + dmc)
               endif

               if (dmc > 0._r8) then
                  dst_num = wght * num_coarse(k,i)*rho(k,i)*1.0e-6_r8
               else 
                  dst_num = 0.0_r8
               end if

               if (pstate_cam%aerosol_dgnum%f(mode_aitken_idx,k,i) > 0._r8) then
                  ! only allow so4 with D>0.1 um in ice nucleation
                  so4_num  = num_aitken(k,i)*rho(k,i)*1.0e-6_r8 &
                     * (0.5_r8 - 0.5_r8*erf(log(0.1e-6_r8/pstate_cam%aerosol_dgnum%f(mode_aitken_idx,k,i))/  &
                     (2._r8**0.5_r8*log(sigmag_aitken))))
               else 
                  so4_num = 0.0_r8 
               end if
               so4_num = max(0.0_r8, so4_num)

             else
 
                if (idxsul > 0) then 
                   so4_num = naer2(idxsul,k,i)/25._r8 *1.0e-6_r8
                end if
                if (idxbcphi > 0) then 
                   soot_num = naer2(idxbcphi,k,i)/25._r8 *1.0e-6_r8
                end if
                if (idxdst1 > 0) then 
                   dst1_num = naer2(idxdst1,k,i)/25._r8 *1.0e-6_r8
                end if
                if (idxdst2 > 0) then 
                   dst2_num = naer2(idxdst2,k,i)/25._r8 *1.0e-6_r8
                end if
                if (idxdst3 > 0) then 
                   dst3_num = naer2(idxdst3,k,i)/25._r8 *1.0e-6_r8
                end if
                if (idxdst4 > 0) then 
                   dst4_num = naer2(idxdst4,k,i)/25._r8 *1.0e-6_r8
                end if
                dst_num = dst1_num + dst2_num + dst3_num + dst4_num
 
             end if
 
             ! *** Turn off soot nucleation ***
             soot_num = 0.0_r8
 
             call nucleati( wsubi(k,i), t(k,i), relhum(k,i), icldm(k,i), qc(k,i), &
                            nfice(k,i), rho(k,i), so4_num, dst_num, soot_num,     &
                            pstate_cam%microp_areo_naai%f(k,i),                       &  
                            nihf(k,i), niimm(k,i), nidep(k,i), nimey(k,i))
 
             pstate_cam%microp_areo_naai_hom%f(k,i) = nihf(k,i)
 
             ! output activated ice (convert from #/kg -> #/m3)
             nihf(k,i)     = nihf(k,i) *rho(k,i)
             niimm(k,i)    = niimm(k,i)*rho(k,i)
             nidep(k,i)    = nidep(k,i)*rho(k,i)
             nimey(k,i)    = nimey(k,i)*rho(k,i)

          end if
       end do
    end do

    if (clim_modal_aero) then
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      !droplet activation for modal aerosol

      ! partition cloud fraction into liquid water part
      lcldn = 0._r8
      lcldo = 0._r8
      do k = top_lev, nlev
         do i = 1, ncol
            qcld = qc(k,i) + qi(k,i)
            if (qcld > qsmall) then
               lcldn(k,i) = cldn(k,i)*qc(k,i)/qcld
               lcldo(k,i) = cldo(k,i)*qc(k,i)/qcld
            end if
         end do
      end do

      call dropmixnuc(  deltatin, ncol, ixnumliq,     &
                        ptend_microp_aero%tend_q%f,  &  
                        wsub,  lcldn, lcldo, nctend_mixnuc)
 
      pstate_cam%microp_areo_npccn%f(:,:ncol) = nctend_mixnuc(:,:ncol)
    else

       !droplet activation for bulk aerosol
    
       do k = top_lev, nlev
          do i = 1, ncol

             if (qc(k,i) >= qsmall) then

                ! get droplet activation rate

                call ndrop_bam_run( &
                   wsub(k,i), t(k,i), rho(k,i), naer2(:,k,i), naer_all, &
                   naer_all, maerosol(:,k,i),  &
                   dum2)
                dum = dum2
             else
                dum = 0._r8
             end if

             ! note: deltatin/2.  accounts for sub step in microphysics
             ! ***** This assumes two sub-steps in microphysics.  It's dangerous to 
             ! ***** make that assumption here.  Should move all coding related to 
             ! ***** microphysics substepping into the microphysics.
             pstate_cam%microp_areo_npccn%f(k,i) = (dum - nc(k,i)/lcldm(k,i))/(deltatin/2._r8)*lcldm(k,i)

          end do
       end do

    end if

    ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
    ! estimate rndst and nanco for 4 dust bins here to pass to MG microphysics

    do k = top_lev, nlev
       do i = 1, ncol

          if (t(k,i) < 269.15_r8) then

             if (clim_modal_aero) then

               ! For modal aerosols:
               !  use size '3' for dust coarse mode...
               !  scale by dust fraction in coarse mode
               
               dmc  = coarse_dust(k,i)
               ssmc = coarse_nacl(k,i)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting 
                  wght = 1._r8
               else
                  ! 3-mode -- needs weighting for dust since dust and seasalt are combined in the "coarse" mode type
                  wght = dmc/(ssmc + dmc)
               endif

               if (dmc > 0.0_r8) then
                  pstate_cam%microp_areo_nacon%f(3,k,i) = wght*num_coarse(k,i)*rho(k,i)
               else
                  pstate_cam%microp_areo_nacon%f(3,k,i) = 0._r8
               end if

               !also redefine parameters based on size...

               pstate_cam%microp_areo_rndst%f(3,k,i) = 0.5_r8*pstate_cam%aerosol_dgnumwet%f(mode_coarse_dst_idx,k,i)
               if (pstate_cam%microp_areo_rndst%f(3,k,i) <= 0._r8) then 
                  pstate_cam%microp_areo_rndst%f(3,k,i) = rn_dst3
               end if

             else

                !For Bulk Aerosols: set equal to aerosol number for dust for bins 2-4 (bin 1=0)

                if (idxdst2 > 0) then 
                   pstate_cam%microp_areo_nacon%f(2,k,i) = naer2(idxdst2,k,i)
                end if
                if (idxdst3 > 0) then 
                   pstate_cam%microp_areo_nacon%f(3,k,i) = naer2(idxdst3,k,i)
                end if
                if (idxdst4 > 0) then 
                   pstate_cam%microp_areo_nacon%f(4,k,i) = naer2(idxdst4,k,i)
                end if
             end if

          end if
       end do
    end do

    !bulk aerosol ccn concentration (modal does it in ndrop, from dropmixnuc)


    if (.not. clim_modal_aero) then

       ! ccn concentration as diagnostic
       call ndrop_bam_ccn(ncol, maerosol, naer2)

       !-----------LiXH 2020-10-21-------------
       if(naer_all .gt. 0)then
       deallocate( naer2,    &
                   maerosol)
       end if 
       !-----------LiXH 2020-10-21-------------

    end if

! output
! wsub, wsubi
! nihf, niimm, nidep, nimey
!---------------LiXH Test----------------
!if(mpi_rank()==0)print*,'LiXH Test microp_aero: rndst:',maxval(pstate_cam%microp_areo_rndst%f),minval(pstate_cam%microp_areo_rndst%f)
!if(mpi_rank()==0)print*,'LiXH Test microp_aero: naai:',maxval(pstate_cam%microp_areo_naai%f),minval(pstate_cam%microp_areo_naai%f)
!if(mpi_rank()==0)print*,'LiXH Test microp_aero: num_coarse:',maxval(num_coarse),minval(num_coarse)
!if(mpi_rank()==0)print*,'LiXH Test microp_aero: nacon:',maxval(pstate_cam%microp_areo_nacon%f(3,:,:)),minval(pstate_cam%microp_areo_nacon%f(3,:,:))
!if(mpi_rank()==0)print*,'LiXH Test microp_aero: npccn:',maxval(pstate_cam%microp_areo_npccn%f),minval(pstate_cam%microp_areo_npccn%f)
!if(mpi_rank()==0)print*,'LiXH Test microp_aero: cld:',maxval(cldn)
!---------------LiXH Test----------------

    end subroutine microp_aero_driver

 end module grist_microp_aero
