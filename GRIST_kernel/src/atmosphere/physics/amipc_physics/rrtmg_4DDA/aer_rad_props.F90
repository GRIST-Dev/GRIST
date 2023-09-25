!========================================================================================
!
!  Created by LiXiaohan on 19/8/13, adopted from CAM5.
!
!  Converts aerosol masses to bulk optical properties for sw and lw radiation
!  computations.  
!========================================================================================

 module aer_rad_props

    use grist_constants,                    only: r8, rga=>gravi
    use grist_nml_module,                   only: nlev, nlevp
    use grist_handle_error,                 only: endrun
    use grist_physics_data_structure,       only: pstate
    use radconstants,                       only: nrh, nswbands, nlwbands, &
                                                  idx_sw_diag, ot_length
    use grist_rad_constituents,             only: rad_cnst_get_info,       &
                                                  rad_cnst_get_aer_mmr,    &
                                                  rad_cnst_get_aer_props
    use grist_wv_saturation,                only: qsat
    use modal_aer_opt,                      only: modal_aero_sw, modal_aero_lw
    use cloud_fraction,                     only: clim_modal_aero_top_lev
    use grist_mpi

 implicit none
 private
 save

 integer :: top_lev = 1
 
 public ::  aer_rad_props_init,        &
            aer_rad_props_sw,          & ! return SW optical props of aerosols
            aer_rad_props_lw             ! return LW optical props of aerosols

 ! Private data
 character(len=256), allocatable :: odv_names(:)  ! outfld names for visible OD

 !-------------------------------------->
 !volcanic_radius should be initiallized and calculated in mozart/mo_sad.F90
 !Should be modified, LiXH
 real(r8), allocatable :: geometric_radius(:,:)
 !<--------------------------------------


 contains

    subroutine aer_rad_props_init(ncol)
    use phys_control, only: phys_getopts
    !io
    integer,    intent(in)  :: ncol
    !local
    integer                    :: i
    integer                    :: numaerosols
    character(len=64), allocatable :: aernames(:)
    logical                    :: prog_modal_aero      ! Prognostic modal aerosols present

    !-------------------------------------->
    !volcanic_radius should be initiallized and calculated in mozart/mo_sad.F90
    !Should be modified, LiXH
    allocate(geometric_radius(nlev,ncol))
    geometric_radius = 1._r8
    !<--------------------------------------


    call phys_getopts( prog_modal_aero_out        = prog_modal_aero )

    ! Limit modal aerosols with top_lev here.
    if (prog_modal_aero) top_lev = clim_modal_aero_top_lev

    ! Contributions to AEROD_v from individual aerosols (climate species).

    ! number of bulk aerosols in climate list
    call rad_cnst_get_info(0, naero=numaerosols)

    ! get names of bulk aerosols
    if(numaerosols .gt. 0)then
        allocate(aernames(numaerosols))
    else
        if(mpi_rank()==0)print*,'no specific bulk aerosols for RRTMG'
        return
    end if

    call rad_cnst_get_info(0, aernames=aernames)

    ! diagnostic output for bulk aerosols
    ! create outfld names for visible OD
    allocate(odv_names(numaerosols))
    do i = 1, numaerosols
       odv_names(i) = 'ODV_'//trim(aernames(i))
    end do

! output:
!    ! Determine default fields
!    if (history_amwg ) then 
!       call add_default ('AEROD_v', 1, ' ')
!    endif   
!    
!    if ( history_aero_optics ) then 
!       call add_default ('AEROD_v', 1, ' ')
!       do i = 1, numaerosols
!          odv_names(i) = 'ODV_'//trim(aernames(i))
!          call add_default (odv_names(i), 1, ' ')
!       end do
!     endif

    if(allocated(aernames))deallocate(aernames)

    end subroutine aer_rad_props_init


!Purpose: Return bulk layer tau, omega, g, f for all spectral intervals.
    subroutine aer_rad_props_sw(ncol, list_idx, nnite, idxnite, &
                                tau, tau_w, tau_w_g, tau_w_f)
    ! io
    integer,             intent(in) :: ncol
    integer,             intent(in) :: list_idx             ! index of the climate or a diagnostic list
    integer,             intent(in) :: nnite                ! number of night columns
    integer,             intent(in) :: idxnite(:)           ! local column indices of night columns
    real(r8), intent(out) :: tau    (nswbands,0:nlev,ncol)  ! aerosol extinction optical depth
    real(r8), intent(out) :: tau_w  (nswbands,0:nlev,ncol)  ! aerosol single scattering albedo * tau
    real(r8), intent(out) :: tau_w_g(nswbands,0:nlev,ncol)  ! aerosol assymetry parameter * tau * w
    real(r8), intent(out) :: tau_w_f(nswbands,0:nlev,ncol)  ! aerosol forward scattered fraction * tau * w
    ! local
    integer :: k, i    ! lev and daycolumn indices
    integer :: iswband ! sw band indices

    ! optical props for each aerosol
    ! hygroscopic
    real(r8) :: h_ext(nrh,nswbands)
    real(r8) :: h_ssa(nrh,nswbands)
    real(r8) :: h_asm(nrh,nswbands)
    ! non-hygroscopic
    real(r8) :: n_ext(nswbands)
    real(r8) :: n_ssa(nswbands)
    real(r8) :: n_asm(nswbands)
    real(r8) :: n_scat(nswbands)
    real(r8) :: n_ascat(nswbands)
    ! radius-dependent
    real(r8), allocatable :: r_ext(:,:)    ! radius-dependent mass-specific extinction
    real(r8), allocatable :: r_scat(:,:)
    real(r8), allocatable :: r_ascat(:,:)
    real(r8), allocatable :: r_mu(:)       ! log(radius) domain variable for r_ext, r_scat, r_ascat

    ! radiative properties for each aerosol
    real(r8) :: ta (nswbands,nlev,ncol)
    real(r8) :: tw (nswbands,nlev,ncol)
    real(r8) :: twf(nswbands,nlev,ncol)
    real(r8) :: twg(nswbands,nlev,ncol)

    ! aerosol masses
    real(r8) :: aermmr(nlev,ncol)      ! mass mixing ratio of aerosols
    real(r8) :: mmr_to_mass(nlev,ncol) ! conversion factor for mmr to mass
    real(r8) :: aermass(nlev,ncol)     ! mass of aerosols

    ! for table lookup into rh grid
    real(r8) :: es(nlev,ncol)     ! saturation vapor pressure
    real(r8) :: qs(nlev,ncol)     ! saturation specific humidity
    real(r8) :: rh(nlev,ncol)
    real(r8) :: rhtrunc(nlev,ncol)
    real(r8) :: wrh(nlev,ncol)
    integer  :: krh(nlev,ncol)
 
    integer  :: numaerosols     ! number of bulk aerosols in climate/diagnostic list
    integer  :: nmodes          ! number of aerosol modes in climate/diagnostic list
    integer  :: iaerosol        ! index into bulk aerosol list

    character(len=ot_length) :: opticstype       ! hygro or nonhygro
 
    integer :: n_mu_samples   !LiXH adds for GRIST


    !-----------LiXH Modified dry pressure to pressure---------->
    !do k = 1, nlev
    !   mmr_to_mass(k,:ncol) = rga * state%pdeldry(k,:ncol)
    !enddo
    do k = 1, nlev
       mmr_to_mass(k,:ncol) = rga * pstate%delp_at_pc_full_level%f(k,:ncol)
    enddo
    !<----------LiXH Modified dry pressure to pressure-----------

    ! top layer (ilev = 0) has no aerosol (ie tau = 0)
    ! also initialize rest of layers to accumulate od's
    tau    (:,:,1:ncol) = 0._r8
    tau_w  (:,:,1:ncol) = 0._r8
    tau_w_g(:,:,1:ncol) = 0._r8
    tau_w_f(:,:,1:ncol) = 0._r8

    ! calculate relative humidity for table lookup into rh grid
    call qsat(pstate%temp_at_pc_full_level%f(1:nlev,1:ncol),     &
              pstate%pressure_at_pc_full_level%f(1:nlev,1:ncol), &
              es(1:nlev,1:ncol), qs(1:nlev,1:ncol))
    rh(1:nlev,1:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(1,1:nlev,1:ncol)   &
                      / qs(1:nlev,1:ncol)

    rhtrunc(1:nlev,1:ncol) = min(rh(1:nlev,1:ncol),1._r8)
    krh(1:nlev,1:ncol) = min(floor( rhtrunc(1:nlev,1:ncol) * nrh ) + 1, nrh - 1) ! index into rh mesh
    wrh(1:nlev,1:ncol) = rhtrunc(1:nlev,1:ncol) * nrh - krh(1:nlev,1:ncol)       ! (-) weighting on left side values

    ! get number of bulk aerosols and number of modes in current list
    call rad_cnst_get_info(list_idx, naero=numaerosols, nmodes=nmodes)

!-------------LiXH Test--------------
! LiXH remove aerosol effect
!numaerosols = 0; nmodes = 0
!-------------LiXH Test--------------

    ! Contributions from modal aerosols.
    if (nmodes > 0) then
       call modal_aero_sw(ncol, list_idx, nnite, idxnite, &
                          tau, tau_w, tau_w_g, tau_w_f)
    else
       tau    (:,:,1:ncol) = 0._r8
       tau_w  (:,:,1:ncol) = 0._r8
       tau_w_g(:,:,1:ncol) = 0._r8
       tau_w_f(:,:,1:ncol) = 0._r8
    end if

    ! Contributions from bulk aerosols.
    do iaerosol = 1, numaerosols

       ! get bulk aerosol mass mixing ratio
       call rad_cnst_get_aer_mmr(list_idx, iaerosol, aermmr)
       aermass(1:top_lev-1,1:ncol) = 0._r8
       aermass(top_lev:nlev,1:ncol) = aermmr(top_lev:nlev,1:ncol) * mmr_to_mass(top_lev:nlev,1:ncol)

       ! get optics type
       call rad_cnst_get_aer_props(list_idx, iaerosol, opticstype=opticstype)

       select case (trim(opticstype))
       case('hygro','hygroscopic','hygroscopi')
          ! get optical properties for hygroscopic aerosols
          call rad_cnst_get_aer_props(list_idx, iaerosol, sw_hygro_ext=h_ext, sw_hygro_ssa=h_ssa, sw_hygro_asm=h_asm)
          call get_hygro_rad_props(ncol, krh, wrh, aermass, h_ext, h_ssa, h_asm, ta, tw, twg, twf)
          tau    (:,1:nlev,1:ncol) = tau    (:,1:nlev,1:ncol) + ta (:,:,1:ncol)
          tau_w  (:,1:nlev,1:ncol) = tau_w  (:,1:nlev,1:ncol) + tw (:,:,1:ncol)
          tau_w_g(:,1:nlev,1:ncol) = tau_w_g(:,1:nlev,1:ncol) + twg(:,:,1:ncol)
          tau_w_f(:,1:nlev,1:ncol) = tau_w_f(:,1:nlev,1:ncol) + twf(:,:,1:ncol)
 
       case('nonhygro','insoluble ')
          ! get optical properties for non-hygroscopic aerosols
          call rad_cnst_get_aer_props(list_idx, iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_ssa=n_ssa, sw_nonhygro_asm=n_asm)
 
          call get_nonhygro_rad_props(ncol, aermass, n_ext, n_ssa, n_asm, ta, tw, twg, twf)
          tau    (:,1:nlev,1:ncol) = tau    (:,1:nlev,1:ncol) + ta (:,:,1:ncol)
          tau_w  (:,1:nlev,1:ncol) = tau_w  (:,1:nlev,1:ncol) + tw (:,:,1:ncol)
          tau_w_g(:,1:nlev,1:ncol) = tau_w_g(:,1:nlev,1:ncol) + twg(:,:,1:ncol)
          tau_w_f(:,1:nlev,1:ncol) = tau_w_f(:,1:nlev,1:ncol) + twf(:,:,1:ncol)
 
       case('volcanic')
          ! get optical properties for volcanic aerosols
          call rad_cnst_get_aer_props(list_idx, iaerosol, sw_nonhygro_ext=n_ext, sw_nonhygro_scat=n_scat, sw_nonhygro_ascat=n_ascat)
 
          call get_volcanic_rad_props(ncol, aermass, n_ext, n_scat, n_ascat, ta, tw, twg, twf)
          tau    (:,1:nlev,1:ncol) = tau    (:,1:nlev,1:ncol) + ta (:,:,1:ncol)
          tau_w  (:,1:nlev,1:ncol) = tau_w  (:,1:nlev,1:ncol) + tw (:,:,1:ncol)
          tau_w_g(:,1:nlev,1:ncol) = tau_w_g(:,1:nlev,1:ncol) + twg(:,:,1:ncol)
          tau_w_f(:,1:nlev,1:ncol) = tau_w_f(:,1:nlev,1:ncol) + twf(:,:,1:ncol)
 
       case('volcanic_radius')
          ! get optical properties for volcanic aerosols
          call rad_cnst_get_aer_props(list_idx, iaerosol, n_mu_samples=n_mu_samples)
          allocate( r_ext(nswbands,n_mu_samples) )
          allocate( r_scat(nswbands,n_mu_samples) ) 
          allocate( r_ascat(nswbands,n_mu_samples) )
          allocate( r_mu(n_mu_samples) )

          call rad_cnst_get_aer_props(list_idx, iaerosol, r_sw_ext=r_ext, r_sw_scat=r_scat, r_sw_ascat=r_ascat, mu=r_mu)
          call get_volcanic_radius_rad_props(ncol, n_mu_samples, aermass, r_ext, r_scat, r_ascat, r_mu, ta, tw, twg, twf)
          tau    (:,1:nlev,1:ncol) = tau    (:,1:nlev,1:ncol) + ta (:,:,1:ncol)
          tau_w  (:,1:nlev,1:ncol) = tau_w  (:,1:nlev,1:ncol) + tw (:,:,1:ncol)
          tau_w_g(:,1:nlev,1:ncol) = tau_w_g(:,1:nlev,1:ncol) + twg(:,:,1:ncol)
          tau_w_f(:,1:nlev,1:ncol) = tau_w_f(:,1:nlev,1:ncol) + twf(:,:,1:ncol)
 
          deallocate( r_ext )
          deallocate( r_scat ) 
          deallocate( r_ascat )
          deallocate( r_mu )

       case('zero')
          ! no effect of "zero" aerosols, so update nothing
       case default
          if(mpi_rank()==0)print*,'aer_rad_props_sw: unsupported opticstype:'//trim(opticstype)//':' 
          call endrun('In aer_rad_props_sw')
       end select
       ! diagnostic output of individual aerosol optical properties
       ! currently implemented for climate list only
       !call aer_vis_diag_out(ncol, nnite, idxnite, iaerosol, ta(idx_sw_diag,:,:), list_idx)

    end do

    ! diagnostic output of total aerosol optical properties
    ! currently implemented for climate list only
    !call aer_vis_diag_out(ncol, nnite, idxnite, 0, tau(idx_sw_diag,:,:), list_idx)


    end subroutine aer_rad_props_sw


! Purpose: Compute aerosol transmissions needed in absorptivity/emissivity calculations
!          lw extinction is the same representation for all 
!          species.  If this changes, this routine will need to do something
!          similar to the sw with routines like get_hygro_lw_abs
    subroutine aer_rad_props_lw(ncol, list_idx, odap_aer)

    use radconstants,  only: ot_length

    ! io
    integer,             intent(in)  :: ncol
    integer,             intent(in)  :: list_idx                      ! index of the climate or a diagnostic list
    real(r8),            intent(out) :: odap_aer(nlwbands,nlev,ncol)  ! [fraction] absorption optical depth, per layer

    ! Local variables

    integer :: bnd_idx     ! LW band index
    integer :: i           ! column index
    integer :: k           ! lev index
    integer :: numaerosols ! number of bulk aerosols in climate/diagnostic list
    integer :: nmodes      ! number of aerosol modes in climate/diagnostic list
    integer :: iaerosol    ! index into bulk aerosol list
    character(len=ot_length) :: opticstype       ! hygro or nonhygro

    ! optical props for each aerosol
    real(r8) :: lw_abs(nlwbands)
    real(r8) :: lw_hygro_abs(nrh,nlwbands)
 
    ! volcanic lookup table
    real(r8), allocatable :: r_lw_abs(:,:)  ! radius dependent mass-specific absorption coefficient
    real(r8), allocatable :: r_mu(:)        ! log(geometric_mean_radius) domain samples of r_lw_abs(:,:)
    integer  :: idx                     ! index to pbuf for geometric radius
    real(r8) :: mu(nlev,ncol)           ! log(geometric_radius)
    real(r8) :: r_mu_min, r_mu_max, wmu, mutrunc
    integer  :: nmu, kmu

    ! for table lookup into rh grid
    real(r8) :: es(nlev,ncol)     ! saturation vapor pressure
    real(r8) :: qs(nlev,ncol)     ! saturation specific humidity
    real(r8) :: rh(nlev,ncol)
    real(r8) :: rhtrunc(nlev,ncol)
    real(r8) :: wrh(nlev,ncol)
    integer  :: krh(nlev,ncol)

    ! aerosol (vertical) mass path and extinction
    ! aerosol masses
    real(r8) :: aermmr(nlev,ncol)    ! mass mixing ratio of aerosols
    real(r8) :: mmr_to_mass(nlev,ncol) ! conversion factor for mmr to mass
    real(r8) :: aermass(nlev,ncol)     ! mass of aerosols

    integer :: n_mu_samples   !LiXH adds for GRIST
    !-----------------------------------------------------------------------------

    ! get number of bulk aerosols and number of modes in current list
    call rad_cnst_get_info(list_idx, naero=numaerosols, nmodes=nmodes)


!----------------LiXH Test---------------
! Lixh remove the effect of aerosol
!nmodes = 0; numaerosols = 0
!----------------LiXH Test---------------

    ! Contributions from modal aerosols.

    if (nmodes > 0) then
       call modal_aero_lw(ncol, list_idx, odap_aer)
    else
       odap_aer = 0._r8
    end if

    ! Contributions from bulk aerosols.
    if (numaerosols > 0) then

       ! compute mixing ratio to mass conversion
       !-----------LiXH Modified dry pressure to pressure---------->
       !do k = 1, nlev
       !   mmr_to_mass(k,:ncol) = rga * state%pdeldry(k,:ncol)
       !enddo
       do k = 1, nlev
          mmr_to_mass(k,:ncol) = rga * pstate%delp_at_pc_full_level%f(k,:ncol)
       enddo
       !<----------LiXH Modified dry pressure to pressure-----------


       ! calculate relative humidity for table lookup into rh grid
       call qsat(pstate%temp_at_pc_full_level%f(1:nlev,1:ncol),     &
                 pstate%pressure_at_pc_full_level%f(1:nlev,1:ncol), &
                 es(1:nlev,1:ncol), qs(1:nlev,1:ncol))
       rh(1:nlev,1:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(1,1:nlev,1:ncol)   &
                      / qs(1:nlev,1:ncol)

       rhtrunc(1:nlev,1:ncol) = min(rh(1:nlev,1:ncol),1._r8)
       krh(1:nlev,1:ncol) = min(floor( rhtrunc(1:nlev,1:ncol) * nrh ) + 1, nrh - 1) ! index into rh mesh
       wrh(1:nlev,1:ncol) = rhtrunc(1:nlev,1:ncol) * nrh - krh(1:nlev,1:ncol)       ! (-) weighting on left side values

    end if

    ! Loop over bulk aerosols in list.
    do iaerosol = 1, numaerosols

       ! get aerosol mass mixing ratio
       call rad_cnst_get_aer_mmr(list_idx, iaerosol, aermmr)
       aermass(1:top_lev-1,1:ncol) = 0._r8
       aermass(top_lev:nlev,1:ncol) = aermmr(top_lev:nlev,1:ncol) * mmr_to_mass(top_lev:nlev,1:ncol)

       ! get optics type
       call rad_cnst_get_aer_props(list_idx, iaerosol, opticstype=opticstype)
       select case (trim(opticstype))
       case('hygroscopic')
           ! get optical properties for hygroscopic aerosols
          call rad_cnst_get_aer_props(list_idx, iaerosol, lw_hygro_ext=lw_hygro_abs)
          do bnd_idx = 1, nlwbands
             do k = 1, nlev
                do i = 1, ncol
                   odap_aer(bnd_idx,k,i) = odap_aer(bnd_idx,k,i) + &
                        aermass(k,i) * &
                        ((1 + wrh(k,i)) * lw_hygro_abs(krh(k,i)+1,bnd_idx) &
                        - wrh(k,i)  * lw_hygro_abs(krh(k,i),  bnd_idx))
                end do
             end do
          end do
       case('insoluble','nonhygro','hygro','volcanic')
           ! get optical properties for hygroscopic aerosols
          call rad_cnst_get_aer_props(list_idx, iaerosol, lw_ext=lw_abs)
          do bnd_idx = 1, nlwbands
             do k = 1, nlev          
                do i = 1, ncol
                   odap_aer(bnd_idx,k,i) = odap_aer(bnd_idx,k,i) + lw_abs(bnd_idx)*aermass(k,i)
                end do
             end do
          end do
          
       case('volcanic_radius')
          ! get optical properties for hygroscopic aerosols
          call rad_cnst_get_aer_props(list_idx, iaerosol, n_mu_samples=n_mu_samples)
          allocate(r_lw_abs(nlwbands,n_mu_samples))
          allocate(r_mu(n_mu_samples))

          call rad_cnst_get_aer_props(list_idx, iaerosol, r_lw_abs=r_lw_abs, mu=r_mu)
          ! get microphysical properties for volcanic aerosols
          !idx = pbuf_get_index('VOLC_RAD_GEOM')
          !call pbuf_get_field(pbuf, idx, geometric_radius )
          
          ! interpolate in radius
          ! caution: clip the table with no warning when outside bounds
          nmu = size(r_mu)
          r_mu_max = r_mu(nmu)
          r_mu_min = r_mu(1)
          do i = 1, ncol
             do k = 1, nlev
                if(geometric_radius(k,i) > 0._r8) then
                   mu(k,i) = log(geometric_radius(k,i))
                else
                   mu(k,i) = 0._r8
                endif
                mutrunc = max(min(mu(k,i),r_mu_max),r_mu_min)
                kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
                wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)
                do bnd_idx = 1, nlwbands
                   odap_aer(bnd_idx,k,i) = odap_aer(bnd_idx,k,i) + &
                      aermass(k,i) * &
                      ((1._r8 - wmu) * r_lw_abs(bnd_idx, kmu  ) + &
                      (wmu) * r_lw_abs(bnd_idx, kmu+1))
                end do
             end do
          end do
 
          deallocate(r_lw_abs)
          deallocate(r_mu)

       case('zero')
          ! zero aerosols types have no optical effect, so do nothing.
       case default
          if(mpi_rank()==0)print*,'aer_rad_props_lw: unsupported opticstype: '//trim(opticstype) 
          call endrun('In aer_rad_props_lw')
       end select
    end do

    end subroutine aer_rad_props_lw


!==============================================================================
! Private methods
!==============================================================================

    subroutine get_hygro_rad_props(ncol, krh, wrh, mass, ext, ssa, asm, &
                               tau, tau_w, tau_w_g, tau_w_f)

    ! Arguments
    integer,  intent(in) :: ncol
    integer,  intent(in) :: krh(nlev,ncol)  ! index for linear interpolation of optics on rh
    real(r8), intent(in) :: wrh(nlev,ncol)  ! weight for linear interpolation of optics on rh
    real(r8), intent(in) :: mass(nlev,ncol)
    real(r8), intent(in) :: ext(nrh,nswbands)
    real(r8), intent(in) :: ssa(nrh,nswbands)
    real(r8), intent(in) :: asm(nrh,nswbands)

    real(r8), intent(out) :: tau    (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w  (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_g(nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_f(nswbands,nlev,ncol)

    ! Local variables
    real(r8) :: ext1, ssa1, asm1
    integer :: icol, ilev, iswband
    !-----------------------------------------------------------------------------

    do iswband = 1, nswbands
       do icol = 1, ncol
          do ilev = 1, nlev
             ext1 = (1 + wrh(ilev,icol)) * ext(krh(ilev,icol)+1,iswband) &
                       - wrh(ilev,icol)  * ext(krh(ilev,icol),  iswband)
             ssa1 = (1 + wrh(ilev,icol)) * ssa(krh(ilev,icol)+1,iswband) &
                       - wrh(ilev,icol)  * ssa(krh(ilev,icol),  iswband)
             asm1 = (1 + wrh(ilev,icol)) * asm(krh(ilev,icol)+1,iswband) &
                       - wrh(ilev,icol)  * asm(krh(ilev,icol),  iswband)
  
             tau    (iswband,ilev,icol) = mass(ilev,icol) * ext1
             tau_w  (iswband,ilev,icol) = mass(ilev,icol) * ext1 * ssa1
             tau_w_g(iswband,ilev,icol) = mass(ilev,icol) * ext1 * ssa1 * asm1
             tau_w_f(iswband,ilev,icol) = mass(ilev,icol) * ext1 * ssa1 * asm1 * asm1
          enddo
       enddo
    enddo

    end subroutine get_hygro_rad_props 


    subroutine get_nonhygro_rad_props(ncol, mass, ext, ssa, asm, &
                                  tau, tau_w, tau_w_g, tau_w_f)

    ! Arguments
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: mass(nlev,ncol)
    real(r8), intent(in) :: ext(nswbands)
    real(r8), intent(in) :: ssa(nswbands)
    real(r8), intent(in) :: asm(nswbands)

    real(r8), intent(out) :: tau    (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w  (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_g(nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_f(nswbands,nlev,ncol)


    ! Local variables
    integer  :: iswband
    real(r8) :: ext1, ssa1, asm1
    !-----------------------------------------------------------------------------
    
    do iswband = 1, nswbands
       ext1 = ext(iswband)
       ssa1 = ssa(iswband)
       asm1 = asm(iswband)
       tau    (iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * ext1
       tau_w  (iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * ext1 * ssa1
       tau_w_g(iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * ext1 * ssa1 * asm1
       tau_w_f(iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * ext1 * ssa1 * asm1 * asm1
    enddo

    end subroutine get_nonhygro_rad_props

    
    subroutine get_volcanic_radius_rad_props(ncol, n_mu_samples, mass, r_ext, r_scat, r_ascat, r_mu, &
                                  tau, tau_w, tau_w_g, tau_w_f)

    ! Arguments
    integer,  intent(in) :: ncol
    integer,  intent(in) :: n_mu_samples
    real(r8), intent(in) :: mass(nlev,ncol)
    real(r8), intent(in) :: r_ext(nswbands,n_mu_samples)
    real(r8), intent(in) :: r_scat(nswbands,n_mu_samples)
    real(r8), intent(in) :: r_ascat(nswbands,n_mu_samples)
    real(r8), intent(in) :: r_mu(n_mu_samples) ! log(radius) domain of mass-specific optics

    real(r8), intent(out) :: tau    (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w  (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_g(nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_f(nswbands,nlev,ncol)


    ! Local variables
    integer  :: iswband
    real(r8) :: g

    integer  :: idx                             ! index to radius in physics buffer
    real(r8) :: mu(nlev,ncol)                   ! log(geometric mean radius of volcanic aerosol)
    integer  :: imu                             ! index into table values of log radius
    integer  :: kmu, nmu
    real(r8) :: wmu, mutrunc, r_mu_max, r_mu_min
 
    ! interpolated values from table
    real(r8) :: ext(nswbands)
    real(r8) :: scat(nswbands)
    real(r8) :: ascat(nswbands)

    integer :: i, k ! column level iterator
    !-----------------------------------------------------------------------------

    tau    =0._r8                 
    tau_w  =0._r8                 
    tau_w_g=0._r8                 
    tau_w_f=0._r8                  

    ! get microphysical properties for volcanic aerosols
    !idx = pbuf_get_index('VOLC_RAD_GEOM')
    !call pbuf_get_field(pbuf, idx, geometric_radius )

    ! interpolate in radius
    ! caution: clip the table with no warning when outside bounds
    nmu = size(r_mu)
    r_mu_max = r_mu(nmu)
    r_mu_min = r_mu(1)
    do i = 1, ncol
       do k = 1, nlev
          if(geometric_radius(k,i) > 0._r8) then
             mu(k,i) = log(geometric_radius(k,i))
          else
             mu(k,i) = 0._r8
          endif
          mutrunc = max(min(mu(k,i),r_mu_max),r_mu_min)
          kmu = max(min(1 + (mutrunc-r_mu_min)/(r_mu_max-r_mu_min)*(nmu-1),nmu-1._r8),1._r8)
          wmu = max(min( (mutrunc -r_mu(kmu)) / (r_mu(kmu+1) - r_mu(kmu)) ,1._r8),0._r8)
          do iswband = 1, nswbands
             ext(iswband) =  &
                ((1._r8 - wmu) * r_ext(iswband, kmu  ) + &
                (wmu) * r_ext(iswband, kmu+1))
             scat(iswband) =  &
                ((1._r8 - wmu) * r_scat(iswband, kmu  ) + &
                (wmu) * r_scat(iswband, kmu+1))
             ascat(iswband) =  &
                ((1._r8 - wmu) * r_ascat(iswband, kmu  ) + &
                (wmu) * r_ascat(iswband, kmu+1))
             if (scat(iswband).gt.0._r8) then
                g = ascat(iswband)/scat(iswband)
             else
                g=0._r8
             endif
             tau    (iswband,k,i) = mass(k,i) * ext(iswband)  
             tau_w  (iswband,k,i) = mass(k,i) * scat(iswband)  
             tau_w_g(iswband,k,i) = mass(k,i) * ascat(iswband)  
             tau_w_f(iswband,k,i) = mass(k,i) * g * ascat(iswband)  
          end do
       enddo
    enddo

    end subroutine get_volcanic_radius_rad_props


    subroutine get_volcanic_rad_props(ncol, mass, ext, scat, ascat, &
                                  tau, tau_w, tau_w_g, tau_w_f)

    ! Arguments
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: mass(nlev,ncol)
    real(r8), intent(in) :: ext(nswbands)
    real(r8), intent(in) :: scat(nswbands)
    real(r8), intent(in) :: ascat(nswbands)

    real(r8), intent(out) :: tau    (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w  (nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_g(nswbands,nlev,ncol)
    real(r8), intent(out) :: tau_w_f(nswbands,nlev,ncol)


    ! Local variables
    integer  :: iswband
    real(r8) :: g
    !-----------------------------------------------------------------------------
    
    do iswband = 1, nswbands
       if (scat(iswband).gt.0._r8) then
          g = ascat(iswband)/scat(iswband)
       else
          g=0._r8
       endif
       tau    (iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * ext(iswband) 
       tau_w  (iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * scat(iswband) 
       tau_w_g(iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * ascat(iswband) 
       tau_w_f(iswband,1:nlev,1:ncol) = mass(1:nlev,1:ncol) * g * ascat(iswband) 
    enddo

    end subroutine get_volcanic_rad_props


 end module aer_rad_props
