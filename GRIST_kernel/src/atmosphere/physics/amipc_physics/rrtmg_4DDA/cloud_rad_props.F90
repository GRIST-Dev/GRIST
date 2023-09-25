!======================================================
!
!  Created by LiXiaohan on 19/8/13.
!======================================================

 module cloud_rad_props

    use grist_constants,                    only: r8
    use grist_nml_module,                   only: nlev, nlevp
    use grist_handle_error,                 only: endrun
    use ebert_curry,                        only: scalefactor
    use radconstants,                       only: nswbands, nlwbands
    use grist_PhysC_data_structure,         only: pstate_cam
    use oldcloud,                           only: oldcloud_lw,              &
                                                  old_liq_get_rad_props_lw, &    
                                                  old_ice_get_rad_props_lw, &
                                                  oldcloud_init
    use grist_wrap_nf
    use grist_mpi
 
    implicit none
    private
    save
    
    public :: cloud_rad_props_init,          &
              cloud_rad_props_get_sw,        & ! return SW optical props of total bulk aerosols
              cloud_rad_props_get_lw,        & ! return LW optical props of total bulk aerosols
              get_ice_optics_sw,             & ! return Mitchell SW ice radiative properties
              ice_cloud_get_rad_props_lw,    & ! Mitchell LW ice rad props
              get_liquid_optics_sw,          & ! return Conley SW rad props
              liquid_cloud_get_rad_props_lw, & ! return Conley LW rad props
              snow_cloud_get_rad_props_lw,   &
              get_snow_optics_sw,            &
              end_of_cloud_rad_props

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
    
    real(r8) :: dmin = 2. *  13. / scalefactor
    real(r8) :: dmax = 2. * 130. / scalefactor

! Minimum cloud amount (as a fraction of the grid-box area) to 
! distinguish from clear sky
 
    real(r8), parameter :: cldmin = 1.0e-80_r8

! Decimal precision of cloud amount (0 -> preserve full resolution;
! 10^-n -> preserve n digits of cloud amount)
 
    real(r8), parameter :: cldeps = 0.0_r8

 contains
    
    subroutine cloud_rad_props_init()

    use grist_rad_constituents, only: iceopticsfile, liqopticsfile
    use slingo,                 only: slingo_rad_props_init
    use ebert_curry,            only: ec_rad_props_init
    use grist_physics_iofile,   only: getfile

    character(len=256) :: liquidfile 
    character(len=256) :: icefile 
    character(len=256) :: locfn
    integer, parameter :: omode = 0  ! nf_nowrite
    integer :: ncid, dimid, mudimid, lambdadimid
    integer :: mu_id, lambda_id, ext_sw_liq_id, ssa_sw_liq_id,  &
               asm_sw_liq_id, abs_lw_liq_id, d_dimid, d_id,     &
               ext_sw_ice_id, ssa_sw_ice_id, asm_sw_ice_id,     &
               abs_lw_ice_id 
    integer :: f_nlwbands, f_nswbands
    ! zhangyiAdds:
    integer :: comm, ierr

    comm = grist_comm()

    liquidfile = liqopticsfile
    icefile    = iceopticsfile

    call slingo_rad_props_init
    call ec_rad_props_init
    call oldcloud_init

    ! read liquid cloud optics
    IF(mpi_rank().eq.0)then
    call getfile(trim(liquidfile), locfn, 0)
    call wrap_open(trim(locfn), omode, ncid)

    call wrap_inq_dimid(ncid, 'lw_band', dimid)
    call wrap_inq_dimlen(ncid, dimid, f_nlwbands)

    call wrap_inq_dimid(ncid, 'sw_band', dimid)
    call wrap_inq_dimlen(ncid, dimid, f_nswbands)

    call wrap_inq_dimid(ncid, 'mu', mudimid)
    call wrap_inq_dimlen(ncid, mudimid, nmu)

    call wrap_inq_dimid(ncid, 'lambda_scale', lambdadimid)
    call wrap_inq_dimlen(ncid, lambdadimid, nlambda)
    END IF
    call mpi_bcast(nmu    , 1, MPI_INTEGER, 0, comm, ierr)
    call mpi_bcast(nlambda, 1, MPI_INTEGER, 0, comm, ierr)

    allocate(g_mu(nmu))
    allocate(g_lambda(nmu,nlambda))
    allocate(ext_sw_liq(nmu,nlambda,nswbands) )
    allocate(ssa_sw_liq(nmu,nlambda,nswbands))
    allocate(asm_sw_liq(nmu,nlambda,nswbands))
    allocate(abs_lw_liq(nmu,nlambda,nlwbands))

    IF(mpi_rank().eq.0)then
    call wrap_inq_varid(ncid, 'mu', mu_id)
    call wrap_get_var_realx(ncid, mu_id, g_mu)

    call wrap_inq_varid(ncid, 'lambda', lambda_id)
    call wrap_get_var_realx(ncid, lambda_id, g_lambda)

    call wrap_inq_varid(ncid, 'k_ext_sw', ext_sw_liq_id)
    call wrap_get_var_realx(ncid, ext_sw_liq_id, ext_sw_liq)

    call wrap_inq_varid(ncid, 'ssa_sw', ssa_sw_liq_id)
    call wrap_get_var_realx(ncid, ssa_sw_liq_id, ssa_sw_liq)

    call wrap_inq_varid(ncid, 'asm_sw', asm_sw_liq_id)
    call wrap_get_var_realx(ncid, asm_sw_liq_id, asm_sw_liq)

    call wrap_inq_varid(ncid, 'k_abs_lw', abs_lw_liq_id)
    call wrap_get_var_realx(ncid, abs_lw_liq_id, abs_lw_liq)
    ENDIF 
    call mpi_bcast(g_mu      , nmu                 , MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(g_lambda  , nmu*nlambda         , MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(ext_sw_liq, nmu*nlambda*nswbands, MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(ssa_sw_liq, nmu*nlambda*nswbands, MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(asm_sw_liq, nmu*nlambda*nswbands, MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(abs_lw_liq, nmu*nlambda*nlwbands, MPI_REAL8, 0, comm, ierr)

    ext_sw_liq = ext_sw_liq / 0.9970449e3_r8 
    abs_lw_liq = abs_lw_liq / 0.9970449e3_r8 

    ! read ice cloud optics
    IF(mpi_rank().eq.0)then
    call getfile(trim(icefile), locfn, 0)
    call wrap_open(trim(locfn), omode, ncid)

    call wrap_inq_dimid(ncid, 'lw_band', dimid)
    call wrap_inq_dimlen(ncid, dimid, f_nlwbands)

    call wrap_inq_dimid(ncid, 'sw_band', dimid)
    call wrap_inq_dimlen(ncid, dimid, f_nswbands)

    call wrap_inq_dimid(ncid, 'd_eff', d_dimid)
    call wrap_inq_dimlen(ncid, d_dimid, n_g_d)
    ENDIF
    call mpi_bcast(n_g_d , 1, MPI_INTEGER, 0, comm, ierr)

    allocate(g_d_eff(n_g_d))
    allocate(ext_sw_ice(n_g_d,nswbands))
    allocate(ssa_sw_ice(n_g_d,nswbands))
    allocate(asm_sw_ice(n_g_d,nswbands))
    allocate(abs_lw_ice(n_g_d,nlwbands))

    IF(mpi_rank().eq.0)then
    call wrap_inq_varid(ncid, 'd_eff', d_id)
    call wrap_get_var_realx(ncid, d_id, g_d_eff)

    call wrap_inq_varid(ncid, 'sw_ext', ext_sw_ice_id)
    call wrap_get_var_realx(ncid, ext_sw_ice_id, ext_sw_ice)

    call wrap_inq_varid(ncid, 'sw_ssa', ssa_sw_ice_id)
    call wrap_get_var_realx(ncid, ssa_sw_ice_id, ssa_sw_ice)

    call wrap_inq_varid(ncid, 'sw_asm', asm_sw_ice_id)
    call wrap_get_var_realx(ncid, asm_sw_ice_id, asm_sw_ice)

    call wrap_inq_varid(ncid, 'lw_abs', abs_lw_ice_id)
    call wrap_get_var_realx(ncid, abs_lw_ice_id, abs_lw_ice)
    ENDIF
    call mpi_bcast(g_d_eff   , n_g_d         , MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(ext_sw_ice, n_g_d*nswbands, MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(ssa_sw_ice, n_g_d*nswbands, MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(asm_sw_ice, n_g_d*nswbands, MPI_REAL8, 0, comm, ierr)
    call mpi_bcast(abs_lw_ice, n_g_d*nlwbands, MPI_REAL8, 0, comm, ierr)

    end subroutine cloud_rad_props_init


    subroutine end_of_cloud_rad_props()
    deallocate(g_mu)
    deallocate(g_lambda)
    deallocate(ext_sw_liq)
    deallocate(ssa_sw_liq)
    deallocate(asm_sw_liq)
    deallocate(abs_lw_liq)
    deallocate(g_d_eff)
    deallocate(ext_sw_ice)
    deallocate(ssa_sw_ice)
    deallocate(asm_sw_ice)
    deallocate(abs_lw_ice)
   
    end subroutine end_of_cloud_rad_props


! Purpose : return totaled (across all species) layer tau, omega, g, f
!           for all spectral interval for aerosols affecting the climate
    subroutine cloud_rad_props_get_sw(ncol, tau, tau_w, tau_w_g, tau_w_f,  &
                                      diagnosticindex, oldliq, oldice)
    ! io
    integer, intent(in)             :: ncol
    integer, optional, intent(in)   :: diagnosticindex    ! index (if present) to radiation diagnostic information
    logical, optional, intent(in)   :: oldliq,oldice
    real(r8), intent(out) :: tau    (nswbands,nlev,ncol) ! aerosol extinction optical depth
    real(r8), intent(out) :: tau_w  (nswbands,nlev,ncol) ! aerosol single scattering albedo * tau
    real(r8), intent(out) :: tau_w_g(nswbands,nlev,ncol) ! aerosol assymetry parameter * tau * w
    real(r8), intent(out) :: tau_w_f(nswbands,nlev,ncol) ! aerosol forward scattered fraction * tau * w
    ! local
    ! rad properties for liquid clouds
    real(r8) :: liq_tau    (nswbands,nlev,ncol) ! aerosol extinction optical depth
    real(r8) :: liq_tau_w  (nswbands,nlev,ncol) ! aerosol single scattering albedo * tau
    real(r8) :: liq_tau_w_g(nswbands,nlev,ncol) ! aerosol assymetry parameter * tau * w
    real(r8) :: liq_tau_w_f(nswbands,nlev,ncol) ! aerosol forward scattered fraction * tau * w

    ! rad properties for ice clouds
    real(r8) :: ice_tau    (nswbands,nlev,ncol) ! aerosol extinction optical depth
    real(r8) :: ice_tau_w  (nswbands,nlev,ncol) ! aerosol single scattering albedo * tau
    real(r8) :: ice_tau_w_g(nswbands,nlev,ncol) ! aerosol assymetry parameter * tau * w
    real(r8) :: ice_tau_w_f(nswbands,nlev,ncol) ! aerosol forward scattered fraction * tau * w

    ! initialize to conditions that would cause failure
    tau     (:,:,:) = -100._r8
    tau_w   (:,:,:) = -100._r8
    tau_w_g (:,:,:) = -100._r8
    tau_w_f (:,:,:) = -100._r8

    ! initialize layers to accumulate od's
    tau    (:,:,1:ncol)  = 0._r8
    tau_w  (:,:,1:ncol)  = 0._r8
    tau_w_g(:,:,1:ncol)  = 0._r8
    tau_w_f(:,:,1:ncol)  = 0._r8

    call get_liquid_optics_sw(ncol, liq_tau, liq_tau_w, liq_tau_w_g, liq_tau_w_f)
    call get_ice_optics_sw   (ncol, ice_tau, ice_tau_w, ice_tau_w_g, ice_tau_w_f)

    tau    (:,:,1:ncol)  =  liq_tau    (:,:,1:ncol)  + ice_tau    (:,:,1:ncol) 
    tau_w  (:,:,1:ncol)  =  liq_tau_w  (:,:,1:ncol)  + ice_tau_w  (:,:,1:ncol) 
    tau_w_g(:,:,1:ncol)  =  liq_tau_w_g(:,:,1:ncol)  + ice_tau_w_g(:,:,1:ncol) 
    tau_w_f(:,:,1:ncol)  =  liq_tau_w_f(:,:,1:ncol)  + ice_tau_w_f(:,:,1:ncol) 

    end subroutine cloud_rad_props_get_sw


! Purpose: Compute cloud longwave absorption optical depth
!          cloud_rad_props_get_lw() is called by radlw() 
    subroutine cloud_rad_props_get_lw(ncol, cld_abs_od, diagnosticindex, oldliq, oldice, oldcloud)
    ! io
    integer,             intent(in)  :: ncol
    real(r8),            intent(out) :: cld_abs_od(nlwbands,nlev,ncol) ! [fraction] absorption optical depth, per layer
    integer, optional,   intent(in)  :: diagnosticindex
    logical, optional,   intent(in)  :: oldliq  ! use old liquid optics
    logical, optional,   intent(in)  :: oldice  ! use old ice optics
    logical, optional,   intent(in)  :: oldcloud  ! use old optics for both (b4b)
    ! local
    integer :: bnd_idx     ! LW band index
    integer :: i           ! column index
    integer :: k           ! lev index
    ! rad properties for liquid clouds
    real(r8) :: liq_tau_abs_od(nlwbands,nlev,ncol) ! liquid cloud absorption optical depth
    ! rad properties for ice clouds
    real(r8) :: ice_tau_abs_od(nlwbands,nlev,ncol) ! ice cloud absorption optical depth

    ! compute optical depths cld_absod 
    cld_abs_od = 0._r8

    if(present(oldcloud))then
       if(oldcloud) then
          ! make diagnostic calls to these first to output ice and liq OD's
          !call old_liq_get_rad_props_lw(state, pbuf, liq_tau_abs_od, oldliqwp=.false.)
          !call old_ice_get_rad_props_lw(state, pbuf, ice_tau_abs_od, oldicewp=.false.)
          ! This affects climate (cld_abs_od)
          call oldcloud_lw(ncol,cld_abs_od,.false.)
          return
       endif
    endif

    if(present(oldliq))then
       if(oldliq) then
          call old_liq_get_rad_props_lw(ncol, liq_tau_abs_od, .false.)
       else
          call liquid_cloud_get_rad_props_lw(ncol, liq_tau_abs_od)
       endif
    else
       call liquid_cloud_get_rad_props_lw(ncol, liq_tau_abs_od)
    endif

    if(present(oldice))then
       if(oldice) then
          call old_ice_get_rad_props_lw(ncol, ice_tau_abs_od, .false.)
       else
          call ice_cloud_get_rad_props_lw(ncol, ice_tau_abs_od)
       endif
    else
       call ice_cloud_get_rad_props_lw(ncol, ice_tau_abs_od)
    endif
      
    cld_abs_od(:,:,1:ncol) = liq_tau_abs_od(:,:,1:ncol) + ice_tau_abs_od(:,:,1:ncol) 

    end subroutine cloud_rad_props_get_lw


    subroutine get_snow_optics_sw (ncol, tau, tau_w, tau_w_g, tau_w_f)
    ! io
    integer,intent(in)   :: ncol
    real(r8),intent(out) :: tau    (nswbands,nlev,ncol) ! extinction optical depth
    real(r8),intent(out) :: tau_w  (nswbands,nlev,ncol) ! single scattering albedo * tau
    real(r8),intent(out) :: tau_w_g(nswbands,nlev,ncol) ! assymetry parameter * tau * w
    real(r8),intent(out) :: tau_w_f(nswbands,nlev,ncol) ! forward scattered fraction * tau * w

    real(r8) :: iciwpth(nlev,ncol), dei(nlev,ncol)
    real(r8) :: dlimited ! d limited to dmin,dmax range

    integer :: i,k,i_d_grid,k_d_eff,i_swband
    real :: wd, onemwd, ext, ssa, asm

    dei(:,1:ncol)     = pstate_cam%microp_des_at_pc_full_level%f(:,1:ncol)
    iciwpth(:,1:ncol) = pstate_cam%microp_icswp_at_pc_full_level%f(:,1:ncol)
    do i = 1,ncol
       do k = 1, nlev
          dlimited = dei(k,i)! linhan 20190918 pstate_cam%microp_dei_at_pc_full_level%f(k,i) ! min(dmax,max(dei(k,i),dmin))
          if( pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) < 1.e-80_r8     &
          .or. dlimited .eq. 0._r8) then
          ! if ice water path is too small, OD := 0
             tau    (:,k,i) = 0._r8
             tau_w  (:,k,i) = 0._r8
             tau_w_g(:,k,i) = 0._r8
             tau_w_f(:,k,i) = 0._r8
          else 
             !if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
                !print*, 'dei from prognostic cldwat2m',dei(k,i)
                !print*, 'grid values of deff ice from optics file',g_d_eff
                !call endrun ('deff of ice exceeds limits')
             !endif
             ! for each cell interpolate to find weights and indices in g_d_eff grid.
             if (dlimited <= g_d_eff(1)) then
                k_d_eff = 2
                wd = 1._r8
                onemwd = 0._r8
             elseif (dlimited >= g_d_eff(n_g_d)) then
                k_d_eff = n_g_d
                wd = 0._r8
                onemwd = 1._r8 
             else
                do i_d_grid = 1, n_g_d
                   k_d_eff = i_d_grid
                   if(g_d_eff(i_d_grid) > dlimited) exit
                enddo
                wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                onemwd = 1._r8 - wd
             endif
             ! interpolate into grid and extract radiative properties
             do i_swband = 1, nswbands
                ext = wd*ext_sw_ice(k_d_eff-1,i_swband) + &
                  onemwd*ext_sw_ice(k_d_eff  ,i_swband) 
                ssa = wd*ssa_sw_ice(k_d_eff-1,i_swband) + &
                  onemwd*ssa_sw_ice(k_d_eff  ,i_swband) 
                asm = wd*asm_sw_ice(k_d_eff-1,i_swband) + &
                  onemwd*asm_sw_ice(k_d_eff  ,i_swband) 
                tau    (i_swband,k,i)=iciwpth(k,i) * ext
                tau_w  (i_swband,k,i)=iciwpth(k,i) * ext * ssa
                tau_w_g(i_swband,k,i)=iciwpth(k,i) * ext * ssa * asm
                tau_w_f(i_swband,k,i)=iciwpth(k,i) * ext * ssa * asm * asm
             enddo
          endif
       enddo
    enddo

    return
    end subroutine get_snow_optics_sw   


    subroutine get_ice_optics_sw (ncol, tau, tau_w, tau_w_g, tau_w_f)
    ! io
    integer,intent(in)   :: ncol
    real(r8),intent(out) :: tau    (nswbands,nlev,ncol) ! extinction optical depth
    real(r8),intent(out) :: tau_w  (nswbands,nlev,ncol) ! single scattering albedo * tau
    real(r8),intent(out) :: tau_w_g(nswbands,nlev,ncol) ! assymetry parameter * tau * w
    real(r8),intent(out) :: tau_w_f(nswbands,nlev,ncol) ! forward scattered fraction * tau * w
    ! local
    real(r8) :: dlimited ! d limited to dmin,dmax range
    integer  :: i,k,i_d_grid,k_d_eff,i_swband
    real(r8) :: wd, onemwd, ext, ssa, asm

    do k = 1, nlev
        do i = 1, ncol
            dlimited = pstate_cam%microp_dei_at_pc_full_level%f(k,i)   ! min(dmax,max(dei(k,i),dmin))
            if( pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) < 1.e-80_r8     &
                .or. dlimited .eq. 0._r8) then
                ! if ice water path is too small, OD := 0
                tau    (:,k,i) = 0._r8
                tau_w  (:,k,i) = 0._r8
                tau_w_g(:,k,i) = 0._r8
                tau_w_f(:,k,i) = 0._r8
            else
                !if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
                   !write(*,*) 'dei from prognostic cldwat2m',dei(k,i)
                   !write(*,*) 'grid values of deff ice from optics file',g_d_eff
                   !call endrun ('deff of ice exceeds limits')
                !endif
                ! for each cell interpolate to find weights and indices in g_d_eff grid.
                if (dlimited <= g_d_eff(1)) then
                   k_d_eff = 2
                   wd = 1._r8
                   onemwd = 0._r8
                elseif (dlimited >= g_d_eff(n_g_d)) then
                   k_d_eff = n_g_d
                   wd = 0._r8
                   onemwd = 1._r8 
                else
                   do i_d_grid = 1, n_g_d
                      k_d_eff = i_d_grid
                      if(g_d_eff(i_d_grid) > dlimited) exit
                   enddo
                   wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                   onemwd = 1._r8 - wd
                endif
                ! interpolate into grid and extract radiative properties
                do i_swband = 1, nswbands
                   ext = wd*ext_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*ext_sw_ice(k_d_eff  ,i_swband) 
                   ssa = wd*ssa_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*ssa_sw_ice(k_d_eff  ,i_swband) 
                   asm = wd*asm_sw_ice(k_d_eff-1,i_swband) + &
                     onemwd*asm_sw_ice(k_d_eff  ,i_swband) 
                   tau    (i_swband,k,i)=pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) * ext
                   tau_w  (i_swband,k,i)=pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) * ext * ssa
                   tau_w_g(i_swband,k,i)=pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) * ext * ssa * asm
                   tau_w_f(i_swband,k,i)=pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) * ext * ssa * asm * asm
                enddo
            endif
        enddo
    enddo

    return
    end subroutine get_ice_optics_sw 


    subroutine get_liquid_optics_sw(ncol, tau, tau_w, tau_w_g, tau_w_f)
    ! io
    integer,intent(in)   :: ncol
    real(r8),intent(out) :: tau    (nswbands,nlev,ncol) ! extinction optical depth
    real(r8),intent(out) :: tau_w  (nswbands,nlev,ncol) ! single scattering albedo * tau
    real(r8),intent(out) :: tau_w_g(nswbands,nlev,ncol) ! assymetry parameter * tau * w
    real(r8),intent(out) :: tau_w_f(nswbands,nlev,ncol) ! forward scattered fraction * tau * w
    ! local
    integer  :: i,k,swband

    do k = 1, nlev
        do i = 1, ncol
            if(pstate_cam%microp_lambdac_at_pc_full_level%f(k,i) > 0._r8) then 
                ! This seems to be clue from microphysics of no cloud
                call gam_liquid_sw(pstate_cam%microp_iclwp_at_pc_full_level%f(k,i),         &
                                   pstate_cam%microp_lambdac_at_pc_full_level%f(k,i),       &
                                   pstate_cam%microp_mu_at_pc_full_level%f(k,i),            &
                                   tau(1:nswbands,k,i), tau_w(1:nswbands,k,i),          &
                                   tau_w_g(1:nswbands,k,i), tau_w_f(1:nswbands,k,i))
            else
                tau(1:nswbands,k,i)     = 0._r8
                tau_w(1:nswbands,k,i)   = 0._r8
                tau_w_g(1:nswbands,k,i) = 0._r8
                tau_w_f(1:nswbands,k,i) = 0._r8
            end if
        end do
    end do

    end subroutine get_liquid_optics_sw


    subroutine liquid_cloud_get_rad_props_lw(ncol, abs_od)
    ! io
    integer,  intent(in)  :: ncol
    real(r8), intent(out) :: abs_od(nlwbands,nlev,ncol)
    ! local
    integer :: i, k

    abs_od = 0._r8

    do k = 1, nlev
       do i = 1,ncol
          if(pstate_cam%microp_lambdac_at_pc_full_level%f(k,i) > 0._r8) then 
             ! This seems to be the clue for no cloud from microphysics formulation
             call gam_liquid_lw(pstate_cam%microp_iclwp_at_pc_full_level%f(k,i),      &  
                                pstate_cam%microp_lambdac_at_pc_full_level%f(k,i),    &  
                                pstate_cam%microp_mu_at_pc_full_level%f(k,i),         &  
                                abs_od(1:nlwbands,k,i))
          else
             abs_od(1:nlwbands,k,i) = 0._r8
          endif
       enddo
    enddo

    end subroutine liquid_cloud_get_rad_props_lw


    subroutine snow_cloud_get_rad_props_lw(ncol, abs_od)
    ! io
    integer,  intent(in)  :: ncol
    real(r8), intent(out) :: abs_od(nlwbands,nlev,ncol)
    ! local
    real(r8) :: dlimited ! d limited to range dmin,dmax
    integer  :: i,k,i_d_grid,k_d_eff,i_lwband
    real(r8) :: wd, onemwd, absor
    real(r8) :: iciwpth(nlev,ncol), dei(nlev,ncol)

    abs_od = 0._r8

    dei(:,1:ncol)     = pstate_cam%microp_des_at_pc_full_level%f(:,1:ncol)
    iciwpth(:,1:ncol) = pstate_cam%microp_icswp_at_pc_full_level%f(:,1:ncol)

    do k = 1, nlev
       do i = 1, ncol
          dlimited = dei(k,i) ! min(dmax,max(dei(k,i),dmin))
          ! if ice water path is too small, OD := 0
          if( iciwpth(k,i) < 1.e-80_r8 .or. dlimited .eq. 0._r8) then
             abs_od (:,k,i) = 0._r8
          !else if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
          !   write(iulog,*) 'dlimited prognostic cldwat2m',dlimited
          !   write(iulog,*) 'grid values of deff ice from optics file',g_d_eff(1),' -> ',g_d_eff(n_g_d)
          !   !call endrun ('deff of ice exceeds limits')
          else
             ! for each cell interpolate to find weights and indices in g_d_eff grid.
             if (dlimited <= g_d_eff(1)) then
                k_d_eff = 2
                wd = 1._r8
                onemwd = 0._r8
             elseif (dlimited >= g_d_eff(n_g_d)) then
                k_d_eff = n_g_d
                wd = 0._r8
                onemwd = 1._r8 
             else
                do i_d_grid = 2, n_g_d
                   k_d_eff = i_d_grid
                   if(g_d_eff(i_d_grid) > dlimited) exit
                enddo
                wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                onemwd = 1._r8 - wd
             endif
             ! interpolate into grid and extract radiative properties
             do i_lwband = 1, nlwbands
                absor = wd*abs_lw_ice(k_d_eff-1,i_lwband) + &
                    onemwd*abs_lw_ice(k_d_eff  ,i_lwband)
                abs_od (i_lwband,k,i)=  iciwpth(k,i) * absor 
             enddo
          endif
       enddo
    enddo

    end subroutine snow_cloud_get_rad_props_lw


    subroutine ice_cloud_get_rad_props_lw(ncol, abs_od)
    ! io
    integer,  intent(in)  :: ncol
    real(r8), intent(out) :: abs_od(nlwbands,nlev,ncol)
    ! local
    real(r8) :: dlimited ! d limited to range dmin,dmax
    integer  :: k,i,i_d_grid,k_d_eff,i_lwband
    real(r8) :: wd, onemwd, absor

    abs_od = 0._r8

    do k = 1, nlev
       do i = 1, ncol
          dlimited = pstate_cam%microp_dei_at_pc_full_level%f(k,i) ! min(dmax,max(dei(k,i),dmin))
          ! if ice water path is too small, OD := 0
          if( pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) < 1.e-80_r8     &
            .or. dlimited .eq. 0._r8) then
             abs_od (:,k,i) = 0._r8
          !else if (dlimited < g_d_eff(1) .or. dlimited > g_d_eff(n_g_d)) then
          !   print*, 'dlimited prognostic cldwat2m',dlimited
          !   print*, 'grid values of deff ice from optics file',g_d_eff(1),' -> ',g_d_eff(n_g_d)
          !   !call endrun ('deff of ice exceeds limits')
          else
             ! for each cell interpolate to find weights and indices in g_d_eff grid.
             if (dlimited <= g_d_eff(1)) then
                k_d_eff = 2
                wd = 1._r8
                onemwd = 0._r8
             elseif (dlimited >= g_d_eff(n_g_d)) then
                k_d_eff = n_g_d
                wd = 0._r8
                onemwd = 1._r8 
             else
                do i_d_grid = 2, n_g_d
                   k_d_eff = i_d_grid
                   if(g_d_eff(i_d_grid) > dlimited) exit
                enddo
                wd = (g_d_eff(k_d_eff) - dlimited)/(g_d_eff(k_d_eff) - g_d_eff(k_d_eff-1))
                onemwd = 1._r8 - wd
             endif
             ! interpolate into grid and extract radiative properties
             do i_lwband = 1, nlwbands
                absor = wd*abs_lw_ice(k_d_eff-1,i_lwband) + &
                    onemwd*abs_lw_ice(k_d_eff  ,i_lwband)
                abs_od (i_lwband,k,i)=  pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) * absor 
             enddo
          endif
       enddo
    enddo

    end subroutine ice_cloud_get_rad_props_lw


    subroutine gam_liquid_lw(clwptn, lamc, pgam, abs_od)
    ! io
    real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
    real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
    real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
    real(r8), intent(out) :: abs_od(1:nlwbands)
    ! for interpolating into mu/lambda
    integer :: imu, kmu, wmu, onemwmu
    integer :: ilambda, klambda, wlambda, onemwlambda, lambdaplus, lambdaminus
    integer :: lwband ! sw band index
    real(r8) :: absc

    if (clwptn < 1.e-80_r8) then
      abs_od = 0._r8
      return
    endif

    if (pgam < g_mu(1) .or. pgam > g_mu(nmu)) then
      print*,'pgam from prognostic cldwat2m in lw',pgam
      print*,'g_mu from file',g_mu(1:nmu)
      call endrun ('pgam exceeds limits')
    endif
    do imu = 1, nmu
      kmu = imu
      if (g_mu(kmu) > pgam) exit
    enddo
    wmu = (g_mu(kmu) - pgam)/(g_mu(kmu) - g_mu(kmu-1))
    onemwmu = 1._r8 - wmu

    do ilambda = 1, nlambda
      klambda = ilambda
      if (wmu*g_lambda(kmu-1,ilambda) + onemwmu*g_lambda(kmu,ilambda) < lamc) exit
    enddo
    if (klambda <= 1 .or. klambda > nlambda) call endrun('lamc exceeds limits')
    lambdaplus = wmu*g_lambda(kmu-1,klambda  ) + onemwmu*g_lambda(kmu,klambda  )
    lambdaminus= wmu*g_lambda(kmu-1,klambda-1) + onemwmu*g_lambda(kmu,klambda-1)
    wlambda = (lambdaplus - lamc) / (lambdaplus - lambdaminus)
    onemwlambda = 1._r8 - wlambda

    do lwband = 1, nlwbands
       absc=     wlambda*    wmu*abs_lw_liq(kmu-1,klambda-1,lwband) + &
             onemwlambda*    wmu*abs_lw_liq(kmu-1,klambda  ,lwband) + &
                 wlambda*onemwmu*abs_lw_liq(kmu  ,klambda-1,lwband) + &
             onemwlambda*onemwmu*abs_lw_liq(kmu  ,klambda  ,lwband)

       abs_od(lwband) = clwptn * absc
    enddo

    return
    end subroutine gam_liquid_lw

    subroutine gam_liquid_sw(clwptn, lamc, pgam, tau, tau_w, tau_w_g, tau_w_f)
    real(r8), intent(in) :: clwptn ! cloud water liquid path new (in cloud) (in g/m^2)?
    real(r8), intent(in) :: lamc   ! prognosed value of lambda for cloud
    real(r8), intent(in) :: pgam   ! prognosed value of mu for cloud
    real(r8), intent(out) :: tau(1:nswbands), tau_w(1:nswbands), tau_w_f(1:nswbands), tau_w_g(1:nswbands)
    ! for interpolating into mu/lambda
    integer :: imu, kmu, wmu, onemwmu
    integer :: ilambda, klambda, wlambda, onemwlambda, lambdaplus, lambdaminus
    integer :: swband ! sw band index
    real(r8) :: ext, ssa, asm

    if (clwptn < 1.e-80_r8) then
      tau = 0._r8
      tau_w = 0._r8
      tau_w_g = 0._r8
      tau_w_f = 0._r8
      return
    endif

    if (pgam < g_mu(1) .or. pgam > g_mu(nmu)) then
      print*, 'pgam from prognostic cldwat2m in sw',pgam
      print*, 'g_mu from file',nmu,g_mu(1:nmu)
      call endrun ('pgam exceeds limits')
    endif
    do imu = 1, nmu
      kmu = imu
      if (g_mu(kmu) > pgam) exit
    enddo
    wmu = (g_mu(kmu) - pgam)/(g_mu(kmu) - g_mu(kmu-1))
    onemwmu = 1._r8 - wmu

    do ilambda = 1, nlambda
       klambda = ilambda
       if (wmu*g_lambda(kmu-1,ilambda) + onemwmu*g_lambda(kmu,ilambda) < lamc) exit
    enddo
    if (klambda <= 1 .or. klambda > nlambda) call endrun('lamc  exceeds limits')
       lambdaplus = wmu*g_lambda(kmu-1,klambda  ) + onemwmu*g_lambda(kmu,klambda  )
       lambdaminus= wmu*g_lambda(kmu-1,klambda-1) + onemwmu*g_lambda(kmu,klambda-1)
       wlambda = (lambdaplus - lamc) / (lambdaplus - lambdaminus)
       onemwlambda = 1._r8 - wlambda

    do swband = 1, nswbands
       ext =     wlambda*    wmu*ext_sw_liq(kmu-1,klambda-1,swband) + &
             onemwlambda*    wmu*ext_sw_liq(kmu-1,klambda  ,swband) + &
                 wlambda*onemwmu*ext_sw_liq(kmu  ,klambda-1,swband) + &
             onemwlambda*onemwmu*ext_sw_liq(kmu  ,klambda  ,swband)
       ! probably should interpolate ext*ssa
       ssa =     wlambda*    wmu*ssa_sw_liq(kmu-1,klambda-1,swband) + &
             onemwlambda*    wmu*ssa_sw_liq(kmu-1,klambda  ,swband) + &
                 wlambda*onemwmu*ssa_sw_liq(kmu  ,klambda-1,swband) + &
             onemwlambda*onemwmu*ssa_sw_liq(kmu  ,klambda  ,swband)
       ! probably should interpolate ext*ssa*asm
       asm =     wlambda*    wmu*asm_sw_liq(kmu-1,klambda-1,swband) + &
             onemwlambda*    wmu*asm_sw_liq(kmu-1,klambda  ,swband) + &
                 wlambda*onemwmu*asm_sw_liq(kmu  ,klambda-1,swband) + &
             onemwlambda*onemwmu*asm_sw_liq(kmu  ,klambda  ,swband)
       ! compute radiative properties
       tau(swband) = clwptn * ext
       tau_w(swband) = clwptn * ext * ssa
       tau_w_g(swband) = clwptn * ext * ssa * asm
       tau_w_f(swband) = clwptn * ext * ssa * asm * asm
    enddo

    return
    end subroutine gam_liquid_sw


 end module cloud_rad_props
