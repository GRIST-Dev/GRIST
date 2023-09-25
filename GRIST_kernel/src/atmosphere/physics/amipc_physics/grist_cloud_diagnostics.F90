!======================================================
!
!  Created by LiXiaohan on 20/5/13.
!  
!  Put cloud physical specifications on the history tape
!  Modified from code that computed cloud optics
!
!  Adopted from CAM5
!======================================================

module grist_cloud_diagnostics

    use grist_constants,                only: i4, r8,   &
                                              gravity
    use grist_nml_module,               only: nlev, nlevp
    use grist_PhysC_data_structure,     only: pstate_cam
    use grist_physics_data_structure,   only: pstate
    use grist_handle_error,             only: endrun
    use grist_mpi

    implicit none
    private
    save

    public :: cloud_diagnostics_init,   &
              cloud_diagnostics_calc,   &
              end_of_cloud_diagnostics, &
              phys_diagnostics_calc,    &
              cldovrlap
          
contains
    
    subroutine cloud_diagnostics_init(ncol)
! io
    integer ,   intent(in)    :: ncol
    
! Initialize pstate
    allocate(pstate_cam%diag_cloud_tot%f(ncol))
    allocate(pstate_cam%diag_cloud_low%f(ncol))
    allocate(pstate_cam%diag_cloud_med%f(ncol))
    allocate(pstate_cam%diag_cloud_hgh%f(ncol))
    allocate(pstate_cam%diag_relhum%f(nlev,ncol))
    allocate(pstate_cam%diag_z_at_500hpa%f(ncol))
    allocate(pstate_cam%diag_u_at_850hpa%f(ncol))
    allocate(pstate_cam%diag_u_at_200hpa%f(ncol))
    allocate(pstate_cam%diag_v_at_850hpa%f(ncol))
    allocate(pstate_cam%diag_tmq%f(ncol))
    allocate(pstate_cam%diag_tgicewp%f(ncol))
    allocate(pstate_cam%diag_tgliqwp%f(ncol))
 
    pstate_cam%diag_cloud_tot%pos   = 0
    pstate_cam%diag_cloud_low%pos   = 0
    pstate_cam%diag_cloud_med%pos   = 0
    pstate_cam%diag_cloud_hgh%pos   = 0
    pstate_cam%diag_relhum%pos      = 0
    pstate_cam%diag_z_at_500hpa%pos = 0
    pstate_cam%diag_u_at_850hpa%pos = 0
    pstate_cam%diag_u_at_200hpa%pos = 0
    pstate_cam%diag_v_at_850hpa%pos = 0
    pstate_cam%diag_tmq%pos         = 0
    pstate_cam%diag_tgicewp%pos     = 0
    pstate_cam%diag_tgliqwp%pos     = 0

    pstate_cam%diag_cloud_tot%f     = 0._r8
    pstate_cam%diag_cloud_low%f     = 0._r8
    pstate_cam%diag_cloud_med%f     = 0._r8
    pstate_cam%diag_cloud_hgh%f     = 0._r8
    pstate_cam%diag_relhum%f        = 0._r8
    pstate_cam%diag_z_at_500hpa%f   = 0._r8
    pstate_cam%diag_u_at_850hpa%f   = 0._r8
    pstate_cam%diag_u_at_200hpa%f   = 0._r8
    pstate_cam%diag_v_at_850hpa%f   = 0._r8
    pstate_cam%diag_tmq%f           = 0._r8
    pstate_cam%diag_tgicewp%f       = 0._r8
    pstate_cam%diag_tgliqwp%f       = 0._r8
 
    end subroutine cloud_diagnostics_init


    subroutine end_of_cloud_diagnostics
    deallocate(pstate_cam%diag_cloud_tot%f)
    deallocate(pstate_cam%diag_cloud_low%f)
    deallocate(pstate_cam%diag_cloud_med%f)
    deallocate(pstate_cam%diag_cloud_hgh%f)
    deallocate(pstate_cam%diag_relhum%f)
    deallocate(pstate_cam%diag_z_at_500hpa%f)
    deallocate(pstate_cam%diag_u_at_850hpa%f)
    deallocate(pstate_cam%diag_u_at_200hpa%f)
    deallocate(pstate_cam%diag_v_at_850hpa%f)
    deallocate(pstate_cam%diag_tmq%f)
    deallocate(pstate_cam%diag_tgicewp%f)
    deallocate(pstate_cam%diag_tgliqwp%f)
 
    end subroutine end_of_cloud_diagnostics


    subroutine phys_diagnostics_calc(ncol)

    use grist_wv_saturation,    only: qsat
! io
    integer ,   intent(in)    :: ncol
! local
    integer  :: k
    real(r8) :: tem(nlev,ncol), ftem(nlev,ncol)
    real(r8) :: z3(nlev,ncol)

    do k = 1, nlev
        z3(k,:ncol) = pstate%z_at_pc_full_level%f(k,:ncol)+pstate%geop_at_pc_surface%f(:ncol)/gravity
    end do

    call vertinterp(ncol, pstate%pressure_at_pc_full_level%f, 50000._r8, z3,                               pstate_cam%diag_z_at_500hpa%f)
    call vertinterp(ncol, pstate%pressure_at_pc_full_level%f, 85000._r8, pstate%u_wind_at_pc_full_level%f, pstate_cam%diag_u_at_850hpa%f)
    call vertinterp(ncol, pstate%pressure_at_pc_full_level%f, 85000._r8, pstate%v_wind_at_pc_full_level%f, pstate_cam%diag_v_at_850hpa%f)
    call vertinterp(ncol, pstate%pressure_at_pc_full_level%f, 20000._r8, pstate%u_wind_at_pc_full_level%f, pstate_cam%diag_u_at_200hpa%f)
    

    call qsat(pstate%temp_at_pc_full_level%f(:,:ncol),      &
              pstate%pressure_at_pc_full_level%f(:,:ncol),  &
              tem, ftem)
    pstate_cam%diag_relhum%f(:,:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(1,:,:ncol)/ftem(:,:ncol)*100._r8

    ftem(:,:ncol) = pstate%tracer_mxrt_at_pc_full_level%f(1,:,:ncol)*pstate%delp_at_pc_full_level%f(:,:ncol)/gravity
    do k = 2, nlev
        ftem(1,:ncol) = ftem(1,:ncol) + ftem(k,:ncol)
    end do
    pstate_cam%diag_tmq%f(:ncol) = ftem(1,:ncol)

    contains

    ! Purpose: Vertically interpolate input array to output pressure level
    subroutine vertinterp(ncol, pmid, pout, arrin, arrout) 
! io
    integer , intent(in)  :: ncol              ! column dimension
    real(r8), intent(in)  :: pmid(nlev,ncol)   ! input level pressure levels 
    real(r8), intent(in)  :: pout              ! output pressure level 
    real(r8), intent(in)  :: arrin(nlev,ncol)  ! input  array
    real(r8), intent(out) :: arrout(ncol)      ! output array (interpolated)
! local
    integer i,k               ! indices
    integer kupper(ncol)      ! Level indices for interpolation
    real(r8) dpu              ! upper level pressure difference
    real(r8) dpl              ! lower level pressure difference
    logical found(ncol)       ! true if input levels found
    logical error             ! error flag 

    do i=1,ncol
       found(i)  = .false.
       kupper(i) = 1
    end do
    error = .false.
    !     
    ! Store level indices for interpolation. 
    ! If all indices for this level have been found, 
    ! do the interpolation 
    !     
    do k=1,nlev-1
       do i=1,ncol
          if ((.not. found(i)) .and. pmid(k,i)<pout .and. pout<=pmid(k+1,i)) then
             found(i) = .true.
             kupper(i) = k
          end if
       end do
    end do
    !
    ! If we've fallen through the k=1,nlev-1 loop, we cannot interpolate and
    ! must extrapolate from the bottom or top data level for at least some
    ! of the longitude points.
    !
    do i=1,ncol
       if (pout <= pmid(1,i)) then
          arrout(i) = arrin(1,i)
       else if (pout >= pmid(nlev,i)) then
          arrout(i) = arrin(nlev,i)
       else if (found(i)) then
          dpu = pout - pmid(kupper(i),i)
          dpl = pmid(kupper(i)+1,i) - pout
          arrout(i) = (arrin(kupper(i),i)*dpl + arrin(kupper(i)+1,i)*dpu)/(dpl + dpu)
       else
          error = .true.
       end if
    end do
    !     
    ! Error check
    !
    if (error) then
       call endrun ('VERTINTERP: ERROR FLAG')
    end if

    return

    end subroutine vertinterp

    end subroutine phys_diagnostics_calc


    subroutine cloud_diagnostics_calc(ncol)

    use grist_physics_update,           only: old_time_level

! io
    integer ,   intent(in)    :: ncol
! local
    integer  :: i,k
    real(r8) :: rgrav
    real(r8) :: cld(nlev,ncol)         ! cloud fraction
    real(r8) :: cwp   (nlev,ncol)      ! in-cloud cloud (total) water path
    real(r8) :: gicewp(nlev,ncol)      ! grid-box cloud ice water path
    real(r8) :: gliqwp(nlev,ncol)      ! grid-box cloud liquid water path
    real(r8) :: gwp   (nlev,ncol)      ! grid-box cloud (total) water path
    real(r8) :: tgwp   (ncol)          ! Vertically integrated (total) cloud water path

    real(r8) :: tpw    (ncol)          ! total precipitable water
    real(r8) :: hl     (ncol)          ! Liquid water scale height
    real(r8) :: clwpold(nlev,ncol)     ! Presribed cloud liq. h2o path

    integer  :: nmxrgn (ncol)          ! Number of maximally overlapped regions
    real(r8) :: pmxrgn (nlevp,ncol)    ! Maximum values of pressure for each


    cld(:,:) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,:)

    ! Compute liquid and ice water paths
    do k = 1, nlev
        do i = 1, ncol
            gicewp(k,i) = pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) * cld(k,i)
            gliqwp(k,i) = pstate_cam%microp_iclwp_at_pc_full_level%f(k,i) * cld(k,i)
        end do
    end do

    ! Determine parameters for maximum/random overlap
    call cldovrlap(ncol, pstate%pressure_at_pc_face_level%f, cld, nmxrgn, pmxrgn)

    ! Cloud cover diagnostics
    call cloud_cover_diags_out(ncol, cld, pstate%pressure_at_pc_full_level%f, nmxrgn, pmxrgn)

    pstate_cam%diag_tgicewp%f(:ncol) = 0._r8
    pstate_cam%diag_tgliqwp%f(:ncol) = 0._r8

    do k = 1, nlev
        pstate_cam%diag_tgicewp%f(:ncol)  = pstate_cam%diag_tgicewp%f(:ncol) + gicewp(k,:ncol)
        pstate_cam%diag_tgliqwp%f(:ncol)  = pstate_cam%diag_tgliqwp%f(:ncol) + gliqwp(k,:ncol)
    end do

    tgwp(:ncol)      = pstate_cam%diag_tgicewp%f(:ncol) + pstate_cam%diag_tgliqwp%f(:ncol)
    gwp(:nlev,:ncol) = gicewp(:nlev,:ncol) + gliqwp(:nlev,:ncol)
    cwp(:nlev,:ncol) = pstate_cam%microp_iciwp_at_pc_full_level%f(:nlev,:ncol)  &
                     + pstate_cam%microp_iclwp_at_pc_full_level%f(:nlev,:ncol)

! output:
! gwp(GCLDLWP), tgwp(TGCLDCWP), gliqwp(TGCLDLWP), tgicewp(TGCLDIWP), cwp(ICLDTWP), 
! pstate_cam%microp_iciwp_at_pc_full_level%f(ICLDIWP)
! pstate_cam%microp_dei_at_pc_full_level%f(dei_cloud)
! pstate_cam%microp_mu_at_pc_full_level%f(mu_cloud)
! pstate_cam%microp_lambdac_at_pc_full_level%f(lambda_cloud)

    !Compute total preciptable water in column (in mm)
    tpw(:ncol) = 0.0_r8
    rgrav = 1.0_r8/gravity
    do k = 1, nlev
        do i = 1, ncol
            tpw(i) = tpw(i) + pstate%delp_at_pc_full_level%f(k,i)*pstate%tracer_mxrt_at_pc_full_level%f(1,k,i)*rgrav
        end do
    end do

    ! Diagnostic liquid water path (old specified form)
    call cldclw(ncol, pstate%z_at_pc_face_level%f, clwpold, tpw, hl)

! output: clwpold(SETLWP), hl(LWSH)

    end subroutine cloud_diagnostics_calc

    
    ! Purpose: 
    ! Partitions each column into regions with clouds in neighboring layers.
    ! This information is used to implement maximum overlap in these regions
    ! with random overlap between them.
    ! On output,
    !    nmxrgn contains the number of regions in each column
    !    pmxrgn contains the interface pressures for the lower boundaries of 
    !    each region! 
    ! Author: W. Collins
    subroutine cldovrlap(ncol, pint, cld, nmxrgn, pmxrgn)
! io
    integer,  intent(in)  :: ncol                ! number of atmospheric columns
    real(r8), intent(in)  :: pint(nlevp,ncol)    ! Interface pressure
    real(r8), intent(in)  :: cld(nlev,ncol)      ! Fractional cloud cover
    integer,  intent(out) :: nmxrgn(ncol)        ! Number of maximally overlapped regions
    real(r8), intent(out) :: pmxrgn(nlevp,ncol)  ! Maximum values of pressure for each
! local 
    integer     :: i, k, n
    real(r8)    :: pnm(nlevp,ncol)               ! Interface pressure
    logical     :: cld_found                     ! Flag for detection of cloud
    logical     :: cld_layer(nlev)               ! Flag for cloud in layer

    do i = 1, ncol
        cld_found = .false.
        cld_layer(:) = cld(:,i) > 0.0_r8
        pmxrgn(:,i) = 0.0_r8
        pnm(:,i)=pint(:,i)*10._r8
        n = 1
        do k = 1, nlev
            if (cld_layer(k) .and.  .not. cld_found) then
                cld_found = .true.
            else if ( .not. cld_layer(k) .and. cld_found) then
                cld_found = .false.
                if (count(cld_layer(k:nlev)) == 0) then
                    exit
                endif
                pmxrgn(n,i) = pnm(k,i)
                n = n + 1
            endif
        end do
        pmxrgn(n,i) = pnm(nlevp,i)
        nmxrgn(i)   = n
    end do

    end subroutine cldovrlap


! Purpose: 
! Evaluate cloud liquid water path clwp (g/m**2)
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: J.T. Kiehl
    subroutine cldclw(ncol, zi, clwp, tpw, hl)
! io
    integer, intent(in)   :: ncol                ! number of atmospheric columns
    real(r8), intent(in)  :: zi(nlevp,ncol)      ! height at layer interfaces(m)
    real(r8), intent(in)  :: tpw(ncol)           ! total precipitable water (mm)
    real(r8), intent(out) :: clwp(nlev,ncol)     ! cloud liquid water path (g/m**2)
    real(r8), intent(out) :: hl(ncol)            ! liquid water scale height
! local
    integer  :: i,k                  ! longitude, level indices
    real(r8) :: clwc0                ! reference liquid water concentration (g/m**3)
    real(r8) :: emziohl(nlevp,ncol) ! exp(-zi/hl)
    real(r8) :: rhl(ncol)           ! 1/hl

! Set reference liquid water concentration

    clwc0 = 0.21_r8

! Diagnose liquid water scale height from precipitable water

    do i=1,ncol
       hl(i)  = 700.0_r8*log(max(tpw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0_r8/hl(i)
    end do

! Evaluate cloud liquid water path (vertical integral of exponential fn)

    do k=1, nlevp
       do i=1, ncol
          emziohl(k,i) = exp(-zi(k,i)*rhl(i))
       end do
    end do
    do k=1, nlev
       do i=1, ncol
          clwp(k,i) = clwc0*hl(i)*(emziohl(k+1,i) - emziohl(k,i))
       end do
    end do

    return
  end subroutine cldclw


    subroutine cloud_cover_diags_out(ncol, cld, pmid, nmxrgn, pmxrgn)
! io
    integer,  intent(in)  :: ncol                ! number of atmospheric columns
    real(r8), intent(in)  :: cld(nlev,ncol)      ! Fractional cloud cover
    real(r8), intent(in)  :: pmid(nlev,ncol)     ! model level pressure
    integer,  intent(in)  :: nmxrgn(ncol)        ! Number of maximally overlapped regions
    real(r8), intent(in)  :: pmxrgn(nlevp,ncol)  ! Maximum values of pressure for each

    call cldsav(ncol, cld, pmid,            & 
                pstate_cam%diag_cloud_tot%f,    &
                pstate_cam%diag_cloud_low%f,    &
                pstate_cam%diag_cloud_med%f,    &
                pstate_cam%diag_cloud_hgh%f,    &
                nmxrgn, pmxrgn)

    end subroutine cloud_cover_diags_out


    ! Purpose: 
    ! Compute total & 3 levels of cloud fraction assuming maximum-random overlap.
    ! Pressure ranges for the 3 cloud levels are specified.
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! Author: W. Collins
    subroutine cldsav(ncol, cld, pmid, cldtot, cldlow, cldmed, cldhgh, nmxrgn, pmxrgn)
! io
    integer, intent(in)   :: ncol                ! number of atmospheric columns
    real(r8), intent(in)  :: cld(nlev,ncol)      ! Fractional cloud cover
    real(r8), intent(in)  :: pmid(nlev,ncol)     ! model level pressure
    integer,  intent(in)  :: nmxrgn(ncol)        ! Number of maximally overlapped regions
    real(r8), intent(in)  :: pmxrgn(nlevp,ncol)  ! Maximum values of pressure for each
    !    maximally overlapped region.
    !    0->pmxrgn(i,1) is range of pressure for
    !    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
    !    2nd region, etc
    real(r8), intent(out) :: cldtot(ncol)        ! Total random overlap cloud cover
    real(r8), intent(out) :: cldlow(ncol)        ! Low random overlap cloud cover
    real(r8), intent(out) :: cldmed(ncol)        ! Middle random overlap cloud cover
    real(r8), intent(out) :: cldhgh(ncol)        ! High random overlap cloud cover
! local 
    integer     :: i,k                           ! Longitude,level indices
    integer     :: irgn(ncol)                    ! Max-overlap region index
    integer     :: max_nmxrgn                    ! maximum value of nmxrgn over columns
    integer     :: ityp                          ! Type counter
    real(r8)    :: clrsky(ncol)                  ! Max-random clear sky fraction
    real(r8)    :: clrskymax(ncol)               ! Maximum overlap clear sky fraction
    real(r8)    :: plowmax                       ! Max prs for low cloud cover range
    real(r8)    :: plowmin                       ! Min prs for low cloud cover range
    real(r8)    :: pmedmax                       ! Max prs for mid cloud cover range
    real(r8)    :: pmedmin                       ! Min prs for mid cloud cover range
    real(r8)    :: phghmax                       ! Max prs for hgh cloud cover range
    real(r8)    :: phghmin                       ! Min prs for hgh cloud cover range

    parameter (plowmax = 120000._r8,plowmin = 70000._r8, &
               pmedmax =  70000._r8,pmedmin = 40000._r8, &
               phghmax =  40000._r8,phghmin =  5000._r8)

    real(r8)    :: ptypmin(4)
    real(r8)    :: ptypmax(4)
    
    data ptypmin /phghmin, plowmin, pmedmin, phghmin/
    data ptypmax /plowmax, plowmax, pmedmax, phghmax/

    ! Initialize region number
    max_nmxrgn = -1
    do i=1,ncol
       max_nmxrgn = max(max_nmxrgn,nmxrgn(i))
    end do

    do ityp = 1, 4
       irgn(1:ncol) = 1
       do k =1,max_nmxrgn-1
          do i=1,ncol
             if (pmxrgn(irgn(i),i) < ptypmin(ityp) .and. irgn(i) < nmxrgn(i)) then
                irgn(i) = irgn(i) + 1
             end if
          end do
       end do

    ! Compute cloud amount by estimating clear-sky amounts

       clrsky(1:ncol)    = 1.0_r8
       clrskymax(1:ncol) = 1.0_r8
       do k = 1, nlev
          do i=1,ncol
             if (pmid(k,i) >= ptypmin(ityp) .and. pmid(k,i) <= ptypmax(ityp)) then
                if (pmxrgn(irgn(i),i) < pmid(k,i) .and. irgn(i) < nmxrgn(i)) then
                   irgn(i) = irgn(i) + 1
                   clrsky(i) = clrsky(i) * clrskymax(i)
                   clrskymax(i) = 1.0_r8
                endif
                clrskymax(i) = min(clrskymax(i),1.0_r8-cld(k,i))
             endif
          end do
       end do
       if (ityp == 1) cldtot(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
       if (ityp == 2) cldlow(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
       if (ityp == 3) cldmed(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
       if (ityp == 4) cldhgh(1:ncol) = 1.0_r8 - (clrsky(1:ncol) * clrskymax(1:ncol))
    end do

    end subroutine cldsav

end module grist_cloud_diagnostics

