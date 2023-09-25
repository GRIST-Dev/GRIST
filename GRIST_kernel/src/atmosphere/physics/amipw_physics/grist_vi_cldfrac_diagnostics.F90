
module grist_vi_cldfrac_diagnostics

    use grist_constants,                only: i4, r8, one, zero

    implicit none
    private
    save

    public :: cloud_diagnostics_calc
          
contains

    subroutine cloud_diagnostics_calc(ncell, nlev, nlevp, &
                                      cldfrac,mpressureFace,mpressureFull, &
                                      cldtot, cldlow, cldmed, cldhgh)
    !use grist_physics_update,           only: old_time_level
! io
    integer ,   intent(in)    :: ncell, nlev, nlevp
    real(r8),   intent(in)    :: cldfrac(ncell,nlev)
    real(r8),   intent(in)    :: mpressureFace(ncell,nlevp), mpressureFull(ncell,nlev)
    real(r8),   intent(out)   :: cldtot(ncell), cldlow(ncell), cldmed(ncell), cldhgh(ncell)
! local
    !integer  :: i,k
    !real(r8) :: rgrav
    !real(r8) :: cld(nlev,ncol)         ! cloud fraction
    !real(r8) :: cwp   (nlev,ncol)      ! in-cloud cloud (total) water path
    !real(r8) :: gicewp(nlev,ncol)      ! grid-box cloud ice water path
    !real(r8) :: gliqwp(nlev,ncol)      ! grid-box cloud liquid water path
    !real(r8) :: gwp   (nlev,ncol)      ! grid-box cloud (total) water path
    !real(r8) :: tgicewp(ncol)          ! Vertically integrated ice water path
    !real(r8) :: tgliqwp(ncol)          ! Vertically integrated liquid water path
    !real(r8) :: tgwp   (ncol)          ! Vertically integrated (total) cloud water path

    !real(r8) :: tpw    (ncol)          ! total precipitable water
    !real(r8) :: hl     (ncol)          ! Liquid water scale height
    !real(r8) :: clwpold(nlev,ncol)     ! Presribed cloud liq. h2o path
    integer  :: nmxrgn (ncell)          ! Number of maximally overlapped regions
    real(r8) :: pmxrgn (ncell,nlevp)   ! Maximum values of pressure for each


!    cld(:,:) = pstate_cam%macrop_cld_at_pc_full_level%f(old_time_level,:,:)
!
! in wrf-physics, there is no "in-cloud" concept for liquid and ice mixing ratio,
! so we may calculate using dynamical vars
!

    ! Compute liquid and ice water paths
    !do k = 1, nlev
    !    do i = 1, ncol
    !        gicewp(k,i) = pstate_cam%microp_iciwp_at_pc_full_level%f(k,i) * cld(k,i)
    !        gliqwp(k,i) = pstate_cam%microp_iclwp_at_pc_full_level%f(k,i) * cld(k,i)
    !    end do
    !end do

    ! Determine parameters for maximum/random overlap
    call cldovrlap(ncell, nlev, nlevp, mpressureFace, cldfrac, nmxrgn, pmxrgn)

    ! Cloud cover diagnostics
    !call cloud_cover_diags_out(ncol, cld, pstate%pressure_at_pc_full_level%f, nmxrgn, pmxrgn)
    call cldsav(ncell, nlev, nlevp, cldfrac, mpressureFull, cldtot, cldlow, cldmed, cldhgh, nmxrgn, pmxrgn)

    !tgicewp(:ncol) = 0._r8
    !tgliqwp(:ncol) = 0._r8

    !do k = 1, nlev
    !    tgicewp(:ncol)  = tgicewp(:ncol) + gicewp(k,:ncol)
    !    tgliqwp(:ncol)  = tgliqwp(:ncol) + gliqwp(k,:ncol)
    !end do

    !tgwp(:ncol)      = tgicewp(:ncol) + tgliqwp(:ncol)
    !gwp(:nlev,:ncol) = gicewp(:nlev,:ncol) + gliqwp(:nlev,:ncol)
    !cwp(:nlev,:ncol) = pstate_cam%microp_iciwp_at_pc_full_level%f(:nlev,:ncol)  &
    !                 + pstate_cam%microp_iclwp_at_pc_full_level%f(:nlev,:ncol)

! output:
! gwp(GCLDLWP), tgwp(TGCLDCWP), gliqwp(TGCLDLWP), tgicewp(TGCLDIWP), cwp(ICLDTWP), 
! pstate_cam%microp_iciwp_at_pc_full_level%f(ICLDIWP)
! pstate_cam%microp_dei_at_pc_full_level%f(dei_cloud)
! pstate_cam%microp_mu_at_pc_full_level%f(mu_cloud)
! pstate_cam%microp_lambdac_at_pc_full_level%f(lambda_cloud)

    !Compute total preciptable water in column (in mm)
    !tpw(:ncol) = 0.0_r8
    !rgrav = 1.0_r8/gravity
    !do k = 1, nlev
        !do i = 1, ncol
        !    tpw(i) = tpw(i) + pstate%delp_at_pc_full_level%f(k,i)*pstate%tracer_mxrt_at_pc_full_level%f(1,k,i)*rgrav
        !end do
    !end do

    ! Diagnostic liquid water path (old specified form)
    !call cldclw(ncol, pstate%z_at_pc_face_level%f, clwpold, tpw, hl)

! output: clwpold(SETLWP), hl(LWSH)

    end subroutine cloud_diagnostics_calc

   subroutine cldovrlap(ncol    ,pver, pverp, pint    ,cld     ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Partitions each column into regions with clouds in neighboring layers.
! This information is used to implement maximum overlap in these regions
! with random overlap between them.
! On output,
!    nmxrgn contains the number of regions in each column
!    pmxrgn contains the interface pressures for the lower boundaries of
!           each region! 
! Method: 

! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------

    implicit none
!
! Input arguments
!
    integer,  intent(in) :: ncol                ! number of atmospheric columns
    integer,  intent(in) :: pver                ! chunk identifier
    integer,  intent(in) :: pverp               ! chunk identifier
    real(r8), intent(in) :: pint(ncol,pverp)   ! Interface pressure
    real(r8), intent(in) :: cld(ncol,pver)     ! Fractional cloud cover
!
! Output arguments
!
    real(r8), intent(out) :: pmxrgn(ncol,pverp)! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
    integer nmxrgn(ncol)                    ! Number of maximally overlapped regions
!
!---------------------------Local variables-----------------------------
!
    integer i                    ! Longitude index
    integer k                    ! Level index
    integer n                    ! Max-overlap region counter

    real(r8) pnm(ncol,pverp)    ! Interface pressure

    logical cld_found            ! Flag for detection of cloud
    logical cld_layer(pver)      ! Flag for cloud in layer
!
!------------------------------------------------------------------------
!

    do i = 1, ncol
       cld_found = .false.
       cld_layer(:) = cld(i,:) > zero
       pmxrgn(i,:) = zero
       pnm(i,:)=pint(i,:)*10._r8
       n = 1
       do k = 1, pver
          if (cld_layer(k) .and.  .not. cld_found) then
             cld_found = .true.
          else if ( .not. cld_layer(k) .and. cld_found) then
             cld_found = .false.
             if (count(cld_layer(k:pver)) == 0) then
                exit
             endif
             pmxrgn(i,n) = pnm(i,k)
             n = n + 1
          endif
       end do
       pmxrgn(i,n) = pnm(i,pverp)
       nmxrgn(i) = n
    end do

    return
  end subroutine cldovrlap
    
    ! Purpose: 
    ! Partitions each column into regions with clouds in neighboring layers.
    ! This information is used to implement maximum overlap in these regions
    ! with random overlap between them.
    ! On output,
    !    nmxrgn contains the number of regions in each column
    !    pmxrgn contains the interface pressures for the lower boundaries of 
    !    each region! 
    ! Author: W. Collins

    subroutine cldovrlap_cam5(ncol, nlev, nlevp, pint, cld, nmxrgn, pmxrgn)
! io
    integer,  intent(in)  :: ncol                ! number of atmospheric columns
    integer,  intent(in)  :: nlev                ! 
    integer,  intent(in)  :: nlevp               ! 
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

    end subroutine cldovrlap_cam5

  subroutine cldsav(ncol  , pver, pverp, cld  , pmid    ,&
                    cldtot  ,cldlow  ,cldmed  , cldhgh  ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute total & 3 levels of cloud fraction assuming maximum-random overlap.
! Pressure ranges for the 3 cloud levels are specified.
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: W. Collins
! 
!-----------------------------------------------------------------------
   !use shr_kind_mod, only: r8 => shr_kind_r8
   !use ppgrid

   implicit none
!------------------------------Parameters-------------------------------
   real(r8) plowmax             ! Max prs for low cloud cover range
   real(r8) plowmin             ! Min prs for low cloud cover range
   real(r8) pmedmax             ! Max prs for mid cloud cover range
   real(r8) pmedmin             ! Min prs for mid cloud cover range
   real(r8) phghmax             ! Max prs for hgh cloud cover range
   real(r8) phghmin             ! Min prs for hgh cloud cover range
!
   parameter (plowmax = 120000._r8,plowmin = 70000._r8, &
              pmedmax =  70000._r8,pmedmin = 40000._r8, &
              phghmax =  40000._r8,phghmin =  5000._r8)

   real(r8) ptypmin(4)
   real(r8) ptypmax(4)

   data ptypmin /phghmin, plowmin, pmedmin, phghmin/
   data ptypmax /plowmax, plowmax, pmedmax, phghmax/
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
!   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in)  :: ncol                 ! number of atmospheric columns
   integer, intent(in)  :: pver, pverp
   real(r8), intent(in) :: cld(ncol,pver)     ! Cloud fraction
   real(r8), intent(in) :: pmid(ncol,pver)    ! Level pressures
   real(r8), intent(in) :: pmxrgn(ncol,pverp) ! Maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc

   integer, intent(in) :: nmxrgn(ncol)        ! Number of maximally overlapped regions
!
! Output arguments
!
   real(r8), intent(out) :: cldtot(ncol)       ! Total random overlap cloud cover
   real(r8), intent(out) :: cldlow(ncol)       ! Low random overlap cloud cover
   real(r8), intent(out) :: cldmed(ncol)       ! Middle random overlap cloud cover
   real(r8), intent(out) :: cldhgh(ncol)       ! High random overlap cloud cover

!
!---------------------------Local workspace-----------------------------
!
   integer i,k                  ! Longitude,level indices
   integer irgn(ncol)          ! Max-overlap region index
   integer max_nmxrgn           ! maximum value of nmxrgn over columns
   integer ityp                 ! Type counter
   real(r8) clrsky(ncol)       ! Max-random clear sky fraction
   real(r8) clrskymax(ncol)    ! Maximum overlap clear sky fraction
!
!-----------------------------------------------------------------------
!
! Initialize region number
!
   max_nmxrgn = -1
   do i=1,ncol
      max_nmxrgn = max(max_nmxrgn,nmxrgn(i))
   end do

   do ityp = 1, 4
      irgn(1:ncol) = 1
      do k =1,max_nmxrgn-1
         do i=1,ncol
            if (pmxrgn(i,irgn(i)) < ptypmin(ityp) .and. irgn(i) < nmxrgn(i)) then
               irgn(i) = irgn(i) + 1
            end if
         end do
      end do
!
! Compute cloud amount by estimating clear-sky amounts
!
      clrsky(1:ncol)    = one
      clrskymax(1:ncol) = one
      do k = 1, pver
         do i=1,ncol
            if (pmid(i,k) >= ptypmin(ityp) .and. pmid(i,k) <= ptypmax(ityp)) then
               if (pmxrgn(i,irgn(i)) < pmid(i,k) .and. irgn(i) < nmxrgn(i)) then
                  irgn(i) = irgn(i) + 1
                  clrsky(i) = clrsky(i) * clrskymax(i)
                  clrskymax(i) = one
               endif
               clrskymax(i) = min(clrskymax(i),one-cld(i,k))
            endif
         end do
      end do
      if (ityp == 1) cldtot(1:ncol) = one - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 2) cldlow(1:ncol) = one - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 3) cldmed(1:ncol) = one - (clrsky(1:ncol) * clrskymax(1:ncol))
      if (ityp == 4) cldhgh(1:ncol) = one - (clrsky(1:ncol) * clrskymax(1:ncol))
   end do

   return
end subroutine cldsav

! Purpose: 
! Compute total & 3 levels of cloud fraction assuming maximum-random overlap.
! Pressure ranges for the 3 cloud levels are specified.
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! Author: W. Collins
    subroutine cldsav_cam5(ncol, nlev, nlevp, cld, pmid, cldtot, cldlow, cldmed, cldhgh, nmxrgn, pmxrgn)
! io
    integer, intent(in)   :: ncol                ! number of atmospheric columns
    integer, intent(in)   :: nlev                ! 
    integer, intent(in)   :: nlevp               ! 
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

    end subroutine cldsav_cam5

end module grist_vi_cldfrac_diagnostics
