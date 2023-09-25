!======================================================
!
!  Created by LiXiaohan on 19/5/13.
!  interface to RRTMG
!  Longwave radiation calculations.
!======================================================

 module radlw
    use grist_constants,                    only: i4, r8
    use grist_nml_module,                   only: nlev, nlevp
    use parrrtm,                            only: nbndlw, ngptlw
    use grist_handle_error,                 only: endrun
    use radconstants,                       only: nlwbands
    use grist_mpi

    implicit none
    
    private
    save

    public :: radlw_init ,   &! initialize constants
              rad_rrtmg_lw    ! driver for longwave radiation code

    ! Private module data
    integer :: ntoplw    ! top level to solve for longwave cooling

   
 contains

    subroutine rad_rrtmg_lw(ncol     ,rrtmg_levs,r_state,                  &
                            pmid     ,aer_lw_abs,cld       ,tauc_lw,       &
                            cicewp_in,csnowp_in ,cliqwp_in ,               &
                            rei_in   ,res_in    ,rel_in    ,               &    
                            qrl      ,qrlc      ,                          &
                            flns     ,flnt      ,flnsc     ,flntc  ,       &
                            flut     ,fldt      ,flus      ,flds   ,       &
                            flutc    ,fldtc     ,flusc     ,fldsc  ,       &
                            fnl      ,fcnl      ,lu        ,ld     ,       &
                            flg_4DDA_lw, coszrs)

    use mcica_subcol_gen_lw, only: mcica_subcol_lw
    use grist_constants,     only: cpair=>cp
    use rrtmg_state,         only: rrtmg_state_t
    use rrtmg_lw_rad,        only: rrtmg_lw, rrtmg_4DDA_lw

! io
    integer, intent(in) :: ncol                  ! number of atmospheric columns
    integer, intent(in) :: rrtmg_levs            ! number of levels rad is applied

! Input arguments which are only passed to other routines
    type(rrtmg_state_t), intent(in) :: r_state

    real(r8), intent(in) :: pmid(nlev,ncol)     ! Level pressure (Pascals)

    real(r8), intent(in) :: aer_lw_abs (nbndlw,nlev,ncol) ! aerosol absorption optics depth (LW)

    real(r8), intent(in) :: cld(nlev,ncol)      ! Cloud cover
    real(r8), intent(in) :: tauc_lw(nbndlw,nlev,ncol)   ! Cloud longwave optical depth by band
    logical, intent(in)  :: FLG_4DDA_LW   ! linhan
    real(r8), intent(in) :: coszrs(ncol) !cosine of solar zenith

    real(r8), intent(in) :: cicewp_in(nlev,ncol)
    real(r8), intent(in) :: csnowp_in(nlev,ncol)
    real(r8), intent(in) :: cliqwp_in(nlev,ncol)
    real(r8), intent(in) :: rei_in(nlev,ncol)
    real(r8), intent(in) :: res_in(nlev,ncol)
    real(r8), intent(in) :: rel_in(nlev,ncol)

    real(r8), intent(out) :: qrl (nlev,ncol)     ! Longwave heating rate
    real(r8), intent(out) :: qrlc(nlev,ncol)     ! Clearsky longwave heating rate
    real(r8), intent(out) :: flns(ncol)          ! Surface cooling flux
    real(r8), intent(out) :: flnt(ncol)          ! Net outgoing flux
    real(r8), intent(out) :: flut(ncol)          ! Upward flux at top of model
    real(r8), intent(out) :: fldt(ncol)          ! Downward flux at top of model
    real(r8), intent(out) :: flus(ncol)          ! Upward flux at surface
    real(r8), intent(out) :: flds(ncol)          ! Downward flux at surface
 
    real(r8), intent(out) :: flnsc(ncol)         ! Clear sky surface cooing
    real(r8), intent(out) :: flntc(ncol)         ! Net clear sky outgoing flux
    real(r8), intent(out) :: flutc(ncol)         ! Upward clear-sky flux at top of model
    real(r8), intent(out) :: fldtc(ncol)         ! Downward clear-sky flux at top of model
    real(r8), intent(out) :: flusc(ncol)         ! Upward longwave clear flux at surface
    real(r8), intent(out) :: fldsc(ncol)         ! Downward longwave clear flux at surface
    real(r8), intent(out) :: fcnl(nlevp,ncol)    ! clear sky net flux at interfaces
    real(r8), intent(out) :: fnl(nlevp,ncol)     ! net flux at interfaces

    real(r8), intent(out) :: lu(nlwbands,nlevp,ncol) ! longwave spectral flux up
    real(r8), intent(out) :: ld(nlwbands,nlevp,ncol) ! longwave spectral flux down

! local
    integer :: i, k, kk, nbnd         ! indices

    real(r8) :: ful(nlevp,ncol)     ! Total upwards longwave flux
    real(r8) :: fsul(nlevp,ncol)    ! Clear sky upwards longwave flux
    real(r8) :: fdl(nlevp,ncol)     ! Total downwards longwave flux
    real(r8) :: fsdl(nlevp,ncol)    ! Clear sky downwards longwv flux

    integer :: inflglw               ! Flag for cloud parameterization method
    integer :: iceflglw              ! Flag for ice cloud param method
    integer :: liqflglw              ! Flag for liquid cloud param method
    integer :: icld                  ! Flag for cloud overlap method
                                     ! 0=clear, 1=random, 2=maximum/random, 3=maximum

    real(r8) :: tsfc(ncol)          ! surface temperature
    real(r8) :: emis(nbndlw,ncol)   ! surface emissivity

    real(r8) :: taua_lw(nbndlw,rrtmg_levs-1,ncol)     ! aerosol optical depth by band

    real(r8), parameter :: dps = 1._r8/86400._r8 ! Inverse of seconds per day

    ! Cloud arrays for McICA 
    integer, parameter :: nsubclw = ngptlw       ! rrtmg_lw g-point (quadrature point) dimension
    integer :: permuteseed                       ! permute seed for sub-column generator

    real(r8) :: cicewp(rrtmg_levs-1,ncol)   ! in-cloud cloud ice water path
    real(r8) :: csnowp(rrtmg_levs-1,ncol)   ! in-cloud cloud snow water path linhan
    real(r8) :: cliqwp(rrtmg_levs-1,ncol)   ! in-cloud cloud liquid water path
    real(r8) :: rei(rrtmg_levs-1,ncol)      ! ice particle effective radius (microns)
    real(r8) :: res(rrtmg_levs-1,ncol)      ! snow particle effective radius (microns)linhan
    real(r8) :: rel(rrtmg_levs-1,ncol)      ! liquid particle radius (micron)

    real(r8) :: cld_stolw(nsubclw,rrtmg_levs-1,ncol)     ! cloud fraction (mcica)
    real(r8) :: cicewp_stolw(nsubclw,rrtmg_levs-1,ncol)  ! cloud ice water path (mcica)
    real(r8) :: csnowp_stolw(nsubclw,rrtmg_levs-1,ncol)  ! cloud sno water path (mcica)linhan
    real(r8) :: cliqwp_stolw(nsubclw,rrtmg_levs-1,ncol)  ! cloud liquid water path (mcica)
    real(r8) :: rei_stolw(rrtmg_levs-1,ncol)               ! ice particle size (mcica)
    real(r8) :: res_stolw(rrtmg_levs-1,ncol)               ! sno particle size (mcica)!linhan
    real(r8) :: rel_stolw(rrtmg_levs-1,ncol)               ! liquid particle size (mcica)
    real(r8) :: tauc_stolw(nsubclw,rrtmg_levs-1,ncol)    ! cloud optical depth (mcica - optional)

    ! Includes extra layer above model top
    real(r8) :: uflx(rrtmg_levs+1,ncol)  ! Total upwards longwave flux
    real(r8) :: uflxc(rrtmg_levs+1,ncol) ! Clear sky upwards longwave flux
    real(r8) :: dflx(rrtmg_levs+1,ncol)  ! Total downwards longwave flux
    real(r8) :: dflxc(rrtmg_levs+1,ncol) ! Clear sky downwards longwv flux
    real(r8) :: hr(rrtmg_levs,ncol)      ! Longwave heating rate (K/d)
    real(r8) :: hrc(rrtmg_levs,ncol)     ! Clear sky longwave heating rate (K/d)
    real(r8) lwuflxs(nbndlw,nlevp+1,ncol)  ! Longwave spectral flux up
    real(r8) lwdflxs(nbndlw,nlevp+1,ncol)  ! Longwave spectral flux down

    ! mji/rrtmg

    ! Calculate cloud optical properties here if using CAM method, or if using one of the
    ! methods in RRTMG_LW, then pass in cloud physical properties and zero out cloud optical 
    ! properties here
    
    ! Zero optional cloud optical depth input array tauc_lw, 
    ! if inputting cloud physical properties into RRTMG_LW
    !          tauc_lw(:,:,:) = 0.
    ! Or, pass in CAM cloud longwave optical depth to RRTMG_LW
    ! do nbnd = 1, nbndlw
    !    tauc_lw(nbnd,:nlev,:ncol) = cldtau(:nlev,:ncol)
    ! end do

    ! Call mcica sub-column generator for RRTMG_LW

    ! Call sub-column generator for McICA in radiation

    ! Select cloud overlap approach (1=random, 2=maximum-random, 3=maximum)
    icld = 2
    ! Set permute seed (must be offset between LW and SW by at least 140 to insure 
    ! effective randomization)
    permuteseed = 150

    ! These fields are no longer supplied by CAM.
    !if (mpi_rank()==0) print*,'linhan0'
    if (.not. FLG_4DDA_LW)then !linhan
        cicewp = 0.0_r8
        csnowp = 0.0_r8
        cliqwp = 0.0_r8
        rei = 0.0_r8
        res = 0.0_r8
        rel = 0.0_r8
    else
        cicewp(1:rrtmg_levs-1,1:ncol) = cicewp_in(nlevp-rrtmg_levs+1:nlevp-1,1:ncol) * 1.E+3_r8 !linhan transfer from kg/m2 to g/m2
        csnowp(1:rrtmg_levs-1,1:ncol) = csnowp_in(nlevp-rrtmg_levs+1:nlevp-1,1:ncol) * 1.E+3_r8 !linhan transfer from kg/m2 to g/m2
        cliqwp(1:rrtmg_levs-1,1:ncol) = cliqwp_in(nlevp-rrtmg_levs+1:nlevp-1,1:ncol) * 1.E+3_r8 !linhan transfer from kg/m2 to g/m2
        rei(1:rrtmg_levs-1,1:ncol) = rei_in(nlevp-rrtmg_levs+1:nlevp-1,1:ncol)
        res(1:rrtmg_levs-1,1:ncol) = res_in(nlevp-rrtmg_levs+1:nlevp-1,1:ncol)  !too huge to cal the scatter
        rel(1:rrtmg_levs-1,1:ncol) = rel_in(nlevp-rrtmg_levs+1:nlevp-1,1:ncol)
    endif

    call mcica_subcol_lw(ncol, rrtmg_levs-1, icld, permuteseed, pmid(nlevp-rrtmg_levs+1:nlevp-1,1:ncol), &
       cld(nlevp-rrtmg_levs+1:nlevp-1,1:ncol), cicewp, csnowp, cliqwp, rei, res, rel, tauc_lw(:, nlevp-rrtmg_levs+1:nlevp-1, :ncol), &
       cld_stolw, cicewp_stolw,csnowp_stolw, cliqwp_stolw, rei_stolw, res_stolw, rel_stolw, tauc_stolw)
    !if (mpi_rank()==0) print*,'linhan1'


    ! Call RRTMG_LW model
    !
    ! Set input flags for cloud parameterizations
    ! Use separate specification of ice and liquid cloud optical depth.
    ! Use either Ebert and Curry ice parameterization (iceflglw = 0 or 1), 
    ! or use Key (Streamer) approach (iceflglw = 2), or use Fu method
    ! (iceflglw = 3), and Hu/Stamnes for liquid (liqflglw = 1).
    ! For use in Fu method (iceflglw = 3), rei is converted in RRTMG_LW
    ! from effective radius to generalized effective size using the
    ! conversion of D. Mitchell, JAS, 2002.  For ice particles outside
    ! the effective range of either the Key or Fu approaches, the 
    ! Ebert and Curry method is applied. 

    ! Input CAM cloud optical depth directly
    if (.not.flg_4DDA_lw)then !linhan
        inflglw = 0
        iceflglw = 0
        liqflglw = 0
    ! Use E&C approach for ice to mimic CAM3
    !   inflglw = 2
    !   iceflglw = 1
    !   liqflglw = 1
    ! Use merged Fu and E&C params for ice
    else  !used when 4DDA is applied linhan
       inflglw = 2
       iceflglw = 3
       liqflglw = 1
    endif

    ! Convert incoming water amounts from specific humidity to vmr as needed;
    ! Convert other incoming molecular amounts from mmr to vmr as needed;
    ! Convert pressures from Pa to hPa;
    ! Set surface emissivity to 1.0 here, this is treated in land surface model;
    ! Set surface temperature
    ! Set aerosol optical depth to zero for now

    emis(:nbndlw,:ncol) = 1._r8
    tsfc(:ncol) = r_state%tlev(rrtmg_levs+1,:ncol)
    taua_lw(:nbndlw,1:rrtmg_levs-1,:ncol) = aer_lw_abs(:nbndlw,nlevp-rrtmg_levs+1:nlevp-1,:ncol)

    lu(:,:,1:ncol) = 0.0_r8
    ld(:,:,1:ncol) = 0.0_r8

    if (FLG_4DDA_LW) then !linhan
    call rrtmg_4DDA_lw(ncol   ,rrtmg_levs   ,icld    ,coszrs ,         &
         r_state%pmidmb  ,r_state%pintmb  ,r_state%tlay    ,r_state%tlev    ,tsfc    ,r_state%h2ovmr, &
         r_state%o3vmr   ,r_state%co2vmr  ,r_state%ch4vmr  ,r_state%o2vmr   ,r_state%n2ovmr  ,r_state%cfc11vmr,r_state%cfc12vmr, &
         r_state%cfc22vmr,r_state%ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
         cld_stolw,tauc_stolw,cicewp_stolw,csnowp_stolw, cliqwp_stolw ,rei, res, rel, &
         taua_lw, &
         uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, &
         lwuflxs, lwdflxs)
    else
    call rrtmg_lw(ncol   ,rrtmg_levs      ,icld            ,                 &
         r_state%pmidmb  ,r_state%pintmb  ,r_state%tlay    ,r_state%tlev    ,tsfc    ,r_state%h2ovmr, &
         r_state%o3vmr   ,r_state%co2vmr  ,r_state%ch4vmr  ,r_state%o2vmr   ,r_state%n2ovmr  ,r_state%cfc11vmr,r_state%cfc12vmr, &
         r_state%cfc22vmr,r_state%ccl4vmr ,emis    ,inflglw ,iceflglw,liqflglw, &
         cld_stolw,tauc_stolw,cicewp_stolw,csnowp_stolw, cliqwp_stolw ,rei, res, rel, &
         taua_lw, &
         uflx    ,dflx    ,hr      ,uflxc   ,dflxc   ,hrc, &
         lwuflxs, lwdflxs)
    endif
    !if (mpi_rank()==0) print*,'linhan2'

    ! All longitudes: store history tape quantities
    ! Flux units are in W/m2 on output from rrtmg_lw and contain output for
    ! extra layer above model top with vertical indexing from bottom to top.
    ! Heating units are in K/d on output from RRTMG and contain output for
    ! extra layer above model top with vertical indexing from bottom to top.
    ! Heating units are converted to J/kg/s below for use in CAM. 

    flut(:ncol)  = uflx (rrtmg_levs,:ncol)
    fldt(:ncol)  = dflx (rrtmg_levs,:ncol)
    flutc(:ncol) = uflxc(rrtmg_levs,:ncol)
    fldtc(:ncol) = dflxc(rrtmg_levs,:ncol)
    flds(:ncol)  = dflx (1,:ncol)
    flus(:ncol)  = uflx (1,:ncol)
    fldsc(:ncol) = dflxc(1,:ncol)
    flusc(:ncol) = uflxc(1,:ncol)

    flns(:ncol)  = uflx (1,:ncol) - dflx (1,:ncol)
    flnsc(:ncol) = uflxc(1,:ncol) - dflxc(1,:ncol)
    flnt(:ncol)  = uflx (rrtmg_levs,:ncol) - dflx (rrtmg_levs,:ncol)
    flntc(:ncol) = uflxc(rrtmg_levs,:ncol) - dflxc(rrtmg_levs,:ncol)

    !
    ! Reverse vertical indexing here for CAM arrays to go from top to bottom.
    !
    ful = 0._r8
    fdl = 0._r8
    fsul = 0._r8
    fsdl = 0._r8
    ful (nlevp-rrtmg_levs+1:nlevp,:ncol)= uflx(rrtmg_levs:1:-1,:ncol)
    fdl (nlevp-rrtmg_levs+1:nlevp,:ncol)= dflx(rrtmg_levs:1:-1,:ncol)
    fsul(nlevp-rrtmg_levs+1:nlevp,:ncol)=uflxc(rrtmg_levs:1:-1,:ncol)
    fsdl(nlevp-rrtmg_levs+1:nlevp,:ncol)=dflxc(rrtmg_levs:1:-1,:ncol)

    !--------------LiXH has not completed scm_crm_mode------------->
    !if (single_column.and.scm_crm_mode) then
    !   call outfld('FUL     ',ful,pcols,lchnk)
    !   call outfld('FDL     ',fdl,pcols,lchnk)
    !   call outfld('FULC    ',fsul,pcols,lchnk)
    !   call outfld('FDLC    ',fsdl,pcols,lchnk)
    !endif
    !<-------------LiXH has not completed scm_crm_mode--------------
    
    fnl(:,:ncol) = ful(:,:ncol) - fdl(:,:ncol)
    ! mji/ cam excluded this?
    fcnl(:,:ncol) = fsul(:,:ncol) - fsdl(:,:ncol)

    ! Pass longwave heating to CAM arrays and convert from K/d to J/kg/s
    qrl = 0._r8
    qrlc = 0._r8
    qrl (nlevp-rrtmg_levs+1:nlev,:ncol)=hr (rrtmg_levs-1:1:-1,:ncol)*cpair*dps
    qrlc(nlevp-rrtmg_levs+1:nlev,:ncol)=hrc(rrtmg_levs-1:1:-1,:ncol)*cpair*dps

    !linhan test
   ! if(mpi_rank()==0)then
   ! print*,'linhan test lw' !linhan
   ! print*,rrtmg_levs+1,uflx(rrtmg_levs+1,1:ncol),dflx(rrtmg_levs+1,1:ncol) !linhan
   ! k=rrtmg_levs
   ! print*,k
   ! print*,hr(k,1:ncol),uflx(k,1:ncol),dflx(k,1:ncol) !linhan
   ! print*,hrc(k,1:ncol),uflxc(k,1:ncol),dflxc(k,1:ncol) !linhan
   !     do k=rrtmg_levs-1,1,-1   !linhan
   !         print*,k !linhan
   !         print*,hr(k,1:ncol),uflx(k,1:ncol),dflx(k,1:ncol) !linhan
   !         !print*,hrc(k,1:ncol),uflxc(k,1:ncol),dflxc(k,1:ncol) !linhan
   !         print*,cicewp(k,1),csnowp(k,1), cliqwp(k,1)
   !         print*,rei(k,1),res(k,1), rel(k,1)
   !     enddo  !linhan
   ! print*,'linhan test lwend' !linhan
   ! end if

    ! Return 0 above solution domain
    if ( ntoplw > 1 )then
       qrl(:ntoplw-1,:ncol) = 0._r8
       qrlc(:ntoplw-1,:ncol) = 0._r8
    end if

    ! Pass spectral fluxes, reverse layering
    ! order=(/3,1,2/) maps the first index of lwuflxs to the third index of lu.
    lu(:,nlevp-rrtmg_levs+1:nlevp,:ncol) = lwuflxs(:,rrtmg_levs:1:-1,:ncol)
    
    ld(:,nlevp-rrtmg_levs+1:nlevp,:ncol) = lwdflxs(:,rrtmg_levs:1:-1,:ncol)
    
    end subroutine rad_rrtmg_lw



! Purpose: Initialize various constants for radiation scheme.
    subroutine radlw_init()
    use rrtmg_lw_init,          only: rrtmg_lw_ini
    use grist_constants,        only: p00
    use grist_hpe_constants,    only: eta_full
 
    implicit none

    !local
    real(r8) :: pref_full(nlev)
    integer  :: k

    pref_full(:) = eta_full(:)*p00 

    ! If the top model level is above ~90 km (0.1 Pa), set the top level to compute
    ! longwave cooling to about 80 km (1 Pa)
    if (pref_full(1) .lt. 0.1_r8) then
       do k = 1, nlev
          if (pref_full(k) .lt. 1._r8) ntoplw  = k
       end do
    else
       ntoplw  = 1
    end if

    call rrtmg_lw_ini

    end subroutine radlw_init





 end module radlw
