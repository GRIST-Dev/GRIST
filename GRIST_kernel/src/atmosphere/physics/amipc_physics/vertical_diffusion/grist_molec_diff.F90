!===================================================================================================
!
!  Created by LiXiaohan on 19/06/01, adopted from CAM5
!
!  to compute molecular diffusivity for various constituents
!
!  init_molec_diff           Initializes time independent coefficients
!  init_timestep_molec_diff  Time-step initialization for molecular diffusivity
!  compute_molec_diff        Computes constituent-independent terms for moleculuar diffusivity
!  vd_lu_qdecomp             Computes constituent-dependent terms for moleculuar diffusivity and
!                            updates terms in the triadiagonal matrix used for the implicit
!                            solution of the diffusion equation
!===================================================================================================

 module grist_molec_diff

    use grist_constants,                    only: r8, i4
    use grist_physics_data_structure,       only: phy_tracer_info
    use grist_handle_error,                 only: endrun


    implicit none
    private
    public  :: init_molec_diff,     &
               compute_molec_diff,  &
               vd_lu_qdecomp
               
! Private:
    real(r8), parameter   :: km_fac = 3.55E-7_r8         ! Molecular viscosity constant [ unit ? ]
    real(r8), parameter   :: pr_num = 1._r8              ! Prandtl number [ no unit ]
    real(r8), parameter   :: pwr    = 2._r8/3._r8        ! Exponentiation factor [ unit ? ]
    real(r8), parameter   :: d0     = 1.52E20_r8         ! Diffusion factor [ m-1 s-1 ] molec sqrt(kg/kmol/K) [ unit ? ]
                                                         ! Aerononmy, Part B, Banks and Kockarts (1973), p39
                                                         ! Note text cites 1.52E18 cm-1 ...

    real(r8)              :: rair                        ! Gas constant for dry air
    real(r8)              :: mw_dry                      ! Molecular weight of dry air
    real(r8)              :: n_avog                      ! Avogadro's number [ molec/kmol ]
    real(r8)              :: gravit     
    real(r8)              :: cpair
    real(r8)              :: kbtz                        ! Boltzman constant

    ! Warning: it seems dangerous to set ntop_molec /= 1, since it is not respected by all relevant modules.
    ! Top interface level to which molecular vertical diffusion is applied ( = 1 )
    integer(i4)           :: ntop_molec
    ! Bottom interface level to which molecular vertical diffusion is applied ( = pver )
    integer(i4)           :: nbot_molec
    real(r8), allocatable :: mw_fac(:)                   ! sqrt(1/M_q + 1/M_d) in constituent diffusivity [  unit ? ]
    real(r8), allocatable :: alphath(:)                  ! Thermal diffusion factor, -0.38 for H, 0 for others
 
contains
    subroutine init_molec_diff(ncnst, rair_in, mw_dry_in, n_avog_in, gravit_in, &
                               cpair_in, kbtz_in, pref_mid, ntop_molec_out, nbot_molec_out)
! io:
    integer(i4),  intent(in)  :: ncnst          ! Number of constituents
    real(r8), intent(in)  :: rair_in
    real(r8), intent(in)  :: mw_dry_in      ! Molecular weight of dry air
    real(r8), intent(in)  :: n_avog_in      ! Avogadro's number [ molec/kmol ]
    real(r8), intent(in)  :: gravit_in
    real(r8), intent(in)  :: cpair_in
    real(r8), intent(in)  :: kbtz_in        ! Boltzman constant
    real(r8), intent(in)  :: pref_mid(:)    ! Reference pressures
 
    integer(i4),  intent(out) :: ntop_molec_out ! Top interface level to which molecular vertical diffusion is applied ( = 1 )
    integer(i4),  intent(out) :: nbot_molec_out ! Bottom interface level to which molecular vertical diffusion is applied.
! local:
    integer(i4)           :: k              ! Level index
    integer(i4)           :: m              ! Constituent index
    integer(i4)           :: indx_H         ! Constituent index for H
    integer(i4)           :: ierr           ! Allocate error check
    real(r8), parameter   :: nbot_molec_pres = 50._r8 ! Pressure above which molecular diffusion is turned off.

    rair       = rair_in
    mw_dry     = mw_dry_in
    n_avog     = n_avog_in
    gravit     = gravit_in
    cpair      = cpair_in
    kbtz       = kbtz_in

! ---------------------------------------------------------------------------------------- !
! Molecular diffusion turned on above ~60 km (50 Pa) if model top is above ~90 km (.1 Pa). !
! ---------------------------------------------------------------------------------------- !

    ntop_molec = 1       ! Should always be 1
    nbot_molec = 0       ! Should be set below about 70 km
    molec_top_loop: do k = 1,size(pref_mid)
       if (pref_mid(k) > nbot_molec_pres) then
          nbot_molec = k-1
          exit molec_top_loop
       end if
    end do molec_top_loop

    ntop_molec_out = ntop_molec
    nbot_molec_out = nbot_molec
    
! Initialize upper boundary condition variables
! CAM: if(waccm) call ubc_init

! Molecular weight factor in constitutent diffusivity
! ***** FAKE THIS FOR NOW USING MOLECULAR WEIGHT OF DRY AIR FOR ALL TRACERS ****
 
    allocate(mw_fac(ncnst))
    do m = 1, ncnst
       mw_fac(m) = d0 * mw_dry * sqrt(1._r8/mw_dry + 1._r8/phy_tracer_info(m)%molec_weight) / n_avog
    end do


    end subroutine init_molec_diff


    integer function compute_molec_diff( pver        , ncnst     , ncol      , t        , pmid   , pint   ,             &
                                         zi          , ztodt     , kvm       , kvt      , tint   , rhoi   , tmpi2     , &
                                         kq_scal     , ubc_mmr   , ubc_flux  , dse_top  , cc_top , cd_top ,             &
                                         cnst_mw_out , cnst_fixed_ubc_out    , cnst_fixed_ubflx_out       , mw_fac_out, &
                                         ntop_molec_out , nbot_molec_out )
    
!    use upper_bc,        only : ubc_get_vals
!    use physconst,       only : cpairv, rairv, kmvis, kmcnd

! io
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
    real(r8), intent(in)    :: t(pver, ncol)             ! Temperature input
    real(r8), intent(in)    :: pmid(pver, ncol)          ! Midpoint pressures
    real(r8), intent(in)    :: pint(pver+1, ncol)        ! Interface pressures
    real(r8), intent(in)    :: zi(pver+1, ncol)          ! Interface heights
    real(r8), intent(in)    :: ztodt                     ! 2 delta-t
    
    real(r8), intent(inout) :: kvm(pver+1, ncol)         ! Viscosity ( diffusivity for momentum )
    real(r8), intent(out)   :: kvt(pver+1, ncol)         ! Kinematic molecular conductivity
    real(r8), intent(inout) :: tint(pver+1, ncol)        ! Interface temperature
    real(r8), intent(inout) :: rhoi(pver+1, ncol)        ! Density ( rho ) at interfaces
    real(r8), intent(inout) :: tmpi2(pver+1, ncol)       ! dt*(g*rho)**2/dp at interfaces

    real(r8), intent(out)   :: kq_scal(pver+1, ncol)     ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8), intent(out)   :: ubc_mmr(ncnst, ncol)      ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(out)   :: ubc_flux(ncnst)           ! Upper boundary flux [ kg/s/m^2 ]
    real(r8), intent(out)   :: cnst_mw_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubc_out(ncnst)
    logical,  intent(out)   :: cnst_fixed_ubflx_out(ncnst)
    real(r8), intent(out)   :: mw_fac_out(ncnst, pver+1, ncol) ! composition dependent mw_fac on interface level
    real(r8), intent(out)   :: dse_top(ncol)             ! dse on top boundary
    real(r8), intent(out)   :: cc_top(ncol)              ! Lower diagonal at top interface
    real(r8), intent(out)   :: cd_top(ncol)              ! cc_top * dse ubc value
    integer,  intent(out)   :: ntop_molec_out   
    integer,  intent(out)   :: nbot_molec_out   

! local 
    integer                 :: m                         ! Constituent index
    integer                 :: i                         ! Column index
    integer                 :: k                         ! Level index

    real(r8)                :: mbarvi                    ! mbarv on interface level
    real(r8)                :: km_top(ncol)              ! molecular conductivity at the top

    real(r8)                :: mkvisc                    ! Molecular kinematic viscosity c*tint**(2/3)/rho
    real(r8)                :: ubc_t(ncol)               ! Upper boundary temperature (K)

! Lixh add variables that should be init in waccmx:
    real(r8)                :: cpairv(pver, ncol)
    real(r8)                :: rairv(pver, ncol)


  ! Get upper boundary values
  !---------------------------------Lixh close this part------------------------------------
  ! if waccmx is avaiable, this part should be modified:
    print*,'ubc_t, ubc_mmr, and ubc_flux can not be set to 0, if do_molec_diff, you should modify waccmx.'
    call endrun("compute_molec_diff")
    ubc_t    = 0._r8
    ubc_mmr  = 0._r8
    ubc_flux = 0._r8
    cpairv(:pver, :ncol) = cpair
    rairv(:pver, :ncol)  = rair
  !  call ubc_get_vals( lchnk, ncol, ntop_molec, pint, zi, ubc_t, ubc_mmr, ubc_flux )
  !---------------------------------Lixh close this part------------------------------------

  ! Below are already computed, just need to be copied for output

    cnst_mw_out(:ncnst)          = phy_tracer_info(:ncnst)%molec_weight
    cnst_fixed_ubc_out(:ncnst)   = phy_tracer_info(:ncnst)%cnst_fixed_ubc
    cnst_fixed_ubflx_out(:ncnst) = phy_tracer_info(:ncnst)%cnst_fixed_ubflx
    ntop_molec_out               = ntop_molec
    nbot_molec_out               = nbot_molec

  !  Need variable mw_fac for kvt and constant otherwise
  !---------------------------------Lixh close this part----------------------------------->
  ! if waccmx is avaiable, this part should be modified:
  !  if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
  !    do m = 1, ncnst
  !      do k = ntop_molec+1, nbot_molec-1
  !        do i = 1, ncol
  !           mbarvi = 0.5_r8 * (mbarv(k,i-1,lchnk)+mbarv(k,i,lchnk))
  !           mw_fac_out(m,k,i) = d0 * mbarvi * sqrt(1._r8/mbarvi + 1._r8/cnst_mw(m)) / n_avog
  !        enddo
  !      enddo
  !      mw_fac_out(m,ntop_molec,:ncol) = 1.5_r8*mw_fac_out(m,ntop_molec+1,:ncol)-.5_r8*mw_fac_out(m,ntop_molec+2,:ncol)
  !      do k = nbot_molec, pver+1
  !        mw_fac_out(m,k,:ncol) = mw_fac_out(m,nbot_molec-1,:ncol)
  !      enddo
  !    end do
  !  else
      do k = 1, pver+1
        do i = 1, ncol
          mw_fac_out(:ncnst,k,i) = mw_fac(:ncnst)
        enddo
      enddo
  !  endif
  !<--------------------------------Lixh close this part------------------------------------
 
  ! Density and related factors for molecular diffusion and ubc.
  ! Always have a fixed upper boundary T if molecular diffusion is active. Why ?
  ! For kvt, set ubc temperature to average of next two lower interface level temperatures

  !---------------------------------Lixh close this part----------------------------------->
  ! if waccmx is avaiable, this part should be modified:
  !  if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
  !    tint(ntop_molec,:ncol) = 1.5_r8*tint(:ncol,ntop_molec+1)-.5_r8*tint(:ncol,ntop_molec+2)
  !  else
      tint (ntop_molec,:ncol) = ubc_t(:ncol)    
  !  endif
  !<--------------------------------Lixh close this part------------------------------------
    
    rhoi (ntop_molec,:ncol) = pint(ntop_molec,:ncol) / ( rairv(ntop_molec,:ncol) * tint(ntop_molec,:ncol) )
    tmpi2(ntop_molec,:ncol) = ztodt * ( gravit * rhoi(ntop_molec,:ncol))**2 &
                                    / ( pmid(ntop_molec,:ncol) - pint(ntop_molec,:ncol) )
    
  ! Compute molecular kinematic viscosity, heat diffusivity and factor for constituent diffusivity
  ! This is a key part of the code.  For WACCM-X, use constituent dependent molecular viscosity and conductivity

    kvt     = 0._r8
    kq_scal = 0._r8
  !---------------------------------Lixh close this part----------------------------------->
  ! if waccmx is avaiable, this part should be modified:
  !  if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
  !    do k = ntop_molec, nbot_molec
  !       do i = 1, ncol
  !         mkvisc  = kmvis(k,i,lchnk) / rhoi(k,i)
  !         kvm(k,i) = kvm(k,i) + mkvisc
  !         mkvisc  = kmcnd(k,i,lchnk) / rhoi(k,i)
  !         kvt(k,i) = mkvisc
  !         kq_scal(k,i) = sqrt(tint(k,i)) / rhoi(k,i)
  !       end do
  !    end do
  !  else
      do k = ntop_molec, nbot_molec
        do i = 1, ncol
          mkvisc   = km_fac * tint(k,i)**pwr / rhoi(k,i)
          kvm(k,i) = kvm(k,i) + mkvisc
          kvt(k,i) = mkvisc * pr_num * cpairv(k,i)  
          kq_scal(k,i) = sqrt(tint(k,i)) / rhoi(k,i)
        end do
      end do
  !  endif
  !<--------------------------------Lixh close this part------------------------------------
    
  ! Top boundary condition for dry static energy

    dse_top(:ncol) = cpairv(ntop_molec,:ncol) * tint(ntop_molec,:ncol) + gravit * zi(ntop_molec,:ncol)

  ! Top value of cc for dry static energy

    do i = 1, ncol
      cc_top(i) = ztodt * gravit**2 * rhoi(ntop_molec,i) * km_fac * ubc_t(i)**pwr / &
                  ( ( pint(2,i) - pint(1,i) ) * ( pmid(1,i) - pint(1,i) ) )
    enddo

    cd_top(:ncol) = cc_top(:ncol) * dse_top(:ncol)
    
    compute_molec_diff = 1
    return
  end function compute_molec_diff


! Purpose: Add the molecular diffusivity to the turbulent diffusivity for a consitutent.
!          Update the superdiagonal (ca(k)), diagonal (cb(k)) and subdiagonal (cc(k))    !
!          coefficients of the tridiagonal diffusion matrix, also ze and denominator.    !
  integer function vd_lu_qdecomp( pver  , ncol   , fixed_ubc  , mw         , ubc_mmr,    &
                                  kv    , kq_scal, mw_facm    , tmpi       , rpdel  ,    &
                                  ca    , cc     , dnom       , ze         , rhoi   ,    &
                                  tint  , ztodt  , ntop_molec , nbot_molec , cd_top ,    &
                                  pmid  , pint   , t          , m      )
! io
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncol                  ! Number of atmospheric columns
    integer,  intent(in)    :: ntop_molec
    integer,  intent(in)    :: nbot_molec
    logical,  intent(in)    :: fixed_ubc             ! Fixed upper boundary condition flag
    real(r8), intent(in)    :: kv(pver+1, ncol)      ! Eddy diffusivity
    real(r8), intent(in)    :: kq_scal(pver+1, ncol) ! Molecular diffusivity ( kq_fac*sqrt(T)*m_d/rho )
    real(r8), intent(in)    :: mw                    ! Molecular weight for this constituent
    real(r8), intent(in)    :: ubc_mmr(ncol)         ! Upper boundary mixing ratios [ kg/kg ]
    real(r8), intent(in)    :: mw_facm(pver+1, ncol) ! composition dependent sqrt(1/M_q + 1/M_d) for this constituent
    real(r8), intent(in)    :: tmpi(pver+1, ncol)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    real(r8), intent(in)    :: rpdel(pver, ncol)     ! 1./pdel ( thickness bet interfaces )
    real(r8), intent(in)    :: rhoi(pver+1, ncol)    ! Density at interfaces [ kg/m3 ]
    real(r8), intent(in)    :: tint(pver+1, ncol)    ! Interface temperature [ K ]
    real(r8), intent(in)    :: ztodt                 ! 2 delta-t [ s ]

    real(r8), intent(inout) :: ca(pver, ncol)        ! -Upper diagonal
    real(r8), intent(inout) :: cc(pver, ncol)        ! -Lower diagonal
    real(r8), intent(inout) :: dnom(pver, ncol)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1)) , 1./(b(k) - c(k)*e(k-1))
    real(r8), intent(inout) :: ze(pver, ncol)        ! Term in tri-diag. matrix system

    real(r8), intent(out)   :: cd_top(ncol)          ! Term for updating top level with ubc

    real(r8), intent(in)    :: pmid(pver, ncol)      ! midpoint pressures
    real(r8), intent(in)    :: pint(pver+1, ncol)    ! interface pressures
    real(r8), intent(in)    :: t(pver, ncol)         ! temperature
    integer,  intent(in)    :: m                     ! cnst index 
! local
    integer                 :: i                     ! Longitude index
    integer                 :: k, kp1                ! Vertical indicies
    real(r8)                :: rghd(pver+1, ncol)    ! (1/H_i - 1/H) * (rho*g)^(-1)
    real(r8)                :: kmq(pver+1, ncol)     ! Molecular diffusivity for constituent
    real(r8)                :: wrk0(ncol)            ! Work variable
    real(r8)                :: wrk1(ncol)            ! Work variable
    real(r8)                :: cb(pver, ncol)        ! - Diagonal
    real(r8)                :: gradm(pver+1, ncol)   ! 1/mbar * d(mbar)/dp *(rho*g)^2 * dt
    real(r8)                :: gradt(pver+1, ncol)   ! alphaTh*(rho*g)^2 1/T * dT/dp * dt, for now alphaTh is non-zero only for H.
    real(r8)                :: mbarvi                ! mbarv at interface

    ! --------------------------------------------------------------------- !
    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the !
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are !
    ! a combination of ca and cc; they are not required by the solver.      !
    !---------------------------------------------------------------------- !

    kmq(:,:)  = 0._r8
    cd_top(:) = 0._r8

  ! Compute difference between scale heights of constituent and dry air
!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
!
!      do k = ntop_molec+1, nbot_molec-1
!    	 do i = 1, ncol
!    	    mbarvi = 0.5_r8 * (mbarv(k-1,i,lchnk)+mbarv(k,i,lchnk))
!    	    rghd(k,i)  = gravit  / (kbtz * n_avog * tint(k,i)) * (mw - mbarvi)
!    	    rghd(k,i)  = ztodt * gravit * rhoi(k,i) * rghd(k,i)
!    	    gradm(k,i) = (mbarv(k,i,lchnk)-mbarv(k-1,i,lchnk))/(pmid(k,i)-pmid(k-1,i))/mbarvi
!    	    gradm(k,i) = ztodt * (rhoi(k,i) * gravit)**2 * gradm(k,i)
!    	    gradt(k,i) = (t(k,i)-t(k-1,i))/(pmid(k,i)-pmid(k-1,i))/tint(k,i)
!    	    gradt(k,i) = ztodt * alphath(m) * (rhoi(k,i) * gravit)**2 * gradt(k,i)
!    	 enddo
!      enddo
!      do k = ntop_molec,ntop_molec
!    	 do i = 1, ncol
!    	    mbarvi = .75_r8*mbarv(k,i,lchnk)+0.5_r8*mbarv(k+1,i,lchnk)-.25_r8*mbarv(k+2,i,lchnk)
!    	    rghd(k,i)  = gravit  / (kbtz * n_avog * tint(k,i)) * (mw - mbarvi)
!    	    rghd(k,i)  = ztodt * gravit * rhoi(k,i) * rghd(k,i)
!    	    gradm(k,i) = (mbarv(k,i,lchnk)-mbarvi)/(pmid(k,i)-pint(k,i))/(mbarv(k,i,lchnk)+mbarvi)*2._r8
!    	    gradm(k,i) = ztodt * (rhoi(k,i) * gravit)**2 * gradm(k,i)
!    	    gradt(k,i) = (t(k,i)-tint(k,i))/(pmid(k,i)-pint(k,i))/(t(k,i)+tint(k,i))*2._r8
!    	    gradt(k,i) = ztodt * alphath(m) * (rhoi(k,i) * gravit)**2 * gradt(k,i)
!    	 enddo
!      enddo
!      do k = nbot_molec,nbot_molec
!    	 do i = 1, ncol
!    	    mbarvi = mbarv(k-1,i,lchnk)
!    	    rghd(k,i)  = gravit  / (kbtz * n_avog * tint(k,i)) * (mw - mbarvi)
!    	    rghd(k,i)  = ztodt * gravit * rhoi(k,i) * rghd(k,i)
!    	    gradm(k,i) = 0._r8
!    	    gradt(k,i) = 0._r8						     ! set to zero because molecular diffusion is small at the lower boundary
!    	 enddo
!      enddo
!
!    else

      do k = ntop_molec, nbot_molec
    	 do i = 1, ncol
    	    rghd(k,i) = gravit / ( kbtz * n_avog * tint(k,i) ) * ( mw - mw_dry )
    	    rghd(k,i) = ztodt * gravit * rhoi(k,i) * rghd(k,i)
    	 enddo
      enddo
      
!    endif
!<----------------------------------Lixh close this part-----------------------------------

    !-------------------- !
    ! Molecular diffusion !
    !-------------------- !

    do k = nbot_molec - 1, ntop_molec, -1
       kp1 = k + 1
       kmq(kp1,:ncol)  = kq_scal(kp1,:ncol) * mw_facm(kp1,:ncol)
       wrk0(:ncol) = ( kv(kp1,:ncol) + kmq(kp1,:ncol) ) * tmpi(kp1,:ncol)
!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
!         wrk1(:ncol) = .5_r8 * (kmq(kp1,:ncol) * rghd(kp1,:ncol) &
!                            - (kv(kp1,:ncol) + kmq(kp1,:ncol)) * gradm(kp1,:ncol) &
!                            - kmq(kp1,:ncol) * gradt(kp1,:ncol))
!       else
         wrk1(:ncol) = kmq(kp1,:ncol) * 0.5_r8 * rghd(kp1,:ncol)
!       endif
!<----------------------------------Lixh close this part-----------------------------------
     ! Add species separation term
       ca(k,:ncol  )  = ( wrk0(:ncol) - wrk1(:ncol) ) * rpdel(k,:ncol)
       cc(kp1,:ncol)  = ( wrk0(:ncol) + wrk1(:ncol) ) * rpdel(kp1,:ncol)
    end do
    
    if( fixed_ubc ) then
!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
!          kmq(ntop_molec,:ncol)  = kq_scal(ntop_molec,:ncol) * mw_facm(ntop_molec,:ncol)
!          wrk0(:ncol) = (kv(ntop_molec,:ncol) + kmq(ntop_molec,:ncol)) * tmpi(ntop_molec,:ncol)/2._r8
!          ! /2. is to extrapolate/(pmid(1)-pint(1)) to /(pmid(1)-pmid(0))
!          wrk1(:ncol) = .5_r8 * (kmq(ntop_molec,:ncol) * rghd(ntop_molec,:ncol) &
!                      - (kv(ntop_molec,:ncol) + kmq(ntop_molec,:ncol)) * gradm(ntop_molec,:ncol) &
!                      - kmq(ntop_molec,:ncol) * gradt(ntop_molec,:ncol))
!          cc(ntop_molec,:ncol)  = (wrk0(:ncol) + wrk1(:ncol)) * rpdel(ntop_molec,:ncol)
!       else
!          cc(ntop_molec,:ncol) = kq_scal(ntop_molec,:ncol) * mw_facm(ntop_molec,:ncol) &
!                               * ( tmpi(ntop_molec,:ncol) + rghd(ntop_molec,:ncol) )   &
!                               * rpdel(ntop_molec,:ncol)
!       endif
!<----------------------------------Lixh close this part-----------------------------------
    end if

  ! Calculate diagonal elements

    do k = nbot_molec - 1, ntop_molec + 1, -1
       kp1 = k + 1
!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
!          cb(k,:ncol) = 1._r8 + ca(k,:ncol) + cc(k,:ncol) &
!                      + rpdel(k,:ncol) * (kmq(kp1,:ncol)*rghd(kp1,:ncol) &
!                      - kmq(k,:ncol)*rghd(k,:ncol) &
!                      -(kv(kp1,:ncol)+kmq(kp1,:ncol)) * gradm(kp1,:ncol)  &
!                      +(kv(k,:ncol)+kmq(k,:ncol)) * gradm(k,:ncol) &
!                      -kmq(kp1,:ncol) *gradt(kp1,:ncol) &
!                      +kmq(k,:ncol) *gradt(k,:ncol))
!       else
          cb(k,:ncol) = 1._r8 + ca(k,:ncol) + cc(k,:ncol)                   &
                      + rpdel(k,:ncol) * ( kmq(kp1,:ncol) * rghd(kp1,:ncol) &
                      - kmq(k,:ncol) * rghd(k,:ncol) )
!       endif
!<----------------------------------Lixh close this part-----------------------------------
    end do

    k   = ntop_molec
    kp1 = k + 1
!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then     
!       if( fixed_ubc ) then
!          cb(k,:ncol) = 1._r8 + ca(k,:ncol) + cc(k,:ncol) &
!                      + rpdel(k,:ncol) * (kmq(kp1,:ncol)*rghd(kp1,:ncol) &
!                      - kmq(k,:ncol)*rghd(k,:ncol) &
!                      -(kv(kp1,:ncol)+kmq(kp1,:ncol)) * gradm(kp1,:ncol)  &
!                      +(kv(k,:ncol)+kmq(k,:ncol)) * gradm(k,:ncol)  &
!                      -kmq(kp1,:ncol) * gradt(kp1,:ncol)  &
!                      +kmq(k,:ncol) * gradt(k,:ncol))
!       else
!          cb(k,:ncol) = 1._r8 + ca(k,:ncol) &
!                      + rpdel(k,:ncol) * (kmq(kp1,:ncol)*rghd(kp1,:ncol) &
!                      - (kv(kp1,:ncol)+kmq(kp1,:ncol)) * gradm(kp1,:ncol) &
!                      -kmq(kp1,:ncol) * gradt(kp1,:ncol))
!       end if
!    else
       if( fixed_ubc ) then
          cb(k,:ncol) = 1._r8 + ca(k,:ncol)                                   &
                      + rpdel(k,:ncol) * kmq(kp1,:ncol) * rghd(kp1,:ncol)     &
                      + kq_scal(ntop_molec,:ncol) * mw_facm(ntop_molec,:ncol) &
                      * ( tmpi(ntop_molec,:ncol) - rghd(ntop_molec,:ncol) )   &
                      * rpdel(ntop_molec,:ncol)
       else
          cb(k,:ncol) = 1._r8 + ca(k,:ncol) &
                      + rpdel(k,:ncol) * kmq(kp1,:ncol) * rghd(kp1,:ncol)
       end if
!    endif
!<----------------------------------Lixh close this part-----------------------------------

    k   = nbot_molec
    cb(k,:ncol) = 1._r8 + cc(k,:ncol) + ca(k,:ncol) &
                - rpdel(k,:ncol) * kmq(k,:ncol)*rghd(k,:ncol)

  ! Compute term for updating top level mixing ratio for ubc
    if( fixed_ubc ) then
        cd_top(:ncol) = cc(ntop_molec,:ncol) * ubc_mmr(:ncol)
    end if

    !-------------------------------------------------------- !
    ! Calculate e(k).                                         !
    ! This term is required in solution of tridiagonal matrix ! 
    ! defined by implicit diffusion equation.                 !
    !-------------------------------------------------------- !

    do k = nbot_molec, ntop_molec + 1, -1
       dnom(k,:ncol) = 1._r8 / ( cb(k,:ncol) - ca(k,:ncol) * ze(k+1,:ncol) )
       ze(k,:ncol)   = cc(k,:ncol) * dnom(k,:ncol)
    end do
    k = ntop_molec
    dnom(k,:ncol) = 1._r8 / ( cb(k,:ncol) - ca(k,:ncol) * ze(k+1,:ncol) )

    vd_lu_qdecomp = 1
    return

  end function vd_lu_qdecomp


 end module grist_molec_diff
