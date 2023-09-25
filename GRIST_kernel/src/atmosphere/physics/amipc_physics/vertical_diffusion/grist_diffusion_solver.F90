!===========================================================================
!  Created by LiXiaohan on 19/06/04, adopted from CAM5
!
!  Module to solve vertical diffusion equations using a tri-diagonal solver
!  The module will also apply countergradient fluxes, and apply molecular 
!  diffusion for constituents.
!
!===========================================================================

module grist_diffusion_solver

    use grist_constants,                    only: r8, i4, latvap, karman
    use grist_mpi,                          only: mpi_rank
    use grist_handle_error,                 only: endrun

    implicit none
    private
    public :: init_vdiff,           &
              compute_vdiff,        &
              new_fieldlist_vdiff,  &
              vdiff_selector,       &
              vdiff_select,         &
              any

  ! logical array of fields to be diffused
    type vdiff_selector
        private
        logical, allocatable, dimension(:) :: fields
    end type vdiff_selector

    ! Below provides functionality of intrinsic any for type vdiff_selector
    interface any
        module procedure my_any
    end interface

    real(r8)   :: cpair                           ! Specific heat of dry air
    real(r8)   :: gravit                          ! Acceleration due to gravity
    real(r8)   :: rair                            ! Gas constant for dry air
    real(r8)   :: zvir                            ! rh2o/rair - 1
    logical    :: do_iss                          ! Use implicit turbulent surface stress computation

  ! Parameters used for Turbulent Mountain Stress

    real(r8), parameter :: z0fac   = 0.025_r8     ! Factor determining z_0 from orographic standard deviation
    real(r8), parameter :: z0max   = 100._r8      ! Max value of z_0 for orography
    real(r8), parameter :: horomin = 10._r8       ! Min value of subgrid orographic height for mountain stress
    real(r8), parameter :: dv2min  = 0.01_r8      ! Minimum shear squared


contains

    subroutine init_vdiff( rair_in, gravit_in, do_iss_in )

    real(r8),             intent(in)  :: rair_in         ! Input gas constant for dry air
    real(r8),             intent(in)  :: gravit_in       ! Input gravitational acceleration
    logical,              intent(in)  :: do_iss_in       ! Input ISS flag
    
    rair   = rair_in     
    gravit = gravit_in 
    do_iss = do_iss_in

  end subroutine init_vdiff

  type(vdiff_selector) pure function new_fieldlist_vdiff(ncnst)

    integer,              intent(in)  :: ncnst           ! Number of constituents

    allocate( new_fieldlist_vdiff%fields( 3 + ncnst ) )
    new_fieldlist_vdiff%fields = .false.

  end function new_fieldlist_vdiff


! Purpose: Driver routine to compute vertical diffusion of momentum, moisture, trace
!          constituents and dry static energy. The new temperature is computed from
!          the diffused dry static energy.
!          Turbulent diffusivities and boundary layer nonlocal transport terms are
!          obtained from the turbulence module.
  subroutine compute_vdiff( pver            , ncnst              , ncol          , pmid         ,               &
                            pint            , rpdel              , t             , ztodt        , taux        , &
                            tauy            , shflx              , cflx          , ntop         , nbot        , &
                            kvh             , kvm                , kvq           , cgs          , cgh         , &
                            zi              , ksrftms            , qmincg        , fieldlist    ,               &
                            u               , v                  , q             , dse          ,               &
                            tautmsx         , tautmsy            , dtk           , topflx       ,               &
                            tauresx         , tauresy            , itaures       , cpairv       , rairi       , &
                            do_molec_diff  , compute_molec_diff, vd_lu_qdecomp, kvt )

! Modification : Ideally, we should diffuse 'liquid-ice static energy' (sl), not the dry static energy.
!                Also, vertical diffusion of cloud droplet number concentration and aerosol number
!                concentration should be done very carefully in the future version.

! io
    integer,  intent(in)    :: pver
    integer,  intent(in)    :: ncnst
    integer,  intent(in)    :: ncol                      ! Number of atmospheric columns
    integer,  intent(in)    :: ntop                      ! Top    interface level to which vertical diffusion is applied ( = 1 ).
    integer,  intent(in)    :: nbot                      ! Bottom interface level to which vertical diffusion is applied ( = pver ).
    integer,  intent(in)    :: itaures                   ! Indicator determining whether 'tauresx,tauresy'
                                                         ! is updated (1) or non-updated (0) in this subroutine.
    real(r8), intent(in)    :: pmid(pver, ncol)          ! Mid-point pressures [ Pa ]
    real(r8), intent(in)    :: pint(pver+1, ncol)        ! Interface pressures [ Pa ]
    real(r8), intent(in)    :: rpdel(pver, ncol)         ! 1./pdel
    real(r8), intent(in)    :: t(pver, ncol)             ! Temperature [ K ]
    real(r8), intent(in)    :: ztodt                     ! 2 delta-t [ s ]
    real(r8), intent(in)    :: taux(ncol)                ! Surface zonal      stress.
                                                         ! Input u-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: tauy(ncol)                ! Surface meridional stress.
                                                         ! Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
                                                         ! Input v-momentum per unit time per unit area into the atmosphere [ N/m2 ]
    real(r8), intent(in)    :: shflx(ncol)               ! Surface sensible heat flux [ W/m2 ]
    real(r8), intent(in)    :: cflx(ncnst, ncol)         ! Surface constituent flux [ kg/m2/s ]
    real(r8), intent(in)    :: zi(pver+1, ncol)          ! Interface heights [ m ]
    real(r8), intent(in)    :: ksrftms(ncol)             ! Surface drag coefficient for turbulent mountain stress. > 0. [ kg/s/m2 ]
    real(r8), intent(in)    :: qmincg(ncnst)             ! Minimum constituent mixing ratios from cg fluxes
    real(r8), intent(in)    :: cpairv(pver, ncol)        ! Specific heat at constant pressure
    real(r8), intent(in)    :: rairi(pver+1, ncol)       ! Gas constant on interface levels
    real(r8), intent(in)    :: kvh(pver+1, ncol)         ! Eddy diffusivity for heat [ m2/s ]
    logical,  intent(in)    :: do_molec_diff             ! Flag indicating multiple constituent diffusivities

    type(vdiff_selector), intent(in) :: fieldlist        ! Array of flags selecting which fields to diffuse

    integer,  external, optional :: compute_molec_diff   ! Constituent-independent moleculuar diffusivity routine
    integer,  external, optional :: vd_lu_qdecomp        ! Constituent-dependent moleculuar diffusivity routine
    real(r8), intent(out), optional :: kvt(pver+1, ncol) ! Kinematic molecular conductivity

    real(r8), intent(inout) :: kvm(pver+1, ncol)         ! Eddy viscosity ( Eddy diffusivity for momentum ) [ m2/s ]
    real(r8), intent(inout) :: kvq(pver+1, ncol)         ! Eddy diffusivity for constituents
    real(r8), intent(inout) :: cgs(pver+1, ncol)         ! Counter-gradient star [ cg/flux ]
    real(r8), intent(inout) :: cgh(pver+1, ncol)         ! Counter-gradient term for heat
    real(r8), intent(inout) :: u(pver, ncol)             ! U wind. This input is the 'raw' input wind to
                                                         ! PBL scheme without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: v(pver, ncol)             ! V wind. This input is the 'raw' input wind to PBL scheme
                                                         ! without iterative provisional update. [ m/s ]
    real(r8), intent(inout) :: q(ncnst, pver, ncol)      ! Moisture and trace constituents [ kg/kg, #/kg ? ]
    real(r8), intent(inout) :: dse(pver, ncol)           ! Dry static energy [ J/kg ]
    real(r8), intent(inout) :: tauresx(ncol)             ! Input  : Reserved surface stress at previous time step
    real(r8), intent(inout) :: tauresy(ncol)             ! Output : Reserved surface stress at current  time step

    real(r8), intent(out)   :: dtk(pver, ncol)           ! T tendency from KE dissipation
    real(r8), intent(out)   :: tautmsx(ncol)             ! Implicit zonal      turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: tautmsy(ncol)             ! Implicit meridional turbulent mountain surface stress
                                                         ! [ N/m2 = kg m/s /s/m2 ]
    real(r8), intent(out)   :: topflx(ncol)              ! Molecular heat flux at the top interface


! local
    integer  :: i, k, m, icol                            ! Longitude, level, constituent indices
    integer  :: status                                   ! Status indicator
    integer  :: nbot_molec                               ! Bottom level where molecular diffusivity is applied
    integer  :: ntop_molec                               ! Top level where molecular diffusivity is applied
    logical  :: lqtst(ncol)                              ! Adjust vertical profiles
    logical  :: need_decomp                              ! Whether to compute a new decomposition
    logical  :: cnst_fixed_ubc(ncnst)                    ! Whether upper boundary condition is fixed
    logical  :: cnst_fixed_ubflx(ncnst)                  ! Whether upper boundary flux is a fixed non-zero value

    real(r8) :: tmpm(pver, ncol)                         ! Potential temperature, ze term in tri-diag sol'n
    real(r8) :: ca(pver, ncol)                           ! - Upper diag of tri-diag matrix
    real(r8) :: cc(pver, ncol)                           ! - Lower diag of tri-diag matrix
    real(r8) :: dnom(pver, ncol)                         ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))

    real(r8) :: tmp1(ncol)                               ! Temporary storage
    real(r8) :: tmpi1(pver+1, ncol)                      ! Interface KE dissipation
    real(r8) :: tint(pver+1, ncol)                       ! Interface temperature
    real(r8) :: rhoi(pver+1, ncol)                       ! rho at interfaces
    real(r8) :: tmpi2(pver+1, ncol)                      ! dt*(g*rho)**2/dp at interfaces
    real(r8) :: rrho(ncol)                               ! 1./bottom level density 

    real(r8) :: zero(ncol)                               ! Zero array for surface heat exchange coefficients 
    real(r8) :: tautotx(ncol)                            ! Total surface stress ( zonal )
    real(r8) :: tautoty(ncol)                            ! Total surface stress ( meridional )

    real(r8) :: dinp_u(pver+1, ncol)                     ! Vertical difference at interfaces, input u
    real(r8) :: dinp_v(pver+1, ncol)                     ! Vertical difference at interfaces, input v
    real(r8) :: dout_u                                   ! Vertical difference at interfaces, output u
    real(r8) :: dout_v                                   ! Vertical difference at interfaces, output v
    real(r8) :: dse_top(ncol)                            ! dse on top boundary
    real(r8) :: cc_top(ncol)                             ! Lower diagonal at top interface
    real(r8) :: cd_top(ncol)                             ! 

    real(r8) :: qtm(pver, ncol)                          ! Temporary copy of q
    real(r8) :: kq_scal(pver+1, ncol)                    ! kq_fac*sqrt(T)*m_d/rho for molecular diffusivity
    real(r8) :: cnst_mw(ncnst)                           ! Molecular weight [ kg/kmole ]
    real(r8) :: ubc_mmr(ncnst, ncol)                     ! Upper boundary mixing ratios [ kg/kg ]
    real(r8) :: ubc_flux(ncnst)                          ! Upper boundary flux [ kg/s/m^2 ]

    real(r8) :: ws(ncol)                                 ! Lowest-level wind speed [ m/s ]
    real(r8) :: tau(ncol)                                ! Turbulent surface stress ( not including mountain stress )
    real(r8) :: ksrfturb(ncol)                           ! Surface drag coefficient of 'normal' stress. > 0.
                                                         ! Virtual mass input per unit time per unit area [ kg/s/m2 ]
    real(r8) :: ksrf(ncol)                               ! Surface drag coefficient of 'normal' stress +
                                                         ! Surface drag coefficient of 'tms' stress.  > 0. [ kg/s/m2 ] 
    real(r8) :: usum_in(ncol)                            ! Vertical integral of input u-momentum. Total zonal
                                                         ! momentum per unit area in column  [ sum of u*dp/g = kg m/s m-2 ]
    real(r8) :: vsum_in(ncol)                            ! Vertical integral of input v-momentum. Total meridional
                                                         ! momentum per unit area in column [ sum of v*dp/g = kg m/s m-2 ]
    real(r8) :: usum_mid(ncol)                           ! Vertical integral of u-momentum after adding explicit residual stress
    real(r8) :: vsum_mid(ncol)                           ! Vertical integral of v-momentum after adding explicit residual stress
    real(r8) :: usum_out(ncol)                           ! Vertical integral of u-momentum after doing implicit diffusion
    real(r8) :: vsum_out(ncol)                           ! Vertical integral of v-momentum after doing implicit diffusion
    real(r8) :: tauimpx(ncol)                            ! Actual net stress added at the current step other than mountain stress
    real(r8) :: tauimpy(ncol)                            ! Actual net stress added at the current step other than mountain stress
    real(r8) :: wsmin                                    ! Minimum sfc wind speed for estimating frictional
                                                         ! transfer velocity ksrf. [ m/s ]
    real(r8) :: ksrfmin                                  ! Minimum surface drag coefficient [ kg/s/m^2 ]
    real(r8) :: timeres                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]
    real(r8) :: ramda                                    ! dt/timeres [ no unit ]
    
    real(r8) :: mw_fac_loc(ncnst, pver+1, ncol)          ! Local sqrt(1/M_q + 1/M_d) for this constituent

    ! Variables needed for WACCM-X
    real(r8) :: ttemp(pver, ncol)                        ! temporary temperature array
    real(r8) :: ttemp0(pver, ncol)                       ! temporary temperature array

    wsmin    = 1._r8                                     ! Minimum wind speed for ksrfturb computation        [ m/s ]
    ksrfmin  = 1.e-4_r8                                  ! Minimum surface drag coefficient                   [ kg/s/m^2 ]
    timeres  = 7200._r8                                  ! Relaxation time scale of residual stress ( >= dt ) [ s ]

    if( (diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) .and. .not. diffuse(fieldlist,'s') )then
        call endrun("diffusion_solver.compute_vdiff: must diffuse s if diffusing u or v" )
    end if

    zero(:) = 0._r8

    ! Compute 'rho' and 'dt*(g*rho)^2/dp' at interfaces
    tint(1,:ncol) = t(1,:ncol)
    rhoi(1,:ncol) = pint(1,:ncol) / (rairi(1,:ncol)*tint(1,:ncol))
    do k = 2, pver
       do i = 1, ncol
          tint(k,i)  = 0.5_r8 * ( t(k,i) + t(k-1,i) )
          rhoi(k,i)  = pint(k,i) / (rairi(k,i)*tint(k,i))
          tmpi2(k,i) = ztodt * ( gravit*rhoi(k,i) )**2 / ( pmid(k,i) - pmid(k-1,i) )
       end do
    end do
    tint(pver+1,:ncol) = t(pver,:ncol)
    rhoi(pver+1,:ncol) = pint(pver+1,:ncol) / ( rair*tint(pver+1,:ncol) )

    rrho(:ncol) = rair  * t(pver,:ncol) / pmid(pver,:ncol)
    tmp1(:ncol) = ztodt * gravit * rpdel(pver,:ncol)

    !--------------------------------------- !
    ! Computation of Molecular Diffusivities !
    !--------------------------------------- !

    ! Modification : Why 'kvq' is not changed by molecular diffusion ? 

    if( do_molec_diff ) then

        if( (.not.present(compute_molec_diff)) .or. (.not.present(vd_lu_qdecomp)) .or. (.not.present(kvt)) ) then
            if(mpi_rank()==0)then
                print*,'compute_vdiff: do_molec_diff true but compute_molec_diff or vd_lu_qdecomp or kvt missing'
            end if
            call endrun("grist_diffusion_solver")
        endif

      ! The next subroutine 'compute_molec_diff' :
      !     Modifies : kvh, kvm, tint, rhoi, and tmpi2
      !     Returns  : kq_scal, ubc_mmr, dse_top, cc_top, cd_top, cnst_mw, 
      !                cnst_fixed_ubc , mw_fac , ntop_molec 

        mw_fac_loc(:,:,:) = 0._r8

        !--------------------------------------------------------------------------------------------------------
        ! In Extended WACCM, kvt is calculated rather than kvh. This is because molecular diffusion operates on 
        ! temperature, while eddy diffusion operates on dse.  Also, pass in constituent dependent "constants"
        !--------------------------------------------------------------------------------------------------------

        status = compute_molec_diff( pver          , ncnst           , ncol       , t         , pmid   , pint   ,           &
                                     zi            , ztodt           , kvm        , kvt       , tint   , rhoi   , tmpi2  ,  &
                                     kq_scal       , ubc_mmr         , ubc_flux   , dse_top   , cc_top , cd_top ,           &
                                     cnst_mw       , cnst_fixed_ubc, cnst_fixed_ubflx, mw_fac_loc , ntop_molec, nbot_molec )

    else

        kq_scal(:,:) = 0._r8
        cd_top(:)    = 0._r8
        cc_top(:)    = 0._r8

    endif

    !---------------------------- !
    ! Diffuse Horizontal Momentum !
    !---------------------------- !
    if( diffuse(fieldlist,'u') .or. diffuse(fieldlist,'v') ) then
        ! Compute the vertical upward differences of the input u,v for KE dissipation
        ! at each interface.
        ! Velocity = 0 at surface, so difference at the bottom interface is -u,v(pver)
        ! These 'dinp_u, dinp_v' are computed using the non-diffused input wind.

        do i = 1, ncol
           dinp_u(1,i) = 0._r8
           dinp_v(1,i) = 0._r8
           dinp_u(pver+1,i) = -u(pver,i)
           dinp_v(pver+1,i) = -v(pver,i)
        end do
        do k = 2, pver
           do i = 1, ncol
              dinp_u(k,i) = u(k,i) - u(k-1,i)
              dinp_v(k,i) = v(k,i) - v(k-1,i)
           end do
        end do

       ! -------------------------------------------------------------- !
       ! Do 'Implicit Surface Stress' treatment for numerical stability !
       ! in the lowest model layer.                                     !
       ! -------------------------------------------------------------- !

       if( do_iss ) then

         ! Compute surface drag coefficient for implicit diffusion 
         ! including turbulent turbulent mountain stress. 

           do i = 1, ncol
              ws(i)       = max( sqrt( u(pver,i)**2._r8 + v(pver,i)**2._r8 ), wsmin )
              tau(i)      = sqrt( taux(i)**2._r8 + tauy(i)**2._r8 )
              ksrfturb(i) = max( tau(i) / ws(i), ksrfmin )
           end do              
           ksrf(:ncol) = ksrfturb(:ncol) + ksrftms(:ncol)  ! Do all surface stress ( normal + tms ) implicitly

         ! Vertical integration of input momentum. 
         ! This is total horizontal momentum per unit area [ kg*m/s/m2 ] in each column.
         ! Note (u,v) are the raw input to the PBL scheme, not the
         ! provisionally-marched ones within the iteration loop of the PBL scheme.  

           do i = 1, ncol
              usum_in(i) = 0._r8
              vsum_in(i) = 0._r8
              do k = 1, pver
                 usum_in(i) = usum_in(i) + (1._r8/gravit)*u(k,i)/rpdel(k,i)
                 vsum_in(i) = vsum_in(i) + (1._r8/gravit)*v(k,i)/rpdel(k,i)
              end do
           end do              

         ! Add residual stress of previous time step explicitly into the lowest
         ! model layer with a relaxation time scale of 'timeres'.

           ramda         = ztodt / timeres
           u(pver,:ncol) = u(pver,:ncol) + tmp1(:ncol)*tauresx(:ncol)*ramda
           v(pver,:ncol) = v(pver,:ncol) + tmp1(:ncol)*tauresy(:ncol)*ramda

         ! Vertical integration of momentum after adding explicit residual stress
         ! into the lowest model layer.

           do i = 1, ncol
              usum_mid(i) = 0._r8
              vsum_mid(i) = 0._r8
              do k = 1, pver
                 usum_mid(i) = usum_mid(i) + (1._r8/gravit)*u(k,i)/rpdel(k,i)
                 vsum_mid(i) = vsum_mid(i) + (1._r8/gravit)*v(k,i)/rpdel(k,i)
              end do
           end do              

         ! Debug 
         ! icol = phys_debug_col(lchnk) 
         ! if ( icol > 0 .and. get_nstep() .ge. 1 ) then
         !      tauresx_in = tauresx(icol)
         !      tauresy_in = tauresy(icol)
         !      u_in  = u(icol,pver) - tmp1(icol) * tauresx(icol) * ramda
         !      v_in  = v(icol,pver) - tmp1(icol) * tauresy(icol) * ramda
         !      u_res = u(icol,pver)
         !      v_res = v(icol,pver)
         ! endif
         ! Debug

       else

         ! In this case, do 'turbulent mountain stress' implicitly, 
         ! but do 'normal turbulent stress' explicitly.
         ! In this case, there is no 'redisual stress' as long as 'tms' is
         ! treated in a fully implicit wway, which is true.

         ! 1. Do 'tms' implicitly

           ksrf(:ncol) = ksrftms(:ncol) 

         ! 2. Do 'normal stress' explicitly

           u(pver,:ncol) = u(pver,:ncol) + tmp1(:ncol)*taux(:ncol)
           v(pver,:ncol) = v(pver,:ncol) + tmp1(:ncol)*tauy(:ncol)

       end if  ! End of 'do iss' ( implicit surface stress )

       ! --------------------------------------------------------------------------------------- !
       ! Diffuse horizontal momentum implicitly using tri-diagnonal matrix.                      !
       ! The 'u,v' are input-output: the output 'u,v' are implicitly diffused winds.             !
       !    For implicit 'normal' stress : ksrf = ksrftms + ksrfturb,                            !
       !                                   u(pver) : explicitly include 'redisual normal' stress !
       !    For explicit 'normal' stress : ksrf = ksrftms                                        !
       !                                   u(pver) : explicitly include 'normal' stress          !
       ! Note that in all the two cases above, 'tms' is fully implicitly treated.                !
       ! --------------------------------------------------------------------------------------- !

       call vd_lu_decomp( pver , ncol  ,                        &
                          ksrf  , kvm  , tmpi2 , rpdel , ztodt , zero , &
                          ca    , cc   , dnom  , tmpm  , ntop  , nbot )

       call vd_lu_solve(  pver , ncol  ,                        &
                          u     , ca   , tmpm  , dnom  , ntop  , nbot , zero )

       call vd_lu_solve(  pver , ncol  ,                        &
                          v     , ca   , tmpm  , dnom  , ntop  , nbot , zero )

       ! ---------------------------------------------------------------------- !
       ! Calculate 'total' ( tautotx ) and 'tms' ( tautmsx ) stresses that      !
       ! have been actually added into the atmosphere at the current time step. ! 
       ! Also, update residual stress, if required.                             !
       ! ---------------------------------------------------------------------- !

       do i = 1, ncol

          ! Compute the implicit 'tms' using the updated winds.
          ! Below 'tautmsx(i),tautmsy(i)' are pure implicit mountain stresses
          ! that has been actually added into the atmosphere both for explicit
          ! and implicit approach. 

          tautmsx(i) = -ksrftms(i)*u(pver,i)
          tautmsy(i) = -ksrftms(i)*v(pver,i)

          if( do_iss ) then

            ! Compute vertical integration of final horizontal momentum

              usum_out(i) = 0._r8
              vsum_out(i) = 0._r8
              do k = 1, pver
                 usum_out(i) = usum_out(i) + (1._r8/gravit)*u(k,i)/rpdel(k,i)
                 vsum_out(i) = vsum_out(i) + (1._r8/gravit)*v(k,i)/rpdel(k,i)
              end do

            ! Compute net stress added into the atmosphere at the current time step.
            ! Note that the difference between 'usum_in' and 'usum_out' are induced
            ! by 'explicit residual stress + implicit total stress' for implicit case, while
            ! by 'explicit normal   stress + implicit tms   stress' for explicit case. 
            ! Here, 'tautotx(i)' is net stress added into the air at the current time step.

              tauimpx(i) = ( usum_out(i) - usum_in(i) ) / ztodt
              tauimpy(i) = ( vsum_out(i) - vsum_in(i) ) / ztodt

              tautotx(i) = tauimpx(i) 
              tautoty(i) = tauimpy(i) 

            ! Compute redisual stress and update if required.
            ! Note that the total stress we should have added at the current step is
            ! the sum of 'taux(i) - ksrftms(i)*u(pver,i) + tauresx(i)'.

              if( itaures .eq. 1 ) then
                  tauresx(i) = taux(i) + tautmsx(i) + tauresx(i) - tauimpx(i)
                  tauresy(i) = tauy(i) + tautmsy(i) + tauresy(i) - tauimpy(i)
              endif

          else

              tautotx(i) = tautmsx(i) + taux(i)
              tautoty(i) = tautmsy(i) + tauy(i)
              tauresx(i) = 0._r8
              tauresy(i) = 0._r8

          end if  ! End of 'do_iss' if

       end do ! End of 'do i = 1, ncol' loop

       k = pver + 1
       do i = 1, ncol
          tmpi1(1,i) = 0._r8
          tmpi1(k,i) = 0.5_r8 * ztodt * gravit * &
                       ( (-u(k-1,i) + dinp_u(k,i))*tautotx(i) + (-v(k-1,i) + dinp_v(k,i))*tautoty(i) )
       end do

       do k = 2, pver
          do i = 1, ncol
             dout_u = u(k,i) - u(k-1,i)
             dout_v = v(k,i) - v(k-1,i)
             tmpi1(k,i) = 0.25_r8 * tmpi2(k,i) * kvm(k,i) * &
                          ( dout_u**2 + dout_v**2 + dout_u*dinp_u(k,i) + dout_v*dinp_v(k,i) )
          end do
       end do

       ! 2. Compute dissipation term at midpoints, add to dry static energy

       do k = 1, pver
          do i = 1, ncol
             dtk(k,i) = ( tmpi1(k+1,i) + tmpi1(k,i) ) * rpdel(k,i)
             dse(k,i) = dse(k,i) + dtk(k,i)
          end do
       end do

    end if ! End of diffuse horizontal momentum, diffuse(fieldlist,'u') routine

    !-------------------------- !
    ! Diffuse Dry Static Energy !
    !-------------------------- !

  ! Modification : In future, we should diffuse the fully conservative 
  !                moist static energy,not the dry static energy.

    if( diffuse(fieldlist,'s') ) then

      ! Add counter-gradient to input static energy profiles

       do k = 1, pver
           dse(k,:ncol) = dse(k,:ncol) + ztodt * rpdel(k,:ncol) * gravit  *                &
                                       ( rhoi(k+1,:ncol) * kvh(k+1,:ncol) * cgh(k+1,:ncol) &
                                       - rhoi(k,:ncol  ) * kvh(k,:ncol  ) * cgh(k,:ncol  ) )
       end do

     ! Add the explicit surface fluxes to the lowest layer

       dse(pver,:ncol) = dse(pver,:ncol) + tmp1(:ncol) * shflx(:ncol)

     ! Diffuse dry static energy
     !----------------------------------------------------------------------------------------------------
     ! In Extended WACCM, kvt is calculated rather kvh. This is because molecular diffusion operates on 
     ! temperature, while eddy diffusion operates on dse.  Also, pass in constituent dependent "constants"
     !----------------------------------------------------------------------------------------------------

!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
!         cc_top(:) = 0._r8      
!         cd_top(:) = 0._r8
!       endif
!<----------------------------------Lixh close this part-----------------------------------

       call vd_lu_decomp( pver , ncol  ,                         &
                          zero  , kvh  , tmpi2 , rpdel , ztodt , cc_top, &
                          ca    , cc   , dnom  , tmpm  , ntop  , nbot    )

       call vd_lu_solve(  pver , ncol  ,                         &
                          dse   , ca   , tmpm  , dnom  , ntop  , nbot  , cd_top )

     ! Calculate flux at top interface
     
     ! Modification : Why molecular diffusion does not work for dry static energy in all layers ?

       if( do_molec_diff ) then
           topflx(:ncol) =  - kvh(ntop_molec,:ncol) * tmpi2(ntop_molec,:ncol) / (ztodt*gravit) * &
                            ( dse(ntop_molec,:ncol) - dse_top(:ncol) )
       end if

       !---------------------------------------------------
       ! Solve for temperature using thermal conductivity 
       !---------------------------------------------------
!-----------------------------------Lixh close this part---------------------------------->
! if waccmx is avaiable, this part should be modified:
!
!       if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
!         call vd_lu_decomp( pver , ncol  ,                           &
!                            zero  , kvt  , tmpi2 , rpdel , ztodt , zero,     &
!                            ca    , cc   , dnom  , tmpm  , ntop  , nbot, cpairv    )
!         do k = 1,pver
!           ttemp0(k,:ncol) = t(k,:ncol)
!           ttemp(k,:ncol) = ttemp0(k,:ncol)
!         enddo
!         call vd_lu_solve(                                                &
!            pver, ncol     ,                                       &
!            ttemp      ,ca       ,tmpm     ,dnom     ,ntop     ,nbot,zero    )   ! upper boundary is zero flux for extended model
!         !  Update dry static energy
!         do k = 1,pver
!           dse(k,:ncol) = dse(k,:ncol) + cpairv(k,:ncol)*(ttemp(k,:ncol) - ttemp0(k,:ncol))
!         enddo
!       endif
!<----------------------------------Lixh close this part-----------------------------------

    endif

    !---------------------------- !
    ! Diffuse Water Vapor Tracers !
    !---------------------------- !

  ! Modification : For aerosols, I need to use separate treatment 
  !                for aerosol mass and aerosol number. 

    ! Loop through constituents

    need_decomp = .true.

    do m = 1, ncnst

       if( diffuse(fieldlist,'q',m) ) then

           ! Add the nonlocal transport terms to constituents in the PBL.
           ! Check for neg q's in each constituent and put the original vertical
           ! profile back if a neg value is found. A neg value implies that the
           ! quasi-equilibrium conditions assumed for the countergradient term are
           ! strongly violated.
 
           qtm(:pver,:ncol) = q(m,:pver,:ncol)

           do k = 1, pver
              q(m,k,:ncol) = q(m,k,:ncol) + &
                             ztodt * rpdel(k,:ncol) * gravit  * ( cflx(m,:ncol) * rrho(:ncol) ) * &
                           ( rhoi(k+1,:ncol) * kvh(k+1,:ncol) * cgs(k+1,:ncol)                    &
                           - rhoi(k,:ncol  ) * kvh(k,:ncol  ) * cgs(k,:ncol  ) )
           end do

           !------------------------------LiXH Test------------------------------->
           !lqtst(:ncol) = all(q(m,:pver,:ncol) >= qmincg(m), 2)
           !do k = 1, pver
           !   q(m,k,:ncol) = merge( q(m,k,:ncol), qtm(k,:ncol), lqtst(:ncol) )
           !end do

           do i = 1, ncol
              lqtst(i) = all(q(m,1:pver,i) >= qmincg(m))
           end do 
           
           do i = 1, ncol
           do k = 1, pver
              q(m,k,i) = merge( q(m,k,i), qtm(k,i), lqtst(i) )
           end do
           end do
           !<-----------------------------LiXH Test-------------------------------

           ! Add the explicit surface fluxes to the lowest layer

           q(m,pver,:ncol) = q(m,pver,:ncol) + tmp1(:ncol) * cflx(m,:ncol)
 
           ! Diffuse constituents.

           if( need_decomp ) then

               call vd_lu_decomp( pver , ncol  ,                         &
                                  zero  , kvq  , tmpi2 , rpdel , ztodt , zero  , &
                                  ca    , cc   , dnom  , tmpm  , ntop  , nbot )

               if( do_molec_diff ) then

                 ! This is for solving molecular diffusion of minor species, thus, for WACCM-X, bypass O and O2 (major species)
                 ! Major species diffusion is calculated separately.  -Hanli Liu
!-----------------------------------Lixh close this part---------------------------------->
if(mpi_rank()==0)print*,'do_molec_diff is not available for now, LiXH has not completed'
!                   if ( diffuse(fieldlistm,'q',m) ) then                  ! decide if diffuse this constituent
!
!                     status = vd_lu_qdecomp( pver   , ncol   , cnst_fixed_ubc(m) , cnst_mw(m), ubc_mmr(m,:), &
!                                             kvq    , kq_scal, mw_fac_loc(m,:,:) , tmpi2     , rpdel       , &
!                                             ca     , cc     , dnom              , tmpm      , rhoi        , &
!                                             tint   , ztodt  , ntop_molec        , nbot_molec, cd_top      , &
!                                             pmid   , pint   , t                 , m         )
!
!                     ! This to calculate the upper boundary flux of H.    -Hanli Liu
!                     if ((cnst_fixed_ubflx(m))) &
!                            q(m,1,:ncol) = q(m,1,:ncol) - ztodt * rpdel(1,:ncol) * gravit * ubc_flux(m)
!
!                   endif
!<----------------------------------Lixh close this part-----------------------------------

               else
                   need_decomp =  .false.
               endif
           end if

           call vd_lu_solve(  pver , ncol  ,                         &
                              q(m,:,1:ncol) , ca, tmpm  , dnom  , ntop  , nbot  , cd_top )

       end if
    end do

    return
  end subroutine compute_vdiff


! Purpose: Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
!          tridiagonal diffusion matrix.
!          The diagonal elements (1+ca(k)+cc(k)) are not required by the solver.
!          Also determine ze factor and denominator for ze and zf (see solver).
  subroutine vd_lu_decomp( pver, ncol ,                               &
                           ksrf , kv  , tmpi , rpdel, ztodt , cc_top, &
                           ca   , cc  , dnom , ze   , ntop  , nbot, cpairv )
! io
    integer,  intent(in)  :: pver                  ! Number of allocated atmospheric levels 
    integer,  intent(in)  :: ncol                  ! Number of computed atmospheric columns
    integer,  intent(in)  :: ntop                  ! Top level to operate on
    integer,  intent(in)  :: nbot                  ! Bottom level to operate on
    real(r8), intent(in)  :: ksrf(ncol)            ! Surface "drag" coefficient [ kg/s/m2 ]
    real(r8), intent(in)  :: kv(pver+1, ncol)      ! Vertical diffusion coefficients [ m2/s ]
    real(r8), intent(in)  :: tmpi(pver+1, ncol)    ! dt*(g/R)**2/dp*pi(k+1)/(.5*(tm(k+1)+tm(k))**2
    real(r8), intent(in)  :: rpdel(pver, ncol)     ! 1./pdel  (thickness bet interfaces)
    real(r8), intent(in)  :: ztodt                 ! 2 delta-t [ s ]
    real(r8), intent(in)  :: cc_top(ncol)          ! Lower diagonal on top interface (for fixed ubc only)

    real(r8), intent(out) :: ca(pver, ncol)        ! Upper diagonal
    real(r8), intent(out) :: cc(pver, ncol)        ! Lower diagonal
    real(r8), intent(out) :: dnom(pver, ncol)      ! 1./(1. + ca(k) + cc(k) - cc(k)*ze(k-1))
    real(r8), intent(out) :: ze(pver, ncol)        ! Term in tri-diag. matrix system

    real(r8), intent(in), optional  :: cpairv(pver, ncol) ! "Variable" specific heat at constant pressure
! local
    integer :: i                                   ! Longitude index
    integer :: k                                   ! Vertical  index

    ! Determine superdiagonal (ca(k)) and subdiagonal (cc(k)) coeffs of the 
    ! tridiagonal diffusion matrix. The diagonal elements  (cb=1+ca+cc) are
    ! a combination of ca and cc; they are not required by the solver.

    ! If switch present and set true, then input kv is kvt for use in diagonal calculations
    if ( present(cpairv) ) then
      do k = nbot-1, ntop, -1
        do i = 1, ncol
           ca(k,i  ) = kv(k+1,i)*tmpi(k+1,i)*rpdel(k,i  ) / cpairv(k,i)
           cc(k+1,i) = kv(k+1,i)*tmpi(k+1,i)*rpdel(k+1,i) / cpairv(k+1,i)
        end do
      end do
    else
      do k = nbot - 1, ntop, -1
        do i = 1, ncol
          ca(k,i  ) = kv(k+1,i) * tmpi(k+1,i) * rpdel(k,i  )
          cc(k+1,i) = kv(k+1,i) * tmpi(k+1,i) * rpdel(k+1,i)
        end do
      end do
    endif

    ! The bottom element of the upper diagonal (ca) is zero (element not used).
    ! The subdiagonal (cc) is not needed in the solver.

    do i = 1, ncol
       ca(nbot,i) = 0._r8
    end do

    ! Calculate e(k).  This term is 
    ! required in solution of tridiagonal matrix defined by implicit diffusion eqn.

    do i = 1, ncol
       dnom(nbot,i) = 1._r8/(1._r8 + cc(nbot,i) + ksrf(i)*ztodt*gravit*rpdel(nbot,i))
       ze(nbot,i)   = cc(nbot,i)*dnom(nbot,i)
    end do

    do k = nbot - 1, ntop + 1, -1
       do i = 1, ncol
          dnom(k,i) = 1._r8/(1._r8 + ca(k,i) + cc(k,i) - ca(k,i)*ze(k+1,i))
          ze(k,i)   = cc(k,i)*dnom(k,i)
       end do
    end do

    do i = 1, ncol
       dnom(ntop,i) = 1._r8/(1._r8 + ca(ntop,i) + cc_top(i) - ca(ntop,i)*ze(ntop+1,i))
    end do

    return
  end subroutine vd_lu_decomp


! Purpose: Solve the implicit vertical diffusion equation with zero flux boundary conditions.
!          Procedure for solution of the implicit equation follows Richtmyer and Morton (1967,pp 198-200).
!          The equation solved is
!
!          -ca(k)*q(k+1) + cb(k)*q(k) - cc(k)*q(k-1) = d(k),
!
!          where d(k) is the input profile and q(k) is the output profile
!
!          The solution has the form
!
!          q(k) = ze(k)*q(k-1) + zf(k)
!
!          ze(k) = cc(k) * dnom(k)
!
!          zf(k) = [d(k) + ca(k)*zf(k+1)] * dnom(k)
!
!          dnom(k) = 1/[cb(k) - ca(k)*ze(k+1)] =  1/[1 + ca(k) + cc(k) - ca(k)*ze(k+1)]
!
!          Note that the same routine is used for temperature, momentum and tracers,
!          and that input variables are replaced.
  subroutine vd_lu_solve( pver  , ncol ,                                    &
                          q     , ca   , ze   , dnom , ntop , nbot , cd_top )
! io
    integer,  intent(in)    :: pver                   ! Number of allocated atmospheric levels 
    integer,  intent(in)    :: ncol                   ! Number of computed atmospheric columns
    integer,  intent(in)    :: ntop                   ! Top level to operate on
    integer,  intent(in)    :: nbot                   ! Bottom level to operate on
    real(r8), intent(in)    :: ca(pver, ncol)         ! -Upper diag coeff.of tri-diag matrix
    real(r8), intent(in)    :: ze(pver, ncol)         ! Term in tri-diag solution
    real(r8), intent(in)    :: dnom(pver, ncol)       ! 1./(1. + ca(k) + cc(k) - ca(k)*ze(k+1))
    real(r8), intent(in)    :: cd_top(ncol)           ! cc_top * ubc value
    real(r8), intent(inout) :: q(pver, ncol)          ! Constituent field
! local
    real(r8)                :: zf(pver, ncol)         ! Term in tri-diag solution
    integer                    i, k                   ! Longitude, vertical indices

    ! Calculate zf(k). Terms zf(k) and ze(k) are required in solution of 
    ! tridiagonal matrix defined by implicit diffusion equation.
    ! Note that only levels ntop through nbot need be solved for.

    do i = 1, ncol
       zf(nbot,i) = q(nbot,i)*dnom(nbot,i)
    end do

    do k = nbot - 1, ntop + 1, -1
       do i = 1, ncol
          zf(k,i) = (q(k,i) + ca(k,i)*zf(k+1,i))*dnom(k,i)
       end do
    end do

    ! Include boundary condition on top element

    k = ntop
    do i = 1, ncol
       zf(k,i) = (q(k,i) + cd_top(i) + ca(k,i)*zf(k+1,i))*dnom(k,i)
    end do

    ! Perform back substitution

    do i = 1, ncol
       q(ntop,i) = zf(ntop,i)
    end do

    do k = ntop + 1, nbot, +1
       do i = 1, ncol
          q(k,i) = zf(k,i) + ze(k,i)*q(k-1,i)
       end do
    end do

    return
  end subroutine vd_lu_solve


  character(128) function vdiff_select( fieldlist, name, qindex )
    ! --------------------------------------------------------------------- !
    ! This function sets the field with incoming name as one to be diffused !
    ! --------------------------------------------------------------------- !
    type(vdiff_selector), intent(inout)        :: fieldlist
    character(*),         intent(in)           :: name
    integer,              intent(in), optional :: qindex
    
    vdiff_select = ''
    select case (name)
    case ('u','U')
       fieldlist%fields(1) = .true.
    case ('v','V')
       fieldlist%fields(2) = .true.
    case ('s','S')
       fieldlist%fields(3) = .true.
    case ('q','Q')
       if( present(qindex) ) then
           fieldlist%fields(3 + qindex) = .true.
       else
           fieldlist%fields(4) = .true.
       endif
    case default
       write(vdiff_select,*) 'Bad argument to vdiff_index: ', name
    end select
    return
    
  end function vdiff_select


  logical function my_any(a)
      ! -------------------------------------------------- !
      ! This function extends the intrinsic function 'any' !
      ! to operate on type vdiff_selector                  ! 
      ! -------------------------------------------------- !
      type(vdiff_selector), intent(in) :: a
      my_any = any(a%fields)
  end function my_any


  logical function diffuse(fieldlist,name,qindex)
    ! ---------------------------------------------------------------------------- !
    ! This function reports whether the field with incoming name is to be diffused !
    ! ---------------------------------------------------------------------------- !
    type(vdiff_selector), intent(in)           :: fieldlist
    character(*),         intent(in)           :: name
    integer,              intent(in), optional :: qindex
    
    select case (name)
    case ('u','U')
       diffuse = fieldlist%fields(1)
    case ('v','V')
       diffuse = fieldlist%fields(2)
    case ('s','S')
       diffuse = fieldlist%fields(3)
    case ('q','Q')
       if( present(qindex) ) then
           diffuse = fieldlist%fields(3 + qindex)
       else
           diffuse = fieldlist%fields(4)
       endif
    case default
       diffuse = .false.
    end select
    return
  end function diffuse

end module grist_diffusion_solver

