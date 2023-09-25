!===================================================================================================
!
!  Created by LiXiaohan on 19/06/01, adopted from CAM5
!  to provide saturation vapor pressure and related calculations      
!
!  Methods:
!  Pure water/ice saturation vapor pressures are calculated on the fly,
!  with the specific method determined by a runtime option.
!  Mixed phase SVP is interpolated from the internal table, estbl,
!  which is created during initialization.
!===================================================================================================

 module grist_wv_saturation

    use grist_constants,            only: i4, r8,               &
                                          epsilo,               &
                                          latvap, latice,       &
                                          rvap, cp, tmelt,      &
                                          h2otrip

    use grist_handle_error,         only: endrun
                                          

    implicit none
    private
    public          ::  read_nml_wv_sat,&
                        wv_sat_init,    &
                        end_of_wv_sat,  &
                        svp_water,      &
                        svp_ice,        &
                        qsat,           &
                        qsat_water,     &
                        qsat_ice,       &
                        findsp_vc,      &
                        qmmr

! This value is slightly high, but it seems to be the value for the
! steam point of water originally (and most frequently) used in the
! Goff & Gratch scheme.
    real(r8),  parameter    :: tboil    = 373.16_r8

! Table of saturation vapor pressure values (estbl) from tmin to
! tmax+1 Kelvin, in one degree increments.  ttrice defines the
! transition region, estbl contains a combination of ice & water
! values.
! Make these public parameters in case another module wants to see the
! extent of the table.
    real(r8),  parameter    :: tmin     = 127.16_r8
    real(r8),  parameter    :: tmax     = 375.16_r8

    real(r8),  parameter    :: ttrice   = 20.00_r8

    integer(i4)             :: plenest          ! length of estbl
    real(r8), allocatable   :: estbl(:)         ! table values of saturation vapor pressure
    real(r8)                :: omeps            ! 1.0_r8 - epsilo
    real(r8)                :: c3

    real(r8)                :: rh2o
    real(r8)                :: cpair

  ! Set coefficients for polynomial approximation of difference
  ! between saturation vapor press over water and saturation pressure
  ! over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial is
  ! valid in the range -40 < t < 0 (degrees C).
    real(r8) :: pcf(5) = (/ 5.04469588506e-01_r8, &
                           -5.47288442819e+00_r8, &
                           -3.67471858735e-01_r8, &
                           -8.95963532403e-03_r8, &
                           -7.78053686625e-05_r8 /)

!   --- Degree 6 approximation ---
!  real(r8) :: pcf(6) = (/ &
!       7.63285250063e-02, &
!       5.86048427932e+00, &
!       4.38660831780e-01, &
!       1.37898276415e-02, &
!       2.14444472424e-04, &
!       1.36639103771e-06 /)

contains

    !subroutine read_nml_wv_sat(nlfile)
    subroutine read_nml_wv_sat()

    use wv_sat_methods,     only: wv_sat_set_default
    !character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input
    ! Local variables
    integer :: unitn, ierr
    character(len=32) :: wv_sat_scheme = "GoffGratch"

! zhangyi, comment to not read this nml

   ! namelist /wv_sat_nl/ wv_sat_scheme
   ! unitn = 111
   ! open(unitn, file=trim(nlfile), status='old')
   ! read(unitn, wv_sat_nl, iostat=ierr)
   ! if (ierr /= 0) call endrun('error reading wv_sat namelist')
   ! close(unitn)

    if(.not. wv_sat_set_default(wv_sat_scheme))then
        call endrun('read_nml_wv_sat :: Invalid wv_sat_scheme.')
    end if

    end subroutine read_nml_wv_sat

    subroutine wv_sat_init

    use wv_sat_methods,             only: wv_sat_methods_init

    integer(i4)        :: i
    character(len=256) :: errstring     ! For wv_sat_methods error reporting.

    !yizhang add
    call read_nml_wv_sat()

    rh2o    = rvap
    cpair   = cp
    omeps   = 1.0_r8 - epsilo

    !---------------------Lixh is not sure for what--------------------------
    ! Transition range method is only valid for transition temperatures at:
    ! -40 deg C < T < 0 deg C
    !call shr_assert_in_domain(ttrice, ge=0._r8, le=40._r8, varname="ttrice",&
    !                msg="wv_sat_init: Invalid transition temperature range.")
    !---------------------Lixh is not sure for what--------------------------

    ! This parameter uses a hardcoded 287.04_r8?
    c3 = 287.04_r8*(7.5_r8*log(10._r8))/cpair

    call wv_sat_methods_init(tmelt, h2otrip, tboil, ttrice, epsilo, errstring)

    if(trim(errstring) .ne. '')then
        call endrun('wv_sat_methods_init:'//trim(errstring))
    end if

    ! Add two to make the table slightly too big, just in case.
    plenest = ceiling(tmax-tmin) + 2
    allocate(estbl(plenest))

    do i = 1, plenest
        estbl(i) = svp_trans(tmin + real(i-1,r8))
    end do

    end subroutine wv_sat_init


    subroutine end_of_wv_sat
    if (allocated(estbl)) deallocate(estbl)
    end subroutine end_of_wv_sat


! Purpose: Compute saturation vapor pressure over water
    elemental function svp_water(t) result(es)

    use wv_sat_methods, only:  wv_sat_svp_water

    real(r8), intent(in) :: t ! Temperature (K)
    real(r8) :: es            ! SVP (Pa)

    es = wv_sat_svp_water(t)

    end function svp_water


! Purpose: Compute saturation vapor pressure over ice
    elemental function svp_ice(t) result(es)

    use wv_sat_methods, only:  wv_sat_svp_ice

    real(r8), intent(in) :: t ! Temperature (K)
    real(r8) :: es            ! SVP (Pa)

    es = wv_sat_svp_ice(t)

    end function svp_ice


! Purpose: Compute saturation vapor pressure with an ice-water transition
    elemental function svp_trans(t) result(es)

    use wv_sat_methods, only:  wv_sat_svp_trans

    real(r8), intent(in) :: t ! Temperature (K)
    real(r8) :: es            ! SVP (Pa)

    es = wv_sat_svp_trans(t)

    end function svp_trans


! Purpose: Does linear interpolation from nearest values found in the table (estbl)
    elemental function estblf(t) result(es)
! io
    real(r8), intent(in) :: t ! Temperature 
    real(r8) :: es            ! SVP (Pa)
! local
    integer  :: i             ! Index for t in the table
    real(r8) :: t_tmp         ! intermediate temperature for es look-up
    real(r8) :: weight        ! Weight for interpolation

    t_tmp = max(min(t,tmax)-tmin, 0._r8)   ! Number of table entries above tmin
    i = int(t_tmp) + 1                     ! Corresponding index.
    weight = t_tmp - aint(t_tmp, r8)       ! Fractional part of t_tmp (for interpolation).
    es = (1._r8 - weight)*estbl(i) + weight*estbl(i+1)

    end function estblf


! Purpose: Get enthalpy based only on temperature and specific humidity.
    elemental function tq_enthalpy(t, q, hltalt) result(enthalpy)
! io
    real(r8), intent(in) :: t      ! Temperature
    real(r8), intent(in) :: q      ! Specific humidity
    real(r8), intent(in) :: hltalt ! Modified hlat for T derivatives
    real(r8) :: enthalpy

    enthalpy = cpair * t + hltalt * q
  
    end function tq_enthalpy


! Purpose: Calculate latent heat of vaporization of pure liquid water at a given temperature.
    elemental subroutine no_ip_hltalt(t, hltalt)
! io
    real(r8), intent(in)  :: t       ! Temperature
    real(r8), intent(out) :: hltalt  ! Appropriately modified hlat

    hltalt = latvap

    ! Account for change of latvap with t above freezing where
    ! constant slope is given by -2369 j/(kg c) = cpv - cw
    if (t >= tmelt) then
        hltalt = hltalt - 2369.0_r8*(t-tmelt)
    end if

    end subroutine no_ip_hltalt


! Purpose: Calculate latent heat of vaporization of water at a given      
!          temperature, taking into account the ice phase if temperature  
!          is below freezing. Optional argument also calculates a term
!          used to calculate d(es)/dT within the water-ice transition range.
    elemental subroutine calc_hltalt(t, hltalt, tterm)
! io
    real(r8), intent(in) :: t        ! Temperature
    real(r8), intent(out) :: hltalt  ! Appropriately modified hlat
    real(r8), intent(out), optional :: tterm

! local
    real(r8) :: tc      ! Temperature in degrees C
    real(r8) :: weight  ! Weight for es transition from water to ice
    integer :: i

    if (present(tterm)) tterm = 0.0_r8

    call no_ip_hltalt(t,hltalt)
    if (t < tmelt) then
        ! Weighting of hlat accounts for transition from water to ice.
        tc = t - tmelt
        if (tc >= -ttrice) then
            weight = -tc/ttrice
            ! polynomial expression approximates difference between es
            ! over water and es over ice from 0 to -ttrice (C) (max of
            ! ttrice is 40): required for accurate estimate of es
            ! derivative in transition range from ice to water
            if (present(tterm)) then
                do i = size(pcf), 1, -1
                    tterm = pcf(i) + tc*tterm
                end do
                tterm = tterm/ttrice
            end if
        else
            weight = 1.0_r8
        end if
        hltalt = hltalt + weight*latice
    end if

    end subroutine calc_hltalt


! Purpose: Temperature derivative outputs, for qsat_*
    elemental subroutine deriv_outputs(t, p, es, qs, hltalt, tterm, gam, dqsdt)
! io
    real(r8), intent(in) :: t      ! Temperature
    real(r8), intent(in) :: p      ! Pressure
    real(r8), intent(in) :: es     ! Saturation vapor pressure
    real(r8), intent(in) :: qs     ! Saturation specific humidity
    real(r8), intent(in) :: hltalt ! Modified latent heat
    real(r8), intent(in) :: tterm  ! Extra term for d(es)/dT in transition region.
    real(r8), intent(out), optional :: gam      ! (hltalt/cpair)*(d(qs)/dt)
    real(r8), intent(out), optional :: dqsdt    ! (d(qs)/dt)
! local
    real(r8) :: desdt        ! d(es)/dt
    real(r8) :: dqsdt_loc    ! local copy of dqsdt

    if (qs == 1.0_r8) then
        dqsdt_loc = 0._r8
    else
        desdt = hltalt*es/(rh2o*t*t) + tterm
        dqsdt_loc = qs*p*desdt/(es*(p-omeps*es))
    end if

    if (present(dqsdt)) dqsdt = dqsdt_loc
    if (present(gam))   gam   = dqsdt_loc * (hltalt/cpair)

    end subroutine deriv_outputs

!---------------------------------------------------------------------
! QSAT (SPECIFIC HUMIDITY) PROCEDURES
!---------------------------------------------------------------------

! Purpose: Provide saturation mass mixing ratio.
!          Note that qmmr is a ratio over dry air, whereas qsat is a ratio
!          over total mass. These differ mainly at low pressures, where SVP
!          may exceed the actual pressure, and therefore qmmr can blow up.
    elemental subroutine qmmr(t, p, es, qm)
    use wv_sat_methods, only: wv_sat_svp_water, &
                              wv_sat_svp_to_qmmr
! io
    real(r8), intent(in) :: t    ! Temperature
    real(r8), intent(in) :: p    ! Pressure
    real(r8), intent(out) :: es  ! Saturation vapor pressure
    real(r8), intent(out) :: qm  ! Saturation mass mixing ratio
                                 ! (vapor mass over dry mass)
    
    es = wv_sat_svp_water(t)
    
    qm = wv_sat_svp_to_qmmr(es, p)
    
    ! Ensures returned es is consistent with limiters on qmmr.
    es = min(es, p)
    
    end subroutine qmmr


! Purpose: Look up and return saturation vapor pressure from precomputed table,
!          then calculate and return saturation specific humidity. Optionally
!          return various temperature derivatives or enthalpy at saturation.
    elemental subroutine qsat(t, p, es, qs, gam, dqsdt, enthalpy)

    use wv_sat_methods,             only: wv_sat_svp_to_qsat
! io
    real(r8), intent(in)            :: t        ! Temperature
    real(r8), intent(in)            :: p        ! Pressure
    real(r8), intent(out)           :: es       ! Saturation vapor pressure
    real(r8), intent(out)           :: qs       ! Saturation specific humidity
    real(r8), intent(out), optional :: gam      ! (l/cpair)*(d(qs)/dt)
    real(r8), intent(out), optional :: dqsdt    ! (d(qs)/dt)
    real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q
! local
    real(r8) :: hltalt                          ! Modified latent heat for T derivatives
    real(r8) :: tterm                           ! Account for d(es)/dT in transition region

    es = estblf(t)
    qs = wv_sat_svp_to_qsat(es, p)

    ! Ensures returned es is consistent with limiters on qs.
    es = min(es, p)
  
    ! Calculate optional arguments.
    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then
        ! "generalized" analytic expression for t derivative of es
        ! accurate to within 1 percent for 173.16 < t < 373.16
        call calc_hltalt(t, hltalt, tterm)

        if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

        call deriv_outputs(t, p, es, qs, hltalt, tterm, &
                           gam=gam, dqsdt=dqsdt)
    end if
    end subroutine qsat


! Purpose: Calculate SVP over water at a given temperature, and then calculate and
!          return saturation specific humidity. 
!          Optionally return various temperature derivatives or enthalpy at saturation.
    elemental subroutine qsat_water(t, p, es, qs, gam, dqsdt, enthalpy)

    use wv_sat_methods, only: wv_sat_qsat_water
! io
    real(r8), intent(in) :: t    ! Temperature
    real(r8), intent(in) :: p    ! Pressure
    real(r8), intent(out) :: es  ! Saturation vapor pressure
    real(r8), intent(out) :: qs  ! Saturation specific humidity
    real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
    real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
    real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q
! local
    real(r8) :: hltalt       ! Modified latent heat for T derivatives

    call wv_sat_qsat_water(t, p, es, qs)

    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then

        ! "generalized" analytic expression for t derivative of es
        ! accurate to within 1 percent for 173.16 < t < 373.16
        call no_ip_hltalt(t, hltalt)

        if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

        ! For pure water/ice transition term is 0.
        call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
                           gam=gam, dqsdt=dqsdt)

    end if

    end subroutine qsat_water


! Purpose: Calculate SVP over ice at a given temperature, and then calculate and 
!          return saturation specific humidity.
!          Optionally return various temperature derivatives or enthalpy at saturation.
    elemental subroutine qsat_ice(t, p, es, qs, gam, dqsdt, enthalpy)
    use wv_sat_methods, only: wv_sat_qsat_ice
! io
    real(r8), intent(in) :: t    ! Temperature
    real(r8), intent(in) :: p    ! Pressure
    real(r8), intent(out) :: es  ! Saturation vapor pressure
    real(r8), intent(out) :: qs  ! Saturation specific humidity
    real(r8), intent(out), optional :: gam    ! (l/cpair)*(d(qs)/dt)
    real(r8), intent(out), optional :: dqsdt  ! (d(qs)/dt)
    real(r8), intent(out), optional :: enthalpy ! cpair*t + hltalt*q
! local
    real(r8) :: hltalt       ! Modified latent heat for T derivatives

    call wv_sat_qsat_ice(t, p, es, qs)

    if (present(gam) .or. present(dqsdt) .or. present(enthalpy)) then
        ! For pure ice, just add latent heats.
        hltalt = latvap + latice

        if (present(enthalpy)) enthalpy = tq_enthalpy(t, qs, hltalt)

        ! For pure water/ice transition term is 0.
        call deriv_outputs(t, p, es, qs, hltalt, 0._r8, &
                           gam=gam, dqsdt=dqsdt)

    end if

    end subroutine qsat_ice

!---------------------------------------------------------------------
! FINDSP (WET BULB TEMPERATURE) PROCEDURES
!---------------------------------------------------------------------

! Purpose: wrapper for findsp which is 1D and handles the output status.
!          Changing findsp to elemental restricted debugging output.
!          If that output is needed again, it's preferable *not* to copy findsp,
!          but to change the existing version.
    subroutine findsp_vc(q, t, p, use_ice, tsp, qsp)
    use grist_mpi
    ! io
    real(r8), intent(in) :: q(:)        ! water vapor (kg/kg)
    real(r8), intent(in) :: t(:)        ! temperature (K)
    real(r8), intent(in) :: p(:)        ! pressure    (Pa)
    logical,  intent(in) :: use_ice     ! flag to include ice phase in calculations
    real(r8), intent(out) :: tsp(:)     ! saturation temp (K)
    real(r8), intent(out) :: qsp(:)     ! saturation mixing ratio (kg/kg)
    ! local
    integer :: status(size(q))     ! flag representing state of output
                                   ! 0 => Successful convergence
                                   ! 1 => No calculation done: pressure or specific
                                   !      humidity not within usable range
                                   ! 2 => Run failed to converge
                                   ! 4 => Temperature fell below minimum
                                   ! 8 => Enthalpy not conserved
    integer :: n, i

    n = size(q)

    call findsp(q, t, p, use_ice, tsp, qsp, status)

    ! Currently, only 2 and 8 seem to be treated as fatal errors.
    do i = 1,n
        if (status(i) == 2) then
            print*, ' rank = ', mpi_rank()
            print*, ' findsp not converging at i = ', i
            print*, ' t, q, p ', t(i), q(i), p(i)
            print*, ' tsp, qsp ', tsp(i), qsp(i)
            call endrun ('wv_saturation::FINDSP -- not converging')
        else if (status(i) == 8) then
            print*, ' rank = ', mpi_rank()
            print*, ' the enthalpy is not conserved at i = ', i
            print*, ' t, q, p ', t(i), q(i), p(i)
            print*, ' tsp, qsp ', tsp(i), qsp(i)
            call endrun ('wv_saturation::FINDSP -- enthalpy is not conserved')
        endif
    end do

    end subroutine findsp_vc


! Purpose: find the wet bulb temperature for a given t and q in a longitude height section
!          wet bulb temp is the temperature and spec humidity that is just saturated and 
!          has the same enthalpy
!          if q > qs(t) then tsp > t and qsp = qs(tsp) < q
!          if q < qs(t) then tsp < t and qsp = qs(tsp) > q
!
! Method:  a Newton method is used
!          first guess uses an algorithm provided by John Petch from the UKMO
!          we exclude points where the physical situation is unrealistic
!          e.g. where the temperature is outside the range of validity for the
!          saturation vapor pressure, or where the water vapor pressure
!          exceeds the ambient pressure, or the saturation specific humidity is 
!          unrealistic
! 
!          Author: P. Rasch
    elemental subroutine findsp (q, t, p, use_ice, tsp, qsp, status)
! io
    real(r8), intent(in)  :: q        ! water vapor (kg/kg)
    real(r8), intent(in)  :: t        ! temperature (K)
    real(r8), intent(in)  :: p        ! pressure    (Pa)
    logical,  intent(in)  :: use_ice  ! flag to include ice phase in calculations
    real(r8), intent(out) :: tsp      ! saturation temp (K)
    real(r8), intent(out) :: qsp      ! saturation mixing ratio (kg/kg)
    integer,  intent(out) :: status   ! flag representing state of output
                                      ! 0 => Successful convergence
                                      ! 1 => No calculation done: pressure or specific
                                      !      humidity not within usable range
                                      ! 2 => Run failed to converge
                                      ! 4 => Temperature fell below minimum
                                      ! 8 => Enthalpy not conserved
! local
    integer, parameter :: iter = 8    ! max number of times to iterate the calculation
    integer :: l                      ! iterator
    real(r8) es                       ! sat. vapor pressure
    real(r8) gam                      ! change in sat spec. hum. wrt temperature (times hltalt/cpair)
    real(r8) dgdt                     ! work variable
    real(r8) g                        ! work variable
    real(r8) hltalt                   ! lat. heat. of vap.
    real(r8) qs                       ! spec. hum. of water vapor

    real(r8) t1, q1, dt, dq
    real(r8) qvd
    real(r8) r1b, c1, c2
    real(r8), parameter :: dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
    real(r8), parameter :: dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
    real(r8) enin, enout

    ! Saturation specific humidity at this temperature
    if (use_ice) then
        call qsat(t, p, es, qs)
    else
        call qsat_water(t, p, es, qs)
    end if

    ! make sure a meaningful calculation is possible
    if (p <= 5._r8*es .or. qs <= 0._r8 .or. qs >= 0.5_r8 &
       .or. t < tmin .or. t > tmax) then
        status = 1
        ! Keep initial parameters when conditions aren't suitable
        tsp = t
        qsp = q
        enin = 1._r8
        enout = 1._r8

        return
    end if

    ! Prepare to iterate
    status = 2

    ! Get initial enthalpy
    if (use_ice) then
        call calc_hltalt(t,hltalt)
    else
        call no_ip_hltalt(t,hltalt)
    end if
    enin = tq_enthalpy(t, q, hltalt)

    ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
    c1  = hltalt*c3
    c2  = (t + 36._r8)**2
    r1b = c2/(c2 + c1*qs)
    qvd = r1b * (q - qs)
    tsp = t + ((hltalt/cpair)*qvd)

    ! Generate qsp, gam, and enout from tsp.
    if (use_ice) then
        call qsat(tsp, p, es, qsp, gam=gam, enthalpy=enout)
    else
        call qsat_water(tsp, p, es, qsp, gam=gam, enthalpy=enout)
    end if

    ! iterate on first guess
    do l = 1, iter

        g = enin - enout
        dgdt = -cpair * (1 + gam)

        ! New tsp
        t1 = tsp - g/dgdt
        dt = abs(t1 - tsp)/t1
        tsp = t1

        ! bail out if past end of temperature range
        if ( tsp < tmin ) then
            tsp = tmin
            ! Get latent heat and set qsp to a value
            ! that preserves enthalpy.
            if (use_ice) then
                call calc_hltalt(tsp,hltalt)
            else
                call no_ip_hltalt(tsp,hltalt)
            end if
                qsp = (enin - cpair*tsp)/hltalt
                enout = tq_enthalpy(tsp, qsp, hltalt)
            status = 4
            exit
        end if

        ! Re-generate qsp, gam, and enout from new tsp.
        if (use_ice) then
            call qsat(tsp, p, es, q1, gam=gam, enthalpy=enout)
        else
            call qsat_water(tsp, p, es, q1, gam=gam, enthalpy=enout)
        end if
        dq = abs(q1 - qsp)/max(q1,1.e-12_r8)
        qsp = q1

        ! if converged at this point, exclude it from more iterations
        if (dt < dttol .and. dq < dqtol) then
            status = 0
            exit
        endif
    end do

    ! Test for enthalpy conservation
    if (abs((enin-enout)/(enin+enout)) > 1.e-4_r8) status = 8

    end subroutine findsp

 end module grist_wv_saturation

