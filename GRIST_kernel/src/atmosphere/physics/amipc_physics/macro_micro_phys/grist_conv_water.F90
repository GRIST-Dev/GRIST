!===================================================================================================
!
!  Created by LiXiaohan on 19/07/20, adopted from CAM5
!
!  Purpose:
!  Computes grid-box average liquid (and ice) from stratus and cumulus
!  Just for the purposes of radiation
!
!  Method:
!  Extract information about deep+shallow liquid and cloud fraction from
!  the physics buffer
!
!===================================================================================================

  module conv_water

    use grist_handle_error,                 only: endrun
    use grist_constants,                    only: r8,      i4,            &
                                                  gravity, latvap,  latice
 
    use grist_nml_module,                   only: nlev


    implicit none
    private
    save

    public :: conv_water_4rad

  contains

    subroutine conv_water_4rad( ncol, conv_water_mode, &
                                rei, pdel, ls_liq, ls_ice, totg_liq, totg_ice )

    use phys_control,                       only: phys_getopts
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_physics_update,               only: old_time_level

    implicit none
! io
    integer,  intent(in) :: ncol
    integer,  intent(in) :: conv_water_mode
    real(r8), intent(in) :: rei(nlev, ncol)        ! Ice effective drop size (microns)
    real(r8), intent(in) :: pdel(nlev, ncol)       ! Moist pressure difference across layer
    real(r8), intent(in) :: ls_liq(nlev, ncol)     ! Large-scale contributions to GBA cloud liq      
    real(r8), intent(in) :: ls_ice(nlev, ncol)     ! Large-scale contributions to GBA cloud ice 

    real(r8), intent(out):: totg_ice(nlev, ncol)   ! Total GBA in-cloud ice
    real(r8), intent(out):: totg_liq(nlev, ncol)   ! Total GBA in-cloud liquid
! local
    real(r8) :: conv_ice(nlev, ncol)               ! Convective contributions to IC cloud ice
    real(r8) :: conv_liq(nlev, ncol)               ! Convective contributions to IC cloud liquid
    real(r8) :: tot_ice(nlev, ncol)                ! Total IC ice
    real(r8) :: tot_liq(nlev, ncol)                ! Total IC liquid
 
    integer  :: i,k                                ! Lon, lev indices buff stuff.
    real(r8) :: cu_icwmr                           ! Convective  water for this grid-box.   
    real(r8) :: ls_icwmr                           ! Large-scale water for this grid-box. 
    real(r8) :: tot_icwmr                          ! Large-scale water for this grid-box.  
    real(r8) :: ls_frac                            ! Large-scale cloud frac for this grid-box. 
    real(r8) :: tot0_frac, cu0_frac, dp0_frac, sh0_frac 
    real(r8) :: kabs, kabsi, kabsl, alpha, dp0, sh0, ic_limit, frac_limit  
    real(r8) :: wrk1         

    ! these calls were already done in convect_shallow...so here I add the same fields to the physics buffer with a "1" at the end
    real(r8) :: sh_cldliq1(nlev, ncol)             ! 'SH_CLDLIQ1' shallow gbm cloud liquid water (kg/kg) 
    real(r8) :: sh_cldice1(nlev, ncol)             ! 'SH_CLDICE1' shallow gbm cloud ice water (kg/kg) 

    parameter( kabsl = 0.090361_r8, frac_limit = 0.01_r8, ic_limit = 1.e-12_r8 )
 
  ! Get microphysics option
 
    character(len=16) :: microp_scheme 

    call phys_getopts( microp_scheme_out = microp_scheme )

    ! --------------------------------------------------------------- !
    ! Loop through grid-boxes and determine:                          !
    ! 1. Effective mean in-cloud convective ice/liquid (deep+shallow) !
    ! 2. Effective mean in-cloud total ice/liquid (ls+convective)     !
    ! --------------------------------------------------------------- !

    do k = 1, nlev
    do i = 1, ncol

       if( pstate_cam%cld_sh_frac_at_pc_full_level%f(k,i) <= frac_limit  &
         .or. pstate_cam%sh_icwmr_at_pc_full_level%f(k,i) <= ic_limit ) then
           sh0_frac = 0._r8
       else
           sh0_frac = pstate_cam%cld_sh_frac_at_pc_full_level%f(k,i)
       endif
       if( pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i) <= frac_limit  &
         .or. pstate_cam%dp_icwmr_at_pc_full_level%f(k,i) <= ic_limit ) then
           dp0_frac = 0._r8
       else
           dp0_frac = pstate_cam%cld_dp_frac_at_pc_full_level%f(k,i)
       endif
       cu0_frac = sh0_frac + dp0_frac

    ! For the moment calculate the emissivity based upon the ls clouds ice fraction

      wrk1 = min(1._r8,max(0._r8, ls_ice(k,i)/(ls_ice(k,i)+ls_liq(k,i)+1.e-36_r8)))

      if( ( cu0_frac < frac_limit )     &
        .or. ( ( pstate_cam%sh_icwmr_at_pc_full_level%f(k,i) + pstate_cam%dp_icwmr_at_pc_full_level%f(k,i) ) < ic_limit ) ) then

            cu0_frac = 0._r8
            cu_icwmr = 0._r8
         
            ls_frac = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,k,i)
            if( ls_frac < frac_limit ) then
                ls_frac  = 0._r8
                ls_icwmr = 0._r8
            else
                ls_icwmr = ( ls_liq(k,i) + ls_ice(k,i) )/max(frac_limit,ls_frac) ! Convert to IC value.
            end if

            tot0_frac = ls_frac
            tot_icwmr = ls_icwmr
           
      else

          ! Select radiation constants (effective radii) for emissivity averaging.
            
            if( microp_scheme == 'RK' ) then
               kabsi = 0.005_r8 + 1._r8/rei(k,i)
            else
               kabsi = 0.005_r8 + 1._r8/min(max(13._r8,rei(k,i)),130._r8)
            endif
            kabs  = kabsl * ( 1._r8 - wrk1 ) + kabsi * wrk1
            alpha = -1.66_r8*kabs*pdel(k,i)/gravity*1000.0_r8

          ! Selecting cumulus in-cloud water.            

            select case (conv_water_mode) ! Type of average
            case (1) ! Area weighted arithmetic average
               cu_icwmr = ( sh0_frac * pstate_cam%sh_icwmr_at_pc_full_level%f(k,i)  &
                          + dp0_frac*pstate_cam%dp_icwmr_at_pc_full_level%f(k,i))/max(frac_limit,cu0_frac)
            case (2)
               sh0 = exp(alpha*pstate_cam%sh_icwmr_at_pc_full_level%f(k,i))
               dp0 = exp(alpha*pstate_cam%dp_icwmr_at_pc_full_level%f(k,i))               
               cu_icwmr = log((sh0_frac*sh0+dp0_frac*dp0)/max(frac_limit,cu0_frac))
               cu_icwmr = cu_icwmr/alpha
            case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
            end select

          ! Selecting total in-cloud water. 
          ! Attribute large-scale/convective area fraction differently from default.

            ls_frac   = pstate_cam%macrop_ast_at_pc_full_level%f(old_time_level,k,i) 
            ls_icwmr  = (ls_liq(k,i) + ls_ice(k,i))/max(frac_limit,ls_frac) ! Convert to IC value.
            tot0_frac = (ls_frac + cu0_frac) 

            select case (conv_water_mode) ! Type of average
            case (1) ! Area weighted 'arithmetic in emissivity' average
               tot_icwmr = (ls_frac*ls_icwmr + cu0_frac*cu_icwmr)/max(frac_limit,tot0_frac)
            case (2)
               tot_icwmr = log((ls_frac*exp(alpha*ls_icwmr)+cu0_frac*exp(alpha*cu_icwmr))/max(frac_limit,tot0_frac))
               tot_icwmr = tot_icwmr/alpha
            case default ! Area weighted 'arithmetic in emissivity' average.
!               call endrun ('CONV_WATER_4_RAD: Unknown option for conv_water_in_rad - exiting')
            end select

      end if

    ! Repartition convective cloud water into liquid and ice phase.
    ! Currently, this partition is made using the ice fraction of stratus condensate.
    ! In future, we should use ice fraction explicitly computed from the convection scheme.

      conv_ice(k,i) = cu_icwmr * wrk1
      conv_liq(k,i) = cu_icwmr * (1._r8-wrk1)

      tot_ice(k,i)  = tot_icwmr * wrk1
      tot_liq(k,i)  = tot_icwmr * (1._r8-wrk1)

      totg_ice(k,i) = tot0_frac * tot_icwmr * wrk1
      totg_liq(k,i) = tot0_frac * tot_icwmr * (1._r8-wrk1)

    end do
    end do

    !Note: if use double_plume scheme, sh_icwmr and cld_sh includ both shallow and deep convections, LiXH.
    !sh_cldliq1(:nlev,:ncol) = pstate_cam%sh_icwmr_at_pc_full_level%f(:nlev,:ncol)      &
    !                         *(1-pstate_cam%macrop_fice_at_pc_full_level%f(:nlev,:ncol))*pstate_cam%cld_sh_frac_at_pc_full_level%f(:nlev,:ncol)
    !sh_cldice1(:nlev,:ncol) = pstate_cam%sh_icwmr_at_pc_full_level%f(:nlev,:ncol)      &
    !                         *pstate_cam%macrop_fice_at_pc_full_level%f(:nlev,:ncol)*pstate_cam%cld_sh_frac_at_pc_full_level%f(:nlev,:ncol)

    ! output:
    ! conv_liq, conv_ice, sh_icwmr, dp_icwmr, tot_liq, tot_ice, sh_frac, dp_frac

    end subroutine conv_water_4rad

  end module conv_water
