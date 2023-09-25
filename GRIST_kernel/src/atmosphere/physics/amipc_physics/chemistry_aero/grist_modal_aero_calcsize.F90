module grist_modal_aero_calcsize
    use grist_constants,                    only: r8, pi
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_nml_module,                   only: nlev
    use grist_rad_constituents,             only: rad_cnst_get_info,        &
                                                  rad_cnst_get_mode_props,  &
                                                  rad_cnst_get_mode_num,    &
                                                  rad_cnst_get_aer_mmr,     &
                                                  rad_cnst_get_aer_props
    use cloud_fraction,                     only: top_lev => clim_modal_aero_top_lev
    use grist_handle_error,                 only: endrun
    use grist_mpi

implicit none
private
save


public ::   modal_aero_calcsize_diag
 
contains


    !-----------------------------------------------------------------------
    !
    ! Calculate aerosol size distribution parameters 
    !
    ! ***N.B.*** DGNUM for the modes in the climate list are put directly into
    !            the physics buffer.  For diagnostic list calculations use the
    !            optional list_idx and dgnum args.
    !-----------------------------------------------------------------------
    subroutine modal_aero_calcsize_diag(ncol, list_idx_in, dgnum_m)

    ! arguments
    integer,                         intent(in)     :: ncol 
    integer,  optional,              intent(in)     :: list_idx_in    ! diagnostic list index
    real(r8), optional, allocatable, intent(out)    :: dgnum_m(:,:,:) ! interstital aerosol dry number mode radius (m)

    ! local
    integer  :: i, k, l1, n
    integer  :: list_idx
    integer  :: nmodes
    integer  :: nspec

    real(r8), parameter :: third = 1.0_r8/3.0_r8

    real(r8) :: mode_num(nlev,ncol) ! mode number mixing ratio
    real(r8) :: specmmr(nlev,ncol)  ! specie mmr
    real(r8) :: specdens      ! specie density

    real(r8) :: dryvol_a(nlev, ncol)   ! interstital aerosol dry volume (cm^3/mol_air)
    real(r8) :: dgncur_a(nlev, ncol) 

    real(r8) :: dgnum, dgnumhi, dgnumlo
    real(r8) :: dgnyy, dgnxx           ! dgnumlo/hi of current mode
    real(r8) :: drv_a                  ! dry volume (cm3/mol_air)
    real(r8) :: dumfac, dummwdens      ! work variables
    real(r8) :: num_a0                 ! initial number (#/mol_air)
    real(r8) :: num_a                  ! final number (#/mol_air)
    real(r8) :: voltonumbhi, voltonumblo
    real(r8) :: v2nyy, v2nxx           ! voltonumblo/hi of current mode
    real(r8) :: sigmag, alnsg
 
    list_idx = 0  ! climate list by default
    if (present(list_idx_in)) list_idx = list_idx_in

    call rad_cnst_get_info(list_idx, nmodes=nmodes)

    if (list_idx /= 0) then
       if (.not. present(dgnum_m)) then
          call endrun('modal_aero_calcsize_diag called for'// &
                      'diagnostic list but dgnum_m pointer not present')
       end if
       allocate(dgnum_m(nmodes, nlev, ncol))
    end if

    do n = 1, nmodes

       ! get mode properties
       call rad_cnst_get_mode_props(list_idx, n, dgnum=dgnum, dgnumhi=dgnumhi, dgnumlo=dgnumlo, &
                                    sigmag=sigmag)

       ! get mode number mixing ratio
       call rad_cnst_get_mode_num(list_idx, n, 'a',  mode_num)

       dgncur_a(:,:) = dgnum
       dryvol_a(:,:) = 0.0_r8

       ! compute dry volume mixrats = 
       !      sum_over_components{ component_mass mixrat / density }
       call rad_cnst_get_info(list_idx, n, nspec=nspec)
       do l1 = 1, nspec

          call rad_cnst_get_aer_mmr(list_idx, n, l1, 'a', specmmr)
          call rad_cnst_get_aer_props(list_idx, n, l1, density_aer=specdens)

          ! need qmass*dummwdens = (kg/kg-air) * [1/(kg/m3)] = m3/kg-air
          dummwdens = 1.0_r8 / specdens

          do i=1,ncol
             do k=top_lev,nlev
                dryvol_a(k,i) = dryvol_a(k,i)    &
                   + max(0.0_r8, specmmr(k,i))*dummwdens
             end do
          end do
       end do

       alnsg  = log( sigmag )
       dumfac = exp(4.5_r8*alnsg**2)*pi/6.0_r8
       voltonumblo = 1._r8 / ( (pi/6._r8)*(dgnumlo**3)*exp(4.5_r8*alnsg**2) )
       voltonumbhi = 1._r8 / ( (pi/6._r8)*(dgnumhi**3)*exp(4.5_r8*alnsg**2) )
       v2nxx = voltonumbhi
       v2nyy = voltonumblo
       dgnxx = dgnumhi
       dgnyy = dgnumlo

       do i = 1, ncol
          do k = top_lev, nlev

             drv_a = dryvol_a(k,i)
             num_a0 = mode_num(k,i)
             num_a = max( 0.0_r8, num_a0 )

             if (drv_a > 0.0_r8) then
                if (num_a <= drv_a*v2nxx) then
                   dgncur_a(k,i) = dgnxx
                else if (num_a >= drv_a*v2nyy) then
                   dgncur_a(k,i) = dgnyy
                else
                   dgncur_a(k,i) = (drv_a/(dumfac*num_a))**third
                end if
             end if

          end do
       end do

       if (list_idx == 0) then
          pstate_cam%aerosol_dgnum%f(n,:,1:ncol) = dgncur_a
       else
          dgnum_m(n,:,1:ncol) = dgncur_a
       end if
 
    end do ! nmodes

    end subroutine modal_aero_calcsize_diag


end module grist_modal_aero_calcsize
