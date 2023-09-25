!======================================================
!
!  Created by LiXiaohan on 19/8/19.
!
!  parameterizes aerosol coefficients using chebychev polynomial
!  parameterize aerosol radiative properties in terms of
!  surface mode wet radius and wet refractive index
!
!  Ghan and Zaveri, JGR 2007.
!  uses Wiscombe's (1979) mie scattering code
!======================================================

 module modal_aer_opt

    use grist_constants,                    only: r8, rhoh2o, rga=>gravi, rdry
    use grist_nml_module,                   only: nlev, nlevp, ntracer
    use grist_handle_error,                 only: endrun
    use radconstants,                       only: nswbands, nlwbands,                   &
                                                  idx_sw_diag, idx_uv_diag, idx_nir_diag
    use grist_rad_constituents,             only: n_diag, rad_cnst_get_info,            &
                                                  rad_cnst_get_call_list,               &
                                                  rad_cnst_get_aer_mmr,                 &
                                                  rad_cnst_get_aer_props,               &
                                                  rad_cnst_get_mode_props  
    use cloud_fraction,                     only: top_lev => clim_modal_aero_top_lev
    use grist_physics_data_structure,       only: pstate
    use grist_PhysC_data_structure,         only: pstate_cam
    use grist_modal_aero_wateruptake,       only: modal_aero_wateruptake_dr
    use grist_modal_aero_calcsize,          only: modal_aero_calcsize_diag
    use grist_mpi

    implicit none
    private
    save
    
    public :: read_nml_modal_aer_opt,   &
              modal_aer_opt_init,       &
              modal_aero_sw,            &
              modal_aero_lw
    
    character(len=*), parameter :: unset_str = 'UNSET'
    
    ! Namelist variables:
    character(256)      :: modal_optics_file = unset_str   ! full pathname for modal optics dataset
    character(256)      :: water_refindex_file = unset_str ! full pathname for water refractive index dataset
    
    ! Dimension sizes in coefficient arrays used to parameterize aerosol radiative properties
    integer, parameter :: ncoef=5, prefr=7, prefi=10
    
    real(r8) :: xrmin, xrmax
    
    ! refractive index for water read in read_water_refindex
    complex(r8) :: crefwsw(nswbands) ! complex refractive index for water visible
    complex(r8) :: crefwlw(nlwbands) ! complex refractive index for water infrared
    
    character(len=4) :: diag(0:n_diag) = (/'    ','_d1 ','_d2 ','_d3 ','_d4 ','_d5 ', &
                                           '_d6 ','_d7 ','_d8 ','_d9 ','_d10'/)

    real(r8) :: fillvalue = 0._r8

 contains

    subroutine read_nml_modal_aer_opt(nlfile)
! io
    character(len=*), intent(in) :: nlfile
! local
    integer :: unitn, ierr

    namelist /modal_aer_opt_nl/ water_refindex_file

    unitn = 111 
    open( unitn, file=trim(nlfile), status='old' )
    read(unitn, modal_aer_opt_nl, iostat=ierr)
    if (ierr /= 0) call endrun('error reading macro_park namelist')
    close(unitn)
 
    end subroutine read_nml_modal_aer_opt


    subroutine modal_aer_opt_init()

    use grist_physics_iofile,       only: getfile

! local
    integer  :: i, m
    real(r8) :: rmmin, rmmax       ! min, max aerosol surface mode radius treated (m)
    character(len=256) :: locfile
    

    logical :: call_list(0:n_diag)
    integer :: ilist, nmodes, m_ncoef, m_prefr, m_prefi
    integer :: errcode

    character(len=*), parameter :: routine='modal_aer_opt_init'
    !----------------------------------------------------------------------------

    rmmin = 0.01e-6_r8
    rmmax = 25.e-6_r8
    xrmin = log(rmmin)
    xrmax = log(rmmax)

    ! Check that dimension sizes in the coefficient arrays used to
    ! parameterize aerosol radiative properties are consistent between this
    ! module and the mode physprop files.
    call rad_cnst_get_call_list(call_list)
    do ilist = 0, n_diag
       if (call_list(ilist)) then
          call rad_cnst_get_info(ilist, nmodes=nmodes)
          do m = 1, nmodes
             call rad_cnst_get_mode_props(ilist, m, ncoef=m_ncoef, prefr=m_prefr, prefi=m_prefi)
             if (m_ncoef /= ncoef .or. m_prefr /= prefr .or. m_prefi /= prefi) then
                print*, routine//': ERROR - file and module values do not match:'
                print*, '   ncoef:', ncoef, m_ncoef
                print*, '   prefr:', prefr, m_prefr
                print*, '   prefi:', prefi, m_prefi
                call endrun(routine//': ERROR - file and module values do not match')
             end if
          end do
       end if
    end do

    call getfile(water_refindex_file, locfile)
    call read_water_refindex(locfile)
    if(mpi_rank()==0)then
        print*,'modal_aer_opt_init: read water refractive index file:',trim(locfile)
    end if

    end subroutine modal_aer_opt_init


! Purpose : calculates aerosol sw radiative properties
    subroutine modal_aero_sw(ncol, list_idx, nnite, idxnite, &
                             tauxar, wa, ga, fa)
    ! io
    integer,             intent(in) :: ncol
    integer,             intent(in) :: list_idx           ! index of the climate or a diagnostic list
    integer,             intent(in) :: nnite              ! number of night columns
    integer,             intent(in) :: idxnite(nnite)     ! local column indices of night columns
    real(r8), intent(out) :: tauxar(nswbands,0:nlev,ncol) ! layer extinction optical depth
    real(r8), intent(out) :: wa(nswbands,0:nlev,ncol)     ! layer single-scatter albedo
    real(r8), intent(out) :: ga(nswbands,0:nlev,ncol)     ! asymmetry factor
    real(r8), intent(out) :: fa(nswbands,0:nlev,ncol)     ! forward scattered fraction
    ! local
    integer :: i, ifld, isw, k, l, m, nc, ns
    integer :: nmodes
    integer :: nspec

    real(r8) :: mass(nlev,ncol)        ! layer mass
    real(r8) :: air_density(nlev,ncol) ! (kg/m3)

    real(r8)             :: specmmr(nlev,ncol)            ! species mass mixing ratio
    real(r8)             :: specdens                      ! species density (kg/m3)
    complex(r8)          :: specrefindex(nswbands)        ! species refractive index
    character*32         :: spectype                      ! species type
    real(r8)             :: hygro_aer 

    real(r8) :: dgnumwet(nlev,ncol)    ! number mode wet diameter
    real(r8) :: qaerwat(nlev,ncol)     ! aerosol water (g/g)
 
    real(r8), allocatable :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
    real(r8), allocatable :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
    real(r8), allocatable :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
    real(r8), allocatable :: wetdens_m(:,:,:)  ! 

    real(r8) :: sigma_logr_aer         ! geometric standard deviation of number distribution
    real(r8) :: radsurf(nlev,ncol)    ! aerosol surface mode radius
    real(r8) :: logradsurf(nlev,ncol) ! log(aerosol surface mode radius)
    real(r8) :: cheb(ncoef,nlev,ncol)

    real(r8)    :: refr(ncol)         ! real part of refractive index
    real(r8)    :: refi(ncol)         ! imaginary part of refractive index
    complex(r8) :: crefin(ncol)       ! complex refractive index
    real(r8) :: refrtabsw(prefr,nswbands) ! table of real refractive indices for aerosols
    real(r8) :: refitabsw(prefi,nswbands) ! table of imag refractive indices for aerosols
    real(r8) :: extpsw(ncoef,prefr,prefi,nswbands) ! specific extinction
    real(r8) :: abspsw(ncoef,prefr,prefi,nswbands) ! specific absorption
    real(r8) :: asmpsw(ncoef,prefr,prefi,nswbands) ! asymmetry factor
 
    real(r8) :: vol(ncol)             ! volume concentration of aerosol specie (m3/kg)
    real(r8) :: dryvol(ncol)          ! volume concentration of aerosol mode (m3/kg)
    real(r8) :: watervol(ncol)        ! volume concentration of water in each mode (m3/kg)
    real(r8) :: wetvol(ncol)          ! volume concentration of wet mode (m3/kg)

    integer  :: itab(ncol), jtab(ncol)
    real(r8) :: ttab(ncol), utab(ncol)
    real(r8) :: cext(ncoef,ncol), cabs(ncoef,ncol), casm(ncoef,ncol)
    real(r8) :: pext(ncol)            ! parameterized specific extinction (m2/kg)
    real(r8) :: specpext(ncol)        ! specific extinction (m2/kg)
    real(r8) :: dopaer(ncol)          ! aerosol optical depth in layer
    real(r8) :: pabs(ncol)            ! parameterized specific absorption (m2/kg)
    real(r8) :: pasm(ncol)            ! parameterized asymmetry factor
    real(r8) :: palb(ncol)            ! parameterized single scattering albedo

    ! Diagnostics
    real(r8) :: extinct(nlev,ncol)
    real(r8) :: absorb(nlev,ncol)
    real(r8) :: aodvis(ncol)          ! extinction optical depth
    real(r8) :: aodabs(ncol)          ! absorption optical depth

    real(r8) :: aodabsbc(ncol)        ! absorption optical depth of BC

    real(r8) :: ssavis(ncol)
    real(r8) :: dustvol(ncol)         ! volume concentration of dust in aerosol mode (m3/kg)

    real(r8) :: burden(ncol)
    real(r8) :: burdendust(ncol), burdenso4(ncol), burdenbc(ncol), &
                burdenpom(ncol), burdensoa(ncol), burdenseasalt(ncol)

    real(r8) :: aodmode(ncol)
    real(r8) :: dustaodmode(ncol)     ! dust aod in aerosol mode

    real(r8) :: specrefr, specrefi
    real(r8) :: scatdust(ncol), scatso4(ncol), scatbc(ncol), &
                scatpom(ncol), scatsoa(ncol), scatseasalt(ncol)
    real(r8) :: absdust(ncol), absso4(ncol), absbc(ncol), &
                abspom(ncol), abssoa(ncol), absseasalt(ncol)
    real(r8) :: hygrodust(ncol), hygroso4(ncol), hygrobc(ncol), &
                hygropom(ncol), hygrosoa(ncol), hygroseasalt(ncol)

    real(r8) :: scath2o, absh2o, sumscat, sumabs, sumhygro
    real(r8) :: aodc                   ! aod of component

    ! total species AOD
    real(r8) :: dustaod(ncol), so4aod(ncol), bcaod(ncol), &
                pomaod(ncol), soaaod(ncol), seasaltaod(ncol)

    logical :: savaervis               ! true if visible wavelength (0.55 micron)
    logical :: savaernir               ! true if near ir wavelength (~0.88 micron)
    logical :: savaeruv                ! true if uv wavelength (~0.35 micron)

    real(r8) :: aoduv(ncol)            ! extinction optical depth in uv
    real(r8) :: aodnir(ncol)           ! extinction optical depth in nir

    ! debug output
    integer, parameter :: nerrmax_dopaer=1000
    integer  :: nerr_dopaer = 0
    real(r8) :: volf            ! volume fraction of insoluble aerosol
    character(len=*), parameter :: subname = 'modal_aero_sw'
 
    ! initialize output variables
    tauxar(:,:,:ncol) = 0._r8
    wa(:,:,:ncol)     = 0._r8
    ga(:,:,:ncol)     = 0._r8
    fa(:,:,:ncol)     = 0._r8

    ! zero'th layer does not contain aerosol
    tauxar(:,0,1:ncol)  = 0._r8
    wa(:,0,1:ncol)      = 0.925_r8
    ga(:,0,1:ncol)      = 0.850_r8
    fa(:,0,1:ncol)      = 0.7225_r8

    !-----------LiXH Modified dry pressure to pressure---------->
    !mass(:ncol,:)         = state%pdeldry(:ncol,:)*rga
    !air_density(:ncol,:) = state%pmid(:ncol,:)/(rair*state%t(:ncol,:))
    ! air density includes water vapor?
    mass(:,:ncol)         = pstate%delp_at_pc_full_level%f(:,:ncol)*rga
    air_density(:,:ncol)  = pstate%pressure_at_pc_full_level%f(:,:ncol)      &
                          /(rdry*pstate%temp_at_pc_full_level%f(:,:ncol))
    !<----------LiXH Modified dry pressure to pressure-----------

    ! diagnostics for visible band summed over modes
    extinct(:,1:ncol)     = 0.0_r8
    absorb(:,1:ncol)      = 0.0_r8
    aodvis(1:ncol)        = 0.0_r8
    aodabs(1:ncol)        = 0.0_r8
    burdendust(:ncol)     = 0.0_r8
    burdenso4(:ncol)      = 0.0_r8
    burdenpom(:ncol)      = 0.0_r8
    burdensoa(:ncol)      = 0.0_r8
    burdenbc(:ncol)       = 0.0_r8
    burdenseasalt(:ncol)  = 0.0_r8
    ssavis(1:ncol)        = 0.0_r8

    aodabsbc(:ncol)       = 0.0_r8 
    dustaod(:ncol)        = 0.0_r8
    so4aod(:ncol)         = 0.0_r8
    pomaod(:ncol)         = 0.0_r8
    soaaod(:ncol)         = 0.0_r8
    bcaod(:ncol)          = 0.0_r8
    seasaltaod(:ncol)     = 0.0_r8

    ! diags for other bands
    aoduv(:ncol)          = 0.0_r8
    aodnir(:ncol)         = 0.0_r8

    ! loop over all aerosol modes
    call rad_cnst_get_info(list_idx, nmodes=nmodes)
    
    if (list_idx == 0) then
       allocate(dgnumwet_m(nmodes,nlev,ncol))
       allocate(dgnumdry_m(nmodes,nlev,ncol))
       allocate(qaerwat_m(nmodes,nlev,ncol))
       allocate(wetdens_m(nmodes,nlev,ncol))

       ! water uptake and wet radius for the climate list has already been calculated
       dgnumwet_m = pstate_cam%aerosol_dgnumwet%f
       qaerwat_m  = pstate_cam%aerosol_qaerwat%f
    else
       ! If doing a diagnostic calculation then need to calculate the wet radius
       ! and water uptake for the diagnostic modes
       call modal_aero_calcsize_diag(ncol,list_idx, dgnumdry_m)  
       call modal_aero_wateruptake_dr(ncol,list_idx, dgnumdry_m, dgnumwet_m, &
                                      qaerwat_m, wetdens_m)
    endif

    do m = 1, nmodes

       ! diagnostics for visible band for each mode
       burden(:ncol)       = 0._r8
       aodmode(1:ncol)     = 0.0_r8
       dustaodmode(1:ncol) = 0.0_r8

       dgnumwet = dgnumwet_m(m,:,:)
       qaerwat  = qaerwat_m(m,:,:)

       ! get mode properties
       call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtabsw=refrtabsw , &
          refitabsw=refitabsw, extpsw=extpsw, abspsw=abspsw, asmpsw=asmpsw)

       ! get mode info
       call rad_cnst_get_info(list_idx, m, nspec=nspec)

       ! calc size parameter for all columns
       call modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

       do isw = 1, nswbands
          savaervis = (isw .eq. idx_sw_diag)
          savaeruv  = (isw .eq. idx_uv_diag)
          savaernir = (isw .eq. idx_nir_diag)

          do k = top_lev, nlev

             ! form bulk refractive index
             crefin(:ncol) = (0._r8, 0._r8)
             dryvol(:ncol) = 0._r8
             dustvol(:ncol) = 0._r8

             scatdust(:ncol)     = 0._r8
             absdust(:ncol)      = 0._r8
             hygrodust(:ncol)    = 0._r8
             scatso4(:ncol)      = 0._r8
             absso4(:ncol)       = 0._r8
             hygroso4(:ncol)     = 0._r8
             scatbc(:ncol)       = 0._r8
             absbc(:ncol)        = 0._r8
             hygrobc(:ncol)      = 0._r8
             scatpom(:ncol)      = 0._r8
             abspom(:ncol)       = 0._r8
             hygropom(:ncol)     = 0._r8
             scatsoa(:ncol)      = 0._r8
             abssoa(:ncol)       = 0._r8
             hygrosoa(:ncol)     = 0._r8
             scatseasalt(:ncol)  = 0._r8
             absseasalt(:ncol)   = 0._r8
             hygroseasalt(:ncol) = 0._r8

             ! aerosol species loop
             do l = 1, nspec
                call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', specmmr)

                call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                            refindex_aer_sw=specrefindex, spectype=spectype, &
                                            hygro_aer=hygro_aer)
 
                do i = 1, ncol
                   vol(i)      = specmmr(k,i)/specdens
                   dryvol(i)   = dryvol(i) + vol(i)
                   crefin(i)   = crefin(i) + vol(i)*specrefindex(isw)
                end do
 
                ! compute some diagnostics for visible band only
                if (savaervis) then
 
                   specrefr = real(specrefindex(isw))
                   specrefi = aimag(specrefindex(isw))
 
                   do i = 1, ncol
                      burden(i) = burden(i) + specmmr(k,i)*mass(k,i)
                   end do
 
                   if (trim(spectype) == 'dust') then
                      do i = 1, ncol
                         burdendust(i) = burdendust(i) + specmmr(k,i)*mass(k,i)
                         dustvol(i)    = vol(i)
                         scatdust(i)   = vol(i)*specrefr
                         absdust(i)    = -vol(i)*specrefi
                         hygrodust(i)  = vol(i)*hygro_aer
                      end do
                   end if
 
                   if (trim(spectype) == 'sulfate') then
                      do i = 1, ncol
                         burdenso4(i) = burdenso4(i) + specmmr(k,i)*mass(k,i)
                         scatso4(i)   = vol(i)*specrefr
                         absso4(i)    = -vol(i)*specrefi
                         hygroso4(i)  = vol(i)*hygro_aer
                      end do
                   end if
                   if (trim(spectype) == 'black-c') then
                      do i = 1, ncol
                         burdenbc(i) = burdenbc(i) + specmmr(k,i)*mass(k,i)
                         scatbc(i)   = vol(i)*specrefr
                         absbc(i)    = -vol(i)*specrefi
                         hygrobc(i)  = vol(i)*hygro_aer
                    end do
                   end if
                   if (trim(spectype) == 'p-organic') then
                      do i = 1, ncol
                         burdenpom(i) = burdenpom(i) + specmmr(k,i)*mass(k,i)
                         scatpom(i)   = vol(i)*specrefr
                         abspom(i)    = -vol(i)*specrefi
                         hygropom(i)  = vol(i)*hygro_aer
                       end do
                   end if
                   if (trim(spectype) == 's-organic') then
                      do i = 1, ncol
                         burdensoa(i) = burdensoa(i) + specmmr(k,i)*mass(k,i)
                         scatsoa(i)   = vol(i)*specrefr
                         abssoa(i)    = -vol(i)*specrefi
                         hygrosoa(i)  = vol(i)*hygro_aer
                      end do
                   end if
                   if (trim(spectype) == 'seasalt') then
                      do i = 1, ncol
                         burdenseasalt(i) = burdenseasalt(i) + specmmr(k,i)*mass(k,i)
                         scatseasalt(i)   = vol(i)*specrefr
                         absseasalt(i)    = -vol(i)*specrefi
                         hygroseasalt(i)  = vol(i)*hygro_aer
                       end do
                   end if
 
                end if
             end do ! species loop
 
             do i = 1, ncol
                watervol(i) = qaerwat(k,i)/rhoh2o
                wetvol(i) = watervol(i) + dryvol(i)
                if (watervol(i) < 0._r8) then
                   if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                      write(*,'(a,2e10.2,a)') 'watervol,wetvol=', &
                         watervol(i), wetvol(i), ' in '//subname
                   end if
                   watervol(i) = 0._r8
                   wetvol(i) = dryvol(i)
                end if
 
                ! volume mixing
                crefin(i) = crefin(i) + watervol(i)*crefwsw(isw)
                crefin(i) = crefin(i)/max(wetvol(i),1.e-60_r8)
                refr(i)   = real(crefin(i))
                refi(i)   = abs(aimag(crefin(i)))
             end do
 
             ! call t_startf('binterp')
 
             ! interpolate coefficients linear in refractive index
             ! first call calcs itab,jtab,ttab,utab
             itab(:ncol) = 0
             call binterp(extpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                          refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                          itab, jtab, ttab, utab, cext)
             call binterp(abspsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                          refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                          itab, jtab, ttab, utab, cabs)
             call binterp(asmpsw(:,:,:,isw), ncol, ncoef, prefr, prefi, &
                          refr, refi, refrtabsw(:,isw), refitabsw(:,isw), &
                          itab, jtab, ttab, utab, casm)
 
             ! call t_stopf('binterp')
 
             ! parameterized optical properties
             do i=1,ncol
 
                if (logradsurf(k,i) .le. xrmax) then
                   pext(i) = 0.5_r8*cext(1,i)
                   do nc = 2, ncoef
                      pext(i) = pext(i) + cheb(nc,k,i)*cext(nc,i)
                   enddo
                   pext(i) = exp(pext(i))
                else
                   pext(i) = 1.5_r8/(radsurf(k,i)*rhoh2o) ! geometric optics
                endif
 
                ! convert from m2/kg water to m2/kg aerosol
                specpext(i) = pext(i)
                pext(i) = pext(i)*wetvol(i)*rhoh2o
                pabs(i) = 0.5_r8*cabs(1,i)
                pasm(i) = 0.5_r8*casm(1,i)
                do nc = 2, ncoef
                   pabs(i) = pabs(i) + cheb(nc,k,i)*cabs(nc,i)
                   pasm(i) = pasm(i) + cheb(nc,k,i)*casm(nc,i)
                enddo
                pabs(i) = pabs(i)*wetvol(i)*rhoh2o
                pabs(i) = max(0._r8,pabs(i))
                pabs(i) = min(pext(i),pabs(i))
 
                palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
                palb(i) = 1._r8-pabs(i)/max(pext(i),1.e-40_r8)
 
                dopaer(i) = pext(i)*mass(k,i)
             end do
 
             if (savaeruv) then
                do i = 1, ncol
                   aoduv(i) = aoduv(i) + dopaer(i)
                end do
             end if
 
             if (savaernir) then
                do i = 1, ncol
                   aodnir(i) = aodnir(i) + dopaer(i)
                end do
             endif
 
             ! Save aerosol optical depth at longest visible wavelength
             ! sum over layers
             if (savaervis) then
                ! aerosol extinction (/m)
                do i = 1, ncol
                   extinct(k,i) = extinct(k,i) + dopaer(i)*air_density(k,i)/mass(k,i)
                   absorb(k,i)  = absorb(k,i) + pabs(i)*air_density(k,i)
                   aodvis(i)    = aodvis(i) + dopaer(i)
                   aodabs(i)    = aodabs(i) + pabs(i)*mass(k,i)
                   aodmode(i)   = aodmode(i) + dopaer(i)
                   ssavis(i)    = ssavis(i) + dopaer(i)*palb(i)
 
                   if (wetvol(i) > 1.e-40_r8) then
 
                      dustaodmode(i) = dustaodmode(i) + dopaer(i)*dustvol(i)/wetvol(i)
 
                      ! partition optical depth into contributions from each constituent
                      ! assume contribution is proportional to refractive index X volume
 
                      scath2o        = watervol(i)*real(crefwsw(isw))
                      absh2o         = -watervol(i)*aimag(crefwsw(isw))
                      sumscat        = scatso4(i) + scatpom(i) + scatsoa(i) + scatbc(i) + &
                                       scatdust(i) + scatseasalt(i) + scath2o
                      sumabs         = absso4(i) + abspom(i) + abssoa(i) + absbc(i) + &
                                       absdust(i) + absseasalt(i) + absh2o
                      sumhygro       = hygroso4(i) + hygropom(i) + hygrosoa(i) + hygrobc(i) + &
                                       hygrodust(i) + hygroseasalt(i)
 
                      scatdust(i)    = (scatdust(i) + scath2o*hygrodust(i)/sumhygro)/sumscat
                      absdust(i)     = (absdust(i) + absh2o*hygrodust(i)/sumhygro)/sumabs
 
                      scatso4(i)     = (scatso4(i) + scath2o*hygroso4(i)/sumhygro)/sumscat
                      absso4(i)      = (absso4(i) + absh2o*hygroso4(i)/sumhygro)/sumabs
 
                      scatpom(i)     = (scatpom(i) + scath2o*hygropom(i)/sumhygro)/sumscat
                      abspom(i)      = (abspom(i) + absh2o*hygropom(i)/sumhygro)/sumabs
 
                      scatsoa(i)     = (scatsoa(i) + scath2o*hygrosoa(i)/sumhygro)/sumscat
                      abssoa(i)      = (abssoa(i) + absh2o*hygrosoa(i)/sumhygro)/sumabs
 
                      scatbc(i)      = (scatbc(i) + scath2o*hygrobc(i)/sumhygro)/sumscat
                      absbc(i)       = (absbc(i) + absh2o*hygrobc(i)/sumhygro)/sumabs
 
                      scatseasalt(i) = (scatseasalt(i) + scath2o*hygroseasalt(i)/sumhygro)/sumscat
                      absseasalt(i)  = (absseasalt(i) + absh2o*hygroseasalt(i)/sumhygro)/sumabs
                      
                      aodabsbc(i)    = aodabsbc(i) + absbc(i)*dopaer(i)*(1.0_r8-palb(i))
 
                      aodc           = (absdust(i)*(1.0_r8 - palb(i)) + palb(i)*scatdust(i))*dopaer(i)
                      dustaod(i)     = dustaod(i) + aodc
 
                      aodc           = (absso4(i)*(1.0_r8 - palb(i)) + palb(i)*scatso4(i))*dopaer(i)
                      so4aod(i)      = so4aod(i) + aodc
 
                      aodc           = (abspom(i)*(1.0_r8 - palb(i)) + palb(i)*scatpom(i))*dopaer(i)
                      pomaod(i)      = pomaod(i) + aodc
 
                      aodc           = (abssoa(i)*(1.0_r8 - palb(i)) + palb(i)*scatsoa(i))*dopaer(i)
                      soaaod(i)      = soaaod(i) + aodc
 
                      aodc           = (absbc(i)*(1.0_r8 - palb(i)) + palb(i)*scatbc(i))*dopaer(i)
                      bcaod(i)       = bcaod(i) + aodc
 
                      aodc           = (absseasalt(i)*(1.0_r8 - palb(i)) + palb(i)*scatseasalt(i))*dopaer(i)
                      seasaltaod(i)  = seasaltaod(i) + aodc
 
                   endif
 
                end do
             endif
 
             do i = 1, ncol
 
                if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 30._r8)) then
                   print*, 'My Rank=',mpi_rank(),'mass=',mass(k,i) 
                   print*, 'dopaer(', i, ',', k, ',', m, ')=', dopaer(i)
                   ! print*, 'itab,jtab,ttab,utab=',itab(i),jtab(i),ttab(i),utab(i)
                   print*, 'k=', k, ' pext=', pext(i), ' specext=', specpext(i)
                   print*, 'wetvol=', wetvol(i), ' dryvol=', dryvol(i), ' watervol=', watervol(i)
                   print*, 'cext=',(cext(l,i),l=1,ncoef)
                   print*, 'crefin=',crefin(i)
                   !print*, 'nspec=', nspec
                   print*, 'cheb=', (cheb(nc,k,i),nc=2,ncoef)
                   do l = 1, nspec
                      call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', specmmr)
                      call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                  refindex_aer_sw=specrefindex)
                   !   volf = specmmr(k,i)/specdens
                   !   print*, 'l=', l, 'vol(l)=', volf
                      print*, 'l=', l, 'specmmr=',specmmr(k,i),'specdens=',specdens
                      print*, 'isw=', isw, 'specrefindex(isw)=', specrefindex(isw)
                   !   print*, 'specdens=', specdens
                   end do
 
                   nerr_dopaer = nerr_dopaer + 1
                   if (nerr_dopaer >= nerrmax_dopaer) then
                      ! print*, '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                      ! call endrun('exit from '//subname)
                   end if
                end if

             end do
 
             do i=1,ncol
                tauxar(isw,k,i) = tauxar(isw,k,i) + dopaer(i)
                wa(isw,k,i)     = wa(isw,k,i)     + dopaer(i)*palb(i)
                ga(isw,k,i)     = ga(isw,k,i)     + dopaer(i)*palb(i)*pasm(i)
                fa(isw,k,i)     = fa(isw,k,i)     + dopaer(i)*palb(i)*pasm(i)*pasm(i)
             end do
 
          end do ! nlev 
       end do ! sw bands
 
       ! mode diagnostics
       ! The diagnostics are currently only output for the climate list.  Code mods will
       ! be necessary to provide output for the rad_diag lists.
       if (list_idx == 0) then
          do i = 1, nnite
             burden(idxnite(i))  = fillvalue
             aodmode(idxnite(i)) = fillvalue
             dustaodmode(idxnite(i)) = fillvalue
          end do
       end if
 
    end do ! nmodes

    ! Output visible band diagnostics for quantities summed over the modes
    ! These fields are put out for diagnostic lists as well as the climate list.
    ! LiXH did not added output code.
    do i = 1, nnite
       extinct(:,idxnite(i)) = fillvalue
       absorb(:,idxnite(i))  = fillvalue
       aodvis(idxnite(i))    = fillvalue
       aodabs(idxnite(i))    = fillvalue
    end do

    ! These diagnostics are output only for climate list
    if (list_idx == 0) then
       do i = 1, ncol
          if (aodvis(i) > 1.e-10_r8) then
             ssavis(i) = ssavis(i)/aodvis(i)
          else
             ssavis(i) = 0.925_r8
          endif
       end do

       do i = 1, nnite
          ssavis(idxnite(i))     = fillvalue

          aoduv(idxnite(i))      = fillvalue
          aodnir(idxnite(i))     = fillvalue

          burdendust(idxnite(i)) = fillvalue
          burdenso4(idxnite(i))  = fillvalue
          burdenpom(idxnite(i))  = fillvalue
          burdensoa(idxnite(i))  = fillvalue
          burdenbc(idxnite(i))   = fillvalue
          burdenseasalt(idxnite(i)) = fillvalue

          aodabsbc(idxnite(i))   = fillvalue

          dustaod(idxnite(i))    = fillvalue
          so4aod(idxnite(i))     = fillvalue
          pomaod(idxnite(i))     = fillvalue
          soaaod(idxnite(i))     = fillvalue
          bcaod(idxnite(i))      = fillvalue
          seasaltaod(idxnite(i)) = fillvalue
        end do

    end if

    deallocate(dgnumwet_m)
    deallocate(dgnumdry_m)
    deallocate(qaerwat_m)
    deallocate(wetdens_m)

    end subroutine modal_aero_sw



! Purpose : calculates aerosol lw radiative properties
    subroutine modal_aero_lw(ncol, list_idx, tauxar)
    ! io
    integer,             intent(in) :: ncol
    integer,             intent(in)  :: list_idx ! index of the climate or a diagnostic list
    real(r8),            intent(out) :: tauxar(nlwbands,nlev,ncol) ! layer absorption optical depth

    ! Local variables
    integer :: i, ifld, ilw, k, l, m, nc, ns
    integer :: nmodes
    integer :: nspec

    real(r8) :: dgnumwet(nlev,ncol)    ! number mode wet diameter
    real(r8) :: qaerwat(nlev,ncol)     ! aerosol water (g/g)
 
    real(r8), allocatable :: dgnumdry_m(:,:,:) ! number mode dry diameter for all modes
    real(r8), allocatable :: dgnumwet_m(:,:,:) ! number mode wet diameter for all modes
    real(r8), allocatable :: qaerwat_m(:,:,:)  ! aerosol water (g/g) for all modes
    real(r8), allocatable :: wetdens_m(:,:,:)  ! 


    real(r8) :: sigma_logr_aer          ! geometric standard deviation of number distribution
    real(r8) :: alnsg_amode
    real(r8) :: xrad(ncol)
    real(r8) :: cheby(ncoef,nlev,ncol)  ! chebychef polynomials

    real(r8) :: mass(nlev,ncol) ! layer mass

    real(r8)             :: specmmr(nlev,ncol)            ! species mass mixing ratio
    real(r8)             :: specdens                      ! species density (kg/m3)
    complex(r8)          :: specrefindex(nlwbands)        ! species refractive index

    real(r8) :: vol(ncol)       ! volume concentration of aerosol specie (m3/kg)
    real(r8) :: dryvol(ncol)    ! volume concentration of aerosol mode (m3/kg)
    real(r8) :: wetvol(ncol)    ! volume concentration of wet mode (m3/kg)
    real(r8) :: watervol(ncol)  ! volume concentration of water in each mode (m3/kg)
    real(r8) :: refr(ncol)      ! real part of refractive index
    real(r8) :: refi(ncol)      ! imaginary part of refractive index
    complex(r8) :: crefin(ncol) ! complex refractive index
    real(r8) :: refrtablw(prefr,nlwbands) ! table of real refractive indices for aerosols
    real(r8) :: refitablw(prefi,nlwbands) ! table of imag refractive indices for aerosols
    real(r8) :: absplw(ncoef,prefr,prefi,nlwbands) ! specific absorption

    integer  :: itab(ncol), jtab(ncol)
    real(r8) :: ttab(ncol), utab(ncol)
    real(r8) :: cabs(ncoef,ncol)
    real(r8) :: pabs(ncol)      ! parameterized specific absorption (m2/kg)
    real(r8) :: dopaer(ncol)    ! aerosol optical depth in layer

    integer, parameter :: nerrmax_dopaer=1000
    integer  :: nerr_dopaer = 0
    real(r8) :: volf             ! volume fraction of insoluble aerosol

    character(len=*), parameter :: subname = 'modal_aero_lw'
    !----------------------------------------------------------------------------

    ! initialize output variables
    tauxar(:,:,:ncol) = 0._r8

    !-----------LiXH Modified dry pressure to pressure---------->
    ! dry mass in each cell
    !mass(:ncol,:) = state%pdeldry(:ncol,:)*rga
    mass(:,:ncol)         = pstate%delp_at_pc_full_level%f(:,:ncol)*rga
    !<----------LiXH Modified dry pressure to pressure-----------

    ! loop over all aerosol modes
    call rad_cnst_get_info(list_idx, nmodes=nmodes)

    if (list_idx == 0) then
       allocate(dgnumwet_m(nmodes,nlev,ncol))
       allocate(dgnumdry_m(nmodes,nlev,ncol))
       allocate(qaerwat_m(nmodes,nlev,ncol))
       allocate(wetdens_m(nmodes,nlev,ncol))
       ! water uptake and wet radius for the climate list has already been calculated
       dgnumwet_m = pstate_cam%aerosol_dgnumwet%f
       qaerwat_m  = pstate_cam%aerosol_qaerwat%f
 
    else
       ! If doing a diagnostic calculation then need to calculate the wet radius
       ! and water uptake for the diagnostic modes
       call modal_aero_calcsize_diag(ncol,list_idx, dgnumdry_m)  
       call modal_aero_wateruptake_dr(ncol,list_idx, dgnumdry_m, dgnumwet_m, &
                                      qaerwat_m, wetdens_m)
    endif

    do m = 1, nmodes

       dgnumwet = dgnumwet_m(m,:,:)
       qaerwat  = qaerwat_m(m,:,:)

       ! get mode properties
       call rad_cnst_get_mode_props(list_idx, m, sigmag=sigma_logr_aer, refrtablw=refrtablw , &
          refitablw=refitablw, absplw=absplw)

       ! get mode info
       call rad_cnst_get_info(list_idx, m, nspec=nspec)

       ! calc size parameter for all columns
       ! this is the same calculation that's done in modal_size_parameters, but there
       ! some intermediate results are saved and the chebyshev polynomials are stored
       ! in a array with different index order.  Could be unified.
       do k = top_lev, nlev
          do i = 1, ncol
             alnsg_amode = log( sigma_logr_aer )
             ! convert from number diameter to surface area
             xrad(i) = log(0.5_r8*dgnumwet(k,i)) + 2.0_r8*alnsg_amode*alnsg_amode
             ! normalize size parameter
             xrad(i) = max(xrad(i), xrmin)
             xrad(i) = min(xrad(i), xrmax)
             xrad(i) = (2*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
             ! chebyshev polynomials
             cheby(1,k,i) = 1.0_r8
             cheby(2,k,i) = xrad(i)
             do nc = 3, ncoef
                cheby(nc,k,i) = 2.0_r8*xrad(i)*cheby(nc-1,k,i)-cheby(nc-2,k,i)
             end do
          end do
       end do
 
       do ilw = 1, nlwbands
 
          do k = top_lev, nlev
 
             ! form bulk refractive index. Use volume mixing for infrared
             crefin(:ncol) = (0._r8, 0._r8)
             dryvol(:ncol) = 0._r8
 
             ! aerosol species loop
             do l = 1, nspec
                call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', specmmr)
                call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                            refindex_aer_lw=specrefindex)
 
                do i = 1, ncol
                   vol(i)    = specmmr(k,i)/specdens
                   dryvol(i) = dryvol(i) + vol(i)
                   crefin(i) = crefin(i) + vol(i)*specrefindex(ilw)
                end do
             end do
 
             do i = 1, ncol
                watervol(i) = qaerwat(k,i)/rhoh2o
                wetvol(i)   = watervol(i) + dryvol(i)
                if (watervol(i) < 0.0_r8) then
                   if (abs(watervol(i)) .gt. 1.e-1_r8*wetvol(i)) then
                      print*, 'watervol,wetvol,dryvol=',watervol(i),wetvol(i),dryvol(i),' in '//subname
                   end if
                   watervol(i) = 0._r8
                   wetvol(i)   = dryvol(i)
                end if

                crefin(i) = crefin(i) + watervol(i)*crefwlw(ilw)
                if (wetvol(i) > 1.e-40_r8) crefin(i) = crefin(i)/wetvol(i)
                refr(i) = real(crefin(i))
                refi(i) = aimag(crefin(i))
             end do
 
             ! interpolate coefficients linear in refractive index
             ! first call calcs itab,jtab,ttab,utab
             itab(:ncol) = 0
             call binterp(absplw(:,:,:,ilw), ncol, ncoef, prefr, prefi, &
                          refr, refi, refrtablw(:,ilw), refitablw(:,ilw), &
                          itab, jtab, ttab, utab, cabs)
 
             ! parameterized optical properties
             do i = 1, ncol
                pabs(i) = 0.5_r8*cabs(1,i)
                do nc = 2, ncoef
                   pabs(i) = pabs(i) + cheby(nc,k,i)*cabs(nc,i)
                end do
                pabs(i)   = pabs(i)*wetvol(i)*rhoh2o
                pabs(i)   = max(0._r8,pabs(i))
                dopaer(i) = pabs(i)*mass(k,i)
             end do
 
             do i = 1, ncol
 
                if ((dopaer(i) <= -1.e-10_r8) .or. (dopaer(i) >= 20._r8)) then
 
                   print*, 'dopaer(',i,',',k,',',m,')=', dopaer(i)
                   print*, 'k=',k,' pabs=', pabs(i)
                   print*, 'wetvol=',wetvol(i),' dryvol=',dryvol(i),' watervol=',watervol(i)
                   print*, 'cabs=', (cabs(l,i),l=1,ncoef)
                   print*, 'crefin=', crefin(i)
                   print*, 'nspec=', nspec
                   do l = 1,nspec
                      call rad_cnst_get_aer_mmr(list_idx, m, l, 'a', specmmr)
                      call rad_cnst_get_aer_props(list_idx, m, l, density_aer=specdens, &
                                                  refindex_aer_lw=specrefindex)
                      volf = specmmr(k,i)/specdens
                      print*, 'l=',l,'vol(l)=',volf
                      print*, 'ilw=',ilw,' specrefindex(ilw)=',specrefindex(ilw)
                      print*, 'specdens=',specdens
                   end do
 
                   nerr_dopaer = nerr_dopaer + 1
                   if (nerr_dopaer >= nerrmax_dopaer) then
                      print*, '*** halting in '//subname//' after nerr_dopaer =', nerr_dopaer
                      call endrun()
                   end if

                end if
             end do
 
             do i = 1, ncol
                tauxar(ilw,k,i) = tauxar(ilw,k,i) + dopaer(i)
             end do
 
          end do ! k = top_lev, nlev
 
       end do  ! nlwbands
 
    end do ! m = 1, nmodes
 
    deallocate(dgnumdry_m)
    deallocate(dgnumwet_m)
    deallocate(qaerwat_m)
    deallocate(wetdens_m)

    end subroutine modal_aero_lw


!===============================================================================
! Private routines
!===============================================================================


! Purpose : read water refractive index file and set module data
    subroutine read_water_refindex(infilename)

    use grist_wrap_nf

    character(*), intent(in) :: infilename   ! modal optics filename

    ! local
    integer, parameter :: omode = 0         ! nf_nowrite
    integer            :: i, ierr
    integer            :: ncid              ! pio file handle
    integer            :: did               ! dimension ids
    integer            :: dimlen            ! dimension lengths
    integer            :: vid               ! variable ids
    real(r8) :: refrwsw(nswbands), refiwsw(nswbands) ! real, imaginary ref index for water visible
    real(r8) :: refrwlw(nlwbands), refiwlw(nlwbands) ! real, imaginary ref index for water infrared

    ! open file
    call wrap_open(trim(infilename), omode, ncid)

    ! inquire dimensions.  Check that file values match parameter values.
    call wrap_inq_dimid(ncid, 'lw_band', did)
    call wrap_inq_dimlen(ncid, did, dimlen)

    if (dimlen .ne. nlwbands) then
       print*, 'lw_band len=', dimlen, ' from ', infilename, ' ne nlwbands=', nlwbands
       call endrun('read_modal_optics: bad lw_band value')
    endif

    call wrap_inq_dimid(ncid, 'sw_band', did)
    call wrap_inq_dimlen(ncid, did, dimlen)

    if (dimlen .ne. nswbands) then
       print*, 'sw_band len=', dimlen, ' from ', infilename, ' ne nswbands=', nswbands
       call endrun('read_modal_optics: bad sw_band value')
    endif

    ! read variables
    call wrap_inq_varid(ncid, 'refindex_real_water_sw', vid)
    call wrap_get_var_realx(ncid, vid, refrwsw)

    call wrap_inq_varid(ncid, 'refindex_im_water_sw', vid)
    call wrap_get_var_realx(ncid, vid, refiwsw)

    call wrap_inq_varid(ncid, 'refindex_real_water_lw', vid)
    call wrap_get_var_realx(ncid, vid, refrwlw)

    call wrap_inq_varid(ncid, 'refindex_im_water_lw', vid)
    call wrap_get_var_realx(ncid, vid, refiwlw)

    ! set complex representation of refractive indices as module data
    do i = 1, nswbands
       crefwsw(i)  = cmplx(refrwsw(i), abs(refiwsw(i)))
    end do
    do i = 1, nlwbands
       crefwlw(i)  = cmplx(refrwlw(i), abs(refiwlw(i)))
    end do

    call wrap_close(ncid)

    end subroutine read_water_refindex


    subroutine modal_size_parameters(ncol, sigma_logr_aer, dgnumwet, radsurf, logradsurf, cheb)

    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: sigma_logr_aer  ! geometric standard deviation of number distribution
    real(r8), intent(in)  :: dgnumwet(:,:)   ! aerosol wet number mode diameter (m)
    real(r8), intent(out) :: radsurf(:,:)    ! aerosol surface mode radius
    real(r8), intent(out) :: logradsurf(:,:) ! log(aerosol surface mode radius)
    real(r8), intent(out) :: cheb(:,:,:)

    integer  :: i, k, nc
    real(r8) :: alnsg_amode
    real(r8) :: explnsigma
    real(r8) :: xrad(ncol) ! normalized aerosol radius
    !-------------------------------------------------------------------------------

    alnsg_amode = log(sigma_logr_aer)
    explnsigma = exp(2.0_r8*alnsg_amode*alnsg_amode)

    do k = top_lev, nlev
       do i = 1, ncol
          ! convert from number mode diameter to surface area
          radsurf(k,i) = 0.5_r8*dgnumwet(k,i)*explnsigma
          logradsurf(k,i) = log(radsurf(k,i))
          ! normalize size parameter
          xrad(i) = max(logradsurf(k,i),xrmin)
          xrad(i) = min(xrad(i),xrmax)
          xrad(i) = (2._r8*xrad(i)-xrmax-xrmin)/(xrmax-xrmin)
          ! chebyshev polynomials
          cheb(1,k,i) = 1._r8
          cheb(2,k,i) = xrad(i)
          do nc = 3, ncoef
             cheb(nc,k,i) = 2._r8*xrad(i)*cheb(nc-1,k,i)-cheb(nc-2,k,i)
          end do
       end do
    end do

    end subroutine modal_size_parameters


! Purpose :    bilinear interpolation of table
    subroutine binterp(table,ncol,km,im,jm,x,y,xtab,ytab,ix,jy,t,u,out)

    implicit none
    integer im,jm,km,ncol
    real(r8) table(km,im,jm),xtab(im),ytab(jm),out(km,ncol)
    integer i,ix(ncol),ip1,j,jy(ncol),jp1,k,ic
    real(r8) x(ncol),dx,t(ncol),y(ncol),dy,u(ncol), &
           tu(ncol),tuc(ncol),tcu(ncol),tcuc(ncol)

    if(ix(1).gt.0)go to 30
    if(im.gt.1)then
      do ic=1,ncol
        do i=1,im
          if(x(ic).lt.xtab(i))go to 10
        enddo
 10     ix(ic)=max0(i-1,1)
        ip1=min(ix(ic)+1,im)
        dx=(xtab(ip1)-xtab(ix(ic)))
        if(abs(dx).gt.1.e-20_r8)then
           t(ic)=(x(ic)-xtab(ix(ic)))/dx
        else
           t(ic)=0._r8
        endif
      end do
    else
      ix(:ncol)=1
      t(:ncol)=0._r8
    endif
    if(jm.gt.1)then
      do ic=1,ncol
        do j=1,jm
          if(y(ic).lt.ytab(j))go to 20
        enddo
 20     jy(ic)=max0(j-1,1)
        jp1=min(jy(ic)+1,jm)
        dy=(ytab(jp1)-ytab(jy(ic)))
        if(abs(dy).gt.1.e-20_r8)then
           u(ic)=(y(ic)-ytab(jy(ic)))/dy
           if(u(ic).lt.0._r8.or.u(ic).gt.1._r8)then
              print*, 'u,y,jy,ytab,dy=',u(ic),y(ic),jy(ic),ytab(jy(ic)),dy
           endif
        else
          u(ic)=0._r8
        endif
      end do
    else
      jy(:ncol)=1
      u(:ncol)=0._r8
    endif
 30 continue
    do ic=1,ncol
       tu(ic)=t(ic)*u(ic)
       tuc(ic)=t(ic)-tu(ic)
       tcuc(ic)=1._r8-tuc(ic)-u(ic)
       tcu(ic)=u(ic)-tu(ic)
       jp1=min(jy(ic)+1,jm)
       ip1=min(ix(ic)+1,im)
       do k=1,km
          out(k,ic)=tcuc(ic)*table(k,ix(ic),jy(ic))+tuc(ic)*table(k,ip1,jy(ic))   &
             +tu(ic)*table(k,ip1,jp1)+tcu(ic)*table(k,ix(ic),jp1)
       end do
    enddo
    return
    end subroutine binterp


 end module modal_aer_opt
