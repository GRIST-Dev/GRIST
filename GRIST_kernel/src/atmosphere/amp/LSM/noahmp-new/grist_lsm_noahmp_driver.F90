
!-----------------------------------------------------------------------
! Purpose: Main driver of Noahmp LSM
! Intially ported and tested to GRIST by cheyz, lixh at 2020/04 (from 
! WRF-V3.8.1?)
! Rev: 1. Seperate for CAM and WRF physics. Layout format changed.
!      2. Add some options by AMIPW_PHYSICS, AMIPW_CLIMATE. The latter is 
!         for climate running with WRFphys, to increase surface albedo 
!         over land ice/snow.
!      3. Some bug fixed and minor changes for consistency. Regress old 
!         codes using REGLSM; WRF_PHYSICS has some inline changes also
!-----------------------------------------------------------------------

 module grist_lsm_noahmp_driver

    use noahmp_tables,                    only: read_landuse_tbl, albd, slmo, sfem, sfz0, therin, sfhc, scfx, li
    use module_sf_noahmpdrv
#ifdef AMIPC_PHYSICS
    use phys_control,                     only: lsm_scheme
#endif
#ifdef AMIPW_PHYSICS
    use grist_PhysW_nml_module,           only: lsm_scheme=>PhysW_lm_scheme
#endif
    use grist_lsm_noahmp_vars
    use noahmp_nml_control   
    use module_sf_noahmp_init,            only:noahmp_init       
!--cheyz
    use grist_mpi
    use grist_data_types
    use grist_lnd_static_vars_module
    use grist_zenith,                      only: zenith, orb_decl
    use grist_physics_data_structure,      only: pstate !%z_at_pc_face_level  ! change to input var later
    use grist_time_manager,                only: get_curr_date  
    use grist_nml_module,                  only: start_ymd, start_tod, model_timestep

    implicit none
       
    public:: grist_noahmp_run 
    
    ! logical,parameter  :: rdmaxalb = .false. !use snow albedo from geogrid;false use table values
    ! logical,parameter  :: myj      = .false. !true if using Mellor-Yamada PBL scheme.
    ! logical,parameter  :: frpcpn   = .false.
    ! logical,parameter  :: rdlai2d  = .false.
       
    LOGICAL              :: FNDSOILW , FNDSNOWH  !  
    integer              :: EM_CORE=2
    integer,parameter    :: sf_urban_physics = 0 !activate urban canopy model (=0: no urban canopy) 

    CONTAINS

    subroutine grist_noahmp_run(ncol, itimestep,  nlev, nsoil,   &
                                uu3d, vv3d, tt3d, qq3d, pp3d,    &
#ifdef AMIPW_PHYSICS
                                dz3d, qsfc_lsm, snowh_lsm, xice_lsm, emiss_lsm, coszrs, julian_in, &
 !xinyf test                               dz3d,  snowh_lsm, xice_lsm, emiss_lsm, coszrs, julian_in, &
#endif
                                swd2d, glw2d, preci2d,           &   !in
#ifdef AMIPC_PHYSICS
                                snow_lsm, xice_lsm,              &
#endif
                                lwup, shflx, qflx, taux, tauy,   &   !out
                                albdvis,albivis, albdnir, albinir, ts )  !out

    integer,intent(in)                            :: ncol, itimestep
    integer,intent(in)                            :: nlev, nsoil  
#ifdef AMIPC_PHYSICS
    real(r8),dimension(nlev,ncol),intent(in)      :: uu3d, vv3d, tt3d, qq3d
    real(r8),dimension(nlev+1,ncol),intent(in)    :: pp3d
#endif
#ifdef AMIPW_PHYSICS
    real(r8),dimension(ncol,nlev,1),intent(in)    :: uu3d, vv3d, tt3d, qq3d, pp3d, dz3d
    real(r8),dimension(ncol)       ,intent(inout) :: qsfc_lsm, snowh_lsm, xice_lsm, emiss_lsm ! snowh: mm->m, xice: fractional, emiss: >0.9
    real(r8),dimension(ncol)       ,intent(in)    :: coszrs
    real(r8)                       ,intent(in)    :: julian_in
#endif
    real(r8),dimension(ncol),intent(in)           :: swd2d, glw2d, preci2d
#ifdef AMIPC_PHYSICS
    real(r8),dimension(ncol),intent(inout)        :: snow_lsm, xice_lsm
#endif
    real(r8),dimension(ncol),intent(inout)        :: lwup,shflx,qflx,taux,tauy
    real(r8),dimension(ncol),intent(out)          :: albdvis,albivis, albdnir, albinir
    real(r8),dimension(ncol),intent(inout)        :: ts
! local
    integer                                       :: yr, mn, dy, sc, hour, minute, second
    real(r8)                                      :: gmt
    character(len=300)                            :: filename
    character(len=24)                             :: currentdate 
    integer                                       :: julyr, i, j, k,ii
    !----------------LiXH Test----------------
    real(r8)                                    :: julian
    !----------------LiXH Test----------------
    !--temporary vars
    character(len=300)                          :: filename_out, sss
    !--local vars
   ! integer, parameter                      :: max_cats = 100 , max_seas = 12, luseas=2, lucats=33
   ! real(r8), dimension(max_cats, max_seas) :: albd, slmo, sfem, sfz0, therin, sfhc
   ! real(r8), dimension(max_cats)           :: scfx, li
    real(r8)                                :: diff_Height 
    integer , parameter                     :: iswater=17, isice=15
    integer                                 :: ls, lc,  isn,isn_mid, iskip, iluvv, ierr, ilu

    ! calcultion of the current_date  with specific form--
    call get_curr_date(start_ymd, start_tod, itimestep, model_timestep, yr, mn, dy, sc) ! yizhang: no minus-1 ensure same julian (get_julgmt should -1)

    hour   = int(sc/3600.)
    minute = int((sc-hour*3600)/60.)
    second = sc-hour*3600-minute*60
    currentdate =''
    write(currentdate(1:19), fmt='(i4.4, 5(a1, i2.2))') yr,'-', mn,'-', dy,'_', hour, ':',minute,':', second

    !********************************************************************************** 
    !  2020/5/27 add
    !----------------------------------------------------------------------------------
    call monthly_min_max ( greenfrac_2d , vegmax_2d , &
            ids , ide , jds , jde , kds , kde , &
            ims , ime , jms , jme , kms , kme , &
            its , ite , jts , jte , kts , kte )
    call monthly_interp_to_date ( greenfrac_2d , currentdate, vegfra_2d , &
            ids , ide , jds , jde , kds , kde , &
            ims , ime , jms , jme , kms , kme , &
            its , ite , jts , jte , kts , kte )

            ! Determine season 1: summer  2: winter 
            !call get_julgmt(currentdate,julyr,julian,gmt)
    call get_julgmt_tod(currentdate,julyr,julian,gmt)

#ifdef AMIPW_PHYSICS
! do a replacement for bit consistency, 
! these two julians are same (using itimestep),
! but have machine-level difference due to different calculation
! lead to zenith-generated coszin difference, lat&lon are exactly the same
! Prefer using this for WRF physics (can be used for CAMphys if needed)
        julian = julian_in
#endif

!#ifdef AMIPC_PHYSICS
    isn_mid=1
    if(julian.lt.105. .or. julian .gt. 288. )isn_mid=2
!#else
!    isn=1
    !    if(julian.lt.105.or.julian.gt.288)isn=2
        ! if(cen_lat.lt.0.0)isn=3-isns
        
!#endif

! this below are uncessary lines
!    open(299, file='LANDUSE.TBL',form='formatted',status='old',iostat=ierr)
    
!    do iskip=1, 103 ! 215: modis ! 103: mminlu='modified_igbp_modis_noah'
!        read (299,*)
!    end do
!    ls=1    !--summer
!    do lc=1,lucats
!        read (299,*)li(lc),albd(lc,ls),slmo(lc,ls),sfem(lc,ls),        &
!                    sfz0(lc,ls),therin(lc,ls),scfx(lc),sfhc(lc,ls)
!    end do
!    read (299,*)
 
!    ls=2   !--winter
!    do lc=1,lucats
!        read (299,*)li(lc),albd(lc,ls),slmo(lc,ls),sfem(lc,ls),        &
!                    sfz0(lc,ls),therin(lc,ls),scfx(lc),sfhc(lc,ls)
!    end do
!    close(299)
        

    do j = jts, jte
        do i = its, ite
            ! if ( snow_2d(i,j) .ge. 10. ) then
            !     snowc_2d(i,j) = 1.
            !     else
            !     snowc_2d(i,j) = 0.0
            ! end if
                ilu= ivgtyp_2d(i,j)
                !  set no-data points to water
                if(ilu==0)then
                ilu=iswater
                end if
                
#ifdef AMIPC_PHYSICS
                if (xlatin_2d(i,j) .lt. 0.0) then
                    isn=3-isn_mid  !southern hemisphere
                else
                    isn=isn_mid !northern hemisphere 
                end if

                albedo_2d(i,j)= albd(ilu,isn)/100._r8  
                if(snowc_2d(i,j) .gt. 0.5) then             
                albedo_2d(i,j)= albedo_2d(i,j)*(1.+scfx(ilu))
                end if
    
                if(xice_2d(i,j).ge.xice_thres_1d)then
                xland_2d(i,j)=1.0  !change to land?
                albedo_2d(i,j)=albd(isice,isn)/100._r8  
                endif
#endif
        enddo
    enddo

    ! call land-surface scheme: might add 1D-version as another case later
    lsm_select: select case (trim(lsm_scheme))
    case("noahmp")
!        filename='grist_lsm_noahmp.nml'
!        call read_nml_noahmp(filename)
#ifdef AMIPC_PHYSICS 
        do ii=1, nlev
            u_phy(ii,:)    =  uu3d(nlev+1-ii,:)
            v_phy(ii,:)    =  vv3d(nlev+1-ii,:)
            t3d(ii,:)      =  tt3d(nlev+1-ii,:)
            qv3d(ii,:)     =  qq3d(nlev+1-ii,:)    ! full level
            p8w3d(ii,:)    =  pp3d(nlev+1+1-ii,:)  ! face level
        end do
        do ii=2, nlev
            dz8w(ii,:)= -pstate%z_at_pc_full_level%f(nlev+2-ii,1:ncol) &      ! not right?
                        +pstate%z_at_pc_full_level%f(nlev+1-ii,1:ncol)
        end do
        dz8w(1,:)=2.0_r8*(pstate%z_at_pc_full_level%f(nlev, 1:ncol)-pstate%z_at_pc_face_level%f(nlev+1,1:ncol))
#endif

#ifdef AMIPW_PHYSICS
! Added by YZ for WRF-physics
! note that, the orignal cheyz implementation of noahmp assumes the CAM vertical structure, nlev=GRIST's full level
! so all these six 3D vars are GRIST-nlev dimension; 
! for WRF, we initialize noahmp using nlevel (grist-nlevp)
! Also note that for these 3D variables, noahmp only uses the lowest level; 3D var was used for clean iostream
        u_phy(1:nlev,1:ncol)   =  transpose(uu3d(1:ncol,1:nlev,1))
        v_phy(1:nlev,1:ncol)   =  transpose(vv3d(1:ncol,1:nlev,1))
        t3d  (1:nlev,1:ncol)   =  transpose(tt3d(1:ncol,1:nlev,1))
        qv3d (1:nlev,1:ncol)   =  transpose(qq3d(1:ncol,1:nlev,1))
        p8w3d(1:nlev,1:ncol)   =  transpose(pp3d(1:ncol,1:nlev,1))
        dz8w (1:nlev,1:ncol)   =  transpose(dz3d(1:ncol,1:nlev,1))
#endif
        swdown   =  swd2d    !--all zero ! corrected: 04/27
        glw      =  glw2d
#ifdef AMIPC_PHYSICS
        precip_in=  preci2d*model_timestep*1000.0_r8  ! total precipitaion [mm]
        snow     =  snow_lsm                          ! ![mm] = [kg/m2], CheYZ, LiXH
#endif
#ifdef AMIPW_PHYSICS
        precip_in=  preci2d                           ! total precipitaion [mm]
        snowh    =  snowh_lsm*1e-3_r8                 ! mm->m
        emiss    =  emiss_lsm
        qsfc     =  qsfc_lsm ! ocean point is from sfclay, land points overwriteen here
#endif
        !call get_julgmt(currentdate,julyr,julian,gmt)
        julian_1d = julian
        call zenith(julian_1d, coszin, ncol, xlongin, xlatin)

        call lsm_noahmp_from_grist()   ! convert grsit vars to noahmp vars

        !sfcrunoff, udrunoff  !   do not consider now??

        call noahmplsm (itimestep = itimestep,  yr=yr,    dt=model_timestep,  julian  = julian_1d,   coszin=coszin_2d,  xlatin  = xlatin_2d,   & 
                        dz8w      = dz8w_2d,       dzs      = dzs_2d,                dx=dx_1d,                                      & ! IN : Model configuration 
                        ivgtyp    = ivgtyp_2d,     isltyp   = isltyp_2d,      vegfra   = vegfra_2d,                                 & ! IN : Vegetation/Soil characteristics
                        vegmax    = vegmax_2d,     tmn      = tmn_2d,         xland    = xland_2d,                                  & ! IN : Vegetation/Soil characteristics
                        xice      = xice_2d,       xice_thres= xice_thres_1d,                                                       & ! IN : Vegetation/Soil characteristics
                        !!!!!! 
                        idveg    =	    idveg,     iopt_crs =	iopt_crs, iopt_btr  =	iopt_btr,  iopt_run =	iopt_run,           & ! IN : User options
                        iopt_sfc =  	iopt_sfc,  iopt_frz = 	iopt_frz, iopt_inf  = 	iopt_inf,  iopt_rad =  	iopt_rad,           & ! IN : User options
                        iopt_alb =   	iopt_alb,  iopt_snf =  	iopt_snf, iopt_tbot = 	iopt_tbot, iopt_stc = 	iopt_stc,           & ! IN : User options
                        iopt_gla = 	    iopt_gla,  iopt_rsf =  	iopt_rsf, iz0tlnd   =	iz0tlnd,   nsoil    =   nsoil,              & ! IN : User options
                        !!!!!! 
                        t3d       = t3d_2d,        qv3d     = qv3d_2d,        u_phy=u_phy_2d,  v_phy=v_phy_2d,  swdown=swdown_2d,   & ! IN : Forcing
                        glw       = glw_2d,        p8w3d    = p8w3d_2d,       precip_in=precip_in_2d,    sr       = sr_2d,          & ! IN : Forcing
                        !!!!!!     
                        tsk       = tsk_2d,        hfx      = hfx_2d,         qfx = qfx_2d, lh  = lh_2d, grdflx = grdflx_2d,        & !IN/OUT : flux related to surface_layer?  
                        smstav   = smstav_2d,      smstot   = smstot_2d,                                                            & !IN/OUT
                        sfcrunoff = sfcrunoff_2d,  udrunoff = udrunoff_2d,    albedo   = albedo_2d,      snowc    = snowc_2d,       & !IN/OUT
                        smois     = smois_2d,      sh2o     = sh2o_2d,        tslb     = tslb_2d,        snow     = snow_2d,        & !IN/OUT
                        snowh     = snowh_2d,      canwat   = canwat_2d,      acsnom   = acsnom_2d,      acsnow   = acsnow_2d,      & !IN/OUT
                        emiss     = emiss_2d,      qsfc     = qsfc_2d,        z0       = z0_2d,          znt      = znt_2d,         & !IN/OUT  !
                        !!!!!!!!!!!! 
                        isnowxy =  	isnowxy_2d,     tvxy =      tvxy_2d,        tgxy  =     tgxy_2d,        canicexy =  canicexy_2d,        & 
                        canliqxy =  canliqxy_2d,    eahxy = 	eahxy_2d,       tahxy=      tahxy_2d,       cmxy     =  cmxy_2d,            &
                        chxy    =   chxy_2d,        fwetxy =    fwetxy_2d,      sneqvoxy =  sneqvoxy_2d,    alboldxy =  alboldxy_2d,        & 
                        qsnowxy =  	qsnowxy_2d,     wslakexy =  wslakexy_2d,                                                                & 
                        zwtxy =     zwtxy_2d,       waxy =      waxy_2d,        wtxy =      wtxy_2d,        tsnoxy = 	tsnoxy_2d,          &    
                        zsnsoxy =  	zsnsoxy_2d,     snicexy =   snicexy_2d,     snliqxy =   snliqxy_2d,     lfmassxy =  lfmassxy_2d,        &   
                        rtmassxy =  rtmassxy_2d,    stmassxy = 	stmassxy_2d,    woodxy =   	woodxy_2d,      stblcpxy =  stblcpxy_2d,        & 
                        fastcpxy =  fastcpxy_2d,    xlaixy =    xlaixy_2d,      xsaixy =    xsaixy_2d,      taussxy = 	taussxy_2d,         &   
                        smoiseq =  	smoiseq_2d,     smcwtdxy =  smcwtdxy_2d,    deeprechxy =deeprechxy_2d,  rechxy =    rechxy_2d,          &     
                        grainxy =   grainxy_2d,     gddxy =  	gddxy_2d,                                                                   &  
                        !!!!!!!!!!! 
                        t2mvxy =   	t2mvxy_2d,      t2mbxy =    t2mbxy_2d,                                                                  & !! OUT Noah MP only 
                        q2mvxy =    q2mvxy_2d,      q2mbxy =    q2mbxy_2d,      tradxy =   	tradxy_2d,                                      &     
                        neexy =     neexy_2d,       gppxy =     gppxy_2d,       nppxy =     nppxy_2d,       fvegxy =    fvegxy_2d,          &       
                        runsfxy = 	runsfxy_2d,     runsbxy =  	runsbxy_2d,     ecanxy =    ecanxy_2d,      edirxy =    edirxy_2d,          &     
                        etranxy =   etranxy_2d,     fsaxy =     fsaxy_2d,       firaxy = 	firaxy_2d,      aparxy =   	aparxy_2d,          &      
                        psnxy =     psnxy_2d,       savxy =     savxy_2d,       sagxy =     sagxy_2d,       rssunxy =   rssunxy_2d,         &     
                        rsshaxy = 	rsshaxy_2d,     bgapxy =   	bgapxy_2d,      wgapxy =    wgapxy_2d,      tgvxy =     tgvxy_2d,           & 
                        tgbxy =     tgbxy_2d,       chvxy =     chvxy_2d,       chbxy = 	chbxy_2d,       shgxy =    	shgxy_2d,           & 
                        shcxy =     shcxy_2d,       shbxy =     shbxy_2d,       evgxy =     evgxy_2d,       evbxy =     evbxy_2d,           &       
                        ghvxy = 	ghvxy_2d,       ghbxy =    	ghbxy_2d,       irgxy =     irgxy_2d,       ircxy =     ircxy_2d,           & 
                        irbxy =     irbxy_2d,       trxy =      trxy_2d,        evcxy = 	evcxy_2d,       chleafxy = 	chleafxy_2d,        & 
                        chucxy =   	chucxy_2d,      chv2xy =    chv2xy_2d,      chb2xy =   	chb2xy_2d,                                      &
                        taux2d =    taux2d_2d,      tauy2d =    tauy2d_2d,      fire2d =    fire2d_2d,&    !cheyz 20200427
                        albd2d =    albd2d_2d,      albi2d =    albi2d_2d, &  !--04/27
                        !!!!!!!!!!!                                   
                        ids = ids, ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,                                              &
                        ims = ims, ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,                                              &
                        its = its, ite = ite , jts = jts , jte = jte , kts = kts , kte = kte )
                            !--optional--
                            !mp_rainc=mp_rainc_2d, mp_rainnc=mp_rainnc_2d, mp_shcv=mp_shcv_2d, mp_snow= mp_snow_2d,                             &
                            !mp_graup=mp_graup_2d,        mp_hail=mp_hail_2d   )
        case default   
    end select lsm_select
         
! copy local arrays to grist grid:
    call lsm_noahmp_to_grist()  
         
! here back to pstate: out
! noahmp only update land points, LiXH
    lwup      = fire2d
    shflx     = hfx
    qflx      = qfx
    taux      = taux2d
    tauy      = tauy2d
!--------------LiXH Test------------------
!    albdvis   = albd2d(1,1:ncol)
!    albdnir   = albd2d(2,1:ncol)
!    albivis   = albi2d(1,1:ncol)
!    albinir   = albi2d(2,1:ncol)

    albdvis   = albedo
    albdnir   = albedo
    albivis   = albedo
    albinir   = albedo
!--------------LiXH Test------------------
 
    ts        = tsk
#ifdef AMIPC_PHYSICS
    snow_lsm(1:ncol)  = snow(1:ncol) ![mm] = [kg/m2], CheYZ, LiXH
    xice_lsm(1:ncol)  = xice(1:ncol)
#endif
#ifdef AMIPW_PHYSICS
! use albedo as in WRF physics
    albdvis   = albedo
    albdnir   = albedo
    albivis   = albedo
    albinir   = albedo
    qsfc_lsm  = qsfc
    emiss_lsm = emiss
    snowh_lsm = snowh*1e3_r8
    xice_lsm  = xice
#endif
    return
 end subroutine grist_noahmp_run
  
 subroutine calc_coszen(ims,ime,jms,jme,its,ite,jts,jte,  &
                        julian,xtime,gmt, &
                        declin,degrad,xlon,xlat,coszen,hrang)
    implicit none
    integer, intent(in) :: ims,ime,jms,jme,its,ite,jts,jte
    real(r8)  , intent(in)    :: julian,declin,xtime,gmt,degrad
    real(r8)  , dimension(ims:ime,jms:jme), intent(in)    :: xlat,xlon
    real(r8)  , dimension(ims:ime,jms:jme), intent(inout) :: coszen,hrang

    integer :: i,j
    real(r8)  :: da,eot,xt24,tloctm,xxlat

    da=6.2831853071795862*(julian-1)/365.
    eot=(0.000075+0.001868*cos(da)-0.032077*sin(da) &
            -0.014615*cos(2*da)-0.04089*sin(2*da))*(229.18)
    xt24=mod(xtime,1440.)+eot
    do j=jts,jte
        do i=its,ite
            tloctm=gmt+xt24/60.+xlon(i,j)/15.
            hrang(i,j)=15.*(tloctm-12.)*degrad
            xxlat=xlat(i,j)*degrad
            coszen(i,j)=sin(xxlat)*sin(declin) &
                        +cos(xxlat)*cos(declin) *cos(hrang(i,j))
        enddo
    enddo
    return
 end subroutine calc_coszen

 subroutine monthly_interp_to_date ( field_in , date_str , field_out , &
               ids , ide , jds , jde , kds , kde , &
               ims , ime , jms , jme , kms , kme , &
               its , ite , jts , jte , kts , kte )
     implicit none

     integer , intent(in)        :: ids , ide , jds , jde , kds , kde , &
                                     ims , ime , jms , jme , kms , kme , &
                                     its , ite , jts , jte , kts , kte

     character (len=24) , intent(in) :: date_str  !2019-07-10_00:00:00
     real(r8), dimension(ims:ime,12,jms:jme) , intent(in)  :: field_in
     real(r8), dimension(ims:ime,   jms:jme) , intent(out) :: field_out
     integer :: i , j , l
     integer , dimension(0:13) :: middle
     integer :: target_julyr , target_julday , target_date
     integer :: julyr , julday , int_month , month1 , month2
     real(r8)  :: gmt
     character (len=4) :: yr
     character (len=2) :: mon , day15

     write(day15,fmt='(i2.2)') 15
     do l = 1 , 12
         write(mon,fmt='(i2.2)') l
         call get_julgmt ( date_str(1:4)//'-'//mon//'-'//day15//'_'//'00:00:00.0000' , julyr , julday , gmt )
         middle(l) = julyr*1000 + julday
     end do

     l = 0
     middle(l) = middle( 1) - 31

     l = 13
     middle(l) = middle(12) + 31

     call get_julgmt ( date_str , target_julyr , target_julday , gmt )
     target_date = target_julyr * 1000 + target_julday
     find_month : do l = 0 , 12
         if ( ( middle(l) .lt. target_date ) .and. ( middle(l+1) .ge. target_date ) ) then
             do j = jts , min ( jde-1 , jte )
             do i = its , min (ide-1 , ite )
                 !if ( skip_middle_points_t ( ids , ide , jds , jde , i , j , em_width , hold_ups ) ) cycle   ?? for wrf boundary?
                 int_month = l
                 if ( ( int_month .eq. 0 ) .or. ( int_month .eq. 12 ) ) then
                     month1 = 12
                     month2 =  1
                 else
                     month1 = int_month
                     month2 = month1 + 1
                 end if
                 field_out(i,j) =  ( field_in(i,month2,j) * ( target_date - middle(l)   ) + &
                                     field_in(i,month1,j) * ( middle(l+1) - target_date ) ) / &
                                     ( middle(l+1) - middle(l) )
             end do
             end do
             exit find_month
         end if
     end do find_month

 end subroutine monthly_interp_to_date
  
! to calculate vegmax 
  subroutine monthly_min_max ( field_in , field_max , &
              ids , ide , jds , jde , kds , kde , &
              ims , ime , jms , jme , kms , kme , &
              its , ite , jts , jte , kts , kte )

  implicit none

  integer , intent(in) :: ids , ide , jds , jde , kds , kde , &
                          ims , ime , jms , jme , kms , kme , &
                          its , ite , jts , jte , kts , kte

  real(r8), dimension(ims:ime,12,jms:jme) , intent(in)  :: field_in
  real(r8), dimension(ims:ime,   jms:jme) , intent(out) :: field_max
  real(r8), dimension(ims:ime,   jms:jme)               :: field_min

  integer :: i , j , l
  real(r8) :: minner , maxxer

  do j = jts , min(jde-1,jte)
  do i = its , min(ide-1,ite)
  minner = field_in(i,1,j)
  maxxer = field_in(i,1,j)
  do l = 2 , 12
      if ( field_in(i,l,j) .lt. minner ) then
          minner = field_in(i,l,j)
      end if
      if ( field_in(i,l,j) .gt. maxxer ) then
          maxxer = field_in(i,l,j)
      end if
  end do
  field_min(i,j) = minner
  field_max(i,j) = maxxer
  end do
  end do

  end subroutine monthly_min_max

  subroutine get_julgmt(date_str,julyr,julday,gmt)
    implicit none
  
    character (len=24) , intent(in) :: date_str
    integer, intent(out  ) :: julyr
    integer, intent(out  ) :: julday
    real(r8), intent(out  ) :: gmt
  
    integer :: ny , nm , nd , nh , ni , ns , nt
    integer :: my1, my2, my3, monss
    integer, dimension(12) :: mmd
    data mmd/31,28,31,30,31,30,31,31,30,31,30,31/
    call split_date_char ( date_str , ny , nm , nd , nh , ni , ns , nt )
    gmt=nh+float(ni)/60.+float(ns)/3600.
    my1=mod(ny,4)
    my2=mod(ny,100)
    my3=mod(ny,400)
    if(my1.eq.0.and.my2.ne.0.or.my3.eq.0)mmd(2)=29
    julday=nd
    julyr=ny
    do monss=1,nm-1
      julday=julday+mmd(monss)
    enddo
  end subroutine get_julgmt

!----------------------LiXH Test-------------------------
  subroutine get_julgmt_tod(date_str,julyr,julday,gmt)
      implicit none
  
      character (len=24) , intent(in) :: date_str
      integer, intent(out  ) :: julyr
      !-----------LiXH Test--------------
      !integer, intent(out  ) :: julday
      real(r8), intent(out  ) :: julday
      !-----------LiXH Test--------------
      real(r8), intent(out  ) :: gmt
  
      integer :: ny , nm , nd , nh , ni , ns , nt
      integer :: my1, my2, my3, monss
      integer, dimension(12) :: mmd
      data mmd/31,28,31,30,31,30,31,31,30,31,30,31/
      call split_date_char ( date_str , ny , nm , nd , nh , ni , ns , nt )
      gmt=nh+float(ni)/60.+float(ns)/3600.
      my1=mod(ny,4)
      my2=mod(ny,100)
      my3=mod(ny,400)
      if(my1.eq.0.and.my2.ne.0.or.my3.eq.0)mmd(2)=29
      julday=nd
      julyr=ny
      do monss=1,nm-1
        julday=julday+mmd(monss)
      enddo
      !-----------LiXH Test--------------
      julday=julday+(nh*3600.+ni*60.+ns)/86400.
      !-----------LiXH Test--------------
#ifndef REGLSM
            julday=julday-1 !yizhang modified for consisten julian
#endif
   end subroutine get_julgmt_tod
!----------------------LiXH Test-------------------------


   subroutine split_date_char ( date , century_year , month , day , hour , minute , second , ten_thousandth)
       implicit none
    
       character(len=24) , intent(in) :: date 
    
       integer , intent(out) :: century_year , month , day , hour , minute , second , ten_thousandth
       
       read(date,fmt='(    i4.4)') century_year
       read(date,fmt='( 5x,i2.2)') month
       read(date,fmt='( 8x,i2.2)') day
       read(date,fmt='(11x,i2.2)') hour
       read(date,fmt='(14x,i2.2)') minute
       read(date,fmt='(17x,i2.2)') second
       read(date,fmt='(20x,i4.4)') ten_thousandth
    
   end subroutine split_date_char

   subroutine current_time(start_date, start_sec, current_date, current_sec, inc_s,iday)
       implicit none
 
       ! io
           integer(i4), intent(in)           :: start_date, start_sec
           integer(i4), intent(in)           :: inc_s
           integer(i4), intent(out)          :: current_date, current_sec, iday
       ! local
           integer(i4)                       :: i
           integer(i4)                       :: year, mon, day, julday
           integer(i4)                       :: month(12),             &
                                                leap_year_month(12)
           logical                           :: leap_year
       
           month = (/31,28,31,30,31,30,31,31,30,31,30,31/)
           leap_year_month = (/31,29,31,30,31,30,31,31,30,31,30,31/)
       
           year  = start_date/10000
           mon   = (start_date-year*10000)/100
           day   = start_date-year*10000-mon*100
       
           current_sec = start_sec+inc_s
           do while(current_sec .ge. 86400)
               day = day+1
               current_sec = current_sec-86400
           end do
       
           julday = day
           if((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or.     &
              (mod(year,400) .eq. 0))then
               do i = 1, mon-1
                   julday = julday+leap_year_month(i)
               end do
           else
               do i = 1, mon-1
                   julday = julday+month(i)
               end do
           end if
       
           do while(julday .gt. 365)
               if((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or.     &
                  (mod(year,400) .eq. 0))then
                   leap_year = .true.
               else
                   leap_year = .false.
               end if
       
               if(leap_year .and. julday .gt. 366)then
                   year = year+1
                   julday = julday-366
               else if(leap_year .and. julday .eq. 366)then
                   exit
               else
                   year = year+1
                   julday = julday-365
               end if
           end do
          iday =julday
          
           mon = 1
           if((mod(year,4) .eq. 0 .and. mod(year,100) .ne. 0) .or.     &
              (mod(year,400) .eq. 0))then
               do i = 1, 12
                   if(julday .gt. leap_year_month(i))then
                       julday = julday-leap_year_month(i)
                       mon    = mon+1
                   else
                       day    = julday
                       exit
                   end if
               end do
           else
               do i = 1, 12
                   if(julday .gt. month(i))then
                       julday = julday-month(i)
                       mon    = mon+1
                   else
                       day    = julday
                       exit
                   end if
               end do
           end if
       
           current_date = year*10000+mon*100+day
       
   end subroutine current_time 
       
   subroutine lsm_noahmp_from_grist() !( its,ite,jts,jte,kts,kte )
   implicit none
   
   integer:: i,j,k
       
   do j = jts,jte
       do  k = kts, kme !kte  
           do i = its,ite
               !arrays for 3D:
               dz8w_2d     (i, k, j)   = dz8w    (k,i)
               t3d_2d  	(i, k, j)	= t3d     (k,i)
               qv3d_2d  	(i, k, j)	= qv3d    (k,i)
               u_phy_2d   	(i, k, j)	= u_phy   (k,i)
               v_phy_2d   	(i, k, j)	= v_phy   (k,i)
               p8w3d_2d	(i, k, j)	= p8w3d   (k,i)
                
           end do    
       end do
   end do
 
   !!!!!!!!!!!!!!!!!!!!!!!!!
   do k = 1,nsoil
      dzs_2d(k) = dzs1(k)
   end do

   do j = jts,jte
       do k = 1, nsoil
           do i = its,ite
   
               sh2o_2d (i,k,j)    = sh2o (k,i)   
               smoiseq_2d(i,k,j) = smoiseq(k,i)    
               smois_2d (i,k,j)   = smois (k,i)
               tslb_2d (i,k,j)    = tslb (k,i)
           
           end do        
       end do
   end do
   
   !!newly added

   do j = jts,jte
       do k = 1, 12 ! no. of month
           do i = its,ite
   
               greenfrac_2d (i,k,j)    = greenfrac (k,i)   
               
           end do        
       end do
   end do

   !!!!!!!!!!!!!!!!!!!!!!!!!
   do j = jts,jte
   do i = its,ite
       ivgtyp_2d   (i, j)	=   ivgtyp     	(i)    
       isltyp_2d   (i, j)	=   isltyp     	(i)
       xice_2d    	(i, j)	=   xice      	(i)
       tsk_2d    	(i, j)	=   tsk     	(i)
       tmn_2d 	    (i, j)	=   tmn  	    (i)
       
       vegfra_2d 	(i, j)	=	vegfra  	(i)
       coszin_2d 	(i, j)	=	coszin  	(i)
       xlatin_2d   (i, j)	=	xlatin    	(i)
       vegmax_2d  	(i, j)	=	vegmax   	(i)
   
       xland_2d  	(i, j)	=	xland   	(i)
       swdown_2d   (i, j)	=	swdown      (i)
       glw_2d 	    (i, j)	=	glw  	    (i)
       precip_in_2d(i, j)	=	precip_in   (i)
       sr_2d  	    (i, j)	=	sr   	    (i)
       hfx_2d      (i, j)	=	hfx       	(i)
       qfx_2d     	(i, j)	=	qfx      	(i)
       lh_2d   	(i, j)	=	lh      	(i)
       grdflx_2d   (i, j)	=	grdflx     	(i)
       smstav_2d  	(i, j)	=	smstav   	(i)
       smstot_2d	(i, j)	=	smstot 	    (i)
       sfcrunoff_2d(i, j)	=	sfcrunoff   (i)
       udrunoff_2d (i, j)	=	udrunoff    (i)
       albedo_2d   (i, j)	=	albedo     	(i)
       snowc_2d   	(i, j)	=	snowc    	(i)
       snow_2d   	(i, j)	=	snow    	(i)
       snowh_2d   	(i, j)	=	snowh    	(i)
       canwat_2d   (i, j)	=	canwat    	(i)
       acsnom_2d 	(i, j)	=	acsnom  	(i)
       acsnow_2d  	(i, j)	=	acsnow   	(i)
       emiss_2d  	(i, j)	=	emiss   	(i)
       qsfc_2d 	(i, j)	=	qsfc  	    (i)
       z0_2d   	(i, j)	=	z0    	    (i)
       znt_2d	    (i, j)	=	znt 	    (i)
       
   end do
   end do 

   do j = jts,jte
   do i = its,ite

       tsnoxy_2d 	    (i, :, j)	=	tsnoxy  	(:,i)
       zsnsoxy_2d  	(i, :, j)	=	zsnsoxy   	(:,i)
       snicexy_2d   	(i, :, j)	=	snicexy    	(:,i)
       snliqxy_2d    	(i, :, j)	=	snliqxy     (:,i)

       !!!!!!!!!!!!!!!!
       isnowxy_2d  	(i, j)	=	isnowxy   	(i)
       tvxy_2d      	(i, j)	=	tvxy       	(i)
       tgxy_2d       	(i, j)	=	tgxy        (i)
       canicexy_2d  	(i, j)	=	canicexy   	(i)
       canliqxy_2d  	(i, j)	=	canliqxy   	(i)
       eahxy_2d 	    (i, j)	=	eahxy   	(i)
       tahxy_2d	    (i, j)	=	tahxy 	    (i)
       cmxy_2d	        (i, j)	=	cmxy 	    (i)
       chxy_2d	        (i, j)	=	chxy 	    (i)
       fwetxy_2d	    (i, j)	=	fwetxy  	(i)
       sneqvoxy_2d	    (i, j)	=	sneqvoxy 	(i)
       alboldxy_2d  	(i, j)	=	alboldxy   	(i)
       qsnowxy_2d  	(i, j)	=	qsnowxy   	(i)
       wslakexy_2d  	(i, j)	=	wslakexy   	(i)
       zwtxy_2d      	(i, j)	=	zwtxy       (i)
       waxy_2d      	(i, j)	=	waxy       	(i)
       wtxy_2d      	(i, j)	=	wtxy       	(i)
       lfmassxy_2d  	(i, j)	=	lfmassxy   	(i)
       rtmassxy_2d  	(i, j)	=	rtmassxy   	(i)
       stmassxy_2d 	(i, j)	=	stmassxy  	(i)
       woodxy_2d   	(i, j)	=	woodxy    	(i)
       grainxy_2d   	(i, j)	=	grainxy    	(i)
       gddxy_2d  	    (i, j)	=	gddxy   	(i)
       stblcpxy_2d  	(i, j)	=	stblcpxy   	(i)
       fastcpxy_2d   	(i, j)	=	fastcpxy    (i)
       xsaixy_2d    	(i, j)	=	xsaixy     	(i)
       xlaixy_2d    	(i, j)	=	lai         (i)  !???  xlaixy     	(i)
   
   
       !======inout vars=====  
       taussxy_2d   	(i, j)	=	taussxy  	(i)
       smcwtdxy_2d  	(i, j)	=	smcwtdxy   	(i)
       deeprechxy_2d   (i, j)	=	deeprechxy  (i)
       rechxy_2d    	(i, j)	=	rechxy     	(i)

       !======only out vars======just inital 0 is enough==========
       t2mvxy_2d   	(i, j)	=	t2mvxy    	(i)
       t2mbxy_2d    	(i, j)	=	t2mbxy     	(i)
       q2mvxy_2d       (i, j)	=	q2mvxy      (i)
       q2mbxy_2d       (i, j)	=	q2mbxy      (i)
       tradxy_2d   	(i, j)	=	tradxy    	(i)
       neexy_2d     	(i, j)	=	neexy      	(i)
       gppxy_2d        (i, j)	=	gppxy       (i)
       nppxy_2d     	(i, j)	=	nppxy      	(i)
       fvegxy_2d    	(i, j)	=	fvegxy     	(i)
       runsfxy_2d 	    (i, j)	=	runsfxy  	(i)
       runsbxy_2d  	(i, j)	=	runsbxy   	(i)
       ecanxy_2d    	(i, j)	=	ecanxy     	(i)
       edirxy_2d       (i, j)	=	edirxy      (i)
       etranxy_2d   	(i, j)	=	etranxy    	(i)
       fsaxy_2d     	(i, j)	=	fsaxy      	(i)
       firaxy_2d 	    (i, j)	=	firaxy  	(i)
       aparxy_2d   	(i, j)	=	aparxy    	(i)
       psnxy_2d     	(i, j)	=	psnxy      	(i)
       savxy_2d        (i, j)	=	savxy       (i)
       sagxy_2d     	(i, j)	=	sagxy      	(i)
       rssunxy_2d   	(i, j)	=	rssunxy    	(i)
       rsshaxy_2d 	    (i, j)	=	rsshaxy  	(i)
       bgapxy_2d   	(i, j)	=	bgapxy    	(i)
       wgapxy_2d    	(i, j)	=	wgapxy     	(i)
       tgvxy_2d        (i, j)	=	tgvxy       (i)
       tgbxy_2d     	(i, j)	=	tgbxy      	(i)
       chvxy_2d     	(i, j)	=	chvxy      	(i)
       chbxy_2d 	    (i, j)	=	chbxy  	    (i)
       shgxy_2d    	(i, j)	=	shgxy     	(i)
       shcxy_2d     	(i, j)	=	shcxy      	(i)
       shbxy_2d        (i, j)	=	shbxy       (i)
       evgxy_2d     	(i, j)	=	evgxy      	(i)
       evbxy_2d     	(i, j)	=	evbxy      	(i)
       ghvxy_2d 	    (i, j)	=	ghvxy    	(i)
       ghbxy_2d    	(i, j)	=	ghbxy     	(i)
       irgxy_2d     	(i, j)	=	irgxy      	(i)
       ircxy_2d        (i, j)	=	ircxy       (i)
       irbxy_2d     	(i, j)	=	irbxy      	(i)
       trxy_2d      	(i, j)	=	trxy       	(i)
       evcxy_2d 	    (i, j)	=	evcxy  	    (i)
       chleafxy_2d 	(i, j)	=	chleafxy  	(i)
       chucxy_2d   	(i, j)	=	chucxy    	(i)
       chv2xy_2d       (i, j)	=	chv2xy      (i)
       chb2xy_2d   	(i, j)	=	chb2xy    	(i)
       
   enddo
   enddo

   end subroutine lsm_noahmp_from_grist

   subroutine lsm_noahmp_to_grist()  
       implicit none
    
       integer:: i,j,k
   
       !-----------     
       
       !do not needed anymore
       ! do j = jts,jte
       !     do  k = kts, kte  
       !         do i = its,ite
       !             !arrays for 3D:
       !             dz8w    (k,i)          =       dz8w_2d     (i, k, j)   
       !             t3d     (k,i)          =       t3d_2d      (i, k, j)
       !             qv3d    (k,i)          =       qv3d_2d     (i, k, j)
       !             u_phy   (k,i)          =       u_phy_2d    (i, k, j)
       !             v_phy   (k,i)          =       v_phy_2d    (i, k, j)
       !             p8w3d   (k,i)          =       p8w3d_2d    (i, k, j)
       !         end do    
       !     end do
       ! end do
   
       !!!!!!!!!!!!!!!!!!!!!!!!!
       do k = 1,nsoil
          dzs_2d(k) = dzs1(k)
       end do
   
       do j = jts,jte
           do k = 1, nsoil
               do i = its,ite
   
                    sh2o (k,i)    =       sh2o_2d (i,k,j)    
                    smoiseq(k,i)  =       smoiseq_2d(i,k,j)  
                    smois (k,i)   =       smois_2d (i,k,j)   
                    tslb (k,i)    =       tslb_2d (i,k,j)    
   
               end do        
           end do
       end do
       
      !!!newly added
   
       do j = jts,jte
           do k = 1, 12 ! no. of month
               do i = its,ite
   
                   greenfrac (k,i)  =  greenfrac_2d (i,k,j)    
                   
               end do        
           end do
       end do
       !---04/27
       do j = jts,jte
           do k = 1, 2 ! no. of band
               do i = its,ite
                  albd2d(k,i)=albd2d_2d(i,k,j)
                  albi2d(k,i)=albi2d_2d(i,k,j)
               end do        
           end do
       end do
       !!!!!!!!!!!!!!!!!!!!!!!!!
       do j = jts,jte
       do i = its,ite

           taux2d      (i) =               taux2d_2d   (i,j)
           tauy2d      (i) =               tauy2d_2d   (i,j)
           fire2d      (i) =               fire2d_2d   (i,j)  !--04/27

           ivgtyp      (i) =               ivgtyp_2d   (i,j)
           isltyp      (i) =               isltyp_2d   (i,j)
           xice        (i) =               xice_2d     (i,j)
           tsk         (i) =               tsk_2d      (i,j)
           tmn         (i) =               tmn_2d      (i,j)
                       
           vegfra      (i) =               vegfra_2d   (i,j)
           coszin      (i) =               coszin_2d   (i,j)
           xlatin      (i) =               xlatin_2d   (i,j)
           vegmax      (i) =               vegmax_2d   (i,j)
                       
           xland       (i) =               xland_2d    (i,j)
           swdown      (i) =               swdown_2d   (i,j)
           glw         (i) =               glw_2d      (i,j)
           precip_in   (i) =               precip_in_2d(i,j)
           sr          (i) =               sr_2d       (i,j)
           hfx         (i) =               hfx_2d      (i,j)
           qfx         (i) =               qfx_2d      (i,j)
           lh          (i) =               lh_2d       (i,j)
           grdflx      (i) =               grdflx_2d   (i,j)
           smstav      (i) =               smstav_2d   (i,j)
           smstot      (i) =               smstot_2d   (i,j)
           sfcrunoff   (i) =               sfcrunoff_2d(i,j)
           udrunoff    (i) =               udrunoff_2d (i,j)
           albedo      (i) =               albedo_2d   (i,j)
           snowc       (i) =               snowc_2d    (i,j)
           snow        (i) =               snow_2d     (i,j)
           snowh       (i) =               snowh_2d    (i,j)
           canwat      (i) =               canwat_2d   (i,j)
           acsnom      (i) =               acsnom_2d   (i,j)
           acsnow      (i) =               acsnow_2d   (i,j)
           emiss       (i) =               emiss_2d    (i,j)
           qsfc        (i) =               qsfc_2d     (i,j)
           z0          (i) =               z0_2d       (i,j)
           znt         (i) =               znt_2d      (i,j)
   
   
       end do
       end do
   
       do j = jts,jte
       do i = its,ite
   
           tsnoxy      (:,i)   =   tsnoxy_2d       (i, :, j)
           zsnsoxy     (:,i)   =   zsnsoxy_2d      (i, :, j)
           snicexy     (:,i)   =   snicexy_2d      (i, :, j)
           snliqxy     (:,i)   =   snliqxy_2d      (i, :, j)
   
           !!!!!!!!!!!!!!!!
           isnowxy   (i)   =   isnowxy_2d      (i, j)
           tvxy      (i)   =   tvxy_2d         (i, j)
           tgxy      (i)   =   tgxy_2d         (i, j)
           canicexy  (i)   =   canicexy_2d     (i, j)
           canliqxy  (i)   =   canliqxy_2d     (i, j)
           eahxy     (i)   =   eahxy_2d        (i, j)
           tahxy     (i)   =   tahxy_2d        (i, j)
           cmxy      (i)   =   cmxy_2d         (i, j)
           chxy      (i)   =   chxy_2d         (i, j)
           fwetxy    (i)   =   fwetxy_2d       (i, j)
           sneqvoxy  (i)   =   sneqvoxy_2d     (i, j)
           alboldxy  (i)   =   alboldxy_2d     (i, j)
           qsnowxy   (i)   =   qsnowxy_2d      (i, j)
           wslakexy  (i)   =   wslakexy_2d     (i, j)
           zwtxy     (i)   =   zwtxy_2d        (i, j)
           waxy      (i)   =   waxy_2d         (i, j)
           wtxy      (i)   =   wtxy_2d         (i, j)
           lfmassxy  (i)   =   lfmassxy_2d     (i, j)
           rtmassxy  (i)   =   rtmassxy_2d     (i, j)
           stmassxy  (i)   =   stmassxy_2d     (i, j)
           woodxy    (i)   =   woodxy_2d       (i, j)
           grainxy   (i)   =   grainxy_2d      (i, j)
           gddxy     (i)   =   gddxy_2d        (i, j)
           stblcpxy  (i)   =   stblcpxy_2d     (i, j)
           fastcpxy  (i)   =   fastcpxy_2d     (i, j)
           xsaixy    (i)   =   xsaixy_2d       (i, j)
           lai       (i)   =   xlaixy_2d       (i, j)
           taussxy   (i)   =   taussxy_2d      (i, j)
           smcwtdxy  (i)   =   smcwtdxy_2d     (i, j)
           deeprechxy(i)   =   deeprechxy_2d   (i, j)
           rechxy    (i)   =   rechxy_2d       (i, j)
   
           !======only out vars======just inital 0 is enough??
           t2mvxy      (i) =   t2mvxy_2d       (i, j)
           t2mbxy      (i) =   t2mbxy_2d       (i, j)
           q2mvxy      (i) =   q2mvxy_2d       (i, j)
           q2mbxy      (i) =   q2mbxy_2d       (i, j)
           tradxy      (i) =   tradxy_2d       (i, j)
           neexy       (i) =   neexy_2d        (i, j)
           gppxy       (i) =   gppxy_2d        (i, j)
           nppxy       (i) =   nppxy_2d        (i, j)
           fvegxy      (i) =   fvegxy_2d       (i, j)
           runsfxy     (i) =   runsfxy_2d      (i, j)
           runsbxy     (i) =   runsbxy_2d      (i, j)
           ecanxy      (i) =   ecanxy_2d       (i, j)
           edirxy      (i) =   edirxy_2d       (i, j)
           etranxy     (i) =   etranxy_2d      (i, j)
           fsaxy       (i) =   fsaxy_2d        (i, j)
           firaxy      (i) =   firaxy_2d       (i, j)
           aparxy      (i) =   aparxy_2d       (i, j)
           psnxy       (i) =   psnxy_2d        (i, j)
           savxy       (i) =   savxy_2d        (i, j)
           sagxy       (i) =   sagxy_2d        (i, j)
           rssunxy     (i) =   rssunxy_2d      (i, j)
           rsshaxy     (i) =   rsshaxy_2d      (i, j)
           bgapxy      (i) =   bgapxy_2d       (i, j)
           wgapxy      (i) =   wgapxy_2d       (i, j)
           tgvxy       (i) =   tgvxy_2d        (i, j)
           tgbxy       (i) =   tgbxy_2d        (i, j)
           chvxy       (i) =   chvxy_2d        (i, j)
           chbxy       (i) =   chbxy_2d        (i, j)
           shgxy       (i) =   shgxy_2d        (i, j)
           shcxy       (i) =   shcxy_2d        (i, j)
           shbxy       (i) =   shbxy_2d        (i, j)
           evgxy       (i) =   evgxy_2d        (i, j)
           evbxy       (i) =   evbxy_2d        (i, j)
           ghvxy       (i) =   ghvxy_2d        (i, j)
           ghbxy       (i) =   ghbxy_2d        (i, j)
           irgxy       (i) =   irgxy_2d        (i, j)
           ircxy       (i) =   ircxy_2d        (i, j)
           irbxy       (i) =   irbxy_2d        (i, j)
           trxy        (i) =   trxy_2d         (i, j)
           evcxy       (i) =   evcxy_2d        (i, j)
           chleafxy    (i) =   chleafxy_2d     (i, j)
           chucxy      (i) =   chucxy_2d       (i, j)
           chv2xy      (i) =   chv2xy_2d       (i, j)
           chb2xy      (i) =   chb2xy_2d       (i, j)
       enddo
       enddo
    
    end subroutine lsm_noahmp_to_grist
  end module grist_lsm_noahmp_driver
