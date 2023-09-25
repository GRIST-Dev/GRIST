!--cheyz- 2019/11
!*********************************************

module grist_lsm_noahmp_vars
    use grist_constants,                    only: i4, r8, deg2rad, zero
    use noahmp_nml_control,                 only: nsoil
    
    implicit none
    public 
    !save
     
    integer                 :: ids,ide,jds,jde,kds,kde
    integer                 :: ims,ime,jms,jme,kms,kme
    integer                 :: its,ite,jts,jte,kts,kte
   
    !=================================================================================================================
    !.. variables and arrays related to land-surface model: noahmp
    !=================================================================================================================

    real(r8),dimension(:,:),allocatable :: &
                    tvxy_2d,      tgxy_2d,       canicexy_2d,  canliqxy_2d, stmassxy_2d,    &
        stblcpxy_2d,  fastcpxy_2d,   xlaixy_2d,    xsaixy_2d,    taussxy_2d,  woodxy_2d,     &
        eahxy_2d, qsnowxy_2d,  wslakexy_2d,  zwtxy_2d,      waxy_2d,      wtxy_2d,           & 
            lfmassxy_2d,  rtmassxy_2d,      & 
        smcwtdxy_2d,  deeprechxy_2d, rechxy_2d,    grainxy_2d,   gddxy_2d,      &
        t2mvxy_2d,   t2mbxy_2d,    q2mvxy_2d,     q2mbxy_2d,  tradxy_2d,   neexy_2d,         &
        gppxy_2d,      nppxy_2d,     fvegxy_2d,    runsfxy_2d, runsbxy_2d,  ecanxy_2d,       &  
        edirxy_2d,     etranxy_2d,   fsaxy_2d,     firaxy_2d, aparxy_2d,   psnxy_2d,         & 
        savxy_2d,      sagxy_2d,     rssunxy_2d,   rsshaxy_2d, bgapxy_2d,   wgapxy_2d,       &  
        tgvxy_2d,      tgbxy_2d,     chvxy_2d,     chbxy_2d, shgxy_2d,    shcxy_2d,          &
        shbxy_2d,      evgxy_2d,     evbxy_2d,     ghvxy_2d, ghbxy_2d,    irgxy_2d,          &
        ircxy_2d,      irbxy_2d,     trxy_2d,      evcxy_2d, chleafxy_2d, chucxy_2d,         &
        chv2xy_2d,     chb2xy_2d,    tahxy_2d,     cmxy_2d,  chxy_2d,     fwetxy_2d,         &
        sneqvoxy_2d, alboldxy_2d, &
        !--newly addedd
        greenfrac, taux2d_2d, tauy2d_2d, fire2d_2d, albd2d, albi2d 
   

    real(r8),dimension(:) ,    allocatable :: dzs_2d
    real(r8),dimension(:,:,:),allocatable :: dz8w_2d, t3d_2d,qv3d_2d, u_phy_2d,v_phy_2d,p8w3d_2d,  &
        smois_2d, sh2o_2d, tslb_2d, tsnoxy_2d, zsnsoxy_2d, snicexy_2d, snliqxy_2d, smoiseq_2d, &
        greenfrac_2d,st_albedo_2d, albd2d_2d, albi2d_2d


    real(r8),dimension(:,:),allocatable :: & 
        dx_2d,   vegfra_2d, coszin_2d, xlatin_2d, xlongin_2d,   &
        vegmax_2d,  tmn_2d, xland_2d,  xice_2d,   &   
        swdown_2d,      glw_2d, precip_in_2d,        sr_2d,  tsk_2d,    &
        hfx_2d,      qfx_2d,     lh_2d,   grdflx_2d,    smstav_2d,  smstot_2d,sfcrunoff_2d,  &
        udrunoff_2d,    albedo_2d,    snowc_2d,  snow_2d,   &
        snowh_2d,   canwat_2d,   acsnom_2d, acsnow_2d,  emiss_2d,  qsfc_2d, z0_2d,   znt_2d

    real(r8),   dimension(:,:),allocatable :: mp_rainc_2d, mp_rainnc_2d, mp_shcv_2d, mp_snow_2d,mp_graup_2d, mp_hail_2d   

    integer ,dimension(:,:),allocatable :: ivgtyp_2d, isltyp_2d,  &  ! vegetation type  ! soil type
                                           isnowxy_2d 

    real(r8)           :: dt_1d, dx_1d, xice_thres_1d
    real(r8)           :: julian_1d     

    !--------------1D vars---grist vars----------------
    real(r8),dimension(:),allocatable :: &
                    tvxy,      tgxy,       canicexy,  canliqxy, stmassxy,    &
        stblcpxy,  fastcpxy,   xlaixy,    xsaixy,    taussxy,  woodxy,     &
        eahxy, qsnowxy,  wslakexy,  zwtxy,      waxy,      wtxy,           & 
            lfmassxy,  rtmassxy,      & 
        smcwtdxy,  deeprechxy, rechxy,    grainxy,   gddxy,      &
        t2mvxy,   t2mbxy,    q2mvxy,     q2mbxy,  tradxy,   neexy,         &
        gppxy,      nppxy,     fvegxy,    runsfxy, runsbxy,  ecanxy,       &  
        edirxy,     etranxy,   fsaxy,     firaxy, aparxy,   psnxy,         & 
        savxy,      sagxy,     rssunxy,   rsshaxy, bgapxy,   wgapxy,       &  
        tgvxy,      tgbxy,     chvxy,     chbxy, shgxy,    shcxy,          &
        shbxy,      evgxy,     evbxy,     ghvxy, ghbxy,    irgxy,          &
        ircxy,      irbxy,     trxy,      evcxy, chleafxy, chucxy,         &
        chv2xy,     chb2xy,    tahxy,     cmxy,  chxy,     fwetxy,         &
        sneqvoxy, alboldxy, taux2d, tauy2d, fire2d  
    

    real(r8),dimension(:) ,    allocatable :: dzs
    real(r8),dimension(:,:),allocatable :: dz8w, t3d,qv3d, u_phy,v_phy,p8w3d,  &
        smois, sh2o, tslb, tsnoxy, zsnsoxy, snicexy, snliqxy, smoiseq, albbck  !newly add


    real(r8),dimension(:),allocatable ::  dx,   vegfra, coszin, xlatin,xlongin,              &
        vegmax,  tmn, xland,  xice,   swdown,   glw, precip_in,    sr,  tsk,                      &
        hfx,      qfx,     lh,   grdflx,    smstav,  smstot,sfcrunoff,  udrunoff,                 &
        albedo,    snowc,  snow, snowh,   canwat,   acsnom, acsnow,  emiss,  qsfc, z0,   znt

    real(r8),   dimension(:),allocatable :: mp_rainc, mp_rainnc, mp_shcv, mp_snow,mp_graup, mp_hail   

    integer ,dimension(:),allocatable :: ivgtyp, isltyp,  &  ! vegetation type  ! soil type
                                         isnowxy 
    real(r8),   dimension(:),allocatable :: lai, chstarxy  !for init                                         

    contains


        !=================================================================================================================
        subroutine lsm_noahmp_init_domain(ncol,nlev)
         implicit none
         integer,intent(in):: ncol, nlev
        !=================================================================================================================
            !initialization of physics dimensions to mimic a rectangular grid for noahmp:
            ims=1   ; ime = ncol
            jms=1   ; jme=1
            kms=1   ; kme = nlev   !+1
        
            ids=ims ; ide=ime + 1
            jds=jms ; jde=jme + 1
            kds=kms ; kde=kme 
        
            its=ims ; ite = ime 
            jts=jms ; jte = jme
            kts=kms ; kte = kme    !-1
        
        end subroutine lsm_noahmp_init_domain 
       !=================================================================================================================


       !=================================================================================================================
        subroutine allocate_lsm_noahmp_2d(ncol, nlev, nsoil)  
       !=================================================================================================================
        implicit none
        integer, intent(in)  ::ncol,nlev,nsoil
        !-----------------------------------------------------------------------------------------------------------------
      
        if(.not.allocated(	taux2d_2d 	    ))		allocate(	taux2d_2d 	    (ims:ime,jms:jme), source = zero  ) ! must initialized, or diagnoised taux has inconsistency at different core numbers
        if(.not.allocated(	tauy2d_2d 	    ))		allocate(	tauy2d_2d 	    (ims:ime,jms:jme), source = zero  )
        if(.not.allocated(	fire2d_2d 	    ))		allocate(	fire2d_2d 	    (ims:ime,jms:jme), source = zero     )
        if(.not.allocated(	 albbck 	    ))		allocate(	 albbck 	    (ims:ime,jms:jme), source = zero     )
     

        if(.not.allocated(	albd2d_2d 	    ))		allocate(	albd2d_2d 	    (ims:ime,1:2,jms:jme), source = zero )
        if(.not.allocated(	albi2d_2d 	    ))		allocate(	albi2d_2d 	    (ims:ime,1:2,jms:jme), source = zero )

        !arrays for soil layer properties:
        if(.not.allocated(	dzs_2d 	        ))		allocate(	dzs_2d 	        (1:nsoil)    )
        if(.not.allocated(	smois_2d 	    ))		allocate(	smois_2d 	    (ims:ime, 1:nsoil,jms:jme)    )
        if(.not.allocated(	sh2o_2d  	    ))		allocate(	sh2o_2d  	    (ims:ime, 1:nsoil,jms:jme)    )
        if(.not.allocated(	tslb_2d    	    ))		allocate(	tslb_2d    	    (ims:ime, 1:nsoil,jms:jme)    )
        if(.not.allocated(	tsnoxy_2d 	    ))		allocate(	tsnoxy_2d 	    (ims:ime,-2:0,jms:jme)    )
        if(.not.allocated(	zsnsoxy_2d  	))		allocate(	zsnsoxy_2d  	(ims:ime,-2:nsoil,jms:jme)    )
        if(.not.allocated(	snicexy_2d   	))		allocate(	snicexy_2d   	(ims:ime,-2:0,jms:jme)    )
        if(.not.allocated(	snliqxy_2d    	))		allocate(	snliqxy_2d    	(ims:ime,-2:0,jms:jme)    )
        if(.not.allocated(	smoiseq_2d  	))		allocate(	smoiseq_2d (ims:ime, 1:nsoil, jms:jme)    )
        !--newly added
        if(.not.allocated(	greenfrac_2d  	))		allocate(	greenfrac_2d  (ims:ime, 1:12, jms:jme)    ) ! 12 month
        if(.not.allocated(	st_albedo_2d  	))		allocate(	st_albedo_2d  (ims:ime, 1:12, jms:jme)    ) ! 12 month
        
        
        !arrays for 3D:
        if(.not.allocated(	dz8w_2d  	    ))		allocate(	dz8w_2d  (ims:ime, kms:kme,jms:jme)    )
        if(.not.allocated(	t3d_2d  	    ))		allocate(	t3d_2d   (ims:ime, kms:kme,jms:jme)    )
        if(.not.allocated(	qv3d_2d  	    ))		allocate(	qv3d_2d  (ims:ime, kms:kme,jms:jme)    )
        if(.not.allocated(	u_phy_2d   	    ))		allocate(	u_phy_2d (ims:ime, kms:kme,jms:jme)    )
        if(.not.allocated(	v_phy_2d   	    ))		allocate(	v_phy_2d (ims:ime, kms:kme,jms:jme)    )
        if(.not.allocated(	p8w3d_2d	    ))		allocate(	p8w3d_2d (ims:ime, kms:kme,jms:jme)    )



        !other arrays:
        if(.not.allocated(	isnowxy_2d  	))		allocate(	isnowxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	tvxy_2d      	))		allocate(	tvxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	tgxy_2d       	))		allocate(	tgxy_2d       	(ims:ime,jms:jme)    )
        if(.not.allocated(	canicexy_2d  	))		allocate(	canicexy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	canliqxy_2d  	))		allocate(	canliqxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	eahxy_2d 	    ))		allocate(	eahxy_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	qsnowxy_2d   	))		allocate(	qsnowxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	wslakexy_2d  	))		allocate(	wslakexy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	zwtxy_2d      	))		allocate(	zwtxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	waxy_2d      	))		allocate(	waxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	wtxy_2d      	))		allocate(	wtxy_2d      	(ims:ime,jms:jme)    )
        
        if(.not.allocated(	lfmassxy_2d  	))		allocate(	lfmassxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	rtmassxy_2d  	))		allocate(	rtmassxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	stmassxy_2d 	))		allocate(	stmassxy_2d 	(ims:ime,jms:jme)    )
        if(.not.allocated(	woodxy_2d   	))		allocate(	woodxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	stblcpxy_2d  	))		allocate(	stblcpxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	fastcpxy_2d   	))		allocate(	fastcpxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	xlaixy_2d    	))		allocate(	xlaixy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	xsaixy_2d    	))		allocate(	xsaixy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	taussxy_2d 	    ))		allocate(	taussxy_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	smcwtdxy_2d  	))		allocate(	smcwtdxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	deeprechxy_2d 	))		allocate(	deeprechxy_2d 	(ims:ime,jms:jme)    )
        if(.not.allocated(	rechxy_2d    	))		allocate(	rechxy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	grainxy_2d   	))		allocate(	grainxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	gddxy_2d  	    ))		allocate(	gddxy_2d  	    (ims:ime,jms:jme)    )
        if(.not.allocated(	t2mvxy_2d   	))		allocate(	t2mvxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	t2mbxy_2d    	))		allocate(	t2mbxy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	q2mvxy_2d     	))		allocate(	q2mvxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	q2mbxy_2d       ))		allocate(	q2mbxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	tradxy_2d   	))		allocate(	tradxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	neexy_2d     	))		allocate(	neexy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	gppxy_2d      	))		allocate(	gppxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	nppxy_2d     	))		allocate(	nppxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	fvegxy_2d    	))		allocate(	fvegxy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	runsfxy_2d   	))		allocate(	runsfxy_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	runsbxy_2d  	))		allocate(	runsbxy_2d  	(ims:ime,jms:jme)    )
        if(.not.allocated(	ecanxy_2d    	))		allocate(	ecanxy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	edirxy_2d     	))		allocate(	edirxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	etranxy_2d   	))		allocate(	etranxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	fsaxy_2d     	))		allocate(	fsaxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	firaxy_2d 	    ))		allocate(	firaxy_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	aparxy_2d   	))		allocate(	aparxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	psnxy_2d     	))		allocate(	psnxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	savxy_2d      	))		allocate(	savxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	sagxy_2d     	))		allocate(	sagxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	rssunxy_2d   	))		allocate(	rssunxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	rsshaxy_2d 	    ))		allocate(	rsshaxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	bgapxy_2d   	))		allocate(	bgapxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	wgapxy_2d    	))		allocate(	wgapxy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	tgvxy_2d      	))		allocate(	tgvxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	tgbxy_2d     	))		allocate(	tgbxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	chvxy_2d     	))		allocate(	chvxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	chbxy_2d 	    ))		allocate(	chbxy_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	shgxy_2d    	))		allocate(	shgxy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	shcxy_2d     	))		allocate(	shcxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	shbxy_2d      	))		allocate(	shbxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	evgxy_2d     	))		allocate(	evgxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	evbxy_2d     	))		allocate(	evbxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	ghvxy_2d 	    ))		allocate(	ghvxy_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	ghbxy_2d    	))		allocate(	ghbxy_2d    	(ims:ime,jms:jme)    )
        if(.not.allocated(	irgxy_2d     	))		allocate(	irgxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	ircxy_2d      	))		allocate(	ircxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	irbxy_2d     	))		allocate(	irbxy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	trxy_2d      	))		allocate(	trxy_2d      	(ims:ime,jms:jme)    )
        if(.not.allocated(	evcxy_2d 	    ))		allocate(	evcxy_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	chleafxy_2d 	))		allocate(	chleafxy_2d 	(ims:ime,jms:jme)    )
        if(.not.allocated(	chucxy_2d   	))		allocate(	chucxy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	chv2xy_2d     	))		allocate(	chv2xy_2d     	(ims:ime,jms:jme)    )
        if(.not.allocated(	chb2xy_2d   	))		allocate(	chb2xy_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	tahxy_2d	    ))		allocate(	tahxy_2d	    (ims:ime,jms:jme)    )	
        if(.not.allocated(	cmxy_2d	        ))		allocate(	cmxy_2d        	(ims:ime,jms:jme)    )	
        if(.not.allocated(	chxy_2d	        ))		allocate(	chxy_2d	        (ims:ime,jms:jme)    )	
        if(.not.allocated(	fwetxy_2d	    ))		allocate(	fwetxy_2d	    (ims:ime,jms:jme)    )	
        if(.not.allocated(	sneqvoxy_2d	    ))		allocate(	sneqvoxy_2d	    (ims:ime,jms:jme)    )	
        if(.not.allocated(	alboldxy_2d  	))		allocate(	alboldxy_2d  	(ims:ime,jms:jme)    )	
    
        !!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!                      
        
       
        if(.not.allocated(	dx_2d   	    ))		allocate(	dx_2d   	    (ims:ime,jms:jme)    )
        if(.not.allocated(	vegfra_2d 	    ))		allocate(	vegfra_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	coszin_2d 	    ))		allocate(	coszin_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	xlatin_2d       ))		allocate(	xlatin_2d       (ims:ime,jms:jme)    )
        if(.not.allocated(	xlongin_2d       ))		allocate(	xlongin_2d       (ims:ime,jms:jme)    )
        if(.not.allocated(	vegmax_2d  	    ))		allocate(	vegmax_2d  	    (ims:ime,jms:jme)    )
        if(.not.allocated(	tmn_2d 		    ))		allocate(	tmn_2d   	    (ims:ime,jms:jme)    )
        if(.not.allocated(	xland_2d  	    ))		allocate(	xland_2d  	    (ims:ime,jms:jme)    )
        if(.not.allocated(	xice_2d    	    ))		allocate(	xice_2d    	    (ims:ime,jms:jme)    )
        if(.not.allocated(	swdown_2d       ))		allocate(	swdown_2d       (ims:ime,jms:jme)    )
        if(.not.allocated(	glw_2d 	        ))		allocate(	glw_2d 	        (ims:ime,jms:jme)    )
        if(.not.allocated(	precip_in_2d    ))		allocate(	precip_in_2d    (ims:ime,jms:jme)    )
        if(.not.allocated(	sr_2d  	        ))		allocate(	sr_2d    	    (ims:ime,jms:jme)    )
        if(.not.allocated(	tsk_2d    	    ))		allocate(	tsk_2d    	    (ims:ime,jms:jme)    )
        if(.not.allocated(	hfx_2d          ))		allocate(	hfx_2d          (ims:ime,jms:jme)    )
        if(.not.allocated(	qfx_2d     	    ))		allocate(	qfx_2d     	    (ims:ime,jms:jme)    )
        if(.not.allocated(	lh_2d   	    ))		allocate(	lh_2d   	    (ims:ime,jms:jme)    )
        if(.not.allocated(	grdflx_2d       ))		allocate(	grdflx_2d       (ims:ime,jms:jme)    )
        if(.not.allocated(	smstav_2d  	    ))		allocate(	smstav_2d  	    (ims:ime,jms:jme)    )
        if(.not.allocated(	smstot_2d	    ))		allocate(	smstot_2d	    (ims:ime,jms:jme)    )
        if(.not.allocated(	sfcrunoff_2d    ))		allocate(	sfcrunoff_2d    (ims:ime,jms:jme)    )
        if(.not.allocated(	udrunoff_2d     ))		allocate(	udrunoff_2d     (ims:ime,jms:jme)    )
        if(.not.allocated(	albedo_2d       ))		allocate(	albedo_2d       (ims:ime,jms:jme)    )
        if(.not.allocated(	snowc_2d   	    ))		allocate(	snowc_2d   	    (ims:ime,jms:jme)    )
        if(.not.allocated(	snow_2d   	    ))		allocate(	snow_2d   	    (ims:ime,jms:jme)    )
        if(.not.allocated(	snowh_2d   	    ))		allocate(	snowh_2d   	    (ims:ime,jms:jme)    )
        if(.not.allocated(	canwat_2d       ))		allocate(	canwat_2d       (ims:ime,jms:jme)    )
        if(.not.allocated(	acsnom_2d 	    ))		allocate(	acsnom_2d   	(ims:ime,jms:jme)    )
        if(.not.allocated(	acsnow_2d  	    ))		allocate(	acsnow_2d  	    (ims:ime,jms:jme)    )
        if(.not.allocated(	emiss_2d  	    ))		allocate(	emiss_2d  	    (ims:ime,jms:jme)    )
        if(.not.allocated(	qsfc_2d 	    ))		allocate(	qsfc_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	z0_2d   	    ))		allocate(	z0_2d   	    (ims:ime,jms:jme)    )
        if(.not.allocated(	znt_2d	        ))		allocate(	znt_2d	        (ims:ime,jms:jme)    )
                
        !!!!!!!!!!!!!!!
        if(.not.allocated(	mp_rainc_2d 	))		allocate(	mp_rainc_2d 	(ims:ime,jms:jme)    )
        if(.not.allocated(	mp_rainnc_2d 	))		allocate(	mp_rainnc_2d 	(ims:ime,jms:jme)    )
        if(.not.allocated(	mp_shcv_2d 	    ))		allocate(	mp_shcv_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	mp_snow_2d	    ))		allocate(	mp_snow_2d	    (ims:ime,jms:jme)    )
        if(.not.allocated(	mp_graup_2d 	))		allocate(	mp_graup_2d 	(ims:ime,jms:jme)    )
        if(.not.allocated(	mp_hail_2d   	))		allocate(	mp_hail_2d   	(ims:ime,jms:jme)    )

        !!!!!!integer
        if(.not.allocated(	ivgtyp_2d 	    ))		allocate(	ivgtyp_2d 	    (ims:ime,jms:jme)    )
        if(.not.allocated(	isltyp_2d  	    ))		allocate(	isltyp_2d  	    (ims:ime,jms:jme)    )
        if(.not.allocated(	isnowxy_2d   	))		allocate(	isnowxy_2d 	    (ims:ime,jms:jme)    )

   
    end subroutine allocate_lsm_noahmp_2d
    !=================================================================================================================



    !=================================================================================================================
    subroutine allocate_lsm_noahmp_1d(ncol, nlev, nsoil) ! allocate vars used in Grist
    !=================================================================================================================
    implicit none
    integer, intent(in)         :: ncol, nlev, nsoil

      
        if(.not.allocated(	taux2d 	    ))		allocate(	taux2d 	    (1:ncol)           )  !--20200427
        if(.not.allocated(	tauy2d 	    ))		allocate(	tauy2d 	    (1:ncol)           )
        if(.not.allocated(	fire2d 	    ))		allocate(	fire2d 	    (1:ncol)           )


        if(.not.allocated(	albd2d 	    ))		allocate(	albd2d 	    (1:2,1:ncol)           )
        if(.not.allocated(	albi2d 	    ))		allocate(	albi2d 	    (1:2,1:ncol)           )



        if(.not.allocated(	dzs 	    ))		allocate(	dzs 	    (1:nsoil)           )
        if(.not.allocated(	smois 	    ))		allocate(	smois 	    (1:nsoil,ncol)      )
        if(.not.allocated(	sh2o  	    ))		allocate(	sh2o  	    (1:nsoil,ncol)      )
        if(.not.allocated(	tslb    	))		allocate(	tslb    	(1:nsoil,ncol)      )
        if(.not.allocated(	tsnoxy 	    ))		allocate(	tsnoxy 	    (-2:0,ncol)         )
        if(.not.allocated(	zsnsoxy  	))		allocate(	zsnsoxy  	(-2:nsoil,ncol)     )
        if(.not.allocated(	snicexy   	))		allocate(	snicexy   	(-2:0,ncol)         )
        if(.not.allocated(	snliqxy    	))		allocate(	snliqxy    	(-2:0,ncol)         )
        if(.not.allocated(	smoiseq  	))		allocate(	smoiseq     (1:nsoil, ncol)     )
        !--newly added
        if(.not.allocated(	greenfrac  	))		allocate(	greenfrac   (1:12, ncol)        )
        
        !arrays for 3D:
        if(.not.allocated(	dz8w  	    ))		allocate(	dz8w  (nlev,ncol)    )
        if(.not.allocated(	t3d  	    ))		allocate(	t3d   (nlev,ncol)    )
        if(.not.allocated(	qv3d  	    ))		allocate(	qv3d  (nlev,ncol)    )
        if(.not.allocated(	u_phy   	))		allocate(	u_phy (nlev,ncol)    )
        if(.not.allocated(	v_phy   	))		allocate(	v_phy (nlev,ncol)    )
        if(.not.allocated(	p8w3d	    ))		allocate(	p8w3d (nlev,ncol)    )



        !other arrays:
        if(.not.allocated(	isnowxy  	))		allocate(	isnowxy  	(ncol)    )
        if(.not.allocated(	tvxy      	))		allocate(	tvxy      	(ncol)    )
        if(.not.allocated(	tgxy       	))		allocate(	tgxy       	(ncol)    )
        if(.not.allocated(	canicexy  	))		allocate(	canicexy  	(ncol)    )
        if(.not.allocated(	canliqxy  	))		allocate(	canliqxy  	(ncol)    )
        if(.not.allocated(	eahxy 	    ))		allocate(	eahxy 	    (ncol)    )
        if(.not.allocated(	qsnowxy   	))		allocate(	qsnowxy  	(ncol)    )
        if(.not.allocated(	wslakexy  	))		allocate(	wslakexy  	(ncol)    )
        if(.not.allocated(	zwtxy      	))		allocate(	zwtxy      	(ncol)    )
        if(.not.allocated(	waxy      	))		allocate(	waxy      	(ncol)    )
        if(.not.allocated(	wtxy      	))		allocate(	wtxy      	(ncol)    )
        
        if(.not.allocated(	lfmassxy  	))		allocate(	lfmassxy  	(ncol)    )
        if(.not.allocated(	rtmassxy  	))		allocate(	rtmassxy  	(ncol)    )
        if(.not.allocated(	stmassxy 	))		allocate(	stmassxy 	(ncol)    )
        if(.not.allocated(	woodxy   	))		allocate(	woodxy   	(ncol)    )
        if(.not.allocated(	stblcpxy  	))		allocate(	stblcpxy  	(ncol)    )
        if(.not.allocated(	fastcpxy   	))		allocate(	fastcpxy   	(ncol)    )
        if(.not.allocated(	xlaixy    	))		allocate(	xlaixy    	(ncol)    )
        if(.not.allocated(	xsaixy    	))		allocate(	xsaixy    	(ncol)    )
        if(.not.allocated(	taussxy 	))		allocate(	taussxy 	(ncol)    )
        if(.not.allocated(	smcwtdxy  	))		allocate(	smcwtdxy  	(ncol)    )
        if(.not.allocated(	deeprechxy 	))		allocate(	deeprechxy 	(ncol)    )
        if(.not.allocated(	rechxy    	))		allocate(	rechxy    	(ncol)    )
        if(.not.allocated(	grainxy   	))		allocate(	grainxy   	(ncol)    )
        if(.not.allocated(	gddxy  	    ))		allocate(	gddxy  	    (ncol)    )
        if(.not.allocated(	t2mvxy   	))		allocate(	t2mvxy   	(ncol)    )
        if(.not.allocated(	t2mbxy    	))		allocate(	t2mbxy    	(ncol)    )
        if(.not.allocated(	q2mvxy     	))		allocate(	q2mvxy     	(ncol)    )
        if(.not.allocated(	q2mbxy      ))		allocate(	q2mbxy   	(ncol)    )
        if(.not.allocated(	tradxy   	))		allocate(	tradxy   	(ncol)    )
        if(.not.allocated(	neexy     	))		allocate(	neexy     	(ncol)    )
        if(.not.allocated(	gppxy      	))		allocate(	gppxy      	(ncol)    )
        if(.not.allocated(	nppxy     	))		allocate(	nppxy     	(ncol)    )
        if(.not.allocated(	fvegxy    	))		allocate(	fvegxy    	(ncol)    )
        if(.not.allocated(	runsfxy   	))		allocate(	runsfxy 	(ncol)    )
        if(.not.allocated(	runsbxy  	))		allocate(	runsbxy  	(ncol)    )
        if(.not.allocated(	ecanxy    	))		allocate(	ecanxy    	(ncol)    )
        if(.not.allocated(	edirxy     	))		allocate(	edirxy     	(ncol)    )
        if(.not.allocated(	etranxy   	))		allocate(	etranxy   	(ncol)    )
        if(.not.allocated(	fsaxy     	))		allocate(	fsaxy     	(ncol)    )
        if(.not.allocated(	firaxy 	    ))		allocate(	firaxy 	    (ncol)    )
        if(.not.allocated(	aparxy   	))		allocate(	aparxy   	(ncol)    )
        if(.not.allocated(	psnxy     	))		allocate(	psnxy     	(ncol)    )
        if(.not.allocated(	savxy      	))		allocate(	savxy      	(ncol)    )
        if(.not.allocated(	sagxy     	))		allocate(	sagxy     	(ncol)    )
        if(.not.allocated(	rssunxy   	))		allocate(	rssunxy   	(ncol)    )
        if(.not.allocated(	rsshaxy 	))		allocate(	rsshaxy   	(ncol)    )
        if(.not.allocated(	bgapxy   	))		allocate(	bgapxy   	(ncol)    )
        if(.not.allocated(	wgapxy    	))		allocate(	wgapxy    	(ncol)    )
        if(.not.allocated(	tgvxy      	))		allocate(	tgvxy      	(ncol)    )
        if(.not.allocated(	tgbxy     	))		allocate(	tgbxy     	(ncol)    )
        if(.not.allocated(	chvxy     	))		allocate(	chvxy     	(ncol)    )
        if(.not.allocated(	chbxy 	    ))		allocate(	chbxy 	    (ncol)    )
        if(.not.allocated(	shgxy    	))		allocate(	shgxy    	(ncol)    )
        if(.not.allocated(	shcxy     	))		allocate(	shcxy     	(ncol)    )
        if(.not.allocated(	shbxy      	))		allocate(	shbxy      	(ncol)    )
        if(.not.allocated(	evgxy     	))		allocate(	evgxy     	(ncol)    )
        if(.not.allocated(	evbxy     	))		allocate(	evbxy     	(ncol)    )
        if(.not.allocated(	ghvxy 	    ))		allocate(	ghvxy 	    (ncol)    )
        if(.not.allocated(	ghbxy    	))		allocate(	ghbxy    	(ncol)    )
        if(.not.allocated(	irgxy     	))		allocate(	irgxy     	(ncol)    )
        if(.not.allocated(	ircxy      	))		allocate(	ircxy      	(ncol)    )
        if(.not.allocated(	irbxy     	))		allocate(	irbxy     	(ncol)    )
        if(.not.allocated(	trxy      	))		allocate(	trxy      	(ncol)    )
        if(.not.allocated(	evcxy 	    ))		allocate(	evcxy 	    (ncol)    )
        if(.not.allocated(	chleafxy 	))		allocate(	chleafxy 	(ncol)    )
        if(.not.allocated(	chucxy   	))		allocate(	chucxy   	(ncol)    )
        if(.not.allocated(	chv2xy     	))		allocate(	chv2xy     	(ncol)    )
        if(.not.allocated(	chb2xy   	))		allocate(	chb2xy   	(ncol)    )
        if(.not.allocated(	tahxy	    ))		allocate(	tahxy	    (ncol)    )	
        if(.not.allocated(	cmxy	    ))		allocate(	cmxy        (ncol)    )	
        if(.not.allocated(	chxy	    ))		allocate(	chxy	    (ncol)    )	
        if(.not.allocated(	fwetxy	    ))		allocate(	fwetxy	    (ncol)    )	
        if(.not.allocated(	sneqvoxy	))		allocate(	sneqvoxy	(ncol)    )	
        if(.not.allocated(	alboldxy  	))		allocate(	alboldxy  	(ncol)    )	

        !!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!                      
        
    
        if(.not.allocated(	dx   	    ))		allocate(	dx   	    (ncol)    )
        if(.not.allocated(	vegfra 	    ))		allocate(	vegfra 	    (ncol)    )
        if(.not.allocated(	coszin 	    ))		allocate(	coszin 	    (ncol)    )
        if(.not.allocated(	xlatin      ))		allocate(	xlatin      (ncol)    )
        if(.not.allocated(	xlongin      ))		allocate(	xlongin      (ncol)    )
        if(.not.allocated(	vegmax  	))		allocate(	vegmax  	(ncol)    )
        if(.not.allocated(	tmn 		))		allocate(	tmn   	    (ncol)    )
        if(.not.allocated(	xland  	    ))		allocate(	xland  	    (ncol)    )
        if(.not.allocated(	xice    	))		allocate(	xice    	(ncol)    )
        if(.not.allocated(	swdown      ))		allocate(	swdown      (ncol)    )
        if(.not.allocated(	glw 	    ))		allocate(	glw 	    (ncol)    )
        if(.not.allocated(	precip_in   ))		allocate(	precip_in   (ncol)    )
        if(.not.allocated(	sr  	    ))		allocate(	sr    	    (ncol)    )
        if(.not.allocated(	tsk    	    ))		allocate(	tsk    	    (ncol)    )
        if(.not.allocated(	hfx         ))		allocate(	hfx         (ncol)    )
        if(.not.allocated(	qfx     	))		allocate(	qfx     	(ncol)    )
        if(.not.allocated(	lh   	    ))		allocate(	lh   	    (ncol)    )
        if(.not.allocated(	grdflx      ))		allocate(	grdflx      (ncol)    )
        if(.not.allocated(	smstav  	))		allocate(	smstav  	(ncol)    )
        if(.not.allocated(	smstot	    ))		allocate(	smstot	    (ncol)    )
        if(.not.allocated(	sfcrunoff   ))		allocate(	sfcrunoff   (ncol)    )
        if(.not.allocated(	udrunoff    ))		allocate(	udrunoff    (ncol)    )
        if(.not.allocated(	albedo      ))		allocate(	albedo      (ncol)    )
        if(.not.allocated(	snowc   	))		allocate(	snowc   	(ncol)    )
        if(.not.allocated(	snow   	    ))		allocate(	snow   	    (ncol)    )
        if(.not.allocated(	snowh   	))		allocate(	snowh   	(ncol)    )
        if(.not.allocated(	canwat      ))		allocate(	canwat      (ncol)    )
        if(.not.allocated(	acsnom 	    ))		allocate(	acsnom   	(ncol)    )
        if(.not.allocated(	acsnow  	))		allocate(	acsnow  	(ncol)    )
        if(.not.allocated(	emiss  	    ))		allocate(	emiss  	    (ncol)    )
        if(.not.allocated(	qsfc 	    ))		allocate(	qsfc 	    (ncol)    )
        if(.not.allocated(	z0   	    ))		allocate(	z0   	    (ncol)    )
        if(.not.allocated(	znt	        ))		allocate(	znt	        (ncol)    )
          
        ! just for init ??
        if(.not.allocated(	lai	        ))		allocate(	lai	        (ncol)    )
        if(.not.allocated(	chstarxy	))		allocate(	chstarxy    (ncol)    )

        ! !!!!!!!!!!!!!!! not needed currently
        if(.not.allocated(	mp_rainc 	))		allocate(	mp_rainc 	(ncol)    )
        if(.not.allocated(	mp_rainnc 	))		allocate(	mp_rainnc 	(ncol)    )
        if(.not.allocated(	mp_shcv 	))		allocate(	mp_shcv 	(ncol)    )
        if(.not.allocated(	mp_snow	    ))		allocate(	mp_snow	    (ncol)    )
        if(.not.allocated(	mp_graup 	))		allocate(	mp_graup 	(ncol)    )
        if(.not.allocated(	mp_hail   	))		allocate(	mp_hail   	(ncol)    )

        !!!!!!integer
        if(.not.allocated(	ivgtyp 	    ))		allocate(	ivgtyp 	    (ncol)    )
        if(.not.allocated(	isltyp  	))		allocate(	isltyp  	(ncol)    )
        if(.not.allocated(	isnowxy   	))		allocate(	isnowxy     (ncol)    )



        ! initialization--5/7
        dx   	 =0.0;	xice    	 =0.0;	lh   	 =0.0;	snowh   	 =0.0;	lai	 =0.0;	mp_rainc 	 =0.0;
        vegfra 	 =0.0;	swdown    	 =0.0;	grdflx      	 =0.0;	canwat     	 =0.0;	chstarxy	 =0.0;	mp_rainnc 	 =0.0;
        coszin 	 =0.0;	glw 	 =0.0;	smstav  	 =0.0;	acsnom 	 =0.0;	ivgtyp 	 =0;	mp_shcv 	 =0.0;
        xlatin     	 =0.0;	precip_in  	 =0.0;	smstot	 =0.0;	acsnow  	 =0.0;	isltyp  	 =0;	mp_snow	 =0.0;
        xlongin     	 =0.0;	sr  	 =0.0;	sfcrunoff  	 =0.0;	emiss  	 =0.0;	isnowxy   	 =0;	mp_graup 	 =0.0;
        vegmax  	 =0.0;	tsk    	 =0.0;	udrunoff  	 =0.0;	qsfc 	 =0.0;			mp_hail   	 =0.0;
        tmn 	 =0.0;	hfx     	 =0.0;	albedo  	 =0.0;	z0   	 =0.0;				
        xland  	 =0;	qfx     	 =0.0;	snowc   	 =0.0;	znt	 =0.0;	snow   	 =0.0;

        taux2d 	 =0.0;	albd2d 	 =0.0;	dzs 	 =0.0;	dz8w  	 =0.0;	lfmassxy  	 =0.0;	deeprechxy 	 =0.0;
        tauy2d 	 =0.0;	albi2d 	 =0.0;	smois 	 =0.0;	t3d  	 =0.0;	rtmassxy  	 =0.0;	rechxy    	 =0.0;
        fire2d 	 =0.0;	greenfrac  	 =0.0;	sh2o  	 =0.0;	qv3d  	 =0.0;	stmassxy 	 =0.0;	grainxy   	 =0.0;
        isnowxy  	 =0;	qsnowxy   	 =0.0;	tslb    	 =0.0;	u_phy   	 =0.0;	woodxy   	 =0.0;	gddxy  	 =0.0;
        tvxy      	 =0.0;	wslakexy  	 =0.0;	tsnoxy 	 =0.0;	v_phy   	 =0.0;	stblcpxy  	 =0.0;	t2mvxy   	 =0.0;
        tgxy       	 =0.0;	zwtxy      	 =0.0;  p8w3d	 =0.0;	fastcpxy   	 =0.0;	t2mbxy    	 =0.0;
        canicexy  	 =0.0;	waxy      	 =0.0;	snicexy   	 =0.0;	fvegxy    	 =0.0;	xlaixy    	 =0.0;	q2mvxy     	 =0.0;
        canliqxy  	 =0.0;	wtxy      	 =0.0;	snliqxy    	 =0.0;	runsfxy   	 =0.0;	xsaixy    	 =0.0;	q2mbxy	 =0.0;
        eahxy 	 =0.0;	shbxy      	 =0.0;	smoiseq  	 =0.0;	runsbxy  	 =0.0;	taussxy 	 =0.0;	tradxy   	 =0.0;
        chucxy   	 =0.0;	evgxy     	 =0.0;	psnxy     	 =0.0;	ecanxy    	 =0.0;	smcwtdxy  	 =0.0;	neexy     	 =0.0;
        chv2xy     	 =0.0;	evbxy     	 =0.0;	savxy      	 =0.0;	edirxy     	 =0.0;	wgapxy    	 =0.0;	gppxy      	 =0.0;
        chb2xy   	 =0.0;	ghvxy 	 =0.0;	sagxy     	 =0.0;	etranxy   	 =0.0;	tgvxy      	 =0.0;	nppxy     	 =0.0;
        tahxy	 =0.0;	ghbxy    	 =0.0;	rssunxy   	 =0.0;	fsaxy     	 =0.0;	tgbxy     	 =0.0;		
        cmxy	 =0.0;	irgxy     	 =0.0;	rsshaxy 	 =0.0;	firaxy 	 =0.0;	chvxy     	 =0.0;		
        chxy	 =0.0;	ircxy      	 =0.0;	bgapxy   	 =0.0;	aparxy   	 =0.0;	chbxy 	 =0.0;		
        fwetxy	 =0.0;	irbxy     	 =0.0;					shgxy    	 =0.0;		
        sneqvoxy =0.0;	trxy      	 =0.0;					shcxy     	 =0.0;		
        alboldxy =0.0;	evcxy 	     =0.0;	chleafxy 	 =0.0;	
        
        zsnsoxy(-2:nsoil,1:ncol)  	 =0; !??????
    end subroutine allocate_lsm_noahmp_1d
     


!=================================================================================================================
    subroutine deallocate_lsm_noahmp_2d()
        !=================================================================================================================
        
         !-----------------------------------------------------------------------------------------------------------------
        
         !arrays for soil layer properties:
         if(allocated(	dzs_2d 	        ))		deallocate(	dzs_2d 	            )
         if(allocated(	smois_2d 	    ))		deallocate(	smois_2d 	        )
         if(allocated(	sh2o_2d  	    ))		deallocate(	sh2o_2d  	        )
         if(allocated(	tslb_2d    	    ))		deallocate(	tslb_2d    	        )
         if(allocated(	tsnoxy_2d 	    ))		deallocate(	tsnoxy_2d 	        )
         if(allocated(	zsnsoxy_2d  	))		deallocate(	zsnsoxy_2d  	    )
         if(allocated(	snicexy_2d   	))		deallocate(	snicexy_2d   	    )
         if(allocated(	snliqxy_2d    	))		deallocate(	snliqxy_2d    	    )
         if(allocated(	smoiseq_2d  	))		deallocate(	smoiseq_2d          )
 
         !arrays for 3D:
         if(allocated(	dz8w_2d  	    ))		deallocate(	dz8w_2d      )
         if(allocated(	t3d_2d  	    ))		deallocate(	t3d_2d       )
         if(allocated(	qv3d_2d  	    ))		deallocate(	qv3d_2d      )
         if(allocated(	u_phy_2d   	    ))		deallocate(	u_phy_2d     )
         if(allocated(	v_phy_2d   	    ))		deallocate(	v_phy_2d     )
         if(allocated(	p8w3d_2d	    ))		deallocate(	p8w3d_2d     )
 
 
 
         !other arrays:
         if(allocated(	isnowxy_2d  	))		deallocate(	isnowxy_2d  	    )
         if(allocated(	tvxy_2d      	))		deallocate(	tvxy_2d      	    )
         if(allocated(	tgxy_2d       	))		deallocate(	tgxy_2d       	    )
         if(allocated(	canicexy_2d  	))		deallocate(	canicexy_2d  	    )
         if(allocated(	canliqxy_2d  	))		deallocate(	canliqxy_2d  	    )
         if(allocated(	eahxy_2d 	    ))		deallocate(	eahxy_2d 	        )
         if(allocated(	qsnowxy_2d   	))		deallocate(	qsnowxy_2d  	    )
         if(allocated(	wslakexy_2d  	))		deallocate(	wslakexy_2d  	    )
         if(allocated(	zwtxy_2d      	))		deallocate(	zwtxy_2d      	    )
         if(allocated(	waxy_2d      	))		deallocate(	waxy_2d      	    )
         if(allocated(	wtxy_2d      	))		deallocate(	wtxy_2d      	    )
         
         if(allocated(	lfmassxy_2d  	))		deallocate(	lfmassxy_2d  	    )
         if(allocated(	rtmassxy_2d  	))		deallocate(	rtmassxy_2d  	    )
         if(allocated(	stmassxy_2d 	))		deallocate(	stmassxy_2d 	    )
         if(allocated(	woodxy_2d   	))		deallocate(	woodxy_2d   	    )
         if(allocated(	stblcpxy_2d  	))		deallocate(	stblcpxy_2d  	    )
         if(allocated(	fastcpxy_2d   	))		deallocate(	fastcpxy_2d   	    )
         if(allocated(	xlaixy_2d    	))		deallocate(	xlaixy_2d    	    )
         if(allocated(	xsaixy_2d    	))		deallocate(	xsaixy_2d    	    )
         if(allocated(	taussxy_2d 	    ))		deallocate(	taussxy_2d 	        )
         if(allocated(	smcwtdxy_2d  	))		deallocate(	smcwtdxy_2d  	    )
         if(allocated(	deeprechxy_2d 	))		deallocate(	deeprechxy_2d 	    )
         if(allocated(	rechxy_2d    	))		deallocate(	rechxy_2d    	    )
         if(allocated(	grainxy_2d   	))		deallocate(	grainxy_2d   	    )
         if(allocated(	gddxy_2d  	    ))		deallocate(	gddxy_2d  	        )
         if(allocated(	t2mvxy_2d   	))		deallocate(	t2mvxy_2d   	    )
         if(allocated(	t2mbxy_2d    	))		deallocate(	t2mbxy_2d    	    )
         if(allocated(	q2mvxy_2d     	))		deallocate(	q2mvxy_2d     	    )
         if(allocated(	q2mbxy_2d       ))		deallocate(	q2mbxy_2d   	    )
         if(allocated(	tradxy_2d   	))		deallocate(	tradxy_2d   	    )
         if(allocated(	neexy_2d     	))		deallocate(	neexy_2d     	    )
         if(allocated(	gppxy_2d      	))		deallocate(	gppxy_2d      	    )
         if(allocated(	nppxy_2d     	))		deallocate(	nppxy_2d     	    )
         if(allocated(	fvegxy_2d    	))		deallocate(	fvegxy_2d    	    )
         if(allocated(	runsfxy_2d   	))		deallocate(	runsfxy_2d 	        )
         if(allocated(	runsbxy_2d  	))		deallocate(	runsbxy_2d  	    )
         if(allocated(	ecanxy_2d    	))		deallocate(	ecanxy_2d    	    )
         if(allocated(	edirxy_2d     	))		deallocate(	edirxy_2d     	    )
         if(allocated(	etranxy_2d   	))		deallocate(	etranxy_2d   	    )
         if(allocated(	fsaxy_2d     	))		deallocate(	fsaxy_2d     	    )
         if(allocated(	firaxy_2d 	    ))		deallocate(	firaxy_2d 	        )
         if(allocated(	aparxy_2d   	))		deallocate(	aparxy_2d   	    )
         if(allocated(	psnxy_2d     	))		deallocate(	psnxy_2d     	    )
         if(allocated(	savxy_2d      	))		deallocate(	savxy_2d      	    )
         if(allocated(	sagxy_2d     	))		deallocate(	sagxy_2d     	    )
         if(allocated(	rssunxy_2d   	))		deallocate(	rssunxy_2d   	    )
         if(allocated(	rsshaxy_2d 	    ))		deallocate(	rsshaxy_2d   	    )
         if(allocated(	bgapxy_2d   	))		deallocate(	bgapxy_2d   	    )
         if(allocated(	wgapxy_2d    	))		deallocate(	wgapxy_2d    	    )
         if(allocated(	tgvxy_2d      	))		deallocate(	tgvxy_2d      	    )
         if(allocated(	tgbxy_2d     	))		deallocate(	tgbxy_2d     	    )
         if(allocated(	chvxy_2d     	))		deallocate(	chvxy_2d     	    )
         if(allocated(	chbxy_2d 	    ))		deallocate(	chbxy_2d 	        )
         if(allocated(	shgxy_2d    	))		deallocate(	shgxy_2d    	    )
         if(allocated(	shcxy_2d     	))		deallocate(	shcxy_2d     	    )
         if(allocated(	shbxy_2d      	))		deallocate(	shbxy_2d      	    )
         if(allocated(	evgxy_2d     	))		deallocate(	evgxy_2d     	    )
         if(allocated(	evbxy_2d     	))		deallocate(	evbxy_2d     	    )
         if(allocated(	ghvxy_2d 	    ))		deallocate(	ghvxy_2d 	        )
         if(allocated(	ghbxy_2d    	))		deallocate(	ghbxy_2d    	    )
         if(allocated(	irgxy_2d     	))		deallocate(	irgxy_2d     	    )
         if(allocated(	ircxy_2d      	))		deallocate(	ircxy_2d      	    )
         if(allocated(	irbxy_2d     	))		deallocate(	irbxy_2d     	    )
         if(allocated(	trxy_2d      	))		deallocate(	trxy_2d      	    )
         if(allocated(	evcxy_2d 	    ))		deallocate(	evcxy_2d 	        )
         if(allocated(	chleafxy_2d 	))		deallocate(	chleafxy_2d 	    )
         if(allocated(	chucxy_2d   	))		deallocate(	chucxy_2d   	    )
         if(allocated(	chv2xy_2d     	))		deallocate(	chv2xy_2d     	    )
         if(allocated(	chb2xy_2d   	))		deallocate(	chb2xy_2d   	    )
         if(allocated(	tahxy_2d	    ))		deallocate(	tahxy_2d	        )	
         if(allocated(	cmxy_2d	        ))		deallocate(	cmxy_2d        	    )	
         if(allocated(	chxy_2d	        ))		deallocate(	chxy_2d	            )	
         if(allocated(	fwetxy_2d	    ))		deallocate(	fwetxy_2d	        )	
         if(allocated(	sneqvoxy_2d	    ))		deallocate(	sneqvoxy_2d	        )	
         if(allocated(	alboldxy_2d  	))		deallocate(	alboldxy_2d  	    )	
     
         !!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!                      
         
        
         if(allocated(	dx_2d   	    ))		deallocate(	dx_2d   	        )
         if(allocated(	vegfra_2d 	    ))		deallocate(	vegfra_2d 	        )
         if(allocated(	coszin_2d 	    ))		deallocate(	coszin_2d 	        )
         if(allocated(	xlatin_2d       ))		deallocate(	xlatin_2d           )
         if(allocated(	xlongin_2d       ))		deallocate(	xlongin_2d           )
         if(allocated(	vegmax_2d  	    ))		deallocate(	vegmax_2d  	        )
         if(allocated(	tmn_2d 		    ))		deallocate(	tmn_2d   	        )
         if(allocated(	xland_2d  	    ))		deallocate(	xland_2d  	        )
         if(allocated(	xice_2d    	    ))		deallocate(	xice_2d    	        )
         if(allocated(	swdown_2d       ))		deallocate(	swdown_2d           )
         if(allocated(	glw_2d 	        ))		deallocate(	glw_2d 	            )
         if(allocated(	precip_in_2d    ))		deallocate(	precip_in_2d        )
         if(allocated(	sr_2d  	        ))		deallocate(	sr_2d    	        )
         if(allocated(	tsk_2d    	    ))		deallocate(	tsk_2d    	        )
         if(allocated(	hfx_2d          ))		deallocate(	hfx_2d              )
         if(allocated(	qfx_2d     	    ))		deallocate(	qfx_2d     	        )
         if(allocated(	lh_2d   	    ))		deallocate(	lh_2d   	        )
         if(allocated(	grdflx_2d       ))		deallocate(	grdflx_2d           )
         if(allocated(	smstav_2d  	    ))		deallocate(	smstav_2d  	        )
         if(allocated(	smstot_2d	    ))		deallocate(	smstot_2d	        )
         if(allocated(	sfcrunoff_2d    ))		deallocate(	sfcrunoff_2d        )
         if(allocated(	udrunoff_2d     ))		deallocate(	udrunoff_2d         )
         if(allocated(	albedo_2d       ))		deallocate(	albedo_2d           )
         if(allocated(	snowc_2d   	    ))		deallocate(	snowc_2d   	        )
         if(allocated(	snow_2d   	    ))		deallocate(	snow_2d   	        )
         if(allocated(	snowh_2d   	    ))		deallocate(	snowh_2d   	        )
         if(allocated(	canwat_2d       ))		deallocate(	canwat_2d           )
         if(allocated(	acsnom_2d 	    ))		deallocate(	acsnom_2d   	    )
         if(allocated(	acsnow_2d  	    ))		deallocate(	acsnow_2d  	        )
         if(allocated(	emiss_2d  	    ))		deallocate(	emiss_2d  	        )
         if(allocated(	qsfc_2d 	    ))		deallocate(	qsfc_2d 	        )
         if(allocated(	z0_2d   	    ))		deallocate(	z0_2d   	        )
         if(allocated(	znt_2d	        ))		deallocate(	znt_2d	            )
                 
         !!!!!!!!!!!!!!!
         if(allocated(	mp_rainc_2d 	))		deallocate(	mp_rainc_2d 	    )
         if(allocated(	mp_rainnc_2d 	))		deallocate(	mp_rainnc_2d 	    )
         if(allocated(	mp_shcv_2d 	    ))		deallocate(	mp_shcv_2d 	        )
         if(allocated(	mp_snow_2d	    ))		deallocate(	mp_snow_2d	        )
         if(allocated(	mp_graup_2d 	))		deallocate(	mp_graup_2d 	    )
         if(allocated(	mp_hail_2d   	))		deallocate(	mp_hail_2d   	    )
 
         !!!!!!integer
         if(allocated(	ivgtyp_2d 	    ))		deallocate(	ivgtyp_2d 	        )
         if(allocated(	isltyp_2d  	    ))		deallocate(	isltyp_2d  	        )
         if(allocated(	isnowxy_2d   	))		deallocate(	isnowxy_2d 	        )
 
        ! !!!!
 

     end subroutine deallocate_lsm_noahmp_2d
    !=================================================================================================================


    !=================================================================================================================
     subroutine deallocate_lsm_noahmp_1d()
    !=================================================================================================================
        
         !arrays for soil layer properties:
         if(allocated(	dzs 	    ))		deallocate(	dzs 	        )
         if(allocated(	smois 	    ))		deallocate(	smois 	        )
         if(allocated(	sh2o  	    ))		deallocate(	sh2o  	        )
         if(allocated(	tslb    	))		deallocate(	tslb    	    )
         if(allocated(	tsnoxy 	    ))		deallocate(	tsnoxy 	        )
         if(allocated(	zsnsoxy  	))		deallocate(	zsnsoxy  	    )
         if(allocated(	snicexy   	))		deallocate(	snicexy   	    )
         if(allocated(	snliqxy    	))		deallocate(	snliqxy    	    )
         if(allocated(	smoiseq  	))		deallocate(	smoiseq         )
       
        !  !arrays for 3D:
         if(allocated(	dz8w  	    ))		deallocate(	dz8w      )
         if(allocated(	t3d  	    ))		deallocate(	t3d       )
         if(allocated(	qv3d  	    ))		deallocate(	qv3d      )
         if(allocated(	u_phy   	))		deallocate(	u_phy     )
         if(allocated(	v_phy   	))		deallocate(	v_phy     )
         if(allocated(	p8w3d	    ))		deallocate(	p8w3d     )
 
          
 
         !other arrays:
         if(allocated(	isnowxy  	))		deallocate(	isnowxy  	    )
         if(allocated(	tvxy      	))		deallocate(	tvxy      	    )
         if(allocated(	tgxy       	))		deallocate(	tgxy       	    )
         if(allocated(	canicexy  	))		deallocate(	canicexy  	    )
         if(allocated(	canliqxy  	))		deallocate(	canliqxy  	    )
         if(allocated(	eahxy 	    ))		deallocate(	eahxy 	        )
         if(allocated(	qsnowxy   	))		deallocate(	qsnowxy  	    )
         if(allocated(	wslakexy  	))		deallocate(	wslakexy  	    )
         if(allocated(	zwtxy      	))		deallocate(	zwtxy      	    )
         if(allocated(	waxy      	))		deallocate(	waxy      	    )
         if(allocated(	wtxy      	))		deallocate(	wtxy      	    )
         
         if(allocated(	lfmassxy  	))		deallocate(	lfmassxy  	    )
         if(allocated(	rtmassxy  	))		deallocate(	rtmassxy  	    )
         if(allocated(	stmassxy 	))		deallocate(	stmassxy 	    )
         if(allocated(	woodxy   	))		deallocate(	woodxy   	    )
         if(allocated(	stblcpxy  	))		deallocate(	stblcpxy  	    )
         if(allocated(	fastcpxy   	))		deallocate(	fastcpxy   	    )
         if(allocated(	xlaixy    	))		deallocate(	xlaixy    	    )
         if(allocated(	xsaixy    	))		deallocate(	xsaixy    	    )
         if(allocated(	taussxy 	))		deallocate(	taussxy 	    )
         if(allocated(	smcwtdxy  	))		deallocate(	smcwtdxy  	    )
         if(allocated(	deeprechxy 	))		deallocate(	deeprechxy 	    )
         if(allocated(	rechxy    	))		deallocate(	rechxy    	    )
         if(allocated(	grainxy   	))		deallocate(	grainxy   	    )
         if(allocated(	gddxy  	    ))		deallocate(	gddxy  	        )
         if(allocated(	t2mvxy   	))		deallocate(	t2mvxy   	    )
         if(allocated(	t2mbxy    	))		deallocate(	t2mbxy    	    )
         if(allocated(	q2mvxy     	))		deallocate(	q2mvxy     	    )
         if(allocated(	q2mbxy      ))		deallocate(	q2mbxy   	    )
         if(allocated(	tradxy   	))		deallocate(	tradxy   	    )
         if(allocated(	neexy     	))		deallocate(	neexy     	    )
         if(allocated(	gppxy      	))		deallocate(	gppxy      	    )
         if(allocated(	nppxy     	))		deallocate(	nppxy     	    )
         if(allocated(	fvegxy    	))		deallocate(	fvegxy    	    )
         if(allocated(	runsfxy   	))		deallocate(	runsfxy 	    )
         if(allocated(	runsbxy  	))		deallocate(	runsbxy  	    )
         if(allocated(	ecanxy    	))		deallocate(	ecanxy    	    )
         if(allocated(	edirxy     	))		deallocate(	edirxy     	    )
         if(allocated(	etranxy   	))		deallocate(	etranxy   	    )
         if(allocated(	fsaxy     	))		deallocate(	fsaxy     	    )
         if(allocated(	firaxy 	    ))		deallocate(	firaxy 	        )
         if(allocated(	aparxy   	))		deallocate(	aparxy   	    )
         if(allocated(	psnxy     	))		deallocate(	psnxy     	    )
         if(allocated(	savxy      	))		deallocate(	savxy      	    )
         if(allocated(	sagxy     	))		deallocate(	sagxy     	    )
         if(allocated(	rssunxy   	))		deallocate(	rssunxy   	    )
         if(allocated(	rsshaxy 	))		deallocate(	rsshaxy   	    )
         if(allocated(	bgapxy   	))		deallocate(	bgapxy   	    )
         if(allocated(	wgapxy    	))		deallocate(	wgapxy    	    )
         if(allocated(	tgvxy      	))		deallocate(	tgvxy      	    )
         if(allocated(	tgbxy     	))		deallocate(	tgbxy     	    )
         if(allocated(	chvxy     	))		deallocate(	chvxy     	    )
         if(allocated(	chbxy 	    ))		deallocate(	chbxy 	        )
         if(allocated(	shgxy    	))		deallocate(	shgxy    	    )
         if(allocated(	shcxy     	))		deallocate(	shcxy     	    )
         if(allocated(	shbxy      	))		deallocate(	shbxy      	    )
         if(allocated(	evgxy     	))		deallocate(	evgxy     	    )
         if(allocated(	evbxy     	))		deallocate(	evbxy     	    )
         if(allocated(	ghvxy 	    ))		deallocate(	ghvxy 	        )
         if(allocated(	ghbxy    	))		deallocate(	ghbxy    	    )
         if(allocated(	irgxy     	))		deallocate(	irgxy     	    )
         if(allocated(	ircxy      	))		deallocate(	ircxy      	    )
         if(allocated(	irbxy     	))		deallocate(	irbxy     	    )
         if(allocated(	trxy      	))		deallocate(	trxy      	    )
         if(allocated(	evcxy 	    ))		deallocate(	evcxy 	        )
         if(allocated(	chleafxy 	))		deallocate(	chleafxy 	    )
         if(allocated(	chucxy   	))		deallocate(	chucxy   	    )
         if(allocated(	chv2xy     	))		deallocate(	chv2xy     	    )
         if(allocated(	chb2xy   	))		deallocate(	chb2xy   	    )
         if(allocated(	tahxy	    ))		deallocate(	tahxy	        )	
         if(allocated(	cmxy	    ))		deallocate(	cmxy        	)	
         if(allocated(	chxy	    ))		deallocate(	chxy	        )	
         if(allocated(	fwetxy	    ))		deallocate(	fwetxy	        )	
         if(allocated(	sneqvoxy	))		deallocate(	sneqvoxy	    )	
         if(allocated(	alboldxy  	))		deallocate(	alboldxy  	    )	
     
         !!!!!!!!!!!!!!!!!!!!!! !!!!!!!!!!!!!!!!!!                      
         
         if(allocated(	dx   	    ))		deallocate(	dx   	        )
         if(allocated(	vegfra 	    ))		deallocate(	vegfra 	        )
         if(allocated(	coszin 	    ))		deallocate(	coszin 	        )
         if(allocated(	xlatin      ))		deallocate(	xlatin          )
         if(allocated(	xlongin      ))		deallocate(	xlongin          )
         if(allocated(	vegmax  	))		deallocate(	vegmax  	    )
         if(allocated(	tmn 		))		deallocate(	tmn   	        )
         if(allocated(	xland  	    ))		deallocate(	xland  	        )
         if(allocated(	xice    	))		deallocate(	xice    	    )
         if(allocated(	swdown      ))		deallocate(	swdown          )
         if(allocated(	glw 	    ))		deallocate(	glw 	        )
         if(allocated(	precip_in   ))		deallocate(	precip_in       )
         if(allocated(	sr  	    ))		deallocate(	sr    	        )
         if(allocated(	tsk    	    ))		deallocate(	tsk    	        )
         if(allocated(	hfx         ))		deallocate(	hfx             )
         if(allocated(	qfx     	))		deallocate(	qfx     	    )
         if(allocated(	lh   	    ))		deallocate(	lh   	        )
         if(allocated(	grdflx      ))		deallocate(	grdflx          )
         if(allocated(	smstav  	))		deallocate(	smstav  	    )
         if(allocated(	smstot	    ))		deallocate(	smstot	        )
         if(allocated(	sfcrunoff   ))		deallocate(	sfcrunoff       )
         if(allocated(	udrunoff    ))		deallocate(	udrunoff        )
         if(allocated(	albedo      ))		deallocate(	albedo          )
         if(allocated(	snowc   	))		deallocate(	snowc   	    )
         if(allocated(	snow   	    ))		deallocate(	snow   	        )
         if(allocated(	snowh   	))		deallocate(	snowh   	    )
         if(allocated(	canwat      ))		deallocate(	canwat          )
         if(allocated(	acsnom 	    ))		deallocate(	acsnom   	    )
         if(allocated(	acsnow  	))		deallocate(	acsnow  	    )
         if(allocated(	emiss  	    ))		deallocate(	emiss  	        )
         if(allocated(	qsfc 	    ))		deallocate(	qsfc 	        )
         if(allocated(	z0   	    ))		deallocate(	z0   	        )
         if(allocated(	znt	        ))		deallocate(	znt	            )
                 
         !!!!!!!!!!!!!!!
         if(allocated(	mp_rainc 	))		deallocate(	mp_rainc 	    )
         if(allocated(	mp_rainnc 	))		deallocate(	mp_rainnc 	    )
         if(allocated(	mp_shcv 	))		deallocate(	mp_shcv 	    )
         if(allocated(	mp_snow	    ))		deallocate(	mp_snow	        )
         if(allocated(	mp_graup 	))		deallocate(	mp_graup 	    )
         if(allocated(	mp_hail   	))		deallocate(	mp_hail   	    )
 
         !!!!!!integer
         if(allocated(	ivgtyp 	    ))		deallocate(	ivgtyp 	        )
         if(allocated(	isltyp  	))		deallocate(	isltyp  	    )
         if(allocated(	isnowxy   	))		deallocate(	isnowxy 	    )
 
        
     end subroutine deallocate_lsm_noahmp_1d



end module grist_lsm_noahmp_vars
